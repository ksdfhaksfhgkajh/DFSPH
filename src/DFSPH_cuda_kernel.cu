#include "DFSPH_cuda_api.cuh"

__device__ double poly6(double distance, double radius)
{
    if (distance < DFSPHSolver::epsilon)
    {
        distance = DFSPHSolver::epsilon;
    }
    if (distance >= radius)
    {
        return 0.0;
    }
    else
    {
        double kernel = 1.0 - pow(distance, 2) / pow(radius, 2);
        return 315.0 / (64.0 * DFSPHSolver::PI * pow(radius, 3)) * pow(kernel, 3);
    }
}

__device__ Vector3D spiky_first_derivative(const Vector3D& diff, double radius)
{
    double distance = sqrt(diff.dot(diff));
    if (distance < DFSPHSolver::epsilon)
    {
        distance = DFSPHSolver::epsilon;
    }
    if (distance >= radius)
    {
        return { 0.0, 0.0, 0.0 };
    }
    else
    {
        Vector3D kernel = 1.0 - (diff / radius);
        return (diff * (- 45.0 / (DFSPHSolver::PI * pow(radius, 4) * distance)))
               * (kernel.dot(kernel)) ;
    }
}

__device__ double viscosity_laplacian(double distance, double radius)
{
    if (distance < DFSPHSolver::epsilon)
    {
        distance = DFSPHSolver::epsilon;
    }
    if (distance >= radius)
    {
        return 0.0;
    }
    else
    {
        return (45.0 / (DFSPHSolver::PI * pow(radius, 6))) * (radius - distance);
    }
}

__global__ void compute_density_kernel(Particle* particles,
                                       const int* neighbor_indices,
                                       const int* neighbor_counts,
                                       int num_particles,
                                       double particle_mass,
                                       double density0,
                                       double radius,
                                       int max_neighbors)
{
    size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];
    if (p_i.type == Particle::Boundary)
    {
        p_i.density = density0;
        return;
    }
    p_i.density = 0.0;

    size_t neighbor_start = idx * max_neighbors;
    int neighbor_count = neighbor_counts[idx];

    for (int n = 0; n < neighbor_count; ++n)
    {
        int j = neighbor_indices[neighbor_start + n];
        if (j == -1) continue;
        Particle& p_j = particles[j];
        if (p_j.type == Particle::Boundary)
        {
            continue;
        }
        double distance = (p_i.position - p_j.position).length();
        p_i.density += particle_mass * poly6(distance, radius);
    }

    if (p_i.density < DFSPHSolver::epsilon)
    {
        p_i.density = DFSPHSolver::epsilon;
    }
}

__global__ void compute_factor_kernel(Particle* particles,
                                       const int* neighbor_indices,
                                       const int* neighbor_counts,
                                       int num_particles,
                                       double particle_mass,
                                       double radius,
                                       int max_neighbors)
{
    size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];
    size_t neighbor_start = idx * max_neighbors;
    int neighbor_count = neighbor_counts[idx];

    Vector3D grad_p_i(0.0, 0.0, 0.0);
    double sum_grad_p_j2 = 0.0;

    for (int n = 0; n < neighbor_count; ++n)
    {
        int j = neighbor_indices[neighbor_start + n];
        if (j == -1) continue;
        Particle &p_j = particles[j];
        if (p_j.type == Particle::Boundary)
        {
            continue;
        }
        Vector3D D_Wij = spiky_first_derivative(p_i.position - p_j.position, radius);
        grad_p_i += D_Wij * particle_mass;
        sum_grad_p_j2 += pow(particle_mass, 2) * (D_Wij.dot(D_Wij));
    }

    double aii = grad_p_i.dot(grad_p_i) + sum_grad_p_j2;
    p_i.factor = 1.0 / (aii + DFSPHSolver::epsilon);

    if (p_i.factor < DFSPHSolver::epsilon)
    {
        p_i.factor = DFSPHSolver::epsilon;
    }
}

__global__ void no_pressure_predict_kernel(Particle* particles,
                                           const int* neighbor_indices,
                                           const int* neighbor_counts,
                                           int num_particles,
                                           double timestep,
                                           double viscosity,
                                           double particle_mass,
                                           Vector3D gravity,
                                           double radius,
                                           int max_neighbors)
{
    int idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];
    if (p_i.type == Particle::Boundary)
    {
        return;
    }
    Vector3D acceleration = gravity;

    int neighbor_start = idx * max_neighbors;
    int neighbor_count = neighbor_counts[idx];

    for (int n = 0; n < neighbor_count; ++n)
    {
        int j = neighbor_indices[neighbor_start + n];
        if (j == -1) continue;
        Particle& p_j = particles[j];

        Vector3D r_ij = p_i.position - p_j.position;
        double distance = r_ij.length();

        if (p_j.type == Particle::Fluid)
        {
            // Viscosity force
            Vector3D v_ij = p_j.velocity - p_i.velocity;
            double W_ij = viscosity_laplacian(distance, radius);
            double avg_density = 0.5 * (p_i.density + p_j.density);
            if (p_j.density > DFSPHSolver::epsilon)
            {
                acceleration += (v_ij / avg_density) * W_ij * viscosity * particle_mass;
            }
        }
    }
    p_i.velocity += acceleration * timestep;
}

__global__ void compute_kappa_kernel(Particle* particles,
                                     const int* neighbor_indices,
                                     const int* neighbor_counts,
                                     int num_particles,
                                     double timestep,
                                     double particle_mass,
                                     double density0,
                                     double radius,
                                     int max_neighbors,
                                     Kappa_t kappa_t)
{
    int idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];

    int neighbor_start = idx * max_neighbors;
    int neighbor_count = neighbor_counts[idx];

    double Dpho_Dt = 0.0;

    for (int n = 0; n < neighbor_count; ++n)
    {
        int j = neighbor_indices[neighbor_start + n];
        if (j == -1) continue;
        Particle& p_j = particles[j];
        if (p_j.type == Particle::Boundary)
        {
            continue;
        }
        Vector3D DW_ij = spiky_first_derivative(p_i.position - p_j.position, radius);
        Dpho_Dt += (p_i.velocity - p_j.velocity).dot(DW_ij) *
                particle_mass;
    }

    if (kappa_t == Is_vel)
    {
        p_i.kappa = Dpho_Dt * p_i.factor / timestep;
    }
    else if (kappa_t == Is_pho)
    {
        double pho_star = p_i.density + timestep * Dpho_Dt;
        p_i.kappa = (pho_star - density0) * p_i.factor / pow(timestep,2);
    }
}

__global__ void adapt_velocities_kernel(Particle* particles,
                                        const int* neighbor_indices,
                                        const int* neighbor_counts,
                                        int num_particles,
                                        double timestep,
                                        double particle_mass,
                                        double radius,
                                        int max_neighbors,
                                        Kappa_t kappa_t)
{
    int idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];
    if (p_i.type == Particle::Boundary)
    {
        return;
    }

    int neighbor_start = idx * max_neighbors;
    int neighbor_count = neighbor_counts[idx];

    Vector3D correction(0.0, 0.0, 0.0);

    for (int n = 0; n < neighbor_count; ++n)
    {
        int j = neighbor_indices[neighbor_start + n];
        if (j == -1) continue;
        Particle& p_j = particles[j];
        if (p_j.type == Particle::Boundary)
        {
            continue;
        }
        correction += spiky_first_derivative(p_i.position - p_j.position, radius)*
                (p_i.kappa / p_i.density + p_j.kappa / p_j.density);
    }
    p_i.velocity = p_i.velocity - correction * (timestep * particle_mass);

    if (kappa_t == Is_pho)
    {
        p_i.position += p_i.velocity * timestep;    // update position after fullfill pho_star == pho_0
    }
}

__global__ void compute_density_error_kernel(Particle* particles,
                                             int num_particles,
                                             double density0,
                                             double* density_errors)
{
    int idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];

    double densityError = p_i.density - density0;
    density_errors[idx] = fabs(densityError);
}

__global__ void apply_boundary_conditions_kernel(Particle* particles,
                                                 const int* neighbor_indices,
                                                 const int* neighbor_counts,
                                                 int num_particles,
                                                 double particle_mass,
                                                 double particle_radius,
                                                 int max_neighbors)
{
    int idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];
    if (p_i.type == Particle::Boundary)
    {
        return;
    }

    int neighbor_start = idx * max_neighbors;
    int neighbor_count = neighbor_counts[idx];

    for (int n = 0; n < neighbor_count; ++n)
    {
        int j = neighbor_indices[neighbor_start + n];
        if (j == -1) continue;
        Particle &p_j = particles[j];

        if (p_j.type == Particle::Fluid)
        {
            continue;
        }

        Vector3D r_ij = p_i.position - p_j.position;
        double distance = r_ij.length();

        double min_distance = 2 * particle_radius;
        if (distance < min_distance && distance > DFSPHSolver::epsilon)
        {
            Vector3D direction = r_ij / distance;
            double penetration_depth = min_distance - distance;
            p_i.position += direction * penetration_depth;

            // Adjust velocities
            double restitution = 0.6; // Coefficient of restitution
            double vn = p_i.velocity.dot(direction);

            if (vn < 0.0)
            {
                double impulse = -(1.0 + restitution) * vn / (1.0 / particle_mass);
                Vector3D impulse_vector = direction * impulse;

                p_i.velocity += impulse_vector / particle_mass;
            }
        }
    }
}

//__global__ void apply_boundary_conditions_kernel(Particle* particles)
//{
//    size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
//
//    double x_min = -1.25, x_max = 1.25;
//    double y_min = -2.0, y_max = 5.0;
//    double z_min = -1.25, z_max = 1.25;
//    double restitution = 0.7;
//
//
//    if (particles[idx].position.x < x_min) {
//        particles[idx].position.x = x_min;
//        particles[idx].velocity.x *= -restitution;
//    }
//    if (particles[idx].position.x > x_max) {
//        particles[idx].position.x = x_max;
//        particles[idx].velocity.x *= -restitution;
//    }
//
//    if (particles[idx].position.y < y_min) {
//        particles[idx].position.y = y_min;
//        particles[idx].velocity.y *= -restitution;
//    }
//    if (particles[idx].position.y > y_max) {
//        particles[idx].position.y = y_max;
//        particles[idx].velocity.y *= -restitution;
//    }
//
//    if (particles[idx].position.z < z_min) {
//        particles[idx].position.z = z_min;
//        particles[idx].velocity.z *= -restitution;
//    }
//    if (particles[idx].position.z > z_max) {
//        particles[idx].position.z = z_max;
//        particles[idx].velocity.z *= -restitution;
//    }
//}