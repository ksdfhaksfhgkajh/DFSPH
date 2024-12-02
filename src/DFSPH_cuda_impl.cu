#include "DFSPH_cuda_api.cuh"

__device__ double poly6(double distance, double radius)
{
    if (distance >= radius)
    {
        return 0.0;
    }
    else if (distance < DFSPHSolver::epsilon)
    {
        return 1.0 / DFSPHSolver::epsilon;
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
    if (distance >= radius)
    {
        return { 0.0, 0.0, 0.0 };
    }
    else if (distance < DFSPHSolver::epsilon)
    {
        return {1.0 / DFSPHSolver::epsilon,
                1.0 / DFSPHSolver::epsilon,
                1.0 / DFSPHSolver::epsilon};
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
    if (distance >= radius)
    {
        return 0.0;
    }
    else if (distance < DFSPHSolver::epsilon)
    {
        return 1.0 / DFSPHSolver::epsilon;
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
                                       double radius,
                                       int max_neighbors)
{
    size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];
    p_i.density = 0.0;

    size_t neighbor_start = idx * max_neighbors;
    int neighbor_count = neighbor_counts[idx];

    for (int n = 0; n < neighbor_count; ++n)
    {
        int j = neighbor_indices[neighbor_start + n];
        if (j == -1) continue;
        Particle& p_j = particles[j];
        double distance = (p_i.position - p_j.position).length();
        p_i.density += particle_mass * poly6(distance, radius);
    }

    if (p_i.density < DFSPHSolver::epsilon)
    {
        p_i.density = DFSPHSolver::epsilon;
    }
}

__global__ void no_pressure_predict_kernel(Particle* particles,
                                           const int* neighbor_indices,
                                           const int* neighbor_counts,
                                           int num_particles,
                                           double timestep,
                                           double viscosity,
                                           double stiffness,
                                           double particle_mass,
                                           double density0,
                                           Vector3D gravity,
                                           double radius,
                                           int max_neighbors)
{
    int idx = static_cast<int>(blockDim.x * blockIdx.x + threadIdx.x);
    if (idx >= num_particles) return;

    Particle& p_i = particles[idx];
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

        // Viscosity force
        Vector3D v_ij = p_j.velocity - p_i.velocity;
        double W_ij = viscosity_laplacian(distance, radius);
        if (p_j.density > DFSPHSolver::epsilon && W_ij > DFSPHSolver::epsilon)
        {
            acceleration += (v_ij / p_j.density) * W_ij * viscosity * particle_mass;
        }

        // Pressure
        double pressure_i = stiffness * (p_i.density - density0);
        double pressure_j = stiffness * (p_j.density - density0);

        if (p_j.density < DFSPHSolver::epsilon)
        {
            continue;
        }

        Vector3D D_Wij = spiky_first_derivative(p_i.position - p_j.position, radius);

        acceleration += D_Wij * (-particle_mass) * (
                (pressure_i / (p_i.density * p_i.density)) +
                (pressure_j / (p_j.density * p_j.density))
        );
    }

    p_i.velocity += acceleration * timestep;
    p_i.position += p_i.velocity * timestep;
}

__global__ void apply_boundary_conditions_kernel(Particle* particles)
{
    size_t idx = blockDim.x * blockIdx.x + threadIdx.x;

    double x_min = -2.5, x_max = 5.0;
    double y_min = -2.5, y_max = 5.0;
    double z_min = -2.5, z_max = 5.0;
    double restitution = 0.5;


    if (particles[idx].position.x < x_min) {
        particles[idx].position.x = x_min;
        particles[idx].velocity.x *= -restitution;
    }
    if (particles[idx].position.x > x_max) {
        particles[idx].position.x = x_max;
        particles[idx].velocity.x *= -restitution;
    }

    if (particles[idx].position.y < y_min) {
        particles[idx].position.y = y_min;
        particles[idx].velocity.y *= -restitution;
    }
    if (particles[idx].position.y > y_max) {
        particles[idx].position.y = y_max;
        particles[idx].velocity.y *= -restitution;
    }

    if (particles[idx].position.z < z_min) {
        particles[idx].position.z = z_min;
        particles[idx].velocity.z *= -restitution;
    }
    if (particles[idx].position.z > z_max) {
        particles[idx].position.z = z_max;
        particles[idx].velocity.z *= -restitution;
    }
}