#include "DFSPH_cuda_api.cuh"

void DFSPHSolver::allocate_device_memory()
{
    size_t num_particles = _particles.size();
    cudaMalloc((void**)&d_particles,
               num_particles * sizeof(Particle));

    int* h_neighbor_counts = new int[num_particles];
    int max_neighbors = 0;
    for (int i = 0; i < num_particles; ++i)
    {
        h_neighbor_counts[i] = static_cast<int>(_neighbor_grid->neighbor_indices[i].size());
        if (h_neighbor_counts[i] > max_neighbors)
            max_neighbors = h_neighbor_counts[i];
    }

    int* h_neighbor_indices = new int[num_particles * max_neighbors];
    for (int i = 0; i < num_particles; ++i)
    {
        int count = h_neighbor_counts[i];
        for (int n = 0; n < count; ++n)
        {
            h_neighbor_indices[i * max_neighbors + n] = _neighbor_grid->neighbor_indices[i][n];
        }
        for (int n = count; n < max_neighbors; ++n)
        {
            h_neighbor_indices[i * max_neighbors + n] = -1;
        }
    }

    cudaMalloc((void**)&d_neighbor_indices,
               num_particles * max_neighbors * sizeof(int));
    cudaMalloc((void**)&d_neighbor_counts,
               num_particles * sizeof(int));

    cudaMemcpy(d_neighbor_indices,
               h_neighbor_indices,
               num_particles * max_neighbors * sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbor_counts,
               h_neighbor_counts,
               num_particles * sizeof(int),
               cudaMemcpyHostToDevice);

    delete[] h_neighbor_indices;
    delete[] h_neighbor_counts;
}

void DFSPHSolver::copy_data_to_device()
{
    size_t num_particles = _particles.size();
    cudaMemcpy(d_particles,
               _particles.data(),
               num_particles * sizeof(Particle),
               cudaMemcpyHostToDevice);
}

void DFSPHSolver::copy_data_to_host()
{
    size_t num_particles = _particles.size();
    cudaMemcpy(_particles.data(),
               d_particles,
               num_particles * sizeof(Particle),
               cudaMemcpyDeviceToHost);
}

void DFSPHSolver::free_device_memory()
{

}

__device__ double poly6(double distance, double radius)
{
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

__global__ void compute_density(Particle* Particles)
{

}

//__global__ void no_pressure_predict(Particle* Particles)
//{
//
//}
//
//__global__ void apply_boundary_conditions(Particle* Particles)
//{
//
//}