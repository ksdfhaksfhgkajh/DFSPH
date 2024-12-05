#include "DFSPHSolver.h"

void DFSPHSolver::allocate_device_memory()
{
    cudaMalloc((void**)&d_particles, _particle_num * sizeof(Particle));

    int* h_neighbor_counts = new int[_particle_num];
    for (int i = 0; i < _particle_num; ++i)
    {
        h_neighbor_counts[i] = static_cast<int>(_neighbor_grid->neighbor_indices[i].size());
        if (h_neighbor_counts[i] > _max_neighbors)
            _max_neighbors = h_neighbor_counts[i];
    }

    int* h_neighbor_indices = new int[_particle_num * _max_neighbors];
    for (int i = 0; i < _particle_num; ++i)
    {
        int count = h_neighbor_counts[i];
        for (int n = 0; n < count; ++n)
        {
            h_neighbor_indices[i * _max_neighbors + n] = _neighbor_grid->neighbor_indices[i][n];
        }
        for (int n = count; n < _max_neighbors; ++n)
        {
            h_neighbor_indices[i * _max_neighbors + n] = -1;
        }
    }

    cudaMalloc((void**)&d_neighbor_indices, _particle_num * _max_neighbors * sizeof(int));
    cudaMalloc((void**)&d_neighbor_counts, _particle_num * sizeof(int));

    cudaMemcpy(d_neighbor_indices, h_neighbor_indices, _particle_num * _max_neighbors * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbor_counts, h_neighbor_counts, _particle_num * sizeof(int), cudaMemcpyHostToDevice);

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
    cudaFree(d_particles);
    cudaFree(d_neighbor_indices);
    cudaFree(d_neighbor_counts);
}

void DFSPHSolver::compute_density_alpha()
{
    int threadsPerBlock = 256;
    int blocksPerGrid = (static_cast<int>(_particle_num) + threadsPerBlock - 1) / threadsPerBlock;

    compute_density_kernel <<<blocksPerGrid, threadsPerBlock >>> (d_particles,
                                                                d_neighbor_indices,
                                                                d_neighbor_counts,
                                                                static_cast<int>(_particle_num),
                                                                _particle_mass,
                                                                _density0,
                                                                _radius,
                                                                _max_neighbors);
    cudaDeviceSynchronize();
    compute_factor_kernel<<<blocksPerGrid, threadsPerBlock >>> (d_particles,
                                                                d_neighbor_indices,
                                                                d_neighbor_counts,
                                                                static_cast<int>(_particle_num),
                                                                _particle_mass,
                                                                _radius,
                                                                _max_neighbors);
    cudaDeviceSynchronize();
}

void DFSPHSolver::no_pressure_predict()
{
    int threadsPerBlock = 256;
    int blocksPerGrid = (static_cast<int>(_particle_num) + threadsPerBlock - 1) / threadsPerBlock;

    no_pressure_predict_kernel <<< blocksPerGrid, threadsPerBlock >>>(d_particles,
                                                                  d_neighbor_indices,
                                                                  d_neighbor_counts,
                                                                  static_cast<int>(_particle_num),
                                                                  _timestep,
                                                                  _viscosity,
                                                                  _particle_mass,
                                                                  _gravity,
                                                                  _radius,
                                                                  _max_neighbors);
    cudaDeviceSynchronize();
}

void DFSPHSolver::adapt_velocities(Kappa_t kappa_t)
{
    int threadsPerBlock = 256;
    int blocksPerGrid = (static_cast<int>(_particle_num) + threadsPerBlock - 1) / threadsPerBlock;
    for (int iter = 0; iter < _vel_ajust_max_iter; ++iter)
    {
        compute_kappa_kernel<<< blocksPerGrid, threadsPerBlock >>>(d_particles,
                                                                   d_neighbor_indices,
                                                                   d_neighbor_counts,
                                                                   static_cast<int>(_particle_num),
                                                                   _timestep,
                                                                   _particle_mass,
                                                                   _density0,
                                                                   _radius,
                                                                   _max_neighbors,
                                                                   kappa_t);
        cudaDeviceSynchronize();
        adapt_velocities_kernel<<< blocksPerGrid, threadsPerBlock >>>(d_particles,
                                                                      d_neighbor_indices,
                                                                      d_neighbor_counts,
                                                                      static_cast<int>(_particle_num),
                                                                      _timestep,
                                                                      _particle_mass,
                                                                      _radius,
                                                                      _max_neighbors,
                                                                      kappa_t);
        cudaDeviceSynchronize();

//        double
//        compute_density_error_kernel<<< blocksPerGrid, threadsPerBlock >>>(d_particles,
//                                                                           num_particles,
//                                                                           _density0,
//                                                     double* density_errors);
    }
}

void DFSPHSolver::apply_boundary_conditions()
{
    int num_particles = static_cast<int>(_particles.size());
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_particles + threadsPerBlock - 1) / threadsPerBlock;
    apply_boundary_conditions_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles,
                                                                    d_neighbor_indices,
                                                                    d_neighbor_counts,
                                                                    static_cast<int>(_particle_num),
                                                                    _particle_mass,
                                                                    _boundary_particle_radius,
                                                                    _fluid_particle_radius,
                                                                    _max_neighbors);
    cudaDeviceSynchronize();
}
