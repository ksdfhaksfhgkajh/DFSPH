#ifndef DFSPH_CUDA_API_CUH
#define DFSPH_CUDA_API_CUH

#ifndef __CUDACC__
#define __host__
#define __device__
#define __global__
#endif

#include "DFSPHSolver.cuh"

__device__ double poly6(double distance, double radius);

__device__ Vector3D spiky_first_derivative(const Vector3D& diff, double radius);

__device__ double viscosity_laplacian(double distance, double radius);

__global__ void compute_density_kernel(Particle* particles,
                                       const int* neighbor_indices,
                                       const int* neighbor_counts,
                                       int num_particles,
                                       double particle_mass,
                                       double radius,
                                       int max_neighbors);

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
                                           int max_neighbors);

__global__ void apply_boundary_conditions_kernel(Particle* particles);

#endif