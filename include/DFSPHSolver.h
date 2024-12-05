#ifndef DFSPHSOLVER_H
#define DFSPHSOLVER_H

#include "Particle.h"
#include "Vector3D.cuh"
#include "NeighborSearch.h"
#include "DFSPH_cuda_api.cuh"
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <random>

class DFSPHSolver
{
public:
    DFSPHSolver(const char* fluid_input_path,
                const char* boundary_input_path,
                const char* output_path);

    void compute_density_alpha();
    void no_pressure_predict();
    void adapt_velocities(Kappa_t kappa_t);
    void simulate();
    void apply_boundary_conditions();

    static constexpr double PI = 3.1415926535897931;
    static constexpr double epsilon = 1e-12;

    void export_to_ply(const std::string& filename) const;


private:
	std::vector<Particle> _particles;
    size_t _particle_num{};
    size_t _fluid_particle_num{};

	const double _timestep{ 0.005 };
	const int _framenum{ 500 };
	const double _radius{ 0.4 };
	const double _particle_mass{ 1.0 };
    const double _boundary_particle_radius{ 0.1 };
    const double _fluid_particle_radius{ 0.005 };
	const double _viscosity{ 8.0e-2 };
	const double _density0{ 400.0 };
	const Vector3D _gravity{0.0, -9.8, 0.0};
    const int _vel_ajust_max_iter{ 4 };
    double _density_error_threshold{ 0.01 * _density0 };
	const char* _output_path = "..";

    std::unique_ptr<NeighborSearch> _neighbor_grid;
    int _max_neighbors{};

    void _load_particles(const char* filepath, Particle::Particle_t particle_type);

public:
    void allocate_device_memory();
    void copy_data_to_device();
    void copy_data_to_host();
    void free_device_memory();

private:
    Particle* d_particles{};
    int* d_neighbor_indices{};
    int* d_neighbor_counts{};
};

#endif
