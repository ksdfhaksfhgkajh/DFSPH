#ifndef DFSPHSOLVER_CUH
#define DFSPHSOLVER_CUH

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
	DFSPHSolver(int particle_num, int  region_length);
	explicit DFSPHSolver(const char* file_path);
    DFSPHSolver(const char* input_path, const char* output_path);

	void compute_density();
	void no_pressure_predict();
    void adapt_velocities(Kappa_t kappa_t);
	void simulate();
	void apply_boundary_conditions();

	static constexpr double PI = 3.1415926535897931;
    static constexpr double epsilon = 1e-8;

	void export_to_ply(const std::string& filename) const;


private:
	std::vector<Particle> _particles;
	const double _timestep{ 0.05 };
	const int _framenum{ 200 };
	const double _radius{ 0.4 };
	const double _particle_mass{ 1.0 };
	const double _viscosity{ 1.0e-3 };
	const double _density0{ 1.0 };
	const Vector3D _gravity{0.0, -9.8, 0.0};
    const int _vel_ajust_max_iter{3};
    double _density_error_threshold{};
	const char* _output_path = "..";

    std::unique_ptr<NeighborSearch> _neighbor_grid;
    int _max_neighbors{0};

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
