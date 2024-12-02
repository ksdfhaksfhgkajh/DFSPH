#ifndef DFSPHSOLVER_H
#define DFSPHSOLVER_H

#include "Particle.h"
#include "NeighborSearch.h"
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

	void compute_density();
	void no_pressure_predict();
	void solve_divergence_free();
	void solve_density_constant();
	void simulate();
	void apply_boundary_conditions();

	static constexpr double PI = 3.1415926535897931;
    static constexpr double epsilon = 1e-8;

	void export_to_ply(const std::string& filename) const;


private:
	std::vector<Particle> _particles;
	const double _timestep{ 0.05 };
	const int _framenum{ 100 };
	const double _radius{ 2 };
	const double _particle_mass{ 1.0 };
	const double _viscosity{ 0.1 };
	const double _density0{ 1.0 };
	const double _stiffness{ 0.01 };
	const Particle::Vector3D _gravity{ 0.0, -9.8, 0.0 };
	const char* _output_path = "../result/";

	NeighborSearch* _neighbor_grid{};

	double _poly6(double distance) const;
	Particle::Vector3D _spiky_first_derivative(const Particle::Vector3D& distance) const;
	double _viscosity_laplacian(double distance) const;

	void _calcu_a_gra(size_t index);
	void _calcu_a_vis(size_t index, const std::vector<int>& neighbors_index_list);
	void _calcu_a_pre(size_t index, const std::vector<int>& neighbors_index_list);
};

#endif
