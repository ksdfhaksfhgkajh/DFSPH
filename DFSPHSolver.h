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
	DFSPHSolver(const char* file_path);

	void compute_density() const;
	void no_pressure_predict() const;
	void solve_divergence_free() const;
	void solve_density_constant() const;
	void simulate() const;
	void apply_boundary_conditions() const;

	static constexpr double PI = 3.1415926535897931;

	void export_to_ply(const std::string& filename) const;

private:
	mutable std::vector<Particle> _particles;
	const double _timestep{ 0.001 };
	const int _framenum{ 100 };
	const double _radius{ 2 };
	const double _particle_mass{ 1.0 };
	const double _viscosity{ 0.01 };
	const double _density0{ 0.001 };
	const double _stiffness{ 10.0 };
	const Particle::Vector3D _gravity{ 0.0, -9.8, 0.0 };

	std::unique_ptr<NeighborSearch> _neighbor_grid;

	double _poly6(double distance) const;
	Particle::Vector3D _spiky_first_derivative(Particle::Vector3D distance) const;
	double _viscosity_laplacian(double distance) const;

	void _calcu_a_gra(size_t index) const;
	void _calcu_a_vis(size_t index, const std::vector<int>& neighbors_index_list) const;
	void _calcu_a_pre(size_t index, const std::vector<int>& neighbors_index_list) const;
};

#endif
