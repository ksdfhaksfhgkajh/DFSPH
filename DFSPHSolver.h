#ifndef DFSPHSOLVER_H
#define DFSPHSOLVER_H

#include "Particle.h"
#include "NeighborSearch.h"
#include <fstream>
#include <sstream>
#include <memory>
#include <random>

class DFSPHSolver
{
public:
	DFSPHSolver(int particle_num, int  region_length);
	//DFSPHSolver(const char* file_path, int timestep);

	void predict_density();

	static constexpr double PI = 3.1415926535897931;

	void export_to_ply(const std::string& filename) const;

private:
	std::vector<Particle> _particles;
	double _timestep{ 0.25 };
	double _radius{ 2 };
	double _particle_mass{ 1 };

	std::unique_ptr<NeighborSearch> _neighbor_grid;

	double _poly6(double distance) const;
	double _poly6_first_derivative(double distance) const;
};

#endif
