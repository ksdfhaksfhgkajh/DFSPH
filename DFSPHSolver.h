#ifndef DFSPHSOLVER_H
#define DFSPHSOLVER_H

#include "Particle.h"
#include "NeighborSearch.h"

class DFSPHSolver
{
private:
	std::vector<Particle> particles;
	double timestep;

	NeighborSearch neighbour_grid;

public:

};

#endif
