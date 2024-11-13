#ifndef NEIGHBORSEARCH_H
#define NEIGHBORSEARCH_H

#include <unordered_map>
#include "Particle.h"

class NeighborSearch
{
private:
	double grid_size;
	std::unordered_map<int, std::vector<int>> grid;

	int _compute_index(const std::vector<double>& position);

public:
	NeighborSearch(const std::vector<Particle>& particles, double grid_size);
	void find_neighbours();
};

#endif