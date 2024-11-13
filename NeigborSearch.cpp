#include "NeighborSearch.h"

int NeighborSearch::_compute_index(const std::vector<double>& position)
{
	return static_cast<int>(position[0] / grid_size) +
		static_cast<int>(position[1] / grid_size) * 73856093 +
		static_cast<int>(position[2] / grid_size) * 19349663;
}

NeighborSearch::NeighborSearch(const std::vector<Particle>& particles, double grid_size)
	:grid_size(grid_size) 
{
	for (int i = 0; i < particles.size(); ++i)
	{
		int index = _compute_index(particles[i].position);
		grid.emplace(index, i);
	}
}

void NeighborSearch::find_neighbours(std::vector<int>& neighbor_index, )
{

}
