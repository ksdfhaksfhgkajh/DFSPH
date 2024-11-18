#include "NeighborSearch.h"

int NeighborSearch::_compute_index(const std::vector<double>& position) const
{
	int hash_index = static_cast<int>(position[X] / _grid_size) * _hash_prime[0] +
		static_cast<int>(position[Y] / _grid_size) * _hash_prime[1] +
		static_cast<int>(position[Z] / _grid_size) * _hash_prime[2];

	return hash_index;
}

NeighborSearch::NeighborSearch(const std::vector<Particle>& particles, double grid_size)
	:_grid_size(grid_size)
{
	for (int i = 0; i < particles.size(); ++i)
	{
		int index = _compute_index(particles[i].position);
		_grid[index].push_back(i);
	}
}

void NeighborSearch::_get_hash_search_range(int(&search_space)[27], int self_hash_index) const
{
	int index = 0;
	for (int dx = -1; dx <= 1; ++dx)
	{
		for (int dy = -1; dy <= 1; ++dy)
		{
			for (int dz = -1; dz <= 1; ++dz)
			{
				search_space[index] = self_hash_index +
					dx * _hash_prime[0] +
					dy * _hash_prime[1] +
					dz * _hash_prime[2];
				++index;
			}
		}
	}

}

std::vector<int> NeighborSearch::find_neighbours(const std::vector<Particle>& particles, int self_index) const
{
	std::vector<int> neighbors_index;
	int self_hash_index = _compute_index(particles[self_index].position);
	int search_space[27]{};
	_get_hash_search_range(search_space, self_hash_index);
	for (auto hash_index : search_space)
	{
		if (_grid.find(hash_index) != _grid.end())
		{
			for (auto neighbor_i : _grid.at(hash_index))
			{
				if (particles[neighbor_i].cal_distance(particles[self_index]) <= _grid_size)
				{
					neighbors_index.push_back(neighbor_i);
				}
			}
		}
	}
	return neighbors_index;
}
