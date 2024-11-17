#ifndef NEIGHBORSEARCH_H
#define NEIGHBORSEARCH_H

#include <unordered_map>
#include "Particle.h"

class NeighborSearch
{
public:
	NeighborSearch(const std::vector<Particle>& particles, double grid_size);
	std::vector<int> find_neighbours(const std::vector<Particle>& particles, int self_index) const;

private:
	double _grid_size;
	static constexpr int _hash_prime[3] = { 1, 73856093, 19349663 };
	std::unordered_map<int, std::vector<int>> _grid;

	int _compute_index(const std::vector<double>& position) const;
	void _get_hash_search_range(int(&search_space)[27], int self_hash_index) const;
};

#endif