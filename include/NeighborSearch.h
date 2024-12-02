#ifndef NEIGHBORSEARCH_H
#define NEIGHBORSEARCH_H

#include "Particle.h"
#include <vector>
#include <unordered_map>

class NeighborSearch
{
public:
    NeighborSearch(const std::vector<Particle>& particles, double grid_size);
    std::vector<std::vector<int>> neighbor_indices; // 存储每个粒子的邻居索引

private:
    const double _grid_size;
    static constexpr int _hash_prime[3] = { 1, 73856093, 19349663 };
    std::unordered_map<int, std::vector<int>> _grid;

    int _compute_index(const Particle::Vector3D& position) const;
    void _get_hash_search_range(int(&search_space)[27], int self_hash_index) const;
};

#endif