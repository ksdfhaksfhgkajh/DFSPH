#ifndef NEIGHBORSEARCH_H
#define NEIGHBORSEARCH_H

#include "Particle.h"
#include <vector>
#include <unordered_map>
#include <cmath>
#include <array>

struct CellCoord {
    int ix;
    int iy;
    int iz;

    bool operator==(const CellCoord &other) const {
        return ix == other.ix && iy == other.iy && iz == other.iz;
    }
};

struct CellCoordHash {
    std::size_t operator()(const CellCoord &c) const {
        // 简单的hash组合, 可以根据需要改进
        std::size_t h1 = std::hash<int>()(c.ix);
        std::size_t h2 = std::hash<int>()(c.iy);
        std::size_t h3 = std::hash<int>()(c.iz);
        // 经典 hash combine 方法
        return h1 ^ (h2 + 0x9e3779b9 + (h1<<6) + (h1>>2)) ^
               (h3 + 0x9e3779b9 + (h2<<6) + (h2>>2));
    }
};

class NeighborSearch
{
public:
    NeighborSearch(const std::vector<Particle>& particles, double grid_size);
    std::vector<std::vector<int>> neighbor_indices;

private:
    const double _grid_size;
    std::unordered_map<CellCoord, std::vector<int>, CellCoordHash> _grid;

    CellCoord _compute_index(const Vector3D& position) const;
    void _get_hash_search_range(std::array<CellCoord, 27>& search_space, const CellCoord& self_cell) const;
};

#endif
