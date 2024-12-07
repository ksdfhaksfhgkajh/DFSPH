#include "NeighborSearch.h"

CellCoord NeighborSearch::_compute_index(const Vector3D& position) const
{
    int ix = static_cast<int>(std::floor(position.x / _grid_size));
    int iy = static_cast<int>(std::floor(position.y / _grid_size));
    int iz = static_cast<int>(std::floor(position.z / _grid_size));
    return CellCoord{ix, iy, iz};
}

void NeighborSearch::_get_hash_search_range(std::array<CellCoord, 27>& search_space, const CellCoord& self_cell) const
{
    int index = 0;
    for (int dx = -1; dx <= 1; ++dx)
    {
        for (int dy = -1; dy <= 1; ++dy)
        {
            for (int dz = -1; dz <= 1; ++dz)
            {
                search_space[index++] = CellCoord{
                        self_cell.ix + dx,
                        self_cell.iy + dy,
                        self_cell.iz + dz
                };
            }
        }
    }
}

NeighborSearch::NeighborSearch(const std::vector<Particle>& particles, double grid_size)
        : _grid_size(grid_size)
{
    for (int i = 0; i < static_cast<int>(particles.size()); ++i)
    {
        CellCoord cell = _compute_index(particles[i].position);
        _grid[cell].push_back(i);
    }

    neighbor_indices.resize(particles.size());

    for (int i = 0; i < static_cast<int>(particles.size()); ++i)
    {
        if (particles[i].type == Particle::Boundary)
        {
            continue;
        }

        CellCoord self_cell = _compute_index(particles[i].position);
        std::array<CellCoord, 27> search_space{};
        _get_hash_search_range(search_space, self_cell);

        std::vector<int> neighbors;
        neighbors.reserve(64); // 预留空间可优化性能

        for (auto &sc : search_space)
        {
            auto it = _grid.find(sc);
            if (it != _grid.end())
            {
                for (auto neighbor_i : it->second)
                {
                    if (neighbor_i != i)
                    {
                        double distance = particles[i].cal_distance(particles[neighbor_i]);
                        if (distance <= _grid_size)
                        {
                            neighbors.push_back(neighbor_i);
                        }
                    }
                }
            }
        }

        neighbor_indices[i] = std::move(neighbors);
    }
}
