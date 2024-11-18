#include "DFSPHSolver.h"

DFSPHSolver::DFSPHSolver(int particle_num, int  region_length)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, region_length);

    _particles.reserve(particle_num);
    for (int i = 0; i < particle_num; ++i) 
    {
        Particle particle;
        particle.position = { dis(gen), dis(gen), dis(gen) };
        particle.velocity = { 0.0, 0.0, 0.0 };
        particle.acceleration = { 0.0, 0.0, 0.0 };

        _particles.emplace_back(std::move(particle));
    }

    _neighbor_grid = std::make_unique<NeighborSearch>(_particles, _radius);
}

inline double DFSPHSolver::_poly6(double distance) const
{
    if (distance >= _radius)
    {
        return 0.0;
    }
    else
    {
        double kernel = 1.0 - pow(distance, 2) / pow(_radius, 2);
        return 315.0 / (64.0 * PI * pow(_radius, 3)) * pow(kernel, 3);
    }
}

inline double DFSPHSolver::_poly6_first_derivative(double distance) const
{
    if (distance >= _radius)
    {
        return 0.0;
    }
    else
    {
        double kernel = 1.0 - pow(distance, 2) / pow(_radius, 2);
        return (-945.0 * distance) / (32.0 * PI * pow(_radius, 5)) * pow(kernel, 2);
    }
}

void DFSPHSolver::predict_density()
{
    for (int i = 0; i < _particles.size(); ++i)
    {
        _particles[i].density = 0.0;
        double sum = 0.0;
        double square_sum = 0.0;

        auto neighbors_index_list = _neighbor_grid->find_neighbours(_particles, i);
        for (const int neighbor_i : neighbors_index_list)
        {
            double distance = _particles[i].cal_distance(_particles[neighbor_i]);
            _particles[i].density += _particle_mass * _poly6(distance);

            sum += _particle_mass * _poly6_first_derivative(distance);
            square_sum += pow(_particle_mass * _poly6_first_derivative(distance), 2);
        }

        _particles[i].factor = _particles[i].density / (pow(sum, 2) + square_sum);
    }
}

void DFSPHSolver::export_to_ply(const std::string& filename) const
{
    std::ofstream ofs(filename);
    if (!ofs.is_open()) 
    {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    ofs << "ply\n";
    ofs << "format ascii 1.0\n";
    ofs << "element vertex " << _particles.size() << "\n"; // 点的数量
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "end_header\n";

    for (const auto& particle : _particles) 
    {
        ofs 
            << particle.position[X] << " "
            << particle.position[Y] << " "
            << particle.position[Z] << "\n";
    }

    ofs.close();
}