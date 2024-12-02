#include "../include/DFSPHSolver.h"

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

DFSPHSolver::DFSPHSolver(const char* file_path)
{
    std::ifstream infile(file_path);
    if (!infile.is_open()) 
    {
        throw std::runtime_error("Failed to open file: " + std::string(file_path));
    }

    std::string line;
    int vertex_count = 0;
    bool header_ended = false;

    // 解析文件头部
    while (std::getline(infile, line)) 
    {
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "element" && line.find("vertex") != std::string::npos) 
        {
            iss >> token >> vertex_count; // 获取顶点数量
        }
        else if (token == "end_header") 
        {
            header_ended = true;
            break;
        }
    }

    if (!header_ended) 
    {
        throw std::runtime_error("PLY file header is invalid or missing 'end_header'.");
    }

    // 读取粒子数据
    for (int i = 0; i < vertex_count; ++i) 
    {
        if (!std::getline(infile, line)) 
        {
            throw std::runtime_error("Unexpected end of file when reading vertex data.");
        }

        std::istringstream iss(line);
        double x, y, z;
        iss >> x >> y >> z;

        Particle particle;
        particle.position = { x, y, z };
        particle.velocity = { 0.0, 0.0, 0.0 };
        particle.acceleration = { 0.0, 0.0, 0.0 };
        // 将粒子位置存储
        _particles.emplace_back(std::move(particle));
    }

    infile.close();
    std::cout << "Loaded " << _particles.size() << " particles from " << file_path << ".\n";
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

Particle::Vector3D DFSPHSolver::_spiky_first_derivative(Particle::Vector3D diff) const
{
    double distance = sqrt(diff * diff);
    if (distance >= _radius)
    {
        return { 0.0, 0.0, 0.0 };
    }
    else
    {
        Particle::Vector3D kernel = 1.0 - (diff / _radius);
        return ((- 45.0 / (PI * pow(_radius, 4) * distance)) * (kernel * kernel) * diff);
    }
}

double DFSPHSolver::_viscosity_laplacian(double distance) const
{
    if (distance >= _radius)
    {
        return 0.0;
    }
    else 
    {
        return (45.0 / (PI * pow(_radius, 6))) * (_radius - distance);
    }
}

void DFSPHSolver::compute_density() const
{
    for (int i = 0; i < _particles.size(); ++i)
    {
        _particles[i].density = 0.0;
        Particle::Vector3D sum(3, 0.0);
        double square_sum = 0.0;

        auto neighbors_index_list = _neighbor_grid->find_neighbours(_particles, i);
        for (const int neighbor_i : neighbors_index_list)
        {
            double distance = _particles[i].cal_distance(_particles[neighbor_i]);
            _particles[i].density += _particle_mass * _poly6(distance);

            Particle::Vector3D D_Wij = _spiky_first_derivative(_particles[i].position - _particles[neighbor_i].position);
            sum += _particle_mass * D_Wij;
            square_sum += pow(_particle_mass, 2) * (D_Wij * D_Wij);
        }

        _particles[i].factor = _particles[i].density / (sum * sum + square_sum);
    }
}

void DFSPHSolver::_calcu_a_gra(size_t index) const
{
    _particles[index].acceleration += _gravity;
}

void DFSPHSolver::_calcu_a_vis(size_t index, const std::vector<int>& neighbors_index_list) const
{
    for (const int neighbor_i : neighbors_index_list) 
    {
        Particle::Vector3D v_ij = _particles[neighbor_i].velocity - _particles[index].velocity;
        double W_ij = _viscosity_laplacian(_particles[index].cal_distance(_particles[neighbor_i]));

        _particles[index].acceleration += W_ij * _viscosity * _particle_mass * (v_ij / _particles[neighbor_i].density);
    }
}

void DFSPHSolver::_calcu_a_pre(size_t index, const std::vector<int>& neighbors_index_list) const
{
    for (const int neighbor_i : neighbors_index_list) 
    {
        double p_i = _stiffness * (_particles[index].density - _density0);
        double p_j = _stiffness * (_particles[neighbor_i].density - _density0);

        Particle::Vector3D D_Wij = _spiky_first_derivative(
            _particles[index].position - _particles[neighbor_i].position
        );

        _particles[index].acceleration += -_particle_mass * (
            (p_i / (_particles[index].density * _particles[index].density)) +
            (p_j / (_particles[neighbor_i].density * _particles[neighbor_i].density))
            ) * D_Wij;
    }
}

void DFSPHSolver::no_pressure_predict() const
{
    for (size_t i = 0; i < _particles.size(); ++i)
    {
        _calcu_a_gra(i);
        auto neighbors_index_list = _neighbor_grid->find_neighbours(_particles, i);
        _calcu_a_vis(i, neighbors_index_list);
        _calcu_a_pre(i, neighbors_index_list);

        _particles[i].velocity += _timestep * _particles[i].acceleration;
        _particles[i].position += _timestep * _particles[i].velocity;
    }
}

//void DFSPHSolver::solve_divergence_free() const
//{
//
//}

void DFSPHSolver::simulate() const
{
    for (int i = 0; i < _framenum; ++i)
    {
        std::ostringstream filename;
        compute_density();
        no_pressure_predict();
        apply_boundary_conditions();

        filename << "frame_" << i + 1 << ".ply";
        export_to_ply(filename.str());
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

void DFSPHSolver::apply_boundary_conditions() const
{
    // 假设边界范围
    double x_min = -2.5, x_max = 5.0;
    double y_min = -2.5, y_max = 5.0;
    double z_min = -2.5, z_max = 5.0;
    double restitution = 0.5; // 反弹系数 (0: 吸收全部能量，1: 完全反弹)

    for (int i = 0; i < _particles.size(); ++i) {
        // 检查 X 方向边界
        if (_particles[i].position[0] < x_min) {
            _particles[i].position[0] = x_min;
            _particles[i].velocity[0] *= -restitution;
        }
        if (_particles[i].position[0] > x_max) {
            _particles[i].position[0] = x_max;
            _particles[i].velocity[0] *= -restitution;
        }

        // 检查 Y 方向边界
        if (_particles[i].position[1] < y_min) {
            _particles[i].position[1] = y_min;
            _particles[i].velocity[1] *= -restitution;
        }
        if (_particles[i].position[1] > y_max) {
            _particles[i].position[1] = y_max;
            _particles[i].velocity[1] *= -restitution;
        }

        // 检查 Z 方向边界
        if (_particles[i].position[2] < z_min) {
            _particles[i].position[2] = z_min;
            _particles[i].velocity[2] *= -restitution;
        }
        if (_particles[i].position[2] > z_max) {
            _particles[i].position[2] = z_max;
            _particles[i].velocity[2] *= -restitution;
        }
    }
}

