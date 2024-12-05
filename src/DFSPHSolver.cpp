#include "DFSPHSolver.h"

DFSPHSolver::DFSPHSolver(const char* fluid_input_path,
                         const char* boundary_input_path,
                         const char* output_path)
        :_output_path(output_path)
{
    _load_particles(fluid_input_path, Particle::Fluid);
    _load_particles(boundary_input_path, Particle::Boundary);
}

void DFSPHSolver::_load_particles(const char* filepath, Particle::Particle_t particle_type)
{
    std::ifstream infile(filepath);
    if (!infile.is_open())
    {
        throw std::runtime_error("Failed to open file: " + std::string(filepath));
    }

    std::string line;
    int vertex_count = 0;
    bool header_ended = false;

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "element" && line.find("vertex") != std::string::npos)
        {
            iss >> token >> vertex_count;
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
        particle.type = particle_type;
        _particles.emplace_back(particle);
    }
    if (particle_type == Particle::Fluid)
    {
        _fluid_particle_num = _particles.size();
    }
    else if(particle_type==Particle::Boundary)
    {
        _particle_num = _particles.size();
    }
    infile.close();
    std::cout << "Loaded " << _particles.size() << " particles from " << filepath << ".\n";
}

void DFSPHSolver::export_to_ply(const std::string& filename) const
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
    {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::cout << "generating " << filename << std::endl;

    ofs << "ply\n";
    ofs << "format ascii 1.0\n";
    ofs << "element vertex " << _fluid_particle_num << "\n";
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "end_header\n";

    for (const auto& particle : _particles)
    {
        if (particle.type == Particle::Boundary)
        {
            continue;
        }
        ofs
                << particle.position.x << " "
                << particle.position.y << " "
                << particle.position.z << "\n";
    }

    ofs.close();
}

void DFSPHSolver::simulate()
{
    _neighbor_grid = std::make_unique<NeighborSearch>(_particles, _radius);
    allocate_device_memory();
    copy_data_to_device();
    compute_density_alpha();

    for (int i = 0; i < _framenum; ++i)
    {
        no_pressure_predict();
        adapt_velocities(Is_pho);
        apply_boundary_conditions();

        copy_data_to_host();
        free_device_memory();

        std::ostringstream filename;
        filename << _output_path << "/frame_" << i + 1 << ".ply";
        export_to_ply(filename.str());

        _neighbor_grid = std::make_unique<NeighborSearch>(_particles, _radius);
        allocate_device_memory();
        copy_data_to_device();
        compute_density_alpha();
        adapt_velocities(Is_vel);
    }
    free_device_memory();
}