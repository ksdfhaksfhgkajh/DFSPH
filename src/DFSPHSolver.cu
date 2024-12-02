#include "DFSPHSolver.cuh"

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

        _particles.emplace_back(particle);

        _neighbor_grid = std::make_unique<NeighborSearch>(_particles, _radius);
    }
}

DFSPHSolver::DFSPHSolver(const char* input_path, const char* output_path)
    :_output_path(output_path)
{
    std::ifstream infile(input_path);
    if (!infile.is_open()) 
    {
        throw std::runtime_error("Failed to open file: " + std::string(input_path));
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
        particle.acceleration = { 0.0, 0.0, 0.0 };
        _particles.emplace_back(particle);
    }

    infile.close();
    std::cout << "Loaded " << _particles.size() << " particles from " << input_path << ".\n";

    _neighbor_grid = std::make_unique<NeighborSearch>(_particles, _radius);
}

DFSPHSolver::DFSPHSolver(const char* file_path)
    :DFSPHSolver(file_path, "../") {}

void DFSPHSolver::no_pressure_predict()
{
    int num_particles = static_cast<int>(_particles.size());
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_particles + threadsPerBlock - 1) / threadsPerBlock;

    no_pressure_predict_kernel <<< blocksPerGrid, threadsPerBlock >>>(d_particles,
                                                                      d_neighbor_indices,
                                                                      d_neighbor_counts,
                                                                      num_particles,
                                                                      _timestep,
                                                                      _viscosity,
                                                                      _stiffness,
                                                                      _particle_mass,
                                                                      _density0,
                                                                      _gravity,
                                                                      _radius,
                                                                      _max_neighbors);
    cudaDeviceSynchronize();
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
    ofs << "element vertex " << _particles.size() << "\n";
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "end_header\n";

    for (const auto& particle : _particles) 
    {
        ofs 
            << particle.position.x << " "
            << particle.position.y << " "
            << particle.position.z << "\n";
    }

    ofs.close();
}

void DFSPHSolver::apply_boundary_conditions()
{
    int num_particles = static_cast<int>(_particles.size());
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_particles + threadsPerBlock - 1) / threadsPerBlock;
    apply_boundary_conditions_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_particles);
    cudaDeviceSynchronize();
}

void DFSPHSolver::allocate_device_memory()
{
    size_t num_particles = _particles.size();
    cudaMalloc((void**)&d_particles, num_particles * sizeof(Particle));

    int* h_neighbor_counts = new int[num_particles];
    for (int i = 0; i < num_particles; ++i)
    {
        h_neighbor_counts[i] = static_cast<int>(_neighbor_grid->neighbor_indices[i].size());
        if (h_neighbor_counts[i] > _max_neighbors)
            _max_neighbors = h_neighbor_counts[i];
    }

    int* h_neighbor_indices = new int[num_particles * _max_neighbors];
    for (int i = 0; i < num_particles; ++i)
    {
        int count = h_neighbor_counts[i];
        for (int n = 0; n < count; ++n)
        {
            h_neighbor_indices[i * _max_neighbors + n] = _neighbor_grid->neighbor_indices[i][n];
        }
        for (int n = count; n < _max_neighbors; ++n)
        {
            h_neighbor_indices[i * _max_neighbors + n] = -1;
        }
    }

    cudaMalloc((void**)&d_neighbor_indices, num_particles * _max_neighbors * sizeof(int));
    cudaMalloc((void**)&d_neighbor_counts, num_particles * sizeof(int));

    cudaMemcpy(d_neighbor_indices, h_neighbor_indices, num_particles * _max_neighbors * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbor_counts, h_neighbor_counts, num_particles * sizeof(int), cudaMemcpyHostToDevice);

    delete[] h_neighbor_indices;
    delete[] h_neighbor_counts;
}


void DFSPHSolver::copy_data_to_device()
{
    size_t num_particles = _particles.size();
    cudaMemcpy(d_particles,
               _particles.data(),
               num_particles * sizeof(Particle),
               cudaMemcpyHostToDevice);
}

void DFSPHSolver::copy_data_to_host()
{
    size_t num_particles = _particles.size();
    cudaMemcpy(_particles.data(),
               d_particles,
               num_particles * sizeof(Particle),
               cudaMemcpyDeviceToHost);
}

void DFSPHSolver::free_device_memory()
{
    cudaFree(d_particles);
    cudaFree(d_neighbor_indices);
    cudaFree(d_neighbor_counts);
}

void DFSPHSolver::compute_density()
{
    int num_particles = static_cast<int>(_particles.size());
    int threadsPerBlock = 256;
    int blocksPerGrid = (num_particles + threadsPerBlock - 1) / threadsPerBlock;

    compute_density_kernel <<<blocksPerGrid, threadsPerBlock >>> (d_particles,
            d_neighbor_indices,
            d_neighbor_counts,
            num_particles,
            _particle_mass,
            _radius,
            _max_neighbors);
    cudaDeviceSynchronize();
}

void DFSPHSolver::simulate()
{
    allocate_device_memory();
    copy_data_to_device();

    for (int i = 0; i < _framenum; ++i)
    {
        _neighbor_grid = std::make_unique<NeighborSearch>(_particles, _radius);

        free_device_memory();
        allocate_device_memory();
        copy_data_to_device();

        compute_density();
        no_pressure_predict();
        apply_boundary_conditions();
        copy_data_to_host();

        std::ostringstream filename;
        filename << _output_path << "/frame_" << i + 1 << ".ply";
        export_to_ply(filename.str());
    }

    free_device_memory();
}
