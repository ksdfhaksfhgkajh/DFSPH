#include "DFSPHSolver.h"
#include <chrono>

int main()
{
    auto start = std::chrono::high_resolution_clock::now();

    DFSPHSolver solver("../sphere.ply",
                       "../bound.ply",
                       "../result");
	//DFSPHSolver solver(1000.0, 1000.0);
	//solver.export_to_ply(".\\particles_1.ply");
	solver.simulate();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << "Elapsed time: " << duration << " ms\n";
    return 0;
}