#include "DFSPHSolver.h"
#include <chrono>

int main()
{
    auto start = std::chrono::high_resolution_clock::now();

    DFSPHSolver solver("../particles/sphere.ply",
                       "../particles/bound.ply",
                       "../result");
	solver.simulate();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cout << "Elapsed time: " << duration << " ms\n";
    return 0;
}