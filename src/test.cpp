#include "DFSPHSolver.cuh"

int main()
{
	DFSPHSolver solver("../sphere.ply");
	//DFSPHSolver solver(1000.0, 1000.0);
	//solver.export_to_ply(".\\particles_1.ply");
	solver.simulate();
}