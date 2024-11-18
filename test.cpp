#include "DFSPHSolver.h"

int main()
{
	DFSPHSolver solver(1000, 1000);
	solver.export_to_ply(".\\particles_1.ply");
}