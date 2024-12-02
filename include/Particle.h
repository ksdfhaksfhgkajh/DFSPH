#ifndef PARTICLE_CUH
#define PARTICLE_CUH

#include <vector>
#include <cmath>
#include <stdexcept>
#include "Vector3D.cuh"

class Particle
{
public:
	Vector3D position{};
	Vector3D velocity{};
	Vector3D acceleration{};
	double density{};
	double factor{};

	double cal_distance(const Particle& other) const;
};

enum Axle :size_t { X, Y, Z };

#endif