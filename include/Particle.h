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
	double density{};
	double factor{};
    double kappa{};

	double cal_distance(const Particle& other) const;
};

enum Axle :size_t { X, Y, Z };
enum Kappa_t :size_t {Is_vel, Is_pho};

#endif