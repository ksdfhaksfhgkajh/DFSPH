#ifndef PARTICLE_CUH
#define PARTICLE_CUH

#include <vector>
#include <cmath>
#include <stdexcept>
#include "Vector3D.cuh"

enum Axle :size_t { X, Y, Z };
enum Kappa_t :size_t {Is_vel, Is_pho};

class Particle
{
public:
    enum Particle_t :size_t {Fluid, Boundary};

	Vector3D position{};
	Vector3D velocity{};
	double density{};
	double factor{};
    double kappa{};
    Particle_t type{};

	double cal_distance(const Particle& other) const;
};

#endif