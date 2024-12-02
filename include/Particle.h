#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <cmath>

class Particle
{
public:
	using Vector3D = std::vector<double>;

	Vector3D position{};
	Vector3D velocity{};
	Vector3D acceleration{};
	double density{};
	double factor{};

	double cal_distance(const Particle& other) const;
};

enum Axle :size_t { X, Y, Z };

Particle::Vector3D operator-(const Particle::Vector3D& self, const Particle::Vector3D& other);

Particle::Vector3D operator-(double self, const Particle::Vector3D& other);

Particle::Vector3D operator/(const Particle::Vector3D& other, double self);

Particle::Vector3D operator*(double self, const Particle::Vector3D& other);

double operator*(const Particle::Vector3D& self, const Particle::Vector3D& other);

void operator+=(Particle::Vector3D& self, const Particle::Vector3D& other);

#endif