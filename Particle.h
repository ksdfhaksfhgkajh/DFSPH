#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <cmath>

class Particle
{
public:
	std::vector<double> position{};
	std::vector<double> velocity{};
	std::vector<double> acceleration{};
	double density{};
	double factor{};

	double cal_distance(const Particle& other) const;
};

enum Axle :size_t { X, Y, Z };

#endif