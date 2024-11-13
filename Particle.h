#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

class Particle
{
public:
	enum axle :size_t {x, y ,z};
	std::vector<double> position;
	std::vector<double> velocity;
	std::vector<double> acceleration;
	double density;
	double factor;

	Particle(std::vector<double> position, std::vector<double> velocity);
};

#endif