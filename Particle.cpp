#include "Particle.h"

double Particle::cal_distance(const Particle& other) const
{
	double distance = 0;
	if (this->position.size() == other.position.size())
	{
		for (size_t i = 0; i < position.size(); ++i)
		{
			distance += pow((this->position[i] - other.position[i]), 2);
		}
	}
	else
	{
		throw "vector size mismatch";
	}
	return sqrt(distance);
}