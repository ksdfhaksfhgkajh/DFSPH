#include "Particle.h"

double Particle::cal_distance(const Particle& other) const
{
    return (position - other.position).length();
}