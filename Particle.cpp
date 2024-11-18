#include "Particle.h"

Particle::Vector3D operator-(const Particle::Vector3D& self, const Particle::Vector3D& other)
{
	Particle::Vector3D result(3, 0.0);

	if (self.size() == other.size() && self.size() == 3)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			result[i] = self[i] - other[i];
		}
	}
	else
	{
		throw "vector size mismatch";
	}
	return result;
}

Particle::Vector3D operator-(double self, const Particle::Vector3D& other)
{
	Particle::Vector3D result(3, 0.0);

	if (other.size() == 3)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			result[i] = self - other[i];
		}
	}
	else
	{
		throw "vector size mismatch";
	}
	return result;
}

Particle::Vector3D operator/(const Particle::Vector3D& other, double self)
{
	Particle::Vector3D result(3, 0.0);

	if (other.size() == 3)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			if (self != 0)
			{
				result[i] = other[i] / self;
			}
			else
			{
				result[i] = other[i] / 1e-6;
			}
		}
	}
	else
	{
		throw "vector size mismatch";
	}
	return result;
}

Particle::Vector3D operator*(double self, const Particle::Vector3D& other)
{
	Particle::Vector3D result(3, 0.0);

	if (other.size() == 3)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			result[i] = other[i] * self;
		}
	}
	else
	{
		throw "vector size mismatch";
	}
	return result;
}

double operator*(const Particle::Vector3D& self, const Particle::Vector3D& other)
{
	double result = 0.0;

	if (self.size() == other.size() && self.size() == 3)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			result += self[i] * other[i];
		}
	}
	else
	{
		throw "vector size mismatch";
	}
	return result;
}

void operator+=(Particle::Vector3D& self, const Particle::Vector3D& other)
{
	if (self.size() == other.size() && self.size() == 3)
	{
		for (size_t i = 0; i < 3; ++i)
		{
			self[i] += other[i];
		}
	}
	else
	{
		throw "vector size mismatch";
	}
}

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
