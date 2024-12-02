#include "Vector3D.cuh"

__host__ __device__ Vector3D::Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

__host__ __device__ Vector3D::Vector3D() : Vector3D(0.0, 0.0, 0.0) {}

__host__ __device__ Vector3D Vector3D::operator-(const Vector3D& other) const
{
    return {x - other.x, y - other.y, z - other.z};
}

__host__ __device__ Vector3D Vector3D::operator*(double scalar) const
{
    return {x * scalar, y * scalar, z * scalar};
}

__host__ __device__ Vector3D Vector3D::operator/(double scalar) const
{
    return {x / scalar, y / scalar, z / scalar};
}

__host__ __device__ Vector3D& Vector3D::operator+=(const Vector3D& other)
{
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

__host__ __device__ double Vector3D::dot(const Vector3D& other) const
{
    return x * other.x + y * other.y + z * other.z;
}

__host__ __device__ double Vector3D::length() const
{
    return sqrt(x * x + y * y + z * z);
}

__host__ __device__ Vector3D operator-(double other, const Vector3D& self)
{
    return {other-self.x, other - self.y, other - self.z};
}