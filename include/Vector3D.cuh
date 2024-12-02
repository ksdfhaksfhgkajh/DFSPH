#ifndef VECTOR3D_CUH
#define VECTOR3D_CUH

#ifndef __CUDACC__
#define __host__
#define __device__
#endif

struct Vector3D
{
    double x;
    double y;
    double z;

    __host__ __device__ Vector3D(double x_, double y_, double z_);

    __host__ __device__ Vector3D();

    __host__ __device__ Vector3D operator-(const Vector3D &other) const;

    __host__ __device__ Vector3D operator*(double scalar) const;

    __host__ __device__ Vector3D operator/(double scalar) const;

    __host__ __device__ Vector3D& operator+=(const Vector3D& other);

    __host__ __device__ double dot(const Vector3D& other) const;

    __host__ __device__ double length() const;
};

__host__ __device__ Vector3D operator-(double other, const Vector3D& self);

#endif