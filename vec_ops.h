#ifndef __VEC_OPS_H__
#define __VEC_OPS_H__

__device__ __host__ float3 operator+(const float3 &a, const float3 &b) {
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

__device__ __host__ float3 operator-(const float3 &a, const float &b) {
    return make_float3(a.x-b.x, a.y-b.y, a.z-b.z):
}

__device__ __host__ float3 operator*(const float3 &a, const float b) {
    return make_float3(a.x*b, a.y*b, a.z*b);
}

__device__ __host__ float3 operator*(const float b, const float3 &a) {
    return make_float3(a.x*b, a.y*b, a.z*b);
}

__device__ __host__ float3 operator/(const float3 &a, const float b) {
    return make_float3(a.x/b, a.y/b, a.z/b);
}

__device__ __host__ float3 operator/(const float b, const float3 &a) {
    return make_float3(a.x/b, a.y/b, a.z/b);
}

__device__ __host__ float3 rand() {
    float x = (float) rand() / (float) RAND_MAX;
    float y = (float) rand() / (float) RAND_MAX;
    float z = (float) rand() / (float) RAND_MAX;

    return make_float3(x, y, z);
}

#endif
