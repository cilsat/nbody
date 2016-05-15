#ifndef __VEC_OPS_H__
#define __VEC_OPS_H__

#include <math.h>
#include <time.h>
#include <stdlib.h>

struct float3 {
    float x;
    float y;
    float z;
};

inline float3 make_float3(const float &x, const float &y, const float &z) {
    float3 temp = {x, y, z};
    return temp;
}

// ADDITION
inline float3 operator+(const float3 &a, const float &b) {
    return make_float3(a.x+b, a.y+b, a.z+b);
}
inline float3 operator+(const float &b, const float3 &a) {
    return make_float3(a.x+b, a.y+b, a.z+b);
}
inline float3 operator+(const float3 &a, const float3 &b) {
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}
inline void operator+=(float3 &a, const float3 &b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
}
inline void operator+=(float3 &a, const float &b) {
    a.x += b; a.y += b; a.z += b;
}

// SUBTRACTION
inline float3 operator-(const float3 &a, const float &b) {
    return make_float3(a.x-b, a.y-b, a.z-b);
}
inline float3 operator-(const float &b, const float3 &a) {
    return make_float3(a.x-b, a.y-b, a.z-b);
}
inline float3 operator-(const float3 &a, const float3 &b) {
    return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
}

// MULTIPLICATION
inline float3 operator*(const float3 &a, const float &b) {
    return make_float3(a.x*b, a.y*b, a.z*b);
}
inline float3 operator*(const float &b, const float3 &a) {
    return make_float3(a.x*b, a.y*b, a.z*b);
}
inline float3 operator*(const float3 &b, const float3 &a) {
    return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);
}

// DIVISION
inline float3 operator/(const float3 &a, const float &b) {
    return make_float3(a.x/b, a.y/b, a.z/b);
}
inline float3 operator/(const float &b, const float3 &a) {
    return make_float3(a.x/b, a.y/b, a.z/b);
}
inline float3 operator/(const float3 &a, const float3 &b) {
    return make_float3(a.x/b.x, a.y/b.x, a.z/b.x);
}

// MAGNITUDE
inline float len_float3(const float3 &a) {
    return pow(a.x*a.x + a.y*a.y + a.z*a.z, 0.5);
}

// INSTANSIATION
float3 const_float3(const float c) {
    return make_float3(c, c, c);
}
float3 rand_float3(float max) {
    float x = max*((float) rand() / (float) RAND_MAX);
    float y = max*((float) rand() / (float) RAND_MAX);
    float z = max*((float) rand() / (float) RAND_MAX);
    return make_float3(x, y, z);
}

#endif
