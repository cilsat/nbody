#include "nbodysys.h"

#define MAX_VAL_X 100   // maximum positional value
#define MAX_VAL_V 10    // maximum (initial) velocity
#define MAX_VAL_M 10    // maximum mass of body
#define G 6.67408e-11   // universal gravitational constant

NBodySys::NBodySys(int n) {
    int i;
    pos = new float3[n];
    mass = new float[n];
    vel = new float3[n];
    force = new float3[n];
    num = n;
    g = 1.f;

    srand(time(NULL));
    for (i = 0; i < n; i++) {
        pos[i] = rand_float3(MAX_VAL_X) - (float)0.5*MAX_VAL_X;
        mass[i] = MAX_VAL_M*((float) rand() / (float) RAND_MAX);
        vel[i] = rand_float3(MAX_VAL_V) - (float)0.5*MAX_VAL_V;
        force[i] = const_float3(0.f);
    }
}

NBodySys::NBodySys(const NBodySys &obj) {
    int n = obj.num;

    memcpy(pos, obj.pos, n*sizeof(float3));
    memcpy(mass, obj.mass, n*sizeof(float));
    memcpy(vel, obj.vel, n*sizeof(float3));
    num = n;
    g = obj.g;
}

NBodySys::~NBodySys() {
    delete[] pos;
    delete[] mass;
    delete[] vel;
    delete[] force;
}

void NBodySys::print() {
    int i;
    int n = num;

    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%.3f %.3f %.3f\n", pos[i].x, pos[i].y, pos[i].z);
    }
    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%.3f %.3f %.3f\n", vel[i].x, vel[i].y, vel[i].z);
    }
    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%.3f ", mass[i]);
    }
    printf("\n");
}

void NBodySys::update_full_seq(float t) {
    int n = num;
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                float3 xr = pos[j] - pos[i];
                float r = pow(len_float3(xr), 3);
                force[i] += xr*mass[j]/r;
            }
        }
    }

    for (i = 0; i < n; i++) {
        float3 vu = vel[i];
        vel[i] += g*force[i]*t;
        pos[i] += 0.5*(vel[i] - vu)*t;
        force[i] = const_float3(0);
    }
}

void NBodySys::update_full_par(float t) {
    int n = num;
    int i, j;

#pragma omp parallel for
    for (i = 0; i < n; i++) {
#pragma omp parallel for
        for (j = 0; j < n; j++) {
            if (i != j) {
                float3 xr = pos[j] - pos[i];
                float r = pow(len_float3(xr), 3);
                force[i] += xr*mass[j]/r;
            }
        }
    }

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        float3 vu = vel[i];
        vel[i] += g*force[i]*t;
        pos[i] += 0.5*(vel[i] - vu)*t;
        force[i] = const_float3(0);
    }
}
