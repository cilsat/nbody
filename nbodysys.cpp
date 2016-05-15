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
    num = n;
    g = 1.f;

    srand(time(NULL));
    for (i = 0; i < n; i++) {
        pos[i] = rand_float3(MAX_VAL_X) - (float)0.5*MAX_VAL_X;
        mass[i] = MAX_VAL_M*((float) rand() / (float) RAND_MAX);
        vel[i] = rand_float3(MAX_VAL_V) - (float)0.5*MAX_VAL_V;
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
    float3 *force = new float3[n];

    for (i = 0; i < n; i++) {
        force[i] = const_float3(0.f);
        for (j = 0; j < n; j++) {
            if (i != j) {
                float3 xr = pos[j] - pos[i];
                float r = len_float3(xr);
                force[i] += xr*(mass[j]/(r*r*r));
            }
        }
    }

    for (i = 0; i < n; i++) {
        float3 vu = vel[i];
        vel[i] += force[i]*(g*t);
        pos[i] += (vel[i] - vu)*(0.5*t);
    }

    delete[] force;
}

void NBodySys::update_full_par(float t) {
    int n = num;
    int i, j;
    float3 *force = new float3[n];
    float gt = g*t;
    float ht = 0.5f*t;

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        force[i] = const_float3(0.f);
#pragma omp parallel for
        for (j = 0; j < n; j++) {
            if (i != j) {
                float3 xr = {pos[j].x-pos[i].x, pos[j].y-pos[i].y, pos[j].z-pos[i].z};
                float r = len_float3(xr);
                float gr = mass[j]/(r*r*r);
                float3 f = {xr.x*gr, xr.y*gr, xr.z*gr};
                force[i] += f;
            }
        }
    }

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        float3 vu = vel[i];
        float3 vt = {force[i].x*gt, force[i].y*gt, force[i].z*gt};
        float3 vr = {(vt.x-vu.x)*ht, (vt.y-vu.y)*ht, (vt.z-vu.z)*ht};
        vel[i] += vt;
        pos[i] += vr;
    }

    delete[] force;
}
