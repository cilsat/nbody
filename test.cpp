#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "vec_ops.h"

#define NB 1000  // number of bodies
#define ITERS 100   // number of iterations
#define DEL 1   // (integer) length of time steps in seconds
#define MAX_VAL_X 100   // maximum positional value
#define MAX_VAL_V 10    // maximum (initial) velocity
#define MAX_VAL_M 10    // maximum mass of body
#define G 6.67408e-11   // universal gravitational constant
#define DEBUG 0

typedef struct {
    float3 *x;
    float *m;
    float3 *v;
    int num;
    float g;
} nbodysys;

nbodysys *init(int n) {
    int i;
    nbodysys *s = new nbodysys;

    srand(time(NULL));
    s->x = new float3[n];
    s->m = new float[n];
    s->v = new float3[n];
    s->num = n;
    s->g = 1.f;

    for (i = 0; i < n; i++) {
        s->x[i] = rand_float3(MAX_VAL_X) - (float)0.5*MAX_VAL_X;
        s->m[i] = MAX_VAL_M*((float) rand() / (float) RAND_MAX);
        s->v[i] = rand_float3(MAX_VAL_V) - (float)0.5*MAX_VAL_V;
    }

    return s;
}

void fin(nbodysys *s) {
    delete[] s->x;
    delete[] s->m;
    delete[] s->v;
    delete s;
}

nbodysys *copy(nbodysys *s) {
    int n = s->num;
    nbodysys *temp = init(n);

    memcpy(temp->x, s->x, n*sizeof(float3));
    memcpy(temp->m, s->m, n*sizeof(float));
    memcpy(temp->v, s->v, n*sizeof(float3));
    temp->num = n;
    temp->g = s->g;

    return temp;
}

void print(nbodysys *s) {
    int i;
    int n = s->num;

    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%.3f %.3f %.3f\n", s->x[i].x, s->x[i].y, s->x[i].z);
    }
    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%.3f %.3f %.3f\n", s->v[i].x, s->v[i].y, s->v[i].z);
    }
    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%.3f ", s->m[i]);
    }
    printf("\n");
}

void update_full_seq(nbodysys *s, float t) {
    int n = s->num;
    int i, j;
    float g = s->g;
    float3 *F = (float3 *) malloc(s->num * sizeof(float3));

    for (i = 0; i < n; i++) {
        F[i] = const_float3(0);
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                // r12 = r1 - r2
                float3 xr = s->x[j] - s->x[i];
                // |r12| = (r12.x^2 + r12.y^2 + r12.z^2)^0.5
                float r = len_float3(xr);
                // ^r12 = r12/|r12|
                float3 u = xr/r;
                // F = G*m1*m2*^r12/|r12|
                float3 f = u*g*s->m[i]*s->m[j]/r;
                F[i] = F[i] + f;
                F[j] = F[j] - f;
            }
        }
    }

    for (i = 0; i < n; i++) {
        // a = F/m
        float3 a = F[i]/s->m[i];
        //printf("%.3f %.3f %.3f\n", a.x, a.y, a.z);
        // v = u + a*t
        float3 vf = s->v[i] + a*t;
        // s = u*t + 0.5*a*t^2
        s->x[i] = s->x[i] + 0.5*(vf - s->v[i])*t;
        s->v[i] = vf;
    }

    free(F);
}

void update_full_par(nbodysys *s, float t) {
    int n = s->num;
    int i, j;
    float g = s->g;
    float3 *F = (float3 *) malloc(s->num * sizeof(float3));

    for (i = 0; i < n; i++) {
        F[i] = const_float3(0);
    }

#pragma omp parallel for
    for (i = 0; i < n; i++) {
#pragma omp parallel for
        for (j = 0; j < n; j++) {
            if (i != j) {
                // r12 = r1 - r2
                float3 xr = s->x[j] - s->x[i];
                // |r12| = (r12.x^2 + r12.y^2 + r12.z^2)^0.5
                float r = len_float3(xr);
                // ^r12 = r12/|r12|
                float3 u = xr/r;
                // F = G*m1*m2*^r12/|r12|
                float3 f = u*g*s->m[i]*s->m[j]/r;
                F[i] = F[i] + f;
                F[j] = F[j] - f;
            }
        }
    }

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        // a = F/m
        float3 a = F[i]/s->m[i];
        //printf("%.3f %.3f %.3f\n", a.x, a.y, a.z);
        // v = u + a*t
        float3 vf = s->v[i] + a*t;
        // s = u*t + 0.5*a*t^2
        s->x[i] = s->x[i] + 0.5*(vf - s->v[i])*t;
        s->v[i] = vf;
    }

    free(F);
}

int main(int argc, char **argv) {
    double dstart, dstop;
    int nb, iters;

    if (argc > 1) {
        nb = strtol(argv[1], NULL, 10);
        iters = strtol(argv[2], NULL, 10);
    } else {
        nb = NB; iters = ITERS;
    }
    nbodysys *sys = init(nb);
    nbodysys *temp = copy(sys);

    if (DEBUG == 1) {
    print(sys);
    print(temp);
    }

    dstart = omp_get_wtime();
    for (int i = 0; i < iters; i++) {
        update_full_seq(sys, (float)DEL);
    }
    dstop = omp_get_wtime();
    printf("%.3f\n", dstop-dstart);

    dstart = omp_get_wtime();
    for (int i = 0; i < iters; i++) {
        update_full_par(temp, (float)DEL);
    }
    dstop = omp_get_wtime();
    printf("%.3f\n", dstop-dstart);

    if (DEBUG == 1) {
        print(sys);
        print(temp);
    }

    fin(sys);
    fin(temp);

    return 0;
}
