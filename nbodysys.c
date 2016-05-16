#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define G 6.67408e-11   // universal gravitational constant

typedef float del_t;

typedef struct body {
    float px, py, pz;
    float vx, vy, vz;
    float m;
} body_t;

typedef struct node {
    float px, py, pz;
    float cx, cy, cz;
    float tm;

    int length;
    int depth;

    body_t *bodies;
    int num_bodies;
    struct node *child;
    int num_child;
} node_t;

typedef struct nbodysys {
    body_t *bodies;
    int num_bodies;
    float *p;
    float *v;
    float *m;
    float maxp, maxv, maxm;
} nbodysys_t;

body_t *init_rand_body(float max_p, float max_v, float max_m) {
    body_t *temp = (body_t *) malloc(sizeof(body_t));
    temp->px = max_p*((float) rand() / (float) RAND_MAX) - 0.5*max_p;
    temp->py = max_p*((float) rand() / (float) RAND_MAX) - 0.5*max_p;
    temp->pz = max_p*((float) rand() / (float) RAND_MAX) - 0.5*max_p;
    temp->vx = max_v*((float) rand() / (float) RAND_MAX) - 0.5*max_v;
    temp->vy = max_v*((float) rand() / (float) RAND_MAX) - 0.5*max_v;
    temp->vz = max_v*((float) rand() / (float) RAND_MAX) - 0.5*max_v;
    temp->m = max_m*((float) rand() / (float) RAND_MAX);

    return temp;
}

nbodysys_t *init_rand_nbodysys(int n, float max_p, float max_v, float max_m) {
    int i;

    nbodysys_t *temp = (nbodysys_t *) malloc(sizeof(nbodysys_t));
    temp->bodies = (body_t *) malloc(n*sizeof(body_t));
    temp->p = (float *) malloc(n*3*sizeof(float));
    temp->v = (float *) malloc(n*3*sizeof(float));
    temp->m = (float *) malloc(n*sizeof(float));
    temp->maxp = max_p; temp->maxv = max_v; temp->maxm = max_m;
    temp->num_bodies = n;

    srand(time(NULL));
    for (i = 0; i < n; i++) {
        temp->bodies[i] = *init_rand_body(max_p, max_v, max_m);
        temp->p[i*3] = temp->bodies[i].px;
        temp->p[i*3 + 1] = temp->bodies[i].py;
        temp->p[i*3 + 2] = temp->bodies[i].pz;
        temp->v[i] = temp->bodies[i].vx;
        temp->v[i*3 + 1] = temp->bodies[i].vy;
        temp->v[i*3 + 2] = temp->bodies[i].vz;
        temp->m[i] = temp->bodies[i].m;
    }

    return temp;
}

void fin_nbodysys(nbodysys_t *nb) {
    free(nb->bodies);
    free(nb->p);
    free(nb->v);
    free(nb->m);
    free(nb);
}

void print_nbodysys(nbodysys_t *nb) {
    for (int i = 0; i < nb->num_bodies; i++) {
        for (int j = 0 ; j < 3; j++) {
            printf("%.3f ", nb->p[i*3 + j]);
        }
        printf("\n");
    }
    printf("\n");
}

typedef void (*update) (nbodysys_t *nb, del_t time);

void all_seq(nbodysys_t *nb, del_t time) {
    int n = nb->num_bodies;
    int i, j;
    float *f = (float *) malloc(n*3*sizeof(float));
    float *p = (float *) malloc(n*3*sizeof(float));
    float *v = (float *) malloc(n*3*sizeof(float));
    memcpy(p, nb->p, n*3*sizeof(float));
    memcpy(v, nb->v, n*3*sizeof(float));
    float gt = 1.f*time;
    float ht = 0.5f*time;

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        f[i*3] = 0;
        f[i*3+1] = 0;
        f[i*3+2] = 0;
#pragma omp parallel for
        for (j = 0; j < n; j++) {
            if (i != j) {
                float rx = p[j*3] - p[i*3];
                float ry = p[j*3+1] - p[i*3+1];
                float rz = p[j*3+2] - p[i*3+2];
                float r = 1.f/sqrt(rx*rx + ry*ry + rz*rz);
                float mr = nb->m[j]*r*r*r;
                f[i*3] = mr*rx;
                f[i*3+1] = mr*ry;
                f[i*3+2] = mr*rz;
            }
        }
        float ux = v[i*3];
        float uy = v[i*3+1];
        float uz = v[i*3+2];
        v[i*3] += f[i*3]*gt;
        v[i*3+1] += f[i*3+1]*gt;
        v[i*3+2] += f[i*3+2]*gt;
        p[i*3] += (v[i*3] - ux)*ht;
        p[i*3+1] += (v[i*3+1] - uy)*ht;
        p[i*3+2] += (v[i*3+2] - uz)*ht;
    }

    free(f);
    free(p);
    free(v);
}
