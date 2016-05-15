#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define NB 100  // number of bodies
#define ND 3    // number of dimensions
#define DEL 1   // (integer) length of time steps in seconds
#define MAX_VAL_X 100   // maximum positional value
#define MAX_VAL_V 10    // maximum (initial) velocity
#define MAX_VAL_M 10    // maximum mass of body

typedef struct {
    float3 *x;
    float *m;
    float3 *v;
    float3 *F;
    int num;
    float g;
} nbodysys;

nbodysys *init(int n) {
    int i, j;
    nbodysys *s = (nbodysys *)malloc(sizeof(nbodysys));

    srand(time(NULL));
    s->x = (float3 *)malloc(n*sizeof(float3));
    s->m = (float *)malloc(n*sizeof(float));
    s->v = (float *)malloc(n*ND*sizeof(float));
    s->F = (float *)malloc(n*ND*sizeof(float));
    s->num = n;
    s->g = 6.67408*pow(10, -11);

    for (i = 0; i < n; i++) {
        for (j = 0; j < ND; j++) { 
            s->x[i*ND + j] = MAX_VAL_X*((float) rand()/(float) RAND_MAX) - 0.5*MAX_VAL_X;
            s->v[i*ND + j] = MAX_VAL_V*((float) rand()/(float) RAND_MAX);
            s->F[i*ND + j] = 0;
        }
        s->m[i] = MAX_VAL_M*((float) rand()/(float) RAND_MAX);
    }

    return s;
}

void fin(nbodysys *s) {
    free(s->x);
    free(s->m);
    free(s->v);
    free(s);
}

void update(nbodysys *s) {
    int i, j, k;
    float d, f, a;
    float *xr = (float *) malloc(ND*sizeof(float));
    float *a = (float *) malloc(ND*sizeof(float));

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            d = 0;
            for (k = 0; k < ND; k++) {
                xr[k] = s->x[j*ND + k] - s->x[i*ND + k];
                d += powf(xr[k], 2);
            }
            d = powf(d, 0.5);
            f = s->g*(s->m[j] - s->m[i])/d;
            for (k = 0; k < ND; k++) {
                a[k] = 

                

void print(nbodysys *s) {
    int i, j;
    int n = s->num;

    for (i = 0; i < n; i++) {
        for (j = 0; j < ND; j++) { 
            printf("%.3f ", s->x[i*ND + j]);
        }
        printf("\n");
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < ND; j++) { 
            printf("%.3f ", s->v[i*ND + j]);
        }
        printf("\n");
    }
    for (i = 0; i < n; i++) {
        printf("%.3f ", s->m[i]);
    }
}

int main(int argc, char **argv) {
    nbodysys *sys = init(NB);
    print(sys, NB);
    fin(sys);

    return 0;
}
