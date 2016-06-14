#ifndef __NBSYS_H__
#define __NBSYS_H__

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

typedef float del_t;

typedef struct body {
    float px, py, pz;
    float vx, vy, vz;
    float ax, ay, az;
    float m;
} body_t;

typedef struct node node_t;
struct node {
    float cx, cy, cz;
    float gm;
    uint8_t depth, num_child;
    node_t *child;
};

typedef struct nbodysys {
    body_t *bodies;
    uint32_t num_bodies;
    float maxp, maxv, maxm;
} nbodysys_t;

typedef void (*update) (nbodysys_t *nb, uint32_t iters, del_t time, uint32_t* opt);

// body methods
body_t *init_rand_body(float max_p, float max_v, float max_m);
void update_body(body_t* b, del_t t);

// node methods
void init_node(node_t *node, body_t **bodies, uint32_t num_bodies, uint8_t depth, uint8_t max_depth, float px, float py, float pz, float *length, float g);
void free_node(node_t *node);
void check_node(node_t *node, body_t *body, float length[]);
void print_node(node_t*, body_t*);

// nbodysys methods
nbodysys_t *init_rand_bodysys(uint32_t, float, float, float);
nbodysys_t *copy_nbodysys(nbodysys_t*);
void free_nbodysys(nbodysys_t*);
void print_nbodysys(nbodysys_t*);

// update methods
void brute(nbodysys_t*, uint32_t iter, del_t);
void barnes(nbodysys_t*, uint32_t iter, del_t);

#endif
