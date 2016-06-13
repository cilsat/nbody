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

typedef struct node {
    uint32_t id;
    float px, py, pz;
    float cx, cy, cz;
    float tm;

    float length;
    uint16_t depth, max_depth;

    uint32_t *bodies;
    uint32_t num_bodies;

    uint32_t *quad[8];
    uint32_t num_quad[8];

    struct node *child[8];
    uint32_t num_child;
} node_t;

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
void init_node(node_t *node, uint8_t dep, uint8_t max_dep, body_t *np_bodies, uint32_t* p_bodies, uint32_t p_nbodies, float p_x, float p_y, float p_z, float p_length);
void free_node(node_t *n);
void check_node(node_t *node, body_t *body);
void check_node_ordered(node_t*, body_t*, uint32_t*, uint32_t*);
void print_node(node_t*, body_t*);

// nbodysys methods
nbodysys_t *init_rand_bodysys(uint32_t, float, float, float);
nbodysys_t *copy_nbodysys(nbodysys_t*);
void free_nbodysys(nbodysys_t*);
void print_nbodysys(nbodysys_t*);

// update methods
void brute(nbodysys_t*, uint32_t iter, del_t);
void barnes(nbodysys_t*, uint32_t iter, del_t);
void brute_ordered(nbodysys_t*, uint32_t iter, del_t, uint32_t**);
void barnes_ordered(nbodysys_t*, uint32_t iter, del_t, uint32_t**);

#endif
