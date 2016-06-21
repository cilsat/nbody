#ifndef __NBSYS_H__
#define __NBSYS_H__

#include <jemalloc/jemalloc.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

typedef float del_t;

typedef struct {
    float px, py, pz;
    float ax, ay, az;
    float m;
    float vx, vy, vz;
} body_t;

typedef struct node node_t;
//#pragma pack(push)
//#pragma pack(4)
struct node {
    float cx, cy, cz;
    float gm;
    node_t *child;
    uint8_t depth, num_child;
};
//#pragma pack(pop)

typedef struct {
    signed char x[8], y[8], z[8];
} ttable_t;

typedef struct {
    body_t *bodies;
    uint32_t num_bodies;
    float maxp, maxv, maxm;
} nbodysys_t;

typedef void (*update) (nbodysys_t *nb, uint32_t iters, del_t time, uint32_t* opt);

// body methods
body_t *init_rand_body(float max_p, float max_v, float max_m);
void update_body(body_t* b, del_t t);

// node methods
void init_node(node_t *node, body_t **bodies, uint32_t num_bodies, uint8_t depth, uint8_t max_depth, float pos[], ttable_t magic, float *length, float g);
void free_node(node_t *node);
void check_node(node_t *node, body_t *body, float *length);
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
