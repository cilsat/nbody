#ifndef __NBSYS_H__
#define __NBSYS_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

typedef float del_t;
typedef float dist_t;

typedef struct body {
    float px, py, pz;
    float vx, vy, vz;
    float ax, ay, az;
    float m;
} body_t;

typedef struct node {
    int id;
    float px, py, pz;
    float cx, cy, cz;
    float tm;

    float length;
    int depth, max_depth;

    body_t *bodies;
    int num_bodies;
    struct node *child;
    int num_child;
} node_t;

typedef struct nbodysys {
    body_t *bodies;
    int num_bodies;
    float maxp, maxv, maxm;
} nbodysys_t;

typedef void (*update) (nbodysys_t *nb, int iters, del_t time);

// body methods
body_t *init_rand_body(float, float, float);
void update_p(body_t*, del_t);

// node methods
node_t *init_node(int, int);
void fin_node(node_t*);
void set_node_members(node_t*, body_t*, float, float, float, float, int);
void set_node_children(node_t*);
void print_node_members(node_t*);
void check_node(node_t*, body_t*);

// nbodysys methods
nbodysys_t *init_rand_bodysys(int, float, float, float);
nbodysys_t *copy_bodysys(nbodysys_t*);
void fin_nbodysys(nbodysys_t*);
void print_nbodysys(nbodysys_t*);

// update methods
void brute(nbodysys_t*, int iter, del_t);
void all_seq(nbodysys_t*, int iter, del_t);
void barnes(nbodysys_t*, int iter, del_t);

#endif
