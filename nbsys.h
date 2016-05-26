#ifndef __NBSYS_H__
#define __NBSYS_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

typedef double del_t;
typedef double dist_t;

typedef struct body {
    double px, py, pz;
    double vx, vy, vz;
    double ax, ay, az;
    double m;
    double gmr;
} body_t;

typedef struct node {
    int id;
    double px, py, pz;
    double cx, cy, cz;
    double tm;

    double length;
    int depth, max_depth;

    body_t *bodies;
    int num_bodies;
    struct node *child;
    int num_child;
} node_t;

typedef struct nbodysys {
    body_t *bodies;
    int num_bodies;
    double maxp, maxv, maxm;
} nbodysys_t;

typedef void (*update) (nbodysys_t *nb, int iters, del_t time);

// body methods
body_t *init_rand_body(double, double, double);
void update_p(body_t*, del_t);

// node methods
node_t *init_node(int, int);
void fin_node(node_t*);
void set_node_members(node_t*, body_t*, double, double, double, double, int);
void set_node_children(node_t*);
void print_node_members(node_t*);
void check_node(node_t*, body_t*);

// nbodysys methods
nbodysys_t *init_rand_bodysys(int, double, double, double);
nbodysys_t *copy_bodysys(nbodysys_t*);
void fin_nbodysys(nbodysys_t*);
void print_nbodysys(nbodysys_t*);

// update methods
void brute(nbodysys_t*, int iter, del_t);
void all_seq(nbodysys_t*, int iter, del_t);
void barnes(nbodysys_t*, int iter, del_t);

#endif
