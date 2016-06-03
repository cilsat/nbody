#include "nbsys.h"
#include <assert.h>

#define G 0.1f // universal gravitational constant should be 6.67408e-11
#define MAX_P 100.f
#define MAX_V 10.f
#define MAX_M 100.f
#define E 50.f
#define E_SQR 2500.f // softening factor
#define DEBUG 2
#define DIST_THRES 0.1f

// debugging
static int node_id;
static int debug;
//static int node_counter;

/* carmack method for fast inverse square root (basically newton's method)
static inline float rsqrt(float x) {
    float hx = 0.5f*x;
    int i = *(int*)&x;
    i = 0x5f3759df - (i>>1);
    x = *(float*)&i;
    x *= (1.5f - hx*x*x);
    return x;
}*/

body_t *init_rand_body(float max_p, float max_v, float max_m) {
    body_t *temp = (body_t*)malloc(sizeof(body_t));
    temp->px = max_p*((float)rand() / (float)RAND_MAX) - 0.5f*max_p;
    temp->py = max_p*((float)rand() / (float)RAND_MAX) - 0.5f*max_p;
    temp->pz = max_p*((float)rand() / (float)RAND_MAX) - 0.5f*max_p;
    temp->vx = max_v*((float)rand() / (float)RAND_MAX) - 0.5f*max_v;
    temp->vy = max_v*((float)rand() / (float)RAND_MAX) - 0.5f*max_v;
    temp->vz = max_v*((float)rand() / (float)RAND_MAX) - 0.5f*max_v;
    temp->ax = 0.f;
    temp->ay = 0.f;
    temp->az = 0.f;
    temp->m = max_m*((float) rand() / (float) RAND_MAX);

    return temp;
}

inline void update_p(body_t* b, del_t t) {
    float ht = 0.5f*t;
    /*
    float hmaxp = MAX_P;
    if ((fabs(b->px) > hmaxp) && (b->vx*b->px > 0.f)) b->vx = -b->vx;
    if ((fabs(b->py) > hmaxp) && (b->vy*b->py > 0.f)) b->vy = -b->vy;
    if ((fabs(b->pz) > hmaxp) && (b->vz*b->pz > 0.f)) b->vz = -b->vz;
    */
    b->vx += b->ax*t;
    b->vy += b->ay*t;
    b->vz += b->az*t;
    b->px += b->vx*ht;
    b->py += b->vy*ht;
    b->pz += b->vz*ht;
}

void free_body(body_t *body) {
    free(body);
}

node_t *init_node(int dep, int max_dep, int* p_bodies, int p_nbodies) {
    node_t *temp = (node_t*)malloc(sizeof(node_t));
    temp->depth = dep;
    temp->max_depth = max_dep;

    temp->num_bodies = p_nbodies;
    temp->bodies = (int*)malloc(p_nbodies*sizeof(int));
    memcpy(temp->bodies, p_bodies, p_nbodies*sizeof(int));

    // quad keeps tracks of which bodies are in which quadrants by
    // index. num_quad keeps track of the number of bodies in each of
    // quad's quadrants.
    temp->num_quad = (int*)malloc(8*sizeof(int));
    temp->quad = (int**)malloc(8*sizeof(int*));
    for (int i = 0; i < 8; i++) {
        temp->quad[i] = (int*)malloc(p_nbodies*sizeof(int));
        for (int j = 0; j < p_nbodies; j++) {
            temp->quad[i][j] = 0;
        }
        temp->num_quad[i] = 0;
    }

    return temp;
}

void set_node(node_t *node, body_t *np_bodies, float p_x, float p_y, float p_z, float p_length) {
    int i;
    int n = node->num_bodies;
    node->id = node_id++;
    node->px = p_x;
    node->py = p_y;
    node->pz = p_z;
    node->length = p_length;
    node->num_child = 0;
    node->child = 0;

    if (n == 1) {
        node->cx = np_bodies[node->bodies[0]].px;
        node->cy = np_bodies[node->bodies[0]].py;
        node->cz = np_bodies[node->bodies[0]].pz;
        node->tm = np_bodies[node->bodies[0]].m;
    }
    else if (n > 1) {
        float t_tm = 0.f;
        float t_cx = 0.f;
        float t_cy = 0.f;
        float t_cz = 0.f;
        for (i = 0; i < n; i++) {
            float t_m = np_bodies[node->bodies[i]].m;
            t_cx += np_bodies[node->bodies[i]].px*t_m;
            t_cy += np_bodies[node->bodies[i]].py*t_m;
            t_cz += np_bodies[node->bodies[i]].pz*t_m;
            t_tm += t_m;
        }
        node->tm = t_tm;
        float invtm = 1.f/t_tm;
        node->cx = t_cx*invtm;
        node->cy = t_cy*invtm;
        node->cz = t_cz*invtm;

        if (node->depth < node->max_depth) {
            for (i = 0; i < n; i++) {
                // derive quadrant q=0..7 from relative position on each axis
                // of body to center of current node
                int b_x = np_bodies[node->bodies[i]].px < node->px ? 0 : 1;
                int b_y = np_bodies[node->bodies[i]].py < node->py ? 0 : 1;
                int b_z = np_bodies[node->bodies[i]].pz < node->pz ? 0 : 1;
                // each quadrant is mapped to a number between 0..7
                int q = b_x*4 + b_y*2 + b_z;

                // keep track of all bodies in each quadrant q
                node->quad[q][node->num_quad[q]] = i;
                // keep track of number of bodies in each quadrant q
                node->num_quad[q]++;
            }

            for (i = 0; i < 8; i++) {
                if (node->num_quad[i] > 0) {
                    // magic functions to determine new position based on
                    // quadrant and length/size of node
                    int m_x = i/4 == 0 ? -1 : 1;
                    int m_y = (i/2)%2 == 0 ? -1 : 1;
                    int m_z = i%2 == 0 ? -1 : 1;
                    float new_px = node->px + m_x*0.25f*node->length;
                    float new_py = node->py + m_y*0.25f*node->length;
                    float new_pz = node->pz + m_z*0.25f*node->length;

                    node->num_child++;
                    node->child = (node_t*)realloc(node->child, node->num_child*sizeof(node_t));
                    node_t *new_node = init_node(node->depth+1, node->max_depth, node->quad[i], node->num_quad[i]);
                    node->child[node->num_child-1] = *new_node;
                    free(new_node);
                    set_node(&node->child[node->num_child-1], np_bodies, new_px, new_py, new_pz, 0.5f*node->length);
                }
            }
        }
    }
}

void free_node(node_t *n) {
    free(n->bodies);
    for (int i = 0; i < 8; i++) {
        free(n->quad[i]);
    }
    free(n->quad);
    free(n->num_quad);
    if (n->num_child > 0) {
        for (int i = 0; i < n->num_child; i++) {
            free_node(&n->child[i]);
        }
    }
    free(n->child);
}

void set_node_children(node_t *node, body_t *np_bodies) {
    int i, j;
    // c_quad keeps tracks of which bodies are in which quadrants by index
    // n_quad keeps track of the number of bodies in each of c_quad's quadrants
    // c_quad and n_quad are statically defined, hence they clean themselves
    int c_quad[8][node->num_bodies];
    int n_quad[8];

    for (i = 0; i < 8; i++) {
        n_quad[i]  = 0;
        for (j = 0; j < node->num_bodies; j++) {
            c_quad[i][j] = 0;
        }
    }

    for (i = 0; i < node->num_bodies; i++) {
        // derive quadrant q=0..7 from relative position on each axis of body
        // to center of current node
        int b_x = np_bodies[node->bodies[i]].px < node->px ? 0 : 1;
        int b_y = np_bodies[node->bodies[i]].py < node->py ? 0 : 1;
        int b_z = np_bodies[node->bodies[i]].pz < node->pz ? 0 : 1;
        // each quadrant is mapped to a number between 0..7
        int q = b_x*4 + b_y*2 + b_z;

        // keep track of all bodies in each quadrant q
        c_quad[q][n_quad[q]] = i;
        // keep track of number of bodies in each quadrant q
        n_quad[q]++;
    }

    for (i = 0; i < 8; i++) {
        if (n_quad[i] > 0) {
            // magic functions to determine new position based on quadrant and
            // length/size of node
            int m_x = i/4 == 0 ? -1 : 1;
            int m_y = (i/2)%2 == 0 ? -1 : 1;
            int m_z = i%2 == 0 ? -1 : 1;
            float new_px = node->px + m_x*0.25f*node->length;
            float new_py = node->py + m_y*0.25f*node->length;
            float new_pz = node->pz + m_z*0.25f*node->length;

            node->num_child++;
            node->child = (node_t*)realloc(node->child, node->num_child*sizeof(node_t));
            node_t *new_node = init_node(node->depth+1, node->max_depth, c_quad[i], n_quad[i]);
            node->child[node->num_child-1] = *new_node;
            free(new_node);
            set_node(&node->child[node->num_child-1], np_bodies, new_px, new_py, new_pz, 0.5f*node->length);
        }
    }
}

/* Recursively checks nodes and decides whether to calculate acceleration based
 * on given node or to check its children instead. Tail call recursion is
 * guaranteed with -O2 as long as check_node is called at the end and no lines
 * follow it. 99.8% of the time is spent in this subroutine.
 */
void check_node(node_t* node, body_t* body) {
    float rx = node->cx - body->px;
    float ry = node->cy - body->py;
    float rz = node->cz - body->pz;

    if (rx+ry+rz != 0.f) {
        float r = 1.f/sqrtf(rx*rx + ry*ry + rz*rz + E_SQR);
        if (node->num_child == 0) {
            float gmr = (G*node->tm)*(r*r*r);
            body->ax += gmr*rx;
            body->ay += gmr*ry;
            body->az += gmr*rz;
        }
        else if (node->length*r < DIST_THRES) {
            float gmr = (G*node->tm)*(r*r*r);
            body->ax += gmr*rx;
            body->ay += gmr*ry;
            body->az += gmr*rz;
        }
        else {// r != E | (ratio > DIST_THRES & node->num_child > 0)
            for (int i = 0; i < node->num_child; i++) {
                check_node(&node->child[i], body);
            }
        }
    }
}

void print_node_members(node_t *node, body_t *np_bodies) {
    for (int j = 0; j < node->depth; j++) {
        printf("\t");
    }
    if (node->num_child == 0) {
        printf("%d %d %.3f %.3f %.3f %.3f\n", node->id, node->num_bodies, node->cx, node->cy, node->cz, node->length);
        for (int i = 0; i < node->num_bodies; i++) {
            for (int j = 0; j < node->depth; j++) {
                printf("\t");
            }
            printf("  %.3f %.3f %.3f %.3f\n", np_bodies[node->bodies[i]].px, np_bodies[node->bodies[i]].py, np_bodies[node->bodies[i]].pz, node->tm);
        }
    } else {
        printf("%d %d %.3f %.3f %.3f %.3f\n", node->id, node->num_bodies, node->cx, node->cy, node->cz, node->length);
        for (int j = 0; j < node->depth; j++) {
            printf("\t");
        }
        printf("%d %.3f %.3f %.3f %.3f\n", node->depth, node->cx, node->cy, node->cz, node->tm);
        for (int i = 0; i < node->num_child; i++) {
            print_node_members(&node->child[i], np_bodies);
        }
    }
}

nbodysys_t *init_rand_nbodysys(int n, float max_p, float max_v, float max_m) {
    int i;

    nbodysys_t *temp = (nbodysys_t*)malloc(sizeof(nbodysys_t));
    temp->bodies = (body_t*)malloc(n*sizeof(body_t));
    temp->maxp = max_p; temp->maxv = max_v; temp->maxm = max_m;
    temp->num_bodies = n;

    srand(time(NULL));
    for (i = 0; i < n; i++) {
        body_t *temp_body = init_rand_body(max_p, max_v, max_m);
        temp->bodies[i] = *temp_body;
        free(temp_body);
    }

    return temp;
}

nbodysys_t *copy_bodysys(nbodysys_t *src) {
    int n = src->num_bodies;
    nbodysys_t *dest = init_rand_nbodysys(n, src->maxp, src->maxv, src->maxm);

    dest->bodies = (body_t*)realloc(dest->bodies, n*sizeof(body_t));
    memcpy(dest->bodies, src->bodies, n*sizeof(body_t));

    return dest;
}

void free_nbodysys(nbodysys_t *nb) {
    free(nb->bodies);
    free(nb);
}

void print_nbodysys(nbodysys_t *nb) {
    printf("\n");
    for (int i = 0; i < nb->num_bodies; i++) {
        printf("%.4f\t%.4f\t%.4f\n", nb->bodies[i].px, nb->bodies[i].py, nb->bodies[i].pz);
    }
}

void brute(nbodysys_t *nb, int iters, del_t time) {
    int iter, i, j;
    for (iter = 0; iter < iters*time; iter += time) {
#pragma omp parallel for private(j)
        for (i = 0; i < nb->num_bodies; i++) {
#pragma omp parallel for
            for (j = 0; j < nb->num_bodies; j++) {
                if (i != j) {
                    float rx = nb->bodies[j].px - nb->bodies[i].px;
                    float ry = nb->bodies[j].py - nb->bodies[i].py;
                    float rz = nb->bodies[j].pz - nb->bodies[i].pz;
                    float r = 1.f/sqrtf(rx*rx + ry*ry + rz*rz + E_SQR);
                    float gmr = (G*nb->bodies[j].m)*(r*r*r);
                    nb->bodies[i].ax += gmr*rx;
                    nb->bodies[i].ay += gmr*ry;
                    nb->bodies[i].az += gmr*rz;
                }
            }
            update_p(&nb->bodies[i], time);
        }
    }
}

void all_seq(nbodysys_t *nb, int iters, del_t time) {
    int n = nb->num_bodies;
    int iter, i, j;
    float ht = 0.5f*time;
    float *p = (float*)malloc(n*3*sizeof(float));
    float *v = (float*)malloc(n*3*sizeof(float));
    float *a = (float*)malloc(n*3*sizeof(float));
    float *gm = (float*)malloc(n*sizeof(float));

    for (i = 0; i < nb->num_bodies; i++) {
        int in = i*3;
        p[in] = nb->bodies[i].px;
        p[in+1] = nb->bodies[i].py;
        p[in+2] = nb->bodies[i].pz;
        v[in] = nb->bodies[i].vx;
        v[in+1] = nb->bodies[i].vy;
        v[in+2] = nb->bodies[i].vz;
        a[in] = nb->bodies[i].ax;
        a[in+1] = nb->bodies[i].ay;
        a[in+2] = nb->bodies[i].az;
        gm[i] = G*nb->bodies[i].m;
    }

    for (iter = 0; iter < iters*time; iter += time) {
#pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            int x = i*3;
            int y = i*3+1;
            int z = i*3+2;
#pragma omp parallel for
            for (j = 0; j < n; j++) {
                if (i != j) {
                    float rx = p[j*3] - p[x];
                    float ry = p[j*3+1] - p[y];
                    float rz = p[j*3+2] - p[z];
                    float r = 1.f/sqrtf(rx*rx + ry*ry + rz*rz + E_SQR);
                    float gmr = gm[j]*(r*r*r);
                    a[x] += gmr*rx;
                    a[y] += gmr*ry;
                    a[z] += gmr*rz;
                }
            }
            v[x] += a[x]*time;
            v[y] += a[y]*time;
            v[z] += a[z]*time;
            p[x] += v[x]*ht;
            p[y] += v[y]*ht;
            p[z] += v[z]*ht;
        }
    }

    for (i = 0; i < n; i++) {
        int in = i*3;
        nb->bodies[i].px = p[in];
        nb->bodies[i].py = p[in+1];
        nb->bodies[i].pz = p[in+2];
        nb->bodies[i].vx = v[in];
        nb->bodies[i].vy = v[in+1];
        nb->bodies[i].vz = v[in+2];
        nb->bodies[i].ax = a[in];
        nb->bodies[i].ay = a[in+1];
        nb->bodies[i].az = a[in+2];
    }

    free(p);
    free(v);
    free(a);
    free(gm);
}

void barnes(nbodysys_t *nb, int iters, del_t time) {
    int iter, i;
    int n = nb->num_bodies;
    node_t *root = 0;
    int *root_bodies = (int*)malloc(n*sizeof(int));
    float max;

    for (i = 0; i < n; i++) root_bodies[i] = i;

    for (iter = 0; iter < iters*time; iter += time) {
        max = 0.f;
        node_id = 1;
        for (i = 0; i < n; i++) {
            max = fabs(nb->bodies[i].px) > max ? fabs(nb->bodies[i].px) : max;
            max = fabs(nb->bodies[i].py) > max ? fabs(nb->bodies[i].py) : max;
            max = fabs(nb->bodies[i].pz) > max ? fabs(nb->bodies[i].pz) : max;
        }
        root = init_node(0, 16, root_bodies, n);
        set_node(root, nb->bodies, 0.f, 0.f, 0.f, max*2.);

#pragma omp parallel for
        for (i = 0; i < n; i++) {
            check_node(root, &nb->bodies[i]);
            update_p(&nb->bodies[i], time);
        }
        if (debug == 2) {
            print_node_members(root, nb->bodies);
        }
        free_node(root);
    }
    free(root_bodies);
    free(root);
}


