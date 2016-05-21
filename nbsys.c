#include "nbsys.h"
#include <assert.h>

#define G 1.f     // universal gravitational constant should be 6.67408e-11
#define E_SQR 1.f // softening factor
#define DEBUG 0
#define DIST_THRES 0.5f

// debugging
static int node_id;
static int debug;

body_t *init_rand_body(float max_p, float max_v, float max_m) {
    body_t *temp = (body_t*)malloc(sizeof(body_t));
    temp->px = max_p*((float) rand() / (float) RAND_MAX) - 0.5*max_p;
    temp->py = max_p*((float) rand() / (float) RAND_MAX) - 0.5*max_p;
    temp->pz = max_p*((float) rand() / (float) RAND_MAX) - 0.5*max_p;
    temp->vx = max_v*((float) rand() / (float) RAND_MAX) - 0.5*max_v;
    temp->vy = max_v*((float) rand() / (float) RAND_MAX) - 0.5*max_v;
    temp->vz = max_v*((float) rand() / (float) RAND_MAX) - 0.5*max_v;
    temp->ax = 0;
    temp->ay = 0;
    temp->az = 0;
    temp->m = max_m*((float) rand() / (float) RAND_MAX);

    return temp;
}

inline void update_p(body_t* b, del_t t) {
    float ht = 0.5*t;
    float ux = b->vx;
    float uy = b->vy;
    float uz = b->vz;
    b->vx += b->ax*t;
    b->vy += b->ay*t;
    b->vz += b->az*t;
    b->px += (b->vx - ux)*ht;
    b->py += (b->vy - uy)*ht;
    b->pz += (b->vz - uz)*ht;
}

void free_body(body_t *body) {
    free(body);
}

node_t *init_node(int dep, int max_dep) {
    node_t *temp = (node_t*)malloc(sizeof(node_t));
    temp->depth = dep;
    temp->max_depth = max_dep;
    temp->id = node_id++;
    return temp;
}

void free_node(node_t *n) {
    free(n->bodies);
    if (n->num_child > 0) {
        for (int i = 0; i < n->num_child; i++) {
            free_node(&n->child[i]);
        }
    }
    free(n->child);
}

void set_node_members(node_t *node, body_t *p_bodies, float p_x, float p_y, float p_z, float p_length, int p_nbodies) {
    if (p_nbodies > 1) {
        int i;
        node->px = p_x;
        node->py = p_y;
        node->pz = p_z;
        node->length = p_length;
        node->num_bodies = p_nbodies;
        node->bodies = (body_t*)malloc(p_nbodies*sizeof(body_t));
        memcpy(node->bodies, p_bodies, p_nbodies*sizeof(body_t));

        float t_tm = 0;
        float t_cx = 0;
        float t_cy = 0;
        float t_cz = 0;
        for (i = 0; i < p_nbodies; i++) {
            float t_m = node->bodies[i].m;
            t_cx += p_x*t_m;
            t_cy += p_y*t_m;
            t_cz += p_z*t_m;
            t_tm += t_m;
        }
        node->tm = t_tm;
        float invtm = 1.f/t_tm;
        node->cx = t_cx*invtm;
        node->cy = t_cy*invtm;
        node->cz = t_cz*invtm;

        node->num_child = 0;
        node->child = 0;
        if (node->depth < node->max_depth) {
            set_node_children(node);
        }
    }
    else if (p_nbodies == 1) {
        node->px = p_x;
        node->py = p_y;
        node->pz = p_z;
        node->length = p_length;
        node->num_bodies = 1;
        node->bodies = (body_t*)malloc(sizeof(body_t));
        memcpy(node->bodies, &p_bodies[0], sizeof(body_t));
        node->cx = node->bodies[0].px;
        node->cy = node->bodies[0].py;
        node->cz = node->bodies[0].pz;
        node->tm = p_bodies[0].m;
        node->num_child = 0;
        node->child = 0;
    }
}

void set_node_children(node_t *node) {
    int i, j;
    body_t **quad = (body_t**)malloc(8*sizeof(body_t*));
    int c_quad[8][node->num_bodies];
    int n_quad[8];

    for (i = 0; i < 8; i++) {
        n_quad[i]  = 0;
        quad[i] = 0;
        for (j = 0; j < node->num_bodies; j++) {
            c_quad[i][j] = 0;
        }
    }

//#pragma omp parallel for
    for (i = 0; i < node->num_bodies; i++) {
        // derive quadrant q=0..7 from relative position of body to center on each axis
        int b_x = node->bodies[i].px < node->px ? 0 : 1;
        int b_y = node->bodies[i].py < node->py ? 0 : 1;
        int b_z = node->bodies[i].pz < node->pz ? 0 : 1;
        int q = b_x*4 + b_y*2 + b_z;

//#pragma omp critical
        {
            // keep track of number of bodies in each quadrant q
            n_quad[q]++;
            // keep track of which quadrant each body belongs to
            c_quad[q][n_quad[q]-1] = i;
        }
    }

    for (i = 0; i < 8; i++) {
        if (n_quad[i] > 0) {
            // make new list(s) of each quadrant and its inhabitants
            quad[i] = (body_t*)realloc(quad[i], (n_quad[i])*sizeof(body_t));
            for (j = 0; j < n_quad[i]; j++) {
                quad[i][j] = node->bodies[c_quad[i][j]];
            }

            // magic functions to determine new position based on quadrant and length/size of node
            int m_x = i/4 == 0 ? -1 : 1;
            int m_y = (i/2)%2 == 0 ? -1 : 1;
            int m_z = i%2 == 0 ? -1 : 1;
            float new_px = node->px + m_x*0.25f*node->length;
            float new_py = node->py + m_y*0.25f*node->length;
            float new_pz = node->pz + m_z*0.25f*node->length;

            node_t *t_node = init_node(node->depth+1, node->max_depth);
            set_node_members(t_node, quad[i], new_px, new_py, new_pz, node->length*0.5, n_quad[i]);
            node->num_child++;
            node->child = (node_t*)realloc(node->child, node->num_child*sizeof(node_t));
            node->child[node->num_child-1] = *t_node;

            free(t_node);
        }
        free(quad[i]);
    }
    free(quad);
}

void check_node(node_t* node, body_t* body) {
    if (node->num_bodies > 0) {
        float rx = node->cx - body->px;
        float ry = node->cy - body->py;
        float rz = node->cz - body->pz;
        float dist = rx*rx + ry+ry + rz*rz;
        float ratio = node->length/sqrtf(dist);

        if (dist <= 0.f) {}
        else if ((ratio < DIST_THRES) | (node->num_child == 0)) {
            //assert((node->bodies[0].px == node->cx) | (node->num_bodies > 1));
            float r = dist + E_SQR;
            float gmr = G*node->tm/sqrtf(r*r*r);
            body->ax += gmr*rx;
            body->ay += gmr*ry;
            body->az += gmr*rz;
        }
        else {
            for (int i = 0; i < node->num_child; i++) {
                check_node(&node->child[i], body);
            }
        }
    }
}

void print_node_members(node_t *node) {
    for (int j = 0; j < node->depth; j++) {
        printf("\t");
    }
    if (node->num_child == 0) {
        printf("%d %d %.3f %.3f %.3f %.3f\n", node->id, node->num_bodies, node->px, node->py, node->pz, node->length);
        for (int i = 0; i < node->num_bodies; i++) {
            for (int j = 0; j < node->depth; j++) {
                printf("\t");
            }
            printf("  %.3f %.3f %.3f\n", node->bodies[i].px, node->bodies[i].py, node->bodies[i].pz);
        }
    } else {
        printf("%d %d %.3f %.3f %.3f %.3f\n", node->id, node->num_bodies, node->px, node->py, node->pz, node->length);
        for (int i = 0; i < node->num_child; i++) {
            print_node_members(&node->child[i]);
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
    for (int i = 0; i < nb->num_bodies; i++) {
        printf("%.3f\t%.3f\t%.3f\n", nb->bodies[i].px, nb->bodies[i].py, nb->bodies[i].pz);
    }
    printf("\n");
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

    for (iter = 0; iter < iters; iter++) {
#pragma omp parallel for
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
                    float r = rx*rx + ry*ry + rz*rz + E_SQR;
                    float gmr = gm[j]/sqrtf(r*r*r);
                    a[x] += gmr*rx;
                    a[y] += gmr*ry;
                    a[z] += gmr*rz;
                }
            }
            float ux = v[x];
            float uy = v[y];
            float uz = v[z];
            v[x] += a[x]*time;
            v[y] += a[y]*time;
            v[z] += a[z]*time;
            p[x] += (v[x] - ux)*ht;
            p[y] += (v[y] - uy)*ht;
            p[z] += (v[z] - uz)*ht;
        }
    }

    for (i = 0; i < nb->num_bodies; i++) {
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
    int i, iter;
    int n = nb->num_bodies;
    node_t *root = init_node(0, 16);
    float max;
    omp_set_nested(1);

    for (iter = 0; iter < iters*time; iter += time) {
        max = 0;
        node_id = 1;
        for (i = 0; i < n; i++) {
            max = fabs(nb->bodies[i].px) > max ? fabs(nb->bodies[i].px) : max;
            max = fabs(nb->bodies[i].py) > max ? fabs(nb->bodies[i].py) : max;
            max = fabs(nb->bodies[i].pz) > max ? fabs(nb->bodies[i].pz) : max;
        }
        set_node_members(root, nb->bodies, 0, 0, 0, max*2, nb->num_bodies);
#pragma omp parallel for
        for (i = 0; i < nb->num_bodies; i++) {
            check_node(root, &nb->bodies[i]);
            update_p(&nb->bodies[i], time);
        }
        if (debug == 2) {
            print_node_members(root);
        }

        free_node(root);
    }
    free(root);
}

