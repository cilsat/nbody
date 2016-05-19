#include "nbsys.h"
#include <assert.h>

#define G 6.67408e-11   // universal gravitational constant
#define E               // softening factor

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

void free_body(body_t *body) {
    free(body);
}

node_t *init_node(int dep, int max_dep) {
    node_t *temp = (node_t *) malloc(sizeof(node_t));
    temp->depth = dep;
    temp->max_depth = max_dep;
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
    int i;

    node->px = p_x;
    node->py = p_y;
    node->pz = p_z;
    node->length = p_length;
    node->num_bodies = p_nbodies;
    if (p_nbodies > 0) {
        node->bodies = (body_t *) malloc(p_nbodies*sizeof(body_t));
        memcpy(node->bodies, p_bodies, p_nbodies*sizeof(body_t));
    }
    else {
        node->bodies = 0;
    }

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
    float invtm = p_nbodies > 0 ? 1.f/t_tm : 0.f;
    node->cx = t_cx*invtm;
    node->cy = t_cy*invtm;
    node->cz = t_cz*invtm;

    node->num_child = 0;
    node->child = 0;
    if ((p_nbodies > 1) && (node->depth < node->max_depth)) {
        set_node_children(node);
    }
}

void set_node_children(node_t *node) {
    int i, j, q;
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

    for (i = 0; i < node->num_bodies; i++) {
        // derive quadrant q=0..7 from relative position of body to center on each axis
        int b_x = node->bodies[i].px < node->px ? 0 : 1;
        int b_y = node->bodies[i].py < node->py ? 0 : 1;
        int b_z = node->bodies[i].pz < node->pz ? 0 : 1;
        q = b_x*4 + b_y*2 + b_z;

        // keep track of which quadrant each body belongs to
        c_quad[q][n_quad[q]] = i;
        // keep track of number of bodies in each quadrant q
        n_quad[q]++;
    }

    for (i = 0; i < 8; i++) {
        // make new list(s) of each quadrant and its inhabitants
        quad[i] = (body_t*)realloc(quad[i], (n_quad[i])*sizeof(body_t));
        for (j = 0; j < n_quad[i]; j++) {
            quad[i][j] = node->bodies[c_quad[i][j]];
        }

        // magic functions to determine new position based on quadrant and length/size of node
        int m_x = i/4 > 0 ? -1 : 1;
        int m_y = (i/2)%2 > 0 ? -1 : 1;
        int m_z = i%2 > 0 ? -1 : 1;
        float new_px = node->px + m_x*0.25f*node->length;
        float new_py = node->py + m_y*0.25f*node->length;
        float new_pz = node->pz + m_z*0.25f*node->length;

        node_t *t_node = init_node(node->depth+1, node->max_depth);
        set_node_members(t_node, quad[i], new_px, new_py, new_pz, node->length*0.5, n_quad[i]);
        node->num_child++;
        node->child = (node_t*)realloc(node->child, node->num_child*sizeof(node_t));
        node->child[node->num_child-1] = *t_node;

        free(t_node);
        free(quad[i]);
    }
    free(quad);
}

void print_node_members(node_t *node) {
    for (int i = 0; i < node->num_child; i++) {
        printf("%d %d %.3f %.3f %.3f %.3f\n", node->num_bodies, node->depth, node->px, node->py, node->pz, node->length);
        print_node_members(&node->child[i]);
    }
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
        body_t *temp_body = init_rand_body(max_p, max_v, max_m);
        temp->bodies[i] = *temp_body;
        free(temp_body);

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

void free_nbodysys(nbodysys_t *nb) {
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
    float *p = (float *) malloc(n*3*sizeof(float));
    float *v = (float *) malloc(n*3*sizeof(float));
    memcpy(p, nb->p, n*3*sizeof(float));
    memcpy(v, nb->v, n*3*sizeof(float));
    float gt = 1.f*time;
    float ht = 0.5f*time;

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        float fx = 0;
        float fy = 0;
        float fz = 0;
        int x = i*3;
        int y = i*3+1;
        int z = i*3+2;
#pragma omp parallel for
        for (j = 0; j < n; j++) {
            if (i != j) {
                float rx = p[j*3] - p[x];
                float ry = p[j*3+1] - p[y];
                float rz = p[j*3+2] - p[z];
                float r = 1.f/sqrt(rx*rx + ry*ry + rz*rz);
                float mr = nb->m[j]*r*r*r;
                fx += mr*rx;
                fy += mr*ry;
                fz += mr*rz;
            }
        }
        float ux = v[x];
        float uy = v[y];
        float uz = v[z];
        v[x] += fx*gt;
        v[y] += fy*gt;
        v[z] += fz*gt;
        p[x] += (v[x] - ux)*ht;
        p[y] += (v[y] - uy)*ht;
        p[z] += (v[z] - uz)*ht;
    }

    free(p);
    free(v);
}

void barnes(nbodysys_t *nb, del_t time) {
    int i;
    node_t *root = init_node(0, 16);

    float max = 0;
    for (i = 0; i < nb->num_bodies; i++) {
        max = nb->p[i*3] > max ? nb->p[i*3] : max; 
        max = nb->p[i*3+1] > max ? nb->p[i*3+1] : max; 
        max = nb->p[i*3+2] > max ? nb->p[i*3+2] : max; 
    }
    set_node_members(root, nb->bodies, 0, 0, 0, max*2, nb->num_bodies);

    print_node_members(root);
    free_node(root);
    free(root);
}
