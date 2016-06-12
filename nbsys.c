#include "nbsys.h"

#define G 0.1f // universal gravitational constant should be 6.67408e-11
#define MAX_P 100.f
#define MAX_V 10.f
#define MAX_M 100.f
#define E 50.f
#define E_SQR 2500.f // softening factor
#define DEBUG 2
#define DIST_THRES 0.25f

// debugging
static uint32_t node_id;
static uint8_t debug;
//static uint32_t node_counter;

// carmack method for fast inverse square root (basically newton's method)
static inline float rsqrt(float x) {
    float hx = 0.5f*x;
    uint32_t i = *(uint32_t*)&x;
    i = 0x5f3759df - (i>>1);
    x = *(float*)&i;
    x *= (1.5f - hx*x*x);
    return x;
}

// initialize a random body, specifying the maximum values for each attribute,
// and returning a pointer to it.
body_t *init_rand_body(float max_p, float max_v, float max_m) {
    body_t *temp = (body_t*)malloc(sizeof(body_t));
    temp->px = max_p*(((float)rand() / (float)RAND_MAX) - 0.5f);
    temp->py = max_p*(((float)rand() / (float)RAND_MAX) - 0.5f);
    temp->pz = max_p*(((float)rand() / (float)RAND_MAX) - 0.5f);
    temp->vx = max_v*(((float)rand() / (float)RAND_MAX) - 0.5f);
    temp->vy = max_v*(((float)rand() / (float)RAND_MAX) - 0.5f);
    temp->vz = max_v*(((float)rand() / (float)RAND_MAX) - 0.5f);
    temp->ax = 0.f;
    temp->ay = 0.f;
    temp->az = 0.f;
    temp->m = max_m*((float) rand() / (float) RAND_MAX);

    return temp;
}

// updates the velocity and position of a body from its acceleration.
inline void update_body(body_t* b, del_t t) {
    float ht = 0.5f*t;
    b->vx += b->ax*t;
    b->vy += b->ay*t;
    b->vz += b->az*t;
    b->px += b->vx*ht;
    b->py += b->vy*ht;
    b->pz += b->vz*ht;
}

// function that initializes a node and returns a pointer to it. this function
// is kept separate from the recursive call to ensure the pointer can be freed
// before the call takes place.
// recursively divide bodies uint32_to quadrants and create new child nodes from
// each of these quadrants. stops when maximum node depth is reached or when
// exactly one body is found in quadrant. quadrants with no bodies are not
// initialized as nodes.
void init_node(node_t *node, uint8_t dep, uint8_t max_dep, body_t *np_bodies, uint32_t* p_bodies, uint32_t p_nbodies, float p_x, float p_y, float p_z, float p_length) {
    uint32_t i;
    uint32_t n = p_nbodies;

    node->id = node_id++;
    node->depth = dep;
    node->max_depth = max_dep;

    node->num_bodies = n;
    node->bodies = (uint32_t*)malloc(p_nbodies*sizeof(uint32_t));
    memcpy(node->bodies, p_bodies, p_nbodies*sizeof(uint32_t));

    memset(&node->quad, 0, sizeof(node->quad));
    memset(&node->num_quad, 0, sizeof(node->num_quad));
    memset(&node->child, 0, sizeof(node->child));

    node->px = p_x;
    node->py = p_y;
    node->pz = p_z;
    node->length = p_length;
    node->num_child = 0;

    if (n == 1) {
        uint32_t b = node->bodies[0];
        node->tm = np_bodies[b].m;
        node->cx = np_bodies[b].px;
        node->cy = np_bodies[b].py;
        node->cz = np_bodies[b].pz;
    }
    else if (n > 1) {
        // quad keeps tracks of which bodies are in which quadrants by
        // index. num_quad keeps track of the number of bodies in each
        // of quad's quadrants.
        uint32_t *quad_data = (uint32_t*)malloc(8*n*sizeof(uint32_t));
        for (uint8_t i = 0; i < 8; i++) {
            node->quad[i] = &quad_data[n*i];
        }
        // calculate center of mass for each axis and total mass
        float t_tm = 0.f;
        float t_cx = 0.f;
        float t_cy = 0.f;
        float t_cz = 0.f;

        for (i = 0; i < n; i++) {
            uint32_t b = node->bodies[i];

            float t_m = np_bodies[b].m;
            t_cx += np_bodies[b].px*t_m;
            t_cy += np_bodies[b].py*t_m;
            t_cz += np_bodies[b].pz*t_m;
            t_tm += t_m;

            // derive quadrant q=0..7 from relative position on each axis
            // of body to center of current node
            uint32_t b_x = np_bodies[b].px > node->px;
            uint32_t b_y = np_bodies[b].py > node->py;
            uint32_t b_z = np_bodies[b].pz > node->pz;
            // each quadrant is mapped to a number between 0..7
            uint32_t q = b_x*4 + b_y*2 + b_z;

            // keep track of all bodies in each quadrant q
            node->quad[q][node->num_quad[q]] = b;
            // keep track of number of bodies in each quadrant q
            node->num_quad[q]++;
        }
        node->tm = t_tm;
        float invtm = 1.f/t_tm;
        node->cx = t_cx*invtm;
        node->cy = t_cy*invtm;
        node->cz = t_cz*invtm;
#pragma omp parallel for
        for (uint8_t i = 0; i < 8; i++) {
            if (node->num_quad[i] > 0) {
                // dynamically add child
#pragma omp atomic
                node->num_child++;
                // magic functions to determine new position based on
                // quadrant and length/size of node
                int8_t m_x = i/4 == 0 ? -1 : 1;
                int8_t m_y = (i/2)%2 == 0 ? -1 : 1;
                int8_t m_z = i%2 == 0 ? -1 : 1;
                float new_px = node->px + m_x*0.25f*node->length;
                float new_py = node->py + m_y*0.25f*node->length;
                float new_pz = node->pz + m_z*0.25f*node->length;
                // initialize new node and push uint32_to child array
                // the recursive call *MUST BE PLACED LAST* to ensure the
                // tail call optimization
                node->child[i] = (node_t*)malloc(sizeof(node_t));
                init_node(node->child[i], node->depth+1, node->max_depth, np_bodies, node->quad[i], node->num_quad[i], new_px, new_py, new_pz, 0.5f*node->length);
            }
        }
    }
}

// recursively free nodes
void free_node(node_t *n) {
    free(n->bodies);
    if (n->num_child > 0) {
        free(n->quad[0]);
#pragma omp parallel for
        for (uint8_t i = 0; i < 8; i++) {
            if (n->child[i] != 0) {
                free_node(n->child[i]);
                free(n->child[i]);
            }
        }
    }
}

/* Recursively checks nodes and decides whether to calculate acceleration based
 * on given node or to check its children instead. Tail call recursion is
 * guaranteed with -O2 as long as check_node is called at the end and no lines
 * follow it. >90% of the time is spent in this subroutine.
 */
void check_node(node_t* node, body_t* body) {
    float rx = node->cx - body->px;
    float ry = node->cy - body->py;
    float rz = node->cz - body->pz;
    float r = rx*rx + ry*ry + rz*rz;
    if (r > 0.f) {
        r = 1.f/sqrtf(r + E_SQR);
        if (node->length*r > DIST_THRES && node->num_child != 0) {
            for (uint8_t i = 0; i < 8; i++) {
                if (node->child[i] != 0)
                    check_node(node->child[i], body);
            }
        }
        else {// r != E | (ratio > DIST_THRES & node->num_child > 0)
            float gmr = G*node->tm*r*r*r;
            body->ax += gmr*rx;
            body->ay += gmr*ry;
            body->az += gmr*rz;
        }
    }
}

void check_node_ordered(node_t* node, body_t* body, uint32_t *order, uint32_t *size) {
    double rx = node->cx - body->px;
    double ry = node->cy - body->py;
    double rz = node->cz - body->pz;
    double r = rx*rx + ry*ry + rz*rz;
    if (r > 0.) {
        r = 1.f/sqrtf(r + E_SQR);
        if (node->length*r > DIST_THRES && node->num_child != 0) {
            for (uint8_t i = 0; i < 8; i++) {
                if (node->child[i] != 0)
                    check_node_ordered(node->child[i], body, order, size);
            }
        }
        else {// r != E | (ratio > DIST_THRES & node->num_child > 0)
            order[*size] = node->bodies[0];
            *size += 1;
            double gmr = G*node->tm*r*r*r;
            body->ax += gmr*rx;
            body->ay += gmr*ry;
            body->az += gmr*rz;
        }
    }
}

void print_node(node_t *node, body_t *np_bodies) {
    for (uint32_t j = 0; j < node->depth; j++) {
        printf("\t");
    }
    if (node->num_child == 0) {
        printf("%d %d %.3f %.3f %.3f %.3f\n", node->id, node->num_bodies, node->px, node->py, node->pz, node->length);
        for (uint32_t i = 0; i < node->num_bodies; i++) {
            for (uint32_t j = 0; j < node->depth; j++) {
                printf("\t");
            }
            printf("  %.3f %.3f %.3f %.3f\n", np_bodies[node->bodies[i]].px, np_bodies[node->bodies[i]].py, np_bodies[node->bodies[i]].pz, node->tm);
        }
    } else {
        printf("%d %d %.3f %.3f %.3f %.3f\n", node->id, node->num_bodies, node->cx, node->cy, node->cz, node->length);
        for (uint32_t j = 0; j < node->depth; j++) {
            printf("\t");
        }
        printf("%d %.3f %.3f %.3f %.3f\n", node->depth, node->cx, node->cy, node->cz, node->tm);
        for (uint32_t i = 0; i < 8; i++) {
            if (node->child[i] != 0)
                print_node(node->child[i], np_bodies);
        }
    }
}

nbodysys_t *init_rand_nbodysys(uint32_t n, float max_p, float max_v, float max_m) {
    nbodysys_t *temp = (nbodysys_t*)malloc(sizeof(nbodysys_t));
    temp->bodies = (body_t*)malloc(n*sizeof(body_t));
    temp->maxp = max_p; temp->maxv = max_v; temp->maxm = max_m;
    temp->num_bodies = n;

    srand(time(NULL));
    for (uint32_t i = 0; i < n; i++) {
        body_t *temp_body = init_rand_body(max_p, max_v, max_m);
        temp->bodies[i] = *temp_body;
        free(temp_body);
    }
    return temp;
}

nbodysys_t *copy_nbodysys(nbodysys_t *src) {
    uint32_t n = src->num_bodies;
    nbodysys_t *dest = (nbodysys_t*)malloc(sizeof(nbodysys_t));
    dest->num_bodies = n;
    dest->maxp = src->maxp;
    dest->maxv = src->maxv;
    dest->maxm = src->maxm;
    dest->bodies = (body_t*)malloc(n*sizeof(body_t));
    memcpy(dest->bodies, src->bodies, n*sizeof(body_t));

    return dest;
}

void free_nbodysys(nbodysys_t *nb) {
    free(nb->bodies);
    free(nb);
}

void print_nbodysys(nbodysys_t *nb) {
    printf("\n");
    for (uint32_t i = 0; i < nb->num_bodies; i+=256) {
        printf("%.9f\t%.9f\t%.9f\n", nb->bodies[i].px, nb->bodies[i].py, nb->bodies[i].pz);
    }
}

void brute(nbodysys_t *nb, uint32_t iters, del_t time) {
    uint32_t iter, i, j, n;
    n = nb->num_bodies;
    for (iter = 0; iter < iters*time; iter += time) {
#pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
#pragma omp parallel for
            for (j = 0; j < n; j++) {
                if (i != j) {
                    float rx = nb->bodies[j].px - nb->bodies[i].px;
                    float ry = nb->bodies[j].py - nb->bodies[i].py;
                    float rz = nb->bodies[j].pz - nb->bodies[i].pz;
                    float r = 1.f/sqrtf(rx*rx + ry*ry + rz*rz + E_SQR);
                    float gmr = G*nb->bodies[j].m*r*r*r;
                    nb->bodies[i].ax += gmr*rx;
                    nb->bodies[i].ay += gmr*ry;
                    nb->bodies[i].az += gmr*rz;
                }
            }
        }
#pragma omp parallel for
        for (i = 0; i < n; i++)
            update_body(&nb->bodies[i], time);
    }
}

void brute_ordered(nbodysys_t *nb, uint32_t iters, del_t time, uint32_t **check_order) {
    uint32_t iter, i, j, n;
    n = nb->num_bodies;
    for (iter = 0; iter < iters*time; iter += time) {
#pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
#pragma omp parallel for
            for (j = 0; j < n-1; j++) {
                uint32_t b = check_order[i][j];
                double rx = nb->bodies[b].px - nb->bodies[i].px;
                double ry = nb->bodies[b].py - nb->bodies[i].py;
                double rz = nb->bodies[b].pz - nb->bodies[i].pz;
                double r = 1.f/sqrtf(rx*rx + ry*ry + rz*rz + E_SQR);
                double gmr = G*nb->bodies[b].m*r*r*r;
                nb->bodies[i].ax += gmr*rx;
                nb->bodies[i].ay += gmr*ry;
                nb->bodies[i].az += gmr*rz;
            }
        }
#pragma omp parallel for
        for (i = 0; i < n; i++)
            update_body(&nb->bodies[i], time);
    }
}

void barnes(nbodysys_t *nb, uint32_t iters, del_t time) {
    uint32_t iter, i;
    uint32_t n = nb->num_bodies;
    node_t *root = 0;
    uint32_t *root_bodies = (uint32_t*)malloc(n*sizeof(uint32_t));
    float max;
    //omp_set_nested(1);

    for (i = 0; i < n; i++) root_bodies[i] = i;

    for (iter = 0; iter < iters*time; iter += time) {
        max = 0.f;
        node_id = 1;
        for (i = 0; i < n; i++) {
            max = fabs(nb->bodies[i].px) > max ? fabs(nb->bodies[i].px) : max;
            max = fabs(nb->bodies[i].py) > max ? fabs(nb->bodies[i].py) : max;
            max = fabs(nb->bodies[i].pz) > max ? fabs(nb->bodies[i].pz) : max;
        }
        root = (node_t*)malloc(sizeof(node_t));
        init_node(root, 0, 16, nb->bodies, root_bodies, n, 0.f, 0.f, 0.f, max*2.f);

#pragma omp parallel for
        for (i = 0; i < n; i++) {
            check_node(root, &nb->bodies[i]);
        }
#pragma omp parallel for
        for (i = 0; i < n; i++)
            update_body(&nb->bodies[i], time);
        if (debug == 2) {
            print_node(root, nb->bodies);
        }
        free_node(root);
    }
    free(root_bodies);
    free(root);
}

void barnes_ordered(nbodysys_t *nb, uint32_t iters, del_t time, uint32_t **check_order) {
    uint32_t iter, i;
    uint32_t n = nb->num_bodies;
    node_t *root = 0;
    uint32_t *root_bodies = (uint32_t*)malloc(n*sizeof(uint32_t));
    float max;
    //omp_set_nested(1);

    for (i = 0; i < n; i++) root_bodies[i] = i;

    for (iter = 0; iter < iters*time; iter += time) {
        max = 0.f;
        node_id = 1;
        for (i = 0; i < n; i++) {
            max = fabs(nb->bodies[i].px) > max ? fabs(nb->bodies[i].px) : max;
            max = fabs(nb->bodies[i].py) > max ? fabs(nb->bodies[i].py) : max;
            max = fabs(nb->bodies[i].pz) > max ? fabs(nb->bodies[i].pz) : max;
        }
        root = (node_t*)malloc(sizeof(node_t));
        init_node(root, 0, 16, nb->bodies, root_bodies, n, 0.f, 0.f, 0.f, max*2.f);

#pragma omp parallel for
        for (i = 0; i < n; i++) {
            uint32_t size = 0;
            check_node_ordered(root, &nb->bodies[i], check_order[i], &size);
        }
#pragma omp parallel for
        for (i = 0; i < n; i++)
            update_body(&nb->bodies[i], time);
        if (debug == 2) {
            print_node(root, nb->bodies);
        }
        free_node(root);
    }
    free(root_bodies);
    free(root);
}

