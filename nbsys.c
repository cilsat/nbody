#include "nbsys.h"

#define TIME 1.f
#define MAX_P 1000.f
#define MAX_V MAX_P*0.01
#define MAX_M 1000.f
#define LEN_MAX MAX_P*0.5f
#define LEN_MIN MAX_P*-0.5f
#define G 0.01f // universal gravitational constant should be 6.67408e-11
#define E_SQR 100.f // softening factor
#define DEBUG 0
#define DIST_THRES 0.5f
#define DEPTH 3
#define DBG_DISPLAY 32

// debugging
static uint32_t node_id = 0;
//static uint32_t node_check = 0;
static uint8_t debug = DEBUG;
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

static inline float rand_range(float max) {
    return max*(((float)rand() / (float)RAND_MAX) - 0.5f);
}

// initialize a random body, specifying the maximum values for each attribute,
// and returning a pointer to it.
body_t *init_rand_body(float max_p, float max_v, float max_m) {
    body_t *temp = malloc(sizeof(body_t));
    temp->px = rand_range(max_p);
    temp->py = rand_range(max_p);
    temp->pz = rand_range(max_p);
    temp->vx = rand_range(max_v);
    temp->vy = rand_range(max_v);
    temp->vz = rand_range(max_v);
    temp->ax = 0.f; temp->ay = 0.f; temp->az = 0.f;
    temp->m = max_m*((float) rand() / (float) RAND_MAX);

    return temp;
}

// updates the velocity and position of a body from its acceleration.
inline void update_body(body_t* b, del_t t) {
    b->vx += b->ax*t;
    b->vy += b->ay*t;
    b->vz += b->az*t;
    b->px += b->vx*t;
    b->py += b->vy*t;
    b->pz += b->vz*t;
    b->ax = 0.f;
    b->ay = 0.f;
    b->az = 0.f;

    if (b->px > LEN_MAX) {
        b->px = LEN_MAX;
        if (b->vx > 0.f) b->vx*=-1.f;
    }
    else if (b->px < LEN_MIN) {
        b->px = LEN_MIN;
        if (b->vx < 0.f) b->vx*=-1.f;
    }
    if (b->py > LEN_MAX) {
        b->py = LEN_MAX;
        if (b->vy > 0.f) b->vy*=-1.f;
    }
    else if (b->py < LEN_MIN) {
        b->py = LEN_MIN;
        if (b->vy < 0.f) b->vy*=-1.f;
    }
    if (b->pz > LEN_MAX) {
        b->pz = LEN_MAX;
        if (b->vz > 0.f) b->vz*=-1.f;
    }
    else if (b->pz < LEN_MIN){
        b->pz = LEN_MIN;
        if (b->vz < 0.f) b->vz*=-1.f;
    }
}

// function that initializes a node and returns a pointer to it. this function
// is kept separate from the recursive call to ensure the pointer can be freed
// before the call takes place.
// recursively divide bodies uint32_to quadrants and create new child nodes from
// each of these quadrants. stops when maximum node depth is reached or when
// exactly one body is found in quadrant. quadrants with no bodies are not
// initialized as nodes.
void init_node(node_t *node, body_t **p_bodies, uint32_t p_nbodies, uint8_t depth, uint8_t max_depth, float pos[], ttable_t magic, float *length, float g) {
    uint32_t i;             // bodies iterator
    uint8_t q;              // quadrants iterator
    uint32_t n = p_nbodies; // number of bodies in node

    //node_id++;
    node->depth = depth;    // depth determines size and position of node
    node->child = 0;        // pointer to child nodes initialized to null
    node->num_child = 0;    // number of child nodes

    if (n == 1) {           // only 1 body in node; recursion stop condition
        body_t *b = p_bodies[0];
        node->cx = b->px;
        node->cy = b->py;
        node->cz = b->pz;
        node->gm = g*b->m;
        free(p_bodies);
        return;
    }
    else if (n > 1) {
        // calculate center of mass for each axis and total mass
        float t_tm = 0.f;
        float t_cx = 0.f;
        float t_cy = 0.f;
        float t_cz = 0.f;

        for (i = 0; i < n; i++) {
            float t_m = p_bodies[i]->m;
            t_cx += p_bodies[i]->px*t_m;
            t_cy += p_bodies[i]->py*t_m;
            t_cz += p_bodies[i]->pz*t_m;
            t_tm += t_m;
        }
        float invtm = 1.f/t_tm;
        node->cx = t_cx*invtm;
        node->cy = t_cy*invtm;
        node->cz = t_cz*invtm;
        node->gm = g*t_tm;

        if (depth >= max_depth) {
            free(p_bodies);
            return;
        }

        // quad keeps tracks of which bodies are in which quadrants by
        // index. num_quad keeps track of the number of bodies in each
        // of quad's quadrants.
        body_t **quad[8];
        uint32_t num_quad[8] = { 0 };
        for (q = 0; q < 8; q++) {
            quad[q] = malloc(n*sizeof(body_t *));
        }
        for (i = 0; i < n; i++) {
            body_t *b = p_bodies[i];
            // derive quadrant q=0..7 from relative position on each axis
            // of body to center of current node
            bool b_x = b->px > pos[0];
            bool b_y = b->py > pos[1];
            bool b_z = b->pz > pos[2];
            // each quadrant is mapped to a number between 0..7
            uint8_t q = b_x*4 + b_y*2 + b_z;
            // keep track of all bodies in each quadrant q
            quad[q][num_quad[q]] = b;
            // keep track of number of bodies in each quadrant q
            num_quad[q]++;
        }
        free(p_bodies);

        uint8_t c = 0, nq[8];
        float n_pos[8][3];
        float len = 0.25f*length[depth];
        for (q = 0; q < 8; q++) {
            if (num_quad[q] > 0) {
                // magic functions to determine new position based on
                // quadrant and length/size of node
                n_pos[c][0] = pos[0] + len*magic.x[q];
                n_pos[c][1] = pos[1] + len*magic.y[q];
                n_pos[c][2] = pos[2] + len*magic.z[q];
                nq[c] = q;
                c++;
            }
            else free(quad[q]);
        }
        // allocate new nodes
        node->child = malloc(c*sizeof(node_t));
        node->num_child = c;
#pragma omp parallel for
        for (q = 0; q < c; q++) {
            // initialize new nodes the recursive call *MUST BE PLACED LAST* to
            // ensure the tail call optimization
            init_node(&node->child[q], quad[nq[q]], num_quad[nq[q]], depth+1, max_depth, n_pos[q], magic, length, g);
        }
    }
}

// recursively free nodes
void free_node(node_t *node) {
    uint8_t i;
    if (node->num_child > 0) {
#pragma omp parallel for
        for (i = 0; i < node->num_child; i++) {
            free_node(&node->child[i]);
        }
        free(node->child);
    }
}

/* Recursively checks nodes and decides whether to calculate acceleration based
 * on given node or to check its children instead. Tail call recursion is
 * guaranteed with -O2 as long as check_node is called at the end and no lines
 * follow it. >90% of the time is spent in this subroutine.
 */
void check_node(node_t *node, body_t *body, float *length) {
    node_t temp = *node;
    float rx = temp.cx - body->px;
    float ry = temp.cy - body->py;
    float rz = temp.cz - body->pz;
    float r = rx*rx + ry*ry + rz*rz;
    if (r > 0.f) {
        r = 1.f/sqrtf(r + E_SQR);
        if (temp.num_child != 0 && length[temp.depth]*r > DIST_THRES) {
            for (uint8_t i = 0; i < temp.num_child; i++) {
                check_node(&temp.child[i], body, length);
            }
        }
        else {// r != E | (ratio > DIST_THRES & temp->num_child > 0)
            float gmr = temp.gm*(r*r*r);
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
        printf("%.3f %.3f %.3f\n", node->cx, node->cy, node->cz);
    } else {
        printf("%.3f %.3f %.3f\n", node->cx, node->cy, node->cz);
        for (uint32_t j = 0; j < node->depth; j++) {
            printf("\t");
        }
        printf("%d %.3f %.3f %.3f %.3f\n", node->depth, node->cx, node->cy, node->cz, node->gm);
        for (uint32_t i = 0; i < node->num_child; i++) {
            print_node(&node->child[i], np_bodies);
        }
    }
}

nbodysys_t *init_rand_nbodysys(uint32_t n, float max_p, float max_v, float max_m) {
    nbodysys_t *temp = malloc(sizeof(nbodysys_t));
    temp->bodies = malloc(n*sizeof(body_t));
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
    nbodysys_t *dest = malloc(sizeof(nbodysys_t));
    dest->num_bodies = n;
    dest->maxp = src->maxp;
    dest->maxv = src->maxv;
    dest->maxm = src->maxm;
    dest->bodies = malloc(n*sizeof(body_t));
    memcpy(dest->bodies, src->bodies, n*sizeof(body_t));

    return dest;
}

void free_nbodysys(nbodysys_t *nb) {
    free(nb->bodies);
    free(nb);
}

void print_nbodysys(nbodysys_t *nb) {
    printf("\n");
    uint32_t a = nb->num_bodies > DBG_DISPLAY ? nb->num_bodies/DBG_DISPLAY : 1;
    for (uint32_t i = 0; i < nb->num_bodies; i+=a) {
        printf("%.9f\t%.9f\t%.9f\n", nb->bodies[i].px, nb->bodies[i].py, nb->bodies[i].pz);
    }
}

void brute(nbodysys_t *nb, uint32_t iters, del_t time) {
    uint32_t iter=0, i, j, k, n;
    n = nb->num_bodies;
    uint32_t n1 = n-1;
    uint32_t n_iter = n*n1/2;
    float *gm = malloc(n*sizeof(float));
    for (i = 0; i < n; i++) {
        gm[i] = G*nb->bodies[i].m;
    }
    while (iter++ < iters) {
#pragma omp parallel for private(i,j)
        for (k = 0; k < n_iter; k++) {
            i = k/n1; j = k%n1;
            if (j >= i) {
                i = n1 - i;
                j = n1 - j - 1;
            }
            float rx = nb->bodies[i].px - nb->bodies[j].px;
            float ry = nb->bodies[i].py - nb->bodies[j].py;
            float rz = nb->bodies[i].pz - nb->bodies[j].pz;
            float r = 1.f/sqrtf(rx*rx + ry*ry + rz*rz + E_SQR);
            r *= r*r;
            float gmr = nb->bodies[j].m*G*r;
            nb->bodies[i].ax -= gmr*rx;
            nb->bodies[i].ay -= gmr*ry;
            nb->bodies[i].az -= gmr*rz;
            gmr = nb->bodies[i].m*G*r;
            nb->bodies[j].ax += gmr*rx;
            nb->bodies[j].ay += gmr*ry;
            nb->bodies[j].az += gmr*rz;
        }
#pragma omp parallel for
        for (i = 0; i < n; i++)
            update_body(&nb->bodies[i], time);
    }
    free(gm);
}

void barnes(nbodysys_t *nb, uint32_t iters, del_t time) {
    uint32_t iter=0, i;
    uint32_t n = nb->num_bodies;
    node_t *root_node = malloc(sizeof(node_t));
    float g = G;
    //omp_set_nested(1);

    // maximum tree depth should be proportional to the number of nodes to the
    // inverse power of eight. the size of a node is proportional to its depth
    uint8_t max_dep = (uint8_t)ceil(DEPTH*pow(nb->num_bodies, 0.125));
    float *length = malloc(max_dep*sizeof(float));
    for (uint8_t d = 0; d < max_dep; d++) {
        length[d] = pow(0.5f, d)*MAX_P;
    }
    if (debug == 2) printf("max depth: %d\n", max_dep);

    // acceleration is calculated each iteration and its value is assumed to be
    // independent of its previous value (starts from 0).
    float pos[3] = { 0.f };

    // boolean truth table to determine a quadrant's position relative to the
    // center of the node, with a unique 3-bit value (in 3 dimensions) for each
    // quadrant. for instance, quadrant 0 is x negative, y negative, and z
    // negative relative to the center of the node.
    ttable_t magic;
    for (uint8_t q = 0; q < 8; q++) {
        magic.x[q] = q/4 == 0 ? -1 : 1;
        magic.y[q] = (q/2)%2 == 0 ? -1 : 1;
        magic.z[q] = q%2 == 0 ? -1 : 1;
    }

    printf("pow!\n");
    while (iter++ < iters) {
        node_id = 0;
        body_t **root_bodies = malloc(n*sizeof(body_t *));
        for (i = 0; i < n; i++){
            root_bodies[i] = &nb->bodies[i];
        }
        init_node(root_node, root_bodies, n, 0, max_dep, pos, magic, length, g);
#pragma omp parallel for
        for (i = 0; i < n; i++)
            check_node(root_node, &nb->bodies[i], length);
#pragma omp parallel for
        for (i = 0; i < n; i++)
            update_body(&nb->bodies[i], time);
        free_node(root_node);
    }
    free(length);
    free(root_node);
}

