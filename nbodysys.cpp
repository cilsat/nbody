#include "nbodysys.h"
#include <string.h>

#define MAX_VAL_X 100   // maximum positional value
#define MAX_VAL_V 10    // maximum (initial) velocity
#define MAX_VAL_M 100   // maximum mass of body
#define G 6.67408e-11   // universal gravitational constant
#define MAX_DEPTH 64    // maximum depth barnes hut oct tree
#define MGC 0.1767766952966369
#define DEBUG 2

inline void min(float &a, const float &b) {
    a = a < b ? a : b;
}

inline void max(float &a, const float &b) {
    a = a > b ? a : b;
}

Body::Body() {
}

Body::Body(float px, float py, float pz, float vx, float vy, float vz, float m) :
    posX(px), posY(py), posZ(pz),
    velX(vx), velY(vy), velZ(vz),
    mass(m)
{}

Body::~Body() {
}

void Body::updatePos() {
}

void Body::print_params() {
    printf("%.3f %.3f %.3f", posX, posY, posZ);
}

Node::Node() {}

Node::Node(unsigned int nDepth) {
    depth = nDepth;
}

Node::~Node() {
}

void Node::reset() {
    bodies.clear();

    for (unsigned int i = 0; i < children.size(); i++) {
        children[i]->reset();
        delete children[i];
    }
    children.clear();
    hasChildren = false;
}

void Node::print_params() {
    float len = 0.5*size/sqrt(3);
    float x0 = posX - len;
    float x1 = posX + len;
    float y0 = posY - len;
    float y1 = posY + len;
    float z0 = posZ - len;
    float z1 = posZ + len;
    printf("d=%d n=%ld %.3f,%.3f %.3f,%.3f %.3f,%.3f %.3f\n", depth, bodies.size(), x0, x1, y0, y1, z0, z1, size);
}

void Node::set_params(std::vector<Body*> nbodies, float px, float py, float pz, float sz) {
    unsigned int i;
    bodies = nbodies;
    hasChildren = false;
    posX = px; posY = py; posZ = pz; size = sz;

    float totX = 0, totY = 0, totZ = 0;

    if (bodies.size() > 0) {
        tm = 0;
        for (i = 0; i < bodies.size(); i++) {
            bodies[i] = bodies[i];
            float m = bodies[i]->mass;
            totX += m*bodies[i]->posX;
            totY += m*bodies[i]->posY;
            totZ += m*bodies[i]->posZ;
            tm += m;
        }
        float invTm = 1.f/tm;
        cmX = totX*invTm;
        cmY = totY*invTm;
        cmZ = totZ*invTm;
    } else {
        tm = 0;
        cmX = 0;
        cmY = 0;
        cmZ = 0;
    }

    if ((DEBUG == 2) && (depth == MAX_DEPTH - 1)) {
        for (i = 0; i < bodies.size(); i++) 
            printf("not found: %d %.3f %.3f %.3f\t%f %f %f\n", depth, bodies[i]->posX, bodies[i]->posY, bodies[i]->posZ, posX, posY, posZ);
    }
    if (bodies.size() > 1 && depth < MAX_DEPTH) {
        set_children();
    }
}

void Node::set_children() {
    unsigned int i;
    std::vector<Body*> q0, q1, q2, q3, q4, q5, q6, q7;
    std::vector<Body*> quad[8];

    for (i = 0; i < bodies.size(); i++) {
        float x = bodies[i]->posX;
        float y = bodies[i]->posY;
        float z = bodies[i]->posZ;

        if (x < posX) {
            if (y < posY) {
                if (z < posZ) {
                    q0.push_back(bodies[i]);
                    quad[0].push_back(bodies[i]);
                }
                else {
                    q1.push_back(bodies[i]);
                    quad[1].push_back(bodies[i]);
                }
            }
            else {
                if (z < posZ) {
                    q2.push_back(bodies[i]);
                    quad[2].push_back(bodies[i]);
                }
                else {
                    q3.push_back(bodies[i]);
                    quad[3].push_back(bodies[i]);
                }
            }
        }
        else {
            if (y < posY) {
                if (z < posZ) {
                    q4.push_back(bodies[i]);
                    quad[4].push_back(bodies[i]);
                }
                else {
                    q5.push_back(bodies[i]);
                    quad[5].push_back(bodies[i]);
                }
            }
            else {
                if (z < posZ) {
                    q6.push_back(bodies[i]);
                    quad[6].push_back(bodies[i]);
                }
                else {
                    q7.push_back(bodies[i]);
                    quad[7].push_back(bodies[i]);
                }
            }
        }
    }

    const double mgc = MGC*size;

    for (i = 0; i < 8; i++) {
        if (quad[i].size() > 0) {
            double mgcX = i < 4 ? -mgc : mgc;
            double mgcY = (i/2)%2 == 0 ? -mgc : mgc;
            double mgcZ = i%2 == 0 ? -mgc : mgc;

            Node *nq = new Node(depth+1);
            if ((DEBUG == 2) && (quad[i].size() == 1)) {
                Body *temp = quad[i][0];
                printf("found: %d %.3f %.3f %.3f\n", depth, temp->posX, temp->posY, temp->posZ);
            }
            nq->set_params(quad[i], posX+mgcX, posY+mgcY, posZ+mgcZ, 0.5*size);

            children.push_back(nq);
            hasChildren = true;
        }
    }

    /*
    Node *nq0 = new Node(depth + 1);
    Node *nq1 = new Node(depth + 1);
    Node *nq2 = new Node(depth + 1);
    Node *nq3 = new Node(depth + 1);
    Node *nq4 = new Node(depth + 1);
    Node *nq5 = new Node(depth + 1);
    Node *nq6 = new Node(depth + 1);
    Node *nq7 = new Node(depth + 1);

    nq0->set_params(q0, posX-mgc, posY-mgc, posZ-mgc, 0.5*size);
    nq1->set_params(q1, posX-mgc, posY-mgc, posZ+mgc, 0.5*size);
    nq2->set_params(q2, posX-mgc, posY+mgc, posZ-mgc, 0.5*size);
    nq3->set_params(q3, posX-mgc, posY+mgc, posZ+mgc, 0.5*size);
    nq4->set_params(q4, posX+mgc, posY-mgc, posZ-mgc, 0.5*size);
    nq5->set_params(q5, posX+mgc, posY-mgc, posZ+mgc, 0.5*size);
    nq6->set_params(q6, posX+mgc, posY+mgc, posZ-mgc, 0.5*size);
    nq7->set_params(q7, posX+mgc, posY+mgc, posZ+mgc, 0.5*size);

    children.push_back(nq0);
    children.push_back(nq1);
    children.push_back(nq2);
    children.push_back(nq3);
    children.push_back(nq4);
    children.push_back(nq5);
    children.push_back(nq6);
    children.push_back(nq7);

    hasChildren = true;
    */
}

const std::vector<Node*>& Node::get_children() {
    return children;
}

Node *Node::get_children(int n) {
    return children[n];
}

const std::vector<Body*>& Node::get_bodies() {
    return bodies;
}

Body *Node::get_bodies(int n) {
    return bodies[n];
}

NBodySys::NBodySys(int n) {
    int i;
    pos = new float3[n];
    mass = new float[n];
    vel = new float3[n];
    num = n;
    g = 1.f;

    srand(time(NULL));
    for (i = 0; i < n; i++) {
        pos[i] = rand_float3(MAX_VAL_X) - (float)0.5*MAX_VAL_X;
        mass[i] = MAX_VAL_M*((float) rand() / (float) RAND_MAX);
        vel[i] = {};
    }
}

void NBodySys::copy(NBodySys *obj) {
    int n = obj->num;

    memcpy(pos, obj->pos, n*sizeof(float3));
    memcpy(mass, obj->mass, n*sizeof(float));
    memcpy(vel, obj->vel, n*sizeof(float3));
    num = n;
    g = obj->g;
}

NBodySys::~NBodySys() {
    delete[] pos;
    delete[] mass;
    delete[] vel;
}

void NBodySys::print() {
    int i;
    int n = num;

    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%.3f %.3f %.3f\n", pos[i].x, pos[i].y, pos[i].z);
    }
    /*
    printf("\n");
    for (i = 0; i < n; i++) {
        printf("%.3f %.3f %.3f\n", vel[i].x, vel[i].y, vel[i].z);
    }
       printf("\n");
       for (i = 0; i < n; i++) {
       printf("%.3f ", mass[i]);
       }*/
    printf("\n");
}

void NBodySys::allpairs_seq(float t) {
    int n = num;
    int i, j;
    float3 *force = new float3[n];

    for (i = 0; i < n; i++) {
        force[i] = {};
        for (j = 0; j < n; j++) {
            if (i != j) {
                float3 xr = pos[j] - pos[i];
                float r = invsqrt_float3(xr);
                force[i] += xr*(mass[j]*(r*r*r));
            }
        }
        float3 vu = vel[i];
        vel[i] += force[i]*(g*t);
        pos[i] += (vel[i] - vu)*(0.5*t);
    }

    delete[] force;
}

void NBodySys::allpairs_par(float t) {
    int n = num;
    int i, j;
    float3 *force = new float3[n];
    float gt = g*t;
    float ht = 0.5f*t;

#pragma omp parallel for
    for (i = 0; i < n; i++) {
        force[i] = {};
#pragma omp parallel for
        for (j = 0; j < n; j++) {
            if (i != j) {
                float3 xr = {pos[j].x-pos[i].x, pos[j].y-pos[i].y, pos[j].z-pos[i].z};
                float r = invsqrt_float3(xr);
                float gr = mass[j]*(r*r*r);
                float3 f = {xr.x*gr, xr.y*gr, xr.z*gr};
                force[i] += f;
            }
        }
        float3 u = vel[i];
        float3 v = {force[i].x*gt, force[i].y*gt, force[i].z*gt};
        float3 s = {(v.x-u.x)*ht, (v.y-u.y)*ht, (v.z-u.z)*ht};
        vel[i] += v;
        pos[i] += s;
    }

    delete[] force;
}

void NBodySys::print_tree(Node *tree) {
    std::vector<Node*> temp = tree->get_children();
    for (unsigned int i = 0; i < temp.size(); i++) {
        temp[i]->print_params();
        printf("%.3f\n", temp[i]->size);
        print_tree(temp[i]);
    }
    temp.clear();
}

void NBodySys::barneshut_seq(float t) {
    int i;

    // BUILD TREE
    std::vector<Body*> nbodies;

    // get bounding box of simulation
    float maxpos = 0;

    for (i = 0; i < num; i++) {
        max(maxpos, fabs(pos[i].x));
        max(maxpos, fabs(pos[i].y));
        max(maxpos, fabs(pos[i].z));
        Body *temp = new Body(pos[i].x, pos[i].y, pos[i].z, vel[i].x, vel[i].y, vel[i].z, mass[i]);
        nbodies.push_back(temp);
    }
    float size = sqrtf(3*maxpos*maxpos);
    printf("%.3f %.3f\n", maxpos, size);

    Node *root = new Node(0);
    root->set_params(nbodies, 0, 0, 0, size);

    // UPDATE
    print_tree(root);

    root->reset();
}
