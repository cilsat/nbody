#ifndef __NBODYSYS_H__
#define __NBODYSYS_H__

#include <stdio.h>
#include <vector>
#include <omp.h>
#include "vec_ops.h"

class Body {
    public:
        Body();
        Body(float px, float py, float pz, float vx, float vy, float vz, float m);
        ~Body();

        void updatePos();
        void print_params();

        float posX, posY, posZ;
        float velX, velY, velZ;
        float mass;
};

class Node {
    public:
        Node();
        Node(unsigned int nDepth);
        ~Node();
        void reset();

        void print_params();
        void set_params(std::vector<Body*> nBodies, float px, float py, float pz, float sz);
        void set_children();
        const std::vector<Node*>& get_children();
        Node *get_children(int n);

        const std::vector<Body*>& get_bodies();
        Body *get_bodies(int n);
        float size;

    protected:
        std::vector<Node*> children;
        std::vector<Body*> bodies;

        bool hasChildren;
        unsigned int depth;

        float posX, posY, posZ;
        float cmX, cmY, cmZ;
        float tm;
}; 

class NBodySys {
    public:
        NBodySys(int num);
        ~NBodySys();
        void copy(NBodySys *obj);

        void print();

        // all pairs
        void allpairs_seq(float t);
        void allpairs_par(float t);

        // heirarchical
        void print_tree(Node *tree);
        void barneshut_seq(float t);

    protected:
        float3 *pos;
        float *mass;
        float3 *vel;
        int num;
        float g;
};

#endif
