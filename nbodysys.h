#ifndef __NBODYSYS_H__
#define __NBODYSYS_H__

#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "vec_ops.h"

class NBodySys {
    public:
        NBodySys(int num);
        NBodySys(const NBodySys &obj);
        ~NBodySys();

        void update_full_par(float t);
        void update_full_seq(float t);
        void print();

    protected:
        float3 *pos;
        float *mass;
        float3 *vel;
        int num;
        float g;
};

#endif
