#include "nbodysys.cpp"

#define NB 1000  // number of bodies
#define ITERS 100   // number of iterations
#define DEL 1   // (integer) length of time steps in seconds

int main(int argc, char **argv) {
    double dstart, dstop;
    int nb, iters;

    if (argc > 1) {
        nb = strtol(argv[1], NULL, 10);
        iters = strtol(argv[2], NULL, 10);
    } else {
        nb = NB; iters = ITERS;
    }

    NBodySys *sys = new NBodySys(nb);
    NBodySys *temp = new NBodySys(nb);
    temp->copy(sys);

    if (DEBUG == 1) {
        //sys->print();
        temp->print();
    }

    /*
    dstart = omp_get_wtime();
    for (int i = 0; i < iters; i++) {
        sys->allpairs_seq((float)DEL);
    }
    dstop = omp_get_wtime();
    printf("%.3f\n", dstop-dstart);

    dstart = omp_get_wtime();
    for (int i = 0; i < iters; i++) {
        temp->allpairs_par((float)DEL);
    }
    dstop = omp_get_wtime();
    printf("%.3f\n", dstop-dstart);
    */

    dstart = omp_get_wtime();
    for (int i = 0; i < iters; i++) {
        temp->barneshut_seq((float)DEL);
    }
    dstop = omp_get_wtime();
    printf("%.3f\n", dstop-dstart);

    if (DEBUG == 1) {
        //sys->print();
        temp->print();
    }

    return 0;
}
