#include "nbodysys.c"

#define NUM_BODIES 8192
#define NUM_ITERS 1
#define MAX_P 100
#define MAX_V 10
#define MAX_M 100
#define DEBUG 0

int main(int argc, char **argv) {
    int num_bodies, num_iters;
    double dstart, dstop;
    update f = all_seq;

    if (argc > 1) {
        num_bodies = strtol(argv[1], NULL, 10);
        num_iters = strtol(argv[2], NULL, 10);
    } else {
        num_bodies = NUM_BODIES;
        num_iters = NUM_ITERS;
    }

    nbodysys_t *nbody = init_rand_nbodysys(num_bodies, MAX_P, MAX_V, MAX_M);
    
    if (DEBUG == 1) print_nbodysys(nbody);

    dstart = omp_get_wtime();
    for (int i = 0; i < num_iters; i++) {
        f(nbody, 1.f);
    }
    dstop = omp_get_wtime();
    if (DEBUG == 1) print_nbodysys(nbody);
    printf("\n%.5f\n", dstop-dstart);

    fin_nbodysys(nbody);

    return 0;
}
