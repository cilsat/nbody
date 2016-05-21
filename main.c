#include "nbsys.c"
#include <getopt.h>

#define NUM_BODIES 16
#define NUM_ITERS 10
#define MAX_P 100
#define MAX_V 10
#define MAX_M 100
#define UPDATE barnes

static struct option long_opt[] = {
    {"file", required_argument, NULL, 'f'},
    {"bodies", required_argument, NULL, 'n'},
    {"iterations", required_argument, NULL, 'i'},
    {"max-pos", required_argument, NULL, 'P'},
    {"max-vel", required_argument, NULL, 'V'},
    {"max-mass", required_argument, NULL, 'M'},
    {"verbose", required_argument, NULL, 'v'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}
};

void read_file(nbodysys_t* nb, char* file) {
    FILE *fstream = fopen(file, "r");
    if (fstream == NULL) {
        printf("no such file found");
        free(nb);
        exit(1);
    }
    char buffer[1024];
    char *line, *newline, *newnum;
    int i = 0;
    while ((line = fgets(buffer, sizeof(buffer), fstream)) != NULL) {
        newline = strtok(line, "\n");
        while (newline != NULL) {
            newnum = strtok(newline, ",");
            nb->bodies[i].px = strtof(newnum, NULL);
            newnum = strtok(NULL, ",");
            nb->bodies[i].py = strtof(newnum, NULL);
            newnum = strtok(NULL, ",");
            nb->bodies[i].pz = strtof(newnum, NULL);
            newnum = strtok(NULL, ",");
            nb->bodies[i].vx = strtof(newnum, NULL);
            newnum = strtok(NULL, ",");
            nb->bodies[i].vy = strtof(newnum, NULL);
            newnum = strtok(NULL, ",");
            nb->bodies[i].vz = strtof(newnum, NULL);
            newnum = strtok(NULL, ",");
            nb->bodies[i].m = strtof(newnum, NULL);
            newline = strtok(NULL, "\n");
            i++;
        }
    }
    fclose(fstream);
    free(line);
    free(newline);
}

int main(int argc, char **argv) {
    const char *short_opt = "f:n:i:P:V:M:v:h";
    int c;
    char *filename=0;
    int num_bodies=NUM_BODIES, num_iters=NUM_ITERS;
    float maxp=MAX_P, maxv=MAX_V, maxm=MAX_M;
    double dstart, dstop;
    nbodysys_t *nbody;

    while ((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1) {
        switch(c) {
            case -1:
            case 0:
                break;
            case 'f':
                filename = optarg;
                break;
            case 'n':
                num_bodies = strtol(optarg, NULL, 10);
                break;
            case 'i':
                num_iters = strtol(optarg, NULL, 10);
                break;
            case 'P':
                maxp = strtof(optarg, NULL);
                break;
            case 'V':
                maxv = strtof(optarg, NULL);
                break;
            case 'M':
                maxm = strtof(optarg, NULL);
                break;
            case 'v':
                debug = strtol(optarg, NULL, 10);
                break;
            case 'h':
                break;
        }
    }

    nbody = init_rand_nbodysys(num_bodies, maxp, maxv, maxm);
    if (filename != 0) {
        read_file(nbody, filename);
    }

    nbodysys_t *copy = copy_bodysys(nbody);
    
    if (debug == 1) {
        print_nbodysys(nbody);
    }

    dstart = omp_get_wtime();
    all_seq(nbody, num_iters, 1.f);
    dstop = omp_get_wtime();
    if (debug == 1) print_nbodysys(nbody);
    printf("\n%.5f\n", dstop-dstart);

    dstart = omp_get_wtime();
    barnes(copy, num_iters, 1.f);
    dstop = omp_get_wtime();
    if (debug == 1) print_nbodysys(copy);
    printf("\n%.5f\n", dstop-dstart);

    free_nbodysys(nbody);
    free_nbodysys(copy);

    return 0;
}
