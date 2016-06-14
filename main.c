#include "nbsys.c"
#include <getopt.h>

#define NUM_BODIES 16
#define NUM_ITERS 10
#define UPDATE barnes

static struct option long_opt[] = {
    {"file", required_argument, NULL, 'f'},
    {"bodies", required_argument, NULL, 'n'},
    {"iterations", required_argument, NULL, 'i'},
    {"time", required_argument, NULL, 't'},
    {"method", required_argument, NULL, 'm'},
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

static inline float min (float a, float b) {
    return a < b ? a : b;
}

static inline float err (float a, float b) {
    return fabs(a - b)/min(fabs(a), fabs(b));
}

int main(int argc, char **argv) {
    const char *short_opt = "f:n:i:P:V:M:v:h";
    int c;
    char *filename=0;
    uint32_t num_bodies=NUM_BODIES, num_iters=NUM_ITERS;
    float maxp=MAX_P, maxv=MAX_V, maxm=MAX_M, t=TIME;
    double dstart, dstop;
    nbodysys_t *nb;

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
            case 't':
                t = strtol(optarg, NULL, 10);
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

    nb = init_rand_nbodysys(num_bodies, maxp, maxv, maxm);
    if (filename != 0) {
        read_file(nb, filename);
    }

    //nbodysys_t *nb1 = copy_nbodysys(nb);
    nbodysys_t *nb2 = copy_nbodysys(nb);
    
    if (debug == 1) {
        print_nbodysys(nb);
    }
    printf("bodies: %d\n", num_bodies);
    printf("iterations: %d\n", num_iters);

    if (debug == 4) {
        uint32_t n = (uint32_t)num_bodies;
        body_t ***check_order = malloc(n*sizeof(body_t **));
        body_t **check_data = malloc(n*(n-1)*sizeof(body_t *));
        for (uint32_t i = 0; i < n; i++) {
            check_order[i] = &check_data[i*(n-1)];
        }
        dstart = omp_get_wtime();
        barnes_ordered(nb2, num_iters, t, check_order);
        dstop = omp_get_wtime();
        printf("%.5f\n", dstop-dstart);
        dstart = omp_get_wtime();
        brute_ordered(nb, num_iters, t, check_order);
        dstop = omp_get_wtime();
        printf("%.5f\n", dstop-dstart);

        free(check_order[0]);
        free(check_order);
        float dx = 0;
        float dy = 0;
        float dz = 0;
#pragma omp parallel for
        for (uint32_t i = 0; i < n; i++ ) {
            dx += err(nb->bodies[i].px, nb2->bodies[i].px);
            dy += err(nb->bodies[i].py, nb2->bodies[i].py);
            dz += err(nb->bodies[i].pz, nb2->bodies[i].pz);
        }
        printf("%.12f %.12f %.12f\n", dx/n, dy/n, dz/n);
        free_nbodysys(nb);
        free_nbodysys(nb2);

        return 0;
    }

    dstart = omp_get_wtime();
    barnes(nb2, num_iters, t);
    dstop = omp_get_wtime();
    if (debug == 1) print_nbodysys(nb2);
    printf("%.5f\n", dstop-dstart);

    /*
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n-1; j++) {
            printf("%d ", check_order[i][j]);
        }
        printf("\n");
    }*/

    if ((debug == 3) || (debug == 1)) {
        dstart = omp_get_wtime();
        brute(nb, num_iters, t);
        dstop = omp_get_wtime();
        if (debug == 1) print_nbodysys(nb);
        printf("%.5f\n", dstop-dstart);
    }

    if (debug == 3) {
        float dx = 0;
        float dy = 0;
        float dz = 0;
        int n = nb->num_bodies;
#pragma omp parallel for
        for (int i = 0; i < n; i++ ) {
            dx += err(nb->bodies[i].px, nb2->bodies[i].px);
            dy += err(nb->bodies[i].py, nb2->bodies[i].py);
            dz += err(nb->bodies[i].pz, nb2->bodies[i].pz);
        }
        printf("%.12f %.12f %.12f\n", dx/n, dy/n, dz/n);
    }

    free_nbodysys(nb);
    free_nbodysys(nb2);

    return 0;
}
