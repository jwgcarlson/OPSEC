/* mkmatrix
 *
 * Utility to generate a matrix stored in a .abn file, for testing. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "abn.h"
#include "cfg.h"

int main(int argc, char* argv[]) {
    int n, a, b;
    double j, m;
    double *Jx;
    char* filename;
    FILE* fp;
    Config opts = cfg_new();

    if(argc != 3) {
        fprintf(stderr, "Usage: %s <n> <filename>\n", argv[0]);
        return 1;
    }

    n = atoi(argv[1]);
    filename = argv[2];

    /* Construct J_x angular momentum matrix for j = (n-1)/2. */
    j = (n - 1)/2.;
    Jx = (double*) calloc(n*n, sizeof(double));
    for(a = 0; a < n-1; a++) {
        m = -j + a;
        Jx[(a+1)*n + a] = Jx[a*n + a+1] = 0.5 * sqrt(j*(j+1) - m*(m+1));
    }
//    for(a = 0; a < n; a++) {
//        Jx[a*n + a] = 1;
//    }

    /* Print matrix to stdout, just to check */
    for(a = 0; a < n; a++) {
        for(b = 0; b < n-1; b++)
            printf("%5.3f ", Jx[a*n + b]);
        printf("%5.3f\n", Jx[a*n + n-1]);
    }

    cfg_set_int(opts, "M", n);
    cfg_set_int(opts, "N", n);
    fp = fopen(filename, "w");
    fprintf(fp, "# J_x angular momentum matrix for j = %g\n", j);
    abn_write(fp, Jx, n*n, "d", opts);
    fclose(fp);

    return 0;
}
