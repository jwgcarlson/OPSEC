#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <vector>
using std::vector;

#include "Model.h"
#include "cfg.h"

void print_usage(FILE* f, char* progname) {
    int myid;
    fprintf(f, "Usage: %s [SWITCHES] [name=value ...]\n", progname);
    fprintf(f,
"Switches:\n"
"  -h            Display this help\n"
"  -c FILE       Read additional configuration options from FILE\n"
"Configuration options:\n"
"  model=NAME    Read list of cells from FILE (required)\n"
    );
}

int main(int argc, char* argv[]) {
    Config cfg = cfg_new();

    /* Parse command line switches */
    int opt;
    const char* optstring = "hc:";
    while((opt = getopt(argc, argv, optstring)) != -1) {
        switch(opt) {
        case 'h':
            print_usage(stdout, argv[0]);
            return 0;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            print_usage(stderr, argv[0]);
            return 1;
        }
    }

    /* Parse additional command line options */
    for(int i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);

    Model* model = InitializeModel(cfg);
    if(!model) {
        fprintf(stderr, "Error loading model. Aborting.\n");
        return 1;
    }
    XiFunc xi = model->GetXi();
    int Nparams = model->NumParams();
    vector<XiFunc> xin(Nparams);
    for(int n = 0; n < Nparams; n++)
        xin[n] = model->GetXiDeriv(n);

    printf("# r  xi(r)  xi,1(r)  xi,2(r) ... xi,N(r)\n");
    double rmin = 1e-1;
    double rmax = 1e3;
    double Nr = 16384;
    for(int i = 0; i < Nr; i++) {
        double r = rmin * exp(i*log(rmax/rmin)/(Nr-1));
        double x1[3] = { 100, 0, 0 };
        double x2[3] = { 100+r, 0, 0 };
        printf("%e %e", r, xi(x1, x2));
        for(int n = 0; n < Nparams; n++)
            printf(" %e", xin[n](x1, x2));
        printf("\n");
    }

    return 0;
}
