#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Cell.h"
#include "Model.h"
#include "TestModel.h"
#include "cfg.h"
#include "sig.h"

int main(int argc, char* argv[]) {
    /* Create test model with constant \xi = 1 */
    Config cfg = cfg_new();
    cfg_set(cfg, "model", "TestModel");
    cfg_set(cfg, "TestModel.form", "exponential");
    cfg_set_double(cfg, "TestModel.A", 1.0);
    cfg_set_double(cfg, "TestModel.b", 0.0);
    TestModel* model = new TestModel(cfg);
    XiFunc xi = model->GetXi();

    /* Create cubical cells */
    Cell c1 = { 0, 0, 0., 1., 0., 1., 0., 1., 1., 1. };
    Cell c2 = { 0, 0, 1., 1.5, 1., 1.5, 1., 1.5, 1., 0.25 };

    double S = ComputeSignalC(c1, c2, xi);
    printf("S = %g\n", S);

    return 0;
}
