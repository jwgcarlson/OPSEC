#ifndef TESTMODEL_H
#define TESTMODEL_H

#include "Model.h"

class TestModel : public Model {
public:
    TestModel(Config cfg);
    ~TestModel();

    int NumParams();
    double GetParam(int n);
    XiFunc GetXi();
    XiFunc GetXiDeriv(int n);

protected:
    int form;   // 0 - exponential, 1 - gaussian
    double A;   // p_0
    double b;   // p_1
};

#endif // TESTMODEL_H
