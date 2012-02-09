#ifndef XIMODEL_H
#define XIMODEL_H

#include <vector>

#include "Model.h"
#include "Spline.h"

class XiModel : public Model {
public:
    XiModel(Config cfg);
    ~XiModel();

    int NumParams();
    XiFunc GetXi();
    XiFunc GetXiDeriv(int n);

protected:
    std::vector<double> rbands;

    Spline pk;
};

#endif // XIMODEL_H
