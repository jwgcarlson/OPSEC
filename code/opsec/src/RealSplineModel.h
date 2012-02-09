#ifndef REALSPLINEMODEL_H
#define REALSPLINEMODEL_H

#include <vector>

#include "Model.h"
#include "Spline.h"

class RealSplineModel : public Model {
public:
    RealSplineModel(Config cfg);
    ~RealSplineModel();

    int NumParams();
    XiFunc GetXi();
    XiFunc GetXiDeriv(int n);

protected:
    int Nknots;                 // number of (k,p) points in P(k) template, not including (0,0)
    std::vector<double> knots;  // k values of knots, not including k=0

    Spline pk;
};

#endif // REALSPLINEMODEL_H
