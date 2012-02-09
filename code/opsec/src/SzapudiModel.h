#ifndef SZAPUDIMODEL_H
#define SZAPUDIMODEL_H

#include <vector>

#include "Model.h"
#include "Spline.h"

/* Parameters:
 *  n = 0 - Nbands-1: band parameters p_n
 *  n = Nbands: anisotropy parameter f
 */

class SzapudiModel : public Model {
public:
    SzapudiModel(Config cfg);
    ~SzapudiModel();

    int NumParams();
    double GetParam(int n);
    XiFunc GetXi();
    XiFunc GetXiDeriv(int n);

protected:
    struct Band { double min, max; };
    int Nbands;
    std::vector<Band> bands;

    /* Linear power spectrum */
    Spline pk;

    /* Logarithmic growth rate $f = d\ln D/d\ln a$ */
    double f;
};

#endif // SZAPUDIMODEL_H
