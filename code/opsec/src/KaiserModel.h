#ifndef KAISERMODEL_H
#define KAISERMODEL_H

#include <vector>

#include "Model.h"
#include "Spline.h"

/* Required configuration options:
 *
 * Parameters:
 *  n = 0 through Nbands-1: band power p_n
 *  n = Nbands: bias b
 *  n = Nbands+1: logarithmic growth rate f
 */

class KaiserModel : public Model {
public:
    KaiserModel(Config cfg);
    ~KaiserModel();

    int NumParams();
    double GetParam(int n);
    XiFunc GetXi();
    XiFunc GetXiDeriv(int n);

protected:
    struct Band { double min, max; };
    int Nbands;
    std::vector<Band> bands;

    double f;           // logarithmic growth rate $f = d\ln D/d\ln a$
    double b;           // linear, scale-independent bias $b$
    Spline pk;          // splined prior matter power spectrum $P(k)$
};

#endif // KAISERMODEL_H
