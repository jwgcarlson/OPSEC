#ifndef REALMODEL_H
#define REALMODEL_H

#include <vector>

#include "Model.h"
#include "Spline.h"

/* Required configuration options:
 *   pkfile:
 *   bands:
 *
 * Parameters:
 *  n = 0 through Nbands-1: band power p_n
 */

class RealModel : public Model {
public:
    RealModel(Config cfg);
    ~RealModel();

    int NumParams();
    double GetParam(int n);
    XiFunc GetXi();
    XiFunc GetXiDeriv(int n);

protected:
    struct Band { double min, max; };
    int Nbands;
    std::vector<Band> bands;

    Spline pk;          // splined power spectrum $P(k)$
};

#endif // REALMODEL_H
