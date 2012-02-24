#ifndef MODEL_H
#define MODEL_H


/* A "model" in OPSEC is a description of how the 2-point clustering of tracers
 * depends on a set of parameters { p_1, p_2, ..., p_N }.  For our purposes, we
 * need to be able to evaluate both the full 2-point correlation function
 *   \xi(\vec{r}_1, \vec{r}_2) ,
 * as well as its partial derivatives with respect to the parameters,
 *   \xi_{,n}(\vec{r}_1, \vec{r}_2) = \partial \xi(\vec{r}_1, \vec{r}_2) / \partial p_n .
 * Both these quantities will, in general, depend on the very parameter values
 * $p_n$ that we are attempting to estimate, so in order to make progress we
 * also need to have prior estimates $p_n^{(0)}$ for these parameters.
 * Ultimately we will compute estimators $\hat{p}_n$ that are unbiased (to
 * second order in $p - p^{(0)}$), so the choice of priors mostly affects the
 * variance of these estimates; the variance is minimized when the prior
 * guesses correspond to truth.
 *
 * The 2-point correlation function $\xi$ (and its derivatives $\xi_{,n}$) in
 * general depends on the locations of the two points $\vec{r}_1$ and
 * $\vec{r}_2$.  For a statistically isotropic field, however, $\xi$ depends
 * only on the geometry of the triangle formed by the two points and the
 * observer.  We parameterize this geometry by the three quantities $r$, $a$,
 * and $b$, illustrated below:
 *                           r                                                    *
 *                   1  ___________  2                                            *
 *                     /_/    A \,'                                               *
 *                    /  B     ,'                                                 *
 *                   /       ,'              T = \theta                           *
 *                  /      ,'                A = \alpha                           *
 *               a /     ,'                  B = \beta                            *
 *                /    ,' b                                                       *
 *               / 2T,'                                                           *
 *              /-.,'                                                             *
 *             / ,'                                                               *
 *            /,'                                                                 *
 *           /'                                                                   *
 *         O                                                                      *
 *
 * The Model class below defines the common interface that all OPSEC models
 * must implement.  Correlation functions are abstracted by XiFunc objects, for
 * which each model should provide its own XiFuncImpl implementation (see
 * RealModel.cpp for an example).
 *
 * Note that if you add your own model subclass, you should update the
 * InitializeModel() method in Model.cpp so that it can be loaded automatically
 * via your config file! */


class Spline;
class XiFunc;

/*****************************************************************************
 * Model
 *
 * Conceptually, a model in OPSEC describes how the 2-point correlation
 * function $\xi$ depends on a set of parameters $p_n$.  A concrete
 * implementation of the Model class must define the following methods:
 *  int NumParams()
 *    - returns the number of parameters in the model
 *  double GetParam(int n);
 *    - returns a prior estimate of parameter $p_n$
 *  XiFunc GetXi();
 *    - returns a XiFunc object that represents a prior estimate of 2-point
 *      correlation function $\xi$ for the model
 *  XiFunc GetXiDeriv(int n);
 *    - returns a XiFunc object that represents a prior estimate for the
 *      derivative $\xi_{,n} = \partial \xi/\partial p_n$ of the 2-point
 *      correlation function with respect to parameter $p_n$
 ******************************************************************************/
class Model {
public:
    Model();
    virtual ~Model();

    /* All Model subclasses must implement the following methods */
    virtual int NumParams() = 0;                // return the number of parameters
    virtual double GetParam(int n) = 0;         // get the prior value for p_n
    virtual XiFunc GetXi() = 0;                 // get the model 2-point function \xi
    virtual XiFunc GetXiDeriv(int n) = 0;       // get the partial derivative d\xi/dp_n
};


#include "cfg.h"

/* Global function to construct a concrete Model instance based on the 'model'
 * option in a config file. */
Model* InitializeModel(Config cfg);


/* TODO: move this somewhere more appropriate */
/* Convenience function for computing $\xi_l^m(r)$ from $P(k)$, where
 *   \xi_l^m(r) = \frac{1}{2\pi^2} \int_0^\infty dk k^m P(k) j_l(kr)
 * In practice the integral is performed over [kmin,kmax] using Simpson's rule,
 * with a damping factor $e^{-(k/kdamp)^2}$ applied to reduce ringing. */
void ComputeXiLM(int l, int m, Spline P, int Nr, const double r[], double xi[],
                 int Nk = 32768, double kmin = 0., double kmax = 10., double kdamp = 1.);

#endif // MODEL_H
