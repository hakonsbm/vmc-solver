#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

#include <armadillo>

using namespace arma;

class VMCsolver;

class BasisFunctions{
public:
    BasisFunctions();
    ~BasisFunctions(){}

    virtual double eval(double k, double k2, double *exp_factor, vec pos) = 0;
};

#endif // BASISFUNCTIONS_H
