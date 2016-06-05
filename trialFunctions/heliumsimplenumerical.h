#ifndef HELIUMSIMPLENUMERICALLY_H
#define HELIUMSIMPLENUMERICALLY_H

#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class HeliumSimpleNumerical : public TrialFunction
{

public:
    HeliumSimpleNumerical(VMCSolver *solver);
    virtual void setSpin(VMCSolver *solver);
    virtual void updateSlaterDeterminant(VMCSolver *solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );
    virtual double lnDerivativeWaveFunction(const mat &r, VMCSolver *solver);
    virtual double lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver);

private:
    double SD;
};

#endif // HELIUMSIMPLENUMERICALLY_H
