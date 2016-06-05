#ifndef HYDROGENTWO_H
#define HYDROGENTWO_H

#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class HydrogenTwo : public TrialFunction
{
public:

    ~HydrogenTwo();

    HydrogenTwo(VMCSolver* solver);
    virtual void setSpin(VMCSolver *solver);
    virtual void updateSlaterDeterminant(VMCSolver *solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver );
    virtual double lnDerivativeWaveFunction(const mat &r, VMCSolver *solver);
    virtual double lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver);
    virtual void calculateAlpha(VMCSolver *solver);


private:
    double SD;


};

#endif // HYDROGENTWO_H
