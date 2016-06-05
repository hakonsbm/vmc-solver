#ifndef QUANTUMDOTS_H
#define QUANTUMDOTS_H


#include "trialfunction.h"
#include "../vmcsolver.h"

#include <armadillo>

using namespace arma;

class VMCSolver;

class QuantumDots : public TrialFunction
{

public:
    QuantumDots(VMCSolver* solver);
    ~QuantumDots(){}
    virtual void updateSlaterDeterminant(VMCSolver *solver);
    virtual double waveFunction(const mat &r, VMCSolver*  solver);
    virtual double localEnergy(const mat &r, VMCSolver *solver);
    virtual double lnDerivativeWaveFunction(const mat &r, VMCSolver *solver);
    virtual double lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver);

private:
    double SD;





};

#endif // QUANTUMDOTS_H
