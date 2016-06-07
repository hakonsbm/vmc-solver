#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H
#include <armadillo>
#include "vmcsolver.h"
#include "derivatives.h"

#include "slaterdeterminantmaster.h"

using namespace arma;

class VMCSolver;

class SlaterDeterminant : public slaterDeterminantMaster
{
public:
    SlaterDeterminant();
    ~SlaterDeterminant();

    virtual void updateSlaterMatrices(const mat &r, VMCSolver *solver);
    virtual double calculateDeterminant(VMCSolver *solver);
    double phi(const mat &r, double alpha, int i, int j, VMCSolver *solver);
//    double SlaterDeterminant::determinantRatioUp(const mat &r, VMCSolver *solver, Derivatives *der);
//    double SlaterDeterminant::determinantRatioDown(const mat &r, VMCSolver *solver, Derivatives *der);
    virtual vec gradientPhi(const mat &r, int i, int j, VMCSolver *solver);
    virtual double laplacianPhi(const mat &r, int i, int j, VMCSolver *solver);

    virtual mat gradientSlaterDeterminant(const mat &r , VMCSolver *solver);
    virtual double laplacianSlaterDeterminant(const mat &r , VMCSolver *solver);

    //For the molecule we need a different wavefunctions
    double phiMolecule(const mat &r, double alpha, int i, int j, VMCSolver *solver);

    void setGTO(bool usingGTO) {useGTO = usingGTO;}
    bool useGTO;
    mat detUpOld;
    mat detDownOld;
    mat detUpInverseOld;
    mat detDownInverseOld;

private:



};

#endif // SLATERDETERMINANT_H
