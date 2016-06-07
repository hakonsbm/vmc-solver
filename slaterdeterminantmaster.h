#ifndef SLATERDETERMINANTMASTER_H
#define SLATERDETERMINANTMASTER_H

#include "vmcsolver.h"
#include "derivatives.h"
//#include "slaterdeterminant.h"
//#include "slaterdeterminantHO.h"
#include <armadillo>

using namespace arma;

class VMCSolver;

class slaterDeterminantMaster
{
public:
    slaterDeterminantMaster();
    ~slaterDeterminantMaster();

    virtual void updateSlaterMatrices(const mat &r, VMCSolver *solver) = 0;
    virtual double calculateDeterminant(VMCSolver *solver) = 0;
    //virtual double phi(const mat &r, double alpha, int i, int j, VMCSolver *solver) = 0;

    virtual vec gradientPhi(const mat &r, int i, int j, VMCSolver *solver) = 0;
    virtual double laplacianPhi(const mat &r, int i, int j, VMCSolver *solver) = 0;

    virtual mat gradientSlaterDeterminant(const mat &r , VMCSolver *solver) = 0;
    virtual double laplacianSlaterDeterminant(const mat &r , VMCSolver *solver) = 0;

    //For the molecule we need a different wavefunctions
    //virtual double phiMolecule(const mat &r, double alpha, int i, int j, VMCSolver *solver) = 0;

    //void setGTO(bool usingGTO) {useGTO = usingGTO;}

    mat detUpOld;
    mat detDownOld;
    mat detUpInverseOld;
    mat detDownInverseOld;

    int updateSlaterCounter = 0;
    int updateDeterminantCounter = 0;

    //bool useGTO;
};

#endif // SLATERDETERMINANTMASTER_H
