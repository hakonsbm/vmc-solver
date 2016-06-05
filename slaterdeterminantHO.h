#ifndef SLATERDETERMINANTHO_H
#define SLATERDETERMINANTHO_H
#include <armadillo>
#include "vmcsolver.h"
#include "derivatives.h"


using namespace arma;

class VMCSolver;

class SlaterDeterminantHO
{

public:
    SlaterDeterminantHO();
    ~SlaterDeterminantHO();



    void updateSlaterMatrices(const mat &r, VMCSolver *solver);
    double calculateDeterminant(VMCSolver *solver);
    //double phi(const mat &r, double alpha, int i, int j, VMCSolver *solver);
//    double SlaterDeterminant::determinantRatioUp(const mat &r, VMCSolver *solver, Derivatives *der);
//    double SlaterDeterminant::determinantRatioDown(const mat &r, VMCSolver *solver, Derivatives *der);
    vec gradientPhi(const mat &r, int i, int j, VMCSolver *solver);
    double laplacianPhi(const mat &r, int i, int j, VMCSolver *solver);

    mat gradientSlaterDeterminant(const mat &r , VMCSolver *solver);
    double laplacianSlaterDeterminant(const mat &r , VMCSolver *solver);

    //For the molecule we need a different wavefunctions
    double phiMolecule(const mat &r, double alpha, int i, int j, VMCSolver *solver);

    mat detUpOld;
    mat detDownOld;
    mat detUpInverseOld;
    mat detDownInverseOld;

    int updateSlaterCounter = 0;
    int updateDeterminantCounter = 0;

private:
    bool useGTO;



};

#endif // SLATERDETERMINANTHO_H
