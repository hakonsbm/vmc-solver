#ifndef ORBITALS_H
#define ORBITALS_H 

#include <armadillo>
#include "vmcsolver.h"
#include "HarmonicOscillator.h"


class BasisFunctions;

class Orbitals {
private:
    int N_dim ;
    BasisFunctions** basis_functions;
    BasisFunctions*** dell_basis_functions;
    BasisFunctions** lapl_basis_functions;
    
public:
    Orbitals();

    void setup_basis();
    double get_basis_function(int N_orbitals, const mat &r, int i, int j,  VMCSolver *solver);
    vec get_gradient(const mat &r, int i, int j, VMCSolver *solver);
    double get_laplacian(const mat &r, int i, int j, VMCSolver *solver);
    void destroy();
};


#endif
