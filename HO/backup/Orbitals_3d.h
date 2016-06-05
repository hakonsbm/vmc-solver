#ifndef ORBITALS_3D_H
#define ORBITALS_3D_H 

#include <armadillo>
#include "vmcsolver.h"
#include "HarmonicOscillator_3d.h"


class BasisFunctions;

class Orbitals_3d {
private:
    int N_dim ;
    BasisFunctions** basis_functions;
    BasisFunctions*** dell_basis_functions;
    BasisFunctions** lapl_basis_functions;
    
public:
    Orbitals_3d();

    void setup_basis();
    double get_basis_function(int N_orbitals, const mat &r, int i, int j,  VMCSolver *solver);
    vec get_gradient(const mat &r, int i, int j, VMCSolver *solver);
    double get_laplacian(const mat &r, int i, int j, VMCSolver *solver);
    void destroy();
};


#endif
