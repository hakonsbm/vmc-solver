#include "Orbitals.h"

using namespace arma;

Orbitals::Orbitals() {
    int implemented = 28; // number of implemented orbitals
    int N_dim = 2; // number of dimensions
    
    basis_functions = new BasisFunctions*[implemented];
    
    dell_basis_functions = new BasisFunctions**[N_dim];
    for (int i = 0; i < N_dim; ++i) {
        dell_basis_functions[i] = new BasisFunctions*[implemented];
    };
    
    lapl_basis_functions = new BasisFunctions*[implemented];
    
    setup_basis();
}

double Orbitals::get_basis_function(int N_orbitals, const mat &r, int i, int j, VMCSolver *solver)
{
    double k, k2, exp_factor;
    
    vec pos(2);
    pos(0)=r(i,0); pos(1)=r(i,1);
    
    k = sqrt(solver->getAlpha() * solver->getOmega());
    k2 = k*k;
    exp_factor = exp(-0.5*k2*(pos(0)*pos(0) + pos(1)*pos(1)));
    
    return basis_functions[j]->eval(k, k2, &exp_factor, pos);
}
vec Orbitals::get_gradient(const mat &r, int i, int j, VMCSolver *solver)
{
    double k, k2, exp_factor;
    
    vec pos(2), derivative(2);
    pos(0)=r(i,0); pos(1)=r(i,1);
    
    k = sqrt(solver->getAlpha() * solver->getOmega());
    k2 = k*k;
    exp_factor = exp(-0.5*k2*(pos(0)*pos(0) + pos(1)*pos(1)));
    derivative(0) = dell_basis_functions[0][j]->eval(k, k2, &exp_factor, pos);
    derivative(1) = dell_basis_functions[1][j]->eval(k, k2, &exp_factor, pos);
    
    return derivative;
}
double Orbitals::get_laplacian(const mat &r, int i, int j, VMCSolver *solver)
{
    double k, k2, exp_factor;
    
    vec pos(2);
    pos(0)=r(i,0); pos(1)=r(i,1);
    
    k = sqrt(solver->getAlpha() * solver->getOmega());
    k2 = k*k;
    exp_factor = exp(-0.5*k2*(pos(0)*pos(0) + pos(1)*pos(1)));
    return lapl_basis_functions[j]->eval(k, k2, &exp_factor, pos);
    
}
    
    
void Orbitals::destroy()
{
    int implemented =  28;
    for (int i = 0; i < implemented; ++i) {
        delete basis_functions[i];
        delete dell_basis_functions[0][i];
        delete dell_basis_functions[1][i];
        delete lapl_basis_functions[i];
    }
    delete basis_functions;
    delete dell_basis_functions;
    delete lapl_basis_functions;
}
    
    
void Orbitals::setup_basis()
    /* REMEMBER: Delete for each "new"!!!! */
{
    basis_functions[0] = new HarmonicOscillator_0();
    basis_functions[1] = new HarmonicOscillator_1();
    basis_functions[2] = new HarmonicOscillator_2();
    basis_functions[3] = new HarmonicOscillator_3();
    basis_functions[4] = new HarmonicOscillator_4();
    basis_functions[5] = new HarmonicOscillator_5();
    basis_functions[6] = new HarmonicOscillator_6();
    basis_functions[7] = new HarmonicOscillator_7();
    basis_functions[8] = new HarmonicOscillator_8();
    basis_functions[9] = new HarmonicOscillator_9();
    basis_functions[10] = new HarmonicOscillator_10();
    basis_functions[11] = new HarmonicOscillator_11();
    basis_functions[12] = new HarmonicOscillator_12();
    basis_functions[13] = new HarmonicOscillator_13();
    basis_functions[14] = new HarmonicOscillator_14();
    basis_functions[15] = new HarmonicOscillator_15();
    basis_functions[16] = new HarmonicOscillator_16();
    basis_functions[17] = new HarmonicOscillator_17();
    basis_functions[18] = new HarmonicOscillator_18();
    basis_functions[19] = new HarmonicOscillator_19();
    basis_functions[20] = new HarmonicOscillator_20();
    basis_functions[21] = new HarmonicOscillator_21();
    basis_functions[22] = new HarmonicOscillator_22();
    basis_functions[23] = new HarmonicOscillator_23();
    basis_functions[24] = new HarmonicOscillator_24();
    basis_functions[25] = new HarmonicOscillator_25();
    basis_functions[26] = new HarmonicOscillator_26();
    basis_functions[27] = new HarmonicOscillator_27();

    dell_basis_functions[0][0] = new dell_HarmonicOscillator_0_x();
    dell_basis_functions[1][0] = new dell_HarmonicOscillator_0_y();
    dell_basis_functions[0][1] = new dell_HarmonicOscillator_1_x();
    dell_basis_functions[1][1] = new dell_HarmonicOscillator_1_y();
    dell_basis_functions[0][2] = new dell_HarmonicOscillator_2_x();
    dell_basis_functions[1][2] = new dell_HarmonicOscillator_2_y();
    dell_basis_functions[0][3] = new dell_HarmonicOscillator_3_x();
    dell_basis_functions[1][3] = new dell_HarmonicOscillator_3_y();
    dell_basis_functions[0][4] = new dell_HarmonicOscillator_4_x();
    dell_basis_functions[1][4] = new dell_HarmonicOscillator_4_y();
    dell_basis_functions[0][5] = new dell_HarmonicOscillator_5_x();
    dell_basis_functions[1][5] = new dell_HarmonicOscillator_5_y();
    dell_basis_functions[0][6] = new dell_HarmonicOscillator_6_x();
    dell_basis_functions[1][6] = new dell_HarmonicOscillator_6_y();
    dell_basis_functions[0][7] = new dell_HarmonicOscillator_7_x();
    dell_basis_functions[1][7] = new dell_HarmonicOscillator_7_y();
    dell_basis_functions[0][8] = new dell_HarmonicOscillator_8_x();
    dell_basis_functions[1][8] = new dell_HarmonicOscillator_8_y();
    dell_basis_functions[0][9] = new dell_HarmonicOscillator_9_x();
    dell_basis_functions[1][9] = new dell_HarmonicOscillator_9_y();
    dell_basis_functions[0][10] = new dell_HarmonicOscillator_10_x();
    dell_basis_functions[1][10] = new dell_HarmonicOscillator_10_y();
    dell_basis_functions[0][11] = new dell_HarmonicOscillator_11_x();
    dell_basis_functions[1][11] = new dell_HarmonicOscillator_11_y();
    dell_basis_functions[0][12] = new dell_HarmonicOscillator_12_x();
    dell_basis_functions[1][12] = new dell_HarmonicOscillator_12_y();
    dell_basis_functions[0][13] = new dell_HarmonicOscillator_13_x();
    dell_basis_functions[1][13] = new dell_HarmonicOscillator_13_y();
    dell_basis_functions[0][14] = new dell_HarmonicOscillator_14_x();
    dell_basis_functions[1][14] = new dell_HarmonicOscillator_14_y();
    dell_basis_functions[0][15] = new dell_HarmonicOscillator_15_x();
    dell_basis_functions[1][15] = new dell_HarmonicOscillator_15_y();
    dell_basis_functions[0][16] = new dell_HarmonicOscillator_16_x();
    dell_basis_functions[1][16] = new dell_HarmonicOscillator_16_y();
    dell_basis_functions[0][17] = new dell_HarmonicOscillator_17_x();
    dell_basis_functions[1][17] = new dell_HarmonicOscillator_17_y();
    dell_basis_functions[0][18] = new dell_HarmonicOscillator_18_x();
    dell_basis_functions[1][18] = new dell_HarmonicOscillator_18_y();
    dell_basis_functions[0][19] = new dell_HarmonicOscillator_19_x();
    dell_basis_functions[1][19] = new dell_HarmonicOscillator_19_y();
    dell_basis_functions[0][20] = new dell_HarmonicOscillator_20_x();
    dell_basis_functions[1][20] = new dell_HarmonicOscillator_20_y();
    dell_basis_functions[0][21] = new dell_HarmonicOscillator_21_x();
    dell_basis_functions[1][21] = new dell_HarmonicOscillator_21_y();
    dell_basis_functions[0][22] = new dell_HarmonicOscillator_22_x();
    dell_basis_functions[1][22] = new dell_HarmonicOscillator_22_y();
    dell_basis_functions[0][23] = new dell_HarmonicOscillator_23_x();
    dell_basis_functions[1][23] = new dell_HarmonicOscillator_23_y();
    dell_basis_functions[0][24] = new dell_HarmonicOscillator_24_x();
    dell_basis_functions[1][24] = new dell_HarmonicOscillator_24_y();
    dell_basis_functions[0][25] = new dell_HarmonicOscillator_25_x();
    dell_basis_functions[1][25] = new dell_HarmonicOscillator_25_y();
    dell_basis_functions[0][26] = new dell_HarmonicOscillator_26_x();
    dell_basis_functions[1][26] = new dell_HarmonicOscillator_26_y();
    dell_basis_functions[0][27] = new dell_HarmonicOscillator_27_x();
    dell_basis_functions[1][27] = new dell_HarmonicOscillator_27_y();

    lapl_basis_functions[0] = new lapl_HarmonicOscillator_0();
    lapl_basis_functions[1] = new lapl_HarmonicOscillator_1();
    lapl_basis_functions[2] = new lapl_HarmonicOscillator_2();
    lapl_basis_functions[3] = new lapl_HarmonicOscillator_3();
    lapl_basis_functions[4] = new lapl_HarmonicOscillator_4();
    lapl_basis_functions[5] = new lapl_HarmonicOscillator_5();
    lapl_basis_functions[6] = new lapl_HarmonicOscillator_6();
    lapl_basis_functions[7] = new lapl_HarmonicOscillator_7();
    lapl_basis_functions[8] = new lapl_HarmonicOscillator_8();
    lapl_basis_functions[9] = new lapl_HarmonicOscillator_9();
    lapl_basis_functions[10] = new lapl_HarmonicOscillator_10();
    lapl_basis_functions[11] = new lapl_HarmonicOscillator_11();
    lapl_basis_functions[12] = new lapl_HarmonicOscillator_12();
    lapl_basis_functions[13] = new lapl_HarmonicOscillator_13();
    lapl_basis_functions[14] = new lapl_HarmonicOscillator_14();
    lapl_basis_functions[15] = new lapl_HarmonicOscillator_15();
    lapl_basis_functions[16] = new lapl_HarmonicOscillator_16();
    lapl_basis_functions[17] = new lapl_HarmonicOscillator_17();
    lapl_basis_functions[18] = new lapl_HarmonicOscillator_18();
    lapl_basis_functions[19] = new lapl_HarmonicOscillator_19();
    lapl_basis_functions[20] = new lapl_HarmonicOscillator_20();
    lapl_basis_functions[21] = new lapl_HarmonicOscillator_21();
    lapl_basis_functions[22] = new lapl_HarmonicOscillator_22();
    lapl_basis_functions[23] = new lapl_HarmonicOscillator_23();
    lapl_basis_functions[24] = new lapl_HarmonicOscillator_24();
    lapl_basis_functions[25] = new lapl_HarmonicOscillator_25();
    lapl_basis_functions[26] = new lapl_HarmonicOscillator_26();
    lapl_basis_functions[27] = new lapl_HarmonicOscillator_27();

}










