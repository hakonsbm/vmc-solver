#include "Orbitals_3d.h"

using namespace arma;

Orbitals_3d::Orbitals_3d() {
    int implemented = 84; // number of implemented orbitals
    int N_dim = 3; // number of dimensions
    
    basis_functions = new BasisFunctions*[implemented];
    
    dell_basis_functions = new BasisFunctions**[N_dim];
    for (int i = 0; i < N_dim; ++i) {
        dell_basis_functions[i] = new BasisFunctions*[implemented];
    };
    
    lapl_basis_functions = new BasisFunctions*[implemented];
    
    setup_basis();
}

double Orbitals_3d::get_basis_function(int N_orbitals, const mat &r, int i, int j, VMCSolver *solver)
{
    double k, k2, exp_factor;
    
    vec pos(3);
    pos(0)=r(i,0); pos(1)=r(i,1); pos(2)=r(i,2);
    
    k = sqrt(solver->getAlpha() * solver->getOmega());
    k2 = k*k;
    exp_factor = exp(-0.5*k2*(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)));
    
    return basis_functions[j]->eval(k, k2, &exp_factor, pos);
}
vec Orbitals_3d::get_gradient(const mat &r, int i, int j, VMCSolver *solver)
{
    double k, k2, exp_factor;
    
    vec pos(3), derivative(3);
    pos(0)=r(i,0); pos(1)=r(i,1); pos(2)=r(i,2);
    
    k = sqrt(solver->getAlpha() * solver->getOmega());
    k2 = k*k;
    exp_factor = exp(-0.5*k2*(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)));
    derivative(0) = dell_basis_functions[0][j]->eval(k, k2, &exp_factor, pos);
    derivative(1) = dell_basis_functions[1][j]->eval(k, k2, &exp_factor, pos);
    derivative(2) = dell_basis_functions[2][j]->eval(k, k2, &exp_factor, pos);
    
    return derivative;
}
double Orbitals_3d::get_laplacian(const mat &r, int i, int j, VMCSolver *solver)
{
    double k, k2, exp_factor;
    
    vec pos(3);
    pos(0)=r(i,0); pos(1)=r(i,1); pos(2)=r(i,2);
    
    k = sqrt(solver->getAlpha() * solver->getOmega());
    k2 = k*k;
    exp_factor = exp(-0.5*k2*(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)));
    return lapl_basis_functions[j]->eval(k, k2, &exp_factor, pos);
    
}
    
    
void Orbitals_3d::destroy()
{
    int implemented =  84;
    for (int i = 0; i < implemented; ++i) {
        delete basis_functions[i];
        delete dell_basis_functions[0][i];
        delete dell_basis_functions[1][i];
        delete dell_basis_functions[2][i];
        delete lapl_basis_functions[i];
    }
    delete basis_functions;
    delete dell_basis_functions;
    delete lapl_basis_functions;
}
    
    
void Orbitals_3d::setup_basis()
    /* REMEMBER: Delete for each "new"!!!! */
{
    basis_functions[0] = new HarmonicOscillator_3d_0();
    basis_functions[1] = new HarmonicOscillator_3d_1();
    basis_functions[2] = new HarmonicOscillator_3d_2();
    basis_functions[3] = new HarmonicOscillator_3d_3();
    basis_functions[4] = new HarmonicOscillator_3d_4();
    basis_functions[5] = new HarmonicOscillator_3d_5();
    basis_functions[6] = new HarmonicOscillator_3d_6();
    basis_functions[7] = new HarmonicOscillator_3d_7();
    basis_functions[8] = new HarmonicOscillator_3d_8();
    basis_functions[9] = new HarmonicOscillator_3d_9();
    basis_functions[10] = new HarmonicOscillator_3d_10();
    basis_functions[11] = new HarmonicOscillator_3d_11();
    basis_functions[12] = new HarmonicOscillator_3d_12();
    basis_functions[13] = new HarmonicOscillator_3d_13();
    basis_functions[14] = new HarmonicOscillator_3d_14();
    basis_functions[15] = new HarmonicOscillator_3d_15();
    basis_functions[16] = new HarmonicOscillator_3d_16();
    basis_functions[17] = new HarmonicOscillator_3d_17();
    basis_functions[18] = new HarmonicOscillator_3d_18();
    basis_functions[19] = new HarmonicOscillator_3d_19();
    basis_functions[20] = new HarmonicOscillator_3d_20();
    basis_functions[21] = new HarmonicOscillator_3d_21();
    basis_functions[22] = new HarmonicOscillator_3d_22();
    basis_functions[23] = new HarmonicOscillator_3d_23();
    basis_functions[24] = new HarmonicOscillator_3d_24();
    basis_functions[25] = new HarmonicOscillator_3d_25();
    basis_functions[26] = new HarmonicOscillator_3d_26();
    basis_functions[27] = new HarmonicOscillator_3d_27();
    basis_functions[28] = new HarmonicOscillator_3d_28();
    basis_functions[29] = new HarmonicOscillator_3d_29();
    basis_functions[30] = new HarmonicOscillator_3d_30();
    basis_functions[31] = new HarmonicOscillator_3d_31();
    basis_functions[32] = new HarmonicOscillator_3d_32();
    basis_functions[33] = new HarmonicOscillator_3d_33();
    basis_functions[34] = new HarmonicOscillator_3d_34();
    basis_functions[35] = new HarmonicOscillator_3d_35();
    basis_functions[36] = new HarmonicOscillator_3d_36();
    basis_functions[37] = new HarmonicOscillator_3d_37();
    basis_functions[38] = new HarmonicOscillator_3d_38();
    basis_functions[39] = new HarmonicOscillator_3d_39();
    basis_functions[40] = new HarmonicOscillator_3d_40();
    basis_functions[41] = new HarmonicOscillator_3d_41();
    basis_functions[42] = new HarmonicOscillator_3d_42();
    basis_functions[43] = new HarmonicOscillator_3d_43();
    basis_functions[44] = new HarmonicOscillator_3d_44();
    basis_functions[45] = new HarmonicOscillator_3d_45();
    basis_functions[46] = new HarmonicOscillator_3d_46();
    basis_functions[47] = new HarmonicOscillator_3d_47();
    basis_functions[48] = new HarmonicOscillator_3d_48();
    basis_functions[49] = new HarmonicOscillator_3d_49();
    basis_functions[50] = new HarmonicOscillator_3d_50();
    basis_functions[51] = new HarmonicOscillator_3d_51();
    basis_functions[52] = new HarmonicOscillator_3d_52();
    basis_functions[53] = new HarmonicOscillator_3d_53();
    basis_functions[54] = new HarmonicOscillator_3d_54();
    basis_functions[55] = new HarmonicOscillator_3d_55();
    basis_functions[56] = new HarmonicOscillator_3d_56();
    basis_functions[57] = new HarmonicOscillator_3d_57();
    basis_functions[58] = new HarmonicOscillator_3d_58();
    basis_functions[59] = new HarmonicOscillator_3d_59();
    basis_functions[60] = new HarmonicOscillator_3d_60();
    basis_functions[61] = new HarmonicOscillator_3d_61();
    basis_functions[62] = new HarmonicOscillator_3d_62();
    basis_functions[63] = new HarmonicOscillator_3d_63();
    basis_functions[64] = new HarmonicOscillator_3d_64();
    basis_functions[65] = new HarmonicOscillator_3d_65();
    basis_functions[66] = new HarmonicOscillator_3d_66();
    basis_functions[67] = new HarmonicOscillator_3d_67();
    basis_functions[68] = new HarmonicOscillator_3d_68();
    basis_functions[69] = new HarmonicOscillator_3d_69();
    basis_functions[70] = new HarmonicOscillator_3d_70();
    basis_functions[71] = new HarmonicOscillator_3d_71();
    basis_functions[72] = new HarmonicOscillator_3d_72();
    basis_functions[73] = new HarmonicOscillator_3d_73();
    basis_functions[74] = new HarmonicOscillator_3d_74();
    basis_functions[75] = new HarmonicOscillator_3d_75();
    basis_functions[76] = new HarmonicOscillator_3d_76();
    basis_functions[77] = new HarmonicOscillator_3d_77();
    basis_functions[78] = new HarmonicOscillator_3d_78();
    basis_functions[79] = new HarmonicOscillator_3d_79();
    basis_functions[80] = new HarmonicOscillator_3d_80();
    basis_functions[81] = new HarmonicOscillator_3d_81();
    basis_functions[82] = new HarmonicOscillator_3d_82();
    basis_functions[83] = new HarmonicOscillator_3d_83();

    dell_basis_functions[0][0] = new dell_HarmonicOscillator_3d_0_x();
    dell_basis_functions[1][0] = new dell_HarmonicOscillator_3d_0_y();
    dell_basis_functions[2][0] = new dell_HarmonicOscillator_3d_0_z();
    dell_basis_functions[0][1] = new dell_HarmonicOscillator_3d_1_x();
    dell_basis_functions[1][1] = new dell_HarmonicOscillator_3d_1_y();
    dell_basis_functions[2][1] = new dell_HarmonicOscillator_3d_1_z();
    dell_basis_functions[0][2] = new dell_HarmonicOscillator_3d_2_x();
    dell_basis_functions[1][2] = new dell_HarmonicOscillator_3d_2_y();
    dell_basis_functions[2][2] = new dell_HarmonicOscillator_3d_2_z();
    dell_basis_functions[0][3] = new dell_HarmonicOscillator_3d_3_x();
    dell_basis_functions[1][3] = new dell_HarmonicOscillator_3d_3_y();
    dell_basis_functions[2][3] = new dell_HarmonicOscillator_3d_3_z();
    dell_basis_functions[0][4] = new dell_HarmonicOscillator_3d_4_x();
    dell_basis_functions[1][4] = new dell_HarmonicOscillator_3d_4_y();
    dell_basis_functions[2][4] = new dell_HarmonicOscillator_3d_4_z();
    dell_basis_functions[0][5] = new dell_HarmonicOscillator_3d_5_x();
    dell_basis_functions[1][5] = new dell_HarmonicOscillator_3d_5_y();
    dell_basis_functions[2][5] = new dell_HarmonicOscillator_3d_5_z();
    dell_basis_functions[0][6] = new dell_HarmonicOscillator_3d_6_x();
    dell_basis_functions[1][6] = new dell_HarmonicOscillator_3d_6_y();
    dell_basis_functions[2][6] = new dell_HarmonicOscillator_3d_6_z();
    dell_basis_functions[0][7] = new dell_HarmonicOscillator_3d_7_x();
    dell_basis_functions[1][7] = new dell_HarmonicOscillator_3d_7_y();
    dell_basis_functions[2][7] = new dell_HarmonicOscillator_3d_7_z();
    dell_basis_functions[0][8] = new dell_HarmonicOscillator_3d_8_x();
    dell_basis_functions[1][8] = new dell_HarmonicOscillator_3d_8_y();
    dell_basis_functions[2][8] = new dell_HarmonicOscillator_3d_8_z();
    dell_basis_functions[0][9] = new dell_HarmonicOscillator_3d_9_x();
    dell_basis_functions[1][9] = new dell_HarmonicOscillator_3d_9_y();
    dell_basis_functions[2][9] = new dell_HarmonicOscillator_3d_9_z();
    dell_basis_functions[0][10] = new dell_HarmonicOscillator_3d_10_x();
    dell_basis_functions[1][10] = new dell_HarmonicOscillator_3d_10_y();
    dell_basis_functions[2][10] = new dell_HarmonicOscillator_3d_10_z();
    dell_basis_functions[0][11] = new dell_HarmonicOscillator_3d_11_x();
    dell_basis_functions[1][11] = new dell_HarmonicOscillator_3d_11_y();
    dell_basis_functions[2][11] = new dell_HarmonicOscillator_3d_11_z();
    dell_basis_functions[0][12] = new dell_HarmonicOscillator_3d_12_x();
    dell_basis_functions[1][12] = new dell_HarmonicOscillator_3d_12_y();
    dell_basis_functions[2][12] = new dell_HarmonicOscillator_3d_12_z();
    dell_basis_functions[0][13] = new dell_HarmonicOscillator_3d_13_x();
    dell_basis_functions[1][13] = new dell_HarmonicOscillator_3d_13_y();
    dell_basis_functions[2][13] = new dell_HarmonicOscillator_3d_13_z();
    dell_basis_functions[0][14] = new dell_HarmonicOscillator_3d_14_x();
    dell_basis_functions[1][14] = new dell_HarmonicOscillator_3d_14_y();
    dell_basis_functions[2][14] = new dell_HarmonicOscillator_3d_14_z();
    dell_basis_functions[0][15] = new dell_HarmonicOscillator_3d_15_x();
    dell_basis_functions[1][15] = new dell_HarmonicOscillator_3d_15_y();
    dell_basis_functions[2][15] = new dell_HarmonicOscillator_3d_15_z();
    dell_basis_functions[0][16] = new dell_HarmonicOscillator_3d_16_x();
    dell_basis_functions[1][16] = new dell_HarmonicOscillator_3d_16_y();
    dell_basis_functions[2][16] = new dell_HarmonicOscillator_3d_16_z();
    dell_basis_functions[0][17] = new dell_HarmonicOscillator_3d_17_x();
    dell_basis_functions[1][17] = new dell_HarmonicOscillator_3d_17_y();
    dell_basis_functions[2][17] = new dell_HarmonicOscillator_3d_17_z();
    dell_basis_functions[0][18] = new dell_HarmonicOscillator_3d_18_x();
    dell_basis_functions[1][18] = new dell_HarmonicOscillator_3d_18_y();
    dell_basis_functions[2][18] = new dell_HarmonicOscillator_3d_18_z();
    dell_basis_functions[0][19] = new dell_HarmonicOscillator_3d_19_x();
    dell_basis_functions[1][19] = new dell_HarmonicOscillator_3d_19_y();
    dell_basis_functions[2][19] = new dell_HarmonicOscillator_3d_19_z();
    dell_basis_functions[0][20] = new dell_HarmonicOscillator_3d_20_x();
    dell_basis_functions[1][20] = new dell_HarmonicOscillator_3d_20_y();
    dell_basis_functions[2][20] = new dell_HarmonicOscillator_3d_20_z();
    dell_basis_functions[0][21] = new dell_HarmonicOscillator_3d_21_x();
    dell_basis_functions[1][21] = new dell_HarmonicOscillator_3d_21_y();
    dell_basis_functions[2][21] = new dell_HarmonicOscillator_3d_21_z();
    dell_basis_functions[0][22] = new dell_HarmonicOscillator_3d_22_x();
    dell_basis_functions[1][22] = new dell_HarmonicOscillator_3d_22_y();
    dell_basis_functions[2][22] = new dell_HarmonicOscillator_3d_22_z();
    dell_basis_functions[0][23] = new dell_HarmonicOscillator_3d_23_x();
    dell_basis_functions[1][23] = new dell_HarmonicOscillator_3d_23_y();
    dell_basis_functions[2][23] = new dell_HarmonicOscillator_3d_23_z();
    dell_basis_functions[0][24] = new dell_HarmonicOscillator_3d_24_x();
    dell_basis_functions[1][24] = new dell_HarmonicOscillator_3d_24_y();
    dell_basis_functions[2][24] = new dell_HarmonicOscillator_3d_24_z();
    dell_basis_functions[0][25] = new dell_HarmonicOscillator_3d_25_x();
    dell_basis_functions[1][25] = new dell_HarmonicOscillator_3d_25_y();
    dell_basis_functions[2][25] = new dell_HarmonicOscillator_3d_25_z();
    dell_basis_functions[0][26] = new dell_HarmonicOscillator_3d_26_x();
    dell_basis_functions[1][26] = new dell_HarmonicOscillator_3d_26_y();
    dell_basis_functions[2][26] = new dell_HarmonicOscillator_3d_26_z();
    dell_basis_functions[0][27] = new dell_HarmonicOscillator_3d_27_x();
    dell_basis_functions[1][27] = new dell_HarmonicOscillator_3d_27_y();
    dell_basis_functions[2][27] = new dell_HarmonicOscillator_3d_27_z();
    dell_basis_functions[0][28] = new dell_HarmonicOscillator_3d_28_x();
    dell_basis_functions[1][28] = new dell_HarmonicOscillator_3d_28_y();
    dell_basis_functions[2][28] = new dell_HarmonicOscillator_3d_28_z();
    dell_basis_functions[0][29] = new dell_HarmonicOscillator_3d_29_x();
    dell_basis_functions[1][29] = new dell_HarmonicOscillator_3d_29_y();
    dell_basis_functions[2][29] = new dell_HarmonicOscillator_3d_29_z();
    dell_basis_functions[0][30] = new dell_HarmonicOscillator_3d_30_x();
    dell_basis_functions[1][30] = new dell_HarmonicOscillator_3d_30_y();
    dell_basis_functions[2][30] = new dell_HarmonicOscillator_3d_30_z();
    dell_basis_functions[0][31] = new dell_HarmonicOscillator_3d_31_x();
    dell_basis_functions[1][31] = new dell_HarmonicOscillator_3d_31_y();
    dell_basis_functions[2][31] = new dell_HarmonicOscillator_3d_31_z();
    dell_basis_functions[0][32] = new dell_HarmonicOscillator_3d_32_x();
    dell_basis_functions[1][32] = new dell_HarmonicOscillator_3d_32_y();
    dell_basis_functions[2][32] = new dell_HarmonicOscillator_3d_32_z();
    dell_basis_functions[0][33] = new dell_HarmonicOscillator_3d_33_x();
    dell_basis_functions[1][33] = new dell_HarmonicOscillator_3d_33_y();
    dell_basis_functions[2][33] = new dell_HarmonicOscillator_3d_33_z();
    dell_basis_functions[0][34] = new dell_HarmonicOscillator_3d_34_x();
    dell_basis_functions[1][34] = new dell_HarmonicOscillator_3d_34_y();
    dell_basis_functions[2][34] = new dell_HarmonicOscillator_3d_34_z();
    dell_basis_functions[0][35] = new dell_HarmonicOscillator_3d_35_x();
    dell_basis_functions[1][35] = new dell_HarmonicOscillator_3d_35_y();
    dell_basis_functions[2][35] = new dell_HarmonicOscillator_3d_35_z();
    dell_basis_functions[0][36] = new dell_HarmonicOscillator_3d_36_x();
    dell_basis_functions[1][36] = new dell_HarmonicOscillator_3d_36_y();
    dell_basis_functions[2][36] = new dell_HarmonicOscillator_3d_36_z();
    dell_basis_functions[0][37] = new dell_HarmonicOscillator_3d_37_x();
    dell_basis_functions[1][37] = new dell_HarmonicOscillator_3d_37_y();
    dell_basis_functions[2][37] = new dell_HarmonicOscillator_3d_37_z();
    dell_basis_functions[0][38] = new dell_HarmonicOscillator_3d_38_x();
    dell_basis_functions[1][38] = new dell_HarmonicOscillator_3d_38_y();
    dell_basis_functions[2][38] = new dell_HarmonicOscillator_3d_38_z();
    dell_basis_functions[0][39] = new dell_HarmonicOscillator_3d_39_x();
    dell_basis_functions[1][39] = new dell_HarmonicOscillator_3d_39_y();
    dell_basis_functions[2][39] = new dell_HarmonicOscillator_3d_39_z();
    dell_basis_functions[0][40] = new dell_HarmonicOscillator_3d_40_x();
    dell_basis_functions[1][40] = new dell_HarmonicOscillator_3d_40_y();
    dell_basis_functions[2][40] = new dell_HarmonicOscillator_3d_40_z();
    dell_basis_functions[0][41] = new dell_HarmonicOscillator_3d_41_x();
    dell_basis_functions[1][41] = new dell_HarmonicOscillator_3d_41_y();
    dell_basis_functions[2][41] = new dell_HarmonicOscillator_3d_41_z();
    dell_basis_functions[0][42] = new dell_HarmonicOscillator_3d_42_x();
    dell_basis_functions[1][42] = new dell_HarmonicOscillator_3d_42_y();
    dell_basis_functions[2][42] = new dell_HarmonicOscillator_3d_42_z();
    dell_basis_functions[0][43] = new dell_HarmonicOscillator_3d_43_x();
    dell_basis_functions[1][43] = new dell_HarmonicOscillator_3d_43_y();
    dell_basis_functions[2][43] = new dell_HarmonicOscillator_3d_43_z();
    dell_basis_functions[0][44] = new dell_HarmonicOscillator_3d_44_x();
    dell_basis_functions[1][44] = new dell_HarmonicOscillator_3d_44_y();
    dell_basis_functions[2][44] = new dell_HarmonicOscillator_3d_44_z();
    dell_basis_functions[0][45] = new dell_HarmonicOscillator_3d_45_x();
    dell_basis_functions[1][45] = new dell_HarmonicOscillator_3d_45_y();
    dell_basis_functions[2][45] = new dell_HarmonicOscillator_3d_45_z();
    dell_basis_functions[0][46] = new dell_HarmonicOscillator_3d_46_x();
    dell_basis_functions[1][46] = new dell_HarmonicOscillator_3d_46_y();
    dell_basis_functions[2][46] = new dell_HarmonicOscillator_3d_46_z();
    dell_basis_functions[0][47] = new dell_HarmonicOscillator_3d_47_x();
    dell_basis_functions[1][47] = new dell_HarmonicOscillator_3d_47_y();
    dell_basis_functions[2][47] = new dell_HarmonicOscillator_3d_47_z();
    dell_basis_functions[0][48] = new dell_HarmonicOscillator_3d_48_x();
    dell_basis_functions[1][48] = new dell_HarmonicOscillator_3d_48_y();
    dell_basis_functions[2][48] = new dell_HarmonicOscillator_3d_48_z();
    dell_basis_functions[0][49] = new dell_HarmonicOscillator_3d_49_x();
    dell_basis_functions[1][49] = new dell_HarmonicOscillator_3d_49_y();
    dell_basis_functions[2][49] = new dell_HarmonicOscillator_3d_49_z();
    dell_basis_functions[0][50] = new dell_HarmonicOscillator_3d_50_x();
    dell_basis_functions[1][50] = new dell_HarmonicOscillator_3d_50_y();
    dell_basis_functions[2][50] = new dell_HarmonicOscillator_3d_50_z();
    dell_basis_functions[0][51] = new dell_HarmonicOscillator_3d_51_x();
    dell_basis_functions[1][51] = new dell_HarmonicOscillator_3d_51_y();
    dell_basis_functions[2][51] = new dell_HarmonicOscillator_3d_51_z();
    dell_basis_functions[0][52] = new dell_HarmonicOscillator_3d_52_x();
    dell_basis_functions[1][52] = new dell_HarmonicOscillator_3d_52_y();
    dell_basis_functions[2][52] = new dell_HarmonicOscillator_3d_52_z();
    dell_basis_functions[0][53] = new dell_HarmonicOscillator_3d_53_x();
    dell_basis_functions[1][53] = new dell_HarmonicOscillator_3d_53_y();
    dell_basis_functions[2][53] = new dell_HarmonicOscillator_3d_53_z();
    dell_basis_functions[0][54] = new dell_HarmonicOscillator_3d_54_x();
    dell_basis_functions[1][54] = new dell_HarmonicOscillator_3d_54_y();
    dell_basis_functions[2][54] = new dell_HarmonicOscillator_3d_54_z();
    dell_basis_functions[0][55] = new dell_HarmonicOscillator_3d_55_x();
    dell_basis_functions[1][55] = new dell_HarmonicOscillator_3d_55_y();
    dell_basis_functions[2][55] = new dell_HarmonicOscillator_3d_55_z();
    dell_basis_functions[0][56] = new dell_HarmonicOscillator_3d_56_x();
    dell_basis_functions[1][56] = new dell_HarmonicOscillator_3d_56_y();
    dell_basis_functions[2][56] = new dell_HarmonicOscillator_3d_56_z();
    dell_basis_functions[0][57] = new dell_HarmonicOscillator_3d_57_x();
    dell_basis_functions[1][57] = new dell_HarmonicOscillator_3d_57_y();
    dell_basis_functions[2][57] = new dell_HarmonicOscillator_3d_57_z();
    dell_basis_functions[0][58] = new dell_HarmonicOscillator_3d_58_x();
    dell_basis_functions[1][58] = new dell_HarmonicOscillator_3d_58_y();
    dell_basis_functions[2][58] = new dell_HarmonicOscillator_3d_58_z();
    dell_basis_functions[0][59] = new dell_HarmonicOscillator_3d_59_x();
    dell_basis_functions[1][59] = new dell_HarmonicOscillator_3d_59_y();
    dell_basis_functions[2][59] = new dell_HarmonicOscillator_3d_59_z();
    dell_basis_functions[0][60] = new dell_HarmonicOscillator_3d_60_x();
    dell_basis_functions[1][60] = new dell_HarmonicOscillator_3d_60_y();
    dell_basis_functions[2][60] = new dell_HarmonicOscillator_3d_60_z();
    dell_basis_functions[0][61] = new dell_HarmonicOscillator_3d_61_x();
    dell_basis_functions[1][61] = new dell_HarmonicOscillator_3d_61_y();
    dell_basis_functions[2][61] = new dell_HarmonicOscillator_3d_61_z();
    dell_basis_functions[0][62] = new dell_HarmonicOscillator_3d_62_x();
    dell_basis_functions[1][62] = new dell_HarmonicOscillator_3d_62_y();
    dell_basis_functions[2][62] = new dell_HarmonicOscillator_3d_62_z();
    dell_basis_functions[0][63] = new dell_HarmonicOscillator_3d_63_x();
    dell_basis_functions[1][63] = new dell_HarmonicOscillator_3d_63_y();
    dell_basis_functions[2][63] = new dell_HarmonicOscillator_3d_63_z();
    dell_basis_functions[0][64] = new dell_HarmonicOscillator_3d_64_x();
    dell_basis_functions[1][64] = new dell_HarmonicOscillator_3d_64_y();
    dell_basis_functions[2][64] = new dell_HarmonicOscillator_3d_64_z();
    dell_basis_functions[0][65] = new dell_HarmonicOscillator_3d_65_x();
    dell_basis_functions[1][65] = new dell_HarmonicOscillator_3d_65_y();
    dell_basis_functions[2][65] = new dell_HarmonicOscillator_3d_65_z();
    dell_basis_functions[0][66] = new dell_HarmonicOscillator_3d_66_x();
    dell_basis_functions[1][66] = new dell_HarmonicOscillator_3d_66_y();
    dell_basis_functions[2][66] = new dell_HarmonicOscillator_3d_66_z();
    dell_basis_functions[0][67] = new dell_HarmonicOscillator_3d_67_x();
    dell_basis_functions[1][67] = new dell_HarmonicOscillator_3d_67_y();
    dell_basis_functions[2][67] = new dell_HarmonicOscillator_3d_67_z();
    dell_basis_functions[0][68] = new dell_HarmonicOscillator_3d_68_x();
    dell_basis_functions[1][68] = new dell_HarmonicOscillator_3d_68_y();
    dell_basis_functions[2][68] = new dell_HarmonicOscillator_3d_68_z();
    dell_basis_functions[0][69] = new dell_HarmonicOscillator_3d_69_x();
    dell_basis_functions[1][69] = new dell_HarmonicOscillator_3d_69_y();
    dell_basis_functions[2][69] = new dell_HarmonicOscillator_3d_69_z();
    dell_basis_functions[0][70] = new dell_HarmonicOscillator_3d_70_x();
    dell_basis_functions[1][70] = new dell_HarmonicOscillator_3d_70_y();
    dell_basis_functions[2][70] = new dell_HarmonicOscillator_3d_70_z();
    dell_basis_functions[0][71] = new dell_HarmonicOscillator_3d_71_x();
    dell_basis_functions[1][71] = new dell_HarmonicOscillator_3d_71_y();
    dell_basis_functions[2][71] = new dell_HarmonicOscillator_3d_71_z();
    dell_basis_functions[0][72] = new dell_HarmonicOscillator_3d_72_x();
    dell_basis_functions[1][72] = new dell_HarmonicOscillator_3d_72_y();
    dell_basis_functions[2][72] = new dell_HarmonicOscillator_3d_72_z();
    dell_basis_functions[0][73] = new dell_HarmonicOscillator_3d_73_x();
    dell_basis_functions[1][73] = new dell_HarmonicOscillator_3d_73_y();
    dell_basis_functions[2][73] = new dell_HarmonicOscillator_3d_73_z();
    dell_basis_functions[0][74] = new dell_HarmonicOscillator_3d_74_x();
    dell_basis_functions[1][74] = new dell_HarmonicOscillator_3d_74_y();
    dell_basis_functions[2][74] = new dell_HarmonicOscillator_3d_74_z();
    dell_basis_functions[0][75] = new dell_HarmonicOscillator_3d_75_x();
    dell_basis_functions[1][75] = new dell_HarmonicOscillator_3d_75_y();
    dell_basis_functions[2][75] = new dell_HarmonicOscillator_3d_75_z();
    dell_basis_functions[0][76] = new dell_HarmonicOscillator_3d_76_x();
    dell_basis_functions[1][76] = new dell_HarmonicOscillator_3d_76_y();
    dell_basis_functions[2][76] = new dell_HarmonicOscillator_3d_76_z();
    dell_basis_functions[0][77] = new dell_HarmonicOscillator_3d_77_x();
    dell_basis_functions[1][77] = new dell_HarmonicOscillator_3d_77_y();
    dell_basis_functions[2][77] = new dell_HarmonicOscillator_3d_77_z();
    dell_basis_functions[0][78] = new dell_HarmonicOscillator_3d_78_x();
    dell_basis_functions[1][78] = new dell_HarmonicOscillator_3d_78_y();
    dell_basis_functions[2][78] = new dell_HarmonicOscillator_3d_78_z();
    dell_basis_functions[0][79] = new dell_HarmonicOscillator_3d_79_x();
    dell_basis_functions[1][79] = new dell_HarmonicOscillator_3d_79_y();
    dell_basis_functions[2][79] = new dell_HarmonicOscillator_3d_79_z();
    dell_basis_functions[0][80] = new dell_HarmonicOscillator_3d_80_x();
    dell_basis_functions[1][80] = new dell_HarmonicOscillator_3d_80_y();
    dell_basis_functions[2][80] = new dell_HarmonicOscillator_3d_80_z();
    dell_basis_functions[0][81] = new dell_HarmonicOscillator_3d_81_x();
    dell_basis_functions[1][81] = new dell_HarmonicOscillator_3d_81_y();
    dell_basis_functions[2][81] = new dell_HarmonicOscillator_3d_81_z();
    dell_basis_functions[0][82] = new dell_HarmonicOscillator_3d_82_x();
    dell_basis_functions[1][82] = new dell_HarmonicOscillator_3d_82_y();
    dell_basis_functions[2][82] = new dell_HarmonicOscillator_3d_82_z();
    dell_basis_functions[0][83] = new dell_HarmonicOscillator_3d_83_x();
    dell_basis_functions[1][83] = new dell_HarmonicOscillator_3d_83_y();
    dell_basis_functions[2][83] = new dell_HarmonicOscillator_3d_83_z();

    lapl_basis_functions[0] = new lapl_HarmonicOscillator_3d_0();
    lapl_basis_functions[1] = new lapl_HarmonicOscillator_3d_1();
    lapl_basis_functions[2] = new lapl_HarmonicOscillator_3d_2();
    lapl_basis_functions[3] = new lapl_HarmonicOscillator_3d_3();
    lapl_basis_functions[4] = new lapl_HarmonicOscillator_3d_4();
    lapl_basis_functions[5] = new lapl_HarmonicOscillator_3d_5();
    lapl_basis_functions[6] = new lapl_HarmonicOscillator_3d_6();
    lapl_basis_functions[7] = new lapl_HarmonicOscillator_3d_7();
    lapl_basis_functions[8] = new lapl_HarmonicOscillator_3d_8();
    lapl_basis_functions[9] = new lapl_HarmonicOscillator_3d_9();
    lapl_basis_functions[10] = new lapl_HarmonicOscillator_3d_10();
    lapl_basis_functions[11] = new lapl_HarmonicOscillator_3d_11();
    lapl_basis_functions[12] = new lapl_HarmonicOscillator_3d_12();
    lapl_basis_functions[13] = new lapl_HarmonicOscillator_3d_13();
    lapl_basis_functions[14] = new lapl_HarmonicOscillator_3d_14();
    lapl_basis_functions[15] = new lapl_HarmonicOscillator_3d_15();
    lapl_basis_functions[16] = new lapl_HarmonicOscillator_3d_16();
    lapl_basis_functions[17] = new lapl_HarmonicOscillator_3d_17();
    lapl_basis_functions[18] = new lapl_HarmonicOscillator_3d_18();
    lapl_basis_functions[19] = new lapl_HarmonicOscillator_3d_19();
    lapl_basis_functions[20] = new lapl_HarmonicOscillator_3d_20();
    lapl_basis_functions[21] = new lapl_HarmonicOscillator_3d_21();
    lapl_basis_functions[22] = new lapl_HarmonicOscillator_3d_22();
    lapl_basis_functions[23] = new lapl_HarmonicOscillator_3d_23();
    lapl_basis_functions[24] = new lapl_HarmonicOscillator_3d_24();
    lapl_basis_functions[25] = new lapl_HarmonicOscillator_3d_25();
    lapl_basis_functions[26] = new lapl_HarmonicOscillator_3d_26();
    lapl_basis_functions[27] = new lapl_HarmonicOscillator_3d_27();
    lapl_basis_functions[28] = new lapl_HarmonicOscillator_3d_28();
    lapl_basis_functions[29] = new lapl_HarmonicOscillator_3d_29();
    lapl_basis_functions[30] = new lapl_HarmonicOscillator_3d_30();
    lapl_basis_functions[31] = new lapl_HarmonicOscillator_3d_31();
    lapl_basis_functions[32] = new lapl_HarmonicOscillator_3d_32();
    lapl_basis_functions[33] = new lapl_HarmonicOscillator_3d_33();
    lapl_basis_functions[34] = new lapl_HarmonicOscillator_3d_34();
    lapl_basis_functions[35] = new lapl_HarmonicOscillator_3d_35();
    lapl_basis_functions[36] = new lapl_HarmonicOscillator_3d_36();
    lapl_basis_functions[37] = new lapl_HarmonicOscillator_3d_37();
    lapl_basis_functions[38] = new lapl_HarmonicOscillator_3d_38();
    lapl_basis_functions[39] = new lapl_HarmonicOscillator_3d_39();
    lapl_basis_functions[40] = new lapl_HarmonicOscillator_3d_40();
    lapl_basis_functions[41] = new lapl_HarmonicOscillator_3d_41();
    lapl_basis_functions[42] = new lapl_HarmonicOscillator_3d_42();
    lapl_basis_functions[43] = new lapl_HarmonicOscillator_3d_43();
    lapl_basis_functions[44] = new lapl_HarmonicOscillator_3d_44();
    lapl_basis_functions[45] = new lapl_HarmonicOscillator_3d_45();
    lapl_basis_functions[46] = new lapl_HarmonicOscillator_3d_46();
    lapl_basis_functions[47] = new lapl_HarmonicOscillator_3d_47();
    lapl_basis_functions[48] = new lapl_HarmonicOscillator_3d_48();
    lapl_basis_functions[49] = new lapl_HarmonicOscillator_3d_49();
    lapl_basis_functions[50] = new lapl_HarmonicOscillator_3d_50();
    lapl_basis_functions[51] = new lapl_HarmonicOscillator_3d_51();
    lapl_basis_functions[52] = new lapl_HarmonicOscillator_3d_52();
    lapl_basis_functions[53] = new lapl_HarmonicOscillator_3d_53();
    lapl_basis_functions[54] = new lapl_HarmonicOscillator_3d_54();
    lapl_basis_functions[55] = new lapl_HarmonicOscillator_3d_55();
    lapl_basis_functions[56] = new lapl_HarmonicOscillator_3d_56();
    lapl_basis_functions[57] = new lapl_HarmonicOscillator_3d_57();
    lapl_basis_functions[58] = new lapl_HarmonicOscillator_3d_58();
    lapl_basis_functions[59] = new lapl_HarmonicOscillator_3d_59();
    lapl_basis_functions[60] = new lapl_HarmonicOscillator_3d_60();
    lapl_basis_functions[61] = new lapl_HarmonicOscillator_3d_61();
    lapl_basis_functions[62] = new lapl_HarmonicOscillator_3d_62();
    lapl_basis_functions[63] = new lapl_HarmonicOscillator_3d_63();
    lapl_basis_functions[64] = new lapl_HarmonicOscillator_3d_64();
    lapl_basis_functions[65] = new lapl_HarmonicOscillator_3d_65();
    lapl_basis_functions[66] = new lapl_HarmonicOscillator_3d_66();
    lapl_basis_functions[67] = new lapl_HarmonicOscillator_3d_67();
    lapl_basis_functions[68] = new lapl_HarmonicOscillator_3d_68();
    lapl_basis_functions[69] = new lapl_HarmonicOscillator_3d_69();
    lapl_basis_functions[70] = new lapl_HarmonicOscillator_3d_70();
    lapl_basis_functions[71] = new lapl_HarmonicOscillator_3d_71();
    lapl_basis_functions[72] = new lapl_HarmonicOscillator_3d_72();
    lapl_basis_functions[73] = new lapl_HarmonicOscillator_3d_73();
    lapl_basis_functions[74] = new lapl_HarmonicOscillator_3d_74();
    lapl_basis_functions[75] = new lapl_HarmonicOscillator_3d_75();
    lapl_basis_functions[76] = new lapl_HarmonicOscillator_3d_76();
    lapl_basis_functions[77] = new lapl_HarmonicOscillator_3d_77();
    lapl_basis_functions[78] = new lapl_HarmonicOscillator_3d_78();
    lapl_basis_functions[79] = new lapl_HarmonicOscillator_3d_79();
    lapl_basis_functions[80] = new lapl_HarmonicOscillator_3d_80();
    lapl_basis_functions[81] = new lapl_HarmonicOscillator_3d_81();
    lapl_basis_functions[82] = new lapl_HarmonicOscillator_3d_82();
    lapl_basis_functions[83] = new lapl_HarmonicOscillator_3d_83();

}










