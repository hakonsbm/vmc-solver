#ifndef HARMONICOSCILLATOR_3D_H
#define HARMONICOSCILLATOR_3D_H 
#include  <armadillo>
#include  "vmcsolver.h"
#include  "BasisFunctions.h"


//Superclass 
class HarmonicOscillator_3d : public BasisFunctions {
protected:
    double* k;
    double* k2;
    double* exp_factor;


public:
    HarmonicOscillator_3d();
    ~HarmonicOscillator_3d();
};


/*
    Subclasses
*/

class HarmonicOscillator_3d_0 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_0();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_0_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_0_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_0_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_0_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_0_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_0_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_0 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_0();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 0 -----------------------------

*/
class HarmonicOscillator_3d_1 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_1();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_1_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_1_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_1_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_1_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_1_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_1_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_1 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_1();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 1 -----------------------------

*/
class HarmonicOscillator_3d_2 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_2();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_2_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_2_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_2_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_2_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_2_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_2_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_2 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_2();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 2 -----------------------------

*/
class HarmonicOscillator_3d_3 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_3();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_3_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_3_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_3_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_3_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_3_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_3_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_3 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_3();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 3 -----------------------------

*/
class HarmonicOscillator_3d_4 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_4();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_4_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_4_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_4_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_4_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_4_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_4_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_4 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_4();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 4 -----------------------------

*/
class HarmonicOscillator_3d_5 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_5();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_5_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_5_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_5_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_5_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_5_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_5_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_5 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_5();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 5 -----------------------------

*/
class HarmonicOscillator_3d_6 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_6();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_6_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_6_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_6_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_6_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_6_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_6_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_6 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_6();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 6 -----------------------------

*/
class HarmonicOscillator_3d_7 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_7();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_7_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_7_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_7_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_7_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_7_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_7_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_7 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_7();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 7 -----------------------------

*/
class HarmonicOscillator_3d_8 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_8();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_8_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_8_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_8_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_8_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_8_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_8_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_8 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_8();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 8 -----------------------------

*/
class HarmonicOscillator_3d_9 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_9();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_9_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_9_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_9_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_9_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_9_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_9_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_9 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_9();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 9 -----------------------------

*/
class HarmonicOscillator_3d_10 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_10();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_10_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_10_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_10_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_10_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_10_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_10_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_10 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_10();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 10 -----------------------------

*/
class HarmonicOscillator_3d_11 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_11();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_11_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_11_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_11_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_11_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_11_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_11_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_11 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_11();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 11 -----------------------------

*/
class HarmonicOscillator_3d_12 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_12();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_12_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_12_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_12_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_12_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_12_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_12_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_12 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_12();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 12 -----------------------------

*/
class HarmonicOscillator_3d_13 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_13();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_13_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_13_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_13_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_13_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_13_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_13_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_13 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_13();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 13 -----------------------------

*/
class HarmonicOscillator_3d_14 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_14();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_14_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_14_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_14_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_14_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_14_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_14_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_14 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_14();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 14 -----------------------------

*/
class HarmonicOscillator_3d_15 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_15();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_15_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_15_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_15_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_15_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_15_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_15_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_15 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_15();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 15 -----------------------------

*/
class HarmonicOscillator_3d_16 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_16();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_16_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_16_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_16_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_16_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_16_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_16_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_16 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_16();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 16 -----------------------------

*/
class HarmonicOscillator_3d_17 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_17();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_17_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_17_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_17_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_17_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_17_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_17_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_17 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_17();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 17 -----------------------------

*/
class HarmonicOscillator_3d_18 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_18();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_18_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_18_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_18_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_18_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_18_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_18_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_18 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_18();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 18 -----------------------------

*/
class HarmonicOscillator_3d_19 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_19();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_19_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_19_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_19_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_19_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_19_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_19_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_19 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_19();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 19 -----------------------------

*/
class HarmonicOscillator_3d_20 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_20();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_20_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_20_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_20_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_20_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_20_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_20_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_20 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_20();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 20 -----------------------------

*/
class HarmonicOscillator_3d_21 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_21();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_21_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_21_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_21_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_21_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_21_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_21_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_21 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_21();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 21 -----------------------------

*/
class HarmonicOscillator_3d_22 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_22();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_22_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_22_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_22_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_22_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_22_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_22_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_22 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_22();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 22 -----------------------------

*/
class HarmonicOscillator_3d_23 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_23();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_23_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_23_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_23_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_23_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_23_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_23_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_23 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_23();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 23 -----------------------------

*/
class HarmonicOscillator_3d_24 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_24();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_24_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_24_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_24_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_24_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_24_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_24_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_24 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_24();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 24 -----------------------------

*/
class HarmonicOscillator_3d_25 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_25();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_25_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_25_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_25_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_25_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_25_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_25_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_25 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_25();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 25 -----------------------------

*/
class HarmonicOscillator_3d_26 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_26();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_26_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_26_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_26_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_26_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_26_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_26_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_26 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_26();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 26 -----------------------------

*/
class HarmonicOscillator_3d_27 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_27();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_27_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_27_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_27_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_27_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_27_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_27_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_27 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_27();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 27 -----------------------------

*/
class HarmonicOscillator_3d_28 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_28();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_28_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_28_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_28_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_28_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_28_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_28_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_28 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_28();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 28 -----------------------------

*/
class HarmonicOscillator_3d_29 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_29();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_29_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_29_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_29_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_29_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_29_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_29_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_29 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_29();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 29 -----------------------------

*/
class HarmonicOscillator_3d_30 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_30();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_30_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_30_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_30_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_30_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_30_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_30_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_30 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_30();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 30 -----------------------------

*/
class HarmonicOscillator_3d_31 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_31();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_31_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_31_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_31_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_31_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_31_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_31_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_31 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_31();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 31 -----------------------------

*/
class HarmonicOscillator_3d_32 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_32();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_32_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_32_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_32_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_32_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_32_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_32_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_32 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_32();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 32 -----------------------------

*/
class HarmonicOscillator_3d_33 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_33();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_33_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_33_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_33_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_33_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_33_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_33_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_33 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_33();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 33 -----------------------------

*/
class HarmonicOscillator_3d_34 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_34();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_34_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_34_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_34_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_34_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_34_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_34_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_34 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_34();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 34 -----------------------------

*/
class HarmonicOscillator_3d_35 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_35();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_35_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_35_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_35_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_35_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_35_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_35_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_35 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_35();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 35 -----------------------------

*/
class HarmonicOscillator_3d_36 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_36();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_36_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_36_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_36_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_36_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_36_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_36_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_36 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_36();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 36 -----------------------------

*/
class HarmonicOscillator_3d_37 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_37();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_37_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_37_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_37_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_37_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_37_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_37_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_37 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_37();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 37 -----------------------------

*/
class HarmonicOscillator_3d_38 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_38();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_38_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_38_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_38_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_38_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_38_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_38_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_38 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_38();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 38 -----------------------------

*/
class HarmonicOscillator_3d_39 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_39();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_39_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_39_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_39_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_39_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_39_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_39_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_39 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_39();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 39 -----------------------------

*/
class HarmonicOscillator_3d_40 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_40();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_40_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_40_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_40_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_40_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_40_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_40_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_40 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_40();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 40 -----------------------------

*/
class HarmonicOscillator_3d_41 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_41();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_41_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_41_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_41_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_41_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_41_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_41_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_41 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_41();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 41 -----------------------------

*/
class HarmonicOscillator_3d_42 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_42();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_42_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_42_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_42_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_42_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_42_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_42_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_42 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_42();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 42 -----------------------------

*/
class HarmonicOscillator_3d_43 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_43();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_43_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_43_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_43_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_43_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_43_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_43_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_43 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_43();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 43 -----------------------------

*/
class HarmonicOscillator_3d_44 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_44();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_44_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_44_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_44_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_44_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_44_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_44_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_44 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_44();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 44 -----------------------------

*/
class HarmonicOscillator_3d_45 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_45();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_45_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_45_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_45_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_45_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_45_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_45_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_45 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_45();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 45 -----------------------------

*/
class HarmonicOscillator_3d_46 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_46();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_46_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_46_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_46_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_46_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_46_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_46_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_46 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_46();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 46 -----------------------------

*/
class HarmonicOscillator_3d_47 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_47();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_47_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_47_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_47_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_47_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_47_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_47_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_47 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_47();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 47 -----------------------------

*/
class HarmonicOscillator_3d_48 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_48();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_48_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_48_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_48_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_48_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_48_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_48_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_48 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_48();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 48 -----------------------------

*/
class HarmonicOscillator_3d_49 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_49();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_49_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_49_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_49_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_49_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_49_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_49_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_49 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_49();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 49 -----------------------------

*/
class HarmonicOscillator_3d_50 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_50();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_50_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_50_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_50_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_50_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_50_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_50_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_50 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_50();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 50 -----------------------------

*/
class HarmonicOscillator_3d_51 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_51();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_51_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_51_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_51_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_51_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_51_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_51_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_51 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_51();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 51 -----------------------------

*/
class HarmonicOscillator_3d_52 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_52();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_52_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_52_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_52_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_52_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_52_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_52_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_52 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_52();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 52 -----------------------------

*/
class HarmonicOscillator_3d_53 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_53();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_53_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_53_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_53_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_53_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_53_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_53_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_53 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_53();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 53 -----------------------------

*/
class HarmonicOscillator_3d_54 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_54();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_54_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_54_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_54_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_54_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_54_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_54_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_54 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_54();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 54 -----------------------------

*/
class HarmonicOscillator_3d_55 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_55();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_55_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_55_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_55_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_55_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_55_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_55_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_55 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_55();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 55 -----------------------------

*/
class HarmonicOscillator_3d_56 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_56();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_56_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_56_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_56_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_56_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_56_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_56_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_56 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_56();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 56 -----------------------------

*/
class HarmonicOscillator_3d_57 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_57();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_57_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_57_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_57_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_57_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_57_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_57_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_57 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_57();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 57 -----------------------------

*/
class HarmonicOscillator_3d_58 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_58();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_58_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_58_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_58_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_58_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_58_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_58_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_58 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_58();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 58 -----------------------------

*/
class HarmonicOscillator_3d_59 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_59();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_59_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_59_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_59_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_59_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_59_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_59_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_59 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_59();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 59 -----------------------------

*/
class HarmonicOscillator_3d_60 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_60();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_60_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_60_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_60_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_60_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_60_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_60_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_60 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_60();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 60 -----------------------------

*/
class HarmonicOscillator_3d_61 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_61();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_61_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_61_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_61_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_61_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_61_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_61_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_61 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_61();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 61 -----------------------------

*/
class HarmonicOscillator_3d_62 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_62();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_62_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_62_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_62_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_62_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_62_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_62_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_62 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_62();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 62 -----------------------------

*/
class HarmonicOscillator_3d_63 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_63();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_63_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_63_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_63_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_63_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_63_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_63_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_63 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_63();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 63 -----------------------------

*/
class HarmonicOscillator_3d_64 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_64();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_64_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_64_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_64_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_64_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_64_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_64_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_64 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_64();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 64 -----------------------------

*/
class HarmonicOscillator_3d_65 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_65();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_65_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_65_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_65_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_65_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_65_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_65_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_65 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_65();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 65 -----------------------------

*/
class HarmonicOscillator_3d_66 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_66();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_66_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_66_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_66_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_66_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_66_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_66_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_66 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_66();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 66 -----------------------------

*/
class HarmonicOscillator_3d_67 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_67();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_67_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_67_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_67_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_67_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_67_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_67_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_67 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_67();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 67 -----------------------------

*/
class HarmonicOscillator_3d_68 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_68();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_68_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_68_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_68_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_68_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_68_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_68_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_68 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_68();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 68 -----------------------------

*/
class HarmonicOscillator_3d_69 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_69();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_69_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_69_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_69_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_69_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_69_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_69_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_69 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_69();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 69 -----------------------------

*/
class HarmonicOscillator_3d_70 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_70();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_70_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_70_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_70_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_70_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_70_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_70_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_70 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_70();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 70 -----------------------------

*/
class HarmonicOscillator_3d_71 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_71();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_71_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_71_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_71_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_71_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_71_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_71_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_71 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_71();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 71 -----------------------------

*/
class HarmonicOscillator_3d_72 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_72();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_72_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_72_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_72_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_72_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_72_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_72_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_72 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_72();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 72 -----------------------------

*/
class HarmonicOscillator_3d_73 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_73();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_73_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_73_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_73_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_73_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_73_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_73_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_73 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_73();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 73 -----------------------------

*/
class HarmonicOscillator_3d_74 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_74();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_74_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_74_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_74_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_74_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_74_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_74_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_74 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_74();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 74 -----------------------------

*/
class HarmonicOscillator_3d_75 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_75();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_75_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_75_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_75_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_75_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_75_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_75_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_75 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_75();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 75 -----------------------------

*/
class HarmonicOscillator_3d_76 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_76();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_76_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_76_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_76_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_76_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_76_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_76_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_76 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_76();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 76 -----------------------------

*/
class HarmonicOscillator_3d_77 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_77();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_77_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_77_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_77_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_77_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_77_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_77_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_77 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_77();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 77 -----------------------------

*/
class HarmonicOscillator_3d_78 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_78();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_78_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_78_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_78_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_78_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_78_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_78_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_78 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_78();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 78 -----------------------------

*/
class HarmonicOscillator_3d_79 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_79();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_79_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_79_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_79_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_79_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_79_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_79_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_79 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_79();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 79 -----------------------------

*/
class HarmonicOscillator_3d_80 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_80();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_80_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_80_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_80_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_80_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_80_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_80_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_80 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_80();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 80 -----------------------------

*/
class HarmonicOscillator_3d_81 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_81();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_81_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_81_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_81_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_81_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_81_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_81_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_81 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_81();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 81 -----------------------------

*/
class HarmonicOscillator_3d_82 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_82();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_82_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_82_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_82_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_82_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_82_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_82_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_82 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_82();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 82 -----------------------------

*/
class HarmonicOscillator_3d_83 : public HarmonicOscillator_3d {
public:

    HarmonicOscillator_3d_83();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_83_x : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_83_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_83_y : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_83_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3d_83_z : public HarmonicOscillator_3d {
public:

    dell_HarmonicOscillator_3d_83_z();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3d_83 : public HarmonicOscillator_3d {
public:

    lapl_HarmonicOscillator_3d_83();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 83 -----------------------------

*/
#endif
