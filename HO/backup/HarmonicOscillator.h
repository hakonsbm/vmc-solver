#ifndef HARMONICOSCILLATOR_H
#define HARMONICOSCILLATOR_H 
#include  <armadillo>
#include  "vmcsolver.h"
#include  "BasisFunctions.h"


//Superclass 
class HarmonicOscillator : public BasisFunctions {
protected:
    double* k;
    double* k2;
    double* exp_factor;


public:
    HarmonicOscillator();
    ~HarmonicOscillator();
};


/*
    Subclasses
*/

class HarmonicOscillator_0 : public HarmonicOscillator {
public:

    HarmonicOscillator_0();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_0_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_0_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_0_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_0_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_0 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_0();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 0 -----------------------------

*/
class HarmonicOscillator_1 : public HarmonicOscillator {
public:

    HarmonicOscillator_1();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_1_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_1_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_1_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_1_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_1 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_1();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 1 -----------------------------

*/
class HarmonicOscillator_2 : public HarmonicOscillator {
public:

    HarmonicOscillator_2();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_2_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_2_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_2_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_2_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_2 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_2();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 2 -----------------------------

*/
class HarmonicOscillator_3 : public HarmonicOscillator {
public:

    HarmonicOscillator_3();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_3_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_3_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_3_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_3 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_3();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 3 -----------------------------

*/
class HarmonicOscillator_4 : public HarmonicOscillator {
public:

    HarmonicOscillator_4();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_4_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_4_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_4_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_4_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_4 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_4();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 4 -----------------------------

*/
class HarmonicOscillator_5 : public HarmonicOscillator {
public:

    HarmonicOscillator_5();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_5_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_5_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_5_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_5_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_5 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_5();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 5 -----------------------------

*/
class HarmonicOscillator_6 : public HarmonicOscillator {
public:

    HarmonicOscillator_6();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_6_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_6_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_6_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_6_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_6 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_6();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 6 -----------------------------

*/
class HarmonicOscillator_7 : public HarmonicOscillator {
public:

    HarmonicOscillator_7();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_7_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_7_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_7_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_7_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_7 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_7();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 7 -----------------------------

*/
class HarmonicOscillator_8 : public HarmonicOscillator {
public:

    HarmonicOscillator_8();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_8_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_8_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_8_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_8_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_8 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_8();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 8 -----------------------------

*/
class HarmonicOscillator_9 : public HarmonicOscillator {
public:

    HarmonicOscillator_9();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_9_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_9_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_9_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_9_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_9 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_9();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 9 -----------------------------

*/
class HarmonicOscillator_10 : public HarmonicOscillator {
public:

    HarmonicOscillator_10();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_10_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_10_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_10_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_10_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_10 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_10();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 10 -----------------------------

*/
class HarmonicOscillator_11 : public HarmonicOscillator {
public:

    HarmonicOscillator_11();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_11_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_11_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_11_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_11_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_11 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_11();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 11 -----------------------------

*/
class HarmonicOscillator_12 : public HarmonicOscillator {
public:

    HarmonicOscillator_12();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_12_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_12_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_12_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_12_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_12 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_12();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 12 -----------------------------

*/
class HarmonicOscillator_13 : public HarmonicOscillator {
public:

    HarmonicOscillator_13();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_13_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_13_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_13_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_13_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_13 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_13();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 13 -----------------------------

*/
class HarmonicOscillator_14 : public HarmonicOscillator {
public:

    HarmonicOscillator_14();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_14_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_14_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_14_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_14_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_14 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_14();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 14 -----------------------------

*/
class HarmonicOscillator_15 : public HarmonicOscillator {
public:

    HarmonicOscillator_15();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_15_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_15_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_15_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_15_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_15 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_15();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 15 -----------------------------

*/
class HarmonicOscillator_16 : public HarmonicOscillator {
public:

    HarmonicOscillator_16();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_16_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_16_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_16_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_16_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_16 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_16();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 16 -----------------------------

*/
class HarmonicOscillator_17 : public HarmonicOscillator {
public:

    HarmonicOscillator_17();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_17_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_17_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_17_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_17_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_17 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_17();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 17 -----------------------------

*/
class HarmonicOscillator_18 : public HarmonicOscillator {
public:

    HarmonicOscillator_18();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_18_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_18_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_18_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_18_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_18 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_18();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 18 -----------------------------

*/
class HarmonicOscillator_19 : public HarmonicOscillator {
public:

    HarmonicOscillator_19();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_19_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_19_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_19_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_19_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_19 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_19();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 19 -----------------------------

*/
class HarmonicOscillator_20 : public HarmonicOscillator {
public:

    HarmonicOscillator_20();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_20_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_20_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_20_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_20_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_20 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_20();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 20 -----------------------------

*/
class HarmonicOscillator_21 : public HarmonicOscillator {
public:

    HarmonicOscillator_21();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_21_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_21_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_21_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_21_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_21 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_21();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 21 -----------------------------

*/
class HarmonicOscillator_22 : public HarmonicOscillator {
public:

    HarmonicOscillator_22();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_22_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_22_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_22_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_22_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_22 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_22();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 22 -----------------------------

*/
class HarmonicOscillator_23 : public HarmonicOscillator {
public:

    HarmonicOscillator_23();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_23_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_23_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_23_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_23_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_23 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_23();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 23 -----------------------------

*/
class HarmonicOscillator_24 : public HarmonicOscillator {
public:

    HarmonicOscillator_24();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_24_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_24_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_24_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_24_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_24 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_24();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 24 -----------------------------

*/
class HarmonicOscillator_25 : public HarmonicOscillator {
public:

    HarmonicOscillator_25();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_25_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_25_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_25_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_25_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_25 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_25();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 25 -----------------------------

*/
class HarmonicOscillator_26 : public HarmonicOscillator {
public:

    HarmonicOscillator_26();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_26_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_26_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_26_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_26_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_26 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_26();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 26 -----------------------------

*/
class HarmonicOscillator_27 : public HarmonicOscillator {
public:

    HarmonicOscillator_27();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_27_x : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_27_x();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class dell_HarmonicOscillator_27_y : public HarmonicOscillator {
public:

    dell_HarmonicOscillator_27_y();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

class lapl_HarmonicOscillator_27 : public HarmonicOscillator {
public:

    lapl_HarmonicOscillator_27();
    virtual double eval(double k, double k2, double *exp_factor, vec pos);

};

/*
   ----------------------------- END 27 -----------------------------

*/
#endif
