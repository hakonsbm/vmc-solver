#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumjastrowanalytical.h"
#include "derivatives.h"
//#include "slaterdeterminant.h"
#include "slaterdeterminantHO.h"
#include "HO/Orbitals.h"
#include "HO/Orbitals_3d.h"

using namespace arma;

class TrialFunction; class Derivatives; class SlaterDeterminantHO; class Orbitals; class Orbitals_3d;

class VMCSolver
{
public:
    VMCSolver();

    void runMasterIntegration();
    //void runMonteCarloIntegration();
    void runMonteCarloIntegrationIS();

    //void calculateOptimalSteplength();
    void QuantumForce(const mat &r, mat &QForce);
    void setTrialFunction(TrialFunction *trialFunction) { m_trialFunction = trialFunction;}
    void setAlpha(double alpha) {m_alpha = alpha;}
    void setBeta(double beta) {m_beta = beta;}
    void setOmega(double omega) {m_omega = omega;}
    void setCharge(int C) {charge = C;}
    void setStepLength(double inStepLength) {stepLength = inStepLength;}
    void setDimensions(int ndimensions) {nDimensions = ndimensions;}
    void setCycles(double cycles);
    void setNParticles(int P);
    void resetSolver();

    int checknCycles(int nCycles, int numprocs);


    void initiateDerivatives(Derivatives *derivatives) {m_derivatives = derivatives; }
    Derivatives *derivatives(){return m_derivatives;}

    //void initiateSlaterDeterminant(SlaterDeterminant *determinant) {m_slaterDeterminant = determinant;}
    //SlaterDeterminant *determinant() {return m_slaterDeterminant;}

    void initiateHO(Orbitals *HOorbitals) {m_orbitals = HOorbitals;}
    Orbitals *orbitals() {return m_orbitals;}

    void initiateHO_3d(Orbitals_3d *HOorbitals) {m_orbitals_3d = HOorbitals;}
    Orbitals_3d *orbitals_3d() {return m_orbitals_3d;}

    void initiateSlaterDeterminantHO(SlaterDeterminantHO *determinant) {m_slaterDeterminant = determinant;}
    SlaterDeterminantHO *determinant() {return m_slaterDeterminant;}

    TrialFunction *trialFunction(){return m_trialFunction;}

    int getNParticles() {return nParticles; }
    int getNDimensions() {return nDimensions; }
    double getAlpha() {return m_alpha; }
    double getBeta() {return m_beta; }
    double getOmega() {return m_omega; }
    int getCharge() {return charge; }
    int getCycles() {return nCycles; }
    double getH()   {return h;}
    double getH2()  {return h2;}
    double getStepLength() {return stepLength;}
    double getMHR() {return MHRatio;}
    int getMy_Rank() {return my_rank;}
    void switchbBlockSampling(bool onOff) { m_blockSampling = onOff;}
    void switchElectronInteraction(bool onOff) {m_electronInteraction = onOff;}
    void setPrinting(int printcode) {printing = printcode;}
    bool getElectronInteration () {return m_electronInteraction; }
    void switchJastrow(bool onOff) {m_Jastrow = onOff;}
    bool getJastrow() {return m_Jastrow; }
    bool getRank() {return my_rank;}
    double getEnergyDerivative() {return energyDerivative;}



    double getEnergyVar() {return m_energyVar;}
    double getEnergy() {return totalEnergy;}
    double getKin() {return totalKin;}
    double getPot() {return totalPot;}
    void storeEnergy(double energy) {m_energy = energy;}
    void storeVariance(double variance) {m_energyVar = variance;}

    void mpiArguments( int nargs, char* args[]){ m_nargs = nargs; m_args = args; }

    string getTF() {return trialFunction()->m_outfileName;}


private:
    TrialFunction *m_trialFunction;
    Derivatives *m_derivatives;
    Orbitals *m_orbitals;
    Orbitals_3d *m_orbitals_3d;
    //SlaterDeterminant *m_slaterDeterminant;
    SlaterDeterminantHO *m_slaterDeterminant;


    double waveFunction(const mat &r);
    double localEnergy(const mat &r);

    bool importanceSampling;    //When this flag is true it uses importance sampling instead of regular sampling
    bool m_blockSampling = false;

    bool m_electronInteraction;
    bool m_Jastrow;

    int printing = 3;


    int nDimensions;
    int charge;
    double stepLength;
    double step_min;
    double m_alpha;
    double m_beta;
    double m_omega;
    int nParticles;
    double h;
    double h2;
    long idum;
    int nCycles;

    double totalEnergy = 0;
    double totalEnergySquared = 0;
    double totalEnergyVar = 0;
    double totalMoves = 0;
    double totalRatio = 0;
    double totalAverageR12 = 0;
    double totalKin = 0;
    double totalPot = 0;

    //Variables for the MPI implementation
    int my_rank;
    int numprocs;
    int m_nargs;
    char** m_args;

    //Each thread should store data as common accessible doubles, so they can be merged into the the master thread
    double m_energyVar;
    double m_energy;
    double m_energySquared;
    double m_averageR12;
    double m_ratio;
    double m_moves;
    double m_kinSum;
    double m_potSum;

    //
    cube mat_pos;
    vec vec_deltaE;
    vec vec_energySquaredSum;
    vec vec_sqrtR12;

    cube rec_mat_pos;
    vec rec_vec_deltaE;
    vec rec_vec_energySquaredSum;
    vec rec_vec_sqrtR12;


    double D; // diffusion constant
    // double timestep; // timestep for gaussian deviate (using steplength)
    double GreensFunction;
    double MHRatio;
    double energyDerivative = 1;


    mat rOld;
    mat rNew;
    mat QForceOld;
    mat QForceNew;

    mat detUpOld;
    mat detDownOld;

    mat detUpNew;
    mat detDownNew;

};

#endif // VMCSOLVER_H
