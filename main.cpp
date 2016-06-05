
#include "vmcsolver.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumjastrowanalytical.h"
#include "trialFunctions/heliumjastrownumerical.h"
#include "trialFunctions/heliumsimpleanalytical.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "trialFunctions/helium.h"
#include "trialFunctions/hydrogen.h"
#include "trialFunctions/beryllium.h"
#include "trialFunctions/neon.h"
#include "trialFunctions/hydrogentwo.h"
#include "trialFunctions/berylliumtwo.h"
#include "trialFunctions/QuantumDots.h"
#include "slaterdeterminant.h"

#include <iostream>
#include <time.h>
#include <unittest++/UnitTest++.h>
#include <iomanip>
//#include <omp.h>
#include <mpi.h>

using namespace std;
ofstream outfile;
ofstream samplefile;

void testDiffOmega(VMCSolver *solver);
void runFindAlphaBeta(VMCSolver *solver);
void runAlphaBetaGSS(VMCSolver * solver, char *args[]);
void runTimedSingle(VMCSolver * solver, bool Blocking);
void runNoJastrow(VMCSolver * solver, bool Interaction, bool Blocking);
void chooseTrialFunction(string args1, VMCSolver * solver);
/*
void runWithDiffConstants(VMCSolver *solver);
void runSIWithDiffTimesteps(VMCSolver *solver);
void runBlockingSampledRun(VMCSolver *solver);
void runCompareAnalytical(VMCSolver *solver);

void runDiffBetaAndR(VMCSolver *solver);
void runNewtonsMethod(VMCSolver *solver);
int runTests(VMCSolver *solver);
*/

int called; // counter for number of times vmc is called
int Gnargs;
double args7;
string trialfunction;

int main(int nargs, char* args[])
{
/*
 * Usage:
 * > mpirun -n [numprocs] vmc QuantumDots [runtype] [nParticles] [omega] [nDimensions] [nCycles] {OPTIONS}
 *                                1           2          3          4         5           6
 *
 * > mpirun -n [numprocs] vmc QuantumDots runAlphaBetaGSS [nParticles] [omega] [nDimensions] [nCycles] { [eps] { [a_low] [a_high] [b_low] [a_high] } }
 * > mpirun -n [numprocs] vmc QuantumDots runTimedSingle [nParticles] [omega] [nDimensions] [nCycles] {[alpha] [beta]}
 */

    VMCSolver *solver = new VMCSolver();

    MPI_Init(&nargs, &args);
    solver->mpiArguments(nargs, args);

    // Command line argument parsing
    int my_rank, numprocs, nCycles;
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    trialfunction = (string)args[1];
    chooseTrialFunction(trialfunction, solver);
    /*
    if((string)args[1]=="HeliumSimpleAnalytical") solver->setTrialFunction(new HeliumSimpleAnalytical(solver));
    else if((string)args[1]=="HeliumSimpleNumerical") solver->setTrialFunction(new HeliumSimpleNumerical(solver));
    else if((string)args[1]=="HeliumJastrowAnalytical") solver->setTrialFunction(new HeliumJastrowAnalytical(solver));
    else if((string)args[1]=="HeliumJastrowNumerical") solver->setTrialFunction(new HeliumJastrowNumerical(solver));
    else if((string)args[1]=="Beryllium") solver->setTrialFunction(new Beryllium(solver));
    else if((string)args[1]=="Neon") solver->setTrialFunction(new Neon(solver));
    else if((string)args[1]=="Helium") solver->setTrialFunction(new Helium(solver));
    else if((string)args[1]=="HydrogenTwo") solver->setTrialFunction(new HydrogenTwo(solver));
    else if((string)args[1]=="BerylliumTwo") solver->setTrialFunction(new BerylliumTwo(solver));
    else if((string)args[1]=="QuantumDots") solver->setTrialFunction(new QuantumDots(solver));

    else {if(my_rank==0) cout << args[1] << " is not a valid atom" << endl; exit(1);}
    */

    // if you want to use GTOs (remember to turn off analytical solving)
    //solver->determinant()->setGTO(true);


    /*
    int nCycles = atoi(args[3]);
    cout<<"Atom: "<< args[1]<<endl;
    cout<<"Run: "<< args[2]<<endl;
    cout<<"Cycles: "<< args[3]<<endl;
    */

    solver->setCycles(atoi(args[6]));
    solver->setNParticles(atoi(args[3]));
    solver->setOmega(atof(args[4]));
    if(atoi(args[5]) == 2 || atoi(args[5]) == 3){
        solver->setDimensions(atoi(args[5]));
    }else{
        if(my_rank==0) cout << args[5]  << " is not a valid number of dimensions" << endl; exit(1);
    }
    Gnargs = nargs;
    if(nargs == 8){
        args7 = atof(args[7]);
    }
    if(nargs == 9){
        args7 = atof(args[7]);
        solver->setAlpha(atof(args[7]));
        solver->setBeta(atof(args[8]));
    }

    if(my_rank==0)
    {
        cout<<"Atom: "<< args[1]<<endl;
        cout<<"Run: "<< args[2]<<endl;
        cout<<"Cycles: "<< args[6]<<endl;
    }

    //Couldn't set this in the subclass of trialfunction, so putting it here for the time being
    if((string)args[1]=="HydrogenTwo") solver->trialFunction()->setNucleusDistance(1.4);
    else if((string)args[1]=="BerylliumTwo") solver->trialFunction()->setNucleusDistance(4.63);

    if((string)args[2]== "testDiffOmega") testDiffOmega(solver);
    else if((string)args[2]=="runFindAlphaBeta") runFindAlphaBeta(solver);
    else if((string)args[2]=="runAlphaBetaGSS") runAlphaBetaGSS(solver, args);
    else if((string)args[2]=="runTimedSingle") runTimedSingle(solver, true);
    else if((string)args[2]=="runTimedSingleNoBlock") runTimedSingle(solver, false);
    else if((string)args[2]=="runNoJastrow") runNoJastrow(solver, true, false);
    else if((string)args[2]=="runNoJastrowBlocking") runNoJastrow(solver, true, true);
    else if((string)args[2]=="runNoInterBlocking") runNoJastrow(solver, false, true);
    /*
    else if((string)args[2]=="runSIWithDiffTimesteps") runSIWithDiffTimesteps(solver);
    else if((string)args[2]=="runBlockingSampledRun") runBlockingSampledRun(solver);
    else if((string)args[2]=="runCompareAnalytical") runCompareAnalytical(solver);
    else if((string)args[2]=="runWithDiffConstants") runWithDiffConstants(solver) ;

    else if((string)args[2]=="runTests") runTests(solver);
    else if((string)args[2]=="runDiffBetaAndR") runDiffBetaAndR(solver);
    else if((string)args[2]=="runNewtonsMethod") runNewtonsMethod(solver);
    */
    else {if(my_rank==0) cout << args[2]  << " is not a valid runtype" << endl; exit(1);}

    // End MPI
    MPI_Finalize ();
    return 0;
}


// Brute force search to minimize energy with respect to alpha and beta.
void runFindAlphaBeta(VMCSolver *solver)
{

    double alphaMin = 0.5;
    double alphaMax = 3.5;
    double betaMin = 0.5;
    double betaMax = 10;//2.0;

    double bestAlpha = 0;
    double bestBeta = 0;

    int nMeshPoints = 5; // number of points in the search-mesh
    //int nSteps = 6; // number of times to improve accuracy

    vec alpha, beta;

    mat eMesh = zeros<mat>(nMeshPoints, nMeshPoints);
    mat vMesh = zeros<mat>(nMeshPoints, nMeshPoints);
    uword eMinx, eMiny; //

    solver->switchElectronInteraction(true);
    solver->switchJastrow(true);
    solver->trialFunction()->setAnalytical(true);

    solver->setPrinting(0);

    int numprocs;
    double start, end;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    if (solver->getMy_Rank() == 0){
        ostringstream ossP, ossD, ossW;
        ossP << solver->getNParticles();
        ossD << solver->getNDimensions();
        ossW << solver->getOmega();
        string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;
        char const * outfilePath = (pathString + ossP.str() + string("_") + ossD.str() + string("D_w") + ossW.str() + string("_AlphaBetaValues.txt")).c_str();
        outfile.open(outfilePath);
        //outfile << "\t\t e \t\t\t tot. e_sq \t\t var \t\t\t alpha \t\t\t beta \t\t avg dist \t\t\t stepl. \t\t cycles \n";
    }
    start = MPI_Wtime();
    //for (int step = 0; step < nSteps; step++){ // step loop start
    /*
    if (solver->getMy_Rank() == 0){
        cout << "Alpha: " << alphaMin << " to " <<  alphaMax << endl;
        cout << "Beta: " <<  betaMin << " to " <<  betaMax << endl;
    }
    */

    //if (solver->getMy_Rank() == 0 && numprocs<=8){cout << "Step " << step << ":" << endl;}

    // Share boundaries among processes
    MPI_Bcast(&alphaMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alphaMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&betaMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&betaMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    alpha = linspace(alphaMin, alphaMax, nMeshPoints);
    beta = linspace(betaMin, betaMax, nMeshPoints);


    for (double alphaN=0; alphaN < nMeshPoints; alphaN++){
        if (solver->getMy_Rank() == 0 && numprocs<=8){cout << alphaN+1 << "/" << nMeshPoints << endl;}
        for (double betaN=0; betaN < nMeshPoints; betaN++){
            //if (solver->getMy_Rank() == 0 && numprocs<=8){cout << alphaN*betaN + betaN << " / " << nMeshPoints*nMeshPoints << "   \r";}
            // run MC
            solver->setAlpha(alpha(alphaN));
            solver->setBeta(beta(betaN));
            solver->runMasterIntegration();
            MPI_Barrier(MPI_COMM_WORLD);

            if (solver->getMy_Rank() == 0){
                vMesh(alphaN, betaN) = solver->getEnergyVar();
                eMesh(alphaN, betaN) = solver->getEnergy();
            }
        }
    }

    if (solver->getMy_Rank() == 0){

        //cout << endl;
        //cout << alpha << endl;
        //cout << beta << endl;
        //cout << eMesh << endl;
        eMesh.min(eMinx, eMiny); // store indices of minimum-energy
        if (numprocs<=8) {cout << "Minx Miny eMesh(Minx,Miny) = " <<  eMinx << " " << eMiny << " " << eMesh(eMinx,eMiny) << endl;}
        bestAlpha = alpha(eMinx);
        bestBeta = beta(eMiny);
        if (numprocs<=8) {cout << "Alpha Beta = " << bestAlpha << " " << bestBeta << endl << endl;}

        if (alpha(eMinx) > alphaMin + (alphaMax - alphaMin)/2){alphaMin = alphaMin + (alphaMax - alphaMin)/2; }
        else {alphaMax = alphaMin + (alphaMax - alphaMin)/2; }

        if (beta(eMiny) > betaMin + (betaMax - betaMin)/2){ betaMin = betaMin + (betaMax - betaMin)/2; }
        else { betaMax = betaMin + (betaMax - betaMin)/2; }

    }
    //} // step loop end


    end = MPI_Wtime();
    if (solver->getMy_Rank() == 0){
        for (double alphaN=0; alphaN < nMeshPoints; alphaN++){
            for (double betaN=0; betaN < nMeshPoints; betaN++){
                outfile << setw(15) << setprecision(8) << alpha(alphaN);
                outfile << setw(15) << setprecision(8) << beta(betaN);
                outfile << setw(15) << setprecision(8) << eMesh(alphaN, betaN);
                outfile << setw(15) << setprecision(8) << vMesh(alphaN, betaN);
                outfile << endl;
            }
        }


        cout << "Found optimal values:" << endl;
        cout << "Alpha = " << bestAlpha << endl;
        cout << "Beta = " << bestBeta << endl;
        cout << "Time used with "<< numprocs << " processes: " << end - start << endl;

        /*
        outfile << "Particles = " << solver->getNParticles()<< endl;
        outfile << "Dims = " << solver->getNDimensions() << endl;
        outfile << "Cycles = " << solver->getCycles() << endl;

        outfile << "E = " << eMesh(eMinx,eMiny) << endl;
        outfile << "Found optimal values:" << endl;
        outfile << "Alpha = " << bestAlpha << endl;
        outfile << "Beta = " << bestBeta << endl;
        */
        outfile.close();
    }
}

// helper function for gcc
double get_solver_energy(double alpha, double beta, VMCSolver *solver)
{
    double E1, V1;

    //cout << "refresh" << endl;
    /*
    delete solver;
    solver = new VMCSolver();
    chooseTrialFunction(trialfunction, solver);
    solver->setCycles(Ncycles);
    solver->setNParticles(Nparticles);
    solver->setOmega(Omega);
    solver->setDimensions(Dimensions);
    solver->switchElectronInteraction(true);
    solver->trialFunction()->setAnalytical(true);
    solver->setPrinting(0);
    */
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //cout << "done" << endl;
    solver->setAlpha(alpha);
    //cout << "alpha" << endl;
    solver->setBeta(beta);
    //cout << "beta" << endl;
    solver->runMasterIntegration();
    MPI_Barrier(MPI_COMM_WORLD);
    //cout << "master" << endl;
    E1 = solver->getEnergy();
    V1 = solver->getEnergyVar();
    //cout << "complete" << endl;

    //return (E1+E2+E3)/3.0; //solver->getEnergy();
    called++;
    if (solver->getNParticles() >= 20)
    {
        return V1;
    }else{
        return E1;
    }

}

// Golden section search to minimize energy with respect to alpha and beta.
void runAlphaBetaGSS(VMCSolver *solver, char *args[])
{
    solver->switchElectronInteraction(true);
    solver->switchJastrow(true);
    solver->trialFunction()->setAnalytical(true);
    solver->setPrinting(0);

    bool findEnergy = false;
    int printprocs = 4;

    double low_value_x, high_value_x, probe_point_x, probe_value_x, gr_point_x, gr_value_x, x_guess;
    double low_value_y, high_value_y, probe_point_y, probe_value_y, gr_point_y, gr_value_y, y_guess;
    double high_bound_x, low_bound_x, high_bound_y, low_bound_y;
    double prev_x_guess, prev_y_guess;


    double init_low_bound_x = 0.5;
    double init_high_bound_x = 2.0;

    double init_low_bound_y = 0.001;
    double init_high_bound_y = 1.0;

    double gr = (sqrt((double) 5) - 1)/2.0;
    double eps = 0.1;
    //double eps2 = 1e-3;

    double temp_E, start, end;


    int numprocs;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    if (Gnargs == 8){ // set a new accuracy
        if (args7 < 0.1 && args7 > 0){
            eps = args7;
        }
    }
    if (Gnargs == 12){ // set new accuracy and search domain
        eps = atof(args[7]);
        init_low_bound_x = atof(args[8]);
        init_high_bound_x = atof(args[9]);

        init_low_bound_y = atof(args[10]);
        init_high_bound_y = atof(args[11]);
    }


    if (solver->getMy_Rank() == 0){
        ostringstream ossP, ossD, ossW;
        ossP << solver->getNParticles();
        ossD << solver->getNDimensions();
        ossW << solver->getOmega();
        string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;
        char const * outfilePath = (pathString + ossP.str() + string("_") + ossD.str() + string("D_w") + ossW.str() + string("_findAlphaBetaGSS.txt")).c_str();
        outfile.open(outfilePath);
        outfile << "\t\t e \t\t\t tot. e_sq \t\t var \t\t\t alpha \t\t\t beta \t\t avg dist \t\t\t stepl. \t\t cycles \n";
        //cout << "eps = " << eps << endl;
    }


    start = MPI_Wtime();
    //cout << gr << endl;
    //cout << gr_point_x << endl;
    //cout << probe_point_x << endl;
    x_guess = (init_low_bound_x + init_high_bound_x)/2.0;
    y_guess = (init_low_bound_y + init_high_bound_y)/2.0;
    prev_x_guess = x_guess + 100;
    prev_y_guess = y_guess + 100;


    /*
    temp_E = get_solver_energy(1.07548, 0.0270139, probe_distance, solver);
    if (solver->getMy_Rank() == 0){
    cout << "To beat:" << temp_E << endl << endl;
    }
    */

    // find best values
    //while (abs(x_guess + y_guess - (prev_x_guess + prev_y_guess))/2.0 > eps2) {

        low_bound_x = init_low_bound_x;
        high_bound_x = init_high_bound_x;
        low_bound_y = init_low_bound_y;
        high_bound_y = init_high_bound_y;

        gr_point_x = high_bound_x + gr * (low_bound_x - high_bound_x);
        probe_point_x = low_bound_x + gr * (high_bound_x - low_bound_x);
        gr_point_y = high_bound_y + gr * (low_bound_y - high_bound_y);
        probe_point_y = low_bound_y + gr * (high_bound_y - low_bound_y);

        gr_value_x = get_solver_energy(gr_point_x, y_guess, solver);
        probe_value_x = get_solver_energy(probe_point_x, y_guess, solver);



        if (solver->getMy_Rank() == 0 && numprocs <=printprocs){
          cout << "Running x with y_guess:" << y_guess << endl;
        }
        // find best current x-value
        while (abs(gr_point_x - probe_point_x)/2.0 > eps) {
          if (solver->getMy_Rank() == 0 && numprocs <=printprocs){
            cout << "gr_p p_p gr_v p_v " << gr_point_x << "   " << probe_point_x << "   " << gr_value_x << "   " << probe_value_x << endl;
          }
          if (gr_value_x < probe_value_x) { // is lowest point under the probe?
            high_bound_x = probe_point_x;
            probe_point_x = gr_point_x;
            gr_point_x = high_bound_x + gr * (low_bound_x - high_bound_x);
            probe_value_x = gr_value_x;
            gr_value_x = get_solver_energy(gr_point_x, y_guess, solver);
          } else{
              low_bound_x = gr_point_x;
              gr_point_x = probe_point_x;
              probe_point_x = low_bound_x + gr * (high_bound_x - low_bound_x);
              gr_value_x = probe_value_x;
              probe_value_x = get_solver_energy(probe_point_x, y_guess, solver);
          }
        }
        prev_x_guess = x_guess;
        x_guess = (low_bound_x + high_bound_x)/2.0;
        gr_value_y = get_solver_energy(x_guess, gr_point_y, solver);
        probe_value_y = get_solver_energy(x_guess, probe_point_y, solver);

        if (solver->getMy_Rank() == 0 && numprocs <=printprocs){
          cout << "New x_guess " << x_guess << endl;
        }

        //find best y-value
        while (abs(gr_point_y - probe_point_y)/2.0 > eps) {
          if (solver->getMy_Rank() == 0 && numprocs <=printprocs){
            cout << "gr_p p_p gr_v p_v " << gr_point_y << "   " << probe_point_y << "   " << gr_value_y << "   " << probe_value_y << endl;
          }
          if (gr_value_y < probe_value_y) {
            high_bound_y = probe_point_y;
            probe_point_y = gr_point_y;
            gr_point_y = high_bound_y + gr * (low_bound_y - high_bound_y);
            probe_value_y = gr_value_y;
            gr_value_y = get_solver_energy(x_guess, gr_point_y, solver);
          } else{
              low_bound_y = gr_point_y;
              gr_point_y = probe_point_y;
              probe_point_y = low_bound_y + gr * (high_bound_y - low_bound_y);
              gr_value_y = probe_value_y;
              probe_value_y = get_solver_energy(x_guess, probe_point_y, solver);
          }
        }
        prev_y_guess = y_guess;
        y_guess = (low_bound_y + high_bound_y)/2.0;

        if (solver->getMy_Rank() == 0 && numprocs <=printprocs){
            cout << "New y_guess " << y_guess << endl << endl;
        }

        temp_E = get_solver_energy(x_guess, y_guess, solver);
        if (solver->getMy_Rank() == 0 && numprocs <=printprocs){
            cout << "a b E:" << x_guess << " " << y_guess << " " << temp_E << endl;
        }

    //}
    end = MPI_Wtime();
    //solver->resetSolver();
    //cout << solver->getEnergy() << endl;
    //double E = get_solver_energy(x_guess, y_guess, solver);
    if (solver->getMy_Rank() == 0){
        cout << "Nparticles: " << solver->getNParticles() << endl;
        cout << "Omega: " << solver->getOmega() << endl;
        cout << "Minimum at: [" << x_guess << ", " << y_guess << "]" << endl;
        cout << "Energy: " << solver->getEnergy() << endl;
        cout << "Var: " << solver->getEnergyVar() << endl;
        cout << "Kin:" << solver->getKin() << endl;
        cout << "Pot:" << solver->getPot() << endl;
        cout << "vmc called " << called << " times." << endl;
        cout << "Time used with "<< numprocs << " processes: " << end - start << endl << endl;

        outfile << "Particles = " << solver->getNParticles()<< endl;
        outfile << "Dims = " << solver->getNDimensions() << endl;
        outfile << "Cycles = " << solver->getCycles()*numprocs << endl;

        outfile << "E = " << temp_E << endl;
        outfile << "Found optimal values:" << endl;
        outfile << "Alpha = " << x_guess << endl;
        outfile << "Beta = " << y_guess << endl;
        outfile.close();

    }
    if (findEnergy){
        solver->setAlpha(x_guess);
        solver->setBeta(y_guess);
        solver->setCycles(10*solver->getCycles()*numprocs);
        runTimedSingle(solver, false);
    }
}

void testDiffOmega(VMCSolver *solver)
{
    //TestSettings
    solver->switchElectronInteraction(true);
    solver->switchJastrow(true);
    solver->trialFunction()->setAnalytical(true);
    solver->switchbBlockSampling(false);

    solver->setPrinting(1);

    //    solver->setAlpha(solver->getCharge());
    string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;
    char const * outfilePath = (pathString + string("_diff_omega.txt")).c_str();
    outfile.open(outfilePath);
    outfile << "\t\t e \t\t\t tot. e_sq \t\t var \t\t\t alpha \t\t\t beta \t\t avg dist \t\t\t stepl. \t\t cycles \n";

    vec nP(7);
    nP(0) = 2;
    nP(1) = 6;
    nP(2) = 12;
    nP(3) = 20;
    nP(4) = 30;
    nP(5) = 42;
    nP(6) = 56;

    int my_rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    for (int index = 0; index < 1 ; index++){
        if (my_rank == 0){
            cout << "################# " << nP(index) << " PARTICLES #################" << endl;
        }
        solver->setNParticles(nP(index));
        solver->setAlpha(1.0);
        //solver->setAlpha(0.1);
        solver->setOmega(10.0);
        solver->runMasterIntegration();
        if (my_rank == 0){cout << endl;}

        //solver->setAlpha(0.2);
        solver->setOmega(5.0);
        solver->runMasterIntegration();
        if (my_rank == 0){cout << endl;}

        //solver->setAlpha(0.5);
        solver->setOmega(2.0);
        solver->runMasterIntegration();
        if (my_rank == 0){cout << endl;}

        //solver->setAlpha(1.0);
        solver->setOmega(1.0);
        solver->runMasterIntegration();
        if (my_rank == 0){cout << endl;}

        //solver->setAlpha(2.0);
        solver->setOmega(0.5);
        solver->runMasterIntegration();
        if (my_rank == 0){cout << endl;}

        //solver->setAlpha(5.0);
        solver->setOmega(0.2);
        solver->runMasterIntegration();
        if (my_rank == 0){cout << endl;}

        //solver->setAlpha(10.0);
        solver->setOmega(0.1);
        solver->runMasterIntegration();

        if (my_rank == 0){
            cout << endl << endl;
        }
    }
    outfile.close();
    return;
}

void runTimedSingle(VMCSolver * solver, bool Blocking)
{
    double start, end;
    int numprocs;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    //TestSettings
    //bool Blocking = true;

    solver->switchElectronInteraction(true);
    solver->switchJastrow(true);
    solver->trialFunction()->setAnalytical(true);
    solver->switchbBlockSampling(Blocking);

    if (numprocs > 8){
        solver->setPrinting(1);
    }else{
        solver->setPrinting(2);
    }


    //solver->setAlpha(args[7]);
    //solver->setBeta(args[]);

    //    solver->setAlpha(solver->getCharge());

    if (solver->getMy_Rank() == 0){
        //cout << "Blocking = " << Blocking << endl;
        // convert int to string for outfiles
        ostringstream ossP, ossD, ossW;
        ossP << solver->getNParticles();
        ossD << solver->getNDimensions();
        ossW << solver->getOmega();

        //string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;
        // for BIG files:
        string pathString = "/home/hakon/Documents/bigDatafiles/" +  solver->trialFunction()->m_outfileName;
        //cout << "Saving to " <<  pathString + ossP.str() + string("_") + ossD.str() + string("D_w") + ossW.str() << endl;
        char const * outfilePath = (pathString + ossP.str() + string("_") + ossD.str() + string("D_w") + ossW.str() + string("_blocking_final.txt")).c_str();
        outfile.open(outfilePath);
        outfile << "\t\t e \t\t\t tot. e_sq \t\t var \t\t\t alpha \t\t\t beta  \t\t\t w \t\t avg dist \t\t\t stepl. \t\t cycles \t\t totKin \t\t totPot"<<endl;
        if (Blocking){
            char const * samplefilePath = (pathString + ossP.str() + string("_") + ossD.str() + string("D_w") + ossW.str() + string("_blocking_final")).c_str();
            samplefile.open(samplefilePath);
        }
    }



    start = MPI_Wtime();
    solver->runMasterIntegration();
    end = MPI_Wtime();

    if (solver->getMy_Rank()==0)
    {
        cout << "Time used with "<< numprocs << " processes: " << end - start << endl;
        outfile.close();
        if (Blocking){
            samplefile.close();
        }
    }

    return;
}


void runNoJastrow(VMCSolver * solver, bool Interaction, bool Blocking)
{
    double start, end;
    int numprocs;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

    //TestSettings

    solver->switchElectronInteraction(Interaction);
    solver->switchJastrow(false);
    solver->trialFunction()->setAnalytical(true);
    solver->switchbBlockSampling(Blocking);


    if (numprocs > 32){
        solver->setPrinting(1);
    }else{
        solver->setPrinting(2);
    }


    //solver->setAlpha(args[7]);
    //solver->setBeta(args[]);

    //    solver->setAlpha(solver->getCharge());

    if (solver->getMy_Rank() == 0){
        //cout << "Blocking = " << Blocking << endl;
        // convert int to string for outfiles
        ostringstream ossP, ossD, ossW;
        ossP << solver->getNParticles();
        ossD << solver->getNDimensions();
        ossW << solver->getOmega();

        //string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;
        // for BIG files:
        string pathString = "/home/hakon/Documents/bigDatafiles/" +  solver->trialFunction()->m_outfileName;
        //cout << "Saving to " <<  pathString + ossP.str() + string("_") + ossD.str() + string("D_w") + ossW.str() << endl;
        char const * outfilePath = (pathString + ossP.str() + string("_") + ossD.str() + string("D_w") + ossW.str() + string("_noInter.txt")).c_str();
        outfile.open(outfilePath);
        outfile << "\t\t e \t\t\t tot. e_sq \t\t var \t\t\t alpha \t\t\t beta  \t\t\t w \t\t avg dist \t\t\t stepl. \t\t cycles \t\t totKin \t\t totPot"<<endl;
        if (Blocking){
            char const * samplefilePath = (pathString + ossP.str() + string("_") + ossD.str() + string("D_w") + ossW.str() + string("_blocking_noInter")).c_str();
            samplefile.open(samplefilePath);
        }
    }



    start = MPI_Wtime();
    solver->runMasterIntegration();
    end = MPI_Wtime();

    if (solver->getMy_Rank()==0)
    {
        cout << "Time used with "<< numprocs << " processes: " << end - start << endl;
        outfile.close();
    }
    if (Blocking){
        samplefile.close();
    }

    return;
}




void chooseTrialFunction(string args1, VMCSolver *solver)
{
    int my_rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    if(args1=="HeliumSimpleAnalytical") solver->setTrialFunction(new HeliumSimpleAnalytical(solver));
    else if(args1=="HeliumSimpleNumerical") solver->setTrialFunction(new HeliumSimpleNumerical(solver));
    else if(args1=="HeliumJastrowAnalytical") solver->setTrialFunction(new HeliumJastrowAnalytical(solver));
    else if(args1=="HeliumJastrowNumerical") solver->setTrialFunction(new HeliumJastrowNumerical(solver));
    else if(args1=="Beryllium") solver->setTrialFunction(new Beryllium(solver));
    else if(args1=="Neon") solver->setTrialFunction(new Neon(solver));
    else if(args1=="Helium") solver->setTrialFunction(new Helium(solver));
    else if(args1=="HydrogenTwo") solver->setTrialFunction(new HydrogenTwo(solver));
    else if(args1=="BerylliumTwo") solver->setTrialFunction(new BerylliumTwo(solver));
    else if(args1=="QuantumDots") solver->setTrialFunction(new QuantumDots(solver));

    else {if(my_rank==0) cout << args1 << " is not a valid atom" << endl; exit(1);}
}

/*
void runDiffBetaAndR(VMCSolver *solver)
{
    //This is for use with the GTO orbitals where alpha is already given, and for two atom cores molecules where we will use the variation of the distance
    // between the nucleuses as a variational parameter along with beta.

    //Can only be run with HydrogenTwo and BerylliumTwo

    solver->trialFunction()->setNucleusDistance(1.4);

    solver->runMasterIntegration();


}

void runWithDiffConstants(VMCSolver *solver)
{
    //Settings for which values it should be cycled over and if we want to use importance sampling or now
    int my_rank;
    //MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    double alpha_min = 1.3; //0.7*solver->getCharge();
    double alpha_max = 2.4; //1.5* solver->getCharge();

    int nSteps = 45;

    double beta_min = 0.01;
    double beta_max = 1.2;
    double d_alpha = (alpha_max-alpha_min)/ (double) nSteps;
    double d_beta = (beta_max-beta_min)/ (double) nSteps;

    double stepcount;

    solver->switchElectronInteraction(true);
    solver->trialFunction()->setAnalytical(false);

    bool ImportanceSampling = true;    //Set to true if you want to run with importance sampling
    solver->switchbBlockSampling(false);
//    solver->setCycles(1000000);

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_alpha_beta")).c_str();


    outfile.open(outfilePath);
//    samplefile.open(samplefilePath);


    clock_t start, end, tot_start, tot_end;     //To keep track of the time
    tot_start = clock();
    for(double alpha = alpha_min ; alpha <= alpha_max; alpha += d_alpha) {
        solver->setAlpha(alpha);
        if(solver->trialFunction()->simpleFlag) {
            if(ImportanceSampling)
            {
                start = clock();
                    solver->runMasterIntegration();
                end = clock();

                double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;
                cout << "Time to run Monte Carlo: " << timeRunMonte << endl;

            }
            else {
                start = clock();
                    solver->calculateOptimalSteplength();
//                      solver->setStepLength(1.4);
                end = clock();

                double timeOptimalStepLength = 1.0*(end - start)/CLOCKS_PER_SEC;

                start = clock();
                    solver->runMonteCarloIntegration();
                end = clock();

                double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;
                cout << "Time to find Optimal Steplength: " << timeOptimalStepLength << endl;
                cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
            }
        }
        else {
            for(double beta = beta_min ; beta <= beta_max; beta += d_beta) {
                stepcount += 1;
                solver->setBeta(beta);
                if(ImportanceSampling)
                {
                    start = clock();
                    solver->runMasterIntegration();
                    end = clock();

                    if (my_rank == 0)
                    {
                        tot_end = clock();
                        int seconds, minutes, hours;
                        double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;
                        double time_so_far = (tot_end - tot_start)/CLOCKS_PER_SEC;
                        double tot_time = time_so_far / (stepcount/(nSteps*nSteps));
                        int time_remaining = round(tot_time - time_so_far);
                        minutes = time_remaining / 60;
                        seconds = time_remaining % 60;
                        hours = minutes / 60;
                        minutes = minutes % 60;
                        cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
                        cout << "Estimated time remaining: " << hours << "h " << minutes << "m " << seconds << "s" << endl;
                    }
                }
                else {
                    start = clock();
                        solver->calculateOptimalSteplength();
                    end = clock();

                    double timeOptimalStepLength = 1.0*(end - start)/CLOCKS_PER_SEC;

                    start = clock();
                        solver->runMonteCarloIntegration();
                    end = clock();

                    double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;

                    cout << "Time to find Optimal Steplength: " << timeOptimalStepLength << endl;
                    cout << "Time to run Monte Carlo: " << timeRunMonte << endl;
                }
            }
        }
    }

    cout << "\nWriting to " << outfilePath << endl;
    outfile.close();
//    samplefile.close();
}

void runSIWithDiffTimesteps(VMCSolver *solver)
{
    solver->switchbBlockSampling(false);
//    solver->setCycles(1000000);

    int nSteps = 100;
    double time_min = 0.01;
    double time_max = 1.;
    double dt = (time_max-time_min)/ (double) nSteps;

    double timeStep;

    solver->setAlpha(4);
    solver->setBeta(0.31);

    solver->switchbBlockSampling(false);    //This also samples the energies at each cycle to do blocking analysis on the data

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_timeStep")).c_str();

    outfile.open(outfilePath);

    clock_t start, end;     //To keep track of the time


    for(timeStep = time_min ; timeStep < time_max ; timeStep += dt )
    {

        solver->setStepLength(timeStep);

        start = clock();
            solver->runMonteCarloIntegrationIS();
        end = clock();

        double timeRunMonte= 1.0*(end - start)/CLOCKS_PER_SEC;

        cout << "Time to run Monte Carlo: " << timeRunMonte << endl;

    }

    outfile.close();
    cout << "\nWriting to " << outfilePath << endl;

}

void runBlockingSampledRun(VMCSolver *solver)
{
    solver->switchbBlockSampling(true);
    solver->switchElectronInteraction(true);
    solver->trialFunction()->simpleFlag = false;
//    solver->setAlpha(solver->getCharge());
    solver->trialFunction()->setAnalytical(false);



    string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_blockingSamples")).c_str();



    samplefile.open(outfilePath);

    solver->runMasterIntegration();

    samplefile.close();

}

void runCompareAnalytical(VMCSolver *solver)
{
    double timeRunAnalytic;
    double timeRunNumerical;


    solver->switchbBlockSampling(false);
    solver->trialFunction()->setAnalytical(true);
//    solver->setCycles(10000000);

    clock_t start, end;     //To keep track of the time

    start = clock();
        solver->runMonteCarloIntegration();
    end = clock();

    timeRunAnalytic = 1.0*(end - start)/CLOCKS_PER_SEC;

    //solver->setTrialFunction(new HeliumJastrowNumerical(solver));
    solver->trialFunction()->setAnalytical(false);

    solver->trialFunction();

    start = clock();
        solver->runMonteCarloIntegration();
    end = clock();

    timeRunNumerical = 1.0*(end - start)/CLOCKS_PER_SEC;

    cout << "Time to calculate analytic vs numerical " << timeRunAnalytic << " vs " << timeRunNumerical << endl;
    cout << "Time run gain "  <<  (timeRunAnalytic - timeRunNumerical) / timeRunNumerical << endl;
    cout << "Time ratia " << timeRunNumerical/timeRunAnalytic << endl;

}


void runNewtonsMethod(VMCSolver *solver)
{

    solver->trialFunction()->setAnalytical(false);
    solver->trialFunction()->setConjugate(false);

    solver->switchbBlockSampling(false);    //This also samples the energies at each cycle to do blocking analysis on the data

    //Opens the file that the relevant wavefunction should be written to, this file is then written to in the
    //vmcSolver class
    string pathString = "../outfiles/" +  solver->trialFunction()->m_outfileName;

    char const * outfilePath = (pathString + string("_conjugate")).c_str();

    outfile.open(outfilePath);

    //Assuming that alpha is known, or GTO trial functions beta is what we want to minimize E_L against. So we guess a value and goes from there

    int steps = 0;
    double lowerEnd = 0;
    double lowerDerivative = 0;
    double midPoint = 0.25;
    double midDerivative = 0;
    double higherEnd = 0.5;
    double higherDerivative = 0;


    //Helium
    //double lowerEnd = 0;
    //double midPoint = 0.5;
    //double higherEnd = 1;
    //Beryllium
    //double lowerEnd = 0;
    //double midPoint = 0.25;
    //double higherEnd = 0.5;
    //Neon (setAnalytical(false))
    //double lowerEnd = 0.07;
    //double midPoint = 0.1;
    //double higherEnd = 0.13;

    solver->setBeta(lowerEnd);
    solver->runMasterIntegration();
    lowerDerivative = solver->getEnergyDerivative();
    solver->setBeta(higherEnd);
    solver->runMasterIntegration();
    higherDerivative = solver->getEnergyDerivative();
    if(higherDerivative/abs(higherDerivative) == -lowerDerivative/abs(lowerDerivative)) {
        solver->setBeta(midPoint);
        while(((higherEnd-lowerEnd)*0.5 > 0.01) && (steps <= 100)) {
            solver->runMasterIntegration();
            midDerivative = solver->getEnergyDerivative();
            if(midDerivative/abs(midDerivative) == lowerDerivative/abs(lowerDerivative)) {
                lowerEnd = midPoint;
                midPoint += (higherEnd-lowerEnd)*0.5;
            }
            else if(midDerivative/abs(midDerivative) == higherDerivative/abs(higherDerivative)) {
                higherEnd = midPoint;
                midPoint -= (higherEnd-lowerEnd)*0.5;
            }
            solver->setBeta(midPoint);
            steps++;
        }
    }
    else {
        cout << "No roots can be guaranteed to be found in the interval between " << lowerEnd << " and " << higherEnd << endl;
    }




//    double beta = 0.5; //This is the initial seed
//    double h= 0.001;
//    double firstDerivative = 0;

//    solver->setBeta(beta);
//    while((abs(solver->getEnergyDerivative()) > 0.01) && (steps <= 100)) {
//        solver->runMasterIntegration();
//        firstDerivative = solver->getEnergyDerivative();
//        solver->setBeta(beta+h);
//        solver->runMasterIntegration();
//        beta += (firstDerivative*h)/(solver->getEnergyDerivative()-firstDerivative);
//        solver->setBeta(beta);
//        cout << "D=" << firstDerivative << endl;
//        steps++;
//   }


    if(solver->getRank() == 0)
    {
        cout << "\nWriting to " << outfilePath << endl;
    }

    outfile.close();

return;
}

int runTests(VMCSolver *solver)
{
    return  UnitTest::RunAllTests();
}
*/
