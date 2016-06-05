#include "vmcsolver.h"
#include "lib.h"
#include "trialFunctions/trialfunction.h"
#include "trialFunctions/heliumsimplenumerical.h"
#include "HO/Orbitals.h"
#include "HO/Orbitals_3d.h"

#include <armadillo>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <mpi.h>


using namespace arma;
using namespace std;
//ofstream outfile;

extern ofstream outfile;
extern ofstream samplefile;


VMCSolver::VMCSolver():
    nDimensions(2),
    charge(2),
    stepLength(0.005),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(clock()),
    nCycles(10000),
    D(0.5),
    my_rank(0)

{
    switchElectronInteraction(true);
    initiateDerivatives(new Derivatives);
    initiateHO(new Orbitals);
    initiateHO_3d(new Orbitals_3d);
    initiateSlaterDeterminantHO(new SlaterDeterminantHO);

    //determinant()->setGTO(false);
}



void VMCSolver::runMasterIntegration()
{


    //MPI initializations
    int numprocs;
    //        MPI_Init(&m_nargs, &m_args);

    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    //        cout << "Hello world, I have rank " << my_rank << " out of "
    //        << numprocs << endl;
    totalEnergy = 0;
    totalEnergySquared = 0;
    totalEnergyVar = 0;
    totalMoves = 0;
    totalRatio = 0;
    totalAverageR12 = 0;
    totalKin = 0;
    totalPot = 0;

    //nCycles = nCycles/numprocs;


    runMonteCarloIntegrationIS();

    MPI_Reduce(&m_energy, &totalEnergy, 1, MPI_DOUBLE,
                  MPI_SUM, 0 , MPI_COMM_WORLD);
    MPI_Reduce(&m_energyVar, &totalEnergyVar, 1, MPI_DOUBLE,
                MPI_SUM, 0 , MPI_COMM_WORLD);
    MPI_Reduce(&m_energySquared, &totalEnergySquared, 1, MPI_DOUBLE,
               MPI_SUM, 0 , MPI_COMM_WORLD);
    MPI_Reduce(&m_moves, &totalMoves, 1, MPI_DOUBLE,
               MPI_SUM, 0 , MPI_COMM_WORLD);
    MPI_Reduce(&m_ratio, &totalRatio, 1, MPI_DOUBLE,
               MPI_SUM, 0 , MPI_COMM_WORLD);
    MPI_Reduce(&m_averageR12, &totalAverageR12, 1, MPI_DOUBLE,
               MPI_SUM, 0 , MPI_COMM_WORLD);

    MPI_Reduce(&m_kinSum, &totalKin, 1, MPI_DOUBLE,
               MPI_SUM, 0 , MPI_COMM_WORLD);
    MPI_Reduce(&m_potSum, &totalPot, 1, MPI_DOUBLE,
               MPI_SUM, 0 , MPI_COMM_WORLD);

    //A lot of the data should be averaged over the results from the different threads

    totalEnergy /=numprocs;
    totalEnergyVar /= numprocs;
    totalEnergySquared /= numprocs;
    totalAverageR12 /= numprocs;
    totalRatio /= numprocs;

    totalKin /= numprocs;
    totalPot /= numprocs;


    if (m_blockSampling){
        // gather cycle-data from all processes
        int blocknCycles = checknCycles(nCycles, numprocs);

        if (my_rank == 0){
            char samplebuf[100];
            char posbuf[200];
            string samplestring = "";
            if (printing>=2){
                cout << "Writing cycle-data" << endl;
                cout << "p0" << endl;
            }
            rec_mat_pos = zeros<cube>(blocknCycles, nParticles, nDimensions);
            rec_vec_deltaE = zeros<vec>(blocknCycles);
            rec_vec_energySquaredSum = zeros<vec>(blocknCycles);
            rec_vec_sqrtR12 = zeros<vec>(blocknCycles);

            for (int i = 0; i < blocknCycles; i++){
                // for the master process
                //char samplestring [100];
                sprintf(samplebuf, "%15f %15f %15f", vec_deltaE(i), vec_energySquaredSum(i), vec_sqrtR12(i));
                samplestring = samplebuf;
                //samplefile << setw(15) << setprecision(8) << vec_deltaE(i);
                //samplefile << setw(15) << setprecision(8) << vec_energySquaredSum(i);
                //samplefile << setw(15) << setprecision(8) << vec_sqrtR12(i);

                for(int n = 0; n < nParticles; n++){
                    if (nDimensions == 2){
                        sprintf(posbuf, " %15f %15f ", mat_pos(i, n, 0), mat_pos(i, n, 1));
                    }else if (nDimensions == 3){
                        sprintf(posbuf, " %15f %15f %15f ", mat_pos(i, n, 0), mat_pos(i, n, 1), mat_pos(i, n, 2));
                    }
                    //samplefile << setw(15) << setprecision(8) << mat_pos(i, n, 0);
                    //samplefile << setw(15) << setprecision(8) << mat_pos(i, n, 1);
                    samplestring += posbuf;
                }
                samplestring += "\n";
                //strcat(samplestring, posstring);
                samplefile << samplestring;
                //samplefile << endl;
            }
            // for the other processes
            //cout << "Rec" << endl;
            for(int p = 1; p < numprocs ; p++) {
                if (printing>=2){
                    cout << "p" << p << endl;
                }

                MPI_Recv(rec_vec_deltaE.begin(), blocknCycles, MPI_DOUBLE, p, 990, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(rec_vec_energySquaredSum.begin(), blocknCycles, MPI_DOUBLE, p, 991, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(rec_vec_sqrtR12.begin(), blocknCycles, MPI_DOUBLE, p, 992, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(rec_mat_pos.begin(), blocknCycles*nParticles*nDimensions, MPI_DOUBLE, p, 993, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //cout << "Process " << p << " done" << endl;
                for (int i = 0; i < blocknCycles; i++){
                    sprintf(samplebuf, "%15f %15f %15f", rec_vec_deltaE(i), rec_vec_energySquaredSum(i), rec_vec_sqrtR12(i));
                    samplestring = samplebuf;
                    //samplefile << setw(15) << setprecision(8) << rec_vec_deltaE(i);
                    //samplefile << setw(15) << setprecision(8) << rec_vec_energySquaredSum(i);
                    //samplefile << setw(15) << setprecision(8) << rec_vec_sqrtR12(i);
                    for(int n = 0; n < nParticles; n++){
                        if (nDimensions == 2){
                            sprintf(posbuf, " %15f %15f ", rec_mat_pos(i, n, 0), rec_mat_pos(i, n, 1));
                        }else if (nDimensions == 3){
                            sprintf(posbuf, " %15f %15f %15f ", rec_mat_pos(i, n, 0), rec_mat_pos(i, n, 1), rec_mat_pos(i, n, 2));
                        }
                        samplestring += posbuf;
                        //samplefile << setw(15) << setprecision(8) << rec_mat_pos(i, n, 0);
                        //samplefile << setw(15) << setprecision(8) << rec_mat_pos(i, n, 1);
                    }
                    samplestring += "\n";
                    //strcat(samplestring, posstring);
                    samplefile << samplestring;
                    //samplefile << endl;
                }
            }
            //cout << "Write" << endl;



                for(int p = 1; p < numprocs; p++){

                }


        }else{
            MPI_Send(vec_deltaE.begin(), vec_deltaE.size(), MPI_DOUBLE, 0, 990, MPI_COMM_WORLD);
            MPI_Send(vec_energySquaredSum.begin(), vec_energySquaredSum.size(), MPI_DOUBLE, 0, 991, MPI_COMM_WORLD);
            MPI_Send(vec_sqrtR12.begin(), vec_sqrtR12.size(), MPI_DOUBLE, 0, 992, MPI_COMM_WORLD);
            MPI_Send(mat_pos.begin(), mat_pos.size(), MPI_DOUBLE, 0, 993, MPI_COMM_WORLD);

        }
    }




    if(my_rank == 0)
    {
        if (printing >= 1){
            cout << "Energy: " << totalEnergy << " Energy (squared sum): " << totalEnergySquared << endl;
            cout << "Variance: " << totalEnergyVar << endl;
            cout << "Alpha: " << m_alpha << " and beta: " << m_beta << endl;
            cout << "Omega: " << m_omega << endl;

            //Write results to file
            outfile << setw(15) << setprecision(8) << totalEnergy;
            outfile << setw(15) << setprecision(8) << totalEnergySquared;
            outfile << setw(15) << setprecision(8) << totalEnergyVar;
            outfile << setw(15) << setprecision(8) << m_alpha;
            outfile << setw(15) << setprecision(8) << m_beta;
            outfile << setw(15) << setprecision(8) << m_omega;
            outfile << setw(15) << setprecision(8) << totalAverageR12;
            outfile << setw(15) << setprecision(8) << stepLength;
            outfile << setw(15) << setprecision(8) << nCycles*numprocs ;
            outfile << trialFunction()->getNucleusDistance();
            outfile << setw(15) << setprecision(8) << totalKin;
            outfile << setw(15) << setprecision(8) << totalPot << endl;
        }
        if (printing >= 2){
            cout << "Kin E: " << totalKin << endl;
            cout << "Pot E: " << totalPot << endl;
            cout << "Moves: " << totalMoves << endl;
            cout << "Ratio: " <<  totalRatio << endl;
            cout << "nParticles: " <<  nParticles << endl;
            cout << "nDimensions: " <<  nDimensions << endl;
            cout << "Average distance between the electrons: " << totalAverageR12 << endl;


            //cout << determinant()->updateSlaterCounter << endl;
            //cout << determinant()->updateDeterminantCounter << endl;
            //cout << "Steplength: " << stepLength << endl;
            //cout << "Nuclei Distance: " << trialFunction()->getNucleusDistance()  << endl;
        }
    }
}



void VMCSolver::runMonteCarloIntegrationIS() {

    //Make sure that the threads have different seeds
    idum = clock();
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    int acc_moves = 0;
    int moves = 0;
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE = 0;
    double potSum = 0;
    double kinSum = 0;
    double r12 = 0;
    double averageR12 = 0;
    double deltaLnDerivative = 0;
    double deltaLnSecondDerivative = 0;
    double lnDerivative = 0;
    double lnSecondDerivative = 0;
    double lnDerivativeLocalEnergy = 0;
    double lnSecondDerivativeLocalEnergy = 0;
    double normalizationFactor = 0;

    //double rec_deltaE, rec_energySquaredSum, rec_sqrtR12, rec_x, rec_y, sqrtR12;
    //char blocking_buffer [250], pos_buffer [50], pos2_buffer [50];
    int blocknCycles = checknCycles(nCycles, numprocs);
    bool reduced_blocking = (blocknCycles != nCycles);

    if (my_rank == 0){
        //cout << blocknCycles << endl;
        //cout << reduced_blocking << endl;
    }

    if (m_blockSampling){
        // vectors and matrices for saving cycle-data
        mat_pos = zeros<cube>(blocknCycles, nParticles, nDimensions);
        vec_deltaE = zeros<vec>(blocknCycles);
        vec_energySquaredSum = zeros<vec>(blocknCycles);
        vec_sqrtR12 = zeros<vec>(blocknCycles);
    }

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    QForceOld = zeros<mat>(nParticles, nDimensions);
    QForceNew = zeros<mat>(nParticles, nDimensions);


    //initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = GaussianDeviate(&idum)*sqrt(stepLength);
        }
    }

    //Set up the Slater Matrices after the move
    determinant()->updateSlaterMatrices(rOld,this);
    trialFunction()->updateSlaterDeterminant(this);


    //exit(1);

    rNew = rOld;
    //loop over Monte Carlo cycles
    int print_cycle = nCycles / 10;
    for(int cycle = 0; cycle < nCycles; cycle++) {

//        determinant()->updateSlaterMatrices(rOld,this);
        //cout << cycle << endl;
        if(my_rank == 0 && printing >= 2)
        {
            if(cycle % print_cycle == 0)
            {
                cout << (double)cycle*100./nCycles << " %" << endl;
                //print_cycle +=(double) nCycles/10;
            }
        }


        //Store the current value of the wave function
        waveFunctionOld = trialFunction()->waveFunction(rOld, this);

//        QuantumForce(rOld, QForceOld);
//        cout << "Old quantumForce "<< endl << QForceOld << endl;
        derivatives()->numericalGradient(QForceOld,rOld, this);
//         cout << "New quantumForce "<< endl << QForceOld << endl;

        QForceOld = 2.*QForceOld*h/waveFunctionOld;
        //New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + GaussianDeviate(&idum)*sqrt(stepLength)+QForceOld(i,j)*stepLength*D;
            }
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < nParticles; k++) {
                if ( k != i) {
                    for (int j=0; j < nDimensions; j++) {
                        rNew(k,j) = rOld(k,j);
                    }
                }
            }
            //Recalculate the value of the wave function
            determinant()->updateSlaterMatrices(rNew,this);
            trialFunction()->updateSlaterDeterminant(this);

            waveFunctionNew = trialFunction()->waveFunction(rNew, this);
            QuantumForce(rNew, QForceNew);
            QForceNew = QForceNew*h/waveFunctionNew; // possible typo



            //  we compute the log of the ratio of the greens functions to be used in the
            //  Metropolis-Hastings algorithm
            GreensFunction = 0.0;
            for (int j=0; j < nDimensions; j++) {
                GreensFunction += 0.5*(QForceOld(i,j)+QForceNew(i,j))*
                        (D*stepLength*0.5*(QForceOld(i,j)-QForceNew(i,j))-rNew(i,j)+rOld(i,j));
            }
            GreensFunction = exp(GreensFunction);

            // The Metropolis test is performed by moving one particle at the time
            MHRatio = GreensFunction*(waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
            if(ran2(&idum) <= MHRatio) {
                acc_moves += 1;
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    QForceOld(i,j) = QForceNew(i,j);
                    waveFunctionOld = waveFunctionNew;
//                    determinant()->updateSlaterMatrices(rNew,this); //Updating the matrices after moving the particle :)

                }

            }

            else {
                for(int j = 0; j < nDimensions; j++) {
                   rNew(i,j) = rOld(i,j);
                   QForceNew(i,j) = QForceOld(i,j);
                   determinant()->updateSlaterMatrices(rOld,this); //Updating the matrices after moving the particle :)
                   trialFunction()->updateSlaterDeterminant(this);
                }
            }


            moves += 1;
            //update energies

            deltaE = trialFunction()->localEnergy(rNew, this);
            kinSum += trialFunction()->kineticEnergy;
            potSum += trialFunction()->potentialEnergy;
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
            if(trialFunction()->getConjugate()){
                deltaLnDerivative = trialFunction()->lnDerivativeWaveFunction(rNew, this);
                lnDerivative += deltaLnDerivative;
                lnDerivativeLocalEnergy += deltaE * deltaLnDerivative;
            }

        }
        //we need to find the average value of r12
        r12 = 0;

        for(int k = 0; k < nDimensions; k++) {

            r12 += (rNew(0,k) - rNew(1,k)) * (rNew(0,k) - rNew(1,k));
        }
        averageR12 += sqrt(r12);

        if (m_blockSampling){
            // saving data at each cycle or each 10th
            // finalize by sending everything to master process
            if (reduced_blocking)
            {
                if (cycle % 10 == 0){
                    int ncycle = cycle/10;
                    vec_deltaE(ncycle) = deltaE;
                    vec_energySquaredSum(ncycle) = energySquaredSum;
                    vec_sqrtR12(ncycle) = sqrt(r12);
                    for(int i = 0; i < nParticles; i++)
                    {
                        for(int d = 0; d < nDimensions; d++)
                        mat_pos(ncycle, i, d) = rNew(i,d);
                    }
                }
            }else{
                vec_deltaE(cycle) = deltaE;
                vec_energySquaredSum(cycle) = energySquaredSum;
                vec_sqrtR12(cycle) = sqrt(r12);
                for(int i = 0; i < nParticles; i++)
                {
                    for(int d = 0; d < nDimensions; d++)
                    mat_pos(cycle, i, d) = rNew(i,d);
                }
            }


            /*
            if (my_rank == 0){
                //samplefile << setw(15) << setprecision(8) << 0;

                //samplefile << setw(15) << setprecision(8) << deltaE;
                //samplefile << setw(15) << setprecision(8) << energySquaredSum;
                //samplefile << setw(15) << setprecision(8) << sqrt(r12);
                //samplefile << setw(15) << setprecision(8) << nDimensions;

                for(int i = 0; i < nParticles; i++)
                {
                        samplefile << setw(15) << setprecision(8) << rNew(i,0);
                        samplefile << setw(15) << setprecision(8) << rNew(i,1);

                        //samplefile << setw(15) << setprecision(8) << rNew(i,2) ;

                        //sprintf(pos_buffer, "%15f\t%15f", rNew(i,0), rNew(i,1));


                }

                //sprintf(blocking_buffer, "DE %15f\t ES %15f\t SQ %15f\t D %15d\t pos %s\n",
                //        deltaE, energySquaredSum, sqrt(r12), nDimensions, pos_buffer);
                //samplefile << blocking_buffer;

                samplefile << endl;

                // receive data from processes, write to file
                for(int p = 1; p < numprocs ; p++) {

                    MPI_Recv(&rec_deltaE, 1, MPI_DOUBLE, p, 990, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&rec_energySquaredSum, 1, MPI_DOUBLE, p, 991, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&rec_sqrtR12, 1, MPI_DOUBLE, p, 992, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    //samplefile << setw(15) << setprecision(8) << p;

                    //samplefile << setw(15) << setprecision(8) << rec_deltaE;
                    //samplefile << setw(15) << setprecision(8) << rec_energySquaredSum;
                    //samplefile << setw(15) << setprecision(8) << rec_sqrtR12;
                    //samplefile << setw(15) << setprecision(8) << nDimensions;


                    for(int i = 0; i < nParticles; i++){
                        MPI_Recv(&rec_x, 1, MPI_DOUBLE, p, 993, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&rec_y, 1, MPI_DOUBLE, p, 994, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        samplefile << setw(15) << setprecision(8) << rec_x;
                        samplefile << setw(15) << setprecision(8) << rec_y;
                        //sprintf(pos_buffer, "%15f\t%15f", rec_x, rec_y);
                    }
                    //sprintf(blocking_buffer, "DE %15f\t ES %15f\t SQ %15f\t D %15d\t pos %s\n",
                    //        rec_deltaE, rec_energySquaredSum, rec_sqrtR12, nDimensions, pos_buffer);
                    //samplefile << blocking_buffer;
                    samplefile << endl;

                }

            }else
            {
                // send data to process 0
                MPI_Send(&deltaE, 1, MPI_DOUBLE, 0, 990, MPI_COMM_WORLD);
                MPI_Send(&energySquaredSum, 1, MPI_DOUBLE, 0, 991, MPI_COMM_WORLD);
                MPI_Send(&sqrtR12, 1, MPI_DOUBLE, 0, 992, MPI_COMM_WORLD);
                for(int i = 0; i < nParticles; i++){
                    MPI_Send(&rNew(i,0), 1, MPI_DOUBLE, 0, 993, MPI_COMM_WORLD);
                    MPI_Send(&rNew(i,1), 1, MPI_DOUBLE, 0, 994, MPI_COMM_WORLD);
                }

            }*/

        }
    }

    /*
    if(m_blockSampling)
    {
        samplefile << "#Alpha: " << m_alpha << " and beta: " << m_beta << endl;
        cout << "blockSampling" << endl;
    }
    */

    normalizationFactor = (double) nCycles*(double) nParticles;
    m_energy = energySum/normalizationFactor;
    m_energySquared = energySquaredSum/normalizationFactor;
    m_energyVar = sqrt((m_energySquared - m_energy*m_energy) /(double) nCycles);
    m_averageR12 = averageR12 / (double) nCycles;
    m_ratio = (double) acc_moves/(double) moves;
    m_moves = moves;
    m_kinSum = kinSum/normalizationFactor;
    m_potSum = potSum/normalizationFactor;




    if(trialFunction()->getConjugate()) {
        energyDerivative = 2*(lnDerivativeLocalEnergy-m_energy*lnDerivative)/normalizationFactor;
    }

}

void VMCSolver::QuantumForce(const mat &r, mat &QForce)
{
    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);


    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = trialFunction()->waveFunction(rMinus, this);
            waveFunctionPlus = trialFunction()->waveFunction(rPlus, this);
            QForce(i,j) =  (waveFunctionPlus-waveFunctionMinus);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
}
/*
void VMCSolver::calculateOptimalSteplength() {
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double step_min = 1;
    double ratio = 0;
    double old_ratio = 1;
    int moves = 0;
    int acc_moves = 0;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength*(ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // find optimal steplength
    for (stepLength = 0.1; stepLength <= 20.0; stepLength += 0.2){
        moves = 0;
        acc_moves = 0;
        waveFunctionOld = 0;
        waveFunctionNew = 0;

        // initial trial positions
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rOld(i,j) = stepLength*(ran2(&idum) - 0.5);
            }
        }
        rNew = rOld;
        for(int cycle = 0; cycle < 1000; cycle++) {
            waveFunctionOld = trialFunction()->waveFunction(rOld, this);
            for(int i = 0; i < nParticles; i++) {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j)+stepLength*(ran2(&idum) - 0.5);
                }
                waveFunctionNew = trialFunction()->waveFunction(rNew, this);
                if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew)/(waveFunctionOld*waveFunctionOld)) {
                    acc_moves += 1;
                    for(int j = 0; j < nDimensions; j++) {
                        rOld(i,j) = rNew(i,j);
                        waveFunctionOld = waveFunctionNew;
                    }
                } else {
                    //acc_moves -= 1;
                    for(int j = 0; j < nDimensions; j++) {
                        rNew(i,j) = rOld(i,j);
                    }
                }
                moves += 1;
            }
        }
        ratio = (double)acc_moves/(double)moves;
        if(abs(0.5-ratio) < abs(0.5-old_ratio)) {
            step_min = stepLength;
            old_ratio = ratio;
        }
    }
    stepLength = step_min;
    cout << endl << "##############################################################################" << endl << endl;
    cout << endl << "Steplength: " << stepLength << endl;
}

void VMCSolver::runMonteCarloIntegration() {
    int acc_moves = 0;
    int moves = 0;
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE = 0;
    double r12 = 0;
    double averageR12 = 0;

    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    //initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }

    rNew = rOld;
    //loop over Monte Carlo cycles
    int print_cycle = 0;

    for(int cycle = 0; cycle < nCycles; cycle++) {
        if(cycle == print_cycle)
        {
            cout << (double)cycle*100./nCycles << " %" << endl;
            print_cycle += (double) nCycles/4;
        }
        //Store the current value of the wave function
        waveFunctionOld = trialFunction()->waveFunction(rOld, this);

        //New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }
//            determinant()->updateSlaterMatrices(rNew,this); //Updating |D| and |D|⁻1, should only be done once
            //Recalculate the value of the wave function
            waveFunctionNew = trialFunction()->waveFunction(rNew, this);

            //Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                acc_moves += 1;
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                }
            }
            else {
                for(int j = 0; j < nDimensions; j++) {
                   rNew(i,j) = rOld(i,j);
                }
            }
            moves += 1;
            //update energies
//            determinant()->updateSlaterMatrices(rNew,this);

            deltaE = trialFunction()->localEnergy(rNew, this);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
        //we need to find the average value of r12
        r12 = 0;
//        for(int k = 0; k < nDimensions; k++) {
//            r12 += (rNew(0,k) - rNew(1,k)) * (rNew(0,k) - rNew(1,k));
//        }
        averageR12 += sqrt(r12);

        if (m_blockSampling &&  cycle % 10 == 0) {
            samplefile << setw(15) << setprecision(8) << deltaE;
            samplefile << setw(15) << setprecision(8) << deltaE*deltaE;
            samplefile << setw(15) << setprecision(8) << sqrt(r12) << endl;

            for(int i = 0; i < nParticles; i++ )
            {
                    samplefile << setw(15) << setprecision(8) << rNew(i,0);
                    samplefile << setw(15) << setprecision(8) << rNew(i,1);
                    samplefile << setw(15) << setprecision(8) << rNew(i,2);
            }
        }
    }
    if(m_blockSampling)
    {
        samplefile << "#Alpha: " << m_alpha << " and beta: " << m_beta << endl;
        cout << "blockSampling" << endl;
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    m_energyVar = sqrt((energySquared - energy*energy) / nCycles);  //Should we add this /(nCycles * nParticles)
    averageR12 /= (double) nCycles;

    //Storing the energy and variance calculated, used in searching for the best fit for alpha and beta
    storeEnergy(energy);
//    storeVariance(m_energyVar); //Not necessary it was already calculated

    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
    cout << "Variance: " << m_energyVar << endl;
    cout << "Moves: " << moves << endl;
    cout << "Accepted moves: " << acc_moves << endl;
    cout << "Ratio: " << (double) acc_moves/(double) moves << endl;
    cout << "Alpha: " << m_alpha << " and beta: " << m_beta << endl;
    cout << "Average distance between the electrons: " << averageR12 << endl;
    cout << "Steplength: " << stepLength << endl;



    outfile << setw(15) << setprecision(8) << energy;
    outfile << setw(15) << setprecision(8) << energySquared;
    outfile << setw(15) << setprecision(8) << m_energyVar;
    outfile << setw(15) << setprecision(8) << m_alpha;
    outfile << setw(15) << setprecision(8) << m_beta;
    outfile << setw(15) << setprecision(8) << averageR12;
    outfile << setw(15) << setprecision(8) << stepLength;
    outfile << setw(15) << nCycles << endl;
}
*/
void VMCSolver::setCycles(double cycles){
    int numprocs;
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    nCycles = cycles/numprocs;
    //cout << "Set cycles from " << cycles << " to " << nCycles << endl;
}

void VMCSolver::setNParticles(int P){
     nParticles = P;
     /*
     cout << "OK 2v" << endl;
     trialFunction()->setSpin(this);
     cout << "OK 3v" << endl;
     */
}

void VMCSolver::resetSolver(){
    totalEnergy = 0;
    totalEnergySquared = 0;
    totalEnergyVar = 0;
    totalMoves = 0;
    totalRatio = 0;
    totalAverageR12 = 0;
    totalKin = 0;
    totalPot = 0;
    m_energy = 0;
    m_energyVar = 0;
    m_energySquared = 0;
    m_moves = 0;
    m_ratio = 0;
    m_averageR12 = 0;
    m_kinSum = 0;
    m_potSum = 0;
}

int VMCSolver::checknCycles(int nCycles, int numprocs){
    if (nCycles*numprocs > 10000000)
    {
        return nCycles/10.0;
    }else{
        return nCycles;
    }
}
