#include "neon.h"
#include "trialfunction.h"
#include "vmcsolver.h"
#include "lib.h"

#include <iostream>

using namespace std;

Neon::Neon(VMCSolver *solver)
{
    simpleFlag = false;
    m_outfileName = "Neon";

    solver->setCharge(10);
    solver->setNParticles(10);
    solver->setAlpha(10.22);
    solver->setBeta(0.109375);
    spin << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1;
}

void Neon::setSpin(VMCSolver *solver)
{
    spin << 0 << 0 << 0 << 0 << 0 << 1 << 1 << 1 << 1 << 1;
}

void Neon::updateSlaterDeterminant(VMCSolver *solver)
{
    SD = solver->determinant()->calculateDeterminant(solver);
}

double Neon::waveFunction(const mat &r, VMCSolver *solver)
{
    double rSingleParticle, alpha, beta, wf, SD, product, rij, a;
    int spin_count = 0;
    vec spins(10);
    for(int s = 0; s<10; s++)
    {
        if(s<5) spins(s)=0;
        else spins(s)=1;
    }
    product = 1.0;
    alpha = solver -> getAlpha();
    beta = solver -> getBeta();
    //vec argument(solver->getNParticles());
    if(solver->getElectronInteration())
        if(!simpleFlag)
        {
            {
            for(int i = 0; i < solver->getNParticles(); i++) {
                //argument[i] = 0.0;
                rSingleParticle = 0;

                for(int j = i + 1; j < solver->getNParticles(); j++) {
                    rij = 0;
                    for(int k = 0; k < solver->getNDimensions(); k++) {
                        rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                    }
                    rij = sqrt(rij);
        //            if(spins(i)==spins(j)) a = 1./4.;
        //            else a = 1./2.;
                    a = spinFactor(i,j);
                    //cout << i << " & " << j  << ": " << a << endl;
                    product = product * exp(a*rij/(1+beta*rij));
                    spin_count++;
                }
            }
        }

    }


    return SD*product;
}

double Neon::lnDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    double rSingleParticle, alpha, beta, sum, rij, a;
    int spin_count = 0;
    vec spins(10);
    for(int s = 0; s<10; s++)
    {
        if(s<5) spins(s)=0;
        else spins(s)=1;
    }
    sum = 0.;
    alpha = solver -> getAlpha();
    beta = solver -> getBeta();
    //vec argument(solver->getNParticles());
    if(solver->getElectronInteration())
    {
        for(int i = 0; i < solver->getNParticles(); i++) {
            //argument[i] = 0.0;
            rSingleParticle = 0;

            for(int j = i + 1; j < solver->getNParticles(); j++) {
                rij = 0;
                for(int k = 0; k < solver->getNDimensions(); k++) {
                rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                rij = sqrt(rij);
                if(spins(i)==spins(j)) a = 1./4.;
                else a = 1./2.;
                //cout << i << " & " << j  << ": " << a << endl;
                sum -= a*rij*rij/((1+beta*rij)*(1+beta*rij));
                spin_count++;
            }
        }
    }
    return sum;
}

double Neon::lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    double rSingleParticle, alpha, beta, sum, rij, a;
    int spin_count = 0;
    vec spins(10);
    for(int s = 0; s<10; s++)
    {
        if(s<5) spins(s)=0;
        else spins(s)=1;
    }
    sum = 0.;
    alpha = solver -> getAlpha();
    beta = solver -> getBeta();
    //vec argument(solver->getNParticles());
    if(solver->getElectronInteration())
    {
        for(int i = 0; i < solver->getNParticles(); i++) {
            //argument[i] = 0.0;
            rSingleParticle = 0;

            for(int j = i + 1; j < solver->getNParticles(); j++) {
                rij = 0;
                for(int k = 0; k < solver->getNDimensions(); k++) {
                rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                rij = sqrt(rij);
                if(spins(i)==spins(j)) a = 1./4.;
                else a = 1./2.;
                //cout << i << " & " << j  << ": " << a << endl;
                sum += 2*a*pow(rij,3)*beta/pow(1+beta*rij,3);
                spin_count++;
            }
        }
    }
    return sum;
}

double Neon::localEnergy(const mat &r, VMCSolver *solver)
{
    //Grabbing all the necessary constants stored in the solver
    double nParticles = solver->getNParticles();
    double nDimensions = solver->getNDimensions();
    double charge = solver->getCharge();


    // Kinetic energy

    double kineticEnergy = 0;
    if(m_analytical)
    {
    //Calculates the kinetic energy as the ratios of -1/2* ( d²/dx²|D| /|D| + 2 * (d/dx |D|/|D|)*d/dx Psi_C/Psi_C + d²/dx² Psi_C /Psi_C )

        //If we want to compute without electroninteraction with a simplified version ofthe trialfunction containing only the slater determinant
        //then only the slater determinant ratio laplacian is used.
        if(solver->getElectronInteration())
        {
            solver->derivatives()->analyticalLaplacianRatio(kineticEnergy, r, solver);
        }
        else
        {
            kineticEnergy += solver->determinant()->laplacianSlaterDeterminant(r, solver);
        }

            kineticEnergy *= -1./2.;
    } else
    {
        kineticEnergy = solver->derivatives()->numericalDoubleDerivative(r, solver) / (2.*waveFunction(r, solver));
    }


    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    if(solver->getElectronInteration())
    {
        double r12 = 0;
        for(int i = 0; i < nParticles; i++) {
            for(int j = i + 1; j < nParticles; j++) {
                r12 = 0;
                for(int k = 0; k < nDimensions; k++) {
                    r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                potentialEnergy += 1 / sqrt(r12);
            }
        }
    }


    return kineticEnergy + potentialEnergy;
}

//double Neon::psi1s(double ri, double alpha)
//{
//    return exp(-alpha*ri);
//}


//double Neon::psi2s(double ri, double alpha)
//{
//    return (1-alpha*ri/2.0)*exp(-alpha*ri/2.0);
//}

//double Neon::phi(const mat &r, double alpha, int i, int j, VMCSolver *solver)
//{
//    // returns an ansatz based on matrix row, M, in the SD
//    int nDimensions = solver->getNDimensions();
//    double ri = 0;
//    for (int k = 0; k < nDimensions; ++k) ri += r(i,k)*r(i,k);
//    ri = sqrt(ri);

//    if (j == 0)
//    {
//        return exp(-alpha*ri); // 1s
//    }
//    else if (j == 1)
//    {
//        return (1-alpha*ri/2.0)*exp(-alpha*ri/2.0); // 2s
//    }
//    else if (j>=2 && j<=4)
//    {
//        int dimension = j-2;
//        return alpha*r(i,dimension)*exp(-alpha*ri/2.0); // 2p
//    }
//}

//double Neon::SlaterDeterminant(const mat &r,double alpha, VMCSolver *solver)
//{
//    int i, j, Nhalf, *indx;
//    double d1, d2, SD;
//    int nParticles= solver->getNParticles();
//    Nhalf = nParticles/2;
//    indx = new int [Nhalf];
//    detUp = zeros<mat>(Nhalf, Nhalf);
//    detDown = zeros<mat>(Nhalf, Nhalf);
//    // fill matrix detUp and detDown

//    for (int k = 0; k <  Nhalf; ++k)
//    {
//        for (i = 0; i < Nhalf; ++i)
//        {
//            // for detUp

//            detUp(i,k) =  phi(r, alpha, i, k, solver);

//            // for detDown
//            detDown(i,k) =  phi(r, alpha, i+Nhalf, k, solver);
//        }
//    }
//    // decompose A (phi matrix) to B & C
//    /*
//     * End up with
//     *     (c00 c01 c02 c03)
//     * A = (b10 c11 c12 c13)
//     *     (b20 b21 c22 c23)
//     *     (b30 b31 b32 c33)
//     */

//    ludcmp(detUp, Nhalf, indx, &d1);
//    ludcmp(detDown, Nhalf, indx, &d2);

//    // compute SD as c00*c11*..*cnn
//    SD = 1;
//    for (i = 0; i < Nhalf; ++i)
//    {
//        SD *= detUp(i, i)*detDown(i, i);
//    }
//    // return SD
//    return d1*d2*SD;
//}
