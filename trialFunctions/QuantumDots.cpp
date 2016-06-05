#include "QuantumDots.h"

QuantumDots::QuantumDots(VMCSolver *solver)
{
    simpleFlag = false;
    m_analytical = true;

    m_outfileName = "QuantumDots";

    solver->setCharge(6);
    solver->setNParticles(2);
    solver->setAlpha(1.0);
    solver->setBeta(1.0000001);
    solver->setOmega(1.0);
    
    spin.zeros(56);
    for (int i = 0; i < 56; ++i)
    {
        spin(i) = 0.5*(1 - pow(-1, i));
    }
}
// NB: SPIN FOR 3D!

void QuantumDots::updateSlaterDeterminant(VMCSolver *solver)
{
    SD = solver->determinant()->calculateDeterminant(solver);
}

double QuantumDots::waveFunction(const mat &r, VMCSolver *solver)
{
    double rSingleParticle, rij, a;
    double product = 1.0;
    double alpha = solver -> getAlpha();
    double beta = solver -> getBeta();

    ///*
    //Calculate the Jastrow factor
    if(solver->getJastrow())
    {
        for(int i = 0; i < solver->getNParticles(); i++) {
            rSingleParticle = 0;
            for(int j = i + 1; j < solver->getNParticles(); j++) {
                rij = 0;
                for(int k = 0; k < solver->getNDimensions(); k++) {
                    rij += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                rij = sqrt(rij);
                a = spinFactor(i,j);
                product = product * exp(a*rij/(1+beta*rij));
            }
        }
    }
    //*/

//cout << product << endl;

//    cout << "Before SD" << endl;
    //SD = solver->determinant()->calculateDeterminant(r,alpha,solver); //SlaterDeterminant(r, alpha, solver);

    return SD*product;
}

double QuantumDots::localEnergy(const mat &r, VMCSolver *solver)
{

    kineticEnergy = 0;
    potentialEnergy = 0;
    double r2sum = 0;
    int Nparticles = solver->getNParticles();

    double r12 = norm(r.row(0) - r.row(1));


    if(m_analytical)
    {
        // Calculates the kinetic energy as the ratios of
        // -1/2* ( d²/dx²|D| /|D| + 2 * (d/dx |D|/|D|)*d/dx Psi_C/Psi_C + d²/dx² Psi_C /Psi_C )
        // If we want to compute without electroninteraction with a simplified version of the
        // trialfunction containing only the slater determinant, then only the slater determinant
        // ratio laplacian is used.
        if(solver->getJastrow())
        {
            solver->derivatives()->analyticalLaplacianRatio(kineticEnergy, r, solver);
        }
        else
        {
            kineticEnergy += solver->determinant()->laplacianSlaterDeterminant(r, solver);
        }
        kineticEnergy *= -1./2.;
    }
    else
        kineticEnergy =  solver->derivatives()->numericalDoubleDerivative(r, solver) / (2.*waveFunction(r, solver)); //The minus looks to be baked into the
                                                                                                                     //doubleDerivative Function


    for (int i = 0; i < Nparticles; ++i){
        r2sum += norm(r.row(i)) * norm(r.row(i));
    }
    potentialEnergy = 0.5 * r2sum;
    potentialEnergy = potentialEnergy*solver->getOmega()*solver->getOmega();
    
    if(solver->getElectronInteration())
    {
        // Contribution from electron-electron potential
        double r12 = 0;
        double yuk_my = 1;
        for(int i = 0; i < solver->getNParticles(); i++) {
            for(int j = i + 1; j < solver->getNParticles(); j++) {
                r12 = 0;
                for(int k = 0; k < solver->getNDimensions(); k++) {
                    r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                }
                r12 = sqrt(r12);
                potentialEnergy += 1 / r12;
                //potentialEnergy += exp(-yuk_my*r12) / r12; // for yukawa interaction
            }
        }
    }

    //cout << SD << endl;
    //kineticEnergy = 4*kineticEnergy*solver->getOmega()*solver->getOmega();
    return kineticEnergy + potentialEnergy;
}

double QuantumDots::lnDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    /*double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));
    double dotproduct = dot(r.row(0),r.row(1));
    double derivative = 0;
    double beta = solver->getBeta();

    derivative = (-1*solver->getAlpha()*(r1+r2)*(1-dotproduct/(r1*r2))*pow(1+r12*beta,2)+3*r12*beta+r12+3)/pow(1+r12*beta,5);

    return derivative;
    */
    cout << "Using unimplemented derivative" << endl;
    return 0;
}

double QuantumDots::lnSecondDerivativeWaveFunction(const mat &r, VMCSolver *solver)
{
    /*
    double r1 = norm(r.row(0));
    double r2 = norm(r.row(1));
    double r12 = norm(r.row(0) - r.row(1));
    double dotproduct = dot(r.row(0),r.row(1));
    double factor = solver->getAlpha()*(r1+r2)*(1-dotproduct/(r1*r2));
    double secondDerivative = 0;
    double beta = solver->getBeta();

    secondDerivative = (r12*(3*r12*r12*factor*beta*beta+6*r12*factor*beta-12*r12*beta-5*r12+3*factor-12))/pow(1+r12*beta,6);

    return secondDerivative;
    */
    cout << "Using unimplemented derivative" << endl;
    return 0;
}

