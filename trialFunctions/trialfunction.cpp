#include "trialfunction.h"
#include "vmcsolver.h"

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


TrialFunction::TrialFunction()
{
//    m_nucleusDistance = 1.4;
    m_nucleusDistance = 4.63;

    m_conjugateMethod = false ;
}

double TrialFunction::spinFactor(int i, int j)
{
    if(spin(i) == spin(j))
        return 1./3; //1./4.;
    else
        return 1.0; //1./2.;
}

void TrialFunction::setNucleusDistance(double R)
{

    m_nucleusDistance = R;
    if(R < 0.001)
    {
        m_zeroDistance = true;
    }
    else m_zeroDistance = false;
}

void TrialFunction::setSpin(VMCSolver *solver)
{
    cout << spin << endl;

    spin << 0 << 1;
    cout << spin << endl;
    spin.zeros(solver->getNParticles());
    double half = (solver->getNParticles() - 1)/2.0;
    cout << solver->getNParticles() << endl;
    cout << half << endl;
    cout << spin.size() << endl;
    
    for (int i = 0; i < solver->getNParticles(); ++i)
    {
        if (i >= half){
            spin << 1;
        }else{
            spin << 0;
        }
    }
}