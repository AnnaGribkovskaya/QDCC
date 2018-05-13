#ifndef QDOTSPBASIS_H
#define QDOTSPBASIS_H

#include "generalspclass.h"
#include <vector>
#include <cmath>
#include <armadillo>
#include "Coulomb_Functions.hpp"


class qdotspbasis: public generalSPclass
{
public:
    //constr
    qdotspbasis(int NumberOfShellsStochastic, int NumberOfShellsExact, int ParticlesNumber, double HOStrenth);

    //vars
    int m_ShellsExact, m_ShellsStochastic, m_FermiLevel, m_StatesExact, m_StatesStochastic, m_nMax;
    double**** m_twoBodyElements;

    //methods
    double TBME(int, int, int, int);
    std::vector<double> getSPenergies();
    void getQuantumStates();
    void getQuantumStatesNumber();
    virtual std::vector<qstate> getStateVec () {return this->m_shells;}

    virtual int getShellsExact      () {return this->m_ShellsExact;}
    virtual int getShellsStochastic () {return this->m_ShellsStochastic;}
    virtual int getFermiLevel       () {return this->m_FermiLevel;}
    virtual int getStatesExact      () {return this->m_StatesExact;}
    virtual int getStatesStochastic () {return this->m_StatesStochastic;}
    virtual int getnMax             () {return this->m_nMax;}

    virtual qstate* oneState(int);
    virtual qstate* sumState(int, int);
    virtual qstate* substractState(int, int);
    virtual qstate* sumSubstractState(int, int, int);
    virtual bool isEqual(qstate*, qstate*);
    virtual void getQuantumDotStates();
    virtual void getQuantumDotStatesNumber();





private:
    // vars
    double homega;
    std::vector<double> m_HOEnergies;

    // methods
    void CalculateSPenergies();
    void fillTwoBodyElements();
    double**** create4dArray(int, int, int, int);
    void setUpStatesPolarSorted();
    void printSPenergies();

protected:
    std::vector<qstate> m_shells;

};
#endif // QDOTSPBASIS_H





