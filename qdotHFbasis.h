#ifndef QDOTHFBASIS_H
#define QDOTHFBASIS_H

#include "generalspclass.h"
#include <vector>
#include <cmath>
#include <armadillo>
#include "Coulomb_Functions.hpp"

class qdotHFbasis: public generalSPclass

{
public:
    //constr
    qdotHFbasis(int NumberOfShellsStochastic, int NumberOfShellsExact, int ParticlesNumber, double HOStrenth);
    //vars
    int m_ShellsExact, m_ShellsStochastic, m_FermiLevel, m_StatesExact, m_StatesStochastic, m_nMax;
    double**** m_twoBodyElements;
    double**** m_twoBodyElementsInterm1;
    double**** m_twoBodyElementsInterm2;
    double**** m_twoBodyElementsInterm3;


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

    virtual qstate oneState(int);
    virtual qstate sumState(int, int);
    virtual qstate substractState(int, int);
    virtual qstate sumSubstractState(int, int, int);
    virtual bool isEqual(qstate, qstate);
    virtual void getQuantumDotStates();
    virtual void getQuantumDotStatesNumber();


    void test();

private:
    // vars
    double homega;
    std::vector<double> m_HOEnergies;

    arma::vec eigval_previous;
    arma::vec eigval;
    arma::mat eigvec;
    //arma::mat m_HOEnergies;
    arma::mat m_C;
    arma::mat m_HF;

    // methods
    void CalculateSPenergies();
    void fillTwoBodyElements();
    double**** create4dArray(int, int, int, int);
    void setUpStatesPolarSorted();
    void printSPenergies();

    //For HF
    void setCoefficientMatrix(arma::mat);
    arma::mat computeDensityMatrix();
    void CalculateNonIntEnergy();
    void computeHFmatrix(arma::mat);
    double computeHartreeFockEnergyDifference();
    void applyHartreeFockMethod();
    void computeHartreeFockEnergy(arma::mat);
    void convertTBMEtoHFbasis();
    void convertOneBodyElementstoHFbasis();
    double**** delete4dArray(int, int, int);



protected:
    std::vector<qstate> m_shells;

};

#endif // QDOTHFBASIS_H
