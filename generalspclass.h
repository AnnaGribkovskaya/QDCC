#ifndef GENERALSPCLASS_H
#define GENERALSPCLASS_H

#include "qstate.h"
#include <vector>
#include <armadillo>

class generalSPclass
{
public:
    generalSPclass();
    virtual void getQuantumDotStates() = 0;
    virtual void getQuantumDotStatesNumber() = 0;
    virtual double TBME(int, int, int, int) = 0;
    virtual std::vector<double> getSPenergies() = 0;


    virtual std::vector<qstate> getStateVec () = 0;
    virtual int getShellsExact      () = 0;
    virtual int getShellsStochastic () = 0;
    virtual int getFermiLevel       () = 0;
    virtual int getStatesExact      () = 0;
    virtual int getStatesStochastic () = 0;
    virtual int getnMax () = 0;
    virtual qstate* oneState(int) = 0;
    virtual qstate* sumState(int, int) = 0;
    virtual qstate* substractState(int, int) = 0;
    virtual qstate* sumSubstractState(int, int, int) = 0;
    virtual bool isEqual(qstate*, qstate*) = 0;



};

#endif // GENERALSPCLASS_H
