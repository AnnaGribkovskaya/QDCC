#include <iostream>
#include <iomanip>      // std::setprecision
#include "generalspclass.h"
#include "qdotspbasis.h"
#include "qdotchannelset.h"
#include "ccd.h"
#include "qdotspbasis.h"
#include "channelset.h"



//using namespace std;

int main(int numberOfArguments, char **argumentList)
{


    int NumberOfShellsExact = 3;
    int NumberOfShellsStochastic = 3;
    int NumberOfElectrons = 2;
    int MC = 200;
    double homega = 1.0;

    // If a first argument is provided, it is the number of shells _STOCHASTIC_
    if(numberOfArguments > 1) NumberOfShellsStochastic = atoi(argumentList[1]);
    // If a second argument is provided, it is the number of electrons
    if(numberOfArguments > 2) NumberOfElectrons = atoi(argumentList[2]);
    // If a third argument is provided, it is the number of shells _EXACT_
    if(numberOfArguments > 3) NumberOfShellsExact = atoi(argumentList[3]);
    // If a fourth argument is provided, it is the number of MC samples
    if(numberOfArguments > 4) MC = atoi(argumentList[4]);


    generalSPclass * SPbasis;
    channelset * AllChannels;

    SPbasis = new qdotspbasis(NumberOfShellsStochastic, NumberOfShellsExact, NumberOfElectrons, homega);
    AllChannels = new qdotchannelset();
    //SPbasis->getQuantumDotStatesNumber();
    //SPbasis->getQuantumDotStates();
    //std::cout << std::setprecision(16) << SPbasis->TBME(0,1,10,11) << std::endl;

        
    ccd QD(SPbasis, AllChannels);
    std::cout << std::setprecision(16) << "Couple Cluster energy " << QD.iterateCCD(-0.01171875) << std::endl;






    return 0;

}
