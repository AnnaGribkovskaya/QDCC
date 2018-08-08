#ifndef MBPT2_H
#define MBPT2_H

#include "generalspclass.h"
#include <iostream>
#include <iomanip>      // std::setprecision
#include "channelindexpair.h"
#include "channel.h"


class mbpt2
{
public:
    mbpt2(generalSPclass *);
    void test();
    //double calculateEnergyBruteForce();
    double calculateEnergyMatMult();

private:

    //vars
    generalSPclass * qsys;
    std::vector<channel>  ChannelVariety;


    //methods
    void setUpChannels();
};


#endif // MBPT2_H





