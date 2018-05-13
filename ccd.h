#ifndef CCD_H
#define CCD_H


#include "generalspclass.h"
#include "channel.h"
#include <iomanip>      // std::setprecision
#include "channelset.h"


class ccd
{
public:
    ccd(generalSPclass *, channelset *);
    double iterateCCD(double);    

private:

    //const
    const double m_toleranceCCD = 1e-6;

    //vars
    generalSPclass * qsys;
    channelset * channels;

    std::vector<channel>  ChannelVariety;
    std::vector<channel>  ChannelVariety1;
    std::vector<channel>  ChannelVariety2;

    std::vector<symblock> m_hhhhVBlock;
    std::vector<symblock> m_ppppVBlock;
    std::vector<symblock> m_hMphMpVBlock;
    std::vector<symblock> m_hhppVBlock;
    std::vector<symblock> m_hMppMhVBlock;
    std::vector<symblock> m_ppMhVBlock;
    std::vector<symblock> m_hhMpVBlock;
    std::vector<symblock> m_pphhVBlock;

    std::vector<symblock> m_pphhTBlock;         // particle + particle = hole + hole sy block
    std::vector<symblock> m_pphhTBlockPrev;
    std::vector<symblock> m_pMhhMpTBlock;       //particle - hole =hole - particle sym block
    std::vector<symblock> m_pMhhMpTBlockPrev;
    std::vector<symblock> m_ppMhTBlock;
    std::vector<symblock> m_ppMhTBlockPrev;
    std::vector<symblock> m_hhMpTBlock;
    std::vector<symblock> m_hhMpTBlockPrev;


    //methods

    void setUpInterractionMatrixBlocksQ3();
    void setUpInitialAmplitudesQ3();
    void setUpInterractionMatrixBlocksQ4();
    void setUpInitialAmplitudesQ4();
    void setUpInterractionMatrixBlocks();
    void setUpInitialAmplitudes();

    void L3();
    void Q2();
    void Q3();
    void Q4();

    void calculateAmplitudes();
    void initializeVandTandChannels();
    void updateAmplitudes();
    double updateAmplitudesL3(int, int, int, int);
    double updateAmplitudesQ3(int, int, int, int);
    double updateAmplitudesQ4(int, int, int, int);
    double computeCCDCorrEnergy();

    void L1plusL2plusQ1();
    void L3Permutations();
    void Q2Permutations();
    void Q3Permutations();
    void Q4Permutations();

    double L3recoupled(int, int, int, int);
    double Q2recoupled(int, int, int, int);
    double Q3recoupled(int, int, int, int);
    double Q4recoupled(int, int, int, int);

    void printAmplitudes();
    void printAmplitudesPrev();


};
#endif // CCD_H
