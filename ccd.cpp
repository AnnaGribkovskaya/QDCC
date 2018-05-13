#include "ccd.h"

ccd::ccd(generalSPclass * qsystem, channelset * allchannels){
    qsys = qsystem;
    channels = allchannels;


}

void ccd::setUpInterractionMatrixBlocks(){
    for(channel onechannel : ChannelVariety){
            m_ppppVBlock.emplace_back(symblock(onechannel.m_ParticleParticleVec.size(), onechannel.m_ParticleParticleVec.size()));
            for(unsigned int ab = 0; ab < onechannel.m_ParticleParticleVec.size(); ab++){
                channelindexpair AB = onechannel.m_ParticleParticleVec.at(ab);
                for(unsigned int cd = 0; cd < onechannel.m_ParticleParticleVec.size(); cd++){
                    channelindexpair CD = onechannel.m_ParticleParticleVec.at(cd);
                    m_ppppVBlock.back().setElement((int)ab, (int)cd, qsys->TBME(AB.first(), AB.second(), CD.first(), CD.second())) ;
                }
            }

            m_hhhhVBlock.emplace_back(symblock(onechannel.m_HoleHoleVec.size(), onechannel.m_HoleHoleVec.size()));
            for(unsigned int kl = 0; kl < onechannel.m_HoleHoleVec.size(); kl++){
                channelindexpair KL = onechannel.m_HoleHoleVec.at(kl);
                for(unsigned int ij = 0; ij < onechannel.m_HoleHoleVec.size(); ij++){
                    channelindexpair IJ = onechannel.m_HoleHoleVec.at(ij);
                    m_hhhhVBlock.back().setElement((int)kl,(int)ij, qsys->TBME(KL.first(), KL.second(), IJ.first(), IJ.second()));
                }
            }

            m_hMphMpVBlock.emplace_back(symblock(onechannel.m_HoleMinusParticleVec.size(), onechannel.m_HoleMinusParticleVec.size()));
            for(unsigned int kc = 0; kc < onechannel.m_HoleMinusParticleVec.size(); kc++){
                channelindexpair KC = onechannel.m_HoleMinusParticleVec.at(kc);
                for(unsigned int jb = 0; jb < onechannel.m_HoleMinusParticleVec.size(); jb++){
                    channelindexpair JB = onechannel.m_HoleMinusParticleVec.at(jb);
                    m_hMphMpVBlock.back().setElement((int)kc,(int)jb, qsys->TBME(KC.first(), JB.second(), KC.second(), JB.first() ));
                }
            }

            m_hhppVBlock.emplace_back(symblock(onechannel.m_HoleHoleVec.size(), onechannel.m_ParticleParticleVec.size()));
            for(unsigned int kl = 0; kl < onechannel.m_HoleHoleVec.size(); kl++){
                channelindexpair KL = onechannel.m_HoleHoleVec.at(kl);
                for(unsigned int cd = 0; cd < onechannel.m_ParticleParticleVec.size(); cd++){
                    channelindexpair CD = onechannel.m_ParticleParticleVec.at(cd);
                    m_hhppVBlock.back().setElement((int)kl,(int)cd, qsys->TBME(KL.first(), KL.second(), CD.first(), CD.second())) ;
                }
            }

            m_hMppMhVBlock.emplace_back(symblock(onechannel.m_HoleMinusParticleVec.size(), onechannel.m_ParticleMinusHoleVec.size()));
            for(unsigned int kc = 0; kc < onechannel.m_HoleMinusParticleVec.size(); kc++){
                channelindexpair KC = onechannel.m_HoleMinusParticleVec.at(kc);
                for(unsigned int dl = 0; dl < onechannel.m_ParticleMinusHoleVec.size(); dl++){
                    channelindexpair DL = onechannel.m_ParticleMinusHoleVec.at(dl);
                    m_hMppMhVBlock.back().setElement((int)kc,(int)dl, qsys->TBME( KC.first(), DL.second(), KC.second(), DL.first()));
                }
            }

            m_pphhVBlock.emplace_back(symblock(onechannel.m_ParticleParticleVec.size(), onechannel.m_HoleHoleVec.size()));
            for(unsigned int kl = 0; kl < onechannel.m_HoleHoleVec.size(); kl++){
                channelindexpair KL = onechannel.m_HoleHoleVec.at(kl);
                for(unsigned int cd = 0; cd < onechannel.m_ParticleParticleVec.size(); cd++){
                    channelindexpair CD = onechannel.m_ParticleParticleVec.at(cd);
                    m_pphhVBlock.back().setElement((int)cd,(int)kl, qsys->TBME(CD.first(), CD.second(), KL.first(), KL.second())) ;
                }
            }
    }
}

void ccd::setUpInterractionMatrixBlocksQ3(){
    for(channel onechannel : ChannelVariety1){
            m_ppMhVBlock.emplace_back(symblock(1, onechannel.m_ParticlePlusParticleMinusHoleVec.size()));
            for(unsigned int abi = 0; abi < onechannel.m_ParticlePlusParticleMinusHoleVec.size(); abi++){
                channelindexpair ABI = onechannel.m_ParticlePlusParticleMinusHoleVec.at(abi);
                for(unsigned int j = 0; j < onechannel.m_HoleVec.size(); j++){
                    channelindexpair J = onechannel.m_HoleVec.at(j);
                    m_ppMhVBlock.back().setElement(0, (int)abi, qsys->TBME( ABI.third(), J.first(), ABI.first(), ABI.second()));
                }
            }
    }
}

void ccd::setUpInterractionMatrixBlocksQ4(){
    for(channel onechannel : ChannelVariety2){
            m_hhMpVBlock.emplace_back(symblock(onechannel.m_HolePlusHoleMinusParticleVec.size(), 1));
            for(unsigned int ija = 0; ija < onechannel.m_HolePlusHoleMinusParticleVec.size(); ija++){
                channelindexpair IJA = onechannel.m_HolePlusHoleMinusParticleVec.at(ija);
                for(unsigned int b = 0; b < onechannel.m_ParticleVec.size(); b++){
                    channelindexpair B = onechannel.m_ParticleVec.at(b);
                    m_hhMpVBlock.back().setElement((int)ija, 0 , qsys->TBME( IJA.first(), IJA.second(), IJA.third(), B.first()));
                }
            }
    }
}

void ccd::setUpInitialAmplitudes(){
    //arma::vec Epsilon;
    //Epsilon.zeros(qsys->getStateVec().size());
     //Epsilon = qsys->getSPenergies();
    for(channel onechannel : ChannelVariety){
        m_pphhTBlock.emplace_back(symblock(onechannel.m_ParticleParticleVec.size(),onechannel.m_HoleHoleVec.size()));
        m_pphhTBlock.back().setZeros();
        m_pphhTBlockPrev.emplace_back(symblock(onechannel.m_ParticleParticleVec.size(),onechannel.m_HoleHoleVec.size()));
        for(unsigned int ab = 0; ab < onechannel.m_ParticleParticleVec.size(); ab++){
            channelindexpair AB = onechannel.m_ParticleParticleVec.at(ab);
            for(unsigned int ij = 0; ij < onechannel.m_HoleHoleVec.size(); ij++){
                channelindexpair IJ = onechannel.m_HoleHoleVec.at(ij);
                m_pphhTBlockPrev.back().setElement((int)ab, (int)ij, (qsys->TBME(AB.first(), AB.second(), IJ.first(), IJ.second()))
                                                                        /(qsys->getSPenergies().at(IJ.first()) + qsys->getSPenergies().at(IJ.second()) - qsys->getSPenergies().at(AB.first()) - qsys->getSPenergies().at(AB.second())));
                m_pphhTBlock.back().setRowMap((int)ab, AB.first(), AB.second());
                m_pphhTBlock.back().setColMap((int)ij, IJ.first(), IJ.second());
            }
        }

        m_pMhhMpTBlock.emplace_back(symblock(onechannel.m_ParticleMinusHoleVec.size(),onechannel.m_HoleMinusParticleVec.size()));
        m_pMhhMpTBlock.back().setZeros();
        m_pMhhMpTBlockPrev.emplace_back(symblock(onechannel.m_ParticleMinusHoleVec.size(),onechannel.m_HoleMinusParticleVec.size()));
        for(unsigned int ai = 0; ai < onechannel.m_ParticleMinusHoleVec.size(); ai++){
            channelindexpair AI = onechannel.m_ParticleMinusHoleVec.at(ai);
            for(unsigned int kc = 0; kc < onechannel.m_HoleMinusParticleVec.size(); kc++){
                channelindexpair KC = onechannel.m_HoleMinusParticleVec.at(kc);
                m_pMhhMpTBlockPrev.back().setElement((int)ai, (int)kc, qsys->TBME(AI.first(), KC.second(), AI.second(), KC.first())
                        /(qsys->getSPenergies().at(KC.first()) + qsys->getSPenergies().at(AI.second()) - qsys->getSPenergies().at(AI.first()) - qsys->getSPenergies().at(KC.second())));
                m_pMhhMpTBlock.back().setRowMap((int)ai, AI.first(), AI.second());
                m_pMhhMpTBlock.back().setColMap((int)kc, KC.first(), KC.second());
            }
        }
    }
}

void ccd::setUpInitialAmplitudesQ3(){
    //arma::vec Epsilon;
    //Epsilon.zeros(qsys->getStateVec().size());
    //Epsilon = qsys->getSPenergies();
    for(channel onechannel : ChannelVariety1){
        m_ppMhTBlock.emplace_back(symblock(onechannel.m_ParticlePlusParticleMinusHoleVec.size(), 1));
        m_ppMhTBlock.back().setZeros();
        m_ppMhTBlockPrev.emplace_back(symblock(onechannel.m_ParticlePlusParticleMinusHoleVec.size(), 1));
        for(unsigned int abj = 0; abj < onechannel.m_ParticlePlusParticleMinusHoleVec.size(); abj++){
            channelindexpair ABJ = onechannel.m_ParticlePlusParticleMinusHoleVec.at(abj);
            m_ppMhTBlock.back().setRowMap3x1((int)abj, ABJ.first(), ABJ.second(), ABJ.third());
            for(unsigned int l = 0; l < onechannel.m_HoleVec.size(); l++){
                channelindexpair L = onechannel.m_HoleVec.at(l);
                m_ppMhTBlockPrev.back().setElement((int)abj, 0, qsys->TBME(ABJ.first(), ABJ.second(), L.first(), ABJ.third())
                        /(qsys->getSPenergies().at(L.first()) + qsys->getSPenergies().at(ABJ.third()) - qsys->getSPenergies().at(ABJ.first()) - qsys->getSPenergies().at(ABJ.second())));
                m_ppMhTBlock.back().setColMap1x3((int)abj, L.first());     // as write matrix of size[some rows x 1] unique index for each element is it's row index
            }
        }
    }
}

void ccd::setUpInitialAmplitudesQ4(){
    //arma::vec Epsilon;
    //Epsilon.zeros(qsys->getStateVec().size());
    //Epsilon = qsys->getSPenergies();
    for(channel onechannel : ChannelVariety2){
        m_hhMpTBlock.emplace_back(symblock(1, onechannel.m_HolePlusHoleMinusParticleVec.size()));
        m_hhMpTBlock.back().setZeros();
        m_hhMpTBlockPrev.emplace_back(symblock(1, onechannel.m_HolePlusHoleMinusParticleVec.size()));
        for(unsigned int klc = 0; klc < onechannel.m_HolePlusHoleMinusParticleVec.size(); klc++){
            channelindexpair KLC = onechannel.m_HolePlusHoleMinusParticleVec.at(klc);
            m_hhMpTBlock.back().setColMap3x1((int)klc, KLC.first(), KLC.second(), KLC.third());
            for(unsigned int a = 0; a < onechannel.m_ParticleVec.size(); a++){
                channelindexpair A = onechannel.m_ParticleVec.at(a);
                m_hhMpTBlockPrev.back().setElement(0, (int)klc, qsys->TBME(A.first(), KLC.third(),  KLC.first(),  KLC.second())
                        /(qsys->getSPenergies().at(KLC.first()) + qsys->getSPenergies().at(KLC.second()) - qsys->getSPenergies().at(A.first()) - qsys->getSPenergies().at(KLC.third())));
                m_hhMpTBlock.back().setRowMap1x3((int)klc, A.first()); // as write matrix of size[some rows x 1] unique index for each element is it's row index

            }
        }
    }
}

void ccd::calculateAmplitudes(){

    L1plusL2plusQ1();
    L3();
    L3Permutations();
    Q2();
    Q2Permutations();
    Q3();
    Q3Permutations();
    Q4();
    Q4Permutations();

    //arma::vec Epsilon;
    //Epsilon.zeros(qsys->getStateVec().size());
    //Epsilon = qsys->getSPenergies();

    for(unsigned int ch = 0; ch < ChannelVariety.size(); ch++){
        channel onechannel = ChannelVariety.at(ch);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            for(unsigned int ab = 0; ab < ChannelVariety.at(ch).m_ParticleParticleVec.size(); ab++){
                channelindexpair AB = onechannel.m_ParticleParticleVec.at(ab);
                for(unsigned int ij = 0; ij < ChannelVariety.at(ch).m_HoleHoleVec.size(); ij++){
                    channelindexpair IJ = onechannel.m_HoleHoleVec.at(ij);
                    double value = 0.0;
                    value = m_pphhVBlock.at(ch).getElement(ab, ij);
                    m_pphhTBlock.at(ch).setElement(ab, ij, (value + m_pphhTBlock.at(ch).getElement(ab,ij))/((qsys->getSPenergies().at(IJ.first()) + qsys->getSPenergies().at(IJ.second()) - qsys->getSPenergies().at(AB.first()) - qsys->getSPenergies().at(AB.second()))) );
                }
            }
        }
    }
}

void ccd::updateAmplitudes(){
    //update m_pphhTblockPrev
    for(unsigned int chan = 0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            for(unsigned int ab = 0; ab < onechannel.m_ParticleParticleVec.size(); ab++){
                for(unsigned int ij = 0; ij < onechannel.m_HoleHoleVec.size(); ij++){
                    m_pphhTBlockPrev.at(chan).setElement(ab,ij, m_pphhTBlock.at(chan).getElement(ab,ij));
                }
            }
        }
    }
    // Updating ampitudes for L3 and Q2 --- m_pMhhMpTBlock update
    for(unsigned int chan1 = 0; chan1 < ChannelVariety.size(); chan1++){
        channel onechannel1 = ChannelVariety.at(chan1);
        symblock amplitudes = m_pMhhMpTBlock.at(chan1);
        for(unsigned int ai = 0; ai < onechannel1.m_ParticleMinusHoleVec.size(); ai++){
            for(unsigned int jb = 0; jb < onechannel1.m_HoleMinusParticleVec.size(); jb++){
                int a = amplitudes.getRowMap(ai).first();
                int i = amplitudes.getRowMap(ai).second();
                int j1 = amplitudes.getColMap(jb).first();
                int b1 = amplitudes.getColMap(jb).second();
                double value = updateAmplitudesL3(a,b1,i,j1);
                m_pMhhMpTBlockPrev.at(chan1).setElement((int)ai, (int)jb, value);
            }
        }
    }
    // updating ampitudes for Q3
    for(unsigned int chanQ3 = 0; chanQ3 < ChannelVariety1.size(); chanQ3++){
        channel onechannelQ3 = ChannelVariety1.at(chanQ3);
        symblock amplitudesQ3 = m_ppMhTBlock.at(chanQ3);
        for(unsigned int abj = 0; abj < onechannelQ3.m_ParticlePlusParticleMinusHoleVec.size(); abj++){
                int a1 = amplitudesQ3.getRowMap3x1(abj).first();
                int b = amplitudesQ3.getRowMap3x1(abj).second();
                int j = amplitudesQ3.getRowMap3x1(abj).third();
                int i1 = amplitudesQ3.getColMap1x3(abj).first();
                double valueQ3 = updateAmplitudesQ3(a1,b,i1,j);
                m_ppMhTBlockPrev.at(chanQ3).setElement((int)abj , 0, valueQ3);
        }
    }
    // Updating ampitudes for Q4
    for(unsigned int chanQ4 = 0; chanQ4 < ChannelVariety2.size(); chanQ4++){
        channel onechannelQ4 = ChannelVariety2.at(chanQ4);
        symblock amplitudesQ4 = m_hhMpTBlock.at(chanQ4);
        for(unsigned int klc = 0; klc < onechannelQ4.m_HolePlusHoleMinusParticleVec.size(); klc++){
                int k2 = amplitudesQ4.getColMap3x1(klc).first();
                int l = amplitudesQ4.getColMap3x1(klc).second();
                int c2 = amplitudesQ4.getColMap3x1(klc).third();
                int a2 = amplitudesQ4.getRowMap1x3(klc).first();
                double valueQ4 = updateAmplitudesQ4(a2,c2,k2,l);
                m_hhMpTBlockPrev.at(chanQ4).setElement(0, (int)klc, valueQ4);
        }
    }
}

double ccd::updateAmplitudesQ4(int p, int q, int r, int s){

    double value = 0.0;
    for(unsigned int chan=0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            channelindexpair PQ = onechannel.m_ParticleParticleVec.at(0);
            if (qsys->isEqual(qsys->sumState(p,q), qsys->sumState(r,s)) && qsys->isEqual(qsys->sumState(p,q), qsys->sumState(PQ.first(),PQ.second()))){
                value =  m_pphhTBlock.at(chan).getElement(m_pphhTBlock.at(chan).getRowMapInverse(p,q),m_pphhTBlock.at(chan).getColMapInverse(r,s) );
            }
        }
    }
    return value;
}

double ccd::updateAmplitudesQ3(int p, int q, int r, int s){

    double value = 0.0;
    for(unsigned int chan=0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            channelindexpair PQ = onechannel.m_ParticleParticleVec.at(0);
            if ( qsys->isEqual(qsys->sumState(p,q), qsys->sumState(r,s)) && qsys->isEqual(qsys->sumState(p,q), qsys->sumState(PQ.first(),PQ.second()))  ){
                value =  m_pphhTBlock.at(chan).getElement(m_pphhTBlock.at(chan).getRowMapInverse(p,q),m_pphhTBlock.at(chan).getColMapInverse(r,s) );
            }
        }
    }
    return value;
}

double ccd::updateAmplitudesL3(int p, int q, int r, int s){

    double value = 0.0;
    for(unsigned int chan = 0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            channelindexpair PQ = onechannel.m_ParticleParticleVec.at(0);
            if (   qsys->isEqual(qsys->sumState(p,q), qsys->sumState(r,s)) && qsys->isEqual(qsys->sumState(p,q), qsys->sumState(PQ.first(),PQ.second())) ){
                value = m_pphhTBlock.at(chan).getElement(m_pphhTBlock.at(chan).getRowMapInverse(p,q), m_pphhTBlock.at(chan).getColMapInverse(r,s));
            }
        }
    }
    return value;
}

double ccd::iterateCCD(double MBPT2Energy){
    initializeVandTandChannels();
    double EnergyDiff = 100.0;
    double CorrelationEnergy;
    double CorrelationEnergyOld = MBPT2Energy;
    while (EnergyDiff > m_toleranceCCD){
        calculateAmplitudes();
        CorrelationEnergy = computeCCDCorrEnergy();
        std::cout << std::setprecision(16) << "Loop Cor Energy " << CorrelationEnergy << std::endl;
        EnergyDiff = CorrelationEnergy - CorrelationEnergyOld;
        EnergyDiff = std::abs(EnergyDiff);
        CorrelationEnergyOld = CorrelationEnergy;
        updateAmplitudes();
    }
    return CorrelationEnergy;
}

double ccd::computeCCDCorrEnergy(){
    double CorrelationEnergy = 0.0;
    for(unsigned int chan = 0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
           CorrelationEnergy += trace(m_hhppVBlock.at(chan).getMatrix()*m_pphhTBlock.at(chan).getMatrix());
           //std::cout << "Energy per channel " << CorrelationEnergy << std::endl;
        }
    }
    return 0.25*CorrelationEnergy;

}

void ccd::initializeVandTandChannels(){

    channels->setUpChannels(qsys);
    ChannelVariety = channels->getChan();
    setUpInterractionMatrixBlocks();
    setUpInitialAmplitudes();

    channels->setUpChannelsQ3(qsys);
    ChannelVariety1 = channels->getChanQ3();
    setUpInterractionMatrixBlocksQ3();
    setUpInitialAmplitudesQ3();

    channels->setUpChannelsQ4(qsys);
    ChannelVariety2 = channels->getChanQ4();
    setUpInterractionMatrixBlocksQ4();
    setUpInitialAmplitudesQ4();
}

void ccd::L1plusL2plusQ1(){
    for(unsigned int chan=0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            m_pphhTBlock.at(chan).setMatrix(0.5*m_ppppVBlock.at(chan).getMatrix()*m_pphhTBlockPrev.at(chan).getMatrix() +
                                            0.5*m_pphhTBlockPrev.at(chan).getMatrix()*m_hhhhVBlock.at(chan).getMatrix() +
                                            0.25*m_pphhTBlockPrev.at(chan).getMatrix()*m_hhppVBlock.at(chan).getMatrix()*m_pphhTBlockPrev.at(chan).getMatrix());
        }
    }
}

void ccd::L3(){
    for(unsigned int chan=0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleMinusHoleVec.size()*onechannel.m_HoleMinusParticleVec.size() != 0){
            m_pMhhMpTBlock.at(chan).setMatrix(m_pMhhMpTBlockPrev.at(chan).getMatrix()*m_hMphMpVBlock.at(chan).getMatrix());
        }
    }
}

void ccd::Q2(){
    for(unsigned int chan=0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleMinusHoleVec.size()*onechannel.m_HoleMinusParticleVec.size() != 0){
            m_pMhhMpTBlock.at(chan).setMatrix(m_pMhhMpTBlockPrev.at(chan).getMatrix()*m_hMppMhVBlock.at(chan).getMatrix()*m_pMhhMpTBlockPrev.at(chan).getMatrix());
        }
    }
}

void ccd::Q3(){
    for(unsigned int chan=0; chan < ChannelVariety1.size(); chan++){
        channel onechannel = ChannelVariety1.at(chan);
        if (onechannel.m_ParticlePlusParticleMinusHoleVec.size()*onechannel.m_HoleVec.size() != 0){
            m_ppMhTBlock.at(chan).setMatrix(-0.5*m_ppMhTBlockPrev.at(chan).getMatrix()*m_ppMhVBlock.at(chan).getMatrix()*m_ppMhTBlockPrev.at(chan).getMatrix());
        }
    }
}

void ccd::Q4(){
    for(unsigned int chan=0; chan < ChannelVariety2.size(); chan++){
        channel onechannel = ChannelVariety2.at(chan);
        if (onechannel.m_HolePlusHoleMinusParticleVec.size()*onechannel.m_ParticleVec.size() != 0){
            m_hhMpTBlock.at(chan).setMatrix(-0.5*m_hhMpTBlockPrev.at(chan).getMatrix()*m_hhMpVBlock.at(chan).getMatrix()*m_hhMpTBlockPrev.at(chan).getMatrix());
        }
    }
}

double ccd::Q3recoupled(int a, int b, int i, int j){

    double value = 0.0;
    for(unsigned chan = 0; chan < ChannelVariety1.size(); chan++){
        channel onechannel = ChannelVariety1.at(chan);
        if (onechannel.m_ParticlePlusParticleMinusHoleVec.size()*onechannel.m_HoleVec.size() != 0){
            channelindexpair PQS = onechannel.m_ParticlePlusParticleMinusHoleVec.at(0);
            if (  qsys->isEqual(qsys->sumSubstractState(a,b,j),qsys->oneState(i)) && qsys->isEqual(qsys->sumSubstractState(a,b,j), qsys->sumSubstractState(PQS.first(),PQS.second(), PQS.third())) ){
                value =  m_ppMhTBlock.at(chan).getElement(m_ppMhTBlock.at(chan).getRowMapInverse3x1(a,b,j), 0);
            }
        }
    }
    return value;
}

void ccd::Q3Permutations(){
    for(unsigned ch = 0; ch < ChannelVariety.size(); ch++){
        channel onechannel = ChannelVariety.at(ch);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            symblock TOrigBlock = m_pphhTBlock.at(ch);
            for(unsigned ab = 0; ab < ChannelVariety.at(ch).m_ParticleParticleVec.size(); ab++){
                for(unsigned ij = 0; ij < ChannelVariety.at(ch).m_HoleHoleVec.size(); ij++){
                    double value = 0.0;
                    int a = TOrigBlock.getRowMap(ab).first();
                    int b = TOrigBlock.getRowMap(ab).second();
                    int i = TOrigBlock.getColMap(ij).first();
                    int j = TOrigBlock.getColMap(ij).second();
                    value = Q3recoupled(a,b,j,i) - Q3recoupled(a,b,i,j);
                    m_pphhTBlock.at(ch).setElement(ab, ij, value + m_pphhTBlock.at(ch).getElement(ab,ij));
                }
            }
        }
    }
}

double ccd::Q4recoupled(int p, int q, int r, int s){

    double value = 0.0;
    for(unsigned int chan = 0; chan < ChannelVariety2.size(); chan++){
        channel onechannel = ChannelVariety2.at(chan);
        if (onechannel.m_HolePlusHoleMinusParticleVec.size()*onechannel.m_ParticleVec.size() != 0){
            channelindexpair RSQ = onechannel.m_HolePlusHoleMinusParticleVec.at(0);
            if ( qsys->isEqual(qsys->sumSubstractState(r,s,q) , qsys->oneState(p)) && qsys->isEqual( qsys->oneState(p), qsys->sumSubstractState(RSQ.first(),RSQ.second(), RSQ.third())) ){
                value = m_hhMpTBlock.at(chan).getElement(0, m_hhMpTBlock.at(chan).getColMapInverse3x1(r,s,q));
            }
        }
    }
    return value;
}

void ccd::Q4Permutations(){
     for(unsigned int ch = 0; ch < ChannelVariety.size(); ch++){
         channel onechannel = ChannelVariety.at(ch);
         if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
             symblock TOrigBlock = m_pphhTBlock.at(ch);
             for(unsigned int ab = 0; ab < ChannelVariety.at(ch).m_ParticleParticleVec.size(); ab++){
                 for(unsigned int ij = 0; ij < ChannelVariety.at(ch).m_HoleHoleVec.size(); ij++){
                     double value = 0.0;
                     int a = TOrigBlock.getRowMap(ab).first();
                     int b = TOrigBlock.getRowMap(ab).second();
                     int i = TOrigBlock.getColMap(ij).first();
                     int j = TOrigBlock.getColMap(ij).second();
                     value = -Q4recoupled(a,b,i,j) + Q4recoupled(b,a,i,j);
                     m_pphhTBlock.at(ch).setElement(ab, ij, value + m_pphhTBlock.at(ch).getElement(ab,ij));
                 }
             }
         }
     }
 }

void ccd::printAmplitudes(){
    std::cout << std::setprecision(16);
    for(unsigned int blockInd = 0; blockInd < m_pphhTBlock.size(); blockInd++){
        symblock oneblock = m_pphhTBlock.at(blockInd);
        if (oneblock.getMatrix().size() != 0){
            for(unsigned int ab = 0; ab < ChannelVariety.at(blockInd).m_ParticleParticleVec.size(); ab++){
                for(unsigned int ij = 0; ij < ChannelVariety.at(blockInd).m_HoleHoleVec.size(); ij++){
                    std::cout << "a b i j " << oneblock.getRowMap(ab).first() << " " << oneblock.getRowMap(ab).second() << " " << oneblock.getColMap(ij).first() << " " <<  oneblock.getColMap(ij).second()  <<   " " << oneblock.getMatrix()(ab,ij)<< std::endl  ;
                }
            }
        }
    }
}

void ccd::printAmplitudesPrev(){
    std::cout << std::setprecision(16);
    for(unsigned int blockInd = 0; blockInd < m_pphhTBlockPrev.size(); blockInd++){
        symblock oneblock = m_pphhTBlockPrev.at(blockInd);
        if (oneblock.getMatrix().size() != 0){
            for(unsigned int ab = 0; ab < ChannelVariety.at(blockInd).m_ParticleParticleVec.size(); ab++){
                for(unsigned int ij = 0; ij < ChannelVariety.at(blockInd).m_HoleHoleVec.size(); ij++){
                    std::cout << "a b i j " << oneblock.getRowMap(ab).first() << " " << oneblock.getRowMap(ab).second() << " " << oneblock.getColMap(ij).first() << " " <<  oneblock.getColMap(ij).second()  <<   " " << oneblock.getMatrix()(ab,ij)<< std::endl  ;
                }
            }
        }
    }
}

double ccd::L3recoupled(int a, int b, int i, int j){

    double value = 0.0;
    for(unsigned int chan = 0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleMinusHoleVec.size()*onechannel.m_HoleMinusParticleVec.size() != 0){
            channelindexpair PQ = onechannel.m_ParticleMinusHoleVec.at(0);
            if (  qsys->isEqual(qsys->substractState(a,i) , qsys->substractState(j,b)) && qsys->isEqual(qsys->substractState(a,i) ,qsys->substractState(PQ.first(),PQ.second())) ){
                value = m_pMhhMpTBlock.at(chan).getElement(m_pMhhMpTBlock.at(chan).getRowMapInverse(a, i),  m_pMhhMpTBlock.at(chan).getColMapInverse(j,b));
            }
        }
    }
    return value;
}

void ccd::L3Permutations(){
    for(unsigned ch = 0; ch < ChannelVariety.size(); ch++){
        channel onechannel = ChannelVariety.at(ch);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            symblock TOrigBlock = m_pphhTBlock.at(ch);
            for(unsigned ab = 0; ab < ChannelVariety.at(ch).m_ParticleParticleVec.size(); ab++){
                for(unsigned ij = 0; ij < ChannelVariety.at(ch).m_HoleHoleVec.size(); ij++){
                    double value = 0.0;
                    int a = TOrigBlock.getRowMap(ab).first();
                    int b = TOrigBlock.getRowMap(ab).second();
                    int i = TOrigBlock.getColMap(ij).first();
                    int j = TOrigBlock.getColMap(ij).second();
                    value = L3recoupled(a,b,i,j) - L3recoupled(b,a,i,j) - L3recoupled(a,b,j,i) + L3recoupled(b,a,j,i);
                    m_pphhTBlock.at(ch).setElement(ab, ij, value + m_pphhTBlock.at(ch).getElement(ab,ij));
                }
            }
        }
    }
}

double ccd::Q2recoupled(int a, int b, int i, int j){

    double value = 0.0;
    for(unsigned int chan = 0; chan < ChannelVariety.size(); chan++){
        channel onechannel = ChannelVariety.at(chan);
        if (onechannel.m_ParticleMinusHoleVec.size()*onechannel.m_HoleMinusParticleVec.size() != 0){
            channelindexpair PQ = onechannel.m_ParticleMinusHoleVec.at(0);
            if (  qsys->isEqual(qsys->substractState(a,j) , qsys->substractState(i,b) ) && qsys->isEqual(qsys->substractState(a,j) , qsys->substractState(PQ.first(),PQ.second())) ){
                value = m_pMhhMpTBlock.at(chan).getElement(m_pMhhMpTBlock.at(chan).getRowMapInverse(a, j), m_pMhhMpTBlock.at(chan).getColMapInverse(i,b));
            }
        }
    }
    return value;
}

void ccd::Q2Permutations(){
    for(unsigned int ch = 0; ch < ChannelVariety.size(); ch++){
        channel onechannel = ChannelVariety.at(ch);
        if (onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            symblock TOrigBlock = m_pphhTBlock.at(ch);
            for(unsigned ab = 0; ab < ChannelVariety.at(ch).m_ParticleParticleVec.size(); ab++){
                for(unsigned ij = 0; ij < ChannelVariety.at(ch).m_HoleHoleVec.size(); ij++){
                    double value = 0.0;

                    int a = TOrigBlock.getRowMap(ab).first();
                    int b = TOrigBlock.getRowMap(ab).second();
                    int i = TOrigBlock.getColMap(ij).first();
                    int j = TOrigBlock.getColMap(ij).second();

                    value = Q2recoupled(a,b,j,i)- Q2recoupled(a,b,i,j);

                    m_pphhTBlock.at(ch).setElement(ab, ij, value + m_pphhTBlock.at(ch).getElement(ab,ij));
                }
            }
        }
    }
}

