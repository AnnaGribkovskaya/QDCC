#include "mbpt2.h"


mbpt2::mbpt2(generalSPclass * qsystem){
    qsys = qsystem;
}


void mbpt2::setUpChannels(){
    int i = 0;
    int Mlmax = 2*(qsys->getShellsStochastic()-1);
    int Smax = 2;
    for(int Ml = -Mlmax; Ml <= Mlmax; Ml++){
      for(int S = -Smax; S <= Smax; S = S + 2){
        i++;
        ChannelVariety.emplace_back(channel());

        for(int i = 0; i < qsys->getFermiLevel(); i++){
          for(int j = 0; j < qsys->getFermiLevel(); j++){
            qstate *QuantumState = qsys->sumState(i,j);
            if(Ml == QuantumState->m() && S == QuantumState->s() && i != j){
              ChannelVariety.back().m_HoleHoleVec.emplace_back(channelindexpair());
              ChannelVariety.back().m_HoleHoleVec.back().set(i, j);
            } delete QuantumState;
          }
        }
        for(int a = qsys->getFermiLevel(); a < qsys->getStatesStochastic(); a++){
          for(int b = qsys->getFermiLevel(); b < qsys->getStatesStochastic(); b++){
            qstate *QuantumState = qsys->sumState(a,b);
            if(Ml == QuantumState->m() && S == QuantumState->s() && a != b){
              ChannelVariety.back().m_ParticleParticleVec.emplace_back(channelindexpair());
              ChannelVariety.back().m_ParticleParticleVec.back().set(a, b);
            } delete QuantumState;
          }
        }
      }
    }
}

/*double mbpt2::calculateEnergyBruteForce(){
    setUpChannels();
    std::vector<CartesianState> shells = qsys->getStateVec();
    double CorrelationEnergy = 0.0;
    arma::vec Epsilon;
    Epsilon.zeros(shells.size());
    Epsilon = qsys->computeSPenergies(shells);

    for(channel onechannel : ChannelVariety){
      for(size_t i = 0; i < onechannel.m_HoleHoleVec.size(); i++){
        channelindexpair ChannelHoleHoleIndexes = onechannel.m_HoleHoleVec.at(i);
        for(size_t a = 0; a < onechannel.m_ParticleParticleVec.size(); a++){
          channelindexpair ChannelParticleParticleIndexes = onechannel.m_ParticleParticleVec.at(a);
          CorrelationEnergy += ((qsys->TBME( ChannelHoleHoleIndexes.first(), ChannelHoleHoleIndexes.second(), ChannelParticleParticleIndexes.first(), ChannelParticleParticleIndexes.second(), shells)*qsys->TBME(ChannelParticleParticleIndexes.first(), ChannelParticleParticleIndexes.second(), ChannelHoleHoleIndexes.first(), ChannelHoleHoleIndexes.second(), shells))/(Epsilon(ChannelHoleHoleIndexes.first()) + Epsilon(ChannelHoleHoleIndexes.second()) - Epsilon(ChannelParticleParticleIndexes.first()) - Epsilon(ChannelParticleParticleIndexes.second())));

        }
      }
    }
    return 0.25*CorrelationEnergy;
}*/

double mbpt2::calculateEnergyMatMult(){
    setUpChannels();
    double CorrelationEnergy = 0.0;
    arma::vec Epsilon;
    Epsilon.zeros(qsys->getStatesStochastic());
    Epsilon = qsys->getSPenergies();

    for(channel onechannel : ChannelVariety){
        if(onechannel.m_ParticleParticleVec.size()*onechannel.m_HoleHoleVec.size() != 0){
            arma::mat A((onechannel.m_HoleHoleVec.size()), onechannel.m_ParticleParticleVec.size());
            arma::mat B(onechannel.m_ParticleParticleVec.size(), onechannel.m_HoleHoleVec.size());
            for(size_t i = 0; i < onechannel.m_HoleHoleVec.size(); i++){
                channelindexpair ChannelHoleHoleIndexes = onechannel.m_HoleHoleVec.at(i);
                for(size_t a = 0; a < onechannel.m_ParticleParticleVec.size(); a++){
                    channelindexpair ChannelParticleParticleIndexes = onechannel.m_ParticleParticleVec.at(a);
                    A(i,a) = qsys->TBME( ChannelHoleHoleIndexes.first(), ChannelHoleHoleIndexes.second(), ChannelParticleParticleIndexes.first(), ChannelParticleParticleIndexes.second());
                    B(a,i) = qsys->TBME(ChannelParticleParticleIndexes.first(), ChannelParticleParticleIndexes.second(), ChannelHoleHoleIndexes.first(), ChannelHoleHoleIndexes.second())/(Epsilon(ChannelHoleHoleIndexes.first()) + Epsilon(ChannelHoleHoleIndexes.second()) - Epsilon(ChannelParticleParticleIndexes.first()) - Epsilon(ChannelParticleParticleIndexes.second()));
                }
            }
            CorrelationEnergy += 0.25*trace(A*B);
        }
    }
    return CorrelationEnergy;
}

