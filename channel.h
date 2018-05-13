#ifndef CHANNEL_H
#define CHANNEL_H

#include <vector>
#include "channelindexpair.h"
#include <armadillo>
#include "symblock.h"


class channel
{
public:
    channel();
    std::vector<channelindexpair> m_HoleHoleVec;
    std::vector<channelindexpair> m_ParticleParticleVec;
    std::vector<channelindexpair> m_ParticleHoleVec;
    std::vector<channelindexpair> m_ParticleMinusHoleVec;
    std::vector<channelindexpair> m_HoleParticleVec;
    std::vector<channelindexpair> m_HoleMinusParticleVec;
    std::vector<channelindexpair> m_ParticlePlusParticleMinusHoleVec;
    std::vector<channelindexpair> m_HoleVec;
    std::vector<channelindexpair> m_HolePlusHoleMinusParticleVec;
    std::vector<channelindexpair> m_ParticleVec;
};

#endif // CHANNEL_H
