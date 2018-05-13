#include "channelindexpair.h"

channelindexpair::channelindexpair()
{
}

void channelindexpair::set(int p, int q){
    m_p = p;
    m_q = q;

}

void channelindexpair::setOne(int p){
    m_p = p;
}

void channelindexpair::setThree(int p, int q, int r){
    m_p = p;
    m_q = q;
    m_r = r;

}

void channelindexpair::setMap(int index, int p, int q){
    m_p = p;
    m_q = q;
    m_index = index;

}

void channelindexpair::setMapOne(int index, int p){
    m_p = p;
    m_index = index;

}

void channelindexpair::setMapThree(int index, int p, int q, int r){
    m_p = p;
    m_q = q;
    m_r = r;
    m_index = index;

}
