#ifndef QDOTCHANNELSET_H
#define QDOTCHANNELSET_H

#include "channelset.h"
#include "channel.h"
#include "qstate.h"

class qdotchannelset : public channelset
{
public:
    qdotchannelset();
    virtual void setUpChannels(generalSPclass *);
    virtual void setUpChannelsQ3(generalSPclass *);
    virtual void setUpChannelsQ4(generalSPclass *);
    virtual std::vector<channel> getChan () {return this->ChannelVariety;}
    virtual std::vector<channel> getChanQ3 () {return this->ChannelVariety1;}
    virtual std::vector<channel> getChanQ4 () {return this->ChannelVariety2;}

private:
    generalSPclass * qsys;
    std::vector<channel>  ChannelVariety;
    std::vector<channel>  ChannelVariety1;
    std::vector<channel>  ChannelVariety2;
};

#endif // QDOTCHANNELSET_H
