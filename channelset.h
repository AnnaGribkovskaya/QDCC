#ifndef CHANNELSET_H
#define CHANNELSET_H

#include "generalspclass.h"
#include "channel.h"

class channelset
{
public:
    channelset();
    virtual void setUpChannels(generalSPclass *) = 0;
    virtual void setUpChannelsQ3(generalSPclass *) = 0;
    virtual void setUpChannelsQ4(generalSPclass *) = 0;

    virtual std::vector<channel> getChan () = 0;
    virtual std::vector<channel> getChanQ3 () = 0;
    virtual std::vector<channel> getChanQ4 () = 0;
};

#endif // CHANNELSET_H
