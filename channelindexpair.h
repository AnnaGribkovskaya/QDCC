#ifndef CHANNELINDEXPAIR_H
#define CHANNELINDEXPAIR_H


class channelindexpair
{

public:
    channelindexpair();
    int first() { return m_p; }
    int second() { return m_q; }
    int third() { return m_r; }
    int index() { return m_index; }
    void set(int, int);

    void setMap(int, int, int);
    void setMapOne(int, int);
    void setMapThree(int, int, int, int);

    void setOne(int);
    void setThree(int, int, int);

private:
    int m_p;     // first index
    int m_q;     // second index
    int m_r;     // third index
    int m_index;   // for mapping

};

#endif // CHANNELINDEXPAIR_H
