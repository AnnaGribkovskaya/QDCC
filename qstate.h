#ifndef QSTATE_H
#define QSTATE_H


class qstate
{
private:
    int m_n;     // n
    int m_m;     // m
    int m_sm;     // spin projection (1 or -1)

public:
    qstate();
    int n() { return m_n; }
    int m() { return m_m; }
    int s() { return m_sm; }
    void set(int n, int m, int sm);
    void flipSpin();
};

#endif // QSTATE_H
