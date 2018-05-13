#include "qstate.h"

qstate::qstate()
{

}

void qstate::set(int n, int m, int sm){
    m_n = n;
    m_m = m;
    m_sm = sm;
}

void qstate::flipSpin(){
    m_sm = m_sm*(-1);
}
