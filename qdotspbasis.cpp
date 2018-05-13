#include "qdotspbasis.h"


bool myComparison(const std::pair<int,double> &a,const std::pair<int,double> &b)
{
       return a.second<b.second;
}


qdotspbasis::qdotspbasis(int NumberOfShellsStochastic, int NumberOfShellsExact, int ParticlesNumber, double HOStrenth)
{
    this->homega = HOStrenth;
    this->m_FermiLevel = ParticlesNumber;
    this->m_ShellsExact = NumberOfShellsExact;
    this->m_ShellsStochastic = NumberOfShellsStochastic;
    setUpStatesPolarSorted();
}

void qdotspbasis::setUpStatesPolarSorted() {

    qstate m_q_state;
    std::vector<int> oddShells;
    int n, m;
    int m_sm = -1;

    //loop to find all odd shell numbers including energyCutOff shell
    for (int i = 1; i <= this->m_ShellsStochastic; i++){
        if (i % 2 != 0){
            oddShells.push_back(i);
        }
    }
    for (unsigned int index = 0; index < oddShells.size(); ++index){
       //positive m
       for (int j=0; j<=this->m_ShellsStochastic - oddShells[index]; j++){
           n = index;
           m = j;
           m_q_state.set(n, m, m_sm);
           m_shells.push_back(m_q_state);
       }
       //negative m
       for (int j=1; j<=this->m_ShellsStochastic - oddShells[index]; j++){
           n = index;
           m = -1*j;
           m_q_state.set(n, m, m_sm);
           m_shells.push_back(m_q_state);
       }
    }

    //now sort m_shells vector
    std::pair<int,double> mapping;
    std::vector<std::pair<int,double>> vector_to_sort;
    for(unsigned int i = 0; i < m_shells.size(); i++) {
        qstate quantum_state = m_shells.at(i);
        double EnergyOfState = homega*((double)2.0*quantum_state.n() + std::abs(quantum_state.m()) + 1.0);
        mapping = std::make_pair(i, EnergyOfState);
        vector_to_sort.push_back(mapping);
    }

    sort(vector_to_sort.begin(),vector_to_sort.end(),myComparison);
    std::vector<qstate> sorted_states;
    for(unsigned int i = 0; i < vector_to_sort.size(); i++) {
        sorted_states.push_back(m_shells.at(vector_to_sort.at(i).first));
        m_shells.at(vector_to_sort.at(i).first).flipSpin();
        sorted_states.push_back(m_shells.at(vector_to_sort.at(i).first));
    }
    m_shells = sorted_states;
    int dim = m_shells.size();
    CalculateSPenergies();
    m_twoBodyElements = create4dArray(dim, dim, dim, dim);
    fillTwoBodyElements();

}

double**** qdotspbasis::create4dArray(int dim1, int dim2, int dim3, int dim4) {

    double**** h = new double***[dim1];

    for(int i=0; i < dim1; i++) {
        h[i] = new double**[dim2];
        for (int j=0; j < dim2; j++) {
            h[i][j] = new double*[dim3];
            for(int k = 0; k < dim3; k++) {
                h[i][j][k] = new double[dim4];
            }
        }
    }
    for(int i=0; i < dim1; i++) {
        for (int j=0; j < dim2; j++) {
            for(int k = 0; k < dim3; k++) {
                for(int l = 0; l < dim4; l++) {
                    h[i][j][k][l] = 0;
                }

            }
        }
    }
    return h;
}

double qdotspbasis::TBME(int p, int q, int r, int s){
    return (m_twoBodyElements[p][q][r][s] - m_twoBodyElements[p][q][s][r]);
}


void qdotspbasis::fillTwoBodyElements(){
    int NumberOfStates = m_shells.size();

    for(int i = 0; i < NumberOfStates; i++) {
        qstate quantum_state_alpha = m_shells.at(i);
        int alpha_n = quantum_state_alpha.n();
        int alpha_m = quantum_state_alpha.m();
        int alpha_sm = quantum_state_alpha.s();

        for(int j = 0; j < NumberOfStates; j++) {
            qstate quantum_state_beta = m_shells.at(j);
            int beta_n = quantum_state_beta.n();
            int beta_m = quantum_state_beta.m();
            int beta_sm = quantum_state_beta.s();

            for(int k = 0; k < NumberOfStates; k++) {
                qstate quantum_state_gama = m_shells.at(k);
                int gama_n = quantum_state_gama.n();
                int gama_m = quantum_state_gama.m();
                int gama_sm = quantum_state_gama.s();

                for(int l = 0; l < NumberOfStates; l++) {
                    qstate quantum_state_delta = m_shells.at(l);
                    int delta_n = quantum_state_delta.n();
                    int delta_m = quantum_state_delta.m();
                    int delta_sm = quantum_state_delta.s();

                    if ((alpha_sm == beta_sm && gama_sm == delta_sm)){
                        m_twoBodyElements[i][k][j][l] = Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, beta_n, beta_m,  delta_n, delta_m);
                    }
                    if ((alpha_sm == delta_sm && gama_sm == beta_sm)){
                        m_twoBodyElements[i][k][l][j] = Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, delta_n, delta_m, beta_n, beta_m);
                    }
                }
            }
        }
    }
}

void qdotspbasis::CalculateSPenergies(){
    int NumberOfStates = m_shells.size();
    //m_HOEnergies.zeros(NumberOfStates, NumberOfStates);
    for(int i = 0; i < NumberOfStates; i++) {
        qstate quantum_state = m_shells.at(i);
        m_HOEnergies.push_back((2.0*(double)quantum_state.n() + (double)std::abs(quantum_state.m()) + 1.0)*homega);
    }
}

std::vector<double> qdotspbasis::getSPenergies(){
    return m_HOEnergies;
}

qstate* qdotspbasis::oneState(int p){
    int N = getStateVec().at(p).n();
    int M = getStateVec().at(p).m();
    int S =  getStateVec().at(p).s();
    qstate *QuantumState = new qstate();
    QuantumState->set(N, M, S);
    return QuantumState;
}

qstate* qdotspbasis::sumState(int p, int q){
    int N = getStateVec().at(p).n() + getStateVec().at(q).n();
    int M = getStateVec().at(p).m() + getStateVec().at(q).m();
    int S = getStateVec().at(p).s() + getStateVec().at(q).s();
    qstate *QuantumState = new qstate();
    QuantumState->set(N, M, S);
    return QuantumState;
}

qstate* qdotspbasis::substractState(int p, int q){
    int N = getStateVec().at(p).n() - getStateVec().at(q).n();
    int M = getStateVec().at(p).m() - getStateVec().at(q).m();
    int S = getStateVec().at(p).s() - getStateVec().at(q).s();
    qstate *QuantumState = new qstate();
    QuantumState->set(N, M, S);
    return QuantumState;
}

qstate* qdotspbasis::sumSubstractState(int p, int q, int r){
    int N = getStateVec().at(p).n() + getStateVec().at(q).n() - getStateVec().at(r).n();
    int M = getStateVec().at(p).m() + getStateVec().at(q).m() - getStateVec().at(r).m();
    int S = getStateVec().at(p).s() + getStateVec().at(q).s() - getStateVec().at(r).s();
    qstate *QuantumState = new qstate();
    QuantumState->set(N, M, S);
    return QuantumState;
}

bool qdotspbasis::isEqual(qstate* state1, qstate* state2){
    if (   //state1->n() == state2->n() &&
           state1->m() == state2->m()
        && state1->s()  == state2->s()
            ) {
        return true;
    } else {
        return false;
    }
}


void qdotspbasis::getQuantumDotStates(){
    for(qstate quantum_state : m_shells){
        std::cout << "n = " << quantum_state.n() << std::endl;
        std::cout << "m = " <<quantum_state.m() << std::endl;
        std::cout << "s = " <<quantum_state.s() << std::endl;
        std::cout << "-----------------" << std::endl;
    }
    printSPenergies();
}

void qdotspbasis::getQuantumDotStatesNumber(){
    std::cout << "Number of available states of system is " << m_shells.size() << std::endl;
}

void qdotspbasis::printSPenergies(){
    for (unsigned int i = 0; i < m_HOEnergies.size(); i++){
        std::cout << m_HOEnergies.at(i) << std::endl;
    }
}
