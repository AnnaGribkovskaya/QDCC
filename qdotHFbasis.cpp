#include "qdotHFbasis.h"

bool Comparison(const std::pair<int,double> &a,const std::pair<int,double> &b)
{
       return a.second<b.second;
}

qdotHFbasis::qdotHFbasis(int NumberOfShellsStochastic, int NumberOfShellsExact, int ParticlesNumber, double HOStrenth)
{
    this->homega = HOStrenth;
    this->m_FermiLevel = ParticlesNumber;
    this->m_ShellsExact = NumberOfShellsExact;
    this->m_ShellsStochastic = NumberOfShellsStochastic;
    setUpStatesPolarSorted();
    this->m_StatesStochastic = m_shells.size();
    test();
    //getQuantumDotStates();
}

void qdotHFbasis::setUpStatesPolarSorted() {

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

    sort(vector_to_sort.begin(),vector_to_sort.end(),Comparison);
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

double**** qdotHFbasis::create4dArray(int dim1, int dim2, int dim3, int dim4) {

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

// PEREPISAT!!!!!!! TBME!!!!

double qdotHFbasis::TBME(int p, int q, int r, int s){
    return (m_twoBodyElements[p][q][r][s] - m_twoBodyElements[p][q][s][r]);
    //return m_twoBodyElements[p][q][r][s];
}

void qdotHFbasis::fillTwoBodyElements(){
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
            //for(int k = 0; k < i; k++) {
                qstate quantum_state_gama = m_shells.at(k);
                int gama_n = quantum_state_gama.n();
                int gama_m = quantum_state_gama.m();
                int gama_sm = quantum_state_gama.s();

                for(int l = 0; l < NumberOfStates; l++) {
                //for(int l = 0; l < j; l++) {
                    qstate quantum_state_delta = m_shells.at(l);
                    int delta_n = quantum_state_delta.n();
                    int delta_m = quantum_state_delta.m();
                    int delta_sm = quantum_state_delta.s();



                     /*if ( (alpha_sm + beta_sm == gama_sm + delta_sm) ){ //&& (alpha_m + beta_m == gama_m + delta_m)
                        tbme1 = Coulomb_HO(homega, alpha_n, alpha_m, beta_n, beta_m, gama_n, gama_m, delta_n, delta_m);
                        tbme2 = Coulomb_HO(homega, alpha_n, alpha_m, beta_n, beta_m, delta_n, delta_m, gama_n, gama_m);

                        if (alpha_sm == gama_sm && beta_sm  == delta_sm && alpha_sm == delta_sm  &&  beta_sm  == gama_sm){
                            m_twoBodyElements[i][j][k][l] = tbme1 - tbme2;
                        }
                        if (alpha_sm == delta_sm  &&  beta_sm  == gama_sm ){
                            m_twoBodyElements[i][j][k][l] =  - tbme2;
                        }
                        if (alpha_sm == gama_sm && beta_sm  == delta_sm){
                            m_twoBodyElements[i][j][k][l] = tbme1;
                        }

                    }*/

                    if (alpha_sm == beta_sm && gama_sm == delta_sm ){
                        m_twoBodyElements[i][k][j][l] = Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, beta_n, beta_m,  delta_n, delta_m);
                    }
                    if ((alpha_sm == delta_sm && gama_sm == beta_sm )){
                        m_twoBodyElements[i][k][l][j] = Coulomb_HO(homega, alpha_n, alpha_m, gama_n, gama_m, delta_n, delta_m, beta_n, beta_m);
                    }

                }
            }
        }
    }
    std::cout << "Two body matrix elements are filled!" << std::endl;
}

void qdotHFbasis::CalculateSPenergies(){
    int NumberOfStates = m_shells.size();
    //m_HOEnergies.zeros(NumberOfStates, NumberOfStates);
    for(int i = 0; i < NumberOfStates; i++) {
        qstate quantum_state = m_shells.at(i);
        m_HOEnergies.push_back((2.0*(double)quantum_state.n() + (double)std::abs(quantum_state.m()) + 1.0)*homega);
    }
}

std::vector<double> qdotHFbasis::getSPenergies(){
    return m_HOEnergies;
}

qstate* qdotHFbasis::oneState(int p){
    int N = getStateVec().at(p).n();
    int M = getStateVec().at(p).m();
    int S =  getStateVec().at(p).s();
    qstate *QuantumState = new qstate();
    QuantumState->set(N, M, S);
    return QuantumState;
}

qstate* qdotHFbasis::sumState(int p, int q){
    int N = getStateVec().at(p).n() + getStateVec().at(q).n();
    int M = getStateVec().at(p).m() + getStateVec().at(q).m();
    int S = getStateVec().at(p).s() + getStateVec().at(q).s();
    qstate *QuantumState = new qstate();
    QuantumState->set(N, M, S);
    return QuantumState;
}

qstate* qdotHFbasis::substractState(int p, int q){
    int N = getStateVec().at(p).n() - getStateVec().at(q).n();
    int M = getStateVec().at(p).m() - getStateVec().at(q).m();
    int S = getStateVec().at(p).s() - getStateVec().at(q).s();
    qstate *QuantumState = new qstate();
    QuantumState->set(N, M, S);
    return QuantumState;
}

qstate* qdotHFbasis::sumSubstractState(int p, int q, int r){
    int N = getStateVec().at(p).n() + getStateVec().at(q).n() - getStateVec().at(r).n();
    int M = getStateVec().at(p).m() + getStateVec().at(q).m() - getStateVec().at(r).m();
    int S = getStateVec().at(p).s() + getStateVec().at(q).s() - getStateVec().at(r).s();
    qstate *QuantumState = new qstate();
    QuantumState->set(N, M, S);
    return QuantumState;
}

bool qdotHFbasis::isEqual(qstate* state1, qstate* state2){
    if (   //state1->n() == state2->n() &&
           state1->m() == state2->m()
        && state1->s()  == state2->s()
            ) {
        return true;
    } else {
        return false;
    }
}

void qdotHFbasis::getQuantumDotStates(){
    int i = 0;
    for(qstate quantum_state : m_shells){
        std::cout << "state Number = " << i << std::endl;
        std::cout << "n = " << quantum_state.n() << std::endl;
        std::cout << "m = " <<quantum_state.m() << std::endl;
        std::cout << "s = " <<quantum_state.s() << std::endl;
        std::cout << "-----------------" << std::endl;
        i++;
    }
    printSPenergies();
}

void qdotHFbasis::getQuantumDotStatesNumber(){
    std::cout << "Number of available states of system is " << m_shells.size() << std::endl;
}

void qdotHFbasis::printSPenergies(){
    for (unsigned int i = 0; i < m_HOEnergies.size(); i++){
        std::cout << m_HOEnergies.at(i) << std::endl;
    }
}

void qdotHFbasis::test(){
    //int N = m_FermiLevel;
    double val = 0.0;
    for(int i = 0; i < m_FermiLevel; i++) {
        val +=  m_HOEnergies.at(i);
        for(int j = 0; j < m_FermiLevel; j++) {
           val += 0.5*TBME(i,j,i,j);
        }
     }
    std::cout << std::setprecision(16) << " initial energy: " << val << std::endl;
    printSPenergies();
    applyHartreeFockMethod();
}


// HERE START THE HARTREE-FOCK PART

void qdotHFbasis::setCoefficientMatrix(arma::mat CoefficientMatrix){
     m_C = CoefficientMatrix;
}

arma::mat qdotHFbasis::computeDensityMatrix(){
    int NumberOfStates = m_shells.size();
    arma::mat DensityMatrix(NumberOfStates, NumberOfStates);

    for(int k = 0; k < NumberOfStates; k++) {
        for(int l = 0; l < NumberOfStates; l++) {
            double sum = 0.0;
            for (int i=0; i < m_FermiLevel; i++) {
                    sum += m_C(k,i)*m_C(l,i);
                    DensityMatrix(k, l) = sum;
            }
        }
    }
    DensityMatrix.print();
    return DensityMatrix;
}
/* NOT NEEDED FUNCTION GETSPENERGIES DOES THE SAME!

void qdotHFbasis::CalculateNonIntEnergy(){
    int NumberOfStates = m_shells.size();
    m_HOEnergies.zeros(NumberOfStates,NumberOfStates);
    for(int i = 0; i < NumberOfStates; i++) {
        QuantumState quantum_state = m_shells.at(i);
        m_HOEnergies(i, i) = (2.0*(double)quantum_state.n() + (double)abs(quantum_state.m()) + 1.0)*homega;
    }
}*/

void qdotHFbasis::computeHFmatrix(arma::mat DensityMatrix){
    int NumberOfStates = m_shells.size();
    m_HF.zeros(NumberOfStates,NumberOfStates);
    double FockElement = 0;

    for(int i = 0; i < NumberOfStates; i++) {
        //qstate quantum_state_alpha = m_shells.at(i);
        //int alpha_n = quantum_state_alpha.n();
        //int alpha_m = quantum_state_alpha.m();
        //int alpha_sm = quantum_state_alpha.s();

        for(int j = 0; j < NumberOfStates; j++) {
            //qstate quantum_state_beta = m_shells.at(j);
            //int beta_n = quantum_state_beta.n();
            //int beta_m = quantum_state_beta.m();
            //int beta_sm = quantum_state_beta.s();

            for(int k = 0; k < NumberOfStates; k++) {
                //qstate quantum_state_gama = m_shells.at(k);
                //int gama_n = quantum_state_gama.n();
                //int gama_m = quantum_state_gama.m();
                //int gama_sm = quantum_state_gama.s();

                for(int l = 0; l < NumberOfStates; l++) {
                    //qstate quantum_state_delta = m_shells.at(l);
                    //int delta_n = quantum_state_delta.n();
                    //int delta_m = quantum_state_delta.m();
                    //int delta_sm = quantum_state_delta.s();
                    //double TBME = 0.0;
                    //double tbme1 = 0.0;
                    //double tbme2 = 0.0;
/*
                    if ((alpha_sm == beta_sm && gama_sm == delta_sm)){
                        tbme1 = m_twoBodyElements[i][k][j][l];
                    }
                    if ((alpha_sm == delta_sm && gama_sm == beta_sm)){
                        tbme2 = m_twoBodyElements[i][k][l][j];
                    }
                    */
                    //TBME = tbme1 - tbme2;
                    FockElement += DensityMatrix(k,l)*TBME(i,k,j,l);
                    if (FockElement !=0.0){
                        //std::cout << "FOOOKKKK   " << FockElement << "  i " << i << " j  "<< j << "  k " << k << " l  "<< j << std::endl;
                    }
                }
            }
            if (i == j) {
                m_HF(i, i) += m_HOEnergies.at(i);
                //std::cout << "SPENERGY   " << m_HOEnergies.at(i) << "  i " << i << std::endl;
            }
            m_HF(i, j) += FockElement;
            FockElement = 0.0;

        }
    }
}

double qdotHFbasis::computeHartreeFockEnergyDifference(){
    return ((arma::accu(abs(eigval - eigval_previous)))/(double)m_shells.size());
}



void qdotHFbasis::computeHartreeFockEnergy(arma::mat DensityMatrix){
    int NumberOfStates = m_shells.size();
    //int FermiLevel = m_FermiLevel;
    double selfConsistentFIeldIterations = 0.0;
    double ExchangePart = 0.0;
    double SingleParticleEnergies = 0.0;

    for(int f = 0; f < m_FermiLevel; f++){
        SingleParticleEnergies += eigval(f);
        std::cout << "SPEnergies " << SingleParticleEnergies << std::endl;
    }

    for(int i = 0; i < NumberOfStates; i++) {
        //qstate quantum_state_alpha = m_shells.at(i);
        //int alpha_n = quantum_state_alpha.n();
        //int alpha_m = quantum_state_alpha.m();
        //int alpha_sm = quantum_state_alpha.s();

        for(int j = 0; j < NumberOfStates; j++) {
            //qstate quantum_state_beta = m_shells.at(j);
            //int beta_n = quantum_state_beta.n();
            //int beta_m = quantum_state_beta.m();
            //int beta_sm = quantum_state_beta.s();

            for(int k = 0; k < NumberOfStates; k++) {
                //qstate quantum_state_gama = m_shells.at(k);
                //int gama_n = quantum_state_gama.n();
                //int gama_m = quantum_state_gama.m();
                //int gama_sm = quantum_state_gama.s();

                for(int l = 0; l < NumberOfStates; l++) {
                    //qstate quantum_state_delta = m_shells.at(l);
                    //int delta_n = quantum_state_delta.n();
                    //int delta_m = quantum_state_delta.m();
                    //int delta_sm = quantum_state_delta.s();
/*
                    double TBME = 0.0;
                    double tbme1 = 0.0;
                    double tbme2 = 0.0;
                    if ((alpha_sm == beta_sm) && (gama_sm == delta_sm)){
                       tbme1 = m_twoBodyElements[i][k][j][l];
                    }
                    if ((alpha_sm == delta_sm) && (gama_sm == beta_sm)){
                       tbme2 = m_twoBodyElements[i][k][l][j];
                    }
                    TBME = tbme1 - tbme2;
                    */
                    selfConsistentFIeldIterations = DensityMatrix(i,j)*DensityMatrix(k,l)*TBME(i,k,j,l);
                    ExchangePart += selfConsistentFIeldIterations;
                }
            }
        }
    }
    double HF_Energy = SingleParticleEnergies - 0.5*ExchangePart;
    // Uncoment for debug
    std::cout << "SPEnergies " << SingleParticleEnergies << std::endl;
    std::cout << "Exchange " << ExchangePart << std::endl;
    std::cout << "===================================================================" << std::endl;
    std::cout << std::setprecision(12);
    std::cout << "Num of electrons = " << m_FermiLevel << std::endl;
    //std::cout << "Num of shells = " << m_EnergyCutOff << std::endl;
    std::cout << "Omega = " << homega << std::endl;
    std::cout << "Total energy " << HF_Energy << std::endl;
    //writeToFile(HF_Energy, NumberOfParticles, m_EnergyCutOff, homega);
}



void qdotHFbasis::applyHartreeFockMethod(){
    int NumberOfStates = m_shells.size();
    std::cout << "DIMENSION    " << NumberOfStates << std::endl;
    arma::mat C(NumberOfStates, NumberOfStates);

    C.eye();
    setCoefficientMatrix(C);
    double difference = 10; //dummy value to handle first iteration
    double epsilon = 10e-8;

    eigval_previous.zeros(NumberOfStates);
    int i = 0;
    while (epsilon < difference && i < 200){
        arma::mat x_DensityMatrix = computeDensityMatrix();
        computeHFmatrix(x_DensityMatrix);
        //std::cout << "=====***===" <<  std::endl;
        //eigval.print();
        arma::eig_sym(eigval, eigvec, m_HF);
        //std::cout << "=====///===" <<  std::endl;
        //eigval.print();
        setCoefficientMatrix(eigvec);
        difference = computeHartreeFockEnergyDifference();
        eigval_previous = eigval;
        i++;

    }
    m_C.print();
    arma::mat y_DensityMatrix = computeDensityMatrix();

    computeHartreeFockEnergy(y_DensityMatrix);
    std::cout << "Number of iterations " << i << std::endl;


}
