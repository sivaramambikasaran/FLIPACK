/*!	
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
*/	
/*!	\file FLIPACK.hpp
 header file of class FLIPACK 
*/
#ifndef __FLIPACK_hpp__
#define __FLIPACK_hpp__

#include"BBFMM2D.hpp"
#include"environment.hpp"

using namespace Eigen;
using namespace std;

/*! Fast linear inversion package */
template <typename T>
class FLIPACK{
public:
    //! Constructor of FLIPACK.
    FLIPACK(double* const Htranspose, double* const X, double* const measurements, double* const R, const unsigned long N, const unsigned m, const unsigned p, const unsigned nMeasurementSets, H2_2D_Tree*Atree);
    /*! This function obtains the cross covariance.*/
    void compute_QHtranspose();
    /*! This function obtains the measurement operator corrected covariance. */
    void compute_HQHtranspose();
    /*! This function obtains the measurement corrected covariance. */
    void compute_Psi();
    /*! This functiin obtains the measurement operator corrected structure. */
    void compute_Phi();
    /*! This function obtains the correction term.*/
    void compute_Xi();
    /*! This function obtains unknown drift coefficients. */
    void compute_Beta();
    /*! This function obtains the unknown.
        Running this function will get 
        all members defined in this class.
     */
    void compute_Solution();
    void get_Posterior_Variance();
    
    void get_QHtranspose(double* &QHtranspose);
    void get_HQHtranspose(double* &HQHtranspose);
    void get_Psi(double* &Psi);
    void get_Phi(double* &Phi);
    void get_Xi(double* &Xi);
    void get_Beta(double* &Beta);
    void get_Solution(double* &Solution);
    
    
    MatrixXd QHtranspose;   /*!< Cross covariance */
    MatrixXd HQHtranspose;  /*!< Measurement operator corrected covariance */
    MatrixXd Psi;           /*!< Measurement corrected covariance */
    MatrixXd Phi;           /*!< Measurement operator corrected structure */
    MatrixXd Solution;      /*!< Unknown */
    MatrixXd Xi;            /*!< Correction term */
    MatrixXd Beta;          /*!< Unknown drift coefficients */
    MatrixXd vDiag;
    H2_2D_Tree *Atree;      /*!< Pointer to a tree */
    
private:
    MatrixXd Htranspose;
    MatrixXd X;
    MatrixXd measurements;
    MatrixXd R;
    
    MatrixXd mainMatrix;
    
    
    bool computedQHtranspose, computedHQHtranspose, computedPsi, computedPhi, computedSolution, computedXi, computedBeta, computedVDiag, computedMainMatrix, computedIntermediateSolution;
    
    unsigned long N;            /*!<  Number of unknowns */
    unsigned m;                 /*!<  Number of measurements */
    unsigned nMeasurementSets;  /*!<  Number of sets of measruements*/
    unsigned p;                 /*!<  Number of terms in the structure*/
    
    void get_Main_Matrix();
    
    void get_Intermediate_Solution();
    
    FullPivLU<MatrixXd> lu;
    
    MatrixXd intermediateSolution;
};



template <typename T>
FLIPACK<T>::FLIPACK(double* const Htranspose, double* const X, double* const measurements, double* const R, const unsigned long N, const unsigned m, const unsigned p, const unsigned nMeasurementSets, H2_2D_Tree*Atree){
    this->Htranspose        =   Map<MatrixXd>(Htranspose, N, m);
    this->X                 =   Map<MatrixXd>(X, N, p);
    this->measurements      =   Map<MatrixXd>(measurements, m, nMeasurementSets);
    this->R                 =   Map<MatrixXd>(R, m, m);
    this->Atree             =   Atree;
    
    computedQHtranspose             =   false;
    computedHQHtranspose            =   false;
    computedPsi                     =   false;
    computedPhi                     =   false;
    computedMainMatrix              =   false;
    computedIntermediateSolution    =   false;
    computedXi                      =   false;
    computedBeta                    =   false;
    computedVDiag                   =   false;
    computedSolution                =   false;
    
    this->N                       =   N;
    this->m                       =   m;
    this->p                       =   p;
    this->nMeasurementSets        =   nMeasurementSets;
    
    mainMatrix              =   MatrixXd(m+p,m+p);
    QHtranspose             =   MatrixXd(N,m);
}


template <typename T>
void FLIPACK<T>::compute_QHtranspose(){
    if (computedQHtranspose==false) {
        cout << endl << "Performing FMM to obtain QHtranspose..." << endl;
        QHtranspose            =   MatrixXd(N,m);
        T A;
        double* QHtranspose_;
        QHtranspose_ = new double[N*m];
        A.calculate_Potential(*Atree,QHtranspose_);
        QHtranspose =   Map<MatrixXd>(QHtranspose_, N, m);
        delete[] QHtranspose_;
        computedQHtranspose    =   true;
        cout << endl << "Obtained QHtranspose" << endl;
    }
}


template <typename T>
void FLIPACK<T>::compute_HQHtranspose(){
    if (computedHQHtranspose==false) {
        compute_QHtranspose();
        HQHtranspose   =   Htranspose.transpose()*QHtranspose;
        computedHQHtranspose  =   true;
        cout << "Obtained HQHtranspose" << endl;
    }
}

template <typename T>
void FLIPACK<T>::compute_Psi(){
    if (computedPsi==false) {
        compute_HQHtranspose();
        Psi =   HQHtranspose + R;
        computedPsi            =   true;
        cout << "Obtained Psi" << endl;
    }
}

template <typename T>
void FLIPACK<T>::compute_Phi(){
    if (computedPhi==false) {
        Phi =   Htranspose.transpose()*X;
        computedPhi            =   true;
        cout << "Obtained Phi" << endl;
    }
}

template <typename T>
void FLIPACK<T>::get_Main_Matrix(){
    if (computedMainMatrix==false) {
        compute_Psi();
        mainMatrix.block(0,0,m,m)  =   Psi;
        compute_Phi();
        mainMatrix.block(0,m,m,p)  =   Phi;
        mainMatrix.block(m,0,p,m)  =   Phi.transpose();
        mainMatrix.block(m,m,p,p)  =   MatrixXd::Zero(p,p);
        lu.compute(mainMatrix);
        computedMainMatrix    =   true;
        cout << "Obtained Main matrix" << endl;
    }
}

template <typename T>
void FLIPACK<T>::get_Intermediate_Solution(){
    if (computedIntermediateSolution==false) {
        get_Main_Matrix();
        MatrixXd rhs                        =   MatrixXd::Zero(m+p,nMeasurementSets);
        rhs.block(0,0,m,nMeasurementSets)   =   measurements;
        intermediateSolution               =   lu.solve(rhs);
        computedIntermediateSolution      =   true;
        cout << "Obtained Intermediate solution" << endl;
    }
}

template <typename T>
void FLIPACK<T>::compute_Xi(){
    if (computedXi==false) {
        get_Intermediate_Solution();
        Xi              =   intermediateSolution.block(0,0,m,nMeasurementSets);
        computedXi     =   true;
        cout << "Obtained Xi" << endl;
    }
}

template <typename T>
void FLIPACK<T>::compute_Beta(){
    if (computedXi==false) {
        get_Intermediate_Solution();
        Beta            =   intermediateSolution.block(m,0,p,nMeasurementSets);
        computedBeta   =   true;
        cout << "Obtained Beta" << endl;
    }
}

template <typename T>
void FLIPACK<T>::compute_Solution(){
    if (computedSolution==false) {
        compute_QHtranspose();
        compute_Beta();
        compute_Xi();
        Solution            =   X*Beta +   QHtranspose*Xi;
        computedSolution   =   true;
        cout << "Obtained Solution" << endl;
    }
}

template <typename T>
void FLIPACK<T>::get_QHtranspose(double*& QHtranspose_) {
    compute_QHtranspose();
    QHtranspose_    =   QHtranspose.data();
}

template <typename T>
void FLIPACK<T>::get_HQHtranspose(double*& HQHtranspose_) {
    compute_HQHtranspose();
    HQHtranspose_   =   HQHtranspose.data();
}

template <typename T>
void FLIPACK<T>::get_Psi(double*& Psi_) {
    compute_Psi();
    Psi_    =   Psi.data();
}

template <typename T>
void FLIPACK<T>::get_Phi(double*& Phi_) {
    compute_Phi();
    Phi_    =   Phi.data();
}

template <typename T>
void FLIPACK<T>::get_Xi(double*& Xi_) {
    compute_Xi();
    Xi_ =   Xi.data();
}

template <typename T>
void FLIPACK<T>::get_Beta(double*& Beta_) {
    compute_Beta();
    Beta_   =   Beta.data();
}

template <typename T>
void FLIPACK<T>::get_Solution(double*& Solution_) {
    compute_Solution();
    Solution_ = Solution.data();
}




#endif //__FLIPACK_hpp__
