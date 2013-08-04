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
    FLIPACK(vector<Point>& location, MatrixXd& Htranspose, MatrixXd& X, MatrixXd& measurements, MatrixXd& R, unsigned short nChebNode, H2_2D_Tree*Atree);
    /*! This function obtains the cross covariance.*/
    void get_QHtranspose();
    /*! This function obtains the measurement operator corrected covariance. */
    void get_HQHtranspose();
    /*! This function obtains the measurement corrected covariance. */
    void get_Psi();
    /*! This functiin obtains the measurement operator corrected structure. */
    void get_Phi();
    /*! This function obtains the correction term.*/
    void get_Xi();
    /*! This function obtains unknown drift coefficients. */
    void get_Beta();
    /*! This function obtains the unknown.
        Running this function will get 
        all members defined in this class.
     */
    void get_Solution();
    void get_Posterior_Variance();
    
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
    vector<Point> location;
    
    MatrixXd mainMatrix;
    
    
    bool computedQHtranspose, computedHQHtranspose, computedPsi, computedPhi, computedSolution, computedXi, computedBeta, computedVDiag, computedMainMatrix, computedIntermediateSolution;
    
    unsigned short nChebNode;   /*!<  Number of Chebyshev nodes per dimension */
    unsigned long N;            /*!<  Number of unknowns */
    unsigned m;                 /*!<  Number of measurements */
    unsigned nMeasurementSets;  /*!<  Number of sets of measruements*/
    unsigned p;                 /*!<  Number of columns of X*/
    
    void get_Main_Matrix();
    
    void get_Intermediate_Solution();
    
    FullPivLU<MatrixXd> lu;
    
    MatrixXd intermediateSolution;
};



template <typename T>
FLIPACK<T>::FLIPACK(vector<Point>& location, MatrixXd& Htranspose, MatrixXd& X, MatrixXd& measurements, MatrixXd& R, unsigned short nChebNode, H2_2D_Tree *Atree){
    this->location          =   location;
    this->Htranspose        =   Htranspose;
    this->X                 =   X;
    this->measurements      =   measurements;
    this->R                 =   R;
    this->nChebNode         =   nChebNode;
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
    
    N                       =   Htranspose.rows();
    m                       =   Htranspose.cols();
    p                       =   X.cols();
    nMeasurementSets        =   measurements.cols();
    
    mainMatrix              =   MatrixXd(m+p,m+p);
    QHtranspose             =   MatrixXd(N,m);
}


template <typename T>
void FLIPACK<T>::get_QHtranspose(){
    if (computedQHtranspose==false) {
        cout << endl << "Performing FMM to obtain QHtranspose..." << endl;
        QHtranspose            =   MatrixXd(N,m);
        T A;
        A.calculate_Potential(*Atree,QHtranspose);
        computedQHtranspose    =   true;
        cout << endl << "Obtained QHtranspose" << endl;
    }
}


template <typename T>
void FLIPACK<T>::get_HQHtranspose(){
    if (computedHQHtranspose==false) {
        get_QHtranspose();
        HQHtranspose   =   Htranspose.transpose()*QHtranspose;
        computedHQHtranspose  =   true;
        cout << "Obtained HQHtranspose" << endl;
    }
}

template <typename T>
void FLIPACK<T>::get_Psi(){
    if (computedPsi==false) {
        get_HQHtranspose();
        Psi =   HQHtranspose + R;
        computedPsi            =   true;
        cout << "Obtained Psi" << endl;
    }
}

template <typename T>
void FLIPACK<T>::get_Phi(){
    if (computedPhi==false) {
        Phi =   Htranspose.transpose()*X;
        computedPhi            =   true;
        cout << "Obtained Phi" << endl;
    }
}

template <typename T>
void FLIPACK<T>::get_Main_Matrix(){
    if (computedMainMatrix==false) {
        get_Psi();
        mainMatrix.block(0,0,m,m)  =   Psi;
        get_Phi();
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
void FLIPACK<T>::get_Xi(){
    if (computedXi==false) {
        get_Intermediate_Solution();
        Xi              =   intermediateSolution.block(0,0,m,nMeasurementSets);
        computedXi     =   true;
        cout << "Obtained Xi" << endl;
    }
}

template <typename T>
void FLIPACK<T>::get_Beta(){
    if (computedXi==false) {
        get_Intermediate_Solution();
        Beta            =   intermediateSolution.block(m,0,p,nMeasurementSets);
        computedBeta   =   true;
        cout << "Obtained Beta" << endl;
    }
}

template <typename T>
void FLIPACK<T>::get_Solution(){
    if (computedSolution==false) {
        get_QHtranspose();
        get_Beta();
        get_Xi();
        Solution            =   X*Beta +   QHtranspose*Xi;
        computedSolution   =   true;
        cout << "Obtained Solution" << endl;
    }
}



#endif //__FLIPACK_hpp__
