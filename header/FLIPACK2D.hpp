//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran</author>
//	
//	FLIPACK2D.hpp
//
#ifndef __FLIPACK2D_hpp__
#define __FLIPACK2D_hpp__

#include"H2_2D_tree.hpp"
#include"iostream"
#include"Eigen/Dense"
//#include"BBFMM2D/H2_2D"

using namespace Eigen;
using namespace std;


template <typename T>
class FLIPACK2D{
public:
    FLIPACK2D(VectorXd* location, MatrixXd& H_transpose, MatrixXd& X, MatrixXd& measurements, MatrixXd& R, unsigned short nchebnode);
    ~FLIPACK2D();
    void get_QH_transpose();
    void build_fmm_tree();
    void get_HQH_transpose();
    void get_Psi();
    void get_Phi();
    void get_Xi();
    void get_Beta();
    void get_Solution();
    void get_Posterior_Variance();
    
    MatrixXd QH_transpose, HQH_transpose, Psi, Phi, Solution, Xi, Beta, V_diag;
    H2_2D_tree *Atree;
    
private:
    MatrixXd H_transpose;
    MatrixXd X;
    MatrixXd measurements;
    MatrixXd R;
    VectorXd* location;
    
    MatrixXd Main_Matrix;
    
    
    bool computed_QH_transpose, computed_HQH_transpose, computed_Psi, computed_Phi, computed_Solution, computed_Xi, computed_Beta, computed_V_diag, computed_Main_Matrix, computed_Intermediate_Solution,built_tree;
    
    unsigned short nchebnode;
    
    unsigned long N;            //  Number of unknowns
    unsigned m;                 //  Number of measurements
    unsigned nmeasurementsets;  //  Number of sets of measruements
    unsigned p;                 //  Number of columns of X
    
    void get_Main_Matrix();
    
    void get_Intermediate_Solution();
    
    FullPivLU<MatrixXd> lu;
    
    MatrixXd intermediate_solution;
};



template <typename T>
FLIPACK2D<T>::FLIPACK2D(VectorXd* location, MatrixXd& H_transpose, MatrixXd& X, MatrixXd& measurements, MatrixXd& R, unsigned short nchebnode){
    this->location          =   location;
    this->H_transpose       =   H_transpose;
    this->X                 =   X;
    this->measurements      =   measurements;
    this->R                 =   R;
    this->nchebnode         =   nchebnode;
    
    computed_QH_transpose           =   false;
    computed_HQH_transpose          =   false;
    computed_Psi                    =   false;
    computed_Phi                    =   false;
    computed_Main_Matrix            =   false;
    computed_Intermediate_Solution  =   false;
    computed_Xi                     =   false;
    computed_Beta                   =   false;
    computed_V_diag                 =   false;
    computed_Solution               =   false;
    built_tree                      =   false;
    
    N                       =   H_transpose.rows();
    m                       =   H_transpose.cols();
    p                       =   X.cols();
    nmeasurementsets        =   measurements.cols();
    
    Main_Matrix             =   MatrixXd(m+p,m+p);
    QH_transpose            =   MatrixXd(N,m);
}

template <typename T>
FLIPACK2D<T>::~FLIPACK2D() {
    delete Atree;
}

template <typename T>
void FLIPACK2D<T>::get_QH_transpose(){
    if (computed_QH_transpose==false) {
        build_fmm_tree();
        cout << endl << "Performing FMM to obtain QH_transpose..." << endl;
        QH_transpose            =   MatrixXd(N,m);
        T A;
        A.calculatepotential(*Atree,QH_transpose);
        computed_QH_transpose   =   true;
        cout << endl << "Obtained QH_transpose" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::build_fmm_tree() {
    if (!built_tree) {
        Atree =   new H2_2D_tree(nchebnode, H_transpose, location);
        built_tree  =   true;
        cout << "fmm tree built" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::get_HQH_transpose(){
    if (computed_HQH_transpose==false) {
        get_QH_transpose();
        HQH_transpose   =   H_transpose.transpose()*QH_transpose;
        computed_HQH_transpose  =   true;
        cout << "Obtained HQH_transpose" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::get_Psi(){
    if (computed_Psi==false) {
        get_HQH_transpose();
        Psi =   HQH_transpose + R;
        computed_Psi            =   true;
        cout << "Obtained Psi" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::get_Phi(){
    if (computed_Phi==false) {
        Phi =   H_transpose.transpose()*X;
        computed_Phi            =   true;
        cout << "Obtained Phi" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::get_Main_Matrix(){
    if (computed_Main_Matrix==false) {
        get_Psi();
        Main_Matrix.block(0,0,m,m)  =   Psi;
        get_Phi();
        Main_Matrix.block(0,m,m,p)  =   Phi;
        Main_Matrix.block(m,0,p,m)  =   Phi.transpose();
        Main_Matrix.block(m,m,p,p)  =   MatrixXd::Zero(p,p);
        lu.compute(Main_Matrix);
        computed_Main_Matrix    =   true;
        cout << "Obtained Main matrix" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::get_Intermediate_Solution(){
    if (computed_Intermediate_Solution==false) {
        get_Main_Matrix();
        MatrixXd rhs                        =   MatrixXd::Zero(m+p,nmeasurementsets);
        rhs.block(0,0,m,nmeasurementsets)   =   measurements;
        intermediate_solution               =   lu.solve(rhs);
        computed_Intermediate_Solution      =   true;
        cout << "Obtained Intermediate solution" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::get_Xi(){
    if (computed_Xi==false) {
        get_Intermediate_Solution();
        Xi              =   intermediate_solution.block(0,0,m,nmeasurementsets);
        computed_Xi     =   true;
        cout << "Obtained Xi" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::get_Beta(){
    if (computed_Xi==false) {
        get_Intermediate_Solution();
        Beta            =   intermediate_solution.block(m,0,p,nmeasurementsets);
        computed_Beta   =   true;
        cout << "Obtained Beta" << endl;
    }
}

template <typename T>
void FLIPACK2D<T>::get_Solution(){
    if (computed_Solution==false) {
        get_QH_transpose();
        get_Beta();
        get_Xi();
        Solution        =   X*Beta +   QH_transpose*Xi;
        computed_Solution   =   true;
        cout << "Obtained Solution" << endl;
    }
}



#endif //__FLIPACK2D_hpp__
