/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	FLIPACK_get_matrix_through_routine_mykernel.cpp
 Input type: Through matrix generating routine;
 Kernel type: kernel defined by user.
*/
#include"environment.hpp"
#include"FLIPACK_Header.hpp"
#include"BBFMM2D.hpp"



using namespace std;
using namespace Eigen;

/*!  Get the location of the unknowns from the user;*/
void get_Location(unsigned long& N, vector<Point>& location){
	N           =	5000;
	VectorXd tmp1	=	VectorXd::Random(N);
	VectorXd tmp2	=	VectorXd::Random(N);
    for (unsigned long i = 0; i < N; i++) {
        Point newPoint;
        newPoint.x =   tmp1[i];
        newPoint.y =   tmp2[i];
        location.push_back(newPoint);
    }
}

/*! Get the measurement operator from the user; */
void get_Measurement_Operator(const unsigned long N, unsigned& m, unsigned& nMeasurementSets, MatrixXd& Htranspose, MatrixXd& measurements, MatrixXd& R){
	m               =	10;                                     //	Number of measurements;
    nMeasurementSets=   5;                                      //  Number of measurement sets;
	Htranspose      =	MatrixXd::Random(N,m);                  //  Transpose of the measurement operator;
    measurements    =   MatrixXd::Random(m,nMeasurementSets);   //  Set of measurements;
    R               =   MatrixXd::Identity(m,m);                //  Covariance of the measurements;
}

/*!  Get the structure of the mean;*/
void get_X(unsigned const N, unsigned short& p, MatrixXd& X){
    p   =   6;
    X   =   MatrixXd::Random(N,p);
}

/*!  Get the number of Chebyshev nodes in one direction;*/
void get_nChebNode(unsigned short& nChebNode){
    nChebNode    =   8;
}


class myKernel: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
        double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
        return 1.0 + rSquare;
    }
};

int main(){
    /**********************************************************/
    /*                                                        */
    /*              Initializing the problem                  */
    /*                                                        */
    /**********************************************************/
    
    cout << endl << "INITIALIZING THE PROBLEM..." << endl;

    clock_t start   =   clock();

    /*******    Getting the configuration of the grid   *******/

	unsigned long N;            //  Number of unknowns;
    vector<Point> location;     //  Location of the unknowns;
    
    get_Location(N,location);
    
    
    cout << endl << "Number of unknowns is: " << N << endl;
    
    /** Getting measurement related information measurements **/
    
    unsigned m;                 //  Number of measurements;
    unsigned nMeasurementSets;  //  Number of measurement sets;
    MatrixXd Htranspose;        //  Transpose of the measurement operator;
    MatrixXd measurements;      //  Actual measurements;
    MatrixXd R;                 //  Covariance matrix;
    
    get_Measurement_Operator(N, m, nMeasurementSets, Htranspose, measurements, R);
    
    cout << endl << "Number of measurements is: " << m << endl;
    cout << endl << "Number of sets of measurements is: " << nMeasurementSets << endl;

    /***************    Getting the structure   ***************/
    
    MatrixXd X;                 //  Structure of the mean;
    unsigned short p;           //  Number of terms in the structure;
    
    get_X(N, p, X);

    cout << endl << "Number of terms in the structure is: " << p << endl;

    /***  Getting the number of Chebyshev nodes for the fmm  ***/
    
    unsigned short nChebNode;   //  Number of Chebyshev nodes( >= 3)
                                //  per dimension;
    
    get_nChebNode(nChebNode);

    cout << endl << "Number of Chebyshev nodes along one direction is: " << nChebNode << endl;

    clock_t end   =   clock();
    
    double timeInitialize  =   double(end-start)/double(CLOCKS_PER_SEC);
    
    cout << endl << "Time taken to initialize the problem is: " << timeInitialize << endl;
    
    /**********************************************************/
    /*                                                        */
    /*                 Fast Linear inversion                  */
    /*                                                        */
    /**********************************************************/
    
    cout << endl << "PERFORMING FAST LINEAR INVERSION..." << endl;
    
    start   =   clock();
    H2_2D_Tree Atree(nChebNode, Htranspose, location);// Build the fmm tree;
    /* Options of kernel:
     LOGARITHM:          kernel_Logarithm
     ONEOVERR2:          kernel_OneOverR2
     GAUSSIAN:           kernel_Gaussian
     QUADRIC:            kernel_Quadric
     INVERSEQUADRIC:     kernel_InverseQuadric
     THINPLATESPLINE:    kernel_ThinPlateSpline
     */
    
    FLIPACK<myKernel> A(Htranspose, X, measurements, R, &Atree);
    
    A.get_Solution();
    
    end   =   clock();
    
    /****     If you want to use more than one kernels    ****/
    
    /*FLIPACK<kernel_Logarithm> C(Htranspose, X, measurements, R, &Atree);
    
    C.get_Solution();*/
    
    double timeFastMethod =   double(end-start)/double(CLOCKS_PER_SEC);

    cout << endl << "Time taken for the fast method is: " << timeFastMethod << endl;

    /**********************************************************/
    /*                                                        */
    /*              Conventional Linear inversion             */
    /*                                                        */
    /**********************************************************/

    cout << endl << "PERFORMING CONVENTIONAL LINEAR INVERSION..." << endl;
    
    start   =   clock();
    
    MatrixXd Q;
    
    myKernel B; // Make sure the type of B here
                // corresponds to the kernel used
                // to generate Q.
    
    B.kernel_2D(N, location, N, location, Q);
    
    MatrixXd temp(m+p,m+p);
    temp.block(0,0,m,m) =   Htranspose.transpose()*Q*Htranspose+R;
    temp.block(0,m,m,p) =   Htranspose.transpose()*X;
    temp.block(m,0,p,m) =   X.transpose()*Htranspose;
    temp.block(m,m,p,p) =   MatrixXd::Zero(p,p);
    

    MatrixXd temprhs(m+p,nMeasurementSets);
    temprhs.block(0,0,m,nMeasurementSets)   =   measurements;
    temprhs.block(m,0,p,nMeasurementSets)   =   MatrixXd::Zero(p,nMeasurementSets);
    MatrixXd tempSolution   =   temp.fullPivLu().solve(temprhs);

    MatrixXd finalSolution  =   X*tempSolution.block(m,0,p,nMeasurementSets)+Q*Htranspose*tempSolution.block(0,0,m,nMeasurementSets);

    end   =   clock();
    
    double  timeExactMethod    =   double(end-start)/double(CLOCKS_PER_SEC);
    
    cout << endl << "Time taken for the exact method is: " << timeExactMethod << endl;
    
    cout << endl << "Relative difference between the fast and conventional solution is: " << (finalSolution-A.Solution).norm()/finalSolution.norm() << endl << endl;
}