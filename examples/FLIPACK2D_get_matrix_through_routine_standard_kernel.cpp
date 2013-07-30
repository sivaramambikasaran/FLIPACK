//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//
//	FLIPACK2D_get_matrix_through_routine_standard_kernel.cpp
//
#include"environment.hpp"
#include"FLIPACK2D.hpp"
#include"read_X_R_measurements.hpp"


using namespace std;
using namespace Eigen;

//  Function gets the location of the unknowns from the user;
void get_Location(unsigned long& N, vector<Point>& location){
	N           =	5000;
	VectorXd tmp1	=	VectorXd::Random(N);
	VectorXd tmp2	=	VectorXd::Random(N);
    for (unsigned long i = 0; i < N; i++) {
        Point new_Point;
        new_Point.x =   tmp1[i];
        new_Point.y =   tmp2[i];
        location.push_back(new_Point);
    }
}

//  Measurement operator from the user;
void get_Measurement_operator(const unsigned long N, unsigned& m, unsigned& nmeasurementsets, MatrixXd& Htranspose, MatrixXd& measurements, MatrixXd& R){
	m               =	10;                                     //	Number of measurements;
    nmeasurementsets=   5;                                      //  Number of measurement sets;
	Htranspose      =	MatrixXd::Random(N,m);                  //  Transpose of the measurement operator;
    measurements    =   MatrixXd::Random(m,nmeasurementsets);   //  Set of measurements;
    R               =   MatrixXd::Identity(m,m);                //  Covariance of the measurements;
}

//  Get the structure of the mean;
void get_X(unsigned const N, unsigned short& p, MatrixXd& X){
    p   =   6;
    X   =   MatrixXd::Random(N,p);
}

//  Get the number of Chebyshev nodes in one direction;
void get_nchebnode(unsigned short& nchebnode){
    nchebnode    =   8;
}


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
    vector<Point> location;       //  Location of the unknowns;
    
    get_Location(N,location);
    
    
    cout << endl << "Number of unknowns is: " << N << endl;
    
    /** Getting measurement related information measurements **/
    
    unsigned m;                 //  Number of measurements;
    unsigned nmeasurementsets;  //  Number of measurement sets;
    MatrixXd Htranspose;        //  Transpose of the measurement operator;
    MatrixXd measurements;      //  Actual measurements;
    MatrixXd R;                 //  Covariance matrix;
    
    get_Measurement_operator(N, m, nmeasurementsets, Htranspose, measurements, R);
    
    cout << endl << "Number of measurements is: " << m << endl;
    cout << endl << "Number of sets of measurements is: " << nmeasurementsets << endl;

    /***************    Getting the structure   ***************/
    
    MatrixXd X;                 //  Structure of the mean;
    unsigned short p;           //  Number of terms in the structure;
    
    get_X(N, p, X);

    cout << endl << "Number of terms in the structure is: " << p << endl;

    /***  Getting the number of Chebyshev nodes for the fmm  ***/
    
    unsigned short nchebnode;   //  Number of Chebyshev nodes( >= 3)
                                //  per dimension;
    
    get_nchebnode(nchebnode);

    cout << endl << "Number of Chebyshev nodes along one direction is: " << nchebnode << endl;

    clock_t end   =   clock();
    
    double time_Initialize  =   double(end-start)/double(CLOCKS_PER_SEC);
    
    cout << endl << "Time taken to initialize the problem is: " << time_Initialize << endl;
    
    /**********************************************************/
    /*                                                        */
    /*                 Fast Linear inversion                  */
    /*                                                        */
    /**********************************************************/
    
    cout << endl << "PERFORMING FAST LINEAR INVERSION..." << endl;
    
    start   =   clock();
    
    H2_2D_tree Atree(nchebnode, Htranspose, location);// Build the fmm tree;
    
    /* Options of kernel:
     LOGARITHM:          kernel_Logarithm
     ONEOVERR2:          kernel_OneOverR2
     GAUSSIAN:           kernel_Gaussian
     QUADRIC:            kernel_Quadric
     INVERSEQUADRIC:     kernel_InverseQuadric
     THINPLATESPLINE:    kernel_ThinPlateSpline
     */
    
    FLIPACK2D<kernel_Gaussian> A(location, Htranspose, X, measurements, R, nchebnode, &Atree);
    
    A.get_Solution();
        
    end   =   clock();
    
    /****     If you want to use more than one kernels    ****/
    
    /*FLIPACK2D<kernel_Logarithm> C(location, Htranspose, X, measurements, R, nchebnode, &Atree);
     
     C.get_Solution();*/
    
    double time_Fast_method =   double(end-start)/double(CLOCKS_PER_SEC);

    cout << endl << "Time taken for the fast method is: " << time_Fast_method << endl;

    /**********************************************************/
    /*                                                        */
    /*              Conventional Linear inversion             */
    /*                                                        */
    /**********************************************************/

    cout << endl << "PERFORMING CONVENTIONAL LINEAR INVERSION..." << endl;
    
    start   =   clock();
    
    MatrixXd Q;
    
    kernel_Gaussian B; // Make sure the type of B here
                       // corresponds to the kernel used
                       // to generate Q.
    
    B.kernel2D(N, location, N, location, Q);
    
    MatrixXd temp(m+p,m+p);
    temp.block(0,0,m,m) =   Htranspose.transpose()*Q*Htranspose+R;
    temp.block(0,m,m,p) =   Htranspose.transpose()*X;
    temp.block(m,0,p,m) =   X.transpose()*Htranspose;
    temp.block(m,m,p,p) =   MatrixXd::Zero(p,p);
    

    MatrixXd temprhs(m+p,nmeasurementsets);
    temprhs.block(0,0,m,nmeasurementsets)   =   measurements;
    temprhs.block(m,0,p,nmeasurementsets)   =   MatrixXd::Zero(p,nmeasurementsets);
    MatrixXd tempsolution   =   temp.fullPivLu().solve(temprhs);

    MatrixXd finalsolution  =   X*tempsolution.block(m,0,p,nmeasurementsets)+Q*Htranspose*tempsolution.block(0,0,m,nmeasurementsets);

    end   =   clock();
    
    double  time_Exact_method    =   double(end-start)/double(CLOCKS_PER_SEC);
    
    cout << endl << "Time taken for the exact method is: " << time_Exact_method << endl;
    
    cout << endl << "Relative difference between the fast and conventional solution is: " << (finalsolution-A.Solution).norm()/finalsolution.norm() << endl << endl;
}