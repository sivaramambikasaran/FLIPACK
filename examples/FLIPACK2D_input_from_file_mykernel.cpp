//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//
//	FLIPACK2D_input_from_file_standard_kernel.cpp
//
#include"iostream"
#include"ctime"
#include"cmath"
#include"Eigen/Dense"
#include"FLIPACK2D.hpp"
#include"H2_2D_tree.hpp"
#include"kernel_types.hpp"
#include"read_X_R_measurements.hpp"
#include"read_location_H.hpp"


using namespace std;
using namespace Eigen;

class mykernel: public kernel_base {
public:
    virtual double kernel_func(double R_square){
        //define your own kernel here
        return 1.0 + R_square;
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

    /*******    Getting location and Htranspose   *******/

	unsigned long N;            //  Number of unknowns;
    VectorXd location[2];       //  Location of the unknowns;
    unsigned m;                 //  Number of measurements;
    MatrixXd Htranspose;        //  Transpose of the measurement operator;
    
    string filename_location_Htranpose = "./../input/test_location_H.txt";
    
    read_Location_and_Measurement_operator (filename_location_Htranpose.c_str(),  N,location,  m, Htranspose);

    
    cout << endl << "Number of unknowns is: "     << N << endl;
    cout << endl << "Number of measurements is: " << m << endl;

    /*******     Getting X, R and measurements     *******/
    
    unsigned nmeasurementsets;  //  Number of measurement sets;
    MatrixXd measurements;      //  Actual measurements;
    MatrixXd R;                 //  Covariance matrix;
    MatrixXd X;                 //  Structure of the mean;
    unsigned short p;           //  Number of terms in the structure;
    
    string filename_X_R_measurements = "./../input/test_X_R_measurements.txt";
    
    read_X_R_measurements(filename_X_R_measurements.c_str(),N,p,m,nmeasurementsets,X,R,measurements);
    
    cout << endl << "Number of sets of measurements is: " << nmeasurementsets << endl;
    cout << endl << "Number of terms in the structure is: " << p << endl;    
    
    /***  Getting the number of Chebyshev nodes for the fmm  ***/
    
    unsigned short nchebnode = 8;   //  Number of Chebyshev nodes( >= 3)
                                    //  per dimension;
    
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
    
    /* Options of kernel:
     LOGARITHM:          kernel_Logarithm
     ONEOVERR2:          kernel_OneOverR2
     GAUSSIAN:           kernel_Gaussian
     QUADRIC:            kernel_Quadric
     INVERSEQUADRIC:     kernel_InverseQuadric
     THINPLATESPLINE:    kernel_ThinPlateSpline
     */
    
    FLIPACK2D<mykernel> A(location, Htranspose, X, measurements, R, nchebnode);
    
    A.get_Solution();
        
    end   =   clock();
    
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
    
    mykernel B; // Make sure the type of B here
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