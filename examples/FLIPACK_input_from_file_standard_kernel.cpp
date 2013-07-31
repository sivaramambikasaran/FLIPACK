//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//
//	FLIPACK_input_from_file_standard_kernel.cpp
//
#include"environment.hpp"
#include"FLIPACK.hpp"
#include"read_X_R_Measurements.hpp"


using namespace std;
using namespace Eigen;

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
    vector<Point> location;       //  Location of the unknowns;
    unsigned m;                 //  Number of measurements;
    MatrixXd Htranspose;        //  Transpose of the measurement operator;
    
    string filenameLocationHtranpose = "./../input/test_Location_H.txt";
    
    read_Location_And_Measurement_Operator (filenameLocationHtranpose.c_str(),  N,location,  m, Htranspose);

    
    cout << endl << "Number of unknowns is: "     << N << endl;
    cout << endl << "Number of measurements is: " << m << endl;

    /*******     Getting X, R and measurements     *******/
    
    unsigned nMeasurementSets;  //  Number of measurement sets;
    MatrixXd measurements;      //  Actual measurements;
    MatrixXd R;                 //  Covariance matrix;
    MatrixXd X;                 //  Structure of the mean;
    unsigned short p;           //  Number of terms in the structure;
    
    string filenameXRMeasurements = "./../input/test_X_R_Measurements.txt";
    
    read_X_R_Measurements(filenameXRMeasurements.c_str(),N,p,m,nMeasurementSets,X,R,measurements);
    
    cout << endl << "Number of sets of measurements is: " << nMeasurementSets << endl;
    cout << endl << "Number of terms in the structure is: " << p << endl;    
    
    /***  Getting the number of Chebyshev nodes for the fmm  ***/
    
    unsigned short nChebNode = 8;   //  Number of Chebyshev nodes( >= 3)
                                    //  per dimension;
    
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
    
    FLIPACK<kernel_Gaussian> A(location, Htranspose, X, measurements, R, nChebNode, &Atree);
    
    A.get_Solution();
        
    end   =   clock();
    
    /****     If you want to use more than one kernels    ****/
    
    /*FLIPACK<kernel_Logarithm> C(location, Htranspose, X, measurements, R, nChebNode, &Atree);
     
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
    
    kernel_Gaussian B; // Make sure the type of B here
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