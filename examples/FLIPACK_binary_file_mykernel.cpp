/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	FLIPACK_binary_file_mykernel.cpp
 Input type: Input binary file;
 Kernel type: kernel defined by user..

*/

#include"environment.hpp"
#include"FLIPACK_Header.hpp"
#include"BBFMM2D.hpp"



using namespace std;
using namespace Eigen;

/*! Deifne user's own kernel */
class myKernel: public kernel_Base {
public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
        double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
        return exp(-rSquare);
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
    /*******            Getting metadata          *******/
    
    unsigned long N;            //  Number of unknowns;
    unsigned m;                 //  Number of measurements;
    unsigned nMeasurementSets;  //  Number of measurement sets;
    unsigned short p;           //  Number of terms in the structure;
    
    string filenameMetadata =   "../input/metadata.txt";
    
    read_Metadata (filenameMetadata, N, p,m, nMeasurementSets);
    
    cout << endl << "Number of unknowns is: "               << N << endl;
    cout << endl << "Number of measurements is: "           << m << endl;
    cout << endl << "Number of sets of measurements is: "   << nMeasurementSets << endl;
    cout << endl << "Number of terms in the structure is: " << p << endl;


    /*******    Getting location and Htranspose   *******/

    vector<Point> location;     //  Location of the unknowns;
    MatrixXd Htranspose;        //  Transpose of the measurement operator;
    
    string filenameLocation = "../input/test_Location.bin";
    string filenameH   = "../input/test_H.bin";      
    
    read_Location_Htranpose_binary(filenameLocation, N, location,  filenameH, m, Htranspose);


    /*******     Getting X, R and measurements     *******/
    
    MatrixXd measurements;      //  Actual measurements;
    MatrixXd R;                 //  Covariance matrix;
    MatrixXd X;                 //  Structure of the mean;
    
    string filenameX    =   "../input/test_X.bin";
    string filenameR    =   "../input/test_R.bin";
    string filenameMeasurement  =   "../input/test_Measurement.bin";
    
    read_X_R_Measurements_Binary (filenameX, N, p, X, filenameR, m, R, filenameMeasurement, nMeasurementSets, measurements);
    
    
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
    
    FLIPACK<myKernel> A(location, Htranspose, X, measurements, R, nChebNode, &Atree);
    
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