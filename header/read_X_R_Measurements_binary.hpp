/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	read_X_R_Measurements_binary.hpp
Read X, R and Measurements from binary file.
*/
#ifndef __read_X_R_Measurements_binary_hpp__
#define __read_X_R_Measurements_binary_hpp__

#include"environment.hpp"


using namespace Eigen;
using namespace std;


void read_X_R_Measurements_Binary (const string& filenameX, unsigned long N, unsigned short p,MatrixXd& X,const string& filenameR, unsigned m, MatrixXd& R,const string& filenameMeasurement, unsigned nMeasurementSets, MatrixXd& measurements);


#endif //(__read_X_R_Measurements_binary_hpp__)
