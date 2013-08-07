/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	read_X_R_Measurements.hpp
Read X, R and Measurements from text file.
*/
#ifndef __read_X_R_Measurements_hpp__
#define __read_X_R_Measurements_hpp__

#include"environment.hpp"


using namespace Eigen;
using namespace std;


void read_X_R_Measurements (const string& filename, unsigned long &N, unsigned short& p, unsigned& m, unsigned& nMeasurementSets, MatrixXd& X, MatrixXd& R, MatrixXd& measurements);

void read_Matrix_By_Line(const string& s, unsigned long row, MatrixXd& M, unsigned m);

#endif //(__read_X_R_Measurements_hpp__)
