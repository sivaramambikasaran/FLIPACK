//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	read_X_R_measurements.hpp
//
#ifndef __read_X_R_measurements_hpp__
#define __read_X_R_measurements_hpp__

#include"environment.hpp"


using namespace Eigen;
using namespace std;

// X: (N, p); R: (m, m); measurements: (m, measurementsets)
/* file format:
 N, p, m, measurementsets
 (Xelem1,Xelem2,....,Xelemp)
 ....
 (Relem1,Relem2,....,Relemm)
 ....
 (meas_elem1,meas_elem2,...,meas_elemmm)
 ....
 */
void read_X_R_measurements (const string& filename, unsigned long &N, unsigned short& p, unsigned& m, unsigned& nmeasurementsets, MatrixXd& X, MatrixXd& R, MatrixXd& measurements);

void read_matrix_by_line(const string& s, unsigned long row, MatrixXd& M, unsigned m);

#endif //(__read_X_R_measurements_hpp__)
