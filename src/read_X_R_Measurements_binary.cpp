/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file read_X_R_Measurements_binary.cpp
 */

#include"read_X_R_Measurements_binary.hpp"

using namespace std;

void read_X_R_Measurements (const string& filenameX, unsigned long N, unsigned short p,MatrixXd& X,const string& filenameR, unsigned m, MatrixXd& R,const string& filenameMeasurementSets, unsigned nMeasurementSets, MatrixXd& measurements) {
    ifstream fin;
    /* Read X */
	fin.open(filenameX.c_str(),ios::binary);
    
	if (!fin.good()){
		cerr << "Failed to open file " << filenameX << endl;
		throw runtime_error("Failed to open file!");
	}
    X.resize(N,p);
    fin.read((char*) X.data(), N*p*sizeof(double));
    fin.close();
    
    /* Read R */
    fin.open(filenameR.c_str(),ios::binary);
    
	if (!fin.good()){
		cerr << "Failed to open file " << filenameR << endl;
		throw runtime_error("Failed to open file!");
	}
    R.resize(m,m);
    fin.read((char*) R.data(), N*p*sizeof(double));
    fin.close();
    
    /* Read Measurements */
    fin.open(filenameMeasurementSets.c_str(),ios::binary);
    
	if (!fin.good()){
		cerr << "Failed to open file " << filenameMeasurementSets << endl;
		throw runtime_error("Failed to open file!");
	}
    measurements.resize(m,m);
    fin.read((char*) measurements.data(), m*nMeasurementSets*sizeof(double));
    fin.close();
}

