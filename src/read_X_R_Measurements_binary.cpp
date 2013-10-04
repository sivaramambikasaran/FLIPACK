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

void read_X_R_Measurements_Binary (const string& filenameX, unsigned long N, unsigned short p,double*& X,const string& filenameR, unsigned m, double*& R,const string& filenameMeasurement, unsigned nMeasurementSets, double*& measurements) {
    ifstream fin;
    /* Read X */
	fin.open(filenameX.c_str(),ios::binary);
    
	if (!fin.good()){
		cerr << "Failed to open file " << filenameX << endl;
		throw runtime_error("Failed to open file!");
	}
    fin.read((char*) X, N*p*sizeof(double));
    fin.close();
    
    /* Read R */
    fin.open(filenameR.c_str(),ios::binary);
    
	if (!fin.good()){
		cerr << "Failed to open file " << filenameR << endl;
		throw runtime_error("Failed to open file!");
	}
    fin.read((char*) R, m*m*sizeof(double));
    fin.close();
    
    /* Read Measurements */
    fin.open(filenameMeasurement.c_str(),ios::binary);
    
	if (!fin.good()){
		cerr << "Failed to open file " << filenameMeasurement << endl;
		throw runtime_error("Failed to open file!");
	}
    fin.read((char*) measurements, m*nMeasurementSets*sizeof(double));
    fin.close();
}

