/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*!\file read_X_R_Measurements.cpp
 */

#include"read_X_R_Measurements.hpp"


void read_By_Line(const string& s, unsigned long row, double* data, unsigned numRows, unsigned numCols);


void read_X_R_Measurements (const string& filename, unsigned long N, unsigned short p, unsigned m, unsigned nMeasurementSets, double* X, double* R, double* measurements) {
    ifstream fin;
	fin.open(filename.c_str());
    
	if (!fin.good()){
		cerr << "Failed to open file " << filename << endl;
		throw runtime_error("Failed to open file!");
	}
    
    string line;
    unsigned long row = 0;
    while(getline(fin, line)){
        line.erase(remove(line.begin(), line.end(), ' '),
                   line.end());
        if (row < N) {
            read_By_Line(line, row, X, N, p);
        }else if (row < N+m) {
            read_By_Line(line, row-N, R, m, m);
        }else {
            read_By_Line(line, row-N-m, measurements, m, nMeasurementSets);
        }
        row++;
    }
    fin.close();
}

void read_By_Line(const string& s, unsigned long row, double* data, unsigned numRows, unsigned numCols) {
    if (!s.empty()) {
        unsigned k = 0;
        const char* start_pt = NULL;
        for (unsigned i = 0; i < s.length();)
            if ( s[i]==',' || s[i]=='(' ) {
                start_pt=&s[++i];
                while(s[i]!=',' && s[i]!=')') {
                    i++;
                }
                if(start_pt!=&s[i]) {
                    data[k*numRows+row]=(double)atof(start_pt);
                }
                k++;
            }
            else {
                i++;
            }
        if(k!=numCols)
            throw runtime_error("Number of measurement is not consistent with input");
    }
}
