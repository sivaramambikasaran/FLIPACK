//	This Source Code Form is subject to the terms of the Mozilla Public
//	License, v. 2.0. If a copy of the MPL was not distributed with this
//	file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//	<author>Sivaram Ambikasaran, Ruoxi Wang</author>
//	
//	read_X_R_measurements.cpp
//

#include"read_X_R_measurements.hpp"



void read_X_R_measurements (const string& filename, unsigned long &N, unsigned short& p, unsigned& m, unsigned& nmeasurementsets, MatrixXd& X, MatrixXd& R, MatrixXd& measurements) {
    ifstream fin;
	fin.open(filename.c_str());
	
	if (!fin.good()){
		cerr << "Failed to open file " << filename << endl;
		throw runtime_error("Failed to open file!");
	}
    // read first line in the file : N, p, m, measurementsets
    string line;
    getline(fin,line);
    line.erase(remove(line.begin(), line.end(), ' '),
               line.end());
    stringstream ss;
    ss << line;
    char comma;
    ss >> N >> comma >> p >> comma >> m >> comma >> nmeasurementsets;
    X               =   MatrixXd::Zero(N,p);
    R               =   MatrixXd::Zero(m,m);
    measurements    =   MatrixXd::Zero(m,nmeasurementsets);
    unsigned long row = 0;
    while(getline(fin, line)){
        line.erase(remove(line.begin(), line.end(), ' '),
                   line.end());
        if (row < N) {
            read_matrix_by_line(line, row, X, p);
        }else if (row < N+m) {
            read_matrix_by_line(line, row-N, R, m);
        }else {
            read_matrix_by_line(line, row-N-m, measurements, nmeasurementsets);
        }
        row++;
    }
    fin.close();
}

void read_matrix_by_line(const string& s, unsigned long row, MatrixXd& M, unsigned m) {
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
                    M(row,k)=(double)atof(start_pt);
                }
                k++;
            }
            else {
                i++;
            }
        if(k!=m)
            throw runtime_error("Number of measurement is not consistent with input");
    }
}


