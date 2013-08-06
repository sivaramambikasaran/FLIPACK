/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran, Ruoxi Wang
 *  \version 3.1
 */
/*! \file	read_metadata.hpp
Read Number of unknowns, number of terms in structure,number of measurements,number of measurement sets from file
*/
#ifndef __read_X_R_Measurements_binary_hpp__
#define __read_X_R_Measurements_binary_hpp__

#include"environment.hpp"


using namespace Eigen;
using namespace std;


void read_Metadata (const string& filenameX, unsigned long& N, unsigned short& p,unsigned& m, unsigned& nMeasurementSets);


#endif //(__read_metadata_hpp__)
