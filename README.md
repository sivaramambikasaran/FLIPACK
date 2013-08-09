
This is the first public release of the FLIPACK library.

Date: July 26th, 2013

**Version 3.1 - First external release.**

%% Copyleft 2013: Sivaram Ambikasaran, Ruoxi Wang, Peter Kitanidis and Eric Darve  
%% Developed by Sivaram Ambikasaran, Ruoxi Wang  
%% Contact: <siva.1985@gmail.com>(Sivaram) , <ruoxi@stanford.edu> (Ruoxi)  
%%   
%% This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license.  
%% The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not %% distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/>.  


###DIRECTORIES AND FILES:

	./examples/		:	Example input C++ codes; Needed to read input from user or from input file.
	./src/			:	Source code in C++
	./header/		:	Relevant header files
	./exec/			:	Executables for BBFMM2D
	./input/		:	The input file.
	./README.md		:	This file
	./License.md	:	License file
	./Makefile		:	Makefile  

###SETTING THINGS UP:

1. To run this package, you need to have **Eigen** and **BBFMM2D**.
    * Set Eigen:
      
	    1). Download Eigen from here: <http://eigen.tuxfamily.org/index.php?title=Main_Page> 
	     
	    2).  Create a directory named `Codes/` inside the main Eigen folder and copy the directory  `FLIPACK/` into the directory `Codes/`. 
	     
	    3).  Open the Makefile, which is inside the folder FLIPACK. Ensure that you have included the path to Eigen in the line containing `CFLAGS`. For instance, in the above setting, the path `"-I ./../../"` should be included in the Makefile.  
    * Set BBFMM2D:  
    
	    1). Download BBFMM2D from here <http://sivaramambikasaran.github.io/BBFMM2D/>  
	    
	    2).  Copy directory `BBFMM2D/` inside of the directory `Codes/`.

2. Once you have this set up, you should be able to run the code. Check whether the code runs by performing the following action. Go to the directory where Makefile is in, then key in the following three commands in the terminal:

		make get_matrix_through_routine_standard_kernel
		cd exec/
		./FLIPACK_get_matrix_through_routine_standard_kernel

The code should now run.

	
###CHANGING THE INPUTS:

The files you have control over are the files inside the directory `./examples/`, read through the files both must be self explanatory for the most part.

1. If you want to generate matrix through your own routine, and use the standard kernels:  

    Go to `/examples/`, open `"FLIPACK_get_matrix_through_routine_standard_kernel.cpp"`.  
    * To generate matrix through routine:   
      Change functions `get_Location()`, `get_Measurement_Operator()`, `get_X()` and `get_nChebNode()`.
    * To use standard kernels:   
      Select the kernel type in `main()`.  
      Options of kernels:
        
  			LOGARITHM:          kernel_Logarithm  
  			ONEOVERR2:          kernel_OneOverR2  
  			GAUSSIAN:           kernel_Gaussian  
  			QUADRIC:            kernel_Quadric  
     		INVERSEQUADRIC:     kernel_InverseQuadric  
     		THINPLATESPLINE:    kernel_ThinPlateSpline  

	
2. If you want to read matrix from file, and use standard kernels:

    Go to the folder `/input/`, put your input file inside of this folder. 
     
    Go to the folder `/examples/`, open `"FLIPACK_input_from_file_standard_kernel.cpp"`.
    * To change input filename:  
    
      change the two lines in `main()`:  
      
      		string filename_location_Htranpose = "./../input/test_Location_H.txt";  
      		string filename_X_R_Measurements = "./../input/test_X_R_Measurements.txt";
    * To use standard kernels:  
    
      The same step as described in 1.


3. If you want to generate matrix through your own routine, and use your own kernel:

    Go to `/examples/`, open `"FLIPACK_get_matrix_through_routine_mykernel.cpp"`.
    * To define your own kernel:  
      Modify `class myKernel`. 
    * To generate your matrix:  
      The same step as described in 1.

4. If you want to read matrix from file, and use your own kernel:

    Go to `/examples/`, open `"FLIPACK_input_from_file_mykernel.cpp"`.
    * To define your own kernel:  
      Modify `class myKernel`. 
    * To change input filename:  
      The same step as described in 2.  
      
5. If you want to read matrix from binary file, and use standard kernel:
 
    Go to `/examples`, open `"FLIPACK_binary_file_standard_kernel.cpp"`.  
    * To change input filename:  
      change the following lines in `main()`:  
      
      		string filenameMetadata 		=   "../input/metadata.txt";  
      		string filenameLocation 		= 	"../input/test_Location.bin";  
      		string filenameH        		= 	"../input/test_H.bin";  
      		string filenameX        		=   "../input/test_X.bin";  
      		string filenameR        		=   "../input/test_R.bin";  
      		string filenameMeasurement  	=   "../input/test_Measurement.bin";  
    * To use standard kernel:  
    
     The same step as described in 1.
     
6. If you want to read matrix from binary file, and use your own kernel:  

    Go to `/examples/`, open `"FLIPACK_binary_file_mykernel.cpp"`.  
    * To change the input filename:  
      The same step as described in 5.  
    * To define your own kernel:  
      Modify `class myKernel`.      


###INPUT FILES

Go to `/input/`, you should put your own input files in the input folder.  
**Donotes:**  

	N: 					Number of unknowns  
	m: 					Number of measurements  
	p: 					Number of terms in the structure  
	nMeasurementSet:	Number of measurement sets  
####TEXT FILES

The file format is described as follows:

* For input file of location and Htranspose:  
   	
   The first row should be like this 
      
		Number of unknowns, Number of sets of measurements
		
   For example:

    	20000,10

   For the rest of the rows, it should start with locations, followed by a row in Htranspose matrix(elements should be separated using ',') If some element is 0, you can leave it as empty instead of 0. If all the elements in a row is 0, nothing need to be typed after the location.(spaces are allowed)

   The row should look like this: 
     
    	(location[0],location[1]) (elem1,elem2,elem3,elem4,…,elemn)

   For example:

		(-0.999984,-0.676221) (0.480685,0.869803,-0.188232,0.548587,-0.771039,0.73709,0.126494)  
		(0.869386,-0.5408)  
		(0.0655345,0.891162) (-0.193033,,-0.0287383,,-0.520512,0.33891,)  
		(0.342299,-0.246828) (0.0732668,,,,,,0.0951028)  
		(-0.984604,-0.44417) (,0.782447,-0.867924,0.485731,-0.729282,-0.481031,0.541473) 
* For input file of X, R and measurements:  

    The first row should be like this:  
    
    	N, p, m, nMeasurementSet

    For example:
 
		20000,6,10,5

    The rest of rows should be rows of X, and then rows of R, and rows of measurements     (in order).  
    Each row should look like this:    
    `(elem1,elem2,elem3,…,elemn)`  
    If some element is 0, you can leave it as empty instead of 0.

    For example:  
    
		(0.819666,-0.996573,-0.21957,-0.532382,0.491672,0.730317)
		(,0.332605,0.0938555,0.0442606,-0.63103,0.322981)
		(0.126619,,,0.262597,-0.470423,0.432109)
		(-0.00923308,0.0931162,-0.57049,-0.112485,0.283007,)

    The whole file should look like:

   		N, p, m, nMeasurementSet	 
    	Rows of X 
    	Rows of R  
    	Rows of measurements 

####BINARY FILES

You should have seperate binary files for each matrix, and a text file containing metadata:  

1. A text file for meta data: 
 
   The file format is:  
   
   		N, p, m, nMeasurement 
   For example:  
   
   		3245,6,288,6  
2. A binary file for Location:

	Elements are stored this way(row-wise):
		
		loc0.x loc0.y  
		loc1.x loc1.y  
		…
3. Binary files for H matrix, X, R, Measurements respectively:  
   The elements of each matrix should be stored in their own binary file row-wise. 

###RUNNING THE CODE:  

Here we give an example:  
If you want to use `"FLIPACK_input_from_file_standard_kernel.cpp"`

1. As stated earlier to run the code, go to the appropriate directory and key in the following:

		make input_from_file_standard_kernel

2. Make sure you have changed or cleaned the .o files from previous compilation. To clean the irrelevant files, key in:

		make clean

3. To tar the file, key in:

		make tar

4. Read through the makefile for other options.

To run other .cpp files:  

1). `FLIPACK_get_matrix_through_routine_mykernel.cpp`  
   key in:  
   
   		make get_matrix_through_routine_mykernel
   
2). `FLIPACK_get_matrix_through_routine_standard_kernel.cpp`  
   key in:  
   
   		make get_matrix_through_routine_standard_kernel
   
3). `FLIPACK_input_from_file_mykernel.cpp`   
   key in:  
   
   		make input_from_file_mykernel
   	
4). `FLIPACK_binary_file_mykernel.cpp`  
   key in:  
   
        make binary_file_mykernel  
        
5). `FLIPACK_binary_file_standard_kernel.cpp`  
   key in:  
   
        make binary_file_standard_kernel
