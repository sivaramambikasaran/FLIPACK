
# FLIPACK
This is the first public release of the FLIPACK library.

Date: July 26th, 2013

**Version 3.1 - First external release.**

%% Copyleft 2013: Sivaram Ambikasaran, Ruoxi Wang, Peter Kitanidis and Eric Darve  
%% Developed by Sivaram Ambikasaran, Ruoxi Wang  
%% Contact: <siva.1985@gmail.com>(Sivaram) , <ruoxi@stanford.edu> (Ruoxi)  
%%   
%% This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license.  
%% The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not %% distributed with this file, You can obtain one at <http://mozilla.org/MPL/2.0/>.  

###1. INTRODUCTION
FLIPACK (Fast Linear Inversion PACKage) is a library for fast linear inversion as described in <a href="http://link.springer.com/article/10.1007/s10596-013-9364-0">this article</a>. Stochastic linear inversion is the backbone of applied inverse methods. Stochastic inverse modeling deals with the estimation of functions from sparse data, which is a problem with a nonunique solution, with the objective to evaluate best estimates, measures of uncertainty, and sets of solutions that are consistent with the data. As finer resolutions become desirable, the computational requirements increase dramatically when using conventional algorithms. The FLIPACK reduces the computational cost from O(N^2) to O(N) by modeling the large dense convariances arising in these problems as a hierarchical matrix, more specifically as a  <img src="http://latex.codecogs.com/svg.latex? $\mathcal{H}^2$ " border="0"/> matrix. Matrix-vector products for hierarchical matrices are accelerated using the <a href="https://github.com/sivaramambikasaran/BBFMM2D">Black Box Fast Multipole Method</a> to accelerate these matrix-vector products.


###2. LINEAR MODEL
####2.1 Prior
Consider that s(x) is a function to be estimated, The basic model of the function to be estimated is taken as:   
<img src="http://latex.codecogs.com/svg.latex? $s(x) = \sum_{k=1}^p f_k(x)\beta_k + \epsilon(x)$ " border="0"/>   
The first term is the prior mean, where <img src="http://latex.codecogs.com/svg.latex? $f_k(x)$ " border="0"/>  are known functions, typically step functions, polynomials, and <img src="http://latex.codecogs.com/svg.latex? $\beta(k)$ " border="0"/>  are unknown coefficients where k = 1,2,...,p. The second term is a random function with zero mean and characterized through a covariance function. After discretization, s(x) is represented through an m by 1 vector s. The mean of s is <img src="http://latex.codecogs.com/svg.latex? $E[s] = X \beta$ " border="0"/>. where X is a known m x p matrix, and <img src="http://latex.codecogs.com/svg.latex? $\beta$ " border="0"/> are p unknown drift coefficients. The covariance of s is <img src="http://latex.codecogs.com/svg.latex? $E[(s-\mu)(s-\mu)^T] = Q$ " border="0"/>.  

####2.2 Measurement equation
The observation/measurement is related to the unknown by the linear relation  
<img src="http://latex.codecogs.com/svg.latex? $y = Hs+v$ " border="0"/>  
where H is n by m given matrix; v is a random vector of observation error, independent from s, with mean zero and covariance matrix R. Then, the prior statistics of y are:  

The mean:  
<img src="http://latex.codecogs.com/svg.latex? $\mu_y = E[Hs + v] = HE[s]+E[v]=HX\beta = \Phi \beta$ " border="0"/>   
The covariance:  
<img src="http://latex.codecogs.com/svg.latex? $\Psi = E[(H(s-X\beta)+v)(H(s-X\beta)+v)^T] = HQH^T$ " border="0"/>  
The y to s corss-covariance:  
<img src="http://latex.codecogs.com/svg.latex? $C_{ys} = E[(H(s-X\beta)+v)(s-X\beta)^T] = HQ$ " border="0"/>
####2.3 The ξ Form 
Introduce the n × 1 vector <img src="http://latex.codecogs.com/svg.latex? $\xi$ " border="0"/> defined through  <img src="http://latex.codecogs.com/svg.latex? $y-HX\beta = \Psi \xi$ " border="0"/>, here <img src="http://latex.codecogs.com/svg.latex? $\xi$ " border="0"/> is the correction term.  Then we have   
<img src="http://latex.codecogs.com/svg.latex? $s = X\beta + QH^T\xi$ " border="0"/>
To obtain the solution we just need to solve this equations:  
<img src="http://latex.codecogs.com/svg.latex? $\begin{pmatrix}
 \Psi & \Phi \\ 
\Phi^T & 0
 \end{pmatrix} \begin{pmatrix}
 \xi  \\ 
\beta
 \end{pmatrix} = \begin{pmatrix}
 y \\ 
0
 \end{pmatrix}$ " border="0"/>
 
 More information can be found at <a href="http://link.springer.com/article/10.1007/s10596-013-9364-0">this article</a>.
###3. DIRECTORIES AND FILES:

	./examples/		:	Example input C++ codes; Needed to read input from user or from input file.
	./src/			:	Source code in C++
	./header/		:	Relevant header files
	./exec/			:	Executables for FLIPACK
	./input/		:	The input file.
	./README.md		:	This file
	./License.md	:	License file
	./Makefile		:	Makefile  

###4. TUTORIAL
####4.1 To Get Started

1. To run this package, you need to have <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">**Eigen**</a> and <a href="https://github.com/sivaramambikasaran/BBFMM2D">**BBFMM2D**</a>.
    * Set Eigen:
      
	    1). Download Eigen from here: <http://eigen.tuxfamily.org/index.php?title=Main_Page> 
	     
	    2).  Create a directory named `Codes/` inside the main Eigen folder and copy the directory  `FLIPACK/` into the directory `Codes/`. 
	     
	    3).  Open the Makefile, which is inside the folder FLIPACK. Ensure that you have included the path to Eigen in the line containing `CFLAGS`. For instance, in the above setting, the path `"-I ./../../"` should be included in the Makefile.  
    * Set BBFMM2D:  
    
	    1). Download BBFMM2D from here <https://github.com/sivaramambikasaran/BBFMM2D>  
	    
	    2).  Copy directory `BBFMM2D/` inside of the directory `Codes/`.

2. To check whether things are set up correctly, you can perform the following action: Go to the directory where Makefile is in, then key in the following three commands in the terminal:

		make binary_file_mykernel
		cd exec/
		./FLIPACK_binary_file_mykernel
		
####4.2 Basic usage

#####4.2.1 FLIPACK with standard kernel

The basic usage of FLIPACK with standard kernel is as follows: 

	#include"FLIPACK_Header.hpp"
	…
	{
	...
	unsigned long N;            	//  Number of unknowns;
    unsigned m;                 	//  Number of measurements;
    unsigned nMeasurementSets;  	//  Number of measurement sets;
    unsigned short p;           	//  Number of terms in the structure;
    vector<Point> location;     	//  Location of the unknowns;
    double* Htranspose;         	//  Transpose of the measurement operator;
    double* measurements;  		 	//  Actual measurements;
    double* R;             		 	//  Covariance matrix;
    double* X;             		 	//  Structure of the mean;
    unsigned short nChebNode = 8;  //  Number of Chebyshev nodes per dimension
    …
    H2_2D_Tree Atree(nChebNode, Htranspose, location, N, m);// Build the fmm tree;
    …
    FLIPACK<kernel_Gaussian> A(Htranspose, X, measurements, R, N, m, p, nMeasurementSets, &Atree);     
    double* solution;    
    A.get_Solution(solution);
    …
    }

This example first build a fmm tree with this line:  
`H2_2D_Tree Atree(nChebNodes, charges, location, N, m);`  
where H2_2D_Tree is a class of fmm tree, the constructor takes 5 arguments:  

* nChebNodes(unsigned short):   
	Number of Chebyshev nodes per dimension. It should take value $\ge$ 3, and we recommend to take value from 3 to 10. (Larger number of Chebyshev nodes would give better result but with much more time)
* charges(double*):   
	All the different sets of charges. This pointer should point to an array with size $N \times m$, and the data should be stored in column-wise. ( i.e. first set of charges, followd by second set of charges, etc)
* location(vector<Point>):  
	Locations of the charges in $2D$ domain. Here Point is a structure type with $x$ and $y$ coordinate defined.  
* N(unsigned long):  
	Number of charges.  
* m(unsigned):  
	Number of sets of charges.  
	
Once you have built the FMM tree, you can solve the linear inversion problem with as many kernels as you want. This code shows an example using Gaussian kernel:  

	FLIPACK<kernel_Gaussian> A(Htranspose, X, measurements, R, N, m, p, nMeasurementSets, &Atree);
	double* solution;
    A.get_Solution(solution);
     
Here FLIPACK is a template class, and it takes 9 arguments:  
(To see what these values are, you read part **2** of this document, or read <a href="http://link.springer.com/article/10.1007/s10596-013-9364-0">this article</a>)
  
* Htranspose(double* ):  
	A pointer to <img src="http://latex.codecogs.com/svg.latex?  $H^T$ " border="0"/> ( the transpose of the measurement operator ), elements of <img src="http://latex.codecogs.com/svg.latex?  $H^T$ " border="0"/>  is stored column-wise in an array, i.e. first column of <img src="http://latex.codecogs.com/svg.latex?  $H^T$ " border="0"/> , followed by second column of <img src="http://latex.codecogs.com/svg.latex?  $H^T$ " border="0"/> , etc. 
* X(double* ):  
	A pointer to X ( structure of the mean ), elements of X is stored column-wise.  
* measurements(double* ):  
	A pointer to y ( actual measurements ), elements of y is stored column-wise.  
* R(double* ):  
	A pointer to R ( Covariance matrix ), elements of R is stored column-wise. 
* N(unsigned long):  
	Number of charges.  
* m(unsigned):  
	Number of sets of charges.  
* p(unsigned short):   
	Number of terms in the structure.
* nMeasurementSets(unsigned):  
	Number of measurement sets.  
* &Atree(H2_2D_Tree* ):  
	A pointer to a FMM tree.  

In this example, the type of kernel we use is Gaussian kernel (kernel_Gaussian), and we have provided several standard kernels (see **4.2.2**)  
The unknown of linear inversion problem is computed via 

	double* solution;
    A.get_Solution(solution);
    
Here `get_Solution(solution)`is a method of FLIPACK. The unknown is stored column-wise in `solution`. By calling this function we can get all members defined in class FLIPACK (see **4.2.3**)  
 

#####4.2.2 Options of provided kernels  
We have provided several standard kernels:  
The entries of the covariance matrix are given by <img src="http://latex.codecogs.com/svg.latex?  $Q_{ij} = k(x_i, y_i )$ " border="0"/>, where <img src="http://latex.codecogs.com/svg.latex?  $x_i$ " border="0"/> and <img src="http://latex.codecogs.com/svg.latex?  $y_i$ " border="0"/> are locations of points. Below are the details of the kernel functions we have provided:

Options of kernels:  

* LOGARITHM kernel:           
	usage: kernel_Logarithm  
	kernel function:  
    <img src="http://latex.codecogs.com/svg.latex?  $k(x,y) = 0.5 \times log(r^2)\, (r\neq 0);\, k(x,y)= 0 \,(r=0).$ " border="0"/> 
	
	
* ONEOVERR2 kernel:  
	usage: kernel_OneOverR  
	kernel function:  
    <img src="http://latex.codecogs.com/svg.latex?  $k(x,y) = 1 / r^2 \,(r \neq 0);\, k(x,y)= 0 \,(r=0)$. (r = |x-y|)" border="0"/>   
	
* GAUSSIAN kernel:  
	usage: kernel_Gaussian  
	kernel function:  
	<img src="http://latex.codecogs.com/svg.latex? $k(x,y) = exp(-r^2)$. (r = |x-y|)" border="0"/>   
	
* QUADRIC kernel:  
	usage: kernel_Quadric  
	kernel function:  
	 <img src="http://latex.codecogs.com/svg.latex? $ k(x,y) = 1 + r^2$. (r = |x-y|)" border="0"/>   

* INVERSEQUADRIC kernel:  
	usage: kernel_InverseQuadric  
	kernel function:  
	 <img src="http://latex.codecogs.com/svg.latex? $k(x,y) = 1 / (1+r^2)$. (r = |x-y|)" border="0"/> 
	
* THINPLATESPLINE kernel:  
	usage:  kernel_ThinPlateSpline  
	kernel function:   <img src="http://latex.codecogs.com/svg.latex? $k(x,y) =  0.5 \times r^2 \times log(r^2 )\, (r \neq 0);\, k(x,y)=0\,(r=0). (r = |x-y|)$" border="0"/>    		
If you want to define your own kernel, please see **4.2.3**.  

#####4.2.3 FLIPACK with user defined kernels

The basic usage is almost the same as **4.2.1** except that you have to define your own routine of computing kernel. One example code is as follows: 

	#include"FLIPACK_Header.hpp"
	
	class myKernel: public kernel_Base {
	public:
    virtual double kernel_Func(Point r0, Point r1){
        //implement your own kernel here
        double rSquare	=	(r0.x-r1.x)*(r0.x-r1.x) + (r0.y-r1.y)*(r0.y-r1.y);
        return exp(-pow(pow(rSquare,0.5)/900.0;8,0.5));
    	}
	};
	
	{
    …
    H2_2D_Tree Atree(nChebNode, Htranspose, location, N, m);// Build the fmm tree;
    …
    FLIPACK<myKernel> A(Htranspose, X, measurements, R, N, m, p, nMeasurementSets, &Atree);     
    double* solution;    
    A.get_Solution(solution);
    …
    }
You can define your own kernel inside `kernel_Func(Point r0, Point r1)`, it takes two Points as input and returns a double value ( <img src="http://latex.codecogs.com/svg.latex? $Q_{ij}$." border="0"/>  ). 	

#####4.2.4 Usage of multiple kernels

You can also use multiple kernels (user defined kernels and standard kernels) in one file, but make sure to have different class names. 
e.g.  
	
	/* You can define your own kernel here */
	class myKernelA: public kernel_Base{
		public:
    	virtual double kernel_Func(Point r0, Point r1) {
    	...
    	}
	}
	class myKernelB: public kernel_Base{
		public:
    	virtual double kernel_Func(Point r0, Point r1) {
    	...
    	}
	}
	…
	
	{
	…
	/* Build the FMM tree */
	H2_2D_Tree Atree(nChebNodes, charges, location, N, m);
	
	/* You can define as many kernels as you want */
	...
	FLIPACK<myKernelA> A(Htranspose, X, measurements, R, N, m, p, nMeasurementSets, &Atree);     
    double* solutionA;    
    A.get_Solution(solutionA);
    
    FLIPACK<myKernelB> B(Htranspose, X, measurements, R, N, m, p, nMeasurementSets, &Atree);     
    double* solutionB;    
    B.get_Solution(solutionB);
    
    /* You can use as many standard kernels as you want */
    FLIPACK<kernel_Gaussian> C(Htranspose, X, measurements, R, N, m, p, nMeasurementSets, &Atree);     
    double* solutionC;    
    C.get_Solution(solutionC);
    
     FLIPACK<kernel_Logarithm> D(Htranspose, X, measurements, R, N, m, p, nMeasurementSets, &Atree);     
    double* solutionD;    
    D.get_Solution(solutionD);
    ...
    }

The basic usage is already domonstrated in **4.2.1** and **4.2.3**. Once you have built the FMM tree, you can use different kernels to solve the linear inversion problem without rebuilding the tree. You can choose kernel type from standard kernels given by us ( see **4.2.2** ), or you can define your own kernel ( see **4.2.3** )
	


####3.3 Methods of FLIPACK

The followings are methods of FLIPACK to get different values that we are interested in ( see **2** ):  
* `void get_QHtranspose(double* &QHtranspose)`:  
	This function obtains the cross covariance, and the matrix is stored column-wise in QHtranspose.  
* `void get_HQHtranspose(double* &HQHtranspose)`:  
	This function obtains the measurement operator corrected covariance, and the matrix is stored column-wise in HQHtranspose.  
* `void get_Psi(double* &Psi)`:  
	This function obtains the measurement corrected covariance, and the matrix is stored column-wise in Psi.  
* `void get_Phi(double* &Phi)`:  
	This functiin obtains the measurement operator corrected structure, and the matrix is stored column-wise in Phi.  
* `void get_Xi(double* &Xi)`:    
	This function obtains the correction term, and the matrix is stored column-wise in Xi.      
* `void get_Beta(double* &Beta)`:    
	This function obtains unknown drift coefficients, and the matrix is stored column-wise in Beta.   
* `void get_Solution(double* &Solution)`:  
	This function obtains the unknown of linear inversion problem, by calling which all the values listed above are computed, and the matrix is stored column-wise in Solution.
	
	
###5. ROUTINES FOR INPUTING AND OUTPUTING DATA
 
####5.1 Reading meta data from text file

We have provided several routines for reading data from text file and binary file, and writing data into binary file.	

	void read_Metadata (const string& filenameMetadata, unsigned long& N, unsigned short& p,unsigned& m, unsigned& nMeasurementSets);
	
The first argument, filenameMetadata is the filename for your metadata. The number of unknowns is stored in N; the number of columns of X is stored in p; the number of measurements is stored in m; the number of sets of measurements is stored in nMeasurementsSets.

**File format:**  
 
`Number of unknowns, Number of terms in the structure, Number of measurements, Number of sets of measurements`

For example:

	3245,6,288,6

####5.2 Reading from binary file  

#####5.2.1 Read locations and transpose of H(.bin)

	void read_Location_Charges_binary(const string& filenameLocation, unsigned long N, vector<Point>& location, const string& filenameHtranspose,unsigned m, double* Htranspose);

The first argument filenameLocation and the forth argument filenameHtranspose are binary file names for location and <img src="http://latex.codecogs.com/svg.latex? $H^T$" border="0"/> respectively. N is the number of unknowns and m is the number of measurements. The data of locations is stored in `location` and the data of <img src="http://latex.codecogs.com/svg.latex? $H^T$" border="0"/>  is stored in `Htranspose` column-wise.  

**File format:** 

1. Binary file for Location:

	Elements are stored column-wise:
		
		loc0.x
		loc1.x  
		…
		loc0.y
		loc1.y
		...
 
2. Binary file for <img src="http://latex.codecogs.com/svg.latex? $H^T$" border="0"/>: 

	Elements of <img src="http://latex.codecogs.com/svg.latex? $H^T$" border="0"/> is stored this way:  
	It should be stored column-wise, i.e.   
	first column of <img src="http://latex.codecogs.com/svg.latex? $H^T$" border="0"/>, followed by second column of <img src="http://latex.codecogs.com/svg.latex? $H^T$" border="0"/>, etc.


#####5.2.2 Read X, R and measurents(.bin)

	void read_X_R_Measurements_Binary (const string& filenameX, unsigned long N, unsigned short p,double*& X,const string& filenameR, unsigned m, double*& R,const string& filenameMeasurement, unsigned nMeasurementSets, double*& measurements);

X, R and measurements should be stored column-wise in seperate binary files. By column-wise, we mean first column, followed by second column, etc. The names of arguments should be self explanatory.

This function will store X column-wise in `X`, store R column-wise in `R`, and store measurements column-wise in `measurements`.	

####5.3 Read from text file

#####5.3.1 Read locations and transpose of H(.txt)
The prototype of function to read input from text file is:  

	void read_Location_Charges (const string& filename, unsigned long N, vector<Point>& location, unsigned m, double*& charges);

The first argument is the filename of your text file, the second argument N and forth argument m are the number of locations and number of sets of charges respectively.
This function stores location in `location` and stores charges column-wise in `charges`.

**File format:**  
For each row, it should start with locations, and followed by a row in <img src="http://latex.codecogs.com/svg.latex? $H^T$" border="0"/>  ( if we do <img src="http://latex.codecogs.com/svg.latex? $QH^T$" border="0"/>  multiplication, <img src="http://latex.codecogs.com/svg.latex? $H^T$" border="0"/>  is the R.H.S.). Here note that elements should be separated using ','. If some element is 0, you can leave it as empty instead of 0. If all the elements in a row is 0, nothing need to be typed after the location.(spaces are allowed)

The row should look like this:  
  
`(location[0],location[1]) (elem1,elem2,elem3,elem4,…,elemn)`

For example:

	(-0.999984,-0.676221)   	(0.480685,0.869803,-0.188232,0.548587,-0.771039,0.73709,0.126494)    
	(0.869386,-0.5408)  
	(0.0655345,0.891162) (-0.193033,,-0.0287383,,-0.520512,0.33891,)  
	(0.342299,-0.246828) (0.0732668,,,,,,0.0951028)  
	(-0.984604,-0.44417) (,0.782447,-0.867924,0.485731,-0.729282,-0.481031,0.541473)  

#####5.3.2 Read X, R and measurents(.txt)	
	
	void read_X_R_Measurements (const string& filename, unsigned long N, unsigned short p, unsigned m, unsigned nMeasurementSets, double* X, double* R, double* measurements);

The first argument is the filename of input. Other arguments are the same as **5.2.2**  

**File format:**
For each line, it should be a row of one of these matrices(X,R,measuremtns). And the file should start with rows of X, then rows of R, and rows of measurements(in order).  
    Each row should look like this:    
    `(elem1,elem2,elem3,…,elemn)`  
    If some element is 0, you can leave it as empty instead of 0.

 For example:  
    
		(0.819666,-0.996573,-0.21957,-0.532382,0.491672,0.730317)
		(,0.332605,0.0938555,0.0442606,-0.63103,0.322981)
		(0.126619,,,0.262597,-0.470423,0.432109)
		(-0.00923308,0.0931162,-0.57049,-0.112485,0.283007,)

 The whole file should look like:

    	Rows of X 
    	Rows of R  
    	Rows of measurements 
    	
####5.4 Writing into binary file

	void write_Into_Binary_File(const string& filename, double* outdata, int numOfElems);  
	
This first argument is the filename for your output data. The second argument is a pointer to the output data, and the last argument is the number of elements in the array of your output data.

###6. EXAMPLES
####6.1 Chage input of examples

We have provided several examples for FLIPACK. Go to examples/, read through the files both must be self explanatory for the most part.
You can use our examples with your own input.

1. If you want to generate input through your own routine, and use the standard kernels:    
    Go to `/examples/`, open `"FLIPACK_get_matrix_through_routine_standard_kernel.cpp"`.  
    * To generate input through routine:   
      Change functions `get_Location()`, `get_Measurement_Operator()`, `get_X()` and `get_nChebNode()`.
    * To use standard kernels:   
      Choose the kernel type in `main()`, options of kernels are in **3.2.2**

	
2. If you want to read input from text file, and use standard kernels:  
    Go to the folder `/input/`, put your input file inside of this folder. 
    Go to the folder `/examples/`, open `"FLIPACK_textfile_standard_kernel.cpp"`.
     * To change input filename:  
      change the two lines in `main()`:  
      
      		string filename_location_Htranpose = "./../input/test_Location_H.txt";  
      		string filename_X_R_Measurements = "./../input/test_X_R_Measurements.txt";
    * To use standard kernels:   
      The same step as described in 1.


3. If you want to generate input through your own routine, and use your own kernel:

    Go to `/examples/`, open `"FLIPACK_get_matrix_through_routine_mykernel.cpp"`.
    * To define your own kernel:  
      Modify `class myKernel`. 
    * To generate your input:  
      The same step as described in 1.

4. If you want to read input from text file, and use your own kernel:

    Go to `/examples/`, open `"FLIPACK_textfile_mykernel.cpp"`.
    * To define your own kernel:  
      Modify `class myKernel`. 
    * To change input filename:  
      The same step as described in 2.  
      
5. If you want to read input from binary file, and use standard kernel:
 
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
     
6. If you want to read input from binary file, and use your own kernel:  

    Go to `/examples/`, open `"FLIPACK_binary_file_mykernel.cpp"`.  
    * To change the input filename:  
      The same step as described in 5.  
    * To define your own kernel:  
      Modify `class myKernel`.      

####6.2 Run examples  

Here we give an example:  
If you want to use `"FLIPACK_binary_file_standard_kernel.cpp"`

1. As stated earlier to run the code, go to the appropriate directory and key in the following:

		make binary_file_standard_kernel

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
   
3). `FLIPACK_textfile_mykernel.cpp`   
   key in:  
   
   		make input_from_file_mykernel
   		
4). `FLIPACK_textfile_standard_kernel.cpp`  
   key in:  
   
        make textfile_standard_kernel
   	
5). `FLIPACK_binary_file_mykernel.cpp`  
   key in:  
   
        make binary_file_mykernel  
        

