/****************************************************************************************************************************************************************THIS FILE CONTAINS GLOBALLY DEFINED PARAMETERS THAT CAN BE ADJUSTED ON DEMAND ****************************************************************************************************************************************************************/
#ifndef PARAM_H
#define PARAM_H

#include<cfloat>
#include<cmath> 

//Pi
//const double PI = 3.141592654;
//maybe use cmath definition:
const double PI = M_PI;

//terminal height and width
const int TERMWIDTH = 80;
const int TERMHEIGHT = 24;

//value from which on entries are considered as zero - especially for triangular matrices or reduced Hessenberg matrices:
const double SMALL = 1.0e-13;


//partition of a matrix into q_part times p_part processions or qq_part times r_part, etc.:
const int q_part = 3;
const int p_part = 2;
const int r_part = 2;
const int qq_part = 4;
const int s_part = 3;


//Factor for "tol" in the QR-algorithm for the unsymmetric eigenproblem:
const double FUCK = 1.5;

//Factor in the inverse iteration method. The constant should be of order unity, i.e. a small unspecified constant in the same order of magnitude as one.
const double ORDUNI=1.2; 
 
//Precision for ouput of matrices (and vectors)

const int VAL = 9;  //6    //precision value of output
const int BROAD = 3;//can be changed when output char is wider than BROAD.



//Since the Francis QR algorithm for determining eigenvalues has QUADRATIC (even CUBIC) convergence, a for loop may replace the while loop in matrix<T>::eigenqr() with boundary value UP:

const int UP = 60;

//maxiumum iterations in matrix<T>::eigenqr() when using the while-loop
const int MITER = 100;

/*the floating point base, also called "radix". Every floating point number is represented via 
     
        mantissa*base^exponent,
  
 where base is usually 2.0. However on certain IBM architectures base is 16 !!!
  
*/

const double RADIX = FLT_RADIX;

//maximum number of iterations for the power iteration and inverse iteration

const int MAXITER = 30; 


//A tiny number which used, for example, in the Guassian elimination when the pivot is zero because of a singular matrix

const double TINY = 1.0e-20;
const double WINY = 2.78e-20;

//perturbation of  a complex number where SIGNIFIC refers to the digits beyond the "." where the number will be changed 

const int SIGNIFIC = 5;

//this is a bound for various applications....

const double TOL = 1.e-13;



//for output-precision: cout.precision(PREC)

const int PREC = 13;



/**************************************************************************************************  FOR READING IN .txt, .dat TABLES *********************************************************************************************************/
const int MAX_TPOINTS = 1300000; //30   //equals row size of the datamatrix
const int MAX_SPECIES = 100;  //200  //equals column size of the datamatrix


const int STRING_SIZE = 12;   //default length of string that stores table entries

const int i_want_it_looong = 256;

//for multiple data treatment
const int N_FIELD = 4;
const int TIMEPTS = 30;
const int NUMOFSPECIES =10 ;
const int STRSIZE = 100; 
const int IGNO = 1024;  //ignore at most 1024 characters

const int MAX_SIZE = 1000;

#endif
