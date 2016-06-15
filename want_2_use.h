#ifndef WANT_2_USE_H_
#define WANT_2_USE_H_

/**
 * \brief Set flags to 1 if you want to use them, else to 0 
 *
 * USAGE:
 *
 *  \code 
    #if USE_BLA_BLA
    //here the stuff you want to use 
    #endif
 *  \endcode
 *
 * NOTE: when 'want_2_use.h' is to be used in a code snipppet it <BB> MUST </BB>
 * precede as <TT> #include "want_2_use.h" </TT> those headers where it is applied. Therefore, it is best to include it as first header in the <I> main file </I> 
 *
*/

#define USE_LAPACK 1

#define USE_NATIVE_LINUX_TIME_MEASUREMENT 1  //0 = standard, 1 = better coz

#define USE_TUNED_FLOATING_POINT_ARITHMETIC 0  //1 = yes, 0 = no

//plot stuff
#define USE_PLOTTING_OPTION 0 //1 = yes, 0 = no

#define USE_PLOTTING_WITH_ALL_SPECIES 1

//! use some update formulations for calculating Hessians, e.g. BFGS
#define USE_APPROX_HESSIAN 0    //1 = yes, 0 = no 


//3RD PARTY SOFTWARE · 3RD PARTY SOFTWARE · 3RD PARTY SOFTWARE · 3RD PARTY SOFT
//3RD PARTY SOFTWARE · 3RD PARTY SOFTWARE · 3RD PARTY SOFTWARE · 3RD PARTY SOFT

#define USE_CPPAD 1   //att.: used with 0 for the long-run test
#define USE_IPOPT 1   //1 = yes, 0 = no
#define USE_UMFPACK 1 //1 = yes, 0 = no

//3RD PARTY SOFTWARE · 3RD PARTY SOFTWARE · 3RD PARTY SOFTWARE · 3RD PARTY SOFT
//3RD PARTY SOFTWARE · 3RD PARTY SOFTWARE · 3RD PARTY SOFTWARE · 3RD PARTY SOFT


////#define USE_HESSIAN_REGULARIZATION 1  //1 = yes, 0 = no

//! some output can be generated
#define SHOW_ME_SOME_OUTPUT 0  // 0/1 no/a calculation is printed to screen

//#define PRINT_2_SCREEN

//! try to desingularize matrix, not recommended for general usage!
#define DESINGULARIZE_MATRIX

#endif //!include guard
