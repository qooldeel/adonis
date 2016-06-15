#ifndef LINKAGE_CONVENTION_FOR_USING_LAPACK_BLAS_HH
#define LINKAGE_CONVENTION_FOR_USING_LAPACK_BLAS_HH

#include<iostream>
#include<cstdlib>
#include<complex>

#if USE_LAPACK
#ifdef LAPACK_USAGE


namespace Adonis{

  
  /**
   *\brief Frequently used <TT> LAPACK </TT> and <TT> BLAS </TT> routines. 
   *
   * NOTE1: Don't forget the underscore "_", hunk!
   * 
   * NOTE2: Apparently, Lapack doesn't like 'const'. So don't use. e.g. 'const double*'. Instead use 'double*' only. 
   *        This applies, of course, to all of your routines built on Lapack.
   */
  extern "C"{ 
    
    //computes C <-- alpha·AB + beta·C 
    extern void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);

    extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

    extern void cgemm_(char*, char*, int*, int*, int*, std::complex<float>*, std::complex<float>*, int*, std::complex<float>*, int*, std::complex<float>*, std::complex<float>*, int*);
    
    extern void zgemm_(char*, char*, int*, int*, int*, std::complex<double>*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*);
    


    //computes C := alpha*A*B + beta*C, where A is symmetric/hermitian
    extern void ssymm_(char*, char*, int*, int*, float*,  float*, int*,float*,int*, float*, float*, int*);

    extern void dsymm_(char*, char*, int*, int*, double*,  double*, int*,double*,int*, double*, double*, int*);

    extern void csymm_(char*, char*, int*, int*, std::complex<float>*,  std::complex<float>*, int*,std::complex<float>*,int*, std::complex<float>*, std::complex<float>*, int*);

    extern void zsymm_(char*, char*, int*, int*, std::complex<double>*,  std::complex<double>*, int*,std::complex<double>*,int*, std::complex<double>*, std::complex<double>*, int*);



    extern void sgtsv_(int*, int*, float*, float*, float*, float*, int*, int*);

    extern void dgtsv_(int*, int*, double*, double*, double*, double*, int*, int*);
    extern void cgtsv_(int*, int*, std::complex<float>*, std::complex<float>*, std::complex<float>*, std::complex<float>*, int*, int*);

    extern void zgtsv_(int*, int*, std::complex<double>*, std::complex<double>*, std::complex<double>*, std::complex<double>*, int*, int*);


    //QR factorisation of A in R^(m x n) (mainly used when m >= n) 
    
    extern void sgeqrf_(int*, int*, float*, int*, float*, float*, int*, int*);

    extern void dgeqrf_(int*, int*, double*, int*, double*, double*, int*, int*);
    extern void cgeqrf_(int*, int*, std::complex<float>*, int*, std::complex<float>*, std::complex<float>*, int*, int*);

    extern void zgeqrf_(int*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, int*);


    //RQ factorisation of A in R^(m x n) (mainly used when m < n)
    extern void sgerqf_(int*, int*, float*, int*, float*, float*, int*, int*);
    
    extern void dgerqf_(int*, int*, double*, int*, double*, double*, int*, int*);
    extern void cgerqf_(int*, int*, std::complex<float>*, int*, std::complex<float>*, std::complex<float>*, int*, int*);

    extern void zgerqf_(int*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, int*);

    //Needed if you want the orthogonal Q of an QR factorisation.
    //Driver: orthogonal matrix ("or"), generate Q after ("gqr"x)
    extern void sorgqr_(int*, int*, int*, float*, int*, float*, float*, int*, int*);
    extern void dorgqr_(int*, int*, int*, double*, int*, double*, double*, int*, int*);
  
    //in complex precision, you must take the UNITARY (driver: "un") pendant: 
    extern void cungqr_(int*, int*, int*, std::complex<float>*, int*, std::complex<float>*, std::complex<float>*, int*, int*);

    extern void zungqr_(int*, int*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, int*);

    
    //Solving a dense general m x n system using LU factorisation (this needs a a call by _getrs preceded by _getrf for the pivoting.
    extern void sgetrf_(int*, int*, float*, int*, int*, int*);

    extern void sgetrs_(char*, int*, int*, float*, int*, int*, float*, int*, int*);

    extern void dgetrf_(int*, int*, double*, int*, int*, int*);

    extern void dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);

    extern void cgetrf_(int*, int*, std::complex<float>*, int*, int*, int*);

    extern void cgetrs_(char*, int*, int*, std::complex<float>*, int*, int*, std::complex<float>*, int*, int*);

    extern void zgetrf_(int*, int*, std::complex<double>*, int*, int*, int*);

    extern void zgetrs_(char*, int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);


    //if only a SQUARE system is to be solved use the following
    extern void sgesv_(int*, int*, float*, int*, int*, float*, int*, int*);
    
    extern void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);

    extern void cgesv_(int*, int*, std::complex<float>*, int*, int*, std::complex<float>*, int*, int*);

    extern void zgesv_(int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);


    //solve general linear system using QR (might be faster than the SVD)
    //NOTE: ?gels assumes A to have FULL RANK!!
    extern void sgels_(char*, int*, int*, int*, float*, int*, float*, int*, float*, int*, int*);
    
    extern void dgels_(char*, int*, int*, int*, double*, int*, double*, int*, double*, int*, int*);
    
    extern void cgels_(char*, int*, int*, int*, std::complex<float>*, int*, std::complex<float>*, int*, std::complex<float>*, int*, int*);
    
    extern void zgels_(char*, int*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>*, int*, int*);


    //can also be used to calculate the inverse together with Xgetrf_
    extern void sgetri_(int*, float*, int*, int*, float*, int*, int*);
    
    extern void dgetri_(int*, double*, int*, int*, double*, int*, int*);

    extern void cgetri_(int*, std::complex<float>*, int*, int*, std::complex<float>*, int*, int*);

    extern void zgetri_(int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);


    //solve linear systems using SV -- may be RANK-DEFICIENT! (faster routine than xgelss, but needs more workspace)-- Note that the <I>complex</I> version has one aditional argument 'RWORK' (float or double array)
    extern void sgelsd_(int*, int*, int*, float*, int*, float*, int*, float*, float*, int*, float*, int*, int*, int*);

    extern void dgelsd_(int*, int*, int*, double*, int*, double*, int*, double*, double*, int*, double*, int*, int*, int*);

    extern void cgelsd_(int*, int*, int*, std::complex<float>*, int*, std::complex<float>*, int*, float*, float*, int*, std::complex<float>*, int*, float*, int*, int*);
    
    extern void zgelsd_(int*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, double*, double*, int*, std::complex<double>*, int*, double*,int*, int*);

    //solve ls using QR -- may be RANK-DEFICIENT!(faster routine than xgelsx but requires more workspace for block operations
    extern void sgelsy_(int*, int*, int*, float*, int*, float*, int*, int*, float*, int*, float*, int*, int*);

    extern void dgelsy_(int*, int*, int*, double*, int*, double*, int*, int*, double*, int*, double*, int*, int*);
    
    extern void cgelsy_(int*, int*, int*, std::complex<float>*, int*, std::complex<float>*, int*, int*, float*, int*, std::complex<float>*, int*, float*, int*);
    
    extern void zgelsy_(int*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, int*, double*, int*, std::complex<double>*, int*, double*, int*);
    

    //returns value of various norms (always real return values, of course
    extern float slange_(char*, int*, int*, float*, int*, float*);
    
    extern double dlange_(char*, int*, int*, double*, int*, double*);
    
    extern float clange_(char*, int*, int*, std::complex<float>*, int*, float*);

    extern double zlange_(char*, int*, int*, std::complex<double>*, int*, double*);
    

    //compute singular value decomposition of a general matrix 
    extern void sgesvd_(char*, char*, int*, int*, float*, int*, float*, float*, int*, float*, int*, float*, int*, int*);

     extern void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);

    extern void cgesvd_(char*, char*, int*, int*, std::complex<float>*, int*, float*, std::complex<float>*, int*, std::complex<float>*, int*, std::complex<float>*, int*, float*, int*);

    extern void zgesvd_(char*, char*, int*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>*, int*, double*, int*);


    //solve symmetric linear system, i.e. A = A^T
    extern void ssysv_(char*, int*, int*, float*, int*, int*, float*, int*, float*, int*, int*);

    extern void dsysv_(char*, int*, int*, double*, int*, int*, double*, int*, double*, int*, int*);

    extern void csysv_(char*, int*, int*, std::complex<float>*, int*, int*, std::complex<float>*, int*, std::complex<float>*, int*, int*);

    extern void zsysv_(char*, int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, int*);

    //solve a (symmetric) positive definit system, where the coefficient matrix obeys \f$ x^T\cdot Ax > 0, \ \forall x \not= 0. \f$ The Cholesky factorization is used to factor \f$ A\f$
    extern void sposv_(char*, int*, int*, float*, int*, float*, int*, int*);

    extern void dposv_(char*, int*, int*, double*, int*, double*, int*, int*);

    extern void cposv_(char*, int*, int*, std::complex<float>*, int*, std::complex<float>*, int*, int*);

    extern void zposv_(char*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, int*);


    //!computes the Cholesky factorization of a <B> s.p.d.</B> matrix \f$ A\f$
    //!if \f$ A\f$ is not s.p.d., then info < 0
    extern void spotrf_(char*, int*, float*, int*, int*);
    extern void dpotrf_(char*, int*, double*, int*, int*);
    extern void cpotrf_(char*, int*, std::complex<float>*, int*, int*); 
    extern void zpotrf_(char*, int*, std::complex<double>*, int*, int*); 

    //!solves solves sym.(hermit.) system using the Cholesky factorization 
    //! computed by the corresponding procedure <TT>xpotrf_</TT> above
    extern void spotrs_(char*, int*, int*, float*, int*, float*, int*, int*);
    extern void dpotrs_(char*, int*, int*, double*, int*, double*, int*, int*);
    extern void cpotrs_(char*, int*, int*, std::complex<float>*, int*, std::complex<float>*, int*, int*);
    extern void zpotrs_(char*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, int*);

    //computes all eigenvalues and, optionally, evecs of a symmetric matrix 
    extern void ssyev_(char*, char*, int*, float* ,int*, float*, float*, int*, int*);

    extern void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);

    //for hermitian matrices
    extern void cheev_(char*, char*, int*, std::complex<float>*, int*, float*, std::complex<float>*, int*, float*, int*); //has an additional argument

    extern void zheev_(char*, char*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*, double*, int*); //has an additional argument


    //computes all eigenvalues and, optionally, evecs of a general matrix
    extern void sgeev_(char*, char*, int*, float*, int*, float*, float*, float*, int*, float*, int*, float*, int*, int*);

    extern void dgeev_(char*, char*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);

    extern void cgeev_(char*, char*, int*, std::complex<float>*, int*, std::complex<float>*, std::complex<float>*, int*, std::complex<float>*, int*, std::complex<float>*, int*, float*, int*);

    extern void zgeev_(char*, char*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>*, int*, double*, int*);

  }//end extern "C"



}//end namespace

#endif //LAPACK_USAGE
#endif 

#endif
