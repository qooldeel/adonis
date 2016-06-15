#ifndef LAPACK_DRIVER_TRAITS_HH
#define LAPACK_DRIVER_TRAITS_HH

#if USE_LAPACK
#ifdef LAPACK_USAGE

#include <cmath>
#include <typeinfo>

#include "externC.hh"
#include "../misc/misctmps.hh"
#include "../common/isclass.hh"
#include "../common/adonisassert.hh"

namespace Adonis{

  /** 
   *   @brief Lapack driver traits 
   *   Firstly, define the lapack routines in "externC.hh" 
   *   
   *   LAPACK knows 4 basic types (1st letter of each routine), viz.  
   *
   *    s ... single precision              (float)
   *    d ... double precision              (double)                    (1)
   *    c ... single complex precision      (std::complex<float>)  
   *    z ... double complex precision      (std::complex<double>) 
   *
   *  
   *  Corresponding LAPACK routines can be found here:
   *  <a href="http://www.netlib.org/lapack//single/">Lapack - single</a>
   *  <a href="http://www.netlib.org/lapack//double/">Lapack - double</a>
   *  <a href="http://www.netlib.org/lapack//complex/">Lapack - complex</a>
   *  <a href"http://www.netlib.org/lapack//complex16/">Lapack -complex16</a>
   *
   * CONVENTION: skip the 1st letter which denotes the Lapack type
   */
  //traits for choosing Lapack precision and the routine according to the type
  //You can define Lapack routines here
  template<char w> class F77Type;   //w represents the precision, cf. (1)

  template<>
  class F77Type<'s'>{
  public:
    typedef float Type;
    typedef float BaseType;
    static const char Value = 's';
    static const bool IsComplex = false;

    static inline void GEMM(char* transA, char* transB, int* m, int* n, int* k, Type* alpha, Type* A, int* lda, Type* B, int* ldb, Type* beta, Type* C, int* ldc){
      sgemm_(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
    }


    static inline void SYMM(char* side, char* uplo, int* m, int* n, Type* alpha, Type* A, int* lda, Type* B, int* ldb, Type* beta, Type* C, int* ldc){
      ssymm_(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc);
    }

    static inline void GTSV(int* a, int* b, Type* c, Type* d, Type* e, Type* f, int* g, int* h){
      sgtsv_(a,b,c,d,e,f,g,h);
      
    }

     static inline void GEQRF(int* a , int* b, Type* c, int* d, Type* e, Type* f, int* g, int* h){
      sgeqrf_(a,b,c,d,e,f,g,h);
    }

    static inline void GERQF(int* a , int* b, Type* c, int* d, Type* e, Type* f, int* g, int* h){
      sgerqf_(a,b,c,d,e,f,g,h);
    }


    static inline void ORGQR(int* a, int* b, int* c, Type* d, int* e, Type* f, Type* g, int* h, int* i){
      sorgqr_(a,b,c,d,e,f,g,h,i);
    }

    static inline void GETRF(int* a, int* b, Type* c, int* d, int* e, int* f){
      sgetrf_(a,b,c,d,e,f);
    }

    static inline void GETRS(char* a, int* b, int* c, Type* d, int* e, int* f, Type* g, int* h, int* i){
      sgetrs_(a,b,c,d,e,f,g,h,i);
    }

    //solves general square linear system using LU decomposition
    static inline void GESV(int* a, int* b, Type* c, int* d, int* e, Type* f, int* g, int* h){
      sgesv_(a,b,c,d,e,f,g,h);
    }

    
    static inline void GELS(char* trans, int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, Type* work, int* lwork, int* info){
      sgels_(trans,m,n,nrhs,A,lda,B,ldb,work,lwork,info);
    }


    static inline void GETRI(int* n, Type* A, int* lda, int* ipiv, Type* work, int* lwork, int* info){
      sgetri_(n,A,lda,ipiv,work,lwork,info);
    }
  
    static inline void GELSD(int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, BaseType* s, BaseType* rcond, int* rank, Type* work, int* lwork, BaseType* rwork,int* iwork, int* info){
      sgelsd_(m,n,nrhs,A,lda,B,ldb,s,rcond,rank,work,lwork,iwork,info);
    }

    /** \brief Rank-deficient least squares. Based on QR decomposition.
     * Cost are the same as in the case of GELSD, namely \f$ \mathcal{O}(MN^2) \f$. 
     * Therefore xgelsd_ and xgelsy_ are the recommended routines.
     */
    static inline void GELSY(int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* jpvt, BaseType* rcond, int* rank, Type* work, int* lwork, BaseType* rwork, int* info){
      sgelsy_(m,n,nrhs,A,lda,B,ldb,jpvt,rcond,rank,work,lwork,info);
    }


    //gives back various norms 
    Type LANGE(char* norm, int* m, int* n, Type* A, int* lda, BaseType* work){
      return slange_(norm,m,n,A,lda,work);
    }
  
    static inline void GESVD(char* jobu, char* jobvt, int* m, int* n, Type* A, int* lda, BaseType* s, Type* u, int* ldu, Type* vt, int* ldvt, Type* work, int* lwork, BaseType* rwork, int* info){
      sgesvd_(jobu,jobvt,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,info);
    }

    //solve system for symmetric matrix \f$ A = A^T.\f$
    static inline void SYSV(char* uplo, int* n, int* nrhs, Type* A, int* lda, int* ipiv, Type* B, int* ldb, Type* work, int* lwork, int* info){
      ssysv_(uplo,n,nrhs,A,lda,ipiv,B,ldb,work,lwork,info);
    }
  

    //solve symmetric p.d. system. A matrix \$A\f$ is p.d. iff there exists a Cholesky-decomposition of it, i.e. \f$ A =GG^T, \ \det(G) \not= 0. \f$ 
    static inline void POSV(char* uplo, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* info){
      sposv_(uplo,n,nrhs,A,lda,B,ldb,info);
    }

    //! compute Cholesky factorization of symmetric p.d. matrix \f$ A \f$
    static inline void POTRF(char* uplo, int* n, Type* A, int* lda, int* info){
      spotrf_(uplo,n,A,lda,info);
    }

    //! solves p.d. system using the Cholesky factorization obtained by POTRF
    static inline void POTRS(char* uplo, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* info){
      spotrs_(uplo,n,nrhs,A,lda,B,ldb,info);
    }

  
    //!solve eigenvalues (or optionally eigenvecs) of symmetric systems. In the REAL case, rwork is just a dummy argument
    static inline void SYEV(char* jobz, char* uplo, int* n, Type* A, int* lda, Type* w, Type* work, int* lwork, BaseType* rwork, int* info){
      ssyev_(jobz,uplo,n,A,lda,w,work,lwork,info);
    }


    //!compute eigenvalues (and optinally eigenvectors) of a general n x n matrix
    //! w and rwork are only dummies here
    static inline void GEEV(char* jobvl, char* jobvr, int* n, Type* A, int* lda, Type* wr, Type* wi, std::complex<BaseType>* w, Type* vl, int* ldvl, Type* vr, int* ldvr, Type* work, int* lwork, BaseType* rwork, int* info){
      sgeev_(jobvl,jobvr,n,A,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info);
    }
  };

  template<>
  class F77Type<'d'>{
  public:
    typedef double Type;
    typedef double BaseType;
    static const char Value = 'd';
    static const bool IsComplex = false;

    static inline void GEMM(char* transA, char* transB, int* m, int* n, int* k, Type* alpha, Type* A, int* lda, Type* B, int* ldb, Type* beta, Type* C, int* ldc){
      dgemm_(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
    }

    static inline void SYMM(char* side, char* uplo, int* m, int* n, Type* alpha, Type* A, int* lda, Type* B, int* ldb, Type* beta, Type* C, int* ldc){
      dsymm_(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc);
    }


    static inline void GTSV(int* a, int* b, Type* c, Type* d, Type* e, Type* f, int* g, int* h){
      dgtsv_(a,b,c,d,e,f,g,h);
      
    }
  
    static inline void GEQRF(int* a , int* b, Type* c, int* d, Type* e, Type* f, int* g, int* h){
      dgeqrf_(a,b,c,d,e,f,g,h);
    }
    
    static inline void GERQF(int* a , int* b, Type* c, int* d, Type* e, Type* f, int* g, int* h){
      dgerqf_(a,b,c,d,e,f,g,h);
    }
 
    static inline void ORGQR(int* a, int* b, int* c, Type* d, int* e, Type* f, Type* g, int* h, int* i){
      dorgqr_(a,b,c,d,e,f,g,h,i);
    }
    
    static inline void GETRF(int* a, int* b, Type* c, int* d, int* e, int* f){
      dgetrf_(a,b,c,d,e,f);
    }

    static inline void GETRS(char* a, int* b, int* c, Type* d, int* e, int* f, Type* g, int* h, int* i){
      dgetrs_(a,b,c,d,e,f,g,h,i);
    }

    static inline void GESV(int* a, int* b, Type* c, int* d, int* e, Type* f, int* g, int* h){
      dgesv_(a,b,c,d,e,f,g,h);
    }
  
   
    static inline void GELS(char* trans, int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, Type* work, int* lwork, int* info){
      dgels_(trans,m,n,nrhs,A,lda,B,ldb,work,lwork,info);
    }
  
     static inline void GETRI(int* n, Type* A, int* lda, int* ipiv, Type* work, int* lwork, int* info){
      dgetri_(n,A,lda,ipiv,work,lwork,info);
    }

    static inline void GELSD(int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, BaseType* s, BaseType* rcond, int* rank, Type* work, int* lwork, BaseType* rwork, int* iwork, int* info){
      dgelsd_(m,n,nrhs,A,lda,B,ldb,s,rcond,rank,work,lwork,iwork,info);
    }
  
    static inline void GELSY(int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* jpvt, BaseType* rcond, int* rank, Type* work, int* lwork, BaseType* rwork, int* info){
      dgelsy_(m,n,nrhs,A,lda,B,ldb,jpvt,rcond,rank,work,lwork,info);
    }

    Type LANGE(char* norm, int* m, int* n, Type* A, int* lda, BaseType* work){
      return dlange_(norm,m,n,A,lda,work);
    }

    //int the real cases rwork is a dummy
    static inline void GESVD(char* jobu, char* jobvt, int* m, int* n, Type* A, int* lda, Type* s, Type* u, int* ldu, Type* vt, int* ldvt, Type* work, int* lwork, BaseType* rwork, int* info){
      dgesvd_(jobu,jobvt,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,info);
    }
  
    static inline void SYSV(char* uplo, int* n, int* nrhs, Type* A, int* lda, int* ipiv, Type* B, int* ldb, Type* work, int* lwork, int* info){
      dsysv_(uplo,n,nrhs,A,lda,ipiv,B,ldb,work,lwork,info);
    }
  
    static inline void POSV(char* uplo, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* info){
      dposv_(uplo,n,nrhs,A,lda,B,ldb,info);
    }
  
    static inline void POTRF(char* uplo, int* n, Type* A, int* lda, int* info){
      dpotrf_(uplo,n,A,lda,info);
    }
    
    static inline void POTRS(char* uplo, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* info){
      dpotrs_(uplo,n,nrhs,A,lda,B,ldb,info);
    }


    static inline void SYEV(char* jobz, char* uplo, int* n, Type* A, int* lda, Type* w, Type* work, int* lwork, BaseType* rwork, int* info){
      dsyev_(jobz,uplo,n,A,lda,w,work,lwork,info);
    }
  
    static inline void GEEV(char* jobvl, char* jobvr, int* n, Type* A, int* lda, Type* wr, Type* wi, std::complex<BaseType>* w, Type* vl, int* ldvl, Type* vr, int* ldvr, Type* work, int* lwork, BaseType* rwork, int* info){
      dgeev_(jobvl,jobvr,n,A,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info);
    }
  };



  template<>
  class F77Type<'c'>{
  public:
    typedef float BaseType;
    typedef std::complex<BaseType> Type;
    static const char Value = 'c';
    static const bool IsComplex = true;

    static inline void GEMM(char* transA, char* transB, int* m, int* n, int* k, Type* alpha, Type* A, int* lda, Type* B, int* ldb, Type* beta, Type* C, int* ldc){
      cgemm_(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
    }

    static inline void SYMM(char* side, char* uplo, int* m, int* n, Type* alpha, Type* A, int* lda, Type* B, int* ldb, Type* beta, Type* C, int* ldc){
      csymm_(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc);
    }

    static inline void GTSV(int* a, int* b, Type* c, Type* d, Type* e, Type* f, int* g, int* h){
      cgtsv_(a,b,c,d,e,f,g,h);
      
    }
  
    static inline void GEQRF(int* a , int* b, Type* c, int* d, Type* e, Type* f, int* g, int* h){
      cgeqrf_(a,b,c,d,e,f,g,h);
    }
  
    static inline void GERQF(int* a , int* b, Type* c, int* d, Type* e, Type* f, int* g, int* h){
      cgerqf_(a,b,c,d,e,f,g,h);
    }
  
    static inline void ORGQR(int* a, int* b, int* c, Type* d, int* e, Type* f, Type* g, int* h, int* i){
      //ATTENTION: you need the unitary version here ;-)
      cungqr_(a,b,c,d,e,f,g,h,i);
    }
    
    static inline void GETRF(int* a, int* b, Type* c, int* d, int* e, int* f){
      cgetrf_(a,b,c,d,e,f);
    }

    static inline void GETRS(char* a, int* b, int* c, Type* d, int* e, int* f, Type* g, int* h, int* i){
      cgetrs_(a,b,c,d,e,f,g,h,i);
    }
  

    static inline void GESV(int* a, int* b, Type* c, int* d, int* e, Type* f, int* g, int* h){
      cgesv_(a,b,c,d,e,f,g,h);
    }
  
  
    static inline void GELS(char* trans, int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, Type* work, int* lwork, int* info){
      cgels_(trans,m,n,nrhs,A,lda,B,ldb,work,lwork,info);
    }
  
     static inline void GETRI(int* n, Type* A, int* lda, int* ipiv, Type* work, int* lwork, int* info){
      cgetri_(n,A,lda,ipiv,work,lwork,info);
    }
  
    static inline void GELSD(int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, BaseType* s, BaseType* rcond, int* rank, Type* work, int* lwork, BaseType* rwork, int* iwork, int* info){
      cgelsd_(m,n,nrhs,A,lda,B,ldb,s,rcond,rank,work,lwork,rwork,iwork,info);
    }
  
    static inline void GELSY(int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* jpvt, BaseType* rcond, int* rank, Type* work, int* lwork, BaseType* rwork, int* info){
      cgelsy_(m,n,nrhs,A,lda,B,ldb,jpvt,rcond,rank,work,lwork,rwork,info);
    }
  
    BaseType LANGE(char* norm, int* m, int* n, Type* A, int* lda, BaseType* work){
      return clange_(norm,m,n,A,lda,work);
    }
  
    static inline void GESVD(char* jobu, char* jobvt, int* m, int* n, Type* A, int* lda, BaseType* s, Type* u, int* ldu, Type* vt, int* ldvt, Type* work, int* lwork, BaseType* rwork, int* info){
      cgesvd_(jobu,jobvt,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info);
    }
  
    static inline void SYSV(char* uplo, int* n, int* nrhs, Type* A, int* lda, int* ipiv, Type* B, int* ldb, Type* work, int* lwork, int* info){
      csysv_(uplo,n,nrhs,A,lda,ipiv,B,ldb,work,lwork,info);
    }
  
    static inline void POSV(char* uplo, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* info){
      cposv_(uplo,n,nrhs,A,lda,B,ldb,info);
    }
    
    static inline void POTRF(char* uplo, int* n, Type* A, int* lda, int* info){
      cpotrf_(uplo,n,A,lda,info);
    }

    static inline void POTRS(char* uplo, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* info){
      cpotrs_(uplo,n,nrhs,A,lda,B,ldb,info);
    }

    //now rwork isn't a dummy any more
    static inline void SYEV(char* jobz, char* uplo, int* n, Type* A, int* lda, BaseType* w, Type* work, int* lwork, BaseType* rwork, int* info){
      cheev_(jobz,uplo,n,A,lda,w,work,lwork,rwork,info); 
    }
    
    //now wi,wr are dummies
    static inline void GEEV(char* jobvl, char* jobvr, int* n, Type* A, int* lda, Type* wr, Type* wi, Type* w, Type* vl, int* ldvl, Type* vr, int* ldvr, Type* work, int* lwork, BaseType* rwork, int* info){
      cgeev_(jobvl,jobvr,n,A,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info);
    }

  };

  template<>
  class F77Type<'z'>{
  public:
    typedef double BaseType;
    typedef std::complex<BaseType> Type;
    static const char Value = 'z';
    static const bool IsComplex = true;


    static inline void GEMM(char* transA, char* transB, int* m, int* n, int* k, Type* alpha, Type* A, int* lda, Type* B, int* ldb, Type* beta, Type* C, int* ldc){
      zgemm_(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
    }

    static inline void SYMM(char* side, char* uplo, int* m, int* n, Type* alpha, Type* A, int* lda, Type* B, int* ldb, Type* beta, Type* C, int* ldc){
      zsymm_(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc);
    }


    static inline void GTSV(int* a, int* b, Type* c, Type* d, Type* e, Type* f, int* g, int* h){
      zgtsv_(a,b,c,d,e,f,g,h);
      
    }
 
    static inline void GEQRF(int* a , int* b, Type* c, int* d, Type* e, Type* f, int* g, int* h){
      zgeqrf_(a,b,c,d,e,f,g,h);
    }
 
    static inline void GERQF(int* a , int* b, Type* c, int* d, Type* e, Type* f, int* g, int* h){
      zgerqf_(a,b,c,d,e,f,g,h);
    } 

    static inline void ORGQR(int* a, int* b, int* c, Type* d, int* e, Type* f, Type* g, int* h, int* i){
      //ATTENTION: you need the unitary version here ;-)
      zungqr_(a,b,c,d,e,f,g,h,i);
    }
  
    static inline void GETRF(int* a, int* b, Type* c, int* d, int* e, int* f){
      zgetrf_(a,b,c,d,e,f);
    }

    static inline void GETRS(char* a, int* b, int* c, Type* d, int* e, int* f, Type* g, int* h, int* i){
      zgetrs_(a,b,c,d,e,f,g,h,i);
    }
  
    static inline void GESV(int* a, int* b, Type* c, int* d, int* e, Type* f, int* g, int* h){
      zgesv_(a,b,c,d,e,f,g,h);
    }
   
    
    static inline void GELS(char* trans, int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, Type* work, int* lwork, int* info){
      zgels_(trans,m,n,nrhs,A,lda,B,ldb,work,lwork,info);
    }
  
    static inline void GETRI(int* n, Type* A, int* lda, int* ipiv, Type* work, int* lwork, int* info){
      zgetri_(n,A,lda,ipiv,work,lwork,info);
    }
  
    static inline void GELSD(int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, BaseType* s, BaseType* rcond, int* rank, Type* work, int* lwork, BaseType* rwork, int* iwork, int* info){
      zgelsd_(m,n,nrhs,A,lda,B,ldb,s,rcond,rank,work,lwork,rwork,iwork,info);
    }
  
    static inline void GELSY(int* m, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* jpvt, BaseType* rcond, int* rank, Type* work, int* lwork, BaseType* rwork, int* info){
      zgelsy_(m,n,nrhs,A,lda,B,ldb,jpvt,rcond,rank,work,lwork,rwork,info);
    }
  
    BaseType LANGE(char* norm, int* m, int* n, Type* A, int* lda, BaseType* work){
      return zlange_(norm,m,n,A,lda,work);
    }
  
    static inline void GESVD(char* jobu, char* jobvt, int* m, int* n, Type* A, int* lda, BaseType* s, Type* u, int* ldu, Type* vt, int* ldvt, Type* work, int* lwork, BaseType* rwork, int* info){
      zgesvd_(jobu,jobvt,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info);
    }
  
    static inline void SYSV(char* uplo, int* n, int* nrhs, Type* A, int* lda, int* ipiv, Type* B, int* ldb, Type* work, int* lwork, int* info){
      zsysv_(uplo,n,nrhs,A,lda,ipiv,B,ldb,work,lwork,info);
    }
  
    static inline void POSV(char* uplo, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* info){
      zposv_(uplo,n,nrhs,A,lda,B,ldb,info);
    }

    static inline void POTRF(char* uplo, int* n, Type* A, int* lda, int* info){
      zpotrf_(uplo,n,A,lda,info);
    }

    static inline void POTRS(char* uplo, int* n, int* nrhs, Type* A, int* lda, Type* B, int* ldb, int* info){
      zpotrs_(uplo,n,nrhs,A,lda,B,ldb,info);
    }
  
    static inline void SYEV(char* jobz, char* uplo, int* n, Type* A, int* lda, BaseType* w, Type* work, int* lwork, BaseType* rwork, int* info){
        zheev_(jobz,uplo,n,A,lda,w,work,lwork,rwork,info); 
    }
  
    //now wi,wr are dummies
    static inline void GEEV(char* jobvl, char* jobvr, int* n, Type* A, int* lda, Type* wr, Type* wi, Type* w, Type* vl, int* ldvl, Type* vr, int* ldvr, Type* work, int* lwork, BaseType* rwork, int* info){
      zgeev_(jobvl,jobvr,n,A,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info);
    }
  };

  


  /**
   * \brief Sometimes complex LAPACK routines require more arguments (e.g. for
   * additional workspace) than the corresponding real versions. This TMP can be used to define, e.g., additional workspace dimensions etc.
   */
  template<bool B> class LapackComplex;

  template<>
  class LapackComplex<true>{  //COMPLEX
  public:
    template<class INT>
    static inline void lwk(INT& lrwork, const INT& m, const INT& n, const INT& smlsiz, const INT& nlvl, const INT& nrhs){
      (m >= n) ? (lrwork = 10*n + 2*n*smlsiz + 8*n*nlvl + 3*smlsiz*nrhs +
		  std::max( ntimes<2>(smlsiz+1), n*(1+nrhs) + 2*nrhs )) :
	(lrwork = 10*m + 2*m*smlsiz + 8*m*nlvl + 3*smlsiz*nrhs +
	 std::max( ntimes<2>(smlsiz+1), n*(1+nrhs) + 2*nrhs));
    }
    

    //needed for eigenvalues/vector determination
    template<class INT>
    static inline INT lwkev(const INT& n) {return std::max(1,3*n-2);}

    template<class INT>
    static inline INT  lwk(const INT& n) {return 2*n;}
  
    template<class INT>
    static inline INT dim(const INT& n) {return 0;}

    template<class INT>
    static inline INT lwk(const INT& n, const INT& fac) {return fac*n;}

    template<class C, class BT>
    static inline void assign_eigenvalues(C& ew, BT* wr, BT* wi){}

    //calculation of lwork might be done one complex values
    template<class T>
    static inline int real_part(const T& c){
      adonis_assert(is_class(c));  //make sure c is a class 
      adonis_assert(typeid(c) == typeid(std::complex<typename T::value_type>));
      return static_cast<int>(c.real());  //real part of a complex number
    }
  };
  
  template<>
  class  LapackComplex<false>{  //REAL
  public:
    template<class INT>
    static inline void lwk(INT& lrwork, const INT& m, const INT& n, const INT& smlsiz, const INT& nlvl, const INT& nrhs){
      lrwork = 0;
    }
  
    //needed for eigenvalues/vector determination
    template<class INT>
    static inline INT lwkev(const INT& n) {return 0;}

    template<class INT>
    static inline INT dim(const INT& n) {return n;}

    template<class INT>
    static inline INT  lwk(const INT& n) {return 0;}
  
    template<class INT> 
    static inline INT lwk(const INT& n, const INT& fac) {return 0;}
    
    template<class C, class BT>
    static inline void assign_eigenvalues(C& ew, BT* wr, BT* wi){
      for(unsigned i = 0; i < ew.size(); ++i){
	ew[i].real() = wr[i];
	ew[i].imag() = wi[i];
      }
    }

    template<class T>
    static inline int real_part(const T& c){
      adonis_assert(!is_class(c));
      return static_cast<int>(c);    //return just the real number
    } 
  };


  /**
   * \brief Computation of eigenvalues is often accompanied by calculation of eigenvectors. Here you can omit the expensive eigenvector calculation if you are only interested in eigenvalues.
   */
  template<class MatrixType, bool EV> class WithEigenvectorComputation;

  template<class MatrixType>
  class WithEigenvectorComputation<MatrixType,true>{
  public:
    typedef typename MatrixType::value_type ReturnType;
  
    template<class INT>
    static inline INT order(INT n){return n;}

    template<class INT>
    static inline INT evecsize(INT ldv, INT n){
      return ldv*n;
    }

    template<class INT>
    static inline void resize(ReturnType* evec, INT n){
      evec = new ReturnType[n];
    }
  
    static inline void delete_ptr(ReturnType* evec){
      delete[] evec;
    }
  };

  template<class MatrixType>
  class WithEigenvectorComputation<MatrixType,false>{
  public:
    typedef void ReturnType;
  
    template<class INT>
    static inline INT order(INT n){return 0;} 

    template<class INT>
    static inline INT evecsize(INT ldv, INT n){
      return 0;
    }

    template<class INT>
    static inline void resize(ReturnType* evec, INT n){}
  
    static inline void delete_ptr(ReturnType* evec){}
    
  };


}

#endif //..LAPACK_USAGE
#endif //

#endif
