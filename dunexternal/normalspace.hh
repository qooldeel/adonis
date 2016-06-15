#ifndef COMPUTE_THE_NORMAL_SPACE_HERE_HH
#define COMPUTE_THE_NORMAL_SPACE_HERE_HH

#if USE_LAPACK
#ifdef LAPACK_USAGE

#include<complex>
#include<algorithm>

#include "../dunestuff/fmatrix.hh" 
#include "../common/error.hh"
#include "../common/adonisassert.hh"
#include "../lapackwrapper/lapackdrivertraits.hh"
#include "../lapackwrapper/wrappertraits.hh"
#include "../misc/useful.hh"
#include "dunecontainercopy.hh"

namespace Adonis{
  
  //!transpose of FieldMatrix A -- needed in what comes
  template<class K, int M, int N> 
  inline Dune::FieldMatrix<K,N,M> transpose(const Dune::FieldMatrix<K,M,N>& A){
    
    Dune::FieldMatrix<K,N,M> A_T;
    //A_T = K();

    for(int i = 0; i < M; ++i)
      for(int j = 0; j < N; ++j)
       	A_T[j][i] = A[i][j];
    
    return A_T;
  }

  
  /**
   *  \brief QR decomposition of a DENSE matrix, i.e. 
   *  \f$ A = Q\cdot R.\f$
   *  Uses DENSE 'Xgegrf' and 'Xorgqr' routines where 
   *  \f$ \textrm{X} \in \{\textrm{s,d,c,z}\}. \f$
   *  
   *  \tparam C character, either s(ingle), d(ouble), c(omplex), z(complex16)
   *  \tparam M rows of input matrix A
   *  \tparam N columns of input matrix A
   *
   *  \return Dune::FieldMatrix<typename F77Type<C>::Type,M,M> the Q-matrix which is orthonormal! The R matrix is neither stored nor needed for my applications.
   *
   *  EXAMPLE:
   *  \code
   
      Dune::FieldMatrix<double,4,3> A;
      A[0][0] = 1.;   A[0][1] = 2.; A[0][2] = 3.;
      A[1][0] = -3.;  A[1][1] = 2.; A[1][2] = 1.;
      A[2][0] = 2.;   A[2][1] = 0.; A[2][2] = -1.;
      A[3][0] = 3.;  A[3][1] = -1.; A[3][2] = 2.;

      std::cout << orthonormalize(A) << std::endl; 

   * \endcode
   *
   *  RESULT (you may wish to cross-check with Matlab(c) by invoking 
   *   \code
       A = [1,2,3; -3,2,1; 2,0,-1; 3,-1,2]
       [Q,R] = qr(A)
   *  \endcode
   *):
   * \code 
      -0.208514 -0.879191 0.156239 -0.398916
      0.625543 -0.414713 0.146474 0.644402
      -0.417029 -0.232239 -0.766547 0.429601
      -0.625543 0.033177 0.605426 0.490973
   *  \endcode
     
   *  Detail: Compute A = Q*R, where Q'*Q = I and R is upper triagonal.
   *   'dgeqrf' can be applied REGARDLESS of whether m or n is larger 
   *    (though in case of m <= n one may also apply 'dgerqf' instead)
   *    Alternatively:
   *     if M > N: A = Q[R 0]^T
     *     or better:
     *	   A = (Q1 Q2)[R 0]^T,
     *	   where Q1 consists of the first n columns of Q, and Q2 consists 
     *     of the remaining m-n columns.
     *
     *  see:
     *    <a href="http://www.nacad.ufrj.br/sgi/007-4325-001/sgi_html/ch03.html">See Eqs. 3-57. - 3-66.</a>
     *
     *Computation of  \f$ A = Q\cdot R\f$ but return Q in FORTRAN style, i.e. in transposed form.
   * This may be beneficial in some computations that require the transpose of Q only. Moreover, it saves time due to the fact that transposition is avoided. 
   * 
   * \tparam K basic type (no need to specify between 's', 'd', 'c' and 'z' 
   * \code
     orthonormalize<double>(A);   //see example above
   * \endcode
   */
 template<typename K, int M, int N>
  inline Dune::FieldMatrix<K,M,M> orthonormalize(const Dune::FieldMatrix<K,M,N>& A){
   const int Tdim = (M >= N) ? M : N;  //get maximum of M and N 

    int m = M,              
      n = N,                
      lda = std::max(1,M),              
      tdim =  Tdim,         
      lwork = M,                                
      info,                 
      k = N;                 

   

    K* tau = new K [tdim];                 
    K* work = new K [lwork];                
    
    Dune::FieldMatrix<K,M,M> Q;  //!the matrix \f$ Q =[Q1,Q2]\in R^{m\times m} \f$
    Q = K();

    //fill Q -- column-wise storage (FORTRAN style)
    for(int i = 0; i < M; ++i){   
      for(int j = 0; j < N; ++j){   //Assumption: M >= N
	Q[j][i] = A[i][j];     
      }
    }

    //! ATTENTION: don't use FieldTypes here as this results in a seg. fault!

    Copy2DynamicArray<K,M,M> q(Q,'r'); //since we store it column-wisely above!
    
    //std::cout << "q = "<< std::endl; ManipulateDynamicArray<K,int,0,M*M>::output(q.get_ptr());
    
    

    //&Q[0][0]
    F77Type<TypeTraits<K>::Value>().GEQRF(&m, &n, q.get_ptr(), &lda, &tau[0], &work[0], &lwork, &info);

    if(info != 0)
      ADONIS_ERROR(LapackError, "LAPACK-failure occurred: info = "<<info<<" by '"<< TypeTraits<K>::Value <<"geqrf_'.");
    
    //&Q[0][0]
    F77Type<TypeTraits<K>::Value>().ORGQR(&m, &tdim, &k, q.get_ptr(), &m, &tau[0], &work[0], &lwork, &info);
    

    if(info != 0)
      ADONIS_ERROR(LapackError, "LAPACK-failure occurred: info = "<<info<<" by '"<< TypeTraits<K>::Value <<"orgqr_'.");	
     
    FillMatrix<ColumnMajor>::now(Q,q.get_ptr());

    delete[] tau;  //proper clean up
    delete[] work;


    return Q;

}


  template<int N, class K, int M>
   inline Dune::FieldMatrix<K,M,N> extract_Q1(const Dune::FieldMatrix<K,M,M>& Q){
   
     Dune::FieldMatrix<K,M,N> Q1;
    
     for(int i = 0; i < M; ++i){
       for(int j = 0; j < N; ++j){
	Q1[i][j] = Q[i][j]; 
       }
     }
     
     return Q1;
   }
  
  template<int N, class K, int M>
  inline Dune::FieldMatrix<K,N,M> extract_Q1_transposed(const Dune::FieldMatrix<K,M,M>& Q){
   
    Dune::FieldMatrix<K,N,M> Q1_T;
    
    //! this is the version I originally intended to use. No transposition
    //! afterwards needed
    for(int j = 0; j < N; ++j){
      for(int i = 0; i < M; ++i){
	Q1_T[j][i] = Q[i][j];
      }
    }	
    return Q1_T;
  }


  /* @brief Mind that the QR decomposition for a M x N matrix (M >= N) returns
   *    A = (Q1 Q2)[R 0]^T where Q1 contains the first N columns of Q and Q2 the
   *   remaining M-N ones.
   *   
   *   -----------------------------------------------------------------------
   *  | Hence Q2 is a M x (M-N) matrix (and can be seen as a representation   | 
   *  | of the unrepresented species and the normalspace, respectively)       |
   *   -----------------------------------------------------------------------  
   *
   *   Taking this into account, the 'orthonormalize' algo is applied. Herafter 
   *   the matrix Q2, viz. the normal space, is extracted 
    
   *  \tparam R number of reaction progress variables
   *  \tparam K basic type 
   *  \tparam M rows of input matrix Q
   *  \tparam N columns of input matrix Q
   * Suppose you have already obtained matrix Q. Retrieve matrix Q2 now by calling 
   * \code
      extract_normalspace<2>(Q);
   * \endcode
   *where the 1st template argument represents the number of reduced variables (a.k.a reaction progress variables)
   *
   *  \return Dune::FieldMatrix<K,M,M-R> M x M-R matrix (see above)
   */
  template<int R, class K, int M, int N>
  inline Dune::FieldMatrix<K,M,M-R> extract_normalspace(const Dune::FieldMatrix<K,M,N>& Q){
    adonis_assert(R < M);    //wouldn't make sense if R >= M

    Dune::FieldMatrix<K,M,M-R> Q2;
    //Q2 = K();
    
    for(int i = 0; i < M; ++i){
      for(int j = R; j < N; ++j){
	Q2[i][j-R] = Q[i][j]; 
      }
    }
    
    return Q2;
  }

  /**
   *\brief Does the same job as 'extract_normalspace' except that \f$ Q^T\f$ is considered as input. This is the case when Lapack-style is used since FORTRAN stores matrices row by row.
   *
   *  \tparam R number of reaction progress variables
   *  \tparam K basic type 
   *  \tparam M rows of input matrix Q_T
   *  \tparam N columns of input matrix Q_T
   *
   *  \return Dune::FieldMatrix<K,M-R,M> matrix \f$Q_2^T.\f$ 
   */
   template<int R, class K, int M, int N>
   inline Dune::FieldMatrix<K,M-R,M> extract_normalspace_transposed(const Dune::FieldMatrix<K,M,N>& Q){
    adonis_assert(R < M);    //wouldn't make sense if R >= M

    Dune::FieldMatrix<K,M-R,M> Q2_T;
    
    //! this is the version I originally intended to use. No transposition
    //! afterwards needed
    for(int j = R; j < N; ++j){
      for(int i = 0; i < M; ++i){
	Q2_T[j-R][i] = Q[i][j];
      }
    }	
    return Q2_T;
  }



  /**
   *  \brief extract space and retrieve Q2 in one step
   *   Compute normalspace and extract it in one function
   */
  template<int R, char C, int M, int N>
  inline Dune::FieldMatrix<typename F77Type<C>::Type,M,M-R> normalspace(const Dune::FieldMatrix<typename F77Type<C>::Type,M,N>& A){
      
    //that's Q2;
    return extract_normalspace<R>(orthonormalize<F77Type<C>::Value>(A));
  }

  
  /** \brief Overload function 'normalspace'
   *   test stuff without defining the char type explicitely.
   *  Instead of writing, e.g., normalspace<2,'d'>(A)
   *  you can also write normalspace<2,double>(A)
   *  thus omitting the need to state the char type for the Fortran routines 
   *  which a user might not know if he/she has never worked with Lapack before.
   *
   *  \tparam R number of reaction progress variables
   *  \tparam K basic type 
   *  \tparam M rows of matrix A
   *  \tparam N columns of matrix A
   *
   *  \return  Dune::FieldMatrix<typename F77Type<TypeTraits<K>::Value>::Type,M,M-R> 
   */
  template<int R, typename K, int M, int N>
  inline  Dune::FieldMatrix<typename F77Type<TypeTraits<K>::Value>::Type,M,M-R> normalspace(const Dune::FieldMatrix<typename F77Type<TypeTraits<K>::Value>::Type,M,N>& A){
      
    return extract_normalspace<R>(orthonormalize<F77Type<TypeTraits<K>::Value>::Value>(A));
  }
    
  
/**
   * \brief Calculate transposed normalspace and extract it in one function
   * \tparam R, K, M, N see above
   */
   template<int R, typename K, int M, int N>
   inline  Dune::FieldMatrix<K,M-R,M> normalspace_transposed(const Dune::FieldMatrix<K,M,N>& A){
      
     return extract_normalspace_transposed<R>(orthonormalize<K>(A));
     
     //! should be equivalent to
     //return transpose(extract_normalspace<R>(orthonormalize<K>(A)));
   }

  
}

#endif  //LAPACK_USAGE
#endif  //USE_LAPACK


#endif
