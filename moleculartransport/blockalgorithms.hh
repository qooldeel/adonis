#ifndef BLOCK_MATRIX_ALGORITHMS_HH
#define BLOCK_MATRIX_ALGORITHMS_HH

#include <vector>

#include "../containers/staticarray.hh"

#include "../common/elementaryoperations.hh"
#include "../templatemetaprograms/matrixunroller.hh"

#if USE_LAPACK
#ifdef LAPACK_USAGE

#include "../lapackwrapper/wrappertraits.hh"
#include "../lapackwrapper/lapackdrivertraits.hh"

namespace Adonis{


  /**
   * \brief solves <BB>dense</BB> \f$ NQ \times NQ\f$ block diagonal system 
   * \f$ A\cdot X = D\f$ via Thomas algo. 
   *Note that a tridiagonal block matrix is a special sparse  quadratic matrix of order \f$ NQ \f$.
   * If Gaussian elimination were to be applied directly, \f$ \mathcal{O}(N^2)^3/3\f$ operations would be necessary. In contrast, Thomas algorithm only needs
   * \f$ \mathcal{O}(4N^2)^2 \f$ operations.
   * Assume that the diag and the sub- and supdiag are composed of 
   * \f$Q \times Q\f$ blocks each, where diag has N such blocks and sub- and 
   * supdiag are composed of N-1 blocks.
   
   \tparam N number of blocks
   \tparam Q order of block matrices, blocksize for short
   \tparam NRHS number of right hand sides
   \tparam T precision of computation (should be also type of all occurring ITs)
   \tparam C block algorithm under consideration, 't' or 'T' for Thomas algo
  */
  template<int N, int Q, int NRHS, class T, char C> class BlockTridiagonalAlgorithm;

  
  template<int N, int Q, int NRHS, class T>  //!Thomas algorithm for block matrices
  class BlockTridiagonalAlgorithm<N,Q,NRHS,T,'t'>{
  public:
    typedef T value_type;
    

    /**
     * \brief implementation of Algorithm (4.5.4) and (4.5.5) in [1,pp. 174]
     * <BB> CAUTION: We assume that all input matrices are stored in FORTRAN-
     * style, i.e. column-major order</BB>
     *  
     * Note: all matrices are going to be overwritten. <TT> rhs </TT> bears 
     *later the solution of the block system 
     *
     * \param low N-1 subdiagonal blocks of size Q x Q
     * \param diag N diagonal blocks of size Q x Q  
     * \param up N-1 superdiagonal blocks of size Q x Q 
     * \param rhs N right hand side blocks of size Q x NRHS
     * \param useTransposedForm (default: true): solve transposed or normal 
     *        systems
     *
     * CAUTION: Here, we assume that the matrices are stored in F77-style,
     *          i.e. column major (in contrast to the row-major format of 
     *          C/C++!!!!! <BB>Hence, on input, be sure to tranpose the matrices
     *          appropriately! </BB>
     *          Therefore we must perform a compile time transposition of the 
     *          non-vector arrays of the LHS of each lin. equation before
     *          the invocation of this member function
     *
     * References:
     * [1]  [GOLUB and VAN LOAN, "Matrix computations", 3rd ed., 1996, p. 174/75]
     */
    template<class IT1, class IT2, class IT3, class IT4>
    static inline void solve(IT1 low, IT2 diag, IT3 up, IT4 rhs){
      adonis_assert(N > 0 && Q > 0 && NRHS > 0);

      //typedef SubtractBasicElements<T> OpType;

      //! needed for solution of square blocks
      char trans = 'N';
    
      //((trans == 'T') && (F77Type<TypeTraits<value_type>::Value>::IsComplex)) ? trans = 'C' : trans = 'T';
      

      int q = Q,              //blocksize
	nrhs = NRHS,
	lda = Q > 1 ? Q : 1,  //max(1,Q)
	ldb = lda,
	info;
      
      int ipiv[Q];
      
      //the following is needed for xgemm_
      int ldb2 =  ((trans == 'N') || (trans == 'n')) ? Q : NRHS,
      	ldc = Q;
     
      value_type alpha = -1,  //these remain fixed throughout the algorithm
      	beta = 1;             //and are needed for BLAS' xgemm_, namely
                              // C <-- alpha·AB + beta·C
     
    
      typedef StaticArray<value_type,Q*Q> RowType;
      RowType U[N > 0 ? N : 1];  //! store diagonal blocks of upper block 
                                 //! tridiagonal matrix
     
      UnrollLoop<0,Q*Q>::assign_rac2rac(U[0],diag[0]);

      

      value_type L_trans[Q*Q];
      char spectrans = 'T'; //needed for Algo (4.5.4). See comment below 
       
      //! Actually we have to solve \f$ U_{i-1}^T L_{i-1}^T = E_{i-1}^T\f$
      //! which requires the transposition of the rhs \f$E_{i-1}\f$ 
      //! Note that we then have to perform the following GEMM operations
      //! involving \f$ L_{i-1}\f$ with its transpose, i.e. the first argument
      //! to GEMM will be 'T' ;)

      for(int i = 1; i < N; ++i){
	
	MatrixUnroller<0,0,Q,Q>::transpose(&L_trans[0],low[i-1]);
	
	//! perform one LU decomposition on D
	F77Type<TypeTraits<value_type>::Value>::GETRF(&q,&q,diag[i-1],&lda,&ipiv[0],&info);
	adonis_assert(info == 0);

	//!solve system, see. algo. (4.5.4)
	F77Type<TypeTraits<value_type>::Value>::GETRS(&spectrans,&q,&q,diag[i-1],&lda,&ipiv[0],&L_trans[0],&ldb,&info); //note that low[i-1] bears solution now


	adonis_assert(info == 0);

	//! update diag 
	F77Type<TypeTraits<value_type>::Value>::GEMM(&spectrans,&trans,&q,&q,&q,&alpha,&L_trans[0],&lda,up[i-1],&ldb,&beta,diag[i],&ldc);

	//! ok, U[i] stores the LU decomposed and updated version of diag[i]
	UnrollLoop<0,Q*Q>::assign_rac2rac(U[i],diag[i]);

	F77Type<TypeTraits<value_type>::Value>::GEMM(&spectrans,&trans,&q,&nrhs,&q,&alpha,&L_trans[0],&lda,rhs[i-1],&ldb2,&beta,rhs[i],&ldc);
      } //end for-loop
	  
      
      //!backward substitution
      //! perform one LU decomposition on D_N, not yet performed
      F77Type<TypeTraits<value_type>::Value>::GETRF(&q,&q,diag[N-1],&lda,&ipiv[0],&info);
      adonis_assert(info == 0);
      //! solve system
      F77Type<TypeTraits<value_type>::Value>::GETRS(&trans,&q,&nrhs,diag[N-1],&lda,&ipiv[0],rhs[N-1],&ldb,&info);
   

      
      //BACK SUBST.:
      //! Solution
      typedef StaticArray<value_type,Q*NRHS> RhsRowType;
      RhsRowType x[N > 0 ? N : 1];
      
      UnrollLoop<0,NRHS*Q>::assign_rac2rac(x[N-1],rhs[N-1]);

      
      for(int i = N-2; i >= 0; --i){

	UnrollLoop<0,NRHS*Q>::assign_rac2rac(x[i],rhs[i]);

	
	F77Type<TypeTraits<value_type>::Value>::GEMM(&trans,&trans,&q,&nrhs,&q,&alpha,up[i],&q,&x[i+1][0],&ldb2,&beta,&x[i][0],&ldc);

	//!each \f$ U_i\f$ must be factored, cf. [1, p. 175] 
	F77Type<TypeTraits<value_type>::Value>::GETRF(&q,&q,&U[i][0],&lda,&ipiv[0],&info);

	F77Type<TypeTraits<value_type>::Value>::GETRS(&trans,&q,&nrhs,&U[i][0],&lda,&ipiv[0],&x[i][0],&ldb,&info);
	adonis_assert(info == 0);

	//assign solution to rhs, i.e. overwrite rhs with solution
	UnrollLoop<0,NRHS*Q>::assign_rac2rac(rhs[i],x[i]);
	
      }//end backwrd. subst.

    }

  };

}//end namespace

#endif
#endif

#endif
