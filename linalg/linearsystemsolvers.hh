#ifndef DENSE_LINEAR_SYSTEM_SOLVERS_HH
#define DENSE_LINEAR_SYSTEM_SOLVERS_HH

#include <cmath>
#include<typeinfo>

#include "../common/adonisassert.hh"

#if USE_LAPACK
#ifdef LAPACK_USAGE

#include "../lapackwrapper/lapackdrivertraits.hh"
#include "../lapackwrapper/wrappertraits.hh"

#include "../misc/useful.hh"
#include "../misc/operations4randomaccesscontainers.hh"

#define CALCULATE_RECIPROCAL_CONDITION_NUMBER 1


namespace Adonis{

  /**
   * \brief Solves general <B>dense</B> SQUARE system, i.e. \f[ A\cdot x = b, \quad A \in \mathbb{R}^{n \times n}\f] via LU decomposition
   * \tparam VType STL-compliant random access container
   * \param A square matrix
   * \param b right hand side vector
   * 
   * CAUTION: it seems that (my) xgesv does not yield the same results like a call of xgetrf, followed by xgetrs (which should be virtually the same)!!! 
   */
  /*  template<class VType>
  inline void square_solve(VType& A, VType& b){
    adonis_assert(b.size()*b.size() == A.size()); //only square matrices allowed
  
    typedef typename VType::value_type value_type;

    int n = (int)b.size(),   //this is the order of the matrix(= size of rhs)
      nrhs = 1,              //number of rhs is always 1 
      lda = n,
      ldb = n,
      info;

    int ipiv[n];
    
    F77Type<TypeTraits<value_type>::Value >().GESV(&n,&nrhs,&A[0],&lda,&ipiv[0],&b[0],&ldb,&info);
    adonis_assert(info == 0);  // info = -i, the i-th argument has illegal value
                               //info = i  U(i,i) = 0 => U sing. => no solution

  }
*/

  /**
   * \brief More general version to solve a system \f$A\cdot X = B, \f$ with \f$ A \f$ being of order n. 
   */
  /*
  template<class V>
  inline void square_solve(V& A, int n, V& B, int nrhs = 1){
    adonis_assert(static_cast<int>(A.size()) == n*n);
    adonis_assert(static_cast<int>(B.size())/nrhs == n);
    
    typedef typename V::value_type value_type;
    
    int lda = std::max(1,n),
      ldb = std::max(1,n),
      info;
    
    int* ipiv = new int[n];
   
    F77Type<TypeTraits<value_type>::Value >().GESV(&n,&nrhs,&A[0],&lda,&ipiv[0],&B[0],&ldb,&info);
    //adonis_assert(info == 0);
    if(info < 0)
      ADONIS_ERROR(LapackError, "the "<<info<<"th argument had an illegal value");
    else if(info > 0)
      ADONIS_ERROR(LapackError,"U("<<info<<","<<info<<") = 0 (F77 indexing, i.e. info-1), hence no solution could be computed, pal.");

    delete[] ipiv;
  }
  */
  

  /**
   * \brief Solves the SQUARE <B>dense</B> \f$ m \times n \f$ system \f[ A\cdot x = b\f], where the rhs b will be overwritten to store the solution.
   * \tparam VType Random Access Container (STL-compliant)
   * \param A general dense matrix stored row-wise as in random access container
   * \param b the 
   * \param nrhs number of rhs (default: 1)
   * \param nb number of block size (default: 1)
   * \param trans A in transposed (default) or normal F77 style
   */
  /*template<class VType>
  inline void solve(VType& A, VType& b, char trans = 'T', int nrhs = 1, int nb = 1){
    
     //previous version
    adonis_assert(A.size()%(b.size()) == 0 );//only meaningful 4 nrhs = 1
    
    adonis_assert(trans == 'N' || trans == 'T' || trans == 'C');

    typedef typename VType::value_type value_type;

    
    int n = (int)A.size()/(b.size()),  //(int)b.size(),    //columns of A
      m = (int)b.size(),//(int)A.size()/(b.size()),  // rows of A
      lda = m,        //>= max(1,m)
      idim = ((m<=n) ? m : n);  //(m >= n) ? m : n;     //min(n,m)
    
    //std::cout<< "LDA = "<<lda<<std::endl<<"dim(IPIV) = "<<idim <<std::endl;
    
    int ipiv[idim];
    
    //int nrhs = 1;                  // number of rhs (here always 1)
    int  info;
    
    
     
    //! performs LU decomposition of \f$ A = P\cdot L \cdot U\f$. 
    //! ipiv stores the indexes of P 
    F77Type<TypeTraits<value_type>::Value>().GETRF(&m,&n,&A[0],&lda,&ipiv[0],&info);
    adonis_assert(info == 0);    // info = -i illegal i-th argument 
                                 // info = i U(i,i) = 0, i.e. U singular
       
    //std::cout << " After getrf: A = "<<A << std::endl; //check output


    
    //!solves \f$ A \cdot x = b \f$ or \f$ A^T \cdot x = b \f$, 
    //!with \f$A \in \mathbb{R}^{n \times n}\f$ using the LU factorisation 
    //! previously computed by GETRF 
    F77Type<TypeTraits<value_type>::Value>().GETRS(&trans,&m,&nrhs,&A[0],&m,&ipiv[0],&b[0],&m,&info);
    adonis_assert(info == 0);    // info = -i illegal i-th argument 
  }
  */



  /**
   * \brief Solves a <B>square</B> linear system \f[ A\cdot X = B \f] via LU decomposition of A.
   * 
   * CAUTION: Lapack already provides the routine 'xgesv' based on, like here, 'xgetrf' followed by 'xgetrs'; however, the former one does not seem to work properly! Hence the prefix 'good_' ;)
   *
   * \tparam V STL-compliant random access container 
   * \param A of order \f$n\f$
   * \param n order of square matrix
   * \param B right hand side matrix
   * \param nrhs number of right hand sides (a.k.a. columns of \f$ B \f$)
   * \param useTransposedForm use \f$ A^T\cdot X = B \f$ or \f$ A^T\cdot H = B \f$ (default: true)
   */
  template<class V>
  inline void good_square_solve(V& A, int n, V& B, int nrhs, bool useTransposedForm =  true){
    adonis_assert(static_cast<int>(A.size()) == n*n); //a must be square
 
    typedef typename V::value_type value_type;
    char trans = 'N';
    
    //!select description -- either \f$ A \f$ or \f$ A^T \f$ (\f$ A^H \f$, when \f$ A \in \mathbb{C}^{m \times n}\f$)
    if(!F77Type<TypeTraits<value_type>::Value>::IsComplex){ //real
      if(useTransposedForm)
	trans = 'T';
    }
    else{                                                   //complex
      if(useTransposedForm)
	trans = 'C';
    }

    
    int lda = std::max(1,n),
      ldb = std::max(1,n),
      info;

    adonis_assert(static_cast<int>(B.size()) == ldb*nrhs);

    int* ipiv = new int [n];

    //! performs LU decomposition of \f$ A = P\cdot L \cdot U\f$. 
    //! ipiv stores the indexes of P 
    F77Type<TypeTraits<value_type>::Value>().GETRF(&n,&n,&A[0],&lda,&ipiv[0],&info);
    adonis_assert(info == 0);

    //!solves \f$ A \cdot x = b \f$ or \f$ A^T \cdot x = b \f$, 
    //!with \f$A \in \mathbb{R}^{n \times n}\f$ using the LU factorisation 
    //! previously computed by GETRF 
    F77Type<TypeTraits<value_type>::Value>().GETRS(&trans,&n,&nrhs,&A[0],&lda,&ipiv[0],&B[0],&ldb,&info);
    
    if(info != 0)
      ADONIS_ERROR(LapackError,"An error occured within Lapack's "<<typeid(value_type).name() << "getrs");

    delete[] ipiv;
  }



  /**
   *  \brief Solve under/overdetermined <B>dense</B> linear systems. 
   * It can, of course, be used for square systems.
   * Check e.g. 
   * <a href="http://www.netlib.org/lapack/double/dgels.f">under(over)determined linear solver based on QR </a>
   *
   * NOTE1: Since ldb = max(1,m,n) the rhs b and the solution x with which b will be overwritten must be of the same size, e.g. \f$ A \in \mathbb{R}^{r \times n}, b \in \mathbb{R}^r\f$, the rhs/solution b must be of size \f$ n\f$. In particular \f$ b = [ \times,\times,\cdot,\cdot,\cdot,\cdot ]^T \f$
   *
   * NOTE2: Recall that F77 style is column-based while C++-style is row-based 
   */
  template<class VType>
  inline void non_square_solve(VType& A, int m, VType& b, char trans = 'T', bool query = false, int nrhs = 1, int nb = 1){
    
    
    adonis_assert(trans == 'N' || trans == 'T');

    typedef typename VType::value_type value_type;

    int n = (int)A.size()/m,        //columns
      lda = m,   //>= max(1,m)
      ldb = (m >= n) ? m : n,   //>= max(1,m,n)
      mn = (m <= n) ? m : n,    //min(m,n),
      maxMNnrhs = (mn >= nrhs) ? mn : nrhs, //>= max(1, mn+max(mn,nrhs)*nb
      lwork = mn + maxMNnrhs*nb,
      info;

    adonis_assert((int)b.size() == ldb*nrhs); //see www.netlib.org
    
    //std::cout << "M = "<<m << std::endl<<"N = "<<n<<std::endl<<"LDA = "<<lda<<std::endl<<"LDB = "<<ldb << std::endl << "mn = "<<mn<<std::endl<<"maxMNnrhs = "<<maxMNnrhs<<std::endl<<"lwork = "<<lwork<<std::endl; //check output

    VType work;
    
    if(query){ //query and allocate the optimal workspace automatically
      lwork = -1;
      value_type wkopt; //lwork is -1
      dgels_(&trans,&m,&n,&nrhs,&A[0],&lda,&b[0],&ldb,&wkopt,&lwork,&info);
      lwork = (int)wkopt;
      //std::cout<<"lwork = "<<lwork<<std::endl;  //check output
      //value_type work[lwork*sizeof(value_type)];
      work.resize(lwork*sizeof(value_type));
    }
    else{
      //value_type work[lwork];
      work.resize(lwork);
    }
    
    //CAUTION: right hand side must be of dimension (ldb,nrhs)!!
    F77Type<TypeTraits<value_type>::Value>().GELS(&trans,&m,&n,&nrhs,&A[0],&lda,&b[0],&ldb,&work[0],&lwork,&info);
    
    if(nrhs > 1) //bring back in usual form
      b = transpose(b,nrhs);
    
    adonis_assert(info == 0);
  }


  /**
   * \brief Solve general under/overdetermined <B>dense</B> linear system. It is assumed that the coefficient matrix has <B>full</B> rank
   * \tparam V random access container which serves as storage container for the matrix
   * \param A coefficient matrix
   * \param m number of rows in A
   * \param B right hand side matrix (or vector in case nrhs = 1)
   * \param nrhs number of right hand sides (= columns of B)
   * \param descr Storage description of A (default: true, i.e. choose \f$A^T (A^H) \f$)
   */
  template<class V>
  inline void solve(V& A, int m, V& B, int nrhs = 1, bool descr = true){
    adonis_assert((int)A.size()%m == 0 && m != 0);
    
    typedef typename V::value_type value_type;
    char trans = 'N';
    
    //!select description -- either \f$ A \f$ or \f$ A^T \f$ (\f$ A^H \f$, when \f$ A \in \mathbb{C}^{m \times n}\f$)
    if(!F77Type<TypeTraits<value_type>::Value>::IsComplex){ //real
      if(descr)
	trans = 'T';
    }
    else{                                                   //complex
      if(descr)
	trans = 'C';
    }


    int n = (int)A.size()/m,  //number of columns
      lda = std::max(1,m),
      ldb = std::max(m,n),
      lwork = -1,             //automatic workspace detection
      info;
  
    adonis_assert((int)B.size() == ldb*nrhs);

    value_type wkopt;
    //calculate optimal lwork -- nothing else
    F77Type<TypeTraits<value_type>::Value>::GELS(&trans,&m,&n,&nrhs,&A[0],&lda,&B[0],&ldb,&wkopt,&lwork,&info);

    //take only the real part in case you calculate with complex values
    lwork = LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::real_part(wkopt);

    value_type* work = new value_type [lwork*sizeof(value_type)];
    
    //solve now general system 
    F77Type<TypeTraits<value_type>::Value>::GELS(&trans,&m,&n,&nrhs,&A[0],&lda,&B[0],&ldb,&work[0],&lwork,&info);

     delete[] work;
     
     if(info < 0)
       ADONIS_ERROR(LapackError,"The "<<info<<"th has an illegal value.");
     else if(info > 0)
       ADONIS_ERROR(LapackError,"The "<<info<<"th diagonal element of the triangular factor of A is zero, \n   so that A does not have full rank; LS solution cannot be computed, pal.");
     
     //if(nrhs > 1)     //transform back into usual form
     //B = transpose(B,nrhs);

  }

  
  /**
   * \brief computes min.-norm solution to linear least squares problem: 
   *           \f$ \min_x \|Ax-b\| \f$
   *  using an SVD of A (any dense (possibly) rank-deficient matrix)
   *
   * \tparam V container type
   * \param A LHS matrix (may be rank deficient)
   * \param m number of rows of A
   * \param B RHS matrix (input). On output: n x nrhs solution of A·X = B 
   * \param nrhs number of columns of B (and later the solution stored in B)
   */
  template<class V>
  inline V solve_rank_deficient_ls(const V& A, int m,  const V& B, int nrhs){
    adonis_assert(static_cast<int>(B.size()/nrhs) == m);
    adonis_assert((static_cast<int>(A.size()%m == 0)));
    
    
    typedef typename V::value_type value_type;
    typedef typename F77Type<TypeTraits<value_type>::Value>::BaseType BaseType;

    int n = static_cast<int>(A.size()/m),
      lda = std::max(1,m),
      ldb = std::max(1,std::max(m,n)),
      minmn = std::min(m,n),
      smlsiz = 25,
      lwork = -1,
      log2term = (int)(log2(static_cast<BaseType>(minmn)/(smlsiz+1))) + 1, 
      nlvl = std::max(0,log2term),
      liwork = std::max(1,3*minmn*nlvl + 11*minmn),
      rank,
      info;
    
    
    //! actually, this is the dimension of the overall solution array from
    //! which on has to extract the solution 
    V b(ldb*nrhs),   //lda*nrhs
      a = transpose(A,n);
    

    for(int i = 0; i < m; ++i) 
      for(int j = 0; j < nrhs; ++j)
	b[ColumnMajor::offset(i,j,n)] = B[RowMajor::offset(i,j,nrhs)];
    
    
    //std::cout << "solve_rank_deficient_ls: b = "<< b << std::endl;
    


    BaseType rcond = -1;
     
    BaseType* S = new BaseType [minmn]; //the sv in decreasing order
    int* iwork = new int [liwork];
     
    //complex stuff
    int lrwork;
    
    LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::lwk(lrwork,m,n,smlsiz,nlvl,nrhs); //void fct that assigns lrwork a value
    
    BaseType* rwork = new BaseType [std::max(1,lrwork)]; 

    value_type wkopt;  // stores optimal work space later on
    //determine optimal lwork size automatically -- only this is computed
    F77Type<TypeTraits<value_type>::Value>::GELSD(&m,&n,&nrhs,&a[0],&lda,&b[0],&ldb,&S[0],&rcond,&rank,&wkopt,&lwork,&rwork[0],&iwork[0],&info);
   
    lwork = LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::real_part(wkopt);  
    lwork = std::max(1,lwork); //(1 > lwork) ? 1 : lwork;   
  
    value_type* work = new value_type [lwork*sizeof(value_type)]; 
    
    F77Type<TypeTraits<value_type>::Value>::GELSD(&m,&n,&nrhs,&a[0],&lda,&b[0],&ldb,&S[0],&rcond,&rank,&work[0],&lwork,&rwork[0],&iwork[0],&info);
  
#if SHOW_ME_SOME_OUTPUT
  std::cout << "EFFECTIVE RANK = "<<rank<< " (of "<< m<<"x"<<n<<" lhs matrix)"<<std::endl;
  std::cout<< "SV's computed within xgelsd_ (in decreasing order) = " << std::endl; for(int i = 0; i < minmn; ++i) std::cout << S[i] << "  "; std::cout << std::endl;
#endif 

    if(info < 0)
      ADONIS_ERROR(LapackError, "The "<<Abs(info)<<"th argument had an illegal value.");
    if(info > 0)
      ADONIS_ERROR(LapackError, "SVD algorithm diverged: "<<info<<"th off-diagonal elements of an intermediate bidiagonal form did not converge to zero.");
  
    //delete the by 'new' defined quantities
    delete[] rwork;
    delete[] work;
    delete[] S;
    delete[] iwork;
    
    if(nrhs > 1)
      b = transpose(b,n);

    //apparently, the n*nrhs part contains the solution in the m > n case
    if(m > n){
      V sol(n*nrhs);
      for(int i = 0; i < n; ++i)
	for(int j = 0; j < nrhs; ++j)
	  sol[RowMajor::offset(i,j,nrhs)] = b[RowMajor::offset(i,j,nrhs)];
      b = sol;
    }

    return b;
  }
  


  /**
   * \brief Although the expicit inversion of a matrix is often unnecessary due to the computational costs (especially in the case of solving square linear systems), there may be circumstances when an inverse is actually needed and can't be circumvented.
   *
   * The inverse is calculated via LU decomposition, using Lapack's xgetri after xgetrf, cf. [HIGHAM, §14, pp. 267]
   * \tparam V STL-compliant random access container 
   * \param A square matrix
   * \param n order of A
   */
  template<class V>
  inline V inverse(const V& A, int n){
    adonis_assert(n*n == static_cast<int>(A.size())); //A must be square
    
    typedef typename V::value_type value_type;
    
    V Inv(A);  //copy A 
   
    if(n == 1){
      adonis_assert(Abs(A[0]) > value_type());
      Inv[0] = static_cast<value_type>(1)/A[0];
      return Inv;
    }

    if(n == 2){
      value_type recdetA = 1./(A[0]*A[3] - A[2]*A[1]);
      Inv[0] = recdetA*A[3];  Inv[1] = -recdetA*A[1]; 
      Inv[2] = -recdetA*A[2]; Inv[3] = recdetA*A[0];
      return Inv;
    }
    
    
    int lda = n,
      info;

    int* ipiv = new int[n];

    //LU decomposition of A
    F77Type<TypeTraits<value_type>::Value>().GETRF(&n,&n,&Inv[0],&lda,&ipiv[0],&info);
    adonis_assert(info == 0);
    int lwork = n;
    value_type* work = new value_type [lwork];
    
    //calulate inverse of A based on the previous LU decomposition
    F77Type<TypeTraits<value_type>::Value>().GETRI(&n,&Inv[0],&lda,&ipiv[0],&work[0],&lwork,&info);

    if(info < 0)
      ADONIS_ERROR(LapackError,"The " << info <<"th argument has an illegal value.");
	
    if(info > 0)
      ADONIS_ERROR(LapackError,"U(" << info <<", "<<info<<") = "<<0<<" (in C++: U["<<info-1<<"]["<<info-1<<"] = 0). Matrix is singular ==>\n   No Inverse exists!");

    delete[] work;
    delete[] ipiv;
    

    return Inv;
  }


  /**
   *\brief Overload inverse -- no local object. May be beneficial in some cases
   */
  template<class V>
  inline V& inverse(V& Inv, const V& A, int n){
    adonis_assert(n*n == static_cast<int>(A.size())); //A must be square

    typedef typename V::value_type value_type;

    Inv = A; //the inverse; copying to keep A
   
    if(n == 1){
      adonis_assert(Abs(A[0]) > value_type());
      Inv[0] = static_cast<value_type>(1)/A[0];
    }
    
    if(n == 2){
      value_type recdetA = 1./(A[0]*A[3] - A[2]*A[1]);
      Inv[0] = recdetA*A[3];  Inv[1] = -recdetA*A[1]; 
      Inv[2] = -recdetA*A[2]; Inv[3] = recdetA*A[0];
      return Inv;
    }

    int lda = n,
      info;

    int* ipiv = new int[n];

    //LU decomposition of A
    F77Type<TypeTraits<value_type>::Value>().GETRF(&n,&n,&Inv[0],&lda,&ipiv[0],&info);
    adonis_assert(info == 0);
    int lwork = n;
    value_type* work = new value_type [lwork];
    
    //calulate inverse of A based on the previous LU decomposition
    F77Type<TypeTraits<value_type>::Value>().GETRI(&n,&Inv[0],&lda,&ipiv[0],&work[0],&lwork,&info);

    if(info < 0)
      ADONIS_ERROR(LapackError,"The " << info <<"th argument has an illegal value.");
	
    if(info > 0)
      ADONIS_ERROR(LapackError,"U(" << info <<", "<<info<<") = "<<0<<" (in C++: U["<<info-1<<"]["<<info-1<<"] = 0). Matrix is singular ==>\n   No Inverse exists!");

    delete[] work;
    delete[] ipiv;

    return Inv;
  }

  

  /**
   * \brief Solve the system \f$ A\cdot X = B, \f$ where \f$ A\f$ is symmetric or hermitian (self-adjoint), i.e. \f$ A = A^T\f$ or in the complex case \f$ A = A^*\f$.
   */
  template<class V>
  inline void solve_symmetric_ls(V& A, int n, V& B, int nrhs, char uplo = 'U', bool swapsol = false){
    adonis_assert(static_cast<int>(A.size())/n == n && static_cast<int>(B.size())/nrhs == n);
    adonis_assert(uplo == 'U' || uplo == 'L');//store upper or lower triangle

    typedef typename V::value_type value_type;

    int lda = std::max(1,n),
      ldb = std::max(1,n),
      lwork = -1,
      info;

    int* ipiv = new int [n];
    value_type wkopt;
    
    //calculate optimal lwork -- nothing else
    F77Type<TypeTraits<value_type>::Value>().SYSV(&uplo,&n,&nrhs,&A[0],&lda,&ipiv[0],&B[0],&ldb,&wkopt,&lwork,&info);

    //take only the real part in case you calculate with complex values
    lwork = LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::real_part(wkopt); 

    value_type* work = new value_type [lwork*sizeof(value_type)];
    
    F77Type<TypeTraits<value_type>::Value>().SYSV(&uplo,&n,&nrhs,&A[0],&lda,&ipiv[0],&B[0],&ldb,&work[0],&lwork,&info);

    //proper clean-up
    delete[] work;
    delete[] ipiv;

    if(info < 0)
      ADONIS_ERROR(LapackError, "The "<<info<<"th argument had an illegal value.");
    if(info > 0)
      ADONIS_ERROR(LapackError, "D("<<info<<","<<info<<") = 0 ==> block diagonal matrix D is singular, so the solution couldn't be computed, pal.");
  
    if(swapsol){  //stores column-by-column:do you want to transpose it ?
      if(nrhs > 1) 
	B = transpose(B,n); 
    }
  }

  
  /**
   * \brief Solves positive definite system.
   *
   * NOTE: for questions concerning swapsol, cf. 'documentation of solve_symmetric_ls'
   */
  template<class V>
  inline void solve_positive_definite_ls(V& A, int n, V& B, int nrhs, char uplo = 'U', bool swapsol = false){
    adonis_assert(static_cast<int>(A.size())/n == n && static_cast<int>(B.size())/nrhs == n);
    adonis_assert(uplo == 'U' || uplo == 'L');

    typedef typename V::value_type value_type;

    int lda = std::max(1,n),
      ldb = std::max(1,n),
      info;
    
    //solve system
    F77Type<TypeTraits<value_type>::Value>().POSV(&uplo,&n,&nrhs,&A[0],&lda,&B[0],&ldb,&info);
    

    if(info < 0) 
      ADONIS_ERROR(LapackError, "The "<<info<<"th argument has an illegal value." );

    if(info > 0)
      ADONIS_ERROR(LapackError, "Leading minor of order "<<info<<" of A is not p.d. \n   ==> No Cholesky factorization computed ==> no solution.");

    
    if(swapsol){  //stores column-by-column:do you want to transpose it ?
      if(nrhs > 1) 
	B = transpose(B,n); 
    }
  }


  /**
   * \brief Computes eigenvalues of symmetric/hermitian matrix
   * 
   * Overload it to adapt it for calculating eigenvectors as well
   */
  template<class V, class W>
  inline void eig_symmetric_hermitian(V& A, int n, W& ew, char jobz = 'N', char uplo = 'U'){ //'N' compute eigenvalues only
    adonis_assert(static_cast<int>(ew.size()) == n );
    adonis_assert(jobz == 'N' || jobz == 'V');
    adonis_assert(uplo == 'U' || uplo == 'L');
    adonis_assert(static_cast<int>(A.size())/n == n); //must be square
   
    typedef typename V::value_type value_type; //might be complex as well
    typedef typename TypeTraits<value_type>::BaseType BaseType;

    adonis_assert(typeid(typename W::value_type) == typeid(BaseType));
    
   
    int lda = std::max(1,n),
      lwork = -1,           //calculate optimal workspace automatically
      info;
    
    //BaseType* w = new BaseType [n]; 
    
    value_type wkopt;
      
    BaseType* rwork = new BaseType [LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::lwkev(n)];  //if complex then max(1,3*n-2), else 0

    //calculate optimal lwork -- nothing else
    F77Type<TypeTraits<value_type>::Value>().SYEV(&jobz,&uplo,&n,&A[0],&lda,&ew[0],&wkopt,&lwork,&rwork[0],&info);
    
    lwork = LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::real_part(wkopt); 
    
    value_type* work = new value_type [lwork*sizeof(value_type)];
    
    //determine eigenvalues (and eigenvectors) of hermitian system
    F77Type<TypeTraits<value_type>::Value>().SYEV(&jobz,&uplo,&n,&A[0],&lda,&ew[0],&work[0],&lwork,&rwork[0],&info);
   
    //delete[] w;
    delete[] rwork;
    delete[] work;

    if(info < 0) 
      ADONIS_ERROR(LapackError, "The "<<info << "th argument had an illegal value");

    if(info > 0)
      ADONIS_ERROR(LapackError, "The algorithm failed to converge. \n   "<<info<<" off-diagonal elements of an intermediate tridiagonal form did not converge to zero.");
 
  }

  /**
   * \brief eigenvalues of a general non-symmetric matrix
   *  
   * Overload it for calculation of eigenvectors
   *
   * Note: eigenvalues are stored consecutively with increasing real parts. 
   *       Conjungate complex eigenvalues are stored with the positive 
   *       imaginary part first.
   */
   template<class V, class C>
   inline void eig_general(V& A, int n, C& ew, char jobvl ='N', char jobvr = 'N'){
     //typedef typename V::value_type eigenvector_type;
  
     typedef typename V::value_type value_type; //might be complex as well
     typedef typename TypeTraits<value_type>::BaseType BaseType;
     
     typedef WithEigenvectorComputation<V,false> NoEigenVectors;

     adonis_assert(jobvl == 'N' || jobvl == 'V');
     adonis_assert(jobvr == 'N' || jobvr == 'V');
     
     int lda = std::max(1,n),
       ldvl = n,
       ldvr = n, 
       lwork = -1, //automatic dimension computation
       info;

     value_type* wr = new value_type[LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::dim(n)];
     value_type* wi = new value_type[LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::dim(n)];
     
     value_type* vl = new value_type [NoEigenVectors::evecsize(ldvl,n)];
     value_type* vr = new value_type [NoEigenVectors::evecsize(ldvr,n)];
     
     value_type wkopt;
      
     BaseType* rwork = new BaseType [LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::lwk(n)];  //2*n if complex otherwise 0 

     //!calculate optimal lwork -- nothing else
     F77Type<TypeTraits<value_type>::Value>().GEEV(&jobvl,&jobvr,&n,&A[0],&lda,&wr[0],&wi[0],&ew[0],&vl[0],&ldvl,&vr[0],&ldvr,&wkopt,&lwork,&rwork[0],&info);
     //dgeev_(&jobvr, &jobvr, &n, &A[0], &lda, &wr[0], &wi[0], &vl[0], &ldvl, &vr[0], &ldvr,&wkopt, &lwork, &info );

     lwork = LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::real_part(wkopt); 
    
     value_type* work = new value_type [lwork*sizeof(value_type)];

     //!calculate eigenvalues
     F77Type<TypeTraits<value_type>::Value>().GEEV(&jobvl,&jobvr,&n,&A[0],&lda,&wr[0],&wi[0],&ew[0],&vl[0],&ldvl,&vr[0],&ldvr,&work[0],&lwork,&rwork[0],&info);
     
   
     //!output
     //for(int i = 0; i < n; ++i){  std::cout << "("<<wr[i] << ", "<<wi[i] <<")"<<std::endl;}

     LapackComplex<F77Type<TypeTraits<value_type>::Value>::IsComplex>::assign_eigenvalues(ew,&wr[0],&wi[0]);
 
      //proper clean up
      delete[] wr;
      delete[] wi;
      delete[] vl;
      delete[] vr;
      delete[] rwork;
      delete[] work;
   
      if(info < 0)
	ADONIS_ERROR(LapackError, "The "<<info << "th argument had an illegal value");
      if(info > 0)
	ADONIS_ERROR(LapackError, "the QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed");
   }


} //end namespace 

#endif
#endif

#endif 
