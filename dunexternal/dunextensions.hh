#ifndef DUNE_EXTENSIONS_HH
#define DUNE_EXTENSIONS_HH

#include<limits>
#include<typeinfo>

#if USE_LAPACK
#ifdef LAPACK_USAGE
#include "../lapackwrapper/lapackdrivertraits.hh"
#include "../lapackwrapper/wrappertraits.hh"
#endif
#endif

#include "../dunestuff/fmatrix.hh"
#include "../common/error.hh"
#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../common/numerictypechecker.hh"
#include "../common/typeadapter.hh"
#include "../common/smartassign.hh"
#include "../common/typeselector.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../misc/misctmps.hh"
#include "../misc/useful.hh"


#include "dunecontainercopy.hh"
#include "normalspace.hh"

namespace Adonis{

  /**
   * \brief fill diagonal matrix with entries from vector w
   * a similar function has been written to work with a random access container
   * instead of a field matrix
   */
  template<class V, int N>
  inline Dune::FieldMatrix<typename V::value_type,N,N>& diagonal_matrix(Dune::FieldMatrix<typename V::value_type,N,N>&  diag, const V& w){
    typedef typename V::value_type value_type;
    adonis_assert((int)w.size() == N);
    
    for(int i = 0; i < N; ++i){
      for(int j = 0; j < N; ++j){
	if(i==j)
	  diag[i][j] = w[i];
	else
	  diag[i][j] = value_type();
      }
    }
    return diag;
  }

 
  /**
   * \brief Concatenate FieldMatrix and matrix based on random access vector
   */
  template<int M, int N, class V>
  class ConcatenateMatrix{
  public:
    typedef V VType;
    typedef typename V::value_type value_type;
    typedef Dune::FieldMatrix<value_type,M,N> FieldMtxType;

    //! m = rows, n = cols
    ConcatenateMatrix(){}

    //! w is tought to be stored row-wise (C-style)
    template<class W>
    V& columnwise(int m, int n, const W& w, const FieldMtxType& ftx){
 
	check_size(m,n,w);
	adonis_assert(n == N);

	int r = m+M;
	//! only resize when necessary. Otherwise just overwrite entries
	((int)data_.size() != r*n) ? data_.resize(r*n) : do_nothing();
	
	for(int i = 0; i < m; ++i)
	  for(int j = 0; j < n; ++j)
	    data_[RowMajor::offset(i,j,n)] = w[RowMajor::offset(i,j,n)];

	for(int i = 0; i < M; ++i){
	  for(int j = 0; j < N; ++j){
	    data_[RowMajor::offset(i+m,j,n)] = ftx[i][j];
	  }
	}      
	
      
      return data_;
    }


    //fieldmatrix first
     template<class W>
     V& columnwise(const FieldMtxType& ftx, int m, int n, const W& w){
    
	check_size(m,n,w);
	adonis_assert(n == N);
      
	int r = m+M;
	((int)data_.size() != r*n) ? data_.resize(r*n) : do_nothing();
	
	for(int i = 0; i < M; ++i)
	  for(int j = 0; j < N; ++j)
	    data_[RowMajor::offset(i,j,N)] = ftx[i][j];

	for(int i = 0; i < m; ++i){
	  for(int j = 0; j < n; ++j){
	    data_[RowMajor::offset(i+M,j,n)] = w[RowMajor::offset(i,j,n)];
	  }
	}      

      return data_;
    }

   template<class W>
   V& rowwise(int m, int n, const W& w, const FieldMtxType& ftx){
     
       check_size(m,n,w);
       adonis_assert(m==M);
      
       int c = n+N;
       ((int)data_.size() != m*c) ? data_.resize(m*c): do_nothing();
       
       
       for(int i = 0; i < m; ++i){
	 for(int j = 0; j < c; ++j){
	   if(j < n)
	     data_[RowMajor::offset(i,j,c)] = w[RowMajor::offset(i,j,n)];
	   else
	     data_[RowMajor::offset(i,j,c)] = ftx[i][j-n];
	 }
       }
       
      return data_;
   }
    

    template<class W>
    V& rowwise(const FieldMtxType& ftx,int m, int n, const W& w){
    
       check_size(m,n,w);
       adonis_assert(m==M);
      
       int c = n+N;
       ((int)data_.size() != m*c) ? data_.resize(m*c) : do_nothing();
       
       
       for(int i = 0; i < M; ++i){
   	 for(int j = 0; j < c; ++j){
   	   if(j < N)
   	     data_[RowMajor::offset(i,j,c)] = ftx[i][j]; 
   	   else
   	     data_[RowMajor::offset(i,j,c)] = w[RowMajor::offset(i,j-N,n)];
   	 }
       }
      
      return data_;
   }



    const V& data() const {return data_;}

    //void reset_flag(){isConcat_=false;}

    const std::size_t size() const {return data_.size();}

  private:
    V data_;

    // void clear_me(){
    //   if(data_.size() != 0 )
    // 	data_.clear(); //erase(data_.begin(),data_.begin()+data_.size());
    // }

    template<class W>
    void check_size(int m, int n, const W& v){
      if((int)v.size() != m*n)
	ADONIS_ERROR(DimensionError,"argument v is not of dimension m*n.");
    }
  };

 
  template<int M, int N, class V>
  inline Dune::FieldMatrix<typename V::value_type,M,N> diagonal_matrix_multiplication(const V& diag, const Dune::FieldMatrix<typename V::value_type,M,N>& mtx){
    adonis_assert(static_cast<int>(diag.size()) == M);
    Dune::FieldMatrix<typename V::value_type,M,N> Prod;

    for(int i = 0; i < M; ++i){
      for(int j = 0; j < N; ++j){
	Prod[i][j] = diag[i]*mtx[i][j];
      }
    }
    return Prod;
  }

   template<int M, int N, class V>
   inline Dune::FieldMatrix<typename V::value_type,M,N> diagonal_matrix_multiplication(const Dune::FieldMatrix<typename V::value_type,M,N>& mtx, const V& diag){
    adonis_assert(static_cast<int>(diag.size()) == N);
    Dune::FieldMatrix<typename V::value_type,M,N> Prod;

    for(int i = 0; i < M; ++i){
      for(int j = 0; j < N; ++j){
	Prod[i][j] = diag[j]*mtx[i][j];
      }
    }
    return Prod;
  }


  template<int M, int N, class V>
  inline Dune::FieldMatrix<typename V::value_type,M,N> dyadic_product(const V& v, const V& w){
    adonis_assert(static_cast<int>(v.size()) == M && static_cast<int>(w.size()) == N);
    typedef typename V::value_type value_type;
    Dune::FieldMatrix<value_type,M,N> tmp;
    for(int i = 0; i < M; ++i){
      for(int j = 0; j < N; ++j){
	tmp[i][j] = v[i]*w[j];//row major by default
      }
    }
    return tmp;
  }

  template<int M, int N, class V>
  inline Dune::FieldMatrix<typename V::value_type,M,N> fill_from_rac(const V& v){
    adonis_assert(static_cast<int>(v.size()) == M*N);
    
    typedef typename V::value_type value_type;
    Dune::FieldMatrix<value_type,M,N> tmp;

    //row-major filling
    for(int i = 0; i < M; ++i)
      for(int  j= 0; j < N; ++j)
	tmp[i][j] = v[RowMajor::offset(i,j,N)];
  
    return tmp;
  }

   template<int M, int N, class V>
   inline Dune::FieldMatrix<typename V::value_type,M,N>& fill_from_rac(Dune::FieldMatrix<typename V::value_type,M,N>& mtx, const V& v){
    adonis_assert(static_cast<int>(v.size()) == M*N);
   
    //row-major filling
    for(int i = 0; i < M; ++i)
      for(int  j= 0; j < N; ++j)
	mtx[i][j] = v[RowMajor::offset(i,j,N)];
  
    return mtx;
  }

  /**
   * \brief FieldMatrix multiplication 
   */
  template<class K,int L, int M, int N>
  inline Dune::FieldMatrix<K,L,N> operator*(const Dune::FieldMatrix<K,L,M>& A,
					    const Dune::FieldMatrix<K,M,N>& B){
    
    Dune::FieldMatrix<K,L,N> C; 
    C = K();             //must be initialised because of += !

    for(int i=0; i < L; ++i){
      for(int j=0; j < N; ++j){ 
	for(int k=0; k < M; ++k){
	  C[i][j] += A[i][k]*B[k][j];
	}
      }
    }  
    return C;
  }


 
  /**
   * \brief <B> Right multiplication </B>  of a matrix, stored as a STL-compliant container with a Dune FieldMatrix.
   * Note: some template arguments must be provided at compile time anyway!! 
   * 
   * \tparam RC Storage organisation: either RowMajor or ColumnMajor (must be specified by you in any case)
   * \tparam L,M,N dimension of matrices, i.e. (L x M) and (M x N) (L, i.e. the <I>row dimension of left matrix</I> 'a' must be specified by you in any case) 
   
   * USAGE:
   * \code 
   double mi[] = { 0.5, -1,   2.5,
		    -1,   0.,  -3.5,
		    2.,   1.,   2.,
		    4.5,  -0.5, -4};

    MyVec<double> kaitain(mi,mi+12);
    cout << "kaitain = "<<kaitain <<endl;

    double b[3][5] = { {-3, 0.5, -0.5, -6.25, 3.5},
		       {0.75, 2.25, -2.35, 0.235, -1},
		       {-1.25, -5, 3.45, 0.55, -1.75 } };
    Dune::FieldMatrix<double,3,5> B;
    for(int i= 0; i < 3; ++i)
      for(int j = 0; j < 5; ++j)
	B[i][j] = b[i][j];
   
    cout << "ROW-storage (C-style) of left matrix:"<<endl;
    cout << "kaitain·B = "<< endl<< matrix_matrix_multiplication<RowMajor,4>(kaitain,B)<<endl;
    
    cout << "COLUMN-storage (F77-style) of left matrix:"<<endl;
    MyVec<double> niatiak = transpose(kaitain,3);
    cout << endl << "niatiak " << niatiak <<endl;
    cout << "niatnak·B = "<< endl<< matrix_matrix_multiplication<ColumnMajor,4>(niatiak,B)<<endl;
   * \endcode
   */
  template<class RC, int L, int M, int N, class V>
  inline Dune::FieldMatrix<typename V::value_type,L,N> matrix_matrix_multiplication(const V& a, const Dune::FieldMatrix<typename V::value_type,M,N>& B){
    adonis_assert(int(a.size()) == L*M);
    adonis_assert(L >= 0 && M >= 0 && N >= 0);
    typedef typename V::value_type value_type; 

    Dune::FieldMatrix<value_type,L,N> C; 
    C = value_type();             //must be initialised because of += !
    
    for(int i=0; i < L; ++i){
      for(int j=0; j < N; ++j){ 
	for(int k=0; k < M; ++k){
	  C[i][j] += a[RC::offset(i,k,RC::proper_dim(L,M))]*B[k][j];
	}
      }
    }  
    return C;
  }


  /**
   * \brief container + field container. Return fieldvector of same size
   */
  template<class T, int N, class V>
  inline Dune::FieldVector<T,N> vector_vector_addition(const V& v, const Dune::FieldVector<T,N>& fv){
    adonis_assert((int)v.size() == N);
    Dune::FieldVector<T,N> add;
    
    for(int i = 0; i < N; ++i)
      smart_assign(add[i], v[i] + fv[i]);

    return add;
      
  }
  
  /**
   * \brief TMP for <B> left multiplication </B> of a Dune-matrix with any random access container of type <TT> W </TT> providing an []-operator and assignment to any (other) random access container of type <TT> V </TT>
   */
  template<class T, int M, int N, class V, class W, int I>
  class MatrixVectorProduct{
  private:
    enum{index_ = (I+1) != M};

  public:
    typedef Dune::FieldMatrix<T,M,N> MatrixType;
    
    static inline void perform_operation(V& res, const MatrixType& mat, const W& w){
      smart_assign(res[I],dot_product<N>(mat[I],w));//res[I] = dot_product<N>(mat[I],w);

      MatrixVectorProduct<T,M,N,V,W,(index_ ? (I+1) : M)>::perform_operation(res,mat,w);
    }
  };
  
  //!end of unrolling loop
  template<class T, int M, int N, class V, class W>
  class MatrixVectorProduct<T,M,N,V,W,M>{
  public:
    typedef Dune::FieldMatrix<T,M,N> MatrixType;

    static inline void perform_operation(V& res, const MatrixType& mat, const W& w){} //!nothing to be done here since this would lead to segmentation faults
  };
  
  //!convenient function 
  template<class T, int M, int N, class V, class W>
  inline void matrix_vector_product(V& res, const Dune::FieldMatrix<T,M,N>& mat, const W& w){
    adonis_assert(int(container_size(res)) == M && N == int(container_size(w)));
    //std::cout << "TMP program used here."<<std::endl;

    //! always start from 0
    MatrixVectorProduct<T,M,N,V,W,0>::perform_operation(res,mat,w);
  }



  /**
   * \brief FieldMatrix addition
   */
  template<class K, int M, int N>
  inline Dune::FieldMatrix<K,M,N> operator+(const Dune::FieldMatrix<K,M,N>& A,
					    const Dune::FieldMatrix<K,M,N>& B){
    Dune::FieldMatrix<K,M,N> C(A);
    C += B;             //invoke FieldMatrix +=
    return C;
  }

   /**
   * \brief FieldMatrix subtraction
   */
  template<class K, int M, int N>
  inline Dune::FieldMatrix<K,M,N> operator-(const Dune::FieldMatrix<K,M,N>& A,
					    const Dune::FieldMatrix<K,M,N>& B){
    Dune::FieldMatrix<K,M,N> C(A);
    C -= B;
    return C;
  }



  /**
   * \brief FieldMatrix-FieldVector multiplication, i.e. res = A*v. In my opinion it is more intuitive to use the *-operator instead of functions such as 'mult' provided with fmatrix.hh
   */
  template<class K, int M, int N>
  inline Dune::FieldVector<K,M> operator*(const Dune::FieldMatrix<K,M,N>& A, 
					  const Dune::FieldVector<K,N>& v){
    Dune::FieldVector<K,M> res;

    for(int i = 0; i < M; ++i){
      res[i] = K();
      for(int j = 0; j < N; ++j){
	res[i] += A[i][j]*v[j] ;
	
      }
    }
    
    return res;
  }


  /**
   * \brief FieldVector-FieldMatrix multiplication, i.e. res = v*A 
   */
  template<class K, int M, int N>
  inline Dune::FieldVector<K,N> operator*(const Dune::FieldVector<K,M>& v, 
					  const Dune::FieldMatrix<K,M,N> & A){
    Dune::FieldVector<K,N> res;

    for(int j = 0; j < N; ++j){
      res[j] = K();
      for(int i = 0; i < M; ++i){
	res[j] += v[i]*A[i][j];
      }
    }
    
    return res;
  }

  
  /**
   * \brief construct matrix filled with random numbers 
   */
  template<class K, int M, int N>
  inline void random_matrix(Dune::FieldMatrix<K,M,N>& m){
    for(int i = 0; i < M; ++i){
      for(int j = 0; j < N; ++j){
	
	m[i][j] =  static_cast<K>(rand())/RAND_MAX;
      }
    }
  }


  /**
   *\brief Fill 1x M FieldMatrix from STL-compliant random access container
   */
  template<class V, int M>
  inline void fill_matrix(const V& v, Dune::FieldMatrix<typename V::value_type,1,M>& A){
    adonis_assert((int)v.size() == M);

    for(int j = 0; j < M; ++j)
      A[0][j] = v[j];
  }

  /**
   *\brief Fill Mx1 FieldMatrix from STL-compliant random access container
   */
  template<class V, int M>
  inline void fill_matrix(const V& v, Dune::FieldMatrix<typename V::value_type,M,1>& A){
    adonis_assert((int)v.size() == M); 

    for(int i = 0; i < M; ++i)
      A[i][0] = v[i];
  }

   

  /**
   * \brief Special (left) Multiplication with \f$ B^T \f$ 
   */
  template<class K, int L, int M, int N, class IV>
  inline Dune::FieldMatrix<K,L,N> BT_left_multiplication(const IV& index, const Dune::FieldMatrix<K,L,M>& BT, const Dune::FieldMatrix<K,M,N>& a){
    
    int dim = (int)std::distance(index.begin(),index.end());
    adonis_assert(dim = L);
   
    Dune::FieldMatrix<K,L,N> C;
    for(int i = 0; i < dim; ++i){
      for(int j = 0; j < N; ++j){
	C[i][j] = BT[i][index[i]]*a[index[i]][j];
      }
    }
    return C;
  }


  /**
   * \brief Calculate \f$ B^T + H^T \f$.
   * NOTE: H^T will be changed! It's o.k. since we use it only once per iteration
   */
  template<class M, class IV>
  inline void BT_addition(const IV& index, const M& BT, M& HT){
    size_t dim = std::distance(index.begin(), index.end()),
      k = 0; 
    
    adonis_assert((int)dim == BT.rows);

    for(size_t i = 0; i < dim; ++i){
      HT[k][index[i]] += BT[k][index[i]];
      k++;
    }
  }

  
  /**
   * \brief Calculate \f$ B^T + H^T \f$.
   * NOTE: H^T will be changed! It's o.k. since we use it only once per iteration
   */
  template<class M, class IV>
  inline M& BT_plus_HT(const IV& index, const M& BT, M& HT){
    size_t dim = std::distance(index.begin(), index.end()),
      k = 0; 
    
    adonis_assert((int)dim == BT.rows);

    for(size_t i = 0; i < dim; ++i){
      HT[k][index[i]] += BT[k][index[i]];
      k++;
    }
    return HT;
  }

  /**
   * \brief Transform any STL -compliant random access container into a FieldMatrix
   * \param v random access container storing matrix
   * \param m FieldMatrix which is to be filled with v
   * \param r char denoting 'r'ow-wise (default) or 'c'olumn-wise filling
   */
  template<class V, int M, int N>
  inline void transform_2_matrix(const V& v, Dune::FieldMatrix<typename V::value_type,M,N>& m, char r = 'r'){
    
    size_t k = 0; 
    switch(r){
    case 'r':                          //row-wise storage of matrix in v
      for(int i = 0; i < M; ++i){
	for(int j = 0; j < N; ++j){
	  m[i][j] = v[k++];
	}
      }
      break;
    case 'c':                         //column-wise storage of matrix in v 
      for(int j = 0; j < N; ++j){
	for(int i = 0; i < M; ++i){
	  m[i][j] = v[k++];
	}
      }
      break;
    default:
      ADONIS_ERROR(DerivedError, " r == "<<r<<" is illegal (either 'r' or 'c')");
    }
  }
  

  /**
   * \brief Overload transform_2_matrix to work with iterators
   */
  template<class T, int M, int N, class ITER>
  inline void transform_2_matrix(ITER i1, ITER i2, Dune::FieldMatrix<T,M,N>& m, char r = 'r'){
   
     if (r == 'r'){                       //row-wise storage 
       size_t i = 0,
	 j = 0;
       for(ITER it = i1; it != i2; ++it){
	 m[i][j] = *it;
	 if(j < N-1){
	   ++j;
	 }
	 else{
	   j = 0;
	   ++i;
	 }
       }
     }
     else if(r == 'c'){                  //column-wise storage 
       size_t i = 0,
	 j = 0;
       for(ITER it = i1; it != i2; ++it){
	 m[i][j] = *it;
	 if(i < M-1){
	   ++i;
	 }
	 else{
	   i = 0;
	   ++j;
	 }
       }
     
     }
     else
       ADONIS_ERROR(DerivedError, " r == "<<r<<" is illegal (either 'r' or 'c')");
  
  }
  

  /**
   *\brief Identity FieldMatrix stored as <B>dense</B> object. Just for Lapack. 
   */
  template<class T, int N>
  inline void identity(Dune::FieldMatrix<T,N,N>& Id){
    for(int i = 0; i < N; ++i){
      for(int j = 0; j < N; ++j){
	Id[i][j] = kronecker_delta<T>(i,j);
      }
    }
  }


  template<class T, int N>
  inline void identity(){
    Dune::FieldMatrix<T,N,N> Id;
    for(int i = 0; i < N; ++i){
      for(int j = 0; j < N; ++j){
	Id[i][j] = kronecker_delta<T>(i,j);
      }
    }
  }


  /**
   * \brief a least squares solution for underdetermined rank-deficient isn't unique and depends strongly on the threshold you have specified.
   * Ideal threshold: macheps*norm(A)
   * \tparam T value type
   * \tparam NORM character specifying the desired norm (default \f$ \|A\|_{\max}. \f$
   */
    template<class T, char NORM = 'm'>
    class LeastSquaresThreshold{
    public:
      LeastSquaresThreshold(const T& t = std::numeric_limits<T>::epsilon()):thresh_(t){
#ifndef NDEBUG
	if(t < (*this).macheps())
	  ADONIS_WARNING(Warning,"t = "<<t << " < macheps = "<<(*this).macheps()<<".");
#endif
      }

      template<int M, int N>
      T& rcond(const Dune::FieldMatrix<T,M,N>& A, const T& fac = 1, char rc = 'c'){
	thresh_ *= (fac * norm(A,NORM,rc));
	return thresh_;
      }

      T& rcond(){
	return thresh_;
      }

      T macheps(){
	return std::numeric_limits<T>::epsilon();
      }

    private:
      T thresh_;
    };
   



  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////// LAPACK/BLAS dependent stuff ////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

#if USE_LAPACK
#ifdef LAPACK_USAGE


  /**
   * \brief compute the inverse of a FieldMatrix using LU decomposition
   * Remember \f[ AX = I, \f]
   * where \f$ X = A^{-1}.\f$
   * A will be overwritten by L and U 
   * B on entry the identity will store the solution \f$ A^{-1}.\f$
   */
  template<class T, int N>
  inline Dune::FieldMatrix<T,N,N> inverse(Dune::FieldMatrix<T,N,N>& A, bool showdet = false){
    
    Dune::FieldMatrix<T,N,N> B;    //the solution
    
    if(N==1){
      adonis_assert(Abs(A[0][0]) > T());
      B[0][0] = static_cast<T>(1)/A[0][0];
    
      if(showdet)
	std::cout <<"determinant = "<< A[0][0] <<"."<< std::endl;
    }
    else if(N==2){
      T det = A[0][0]*A[1][1]-A[1][0]*A[0][1];
      adonis_assert(Abs(det) > T());
      B[0][0] = A[1][1];  B[0][1] = -A[0][1];
      B[1][0] = -A[1][0]; B[1][1] = A[0][0];
      B /= det;

      if(showdet)
	std::cout <<"determinant = "<< det <<"."<< std::endl;
    }
    else{
      
      int n = N,
	nrhs = n,
	lda = n,         //>= max(1,N)
	ldb = n,         //>= max(1,N)
	info;
      
	int ipiv[N];
	
   
	identity(B);               //the rhs and later solution

	F77Type<TypeTraits<T>::Value >().GESV(&n,&nrhs,&A[0][0],&lda,&ipiv[0],&B[0][0],&ldb,&info);
	
	//show determinant
	if(showdet){   
	  T det = 1;
	  for(int i = 0; i < N; ++i)
	    det *= A[i][i];
	  std::cout <<"determinant = "<< det <<"."<< std::endl;
	}
	

#ifndef NDEBUG
	if(info < 0)
	  ADONIS_ERROR(LapackError,"The " << info <<"th argument has an illegal value.");
	
	if(info > 0)
	  ADONIS_ERROR(LapackError,"U(" << info <<", "<<info<<") = "<<0<<" (in C++: U["<<info-1<<"]["<<info-1<<"] = 0). U is singular ==>\n   No solution could be computed, pal.");
      
#endif
    }

    return B;
  }


  /**\brief Alternative matrix inversion using LAPACK's ?getri after ?getrf. 
   *Use LU decomposition and do not consider any cases, e.g. where N==1 or N==2. This is equivalent to the above routine, cf. [HIGHAM, §14, pp. 267]
*/
  template<class T, int N>
  inline Dune::FieldMatrix<T,N,N>& inversion(Dune::FieldMatrix<T,N,N>& A){
    int n = N,
      lda = n,
      info;
	
    int ipiv[N];
    
    //LU decomposition of A
    F77Type<TypeTraits<T>::Value >().GETRF(&n,&n,&A[0][0],&lda,&ipiv[0],&info);
    if(info > 0){
      ADONIS_ERROR(LapackError," LU factorization: U("<<info<<","<<info<<") ~ 0. \n    Division by zero will occur.");
    }
    if(info < 0){
      ADONIS_ERROR(LapackError,"The " << info <<"th argument has an illegal value.");
    }
    int lwork = N;
    T work[N]; //lwork
    
    //calulate inverse of A based on the previous LU decomposition
    F77Type<TypeTraits<T>::Value >().GETRI(&n,&A[0][0],&lda,&ipiv[0],&work[0],&lwork,&info);
    
#ifndef NDEBUG
    if(info < 0)
      ADONIS_ERROR(LapackError,"The " << info <<"th argument has an illegal value.");
	
    if(info > 0)
      ADONIS_ERROR(LapackError,"U(" << info <<", "<<info<<") = "<<0<<" (in C++: U["<<info-1<<"]["<<info-1<<"] = 0). Matrix is singular ==>\n   No Inverse exists!");
    
#endif
    
    return A;
  }




  /**
   * \brief Copy smaller fmatrix into bigger one.
   */
  template<class T, int M, int NRHS, int S>
  inline void copy(const Dune::FieldMatrix<T,M,NRHS>& B,
		   Dune::FieldMatrix<T,S,NRHS>& Sol){
    adonis_assert(S >= M);
    for(int i = 0; i < M; ++i){
      for(int j = 0; j < NRHS; ++j){
	Sol[i][j] = B[i][j];
      }
    }
  }
  

  /**
   * \brief for a chemical source \f[ S: \mathbb{R}^D  \rightarrow \mathbb{R}^R\] a matrix \f$ A \in \mathbb{R}^{R \times D}\f$ is assigned.
   *
   *  \tparam T value type
   *  \tparam D domain space dimension
   *  \tparam R range space dimension
   *  \tparam MT character designating the matrix form (currently available: 'f'
   *   dense Dune::FieldMatrix, 'd' dense FieldDiagonal)
   */
  template<class T, int D, int R, template<class S> class F, char MT>
  class CanonicalSourceTermIsomorphism{
  public:
    typedef MatrixTypeSelector<T,R,D,MT> MSelector;
    typedef typename MSelector::MatrixType MatrixType;
    typedef F<T> FunType;
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> XVecType;

    CanonicalSourceTermIsomorphism():S_(R){}
    
    //creates a matrix having a (in general) non-zero 1st row, zeros elsewhere
    template<class V>
    MatrixType& fill_1st_row(const V& zM ){
      A_ = T();
      return MSelector::fill_1st_row(S_(zM),A_);
    }

    //creates a matrix havin a (in general) non-zero 1st column, zeros elsewhere
    template<class V>
    MatrixType& fill_1st_column(const V& zM){
      A_ = T();
      return MSelector::fill_1st_column(S_(zM),A_);
    }

  private:
    FunType S_;
    MatrixType A_;
  };

 


  /**
   * \brief Solve a symmetric linear system AX = B. B will overwritten to hold the solution later on.
   */
  template<class T, int N, int NRHS>
  inline void solve_symmetric_ls(const Dune::FieldMatrix<T,N,N>& A, Dune::FieldMatrix<T,N,NRHS>& B, char uplo = 'U', char rc = 'c'){
    adonis_assert(uplo == 'U' || uplo == 'L');
    adonis_assert(rc == 'c' || rc == 'r');

    Copy2DynamicArray<T,N,N> a(A,rc);
    Copy2DynamicArray<T,N,NRHS> b(B,rc);
    

    int n = N,
      nrhs = NRHS,
      lda = std::max(1,n),
      ldb = std::max(1,n),
      lwork = -1,
      info;

    int* ipiv = new int [N];
    T wkopt;

    //calculate optimal lwork -- nothing else
    F77Type<TypeTraits<T>::Value>().SYSV(&uplo,&n,&nrhs,a.get_ptr(),&lda,&ipiv[0],b.get_ptr(),&ldb,&wkopt,&lwork,&info);
    
    //take only the real part in case you calculate with complex values
    lwork = LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::real_part(wkopt); 

    T* work = new T [lwork*sizeof(T)];
    
    //solve the symmetric system
    F77Type<TypeTraits<T>::Value>().SYSV(&uplo,&n,&nrhs,a.get_ptr(),&lda,&ipiv[0],b.get_ptr(),&ldb,&work[0],&lwork,&info);

    delete[] work;
    delete[] ipiv;

    transform_2_matrix(b.begin(),b.end(),B,rc);  //B is overwritten with solut.
  }
  
  /**
   * \brief Solve <B>symmetric positive definit </B> system with the lhs matrix \f$ A\f$ being s.p.d. via Cholesky decomposition.  
   */
  template<class T, int N, int NRHS>
  inline void solve_positive_definite_ls(const Dune::FieldMatrix<T,N,N>& A, Dune::FieldMatrix<T,N,NRHS>& B, char uplo = 'U', char rc = 'c'){
    adonis_assert(uplo == 'U' || uplo == 'L');
    adonis_assert(rc == 'c' || rc == 'r');
    
    Copy2DynamicArray<T,N,N> a(A,rc);
    Copy2DynamicArray<T,N,NRHS> b(B,rc);
   
    int n = N,
      nrhs = NRHS,
      lda = std::max(1,n),
      ldb = std::max(1,n),
      info;

    //solve system
    F77Type<TypeTraits<T>::Value>().POSV(&uplo,&n,&nrhs,a.get_ptr(),&lda,b.get_ptr(),&ldb,&info);
    

    if(info < 0) 
      ADONIS_ERROR(LapackError, "The "<<info<<"th argument has an illegal value." );

    if(info > 0)
      ADONIS_ERROR(LapackError, "Leading minor of order "<<info<<" of A is not p.d. \n   ==> No Cholesky factorization computed ==> no solution.");

  
     transform_2_matrix(b.begin(),b.end(),B,rc);  //B stores solution now
  }

  /**
   *\brief Solve Over-/Underdetermined <B>rank-deficient, dense</B> linear system using LAPACK's ?gelsd.
   * See, e.g.:
   * <a href="http://www.netlib.org/lapack/double/dgelsd.f"> Rank-deficient systems</a>
   * \param A \f$ M \times N \f$ matrix on input which will be destroyed on exit
   * \param B \f$ M \times NRHS \f$ rhs matrix on entry, the \f$ N \times NRHS \f$ solution on exit. Therefore the actual storage of the rhs \f$B \f$ will be an  array of dimension \f$ LDB= \max(M,N) \times NRHS. \f$ 
   * \param rcond for determination of effective rank; if rcond < 0 (default) machine precision is used
   * \param rc 'r'ow-wise or 'c'olumn-wise storage for object 'Copy2DynamicArray<T,M,N>'
   * \return a \f$ LDA:= \max(M,N) \times NRHS \f$ FieldMatrix
   * NOTE: since A is destroyed and B converted into a differently sized object,  it might be beneficial to copy the stuff to dynamic arrays first! 
   *
   * For detailed information, visit
   * <a href="http://www.netlib.org/lapack/lug/node27.html#llsq"> Linear Least Squares Problems</a>
   */
  template<class T, int M, int N, int NRHS> 
  inline Dune::FieldMatrix<T,N,NRHS> solve_rank_deficient_ls(const Dune::FieldMatrix<T,M,N>& A, const Dune::FieldMatrix<T,M,NRHS>& B,
 typename F77Type<TypeTraits<T>::Value>::BaseType rcond = -1, char rc = 'c'){
    adonis_assert(M != 0 || N != 0 || NRHS != 0);

    //copy stuff
    Copy2DynamicArray<T,M,N> a(A,rc);
    Copy2DynamicArray<T,M,NRHS> b(B,rc);
    //print_in_matrix_style(b.begin(),b.end(),NRHS);
    
    const int LDB = (M > N) ? M : N;
    const int MINMN = (M < N) ? M : N;
    b.restructure(LDB); //(N); //(LDB);
    
    
    int m = M,
      n = N,
      nrhs = NRHS,
      lda = (M > 1) ? M : 1,
      ldb = LDB,
      minmn = MINMN, //std::min(M,N),
      smlsiz = 25,   //usual value
      lwork = -1,   //automatically determine optimal worksize by a call of
                    //dgelsd_ which ONLY calculates the optimal worksize 
      
      log2term = (int)(log2(static_cast<double>(minmn)/(smlsiz+1))) + 1, 
      nlvl = std::max(0,log2term),
      liwork = std::max(1,3*minmn*nlvl + 11*minmn),
      rank,
      info;

    //output
    //std::cout<< "M = "<<M <<std::endl<<"N = "<<N<< std::endl << "LDA = "<<lda <<std::endl<< "LDB = "<<ldb<<std::endl<<std::endl<<"NRHS = "<<nrhs<<std::endl<<"minmn = "<<minmn << std::endl<<"log2term = "<<log2term<<std::endl<<"nlvl = "<<nlvl << std::endl<< "liwork = "<<liwork <<std::endl<< std::endl<< "lwork = "<<lwork<<std::endl;

    typedef typename F77Type<TypeTraits<T>::Value>::BaseType BaseType;

    //workspaces 
    BaseType* S = new BaseType [minmn]; //the sv in decreasing order
    int* iwork = new int [liwork];

    //complex stuff
    int lrwork;
    
    LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::lwk(lrwork,M,N,smlsiz,nlvl,NRHS); //void fct that assigns lrwork a value
    
    //std::cout << "lrwork = "<<lrwork << std::endl;
      

    BaseType* rwork = new BaseType [std::max(1,lrwork)]; 
   
    T wkopt;  // stores optimal work space later on
    //determine optimal lwork size automatically -- only this is computed
    F77Type<TypeTraits<T>::Value>::GELSD(&m,&n,&nrhs,a.get_ptr(),&lda,b.get_ptr(),&ldb,&S[0],&rcond,&rank,&wkopt,&lwork,&rwork[0],&iwork[0],&info);
   
    lwork = LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::real_part(wkopt);  
    lwork = (1 > lwork) ? 1 : lwork;      //std::max(1,(int)wkopt);
    //std::cout << "lwork (auto-determined) = "<< lwork << std::endl;

    //std::cout << "lwork_opt = "<<lwork<<std::endl<<std::endl<< "lwork*sizeof(T) = "<<lwork*sizeof(T) <<std::endl;

    T* work = new T [lwork*sizeof(T)]; 
    
    F77Type<TypeTraits<T>::Value>::GELSD(&m,&n,&nrhs,a.get_ptr(),&lda,b.get_ptr(),&ldb,&S[0],&rcond,&rank,&work[0],&lwork,&rwork[0],&iwork[0],&info);
    std::cout << "EFFECTIVE RANK = "<<rank<< " (of "<< M<<"x"<<N<<" lhs matrix)"<<std::endl;

    std::cout<< "SV's computed within xgelsd_ (in decreasing order) = " <<std::endl;
    ManipulateDynamicArray<BaseType,int,0,MINMN>::output(S);
    std::cout << std::endl;

#ifndef NDEBUG
    if(info < 0)
      ADONIS_ERROR(LapackError, "The "<<Abs(info)<<"th argument had an illegal value.");
    if(info > 0)
      ADONIS_ERROR(LapackError, "SVD algorithm diverged: "<<info<<"th off-diagonal elements of an intermediate bidiagonal form did not converge to zero.");
#endif

    //check output as stored in the dynamic vector
    // std::cout << "Solution (calculated by ?gelsd) = "<<std::endl<< b << std::endl;

    //std::cout << std::endl<< "b.size() = "<<b.size()<<std::endl;

    //delete the by 'new' defined quantities
    delete[] rwork;
    delete[] work;
    delete[] S;
    delete[] iwork;
    
    Dune::FieldMatrix<T,N,NRHS> Sol;
    std::cout << "SOLUTION b = "; print_in_matrix_style(b.begin(),b.end(),nrhs,"Solution matrix",5);
    //assign to FieldMatrix of appropriate dimensions (LDB x NRHS)
    transform_2_matrix(b.begin(),b.end(),Sol,rc);  

    return Sol;  
  }
  


  /** 
   * \brief Overload rank deficient routine
   */
  template<class T, int M, int N> 
  inline Dune::FieldVector<T,N> solve_rank_deficient_ls(const Dune::FieldMatrix<T,M,N>& A, const Dune::FieldVector<T,M>& B, typename F77Type<TypeTraits<T>::Value>::BaseType rcond = -1, char rc = 'c'){
     
    //std::cout << "-------------TEST ME, PAL ------------"<<std::endl;

    NumericDataTypeChecker<T>::certify(); //if type is below double, throw error

    typedef typename TypeAdapter<T>::Type ProperType;
    
    typedef typename F77Type<TypeTraits<T>::Value>::BaseType BaseType;

    //copy lhs matrix
    Copy2DynamicArray<T,M,N,typename TypeAdapter<T>::Type> a(A,rc);
    

    const int LDB = (M > N) ? M : N;
    const int MINMN = (M < N) ? M : N;
    
    //copy rhs vector
    ProperType* solution = new ProperType [LDB];                //solution
    for(int i = 0; i < M; ++i){                          //fill with B[i]  
      smart_assign(solution[i],B[i]);
    }
   

    int m = M,
      n = N,
      nrhs = 1,
      lda = M,
      ldb = LDB,
      minmn = MINMN, //std::min(M,N),
      smlsiz = 25,   //usual value
      lwork = -1,   //automatically determine optimal worksize by a call of
                    //dgelsd_ which ONLY calculates the optimal worksize 
      
      log2term = (int)(log2(static_cast<BaseType>(minmn)/(smlsiz+1))) + 1, 
      nlvl = std::max(0,log2term),
      liwork = std::max(1,3*minmn*nlvl + 11*minmn),
      rank,
      info;

    //workspaces 
    BaseType* S = new BaseType [minmn]; //the sv in decreasing order
    int* iwork = new int [liwork];

    //complex stuff
    int lrwork;
    
    LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::lwk(lrwork,M,N,smlsiz,nlvl,nrhs); //void fct that assigns lrwork a value
    
    //std::cout << "lrwork = "<<lrwork << std::endl;
      

    BaseType* rwork = new BaseType [std::max(1,lrwork)]; 
   
    ProperType wkopt;  // stores optimal work space later on
    //determine optimal lwork size automatically 
    F77Type<TypeTraits<T>::Value>().GELSD(&m,&n,&nrhs,a.get_ptr(),&lda,&solution[0],&ldb,&S[0],&rcond,&rank,&wkopt,&lwork,&rwork[0],&iwork[0],&info);
   
    lwork = LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::real_part(wkopt);  
    lwork = (1 > lwork) ? 1 : lwork;      //std::max(1,(int)wkopt);
    //std::cout << "lwork (auto-determined) = "<< lwork << std::endl;

    //std::cout << "lwork_opt = "<<lwork<<std::endl<<std::endl<< "lwork*sizeof(T) = "<<lwork*sizeof(T) <<std::endl;

    ProperType* work = new ProperType [lwork*sizeof(ProperType)]; 
    
    F77Type<TypeTraits<T>::Value>().GELSD(&m,&n,&nrhs,a.get_ptr(),&lda,&solution[0],&ldb,&S[0],&rcond,&rank,&work[0],&lwork,&rwork[0],&iwork[0],&info);
    
    //std::cout << "EFFECTIVE RANK = "<<rank<< " (of "<< M<<"x"<<N<<" lhs matrix)"<<std::endl;

    //std::cout<< "SV's computed within xgelsd_ (in decreasing order) = " <<std::endl;
    //ManipulateDynamicArray<BaseType,int,0,MINMN>::output(S);
    //std::cout << std::endl;

    if(info < 0)
      ADONIS_ERROR(LapackError, "The "<<Abs(info)<<"th argument had an illegal value.");
    if(info > 0)
      ADONIS_ERROR(LapackError, "SVD algorithm diverged: "<<info<<"th off-diagonal elements of an intermediate bidiagonal form did not converge to zero.");

    //assign back to output container
    typedef Dune::FieldVector<T,N> OutContainerType;
    OutContainerType Sol;
    for(int i = 0; i < LDB; ++i) Sol[i] = solution[i];
    //FillLinearContainer<OutContainerType,0,N>::fill_me(Sol.begin(),solution);
    


    //delete the by 'new' defined quantities
    delete[] rwork;
    delete[] work;
    delete[] S;
    delete[] iwork;
    delete[] solution;

    return Sol;  
  }
  
  
  /***
   * \brief Extract up to row k
   */
  template<class M1, class M2>
  inline void extract_up_to_row(int k, M1& A, const M2& G){
    adonis_assert(k >0 && k < G.rows);
    for(int i = 0; i < k; ++i)
      for(int j = 0; j < G.cols; ++j)
	A[i][j] = G[i][j];
  }




  /**
   * \brief Calculate the (Moore-Penrose-)pseudoinverse for a general matrix \f$  A \in \mathbb{R}^{M \times N}.\f$ 
   * Check with Matlab's(r) pinv
   * 
   * NOTE: in the underdetermined case, i.e. \f$ M < N \f$ the result is different from Matlab as it is not using xgelsd nor xgelss nor xgelsy but it's own hack. So don't panic.
   See <a href="http://icl.cs.utk.edu/lapack-forum/viewtopic.php?p=497&">pseudo-inverse </a>
   *
   * NOTE: octave yields the same results. 
  */
  template<class T, int M, int N>
  inline Dune::FieldMatrix<T,N,M> pseudo_inverse(const Dune::FieldMatrix<T,M,N>& A, typename F77Type<TypeTraits<T>::Value>::BaseType rcond = -1){
    //! rhs B is a identity matrix of order \f$ \max(M,N)\f$ 
    const int mx =  M;   //(M > N) ? M : N;
    
    Dune::FieldMatrix<T,mx,mx> B;  //rhs, just appropriate identity
    identity(B);
    
    Dune::FieldMatrix<T,N,M> PI;
    
    const char c = 'c';

    if(M > N)
      extract_up_to_row(N,PI,solve_rank_deficient_ls(A,B,rcond,c));
    else
      PI = solve_rank_deficient_ls(A,B,rcond,c);
      
    return PI;
  } 
  

  /**
   * \brief If the rows of a real matrix \f$ A \in \mathbb{R}^{m \times n}, \ m \leq n\f$ are linearly independent, and explicit formula is given via \f[ A^{+} = A^T(AA^T)^{-1}. \]
   *
   * NOTE: if the rows are orthonormal, like in the case of \f$ N^T\f$, then the above formula simplifies to \f[  A^{+} = A^T. \f]
   *
   * NB: No efforts are made to check whether A is orthonormal or not.
   */
  template<class T, int M, int N>
  inline Dune::FieldMatrix<T,N,M> pseudo_inverse_full_row_rank(const Dune::FieldMatrix<T,M,N>& A){
    Dune::FieldMatrix<T,N,M> AT = transpose(A);
    Dune::FieldMatrix<T,M,M> Iv =  A*AT;
    Iv.invert();  //alt.: inversion(Iv);   
    return AT*Iv;  
  }


  /**
   * \brief Performs calculation of entry \f$ I_{ii} - N_i \cdot N_j^T without actually transposing matrix \f$ N \f$ 
   */
  template<class M>
  inline typename M::field_type N_i_times_N_T_j(int i, int j, const M& A){
    typename M::field_type dp = 0;
    for(int k = 0; k < A.rows; ++k){
      dp += A[k][i]*A[k][j];
    }
    if(i == j)
      dp = 1 - dp;
    return dp;
  }
  
  /**
   * \brief Computes various field matrix norms via LAPACK
   * \param A \f$ M \times N\f$ FieldMatrix
   * \param norm character defining the norm. Available are
   *   \f$ \| A \|_{\max}\f$     'M' or 'm'
   *   \f$ \| A \|_1 \f$         '1', 'O' or 'o'
   *   \f$ \| A \|_{\infty}\f$   'I' or 'i'
   *   \f$ \| A \|_{\operatorname{F}} \f$ 'F', 'f', 'E' or 'e'
   *   \f$ \| A \|_2 \f$ '2' (the largest singular value, computed via xgesvd_)
   * \param rc specifies the storage format ('c'olumn by default, or 'r'ow)
   * \return depending on the precision, either single or double
   */
  template<class T, int M, int N>
  inline typename F77Type<TypeTraits<T>::Value>::BaseType norm(const Dune::FieldMatrix<T,M,N>& A, char norm, char rc = 'c'){
    adonis_assert(rc == 'c' || rc == 'r');
    adonis_assert(norm == 'M' || norm == 'm' ||
		  norm == '1' || norm == 'O' || norm == 'o' ||
		  norm == 'I' || norm == 'i' ||
		  norm == 'F' || norm == 'f' || norm == 'E' || norm == 'e' ||
		  norm == '2'); //the largest sv
     
    Copy2DynamicArray<T,M,N> a(A,rc);
    
    int m = M,
      n = N,
      lda = std::max(m,1),
      lwork = m;  //except for the infinity norm, workspace 'work' 
                  //won't be referenced;

    typedef typename F77Type<TypeTraits<T>::Value>::BaseType BaseType;
    BaseType work[std::max(1,lwork)];
    
    BaseType normOfA = -1;  //norm is never negative
    if(norm == '2'){
      normOfA = svd(A,'A','A',rc)[0]; //the largest value is just the 1st value
    }
    else
      normOfA = F77Type<TypeTraits<T>::Value>().LANGE(&norm,&m,&n,a.get_ptr(),&lda,&work[0]);
    
    adonis_assert(normOfA >= T());  //if norm < 0, something's got wrong!

    return normOfA;
  }


  /**
   *  \brief Solve a possibly rank deficient system using QR  
   * <a href="http://www.netlib.org/lapack//double/dgelsy.f"> Rank-deficient systems using QR</a>
   * xgelsd_ and xgelsy_ are the recommended routines when solving rank-deficient least squares and both routines' cost are \f$ \mathcal{O}(M\cdot N^2)\f$, although the latter one is still a bit faster due to qr decomposition. 
   * NOTE: Seems to have a bug within xgelsy_: Some uninitialized data within usr/lib/sse2/atlas/libblas.so.3gf.0
   */
  template<class T, int M, int N, int NRHS>
  inline Dune::FieldMatrix<T,((M > N)? M : N), NRHS> solve_rank_deficient_ls_via_qr(const Dune::FieldMatrix<T,M,N>& A, const Dune::FieldMatrix<T,M,NRHS>& B, typename F77Type<TypeTraits<T>::Value>::BaseType rcond = -1, char rc = 'c'){
    Copy2DynamicArray<T,M,N> a(A,rc);
    Copy2DynamicArray<T,M,NRHS> b(B,rc);
    const int LDB = (M > N) ? M : N;
    b.restructure(LDB);


    //std::cout << "a ("<<a.size()<<" entries) = "<<a<<std::endl<<std::endl<<"b ("<<b.size()<<" entries) = "<<b<<std::endl;


    int m = M,
      n = N,
      nrhs = NRHS,
      lda = std::max(1,M),
      ldb = LDB,
      rank,
      lwork = -1,
      info;

    // std::cout <<"M = "<<m<<std::endl<<"N = "<<n<<std::endl<<"NRHS = "<<nrhs << std::endl<< "LDA = "<<lda<<std::endl<<"LDB = "<<ldb<<std::endl;

    // ============================ INPUT =====================================
    // Due to the above link, if jpvt[i] != 0 the i-th column of A is permutated
    //to the front of AP, otherwise column i is a free column
    //initialise jpvt to be zero so that all columns are free, see
    //  http://www.nag.co.uk/lapack-ex/examples/source/dgelsy-ex.f
    int* jpvt = new int [N];
    ManipulateDynamicArray<int,int,0,N>::fill(jpvt,T()); //fill with zeros
    //========================================================================

    T wkopt;
    
    //only for complex stuff
    typedef typename F77Type<TypeTraits<T>::Value>::BaseType BaseType;
    BaseType* rwork = new BaseType [LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::lwk(N)];  //this time an integer is return in []
    

    //automatical determination of optimal workspace (only workspace computed)
    F77Type<TypeTraits<T>::Value>().GELSY(&m,&n,&nrhs,a.get_ptr(),&lda,b.get_ptr(),&ldb,&jpvt[0],&rcond,&rank,&wkopt,&lwork,&rwork[0],&info);
     
    lwork =  LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::real_part(wkopt)*sizeof(T);
    //std::cout << "lwork = "<<lwork << std::endl;
    T* work = new T [std::max(1,lwork)];
    
    F77Type<TypeTraits<T>::Value>().GELSY(&m,&n,&nrhs,a.get_ptr(),&lda,b.get_ptr(),&ldb,&jpvt[0],&rcond,&rank,&work[0],&lwork,&rwork[0],&info);

    adonis_assert(info == 0);

    std::cout <<"effective rank = "<<rank << " of "<<M <<"x"<<N <<" matrix"<< std::endl;
    //std::cout << "Solution of rank-deficient system solved by QR = "<<std::endl<< b <<std::endl;

    delete[] rwork;
    delete[] work;
    delete[] jpvt;
  
    Dune::FieldMatrix<T,LDB,NRHS> Sol;
    transform_2_matrix(b.begin(),b.end(),Sol,rc);  

    return Sol;  
  }


  /**
   * \brief Singular value decomposition for a <B>dense</B> general matrix.
   * Singular values are stored in a real ((M < N) ? M : N)-dimensional FieldVector s, such that s[i] >= s[i+1] 
   */
  template<class T, int M, int N>
  inline Dune::FieldVector<typename F77Type<TypeTraits<T>::Value>::BaseType,((M < N) ? M : N)> svd(const Dune::FieldMatrix<T,M,N>& A, char jobu = 'A', char jobvt = 'A', char rc = 'c'){
    adonis_assert((jobu == 'A' || jobu == 'S' || jobu == 'O' || jobu == 'N') &&
		  (jobvt == 'A' || jobvt == 'S' || jobvt == 'O' ||jobvt == 'N'));
    
    Copy2DynamicArray<T,M,N> a(A,rc);
    const int minmn = (M < N) ? M : N;

    int m = M,
      n = N,
      lda = std::max(1,m),
      ldu = 1,
      ldvt = 1,
      lwork = -1,
      info;
  
    if(jobvt == 'A')
      ldvt = N;

    if(jobvt == 'S')
      ldvt = minmn;

    T wkopt;
    
    typedef F77Type<TypeTraits<T>::Value> FORT;
    typedef typename FORT::BaseType BaseType;

    Dune::FieldVector<BaseType,minmn> s;
    if(jobu == 'S' || jobu == 'A')
      ldu = m;
    
    int udim = ldu*m;

    if(jobu == 'S')
      udim = ldu*minmn;
    
    T* u = new T [udim]; 
    T* vt = new T [ldvt*n];
    
    //for the complex version only
    
    T* rwork = new T [LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::lwk(minmn,5)];
 
   
    //determine lwork automatically
    FORT().GESVD(&jobu,&jobvt,&m,&n,a.get_ptr(),&lda,&s[0],&u[0],&ldu,&vt[0],&ldvt,&wkopt,&lwork,&rwork[0],&info);

    lwork = std::max(1, LapackComplex<F77Type<TypeTraits<T>::Value>::IsComplex>::real_part(wkopt));
    
    T* work = new T [lwork*sizeof(T)];

    //compute svd 
    FORT().GESVD(&jobu,&jobvt,&m,&n,a.get_ptr(),&lda,&s[0],&u[0],&ldu,&vt[0],&ldvt,&work[0],&lwork,&rwork[0],&info);
   
    delete[] work;
    delete[] rwork;
    delete[] u;
    delete[] vt;

 
#ifndef NDEBUG
    if(info > 0)
      ADONIS_ERROR(LapackError, "SVD algo did not converge!");

    if(info < 0)
      ADONIS_ERROR(LapackError, "The "<<info<<"th argument has an illegal value.");
#endif

    
    return s;

  }

  
#endif
#endif

  /////////////////////////// END OF LAPACK/BLAS DEPENDENCIES //////////////////

  
  
} //end namespace 

#endif
