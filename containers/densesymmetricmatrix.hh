#ifndef DENSE_SYMMETRIC_MATRIX_HH
#define DENSE_SYMMETRIC_MATRIX_HH

#include <iostream>

#include "../common/adonisassert.hh"
#include "../expressiontemplates/expressionvectortraits.hh"
#include "../misc/useful.hh"
#include "../misc/misctmps.hh"
#include "../linalg/linearsystemsolvers.hh"
#include "../common/elementaryoperations.hh"
#include "../common/globalfunctions.hh"

namespace Adonis{

  /**
   * \brief Storage representation of symmetric dense matrix
   *
   * Reference: [GOLUB/VANLOAN, "Matrix Computations", 3rd ed., § 1.2.7]
   */
  template<class US>
  class SymmetricStorageByColumn{
  public:
    typedef US index_value;

    static const char Value = 'c';
 
    static inline US index(const US& n, const US& i, const US& j){
      adonis_assert(i >= j);
      return ((j*n - (j+1)*j/2 + i));  //i >= j
    }
  
    //!fill from full dense symmetric matrix stored as vector
    template<class SY,class V>
    static inline void fill(SY& symm, const V& mat){
      //ADONIS_ERROR(UndefinedFunctionError,"In order to use this function properly, you have to define its body first ;)");
      unsigned n = symm.order();
      adonis_assert(mat.size() == ntimes<2>(n));
      
      unsigned count = 0;
      for(unsigned j = 0; j < n; ++j){
	for(unsigned i = 0; i < n-j; ++i){
	  symm[count++] = mat[i + j*(n+1)];
	}
      }

    }
    
  };

  template<class US>
  class SymmetricStorageByDiagonal{
  public:
    typedef US index_value;
    
    static const char Value = 'd';
    
    static inline US index(const US& n, const US& i, const US& k){
      adonis_assert(k >= 0);
      return (i + n*k - k*(k-1)/2);      //k >= 0
    }
  
    template<class SY,class V>
    static inline void fill(SY& symm, const V& mat){
      unsigned n = symm.order();
      adonis_assert(mat.size() == ntimes<2>(n));
      
      for(unsigned k = 0; k < n; ++k){
	for(unsigned i = 0; i < n-k; ++i){
	  symm[index(n,i,k)] = mat[i*n + (i+k)];//mat_{i,i+k} are the subdiag entries
	}
      }
    }

  };

  /**
   * \brief Object for storing a <B> dense symmetric </B> matrix. Store only upper(lower) part of symmetric dense matrix. 
   * \tparam T value_type
   * \tparam V any STL compliant random access container that stores the triangular matrix
   * \tparam R Storage Organisation (storage by diagonal by default)
   */
  template<class T, template<class S, class A = std::allocator<S> > class V = ExprTmpl::MyVec, class R = SymmetricStorageByDiagonal<unsigned> >
  class DenseSymmetricMatrix{
  public:
    typedef V<T> VecType;
    
    typedef typename VecType::iterator iterator;
    typedef typename VecType::const_iterator const_iterator;

    DenseSymmetricMatrix(unsigned n = 0, const T& d = T()):n_(n), v_(n*(n+1)/2){
      for(unsigned i = 0; i < v_.size(); ++i)
	v_[i] = d;
    }

    template<class ITER>
    DenseSymmetricMatrix(unsigned n, ITER first, ITER last):n_(n){
      adonis_assert(n*(n+1)/2 == static_cast<unsigned>(std::distance(first,last)));
      for(ITER it = first; it != last; ++it){
	v_.push_back(*it);
      }
    }


    /*//! assign from given vector
      template<class W> 
    DenseSymmetricMatrix& operator=(const W& w){
      adonis_assert(v_.size() == w.size());
      for(unsigned i = 0; i < v_.size(); ++i)
	v_[i] = w[i];

      return *this;
      }*/


    char storage_organisation() {return R::Value;} 

    unsigned size() const{
      return  v_.size();// alternatively return "n_*(n_+1)/2"
    }


    T& operator[](unsigned i){
      adonis_assert(i < size());
      return v_[i];
    }

    const T& operator[](unsigned i) const{
      adonis_assert(i < size());
      return v_[i];
    }

    unsigned order() const {return n_;}

    void change_order(unsigned n1){
      n_ = n1;
    }

    //a tiny iteration interface
    iterator begin(){
      return  v_.begin();
    }
  
    iterator end(){
      return v_.end();
    }

    const_iterator begin() const{
      return v_.begin();
    }
    
    const_iterator end() const{
      return v_.end();
    }

    
    //!matrix-like indexing -- DON'T use it like that! It just gives \f$ a_{i+k,i}, k \geq  j = i+k \f$
    unsigned index(unsigned i, unsigned j) const {
      return R::index(n_,i,j);
    }
    

    const VecType& get_container() const {return v_;}
    VecType& get_container() {return v_;}
    
    /**
     * \brief symmetric matrix stored as <I> full</I> triangular matrix in a random access container.
     * 
     * \tparam RCM row/column major storage. If RCM = RowMajor, then the lower part is filled and with RCM = ColumnMajor the opposite is true
     * \tparam W random access container that holds the matrix
     */
    template<class RCM, class W>
    void fill_me_into_a_full_matrix_given_as_rac(W& w) const{
      adonis_assert(w.size() == ntimes<2>(n_)); //must be square

      for(unsigned k = 0; k < n_; ++k){   //! \f$ k\f$-th (sub)diagonal
	for(unsigned i = 0; i < n_-k; ++i){
       
	  w[RCM::offset(i,i+k,n_)] = v_[(*this).index(i,k)];
	
	}
      }
    }
    
    
    /**
     * \brief this is needed for 'regularization' stuff, i.e. to make the symmetric matrix positive definite. See [FLETCHER, "Practical Methods of Optimization", 2nd ed, p. 30 ?], i.e. \f[ A \longleftarrow (A +\nu I), \f] where \f$ nu \f$ is a small constant if \f$ A \f$ is close to p.d. . Otherwise \f$ \nu geq 0\f$ can be determined as in algorithm (5.2.7), p. 102/03 in the reference above.
     */
    template<class OP>
    void update_diagonal(const T& nu = 0.021568){
      for(unsigned i = 0; i < n_; ++i)
	OP::apply(v_[i],nu);
    }

#if USE_LAPACK
    /**
     * \brief Solves symmetric systme using Lapack's xsysv.f
     * 
     * NOTE: output is column-based storage (F77-style): you can either set 'swapsol' true or, depending on further computations concerning LAPACK, leave it like that and access it via ColumnMajor::offset(i,j,n_). 
     */
    template<class W>
    void solve(W& rhs, int nrhs, bool swapsol = false) const{
      int order = static_cast<int>(n_);
      W a(ntimes<2>(n_));
      (*this).fill_me_into_a_full_matrix_given_as_rac<ColumnMajor>(a);
      //print_in_matrix_style(a,(*this).order());
      //!NOTE: Due to F77: if ColumnMajor then 'U'.Likewise: if RowMajor then 'L'
      solve_symmetric_ls(a,order,rhs,nrhs,'U',swapsol);  
    }
    

    /**
     * \brief Assuming positive definiteness, use cholesky factorization to compute the solution.
     * 
     * The most effective way to check if a given matrix is p. d. is to check whether the Cholesky factors \f$ LL^*\f$ exist with \f$ l_{ii} > 0\f$, cf. [FLETCHER, 2nd ed, §2, p. 15]   
     */
    template<class W>
    void solve_pos_def(W& rhs, int nrhs, bool swapsol = false) const{
      int order = static_cast<int>(n_);
      W a(ntimes<2>(n_));
      (*this).fill_me_into_a_full_matrix_given_as_rac<ColumnMajor>(a);
      
      solve_positive_definite_ls(a,order,rhs,nrhs,'U',swapsol);  
    }


    

    template<class W>
    void regularization(W& rhs, int nrhs, bool isreg = false, const T& nu = 0.021568, bool swapsol = false){
      typedef typename W::value_type value_type;
      if(!isreg)
	(*this).update_diagonal<AddBasicElements<value_type> >(nu);
      
      int order = static_cast<int>(n_);
      W a(ntimes<2>(n_));
      (*this).fill_me_into_a_full_matrix_given_as_rac<ColumnMajor>(a);
       solve_positive_definite_ls(a,order,rhs,nrhs,'U',swapsol);
    }

    
    /**
     * \brief Eigenvalues are always real for symmetric (hermitian) matrices.
     * \param e the vector which will be filled with the eigenvalues later on
     */
    template<class W>
    W& eigenvalues(W& w){
      adonis_assert(w.size() == n_);
      W a(ntimes<2>(n_));  
      (*this).fill_me_into_a_full_matrix_given_as_rac<ColumnMajor>(a);
      
      eig_symmetric_hermitian(a,n_,w,'N','U');
      
      return w;
    }


    //! Check for positive definiteness. If eigenvalues are all strictly positive then matrix is p.d. 
    bool is_positive_definite(){
      typedef typename TypeTraits<T>::BaseType BaseType;//must be a non-complex
      V<BaseType> w(n_);
      (*this).eigenvalues(w);
      
      bool ispd = false;
      for(unsigned i = 0; i < w.size(); ++i){
	if(w[i] < 0){
	  ispd = false;
	  break;    //leave loop ahead of time since the found one is not > 0  
	}
	else{
	  ispd = true;
	}
      }
      return ispd;
    }

  #endif  
    


    /**
       \brief Fills symmetric matrix from full matrix stored in random access container. 
     */
    template<class MX>
    void fill_from_full_matrix(const MX& mat, unsigned newOrder = 0){
      //std::cout << "(*this).size() = "<< (*this).size() << std::endl << "hessian.size() = "<< mat.size() << std::endl;
      if((*this).size() == 0){
	adonis_assert(newOrder > 0);
	
	n_ = newOrder;                 //changes object
	v_.resize(gauss_sum(n_));
      }
      adonis_assert(mat.size() == ntimes<2>(n_));
      R::fill(*this,mat);
    }
    

    /**
     * \brief Given a <I>column </I> based representation of a dense symmetric matrix, 
     *  transform it into one that is based on storing the <I> (sub)diagonals </I>
     *
     * \tparam CSY can be any object providing a []-operator and a size() member
     */
    template<class CSY>
    void transform(const CSY& Av){
      adonis_assert(R::Value == 'd');
      adonis_assert((*this).size() == Av.size());
      
      unsigned c = 0;
      for(unsigned k = 0; k < n_; ++k){
	for(unsigned i = 0; i < n_-k; ++i){
	  v_[c++] = Av[k + n_*i - i*(i-1)/2];
	}
      }
    }

    friend std::ostream& operator<<(std::ostream& os, const DenseSymmetricMatrix& Ds){
      std::cout << "Dense symmetric matrix:"<<std::endl;
      for(const_iterator it = Ds.begin(); it != Ds.end(); ++it)
	os << *it << "  ";
      /*for(unsigned i = 0; i < Ds.v_.size(); ++i)
	os << Ds.v_[i] << "  ";*/
      os << std::endl;
      return os;
    }

   

  private:
    unsigned n_;   //! order of matrix
    VecType v_;    //! stores values of upper (lower resp.) matrix
  };





  /** 
   * \brief Dense symmetric matrix - vector multiplication 
   *
   * Ref.: [GOLUB/VANLOAN, "Matrix Computations", 3 rd ed., § 1.2.8, p.21/22] 
   * \tparam the outputcontainer can be any random access container (not necessarily an expression vector)
   
   */
  template<class W, template<class S, class A = std::allocator<S> > class V, class R>
  inline W operator*(const DenseSymmetricMatrix<typename W::value_type,V,R>& A, const W& x){
      adonis_assert(A.order() == x.size());
      adonis_assert(R::Value == 'd'); //!this is written for diagonal storage which seems to be superior to column-based

      unsigned n = A.order();
      // std::cout << "order n = "<< n << std::endl;

      W y(n);
      
      for(unsigned i = 0; i < n; ++i)
	y[i] += A[i]*x[i];
      
      //unsigned t = 0;
      for(unsigned k = 1; k <= n-1; ++k){
	//t = n*k - k*(k-1)/2; 
	//std::cout << "t = "<< t << std::endl;
	for(unsigned i = 0; i < n-k; ++i){
	  //std::cout << "i = "<< i << "    i+t = "<< i+t << std::endl; 
	  y[i] += A[A.index(i,k)]*x[i+k]; //equiv: ...A[i+t]...
	}
	for(unsigned i = 0; i < n-k; ++i){
	  y[i+k] += A[A.index(i,k)]*x[i];
	}
      }
      return y;
    }

  
  //the same for swapped function arguments
  template<class W, template<class S, class A = std::allocator<S> > class V, class R>
  inline W operator*(const W& x, const DenseSymmetricMatrix<typename W::value_type,V,R>& A){
    return A*x;
  }


} //end namespace 

#endif
