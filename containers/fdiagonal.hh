#ifndef DIAGONAL_MATRIX_FORMAT_HH
#define DIAGONAL_MATRIX_FORMAT_HH

#include<iostream>
#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../dunestuff/fmatrix.hh"

namespace Adonis{

  /**
   * \brief a diagonal matrix is a square matrix of order N whose off-diagonal entries are all zero. Hence it can be stored in a single array of size N, endowed with some operations with general dense matrices.
   */
  template<class T, int N>
  class FieldDiagonal{
  public:
    typedef std::size_t size_type;

    typedef T value_type;
    typedef T field_type;
    typedef FieldDiagonal<T,N> ThisType;
    
    typedef T* iterator;
    typedef const T* const_iterator;

    enum{dimension = N};

    //a little STL-compliance
    size_t size() const {return static_cast<size_t>(N);}
    
    size_t rows() const {return static_cast<size_t>(N);}
    
    size_t cols() const {return static_cast<size_t>(N);}
  
    FieldDiagonal(){}

    FieldDiagonal(const T& d){
      for(int i = 0; i < N; ++i)
	a_[i] = d;
    }

    //construct from iterators
    template<class IT>
    FieldDiagonal(IT it1, IT it2){
      size_t k = 0;
      
      adonis_assert((int)std::distance(it1,it2) == N);
      
      for(IT it = it1; it != it2; ++it)
	a_[k++] = *it;
    }

    const T& operator[](size_type i) const{
      adonis_assert(i < N);
      return a_[i];
    }

    T& operator[](size_type i){
      adonis_assert(i < N);
      return a_[i];
    }

    //copy construction
    FieldDiagonal(const FieldDiagonal& D){
      for(int i = 0; i < N; ++i)
	a_[i] = D[i];
    }
    
    //copy assignment
    FieldDiagonal& operator=(const FieldDiagonal& D){
      if(this != &D){
	for(int i = 0; i < N; ++i)
	  a_[i] = D[i];
      }
      return *this;
    }
    
    FieldDiagonal& operator=(const T& d){
	for(int i = 0; i < N; ++i)
	  a_[i] = d;
      return *this;
    }
    
    
    FieldDiagonal& operator*=(const T& d){
      for(int i = 0; i < N; ++i)
	  a_[i] *= d;
      return *this;
    }

    FieldDiagonal& operator/=(const T& d){
      adonis_assert(Abs(d) > T());
      for(int i = 0; i < N; ++i)
	  a_[i] /= d;
      return *this;
    }

    FieldDiagonal& operator*=(const FieldDiagonal& D){
      for(int i = 0; i < N; ++i)
	  a_[i] *= D[i];
      return *this;
    }

    FieldDiagonal& operator+=(const FieldDiagonal& D){
      for(int i = 0; i < N; ++i)
	  a_[i] += D[i];
      return *this;
    }

    FieldDiagonal& operator-=(const FieldDiagonal& D){
      for(int i = 0; i < N; ++i)
	  a_[i] -= D[i];
      return *this;
    }

    //!unary minus: i.e. u can use -D instead of D *= -1 which faciliates code
    FieldDiagonal& operator-(){
      return ((*this) *= -1);
    }



    iterator begin(){
      return &a_[0];
    }

    iterator end(){
      return &a_[N];
    }
    
    const_iterator begin() const{
      return &a_[0];
    }

    const_iterator end() const{
      return &a_[N];
    }

    //determinant -- just the product of the diagonal entries
    T determinant() const{
      T det = 1;
      for(int i = 0; i < N; ++i)
	det *= a_[i];

      return det;
    }

    inline bool is_singular() const{
      if(Abs(determinant()) < 0)
	return true;
      else 
	return false;
    }
    
    //inverse -- just the reciprocal values of the diagonal 
    void invert(){
      ThisType A(*this);
      *this = T();          //reset for holding the inverse
      for(int i = 0; i < N; ++i){
	adonis_assert(!is_zero(A[i]));  //assert that vals are above e.g. 1.e-13
	a_[i] = static_cast<T>(1)/A[i];
      }
    }

    //a somehow useless operation since the transpose of a diagonal matrix is the diagonal matrix itself.
    FieldDiagonal& transpose(){
      return *this;
    }

    //create a dense field matrix out of a diagonal matrix (seems to be a nuisance but anyway, perhaps someone might find a use for it ;)
    void make_dense_matrix(Dune::FieldMatrix<T,N,N>& DenseDiag) const{
      DenseDiag = T();
      for(int i = 0; i < N; ++i)
	DenseDiag[i][i] = a_[i];
    }

    //output via <<
    friend std::ostream& operator<<(std::ostream& os, const FieldDiagonal& D){
      os << "{ ";
      for(const_iterator it = D.begin(); it != D.end(); ++it)
	os << *it << "  ";
      os << "}"<<std::endl;
      return os;
    }

    inline void identity(){    //the identity matrix
      for(int i = 0; i < N; ++i)
	a_[i] = 1;
    }

  private:
    T a_[N];
  };




  //some operations -- these could be easily reformed as expression templates. 
  //However, this would need a definition of an additional 
  //operator=(Expression e) within the class FieldMatrix in 'fmatrix.hh'....
  
  /**
   *\brief Left addition of FieldDiagonal with <B>dense</B> FieldMatrix
   */
  template<class T, int N>
  inline Dune::FieldMatrix<T,N,N> operator+(const FieldDiagonal<T,N>& D, 
					    const Dune::FieldMatrix<T,N,N>& A){
   
    typedef Dune::FieldMatrix<T,N,N> DenseType;
    DenseType Add(A);
    
    for(int i = 0; i < N; ++i){
      Add[i][i] += D[i];
    }
    return Add;
  }
  
  template<class T, int N>
  inline Dune::FieldMatrix<T,N,N> operator+(const Dune::FieldMatrix<T,N,N>& A,
					    const FieldDiagonal<T,N>& D){
    return D + A;  //same as above   
  }

  //some operations 
  /**
   *\brief Left subtraction of FieldDiagonal with <B>dense</B> FieldMatrix
   */
  template<class T, int N>
  inline Dune::FieldMatrix<T,N,N> operator-(const FieldDiagonal<T,N>& D, 
					    const Dune::FieldMatrix<T,N,N>& A){
   
    typedef Dune::FieldMatrix<T,N,N> DenseType;
    DenseType Subtract(A);

    for(int i = 0; i < N; ++i){
      for(int j = 0; j < N ; ++j){
	if(i == j)
	  Subtract[i][j] = D[i]-A[i][j];
	else
	  Subtract[i][j] *= -1; 
      }
    }
    return Subtract;
  }

  template<class T, int N>
  inline Dune::FieldMatrix<T,N,N> operator-(const Dune::FieldMatrix<T,N,N>& A,
					    const FieldDiagonal<T,N>& D){
   
    typedef Dune::FieldMatrix<T,N,N> DenseType;
    DenseType Subtract(A);
    
    for(int i = 0; i < N; ++i){
      Subtract[i][i] -= D[i];
    }
    return Subtract;
  }
  
  
  /**
   *\brief Left multiplication of FieldDiagonal with <B>dense</B> FieldMatrix
   */
  template<class T, int M, int N>
  inline Dune::FieldMatrix<T,M,N> operator*(const FieldDiagonal<T,M>& D, 
					    const Dune::FieldMatrix<T,M,N>& A){
   
    typedef Dune::FieldMatrix<T,M,N> DenseType;
    DenseType Prod;
    
    for(int i = 0; i < M; ++i){
      for(int j = 0; j < N; ++j){
	Prod[i][j] = D[i]*A[i][j];
      }
    }
    return Prod;
  }
  
  template<class T, int M, int N>
  inline Dune::FieldMatrix<T,M,N> operator*(const Dune::FieldMatrix<T,M,N>& A,
					    const FieldDiagonal<T,N>& D){
   
    typedef Dune::FieldMatrix<T,M,N> DenseType;
    DenseType Prod;
    
    for(int i = 0; i < M; ++i){
      for(int j = 0; j < N; ++j){
	Prod[i][j] = D[j]*A[i][j];
      }
    }
    return Prod;
  }
  


  /**
   * \brief Multiplication with vector from right
   */
  /*template<int N, class V>
  inline V& operator*(const FieldDiagonal<typename V::value_type,N>& D, V& v){
    adonis_assert((int)v.size() == N);
    for(int i = 0; i < N; ++i)
      v[i] *= D[i];
    return v;
  }

  
   //!brief Multiplication with vector from left -- the same as above (it is the transposed vector, strictly speaking
  template<class V,int N>
  inline V& operator*(V& v, const FieldDiagonal<typename V::value_type,N>& D){
    adonis_assert((int)v.size() == N);
    for(int i = 0; i < N; ++i)
      v[i] *= D[i];
    return v;
  }
  */
  /////////////
  
  /**
   * \brief Multiplication Diagonal with vector
   */
  template<int N, class V>
  inline V operator*(const FieldDiagonal<typename V::value_type,N>& D, const V& x){
    V v(N);
    for(int i = 0; i < N; ++i)
      v[i] = D[i]*x[i]; 
    return v;
  }

  template<int N, class V>
  inline V operator*(const V& x, const FieldDiagonal<typename V::value_type,N>& D){
    return D*x;
  }
  
  template<int N, class T>
  inline Dune::FieldVector<T,N> operator*(const FieldDiagonal<T,N>& D, const Dune::FieldVector<T,N>& x){
    Dune::FieldVector<T,N> v;
    for(int i = 0; i < N; ++i)
      v[i] = D[i]*x[i]; 
    return v;
  }

 template<int N, class T>
 inline Dune::FieldVector<T,N> operator*(const Dune::FieldVector<T,N>& x, const FieldDiagonal<T,N>& D){
   return D*x;
 }


  /**
   * \brief Scale system prior to Gaussian elimination (based on LU) decomposition, e.g. 
   * [Higham, 2nd, § 9.8] or [Golub/vanLoan, 3rd, § 3.5.2]
   */
  template<class FD1, class FD2, class LHS, class RHS>
  inline void scaling_prior_to_GE(const FD1& D1, const FD2& D2, LHS& A, RHS& B){
    FD2 D2_Inv(D2);

    D2_Inv.invert();
    
    A = D1*A*D2*D2_Inv;
    B = D1*B;
  }


}//end of namespace 

#endif
