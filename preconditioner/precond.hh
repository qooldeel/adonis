#ifndef PRECONDITIONERS_TO_BE_USED_TO_SPEED_UP_COMPUTATIONS_HH
#define PRECONDITIONERS_TO_BE_USED_TO_SPEED_UP_COMPUTATIONS_HH

#include "../common/globalfunctions.hh"
#include "../misc/useful.hh"
#include "../misc/operations4randomaccesscontainers.hh"

namespace Adonis{

  template<class T>
  class PrecConsts{
  public:
    static const T omega;        //used in conjunction with the SSOR precond.

  };
  template<class T> const T PrecConsts<T>::omega = 1.3; //1.2 //1.3

  /**
   * \brief Some preconditioner
   */
  template<int PC, char MT> class Preconditioner;

  template<>
  class Preconditioner<1,'f'>{  //1 = SSOR for full matrix
  public:
    template<class V>
    static inline V prec(const V& A, const V& r){
      typedef typename V::value_type value_type;
      
      //! CAUTION it is very important to use ints here !!!!!!!!!!!!
      int n = static_cast<int>(r.size());
      V z(n);
      
      value_type sum;
      for(int i = 0; i < n; ++i){
	sum = value_type();
	for(int j = 0; j < i; ++j){
	  sum += A[RowMajor::offset(i,j,n)]*z[j];
	}
	adonis_assert(!is_zero(A[RowMajor::offset(i,i,n)]));
	z[i] = (r[i] - PrecConsts<value_type>::omega*sum)/A[RowMajor::offset(i,i,n)];
      }
      
      for(int i = n-1; i >= 0; --i){
	sum = value_type();
	for(int j = i+1; j < n; ++j){
	  sum += A[RowMajor::offset(i,j,n)]*z[j];
	}
	adonis_assert(!is_zero(A[RowMajor::offset(i,i,n)]));
	z[i] -= PrecConsts<value_type>::omega*sum/A[RowMajor::offset(i,i,n)];
      }
      return z;
    }
  };

  
  //! it is assumed that the dense symmetric matrix A is stored in lower 
  //! triangular form (or upper triangular form).   
  template<>
  class Preconditioner<1,'y'>{  //1 = SSOR for symmetric dense matrix
  public:
    template<class V>
    static inline V prec(const V& A, const V& r){
      //std::cout << "SSOR for symmetric matrices"<<std::endl;
      adonis_assert(A.size() == gauss_sum(r.size()));
      typedef typename V::value_type value_type;
      
      //! CAUTION it is very important to use ints here !!!!!!!!!!!!
      int n = static_cast<int>(r.size());
      V z(n);
      
      value_type sum;
      for(int i = 0; i < n; ++i){
	sum = value_type();
	for(int j = 0; j < i; ++j){
	  sum += A[SymmetricAccess::offset(i,j,n)]*z[j];
	}
	adonis_assert(!is_zero(A[SymmetricAccess::offset(i,i,n)]));
	z[i] = (r[i] - PrecConsts<value_type>::omega*sum)/A[SymmetricAccess::offset(i,i,n)];
      }
      
      for(int i = n-1; i >= 0; --i){
	sum = value_type();
	for(int j = i+1; j < n; ++j){
	  sum += A[SymmetricAccess::offset(i,j,n)]*z[j];
	}
	adonis_assert(!is_zero(A[SymmetricAccess::offset(i,i,n)]));
	z[i] -= PrecConsts<value_type>::omega*sum/A[SymmetricAccess::offset(i,i,n)];
      }
      return z;
    }
  };

  

}

#endif
