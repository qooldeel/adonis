#ifndef UMFPACK_TYPE_TRAITS_HH
#define UMFPACK_TYPE_TRAITS_HH

#include <iostream>
#include <complex>

#include "../common/globalfunctions.hh"
#include "../common/numerictypechecker.hh"

#if USE_UMFPACK
#include "umfpack.h"
#endif

namespace Adonis{

  /**
   * UMFPACK promotes double, int and SuiteSparse_long integers (di, dl, zi, zl).
   * Hence we must guarantee that all other types are up-/or downgraded to the
   * aforementiones types
   */
  template<class T>    //! e.g. when a long double is used
  class UmfpackTypeTraits{  
  public:
    typedef double Type;
    typedef double BaseType;
  };

  //!specializations
  template<class T>
  class UmfpackTypeTraits<std::complex<T> >{
  public:
    typedef std::complex<double> Type;
    typedef double BaseType;
  };

  //!INDICIES: INTEGRAL TYPES, only int and SuiteSparse_long are permitted
  //! note that std::size_t is just an alias for int (32bit) and long (64bit)
  template<>
  class UmfpackTypeTraits<short>{
  public:
    typedef int Type;
    typedef int BaseType;
  };

  template<>
  class UmfpackTypeTraits<int>{
  public:
    typedef int Type;
    typedef int BaseType;
  };
  
  template<>
  class UmfpackTypeTraits<unsigned int>{
  public:
    typedef int Type;
    typedef int BaseType;
  };

  template<>
  class UmfpackTypeTraits<bool>{
  public:
    typedef int Type;
    typedef int BaseType;
  };

  template<>
  class UmfpackTypeTraits<char>{
  public:
    typedef int Type;
    typedef int BaseType;
  };

 
  ////! Apparently, <TT> SuiteSparse_long </TT>  and <TT> long </TT>  are the
  ////! same types 
// #if USE_UMFPACK
//   template<>
//   class UmfpackTypeTraits<SuiteSparse_long>{
//   public:
//     typedef SuiteSparse_long Type;
//     typedef SuiteSparse_long BaseType;
//   };
// #endif
 
  template<>
  class UmfpackTypeTraits<long>{
  public:
#if USE_UMFPACK
    typedef SuiteSparse_long Type;
    typedef Type BaseType;
#else
    typedef int Type;
    typedef int BaseType;
#endif
  };

//! INDICIES: COMPLEX INTEGRAL TYPES
  template<>
  class UmfpackTypeTraits<std::complex<short> >{
  public:
    typedef int Type;
    typedef int BaseType;
  };

  template<>
  class UmfpackTypeTraits<std::complex<int> >{
  public:
    typedef int Type;
    typedef int BaseType;
  };
  
  template<>
  class UmfpackTypeTraits<std::complex<unsigned int> >{
  public:
    typedef int Type;
    typedef int BaseType;
  };

  template<>
  class UmfpackTypeTraits<std::complex<bool> >{
  public:
    typedef int Type;
    typedef int BaseType;
  };

  template<>
  class UmfpackTypeTraits<std::complex<char> >{
  public:
    typedef int Type;
    typedef int BaseType;
  };

 
// #if USE_UMFPACK
//   template<>
//   class UmfpackTypeTraits<std::complex<SuiteSparse_long> >{
//   public:
//     typedef SuiteSparse_long Type;
//     typedef SuiteSparse_long BaseType;
//   };
// #endif


  template<>
  class UmfpackTypeTraits<std::complex<long> >{
  public:
#if USE_UMFPACK
    typedef SuiteSparse_long Type;
    typedef Type BaseType;
#else
    typedef int Type;
    typedef int BaseType;
#endif
  };



  /**
   * \brief When complex numbers are used, UMFPACK requires that the real parts
   * and the complex parts of a complex array A are passed by 2 double arrays, 
   * Ax and Az, respectively, the first storing the real parts, the latter one 
   * the imaginary parts of A.
   *
   * NOTE: no range check is performed, only iterators are used and it is 
   * assumed that a []-operator is available
   */
  template<bool B> class UmfpackUsesComplexNumbers;

  //! specializations
  //! complex numbers can be detected via <TT> NumericDataTypeChecker<T>::IsComplex </TT>: returns true or false whether a complex number is encountered or not 
  template<>
  class UmfpackUsesComplexNumbers<true>{  //COMPLEX numbers are used
  public:
    template<class INT, class V>
    static inline void resize(INT length, V& v, const typename V::value_type& val = typename V::value_type()){  
      INT dim = static_cast<INT>(std::distance(v.begin(),v.end()));
      //std::cout << "----------- RESIZE "<<std::endl;

      (length > 0 && dim != length) ? v.resize(length,val) : do_nothing();
    }

    template<class INT, class ITZ, class ITD1, class ITD2>
    static inline void split_complex_container_into_real_values(INT length, ITZ z, ITD1 re, ITD2 im){
      for(INT i = 0; i < length; ++i){
	re[i] = z[i].real();
	im[i] = z[i].imag();
      }
    }
  
    template<class INT, class ITZ, class ITD1, class ITD2>
    static inline void join(INT length, ITZ z, ITD1 re, ITD2 im){
      for(INT i = 0; i < length; ++i){
	z[i].real() = re[i];
	z[i].imag() = im[i];
      }
    }

    template<class ITZ, class ITRE>
    static inline ITRE proper_iterator(ITZ z, ITRE re){
      return re;
    }
  };


  template<>
  class UmfpackUsesComplexNumbers<false>{ // NON-complex numbers   
  public:
    template<class INT, class V> //do nothing
    static inline void resize(INT length, V& v, const typename V::value_type& val = typename V::value_type() ){}

    template<class INT, class ITZ, class ITD1, class ITD2>
    static inline void split_complex_container_into_real_values(INT length, ITZ z, ITD1 re, ITD2 im){}  //do nothing at all
  
    template<class INT, class ITZ, class ITD1, class ITD2>  
    static inline void join(INT length, ITZ z, ITD1 re, ITD2 im){
        for(INT i = 0; i < length; ++i){
	  z[i] = re[i];
	}
    }

    template<class ITZ, class ITRE>
    static inline ITZ proper_iterator(ITZ z, ITRE re){
      return z;
    }
  };

  
 
} //end namespace 

#endif
