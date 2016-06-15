#ifndef WRAPPER_TRAITS_4_LAPACK_STUFF_HH
#define WRAPPER_TRAITS_4_LAPACK_STUFF_HH

/** \brief These traits prevent you from inputting unsupported BASIC types 
 *    when using FORTRAN 77.
 *   Unsupported types such as 'short' and 'int', for instance, will be 
 *   interpreted  as the next "higher" type which -in the case of the example- 
 *   will be 'float'
 *   Implementation is done via Value Traits.
*/


#if USE_LAPACK
#ifdef LAPACK_USAGE


#include<complex>

#if USE_CPPAD
#include<cppad/cppad.hpp>
#endif

#include "../common/typeadapter.hh" //used in conjunction with CppAD


namespace Adonis{

  template<typename T>
  class TypeTraits{
  public:
    typedef typename TypeAdapter<T>::Type Type;
    typedef typename TypeAdapter<T>::BaseType BaseType;
  };


  //specialisations
  //###### Real ######
  template<>
  class TypeTraits<short>{
  public:
    typedef float Type;
    typedef float BaseType;
    static const char Value = 's'; //the constant static data member
  };                               //NOTE: Assignment of a value does not work 
                                   //when you want to use, say, 'double' since 
                                   //this is not a const static data member 
                                   //like int, bool and char!

  template<>
  class TypeTraits<long>{
  public:
    typedef float Type;
    typedef float BaseType;
    static const char Value = 's';
  };

  template<>
  class TypeTraits<bool>{
  public:
    typedef float Type;
    typedef float BaseType;
    static const char Value = 's';
  };


  template<>
  class TypeTraits<int>{
  public:
    typedef float Type;
    typedef float BaseType;
    static const char Value = 's';
  };

  template<>
  class TypeTraits<unsigned int>{   //this is the size_t case
  public:
    typedef float Type;
    typedef float BaseType;
    static const char Value = 's';
  };

  

  template<>
  class TypeTraits<float>{
  public:
    typedef float Type;
    typedef float BaseType;
    static const char Value = 's';
  };


  template<>
  class TypeTraits<char>{
  public:
    typedef float Type;
    typedef float BaseType;
    static const char Value = 's';
  };



  //the double precision case
  template<>
  class TypeTraits<double>{
  public:
    typedef double Type;
    typedef double BaseType;
    static const char Value = 'd';
  };

  
  template<>
  class TypeTraits<long double>{
  public:
    typedef double Type;
    typedef double BaseType;
    static const char Value = 'd';
  };


  


  //###### Complex ########

  template<>
  class TypeTraits<std::complex<short> >{
  public:
    typedef std::complex<float> Type;
    typedef float BaseType;
    static const char Value = 'c';
  };
 
  template<>
  class TypeTraits<std::complex<bool> >{
  public:
    typedef std::complex<float> Type;
    typedef float BaseType;
    static const char Value = 'c';
  };

  template<>
  class TypeTraits<std::complex<long> >{
  public:
    typedef std::complex<float> Type;
    typedef float BaseType;
    static const char Value = 'c';
  };

  template<>
  class TypeTraits<std::complex<int> >{
  public:
    typedef std::complex<float> Type;
    typedef float BaseType;
    static const char Value = 'c';
  };

  template<>
  class TypeTraits<std::complex<unsigned int> >{ //this is the size_t case
  public:
    typedef std::complex<float> Type;
    typedef float BaseType;
    static const char Value = 'c';
  };

  

  template<>
  class TypeTraits<std::complex<float> >{
  public:
    typedef std::complex<float> Type;
    typedef float BaseType;
    static const char Value = 'c';
  };

  template<>
  class TypeTraits<std::complex<char> >{
  public:
    typedef std::complex<float> Type;
    typedef float BaseType;
    static const char Value = 'c';
   };



  //the double precision case
  template<>
  class TypeTraits<std::complex<double> >{
  public:
    typedef std::complex<double>  Type;
    typedef double BaseType;
    static const char Value = 'z';
  };

  
  template<>
  class TypeTraits<std::complex<long double> >{
  public:
    typedef std::complex<double> Type;
    typedef double BaseType;
    static const char Value = 'z';
  };


#if USE_CPPAD
  /**
   * \brief Specialisation for CppAD types
   */
  template<class T>
  class TypeTraits<CppAD::AD<T> >{
  public:
    typedef typename TypeAdapter<typename CppAD::AD<T>::value_type>::Type GenType; //can still be something like long double or std::complex<long double>
    typedef typename TypeTraits<GenType>::Type Type;
    typedef typename TypeTraits<GenType>::BaseType BaseType;
  
    static const char Value = TypeTraits<T>::Value;
  };
#endif //USE_CPPAD


}//end namespace


#endif //LAPACK_USAGE
#endif 


#endif
