#ifndef A_TYPE_ADAPTER_WHICH_CAN_BE_EASILY_EMPLOYED_HH
#define A_TYPE_ADAPTER_WHICH_CAN_BE_EASILY_EMPLOYED_HH

#include <iostream>
#include <sstream>
#include <complex>
#include <string>
#include <typeinfo>

#include "isclass.hh"
#include "adonisassert.hh"
#include "error.hh"

#if USE_CPPAD
#include <cppad/cppad.hpp>
#endif 

namespace Adonis{
  
  //! Check if typee is a class or not
  template<class X, bool B> class WhatType;

  //! is a class, e.g. std::complex
  template<class X>
  class WhatType<X,true>{
  public:
    typedef X Type;
    typedef typename X::value_type BaseType;
  };

  //! is not class, e.g. long double
  template<class X>
  class WhatType<X,false>{ 
  public:
    typedef X Type;
    typedef X BaseType;
  };
  

  /**
   * \brief Define iterator type
   * USAGE:
   * \code
      typedef std::complex<long double> some_type;
      some_type cl(3.000000000000000000, 4.0000000000000000000000000);

      for(WhatIteratorType<some_type, IsContainer<some_type>::Value>::const_iterator it = &cl; it != &cl+1; ++it){
      cout<< (*it).real() << "   "<<(*it).imag() <<endl;
    }
   * \endcode
   */
  template<class X, bool B> class WhatIteratorType;
  
  //! is container
  template<class X>
  class WhatIteratorType<X,true>{
  public:
    typedef typename X::iterator iterator;
    typedef typename X::const_iterator const_iterator;
  };


  //! non container
  template<class X>
  class WhatIteratorType<X,false>{
  public:
    typedef X* iterator;
    typedef const X* const_iterator;
  };


  /**
   * \brief Every type has now a type and a base type. Those are the same in case that non-class types are used. Otherwise the base type is the <TT> value_type </TT> of the class, assumed of course the class implemention provides a value_type (which is true for every foresightedly implemented class)!!
   */
  template<class T> 
  class TypeAdapter{
  public:
    typedef typename WhatType<T,IsClass<T>::Value>::Type Type;
    typedef typename WhatType<T,IsClass<T>::Value>::BaseType BaseType;
  };
  

#if USE_CPPAD
  template<class T>
  class TypeAdapter<CppAD::AD<T> >{
  public:
    //!OLD SETTING
    typedef T Type;  //!note this is already the AD template type !!
    typedef typename TypeAdapter<T>::BaseType BaseType;
    // typedef typename TypeAdapter<T>::BaseType Type;
    // typedef typename TypeAdapter<Type>::BaseType BaseType;
    //typedef CppAD::AD<T> BaseType;
  };
#endif

  //

  /**
   * \brief typedef for constant refs
   */
  template<class T>
  class ConstReference{
  public:
    typedef const typename  TypeAdapter<T>::Type Type;
    typedef const typename TypeAdapter<T>::BaseType BaseType;
  };
  
  template<class T>
  class ConstReference<const T>{
  public:
    typedef const typename TypeAdapter<T>::Type Type;
    typedef const typename TypeAdapter<T>::BaseType BaseType;
  };
  
  template<class T>
  class ConstReference<T&>{
  public:
    typedef const typename TypeAdapter<T>::Type& Type;
    typedef const typename TypeAdapter<T>::BaseType& BaseType;
  };

  template<class T>
  class ConstReference<const T&>{
  public:
    typedef const typename TypeAdapter<T>::Type& Type;
    typedef const typename TypeAdapter<T>::BaseType& BaseType;
  };


  /**
   * \brief typedefs for (mutable) references
   */
  template<class T>
  class Reference{
  public:
    typedef typename TypeAdapter<T>::Type Type;
    typedef typename TypeAdapter<T>::BaseType BaseType;
  };
  
  template<class T>
  class Reference<const T>{
  public:
    typedef typename TypeAdapter<T>::Type Type;
    typedef typename TypeAdapter<T>::BaseType BaseType;
  };
  
  template<class T>
  class Reference<T&>{
  public:
    typedef typename TypeAdapter<T>::Type& Type;
    typedef typename TypeAdapter<T>::BaseType& BaseType;
  };

  template<class T>
  class Reference<const T&>{
  public:
    typedef typename TypeAdapter<T>::Type& Type;
    typedef typename TypeAdapter<T>::BaseType& BaseType;
  };


  /**
   * \brief Non-container adapter for scalars 
   */
  template<class T>
  class ConstScalarAdapter{
  public:
    typedef T value_type;
    typedef const T* const_iterator;

    ConstScalarAdapter(const T& scal = T()):scal_(scal){}

    const_iterator begin() const{
      return &scal_;
    }

    const_iterator end() const{
      return &scal_+1;  //plus 1
    }
   
  private:
    const T& scal_;

  };


  /**
   * \brief mutable container adapter for scalars
   */
  template<class T>
  class MutableScalarAdapter{
  public:
    typedef T value_type;
    typedef T* iterator;

    MutableScalarAdapter(T& scal = T()):scal_(scal){}

    iterator begin(){
      return &scal_;
    }

    iterator end(){
      return &scal_+1;  //plus 1
    }
   
  private:
     T& scal_;

  };

  /**
   * \brief output types which I use for my class to read in non-compile time parameters, namely with class <TT> ParameterData </TT>
   */
  template<class T> 
  class OutputTypeAdapter{  //! string to number 
  public:
    typedef typename TypeAdapter<T>::Type Type;   //for CppAD this is the num. type, for double it's just double


    static inline Type convert(const std::string& ds){
      adonis_assert(ds.size() > 0); //otherwise you don't need to read in
      Type number;
      std::stringstream ss(ds);
      ss >> number;
      return number;
    }
  };
 

// #if USE_CPPAD
// template<class T> 
// class OutputTypeAdapter<CppAD::AD<T> >{  //! string to number 
//   public:
//     typedef T Type; 


//     static inline Type convert(const std::string& ds){
//       adonis_assert(ds.size() > 0); //otherwise you don't need to read in
//       Type number;
//       std::stringstream ss(ds);
//       ss >> number;
//       return number;
//     }
//   };
// #endif


  template<> 
  class OutputTypeAdapter<std::string>{ //! just return the string
  public:
    static inline const std::string& convert(const std::string& ds){
      adonis_assert(ds.size() > 0); //otherwise you don't need to read in
      return ds;
    }
  };

  
  template<> 
  class OutputTypeAdapter<char>{  //!string only possesses one character
  public:
    static inline char convert(const std::string& ds){
      adonis_assert(ds.size() > 0); //otherwise you don't need to read in
      return ds[0];
    }
  };
  

  template<> 
  class OutputTypeAdapter<bool>{ //!convert string to boolean value
  public:
    static inline bool convert(const std::string& ds){
      adonis_assert(ds.size() > 0); //otherwise you don't need to read in
      return ((ds == "true" || ds == "1") ? true : false);
    }
  };

 

  /**
   * \brief Checks if a given type is a container or not. The first one is by definition no basic data type and an error will be thrown
   */
  template<class C, bool B> class IsBasicDataType;

  template<class C>
  class IsBasicDataType<C,true>{
  public:
    static const bool Value = true;
  
    static inline void check(){}
  };

  template<class C>
  class IsBasicDataType<C,false>{
  public:
    static const bool Value = false;
  
    static inline void check(){
      ADONIS_ERROR(TypeError, "Container \"" << typeid(C).name() <<"\" is not a basic data type.");
    }
  };


} //end namespace

#endif
