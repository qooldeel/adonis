#ifndef UNIVERSAL_CONSTANTS_THAT_CAN_BE_USED_IN_ANY_CONTEXT_HH
#define UNIVERSAL_CONSTANTS_THAT_CAN_BE_USED_IN_ANY_CONTEXT_HH

#include "typeadapter.hh"

namespace Adonis{

  /**
   * \brief Universally applicable constants, including mathematical constants
   * such as $\mathit{e}.$ All constants are given in their decimal representation
   */
  template<class T>
  class UniversalConstants{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef int IntegerType;
    typedef std::size_t SizeType;
    typedef unsigned UnsignedIntegerType;

    static const T enlargeRange;
    static const T aboutZero;
    static const T aboutInfinity;
    static const T smallValue;
    static const T smallLow;
    static const T smallUp;

    //!mathematical constants
    static const T Pi;
    static const T EulersNumber;
    static const T EulerMascheroni;
    static const T goldenRatio;
  };

  //! a given range [a,b] can be enlarged by a prescribed value,
  //! i.e. [a-eps,b+eps]
  template<class T> const T UniversalConstants<T>::enlargeRange = 1.e-07;
  
  //! define non-integral types right here -- and <B> only </B> here
  template<class T> const T UniversalConstants<T>::aboutZero = 1.e-13;
  template<class T> const T UniversalConstants<T>::aboutInfinity = 1.e+28;
  template<class T> const T UniversalConstants<T>::smallValue = 1.e-16;
  template<class T> const T UniversalConstants<T>::smallLow = 1.e-11;
  template<class T> const T UniversalConstants<T>::smallUp = 1.e-10;
  
  template<class T> const T UniversalConstants<T>::Pi = 3.1415926535897932384626433832795028841971;
  template<class T> const T UniversalConstants<T>::EulersNumber = 2.7182818284590452353602874713526624977572;                                   
   template<class T> const T UniversalConstants<T>::EulerMascheroni = 0.5772156649015328606065120900824024310421;
  template<class T> const T UniversalConstants<T>::goldenRatio = 1.6180339887498948482;
  

} //end namespace 

#endif
