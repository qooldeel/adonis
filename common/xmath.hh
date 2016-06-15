#ifndef EXTENDED_MATH_COLLECTION_HH
#define EXTENDED_MATH_COLLECTION_HH

#if USE_CPPAD
#include <cppad/cppad.hpp>
#endif

#include "error.hh"

namespace Adonis{



  /**
   * \brief Throw error when real square root is performed on negative value
   */
  inline double Sqrt(const double& val){
    if(val < 0 )
      ADONIS_ERROR(ValueError,"val is NEGATIVE!");
    return std::sqrt(val);

  }

    

  inline float Sqrt(const float& val){
     if(val < 0 )
      ADONIS_ERROR(ValueError,"val is NEGATIVE!");
    return std::sqrt(val);

  }

    

  inline long double Sqrt(const long double& val){
     if(val < 0 )
      ADONIS_ERROR(ValueError,"val is NEGATIVE!");
    return std::sqrt(val);

  }

    

  //T is any other type

  template<class T>

  inline double Sqrt(const T& val){
     if(val < 0 )
      ADONIS_ERROR(ValueError,"val is NEGATIVE!");
    return std::sqrt(val);

  }

    

  template<class T>

  inline std::complex<T> Sqrt(const std::complex<T>& z){
    return std::sqrt(z);
  }

    

     

#if USE_CPPAD

  inline CppAD::AD<double> Sqrt(const CppAD::AD<double>& val){
    if(val < 0 )
      ADONIS_ERROR(ValueError,"val is NEGATIVE!");
    return CppAD::sqrt(val);

  }

    

  inline CppAD::AD<float> Sqrt(const CppAD::AD<float>& val){
    if(val < 0 )
      ADONIS_ERROR(ValueError,"val is NEGATIVE!");
    return CppAD::sqrt(val);

  }

    

  // inline CppAD::AD<long double> Sqrt(const CppAD::AD<long double>& val){
  //   if(val < 0 )
  //     ADONIS_ERROR(ValueError,"val is NEGATIVE!");
  //   return CppAD::sqrt(val);

  // }

    

  // //!T is any other type -- this is not working when T != double
  // template<class T>
  // inline CppAD::AD<double> Sqrt(const CppAD::AD<T>& val){
  //   if(val < 0 )
  //     ADONIS_ERROR(ValueError,"val is NEGATIVE!");
  //   return CppAD::sqrt(val);

  // }

    

  template<class T>

  inline CppAD::AD<std::complex<T> > Sqrt(const CppAD::AD<std::complex<T> >& z){
    return CppAD::sqrt(z);
  }

    

#endif

    

} //end namespace

    

#endif 
