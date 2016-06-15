#ifndef SMART_ASSIGN_CAN_BE_HELPFUL_HH
#define SMART_ASSIGN_CAN_BE_HELPFUL_HH

#include "../expressiontemplates/exprvec.hh"
#include "adonisassert.hh"

#include "globalfunctions.hh"
#include "typeadapter.hh"


#if USE_CPPAD
#include <cppad/cppad.hpp>
#endif 

namespace Adonis{
  
 
  //! assignment of to unequal types -- assumes that these values can be 
  //! converted by an operator=
  template<class T1, class T2>
  inline void smart_assign(T1& t1, const T2& t2){
    t1 = t2;
  }
 
  /**
   * \brief smart assign functions: if second argument is a CppAD::AD then retrieve its value and assign it to the first value of type T
   */
#if USE_CPPAD
  template<class T>
  inline void smart_assign(T& t1, const CppAD::AD<T>& ad){
    t1 = convert_fancy_num_2_num(ad); //CppAD::Value(ad); 
  }
#endif

  
} //end namespace 

#endif 
