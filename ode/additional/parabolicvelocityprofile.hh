#ifndef SPECIFY_A_PARABOLIC_INLET_VELOCITY_PROFILE_HH
#define SPECIFY_A_PARABOLIC_INLET_VELOCITY_PROFILE_HH

#include "../../common/globalfunctions.hh"
#include "../../common/numericaltraits.hh"

namespace Adonis{


  template<class T1, class T2, class T3>
  inline typename NumericalTypePromotion<T1,typename NumericalTypePromotion<T2,T3>::ReturnType>::ReturnType parabolic_inlet_velocity(const T1& r, const T2& R, const T3& v1max){
    return v1max*(1. - r*r/(R*R));
  }


  /**
   * \brief Compile time selection of various velocity profiles
   */
  template<class T1, class T2, class T3, char C> class InletVelocity;

  //! partial specializations
  template<class T1, class T2, class T3>
  class InletVelocity<T1,T2,T3,'p'>{
  public:
    typedef typename NumericalTypePromotion<T1,typename NumericalTypePromotion<T2,T3>::ReturnType>::ReturnType ReturnType;

    static const char Value = 'p';

    static inline ReturnType velocity(const T1& r, const T2& R, const T3& v1max){
      return zero(parabolic_inlet_velocity(r,R,v1max));
    }

  };

  template<class T1, class T2, class T3>
  class InletVelocity<T1,T2,T3,'P'>{
  public:
    typedef typename NumericalTypePromotion<T1,typename NumericalTypePromotion<T2,T3>::ReturnType>::ReturnType ReturnType;

    static const char Value = 'P';

    static inline ReturnType velocity(const T1& r, const T2& R, const T3& v1max){
      return InletVelocity<T1,T2,T3,'p'>::velocity(r,R,v1max);
    }

  };
  
  //! constant inlet velocity profile
  template<class T1, class T2, class T3>
  class InletVelocity<T1,T2,T3,'c'>{
  public:
    typedef typename NumericalTypePromotion<T1,typename NumericalTypePromotion<T2,T3>::ReturnType>::ReturnType ReturnType;

    static const char Value = 'c';

    static inline const T3& velocity(const T1& r, const T2& R, const T3& v1max){
      return v1max;  //the first two arguments are just dummies
    }

  };
  
  template<class T1, class T2, class T3>
  class InletVelocity<T1,T2,T3,'C'>{
  public:
    typedef typename NumericalTypePromotion<T1,typename NumericalTypePromotion<T2,T3>::ReturnType>::ReturnType ReturnType;

    static const char Value = 'C';

    static inline const T3& velocity(const T1& r, const T2& R, const T3& v1max){
      return v1max;  //the first two arguments are just dummies
    }

  };
  

}  //end namespace 

#endif
