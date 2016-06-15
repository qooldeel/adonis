#ifndef UTILITIES_FOR_DIFFERENTIAL_EQUATIONS_HH
#define UTILITIES_FOR_DIFFERENTIAL_EQUATIONS_HH

#include <cmath>
#include "../common/adonisassert.hh"
#include "../common/isclass.hh"
#include "../misc/operations4randomaccesscontainers.hh"


namespace Adonis{

  /**
   * \brief TMP to select whether one regards containers or types
   */
  template<class X, bool B> class SpecificCalculationOfHNew;

  template<class X>
  class SpecificCalculationOfHNew<X,true>{ //is container
  public:
    static inline typename X::value_type will_be_performed(const X& x, const typename X::value_type& h, unsigned p, const typename X::value_type& tol) {
      return 0.9*h*std::pow(tol/(Norm<'i',typename X::value_type>::norm(x)), 1./(1+p));
    }
  };

  template<class T>
  class SpecificCalculationOfHNew<T,false>{ //a numeric type
  public:
    static inline T will_be_performed(const T& x, const T& h, unsigned p, const T& tol) {
      return 0.9*h*std::pow(tol/(Abs(x)), 1./(1+p));
    }
  };

  
  /**
   * \brief The formula \f[ h_{\textrm{new}} := 0.9\cdot h \left( \frac{\textrm{tol}}{\| e\|_{\infty}}\right)^{\frac{1}{1+p}}\f]
   * is used for both increasing and decreasing the (old) stepsize \f$ h.\f$
   * In the above equation \f$ e \f$ respresents an estimate of the local trunction error, tol the local error tolerance and \f$ p\f$ the order of the method. If two methods with different order are employed, e.g. as in the case of the adaptive RKF45 then \f$ p = 4, \f$ i.e. takes the value of the first formula in the pair of RK formulae (see [KINCAID,CHENEY, <I>"Numerical Analysis"</I>, ch. 8, p. 548]). I think it takes the order of the method, hence 5.
   * This formula can be found in
   * [1] [PRINCE/DORMAND, <I>"Higher order embedded Runge-Kutta formulae"</I>, eq. (1.3)]
   * [2] [SHAMPINE/WATTS, <I>"Global Error Estimation for Ordinary Differential Equations"</I>, sect. 4, p. 178] -- here the negative power is employed
   [3] [HULL/ENRIGHT/FELLEN/SEDGWICK, <I>"COMPARING NUMERICAL METHODS FOR ORDINARY DIFFERENTIAL EQUATIONS"</I>, p. 622] -- here the formula is slightly different, for it takes into account the index of the current iteration 
   *
   * NOTE: this formula is used for <B> both increasing and decreasing </B> the step size. 
  */
  template<class V, class T>
  inline T h_new(const V& x, const T& h, unsigned p, const T& tol){
    return SpecificCalculationOfHNew<V,IsContainer<V>::Value>::will_be_performed(x,h,p,tol);
  }

 
} //end namespace 

#endif
