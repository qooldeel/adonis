#ifndef FUNCTIONS_TO_TEST_CONVERGENCE_OPTIMIZATION_ALGORITHMS_HH
#define FUNCTIONS_TO_TEST_CONVERGENCE_OPTIMIZATION_ALGORITHMS_HH

#include <cmath>

#include "../common/adonisassert.hh"
#include "../misc/misctmps.hh"


namespace Adonis{

  /**
   * \brief The classical Rosenbrock which has a global minimum at \f$ x^* = [1,1]^T.\f$
   */
  template<class X>
  inline typename X::value_type rosenbrock_function(const X& x){
    return  ( 100*ntimes<2>(x[1] - ntimes<2>(x[0])) + ntimes<2>(1 - x[0]) );
  }


  /**
   * \brief Himmelblau's function it has several local minima, four of them with
   * function value 0.
   */
  template<class X>
  inline typename X::value_type himmelblau_function(const X& x){
    return ( ntimes<2>(ntimes<2>(x[0]) + x[1] -11) + ntimes<2>(x[0] + ntimes<2>(x[1]) -7) );
  }


  /**
   * \brief Griewank function of order n for \f$ x_i \in [-600.,600]\f$
   * It has a <B> global minumum </B> at the origin \f$ x = 0 \f$ but
   * has many local minima
  */
  template<class X>
  inline typename X::value_type griewank_function(size_t n, const X& x){
    typedef typename X::value_type value_type;
    value_type sum = value_type(0),
      prod = value_type(1);

    for(size_t i = 0; i < n; ++i){
      sum += x[i]*x[i];
      prod *= cos(x[i]/sqrt(static_cast<value_type>(i)));
    }
    return ( 1 + sum/4000. - prod);
  }
  

} //end namespace 

#endif
