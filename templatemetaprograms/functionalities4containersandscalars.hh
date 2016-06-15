#ifndef FUNCTIONALITIES_4_BOTH_CONTAINERS_AND_SCALARS_HH
#define FUNCTIONALITIES_4_BOTH_CONTAINERS_AND_SCALARS_HH

#include "../common/adonisassert.hh"
#include "../common/isclass.hh"

namespace Adonis{

  /**
   * \brief in some cases it is useful to apply a function of the same name 
   * to a scalar as well as a vector, e.g. we want to apply the abs function to
   * both vectors and scalars without maintaining two differently named 
   * functions for vectors and scalars.
   *
   * This TMP ensures that the right type is considered
   */
  template<class X, bool V> class Functionalities4ContainersAndScalars;

  template<class X>
  class Functionalities4ContainersAndScalars<X,true>{ //!X is a container
  public:
    typedef typename X::value_type value_type;

    static inline X& perturbation(X& x, const typename X::value_type& eps = 1.e-13){
      adonis_assert(eps > 0 && eps < 1.e-9); //don't let become perturbation to large
      for(typename X::iterator it = x.begin(); it != x.end(); ++it){
	*it += eps;
      }
      return x;
    }
  };

  
  template<class X>
  class Functionalities4ContainersAndScalars<X,false>{ //!X is a scalar
  public:
    typedef X value_type;

    static inline X& perturbation(X& x, const X& eps = 1.e-13){
      adonis_assert(eps > 0 && eps < 1.e-9); //don't let become perturbation to large
      x += eps;
      return x;
    }
  };

  
  //!convenient functions
  template<class X, class T>
  X& perturbation(X& x, const T& eps){
    return Functionalities4ContainersAndScalars<X,IsContainer<X>::Value>::perturbation(x,eps);
  }
  
} //end namespace

#endif
