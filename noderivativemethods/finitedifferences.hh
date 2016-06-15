#ifndef FINITE_DIFFERENCE_APPROXIMATIONS_TO_FUNCTIONS_HH
#define FINITE_DIFFERENCE_APPROXIMATIONS_TO_FUNCTIONS_HH

#include <iostream>
#include <limits> 
#include <cmath>
#include <string>

#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../common/error.hh"
#include "../expressiontemplates/expressionvectortraits.hh"
#include "../expressiontemplates/exprvec.hh" 
#include "../expressiontemplates/scalar.hh"
#include "../misc/useful.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../common/isclass.hh"

namespace Adonis{

  /**
   * \brief compute finite difference. Compile-time determination whether vector or double is used.
   */
  template<bool B, class J> class ArgumentTypeDependency;
 
  /**
     * \brief central finite difference for approximating first derivatives of \f$ f: R^n \longrightarrow R^m\f$ , i.e. \f[ \frac{\partial f}{\partial x_i}(x) = \nabla f(x)^Te_i \approx \frac{f(x + \eps e_i) - f(x + \eps e_i)}{2\eps},\f]
     * where $\varepsilon = \sqrt[3]{\mathrm{MEPS}}.$ MEPS denotes the machine precision and $e_i$ the \f$i\f$-th canonical unit vector (cf. ref. [2]). Note that this choice of $\varepsilon$ is optimal w.r.t. central fd (anyway you might try \f$ \sqrt[2]{\mathrm{MEPS}}\f$ as in the forward fd case as well.
     *
     * NOTE 1: To compute a dense Jacobian with central differences, \f$(2n+1) \f$  function evaluations of \f$ f.\f$ 
     *
     *NOTE 2: a construction, such as <TT> col = (fun_(xP + epsithP) - fun_(xP - epsithP))/2*eps </TT> is possible when <TT>fun_</TT>'s operator returns a vector instead of a reference. 
     *
     * REFERENCE:
         [1] [FLETCHER, Practical Methods of Optimization, 2nd ed., p. 20]
	 [2] [NOCEDAL,WRIGHT], Numerical Optimization, 2nd ed., § 8.1, p. 198]
 */
  template<class J>
  class ArgumentTypeDependency<true,J>{  //container
  public:
    typedef J ExprVecType;
    typedef typename J::value_type value_type;
     
    template<class X, class FUN>
    static inline void forward_fd_jacobian(J& jac, const X& x,unsigned n, unsigned m, const value_type& eps, FUN& fun){
      adonis_assert(x.size() == n);  //only plausible for containers
      
      (jac.size() == 0) ? jac.resize(n*m) : do_nothing();

      ExprVecType col(m),  //ith column of Jacobian
	ithunitvector(n),  //size of domain vector
	epsith(n);
      
      for(unsigned i = 0; i < n; ++i){
	epsith = eps*canonical_base_vector(ithunitvector,i);

	col = fun(x + epsith);
	col -= fun(x);
	col /= eps;
	
	//std::cout << "column "<< i << " = "<< col << std::endl;
	
	build_vector(jac,col,n,i);  //create Jacobian -- stored row-wisely
      }
    }

    //when fx is known
    template<class X, class FUN, class FX>
    static inline void forward_fd_jacobian(J& jac, const X& x,unsigned n, unsigned m, const value_type& eps, FUN& fun, const FX& fx){
      adonis_assert(x.size() == n);  //only plausible for containers
      
      (jac.size() == 0) ? jac.resize(n*m) : do_nothing();

      ExprVecType col(m),  //ith column of Jacobian
	ithunitvector(n),  //size of domain vector
	epsith(n);
      
      for(unsigned i = 0; i < n; ++i){
	epsith = eps*canonical_base_vector(ithunitvector,i);

	col = fun(x + epsith);
	col -= fx;
	col /= eps;
	
	//std::cout << "column "<< i << " = "<< col << std::endl;
	
	build_vector(jac,col,n,i);  //create Jacobian -- stored row-wisely
      }
    }

    template<class X, class FUN>
    static inline void central_fd_jacobian(J& jac, const X& x,unsigned n, unsigned m, const value_type& eps, FUN& fun){
      adonis_assert(x.size() == n);  //only plausible for containers
      
      (jac.size() == 0) ? jac.resize(n*m) : do_nothing();

      ExprVecType col(m),  //ith column of Jacobian
	ithunitvector(n),  //size of domain vector
	epsith(n);
      
      for(unsigned i = 0; i < n; ++i){
	epsith = eps*canonical_base_vector(ithunitvector,i);

	col = fun(x + epsith);
	col -= fun(x - epsith);
	col /= (2.*eps);
	
	//std::cout << "column "<< i << " = "<< col << std::endl;
	
	build_vector(jac,col,n,i);  //create Jacobian -- stored row-wisely
      }
    }
  };


template<class J>
class ArgumentTypeDependency<false,J>{  //no container
  public:
   
  typedef J ExprVecType;
  typedef typename J::value_type value_type;
    
  template<class X, class FUN>
  static inline void forward_fd_jacobian(J& jac, const X& x,unsigned n, unsigned m, const value_type& eps, FUN& fun){
    adonis_assert(n*m == 1);
    jac.resize(1);
    jac[0] = (fun(x + eps) - fun(x))/eps;
  }

  //when fx = fun(x) is known
  template<class X, class FUN, class FX>
    static inline void forward_fd_jacobian(J& jac, const X& x,unsigned n, unsigned m, const value_type& eps, FUN& fun, const FX& fx){
    adonis_assert(n*m == 1);
    jac.resize(1);
    jac[0] = (fun(x + eps) - fx)/eps;
  }

  template<class X,class FUN>
    static inline void central_fd_jacobian(J& jac, const X& x,unsigned n, unsigned m, const value_type& eps, FUN& fun){
      adonis_assert(n*m == 1);
      jac.resize(1);

      jac[0] = (fun(x + eps) - fun(x - eps))/(2*eps);
    }
};
  

  template<int N> class WhatFiniteDifference;

  template<> 
  class WhatFiniteDifference<1>{  //! forward
  public:
     template<class RES, class X, class D, class FUN>
    static inline void fd(RES& jac, const X& x, unsigned n, unsigned m, const D& eps, FUN& fun){ 
      ArgumentTypeDependency<IsContainer<X>::Value,RES>::forward_fd_jacobian(jac,x,n,m,eps,fun);
     }

    template<class RES, class X, class FX, class D, class FUN>
    static inline void fd(RES& jac, const X& x, unsigned n, unsigned m, const D& eps, FUN& fun, const FX& fx){ 
      ArgumentTypeDependency<IsContainer<X>::Value,RES>::forward_fd_jacobian(jac,x,n,m,eps,fun,fx);
    }
  };

    
  template<> 
  class WhatFiniteDifference<2>{  //! central
  public:
    template<class RES, class X, class D, class FUN>
    static inline void fd(RES& jac, const X& x, unsigned n, unsigned m, const D& eps, FUN& fun){  
      ArgumentTypeDependency<IsContainer<X>::Value,RES>::central_fd_jacobian(jac,x,n,m,eps,fun);
    }

    template<class RES, class X, class FX, class D, class FUN>
    static inline void fd(RES& jac, const X& x, unsigned n, unsigned m, const D& eps, FUN& fun, const FX& fx){  //! fx is only a dummy here!!
      ArgumentTypeDependency<IsContainer<X>::Value,RES>::central_fd_jacobian(jac,x,n,m,eps,fun);
    }
  };


  /**
   * Finite Difference evaluation of a function 
   *
   * \tparam FUN (any) functor possessing a ()-operator 
   * \tparam N chooses type of finite difference: 1 = forward, 2 = central
   */
  template<class FUN, int N = 2> 
  class FiniteDifference{
  public:
    typedef FUN FunType;
    typedef typename FUN::value_type value_type;
    typedef typename ExprVecTraits<value_type>::ExprVecType ExprVecType;
    typedef ExprVecType ReturnType; //type of calculated Jacobian

    FiniteDifference(FunType& fun):fun_(fun){}

    enum{whatfd = N};

    void info(){
      std::string str;
      if(N == 1)
	str = "FORWARD";
      if(N == 2)
	str = "CENTRAL";

      std::cout << str << " FD used."<< std::endl;
    }
    
    template<class X>
    ExprVecType& fwd(const X& x, const value_type& eps = std::pow(std::numeric_limits<value_type>::epsilon(),1./3)){
      //std::cout << "Central finite difference" << std::endl;

      unsigned n = fun_.domain_dim(),
	m = fun_.dim();

      ArgumentTypeDependency<IsContainer<X>::Value,ExprVecType>::forward_fd_jacobian(jac_,x,n,m,eps,fun_);
    
      return jac_;
    }

    template<class X>
    ExprVecType& central(const X& x, const value_type& eps = std::pow(std::numeric_limits<value_type>::epsilon(),1./3)){
      //std::cout << "Central finite difference" << std::endl;

      unsigned n = fun_.domain_dim(),
	m = fun_.dim();

      ArgumentTypeDependency<IsContainer<X>::Value,ExprVecType>::central_fd_jacobian(jac_,x,n,m,eps,fun_);
    
      return jac_;
    }
     
    //! for convenience.
    template<class X>
    ExprVecType& jacobian(const X& x, const value_type& eps = std::pow(std::numeric_limits<value_type>::epsilon(),1./3)){
      return (*this).central(x,eps);
    }

    //! yet for another convenience -- Note the capital 'J' in the member's name
    template<class X>
    ExprVecType& Jacobian(const X& x, const value_type& eps = std::pow(std::numeric_limits<value_type>::epsilon(),1./3)){
      return (*this).central(x,eps);
    }
    
    //! tuned version to get Jacobian 
    template<class X, class FX>
    ExprVecType& Jacobian(const X& x, const FX& fx, const value_type& eps){
      unsigned n = fun_.domain_dim(),
	m = fun_.dim();

      WhatFiniteDifference<N>::fd(jac_,x,n,m,eps,fun_,fx);

      ArgumentTypeDependency<IsContainer<X>::Value,ExprVecType>::central_fd_jacobian(jac_,x,n,m,eps,fun_);
    
      return jac_;
    }

    template<class X, class FX>
    ExprVecType& jacobian(const X& x, const FX& fx, const value_type& eps){
      //std::cout << "templatized Jacobian"<< std::endl;
      return (*this).Jacobian(x,fx,eps);
    }
 
    //! dummy functions
    template<class X>
    ExprVecType& hessian(const X& x, int i=0){
      ADONIS_ERROR(DerivativeError,"No Hessian implemented so far due to accuracy, pal ;)");
      return jac_;
    }

      template<class X>
    ExprVecType& Hessian(const X& x, int i=0){
      ADONIS_ERROR(DerivativeError,"No Hessian implemented so far due to accuracy, pal ;)");
      return jac_;
    }

  private:
    FunType& fun_;  
    ExprVecType jac_;
  };

}

#endif 
