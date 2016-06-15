#ifndef SOME_SYSTEMS_OF_ODE_FOR_TESTING_REASONS_HH
#define SOME_SYSTEMS_OF_ODE_FOR_TESTING_REASONS_HH

#include<cassert>
#include "../../expressiontemplates/exprvec.hh"
#include "../../common/globalfunctions.hh"

namespace Adonis{

  /**
   * \brief Van der Pol equations as an example of a simple stiff system
   * IC: [2,0], time horizon: [0,3000], tolerances: 1.e-03 and 1.e-06, resp.
   *
   * Source:
   * [1] <a href="http://de.mathworks.com/help/matlab/ref/ode15s.html?searchHighlight=ode15s"> Example of a stiff system in two variables </a> 
   *
   * 
   */
  template<typename T>
  class VanderPolOscillator{
  public:
  
    typedef Adonis::ExprTmpl::MyVec<T> VType;
    typedef T value_type;

    VanderPolOscillator(size_t dim = 0):rhs_(dim){}

    enum{
      domainDim =2
    };

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const{return domainDim;}

    template<class X>
    VType& operator()(const X& y){
    
      //========================================================================
  
      rhs_[0] = y[1];
      rhs_[1] = 1000.*(1-ntimes<2>(y[0]))*y[1] - y[0];

      //========================================================================

      return rhs_;
    }

    //! fake time dependent operator
    template<class V>
    VType& operator()(const T& time, const V& y){
      return (rhs_ = (*this).operator()(y));
    }

  private:
    VType rhs_;
  };

  
  /**
   * \brief A simple nonlinear stiff system with known exact solution at time t
   * \f$ y(t) = (exp(-2t),exp(-t))\f$ evaluted for \f$ y(0) = (1,1).\f$ 
   * The solution is sought on [0,1], cf. [1]
   *
   * Reference:
   *
   *  [1] <a href="http://www.m-hikari.com/imf-password2008/13-16-2008/tahmasbyIMF13-16-2008.pdf"> Nice nonlinear stiff system </a> 
   *
   */
  template<typename T>
  class SimpleStiffSys{
  public:
  
    typedef Adonis::ExprTmpl::MyVec<T> VType;
    typedef T value_type;

    SimpleStiffSys(size_t dim = 0):rhs_(dim){}

    enum{
      domainDim =2
    };

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const{return domainDim;}

    template<class X>
    VType& operator()(const X& y){
    
      //========================================================================
  
      rhs_[0] = -1002*y[0] + 1000*y[1]*y[1]; //no std::sin coz of CppAD!
      rhs_[1] = y[0] - y[1]*(1+y[1]);

      //========================================================================

      return rhs_;
    }

    //! fake time dependent operator
    template<class V>
    VType& operator()(const T& time, const V& y){
      return (rhs_ = (*this).operator()(y));
    }

  private:
    VType rhs_;
  };

  /**
   * \brief Note that operator version with expression are often to strict (especially when using AD via CppAD. <B> Better: Use Abstract template argument </B> instead of special ones and let your compiler do the job ;)
   *
   * NOTE 1: This may not possible when using, e.g. CppAD. For that you can define the reactions based on some STL-compliant vector (see below) 
   *
   * NOTE 2: Mind that size_t and unsigned are different on x86_64 machines!!!
   *         Therefore it is recommended to use either int or size_t 
   *
   */
  //example taken from [KINCAID/CHENEY, "Numerical Analysis: Mathematics of Scientific Computing, p. 524"]
  template<typename T>
  class KincaidCheney{
  public:
    
    typedef  Adonis::ExprTmpl::MyVec<T> VType;
    typedef T value_type;

    KincaidCheney(size_t dim=0):rhs_(dim){} //!invoke constructor of base class
    
    size_t domain_dim() const {return 2;}
    size_t dim() const {return rhs_.size();}

  
    template<class E>
    VType& operator()(const T& t, const E& y){
   
      //============== Define test RHS here, hunk =========================
      
      
      rhs_[0] = y[0] + 4*y[1] - exp(t);
      rhs_[1] = y[0] + y[1] + 2.*exp(t);
      
      //===================================================================
   
      return rhs_;
    }

  
    Adonis::ExprTmpl::MyVec<T>& exact_solution(const T& t, const T& a, const T& b){
    
      rhs_[0] = 2.*a*exp(3.*t) - 2.*b*exp(-t) - 2.*exp(t);
      rhs_[1] = a*exp(3.*t) + b*exp(-t) + 0.25*exp(t);
      
      return rhs_;
    }
  
  private:
    VType rhs_;
  };



  //another equation system of 2 equations as stated in [LAMBERT, "Numerical methods for ODEs"]. Integrate over [0,10] with iv's [-12,6]
  //Dies the computed solution compare favorably with the exact soltution?
  template<typename T>
  class Lambert{
  public:
    
    typedef Adonis::ExprTmpl::MyVec<T> VType;
    typedef T value_type;

    Lambert(size_t dim=0):rhs_(dim){}

    size_t domain_dim() const {return 2;}
    size_t dim() const {return rhs_.size();}
    
    template<class E>
    VType& operator()(const T& t, const E& y){
   
      //============== Define test RHS here, hunk =========================
   
      T c = std::cos(6.*t), s = std::sin(6.*t);
   
      rhs_[0] = (-1.-9.*c*c + 12.*s*c)*y[0] + (12.*c*c + 9.*s*c)*y[1];
      rhs_[1] = (-12.*s*s +9.*s)*y[0] + (-1. - 9.*s*s - 12*s*c)*y[1];
    
      //===================================================================
    
      return rhs_;
    }
    
    
    //the exact solution looks somewhat freaky (and is more or less 0) ;-)
    Adonis::ExprTmpl::MyVec<T>& exact_solution(const T& t){
    
      T c = std::cos(6.*t), s = std::sin(6.*t);

      rhs_[0] = std::exp(-13.*t)*(s-2.*c);
      rhs_[1] = std::exp(-13.*t)*(2.*s+c);

      return rhs_;
    }
 
  private:
    VType rhs_;

  };


  //a system comprising 3 equations drawn from [BURDEN/FAIRES, "Numerical Analysis"]
  template<typename T>
  class BurdenFaires{
  public:
  
    typedef Adonis::ExprTmpl::MyVec<T> VType;
    typedef T value_type;

    BurdenFaires(size_t dim=0):rhs_(dim){}

    size_t domain_dim() const {return 2;}
    size_t dim() const {return rhs_.size();}
  
    template<class E>
    VType& operator()(const T& t, const E& y){
   
      //============== Define test RHS here, hunk =========================
  
   
      rhs_[0] = y[1];
      rhs_[1] = -y[0] - 2.*std::exp(t) +1;
      rhs_[2] = -y[0] - std::exp(t) +1;
     
      //===================================================================
      
      return rhs_;
    }

  private:
    VType rhs_;
  };


//the famous LOTKA-VOLTERRA predator-prey model. It has a PERIODIC solution which does not have an explicit analytical formula. 
template<typename T>
class LotkaVolterra{
public:
  
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  typedef T value_type;
  //always of dim 2, hence set as default
  LotkaVolterra(size_t dim = 0):rhs_(dim){}
  
  template<class E>
  VType& operator()(const T& t, const E& y){
   
    //============== Define Lotka-Volterra equations ======================
  
    //the for growth parameters
    T alpha = 0.03, beta = 0.2, 
      gamma = 0.4, delta = 0.01; 

    rhs_[0] = (alpha - beta*y[1])*y[0];
    rhs_[1] = (-gamma + delta*y[0])*y[1];
 
     
    //====================================================================
   
    return rhs_;
  }

private:
  VType rhs_;
};









//Chemical Oscillation -- the Oregonator reactions (simplyfied version of the famous Zhabotinski-Belousov reaction [cf. DEUFLHARD/BORNEMANN, "Scientific Computing with ODEs"]. Despite its apparent simplicity, this coupled system cannot be solved analytically any more!    
template<typename T>
class Oregonator{
public:
  
  //if you need the dimension of the domain space as well
  enum{
    domainDim = 5  
  };
  
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  typedef T value_type;

  Oregonator(size_t dim=0):rhs_(dim), k1_(1.34),k2_(1.6e9),k3_(8.0e3),k4_(4.0e7),k5_(1.){}

  size_t dim() const {return rhs_.size();}
  size_t domain_dim() const{return domainDim;}

 
  template<class V>
  VType& operator()(const V& c){
   
    //============== Define test RHS here, hunk =========================
                      //STIFF:
   
    rhs_[0] = -k1_*c[0]*c[1] - k3_*c[0]*c[2];
    rhs_[1] = -k1_*c[0]*c[1] - k2_*c[1]*c[2] + k5_*c[4];
    rhs_[2] = k1_*c[0]*c[1] - k2_*c[1]*c[2] + k3_*c[0]*c[2] - 2.*k4_*Adonis::sqr(c[2]);
     
    rhs_[3] = k2_*c[1]*c[2] + k4_*Adonis::sqr(c[2]);
    rhs_[4] = k3_*c[0]*c[2] - k5_*c[4];
      
    //===================================================================
   
    return rhs_;
  }

  template<class V>
  Adonis::ExprTmpl::MyVec<T>& operator()(const T& time, const V& c){
    rhs_ = (*this).operator()(c);
    return rhs_;
  }

private:
  VType rhs_;
  T k1_, k2_, k3_, k4_, k5_;
};

//non-stiff but OSCILLATORY test case 
template<typename T> 
class Brusselator{
public:
   enum{
    domainDim = 2  
  };


  typedef Adonis::ExprTmpl::MyVec<T> VType;
  typedef T value_type;

  Brusselator(size_t dim=0):rhs_(dim){}

  size_t dim() const {return rhs_.size();}
  size_t domain_dim() const{return domainDim;}

  /**
   *\brief The Brusselator (without mass transport). See [Hairer,Wanner, Vol II, §IV, eq. (1.6)]
   * 
   * In the non-diffusion case, the ICs read as \f$ y = [1,3]^T\f$ (altern.: \f$ y = [1,1]^T\f$)
   */
  template<class W>
  VType& operator()(const W& y){
    const T A = 1,
      B = 3;
    
    rhs_[0] = A + Adonis::sqr(y[0])*y[1] -(B+1)*y[0];
    rhs_[1] = B*y[0] - Adonis::sqr(y[0])*y[1]; 

    return rhs_;
  }

  //'time-dependent' version to use with time dependent solvers
  template<class TIME,class W>
  VType& operator()(const TIME& time, const W& y){
    return (rhs_ = (*this).operator()(y));
  }

private:
  VType rhs_;
};


/**
 * \brief A stiff (hopefully) non-oscillatory chemical system. See [HAIRER,WANNER,vol II,p.3] 
 * IC: \f$ y_0 = [1,0,0]^T \f$

 * Reference solution at the end of integration interval:
 * [2.083344015e-08,8.33336077e-14,9.99999999e-01]^T
*/
template<typename T>
class HairerRobertson{
public:
  enum{
    domainDim = 3  
  };
  
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  typedef T value_type;

  HairerRobertson(size_t dim=0):rhs_(dim){}

  size_t dim() const {return rhs_.size();}
  size_t domain_dim() const{return domainDim;}

  template<class X>
  VType& operator()(const X& y){
    //===================== Model after Robertson ============================

    rhs_[0] = -0.04*y[0] + 1e4*y[1]*y[2];
    rhs_[1] = 0.04*y[0] - 1e4*y[1]*y[2] - 3e7*Adonis::sqr(y[1]);
    rhs_[2] = 3.e7*Adonis::sqr(y[1]);
  
    //========================================================================
    return rhs_;
  }


  //! fake time dependent operator
  template<class V>
  VType& operator()(const T& time, const V& y){
    return (rhs_ = (*this).operator()(y));
  }
  
private:
  VType rhs_;
};


/**
 * \brief A stiff and nonlinear example taken from [RALSTON/RABINOWITZ, p.231]
 * \f$ y_0 = [0,0]^T, \quad t \in [0,100]\f$
 *
 * The numerical solution obtained with Gear's implicit 3rd order method are
 *  y(t=100) = [-0.99164187, 0.98333613]^T
 */
template<typename T>
class RalstonRabinowitz{
public:
 
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  typedef T value_type;

  enum{
    domainDim =2
  };

  RalstonRabinowitz(size_t dim=0):rhs_(dim){}
  
  size_t dim() const {return rhs_.size();}
  size_t domain_dim() const{return domainDim;}

  template<class V>
  VType& operator()(const V& y){
    
    //========================================================================
    
    rhs_[0] = .01 - (.01 + y[0] + y[1])*(1+(y[0] + 1000)*(y[0] + 1));
    rhs_[1] = .01 - (.01 + y[0] + y[1])*(1 + Adonis::sqr(y[1]));

    //========================================================================

    return rhs_;
  }
  
  //! fake time dependent operator
  template<class V>
  VType& operator()(const T& time, const V& y){
    return (rhs_ = (*this).operator()(y));
  }

private:
  VType rhs_;
};


/**
 *\brief Stiff non-autonomous (scalar) example \f$ y'(t) = -100(y(t) - sin(t)), t \geq 0, y(0) = 1 \f$ from [Ascher/Petzold, §3.4, p. 49]
 * Note that every non-autonomous \f$n\f$-dim system \f$ y' = f(t,y)\f$ can be transformed to an autonomous one via introducing an additional variable \f$ w'(t) = 1, w(t_0) = t_0 , w(t) = t. \f$ This leads to the new \f$(n+1)\f$-dim sytem \f[ [y'(t),w'(t)]^T = [f(y(t),w(t)),1]^T\f] 
 */
template<typename T>
class Non2Autonomous{
public:
  
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  typedef T value_type;

  Non2Autonomous(size_t dim = 0):rhs_(dim){}

  enum{
    domainDim =2
  };

  size_t dim() const {return rhs_.size();}
  size_t domain_dim() const{return domainDim;}

  template<class X>
  VType& operator()(const X& y){
    
    //========================================================================
  
    rhs_[0] = -100*(y[0] - sin(y[1]) ); //no std::sin coz of CppAD!
    rhs_[1] = 1;

    //========================================================================

    return rhs_;
  }


  //! fake time dependent operator
  template<class V>
  VType& operator()(const T& time, const V& y){
    return (rhs_ = (*this).operator()(y));
  }

private:
  VType rhs_;
};


/**
 * \brief Stiff non-autonomous system taken from [Burden/Faires,§5.11, p. 335] 
 *  \f$ \f$
 */
template<typename T>
class StiffBurdenFaires{
public:
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  typedef T value_type;

  StiffBurdenFaires(size_t dim = 0):rhs_(dim){}

  enum{
    domainDim = 3     //2 species + 1 (time) 
  };

  size_t dim() const {return rhs_.size();}
  size_t domain_dim() const{return domainDim;}

  template<class X>
  VType& operator()(const X& u){
    //============ transformed non-autonomous system to autonomous one =========
    rhs_[0] = 9*u[0] + 24*u[1] + 5*cos(u[2]) - 1./3.*sin(u[2]);
    rhs_[1] = -24*u[0] - 51*u[1] - 9*cos(u[2]) + 1./3.*sin(u[2]);
    rhs_[2] = 1.;       //u[2] ^= time 
    
    //==========================================================================
    
    return rhs_;
  }

  //! fake time dependent operator
  template<class V>
  VType& operator()(const T& time, const V& y){
    return (rhs_ = (*this).operator()(y));
  }

private:
  VType rhs_;
};


} //end namespace

#endif
