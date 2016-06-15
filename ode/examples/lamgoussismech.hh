#ifndef LAM_GOUSSIS_MECHANISM_HH
#define LAM_GOUSSIS_MECHANISM_HH

#include "../common/adonisassert.hh"
#include "../../misc/misctmps.hh"
#include "../../io/readinparameters.hh"
namespace Adonis{

  /**
   * \brief Lam-Goussis combustion mechanism
   *
   *  ODEs are given in [1], the full stuff can be found in [2]
   *  
   *  References:
   *  
   * [1]  Gear, Kaper, Kevrekidis and Zagaris, "Projecting to a slow manifold: Singularly perturbed systems and legacy codes", SIADS, 4, 2005, pp. 711--732 
   * [2]  Lam and Goussis, "Understanding complex chemical kinetics with CSP", 22nd Symp. on Combustion/The Combustion Institute, 1988, pp. 931--941
   *
   */
  template<class T>
  class LamGoussisMech{
  public:
    typedef ExprTmpl::MyVec<T> VType;
    typedef T value_type;
    typedef typename TypeAdapter<T>::BaseType BaseType;

    LamGoussisMech(unsigned n = 0):rhs_(n), mu_(4.5e-4){
   
      k1f_ = 1.0136e12;    k1b_ = 1.1007e13;
      k2f_ = 3.5699e12;    k2b_ = 3.2105e12;
      k3f_ = 4.7430e12;    k3b_ = 1.8240e11;
      k4f_ = 6.0000e13;   
      k5f_ = 6.2868e15;
      k8f_ = 6.5325e12;    k8b_ = 3.1906e11;
      
     
    }

    size_t dim() const {return 7;}
    size_t domain_dim() const{return 7;}

  
    template<class X>
    VType& operator()(const X& y){
      rhs_[0] = -k1f_*y[0]*y[1]  + k1b_*y[2]*y[3] + k4f_*y[2]*y[6] - mu_*k5f_*y[0]*y[1];
     
      rhs_[1] = -k1f_*y[0]*y[1] + k1b_*y[2]*y[3] + k2f_*y[3]*y[4] - k2b_*y[1]*y[2] + k3f_*y[2]*y[4] - k3b_*y[1]*y[5] - mu_*k5f_*y[0]*y[1];

      rhs_[2] = k1f_*y[0]*y[2] - k1b_*y[2]*y[3] + k2f_*y[3]*y[4] + k2b_*y[1]*y[2] 
	-k3f_*y[2]*y[4] + k3b_*y[1]*y[5] - k4f_*y[2]*y[6] -2*k8f_*y[2]*y[2] + 2*k8b_*y[3]*y[5];

      rhs_[3] = k1f_*y[0]*y[2] - k1b_*y[2]*y[3] - k2f_*y[3]*y[4] + k2b_*y[1]*y[2] + k8f_*y[2]*y[2] - k8b_*y[3]*y[5];

      rhs_[4] = -k2f_*y[3]*y[4] + k2b_*y[1]*y[2] - k3f_*y[2]*y[4] + k3b_*y[1]*y[5];
      
      rhs_[5] = k3f_*y[2]*y[4] - k3b_*y[1]*y[5] + k4f_*y[2]*y[6] + k8f_*y[2]*y[2] - k8b_*y[3]*y[5];

      rhs_[6] = -k4f_*y[2]*y[6] + mu_*k5f_*y[0]*y[2];

      return rhs_;
    }

   

  private:
    VType rhs_;
    BaseType mu_,
      k1f_, k2f_, k3f_, k4f_, k5f_, k8f_,
      k1b_, k2b_, k3b_, k8b_;
  };

}

#endif
