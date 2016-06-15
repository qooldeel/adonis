#ifndef CONTINUOUS_STIRRED_TANK_REACTOR_HH
#define CONTINUOUS_STIRRED_TANK_REACTOR_HH

#include "../../common/error.hh"
#include "../../expressiontemplates/exprvec.hh"
#include "../../containers/staticarray.hh"
#include "../../misc/misctmps.hh"

namespace Adonis{

  /**
   * Taken from M. Diehl, PhD thesis, 2001, Heidelberg
   */
  template<class T>
  class CSTR{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;

    CSTR(int n = 0):rhs_(4), q_(2),
		    k10_(1.287e+12), 
		    k20_(1.287e+12), 
		    k30_(9.043e+9), 
		    E1_(-9758.3), 
		    E2_(-9758.3), 
		    E3_(-8560), 
		    H1_(4.2), 
		    H2_(-11.0), 
		    H3_(-41.85), 
		    rho_(0.9342), 
		    Cp_(3.01), 
		    kw_(4032), 
		    AR_(0.215),
		    VR_(10),
		    mK_(5), 
		    CpK_(2.0), 
		    cA0_(5.1), 
		    theta0_(104.9){
      q_[0] = q_[1] = 1.; //default
    }


    const VType& controls() const {return q_;}

    std::size_t dim() const {return 4;}
    std::size_t domain_dim() const {return 4;}

    void set_controls(const T& q1, const T& q2){
      q_[0] = q1;
      q_[1] = q2;
    }

    const T& q1() const {return q_[0];}
    const T& q2() const {return q_[1];}

    template<class X>
    VType& operator()(const X& x){
      
      T k1 = rate(k10_,E1_,x[2]),
	k2 = rate(k20_,E2_,x[2]),
	k3 = rate(k30_,E3_,x[2]);

      rhs_[0] = q_[0]*(cA0_-x[0]) - k1*x[0] - k3*ntimes<2>(x[0]);
      rhs_[1] = -q_[0]*x[1] + k1*x[0] - k2*x[1];
      rhs_[2] = q_[0]*(theta0_-x[2])+(kw_*AR_)/(rho_*Cp_*VR_)*(x[3]-x[2])
	-(1./(rho_*Cp_))*(k1*x[0]*H1_ + k2*x[1]*H2_ + k3*ntimes<2>(x[0])*H3_);
      rhs_[3] = (1./(mK_*CpK_))*(q_[1] + kw_*AR_*(x[2]-x[3]));
    
      return rhs_;
    }

  private:
    VType rhs_, 
      q_;        //controls
    T k10_, k20_, k30_, E1_, E2_, E3_, H1_, H2_, H3_, rho_, Cp_, kw_, AR_,VR_,mK_, CpK_, cA0_, theta0_;

    

    T rate(const T& a, const T& energ, const T& temp){
      return a*exp(energ/(temp+273.15));
    }
  };

}

#endif
