#ifndef ODE_SOLVERS_THE_OTHER_WAY_HH
#define ODE_SOLVERS_THE_OTHER_WAY_HH

namespace Adonis{
   
  template<class MODEL>
  class ClassicalRK{
  public:
    typedef typename MODEL::value_type value_type;
    typedef typename MODEL::VType VType;

    ClassicalRK(MODEL& ode):ode_(ode),next_(ode_.dim()),k1_(ode_.dim()),k2_(ode_.dim()),k3_(ode_.dim()),k4_(ode_.dim()){}

    template<class TSTEP, class V>
    VType& step(const TSTEP& deltaT, const V& x){
      k1_ = deltaT*ode_(x);
      k2_ = deltaT*ode_(x+0.5*k1_);
      k3_ = deltaT*ode_(x+0.5*k2_);
      k4_ = deltaT*ode_(x + k3_);

      next_ = x + 1./6.*(k1_ + 2.*k2_ + 2.*k3_ + k4_);
      return next_;
    }
    

  private:
    MODEL& ode_;  //reference
    VType next_, k1_,k2_, k3_,k4_;
  };


}

#endif
