#ifndef SOME_USEFUL_STUFF_FOR_ODE_HH
#define SOME_USEFUL_STUFF_FOR_ODE_HH

namespace Adonis{

  template<class T>
  class TemperatureDistribution{
  public:
    typedef T value_type;
    
    TemperatureDistribution(const T& hx = T(), const T& hy = T(), const T& x1 = T(), const T& x2 = T(), const T& t1 = T(), const T& t2 = T(), const T& t3 = T()):hx_(hx),hy_(hy),x1_(x1),x2_(x2),t1_(t1),t2_(t2),t3_(t3),t_(T()){}

    T& operator()(const T& x, const T& y){
      if(x < x1_)
	t_ = t1_;
      else if (x1_ <= x && x <= x2_)
	t_ = t2_;
      else if(x > x2_)
	t_ =  t3_;
      return t_;
    }
    
  private:
    const T& hx_, hy_, x1_, x2_,
      t1_,t2_,t3_;
    T t_;
  };
  
} //end namespace

#endif
