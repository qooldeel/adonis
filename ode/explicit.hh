#ifndef EXPLICIT_METHODS_EASY_IMPLEMENTATION_HH
#define EXPLICIT_METHODS_EASY_IMPLEMENTATION_HH


namespace Adonis{

  template<class ODE>
  class ExplicitTrapezoidal{
  public:
    typedef typename ODE::value_type value_type;
    typedef typename ODE::VType VType;
    typedef std::size_t SizeType;

    enum{OOM = 2};  //order of method

    ExplicitTrapezoidal(ODE& ode):model_(ode),predict_(ode.dim()),res_(ode.dim()){}

    template<class T, class X>
    VType& step(const T& kn, const X& x){
      predict_ = x + kn*model_(x); 
      res_ = x + 0.5*kn*(model_(x) + predict_);
      return res_;
    }

    template<class T, class X, class Y>
    VType& step(const T& kn, const X& eval, const Y& u){
      predict_ = u + kn*model_(eval); 
      res_ = u + 0.5*kn*(model_(eval) + predict_);
      return res_;
    }

  private:
    ODE& model_;
    VType predict_,res_;
  };

}

#endif
