#ifndef BASICS_OF_PROBABILITY_HH
#define BASICS_OF_PROBABILITY_HH

#include "../accuracy/floatingpointarithmetic.hh"
#include "../expressiontemplates/exprvec.hh"

namespace Adonis{
  
  /**
   * \brief Calculate sample mean. Note that the sample mean of a scalar is the
   * scalar itself.
   */
  template<class X, bool B> class SampleMean;

  template<class X>
  class SampleMean<X,true>{  //o.k. container detected  
  public:
    typedef typename X::value_type value_type;
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType; //! floating point addition type
#endif

    static inline value_type value(const X& x){
      size_t n = static_cast<size_t>(std::distance(x.begin(),x.end()));
      adonis_assert(n != 0);
      value_type mean = value_type();
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      value_type s = value_type(),
	c = value_type();
#endif
      for(size_t i = 0; i < n; ++i){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	mean = AdditionType::add(x[i],s,c);
#else
	mean += x[i];
#endif
      }

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<value_type,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
      if(Sgn(Abs(c)) != 0){  //! if c == 0 then do not assign s since 
	//! it might be 0 because no addition was involved!
	CorrectRoundingErrorAfterwards<value_type,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(mean,s);
      }
#endif 

      return 1./n*mean;
    }
  };


  template<class X>
  class SampleMean<X,false>{  //o.k. only scalar 
  public:
    typedef X value_type;
    static inline const X& value(const X& x){
      return x;
    }
  };


  //!convenient functino
  template<class T>
  inline typename SampleMean<T,IsContainer<T>::Value>::value_type mean(const T& x){
    return SampleMean<T,IsContainer<T>::Value>::value(x);
  }

}

#endif
