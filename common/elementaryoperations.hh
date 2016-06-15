#ifndef ELEMENTARY_OPERATIONS_HH
#define ELEMENTARY_OPERATIONS_HH

#include "../accuracy/floatingpointarithmetic.hh"
#include "../expressiontemplates/xsettings.hh"

namespace Adonis{
  
  /**
   *  \brief Addition, one element is altered using += 
   */
  template<class T>
  class AddBasicElements{
  public:
    typedef T value_type;

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif

    static inline void apply(T& a, const T& b){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      AdditionType::pluseq(a,b);
#else

      a += b;
#endif         
}
  
    //commutative
    static inline void apply(const T& b, T& a){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      AdditionType::pluseq(a,b);
#else  
      a += b;
#endif    
}
  
    //for const refs
    static inline T apply(const T& a, const T& b){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      return AdditionType::add(a,b);
#else
      return a+b;
#endif    
    }
  };


  /**
   * \brief Subtraction
   */
  template<class T>
  class SubtractBasicElements{
  public:
    typedef T value_type;

    static inline void apply(T& a, const T& b){
      a -= b;
    }
  
    //commutative
    static inline void apply(const T& b, T& a){
      a -= b;
    }
  
    //for const refs
    static inline T apply(const T& a, const T& b){
      return a-b;
    }
  };



  /**
   *  \brief Multiplication, one element is altered using += 
   */
  template<class T>
  class MultiplyBasicElements{
  public:
    typedef T value_type;
    static inline void apply(T& a, const T& b){
      a *= b;
    }
  
    //commutative 
    static inline void apply(const T& b, T& a){
      a *= b;
    }
  
    //for const refs
    static inline T apply(const T& a, const T& b){
      return a*b;
    }
    
  };
  

  /**
   *  \brief Division -- no efforts are made to check whether a is (close to) zero! 
   */ 
  template<class T>
  class DevideBasicElements{
  public:
    typedef T value_type;
    static inline void apply(T& a, const T& b){
      a /= b;
    }
  
    //commutative 
    static inline void apply(const T& b, T& a){
      a /= b;
    }
  
    //for const refs
    static inline T apply(const T& a, const T& b){
      return a/b;
    }
    
  };

}

#endif 
