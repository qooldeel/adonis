#ifndef TUNE_YOUR_ROUNDING_ERRORS_INDUCED_BY_FLOATING_POINT_OPERATIONS_HH
#define TUNE_YOUR_ROUNDING_ERRORS_INDUCED_BY_FLOATING_POINT_OPERATIONS_HH

#include "../common/globalfunctions.hh"
#include "../common/fancymessages.hh"
#include "../common/adonisassert.hh"
#include "../common/numerictypechecker.hh"

namespace Adonis{

/**
 * \brief Kahan summation algorithm 
 *
 *Problem: Summation of many floating-point numbers involves <I> rounding </I>. Balancing algorithms such as Kahan's algorithm try to correct the rounding errors, cf. [1, 2].
 *
 *  Given : \f$ x_1, \ldots, x_n \in R\f$. 
 *  Goal: calculate \f$ S = \sum_{i=1}^n x_i\f$ with smaller rounding errors
 *
 * References:
 *
 *   [1] [A. Klein, "A generalized Kahan-Babu\<ska-summation-algorithm", Computing, 76, 2006, pp. 279 -- 293]
 * 
 *   [2] [D. Goldberg, "What every computer scientist should know about floating-
point arithmetic", ACM Computing Surveys, 23, 1991, p. 45]
*
*   for a slightly different form, cf. 
*
*    [3] [N. Higham, "Accuracy and stability of numerical algorithms", 2nd ed., 2002, § 4.3, pp. 83, esp. Algo 4.2, p. 84]
*
* NOTE1: depending on the optimization options (don't use "-fassociative-math" 
*        when you want to use Kahan(-Babuska) summation) and when using <TT> 
*        valgrind </TT>, different rounding results may occur.
*
* NOTE2: it even seems that an unrolled loop (or equivalently, a TMP or even 
*        a unrolled loop by hand) yields different results.
*/
  template<char C, class A> class TunedFloatingPointAddition;

  //! if you intend to use a direct loop unrolling, use the following in 
  //! conjunction with, e.g. <TT> basicmetaprogramms/unrollloop.hh </TT> 
  template<class A>
  class TunedFloatingPointAddition<'k',A>{
  public:
    typedef typename A::value_type value_type;
    
    static inline value_type addition(size_t i, const A& a, value_type& s, value_type c){
      value_type y = a[i] - c,
	t = s + y;
      c = (t - s) - y;
      return (s = t); 
    }
  };

  template<class A>
  class TunedFloatingPointAddition<'K',A>{
  public:
    typedef typename A::value_type value_type;
    
    static inline value_type addition(size_t i, const A& a, value_type& s, value_type c){
      return TunedFloatingPointAddition<'k',A>::addition(i,a,s,c);
    }
  };
  

  //!Babuska-Kahan
  template<class A>
  class TunedFloatingPointAddition<'b',A>{
  public:
    typedef typename A::value_type value_type;
    
    static inline value_type addition(size_t i, const A& a, value_type& s, value_type c){
      value_type t = s + a[i];
      (Abs(s) >= Abs(a[i])) ?  c += ((s-t)+a[i]) : c += ((a[i]-t)+s);

      return (s = t);
    }
  };
  

   template<class A>
   class TunedFloatingPointAddition<'B',A>{
   public:
     typedef typename A::value_type value_type;
     
     static inline value_type addition(size_t i, const A& a, value_type& s, value_type c){
       return TunedFloatingPointAddition<'b',A>::addition(i,a,s,c);
     }
   };


  //!usual for-loop case
  template<class AT>
  inline typename AT::value_type kahan_algorithm(const AT& a){
    typedef typename AT::value_type value_type;
    
    //! only meaningful for floatingpoints
    NumericDataTypeChecker<value_type>::certify();

    unsigned n = std::distance(a.begin(),a.end());
    
    value_type s = value_type(),       //the sum
      c = value_type(),
      y, t;
    
    for(unsigned k = 0; k < n; ++k){
      y = a[k] - c;
      t = s + y;
      c = (t - s) - y;
      s = t;
    }
    
    return s;
  }


  /**
   * \brief Kahan summation for only 2 numbers
   */
  template<class T>
  inline T kahan_algorithm(const T& a, const T& b){
    NumericDataTypeChecker<T>::certify();

    T s = a,  //a assigned 
      c = T(),
      y, t;

    y = b - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
    
    return s;
  }


  //! <TT> val </TT> is the value to be added (corresponds to <TT>x[k] </TT> 
  //! in ref. [1]
  template<class T>
  inline T& kahan_algorithm(const T& val, T& s, T& c){
    NumericDataTypeChecker<T>::certify();
    
    T y = val - c,
      t = s + y;
    c = (t -s) - y;   //! update c and s
    s = t;
    return s;
  }



  /**
   * \brief For 2 numbers, circumvent a+=b by defining the following assignment where the Kahan summation is applied
   */
  template<class T>
  inline T& kahan_assign(T& a, const T& b){
    NumericDataTypeChecker<T>::certify();

    T c = T(),
      y, t;
       
    y = b - c;
    t = a + y;
    c = (t - a) - y;
    a = t;
  
    return a;
  }
  

  
  
  /**
   * \brief Kahan-Babuska algorithm, cf. [1]
   *  
   * NOTE: since it involves conditional statements such as <TT> if </TT>, it 
   *       might cause problems when being used in conjunction with Algorithmic 
   *       Differentiation (AD).
   */
  template<class AT>
  inline typename AT::value_type kahan_babuska_algorithm(const AT& x){
    typedef typename AT::value_type value_type;
    NumericDataTypeChecker<value_type>::certify();
     
    unsigned n = std::distance(x.begin(),x.end());
    
    value_type s = x[0],//value_type(),       //the sum
      c = value_type(),
      t;
    
    for(unsigned i = 1; i < n; ++i){
      t = s + x[i];
      if(Abs(s) >= Abs(x[i]))
	c += ((s-t)+x[i]);
      else
	c += ((x[i]-t)+s);
      s = t;
    }
    s += c; //correct after loop has been invoked!
    return s;
  }


   /**
   * \brief Babuska-Kahan summation for only 2 numbers
   */
  template<class T>
  inline T kahan_babuska_algorithm(const T& a, const T& b){
    NumericDataTypeChecker<T>::certify();
    
    T s = a,   //a assigned
      c = T(),
      t = s + b;

    if(Abs(s) >= Abs(b))
      	c += ((s-t)+b);
    else
      c += ((b-t)+s);
    s = t;
    
    s += c; //correct after all addition has been done!
    
    return s;
  }
 
  //! note that this compensated addition needs the correction to be done at the very end of all additions!
  template<class T>
  inline T& kahan_babuska_algorithm(const T& val, T& s, T& c){
    NumericDataTypeChecker<T>::certify();

    T t = s + val;

    if(Abs(s) >= Abs(val))
      	c += ((s-t)+val);
    else
      c += ((val-t)+s);
    s = t;
    
    return s;
  }



/**
   * \brief For 2 numbers, circumvent <TT> a+=b </TT> by defining the following assignment where the Kahan-Babuska summation, cf. [1], is applied
   */
  template<class T>
  inline T& kahan_babuska_assign(T& a, const T& b){
    NumericDataTypeChecker<T>::certify();
     
    T c = T(),
      t;

    t = a + b;
    if(Abs(a) >= Abs(b))
      	c += ((a-t)+b);
    else
      c += ((b-t)+a);
    a = t;
    
    a += c; //correct after all addition has been done!
    
    return a;
  }

  /**
   * \brief Only test: conventional sum
   */
   template<class AT>
   inline typename AT::value_type conventional_sum(const AT& a){
    typedef typename AT::value_type value_type;
    NumericDataTypeChecker<value_type>::certify();

    unsigned n = std::distance(a.begin(),a.end());
    
    value_type s = a[0]; // value_type();

    for(unsigned i = 1; i < n; ++i)
      s += a[i];
   
    return s;
   }


  //==================== NON-STAND-ALONE USAGE ================================

  /**
   * \brief Some algorithms to reduce rounding errors from addition.
   * This code snippet can be applied to 2 numbers, given \f$ s\f$ and \f$c\$, e.g. in conjunction with <BB>expression templates</BB> and class Add<T>.
   *
   * \tparam T precision of computation 
   * \tparam C select algorithm (currently: 'k','K' = Kahan, 'b','B' = Babuska-Kahan-algorithm
   * 
   * NOTE1: needs <I>temporary objects </I> to store the correction for the
   *       rounding error
   *
   * NOTE2: apparently, expressions are evaluated from left to right after their
   *        generation. Therefore the left elementary summand \f$ a \f$ already
   *        contains information about the sum \f$s\f$. The rouning error 
   *        correction is passed as reference and is altered in each +-operation
   *
   * NOTE3:  # of KAHAN-summations must agree with # of '+' occuring in 
   *          expression times length of assignment expression vector           
   */
  template<class T, char C> class BalancingAlgorithm;

  template<class T>
  class BalancingAlgorithm<T,'k'>{
  public:
    //! addition of 2 numbers
    static inline T& addition(const T& a, const T& b, T& s, T& c){
      NumericDataTypeChecker<T>::certify();
      
      T y, 
	t; 
      //                             or:
    
      s = a;                           // //s = a; //i.e. uncomment

     
      y = b - c;
      t = s + y;                      //t = a + y; 
      c = (t - s) - y;                // c = (t - a) - y;
      s = t; 
      
      //!apparently, this may affect the result (<TT>std::cout</TT> as well)!!
      //FancyMessages().nice_output("KAHAN-SUMMATION applied", 32,40);
      
      return s;
    }
  };

  template<class T>
  class BalancingAlgorithm<T,'K'>{
  public:
    static inline T& addition(const T& a, const T& b, T& s, T& c){
      return  BalancingAlgorithm<T,'k'>::addition(a,b,s,c);
    }
  };


  //!specialization -- Kahan-Babuska balanced addition 
  //! Needs correction of sum afterwards!, cf. CorrectRoundingErrorAfterwards, when used.
  template<class T>
  class BalancingAlgorithm<T,'b'>{
  public:
    static inline T& addition(const T& a, const T& b, T& s, T& c){
      NumericDataTypeChecker<T>::certify();
      
      T t;
      
      s = a;
     
      t = s + b;
      if(Abs(s) >= Abs(b))
	c += ((s - t) + b);
      else
	c += ((b - t) + s);
      
      s = t;
     
      //!apparently, this may affect the result (<TT>std::cout</TT> as well)!!
      //FancyMessages().nice_output("BABUSKA-KAHAN-SUMMATION applied", 32,40);

      
      return s;
    }
  
  };

  template<class T>
  class BalancingAlgorithm<T,'B'>{
  public:
    static inline T& addition(const T& a, const T& b, T& s, T& c){
      return BalancingAlgorithm<T,'b'>::addition(a,b,s,c);
    }
   
  };

  /**
   * \brief Since the Babuska-Kahan algorithm for addition performs the correction <I> after </I> the summation one has to perform the correction afterwards 
   *
   * default: do nothing during assignment
   */
  template<class T, char C>
  class CorrectRoundingErrorAfterwards{
  public:
    static inline void update_s(T& s, const T& c){} //do nothing
  
    static inline void assign(T& val, const T& s){} //do nothing
  
    template<class V>  //do nothing
    static inline void matvecmult_special(V& v, V& s, const V& c){}  
  
    template<class IT, class V>  //do nothing
    static inline void matvecmult_special_x(IT it, V& s, const V& c){} 

  };

  //! specialization for Kahan-Babuska algorithm
  template<class T>
  class CorrectRoundingErrorAfterwards<T,'b'>{
  public:
    
    static inline void update_s(T& s, const T& c){
      s += c;
    }
  
    static inline void assign(T& val, const T& s){
      val = s;
    }
  
    //! for matrix-vector multiplication which always contains '+'
    template<class V>  
    static inline void matvecmult_special(V& v, V& s, const V& c){
      adonis_assert((v.size() == s.size()) && (s.size() == c.size()));

      typedef typename V::value_type value_type;
      //! should be the same types
      adonis_assert(typeid(T) == typeid(value_type));

      for(size_t j = 0; j < v.size(); ++j){
	CorrectRoundingErrorAfterwards<value_type,'b'>::update_s(s[j],c[j]);
	CorrectRoundingErrorAfterwards<value_type,'b'>::assign(v[j],s[j]);
      }
    }

    template<class IT, class V>  
    static inline void matvecmult_special_x(IT it, V& s, const V& c){
      adonis_assert((s.size() == c.size()));

      typedef typename V::value_type value_type;
      //! should be the same types
      adonis_assert(typeid(T) == typeid(value_type));

      for(size_t j = 0; j < s.size(); ++j, ++it){ //increment iterator as well
	CorrectRoundingErrorAfterwards<value_type,'b'>::update_s(s[j],c[j]);
	CorrectRoundingErrorAfterwards<value_type,'b'>::assign(*it,s[j]);
      }
    }

  };

  //! u can use 'B' as well (instead of 'b')
  template<class T>
  class CorrectRoundingErrorAfterwards<T,'B'>{
  public:
    static inline void update_s(T& s, const T& c){
      CorrectRoundingErrorAfterwards<T,'b'>::update_s(s,c);
    }

    static inline void assign(T& val, const T& s){
      CorrectRoundingErrorAfterwards<T,'b'>::assign(val,s);
    }

    template<class V>  
    static inline void matvecmult_special(V& v, V& s, const V& c){
      typedef typename V::value_type value_type;
      CorrectRoundingErrorAfterwards<value_type,'b'>::matvecmult_special(v,s,c);
    }

    template<class IT, class V>  
    static inline void matvecmult_special_x(IT it, V& s, const V& c){
      typedef typename V::value_type value_type;
      CorrectRoundingErrorAfterwards<value_type,'b'>::matvecmult_special(it,s,c);
    }

  };

  //================================ END =======================================


  /**
   * \brief A smart addition: default case: standard summation
   */
  template<char C>  //!default += 
  class IntelligentArithmeticOperation{
  public:
    template<class T>
    static inline void pluseq(T& a, const T& b){
      a += b;
    }

    template<class T>
    static inline void minuseq(T& a, const T& b){
      a -= b;
    }

    template<class T>
    static inline void timeseq(T& a, const T& b){
      a *= b;
    }
    
    template<class T>
    static inline void divideeq(T& a, const T& b){
      adonis_assert(b != T());
      a /= b;
    }

    template<class T>  //! addition of 2 numbers
    static inline T add(const T& a, const T& b){
      return a+b;
    }

    template<class T>
    static inline T& add(const T& val, T& s, T& c){
      return s;
    }

  };

  //! specializations 
  template<>
  class IntelligentArithmeticOperation<'k'>{
  public: 
    template<class T>
    static inline void pluseq(T& a, const T& b){
      kahan_assign(a,b);
    }
  
    //those won't be affected neither by kahan nor babuska
    template<class T>
    static inline void minuseq(T& a, const T& b){
      a -= b;
    }

    template<class T>
    static inline void timeseq(T& a, const T& b){
      a *= b;
    }
    
    template<class T>
    static inline void divideeq(T& a, const T& b){
      adonis_assert(b != T());
      a /= b;
    }
    
    template<class T>  //! Kahan addition of 2 numbers
    static inline T add(const T& a, const T& b){
      return  kahan_algorithm(a,b);
    }

    template<class T>
    static inline T& add(const T& val, T& s, T& c){
      return kahan_algorithm(val,s,c);
    }
    
  };

  template<>
  class IntelligentArithmeticOperation<'K'>{
  public:
    template<class T>
    static inline void pluseq(T& a, const T& b){
      IntelligentArithmeticOperation<'k'>::pluseq(a,b);
    }
    
    template<class T>
    static inline void minuseq(T& a, const T& b){
      IntelligentArithmeticOperation<'k'>::minuseq(a,b);
    }
    
    template<class T>
    static inline void timeseq(T& a, const T& b){
      IntelligentArithmeticOperation<'k'>::timeseq(a,b);
    }
    
    template<class T>
    static inline void divideeq(T& a, const T& b){
      IntelligentArithmeticOperation<'k'>::divideeq(a,b);
    }

     template<class T>
     static inline T add(const T& a, const T& b){
       return IntelligentArithmeticOperation<'k'>::add(a,b);
     }

    template<class T>
    static inline T& add(const T& val, T& s, T& c){
      return IntelligentArithmeticOperation<'k'>::add(val,s,c);
    }
  };

  template<>
  class IntelligentArithmeticOperation<'b'>{
  public: 
    template<class T>
    static inline void pluseq(T& a, const T& b){
      kahan_babuska_assign(a,b);
    }
  
     //those won't be affected neither by kahan nor babuska
    template<class T>
    static inline void minuseq(T& a, const T& b){
      a -= b;
    }

    template<class T>
    static inline void timeseq(T& a, const T& b){
      a *= b;
    }
    
    template<class T>
    static inline void divideeq(T& a, const T& b){
      adonis_assert(b != T());
      a /= b;
    }

     template<class T>
     static inline T add(const T& a, const T& b){
       return  kahan_babuska_algorithm(a,b);
    }

    template<class T>
    static inline T& add(const T& val, T& s, T& c){
      return kahan_babuska_algorithm(val,s,c);
    }

  };

  template<>
  class IntelligentArithmeticOperation<'B'>{
  public: 
    template<class T>
    static inline void pluseq(T& a, const T& b){
      IntelligentArithmeticOperation<'b'>::pluseq(a,b);
    }

    template<class T>
    static inline void minuseq(T& a, const T& b){
      IntelligentArithmeticOperation<'b'>::minuseq(a,b);
    }

    template<class T>
    static inline void timeseq(T& a, const T& b){
      IntelligentArithmeticOperation<'b'>::timeseq(a,b);
    }

    template<class T>
    static inline void divideeq(T& a, const T& b){
      IntelligentArithmeticOperation<'b'>::divideeq(a,b);
    }
  
    template<class T>
    static inline T add(const T& a, const T& b){
      return IntelligentArithmeticOperation<'b'>::add(a,b);
    }
    
    template<class T>
    static inline T& add(const T& val, T& s, T& c){
      return IntelligentArithmeticOperation<'b'>::add(val,s,c);
    }

  };

} //end namespace 

#endif 
