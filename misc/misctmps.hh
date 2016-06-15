#ifndef MISCELLANEOUS_TMPS_THAT_MAY_AMELIORATE_YOUR_LIFE_HH
#define MISCELLANEOUS_TMPS_THAT_MAY_AMELIORATE_YOUR_LIFE_HH

#include <iostream>
#include <typeinfo>

#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"
#include "../common/numerictypechecker.hh"
#include "../common/isclass.hh"
#include "../common/typeadapter.hh"

#include "../common/smartassign.hh"

namespace Adonis{

  /**
   * \brief Note that for \f$ a \in \mathbb{K}\f$ and for \f$ n \in \mathbb{N}\f$, we have \f$ a^{-n} = \frac{1}{a^n}.\f$
   */
  template<bool B, class T> class InverseExponentiation;

  template<class T>
  class InverseExponentiation<true,T>{
  public:
    static inline T invexp(const T& e){
      adonis_assert(e != T());
      NumericDataTypeChecker<T>::certify(); //make sure it's a numeric value

      return 1./e;
    }
  };

  template<class T>
  class InverseExponentiation<false,T>{
  public:
    static inline T invexp(const T& e){
      return e;
    }
  };
  

/**
 * \brief computes the natural power of a given value, i.e. powers whose exponent is taken from \f$ \mathbb{N}.\f$ 
 *
 * \code
      std::cout << ntimes<3>(-5.75) << std::endl;
 * \endcode
 * should return -190.11
 * \code
      std::complex<double> z(-3.,4.);
      std::cout << ntimes<2>(z) << std::endl;
 * \endcode
 *should return (-7,-24)
 */
  template<int N, class X>
  class NTimes{
  public:
    static inline X raise_to_the_power_of_N(const X& x){
      return x*NTimes<N-1,X>::raise_to_the_power_of_N(x);
    }
  };

  //end of recursion
  template<class X>
  class NTimes<0,X>{
  public:
    static inline X raise_to_the_power_of_N(const X& x){
      return X(1);
    }
  };

  /**
   * \brief Convenient function for using natural exponentiation
   *
   * NOTE: You can use negative N as well here, since it is checked at compile time whether \f$ a^N\f$ or \f$ \frac{1}{a^N}\$ is applied.
   *
   * 
   */
  template<int N, class X>
  inline X ntimes(const X& x){
    //return NTimes<N,X>::raise_to_the_power_of_N(x);
    return InverseExponentiation<(N < 0),X>::invexp(NTimes<((N < 0) ? -N : N),X>::raise_to_the_power_of_N(x));
    
  }
  


  /** \brief Fill a dynamic array created previously by e.g. T* da = new T [n];    
   * NOTE: iterators and pointers don't allow range checks!
   * USAGE EXAMPLE:
   *\code
     const size_t dip = 6;
     double* ptr = new double [dip];
     ManipulateDynamicArray<double,size_t,0,dip>::fill(ptr,42.);

     //output option
     ManipulateDynamicArray<double,size_t,0,dip>::output(ptr);
     
     delete[] ptr;
   *\endcode
   */
  template<class T, class INT, INT I, INT N>
  class ManipulateDynamicArray{
  private:
    enum{index_ = (I+1) != N};
    
  public:
    static inline void fill(T* p, const T& d){
      *p = d;   //
      ++p;      //do the pointer arithmetics
      ManipulateDynamicArray<T,INT,((index_)?(I+1):N),N>::fill(p,d);//recursion
    }
  
    template<class ITER>
    static inline void fill(T* p, ITER i1, ITER i2){
      smart_assign(*p,*i1); //*p = static_cast<T>(*i1);
      //std::cout << "*p = "<<*p<< std::endl;
      ++i1;
      ++p;
      ManipulateDynamicArray<T,INT,((index_)?(I+1):N),N>::fill(p,i1,i2);
    }
  
    
    template<class V>
    static inline void fill(T& cont1, const V& cont2){
      cont1[I] = cont2[I];
      ManipulateDynamicArray<T,INT,((index_)?(I+1):N),N>::fill(cont1,cont2);
    }

    static inline void output(T* p){
      std::cout << *p << "  ";
      ++p;
      ManipulateDynamicArray<T,INT,((index_)?(I+1):N),N>::output(p);//recursion
    }
  
    template<class ITER>
    static inline void output(ITER i1, ITER i2){
      std::cout <<  *i1 << "  ";
      ++i1;
      ManipulateDynamicArray<T,INT,((index_)?(I+1):N),N>::output(i1,i2);
    }
      
  };
  
  //end of recursion
  template<class T, class INT, INT N>
  class ManipulateDynamicArray<T,INT,N,N>{
  public:
    static inline void fill(T* p, const T& d){}  //do nothing
    
    template<class ITER>
    static inline void fill(T* p, ITER i1, ITER i2){
      adonis_assert(i1 == i2);      //make sure iterator i1 has reached i2
    }
  
    template<class V>
    static inline void fill(T& cont1, const V& cont2){}

    static inline void output(T* p){
      std::cout << std::endl;
    }
  
    template<class ITER>
    static inline void output(ITER i1, ITER i2){
      adonis_assert(i1 == i2);
      std::cout <<  std::endl;
      
    }
  };
  

  
 
  /**
   * \brief Choose at compile time the value type of a container. This is very helpful because DUNE's field types don't have a value_type (at least not up to now) 
   */
  template<bool B, class C> class ValueTypeWrapper;

  template<class C>
  class ValueTypeWrapper<true,C>{  //a DUNE field type
  public:
    typedef typename C::field_type value_type;
  };

  template<class C>
  class ValueTypeWrapper<false,C>{  //no field type
  public:
    typedef typename C::value_type value_type;
  };

  
  /**
   * \brief Fill a linear container from, e.g. a given pointer (array) 
   */
  template<class C, int I, int N>
  class FillLinearContainer{
  private:
    enum{index_ = (I+1)!= N};

  public :
    typedef typename TypeAdapter<typename ValueTypeWrapper<IsADuneFieldContainer<C>::Value,C>::value_type>::Type Type;
    typedef typename C::iterator iterator;

    static inline void fill_me(iterator it, Type* ptr){
      *it = *ptr;
      ++it;
      ++ptr;
      
      FillLinearContainer<C,((index_)?(I+1):N),N>::fill_me(it,ptr);
    }
  };

  
  template<class C, int N>
  class FillLinearContainer<C,N,N>{
  public:
    // typedef typename ValueTypeWrapper<IsADuneFieldContainer<C>::Value,C>::value_type value_type;
    typedef typename TypeAdapter<typename ValueTypeWrapper<IsADuneFieldContainer<C>::Value,C>::value_type>::Type Type;
    typedef typename C::iterator iterator;

    //do nothing
    static inline void fill_me(iterator it, Type* ptr){}
  };


  //!no convenient function here
  


  /**
   * \brief Calculate the dot product of two random access containers providing an []-operator. This template meta program might be beneficial, e.g. when a dot product bewteen a static container (size is known at compile time) with a dynamic container is going to be performed or the dimension of the containers is known at compile time anyway.
   * 
   */
  template<int N, class T>
  class DotProduct{
  public:
    template<class V, class W>
    static T result(const V& v, const W& w){
      return v[N-1]*Conj(w[N-1]) + DotProduct<N-1,T>::result(v,w);
    }
  };

  //end of recursion
  template<class T>
  class DotProduct<1,T>{
  public:
    template<class V, class W>
    static T result(const V& v, const W& w){
      return v[0]*Conj(w[0]);
    }
  };
  
  //! specialization in the case of zero size; just return 0 {rather unorthodox
  //! way which satisfies my demands ;) }
  template<class T>
  class DotProduct<0,T>{
  public:
    template<class V, class W>
    static T result(const V& v, const W& w){
      return T();
    }
  };

  //convenient function
  template<int N, class V, class W>
  inline typename  ValueTypeWrapper<IsADuneFieldContainer<W>::Value,W>::value_type dot_product(const V& v, const W& w){
    adonis_assert(std::distance(v.begin(),v.end()) == std::distance(w.begin(),w.end()) && N == static_cast<int>(std::distance(v.begin(),v.end())));

    typedef typename ValueTypeWrapper<IsADuneFieldContainer<W>::Value,W>::value_type value_type;

    return DotProduct<N,value_type>::result(v,w);
  }
  

} //end namespace

#endif
