#ifndef NUMERICAL_TYPE_PROMOTION_HH
#define NUMERICAL_TYPE_PROMOTION_HH

#include <complex>
#include "nihil.hh"

#if USE_CPPAD
#include <cppad/cppad.hpp>
#endif

namespace Adonis{
  
 
  /**
   * \brief The type promotion guarantees that the right type is selected as
   *  return type in binary relation L ° R, e.g.
   *   double ° float --> double,
   *   complex<float> ° double --> complex<double>,
   * i.e. in mixed operations we switch to the next "higher" type.
   *
   * Note: C/C++ rules including integral types haven't been considered 
   * so far, e.g.
   *   int ° double --> double  
   *
   * Decisions about the ReturnType of, e.g. 3 different types can be obtained 
   * via nested application of NumericalTypePromotion:
   *
   * \code 
      typedef NumericalTypePromotion<long double, NumericalTypePromotion<double,std::complex<float> >::ReturnType>::ReturnType value_type;
   * \endcode
   * Here <I> typeid(value_type).name() </I> gives <TT> St7complexIeE</TT>, i.e. std::complex<long double> as expected ;)
   */
  template<class L, class R>
  class NumericalTypePromotion{
  public:
    typedef Nihil ReturnType; //throw this if something isn't specified yet
  };
  

  //! o.k. if two types are equal then the one type is also the return type
  template<class S>
  class NumericalTypePromotion<S,S>{
  public:
    typedef S ReturnType;
  };




  //L = float
  template<>
  class  NumericalTypePromotion<float,double>{
  public:
    typedef double ReturnType;
  };
  
  template<>
  class  NumericalTypePromotion<float,long double>{
  public:
    typedef long double ReturnType;
  };
  

  //L = double
  template<>
  class  NumericalTypePromotion<double,float>{
  public:
    typedef double ReturnType;
  };

  template<>
  class  NumericalTypePromotion<double,long double>{
  public:
    typedef long double ReturnType;
  };


  //L = long double 
  template<>
  class  NumericalTypePromotion<long double,float>{
  public:
    typedef long double ReturnType;
  };

  template<>
  class  NumericalTypePromotion<long double,double>{
  public:
    typedef long double ReturnType;
  };


 

  //////////////////////////////////////////////////////////////////////////
  ///////////////////////// COMPLEX VERSION ////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  
   //L = complex<float>
  template<>
  class  NumericalTypePromotion<std::complex<float>,std::complex<double> >{
  public:
    typedef std::complex<double> ReturnType;
  };
  
  template<>
  class  NumericalTypePromotion<std::complex<float>,std::complex<long double> >{
  public:
    typedef std::complex<long double> ReturnType;
  };
  
  
  //L = std::complex<double>
  template<>
  class  NumericalTypePromotion<std::complex<double>,std::complex<float> >{
  public:
    typedef std::complex<double> ReturnType;
  };

  template<>
  class  NumericalTypePromotion<std::complex<double>,std::complex<long double> >{
  public:
    typedef std::complex<long double> ReturnType;
  };


  //L = std::complex<long double>
  template<>
  class  NumericalTypePromotion<std::complex<long double> ,std::complex<float> >{
  public:
    typedef std::complex<long double> ReturnType;
  };

  template<>
  class  NumericalTypePromotion<std::complex<long double>,std::complex<double> >{
  public:
    typedef std::complex<long double> ReturnType;
  };


  //!Mixed complex and basic numerical types
  template<class T>
  class  NumericalTypePromotion<std::complex<T>, std::complex<T> >{
  public:
    typedef std::complex<T> ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<std::complex<T1>,T2>{
  public:
    typedef std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<T1, std::complex<T2> >{
  public:
    typedef std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<std::complex<T1>, std::complex<T2> >{
  public:
    typedef std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> ReturnType;
  };



  /////////////////////////////////////////////////////////////////////////////
  ////////////////////////// OTHER SPECIALIZATIONS ////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  //! CPPAD specialisation
#if USE_CPPAD
  //! in case we have the same CppAD Type
  template<class T>
  class  NumericalTypePromotion<CppAD::AD<T>,CppAD::AD<T> >{
  public:
    typedef CppAD::AD<T> ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<T1,CppAD::AD<T2> >{
  public:
    typedef CppAD::AD<typename NumericalTypePromotion<T1,T2>::ReturnType> ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<CppAD::AD<T1>,T2>{
  public:
    typedef CppAD::AD<typename NumericalTypePromotion<T1,T2>::ReturnType> ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<CppAD::AD<T1>,CppAD::AD<T2> >{
  public:
    typedef CppAD::AD<typename NumericalTypePromotion<T1,T2>::ReturnType> ReturnType;
  };

  /////////////////////////////////////////////////////////////////////////////
  ///////////////////// WEIRD COMPOSITIONS ////////////////////////////////////
  //!weird type promotions, perhaps sometimes necessary
  template<class T1, class T2>
  class  NumericalTypePromotion<std::complex<CppAD::AD<T1> >,CppAD::AD<T2> >{
  public:
    typedef std::complex<CppAD::AD<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<CppAD::AD<T1>, std::complex<CppAD::AD<T2> > >{
  public:
    typedef std::complex<CppAD::AD<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };

 
  
  template<class T1, class T2>
  class  NumericalTypePromotion<CppAD::AD<std::complex<T1> >,CppAD::AD<T2> >{
  public:
    typedef CppAD::AD<std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };

  
  template<class T1, class T2>
  class  NumericalTypePromotion<CppAD::AD<T1>, CppAD::AD<std::complex<T2> > >{
  public:
    typedef CppAD::AD<std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<CppAD::AD<std::complex<T1> >, CppAD::AD<std::complex<T2> > >{
  public:
    typedef CppAD::AD<std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };

  template<class T>
  class  NumericalTypePromotion<CppAD::AD<std::complex<T> >, CppAD::AD<std::complex<T> > >{
  public:
     typedef CppAD::AD<std::complex<T> > ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<std::complex<T1>, CppAD::AD<std::complex<T2> > >{
  public:
    typedef CppAD::AD<std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };

  template<class T1, class T2>
  class  NumericalTypePromotion<CppAD::AD<std::complex<T1> >, std::complex<T2> >{
  public:
    typedef CppAD::AD<std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };


  template<class T1, class T2>
  class  NumericalTypePromotion<CppAD::AD<T1>, std::complex<T2> >{
  public:
    typedef CppAD::AD<std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };

   template<class T1, class T2>
   class  NumericalTypePromotion<std::complex<T1>, CppAD::AD<T2> >{
  public:
    typedef CppAD::AD<std::complex<typename NumericalTypePromotion<T1,T2>::ReturnType> > ReturnType;
  };


#endif
  

 
}

#endif
