#ifndef BOUNDARY_TEMPERATURE_DISTRIBUTION_HH
#define BOUNDARY_TEMPERATURE_DISTRIBUTION_HH

#include <iostream>
#include <cmath>

#include "../massactionkinetics/physicalconstants.hh"
#include "../../common/globalfunctions.hh"
#include "../../common/adonisassert.hh"
#include "../common/typeadapter.hh"
#include "../../misc/misctmps.hh"

namespace Adonis{

  /**
   * \brief Boundary value profile which has been used in [1] and which 
   *        has been kindly provided by Ch. Frouzakis
   */
  template<class T>
  inline T Frouzakis_profile(const T& Tmin, const T& Tmax, const T& x0, const T& x){
    return ( (x <= x0) ? Tmin*(1.1*cos(PhysicalConstants<T>::Pi*(x/0.5+1)) +2.1) : Tmax );
  }


  /**
   * \brief Natural cubic spline-like temperature distribution in four points
   * (0,Tmin), (x0,Tmax), (x0+b,Tmax-Tmin), (length,Tmin).
   * Values above Tmax and below Tmin, respectively, have been trimmed to 
   * preserve physically meaningful values for the cubic spline approximation.
   */
  template <class T>
  inline T natural_cubic_spline_profile(const T& Tmin, const T& Tmax, const T& x0, const T& x, const T& b, const T& length){
    T h0 = x0,
      h1 = b,
      h2 = length-(x0+b),
      b0 = (Tmax-Tmin)/h0,
      b1 = -Tmin/h1,
      b2 = (2*Tmin - Tmax)/h2,
      v1 = 2*(h0+h1),
      v2 = 2*(h1+h2),
      u1 = 6*(b1-b0),
      u2 = 6*(b2-b1);

    //2x2 'tridiagonal' system: note: z0 = z3 = 0 (natural spline=
    //    v1   h1      z1     u1
    //              ·      =  
    //    h1   v2      z2     u2

    T det = 1./(v1*v2 - ntimes<2>(h1)),
      z0 = 0.,  //natural spline cdt
      z1 = det*(v2*u1 - h1*u2),
      z2 = det*(-h1*u1 + v1*u2),
      z3 = 0.;  //natural spline cdt

   
    // the piecewise cubic polynomial reads as \f$ S_i(x) = \frac{z_{i+1}}{6h_i}(x-x_i)^3 + \frac{z_i}{6hi}(x_{i+1}-x)^3 + (\frac{y_{i+1}}{h_i} - \frac{z_{i+1}h_i}{6})(x-x_i) + (\frac{y_i}{h_i} - \frac{z_ih_i}{6})(x_{i+1}-x), \qquad i = 0,\ldots,n-1\f$
    T spln = 2./3.*Tmax; //just initialize variable
    if (x <= x0){
      spln = z1/(6.*h0)*ntimes<3>(x-0.) + z0/(6.*h0)*ntimes<3>(x0-x) + (Tmax/h0 - z1/6.*h0)*(x-0.) + (Tmin/h0 - h0/6.*z0)*(x0-x);
      //spln = (spln < Tmin) ? Tmin : spln;
      spln = (spln < Tmin) ? Tmin : ((spln > Tmax) ? Tmax : spln);
    }
    else if ((x0 <= x) && (x <= (x0+b))){
      spln = z2/(6.*h1)*ntimes<3>(x-x0) + z1/(6.*h1)*ntimes<3>(x0+b-x) + ((Tmax-Tmin)/h1 - z2/6.*h1)*(x-x0) + (Tmax/h1 - h1/6.*z1)*(x0+b-x);
       //spln = (spln < Tmin) ? Tmin : spln;
      spln = (spln < Tmin) ? Tmin : ((spln > Tmax) ? Tmax : spln);
    }
    else if (((x0+b) <= x) && (x <= length)){
      spln = z3/(6.*h2)*ntimes<3>(x-(x0+b)) + z2/(6.*h2)*ntimes<3>(length-x) + (Tmin/h2 - z3/6.*h2)*(x-(x0+b)) + ((Tmax-Tmin)/h2 - h2/6.*z2)*(length-x);
      //spln = (spln < Tmin) ? Tmin : spln;
      spln = (spln < Tmin) ? Tmin : ((spln > Tmax) ? Tmax : spln);
    }
    else{
      //do nothing
    }
    return spln;
    
  }

  //3 points only 
  template<class T >
  inline T natural_cubic_spline_profile(const T& Tmin, const T& Tmax, const T& x0, const T& x, const T& length){
    T h0 = x0,
      h1 = length - x0,
      b0 = (Tmax-Tmin)/h0,
      b1 = (Tmin-Tmax)/h1,
      v1 = 2*(h0+h1),
      u1 = 6*(b1-b0),
      z0 = 0.,
      z1 = u1/v1,
      z2 = 0.;

    T spln = 2./3.*Tmax; //just initialize variable
    if(x <= x0){
      spln = z1/(6.*h0)*ntimes<3>(x-0.) + z0/(6.*h0)*ntimes<3>(x0-x) + (Tmax/h0 - z1/6.*h0)*(x-0.) + (Tmin/h0 - h0/6.*z0)*(x0-x);
      //spln = (spln < Tmin) ? Tmin : spln;
    }
    else if((x0 <= x) && (x <= length)){
      spln =  z2/(6.*h1)*ntimes<3>(x-x0) + z1/(6.*h1)*ntimes<3>(length-x) + (Tmin/h1 - z2/6.*h1)*(x-x0) + (Tmax/h1 - h1/6.*z1)*(length-x);
      //spln = (spln < Tmin) ? Tmin : spln;
    }
    else{
      //do nothing
    }
    return spln;
  }


  /**
   * \brief Wall temperature functors
   */
  template<char C, class T> class BoundaryWallTemperature;

  //!partial specialization
  template<class T>
  class BoundaryWallTemperature<'h',T>{ //! hyperbolic tangent profile
  public:
    typedef typename TypeAdapter<T>::BaseType BaseType;

    static const char Value = 'h';

    static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){ //length is dummy
      return ( Tmin + 0.5*(Tmax-Tmin)*(1+tanh(a*(x-x0))) );
      //! normalized hyperbolic tangent profiles between Tmin and Tmax centered 
      //! at x0, due to Ch. Frouzakis, ETHZ, private correspondence
      //return ( 1 + (Tmax-Tmin)*(1 + tanh((x-x0)/a)) );
    }
    
  };

   template<class T>
   class BoundaryWallTemperature<'H',T>{ //! hyperbolic tangent profile
   public:
     typedef typename TypeAdapter<T>::BaseType BaseType;
     
     static const char Value = 'H';

     static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){ //length is dummy
       return BoundaryWallTemperature<'h',T>::temperature(Tmin,Tmax,x0,x,a,length);
     }
     
   };


  //!Wall temperature profile due to Ch. Frouzakis, ETHZ, private corresp.
  template<class T>
  class BoundaryWallTemperature<'f',T>{ //! Frouzakis profile
  public:
    typedef typename TypeAdapter<T>::BaseType BaseType;
    
    static const char Value = 'f';
    
    static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){ //a and length are dummies
      return Frouzakis_profile(Tmin,Tmax,x0,x); 
    }
    
  };
  
   template<class T>
   class BoundaryWallTemperature<'F',T>{ //! Frouzakis profile
   public:
     typedef typename TypeAdapter<T>::BaseType BaseType;

     static const char Value = 'F';

     static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){
       return BoundaryWallTemperature<'f',T>::temperature(Tmin,Tmax,x0,x,a,length);
     }
     
   };



  //! my profile, based on Ch. Frouzakis's wall temperature distribution
   template<class T>
   class BoundaryWallTemperature<'m',T>{ 
   public:
    typedef typename TypeAdapter<T>::BaseType BaseType;
    
    static const char Value = 'm';
    
     static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){ //length is a dummy 
      T res = T();
      if(x < x0) 
	res = Frouzakis_profile(Tmin,Tmax,x0,x);
      else if((x0 <= x) && (x <= (x0 + a))) //peak or plateau, depending on a
	res = Tmax;
      else if(x > x0 + a)
	res = Frouzakis_profile(Tmin,Tmax,x0,-x); //minus x here for reflection
      else{
	// do nothing
      }
      return res;
    }
    
  };
  
   template<class T>
   class BoundaryWallTemperature<'M',T>{ //! my profile
   public:
     typedef typename TypeAdapter<T>::BaseType BaseType;

     static const char Value = 'M';

     static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){
       return BoundaryWallTemperature<'m',T>::temperature(Tmin,Tmax,x0,x,a,length);
     }
   };
  
  
   template<class T>
   class BoundaryWallTemperature<'a',T>{  //profile roughly resembling an 'A'
   public:
    typedef typename TypeAdapter<T>::BaseType BaseType;
    
    static const char Value = 'a';
    
     static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){ //length is here now an additional parameter, not the x-length
      T res = T();
      if(x < (x0-length)) 
	res = Tmin;
      else if(((x0-length) <= x) && (x < x0)) //peak or plateau, depending on a
	res = 0.9*Tmax;  //90% of Tmax
      else if((x0 <= x ) && (x <= (x0 + a)))
	res = Tmax;
      else if (((x0 + a) < x) && (x < (x0 + a + length)))
	res =  0.9*Tmax;  //90% of Tmax
      else if(x >= (x0 + a + length))
	res = Tmin;
      else{
	// do nothing
      }
      return res;
     }
   };
  
   template<class T>
   class BoundaryWallTemperature<'A',T>{ 
   public:
     typedef typename TypeAdapter<T>::BaseType BaseType;

     static const char Value = 'A';

     static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){
       return BoundaryWallTemperature<'a',T>::temperature(Tmin,Tmax,x0,x,a,length);
     }
   };
  


  template<class T>
  class BoundaryWallTemperature<'e',T>{ //! exponential profile
  public:
    typedef typename TypeAdapter<T>::BaseType BaseType;

    static const char Value = 'e';

    static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){ //length is a dummy
      return ( Tmin + (Tmax-Tmin)*exp(-1.e+05*ntimes<2>(x-x0)) );
    }
    
  };

  template<class T>
  class BoundaryWallTemperature<'E',T>{ //! exponential  profile
  public:
    typedef typename TypeAdapter<T>::BaseType BaseType;

    static const char Value = 'E';

    static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){
      return BoundaryWallTemperature<'e',T>::temperature(Tmin,Tmax,x0,x,a,length);
    }
     
  };



  template<class T>
  class BoundaryWallTemperature<'s',T>{ //! natural cubic spline-like profile
  public:
    typedef typename TypeAdapter<T>::BaseType BaseType;

    static const char Value = 's';

    static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){ //length is a dummy
      return (natural_cubic_spline_profile(Tmin,Tmax,x0,x,a,length) );
    }
    
  };

   template<class T>
   class BoundaryWallTemperature<'S',T>{ //! exponential  profile
   public:
     typedef typename TypeAdapter<T>::BaseType BaseType;

     static const char Value = 'S';

     static inline T temperature(const T& Tmin, const T& Tmax, const T& x0, const T& x, const BaseType& a, const BaseType& length){
       return BoundaryWallTemperature<'s',T>::temperature(Tmin,Tmax,x0,x,a,length);
     }
     
   };



  /**
   * \brief Thx goes to Malte Braack for sharing his data with me.
   */
  template<class T, char TEMPPROFTYPE = 'e'>
  class BoundaryTemperatureDistribution{
  public:
    typedef typename TypeAdapter<T>::BaseType BaseType;

    enum{
       Value = BoundaryWallTemperature<TEMPPROFTYPE,T>::Value
    };
    
   
    BoundaryTemperatureDistribution(const T& Tmin = T(), const T& Tmax = T(), const T& x0 = T(), const BaseType& a = BaseType(1), const BaseType& length = BaseType(1)):Tmin_(Tmin),Tmax_(Tmax),x0_(x0),a_(a),length_(length){
      adonis_assert(Tmin <= Tmax);
    }

    void initialize(const T& Tmin, const T& Tmax, const T& x0, const BaseType& a = 1, const BaseType& length = 1){
      adonis_assert(Tmin <= Tmax);
      Tmin_ = Tmin;
      Tmax_ = Tmax;
      x0_ = x0;
      a_ = a;
      length_ = length;
    }

    //! actually, the temperature distribution only hinges on spatial coord. x
    T operator()(const T& x){
      return ( BoundaryWallTemperature<TEMPPROFTYPE,T>::temperature(Tmin_,Tmax_, x0_,x,a_,length_) );
    }

    T operator()(const T& x, const T& y){
      return ( BoundaryWallTemperature<TEMPPROFTYPE,T>::temperature(Tmin_,Tmax_, x0_,x,a_,length_) );
    }

  private:
    T Tmin_, Tmax_, x0_;
    BaseType a_, length_; 
  };

} //end namespace 

#endif
