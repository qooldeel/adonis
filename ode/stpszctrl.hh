#ifndef STEPSITE_CONTROL_FOR_1_STEP_METHOD_HH
#define STEPSITE_CONTROL_FOR_1_STEP_METHOD_HH

#include <cmath>

#include "../common/typeadapter.hh"
#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../misc/misctmps.hh"
#include "diffequtilities.hh"

namespace Adonis{

  //! maybe better choices may be made
  template<class T>
  inline void enforce_bound_on_step_size(T& h, const T& hmin, const T& hmax){
    if(h < hmin)
      h = hmin;
    if(h > hmax) 
      h = hmax; 
    //same as h = std::max(hmin,std::min(h,hmax)); 
  }

  /**
   * \brief a nice step size control
   *
   *  \tparam V Container type to compute with
   *  \tparam NORM character specifying the norm ('1', '2' and 'i'(inf))
   *  \tparam B boolean value (true  = use adaptive stepsize, 
   *                           false = use equidistant stepsize)
   */
  template<class V, char NORM, bool B> class StepSizeControl4OneStepMethod;


  /**
   * \brief Equidistant stepsize hestim is used to increment the time. The 
   *  maximum number of steps is computed via \f$ \lceil (T - t_0)/h_{\mathrm{estim}}\f$
   */
  template<class V, char NORM>    //! no stepsize control, use equidistant steps
  class StepSizeControl4OneStepMethod<V,NORM,false>{
  public:
    static const bool Value = false; 
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef Norm<NORM,BaseType> NormType;
   
    
    template<class T, class C>
    static inline value_type& compute_new_step_size(V& err, const V& y_n, const  V& y_prev, const V& y_prev_prev, value_type& h, value_type& h_prev, const value_type& hmin, const value_type& hmax, const T& rtol, const C& Cstszctrl, std::size_t accuracyOfMethod){
      enforce_bound_on_step_size(h,hmin,hmax);
      return h;   //just give back estimate
    }


    template<class T, class H>
    static inline unsigned max_iters(unsigned m, const T& t0, const T& tend, const H& hestim){
      if(Abs(hestim) <= H())
	return m;               //avoid a possible division by zero
      else
	return Ceil((tend - t0)/hestim);
    }

    //! simplistic stepsize control
    static inline const value_type& compute_new_step_size(const value_type& kn, bool& isNewtonConverged, const value_type& hmin, const value_type& hmax, std::size_t accuracyOfMethod = 0){
      enforce_bound_on_step_size(kn,hmin,hmax);
      return kn;
    }


    //! normal stepsize 
    template<class T>
    static inline value_type& compute_standard_step_size(V& err, const V& y_n, const  V& y_prev, const V& y_prev_prev, value_type& h, const value_type& hmin, const value_type& hmax, const T& rtol, std::size_t accuracyOfMethod){
      enforce_bound_on_step_size(h,hmin,hmax);
      return h;
    }
  };


  /**
   * \brief Adaptive step size computation. Max. number of iterations is just
   *  prescribed number (inlike in the equidist. case in which this number will
   *  be overwritten)
   */
  template<class V, char NORM>   //! use stepsizectrl
  class StepSizeControl4OneStepMethod<V,NORM,true>{
  public:
    static const bool Value = true; 
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef Norm<NORM,BaseType> NormType;

    template<class T, class C>
    static inline value_type& compute_new_step_size(V& err, const V& y_n, const  V& y_prev, const V& y_prev_prev, value_type& h, value_type& h_prev, const value_type& hmin, const value_type& hmax, const T& rtol, const C& Cstszctrl, std::size_t accuracyOfMethod){
     
      err = y_prev - y_prev_prev;
      BaseType errnorm = NormType::norm(err);
      if(errnorm < rtol){
	if(errnorm > BaseType()) //avoid division by zero
	  h = (rtol*h_prev)/(Cstszctrl*errnorm); //here: C = 1. (in gen. C > 0)
      }
      else{ 
	h  *= 0.5;
      }
      enforce_bound_on_step_size(h,hmin,hmax);
      h_prev = h;
      
      return h;
    }


    template<class T, class H>
    static inline unsigned max_iters(unsigned m, const T& t0, const T& tend, const H& hestim){
      return m;
    }

    //!simplistic stepsize control
    static inline value_type& compute_new_step_size(value_type& kn, bool& isNewtonConverged,const value_type& hmin, const value_type& hmax, std::size_t accuracyOfMethod){
      if(isNewtonConverged)
	kn *= 2;
      else 
	kn *= 0.5;
      //! don't let stepsize become too small or large
      enforce_bound_on_step_size(kn,hmin,hmax);
      //reset for next iteration
      isNewtonConverged = false;
      return kn;
    }

    //! normal stepsize
    template<class T>
    static inline value_type& compute_standard_step_size(V& err, const V& y_n, const  V& y_prev, const V& y_prev_prev, value_type& h, const value_type& hmin, const value_type& hmax, const T& rtol, std::size_t accuracyOfMethod){
      err = y_prev - y_prev_prev;
      BaseType errnorm = NormType::norm(err);
      //if(errnorm < rtol){
      if(errnorm > BaseType()){ //avoid division by zero
	h = 0.9*h*pow(rtol/errnorm, 1./(accuracyOfMethod+1));
      }
      //}
      else{ 
	h *= 0.25;  //quarter it otherwise
      }
      enforce_bound_on_step_size(h,hmin,hmax);
    return h;
    }
  };


  /**
   * \brief Select a stepsize control. At the moment, two are at choice:
   *  'E' stepsize control due to Johnson et al. and a simplistic one ('s')
   *  
   * NOTE: if SZCTRL::Value is <I>false</I> then equidistant stepsize will be 
   *       used in all cases.
   */
  template<class SZCTRL, char C> class ChooseStepsizeControl;

  //! stepsize control due to Johnson et all 
  template<class SZCTRL> 
  class ChooseStepsizeControl<SZCTRL,'E'>{
  public:
    template<class V, class T, class C>
    static inline typename V::value_type& compute_new_step_size(V& err, const V& y_n, V& y_prev, V& y_prev_prev, typename V::value_type& h, typename V::value_type& h_prev, const typename V::value_type& hmin, const typename V::value_type& hmax, const T& rtol, const C& Cstszctrl, bool& isNewtonConverged, std::size_t accuracyOfMethod){
      return SZCTRL::compute_new_step_size(err,y_n,y_prev,y_prev_prev,h,h_prev,hmin,hmax,rtol,Cstszctrl,accuracyOfMethod);
    }

    template<class V>
    static inline void overwrite_y_prev_prev(const V& yp, V& ypp){
      ypp = yp;
    }
  };

 //! standard stepsize control 
  template<class SZCTRL> 
  class ChooseStepsizeControl<SZCTRL,'n'>{
  public:
    template<class V, class T, class C>
    static inline typename V::value_type& compute_new_step_size(V& err, const V& y_n, V& y_prev, V& y_prev_prev, typename V::value_type& h, typename V::value_type& h_prev, const typename V::value_type& hmin, const typename V::value_type& hmax, const T& rtol, const C& Cstszctrl, bool& isNewtonConverged, std::size_t accuracyOfMethod){
      return SZCTRL::compute_standard_step_size(err,y_n,y_prev,y_prev_prev,h,hmin,hmax,rtol,accuracyOfMethod);
    }
  
    template<class V>
    static inline void overwrite_y_prev_prev(const V& yp, V& ypp){
      ypp = yp;
    }
  };
  
  //! simple stepsize control
  template<class SZCTRL> 
  class ChooseStepsizeControl<SZCTRL,'s'>{
  public:
    template<class V, class T, class C>
    static inline typename V::value_type& compute_new_step_size(V& err, const V& y_n, V& y_prev, V& y_prev_prev, typename V::value_type& h, typename V::value_type& h_prev, const typename V::value_type& hmin, const typename V::value_type& hmax, const T& rtol, const C& Cstszctrl, bool& isNewtonConverged,std::size_t accuracyOfMethod){
      return SZCTRL::compute_new_step_size(h,isNewtonConverged,hmin,hmax,accuracyOfMethod);
    }
   
    template<class V>
    static inline void overwrite_y_prev_prev(const V& yp, V& ypp){} //do nothing 
  };


  /**
   * \brief simple stepsize control for 
   */
  template<class V, char NORM = '2'>
  class SimpleStepSizeAdaption{
  public:
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef Norm<NORM,BaseType> NormType;
    

    SimpleStepSizeAdaption(const BaseType& tol, const V& y1, const V& y2):tol_(tol),y1_(y1),y2_(y2),tmp_(y1.size()){
      adonis_assert(y1.size() == y2.size());
    }
    
    template<class T>  
    value_type new_h(const T& hOld, int m, const BaseType& scal = 1. ){
      tmp_ = y1_-y2_;
      BaseType nm = NormType::norm(tmp_);
      adonis_assert(!is_zero(nm));
      return scal*hOld*NaturalPow(tol_/nm, -(m+1));
    }
    
  private:
    BaseType tol_;
    const V& y1_, y2_; //approx solutions: y1_ obtained by method of 
                             //order m, y2_ by method of order m+1
    V tmp_;
  };


  /**
   * \briefe Store smallest and largest stepsize $k_n$
   *
   * \tparam T stepsize type
   * \tparam B boolean. If 'true' a stepsize control is switched on, else
   *                    equidistant timestepping is performed
   */
  template<class T, bool B> class MinMaxStepsize;

  template<class T>
  class MinMaxStepsize<T,false>{  //equidistant stepsize, all kn are equal
  public:
    MinMaxStepsize(const T& kstart):kn_(kstart){}

    void compute(const T& kn){}  //do nothing 

    const T& min_stepsize() const {return kn_;}
    const T& max_stepsize() const {return kn_;}

  private:
    const T& kn_;
  };

  template<class T>
  class MinMaxStepsize<T,true>{ //stepsize controll switched on
  public:
    MinMaxStepsize(const T& kstart):minkn_(kstart),maxkn_(kstart){}

    void compute(const T& kn){
      minkn_ = Min(minkn_,kn);   //note that std::min, std::max don't work
      maxkn_ = Max(maxkn_,kn);   //with complex numbers. But an order relation
                                 //is z1 < z2 :<==> |z1| < |z2|
    }

    const T& min_stepsize() const {return minkn_;}
    const T& max_stepsize() const {return maxkn_;}

  private:
    T minkn_,
      maxkn_;
  };
  

} //end namespace Adonis

#endif
