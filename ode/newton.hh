#ifndef NEWTON_METHOD_HH
#define NEWTON_METHOD_HH

#include <iostream>

#include "../common/typeselector.hh"
#include "../linalg/linearsystemsolvers.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../common/elementaryoperations.hh"
#include "../misc/misctmps.hh"
#include "../misc/useful.hh"


//============= TODO: (un)comment the follwing line ==================
//#define EVALUATE_FUN_AT_YPREV_IN_FCT_AGAIN_4_TRAPEZOIDAL_METH
//====================================================================

namespace Adonis{

  /**
   * \brief Non Linear Equation which arise when <B>implicit</B> methods are to be solved. Currently only two one-step methods are implemented but it is easy to overload function for general purpose Linear Multistep Methods (LMMs)
   * \tparam N Order of the method
   * \tparam F &quot original function &quot, i.e. the right hand side of the (autonomous) ode \f[ y'(t) = f(y(t))\f]
   */
  template<unsigned N, class F> class NLEQ4ImplicitMethod;
  

  //! nonlinear equation for implicit Euler ( 
  template<class F>
  class NLEQ4ImplicitMethod<1,F>{
  public:
    typedef typename F::value_type value_type;
    typedef NLEQ4ImplicitMethod<1,F> NLEQType;  //ThisType

    template<class V>
    static inline void resize(size_t n, V& fprev){} //do nothing

    //! this form is needed by Newton's method (the RHS)
    template<class Y, class X>
    static inline Y& negative_equation(Y& g, value_type& h, const X& xn, const X& xprev, F& fun, Y& fprev, unsigned& count){ //! the last argument is not needed here
      return (g = -1.*xn + xprev + h*fun(xn));
    }
  
    //! this form is needed by the natual criterion function for the damping strategy
    template<class Y, class X>
    static inline Y& equation(Y& g,  value_type& h, const X& xn, const X& xprev, F& fun, Y& fprev, unsigned& count){  //! the last argument is not needed here
      return (g = xn - xprev - h*fun(xn));
    }
  
    //! check also for stationarity of solution
    template<class Y, class X>
    static inline Y& equation(bool& isStationary, Y& g,  value_type& h, const X& xn, const X& xprev, F& fun, Y& fprev, unsigned& count){ 
      typedef typename TypeAdapter<typename Y::value_type>::BaseType BaseType;
      
      g = fun(xn);  //only one function evaluation needed

      if(is_zero(Norm<'2',BaseType>::norm(g)))
	isStationary = true;
      
      return (g = xn - xprev - h*g);
     }

    static inline value_type a0() {return 1.;}
    static inline value_type b0() {return 1.;}
  };
  
   
  /** 
   * Nonlinear equation for implicit trapezoidal method 
   * Note, despite being A-stable, the trapezoidal method does not have an 
   * overall-stability, because it damps rapidly decaying components only very
   * mildly. In other words the trapezoidal method is not L-stable, whereas
   * the implicit Euler is.
   *
   * Actually, this is a special kind of \f$ \theta \f$-method, i.e. \f[ u_{n+1} = \theta F(u_{n+1}) + (1 - \theta)F(u_n),\f]
   * where \f$ \theta = 0.5.\f$
   * Some practitioners use, e.g., \f$ \theta = 0.51\f$ to make the procedure
   * behave more like the implicit Euler
   */
  template<class F>
  class NLEQ4ImplicitMethod<2,F>{
  public:
    typedef typename F::value_type value_type;
    typedef NLEQ4ImplicitMethod<2,F> NLEQType;  //ThisType 

    template<class V>
    static inline void resize(size_t n, V& fprev){
      (fprev.size() == 0) ? fprev.resize(n) : do_nothing(); 
    }

    //! Recall that xprev is const during Newton iteration, so is fun_(xprev)
    template<class Y, class X>
    static inline Y& negative_equation(Y& g, value_type& h, const X& xn, const X& xprev, F& fun, Y& fprev, unsigned& count){
      g = fun(xn);
#ifndef EVALUATE_FUN_AT_YPREV_IN_FCT_AGAIN_4_TRAPEZOIDAL_METH
      if(count == 0)
	fprev = g;   //xn has been seeded with xprev as starting val for Newton
#endif
      //!                                   ... + fprev does not yield good res.
      g = -1.*xn + xprev + 0.5*h*( g + 
#ifdef EVALUATE_FUN_AT_YPREV_IN_FCT_AGAIN_4_TRAPEZOIDAL_METH
				   fun(xprev)
#else
				   fprev
#endif
				   ); 
      return g;
    }
    
    template<class Y, class X>
    static inline Y& equation(Y& g,  value_type& h, const X& xn, const X& xprev, F& fun, Y& fprev, unsigned& count){
      g = fun(xn);
#ifndef EVALUATE_FUN_AT_YPREV_IN_FCT_AGAIN_4_TRAPEZOIDAL_METH
      if(count == 0)
	fprev = g;
#endif
      //!                                   ... + fprev does not yield good res.
      g = xn - xprev - 0.5*h*( g + 
#ifdef EVALUATE_FUN_AT_YPREV_IN_FCT_AGAIN_4_TRAPEZOIDAL_METH
			       fun(xprev)
#else
			       fprev
#endif
);
      
      //std::cout << "fprev =  "<< std::setprecision(16)<< fprev << "fun(x) = "<< fun(xprev) << std::endl;
      
      return g;
    }

    //! check also for stationarity of solution
    template<class Y, class X>
    static inline Y& equation(bool& isStationary, Y& g,  value_type& h, const X& xn, const X& xprev, F& fun,  Y& fprev, unsigned& count){ 
      typedef typename TypeAdapter<typename Y::value_type>::BaseType BaseType;
      
      g = fun(xn);  //only one function evaluation needed
          
#ifndef EVALUATE_FUN_AT_YPREV_IN_FCT_AGAIN_4_TRAPEZOIDAL_METH
      if(count == 0)
	fprev = g;
#endif
      if(is_zero(Norm<'2',BaseType>::norm(g)))
	isStationary = true;
      
      return (g = xn - xprev - 0.5*h*(g + 
#ifdef EVALUATE_FUN_AT_YPREV_IN_FCT_AGAIN_4_TRAPEZOIDAL_METH
				      fun(xprev)
#else
				      fprev
#endif
));
     }

  
    static inline value_type a0() {return 1.;}
    static inline value_type b0() {return 0.5;}
  };
  

  /**
   * \brief Every <I> implicit </I> one step, lmm as well as  bdf method possesses an iteration matrix of the form \f[ a_0 I - b_0 h f_y(y_n),  \quad a_0, b_0 \not= 0 \f]
   * \tparam NLEQ type of nonlinear equation (only needed for correct assignment of parameters \f$ a_0\f$ and \f$ b_0.\f$ 
   * \tparam JV type of iteration matrix (Jacobian at input) stored as random access container
   */
  template<class NLEQ, class JV> 
  class DyNLEQ{
  public:
    typedef typename JV::value_type value_type;

    //! mind the reference to stepsize h !!
    template<class INT> 
    static inline void form_iteration_matrix(JV& iterationMatrix, value_type& h, const INT& cols, const value_type& a0, const value_type& b0){
      //! implicit euler: \f$ a_0 = b_0 = 1\f$, implicit trapezoidal: \f$ a_0 = 1, b_0 = 0.5 \f$ 
      iterationMatrix *= -(NLEQ::NLEQType::b0()*h);
      update_diagonal<AddBasicElements>(iterationMatrix,cols,NLEQ::NLEQType::a0());
    }

    //! if you've already devided by a0, then this function may be appropriate
    //! mind the reference to stepsize h !!
    template<class INT> //!JV jacobian stored as random access cont.
    static inline void form_iteration_matrix(JV& iterationMatrix, value_type& h, const INT& cols, const value_type& b0deva0){
      iterationMatrix *= -h;
      update_diagonal<AddBasicElements>(iterationMatrix,cols,b0deva0);
    }
  };


  /**
   * \brief This class wraps the above implementation to get something of the form \f$ F(x) = 0, \f$ i.e. we are able to use it in a more 
   * <I>mathematical</I> form.
   *
   * \tparam N accuracy of nonlinear equation
   * \tparam OFU original function entering the nonlinear equation
   * \tparam V the this stores \f$ F(x) \f$
   * \tparam X argument of F
   */
  template<unsigned N, class OFU, class V, class X = V>
  class NonLinearEquationForImplicitONEStepMethods{
  public:
    typedef typename X::value_type value_type;
    typedef NLEQ4ImplicitMethod<N,OFU> NLEQType;
    typedef V VType;
    typedef X XType; 

    enum{ORDER = N}; //!access via <TT> ::ORDER </TT>

    NonLinearEquationForImplicitONEStepMethods(V& g, const X& xprev, value_type& h, OFU& fun, X& fprev, unsigned& count):g_(g),xprev_(xprev),h_(h),fun_(fun),fprev_(fprev),count_(count){}

    V& operator()(const X& yn){
      return NLEQType::equation(g_,h_,yn,xprev_,fun_,fprev_,count_);
    }
  
    //! check for stationarity, i.e. \f$ fun_(xn) ~ 0\f$
    V& operator()(const X& yn, bool& isStationary){
      return NLEQType::equation(isStationary,g_,h_,yn,xprev_,fun_,fprev_,count_);
    }

  private:
    V& g_;                        //store references only!
    const X& xprev_;
    value_type& h_;
    OFU& fun_;
    X& fprev_;
    unsigned& count_;
  };


  
}//end namespace

#endif
