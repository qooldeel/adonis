#ifndef CONSTANTS_WHICH_MAY_BE_USED_BY_VARIOUS_INTERFACES_HH
#define CONSTANTS_WHICH_MAY_BE_USED_BY_VARIOUS_INTERFACES_HH



/**
 * \brief Defines for ODE solver
*/
//!switch 2-stage 2nd order L-stable Rosenbrock method on  
#define USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD

//!Enable globalization (note, can only be used with Newton's method, and when 
//! constant below is set to 'true'. When switched off (commented out)
//! a simple damped Newton is used if dampMe = 'true', otherwise the normal
//! Newton is used
//#define ENABLE_GLOBALIZATION_4_NEWTION_S_METHOD

namespace Adonis{

  /**
   * \brief Constants/Settings for Differential Equation stuff  
   *  Wraps constants of your choice -- define global constants here
   *
   * \tparam T precision of computation 
   *\code
        Constants<double,int> con;
	std::cout << con.theUltimateAnswer << std::endl;

	//also possible
	double blicky = Constants<double>::theUltimateAnswer;
	std::cout << "blicky = "<<blicky << std::endl;
   * \endcode 
   */
  template<class T>
  class Constant{
  public:
    typedef T value_type;
    typedef unsigned index_type;
    
    //!don't forget the 'static const', when defining some new constants, hunk.
    static const T theUltimateAnswer;
   
    //constants for method of lines 
    static const unsigned numberOfSpecies = 2;
    static const unsigned spacePoints = 51; //41; //11; //101; //41; Brusselator with diff           
    static const unsigned spPts = 101; 
    
    static const unsigned discretizationPointsInSpace = 100;
  
    //threshold value below which is considered zero
    static const T ApproximatelyZero;
    static const T VerySmall;
    

    //reserve some default consecutive space
    static const unsigned Reserve = 20;
 
    static const char whatNorm = '2';
    static const char locErrNorm = '2';
     
    static const bool useSTSZctrl = true;   //!true = yes, false = equidist
    
    //! 's' might be not so reliable
    static const char whatSTSZctrl = 'E';   //!'E' = Eriksson, 's' = simple, 'n' = normal

    //!Damped newton (currently for mass fractions and 2D Mol only)
    static const bool dampMe = true;  //true: some damped Newton
    //! only available when ENABLE_GLOBALIZATION_4_NEWTION_S_METHOD is uncommented
    static const bool globalizeMe = true; //true: use globalization
    static const unsigned int dampFac = 2;
    static const T lambdaMin;
    //!Simplified Newton method
    static const bool useSimplifiedNewt = false; //true: G'(y_0), false: G'(y_n)
    static const unsigned int iterFrom = 5; //! recalculate G' from n*iterFrom, n= 1, 2, ...

    //! check for additional convergence of method
    static const bool additionConvergenceCheck = true;
    static const T convergeEps;

    //! sparsity pattern only once computed by UMFPACK
    static const bool CalcPatternOnce = true; 

    //! set error and stop calculations due to assumed severity of numerical
    //! result (false = proceed computation with warning, true = error)
    static const bool setErrorAndStop = false;  

    //contraction security factor ( 0 < fac < 1 )
    static const T contractionSecurityFac;

    // if contraction isn't fulfilled reduce time stepsize within Newton iter
    static const T decreaseTimeStepSizeWithinNewton;

    static const unsigned renpopeparameter = 2; //possibilities: 1,2,3,4 (2 and 4 nonlinear)

    //! less-than-ideal solution for specifying a 1D right bdy ;)
    //! must be equal to specifications in "../modred/dat" !!
    static const int rightbdy1D = 51;

    //! ============ order of (implicit) method ===============================
    static const unsigned accuracyOfMethod = 1; //1 = impl. Eul,
                                                //2 = impl. trapezoid
    //=========================================================================
 

    //! print solution of an ode (when set to true, then 1 against 2 is plotted) 
    static const bool write2FileAndPrint = true;
    static const unsigned xax = 1;
    static const unsigned yax = 2;

    static const unsigned whatApproximation = 1; //1st or 2nd approx.

  };

  //! non-integral types must be defined outside of class!!!
  template<class T> const T Constant<T>::theUltimateAnswer = 42.;
  template<class T> const T Constant<T>::ApproximatelyZero = 1.e-13;
  template<class T> const T Constant<T>::VerySmall = 1.e-10;
  template<class T> const T Constant<T>::contractionSecurityFac = 0.3;
  template<class T> const T Constant<T>::decreaseTimeStepSizeWithinNewton = 0.25;
  //! damped Newton
  template<class T> const T Constant<T>::lambdaMin = 1.e-10;

  //convergence detector
  template<class T> const T Constant<T>::convergeEps = 1.e-08;
}

#endif
