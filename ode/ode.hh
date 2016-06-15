#ifndef ORDINARY_DIFFERENTIAL_EQUATION_INTERFACE_HH
#define ORDINARY_DIFFERENTIAL_EQUATION_INTERFACE_HH

#include <iomanip> 
#include <string>
#include <cmath>

#include "../expressiontemplates/exprvec.hh" 
#include "../linalg/linearsystemsolvers.hh"
#include "../misc/operations4randomaccesscontainers.hh" 
#include "tol.hh"
#include "constants.hh"
#include "diffequtilities.hh"
#include "../expressiontemplates/expressionvectortraits.hh"

#include "../dunestuff/fmatrix.hh"
#include "newton.hh"

#include "../misc/useful.hh"
#include "../graphics/printer.hh"
#include "../graphics/printmoldata.hh"

#include "../common/numerictypechecker.hh"
#include "../marcvecmatrix/myfunctions.h"

#include "stpszctrl.hh"

#include "../sparse/squarelinearsystem.hh"
#include "../sparse/sparseutilities.hh"
#include "../sparse/sparsitypattern.hh"

#include "../misc/xscaling.hh"
#include "../common/tempusfugit.hh"

#include "dampednewton.hh"
#include "simplifiednewton.hh"
#include "globalization4newton.hh"

//===== You can switch off/on AD by uncommenting /commenting the next line ====
//=== as long as you've provided a Jacobian matrix elsewhere ==================
#define APPLY_ALGORITHMIC_DIFF  USE_CPPAD
//=============================================================================


#include "../derivatives/derivativecalculator.hh"

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
#include "../accuracy/floatingpointarithmetic.hh"
#include "../expressiontemplates/xsettings.hh"
#endif 

#include "../sparse/sparsesettings.hh"
#include "../sparse/sparsematrix.hh"

#include "repair.hh"
#include "miscstuff.hh"
#include "convergencedetector.hh"

namespace Adonis{

  /** 
   *\brief General initial value formulation  
   * \f[ Y'(t) = F(t, Y(t))   \mathrm{on}\ [t0, tf] x R^n,  \\ Y(t0) = Y0. \f]
   *                         
   * \tparam K precision type
   * \tparam F functor for the RHS \f$ F(\cdot)\f$
   * \tparam PRINTMOL use special printing option when method of lines is 
   *  considered (default: 0, i.e. no method of lines printing style) 
   * \tparam SPARSELINALG use sparse linear algebra? true = yes, 
   *                                                 false = no (default)
   * \tparam SCALEITERATIONSYS scale iteration system of Newton iteration
   *    true = yes, false = no (default)
   * \tparam SPARSITYPATTERNVIAAD when sparsity pattern is computed, assume AD
   * pattern generation
   * \tparam REPAIR repair quantities inside solver
   * \tparam REPAIRTYPE select how quantities are to be corrected
   * \tparam PRINTPRESSURE print pressure, computed from differential quantities
   *                       via the ideal gas law ( true = print, false = don't
   *                       (default)
   */
  template<typename K, template<typename D> class F, int PRINTMOL = 0, bool SPARSELINALG = false, bool SCALEITERATIONSYS = false, bool SPARSITYPATTERNVIAAD = (SPARSELINALG == true) ? true : false, bool REPAIR = false, char REPAIRTYPE = 'n', bool PRINTPRESSURE = false>
  class ODE{
  
  public:
    enum{molPrintValue = PRINTMOL, isSparse = SPARSELINALG};

    typedef K value_type;
    typedef typename TypeAdapter<K>::BaseType BaseType;
    typedef typename ExprVecTraits<K>::ExprVecType ExprVecType;
    typedef typename ExprVecTraits<BaseType>::ExprVecType VecBaseType;
    
    typedef typename ExprVecTraits<int>::ExprVecType IndexVecType;
    typedef std::string StringType;
    typedef std::size_t SizeType;
    
    typedef F<K> FunType;  

    //matrix storage format of a dense matrix stored contiguous in memory
    typedef RowMajor DenseMatrixStorageOrganizationType;

  private:
    SizeType niter_;
    K t0_,                   //begin of time horizon
      tf_;                   //end of time horizon
    ExprVecType Y_;  //(inital) solution of (1)
    FunType fun_;                //the RHS - default constructor
    SizeType oom_;                 //order of method 
    
  public:
   
    //serves also as default constructor
    ODE(const K& t0 = K(), const K& tf = K(), const ExprVecType& YStart = ExprVecType()):niter_(0),t0_(t0), tf_(tf), Y_(YStart),fun_(YStart.size()),oom_(0){}
    
    
    //construct from any random access container providing iterators
    template<class V>
    ODE(const K& t0, const K& tf, const V& YStart):niter_(0),t0_(t0), tf_(tf),fun_(std::distance(YStart.begin(),YStart.end())),oom_(0){
      SizeType dim = std::distance(YStart.begin(),YStart.end());
      
      Y_.resize(dim);
      for(SizeType i = 0; i  < dim; ++i){
	Y_[i] = YStart[i];
      }
    }

    
    //construct from given iterators
    template<class ITER>
    ODE(const K& t0, const K& tf, ITER i1, ITER i2):niter_(0), t0_(t0), tf_(tf),fun_(std::distance(i1,i2)),oom_(0){
     
      SizeType dim = std::distance(i1,i2);
      Y_.reserve(dim);
      for(ITER it = i1; it  != i2; ++it){
	Y_.push_back(*it);
      }
    }

    
    //constructor for scalar system (1 equation only) 
    ODE(const K& t0, const K& tf, const K& start):niter_(0), t0_(t0), tf_(tf),fun_(1),oom_(0){
      Y_.resize(1); //one equation only 
      Y_[0] = start;
    }


    //copy-construction
    ODE(const ODE& ode):niter_(ode.niter_), t0_(ode.t0_), tf_(ode.tf_), Y_(ode.Y_), fun_(ode.fun_),oom_(0){}

    //copy-assignment
    ODE& operator=(const ODE& ode){
      if(this != &ode){
	niter_ = ode.niter_;
	t0_ = ode.t0_;
	tf_ = ode.tf_;
	Y_ = ode.Y_;
	fun_= ode.fun_;
	oom_ = ode.oom_;
      }
      return *this;
    }


    SizeType dimension() const {return Y_.size();}

    //access to fields
    K& tstart() {return t0_;}
    const K& tstart() const {return t0_;}
    K& tend() {return tf_;}
    const K& tend() const {return tf_;}

    void change_time_horizon(const K& t0, const K& tf){
      t0_ = t0;
      tf_ = tf;
    }

    template<class X>
    void change_initial_value(const X& x){
      for(SizeType i = 0; i < Y_.size(); ++ i)
	Y_[i] = x[i];
    }

    ExprVecType& state() {return Y_;}
    FunType& function()  {return fun_;}
    
    SizeType order_of_method() const {return oom_;}
    const SizeType number_of_iterations() const {return niter_;} 

   
    template<class C>
    ExprVecType& function(const C& x){
      return fun_(x);
    }

    SizeType dim() const{return Y_.size();}

    virtual ~ODE(){}

    virtual void print() const{
      std::cout <<"------------------------------------------------------"<<std::endl;
      std::cout <<" TIME HORIZON:  t_0 = "<<t0_<<"    t_f = "<<tf_<<std::endl<<std::endl;
      std::cout <<" INITIAL VALUE(s) X_0 = "<< Y_ << std::endl;
      std::cout <<" SYTEM DIMENSION (dim(X_0)) =  "<< Y_.size() <<std::endl<<std::endl;
      //std::cout <<"------------------------------"<<std::endl;
    }

    
    /**
   * \brief classical RUNGE-KUTTA 4 scheme just for testing the expression stuff
   *
   * NOTE: equidistant step size -- no 'virtual' can be used here since in ODE we use a function depending on one expression and for the derived ParamODE we use a function defined with two expressions!!
   */
    ExprVecType rk4(SizeType M, const StringType& fname = "rk4.dat"){
      oom_ = 4;

      if(M==0) //in case of only one iteration just return initial vector
	return Y_;

      K h = (tf_ - t0_)/M,
	t1 , t2; //time iterates
    
      SizeType dim = Y_.size();

      ExprVecType Y(Y_);  //start values
   
      //NOTE: Vec dim has to be known at compile time [cf. VELDHUIZEN, Expression Templates (in C++ Gems), 2003 ] 
      ExprVecType R1(dim), R2(dim), R3(dim); //the Runge-Kutta functions
      
     
      //PrintSolution PS(fname);

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif
      niter_ = M+1;
      for(SizeType i = 0; i < M+1; ++i){ //loop over time horizon

	//! floating point arithmetic for adding 2 numbers 
	//! each time step is considered as addition of 2 numbers, because the 
	//! result is used in consecutive operations
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	t1 = AdditionType::add(t0_,i*h);
	t2 = AdditionType::add(t1,0.5*h);
#else
	t1 = t0_ + i*h;
	t2 = t1 + 0.5*h;
#endif

	//std::cout << "i = "<< i << "   t = "<< time << "   h = "<< h <<std::endl;

	//! my expression vectors already posses the ability to perform
	//! compensated and balanced addition
	R1 = h*fun_(t1, Y);
	R2 = h*fun_(t2, Y + 0.5*R1);
	R3 = h*fun_(t2, Y + 0.5*R2);
      
	
	//overwrites Y_ with new iterate
	Y += ( ((R1 + 2.*(R2 + R3) + h*fun_(t1 + h, Y + R3)))/6. );
    
	//PS.write_2_file(t2,Y);
      }
    
      //PS.plot(2,4);

      return Y;
    }

    /*
    //! a 3-stage Runge-Kutta method -- Simpson
    ExprVecType rk3(SizeType M, const StringType& fname = "rk3.dat"){
      adonis_assert(M != 0);
      
      oom_ = 3;
      
      if(M==0) //in case of only one iteration just return initial vector
	return Y_;

       K h = (tf_ - t0_)/M,
	 time; //time iterate

       ExprVecType y(Y_);   //!this is \f$ y(t_0) = y_0 \f$
       
       //PrintSolution PS(fname);

       SizeType dim = Y_.size();
       ExprVecType K1(dim), K2(dim), K3(dim);


#if USE_TUNED_FLOATING_POINT_ARITHMETIC
 typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif

       for(SizeType i = 0; i < M+1; ++i){ //loop over time grid
	 //! current time point
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	 time = AdditionType::add(t0_,i*h);
#else	 

	 time = t0_ + i*h; 
#endif

	 //std::cout << "i = "<< i << "   t = "<< time << "   h = "<< h <<std::endl;

	 //intermediate stages
	 K1 = fun_(time,y);
	 K2 = fun_(
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
		   AdditionType::add(time,0.5*h),		   
#else
		   time + 0.5*h,
#endif		   
		   
		   y + 0.5*K1);
	 
	 K3 = fun_(
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
		   AdditionType::add(time,h),
#else		   
		   time + h,
#endif
		   y - h*K1 + 2.*h*K2);

	 y += h/6.*(K1 + 4.*K2 + K3);

	 //PS.write_2_file(t2,Y);

       }

       //PS.plot(2,4);
       return y;
       
    }
*/

    /**
     * \brief Adaptive RUNGE-KUTTA-FEHLBERG method with reduction extension
     *
     * This is an <I> explicit </I> method for solving <B>non-stiff</B> ODE systems of order (4)5 which uses 6 function evaluations
     * Ref.:
     * [1] [KINCAID, CHENEY, "Numerical Analysis", § 8, pp. 544]
     */
    ExprVecType adaptive_rkf(unsigned M, const K& hEstim, const K& tol, const StringType& fname = "adaptiveRungeKuttaFehlberg.dat", const K& hmin = 1e-12){
      
      NumericDataTypeChecker<K>::certify(); //guarantee numeric data types 

      adonis_assert(M != 0);
      
      oom_ = 5;

      if(M==0) //in case of only one iteration just return initial vector
	return Y_;


      unsigned dim =  Y_.size();

      ExprVecType x(Y_),  //assign from full
	y(dim),
	F1(dim), F2(dim), F3(dim), F4(dim), F5(dim), F6(dim),
	e(dim);  //local truncation error
       
      K h = (hEstim <= K()) ? (tf_ - t0_)/M : hEstim,
	t = t0_,
	d = K(),
	s = K();

      bool iflag = true;
      
      unsigned k = 0;
      
      
      typedef Norm<Constant<K>::whatNorm,typename TypeAdapter<K>::BaseType> NormType;

      typename TypeAdapter<K>::BaseType errnorm;

      unsigned count = 1;
      
      PrintSolution PS(fname);

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif 

      while( k < M ){  //works as well 
	k++;
	d = tf_ - t;

	if(Abs(d) <= Abs(h)){
	  iflag = false;
	  h = d;
	}
	
	s = t;
	y = x;   
	 
	//!compute \f$ F_i, \  i = 1, \ldots, 6 \f$ via function expressions
	F1 = h*fun_(t, x); 
#if USE_TUNED_FLOATING_POINT_ARITHMETIC

	F2 = h*fun_(AdditionType::add(t,1./4*h), x + 1./4*F1);
	F3 = h*fun_(AdditionType::add(t,3./8*h), x + 3./32*F1 + 9./32*F2);
	F4 = h*fun_(AdditionType::add(t,12./13*h), x + 1932./2197*F1 - 7200./2197*F2 + 7296./2197*F3);
	F5 = h*fun_(AdditionType::add(t,h), x + 439./216*F1 -8.*F2 + 3680./513*F3 - 845./4104*F4);
	F6 = h*fun_(AdditionType::add(t,0.5*h), x + -8./27*F1 + 2.*F2 - 3544./2565*F3 + 1859./4104*F4 - 11./40*F5);

#else

	F2 = h*fun_(t + 1./4*h, x + 1./4*F1);
	F3 = h*fun_(t + 3./8*h, x + 3./32*F1 + 9./32*F2);
	F4 = h*fun_(t + 12./13*h, x + 1932./2197*F1 - 7200./2197*F2 + 7296./2197*F3);
	F5 = h*fun_(t + h, x + 439./216*F1 -8.*F2 + 3680./513*F3 - 845./4104*F4);
	F6 = h*fun_(t + 0.5*h, x + -8./27*F1 + 2.*F2 - 3544./2565*F3 + 1859./4104*F4 - 11./40*F5);
#endif	

	//!compute \f$ x(t+h) = x(t) + \sum_{i=1}^6 a_iF_i \f$ -- 5th order
	x += (16./135*F1 + 0.*F2 + 6656./12825*F3 + 28561./56430*F4 - 9./50*F5 + 2./55*F6);
	
	
	//!compute local truncation error \f$ \sum_{i=1}^6 (a_i - b_i)F_i\f$
	e = 1./360*F1 + 0.*F2 - 128./4275*F3 - 2197./75240*F4 + 1./50*F5 + 2./55*F6;

	
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	AdditionType::pluseq(t,h);
#else
	t += h; //time update
#endif
	
	//some output
	std::cout << "k = "<< k << "   t = "<< std::showpoint << t << "   h = "<< std::showpoint << h <<std::endl; 
	//std::cout << "x = "<<x <<std::endl;//std::cout<< "e = "<< e << std::endl;     

	PS.write_2_file(t,x);

	if(!iflag)
	  break;   //stop

	errnorm = NormType::norm(e);
	
     
	if(errnorm >= tol){
	  t = s;
	  h = h_new(e,h,oom_,tol); //h *= 0.5;  //!change to oom_-1 ?
	  x = y;
	  k--;    //really set back k? Is that really needed? 
	}
	else{
	  if(errnorm < tol/128.) 
	    h = h_new(e,h,oom_,tol);  //h *= 2;  //!change to oom_-1
	}
	
	count++;  //count overall iterations


      } //end loop
      

      if(count > M)
	ADONIS_INFO(Information, "Maximum number M = " << M << " reached by count = "<<count<<". \n   Possibly no convergence...");

      niter_ = count;
      std::cout << "Solution found at t_final = "<< std::noshowpoint << t <<" after "<< count << " iterations."<<std::endl;
      
      PS.plot(1,2,//4,
	      "Adaptive RKF -- Order of method: " + my_function_collection::num2str(oom_)
	      //,"set terminal postscript eps color enhanced; set output \"" + PS.name()+".eps\"; " 
	      //,"set xrange[-2.1 : 1]; "
	      //, "set yrange[-14.5 : 14.5]; "
	      ); 
      
      return x;

    }


    /**
     * \brief slightly different version of adaptive RKF45
     *
     * This version is based on the HDNUM library of P. Bastian, cf. [1]
     *
     * [1] <a href="http://conan.iwr.uni-heidelberg.de/teaching/numerik1_ss2010/"> Peter Bastians Heidelberger Numerikbibliothek </a>
     */
    ExprVecType rkf45(unsigned M, const K& hEstim, const K& tol, const StringType& fname = "rkf45.dat", const K& hmin = 1e-12, const K& rho = 0.8, const K& alpha = 0.25, const K& beta = 4.){
       NumericDataTypeChecker<K>::certify(); //guarantee numeric data types 

      adonis_assert(M != 0);
      
      oom_ = 5;

      if(M==0) //in case of only one iteration just return initial vector
	return Y_;


      unsigned dim =  Y_.size();

      ExprVecType x(Y_),  //assign from full
	w(dim), ww(dim),
	F1(dim), F2(dim), F3(dim), F4(dim), F5(dim), F6(dim);
       
      K h = (hEstim <= K()) ? (tf_ - t0_)/M : hEstim,
	t = t0_,
	d = K();

      bool iflag = true,
	printflag = true;
      
      unsigned k = 0, rejected(0),count(0);
      
      
      typedef Norm<Constant<K>::whatNorm,typename TypeAdapter<K>::BaseType> NormType;

      typename TypeAdapter<K>::BaseType errnorm;
      
      PrintSolution PS(fname);

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif 

      SizeType c = 1;

      K c2 = 1./4, c3 = 3./8, c4 = 12./13, c5 = 1., c6 = 1./2.,
	a21 = 1./4., a31 = 3./32, a32 = 9./32,
	a41 = 1932./2197, a42 = -7200./2197, a43 = 7296./2197,
	a51 = 439./216, a52 = -8., a53 = 3680./513, a54 = -845./4104,
	a61 = -8./27, a62 = 2., a63 = -3544./2565, a64 = 1859./4104,
	a65 = -11./40.,
	b1 = 25./216, b2 = 0., b3 = 1408./2565, b4 = 2197./4104, b5 = -1./5,
	bb1 = 16./135, bb2 = 0., bb3 = 6656./12825, bb4 = 28561./56430, 
	bb5 = -9./50., bb6 = 2./55;

      for(k = 0; k < M; k+=c){
	count++;
	d = tf_ - t;

	if(Abs(d) <= Abs(h)){
	  iflag = false;
	  h = d;
	}

	//! COMPUTE stages
	//!compute \f$ F_i, \  i = 1, \ldots, 6 \f$ via function expressions
	F1 = h*fun_(t, x); 
#if USE_TUNED_FLOATING_POINT_ARITHMETIC

	F2 = h*fun_(AdditionType::add(t,c2*h), x + a21*F1);
	F3 = h*fun_(AdditionType::add(t,c3*h), x + a31*F1 + a32*F2);
	F4 = h*fun_(AdditionType::add(t,c4*h), x + a41*F1 + a42*F2 + a43*F3);
	F5 = h*fun_(AdditionType::add(t,c5*h), x + a51*F1 + a52*F2 + a53*F3 + a54*F4);
	F6 = h*fun_(AdditionType::add(t,c6*h), x + a61*F1 + a62*F2 + a63*F3 + a64*F4 + a65*F5);

#else

	F2 = h*fun_(t + c2*h, x + a21*F1);
	F3 = h*fun_(t + c3*h, x + a31*F1 + a32*F2);
	F4 = h*fun_(t + c4*h, x + a41*F1 + a42*F2 + a43*F3);
	F5 = h*fun_(t + c5*h, x + a51*F1 + a52*F2 + a53*F3 + a54*F4);
	F6 = h*fun_(t + c6*h, x + a61*F1 + a62*F2 + a63*F3 + a64*F4 + a65*F5);
#endif	

	//!compute 4th order approximation
	w = x;
	w += (b1*F1 + b2*F2 + b3*F3 + b4*F4 + b5*F5);

	//!compute 5th order approximation
	ww = x;
	ww += (bb1*F1 + bb2*F2 + bb3*F3 + bb4*F4 + bb5*F5 + bb6*F6);

	//! local error estimate
	w -= ww;

	//some output
	std::cout << "k = "<< k << "   t = "<< std::showpoint << t << "   h = "<< std::showpoint << h <<std::endl; 

	if(printflag) //! only print to file when we iterate
	  PS.write_2_file(t,x);

	if(!iflag)
	  break;   //stop


	errnorm = NormType::norm(w);
	
	BaseType h_opt(h*std::pow(rho*tol/errnorm,0.2));
	h_opt = std::min(beta*h,std::max(alpha*h,h_opt));

	if(errnorm <= tol){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	  AdditionType::pluseq(t,h);
#else
	  t += h; //time update
#endif
	  x = ww;
	  h = h_opt;
	  c = 1;
	  printflag = true;
	}
	else{
	  rejected++;
	  c = 0;         //don't increment counter in that case
	  h = h_opt;
	  printflag = false;
	  if(h > hmin) continue;
	}

      }//end for
        
       if(k > M)
	ADONIS_INFO(Information, "Maximum number M = " << M << " reached by count = "<< k <<". \n   Possibly no convergence...");

       niter_ = count;
       std::cout << "Solution found at t_final = "<< std::noshowpoint << t <<". Steps: "<< count << ".  Accepted: "<< k << "  Rejected: "<< rejected <<std::endl;
      
      PS.plot(1,2,//4,
	      "Adaptive RKF45 -- Order of method: " + my_function_collection::num2str(oom_)); 
      
      return x;
    }


    /**
     *\ brief the explicit euler
     */
    ExprVecType explicit_euler(SizeType M, const StringType& fname = "Explicit_Euler.dat", const IndexVecType& index = IndexVecType(), const StringType& additionalSettingsFile = StringType(),int printPrecision = 8){
      
      adonis_assert(M != 0);
      
      oom_ = 1;
      
      if(M==0) //in case of only one iteration just return initial vector
	 return Y_;

       K h = (tf_ - t0_)/M,
	 time; //time iterate

       ExprVecType y_current;

       y_current = Y_; 

       //std::ofstream of("ExplicitEuler.dat", std::ios_base::out);
       //PrintSolution PS(fname);
       PrintSolution PS;  
       PrintMOLData<PRINTMOL,BaseType,PRINTPRESSURE> DMP;
       DMP.init(additionalSettingsFile,fname,PS,printPrecision);

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
 typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif
 niter_ = M+1; 
 for(SizeType i = 0; i < M+1; ++i){ //loop over time grid

	 //! current time point
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	 time = AdditionType::add(t0_,i*h);
#else	 

	 time = t0_ + i*h; 
#endif



	 //normalize(y_current);   //may help !
	 
	 //output
	 std::cout << i << ".)   t = "<< time << "   dt = "<< h <<std::endl;
	 //std::cout << i << ".)    y = "<< y_current << std::endl; 
	 
	 //! Euler's method
	 y_current += h*fun_(y_current); 

       
	 /*
	 // modified Euler -- 2nd order (rk2) , o.k. with Babuska-Kahan 
	 y_current += h*fun_(
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
			     AdditionType::add(time,0.5*h),		     
#else
			     time + 0.5*h,
#endif
			     
			     y_current +  0.5*h*fun_(time,y_current));
	 */
	 

	 //PS.write_2_file(time, y_current);
	 DMP.write_2_file(i,time,y_current,fname,PS,index);
       } //end time-loop
       
       std::cout << "Solution found at t_final = "<< time <<" after "<< M+1 << " iterations."<<std::endl;

#if USE_PLOTTING_OPTION    
       PS.plot(1,2,"Explicit Euler -- Order of method: " + my_function_collection::num2str(oom_));
#endif       

       return y_current;
    }



/**
     *\ brief the Heun's method, a 2nd order 2-stage RK method
     */
    ExprVecType heun(SizeType M, const StringType& fname = "Explicit_Euler.dat"){
      
      adonis_assert(M != 0);
      
      oom_ = 2;
      
      if(M==0) //in case of only one iteration just return initial vector
	 return Y_;

       K h = (tf_ - t0_)/M,
	 time; //time iterate


       ExprVecType y_current;

       y_current = Y_; 

       ExprVecType K1(Y_.size()), K2(Y_.size());

       //std::ofstream of("ExplicitEuler.dat", std::ios_base::out);
       PrintSolution PS(fname);

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
 typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif
 niter_ = M+1;
       for(SizeType i = 0; i < M+1; ++i){ //loop over time grid

	 //! current time point
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	 time = AdditionType::add(t0_,i*h);
#else	 

	 time = t0_ + i*h; 
#endif

	 K1 = fun_(time,y_current);
	 K2 = fun_(
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
		   AdditionType::add(time,h),
#else
		   time+h,
#endif
		   
		   y_current + h*K1);


	 y_current += h*0.5*(K1 + K2);
	 
	 /*	 	 
	 //! Heun's method (RK2) -- 2nd order o.k. with Babuska-Kahan 
	 y_current += h*0.5*(fun_(time,y_current) + fun_(
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
							 AdditionType::add(time,h),							 
#else
							 time+h,
#endif							 
							 y_current + h*fun_(time,y_current)));
  	 */

	 PS.write_2_file(time, y_current);
       
       } //end time-loop
       
       std::cout << "Solution found at t_final = "<< time <<" after "<< M+1 << " iterations."<<std::endl;

       PS.plot(1,2,"Explicit Euler -- Order of method: " + my_function_collection::num2str(oom_));
       

       return y_current;
    }

    


    /**
   * \brief Predictor Corrector Euler
   *  Predictor: \f$ y_new \f$ is an estimate for the future value. Computed via  explicit Euler (i.e. method of same order). 
   *  Corrector: Evaluate rhs of ODE at \f$ y_new. \f$
   * Note: PC methods might be too inaccurate for stiff equations, see [LAMBERT, p. 238]
   */
    ExprVecType pc_euler(SizeType M, const StringType& fname = StringType()){
      oom_ = 1;

      if(M==0) //in case of only one iteration just return initial vector
	return Y_;
      
      K h = (tf_ - t0_)/M,
	tNext; //time iterate

      ExprVecType y_current, 
	y_new(Y_.size());

      y_current = Y_; 

      std::ofstream of("PCEuler.dat", std::ios_base::out);

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
 typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif

 niter_ = M+1;
      for(SizeType i = 0; i < M+1; ++i){ //loop over time grid
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	tNext = AdditionType::add(t0_,i*h);
#else	 

	tNext = t0_ + i*h; 
#endif

	//normalize(y_current);
	
	y_new = y_current + h*fun_(y_current);  //predict value via explicit Euler
	y_current += h*fun_(y_new);  //correct by y_new
      

	of << tNext << "    "<<y_current[0] << "    "<< y_current[1]<<std::endl;
      }
      of.close();

      return y_current;

    }
   
 

    /**
     * \brief Strategy for solving stiff ODEs with explicit methods with standard cG(1) method.
     *
     * The cG(1) (Galerkin's method with piecewise continuous linear functions) for an IVP of the form \f[ \dot u(t) = f(u(t))  t \in ]0,T], \ u(0) = u_0 \f] is derived in [2]. In short, on each subinterval \$ I_n = ]t_{n-1},t_n[ $\ of a partition $T_h$ of the time horizon $[0,T]$, we have \f[ \int_{t_{n-1}}^{t_n} (\dot U,v)dt =   \int_{t_{n-1}}^{t_n} (f,v)dt  \forall v \in P^{q-1}]t_{n_1},t_n[,\f] with \f$ P^{q-1}]t_{n_1},t_n[\f$ being the space of polynomials of degree at most \f$ q\f$ on \f$ I_n.\f$ This space is introduced since it produces more accurate approximations. Since \f$ q = 1 \f$ and testing against polynomials $t^j, \ j = 0, \ldots, q-1$ they vanish, i.e. for all test fcts we have \f$ v = 1.\f$
Therefore, the approximation simplifies to \f[ \int_{t_{n-1}}^{t_n} \dot U dt =   \int_{t_{n-1}}^{t_n} f dt.\f] Applying the fundamental theorem of calculus to the first term, we end up with \f[ U(t_n) - U(t_{n-1}) =  \int_{t_{n-1}}^{t_n} f dt.\f] The RHS-integral can be evaluated by, e.g. the Gauss midpoint rule, \f$ \int_a^b f dt \approx (b-a)f(\frac{a+b}{2})\f$ or, when function evaluations doesn't count much, by \f$ \approx (b-a)\frac{f(a)+f(b)}{2}\f$, which gives in the above notation (only for the midpoint rule): \f$ \int_{t_{n-1}}^{t_n} f dt \approx k_n(f(\frac{U_n + U_{n-1}}{2}), \f$ where $U_n = U(t_n), k_n = t_n - t_{n-1}.$ In this case we call it implicit cG(1) method with midpoint quadrature. This can be solved either by Newton's iteration or a fixed point iteration (as done in the following algo). 
     *
     * \param TOL measures error \f$ e(T) := U(T) -u(t), \f$ where \f$ U\f$ and \f$ u\f$ are the approximated and analytical solution respectively
     * \param tol gives the value at which convergence of the fixed point iteration is accepted (try TOL = tol)
     * \param hEstim an appropriate initial time step
     * \param stabilityFactor (e.g. <= 1)
     * \param c close to 1
     * \param maxFixedIters number of iterations after which divergence of fixed point iteration is taken
     * \param fname name of file writing results
     * \param kmin minimal time step that is allowed
     * \param kmax maximal time step that is allowed
     * \param index, additionalSettingsFile, printPrecision these are arguments
     *        needed for 1D, 2D MOL plotting
     * Reference:
     *
     *     [1] ERIKSSON, JOHNSON and LOGG, "Explicit time-stepping for stiff ODEs", SIAM, J. Sci. Comput., 25, 2003
     *
     *     [2] ERIKSSON, ESTEP, HANSBO and JOHNSON, "Computational Differential Equations", Studentlitteratur, 1996, p. 210, 211 (and 182,183 for introductory remarks)
     */
    ExprVecType explicit_4_stiff(const BaseType& TOL, const BaseType& tol, const K& hEstim, const BaseType& stabilityFactor, const BaseType& c = 0.9675, SizeType maxFixedIters = 13, const StringType& fname = "explicit4stiffDue2ErikssonJohnsonLogg.dat", const K& kmin = 0.0000125, const K& kmax = 0.75,const IndexVecType& index = IndexVecType(), const StringType& additionalSettingsFile = StringType(), int printPrecision = 8) {
      
      adonis_assert(maxFixedIters >= 2);

      oom_ = 1;

      SizeType dim = Y_.size();
      
      ExprVecType u(Y_),   //u0 
	u_prev(Y_),
	residual(dim),    //measures how well u satisfies the ODE
	discreteResidual(dim),
	ul(dim);  //needded for fixed point iteration

      K time = t0_,
	kn = hEstim,   //suitable initial time step
	ktilde;

      //K k;

      BaseType resnorm = BaseType(),  //!\f$r_{n,l}\f$
	resPrevNorm = BaseType();     //!\f$r_{n,l-1} \f$

      //ExprVecType ul_1(dim), err(dim);
     

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif


      
     bool performMsteps = false;

      //=================== NORM ===================
      typedef Norm<Constant<K>::whatNorm,typename TypeAdapter<K>::BaseType> NormType;

      //UserDefinedTolerance<NormType> canIbreakLoop(tol);
      //============================================



      //PrintSolution PS(fname); 
      PrintSolution PS;  
      PrintMOLData<PRINTMOL,BaseType,PRINTPRESSURE> DMP;
      DMP.init(additionalSettingsFile,fname,PS,printPrecision);


      BaseType Lscal;
      SizeType m = 1;


      
      //2.) TIME LOOP
      SizeType countTimeSteps = 1;
      
      while(time < tf_){
	
	//(a)
	if(countTimeSteps > 1){
	  //!"continuous" residual based on the solution u
	  residual = 1./kn*(u-u_prev) - fun_(u_prev);
	  

	  //(e) update 
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	  AdditionType::pluseq(time,m*kn); 
#else
	  time += m*kn;   //iterate time 
#endif	
	

	  //reset m
	  m = 1;

	  u_prev = u; //save old approximation
	  
	  //(b)
	  //std::cout << "||R(U(t_{n-1}))|| = "<<  NormType::norm(residual) <<std::endl;
	  ktilde = TOL/(stabilityFactor*NormType::norm(residual));
	  kn = 2.*ktilde*kn/(
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
			     AdditionType::add(ktilde,kn)
#else
			     ktilde + kn
#endif   

			     );

	}//end if countTimeSteps > 1

	if(kn > kmax)  //! if stepsize get's too big, stub it reasonably
	  kn = kmax;
	if(kn < kmin) //added this here
	  kn = kmin;

	//(c)
	ul = u_prev; //!start fixed point iteration with \f$U_{n-1}.\f$

	
        performMsteps = false; //reset
       

	//! solve discrete equations via FIXED POINT ITERATION on \f$ ]t_{n-1}, t_n[\f$
	//! This is a functional iteration and it is required \f$k_n\|\frac{\partial fun_}{\partial u}\| < 1\f$ for some norm in order to establish a contraction mapping, which in turn forces the step size \f$k_n\f$ to be restricted.
  	for(SizeType l = 1; l <= maxFixedIters; ++l){
	  //ul_1 = ul;
	  //! use midpoint quadrature to compute \f$\int f(u(t)dt \f$
	  ul = u_prev + kn*fun_(0.5*(u_prev + ul));
	  
	  discreteResidual = 1./kn*(ul - u_prev) - fun_(0.5*(u_prev + ul));
	 
	  // err = ul - ul_1;
	  

	  resnorm = NormType::norm(discreteResidual);  //::norm(err);  
	  
	  adonis_assert(is_well_defined_value(resnorm));
	 
	  //! discrete residual serves as stop criterion for fixed pt. iter.
	  if(resnorm <= tol){  
	    u = ul;  //UPDATE u
	    break;
	  }
	  
	  //(d) first step 
	  //! no convergence up to here
	  if(l == maxFixedIters-1){
	    resPrevNorm = resnorm;
	  }
	  
	  //!Damping: stabilize system with a couple of small euler steps
	  if(l == maxFixedIters){ //o.k. now we reached max. number of allowed
                                  //iterations and there is no convergence 
	    Lscal = 2.*resnorm/(kn*resPrevNorm);
	  
	    //! this might be more adequate for parabolic systems
	    //BaseType logkL = std::log(kn*Lscal)/std::log(2); 
	    BaseType logkL = std::log(kn*Lscal);

	    //std::cout << "LOG(kn·L) = "<< logkL << std::endl;

	    adonis_assert(is_well_defined_value(logkL));
	  
	      
	    (logkL < BaseType()) ? logkL = std::floor(logkL) : logkL = std::ceil(logkL);
	    
	    m = Abs(logkL);
	    
	    adonis_assert(m != 0);

	    kn = c/Lscal;            //! c~1 affects damping
	    //k = c/Lscal;           

	    performMsteps = true;
	    
	    if(kn < kmin)  //!guarantee that we don't go below kmin 
	      kn = kmin;
	   
	  }
	  
	}//end FIXED POINT ITERATION


	

	//(d) continued
	//! take m explicit Euler steps with time step k -- in case of fixed
	//! point iteration has failed to converge 
	if(performMsteps){
	  ADONIS_INFO(Information, "Fixed point iteration failed at the " << 	countTimeSteps <<"-th time step. \n   Perform m = "<< m << " explicit Euler steps with of size "<< kn);
	  for(SizeType i = 0; i < m; ++i){
	    u += kn*fun_(u);  
	    //u += k*fun_(u);
	  }
	  //(e)   //now kn is upated 
	  //#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	  //AdditionType::pluseq(time,m*kn); 
	  //#else
	  //time += m*kn;   //iterate time 
	  //#endif     
	}


	//!correct time towards final time T
	if(time >= tf_){
	  break;
	}
	else if (
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
		 AdditionType::add(time,kn)
		 
#else
		 
		 time+kn 
#endif     
		 > tf_)
	  kn = tf_-time;


	//! Print MOL variant
	PrintMeOrNot<Constant<BaseType>::write2FileAndPrint>::w_2_f(DMP,kn,countTimeSteps,time,u,fname,PS,index);

	std::cout << countTimeSteps << ".)    time = "<< time << "     kn = "<< kn << std::endl; 

	countTimeSteps++; //! count number of time steps

      }//end TIME loop
      
     
      

#if USE_PLOTTING_OPTION     
#if USE_PLOTTING_WITH_ALL_SPECIES
      PrintMeOrNot<Constant<BaseType>::write2FileAndPrint>::all_graphics(PS,dim,"Explicit Method for stiff ODEs by ERIKSSON, JOHNSON and LOGG");
#else 
       //! print a single species -- e.g. against time or another one
      PrintMeOrNot<Constant<BaseType>::write2FileAndPrint>::graphics(PS,Constant<BaseType>::xax,Constant<BaseType>::yax,"ErikssonJohnsonLogg-method -- 1st order");
#endif
#endif
      niter_ = countTimeSteps;
       std::cout << "Weird method: Solution found at t_end = "<<time<< " (tf_ = "<< tf_<< ") after "<< countTimeSteps << " time iterations."<<std::endl;


      return u;

    }


     /** 
      *  \brief <B>Implicit (backward) Euler method </B> with stepsize control (*) and a simple damped Newton strategy.
      * \param Mxsteps number of (uniform) steps over \f$ [t_0, t_{\textrm{f}} \f$. 
      * \param maxNewt precaution which guarantees that somehow the NEWTON iteration stops even if the local error criterion, e.g. \f$ \|y^{[\nu+1]}_{n+1} -y^{[\nu]}_{n+1}\| \leq ntol, \f$ isn't met.
      * \param ntol user-defined tolerance for the newton iteration. If 'ntol' is too big (e.g. 1e-2) it might be that more iterations are needed to achieve convergence (then the 'maxNewt' is beneficial to confine the number of Newton iterations (due to its local quadratic convergence no more than a couple of steps,say 3 to 7, should suffice).
      \param rtol user-defined tolerance for controlling the local error (needed for the stepsize determination)
      * \param fname you can specify a file to which data will be written
      *              according to the ODE printflag PRINTMOL
      * \param hEstim an estimated stepsize, if set to <= zero, then hEstim = (tf_ - t0_)/M is taken as starting stepsize
      * \param Cstszctrl constant stepsize control (\f$ C > 0\f$), default C = 1
      * \param hmin, hmax ranges for stepsize h (default values given)
      * \param index for printing, selects only species whose index is specified
      * \param additionalSettingsFile additional settings can be given. This is
      * \param printPrecision precision for printing to file 
      *  mainly needed in conjunction with special <b> printing </b> options 
      *  when the method of lines is considered.
      * \param trans specifies the form of the linear system ('N' no transpose, 'T' transpose, 'C' conjugate transpose)
      *
      *  NOTE: since Fortran is row based you should use 'T' in order to match with C style column-wise notation (set to 'T' by default). 
      *
      * (*) REMARK: The stepsize control is remarkably simple, namely \f[ h_{n+1} = \frac{rtol \cdot h_n}{C|Y_{n-1} - Y_{n-2}|},\f]
      * i.e. <I> no </I> higher order formulas or repeated solution (e.g. by taking two computations, one with halved and one with full stepsize) are used in order to merely calculate a stepsize. Instead, it only relies on previously computed values which, from the computational point of view, makes it very attractive!
      *
      * -- for the formula, see [1, p. 911 (mid)]; for a simple algo employing this formula, cf. [2, p. 287 (bottom)] 
      * 
      * NOTE: Upon uncommenting 'USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD' in 
      *       file 'ode/constants.hh', a 2nd order L-stable Rosenbrock function is used,
      *       cf. [4]
      *
      * References:
      *
      *[1] [C. JOHNSON, "ERRROR ESTIMATES AND ADAPTIVE TIME-STEP CONTROL FOR A CLASS OF ONE-STEP METHODS FOR STIFF ORDINARY DIFFERENTIAL EQUATIONS", SIAM, vol. 25, 1988]
      * 
      * [2] [C.JOHNSON, Y.Y. NIE, V. THOM´EE, "AN A POSTERIORI ERROR ESTIMATE AND ADAPTIVE TIMESTEP CONTROL FOR A BACKWARD EULER DISCRETIZATION OF A PARABOLIC PROBLEM", SIAM, vol. 27, 1990]
      *
      * [3] [U.M.ASCHER, L.R.PETZOLD, "Computer methods for ODEs and DAEs", SIAM, 1998]
      * [4] [VERWER, SPEE, BLOM, HUNDSDORFER], "A SECOND-ORDER ROSENBROCK METHOD APPLIED TO PHOTOCHEMICAL DISPERSION PROBLEMS", SISC, vol. 20, 1999]
      */
    ExprVecType implicit_euler(unsigned Mxsteps, unsigned maxNewt, const BaseType& ntol, const BaseType& rtol, const StringType& fname = "BackwardEuler.dat",
			       const K& hEstim = K(), 
			       const BaseType& Cstszctrl = 1., //stpz param
			       const BaseType& hmin = 1.25e-12,
			       const BaseType& hmax = 1.25e+03,
			       const IndexVecType& index = IndexVecType(),
			       const StringType& additionalSettingsFile = StringType(),
			       int printPrecision = 12, 
			       const typename ExprVecType::value_type* pMolarMasses = 0, int fromLastStepOn = 0, int excessSpecIndex = 0, const K& kthStep = K()){
      
      adonis_assert(hEstim >= hmin && hEstim <= hmax);
      adonis_assert(hEstim <= tf_);

      typedef RepairPhysicalQuantity<K,PRINTMOL,REPAIR,REPAIRTYPE> RepairType;
      RepairType Repairer;
      Repairer.init(additionalSettingsFile,Y_);

 
      //! use or don't use stepsizectrl, cf. [1,2]
      typedef StepSizeControl4OneStepMethod<ExprVecType,Constant<K>::locErrNorm,Constant<K>::useSTSZctrl> SzctrlType;
      //! choose between simple and advanced stepsize controls
      typedef ChooseStepsizeControl<SzctrlType,Constant<BaseType>::whatSTSZctrl> ChooseStpszCtrlType;

      typedef ConvergenceDetector<ExprVecType,Constant<BaseType>::additionConvergenceCheck,PRINTMOL> ConverType;
      //! record minimum and maximum stepsize (adaptive stepsize control
      MinMaxStepsize<K,Constant<K>::useSTSZctrl> SmallestLargestStepSize(hEstim);

      //! you can set the boundary for some specific reasons by applying the
      //! functor's 'set_boundary' member, assumed it is a MOL functor
      //! not necessary when ghost points are used for MOL
      typedef BoundaryFromFunctor<FunType,IsMOLFunctor<FunType>::Value> BdyFromFuncType;
      //! decide whether equidistant or given M is to be used
      unsigned maxit = SzctrlType::max_iters(Mxsteps,t0_,tf_,hEstim);
      
      bool isNewtonConverged = false;
      bool isIntegrationFinished(false);
      
      NumericDataTypeChecker<K>::certify(); //guarantee numeric data types 

      adonis_assert(maxit != 0);
  
#ifdef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
      oom_ = 2;
#else
      oom_ = Constant<SizeType>::accuracyOfMethod;    //order of implicit method
      std::cout<< "Method is of "<< oom_ <<".-order accuracy."<< std::endl;
#endif

      K h_fixed =  (hEstim <= K()) ? (tf_ - t0_)/maxit : hEstim,
	h = h_fixed,                   //uniform stepsize
	time = t0_;                          //time iterate

      unsigned i = 0;  //count timesteps

      BaseType storeMe(0);
      bool readyToPrint(true);
      if(fromLastStepOn > 0){
	std::cout << "ODE: Print from slide "<<fromLastStepOn<<". Repair FIRST; just in case..."<<std::endl;
	Repairer.repair(Y_,excessSpecIndex);
	if(additionalSettingsFile.size() != 0){
	  ParameterData Parm;
	  Parm.read_from_file(additionalSettingsFile);
	  h = Parm.get_datum<K>("t0_last_timestepsize");
	  storeMe = Parm.get_datum<BaseType>("elapsed_time_so_far");
	  //!this is the step right after the last printed slide
	  //i = ((unsigned)fromLastStepOn-1)*Parm.get_datum<unsigned>("nthstep")+1;
	  i = Parm.get_datum<unsigned>("number_lasttimestep"); //+1;
	  readyToPrint = false;
	}
      }
    

      K h_prev = h;                              //STPSZCTRL
      //K t_prev = t0_;                          //STPSZCTRL
      

      

      int dim = (int)Y_.size();

      bool isStationary(false);

      ExprVecType y_prev(dim),         //!CONSTANT during Newton iteration!!
	g(dim),                  //rhs of Newton system
	y_n(dim),                //current solution
      	y_prev_prev(dim),     //STPSZCTRL
	err(dim),             //STPSTCTRL -- local discretization error estimate
      
	fprev;              //for 2nd order trapezoidal  

      typedef NLEQ4ImplicitMethod<Constant<K>::accuracyOfMethod,FunType> MethodType;
      MethodType::resize(dim,fprev); //resize it or not (in case of trapezoidal)


      //maybe an extraction has to be performed here
      y_prev = Y_; //better choices might be selected [Ascher,Petzold]
      
      y_n = y_prev;          //   -"- dito , e.g. y_n = 2.5*y_prev


      //BdyFromFuncType::set_boundary(y_prev,fun_);
      //BdyFromFuncType::set_boundary(y_n,fun_);
      
      //! JACOBIAN of iteration matrix
      //! if sparse, use compressed sparse column scheme, else dense format
      const char MtxFormat = ((SPARSELINALG==true) ? CompressedSparseStorageIdentity<SparseSettings::TypeOfCompression>::Value : 'D');
      
      //! currently not implemented
      if((SparseSettings::TypeOfCompression == 'r') || (SparseSettings::TypeOfCompression == 'R'))
	ADONIS_ERROR(ImplementationError,"Currently only dense (default) or compressed sparse column ('c', 'C') format is allowed.");

      //! in case of a 2-D MOL functor, i.e. when PRINTMOL = 2 the previous 
      //! composition y_prev is stored as iterator in the functor
      //! must be invoked here also
      ReadOutCompositionFromODEScheme<PRINTMOL,typename ExprVecType::iterator>::get_previous(fun_,y_prev.begin(),h);

      typedef JacobianMatrixCalculator<SPARSELINALG,K,F,ExprTmpl::MyVec,MtxFormat,int> JacobianType;  //!last argument must be int here (coz of UMFPACK)
      JacobianType Gprime(true,fun_);   //true = calculate G'
      Gprime.initialize(dim,dim,Y_);    //proper initialization with Y_
     
   
      //! SCALING -- the get_pattern() just returns a double in case of density.
      //!            (likewise, the PatternType is just a fake double)
      typedef ScaleSquareLinearSystem<SCALEITERATIONSYS,SPARSELINALG,ExprVecType,typename JacobianType::PatternType> ScalingType;
      ScalingType ScaleSys(Gprime.get_pattern(),rtol,dim);

      //! that's the right hand side of the Newton system, i.e.
      //! G(y_n) = y_n - y_prev - h0·S(y_n),   if implicit Euler is used
      //! G(y_n) = y_n - y_prev - 0.5·h0·(S(y_n) + S(y_prev)), if the trapezoid.
      //!                                                      is used
#ifndef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
      typedef NonLinearEquationForImplicitONEStepMethods<Constant<K>::accuracyOfMethod,FunType,ExprVecType> NLEQType;
#endif  
   
      //======================= SET NORM HERE ==================================
      typedef Norm<Constant<K>::whatNorm,typename TypeAdapter<K>::BaseType> NormType;
      //typedef Norm<Constant<K>::locErrNorm,BaseType> LocErrNormType;
      
      UserDefinedTolerance<NormType> Tol(ntol);
   
      //========================================================================
  
      
      //P·P·P·P·P·P·P·P·P·P·P·P·P·P· PRINT ·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P
      PrintSolution PS;  
      PrintMOLData<PRINTMOL,BaseType,PRINTPRESSURE> PMD;
      PMD.init(additionalSettingsFile,fname,PS,printPrecision);
      //P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P
    
#ifdef ENABLE_GLOBALIZATION_4_NEWTION_S_METHOD
      Globalization4Newton<ExprVecType,RepairType,Constant<BaseType>::globalizeMe> Globalizer(Repairer,dim);
#else
      SimpleDampedNewtonMassFraction<ExprVecType,PRINTMOL,Constant<K>::dampMe> DampedNewton(1.,Constant<BaseType>::lambdaMin);
#endif
      SimplifiedNewton<JacobianType,Constant<BaseType>::useSimplifiedNewt> SNewt(Gprime,Constant<BaseType>::iterFrom);

#ifndef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD //only for Newton
      unsigned count = 0;                           //count Newton iterations
      bool leaveLoop = false; //check if Newton iterations can be stopped
#endif
      
      //****************** SOLVE SQUARE LIN. SYSTEM ************************
      typedef SquareLinearSystem<typename JacobianType::MatrixType,ExprVecType,MtxFormat,Constant<BaseType>::CalcPatternOnce,int> LinearSystemType;
      //                   A,                  b,nrhs //systemType
      LinearSystemType LSS(Gprime.get_matrix(),g,1,    0); 
      
      
#ifdef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
      ExprVecType ev2(dim), k2(dim); //intermediate results

      typedef typename JacobianType::MatrixType IterationMatrixType;
      IterationMatrixType Gp2(Gprime.get_matrix());
      //two systems to be solved. TODO: ideally with one LU decomposition only
      LinearSystemType LSS2(Gp2,k2,1,    0); //(Gp2,k2,1, 0)
      
      //scaling
      ExprVecType vals;
      //! only resize vals when scaling is desired and in column compressed form
      ((MtxFormat == 'C') || (MtxFormat == 'c')) ? vals.resize(Gprime.dimension_of_values()) : do_nothing();
      //!  ScaleSys can be reused since only a const reference to the pattern is stored (which stays the same)
      typedef SparseMatrixSelector<MtxFormat,K,int> MtxSelectType;
      typedef MatrixTypeDetector<MtxFormat,K> MtxDetectType;
#endif

      //********************************************************************

      
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif


      CPUTime<BaseType> stepCPUtm,
	IntermediateCPUtm;
      StringType toScreen;
      BaseType stepElapsedCpuTime(0), intermedElapsedTime(0),
	 //add time for each step to get totally elapsed CPU time
	addCPUTime(0);            
                          
      if(fromLastStepOn > 0){
	addCPUTime = storeMe;  //now addCPUTime starts from the last recorded elapsed time 
      }

      stepCPUtm.start(false); //print nothing to screen

      unsigned l(0); //meaningful when kthStep != 0

      //TIME ITERATION
      for( ; i < maxit; ++i){       //loop over time points
	stepCPUtm.reset(); //reset -- calculate time of one timestep
	stepCPUtm.start(false); //start measuring time, print nothing to screen

	//P·P·P·P·P·P·P·P·P·P·P·P·P·P· PRINT ·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P
	
	if(readyToPrint){
	  if((!is_zero(Abs(kthStep)))){
	    //also guarantee that final slide is printed
	    if((Abs(time+h) >= l*Abs(kthStep)) || (is_equal(time,tf_))){
	      if((!is_zero(time))){
		//now are close to the value specified in l*k while not
		//letting h become too small
		 
		if(!is_equal(Abs(time+h),l*Abs(kthStep))){
		  std::cout << "kthStep: Adjusting h" << std::endl;
		  h = Max(hmin,Sgn(hEstim)*(Abs(time+h) - Abs(l*kthStep)));
		}
	      }
	      l++; //increment l
	      PMD.write_2_file(h,i,time,y_n,fname,PS,index,pMolarMasses,addCPUTime,fromLastStepOn,kthStep);
	    }
	  }
	  else{ //usual printing
	    PMD.write_2_file(h,i,time,y_n,fname,PS,index,pMolarMasses,addCPUTime,fromLastStepOn);
	  }
	}
	  
	if(fromLastStepOn > 0){  //omit printing the very first slide since this
	  readyToPrint = true;  //will be nothing than the already printed
	}
	
	//P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P·P
			    
	
	std::cout << "Timestep "<< i << ".)    t = "<< time << "    h = "<< h << ((!is_zero(stepElapsedCpuTime)) ? "    CPU: "+time_in_human_readable_form(stepElapsedCpuTime) : "") <<std::endl; //only print when greater than approx. zero

	if(isIntegrationFinished)
	  break; //end time reached or exceeded -- leave loop

	
#ifndef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
	count = 0;       //reset counter for Newton iters
#endif	
	//! increment time step
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	AdditionType::pluseq(time,h); 
#else
	time += h;   //iterate time 
#endif	
	
#ifndef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
	leaveLoop = false;   //reset
	//isNewtonConverged = false; //reset
#endif
	//! in case of a 2-D MOL functor, i.e. when PRINTMOL = 2 the previous 
	//! composition y_prev is stored as iterator in the functor
	ReadOutCompositionFromODEScheme<PRINTMOL,typename ExprVecType::iterator>::get_previous(fun_,y_prev.begin(),h);

#ifdef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
	//1st linear system
	//possible clipping performed
	Repairer.repair(y_n,excessSpecIndex); //repair me or not (obey bounds on species etc.)
	Repairer.info(i,-1, "1st Rosenbrock sys: ");

	//BdyFromFuncType::set_boundary(y_n,fun_);
	  
	g = fun_(y_n);
	Gprime.evaluate_jacobian(y_n);
	Gprime.compute_G_prime(h,1.,1+1./sqrt(2.));
	Gp2 = Gprime.get_matrix();   //use G' for second lin. system 
	
	ScaleSys.scale(Gprime.value_pointer(),g.begin(),y_n,rtol,ntol);
	Gprime.reorder_values(); //needed for solving the sys when sparse
	LSS.solve(); //g (a.k.a. k1) is overwritten now with solution of 1st sys, so might be G'
	ScaleSys.retrieve_solution(g);


	//2nd linear system
	ev2 = y_n + h*g;       //1st order consistend at t = t_n+1 ==> cheap local error estimation ofr step size control [4, p. 1461]
	//possible clipping performed
	Repairer.repair(ev2,excessSpecIndex); //repair me or not (obey bounds on species etc.)
	Repairer.info(i,-1,"2nd Rosenbrock sys: ");

	//BdyFromFuncType::set_boundary(ev2,fun_);

	k2 = fun_(ev2) - 2.*g;
	ScaleSys.scale(MtxDetectType::value_ptr(Gp2),k2.begin(),ev2,rtol,ntol);
	//!reorder values
	MtxSelectType::copy(vals,MtxDetectType::value_ptr(Gp2));
	MtxSelectType::reorder_values(MtxDetectType::value_ptr(Gp2),vals, Gprime.get_mapper_idx());
       
	LSS2.solve();  // solves Gp2 · x = k2 and overwrites k2 with sol x
	ScaleSys.retrieve_solution(k2);

	//! update solution
	y_n += 1.5*h*g + 0.5*h*k2;

#else //NEWTON METHOD for 1-step methods used

	//!the nonlinear equation (only references are stored)
	//!note that y_prev remains fixed during Newton iteration
	//!note that fprev stays fixed during Newton when trapezoidal is used
	NLEQType nleq(g,y_prev,h,fun_,fprev,count); 
	
#ifndef ENABLE_GLOBALIZATION_4_NEWTION_S_METHOD
	DampedNewton.reset_lambda(); //! before each Newton iter, lambda reset
#endif
	y_n = y_prev;    //assign previous value to start Newton iteration
	while((leaveLoop == false)){ //!NEWTON ITERATION LOOP
	  
	  Repairer.repair(y_n,excessSpecIndex); //repair me or not (obey bounds on species etc.)
	  Repairer.info(i,count);
	  //BdyFromFuncType::set_boundary(y_n,fun_);
	  g = -nleq(y_n,isStationary);  //negative nonlin. eq., i.e. -g(y_n), 
	  //'count' is addressed
	  //2nd arg. checks for stationarity
#ifdef ENABLE_GLOBALIZATION_4_NEWTION_S_METHOD
	  Globalizer.get_f0(i,g);
#endif

	  if(isStationary)
	    break;
	   
	  toScreen.clear();
	  IntermediateCPUtm.reset();
	  IntermediateCPUtm.start(false);
	  
	  //! 1.) compute Jacobian matrix S'
	  //Gprime.evaluate_jacobian(y_n);
	  SNewt.evaluate_jacobian(y_n,count); //use 

	  intermedElapsedTime = IntermediateCPUtm.stop(false);

	  toScreen = "   Jac evaluation: "+ time_in_human_readable_form(intermedElapsedTime);
	  IntermediateCPUtm.reset();
	  IntermediateCPUtm.start(false);
	  //! 2.) compute G' = a0·I-b0·h·S'
	  //Gprime.compute_G_prime(h,MethodType::a0(),MethodType::b0());
	  SNewt.compute_G_prime(h,count,MethodType::a0(),MethodType::b0());

	  intermedElapsedTime = IntermediateCPUtm.stop(false);
	  toScreen += "   Set up G': "+ time_in_human_readable_form(intermedElapsedTime);

	  //! SCALING of Newton system. Note that even when the CSC format is
	  //! used, we expect the values to be in CSR format since the pattern
	  //! is also row-oriented.
	  ScaleSys.scale(Gprime.value_pointer(),g.begin(),y_n,rtol,ntol);

	  //! 3.) reorder values before solving the system; 
	  //!     This becomes only necessary when compressed  sparse 'C'olumn 
	  //!     is used, since the default sparse solver is UMFPACK.
	  IntermediateCPUtm.reset();
	  IntermediateCPUtm.start(false);
	  
	  Gprime.reorder_values(); 
	  intermedElapsedTime = IntermediateCPUtm.stop(false);
	  toScreen += "   Reorder (if necessary): "+ time_in_human_readable_form(intermedElapsedTime);

	  //****************** SOLVE LIN. SYSTEM ******************************
	  IntermediateCPUtm.reset();
	  IntermediateCPUtm.start(false);
	  LSS.solve();  //now g is overwritten with (scaled) solution

	  intermedElapsedTime = IntermediateCPUtm.stop(false);
	  toScreen += "   Lin sys: "+ time_in_human_readable_form(intermedElapsedTime);
	  //*******************************************************************

	  if(toScreen.size() != 0)
	    std::cout << count << ".) "<< toScreen <<std::endl;
	  
	  //SCALING
	  ScaleSys.retrieve_solution(g); //scale solution back (or not)

	  //======= UPDATE Newton iterate y_n towards 'Newton direction' g ====
#ifdef ENABLE_GLOBALIZATION_4_NEWTION_S_METHOD
	  Globalizer.armijo_line_search(i,y_n,g,nleq,excessSpecIndex,isStationary);  
	  //y_n += g;
#else
	  y_n += (DampedNewton.damping(Repairer.are_phys_quantities_violated())*g);
#endif
	  //=====================================================================
	
	  // Repairer.repair(y_n,excessSpecIndex); //repair me or not (obey bounds on species etc.)
	  // Repairer.info(i,count);

	  count++;  //increment Newton iteration counter

	  
	  //if ||delta||<= TOL then convergence and exit
	  leaveLoop = Tol.terminate_iteration(g); //cf. [ASCHER/PETZOLD,p.54]
	  if(leaveLoop == true) 
	    isNewtonConverged = true;

	  if(count > maxNewt-1){ //count starts from 0 to maxNewt-1
	    //ADONIS_ERROR(MaxIterationError,
	    ADONIS_INFO(Information, "Maximum number of Newton iterations ("<<maxNewt<<") exceeded (count = "<<count<<") at the \n   "<< i <<"-th time step. \n   Newton iteration does NOT converge. \n   You may increase 'maxNewt' and invoke routine again.");
	    break;
	  }

	}//end while-loop of NEWTON iteration

#endif //USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD	
	
	//repair output again after solution has been found, for instance
	Repairer.repair(y_n,excessSpecIndex); //repair me or not (obey bounds on species etc.)
	//! for correct OUTPUT only since in most MOL functors rhs_|delta_Omega
	//! is Zero, therefore the solution needs to be updated w.r.t.the bdy
	BdyFromFuncType::set_boundary(y_n,fun_);
	
	//! not relevant for Rosenbrock method
	ChooseStpszCtrlType::overwrite_y_prev_prev(y_prev,y_prev_prev); 

	if((isStationary) || (ConverType::is_converged(y_n,y_prev,PMD.num_of_quants(y_n),Constant<BaseType>::convergeEps)==true) ){
	  FancyMessages().nice_output("Stationary solution found (within prescribed tolerance).",33);
	  h = tf_ -time;
	  time = tf_;
	  break;  //leave time loop, since nothing changes any more
	}

	y_prev = y_n; //overwrite previous value with new val from Newton iter

	

	//=========================================================
	//============== stepsize control  ========================
	//=========================================================
	
	//! h will be computed anew and h_prev will be overwritten by h,
	//! depending on the stepsize control method
#ifdef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
	//!ev2 is a 1st order approximation at t_n, y_n is 2nd order and can be
	//!used for the stpszcontrol.
	//!2nd argument y_n is a dummy anyway, 3rd argument corresponds to  
	//!y_prev and 4th to y_prev_prev now
	ChooseStpszCtrlType::compute_new_step_size(err,y_n,y_n,ev2,h,h_prev,hmin,hmax,rtol,Cstszctrl,isNewtonConverged,oom_);
	SmallestLargestStepSize.compute(h);
#else
	ChooseStpszCtrlType::compute_new_step_size(err,y_n,y_prev,y_prev_prev,h,h_prev,hmin,hmax,rtol,Cstszctrl,isNewtonConverged,oom_);
	SmallestLargestStepSize.compute(h); //only record largest and smallest stepsize during the time integration for later output
#endif
	//! adjustment at the end of time horizon 
	if(time >= tf_){
	  //break;
	  
	  isIntegrationFinished = true;
	}
	//else
	if (
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
		 AdditionType::add(time,h)
		 
#else
		 
		 time+h 
#endif     
		 > tf_){
	  h_prev = h;  //won't be used anymore
	  h = ( (!is_zero(tf_-time)) ? (tf_-time) : h_prev ) ;
	  //std::cout << "h_final = "<< h << std::endl;
	  //time = tf_; //reset to tf_
	}
	//==========================================================
	//==========================================================
	//==========================================================

	stepElapsedCpuTime = stepCPUtm.stop(false); //do not print anything to screen
	addCPUTime += stepElapsedCpuTime; //sum up time

      }//end time-loop
      //of.close();

      if(!is_zero(kthStep)){
	if(l == 0)
	  ADONIS_WARNING(Warning, "kthStep != 0 but no slide printed.\n Check kthStep and print criteria in fct 'implicit_euler'");
      }
      
      if(i >= maxit)
      	ADONIS_INFO(Information, "Maximum number of time steps reached (M = "<<maxit<<") \n   This may give you not the real solution! \n   Increase M and/or play around with rtol, pal.");

      //! inform me about some general settings of the implicit method
#if USE_PLOTTING_OPTION      
#if USE_PLOTTING_WITH_ALL_SPECIES
      PrintMeOrNot<Constant<BaseType>::write2FileAndPrint>::all_graphics(PS,dim,"Implicit Euler");
#else
      PrintMeOrNot<Constant<BaseType>::write2FileAndPrint>::graphics(PS,Constant<BaseType>::xax,Constant<BaseType>::yax,"Implicit Euler -- Order of method: " + my_function_collection::num2str(oom_));
      //PS.plot(1,2,". Order of method: " + my_function_collection::num2str(oom_)); //2 = H2, H2O = 6
#endif
#endif
     
      niter_ = i;

      std::cout << "Solution found at t_end = "<<time<< " (tf_ = "<< tf_<< ") after "<< i << " iterations (out of "<<maxit<<"). LS solver category: ";
      ((SPARSELINALG)?FancyMessages().display_in_bold("Sparse."):FancyMessages().display("Dense."));
      std::cout << "("<<LSS.system_to_solve() <<")"<< std::endl;
      std::cout << "Jacobian Type: ";
      ((USE_CPPAD == 1)?FancyMessages().display_in_bold("AD."):FancyMessages().display_in_bold("FD."));

#ifndef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
      std::cout << std::endl<< "SCALED Newton system?: ";
#else
      std::cout << "Rosenbrock system scaled?: ";
#endif
      ((SCALEITERATIONSYS)?FancyMessages().nice_output("yes"):FancyMessages().display("no."));

      std::cout << std::endl << "Stepsize info:"<<std::endl<< " type of stepsize control: " << (Constant<BaseType>::useSTSZctrl == false ? "equidistant" : "adaptive") <<std::endl<< " h_min: "<<SmallestLargestStepSize.min_stepsize()
		<< ",  h_max: "<< SmallestLargestStepSize.max_stepsize() << std::endl;

#ifndef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
      std::cout << "METHOD in use: implicit "<< ((Constant<K>::accuracyOfMethod == 1) ? "EULER." : ((Constant<K>::accuracyOfMethod == 2) ? "TRAPEZOIDAL." : "unknown.")) << std::endl;
      std::cout << ((Constant<BaseType>::useSimplifiedNewt == true) ? "SIMPLIFIED NEWTON" : "Conventional Newton method") << " used." << std::endl;
#ifdef ENABLE_GLOBALIZATION_4_NEWTION_S_METHOD
      std::cout << "  GLOBAL CONVERGENCE enhanced via Armijo LS" << std::endl;
#endif
      if(Constant<K>::accuracyOfMethod == 2){
#ifdef EVALUATE_FUN_AT_YPREV_IN_FCT_AGAIN_4_TRAPEZOIDAL_METH
	std::cout << "Addendum: I used 'fun_(yprev)' explicitly; this might get computationally expensive for complicated fun_..." << std::endl;
#endif
      }
#else
      std::cout << "2nd order ROSENBROCK METHOD used" << std::endl;
#endif
      ConverType::info();  //some additional convergence check applied?
      if(PRINTMOL == 2){
#ifndef NONCONSERVATIVE_FORM
      std::cout << "Calculation in CONSERVATIVE form. Can this be applied to your problem? \n If not uncomment '#define NONCONSERVATIVE_FORM' in ../fdm/gensettings.hh."<<std::endl;
#endif
      }
      std::cout << "------------------------------------------------------------"<< std::endl;
      return y_n;
    }

  }; //end of ODE class


 
 


  //================ PARAMETRIC VERSION -- just an example ==================== 
  /**
   * \brief PARAMETER Version of ODE 
   */
  template<typename T, template<typename D> class RHS>
  class ParamODE: public ODE<T,RHS>{
  
    typedef typename ODE<T,RHS>::value_type value_type;
    typedef typename ODE<T,RHS>::ExprVecType ExprVecType;
    typedef std::size_t SizeType;

  private:
    ExprVecType P0_;                   //parameter vector in R^p

  public:
    //derived constructor -- also default version
    ParamODE(const T& t0 = T(), const T& tf = T(), const ExprVecType& YStart = ExprVecType(), const ExprVecType& P0 = ExprVecType()):ODE<T,RHS>(t0,tf,YStart),P0_(P0){} 

  
   
 //for scalar cases in parameters as well as initial state 
   ParamODE(const T& t0, const T& tf, const T& sx0, const T& sp0):ODE<T,RHS>(t0,tf,sx0){
     P0_.resize(1);
     P0_[0] = sp0;
   }
   

   //mixed version 1
   ParamODE(const T& t0, const T& tf, const ExprVecType& YStart, const T& sp0):ODE<T,RHS>(t0,tf,YStart){
     P0_.resize(1);
     P0_[0] = sp0;
   }
   
   //mixed version 2
   ParamODE(const T& t0, const T& tf, const T& sx0, const ExprVecType& P0):ODE<T,RHS>(t0,tf,sx0), P0_(P0){}



   //copy constructor
   ParamODE(const ParamODE& pode):ODE<T,RHS>(pode),P0_(pode.P0_){}

   //copy assignment
   ParamODE& operator=(const ParamODE& pode){
     if(this != &pode){
       ODE<T,RHS>::operator=(pode);
       P0_ = pode.P0_;
     }
     return *this;
   }


   void print() const{
     ODE<T,RHS>::print();
     std::cout << " PARAMETER INITAL VALUE(s) = "<<P0_<<std::endl;
     std::cout << " PARAMETER DIMENSION (dim(P_0)) = "<<P0_.size()<<std::endl<<std::endl;
   }

    //the classical RK4 for a paramterised ode -- I don't want to integrate with it; it just serves as an example to check my expression templates ;-) 
    ExprVecType rk4(SizeType M){
      typedef ODE<T,RHS> BS;
  
    
      if(M==1) 
	return BS::state();

      T h = (BS::tend() - BS::tstart())/M,
	t1, t2; 
    
      SizeType dim = BS::dim();
    
      ExprVecType x(BS::state()), //initial vector copied
	p(P0_);                   //parameter vector copied

      //needed for wrapping
      ExprVecType R1(dim), R2(dim), R3(dim); //the Runge-Kutta functions
    
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
#endif

     
      for(SizeType i = 0; i < M; ++i){ //loop over time horizon
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	t1 = AdditionType::add(BS::tstart(),i*h);
	t2 = AdditionType::add(t1,0.5*h);
#else
	t1 = BS::tstart() + i*h;
	t2 = t1 + 0.5*h;	
#endif

	//std::cout<< "TIME(S): "<< t1<<"  "<<t2<<std::endl; //if you wanna check

	R1 = h*BS::function()(t1, x, p);
	R2 = h*((BS::function())(t2, x + 0.5*R1, p));
	R3 = h*BS::function()(t2, x + 0.5*R2, p);
      
	//to check output uncomment the next line(s)
	//std::cout<<" R1 =" << R1 <<" R2 =" << R2 <<" R3 =" << R3 <<std::endl;

	//overwrites Y_ with new iterate
	x += ( ((R1 + 2.*(R2 + R3) + h*BS::function()(t1 + h, x + R3, p)))/6. );
	
      }//end loop
    
    

      return x;

    }
 };


} //end namespace 

#endif
