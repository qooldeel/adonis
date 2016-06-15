#ifndef TODO__WRITE_YOUR_SPECIFIC_SBNLPS_HEREAFTER_HH
#define TODO__WRITE_YOUR_SPECIFIC_SBNLPS_HEREAFTER_HH

#include <cmath>

#include "sbeqnlp.hh"
#include "../common/typeadapter.hh"
#if USE_CPPAD
#include "../derivatives/jacobianofsource.hh"
#else 
#include "../noderivativemethods/finitedifferences.hh"
#endif


#include "../misc/operations4randomaccesscontainers.hh"
#include "../misc/misctmps.hh"
#include "../templatemetaprograms/unrollloop.hh"
#include "functions4opttests.hh"
#include "../massactionkinetics/data/thermochemicaldata.hh"

#include "../ode/sourceterms.hh"  //H2C6 by hand

namespace Adonis{


  ////////////////////////////////////////////////////////////////////////////
  ////////// TODO: YOUR EXAMPLES CAN BE SPECIFIED RIGHT HERE ///////////////// 
  //                                                                        //
  //   Solve the problem                                      ////////////////
  //        min{ f(x) | c(x) = 0, h(x) <= 0, l <= x <= u}     //////////////// 
  //                                                                        //
  ////////////////////////////////////////////////////////////////////////////


  template<class T>
  class SBEQProblem1: public GeneralSimplyConstrainedProblem<T,SBEQProblem1<T> >{
  public:
    typedef GeneralSimplyConstrainedProblem<T,SBEQProblem1<T> > Basis;
    typedef typename Basis::VType VType;
    typedef T value_type;

    typedef H2Combustion6Spex<T> FunType;
#if USE_CPPAD 
   typedef JacS<T,H2Combustion6Spex,ExprTmpl::MyVec> JacobianType;
#else
    typedef FiniteDifference<FunType> FDType;
#endif

//typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef Norm<'2',T> NormType;

    /////////////
    typedef ThermoData4Mechanism<T,6> DataType; //here the rpv index is readily
    //available
    ////////////

    enum{
      nvar = 6,
      neq = 2  + DataType::rednspec ///////TODO: uncomment the last term  
    };


    SBEQProblem1(const T& scal = 1.e-10):objScal_(scal),obj_(nvar),fun_(nvar){
      Basis::n_ = nvar;
      Basis::m_ = neq;
      
      Basis::initialize();   //perform resizing once

  #if USE_CPPAD
      jac_.set(Basis::n_,Basis::n_);
#endif
    }

    VType& low() {
      for(size_t i = 0; i < Basis::n_; ++i)
	Basis::low_[i] = 0;   
      
      return Basis::low_;
    }
    

    VType& up() {
      for(size_t i = 0; i < Basis::n_; ++i)
	Basis::up_[i] = 1.e+30; //! de factor: \f$ +\infty, \f$ i.e.     
                                //! bounded from below only   
      return Basis::up_;
    }
  
    template<class E>
    VType& equality_constraints(const E& x){
      //std::cout << "Basis::red_.size() = "<< Basis::red_.size() << std::endl;
      adonis_assert(static_cast<int>(Basis::red_.size()) == DataType::rednspec);
      Basis::eq_[0] = 2.*x[0] + 2.*x[4] + x[1] + x[5] - 2;
      Basis::eq_[1] = 2.*x[2] +    x[4] + x[3] + x[5] - 1; 

      Basis::eq_[2] = x[DataType::rpv_index()[0]] - Basis::red_[0];
      Basis::eq_[3] = x[DataType::rpv_index()[1]] - Basis::red_[1];

      return Basis::eq_;
    }
    

    template<class E>
    T objective(const E& x) {  //return in original value_type
      //unit_ = fun_(x);
      //normalize(unit_);
#if USE_CPPAD
      obj_ = matrix_vector_product(jac_.jacobian(x),Basis::n_,fun_(x));
#else
      FDType fd(fun_);
      obj_ = matrix_vector_product(fd.jacobian(x),Basis::n_,fun_(x));
#endif      
//std::cout << "obj_ ="<< obj_ << std::endl;
      return objScal_*ntimes<2>(NormType::norm(obj_));
    }

  private:
    T objScal_;   //!constant scaling for objective, doesn't change anything
    VType obj_;
    FunType fun_;
#if USE_CPPAD
    JacobianType jac_;
#endif    
    //VType unit_;
  };
  


  /**
   * Example taken from [1, eq. (17.3), p. 499 and eq. (17.40), p. 516 resp.]:
   * 
   *    min x1 +x2
   *
   *    s.t.
   *
   *    x1*x1 + x2*x2 - 2 = 0
   *
   * The optimal point is \f$ x^* = [-1,-1]^T. \f$
   * References:
   *  [1]  [NOCEDAL,WRIGHT, "Numerical Optimization", 2nd ed., 2006]  
   */
  template<class T>
  class EasyNLPProblem: public GeneralSimplyConstrainedProblem<T,EasyNLPProblem<T> >{
  public:
    typedef GeneralSimplyConstrainedProblem<T,EasyNLPProblem<T> > Basis;
    typedef typename Basis::VType VType;
    typedef T value_type;

    enum{
      nvar = 2,
      neq = 1
    };


    EasyNLPProblem(){
      Basis::n_ = nvar;
      Basis::m_ = neq;
      
      Basis::initialize();   //perform resizing once
    }

    VType& low() {
      for(size_t i = 0; i < Basis::n_; ++i)
	Basis::low_[i] = -1.e+40;   //! de facto: \f$ -\infty, \f$ i.e.
                                    //! not bounded from below
      return Basis::low_;
    }
    

    VType& up() {
      for(size_t i = 0; i < Basis::n_; ++i)
	Basis::up_[i] = 1.e+40; //! de facto: \f$ +\infty, \f$ i.e.     
                                //! bounded from above   
      return Basis::up_;
    }
  
    template<class E>
    VType& equality_constraints(const E& x){
      Basis::eq_[0] = x[0]*x[0] + x[1]*x[1] -2;
      return Basis::eq_;
    }
    

    template<class E>
    T objective(const E& x) {  //return in original value_type
      return (x[0] + x[1]);
    }
  };
  

  /**
   * \brief linearly constrained convex problem due to [SNYMAN, "Practical Mathematical Optimization", Springer 2005, Prob. 5.3.8.1, p.202] 
   * Take \f$ \lambda_0 = 0\f$
   * Optimal solution \f$ x^* = [1,4]^T\f$
   */
  template<class T>
  class SnymanNLPProblem: public GeneralSimplyConstrainedProblem<T,SnymanNLPProblem<T> >{
  public:
    typedef GeneralSimplyConstrainedProblem<T,SnymanNLPProblem<T> > Basis;
    typedef typename Basis::VType VType;
    typedef T value_type;

    enum{
      nvar = 2,
      neq = 1
    };


    SnymanNLPProblem(){
      Basis::n_ = nvar;
      Basis::m_ = neq;
      
      Basis::initialize();   //perform resizing once
    }

    VType& low() {
      for(size_t i = 0; i < Basis::n_; ++i)
	Basis::low_[i] = -1.e+40;   //! de facto: \f$ -\infty, \f$ i.e.
                                    //! not bounded from below
      return Basis::low_;
    }
    

    VType& up() {
      for(size_t i = 0; i < Basis::n_; ++i)
	Basis::up_[i] = 1.e+40; //! de facto: \f$ +\infty, \f$ i.e.     
                                //! bounded from above   
      return Basis::up_;
    }
  
    template<class E>
    VType& equality_constraints(const E& x){
      Basis::eq_[0] = x[0] + x[1] -5;
      return Basis::eq_;
    }
    

    template<class E>
    T objective(const E& x) {  //return in original value_type
      return (6*x[0]*x[0] + 4*x[0]*x[1] + 3*x[1]*x[1]);
    }
  };
  

  /**
   * \brief This is just the Rosenbrock function in a more general setting ;)
   * Note that the Rosenbrock function has a <B> global minimum </B> at 
   * \f$ x^* = [1,1]^T. \f$
   */
  template<class T>
  class RosenNLPProblem: public GeneralSimplyConstrainedProblem<T,RosenNLPProblem<T> >{
  public:
    typedef GeneralSimplyConstrainedProblem<T,RosenNLPProblem<T> > Basis;
    typedef typename Basis::VType VType;
    typedef T value_type;

    enum{
      nvar = 2,
      neq = 0    //! NO eq. constraints!! 
    };

    RosenNLPProblem(){
      Basis::n_ = nvar;
      Basis::m_ = neq;
      
      Basis::initialize();   //perform resizing once
    }

     VType& low() {
      Basis::low_[0] = -1.e+40; 
      Basis::low_[1] = -1.e+40;   //!\f$-\inf, \f$, second var isn't bounded from below 

      return Basis::low_;
    }
    

    VType& up() {
      Basis::up_[0] = +1.e+40;
      Basis::up_[1] = +1.e+40;
      return Basis::up_;
    }
  
    //! no constraints at all 
    template<class E>
    VType& equality_constraints(const E& x){
      return Basis::eq_;
    }
    

    template<class E>
    T objective(const E& x) {  //return in original value_type
      return rosenbrock_function(x);
    }
  };


  /**
   * \brief Problem with convex function and linear constraints as well as 
   *  simple bounds. 
   *  Example has been taken from a talk(?) at MIT (2006):
   *
   *  \f[ \min f(x) = 2x_1^2 + x_1x_2 + 2x_2^2 - 6x_1 - 8x_2 +15 \f]
   * subject to 
   * \f[ -2x_1 + 2x_2 +1 = 0, \\ x_1 + 2x_2 -5 \leq 0 \\ x_1 \leq 1.75, x_2 \leq 2.\f]
   *  approximate solution: \f$ x^* \approx [2.7, 2.2]^T \f$
   */
  template<class T>
  class MITEasyNLPProblem: public GeneralSimplyConstrainedProblem<T,MITEasyNLPProblem<T> >{
  public:
    typedef GeneralSimplyConstrainedProblem<T,MITEasyNLPProblem<T> > Basis;
    typedef typename Basis::VType VType;
    typedef T value_type;

    enum{
      nvar = 3,   // 2 vars + 1 slack
      neq = 2    //! 1 eq + 1 ineq 
    };


    MITEasyNLPProblem(){
      Basis::n_ = nvar;
      Basis::m_ = neq;
      
      Basis::initialize();   //perform resizing once
    }

    VType& low() {
      Basis::low_[0] = -1.e+19; 
      Basis::low_[1] = -1.e+19;   //!\f$-\inf, \f$, second var isn't bounded from below 
      Basis::low_[2] = 0.;   //slack 

      return Basis::low_;
    }
    

    VType& up() {
      Basis::up_[0] = 1.75;
      Basis::up_[1] = 2.;
      Basis::up_[2] = 1.e+19; //slack
      return Basis::up_;
    }
  
    //! no constraints at all 
    template<class E>
    VType& equality_constraints(const E& x){
      Basis::eq_[0] = -2*x[0] + 2*x[1] + 1;  //equality constraint
      Basis::eq_[1] = x[0] + 2*x[1] - 5 + x[2];  //ineq transformed to eq. via
      return Basis::eq_;                         //slack variable x[2]
    }
    

    template<class E>
    T objective(const E& x) {  //return in original value_type
      return (2*x[0]*x[0] + x[0]*x[1] + 2*x[1]*x[1] - 6*x[0] - 6*x[1] + 15);
    }
  };
  

  /**
   *  \brief Ipopt problem (problem no. 71 of Hock-Schittkowsky test suite)
   *
   *  See p. 92 of the following link
   *
   *   <a href="http://www.ai7.uni-bayreuth.de/test_problem_coll.pdf">All examples of the Hock-Schittkowsky-Collection</a>
   *
   * x0 = [1,5,5,1]
   *
   * x* = [1.0, 4.74299963, 3.82114998, 1.37940829]
   */
  template<class T>
  class HockSchittkowsky71: public GeneralSimplyConstrainedProblem<T,HockSchittkowsky71<T> >{
  public:
    typedef GeneralSimplyConstrainedProblem<T,HockSchittkowsky71<T> > Basis;
    typedef typename Basis::VType VType;
    typedef T value_type;

    enum{
      nvar = 5, //4 variables + 1 slack
      neq = 2   // 1 eq. + 1 ineq.
    };


    HockSchittkowsky71(){
      Basis::n_ = nvar;
      Basis::m_ = neq;
      
      Basis::initialize();   //perform resizing once
    }

     VType& low() {
       for(int i = 0; i < 4; ++i)
	 Basis::low_[i] = 1.;
       Basis::low_[4] = 0.;   //slack
      return Basis::low_;
    }
    

    VType& up() {
      for(int i = 0; i < 4; ++i)
	 Basis::up_[i] = 5.;
      Basis::up_[4] = 1.e+19;   //slack
      return Basis::up_;
    }
  
    //! constraints -- the last entry of x is the slack
    //! only take sum and product over the first 4 species 
    template<class E>
    VType& equality_constraints(const E& x){
      Basis::eq_[0] = UnrollLoop<0,4>::template sum<2>(x) - 40.;
      Basis::eq_[1] = -1.*UnrollLoop<0,4>::product(x) + 25. + x[4]; //slack  
      // Basis::eq_[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 40.;
      // Basis::eq_[1] = -x[0]*x[1]*x[2]*x[3] + 25. + x[4];  //+ slack 
      return Basis::eq_;
    }
    

    template<class E>
    T objective(const E& x) {  //return in original value_type
      return ( x[0]*x[3]*(x[0] + x[1] + x[2]) + x[2] );
    }
  };


  /**
   *  \brief Purely nonlinear example taken from [NOCEDAL and WRIGHT, "Numerical Optimization", 2nd edition, exercise 18.3, p. 562].
   *
   * Starting point: \f$ x_0 = [-1.71, 1.59, 1.82, -0.763, -0.763]^T \f$
   * Solution: \f$ x^* = [-1.8, 1.7, 1.9, -0.8, -0.8]^T\f$
   */
  template<class T>
  class TestNLPProblem: public GeneralSimplyConstrainedProblem<T,TestNLPProblem<T> >{
  public:
    typedef GeneralSimplyConstrainedProblem<T,TestNLPProblem<T> > Basis;
    typedef typename Basis::VType VType;
    typedef T value_type;

    enum{
      nvar = 5,
      neq = 3
    };


    TestNLPProblem(){
      Basis::n_ = nvar;
      Basis::m_ = neq;
      
      Basis::initialize();   //perform resizing once
    }


    //! no variable is bounded from below or above
    VType& low() {
      Basis::low_[0] = -1.e+40;
      Basis::low_[1] = -1.e+40;   
      Basis::low_[2] = -1.e+40;
      return Basis::low_;
    }
    

    VType& up() {
      Basis::up_[0] = +1.e+40;
      Basis::up_[1] = +1.e+40;
      Basis::up_[2] = +1.e+40;
      return Basis::up_;
    }
  
    //! constraints
    template<class E>
    VType& equality_constraints(const E& x){
      Basis::eq_[0] = UnrollLoop<0,nvar>::template sum<2>(x) - 10.;
      Basis::eq_[1] = x[1]*x[2] - 5.*x[3]*x[4];
      Basis::eq_[2] = ntimes<3>(x[0]) + ntimes<3>(x[1]) + 1.;
      
      return Basis::eq_;
    }
    

    template<class E>
    T objective(const E& x) {  //return in original value_type
      return exp(UnrollLoop<0,nvar>::product(x)) - 0.5*ntimes<2>(ntimes<3>(x[0])+ ntimes<3>(x[1]) + 1.);
    }
  };




  /**
   * \brief The following example is taken from the LANCELOT_SIMPLE user documentation, cf. pp. 8--11 of 
   <a href="http://www.galahad.rl.ac.uk/doc/lancelot_simple.pdf">Galahad -- Lancelot_simple </a> 
   * Note that the inequality constraint can be made an equality constraint by
   * introducing a slack variable \f$ s \geq 0 \f$ via \f$ g(x) + s = 0.\f$
   * 
   * Starting point is given by \f$ x_0 = [-1.2, 1.]^T\f$
   *
   * Solution of LANCELOT_SIMPLE:
     8 iterations. Optimal value = 0.023314
     x* = [ 0.8475, 0.7175]^T
   */
  template<class T>
  class LancSimpleNLPProblem: public GeneralSimplyConstrainedProblem<T,LancSimpleNLPProblem<T> >{
  public:
    typedef GeneralSimplyConstrainedProblem<T,LancSimpleNLPProblem<T> > Basis;
    typedef typename Basis::VType VType;
    typedef T value_type;

    enum{
      nvar = 3,  //2 variable + 1 slack variable
      neq = 2    //1 ec + 1 ec (converted inequality constraint via slack)
    };


    LancSimpleNLPProblem():infinity_(1.e+20){
      Basis::n_ = nvar;
      Basis::m_ = neq;
      
      Basis::initialize();   //perform resizing once
    }

    VType& low() {  //!only \f$ x_1\f$ and \f$s = x_3\f$ are bounded from below 
      Basis::low_[0] = 0.; 
      Basis::low_[1] = -infinity_;  
      Basis::low_[2] = 0.;    //slack variable
      return Basis::low_;
    }
    

    VType& up() {
      Basis::up_[0] = infinity_;
      Basis::up_[1] = 3.;
      Basis::up_[2] = infinity_;  //slack variable
      return Basis::up_;
    }
  
    //! no constraints at all 
    template<class E>
    VType& equality_constraints(const E& x){  //! x[2] = slack variable
      Basis::eq_[0] = x[0] + 3.*x[1] - 3.;
      Basis::eq_[1] = x[0]*x[0] + x[1]*x[1] - 4. + x[2];
      return Basis::eq_;
    }
    

    template<class E>
    T objective(const E& x) {  //actually that's the Rosenbrock function
      return rosenbrock_function(x);
    }

  private:
    T infinity_;
  };
  


} //end namespace 

#endif
