#ifndef OPTIMIZATION_APPROACH_BY_LEBIEDZ_TEST_HH
#define OPTIMIZATION_APPROACH_BY_LEBIEDZ_TEST_HH

#include <iostream>
#include <stdio.h>

#if USE_IPOPT
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#else
#error "U need to switch on IPOPT flag, i.e. USE_IPOPT = 1"
#endif

//#if USE_CPPAD
#include <cppad/cppad.hpp>
// #else
// #error "U need to switch on CPPAD flag, i.e. USE_CPPAD = 1"
// #endif



#include "../common/typeadapter.hh"
#include "../common/globalfunctions.hh"
#include "../common/universalconstants.hh"
#include "../common/smartassign.hh"
#include "../common/error.hh"
#include "../derivatives/jacobianofsource.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../misc/misctmps.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../templatemetaprograms/unrollloop.hh"

#include "simstartingpoint.hh"

#include "../ode/examples/automaticallygeneratedsourceterms/h2c6.hh"

#include "../dunestuff/fmatrix.hh"
#include "../dunexternal/dunextensions.hh"
#include "indexmanipulation.hh"
#include "rvp.hh"

namespace Adonis{

  template<int N, int R, template<class D> class FUN, int ENCOD, char NORM>
  class NLPSIMCalculation: public Ipopt::TNLP{
  public:
    typedef Ipopt::Number Number;
    typedef Ipopt::Index Index;
    typedef Number value_type;

    typedef FUN<Number> FunType;
    typedef JacS<Number,FUN,ExprTmpl::MyVec> JacobianType;
    typedef Norm<NORM,Number> NormType;

    typedef ExprTmpl::MyVec<Number> VType;
    typedef CppAD::AD<Number> ADType;
    typedef CppAD::ADFun<Number> ADFunType;
    typedef ExprTmpl::MyVec<ADType> VecADType;

    typedef SIMStartingPoint<Number,ENCOD> SimStartType;

    enum{
      M = 2*R    //# of constraints
    };   

    
    template<class V>
    NLPSIMCalculation(const Number& objScal, const V& red, 
		      VType& xstar, ADFunType& adseqObj, FunType& fun, SimStartType& simstart, VType& eval, VType& obj, JacobianType& jac, VType& gradf, VType& hessf):objScal_(objScal),xstar_(xstar),adseqObj_(adseqObj),fun_(fun),simstart_(simstart),eval_(eval),obj_(obj),jac_(jac),gradf_(gradf),hessf_(hessf){
      red_.resize(R);
      UnrollLoop<0,R>::assign_rac2rac(red_,red);
    }

   
    ~NLPSIMCalculation(){}

    
     bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
		      Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style)
    {
      n = N;    //number of variables
      m = M;    //number of equality and inequality constraints(here: 3 eq)
      nnz_jac_g = M*N; //assume that jacobian of constraints is dense
      nnz_h_lag = N*(N+1)/2; // we assume the jacobian is dense, too 
                      //(store lower triang only, i.e. 6*7/2 elements)

      // use the C style indexing (0-based)
      index_style = TNLP::C_STYLE;
      
      return true;
    }

    //! returns the variable bounds
    bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
    {
      // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
      // If desired, we could assert to make sure they are what we think they are.
      assert(n == N);
      assert(m == M);
      
       
      for (Index i=0; i<N; ++i) {
	x_l[i] = 0.;          // lower bounds
	x_u[i] = UniversalConstants<Number>::aboutInfinity;  // no upper bounds exist
      }
    
      //constraints g are equality constraints, so we set the
      // upper and lower bound to the same value
      g_l[0] = g_u[0] = 2.;           //this is rhs b of C·x
      g_l[1] = g_u[1] = 1.;
      
      g_l[2] = g_u[2] = red_[0];      //this is rhs r of B^T·x
      g_l[3] = g_u[3] = red_[1];

      // std::cout << "MARC:"<<std::endl;
      // for(Index i = 0; i < m; ++i)
      // 	std::cout << "("<< g_l[i] << ", " << g_u[i] << ")  ";
      // std::cout << std::endl;

      return true;
    }


     bool get_starting_point(Index n, bool init_x, Number* x,
			    bool init_z, Number* z_L, Number* z_U,
			    Index m, bool init_lambda,
			    Number* lambda)
    {
      // Here, we assume we only have starting values for x, if you code
      // your own NLP, you can provide starting values for the dual variables
      // if you wish
      assert(init_x == true);
      assert(init_z == false);
      assert(init_lambda == false);
      
      
      for(Index i = 0; i < N; ++i)
	x[i] = simstart_.get_point()[i];
     
      return true;
    }

    //! evaluate objective function 
    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
       assert(n == N);
       UnrollLoop<0,N>::assign(eval_.begin(),&x[0]);

       obj_ = matrix_vector_product(jac_.jacobian(eval_),N,fun_(eval_));
       obj_value = objScal_*ntimes<2>(NormType::norm(obj_));
       
       return true;
    }


     //! gradient of objective
    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
    {
      assert(n == N);
      UnrollLoop<0,N>::assign(eval_.begin(),&x[0]);
      
      gradf_ = adseqObj_.Jacobian(eval_);
      assert(gradf_.size() == N);
      for(Index i = 0; i < n; ++i)
	grad_f[i] = gradf_[i];
      
      return true;
    }
    

    //! evaluate constraints
    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
      assert(n == N);
      assert(m == M);

      g[0] = 2.*x[0] + 2.*x[4] + x[1] + x[5];   //C·x
      g[1] = 2.*x[2] +    x[4] + x[3] + x[5];
      g[2] = x[0];                              //B^T·x
      g[3] = x[4];
      
      // std::cout << "MARC:"<<std::endl;
      // for(Index i = 0; i < m; ++i)
      // 	std::cout << g[i] << " ";
      // std::cout << std::endl;

      return true;
    }

     //! Jacobian of constraints in triplex/coordinate format
      /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
    bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
		    Number* values){
      if(values == NULL){
	// return the structure of the jacobian, row-wise storage
	// this particular jacobian is dense
	//! store the values of the jacobian row-wise
	iRow[0] = 0; iRow[1] = 0; iRow[2] = 0; iRow[3] = 0; iRow[4] = 0; iRow[5] = 0;

	iRow[6] = 1; iRow[7] = 1; iRow[8] = 1; iRow[9] = 1; iRow[10] = 1; iRow[11] = 1;
	
	iRow[12] = 2; iRow[13] = 2; iRow[14] = 2; iRow[15] = 2; iRow[16] = 2; iRow[17] = 2;

	iRow[18] = 3; iRow[19] = 3; iRow[20] = 3; iRow[21] = 3; iRow[22] = 3; iRow[23] = 3;
	
	jCol[0] = 0; jCol[1] = 1; jCol[2] = 2; jCol[3] = 3; jCol[4] = 4; jCol[5] = 5;
	jCol[6] = 0; jCol[7] = 1; jCol[8] = 2; jCol[9] = 3; jCol[10] = 4; jCol[11] = 5;
	jCol[12] = 0; jCol[13] = 1; jCol[14] = 2; jCol[15] = 3; jCol[16] = 4; jCol[17] = 5;
	jCol[18] = 0; jCol[19] = 1; jCol[20] = 2; jCol[21] = 3; jCol[22] = 4; jCol[23] = 5;
	
      }
      else{
	// return the values of the jacobian of the constraints -- always const
	values[0] = 2; //0,0
	values[1] = 1; //0,1
	values[2] = 0; //0,2
	values[3] = 0; //0,3
	values[4] = 2; //0,4
	values[5] = 1; //0,5
	
	values[6] = 0; //1,0
	values[7] = 0; //1,1
	values[8] = 2; //1,2
	values[9] = 1; //1,3
	values[10] = 1; //1,4
	values[11] = 1; //1,5

	values[12] = 1.; //2,0
	values[13] = 0.; //2,1
	values[14] = 0.; //2,2
	values[15] = 0.; //2,3
	values[16] = 0.; //2,4
	values[17] = 0.; //2,5

	values[18] = 0.; //3,0
	values[19] = 0.; //3,1
	values[20] = 0.; //3,2
	values[21] = 0.; //3,3
	values[22] = 1.; //3,4
	values[23] = 0.; //3,5

      }
      return true;
    }
    
     /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   *
   * \[ \sigma_f \nabla^2 f(x) + \sum_{i=1}^m \nabla^2 \lambda_i g_i(x)\]
   *
   *  NOTE: since the elementary mass balance is always linear, the Hessians
   *        of the constraints, \f$g_i\f$, are alwas zero, which is why we 
   *        need <B> not </B> to compute them when forming the Hessian of 
   *        the Lagrangian :) :)
   */
    bool eval_h(Index n, const Number* x, bool new_x,
		Number obj_factor, Index m, const Number* lambda,
		bool new_lambda, Index nele_hess, Index* iRow,
		Index* jCol, Number* values)
    {
      if (values == NULL) {
	// return the structure. This is a symmetric matrix, fill the lower left
	// triangle only.
	// the hessian for this problem is actually dense
	Index idx=0;
	for (Index row = 0; row < N; row++) {
	  for (Index col = 0; col <= row; col++) {
	    iRow[idx] = row;
	    jCol[idx] = col;
	    idx++;
	  }
	}
	assert(idx == nele_hess);
      }
      else{
	// return the values. This is a symmetric matrix, fill the lower left
	// triangle only
	// fill the objective portion
      
	//!NOTE: since the Hessian of the constraints are zero, these parts can
	//!      be skipped ;)

	UnrollLoop<0,N>::assign(eval_.begin(),&x[0]);

	hessf_ = adseqObj_.Hessian(eval_,0);

	//! fill objective portion, lower left
	Index count = 0;
	for(Index i = 0; i < N; ++i){
	  for(Index j = 0; j <= i; ++j){
	    values[count] = obj_factor*hessf_[i*N + j]; //access lower elements
	    count++;
	  }
	}

      }
      return true;
    }


     /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    void finalize_solution(Ipopt::SolverReturn status,
			   Index n, const Number* x, const Number* z_L, const Number* z_U,
			   Index m, const Number* g, const Number* lambda,
			   Number obj_value,
			   const Ipopt::IpoptData* ip_data,
			   Ipopt::IpoptCalculatedQuantities* ip_cq){
      // here is where we would store the solution to variables, or write to a file, etc
      // so we could use the solution.
      
      // For this example, we write the solution to the console
      adonis_assert(xstar_.size() == N);
      //std::cout << " SOLUTION (primal) x* = " << std::endl;
      for(Index i = 0; i < n; ++i){
	//	std::cout << x[i] << "  ";
	xstar_[i] = x[i]; 
      }
      //std::cout << std::endl;
      
     
      // printf("\n\nObjective value\n");
      // printf("f(x*) = %e\n", obj_value);
      
      // printf("\nFinal value of the constraints:\n");
      // for (Index i=0; i<m ;i++) {
      // 	printf("g(%d) = %e\n", i, g[i]);
      // }
    }


   
  private:
    //! ONLY REFERENCES here 
    const Number& objScal_;
    //const VType& red_;
    VType& xstar_;      //solution vector
    ADFunType& adseqObj_;
    FunType& fun_;
    SimStartType& simstart_;
    VType& eval_;
    VType& obj_;
    JacobianType& jac_;
    VType& gradf_;
    VType& hessf_;
   
    VType red_;

    //make copy constructio and copy assignment private
    NLPSIMCalculation(const NLPSIMCalculation&);
    NLPSIMCalculation& operator=(const NLPSIMCalculation&);
  };







  //===========================================================================
  //===========================================================================
  //===========================================================================
  //===========================================================================
  //===========================================================================

  template<int N, int R, template<class D> class FUN, int ENCOD, char NORM = '2'>
  class LebdiedzOptimizationApproach{
  public:
    typedef NLPSIMCalculation<N,R,FUN,ENCOD,NORM> ProblemType;
  
    typedef typename ProblemType::Number Number;
    typedef typename ProblemType::Index Index;
    typedef typename ProblemType::FunType FunType;
    typedef typename ProblemType::JacobianType JacobianType;
    typedef typename ProblemType::NormType NormType;

    typedef typename ProblemType::VType VType;
    typedef typename ProblemType::ADType ADType;
    typedef typename ProblemType::ADFunType ADFunType;
    typedef typename ProblemType::VecADType VecADType;

    typedef typename ProblemType::SimStartType SimStartType;

    typedef Ipopt::SmartPtr<Ipopt::TNLP> IPTNLPType;
    typedef Ipopt::SmartPtr<Ipopt::IpoptApplication> IPApplType;
    
    typedef Dune::FieldMatrix<Number,R,N> BTType;
    typedef Dune::FieldMatrix<Number,N,N-R> UType;

    LebdiedzOptimizationApproach(const Number objScal = 1.e-10):objScal_(objScal),xstar_(N),obj_(N),eval_(N),fun_(N){
       jac_.set(N,N);

      //! create sequence for objective function
      VecADType Xo(N), Yo(1);
      CppAD::Independent(Xo);
      VecADType mtxvec(N);
      JacS<ADType,FUN,ExprTmpl::MyVec> Jacobi;
      Jacobi.set(N,N);
      FUN<ADType> Func(N);
      mtxvec = matrix_vector_product(Jacobi.jacobian(Xo),N,Func(Xo));
      Yo = objScal*ntimes<2>(Norm<'2',ADType>::norm(mtxvec));
      adseqObj_.Dependent(Xo,Yo);
      adseqObj_.optimize();
    
      red_.resize(R);
     

      simstart_.initialize(N);
      v_.resize(N);
      for(int i = 0; i < N; ++i) //! by default, assign simstart to v_
	smart_assign(v_[i],simstart_.get_point()[i]);
    

      //create B^T and U
      ExprTmpl::MyVec<int> index(2);
      index <<= 0,4;
     
      IndexManipulator<N,R,ExprTmpl::MyVec> IM;
      IM.create_unrepresented_index(index);
      ReactionProgressVariables<Number,N,R> RPV(BT_,U_); //references 
      RPV.create_B_T(index);     //now B^T is filled 
      RPV.create_U(IM.get_index()); //now U is filled
    }

    const BTType& B_T() const {return BT_;}
    const UType& U() const {return U_;}

    template<class V, class W>
    VType& evaluate(const V& red, const W& full){
      UnrollLoop<0,N>::assign_rac2rac(v_,full);

      UnrollLoop<0,R>::assign_rac2rac(red_,red);

      nlp_ = new ProblemType(objScal_,red,xstar_,adseqObj_,fun_,simstart_,eval_,obj_,jac_,gradf_,hessf_);
      app_= new Ipopt::IpoptApplication();

      // // Change some options
      // // Note: The following choices are only examples, they might not be
      // //       suitable for your optimization problem.
      app_->Options()->SetNumericValue("tol", 1e-5); //def: 1e-08
      app_->Options()->SetStringValue("mu_strategy", "monotone");  //"monotone"
                                                              //"adaptive"
      //due to the guys here, the "adaptive" flag is a heuristic

      app_->Options()->SetStringValue("output_file", "ipopt.out");
      app_->Options()->SetIntegerValue("print_level", 0); //don't print anything(0)
                                                     //def.: 5

      ////Dominik's tip:
      app_->Options()->SetNumericValue("mu_init",1.e-06);
      app_->Options()->SetStringValue("nlp_scaling_method", "none");
      app_->Options()->SetNumericValue("constr_viol_tol",1.e-04); //def: 1e-04
       app_->Options()->SetNumericValue("acceptable_tol",1.e-04); //def: 1e-06
      //app_->Options()->SetStringValue("derivative_test", "second-order");

      //app_->Options()->SetStringValue("hessian_approximation","limited-memory");
      // app_->Options()->SetStringValue("limited_memory_update_type","bfgs");
      // app_->Options()->SetStringValue("limited_memory_aug_solver","sherman-morrison");
      // app_->Options()->SetIntegerValue("limited_memory_max_history",8);
      
      //! don't relax simple bounds l-eps <= x <= u + eps
      app_->Options()->SetNumericValue("bound_relax_factor",1.e-04); //def.: 1e-08
      
      //my settings -- similar to those Jochen has chosen
      // app_->Options()->SetIntegerValue("print_level", 3);
      // app_->Options()->SetNumericValue("point_perturbation_radius", 0.0);
      // app_->Options()->SetNumericValue("bound_relax_factor",        0.0);
      // app_->Options()->SetNumericValue("bound_push",                1.0e-99);
      // app_->Options()->SetNumericValue("bound_frac",                1.0e-99);
      // app_->Options()->SetNumericValue("slack_bound_push",          1.0e-99);
      // app_->Options()->SetNumericValue("slack_bound_frac",          1.0e-99);
      // app_->Options()->SetNumericValue("constr_mult_init_max",      1.0e+99);
      // app_->Options()->SetNumericValue("bound_mult_init_val",       1.0);
      // //app_->Options()->SetStringValue("bound_mult_init_method", "mu-based");
      // app_->Options()->SetNumericValue("s_max",                     100);
      // app_->Options()->SetNumericValue("tol",                       1e-06);
      // app_->Options()->SetNumericValue("dual_inf_tol",              1e-01);
      // app_->Options()->SetNumericValue("constr_viol_tol",           2e-6); //-8
      // app_->Options()->SetNumericValue("compl_inf_tol",             1e-1);
      
      // app_->Options()->SetNumericValue("acceptable_tol",            1e-03);
      // app_->Options()->SetIntegerValue("acceptable_iter",           5);
      // app_->Options()->SetNumericValue("acceptable_dual_inf_tol",   1e+5);
      // app_->Options()->SetNumericValue("acceptable_constr_viol_tol",1e-6);
      // app_->Options()->SetNumericValue("acceptable_compl_inf_tol",  1e+7);
      // app_->Options()->SetNumericValue("acceptable_obj_change_tol", 1e+5);
      // app_->Options()->SetIntegerValue("accept_after_max_steps",    10);
      


      //int numiter = app_->Statistics()->IterationCount(); //give back number of iters
 

      // // The following overwrites the default name (ipopt.opt) of the
      // // options file
      // // app_->Options()->SetStringValue("option_file_name", "hs071.opt");
      
      // // Intialize the IpoptApplication and process the options
      Ipopt::ApplicationReturnStatus status;
      status = app_->Initialize();
      if (status != Ipopt::Solve_Succeeded) {
	printf("\n\n*** Error during initialization!\n");
	std::cout<< "IP STATUS: "<< (int) status << std::endl;
      }

      // // Ask Ipopt to solve the problem
      status = app_->OptimizeTNLP(nlp_);

      if (status == Ipopt::Solve_Succeeded) {
	std::cout <<"*** The problem solved :)" <<  std::endl;
      }
      else {
	//printf("\n\n*** The problem FAILED!\n");
	ADONIS_ERROR(DerivedError,"**** The problem FAILED, pal!\n   IP status: "<<(int)status << ".");
      }

      return xstar_;
    }
    

    const VType& get_z() const {return v_;}
    VType& get_z() {return v_;}

  private:
    Number objScal_;
    VType xstar_;      //solution vector
    VType obj_, gradf_, hessf_, eval_;
    FunType fun_;
    JacobianType jac_;
    ADFunType adseqObj_;
    VType red_;
    SimStartType simstart_;
    
    VType v_; 
    
    IPTNLPType nlp_;
    IPApplType app_;

    BTType BT_;
    UType U_;
  };

} //end namespace 

#endif
