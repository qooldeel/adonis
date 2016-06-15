#ifndef GLOBALLY_CONVERGENT_AUGMENTED_LAGRANGIAN_METHOD_HH
#define GLOBALLY_CONVERGENT_AUGMENTED_LAGRANGIAN_METHOD_HH

#include <iostream>
#include <algorithm>
#include <string>

#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../common/error.hh"
#include "../common/typeadapter.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../misc/misctmps.hh"
#include "../misc/useful.hh"
#include "../templatemetaprograms/unrollloop.hh"
#include "../io/readinparameters.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../statistics/probability.hh"

#include "steihaug.hh"
#include "trustregion.hh"
#include "derivauglag.hh"


#include "../regularization/guaranteepositivity.hh"


namespace Adonis{
    /**
   * \brief Implementation of Algorithm 1, p. 551 of Ref. [1]
   *
   * References:
   *
   *    [1] [CONN,GOULD, and TOINT, "A globally convergent augmented Lagrangian Method", SIAM J. Numer. Anal., Vol. 28, No. 2, pp. 545-572, 1991]
   */

  /**
   * \brief Functor of the  <I> augmented Lagrangian function </I> of a  simply bounded problem  \f[ min\{f(x) | c(x) = 0 \wedge l \leq x \leq u \}.\f] Remember that the simple constraints aren't considered in the augmented Lagrangian. The function is therefore given by the expression \f[L_{\mathrm{aug}}(x,\lambda,S,\mu) := f(x) + \sum_{i=1}^m \lambda_i c_i(x) +\frac{1}{2\mu}\sum_{i=1}^m s_{ii}c_i(x)^2.\f]
   */
  template<class T, template<class D> class PROBL>
  class AugmentedLagrangianFunction{
  public:
    typedef PROBL<T> ProbType;
    typedef typename ProbType::VType VType;
 
    AugmentedLagrangianFunction():ec_(ProbType::neq){}

    size_t dim() const {return 1;} //! range space will always be one
    size_t domain_dim() const {return static_cast<size_t>(ProbType::nvar);}

    ProbType& get_problem() {return prob_;}


    template<class X, class L, class DIAG>
    T operator()(const X& x, const L& lambda, const DIAG& diag, const T& mu){
      adonis_assert(mu > T()); //! positivity of \$ mu\f$
      ec_ = prob_.equality_constraints(x); //! 2 evaluations are coming...
      return ( prob_.objective(x) + dot_product<ProbType::neq>(lambda,ec_) + 1./(2.*mu)*UnrollLoop<0,ProbType::neq>::template sum<2>(diag,ec_) );
    }

  private:
    ProbType prob_;  //! create problem instance of the form
                     //! \f$ min\{f(x) | c(x) \wedge l \leq x \leq u \}\f$
    VType ec_;    //! equality constraints vector to store evaluation 
  };

 
  
  /**
   * \brief Projection onto box constraints \f$ l_i \leq x_i \leq u_i, \quad i = 1,\ldots, m\f$. More precisely, in this special case, the projection operator
   maps each component of \f$x\f$ as follows: \f$ l_i,  \ x_i \leq l_i, u_i, \ x_i \geq u_i, x_i \f$ otherwise. 
   */
  template<class X>
  class ProjectionOntoBoxConstraints{
  public:
    typedef X VType;
    typedef typename VType::value_type value_type;

    ProjectionOntoBoxConstraints(const X& l, const X& u):low_(l),up_(u),dim_(static_cast<size_t>(std::distance(l.begin(),l.end()))){
      adonis_assert(std::distance(l.begin(),l.end()) == std::distance(u.begin(),u.end()));
    } 

    const size_t size() const {return dim_;}

    //!componentwise projection 
    template<class E>
    value_type operator()(size_t i, const E& x){
      adonis_assert(low_[i] <= up_[i]);

      return ( (x[i] <= low_[i]) ? low_[i] : ( (x[i] >= up_[i]) ? up_[i] : x[i] ) );
    }

    template<class C1, class C2>
    C1& projection(C1& proj, const C2& x){
      for(size_t i = 0; i < dim_; ++i)
	proj[i] = (*this).operator()(i,x);
      return proj;
    }
      

    //! this might be more convenient in usage
    template<class E>
    VType& projection(const E& x){
      if(proj_.size() != dim_) //only resize when necessary, normally once
	proj_.resize(dim_);
      
      for(size_t i = 0; i < dim_; ++i)
	proj_[i] = (*this).operator()(i,x);
      return proj_;
    }

    const X& low() const {return low_;}
    const X& up() const {return up_;}

    template<class E>
    bool is_feasible_wrt_simple_bounds(const E& x){
      bool ct = true;
      for(size_t i = 0; i < dim_; ++i){
	if(!is_contained(x[i],low_[i],up_[i])){
	  ct = false;
	  break;   //leave loop betimes
	} 
      }
      return ct;
    }

  private:
    const X& low_;  //! these are fixed
    const X& up_;
    size_t dim_;
    VType proj_;
    };
  


  template<char C>          
  class MuStrategy{    //! standard: decrease mu by factor tau
  public:
    template<class T>
    static inline void applied(T& mu, const T& tau){
      mu *= tau;         
    }
  };

  
  template<>
  class MuStrategy<'n'>{  //! don't change mu
  public:
    template<class T>
    static inline void applied(T& mu, const T& tau){}
  };

  template<>
  class MuStrategy<'N'>{  //! don't change mu
  public:
    template<class T>
    static inline void applied(T& mu, const T& tau){}
  };


#if USE_CPPAD
  /**
   * \brief Projection algorithm for inner optimization/iteration of 
   * the augmented Lagrangian subproblem with <I> simple bounds</I>, i.e.
   * \f$ min_x\{L_{\mathrm{aug}}(x) | l \leq x \leq u \}\f$, with \f$ f\f$ denoting the
   * augmented Lagrangian function which involves <I> equality </I> constraints
   * from the original problem.
   *
   * \tparam T value type
   * \tparam PROBLEM equality constrained and simply constrained NLP problem
   * \tparam MU updates mu: if 'n' or 'N' don't update mu, else decrease it
   *
   * References:
   *
   * [1]  [CONN,GOULD and TOINT, "A globally convergent augmented Lagrangian method", SIAM J. Numer. Anal., Vol. 28, pp. 545-572, 1991] 
   */
  template<class T, template <class T2> class PROBLEM, char MU = 'a'> 
  class GloballyConvergentAugmentedLagrangianMethod{
  public:
    typedef PROBLEM<T> ProblemType;
    typedef typename ProblemType::VType VType;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef typename ExprTmpl::MyVec<int> IndexVecType;

    typedef MuStrategy<MU> MuUpdateType;

    //! TR subproblem solvers
    typedef SteihaugPCG<VType> TRSolverType;

    typedef TRRadius<T,'2'> RadiusType;


    //============= any norm ==============================================
    typedef Norm<'2',BaseType> NormType; //! standard norm in \f$ R^n\f$
    //=====================================================================


    enum{nvar = ProblemType::nvar,
	 neq = ProblemType::neq
    };

    GloballyConvergentAugmentedLagrangianMethod():numOfIterations_(0),iterInner_(0),cgMaxiter_(0),trMaxiter_(0),maxit_(0),eta0_(T()), mu0_(T()), omega0_(T()),tau_(T()),gamma1_(T()),omegaStar_(T()),etaStar_(T()),alphaW_(T()),betaW_(T()), alphaE_(T()),betaE_(T()),trtol_(T()),muMin_(T()),aboveEps_(T()),innerprod_(T()),delta_(T()),mk_(T()),ratio_(T()),deltaMin_(T()){}

    //! MyVec allows assignment from different racs having operator[]
    //! \param x, l start iterate for state and Lagrange multipliers, resp.
    //! \param sdiag, mdiag (fixed) p.d. diagonal scaling matrices 
    //!
    template<class V>
    GloballyConvergentAugmentedLagrangianMethod(const V& x, const V&lda, const V& sdiag, const V& mdiag, const std::string& filename, const V& reduced = V()):numOfIterations_(0),iterInner_(0),cgMaxiter_(0), trMaxiter_(0), maxit_(0), x0_(x), lambda0_(lda), Sdiag_(sdiag), Mdiag_(mdiag),red_(reduced),projGrad_(std::distance(x.begin(),x.end())),p_(std::distance(x.begin(),x.end())),eval_(std::distance(x.begin(),x.end())),xprev_(x.size()),ScalD_(ntimes<2>(std::distance(x.begin(),x.end()))){
      ParameterData FreeParam;
      FreeParam.read_from_file(filename);
      
      eta0_ = FreeParam.get_datum<T>("eta0");
      mu0_ = FreeParam.get_datum<T>("mu0");
      omega0_ = FreeParam.get_datum<T>("omega0");
      tau_ = FreeParam.get_datum<T>("tau");
      gamma1_ = FreeParam.get_datum<T>("gamma1");
      omegaStar_ =  FreeParam.get_datum<T>("omegaStar"); 
      etaStar_ = FreeParam.get_datum<T>("etaStar");
      alphaW_ = FreeParam.get_datum<T>("alphaW");
      betaW_ = FreeParam.get_datum<T>("betaW");
      alphaE_ = FreeParam.get_datum<T>("alphaE");
      betaE_ = FreeParam.get_datum<T>("betaE");
      //! for inner iteration -- TR version
      cgMaxiter_ = FreeParam.get_datum<size_t>("cgMaxiter");
      trMaxiter_ = FreeParam.get_datum<size_t>("trMaxiter");
      maxit_ = FreeParam.get_datum<size_t>("maxit");
      trtol_ =  FreeParam.get_datum<T>("trtol");
      muMin_ = FreeParam.get_datum<T>("muMin");
      aboveEps_ = FreeParam.get_datum<T>("aboveEps");

      innerprod_= T();
      delta_ = T();
      mk_ = T();
      ratio_ = T();
      deltaMin_ = FreeParam.get_datum<T>("deltaMin");

      //set up TR problem
      radius_.initialize(FreeParam.get_datum<T>("delta0"),FreeParam.get_datum<T>("deltaMax"),FreeParam.get_datum<T>("trgamma"));
      approxTRProblem_.initialize(x.size());
      
      symm_identity(Cprec_,x.size());  //just identity, i.e. no preconditioning
        
      //! set additional constraints depending on whether red.size() != 0 
      if(red_.size() != 0){
	DerivLaug_.get_problem().set_additional_bounds(red_);
	DerivLaug_.get_AD_problem().set_additional_bounds(red_); 
	prob_.set_additional_bounds(red_);
	Laug_.get_problem().set_additional_bounds(red_);
      }
      //! initialize derivatives
      DerivLaug_.set();

      std::cout << "CONSTR. invoked" <<std::endl;
    } //end constructor


    //! you should check if a solution point satisfies the equality constraints
    template<class X>
    VType& equality_constraints(const X& x){
      return prob_.equality_constraints(x); 
    }

    size_t number_of_iterations() const{return numOfIterations_;}
    size_t number_of_inner_iterations() const{return iterInner_;}
    size_t number_of_additional_constraints() const{return red_.size();}

    size_t max_iterations() const {return maxit_;}

    ProblemType& get_problem() {return prob_;}
    
    //! Find the optimum of problem using augmented Lagrangian to bind the
    //! equality constraints and the simple bounds on state variable x
    VType& optimum(const VType& Cprec = VType()){
      if(Cprec != VType())
	Cprec_ = Cprec;
      
      
      //! stores only const refs
      ProjectionOntoBoxConstraints<VType> Box(prob_.low(),prob_.up());
      
      T mu = mu0_, 
	alpha = std::min(mu,gamma1_), 
	omega = omega0_*std::pow(alpha,alphaW_), 
	eta = eta0_*std::pow(alpha,alphaE_);
      
      //these are for the inner iteration
      x_ = x0_;  lambda_ = lambda0_;  
	
      BaseType w, nmg0;

     
      //! MAIN LOOP
      size_t k;
      //size_t diffcount = 0;
      for(k = 0; k < maxit_; ++k) {  
	grad_ = DerivLaug_.gradient(x_,lambda_,Sdiag_,mu);
	//normalize(grad_);
	B_ = DerivLaug_.hessian(x_,lambda_,Sdiag_,mu);
       

	// //!Regularization
	// symmetric_matrix_2_full_matrix(Bfull_,B_);
	// CholeskyWithAddedMultipleOfIdentity<VType> RegularizeB(Bfull_);
	// RegularizeB.compute(0.001,100.5,'L',10);
	// std::cout << "Hessian made p.d. after "<< RegularizeB.number_of_iterations() << " iterations"<< std::endl;
	// full_matrix_2_symmetric_matrix(B_,Bfull_);
	
	//std::cout << "Mdiag_ = "<< Mdiag_ << std::endl;

	projGrad_ = x_ - Box.projection(x_ - Mdiag_*grad_); 

	std::cout << std::endl<< "ITERATION k = "<< k << std::endl<<
	  "---------------" << std::endl<<  std::setprecision(12) << "x_ = "<< x_ << "lambda_ = "<< lambda_ << std::setprecision(5) << std::endl<< "f(x_) = "<< prob_.objective(x_) << std::endl<<"c(x_) = "<< prob_.equality_constraints(x_) <<std::endl<< "grad_ = "<< grad_ << std::endl << "projGrad_ = "<< projGrad_ << std::endl;


	//===============================================================
	//INNER ITERATION -- mu, S, and lambda fixed
	// Goal: find an \f$ l \leq  x \leq u \f$ 
	//perform Trust Region algorithm here, see [N/W,algo. 4.1 p.69]
	//===============================================================	
	gp_ = projGrad_;
	nmg0 = NormType::norm(gp_);
	delta_ = radius_.get_delta0();
	iterInner_ = 0; //reset

	while(NormType::norm(gp_) > trtol_*std::max(1.,nmg0)){
	  //std::cout << "||gp_|| = "<< NormType::norm(gp_) << "  trtol_*std::max(1.,nmg0) = "<< trtol_*std::max(1.,nmg0) << std::endl;
	  iterInner_++;	    
	  //! compute direction
	  p_ = approxTRProblem_.approximation(B_,gp_,delta_,std::min(0.5,std::sqrt(NormType::norm(gp_))),Cprec_,cgMaxiter_);
	  
	  // std::cout << "# of PCG iterations: "<< approxTRProblem_.number_of_iterations() << std::endl;

	  mk_ = dot(gp_,p_) + 0.5*dot(p_,symm_matrix_vector_multiplication(B_,p_));
	  
	  eval_ = x_ + p_;
	  ratio_ = (Laug_(x_,lambda_,Sdiag_,mu) - Laug_(eval_,lambda_,Sdiag_,mu))/(-mk_);
	  radius_.update(delta_,ratio_,p_);
	  // std::cout << "delta = "<< delta_ << std::endl;
	  if(delta_ < deltaMin_)  //prevent delta_ of becoming too small
	    delta_ = deltaMin_;

	  if(radius_.is_direction_changed() == true){
	    //! when direction is too small, we won't get a sufficient reduction
	    //! of the 2nd order Taylor series model 
	    if(is_zero(NormType::norm(p_),aboveEps_)){
	      std::cout << "p_ too small ==> Too less progress towards direction p_. Break TR iteration" << std::endl;
	      break;
	    }
	    x_ += p_;
      
	    //! check if point is feasible 
	    // if(Box.is_feasible_wrt_simple_bounds(x_) == false){
	    //   x_ = Box.projection(x_);
	    // }
	    
	    grad_ = DerivLaug_.gradient(x_,lambda_,Sdiag_,mu);
	    //normalize(grad_);
	    gp_ = x_ - Box.projection(x_ - Mdiag_*grad_); 
	    
	    B_ = DerivLaug_.hessian(x_,lambda_,Sdiag_,mu);
	  
	    // //!Regularization
	    // symmetric_matrix_2_full_matrix(Bfull_,B_);
	    // CholeskyWithAddedMultipleOfIdentity<VType> RegularizeBIn(Bfull_);
	    // RegularizeBIn.compute(0.001,100.5,'L',10);
	    // std::cout << "Hessian made p.d. after "<< RegularizeBIn.number_of_iterations() << " iterations"<< std::endl;
	    // full_matrix_2_symmetric_matrix(B_,Bfull_);
	    
	  }

	  if(iterInner_ >= trMaxiter_){
	    ADONIS_WARNING(Warning, "Too many iterations in TR algorithm (max. numb. iters = "<< trMaxiter_<<") \n   Taking current iterate as best approximation available so far...");
	  break;
	  }
	  
	}//end while
	//===========================================
	//END IN. ITERATION
	//===========================================
	std::cout << "Number of TR (inner) iterations: "<< iterInner_ << std::endl;

	//! check if point is feasible 
	// if(Box.is_feasible_wrt_simple_bounds(x_) == false){
	//   x_ = Box.projection(x_);
	// }
       
	//std::cout << "x_ = "<< x_ << std::endl;

	//======================================================================
	//======================================================================

	//now we have computed the iterate \f$ x_k\f$
	//projGrad_ = x_ - Box.projection(x_ - Mdiag_*grad_); 
	w = NormType::norm(projGrad_);  //gp_??
	std::cout << "||x -P(x - gradf)|| = "<< w  << std::endl; 
	std::cout<< "penalty mu = "<< mu   <<  "  alpha = "<< alpha << "  omega = "<< omega << "  eta = "<<  eta << std::endl << " lambda = "<< lambda_<< std::endl ;

	//! leave loop if projection is smaller than variable tolerance
	if(w <= omega){
	  //std::cout << "<= omega detected" << std::endl;
	  break;
	}


	//========== OUTER ITERATION ========================================
	eq_ = prob_.equality_constraints(x_);

	//!test for convergence and update Lagrange multiplier estimates
	std::cout << "c(x) = "<< eq_ << std::endl <<"||c(x)|| = "<< NormType::norm(eq_) << std::endl;
	
	if(NormType::norm(eq_) <= eta){
	  if((NormType::norm(projGrad_) <= omegaStar_) && (NormType::norm(eq_) <= etaStar_)){
	    break;  //x_ is approximate solution
	  }
	  //else{
	    //! 1st order Lagrange multipier estimates, don't change mu
	  lambda_ += (Sdiag_*eq_/mu ); 
	  //don't update mu
	  alpha = std::min(mu,gamma1_);
	  omega *= std::pow(alpha,betaW_);
	  eta *= std::pow(alpha,betaE_);
	    //}
	}
	else{ //! decrease penalty parameter and update Lagrange multipl. estim.
	  //don't update lambda
	  //mu *= tau_;
	  MuUpdateType::applied(mu,tau_);
	  alpha = std::min(mu,gamma1_);
	  omega = omega0_*std::pow(alpha,alphaW_);
	  eta = eta0_*std::pow(alpha,alphaE_);
	}

	//! Safeguard
	if(mu < muMin_) {  //note mu is basis for calculation
	  ADONIS_INFO(Information,"mu < muMin_; Apply safeguard");
	  mu = mu0_;  //muMin_;  //MARC: 
	  //!set omega and eta to reasonable values as well
	  omega = omega0_; //omegaStar_;
	  eta = eta0_;  //etaStar_;
	}

	 
      } //end  k=0,1,2,...
     
      numOfIterations_ = k; //performed iterations
      
      
      // if(numOfIterations_ >= maxit_-1)
      // 	ADONIS_INFO(Information, "Maximum number of iterations reached (maxit = "<< maxit<<". \n   Errorneous result might be produced...");
       
      return x_;   //give back alleged local minimum
      //return ( x_ = Box.projection(x_) );
    }//end of function


    

  private:
    size_t numOfIterations_, iterInner_, cgMaxiter_, trMaxiter_, maxit_;
    
    VType x0_, lambda0_, Sdiag_, Mdiag_,x_, lambda_, red_,grad_, B_, projGrad_, p_, gp_, eval_,eq_, xprev_,gradprev_, Bfull_,
      Cprec_;
    T eta0_, mu0_, omega0_,tau_,gamma1_,omegaStar_,etaStar_,alphaW_,betaW_,
      alphaE_,betaE_,
      trtol_, muMin_,aboveEps_,innerprod_, delta_, mk_,ratio_,deltaMin_;
    
    ProblemType prob_;
    AugmentedLagrangianFunction<T,PROBLEM> Laug_;

    VType ScalD_;
    RadiusType radius_;
    TRSolverType approxTRProblem_;

    DerivativesOfAugmentedLagrangianFunction<T,PROBLEM> DerivLaug_;

  };

#endif //USE_CPPAD in 'DerivativesOfAugmentedLagrangianFunction<T,PROBLEM>'

} //end of namespace 

#endif
