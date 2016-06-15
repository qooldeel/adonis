#ifndef  DERIVATIVE_OF_AUGMENTED_LAGRANGIAN_HH
#define  DERIVATIVE_OF_AUGMENTED_LAGRANGIAN_HH

#if USE_CPPAD
#include <cppad/cppad.hpp>
#endif

#include "../common/adonisassert.hh" 
#include "../common/error.hh"

#include "../misc/operations4randomaccesscontainers.hh"
#include "../templatemetaprograms/matrixunroller.hh"
#include "../templatemetaprograms/unrollloop.hh"

namespace Adonis{

  template<int I, int NEQ>
  class CalculateTermsOfHessianOfAugLagr{
  public:
  private:
    enum{go_ = (I+1) != NEQ};  //!sum over NEQ only!!!!
  
  public:
    typedef CalculateTermsOfHessianOfAugLagr<go_ ? (I+1) : NEQ,NEQ> NextIterType; 

    //!NOTE: NVAR has to be given explicitly, denoting the number of variables
    template<int NVAR, class X, class L, class SDIAG, class T, class HCI, class DX, class M1, class V, class JACEQ, class DY, class EQ, class SUM>
    static inline DY& nicely(const X& x, const L& lambda, const SDIAG& sdiag, const T& mu, HCI& hessci, DX& dx, M1& m1, V& v, const JACEQ& jaceq, DY& dy, const EQ& eq, SUM& m2){
      hessci = dx.Hessian(x,I);
      m1 += lambda[I]*hessci;                 //m1 is changed ==> reset
      UnrollLoop<0,NVAR>::template get_row<I>(v,jaceq); 
      MatrixUnroller<0,0,NVAR,NVAR>::dyad(dy,v,v);
      m2 += sdiag[I]*(dy + eq[I]*hessci);
      
      return NextIterType::template nicely<NVAR>(x,lambda,sdiag,mu,hessci,dx,m1,v,jaceq,dy,eq,m2);
    }
  };

  
  //!specialisation -- end of recursion
  template<int NEQ>
  class CalculateTermsOfHessianOfAugLagr<NEQ,NEQ>{
  public:
    template<int NVAR, class X, class L, class SDIAG, class T, class HCI, class DX, class M1, class V, class JACEQ, class DY, class EQ, class SUM>
    static inline DY& nicely(const X& x, const L& lambda, const SDIAG& sdiag, const T& mu, HCI& hessci, DX& dx, M1& m1, V& v, const JACEQ& jaceq, DY& dy, const EQ& eq, SUM& m2){
      return (m2 /= mu);
    }
  };



#if USE_CPPAD
  /**
   * \brief I consider derivatives of the following augmented Lagrangian:
   * \f[ L_{\mathrm{aug}}(x,\lambda,S,\mu):= f(x) + \langle \lambda, c(x) \rangle + \frac{1}{2\mu} \sum_{i=1}^m s_{ii}c_i(x)^2.\f]
   * The derivatives of 1st order  can be readily obtained from \f$ \nabla_x f, \nabla_x c_i\f$ and those of 2nd order need \f$ \nabla^2_xx f, \nabla^2_xx c_i\f$ in addition. 
   */
  template<class T, template <class D> class PROBLEM>
  class DerivativesOfAugmentedLagrangianFunction{
  public:
    typedef T value_type;

    typedef CppAD::AD<T> ADType;
    typedef PROBLEM<ADType> ADProblemType;
    typedef typename ADProblemType::VType ADVecType;
    typedef CppAD::ADFun<T> ADFunType;
    typedef PROBLEM<T> ProblemType;
    typedef typename ProblemType::VType VType;

    DerivativesOfAugmentedLagrangianFunction():isSet_(false){}
    
    //! set up AD 
    void set(){
      if(!isSet_){  //o.k. if we haven't set up problem, do it now
	std::cout << "DLaug::set() ==> prob_.get_red() = "<< prob_.get_red() << std::endl;
	ADVecType Xfun(ADProblemType::nvar), Yfun(1), Xeq, Yeq;
      
	if(ProblemType::neq > 0){ //this makes only sense if we have constraints
	  Xeq.resize(ADProblemType::nvar); 
	  Yeq.resize(ADProblemType::neq);
	}

	//AD sequence for objective function
	CppAD::Independent(Xfun);
	Yfun[0] = adprob_.objective(Xfun);
	seqfun_.Dependent(Xfun,Yfun);
	seqfun_.optimize();
	
	if(ProblemType::neq > 0){ //when constraints are existant
	  //AD sequence for equality constraints
	  CppAD::Independent(Xeq);
	  Yeq = adprob_.equality_constraints(Xeq);
	  seqeq_.Dependent(Xeq,Yeq);
	  seqeq_.optimize();
	}
	
	//set eval
	eval_.resize(ADProblemType::nvar); //nvar vector
	lave_.resize(ADProblemType::nvar);
	dy_.resize(ADProblemType::nvar*ADProblemType::nvar);
	//!1st sum = \f$ \sum lambda_i \nabla^2c_i\f$
	firstsum_.resize(ADProblemType::nvar*ADProblemType::nvar);
	//! 2nd sum = \f$ \frac{1}{\mu}\sum s_{ii}\left( \nabla c_i \nabla c_i^T + c_i\nabla^2 c_i\right) \f$
	secondsum_.resize(ADProblemType::nvar*ADProblemType::nvar);

	isSet_ = true;
      }
    }


    template<class V1, class V2, class V3>
    VType& gradient(const V1& x, const V2& lambda, const V3& Sdiag, const T& mu){
      adonis_assert(ProblemType::neq == static_cast<int>(lambda.size())); 
      adonis_assert(isSet_); //when thrown function .set() hasn't been invoked

      typedef MatrixUnroller<0,0,ProblemType::neq,ProblemType::nvar> TwoLoopsType;

      gradf_ = seqfun_.Jacobian(x); //overwrites gradf_ at each call
      //! each row of the equality constraint jacobian is \f$ \nabla_x c_i^T\f$

      if(ProblemType::neq > 0){ //!when constraints are existant
	jaceq_ = seqeq_.Jacobian(x);
	//! RESET: same as eval_ = T() and lave_ = T() but at compile time;
	UnrollLoop<0,ProblemType::nvar>::assign_value(eval_.begin(),T());
	UnrollLoop<0,ProblemType::nvar>::assign_value(lave_.begin(),T());
      
	eq_ = prob_.equality_constraints(x); //overwrites eq_ again

	gradf_ += ( TwoLoopsType::multiply_vector_entry_with_vector(eval_,lambda,jaceq_) + 1/mu*TwoLoopsType::multiply_vector_entry_with_vector(lave_,eq_,jaceq_,Sdiag) );
	//! alternative evaluation
	//eq_ = lambda + Sdiag*eq_/mu; //actually 1st order Lagr. mlplr update
	//gradf_ += ( TwoLoopsType::multiply_vector_entry_with_vector(eval_,eq_,jaceq_) );
      }
      return gradf_; //!now this is \f$ \nabla L_{\mathrm{aug}} \f$
    }
    
    
    //!Hessian of augmented Lagrangian
    template<class V1, class V2, class V3>
    VType& hessian(const V1& x, const V2& lambda, const V3& Sdiag, const T& mu){
      adonis_assert(ProblemType::neq == static_cast<int>(lambda.size())); 
      adonis_assert(isSet_);  //when thrown function .set() hasn't been invoked
       
      hessFullf_ = seqfun_.Hessian(x,0);
    
      if(ProblemType::neq > 0){ //!when constraints are existant
	const int nvar2 = ProblemType::nvar*ProblemType::nvar;
	jaceq_ = seqeq_.Jacobian(x);
	//std::cout << "A = "<< jaceq_ << std::endl;
	//UnrollLoop<0,nvar2>::assign_value(dy_.begin(),T());   //reset
	UnrollLoop<0,nvar2>::assign_value(firstsum_.begin(),T());    //reset
	UnrollLoop<0,nvar2>::assign_value(secondsum_.begin(),T());    //reset

	eq_ = prob_.equality_constraints(x);

	//overwrites dy with \f[ \frac{1}{\mu}\sum_{i=1}^m s_{ii}\left(\nabla c_i \nabla c_i^T + c_i \nabla^2c_i \right) \]
	CalculateTermsOfHessianOfAugLagr<0,ProblemType::neq>::template nicely<ProblemType::nvar>(x,lambda,Sdiag,mu,hessFullci_,seqeq_,firstsum_,eval_,jaceq_,dy_,eq_,secondsum_); 

	hessFullf_ += (firstsum_ + secondsum_);
      }

      //! only store upper part of Hessian and return stuff
      full_matrix_2_symmetric_matrix(hessSym_,hessFullf_);

      return hessSym_;
    }

   

    ProblemType& get_problem() {return prob_;}

    ADProblemType& get_AD_problem() {return adprob_;}

  private:
    ADProblemType adprob_;   //for derivatives only
    ADFunType seqfun_, seqeq_;
    mutable bool isSet_;
    VType gradf_, jaceq_, hessFullf_, hessFullci_, hessSym_, eval_, lave_,
						   dy_, eq_, firstsum_, 
						   secondsum_;
    ProblemType prob_;

  };
#endif //CPPAD

}

#endif //namespace
