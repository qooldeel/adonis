#ifndef SPECIES_RECONSTRUCTION_HH
#define SPECIES_RECONSTRUCTION_HH

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../common/elementaryoperations.hh"
#include "../common/typeadapter.hh"
#include "../common/smartassign.hh"

#include "../expressiontemplates/exprvec.hh"
#include "../templatemetaprograms/unrollloop.hh"

#include "../derivatives/jacobianofsource.hh"
#include "../linalg/linearsystemsolvers.hh"
#include "../misc/operations4randomaccesscontainers.hh"

#include "massbalanceinvariant.hh"
#include "indexmanipulation.hh"
#include "rvp.hh"
#include "bounds.hh"
#include "../dunestuff/fmatrix.hh"
#include "../dunexternal/dunextensions.hh"
namespace Adonis{

  /**
   * \brief This is only an alternative attempt aiming at replacing the
   * NLP formulation
   *
   * \tparam T value type
   * \tparam N full dimension
   * \tparam R reduced dimension
   * \tparam FUN full source term
   * \tparam ENCO integer encoding for mechanism needed e.g. for bounds on vars 
   */
  template<class T, int N, int R, template<class D> class FUN, int ENCOD, char NORM = '2'> 
  class SpeciesReconstruction{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType; 
    typedef typename TypeAdapter<T>::Type DType;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef Norm<NORM,typename TypeAdapter<T>::BaseType> NormType;
    typedef JacS<T,FUN,ExprTmpl::MyVec> JacobianType;

    typedef MassBalanceInvariant<VType,NORM> LinearInvariantType;
    
    typedef ExprTmpl::MyVec<size_t> IndexVType; 
    typedef Dune::FieldMatrix<T,R,N> BTType;
    typedef Dune::FieldMatrix<T,N,N-R> UType;

    typedef MechanismBounds<DType,ENCOD> MechBdsType;

    //========= maxnumber of newton iters ============
    enum{mxit = 5};
    //================================================

    
    //default construction 
    //!NOTE: since in the default constructor of fun_ we set the stuff for the
    //!      calc. of the source, it is mandatory to initialize it here!!
    SpeciesReconstruction():isConverged_(false),isPerturbed_(false),kn_(T()),kEstim_(T()),newtTol_(T()),fun_(N),indexptr_(0),whichTimeStepper_('e'),maxit_(0),nmtol_(0.),Tol_(0.){}

   
    size_t size() const {return v_.size();}

    
    template<class V, class W>
    void initialize(const V& v, const T& kn, const BaseType& newtTol, size_t* ixit, const W& C, const W& b, char wts = 'e', int maxit= 12, const BaseType& nmtol = 1.e-09, const BaseType& Tol = 1.e-09 ){
      adonis_assert((int)std::distance(v.begin(),v.end()) == N);
      adonis_assert((int)fun_.dim() == N); 
      
      isConverged_ = false;
      isPerturbed_ = false;
      v_.resize(N);
      vl_.resize(N);
      v_prev_.resize(N);
      farg_.resize(N);
      Dg_.resize(N*N);   //Jacobian
      C_ = C;
      b_ = b;
      g_.resize(N);
      //dotomega_.resize(N);
      index_.resize(R);
      
      for(int i = 0; i < N; ++i)
	smart_assign(v_[i],v[i]); //full composition with zM initialized
      
      
      NablaChem_.set_with_init(N,N,v_);
      
      kn_ = kn;
      kEstim_ = kn;
      newtTol_ = newtTol;
      indexptr_ = ixit;
      whichTimeStepper_ = wts;
      maxit_ = maxit;
      nmtol_ = nmtol;
      Tol_ = Tol;

      UnrollLoop<0,R>::assign(index_.begin(),indexptr_);
      IndexManipulator<N,R,ExprTmpl::MyVec> IM;
      IM.create_unrepresented_index(index_);
      ReactionProgressVariables<T,N,R> RPV(BT_,U_); //references 
      RPV.create_B_T(index_);     //now B^T is filled 
      RPV.create_U(IM.get_index()); //now U is filled

      mechbds_.initialize(N);  //! set bounds on variables of mechanism
    }


    const BTType& B_T() const {return BT_;}
    const UType& U() const {return U_;}


    template<class V, class W>
    VType& evaluate(const V& red, const W& full){
      //assign full estimate to v_
      //if(full.size() != 0) //otherwise default construction is used
      UnrollLoop<0,N>::assign_rac2rac(v_,full); //incorporates smart_assign

      //! \f$ B^T u = r\f$, i.e. replace appropriate components of \f$u\f$ 
      //! with those of reduced type \f$ r\f$ 
      UnrollLoop<0,R>::assign_smaller_2_greater_rac(v_,red,indexptr_);
      
      // V lh(R);
      // matrix_vector_product(lh,BT_,v_);
      // lh = -lh + red;
      // UnrollLoop<0,R>::assign_smaller_2_greater_rac(v_,lh,indexptr_);
      // VType A, l, bt(R*N);
      // int k = 0;
      // for(int i = 0; i < R; ++i){
      // 	for(int j = 0; j < N; ++j){
      // 	  bt[k++] = BT_[i][j];
      // 	}
      // }

      // join(A,C_,bt);
      // join(l,b_,red);

      VType r(R);
      LinearInvariantType LinVar(C_,b_); //stores const references only
      //LinearInvariantType LinVar(A,l);

      //! check if mass balance is (approximately) violated
      isPerturbed_ = false;
      r = matrix_vector_product(C_,v_) - b_; //!residual
      if(NormType::norm(r)  > nmtol_){
	//std::cout << "C = "<< C_ << " b = "<< b_ << std::endl; //is const
	std::cout << " r = "<< r << std::endl << " ||r|| = "<< NormType::norm(r) << std::endl;
	//! a few minimum norm iterations should work
	v_ = LinVar.enforce_balance(maxit_,Tol_,v_);

	isPerturbed_ = true;
      }

      //! project back onto box constraints
      v_ = Max(mechbds_.low(),Min(v_,mechbds_.up()));
      
     
      // if(isPerturbed_){
      if(whichTimeStepper_ == 'i' || whichTimeStepper_ == 'I'){
	//!save old approx -- for IMPLICIT TIMESTEPPER
	//!========= maybe better guesses may be possible =====================
	v_prev_ = v_; 
	vl_ = v_prev_;
	//!====================================================================

	size_t count = 0;
	for(int i = 0; i < mxit; ++i){

	  //! calculate Dg·x= -g //IMPLICIT MIDPOINT RULE
	  //! impl. midpoint rule
	  // farg_ = 0.5*(vl_ + v_prev_);  
	  // g_ = -(vl_ - v_prev_ - kn_*fun_(farg_)); 
	  // NablaChem_.jacobian(farg_);   //vl_
	  // Dg_ = -0.5*kn_*NablaChem_.get_jacobian();

	  // //! implicit euler
	  //dotomega_ = fun_(vl_);
	  g_ = -(vl_ - v_prev_ - kn_*fun_(vl_)); //dotomega_);    //this is -g
	  Dg_ = -kn_*NablaChem_.jacobian(vl_);
	
	
	  update_diagonal<AddBasicElements>(Dg_,N,1.); //Now: I -0.5kn·Df
	  good_square_solve(Dg_,N,g_,1);
		
	  //update with newton direction
	  vl_ += g_;

	  // //! DOES NOT SEEM TO IMPROVE THINGS HERE
	  // //! check if mass balance is (approximately) violated
	  // r = matrix_vector_product(C_,vl_) - b_;
	  // if(NormType::norm(r)  > 1.e-05){
	  //   vl_ = LinVar.enforce_balance(30,1.e-05,vl_);
	  // }

	  //early convergence (within tolerance newtTol_)
	  if(NormType::norm(g_) <= newtTol_){
	    //---- ASSIGN -------
	    v_ = vl_;
	    // kn_ *= 2;
	    // if(kn_ > 5.75e+02)
	    //   kn_ = kEstim_;

	    isConverged_ = true;
	    //---------------------
	    break;
	  }
	
	  if(i == mxit){
	    ADONIS_INFO(Information, "Maximum number of Newton iterations ("<<mxit<<") reached");
	    i = 0; //set back i and repeat with kn/2 
	    kn_ *= 0.5;
	  }

	  adonis_assert(Abs(kn_) >= 1.e-12); //don't let kn_ become too small

	  count++;

	  if(count == 4*mxit){
	    ADONIS_INFO(Information, "After 4 iterations with decreasing stepsize, there hasn't been encountered any convergence...");
	    v_ = vl_;
	    break;  //o.k. 4 times with decreasing step size yielded no result
	  }
	} //end Newton iteration
      
      } //end implicit
      //!just perform one explicit Euler step with small kn_
      if(whichTimeStepper_ == 'e' || whichTimeStepper_ == 'E'){ 
	//dotomega_ = fun_(v_);
	v_ += kn_*fun_(v_); 
      }
      // } is perturbed
      // //assign values of rpv to full composition
      //UnrollLoop<0,R>::assign_smaller_2_greater_rac(v_,red,indexptr_);
      
      return v_;
    }

    const VType& get_z() const {return v_;}
    VType& get_z() {return v_;}

    
    //const VType& dot_omega() const {return dotomega_;}

    // const T& operator[](int i) const{
    //   adonis_assert(i>=0 && i<N);
    //   return dotomega_[i];
    // }

  private:
    mutable bool isConverged_, isPerturbed_;
    VType v_,
      vl_,
      v_prev_,
      farg_,
      Dg_,
      g_;
    //VType dotomega_;
    VType C_,  //mass balance matrix
      b_;  //sumed up balances 
    T kn_,
      kEstim_,
      newtTol_;
    FUN<T> fun_;
    size_t* indexptr_;
    char whichTimeStepper_;
    BaseType maxit_,nmtol_,Tol_;
    JacobianType NablaChem_;
    

    IndexVType index_; 
    BTType BT_;
    UType U_;
    MechBdsType mechbds_;
  };



}

#endif
