#ifndef SOME_HAND_REDUCED_SOURCE_TERMS_HH
#define SOME_HAND_REDUCED_SOURCE_TERMS_HH

#include "../common/error.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../linalg/linearsystemsolvers.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../derivatives/jacobianofsource.hh"

#include "../ode/sourceterms.hh"  //full ODE system

namespace Adonis{

 

  /**
   * \brief Reduced version of a toy mechanism. The full mechanism can be found
   * in file '../ode/sourceterms.hh'
   */
  template<class T>
  class ReducedToyMechanism{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type DType;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef ExprTmpl::MyVec<T> VType;
    typedef std::size_t SizeType;
    typedef Norm<'2',BaseType> NormType;
    
    typedef ToyMechanism<DType> FullModelType;
    typedef JacS<DType,ToyMechanism,ExprTmpl::MyVec> JacobianType;

    enum{RPV=2   //index of rpv
    };

    ReducedToyMechanism(SizeType n = 1):maxIt_(13),     count_(1),
					maxNewt_(5),
					tol_(1.e-10),
					newtTol_(1.e-08),
					fullode_(3),
					kf2_(0.01),kr2_(1.e-05),fulldim_(3),rhs_(n),zM_(3), zM0_(3), C_(3*1), b_(1), residual_(1), z_nu_(3), g_(1), defect_(3), G_(3),Gprime_(3*3),up_(3), low_(3){

      //======== TODO: change setting here================= ==================
      zM0_ <<= 0.35, 0.25, 0.4;

      kn_ = 0.7; //1.25e-01;  //should be of approximately the same magnitude as
                       //in the outer scheme 
      whichTimeStepper_ = 'E'; //e, E = explicit Euler, i, I = impl. Euler
      //======================================================================

      //!lower and upper bounds
      low_ <<= 0.,0.,0.;
      up_ <<= 1.,1.,1.;

      //! conservation relation x[0] + x[1] + x[2] = 1.0
      C_ <<= 1., 1., 1.;
      b_ <<= 1.;

      zM_ = zM0_;
      
      NablaChem_.set(fulldim_,fulldim_); 

}

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}


    template<class X>
    VType& operator()(const X& red){

      zM_[RPV] = red[0];  //assign rpv to zM_
      
      residual_ = matrix_vector_product(C_,zM_) - b_;
      
      if(NormType::norm(residual_) > 1.e-09){
	//!perform "Gauss-Newton" iteration
	z_nu_ = zM_;
	int m = (int)b_.size();
	for(SizeType nu = 1; nu <= maxIt_; ++nu){
	  count_ = nu;
	  g_ = -matrix_vector_product(C_,z_nu_) + b_;
	  defect_ = solve_rank_deficient_ls(C_,m,g_,1); 
	  z_nu_ += defect_;
	  if(NormType::norm(defect_) <= tol_)
	    break;  //leave loop since convergence has been achieved
	}
	if(count_ == maxIt_)
	  ADONIS_INFO(Information, "Max. number of iterations reached (maxit = "<< maxIt_ << ").\n   No convergence with tol = "<< tol_ << ".");
      
	zM_ = z_nu_;
      }
      
      zM_ = Max(low_,Min(zM_,up_));

      //! PERFORM 1 IMPLICIT or EXPL. EULER TIME STEP
      if(whichTimeStepper_ == 'i' ||  whichTimeStepper_ == 'I'){
	//std::cout << "IMPL. EULER" << std::endl;
	z_prev_ = zM_;
	z_nu_ = z_prev_; //reuse z_nu_
	count_ = 1; //reset
	for(SizeType l = 1; l <= maxNewt_; ++l){
	  G_ = -(z_nu_ - z_prev_ - kn_*fullode_(z_nu_));
	  Gprime_ = -kn_*NablaChem_.jacobian(z_nu_);
	  //print_in_matrix_style(Gprime_,fulldim_);
	  update_diagonal<AddBasicElements>(Gprime_,fulldim_,1.);
	  good_square_solve(Gprime_,fulldim_,G_,1);
	  z_nu_ += G_;

	  if(NormType::norm(G_) <= newtTol_){
	    //---- ASSIGN -------
	    zM_ = z_nu_;
	    //---------------------
	    break;
	  }

	  if(l == maxNewt_){
	    ADONIS_INFO(Information, "Maximum number of Newton iterations ("<<maxNewt_<<") reached. Decrease kn_....");
	    l = 0; //set back i and repeat with kn/2 
	    kn_ *= 0.5;
	  }

	  adonis_assert(Abs(kn_) >= 1.e-12); //don't let kn_ become too small

	  count_++;

	  if(count_ == 4*maxNewt_){
	    ADONIS_INFO(Information, "After 4 iterations with decreasing stepsize, there hasn't been encountered any convergence...");
	    break;  //o.k. 4 times with decreasing step size yielded no result
	  }
	} //end NEWTON iter
      }
      
      else if (whichTimeStepper_ == 'e' || whichTimeStepper_ == 'E'){
	//perform explicit EULER step
	//std::cout << "EX. EULER" << std::endl;
	zM_ += kn_*fullode_(zM_);
      }
      else
	ADONIS_ERROR(NotImplementedError,"Time stepping method whichTimeStepper_ = "<<whichTimeStepper_ << " undefined.");

      //now zM_ has been determined to the best of my knowlede and belief
	

      //!reduced source term
      rhs_[0] = kf2_*zM_[1] - kr2_*zM_[2];

      return rhs_;
    }

  private:
    SizeType maxIt_, count_, maxNewt_;
    BaseType tol_, newtTol_;
    FullModelType fullode_;
    value_type kf2_, kr2_;
    SizeType fulldim_;
    VType rhs_, zM_, zM0_, C_, b_,residual_, z_nu_, g_, defect_, G_,Gprime_,
	  up_, low_, z_prev_;
  
    value_type kn_;
    char whichTimeStepper_;
    JacobianType NablaChem_;
  
  };

}//end namespace

#endif

