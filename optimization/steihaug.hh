#ifndef STEIHAUG_PRECONDITIONED_CONJUGATE_GRADIENT_TRUST_REGION_HH
#define STEIHAUG_PRECONDITIONED_CONJUGATE_GRADIENT_TRUST_REGION_HH

#include "../common/isclass.hh"
#include "../common/typeselector.hh"
#include "../common/typeadapter.hh"
#include "../common/globalfunctions.hh"
#include "../common/error.hh"
#include "../common/adonisassert.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../lapackwrapper/lapackdrivertraits.hh"
#include "../lapackwrapper/wrappertraits.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../misc/misctmps.hh"
#include "../preconditioner/precond.hh"

namespace Adonis{

  /**
   * \brief Steihaug's preconditioned conjugate gradient method for the
   * approximate solution of the TR subproblem \f[ \min\{m(p):= (\nabla f(x),p) + \frac{1}{2}(Bp,p) | \|p\|_C \leq \Delta\}, \f]
    where
    *  \f$ B\f$ is a real \f$ n \times n \f$ symmetric matrix, and \f$ m(p)\f$ is a linear quadratic approximation to a smooth function \f$f:R^N \rightarrow R \f$
    * Use this method in conjunction with any TR algorithm, cf. e.g. [2, § 5, Algo (5.1.6), p. 96]
    * References:
    *  
    * [1] [STEIHAUG, "The conjugate gradient method and trust regions in large scale optimization", SIAM J. Numer. Anal., Vol. 20, No. 3, 626-637, 1983]
    *
    * [2] [FLETCHER, "Practical methods of optimization", 2nd ed, John Wiley & Sons, 1987]
    *
    * \tparam SYMMTX (dense or sparse) symmetric \f$ n \times n\f$ matrix. In 
    *  the <B>dense symmetric case, only the lower part (column-wise, which is tantamount to storing the upper triangular part row-wise) is stored</B>.
    *
    * \tparam NT norm type for trust region ('1', '2', 'i' or 'I'). Note that
    * \tparam WHATNM  ('o'. 'O' uses an ordinary norm specified by NT 
    *   (see above). Use 's' or 'S' for <I>scaled norm</I>, cf. Steihaug, [1],
    *   uses a scaled norm \f$ \|p\|_C := \sqrt{\langle p, Cp\rangle}\f$, 
    *   where \f$C\f$ is a s.p.d. matrix. In the latter case NT becomes
    *   insignificant) 
    *
    */
  template<class SYMMTX, char NT = '2', char WHATNM = 'o'>
  class SteihaugPCG{
  public:
    typedef typename ValueTypeSelector<SYMMTX,IsADuneFieldContainer<SYMMTX>::Value>::value_type value_type;
    typedef SYMMTX SymmetricMatrixType;

    typedef ExprTmpl::MyVec<value_type> VType; 
    typedef typename TypeAdapter<value_type>::BaseType BaseType; 
    //********* set norm here *****************
    //typedef Norm<NT,BaseType> NormType;
    typedef NormWrapper<WHATNM,NT> NormType;  //'o'rdinary or 's'caled norm
    //*****************************************
    //********* set preconditioner here *******
    typedef Preconditioner<1,'y'> PrecondType;   //1=SSOR, 'y'= dense symmetric
    //*****************************************


    SteihaugPCG(size_t n = 0):det_(value_type()),p_(n),pp1_(n),sh_(n),iter_(0){}

    //! if default constructor is used and you wanna initialize it later on
    void initialize(size_t n){
      det_ = value_type();
      p_.resize(n);
      pp1_.resize(n);
      sh_.resize(n);
      iter_ = 0;
    }

    size_t number_of_iterations() const {return iter_;}

    //! return approximated direction \f$p \f$
    //! the tolerance should be taken similar to \f$ \min(0.5,\sqrt{\|\nabla f_k\|}\f$, in which case one obtains a rapid convergence (superlinear)
    VType& approximation(const SYMMTX& B, const VType& gradf, const value_type& delta, const value_type& tol, const SYMMTX& C = SYMMTX(), size_t maxit = 15){
          
      int n = static_cast<int>(gradf.size()),
	ld = std::max(1,n),
	nrhs = 1,
	info;
      char uplo = 'L';
    
      if(C != SYMMTX()){ 
	//!Cholesky decomposition of C; have to be performed only once
	symmetric_matrix_2_full_matrix(Cfullmtx_,C); //only needed for LAPACK
	F77Type<TypeTraits<value_type>::Value>::POTRF(&uplo,&n,&Cfullmtx_[0],&ld,&info);
      }
      
      
      //!STEP 1:
      p_ = value_type();   //p_ = p0 = 0
      g0_ = gradf;
      r_ = gradf;
      rTilde_ = r_;
      if(C != SYMMTX()){
	//!preconditioning -- B could be preconditioned as well
	F77Type<TypeTraits<value_type>::Value>::POTRS(&uplo,&n,&nrhs,&Cfullmtx_[0],&ld,&rTilde_[0],&ld,&info);   //rTilde is overwritten
	d_ = rTilde_;
      }
      else
	d_ = PrecondType::prec(B,rTilde_);
     
      value_type rrtil = dot(r_,rTilde_),
      rrtil1;
      bool flag = true;
      value_type gamma, alpha,
	a,b,c,tau, //needed in the execution of step (2.1) and (2.2), cf. [1] 
	bek, pred1, pred2;
     
      
      iter_ = 0; //reset counter
      //MAIN LOOP
      while(NormType::norm(r_,C) > tol*std::max(1.e-12,NormType::norm(g0_,C))){
	iter_++;
	//! STEP 2:
	Bd_ = symm_matrix_vector_multiplication(B,d_);//matrix_vector_product(B,n,d_);
	gamma = dot(d_,Bd_);
	
	//! STEP 3:
	if(gamma > 0){         //positive curvature
	  alpha = rrtil/gamma;
	  pp1_ = p_ -  alpha*d_;
	  //! execute (2.2) of [1]  -- iterate leaves TR
	  if((flag == true) && (NormType::norm(pp1_,C) >= delta)){
	    //! perform quadratic interpolation
	    c = dot(p_,p_) - ntimes<2>(delta);
	    b = -2.*dot(p_,d_);
	    a = dot(d_,d_);
	    det_ = ntimes<2>(b) - 4.*a*c;
	    if(det_ < value_type()){ //beware of negative roots
	      ADONIS_ERROR(DerivedError, "det_ < 0; only non-complex roots allowed here");
	    }
	    //! positive root of quadtratic equation \f$ ax^2 +bx +c = 0 \f$
	    tau = (-b + std::sqrt(det_))/(2.*a);
	    sh_ = p_ - tau*d_;
	    flag = false;
	  }
	
	  p_ = pp1_;
	  //!part of STEP 4:
	  r_ -= alpha*Bd_;
	  
	  //! STEP 5:
	  //!preconditioning
	  rTilde_ = r_;
	  if(C != SYMMTX()){
	    F77Type<TypeTraits<value_type>::Value>::POTRS(&uplo,&n,&nrhs,&Cfullmtx_[0],&ld,&rTilde_[0],&ld,&info);  //again rTilde_ is overwritten
	  }
	  else
	    rTilde_ = PrecondType::prec(B,rTilde_);

	  rrtil1 = dot(r_,rTilde_);
	  bek = rrtil1/rrtil;
	  d_ = rTilde_ + bek*d_;
	  rrtil = rrtil1;
	} //end if gamma > 0
	else{  //!direction of negative curvature, i.e. \f$ \gamma = (d,Bd) < 0 \f$
     
	  d_ = static_cast<value_type>(Sgn(dot(g0_,d_)))*d_;
	  //normalize(d_);
	  c = dot(p_,p_) - ntimes<2>(delta);
	  b = -2.*dot(p_,d_);
	  a = dot(d_,d_);
	  det_ = ntimes<2>(b) - 4.*a*c;
	  if(det_ < value_type()){ //beware of negative roots
	    ////! 1st alternative
	    //ADONIS_ERROR(DerivedError, "det_ < 0; only noncomplex roots allowed here");

	    ////! 2nd alternative
	    //std::cout << "a = "<< a << "  b = "<< b << "  c = "<< c << std::endl;
	    //ADONIS_INFO(Information,"det is negative ==> complex solution");
	    std::complex<value_type> z = det_;
	    z = std::sqrt(z);
	    z = (-b + z)/(2*a);
	    //std::cout << "POS. COMPLEX ROOT: "<< z << std::endl;
	    tau = Abs(z); //z.real();
	    p_ -= tau*d_;
	    break;
	  }
	  //if(det_ < value_type()) det_ = value_type(); //just make it zero
	  tau = (-b + std::sqrt(det_))/(2*a);
	  p_ -= tau*d_;
	  //!CG may become unreliable in case of neg. curvature
	  std::cout << "NEGATIVE CURVATURE found in Steihaug (P)CG in "<<iter_<<". iteration"<< std::endl;
	  //ADONIS_INFO(Information,"Negative curvature direction in Steihaug CG-iteration i = " << iter_<<".");  //return;
	  break;
	}

	// //safeguard -- don't perform too many CG iterations
	if(iter_ >= maxit){
	  ADONIS_INFO(Information,"Too many CG-iterations (iter_ = "<< iter_ << ") \n   Accept current iterate as best one");
	  break;
	}

      }//end while
      res_ = p_;
      
      if(!flag){
	resize_me_when_empty(n,res_);
	res_ = (delta/NormType::norm(p_,C))*p_;
	pred1 = -(dot(gradf,res_) + 0.5*dot(res_,symm_matrix_vector_multiplication(B,res_)));
	pred2 = -(dot(gradf,sh_) + 0.5*dot(sh_,symm_matrix_vector_multiplication(B,sh_)));
	if(pred2 > pred1)
	  res_ = sh_;
      }
      return res_;
    }

  private:
    value_type det_;
    VType Cfullmtx_;
    VType p_, pp1_, g0_, r_, rTilde_, d_, sh_, Bd_, res_;
    mutable size_t iter_;
  };

}

#endif
