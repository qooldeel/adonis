#ifndef TRUST_REGION_APPROACH_HH
#define TRUST_REGION_APPROACH_HH

#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"
#include "../common/typeadapter.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "steihaug.hh"
#include "hessianapproximations.hh"

namespace Adonis{


  /**
   * \brief update trust region radius \f$ \Delta_k\f$
   *
   * This algorithm has been  adapted from [1, algo 4.1, p. 69] 
   *  I have implemented the computation of delta and nu at the \f$k\f$-
   *  iteration. This class can be simply used with a TR algorithm. It 
   *  returns the (new) TR radius \f$ \Delta \f$ and updates a Boolean flag 
   *  <TT> nu_ </TT>.  Due to the value of that falg, we alter the direction, 
   *  say \f$ x \leftarrow x +p \f$,  or leave it unchanged in further 
   *  computations. See also [2, algo (5.1.6), p. 96]
   *  
   * References:
   *  
   *  [1] [NOCEDAL & WRIGHT, "Numerical Optimization", 2nd ed., Springer, 2006]
   *  [2] [FLETCHER, "Practical Methods of Optimization", 2nd ed., Wiley, 1987]
   */
  template<class T, char NT = '2'>
  class TRRadius{
  public:
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef Norm<NT,BaseType> NormType;

    TRRadius(const T& d0 = 1.e-04, const T& dmax = 0.45, const T& gamma = T()):delta0_(d0),deltaMax_(dmax),gamma_(gamma),nu_(false){
      check(d0,dmax,gamma);
    }
      
    bool is_direction_changed() const {return nu_;}
    

    const T& get_delta0() const {return delta0_;}

    //! overwrites fields 
    void initialize(const T& d0, const T& dmax, const T& gamma){ 
      check(d0,dmax,gamma);
      delta0_ = d0;
      deltaMax_ = dmax;
      gamma_ = gamma;
      nu_ = false;
    }

    //! update radius (at the kth itertation)
    template<class V>
    T& update(T& delta, const T& ratio, const V& p) const {
      if(ratio < 0.25)
	delta *= 0.25;
      else{
	if((ratio > 0.75) && (NormType::norm(p) == delta))
	  delta = std::min(2.*delta,deltaMax_);
	//the else-case leaves delta unchanged
      }
    
      if(ratio > gamma_){
	nu_ = true;      //!\f$ x \leftarrow x + p\f$
	//std::cout << "DIRECTION MAY BE CHANGED"<<std::endl;
      }
      else 
	nu_ = false;     //! \f$ x \f$ isn't altered
      //std::cout << "nu_ = "<< nu_ << std::endl;
      return delta;
    }

  private:
    T delta0_, deltaMax_, gamma_;
    mutable bool nu_;

    void check(const T& d0, const T& dmax, const T& gamma){
      adonis_assert((T() < d0) && (d0 < dmax));
      adonis_assert(dmax > T());
      adonis_assert((T()<= gamma) && (gamma < 0.25)); 
    }
    
  };
  

  /**
   * \brief Trust region algorithm for the globalization of unconstrained 
   *  optimization problems of the form \f[\min_x\{f(x)\} \f] by using quadratic
   *  model function, \f$ m_k(p) \approx f(x_k) + \langle \nabla f(x_k),p\rangle + 0.5\cdot \langlep,Bp \rangle\f$  of \f$ f\f$
   *
   * Unlike line search methods, TR methods do <I>not</I> require positive 
   * definiteness of the involved Hessians, \f$B\f$,
   * (or approximations thereof). 
   * 
   *
   * \tparam TRSOLVER algorithm for the solution of \f[ \min_p \{m_k(p) | \|p\| \leq \Delta_k \} \f]
   * \tparam NT norm type for trust region ('1','2' and 'i'/'I' are at our disposal)
   * \tparam HAX update formula for Hessian ('b' = BFGS or 'd' = damped BFGS) 
   *
   * References:
   *
   *  [1] [NOCEDAL and WRIGHT, "Numerical Optimization", 2nd edition, Springer 2006]
   *  [2] [STEIHAUG, "The conjugate gradient method and trust regions in large scale optimization", SIAM J. Numer. Anal., 20(3), pp. 626-637, 1983]
   */
  template<class TRSOLVER, char NT = '2', char HAX = 'b'>
  class TrustRegion{
  public:
    typedef typename TRSOLVER::value_type value_type;
    typedef typename TRSOLVER::SymmetricMatrixType SymmetricMatrixType;
    typedef typename TRSOLVER::VType VType;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef std::string StringType;
    typedef TRRadius<value_type,NT> RadiusType;
    typedef typename RadiusType::NormType NormType;
    
    typedef UpdateHessianApproximation<HAX,VType> HessUpdateType;

    TrustRegion(const value_type& delta0 = 0.005, const value_type& deltaMax = 0.45, const value_type& gamma = 0.15, const BaseType& tol = 0.001):count_(0),radius_(delta0,deltaMax,gamma),tol_(tol),nmg0_(BaseType()),delta_(delta0),ratio_(BaseType()){}

 
    size_t number_of_iterations() const {return count_;}

    template<class FUN, class DERIV>
    VType& method(const VType& x0, const SymmetricMatrixType& B0, const SymmetricMatrixType& Cprec, FUN& fun, DERIV& dx, size_t maxit = 75){
      count_ = 0; //reset counter
      size_t dim = static_cast<size_t>(x0.size());
      TRSOLVER approxModel(dim);
      HessUpdateType hessUp(dim);
 
      VType Bp(dim);
      BaseType mk;

      //value_type eps = std::sqrt(std::numeric_limits<value_type>::epsilon());

      x_ = x0;
      B_ = B0;
      gradf_ = dx.Jacobian(x_);
      nmg0_ = NormType::norm(gradf_);
      
      //REPEAT UNTIL SOME STOPPING CRITERION IS FULFILLED 
      while(NormType::norm(gradf_) > tol_*std::max(1.,nmg0_)){
	count_++;
	//! solve \f$ \min_p \{m(p) | \|p\|\leq \Delta\}\f$ approximately
	//!NOTE: The choice of the tolerance (the 4th argument, tol, of the CG 
	//!method) is crucial insofar as keeping the cost of the TR low. 
	//! It may be chosen akin to that in [1, algo 7.1, 3rd line, p. 169], 
	//!i.e. \f$ \min(0.5,\sqrt{\|\nabla f_k\|}\f$ to obtain superlin. cvgce.
	p_ = approxModel.approximation(B_,gradf_,delta_,std::min(0.5,std::sqrt(NormType::norm(gradf_))),Cprec);

	std::cout << "# of iterations of STEIHAUG PCG: "<< approxModel.number_of_iterations() << std::endl;
	

	//! calculate ratio
	Bp = symm_matrix_vector_multiplication(B_,p_);
	mk = dot(gradf_,p_) + 0.5*dot(p_,Bp);
	//mk = fun(x_) +  dot(gradf_,p_) + 0.5*dot(p_,Bp); //due to [1, eq. (4.4), p. 69]
	ratio_ = (fun(x_) - fun(x_ + p_))/(-mk);  //! due to [2, eq. (3.1). p. 631]
	//ratio_ = (fun(x_) - fun(x_ + p_))/(fun(x_)- mk);; //due to [1, eq. (4.4), p. 69]
	
	//! update TR radius \f$ \Delta \f$
	radius_.update(delta_,ratio_,p_); //delta is updated (or not)
	
	
	
	//! update x_ (or not)
	if(radius_.is_direction_changed() == true){
	  x_ += p_;
	  
 
	  //!only compute derivatives only when we direction is changed
#if USE_APPROX_HESSIAN
	  gradprev_ = gradf_;
#endif 	  

	  gradf_ = dx.Jacobian(x_); //! calculate new gradient

#if USE_APPROX_HESSIAN
	  //! update hessian approximation -- note that \f$ p = x_{k+1}-x_k \f$
	  hessUp.hessian_update(B_,p_,gradf_,gradprev_); //B_ is changed 
#else
	  full_ = dx.Hessian(x_,0);   //calculate new exact Hessian here
	  full_matrix_2_symmetric_matrix(B_,full_);

#endif
	
	  //std::cout << "B_ =" <<  B_ << std::endl;

	}

	//std::cout << "CURRENT iterate x_{"<<count_<<"} = "<< x_ << std::endl;
	
	
	
	if(count_ >= maxit){
	  ADONIS_WARNING(Warning, "Too many iterations in TR algorithm! \n   Taking current iterate as best approximation...");
	  break;
	}

      } //end while
     
      return x_; //approximate solution of 
    }


  private:
    size_t count_;
    RadiusType radius_;
    BaseType tol_, nmg0_, delta_, ratio_;
    VType x_, gradf_, p_, full_;
    VType gradprev_; //! used for possible Hessian updates
    SymmetricMatrixType B_;
  };

}

#endif
