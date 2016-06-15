#ifndef LINEAR_INVARIANT_AKA_MASS_BALANCE_HH
#define LINEAR_INVARIANT_AKA_MASS_BALANCE_HH

#include <iostream>

#include "../common/adonisassert.hh"
#include "../common/error.hh"
#include "../common/typeadapter.hh"
#include "../common/globalfunctions.hh"
#include "../linalg/linearsystemsolvers.hh"
#include "../misc/operations4randomaccesscontainers.hh"

namespace Adonis{


	/**
	 * \brief Linear invariants arise in practice as mass balance invarinats, see e.g. [1,2]
	 *
	 * The algorithm is an adaption of the method discussed in [1].
	 *
	 * References:
	 *  [1] [SHAMPINE, "Conservation Laws and the Numerical Solution of ODEs", Comp. & Maths. with Appls. 128, 1986, pp. 1287-1296]
	 *
	 *  [2] [GEAR, "Invariants and numerical methods for ODEs", Physica D 60, 1992, pp. 303-310]
	 *
	 * \tparam V vector for storage of dense matrix and vector \f$C\f$ and \f$b,\f$ respectively
	 * \tparam NM Norm type (since a linear least squares problem is considered, it is chosen to be '2' by default for the Euclidean norm
	 */
	template<class V, char NM = '2'>
	class MassBalanceInvariant{
	public:
	  typedef std::size_t SizeType;
	  typedef typename V::value_type value_type;
	  typedef typename TypeAdapter<value_type>::BaseType BaseType;
	  typedef Norm<NM,BaseType> NormType;
	  
	  MassBalanceInvariant(const V& c, const V& b):count_(0),C_(c),b_(b){}

	  //! \param maxit maximum number of minimum norm iterations
	  //! \param tol tolerance, i.e. convergence is achieved when \f$ \|y_{n+1} - y_{n+1}^* \| <= tol \f$
	  //! \param u current full space estimate
	  V& enforce_balance(SizeType maxit, const BaseType& tol, const V& u){
	    SizeType cols = C_.size()/b_.size();
	    adonis_assert(cols != 0);

	    int m = static_cast<int>(b_.size()),
	      nrhs = 1;   //always one in this application
	      
	    //!o.k., only resize when necessary!
	    (g_.size() != b_.size()) ? g_.resize(b_.size()) : do_nothing();
	    
	    //start with given full space estimate u
	    uk_ = u; 

	    //! Computes \f$ C \delta = -(Cu - b)\f$ iteratively
	    for(SizeType k = 0; k <= maxit; ++k){
	      count_ = k;
	      //std::cout << k << ".)   ";
	      g_ = -matrix_vector_product(C_,uk_) + b_;
	      delta_ = solve_rank_deficient_ls(C_,m,g_,nrhs);
	      //std::cout << "delta = "<< delta_ << std::endl;
	      uk_ += delta_;
	      
	      if(NormType::norm(delta_) <= tol){
		std::cout << "Mass balance granted after "<<count_ << " iters."<<std::endl;
		break;  //leave loop since convergence has been achieved
	      }
	    }
	    
	    if(count_ == maxit){
	      std::cout << " you get: "<< uk_ << std::endl;
	      //ADONIS_INFO(Information, "Max. number of iterations reached (maxit = "<< maxit << ").\n   No convergence with tol = "<< tol << ".");
	      ADONIS_ERROR(IterationError, "Max. number of iterations reached (maxit = "<< maxit << ").\n   No convergence detected with tol = "<< tol << ".");

	    }
	  
	    return uk_;   //this is a perturbation of u that satisfies linear
	                  //mass balances
	  }


	private:
	  SizeType count_;
	  const V& C_;
	  const V& b_;
	  V g_, delta_, uk_;
	  
	};

} //end namespace 




//// JUST SOME PLAY-AROUND THAT WAS NOT FOUND TO WORK PROPERLY
 /**
	   * \brief Overload function to try this "hidden constraint" variant
	   */
	  // template<class JAC, class FUN, class IXIT>
	  // V& enforce_balance(SizeType maxit, const BaseType& tol, const V& u, FUN& fun, JAC& jac, IXIT indexptr, int nred){
	  //   SizeType cols = C_.size()/b_.size();
	  //   adonis_assert(cols != 0);

	  //   int m = static_cast<int>(b_.size()),
	  //     nrhs = 1;   //always one in this application
	      
	  //   //!o.k., only resize when necessary!
	  //   (g_.size() != b_.size()) ? g_.resize(b_.size()) : do_nothing();
	    
	  //   //start with given full space estimate u
	  //   uk_ = u; 

	  //   V f, mtx;

	  //   //! Computes \f$ C \delta = -(Cu - b)\f$ iteratively
	  //   for(SizeType k = 0; k <= maxit; ++k){
	  //     count_ = k;
	  //     std::cout << k << ".)   ";
	  //     f = fun(uk_);
	  //     mtx = matrix_vector_product(C_,matrix_vector_product(jac.jacobian(uk_),f)); //! mtx is R x 1 "matrix"
	  //     g_ = -matrix_vector_product(C_,f); //! is R x 1 vector
	  //     delta_ = solve_rank_deficient_ls(mtx,m,g_,nrhs);
	  //     std::cout << "delta_reduced = "<< delta_ << std::endl;
	  //     //uk_ += delta_; //ERROR delta_ is (R,1) !!!
	      
	  //     for(int i = 0; i < nred; ++i){
	  // 	uk_[indexptr[i]] = delta_[i];
	  //     }

	  //     if(NormType::norm(delta_) <= tol){
	  // 	break;  //leave loop since convergence has been achieved
	  //     }
	  //   }
	    
	  //   if(count_ == maxit){
	  //     //ADONIS_INFO(Information, "Max. number of iterations reached (maxit = "<< maxit << ").\n   No convergence with tol = "<< tol << ".");
	  //     ADONIS_ERROR(IterationError, "Max. number of iterations reached (maxit = "<< maxit << ").\n   No convergence detected with tol = "<< tol << ".");

	  //   }
	  
	  //   return uk_;   //this is a perturbation of u that satisfies linear
	  //                 //mass balances
	  // }


#endif
