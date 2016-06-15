#ifndef EXACT_TRUST_REGION_SUB_PROBLEM_HH
#define EXACT_TRUST_REGION_SUB_PROBLEM_HH

#include <iostream>
#include <cmath>

#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"
#include "../common/typeadapter.hh"
#include "../common/elementaryoperations.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../misc/useful.hh"
#include "../misc/misctmps.hh"
#include "../lapackwrapper/lapackdrivertraits.hh"
#include "../lapackwrapper/wrappertraits.hh"
#include "../linalg/linearsystemsolvers.hh"

namespace Adonis{
  
  /**
   *   \brief This algorithm is similar to that found in [1, §4, algo 4.3, p.87] 
   *    which calculates the direction \f$p(\lambda), \ \lambda > 0\f$, such 
   *    such that \f$ B+\lambda I\f$ is p.d. and \f$ \lambda \f$ is updated
   *    with the trust region radius \f$ \Delta \f$
   *
   * [1] [NOCEDAL/WRIGHT, "Numerical Analysis", 2nd ed., 2006]
  */
  template<class HESS, class GRAD, char NT = '2'>
  class ExactSubTRNewton{
  public:
    typedef typename HESS::value_type value_type;
    typedef GRAD VType;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;

    typedef Norm<NT,BaseType> NormType;

    ExactSubTRNewton(size_t n = 0, const value_type& lambda = value_type()):lambda_(lambda),p_(n),q_(n){}
    
    
    void initialize(size_t n = 0, const value_type& lambda = value_type()){
      lambda_ = lambda_;
      p_.resize(n);
      q_.resize(n);
    }

    VType& compute_direction(const value_type& lambda0, const HESS& B, const GRAD& g, const value_type& delta, const value_type& tol = 1.e-05, size_t maxit = 10, char uplo = 'L'){

      int n = std::max(1,static_cast<int>(g.size())),
	nrhs = 1,
	info = 1;  //false
      
      VType eigvals(n), tmpHess, R(n*n);
      // VType RT(n*n);

      hessFull_ = symmetric_matrix_2_full_matrix(hessFull_,B);

      tmpHess = hessFull_;

      //std::cout << "Hessian (full) = "; print_in_matrix_style(hessFull_,n);
      
      //! eigenvals are stored in ascending order, i.e. 
      //! \f$ \lambda_1 \leq \lambda_2, ..., \leq \lambda_n \f$
      eig_symmetric_hermitian(tmpHess,n,eigvals);

      
      //std::cout << "eigen values of B: "<< eigvals << std::endl;
      lambda_ = lambda0;
      
      //safeguard
      // if(lambda_ < -eigvals[0])
      // 	lambda_ = Sgn(eigvals[0])*eigvals[0];
     
      size_t count = 1, inner = 1;
      

      for(size_t l = 0; l < maxit; ++l){
	count++;
	p_ = -g;
		
	//!perform Cholesky decomposition of 
	inner = 1; //reset
	while(info != 0){
	  inner++;
	  A_ = hessFull_;
	  update_diagonal<AddBasicElements>(A_,n,lambda_);	

	  F77Type<TypeTraits<value_type>::Value>::POTRF(&uplo,&n,&A_[0],&n,&info);
	  lambda_ = std::max(10*lambda_,1.e-03);

	  if(inner >= 7){
	    ADONIS_INFO(Information,"Too many regularization steps");
	    break;
	  }
	  
	}
	adonis_assert(info == 0);  //succesful exit

	//copy 
	if(uplo == 'L'){
	  for(int i = 0; i < n; ++i){
	    for(int j = i; j < n; ++j){
	      // R[RowMajor::offset(i,j,n)] = A_[RowMajor::offset(i,j,n)];
	      //! transpose directly
	      R[RowMajor::offset(j,i,n)] = A_[RowMajor::offset(i,j,n)];
	    }
	  }
	}

	F77Type<TypeTraits<value_type>::Value>::POTRS(&uplo,&n,&nrhs,&A_[0],&n,&p_[0],&n,&info);
	adonis_assert(info == 0); //succesful lin. solve

	q_ = p_;
	non_square_solve(R,n,q_);
	
	//RT = transpose(R,n);
	//non_square_solve(RT,n,q_);


	BaseType nmp = NormType::norm(p_),
	  diff = nmp - delta; 
       
	lambda_ += ntimes<2>(nmp/NormType::norm(q_))*(diff/delta);
	//std::cout << "lambda = " << lambda_ << std::endl;

	// if(Abs(diff) <= tol){
	//   break;
	// }

	//!practical versions are content with an approx. solution that can
	//! be obtained within 3 iterations (see [1, p. 87])
	if(count >= 3){ 
	  break;
	}
	
	
	adonis_assert(count <= maxit);
      }//end loop
      return p_;
    }

  private:
    value_type lambda_;
    VType hessFull_, p_,q_, A_;
    
    
  };

}

#endif
