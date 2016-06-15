#ifndef MAKE_SYMMETRIC_MATRIX_POSITIVE_DEFINITE_HH
#define MAKE_SYMMETRIC_MATRIX_POSITIVE_DEFINITE_HH

#include <algorithm>

#include "../common/globalfunctions.hh"
#include "../common/elementaryoperations.hh"
#include "../common/adonisassert.hh"
#include "../lapackwrapper/lapackdrivertraits.hh"
#include "../lapackwrapper/wrappertraits.hh"
#include "../misc/useful.hh"

namespace Adonis{

  /**
   * \brief Cholesky factorization with added multiple of identity, i.e.
   * perform  Cholesky for \f$ A \leftarrow A + \tau I\f$, cf. [1].
   *
   * Refernces:
   *  [1] [Nocedal and Wright, "Numerical Optimization", 2nd ed., 2006, Algo. 3.3, p. 51]
   *
   * NOTE: no effords have been made to base the whole stuff on symmetric matrix
   *       objects since we use LAPACK here anyway
   */
  template<class V>
  class CholeskyWithAddedMultipleOfIdentity{
  public:
    typedef typename V::value_type value_type;
    
    CholeskyWithAddedMultipleOfIdentity(V& symtx):symmMtx_(symtx),tau_(value_type()),count_(1){}

    const size_t number_of_iterations() const {return count_;} 

    //give back A <-- A + tau·I
    const V& give_back_matrix() const {return res_;}

    /**it might be beneficial, in case of many iterations, to increase prefac e.g. by say a factor of 10. I wouldn't call it a rule of thumb, but it seems to me that the nearer to 0 the diagonal elements are, the greater prefac should be chosen  
     * 
     * \param beta, prefact Parameters affecting tau_ which is used to update the diagonal
     * \param uplo Give back Cholesky decomposition in 'L'ower or 'U'pper triangular form
     * \param mxit Number of maximum allowed regularization steps
     * \param giveBackMtx if 'true' field res_ is overwritten to contain \f$ A + \tau I\f$; default: false
     * \return Reference pointing to the Cholesky decomposition.
     *
     * HINT: the return Cholesky factorization, say \f$ LL^T\f$ can be used later on to solve a linear system via F77Type<TypeTraits<value_type>::Value>::POTRF(·)
     */
    V& compute(const value_type& beta = 0.001, const value_type& prefac = 2, char uplo = 'L', size_t mxit = 10, bool giveBackMtx = false){
      adonis_assert(beta > value_type() && prefac > 1);
      adonis_assert(uplo == 'L' || uplo == 'U');
      int n = static_cast<int>(std::sqrt(std::distance(symmMtx_.begin(),symmMtx_.end()))),
	lda = std::max(1,n),
	info = 1;        //define it > 0!
      
      //std::cout << "symmMtx_ = " << symmMtx_ << std::endl;

      value_type mindiag = symmMtx_[0];
      for(int i = 1; i < n; ++i){
	mindiag = std::min(symmMtx_[RowMajor::offset(i,i,n)],mindiag);
      }

      //std::cout << "min_i{a_ii} = "<<mindiag<< std::endl;

      if(mindiag > value_type())
	tau_ = value_type();
      else
	tau_ = -mindiag + beta;

      count_ = 1;  //reset
      while(info != 0){
	A_ = symmMtx_; //assign original matrix
	//std::cout << "tau_ = "<< tau_ << std::endl;      //OUTPUT
	update_diagonal<AddBasicElements>(A_,n,tau_); //update with new tau
	//std::cout << "A_ + tau_·I = "<< A_ << std::endl; //OUTPUT
	
	if(giveBackMtx)
	  res_ = A_;

	//!perform Choleky factorization on symmMtx_;
	F77Type<TypeTraits<value_type>::Value>::POTRF(&uplo,&n,&A_[0],&lda,&info);
	tau_ = std::max(prefac*tau_,beta);
	count_++;
	if(count_ == mxit)
	  ADONIS_ERROR(DerivedError,"Too many iterations. \n   Increase argument 'prefac' or number of iterations (less efficient), pal!");
      }
      return A_;  //return successful Cholesky decomposition
    }


    //! when it is desired to give back \f$ A \f$ instead of its Cholesky 
    //! decomposition, you may find this routine useful in conjuction with 
    //! setting the 5th argument of compute to 'true'
    const V& updated_matrix() const {return res_;}
    

  private:
    V& symmMtx_;
    V A_, res_;
    value_type tau_;
    size_t count_;
  };

}

#endif
