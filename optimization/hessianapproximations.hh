#ifndef HESSIAN_APPROXIMATION_FORMULAS_HH
#define HESSIAN_APPROXIMATION_FORMULAS_HH

#include "../common/globalfunctions.hh"
#include "../common/elementaryoperations.hh"
#include "../misc/operations4randomaccesscontainers.hh"

namespace Adonis{

  /**
   * \brief Hessian matrix updates
   *
   * NOTE: we assume that the Hessian is <I>dense</I>. Since the Hessian is 
   *       symmetric/hermitian, we only store the lower triangular part (then 
   *       column-wise) which is equivalent to storing the upper triangular 
   *       part (in row-wise fashion). All operations will be performed with 
   *       that storage format!
   */
  template<class V, class UP>
  class MatrixUpdate{
  private:
    UP& r2d(){return static_cast<UP&>(*this);}
    
    V& update(V& B, const V& xnext, const V& x, const V& gradxnext, const V& gradx){
      return r2d().update_formula(B,xnext,x,gradxnext,gradx);
    }

    V& update(V& B, const V& p, const V& gradxnext, const V& gradx){
      return r2d().update_formula(B,p,gradxnext,gradx);
    }

  protected:
    typedef typename V::value_type value_type;
    V sk_, yk_,
      dy1_,dy2_,
      Bs_;

    //! these are for general usage 
  public:
    V& hessian_update(V& B, const V& xnext, const V& x, const V& gradxnext, const V& gradx){
      return update(B,xnext,x,gradxnext,gradx);
    }
    
    V& hessian_update(V& B, const V& p, const V& gradxnext, const V& gradx){
      return update(B,p,gradxnext,gradx);
    }
    

    void init(size_t dim, const value_type& d = value_type()){
      sk_.resize(dim,d);
      yk_.resize(dim,d);
    }
    
  };

  
  /**
   * \tparam C select Hessian update formula ('b','d' or 'i': BFGS, damped BFGS,   *                                           or inverse Hessian, respectively)    
   * \tparam V array type (STL-compliant) to store symmetric matrix
   */
  template<char C, class V> class UpdateHessianApproximation;

  //!partial specialisations
  /**
   * BFGS formula \f[ B_{k+1} \leftarrow B_k - \frac{B_ks_ks_k^TB_k}{s_k^TB_ks_k}-\frac{y_ky_k^T}{y_k^Ts_k} \quad s_k = x_{k+1} -x_k, y_k = \nabla f_{k+1}- \nabla f_k. \f] Note, that \f$ B_ks_k = s_k^TB_k \f$, since \f$B_k\f$ is symmetric.
   *
   * The BFGS update preserves symmetric positive definiteness, i.e. for every
   * s.p.d. start iterate, say \f$ B_0\f$, we get s.p.d. successors, \f$ B_k, \quad k = 1,2,3,\ldots\f$
   * For convenience: \f$ B_0 = I \f$
   */
  template<class V> class UpdateHessianApproximation<'b',V>: public MatrixUpdate<V,UpdateHessianApproximation<'b',V> >{
  public:
    typedef MatrixUpdate<V,UpdateHessianApproximation<'b',V> > BaseClassType;
    typedef typename BaseClassType::value_type value_type;

    UpdateHessianApproximation(size_t dim = 0){
      BaseClassType::init(dim);
    }

    void initialize(size_t dim){
      BaseClassType::init(dim);
    }

    //! note only lower triangle (column-wise) is stored and computed (resp. upper in row-wise fashion) in array!!
    //! assume that vectors admit the usual vector operations -, +
    //! otherwise, use vector_vector_subtraction(u,v,w), etc.
    //! Formula taken from [1, Eq. (2.19), p. 24]
    V& update_formula(V& B,const V& xnext, const V& x, const V& gradxnext, const V& gradx){
      //! BEWARE of cancellation, i.e. when xnext ~ x
      BaseClassType::sk_ = xnext-x;
      BaseClassType::yk_ = gradxnext - gradx;
      //adonis_assert(!nearly_zero(BaseClassType::sk_) || !nearly_zero(BaseClassType::yk_));

      BaseClassType::Bs_ = value_type(); //reset Bs_ to zero!!
      symm_matrix_vector_multiplication(BaseClassType::Bs_,B,BaseClassType::sk_);
      BaseClassType::dy1_ = dyad(BaseClassType::Bs_)/dot(BaseClassType::sk_,BaseClassType::Bs_);
      BaseClassType::dy2_ = dyad(BaseClassType::yk_)/dot(BaseClassType::yk_,BaseClassType::sk_);
      return ( B = B - BaseClassType::dy1_ + BaseClassType::dy2_ );
    }
  
    //! if \f$ p_k := x_{k+1} + x_k \f$ is already computed
      V& update_formula(V& B,const V& pk, const V& gradxnext, const V& gradx){
	BaseClassType::yk_ = gradxnext - gradx;
	adonis_assert(!nearly_zero(BaseClassType::yk_));
	BaseClassType::Bs_ = value_type(); //reset Bs_ to zero!!
	symm_matrix_vector_multiplication(BaseClassType::Bs_,B,pk);
	BaseClassType::dy1_ = dyad(BaseClassType::Bs_)/dot(pk,BaseClassType::Bs_);
      BaseClassType::dy2_ = dyad(BaseClassType::yk_)/dot(BaseClassType::yk_,pk);
      return ( B = B - BaseClassType::dy1_ + BaseClassType::dy2_ );
      }

  };
  

  //! damped BFGS, cf. [NOCEDAL, WRIGHT, "Numerical Optimization",2nd. ed., 2006, algo 18.2, p. 537]
  template<class V> class UpdateHessianApproximation<'d',V>: public MatrixUpdate<V,UpdateHessianApproximation<'d',V> >{
  public:
    typedef MatrixUpdate<V,UpdateHessianApproximation<'d',V> > BaseClassType;
    typedef typename BaseClassType::value_type value_type;

    UpdateHessianApproximation(size_t dim = 0):theta_(value_type()),sTy_(value_type()),sTBs_(value_type()),r_(dim){
      BaseClassType::init(dim);
    }

    //! note only lower triangle (column-wise) is stored and computed (resp. upper in row-wise fashion) in array!!
    V& update_formula(V& B,const V& xnext, const V& x, const V& gradxnext, const V& gradx){
      BaseClassType::sk_ = xnext-x;
      BaseClassType::yk_ = gradxnext - gradx;
      //adonis_assert(!nearly_zero(BaseClassType::sk_) || !nearly_zero(BaseClassType::yk_));
      BaseClassType::Bs_ = value_type(); //reset Bs_ to zero since Bs_ += !!
      symm_matrix_vector_multiplication(BaseClassType::Bs_,B,BaseClassType::sk_);
      sTy_ = dot(BaseClassType::sk_, BaseClassType::yk_);
      sTBs_ = dot(BaseClassType::sk_,BaseClassType::Bs_);
      
      if(sTy_ >= 0.2*sTBs_)
	theta_ = 1.;
      else
	theta_ = (0.8*sTBs_)/(sTBs_ - sTy_);

      r_ = theta_*BaseClassType::yk_ + (1-theta_)*BaseClassType::Bs_;

      BaseClassType::dy1_ = dyad(BaseClassType::Bs_)/sTBs_;
      BaseClassType::dy2_ = dyad(r_)/dot(BaseClassType::sk_,r_);

      return ( B = B - BaseClassType::dy1_ + BaseClassType::dy2_ );      
    }
  
  private:
    value_type theta_, sTy_, sTBs_;
    V r_;
  };
  

  /**
   * \brief Approximation of inverse Hessian (inverse is also symmetric). 
   * Some practical implementations avoid the need to factorize the Hessian at
   * each iteration by updating the inverse of it. Then we only have to perform
   * matrix-vector operations in each step instead of solving a linear system!
   *
   * It is suggested to use \f$ H_0 = \frac{y^Ts}{y^Ty}I\f$ rather than \f$ H_0 = I \f$ as a starting iterate, see [1, eq. (6.20), p. 143]
   *
   * The implementation of the inverse Hessian approximation follows [1, eq. (2.21), p. 25]
   * 
   * References:
   *
   *  [1] [NOCEDAL and WRIGHT, "Numerical Analysis", 2nd ed., Springer, 2006]
   */
   template<class V> class UpdateHessianApproximation<'i',V>: public MatrixUpdate<V,UpdateHessianApproximation<'i',V> >{
  public:
    typedef MatrixUpdate<V,UpdateHessianApproximation<'i',V> > BaseClassType;
    typedef typename BaseClassType::value_type value_type;

     UpdateHessianApproximation(size_t dim = 0):rho_(value_type()){
      BaseClassType::init(dim);
    }

    //! note only lower triangle (column-wise) is stored and computed (resp. upper in row-wise fashion) in array!!
     V& update_formula(V& H, const V& xnext, const V& x, const V& gradxnext, const V& gradx){
       BaseClassType::sk_ = xnext-x;
       BaseClassType::yk_ = gradxnext - gradx;
       //adonis_assert(!nearly_zero(BaseClassType::sk_) || !nearly_zero(BaseClassType::yk_));
       //BaseClassType::dy1_ = dyad(
       rho_ = 1./dot(BaseClassType::yk_,BaseClassType::sk_);
       BaseClassType::dy1_ = dyad(BaseClassType::sk_,BaseClassType::yk_)*(-rho_);
       //BaseClassType::dy1_ *= -rho_;
       update_diagonal<AddBasicElements>(BaseClassType::dy1_,x.size(),1.);
     
       BaseClassType::dy2_ = transpose(BaseClassType::dy1_,x.size());
       
       symmetric_matrix_2_full_matrix(H_,H);
       
       //reset necessary since += !
       Left_ = value_type();
       Right_ = value_type();

       matrix_matrix_multiplication(Left_,BaseClassType::dy1_,x.size(),H_,x.size());
       matrix_matrix_multiplication(Right_,Left_,x.size(),BaseClassType::dy2_,x.size());
     
       //Right_ stores the symmetric matrix. Now overwrite H appropriately
       full_matrix_2_symmetric_matrix(H,Right_);
       H += rho_*dyad(BaseClassType::sk_);
       return H;
     }
   
   private:
     value_type rho_;
     V Left_, //non-symm
       H_,      //stores symmetric matrix in full format
       Right_;
     };

     
} //end namespace 

#endif
