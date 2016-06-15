#ifndef INTERFACE_FOR_SIMPLIFIED_NEWTON_METHOD_HH
#define INTERFACE_FOR_SIMPLIFIED_NEWTON_METHOD_HH

namespace Adonis{

  /**
   * \brief Solves the itertion \f$ G'(y_0)\Delta y_n = -G(y_n), y_{n+1} = y_n + \Delta y_n,\f$ the evaluation of the system Jacobian is only performed once per set of Newton iterations. However, to increase convergence it may be beneficial to compute the Jacobian anew after a few iterations.
   * \tparam JAC Jacobian type of iteration, e.g. Gprime
   * \tparam B boolean to decide whether a simplified Newton shall be used (yes)
   *           or not (false)
   */
  template<class JAC, bool B> class SimplifiedNewton;


  //! simplified Newton method
  template<class JAC>
  class SimplifiedNewton<JAC,true>{
  public:
    typedef typename JAC::VecBaseType VType;
    typedef typename JAC::MatrixType MatrixType;

    SimplifiedNewton(JAC& Gp, unsigned n):Gprime_(Gp),afterIter_(n){}

    template<class X>
    VType& evaluate_jacobian(const X& x, unsigned NewtonIter, bool reassign = true){
      //calculate system Jacobian only at iter = 0 and if multiple of afterIter_
      if((NewtonIter==0) || (NewtonIter%afterIter_ == 0)){
	Jac_ = Gprime_.evaluate_jacobian(x,reassign);
      }
      return Jac_;
    }

    //! must be followed by member evaluate_jacobian
    template<class H, class C>
    MatrixType& compute_G_prime(const H& h, unsigned NewtonIter, const C& a, const C& b){
      //calculate system Jacobian only at iter = 0 and if multiple of afterIter_
      if((NewtonIter==0) || (NewtonIter%afterIter_ == 0)){
	Gp0_ = Gprime_.compute_G_prime(h,a,b); //store once computed G'
      }
      //! this is just for the purpose that we don't lose information about G' computed some steps before
      if(NewtonIter > 0){
	  Gprime_.get_matrix() = Gp0_;  //assign to iteration matrix in case it will be overwritten/modified by a solution procedure
      }
      return Gp0_;
    }

  private:
    JAC& Gprime_;
    unsigned afterIter_;  
    MatrixType Gp0_;
    VType Jac_;
  };

  
  //! Conventional Newton method
  template<class JAC>
  class SimplifiedNewton<JAC,false>{
  public:
    typedef typename JAC::VecBaseType VType;
    typedef typename JAC::MatrixType MatrixType;

    SimplifiedNewton(JAC& Gp, unsigned n):Gprime_(Gp),afterIter_(n){}

    //! 2nd argument has no effect heare
    template<class X>
    VType& evaluate_jacobian(const X& x, unsigned NewtonIter, bool reassign = true){
      return Gprime_.evaluate_jacobian(x,reassign);
    }
   
    //! must be followed by member evaluate_jacobian
    //! 2nd argument has no effect heare
    template<class H, class C>
    MatrixType& compute_G_prime(const H& h, unsigned NewtonIter, const C& a, const C& b){
	return Gprime_.compute_G_prime(h,a,b);
    }

  private:
    JAC& Gprime_;
    unsigned afterIter_;  
  };


  

} //end namespace 

#endif
