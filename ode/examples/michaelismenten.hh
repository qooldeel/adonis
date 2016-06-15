#ifndef MICHAELIS_MENTEN_SYSTEMS_AS_ODES_HH
#define MICHAELIS_MENTEN_SYSTEMS_AS_ODES_HH

#include "../common/adonisassert.hh"
#include "../../misc/misctmps.hh"
#include "../../io/readinparameters.hh"
namespace Adonis{

  /**
   * \brief This example has been taken from [1]
   *
   * References:
   *
   *  [1]  Gear, Kaper, Kevrekidis and Zagaris, "Projecting to a slow manifold: Singularly perturbed systems and legacy codes", SIADS, 4, 2005 
   */
  template<class T>
  class MichaelisMentenOMalley{
  public:
    typedef ExprTmpl::MyVec<T> VType;
    typedef T value_type;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef Constant<BaseType> ConstantType;

    MichaelisMentenOMalley(unsigned n = 0):rhs_(n),jac_(n*n){
         
      ParameterData PD;
      PD.read_from_file("/home/mfein/MARC++/ode/inputfiles/michaelismenten.dat");
      kappa_ = PD.get_datum<T>("kappa"); 
      lambda_ = PD.get_datum<T>("lambda");
      epsilon_ = PD.get_datum<T>("epsilon");
      adonis_assert(epsilon_ > T());
    }

    size_t dim() const {return 2;}
    size_t domain_dim() const{return 2;}

    const T& kappa() const {return kappa_;}
    const T& lambda() const {return lambda_;}
    const T& epsilon() const {return epsilon_;}

    template<class X>
    VType& operator()(const X& z){
      rhs_[0] = -z[0] + (z[0] + kappa_ - lambda_)*z[1];
      rhs_[1] = 1./epsilon_*(z[0] - (z[0] + kappa_)*z[1]);
      return rhs_;
    }

    template<class X>
    VType& exact_SIM(const X& z){
      rhs_[0] = z[0];
      rhs_[1] = z[0]/(z[0]+kappa_);
      return rhs_;
    }

    template<class X>
    VType& exact_Jacobian(const X& z){
      T ieps = 1./epsilon_;
      jac_[0] = -1 + z[1];
      jac_[1] = z[0] + kappa_ - lambda_;
      jac_[2] = ieps*(1 - z[1]);
      jac_[3] = -ieps*(z[0] + kappa_);

      return jac_;
    }

    //! just replace <TT>z[1]</TT> with <TT>h(z[0]) = z[0]/(z[0]+kappa_)</TT>
    template<class X>
    T reduced_dynamics(const X& z){
      return -lambda_*z[0]/(z[0]+kappa_);
    }

    //! use this to check method when computing SIM approximately
    template<class X>
    T redfun(const X& z) const{
      return -z[0] + (z[0] + kappa_ - lambda_)*z[1];
    }

    //!scalar argument version
    T reduced_rhs(const T& z){
      return  -lambda_*z/(z+kappa_);
    }

    template<class X> //! I only cross-checked with ipopt
    T reduced_dynamics_Jacobian(const X& z){
      return -lambda_*kappa_/(ntimes<2>(z[0] + kappa_)); 
      //( -1 + (z[0]*z[0] + 2*kappa_*z[0] + kappa_*(kappa_-lambda_))/ntimes<2>(z[0]+kappa_) );
    }

     //!scalar argument version
    T reduced_rhs_prime(const T& z){
      return -lambda_*kappa_/(ntimes<2>(z + kappa_)); 
	//( -1 + (z*z + 2*kappa_*z + kappa_*(kappa_-lambda_))/ntimes<2>(z+kappa_) );
    }

  private:
    VType rhs_, jac_;
    T kappa_, lambda_, epsilon_;
  };

}

#endif
