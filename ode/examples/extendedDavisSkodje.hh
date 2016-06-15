#ifndef EXTENDED_DAVIS_SKODJE_HH
#define EXTENDED_DAVIS_SKODJE_HH

#include "../../misc/misctmps.hh"
#include "../../io/readinparameters.hh"
#include "../constants.hh"


namespace Adonis{
  
  /**
     \brief the source of the extendend Davis-Skodje model, which can be found in Ren,Pope, Combust. Flame 147, (2006)
     * Unlike the other implementation of this source code, the Parameter file
     * are passed by a file called "Parameter.dat" in the corresponding pwd.
   */
  template<class T>
  class RPDavisSkodje{
  public:
   
    typedef ExprTmpl::MyVec<T> VType;
    typedef T value_type;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef Constant<BaseType> ConstantType;

    RPDavisSkodje(unsigned dim = 0):rhs_(dim),jac_(2*2){
      ParameterData PD;
      PD.read_from_file("inputfiles/renpopeParams.dat");
      a_ = PD.get_datum<BaseType>("a");
      b_ = PD.get_datum<BaseType>("b");
      c_ = PD.get_datum<BaseType>("c");
      d_ = PD.get_datum<BaseType>("d");
      eps_ = PD.get_datum<BaseType>("eps");
    }
  
    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const{return 2;} 

    void parameter(){
      std::cout << "a = "<<a_ << "  b = "<< b_ << "  c = "<<c_ << "  d = "<<d_ << "  eps = "<< eps_ << std::endl;
    }

    template<class X>
    VType& operator()(const X& z){
      T fprime = (*this).fprime(z);

      rhs_[0] = c_/eps_*(z[1] - fprime) - d_*z[0];
      rhs_[1] = -1./eps_*(z[1] - fprime) - z[0]/(ntimes<2>(1+b_*z[0]));

      return rhs_;
    }

    //! see [1, eq. 30 on page 251]
    template<class X>
    VType& exact_SIM(const X& z){
      rhs_[0] = z[0];
      rhs_[1] = (*this).fprime(z);
      return rhs_;
    }

    template<class X>
    VType& exact_Jacobian(const X& z){
      T tmp = (ntimes<-2>(1+a_*z[0]));
      jac_[0] = -c_/eps_*tmp - d_;                   
      jac_[1] = c_/eps_;
      jac_[2] = 1./eps_*tmp - (1-b_*z[0])/ntimes<3>(1+b_*z[0]);     
      jac_[3] = -1./eps_;
  
      return jac_; 
    }

    template<class X>
    VType& dz1(const X& z){
      value_type tmp = 1/ntimes<2>(1+a_*z[0]);
      rhs_[0] = -c_/eps_*tmp -d_;
      rhs_[1] = 1/eps_*tmp - (1-b_*z[0])/ntimes<3>(1+b_*z[0]);
      return rhs_;
    }

    template<class X>
    VType& tangent_space(const X& z){
      rhs_[0] = 1;
      rhs_[1] = 1./ntimes<2>(1+a_*z[0]);
      return rhs_;
    }
    
  private:
    BaseType a_, b_, c_, d_, eps_;
    VType rhs_, jac_;

    template<class X>
    T fprime(const X& z){
      return ( z[0]/(1 + a_*z[0]) );
    }
  };


}

#endif
