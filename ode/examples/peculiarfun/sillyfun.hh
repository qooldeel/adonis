#ifndef JUST_SOME_TEST_ODE_RHS_HH
#define JUST_SOME_TEST_ODE_RHS_HH

#include <string>
#include "../../../expressiontemplates/exprvec.hh"
#include "../../../misc/operations4randomaccesscontainers.hh"
#include "../../../misc/misctmps.hh"
#include "../../../io/readinparameters.hh"

namespace Adonis{

  template<class T>
  class VarEqual{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type Type;
    
    typedef ExprTmpl::MyVec<T> VType; 

    VarEqual(size_t n = 0, const std::string& fname = "TestParameter.dat"):rhs_(n),jac_(n*n){
      ParameterData PD;
      PD.read_from_file(fname);
      a_ = PD.get_datum<Type>("a");
      b_ = PD.get_datum<Type>("b");
      c_ = PD.get_datum<Type>("c");
      d_ = PD.get_datum<Type>("d");
      eps_ = PD.get_datum<Type>("eps");
    }

    Type& get_a() {return a_;}
    Type& get_b() {return b_;}
    Type& get_c() {return c_;}
    Type& get_d() {return d_;}
    Type& get_eps() {return eps_;}

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    template<class X>
    VType& operator()(const X& z){
      T tmp = 1/ntimes<2>(1+a_*z[0]);
      jac_[0] = -c_/eps_*tmp - d_;  jac_[1] = c_/eps_;
      jac_[2] = 1/eps_*tmp - (1-b_*z[0])/ntimes<3>(1+b_*z[0]); jac_[3] = -1/eps_;

      VType matvecprod(dim());
      
      for(size_t i = 0; i < dim(); ++i){
	for(size_t j = 0; j < dim(); ++j){
	  matvecprod[i] += jac_[i*dim()+j]*z[j];
	}
	rhs_[i] = matvecprod[i];
      }
     
      //! add derivatives w.r.t. parameter z_1
      rhs_[0] += jac_[0];
      rhs_[1] += jac_[2];
      

      return rhs_;
    }
    
  private:
    VType rhs_, jac_;
    Type a_,b_,c_,d_,eps_;
  };

}

#endif
