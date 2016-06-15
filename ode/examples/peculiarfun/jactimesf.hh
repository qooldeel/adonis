#ifndef JACOBIAN_TIMES_FUN_HH
#define JACOBIAN_TIMES_FUN_HH

#include <iostream>

#include "../../../expressiontemplates/exprvec.hh"
#include "../../../derivatives/jacobianofsource.hh"
#include "../../../misc/operations4randomaccesscontainers.hh"

namespace Adonis{

  template<class T>
  class SimpleToyMech{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType; 
    
    SimpleToyMech(size_t n = 0):rhs_(n){
      k1_ = 1.;
      ki1_ = 1.e-5;
      k2_ = 10;
      ki2_ = 1.e-5;

    }

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    template<class X>
    VType& operator()(const X& x){
      
      rhs_[0] = -k1_*x[0]*x[0] + ki1_*x[1];
      rhs_[1] =  k1_*x[0]*x[0] -(ki1_ + k2_)*x[1] + ki2_*x[2];
      rhs_[2] = k2_*x[1] - ki2_*x[2];

      return rhs_;
    }

    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
    }


  private:
    VType rhs_;
    T k1_, ki1_, k2_, ki2_;
  };


  //Peculiar function
  template<typename T>
  class JacTimesFun{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType; 
    
    // ============= ORIGINAL SOURCE EMPLOYED ====================
    typedef SimpleToyMech<T> FunType;
    typedef JacS<T,SimpleToyMech,ExprTmpl::MyVec> JacobianType;
    // ===========================================================

    JacTimesFun(size_t n = 0):rhs_(n),fun_(n){
      if(n == 3)
	NablaChem_.set(n,n);
    }


    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    
    template<class X>
    VType operator()(const X& x){
      // NablaChem_.jacobian(x);
      //std::cout << "Jacobian = "<< NablaChem_.get_jacobian() << std::endl;
      rhs_ = matrix_vector_product(NablaChem_.jacobian(x),(*this).dim(),fun_(x));
      return rhs_;
    }


    //! fake operator
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
    }


  private:
    VType rhs_;
    FunType fun_;  
    JacobianType NablaChem_;
  };

}

#endif
