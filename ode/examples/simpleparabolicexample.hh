#ifndef LINEAR_PARABOLIC_EXAMPLE_DUE_TO_JOHNSON_NIE_THOMEE_HH
#define LINEAR_PARABOLIC_EXAMPLE_DUE_TO_JOHNSON_NIE_THOMEE_HH

#include "../../expressiontemplates/exprvec.hh"
#include "../../common/typeadapter.hh"
#include "../constants.hh"
#include "../../misc/misctmps.hh"
#include "../../common/adonisassert.hh"
#include "../../common/globalfunctions.hh"

namespace Adonis{
  /**
   * \brief This example can be found in 
   [C.JOHNSON, Y.Y. NIE, V. THOM´EE, "AN A POSTERIORI ERROR ESTIMATE AND ADAPTIVE TIMESTEP CONTROL FOR A BACKWARD EULER DISCRETIZATION OF A PARABOLIC PROBLEM", SIAM, vol. 27, 1990, ep. (4.1), pp. 287--291]
   */
  template<class T>
  class SimpleLinearParabolicProblem{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<T> VType;
    typedef Constant<DType> ConstantsType;

    enum{
      NumberOfSpacePoints = ConstantsType::discretizationPointsInSpace,
      DomDim = NumberOfSpacePoints+2
    };

    SimpleLinearParabolicProblem(unsigned dim = 0):rhs_(dim),h_(1./(NumberOfSpacePoints+1)){}

    unsigned dim() const {return rhs_.size();}
    unsigned domain_dim() const{return DomDim;}


    template<class X>
    VType& operator()(const X& x){
      T h_2 = ntimes<2>(h_); 

      //BCs
      rhs_[0] = 1.; rhs_[NumberOfSpacePoints+1] = 1.;
      
      

      //interior
      for(unsigned i = 1; i < NumberOfSpacePoints+1; ++i){
	rhs_[i] = (x[i+1] - 2.*x[i] + x[i-1])/h_2;
      }

      return rhs_;
    }

    
    //! fake time-dependent operator()(·)
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
    }

  private:
    VType rhs_;
    T h_;
  };


  /**
   * \brief The three initial conditions as stated in the above reference
   */

  template<class VType>
  inline VType& IC_I(VType& psi, const typename VType::value_type& h){
    for(unsigned j = 1; j < psi.size()-1; ++j){  // 0<x<1
      psi[j] = 5.*(1-2*Abs(j*h - 0.5))+1;
    }
    return psi;
  }
   
  template<class VType>
  inline VType& IC_II(VType& psi, const typename VType::value_type& h){
    for(unsigned j = 1; j < psi.size()-1; ++j){ // 0<x<1
      psi[j] = -3.;
    }
    return psi;
  }
  
  template<class VType>
  inline VType& IC_III(VType& psi, const typename VType::value_type& h){
    for(unsigned j = 1; j < psi.size()-1; ++j){   // 0<x<1
      psi[j] = 1./std::sqrt(j*h);
    }
    return psi;
  }
    
  


  //===================================================================
  template<class X>
  inline X source_term_for_1_D_heat_equation(const X& x){
    return DiracDeltaDistribution<'l'>::approximation(x,1.e-04);
  }




  /**
   * 1D heat equation with source \f$ f(x,t) = f(x) = \Delta_{\epsilon}(x)\f$
   */
 template<class T>
  class HeatEquation1D{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<T> VType;
    typedef Constant<DType> ConstantsType;

    enum{
      NumberOfSpacePoints = ConstantsType::discretizationPointsInSpace,
      DomDim = NumberOfSpacePoints+2
    };

    HeatEquation1D(unsigned dim = 0):rhs_(dim),h_(1./(NumberOfSpacePoints+1)){}

    unsigned dim() const {return rhs_.size();}
    unsigned domain_dim() const{return DomDim;}


    template<class X>
    VType& operator()(const X& x){
      T h_2 = ntimes<2>(h_); 

      //BCs
      rhs_[0] = 0.;  rhs_[NumberOfSpacePoints+1] = 0.;
      
      

      //interior
      for(unsigned i = 1; i < NumberOfSpacePoints+1; ++i){

	rhs_[i] = (x[i+1] - 2.*x[i] + x[i-1])/h_2 +  source_term_for_1_D_heat_equation(i*h_);
      
      }

      return rhs_;
    }

    
    //! fake time-dependent operator()(·)
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
    }

  private:
    VType rhs_;
    T h_;
  };


}

#endif 
