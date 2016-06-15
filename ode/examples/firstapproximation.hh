#ifndef FIRST_APPROXIMATION_MOL_HH
#define FIRST_APPROXIMATION_MOL_HH

#include "../../misc/misctmps.hh"
#include "../../io/readinparameters.hh"
#include "../constants.hh"


namespace Adonis{
  template<class T>
  class FirstApproximation{
  private:
    T d_,        //only these two are needed
      D1_;
  
  public:

    typedef ExprTmpl::MyVec<T> VType;
    typedef T value_type;
    typedef Constant<double> ConstantType;

    enum{
      space = ConstantType::spacePoints            //number of space points
    };

    
    FirstApproximation(size_t dim=0):rhs_(dim){
      Parameter<double,7> P("Parameter.dat");
      d_  = P[3];
      D1_ = P[5];
    
      std::cout<<" d = "<< d_ <<"  D1 = " << D1_ <<std::endl;
    }
    
    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const{return space;} //one species; hence equals space points


    /**
     *  \brief first approximation. Note that \f$ z_2 = \frac{z_1}{1+az_1}\f$ effects that a part of the chemistry is zero  
     */
    template<class X>
    VType& operator()(const X& z){
      //============== the first approximation =================================  
      //adonis_assert((int)z.size() == space);  //only one species in the reduction

      const T h = 1./(space-1); 
      
      //left bdy                                  //right bdy
      rhs_[0] = 0.;             rhs_[space-1] = 0.; 
     
      for(int i = 1; i < space-1; ++i){
	rhs_[i] = -d_*z[i] + D1_*(z[i+1]- 2.*z[i] + z[i-1])/(ntimes<2>(h));
      }
      
      //========================================================================

      return rhs_;
    }


    //! fake time dependent operator
    template<class V>
    VType& operator()(const T& time, const V& y){
      return (rhs_ = (*this).operator()(y));
    }

  private:
    VType rhs_;

  };

}

#endif 
