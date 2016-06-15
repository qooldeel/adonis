#ifndef H2C6_1ST_APPROX_PURE_SOURCE_HH
#define H2C6_1ST_APPROX_PURE_SOURCE_HH

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iomanip>

#include "../../../expressiontemplates/exprvec.hh"
#include "../../constants.hh"
#include "../../../io/readinparameters.hh"
#include "../../../modred/manifold.hh"

#include "../../../common/typeadapter.hh"
#include "../../../common/smartassign.hh"


namespace Adonis{

  template<class T>
  class H2MechPURESOURCE1d1STAPPROX{
    public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    
    typedef typename TypeAdapter<T>::Type DType;
    
    typedef ComputeManifold<DType,6,2,true> ManifoldType;
    typedef typename ManifoldType::CompositionType CompositionType;
    typedef typename ManifoldType::VecD VecD; //std::vector<double>
    typedef typename ManifoldType::VecS VecS; //std::vector<std::string> 
    typedef typename ManifoldType::CoordinateType CoordinateType;
   

    H2MechPURESOURCE1d1STAPPROX(unsigned dim = 0, const DType& rpv1 = 0.3, const DType& rpv2 = 0.6):rhs_(dim), rpvh2_(rpv1), rpvh2o_(rpv2),reduced_(2){

      VType pvvalues(2);
      pvvalues[0] = rpv1;
      pvvalues[1] = rpv2;
      
      VecS names(2);
      names[0] = "H2";  
      names[1] = "H2O";
      
      CoordinateType rpvindex;
      rpvindex[0] = 0; rpvindex[1] = 4;

      Mani_.activate(names,pvvalues,rpvindex); //first(...) invoked
      
      std::cout << "RED -- FIRST (RHS) = " << pvvalues <<std::endl;
      std::cout << "SIM -- FIRST (RHS) = " << Mani_.get_z() << std::endl;
    }


    unsigned dim() const {return rhs_.size();}
    unsigned domain_dim() const{return 2;}


    template<class X>
    VType& operator()(const X& x){
      //=================================================
      std::vector<T> rV(12);

      T pk[12]; //reaction rates 
      pk[0] = 2.0;
      pk[1] = 1.0;
      pk[2] = 1.0;
      pk[3] = 1000.0;
      pk[4] = 1000.0;
      pk[5] = 100.0;
      pk[6] = 216.0;
      pk[7] = 337.5;
      pk[8] = 1400.0;
      pk[9] = 10800.0;
      pk[10] = 33750.0;
      pk[11] = 0.7714;
      

      //-------calculate new point             //only if calculated before!
      smart_assign(reduced_[0], x[0]);       //x[0];
      smart_assign(reduced_[1], x[1]);      //x[1];
      Mani_.evaluate(reduced_);   //evaluate new point here

      //get computed full composition from reduced 
      zM_ = Mani_.get_z(); 

      //std::cout << "RED(in) (RHS) = "; print_all(reduced_,16);
      //TypeDependentOutput<T>::show_me(zM_, "SIM(out) (RHS)=");
      
     

      

      rV[0] =   pk[0] * x[0];             //zM_[0];          
      rV[1]  =  pk[6] * zM_[1]*zM_[1];
      
      rV[4]  =  pk[2] * x[1];       //zM_[4];             
      rV[5]  =  pk[8] * zM_[1]*zM_[5];
      rV[6]  =  pk[3] * x[0]*zM_[3];       // zM_[0]
      rV[7]  =  pk[9] * zM_[1]*zM_[5];
      
      rV[10] =  pk[5] * x[0]*zM_[3];       // zM_[0] 
      rV[11] =  pk[11] * x[1];      // zM_[4]; 
      
     
      rhs_[0] = - rV[0]  + rV[1]              //H2
	  - rV[6]  + rV[7]
	  - rV[10] + rV[11];
      	
	
      rhs_[1] = - rV[4]  + rV[5]               //H2O
	  + rV[10] - rV[11];


      
      
      //---------------------------------------------------------------------

      std::cout << "mass conservation (must be const. throughout!) = "<< 2.*x[0] + zM_[1] + 2.*x[1] + zM_[5] << std::endl;
      //=================================================
      return rhs_;
      
    }


    //fake time-dependent operator -- needed (currently) for adaptive methods)
     template<class TIME, class FV>
     VType& operator()(const TIME& time, const FV& z){   
       rhs_ = (*this).operator()(z);
       return rhs_;
     }

  private:
    VType rhs_;
    T rpvh2_, rpvh2o_;
    ManifoldType Mani_;
    VecD reduced_;
    CompositionType zM_; //thats the full composition having 6 species
  };

} //end namespace 

#endif
