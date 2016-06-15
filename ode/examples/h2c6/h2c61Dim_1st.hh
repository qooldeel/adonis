#ifndef FIRST_APPROX_H2C6_1DIM_WITH_TRANSPORT_HH
#define FIRST_APPROX_H2C6_1DIM_WITH_TRANSPORT_HH

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

#include "../../../misc/operations4randomaccesscontainers.hh"


namespace Adonis{

  template<class T>
  class H2MechWithDiffusion1d1STAPPROX{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    
    typedef typename TypeAdapter<T>::Type DType;
    
    typedef ComputeManifold<DType,6,2,true> ManifoldType;
    typedef typename ManifoldType::CompositionType CompositionType;
    typedef typename ManifoldType::VecD VecD; //std::vector<double>
    typedef typename ManifoldType::VecS VecS; //std::vector<std::string> 
    typedef typename ManifoldType::CoordinateType CoordinateType;
   

  private:
    T D1_,D5_,  //the diffusion coefficients
      rpvh2o_, rpvh2_;

    ManifoldType Mani_;
    std::ofstream of_;

    template<class U>
    T fd1D2nd(const U& u, int i, const T& h){ //const diffusion
      return ( (u[i+1] - 2.*u[i] + u[i-1])/(h*h) );
    }
    
   
    
  public:
    
    enum{ 
      space = Constant<double>::spacePoints,
      numOfSpecies = 2,                      //only 2 now !!
      domainDim = numOfSpecies*space
    };

    
    H2MechWithDiffusion1d1STAPPROX(unsigned dim = 0, const DType& rpvh2 = 0.3, const DType& rpvh2o = 0.6, const DType& tol = 1.e-03):rhs_(dim), rpvh2_(rpvh2), rpvh2o_(rpvh2o),tol_(tol),reduced_(numOfSpecies){
      
      Parameter<double,6> Dc("Diffcoefs.dat");
      D1_ = Dc[0]; D5_ = Dc[4];
      

      
      
      VType pvvalues(numOfSpecies);
      pvvalues[0] = rpvh2;
      pvvalues[1] = rpvh2o;
      
      VecS names(numOfSpecies);
      names[0] = "H2";  
      names[1] = "H2O";
      
      CoordinateType rpvindex;
      rpvindex[0] = 0; rpvindex[1] = 4;

      Mani_.activate(names,pvvalues,rpvindex); //first(...) invoked
      

      //std::cout << "z_sim = "<< Mani_.get_z() << std::endl;
      //std::string simname = "SIM_x_t.dat"; 
      //int isys;
      //isys = system(("rm " + simname).c_str()); //because we append 
      //of_.open(simname.c_str(), std::ios_base::out | std::ios_base::app);
      //of_ << "# Sim as computed by Jochen's tool at the discretisation points" << std::endl;
    }

    ~H2MechWithDiffusion1d1STAPPROX(){
      //of_.close();
    }

    unsigned dim() const {return rhs_.size();}
    unsigned domain_dim() const{return domainDim;}

    /////////////////////////////////////////////////////
    //operator()(·) -- here, warmstart is invoked
    /////////////////////////////////////////////////////
    template<class X>
    VType& operator()(const X& red){
      
      /*   
      T norm = 0.0;
      for(size_t i = 0; i < domain_dim(); ++i){
	norm += ntimes<2>(Abs(red[i])); 
      }

      norm = sqrt(norm);

      std::cout << "norm (input) = "<<norm << std::endl;
      */
      
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
      
      const T h = 1./(space-1);

      //bdy constant? Sim characterized via rpv (no time- and space-dpcy)
      //no flux across bdy
      for(int i = 0; i < numOfSpecies; ++i){
	rhs_[i*space] = 0.;                           //left
	rhs_[(i+1)*(space-1)] = 0.;                   //right
      }
     
     

      //  diff_ = DType();  //initialize [0,0,...,0]

      //interior
      for(int i = 1; i < space-1; ++i){
	
	

	//if(Norm<'2',DType>::norm(diff_) > tol_){
	  //diff_ = zM_;

	//! compute only at selected points 
	//if(i == 1 || i == (space-2)/2 || i == space-2){
	//if(i==1 || i%2 != 0 || i == space-2){  //!selected points //evaluate only at one half of the original spatial points to increase performance
	//if(i == 1){
	smart_assign(reduced_[0],red[i]);        
	smart_assign(reduced_[1],red[space+i]);    
	Mani_.evaluate(reduced_);   //evaluate new point here
	zM_ = Mani_.get_z();  //very first entry is from intialisation
	  
	  //}

	 
	//note that here are species in dependence of others, only z_1 and z_5 occur (represented by corresponding forward and backward reactions)
	rV[0] =   pk[0] * red[i];             //zM[0];          
	rV[1]  =  pk[6] * zM_[1]*zM_[1];

	rV[4]  =  pk[2] * red[space+i];       //zM[4];             
	rV[5]  =  pk[8] * zM_[1]*zM_[5];
	rV[6]  =  pk[3] * red[i]*zM_[3];       // zM[0]
	rV[7]  =  pk[9] * zM_[1]*zM_[5];

	rV[10] =  pk[5] * red[i]*zM_[3];       // zM[0] 
	rV[11] =  pk[11] * red[space+i];      // zM[4]; 


	// ========= Approximation of 2n derivs via FDs ==================
	T diffH2 = D1_*fd1D2nd(red,i,h),
	  diffH2O =  D5_*fd1D2nd(red,space+i,h);

	//================================================================

	//with diffusion -- the reaction progress variables
	rhs_[i] = - rV[0]  + rV[1]              //H2
	  - rV[6]  + rV[7]
	  - rV[10] + rV[11]   + diffH2;
	
	
	//std::cout <<" D1·ð^2_x(H2) = " << diffH2 << std::endl;
	
	
	rhs_[space+i] = - rV[4]  + rV[5]               //H2O
	  + rV[10] - rV[11]     + diffH2O;
      
	//std::cout <<" D5·ð^2_x(H2O) = " << diffH2O << std::endl;

	//std::cout << "MASS BALANCE (spatial) = "<< 2.*red[i] + zM[1] + 2.*red[space+i] + zM[5] << std::endl;

      } //end spatial discretisation

      return rhs_;
    }
    
    //fake time-dependent operator
     template<class TIME, class FV>
     VType& operator()(const TIME& time, const FV& z){   
       rhs_ = (*this).operator()(z);
       return rhs_;
     }

  private:
    VType rhs_;
    CompositionType zM_,   //thats the full composition having 6 species
      diff_;
    VecD reduced_;
  
    DType tol_;
  };

}

#endif 
