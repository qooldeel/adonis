#ifndef H2C6_1DIM_WITH_TRANSPORT_HH
#define H2C6_1DIM_WITH_TRANSPORT_HH


#include "../../../expressiontemplates/exprvec.hh"
#include "../../constants.hh"
#include "../../../io/readinparameters.hh"


namespace Adonis{

  template<class T>
  class H2MechWithDiffusion1D{
  private:
    T D1_,D2_,D3_,D4_,D5_,D6_,  //the diffusion coefficients
      rpvh2o_, rpvh2_;

    template<class U>
    T fd1D2nd(const U& u, int i, const T& h){ //const diffusion
      return ( (u[i+1] - 2.*u[i] + u[i-1])/(h*h) );
    }

  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
 
    enum {
      space = Constant<double>::spacePoints,
      numOfSpecies = 6,
      domainDim = numOfSpecies*space
    };


    H2MechWithDiffusion1D(unsigned dim = 0):rhs_(dim){
      Parameter<double,6> Dc("Diffcoefs.dat");
      D1_ = Dc[0]; D2_ = Dc[1];  D3_ = Dc[2]; D4_ = Dc[3]; 
      D5_ = Dc[4]; D6_ = Dc[5];
    }

    unsigned dim() const {return rhs_.size();}
    unsigned domain_dim() const{return domainDim;}

    template<class X>
    VType& operator()(const X& z){
     
      T rV[12];

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
      //left bdy                                  //right bdy
      //rhs_[0] = 0.;             rhs_[space-1] = 0.; 
      //rhs_[space] = 0.;         rhs_[2*space-1] = 0.;
      //rhs_[2*space] = 0.;       rhs_[3*space-1] = 0.;
      //rhs_[3*space] = 0.;       rhs_[4*space-1] = 0.;
      //rhs_[4*space] = 0.;       rhs_[5*space-1] = 0.;
      //rhs_[5*space] = 0.;       rhs_[6*space-1] = 0.;

      //interior
      for(int i = 1; i < space-1; ++i){
	rV[0] =   pk[0] * z[i];
	rV[1]  =  pk[6] * z[space+i]*z[space+i];
	rV[2]  =  pk[1] * z[2*space+i];
	rV[3]  =  pk[7] * z[3*space+i]*z[3*space+i];
	rV[4]  =  pk[2] * z[4*space+i];
	rV[5]  =  pk[8] * z[space+i]*z[5*space+i];
	rV[6]  =  pk[3] * z[i]*z[3*space+i];
	rV[7]  =  pk[9] * z[space+i]*z[5*space+i];
	rV[8]  =  pk[4] * z[2*space+i]*z[space+i];
	rV[9]  =  pk[10] * z[3*space+i]*z[5*space+i];
	rV[10] =  pk[5] * z[i]*z[3*space+i];
	rV[11] =  pk[11] * z[4*space+i];

	//with diffusion
	rhs_[i] = - rV[0]  + rV[1]              //H2
	  - rV[6]  + rV[7]
	  - rV[10] + rV[11]   + D1_*fd1D2nd(z,i,h);

	rhs_[space+i] = 2.0*rV[0] - 2.0*rV[1]          //H
	  + rV[4] - rV[5]
	  + rV[6] - rV[7]
	  - rV[8] + rV[9]     + D2_*fd1D2nd(z,space+i,h);
      
	rhs_[2*space+i] = - rV[2] + rV[3]                //O2                 
	  - rV[8] + rV[9]   + D3_*fd1D2nd(z,2*space+i,h);
	
	rhs_[3*space+i] = + 2.0*rV[2] - 2.0*rV[3]        //O
	  - rV[6]     + rV[7]
	  + rV[8]     - rV[9]
	  - rV[10]    + rV[11]  + D4_*fd1D2nd(z,3*space+i,h);; 
	
	rhs_[4*space+i] = - rV[4]  + rV[5]               //H2O
	  + rV[10] - rV[11]     + D5_*fd1D2nd(z,4*space+i,h);
      
	rhs_[5*space+i] = rV[4] - rV[5]                  //OH
	  + rV[6] - rV[7]
	  + rV[8] - rV[9]   + D6_*fd1D2nd(z,5*space+i,h);
      }

      return rhs_;
    }
    

     template<class TIME, class FV>
     VType& operator()(const TIME& time, const FV& z){   
       rhs_ = (*this).operator()(z);
       return rhs_;
     }

  private:
    VType rhs_;

  };

}

#endif 
