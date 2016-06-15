#ifndef XDS_FULL_NUMERICAL_TREATMENT_HH
#define XDS_FULL_NUMERICAL_TREATMENT_HH

#include "../../../expressiontemplates/exprvec.hh"
#include "../../../misc/misctmps.hh"

#include "../weirdtransport.hh"

#include "../../../io/readinparameters.hh"

#include "../../constants.hh"

namespace Adonis{
  
  /**
   * \brief Full model
   */
  template<class T>
  class XDSFull{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;

    typedef Constant<double> ConstantType;

  private:
    T a_,
      b_,
      c_,
      d_,
      e_,
      D1_,
      eps_;

      VType rhs_;

    template<class U>
    T central_fd_2nd(const U& u, int i, const T& h){
      return ( (u[i+1] - 2.*u[i] + u[i-1])/(h*h) );
    }
      

  public:
    enum{
      numOfSpecies = 2,          //full
      space = ConstantType::spPts,            //number of space points
      domainDim = numOfSpecies*space
    };

    XDSFull(size_t dim=0):rhs_(dim){
      Parameter<double,7> P("Parameter.dat");
      a_ = P[0]; b_ = P[1]; c_ = P[2]; d_ = P[3];
      e_ = P[4]; D1_ = P[5]; eps_ = P[6];

      std::cout << "The parameters are as follows:"<<std::endl;
      std::cout << "a = "<<a_<<", b = "<<b_<<", c = "<<c_<<", d = "<<d_<<", e = "<<e_<<", D1 = "<<D1_ << " and  eps = "<<eps_<<std::endl;
  
    }
    
    
    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const{return domainDim;}

    /**
     * \brief the discretised system in space 
     */
    template<class X>
    VType& operator()(const X& y){

      //==================== SPACE DISCRETISATION ============================
      //adonis_assert((int)y.size() == domainDim);
      
      
      const T h = 1./(space-1);  //interval [0,1]

      DiffusiveTransport<T> D2(D1_,e_); //fields consist of const refs!
      
      //!NOTE: all boundary values are constant, hence \f$\frac{d}{dt} y_{\operatorname{boundary}} = 0.\f$
      //left bdy                                  //right bdy
      rhs_[0] = 0.;             rhs_[space-1] = 0.; 
      rhs_[space] = 0.;         rhs_[2*space-1] = 0.;   

      //interior
      T z1_z2 = T(), x = T();
      for(int i = 1; i < space-1; ++i){
	z1_z2 = y[space+i] - y[i]/(1+a_*y[i]); //z2 - z1/(1+az1)

	x = i*h;

	rhs_[i] = c_/eps_*z1_z2 - d_*y[i] + D1_*central_fd_2nd(y,i,h); 
	  //D1_*(y[i+1] - 2* y[i] + y[i-1])/(ntimes<2>(h));

	rhs_[space+i] = -1./eps_*z1_z2 - y[i]/(ntimes<2>(1. + b_*y[i])) +   ( D2(x+0.5*h)*(y[space+i+1] - y[space+i]) - D2(x-0.5*h)*(y[space+i] - y[space+i-1]) )/(ntimes<2>(h));
      }
     
      //======================================================================
    

      return rhs_;
    }


    //! fake time dependent operator
    template<class V>
    VType& operator()(const T& time, const V& y){
      return (rhs_ = (*this).operator()(y));
    }
    
  };


}

#endif
