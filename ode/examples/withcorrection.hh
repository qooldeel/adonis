#ifndef WITH_CORRECTION_TERMS_HH
#define WITH_CORRECTION_TERMS_HH

#include "../../misc/misctmps.hh"

#include "weirdtransport.hh"

#include "../../io/readinparameters.hh"

#include "../constants.hh"


namespace Adonis{

  template<class T>
  class WithCorrection{
  private:
    T a_,
      b_,
      c_,
      d_,
      e_,
      D1_;    //eps not needed here

  public:
    
    typedef ExprTmpl::MyVec<T> VType;
    typedef T value_type;
    typedef Constant<double> ConstantType;

    enum{
      space = ConstantType::spacePoints            //number of space points
    };

    
    WithCorrection(size_t dim=0):rhs_(dim){
      Parameter<double,7> P("Parameter.dat");
      a_ = P[0]; b_ = P[1]; c_ = P[2]; d_ = P[3];
      e_ = P[4]; D1_ = P[5];
      
      std::cout << "a = "<<a_<<", b = "<<b_<<", c = "<<c_<<", d = "<<d_<<", e = "<<e_<<", D1 = "<<D1_ <<std::endl;
    }

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const{return space;} //one species(i.e. = #spacepoints)

    
    /**
     *\brief Note that \f[ B^TJU = \frac{c}{\epsilon} \f]
     *  Furthermore, we have 
     * \f[ \delta u = \frac{\epsilon}{cf'+1}\cdot \left(df'z_1 - \frac{z_1}{(1+bz_1)^2} - D_1f'\frac{\partial^2 z_1}{\partial x^2} + \frac{partial}{\partial x} D_2 \frac{\partial z_2}{\partial x} \right). \f]
     */
    template<class X>
    VType& operator()(const X& z){
      //=================== 1st approx. + CORRECTION ========================== 
      //adonis_assert((int)z.size() == space);  //only one species

      const T h = 1./(space-1);
     
      DiffusiveTransport<T> D2(D1_,e_); //fields consist of const refs!
      
      //left bdy                                  //right bdy
      rhs_[0] = 0.;             rhs_[space-1] = 0.; 

      T fprime = T(),
	diffuse1 = T(),
	x = T();
     
      //INTERIOR
      for(int i = 1; i < space-1; ++i){
	fprime = 1./(ntimes<2>(1 + a_*z[i]));
	
	diffuse1 = D1_*(z[i+1]- 2.*z[i] + z[i-1])/(ntimes<2>(h));
	
	x = i*h;

	rhs_[i] = -d_*z[i] +  diffuse1 +     //1st approx
	  c_/(c_*fprime + 1)*( d_*z[i]*fprime - z[i]/(ntimes<2>(1+b_*z[i])) - fprime*diffuse1 + ( D2(x+0.5*h)*(z[i+1]/(1.+a_*z[i+1]) - z[i]/(1.+a_*z[i])) - D2(x-0.5*h)*(z[i]/(1.+a_*z[i]) - z[i-1]/(1.+a_*z[i-1])) )/(ntimes<2>(h)) );
      }


      //=======================================================================
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
