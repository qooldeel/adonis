#ifndef XTENDED_DAVIS_SKODJE_WITH_CPA_HH
#define XTENDED_DAVIS_SKODJE_WITH_CPA_HH

#include<iostream>

#include "../../misc/misctmps.hh"
#include "weirdtransport.hh"
#include "../../io/readinparameters.hh"
#include "../constants.hh"

#include "../../common/typeadapter.hh"

namespace Adonis{

  /**
   * \brief Taken from <I> Ren, Pope, &ldquo The use of slow manifolds in reactive flows &rdquo, Combust. Flame 147 (2006),</I>
   * the 2nd approximation is implemented by using eqs. (37) and (43) on page 255, thus neglecting the specialized formula (44) on the subsequent page.
   * 
   *  Actually it does the same as <TT>withcorrection.hh</TT>. In addition, it forms \f$ z_2(x,t) = z^M(x,t) + \delta z\f$ which can be accessed via &ldquo get_z2() &rdquo .
   */
  template<class T>
  class RPSecApproxCPA{
  public:
    typedef ExprTmpl::MyVec<T> VType;
    typedef T value_type;
    
  private:
    VType z2_;
    
    T a_,
      b_,
      c_,
      d_,
      e_,
      D1_,
      eps_;
 
    T delta_u_;

  public:
   
    typedef typename TypeAdapter<T>::BaseType BaseType; 
    typedef Constant<BaseType> ConstantType;

    enum{
      space = ConstantType::spacePoints            //number of space points
    };


    RPSecApproxCPA(unsigned dim = 0):rhs_(dim), z2_((dim == 0) ? 0 : space),delta_u_(T()){
      Parameter<BaseType,7> P("Parameter.dat");
       a_ = P[0]; b_ = P[1]; c_ = P[2]; d_ = P[3];
       e_ = P[4]; D1_ = P[5]; eps_ = P[6];
    
       std::cout << "a = "<<a_<<", b = "<<b_<<", c = "<<c_<<", d = "<<d_<<", e = "<<e_<<", D1 = "<<D1_ << ", eps = "<< eps_ << std::endl;
    }

    unsigned dim() const {return rhs_.size();}
    unsigned domain_dim() const{return space;} //one species(i.e. = #spacepoints)

    const VType& get_z2() const {return z2_;}

    T& z2_index(unsigned i)  {
      adonis_assert(i < z2_.size());
      return z2_[i];
    }

    T sim(const T& z, const T& a){
      return z/(1 + a*z);         //!the SIM value for z_2
    }

    const T& operator[](unsigned i) const{   //return z_1
      adonis_assert(i < rhs_.size());
      return rhs_[i];
    }

    const T& get_delta_u() const {return delta_u_;} 


    template<class X>
    VType& operator()(const X& z){
      //adonis_assert((int)z.size() == space);  //only one species
      //adonis_assert((int)z2_.size() == space);

       const T h = 1./(space-1);

       DiffusiveTransport<T> D2(D1_,e_); //fields consist of const refs!

       //left bdy                                  //right bdy
       rhs_[0] = 0.;             rhs_[space-1] = 0.;

       //fill z2 at the boundary
       z2_[0] = sim(z[0],a_);    z2_[space-1] = sim(z[space-1],a_);
		     

        T fprime = T(),
	  diffuse1 = T(),
	  x = T();
	
	for(int i = 1; i < space-1; ++i){ //loop over domain (without) bdy
	  fprime = 1./(ntimes<2>(1 + a_*z[i]));
	  diffuse1 = D1_*(z[i+1]- 2.*z[i] + z[i-1])/(ntimes<2>(h));
	  
	  x = i*h;  //equidistant spacing
	  
	  delta_u_ =  eps_/(c_*fprime + 1)*( d_*z[i]*fprime -  z[i]/(ntimes<2>(1+b_*z[i])) - fprime*diffuse1 + ( D2(x+0.5*h)*(z[i+1]/(1.+a_*z[i+1]) - z[i]/(1.+a_*z[i])) - D2(x-0.5*h)*(z[i]/(1.+a_*z[i]) - z[i-1]/(1.+a_*z[i-1])) )/(ntimes<2>(h)) );

	  rhs_[i] = -d_*z[i] + diffuse1 +  //1st approx
	    c_/eps_*              //B^T·J·U  (cf. (37) second line of equation)
	    delta_u_;
	
	  //! Note: SIM is given by \f$ [z_1, z_1/(1 + az_1)]^T\f$
	  //! Note: \f$ z(x,t) = z^M(x,t) + \delta z \f$
	  z2_[i] = sim(rhs_[i],a_) + delta_u_;
	  
	}

	return rhs_;
    }

  private:
    VType rhs_;

  };
  
}

#endif 
