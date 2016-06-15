#ifndef BOUNDS_ON_VARIABLES_FOR_MECHANISM_HH
#define BOUNDS_ON_VARIABLES_FOR_MECHANISM_HH

#include <iostream>
#include "../common/universalconstants.hh"
#include "../expressiontemplates/exprvec.hh"

namespace Adonis{

  template<class T, class CLASS> 
  class BoundsOnVariables{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef std::size_t SizeType;

    void initialize(SizeType n){
      n_ = n;
      low_.resize(n);
      up_.resize(n);
    }
    

  private:
    CLASS& refer_2_derived() {return static_cast<CLASS&>(*this);}
    
    VType& low() {return refer_2_derived().low();}
    VType& up() {return refer_2_derived().up();}

  
  protected:
    SizeType n_;
    VType low_,
      up_;
  };


  //======================================================================
  //===================== Derived Classes ================================
  //======================================================================

  /**
   * \brief bounds on variables
   * \tparam T value type
   * \tparam ENCOD integer that encodes the mechanism of interest
   */
  template<class T, int ENCOD> class MechanismBounds: public BoundsOnVariables<T, MechanismBounds<T,ENCOD> >{};

  //! specializations 
  template<class T>  //H2C6 Mechanism
  class MechanismBounds<T,6>: public BoundsOnVariables<T,MechanismBounds<T,6> >{
  public:
    typedef BoundsOnVariables<T,MechanismBounds<T,6> > BaseClassType;
    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::VType VType;
    typedef typename BaseClassType::SizeType SizeType;

    MechanismBounds(SizeType n = 0){
      if(n > 0)
	BaseClassType::initialize(n);
    }
    
    void initialize(SizeType n){
      //adonis_assert(BaseClassType::n_ == 6);
      BaseClassType::initialize(n);
    }

    VType& low() {
      //adonis_assert(BaseClassType::n_ == 6);
      BaseClassType::low_[0] = 0.;
      BaseClassType::low_[1] = 0.;
      BaseClassType::low_[2] = 0.;
      BaseClassType::low_[3] = 0.;
      BaseClassType::low_[4] = 0.;
      BaseClassType::low_[5] = 0.;
      
      return  BaseClassType::low_;
    }

    VType& up() {
      //adonis_assert(BaseClassType::n_ == 6);
      BaseClassType::up_[0] = UniversalConstants<T>::aboutInfinity;
      BaseClassType::up_[1] = UniversalConstants<T>::aboutInfinity;
      BaseClassType::up_[2] = UniversalConstants<T>::aboutInfinity;
      BaseClassType::up_[3] = UniversalConstants<T>::aboutInfinity;
      BaseClassType::up_[4] = UniversalConstants<T>::aboutInfinity;
      BaseClassType::up_[5] = UniversalConstants<T>::aboutInfinity;
    
      return BaseClassType::up_;
    }

  };



} //end namespace 

#endif
