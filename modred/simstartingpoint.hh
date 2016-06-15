#ifndef SIM_STARTING_POINT_DEFAULT_HH
#define SIM_STARTING_POINT_DEFAULT_HH

#include <iostream>
#include "../common/universalconstants.hh"
#include "../expressiontemplates/exprvec.hh"

namespace Adonis{

  template<class T, class CLASS> 
  class GeneralSimStartingPoint{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef std::size_t SizeType;

    void initialize(SizeType n){
      n_ = n;
      simstart_.resize(n);
    }
    

  private:
    CLASS& refer_2_derived() {return static_cast<CLASS&>(*this);}
    
    VType& get_point() {return refer_2_derived().get_point();}

  
  protected:
    SizeType n_;
    VType simstart_;
  };

  //======================================================================
  //===================== Derived Classes ================================
  //======================================================================

  /**
   * \brief SIM default starting point
   * \tparam T value type
   * \tparam ENCOD integer that encodes the mechanism of interest
   */
  template<class T, int ENCOD> class SIMStartingPoint: public GeneralSimStartingPoint<T, SIMStartingPoint<T,ENCOD> >{};


  //! specialization
  //! H2C6 mechanism
  template<class T>
  class SIMStartingPoint<T,6>: public GeneralSimStartingPoint<T, SIMStartingPoint<T,6> >{
  public:
    typedef GeneralSimStartingPoint<T, SIMStartingPoint<T,6> > BaseClassType;
    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::VType VType;
    typedef typename BaseClassType::SizeType SizeType;

    SIMStartingPoint(SizeType n = 0){
      if(n > 0)
	BaseClassType::initialize(n);
    }
    
    void initialize(SizeType n){
      //adonis_assert(BaseClassType::n_ == 6);
      BaseClassType::initialize(n);
    }

    VType& get_point(){
      BaseClassType::simstart_[0] = 3.0004541410314495e-01;  //H2
      BaseClassType::simstart_[1] = 1.8944296646551151e-01;  //H
      BaseClassType::simstart_[2] = 1.4354979042038563e-01;  //O2
      BaseClassType::simstart_[3] = 1.0238618155887168e-01;  //O
      BaseClassType::simstart_[4] = 5.9995196772784154e-01;  //H2O
      BaseClassType::simstart_[5] = 1.0562269872515486e-02;  //OH
   

      // BaseClassType::simstart_[0] = 0.4549999999992859;  //H2
      // BaseClassType::simstart_[1] = 0.7780585738931507;  //H
      // BaseClassType::simstart_[2] = 0.2366143850825262;  //O2
      // BaseClassType::simstart_[3] = 0.3628298037265891;  //O
      // BaseClassType::simstart_[4] = 0.1479999999999196;  //H2O
      // BaseClassType::simstart_[5] = 0.01594142610843904;  //OH
      

      return BaseClassType::simstart_;
    }


  };

} //end namespace

#endif
