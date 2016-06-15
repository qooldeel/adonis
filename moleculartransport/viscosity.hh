#ifndef DYNAMIC_SHEAR_VISCOSITY_HH
#define DYNAMIC_SHEAR_VISCOSITY_HH

#include "generaltransport.hh"
#include "fancytransporttmps.hh"
#include "../templatemetaprograms/unrollloop.hh"
#include "../templatemetaprograms/matrixunroller.hh"
#include "../common/globalfunctions.hh"

namespace Adonis{

  template<class TPTDATATYPE, bool BOOL>
  class PureSpeciesViscosity:  public GeneralTransport<TPTDATATYPE,PureSpeciesViscosity<TPTDATATYPE,BOOL> >{
  public:
    typedef GeneralTransport<TPTDATATYPE,PureSpeciesViscosity<TPTDATATYPE,BOOL> > BaseClassType;
    typedef typename BaseClassType::DataType DataType;

    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::PointerType PointerType;
    
    enum{NSPEC = DataType::nspec,
	 DIM = SystemDimension<DataType::nspec,DataType::rednspec,BOOL>::Value
    };

    PureSpeciesViscosity(PointerType ptr = 0){
      BaseClassType::ptr_ = ptr;
    }
    
    PureSpeciesViscosity(const PureSpeciesViscosity& Pure){
      BaseClassType::ptr_ = Pure.ptr_;
    }
    
    PureSpeciesViscosity& operator=(const PureSpeciesViscosity& Pure){
      if(this != &Pure){
	BaseClassType::ptr_ = Pure.ptr_;
      }
      return *this;
    }

    void initialize(PointerType ptr){
      BaseClassType::ptr_ = ptr;
    }
    
    PointerType coefficients(){return &vis_[0];}

    value_type& get_viscosity(size_t k){
      adonis_assert(k<DIM && k >= 0);
      return vis_[k];
    }

    const value_type& get_viscosity(size_t k) const{
      adonis_assert(k<DIM && k >= 0);
      return vis_[k];
    }

    template<bool CHS, class T>
    T viscosity_coefficient(size_t k, const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<".");
#endif
      return ChooseRightPrefactor<T,CHS>::visprefac()*sqrt(temp*ChooseRightPrefactor<T,CHS>::in_gram(BaseClassType::molar_mass(k)))/(ntimes<2>(BaseClassType::sigma(k))*BaseClassType::Omega_2_2_star_approximation(k,k,BaseClassType::reduced_temperature(k,temp)));
    }

    template<class T>
    void compute_viscosities(const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<".");
#endif
      ComputeTransportProperties<BOOL>::compute_viscosities(*this,temp);
   }

  private:
    value_type vis_[DIM]; //in kg/(m·s)
  };
  
  

  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  template<class TPTDATATYPE, bool BOOL>
  class MixtureAveragedViscosity: public GeneralTransport<TPTDATATYPE,MixtureAveragedViscosity<TPTDATATYPE,BOOL> >{
  public:
    typedef GeneralTransport<TPTDATATYPE,MixtureAveragedViscosity<TPTDATATYPE,BOOL> > BaseClassType;
    typedef typename BaseClassType::DataType DataType;

    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::PointerType PointerType;
    
    //! pure species properties
    typedef PureSpeciesViscosity<TPTDATATYPE,BOOL> PureViscosityType;

    enum{DIM = SystemDimension<DataType::nspec,DataType::rednspec,BOOL>::Value
    };

    
   MixtureAveragedViscosity(PointerType ptr = 0):PureVis_(ptr),mixVis_(value_type()){
      BaseClassType::ptr_ = ptr;
    }

    
    MixtureAveragedViscosity(const MixtureAveragedViscosity& mav):PureVis_(mav.PureVis_),mixVis_(mav.mixVis_){
      BaseClassType::ptr_ = mav.ptr_;
    }
    
    MixtureAveragedViscosity& operator=(const MixtureAveragedViscosity& mav){
      if(this != &mav){
	PureVis_ = mav.PureVis_;
	mixVis_ = mav.mixVis_;
	BaseClassType::ptr_ = mav.ptr_;
      }
      return *this;
    }
    
    void initialize(PointerType ptr){
      BaseClassType::ptr_ = ptr;
      PureVis_.initialize(ptr);
      mixVis_ = value_type();
    }

 
    template<class T, class X>
    value_type& compute_mixture_averaged_viscosity(const T& temp, const X& Xfrac){
      PureVis_.compute_viscosities(temp);

      return (mixVis_ = MatrixUnroller<0,0,DIM,DIM>::template compute_mixture_averaged_viscosity<DataType,BOOL>(Xfrac,PureVis_.coefficients()));
    }
    

    //! get pure viscosity coefficients 
    PointerType coefficients() {return PureVis_.coefficients();}
    
    //! get mixture-averaged diffusion coefficient
    value_type& viscosity_coeff() {return mixVis_;}
    const value_type& viscosity_coeff() const {return mixVis_;}

  private:
    PureViscosityType PureVis_;
    value_type mixVis_;
  };

} //end namespace 



#endif
