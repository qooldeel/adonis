#ifndef PURE_SPECIES_AND_MIXTURE_AVERAGED_THERMAL_CONDUCTIVITIES_HH
#define PURE_SPECIES_AND_MIXTURE_AVERAGED_THERMAL_CONDUCTIVITIES_HH

#include "../common/globalfunctions.hh"
#include "../common/smartassign.hh"

#include "../massactionkinetics/physicalconstants.hh"
#include "../massactionkinetics/thermochemistry.hh"

#include "generaltransport.hh"
#include "transportsettings.hh"
#include "fancytransporttmps.hh"
#include "prefactors.hh"
#include "diffusion.hh"
#include "viscosity.hh"

#include "../misc/misctmps.hh"

namespace Adonis{
 
  template<class TPTDATATYPE, bool BOOL>
  class PureSpeciesConductivity: GeneralTransport<TPTDATATYPE,PureSpeciesConductivity<TPTDATATYPE,BOOL> >{
  public:
    typedef GeneralTransport<TPTDATATYPE,PureSpeciesConductivity<TPTDATATYPE,BOOL> > BaseClassType;
    typedef typename BaseClassType::DataType DataType;

    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::PointerType PointerType;
    typedef typename BaseClassType::SizeType SizeType;

    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;

    typedef PureSpeciesViscosity<TPTDATATYPE,BOOL> ViscosityType;
    typedef PureSpeciesDiffusion<TPTDATATYPE,BOOL> BinDiffusionType;

    enum{DIM = SystemDimension<DataType::nspec,DataType::rednspec,BOOL>::Value};

    PureSpeciesConductivity(PointerType ptr = 0):nasa_(DataType::nspec,DataType::thermo(),DataType::temperature_bounds()), Vis_(ptr), BinDiff_(ptr),twoDevPi_(2./UniversalConstants<value_type>::Pi){
      BaseClassType::ptr_ = ptr; 
    }
    
    PureSpeciesConductivity(const PureSpeciesConductivity& Pure):nasa_(DataType::nspec,DataType::thermo(),DataType::temperature_bounds()),Vis_(Pure.Vis_),BinDiff_(Pure.BinDiff_),twoDevPi_(2./UniversalConstants<value_type>::Pi){
      BaseClassType::ptr_ = Pure.ptr_;
    }
    
    PureSpeciesConductivity& operator=(const PureSpeciesConductivity& Pure){
      if(this != &Pure){
	nasa_.initialize(DataType::nspec,DataType::thermo(),DataType::temperature_bounds());
	Vis_ = Pure.Vis_;
	BinDiff_ = Pure.BinDiff_;
	twoDevPi_ = 2./UniversalConstants<value_type>::Pi;
	BaseClassType::ptr_ = Pure.ptr_;
      }
      return *this;
    }
    
    void initialize(PointerType ptr){
      nasa_.initialize(DataType::nspec,DataType::thermo(),DataType::temperature_bounds());
      Vis_.initialize(ptr);
      BinDiff_.initialize(ptr);
      twoDevPi_ = 2./UniversalConstants<value_type>::Pi;
      BaseClassType::ptr_ = ptr;
    }

    PointerType coefficients(){return &conduct_[0];}

    value_type& get_conductivity(size_t k){
      adonis_assert(k<DIM && k >= 0);
      return conduct_[k];
    }

    const value_type& get_conductivity(size_t k) const{
      adonis_assert(k<DIM && k >= 0);
      return conduct_[k];
    }


    template<bool CHS, class RHO, class P, class T>
    T conductivity_coefficient(SizeType k, const RHO& rho, const P& p, const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()) || (!is_well_defined_value(p)) || (!is_well_defined_value(rho)) || (rho < RHO()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<" p = "<< p << "  rho = "<< rho << ".");
#endif
      T Zrot = (*this).Z_rot(k,temp),
	fvib = (*this).f_vib(k,rho,p,temp),
	AdevB = (*this).A_dev_B(k,Zrot,fvib);
      return ( Vis_.template viscosity_coefficient<CHS>(k,temp)/DataType::molar_masses()[k]*((*this).f_trans(k,AdevB)*(*this).C_v_trans_dev_R(k)*PhysicalConstants<value_type>::IdealGasConstant + (*this).f_rot(fvib,AdevB)*(*this).C_v_rot_dev_R(k)*PhysicalConstants<value_type>::IdealGasConstant + fvib*(*this).C_v_vib(k,temp)) );
    }


    template<class RHO, class P, class T>
    void compute_conductivities(const RHO& rho, const P& p, const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()) || (!is_well_defined_value(p)) || (!is_well_defined_value(rho)) || (rho < RHO()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<" p = "<< p << "  rho = "<< rho << ".");
#endif
      ComputeTransportProperties<BOOL>::compute_conductivities(*this,rho,p,temp);
    }
    

  private:
    ThermoType nasa_;
    ViscosityType Vis_;
    BinDiffusionType BinDiff_;
    value_type twoDevPi_;         //value that is needed again and again
    value_type conduct_[DIM];

    //! don't let anyone use these member functions outside the class
    //! translational part of \f$C_V\f$ is in all cases the same
    value_type C_v_trans_dev_R(SizeType k) { return 1.5;}

    value_type C_v_rot_dev_R(SizeType k) {
      value_type res = 0.;
      if(BaseClassType::geometry(k) == 0) //single atom
	res = 0.;
      else if(BaseClassType::geometry(k) == 1) //linear molecule
	res = 1.;
      else if(BaseClassType::geometry(k) == 2)  //nonlinear molecule
	res = 1.5;
      return res;
    }

    //! molar heat capacity \f$C_V\f$ of species \f$k\f$ at const. volume which
    //! is given by \f$ C_{V,k} = C_{p,k} - R.\f$
    //! See, [Chemkin-Pro, "Theory Manual", 2008, Eq. 2-27, p. 28]
    template<class T>
    value_type C_V(SizeType k, const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp << ".");
#endif
      return convert_number(nasa_.heat_capacity(k,temp)) - PhysicalConstants<value_type>::IdealGasConstant;
    }

    
    template<class T>
    value_type C_v_vib(SizeType k, const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp << ".");
#endif
      value_type res = 0.;
      if(BaseClassType::geometry(k) == 0) //single atom
	res = 0.;
      else if(BaseClassType::geometry(k) == 1) //linear molecule
	res = (*this).C_V(k,temp) - 2.5*PhysicalConstants<value_type>::IdealGasConstant;
      else if(BaseClassType::geometry(k) == 2)  //nonlinear molecule
	res = (*this).C_V(k,temp) - 3.*PhysicalConstants<value_type>::IdealGasConstant;
      else{} //default
      return res;
    }

    
    template<class T1, class T2, class T3>
    T3 f_vib(SizeType k, const T1& rho, const T2& p, const T3& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T3()) || (!is_well_defined_value(p)) || (!is_well_defined_value(rho)) || (rho < T1()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<" p = "<< p << "  rho = "<< rho << ".");
#endif
      //! incorporates "self-diffusion"
      return ( rho*BinDiff_.template binary_diffusion_coefficient<TransportSettings::molarMassInGrams>(k,k,p,temp)/(Vis_.template viscosity_coefficient<TransportSettings::molarMassInGrams>(k,temp)) );
    }

    template<class T>
    T F_depend(SizeType k, const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp << ".");
#endif
       value_type piPowers1point5 = pow(UniversalConstants<value_type>::Pi,1.5);
       T ekbdevtemp = BaseClassType::epsilon_kB(k)/temp;
       return ( 1 + 0.5*piPowers1point5*sqrt(ekbdevtemp) + (ntimes<2>(UniversalConstants<value_type>::Pi)/4. + 2)*ekbdevtemp + piPowers1point5*pow(ekbdevtemp,1.5) );
    }

    template<class T>
    T Z_rot(SizeType k, const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp << ".");
#endif
      return ( BaseClassType::Z_rot(k)*(*this).F_depend(k,298.)/(*this).F_depend(k,temp) );
    }


    template<class T>
    T A_dev_B(SizeType k, const T& Zrot, const T& fvib){      
      return ( (2.5 - fvib)/(Zrot + twoDevPi_*(5./3.*(*this).C_v_rot_dev_R(k)+fvib)) );
    }
    
    template<class T>
    T f_trans(SizeType k, const T& AdevB){
      return 2.5*(1. - twoDevPi_*(*this).C_v_rot_dev_R(k)/(*this).C_v_trans_dev_R(k)*AdevB);
    }
    
    template<class T1, class T2>
    T1 f_rot(const T1& fvib, const T2& AdevB){
      return fvib*(1. + twoDevPi_*AdevB);
    }

  }; //end class
  


  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  template<class TPTDATATYPE, bool BOOL>
  class MixtureAveragedConductivity: public GeneralTransport<TPTDATATYPE,MixtureAveragedConductivity<TPTDATATYPE,BOOL> >{
  public:
    typedef GeneralTransport<TPTDATATYPE,MixtureAveragedConductivity<TPTDATATYPE,BOOL> > BaseClassType;
    typedef typename BaseClassType::DataType DataType;

    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::PointerType PointerType;

    typedef PureSpeciesConductivity<TPTDATATYPE,BOOL> PureConductivityType;
    
    enum{DIM = SystemDimension<DataType::nspec,DataType::rednspec,BOOL>::Value};

    MixtureAveragedConductivity(PointerType ptr = 0):PureCond_(ptr),mixCond_(value_type()){
      BaseClassType::ptr_ = ptr;
    }

    MixtureAveragedConductivity(const MixtureAveragedConductivity& mac):PureCond_(mac.PureCond_),mixCond_(mac.mixCond_){
      BaseClassType::ptr_ = mac.ptr_;
    }

    MixtureAveragedConductivity& operator=(const MixtureAveragedConductivity& mac){
      if(this != &mac){
	PureCond_ = mac.PureCond_;
	mixCond_ = mac.mixCond_;
	BaseClassType::ptr_ = mac.ptr_;
      }
      return *this;
    }
    
    void initialize(PointerType ptr){
      BaseClassType::ptr_ = ptr;
      PureCond_.initialize(ptr);
      mixCond_ = value_type();
    }

    //! computes mix. avg. thermal conductivity 
    //! \f[ \lambda = \frac{1}{2}\left(\sum_{k=1}^K X_k\lamda_k + \frac{1}{\sum_{k=1}^K X_k/\lambda_k}\right)\f]
    template<class RHO, class P, class T, class X>
    value_type& compute_mixture_averaged_conductivity(const RHO& rho, const P& p, const T& temp, const X& Xfrac){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()) || (!is_well_defined_value(p)) || (!is_well_defined_value(rho)) || (rho < RHO()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<" p = "<< p << "  rho = "<< rho << ".");
#endif
      //! compute full or reduced conductivities and store them in PureCond_
      PureCond_.compute_conductivities(rho,p,temp);
      
      //UnrollLoop<0,DIM>::print(PureCond_.coefficients());

      return (mixCond_ =  convert_number(0.5*(UnrollLoop<0,DIM>::sum_prod(Xfrac,PureCond_.coefficients()) + 1./(UnrollLoop<0,DIM>::template sum_product_x<1,-1>(Xfrac,PureCond_.coefficients())))));
     

      //return (mixCond_ =  convert_number(0.5*(UnrollLoop<0,DIM>::sum_product(Xfrac.begin(),PureCond_.coefficients()) + 1./(UnrollLoop<0,DIM>::template sum_product_xtended<1,-1>(Xfrac.begin(),PureCond_.coefficients())))));
     
    }

  private:
    PureConductivityType PureCond_;
    value_type mixCond_;

  };
  

} //end namespace 

#endif
