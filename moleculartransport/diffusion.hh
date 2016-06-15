#ifndef MOLECULAR_TRANSPORT_DIFFUSION_HH
#define MOLECULAR_TRANSPORT_DIFFUSION_HH

#include "generaltransport.hh"
#include "avoidsingularities.hh"
#include "transportsettings.hh"

#include "../misc/misctmps.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../templatemetaprograms/unrollloop.hh"

#include "../lapackwrapper/lapackdrivertraits.hh"
#include "../lapackwrapper/wrappertraits.hh"
#include "../common/smartassign.hh"
#include "../common/globalfunctions.hh"
#include "../massactionkinetics/physicalconstants.hh"

#include "fancytransporttmps.hh"
#include "../misc/useful.hh"

namespace Adonis{

  /** 
   * \brief stand-alone object for diffusion
   * \tparam TPTDATATYPE transport type. 
   * \tparam BOOL boolean type to switch on/off reduced description
   *         (true = reduced, false = full)
   * \param  ptr Pointer to transport data
   */
  template<class TPTDATATYPE, bool BOOL>
  class PureSpeciesDiffusion: public GeneralTransport<TPTDATATYPE,PureSpeciesDiffusion<TPTDATATYPE,BOOL> >{
  public:
    typedef GeneralTransport<TPTDATATYPE,PureSpeciesDiffusion<TPTDATATYPE,BOOL> > BaseClassType;
    typedef typename BaseClassType::DataType DataType;

    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::PointerType PointerType;
   
    enum{NSPEC = TPTDATATYPE::nspec, //full system dimension
	 DIM = SystemDimension<DataType::nspec,DataType::rednspec,BOOL>::Value
    };
    PureSpeciesDiffusion(PointerType ptr = 0){
      BaseClassType::ptr_ = ptr;
    }
    
    PureSpeciesDiffusion(const PureSpeciesDiffusion& Pure){
      BaseClassType::ptr_ = Pure.ptr_;
    }
    
    PureSpeciesDiffusion& operator=(const PureSpeciesDiffusion& Pure){
      if(this != &Pure){
	BaseClassType::ptr_ = Pure.ptr_;
      }
      return *this;
    }

    void initialize(PointerType ptr){
      BaseClassType::ptr_ = ptr;
    }
    
    //! here offset denotes the index associated with storage format, i.e.
    //! when storage by rows is used, take offset = I*N - (I+1)*I/2 + J
    value_type& get_binary_diffusion_coefficient(int offset) {return BinDiff_[offset];}
    const value_type& get_binary_diffusion_coefficient(int offset) const {return BinDiff_[offset];}


    PointerType get_matrix_ptr() {return &BinDiff_[0];}
    
    //!returns pointer to crucial coefficients, i.e. here to BinDiff_
    PointerType coefficients(){return &BinDiff_[0];}

    //!binary diffusion coefficient. 
    //!USAGE: <TT> .template binary_diffusion_coefficient<true>(i,j,p,temp)</TT>
    //! true for g/mol and pressure in bar scaling
    template<bool CHS, class T1, class T2>
    T2 binary_diffusion_coefficient(size_t i, size_t j, const T1& p, const T2& temp) {
#ifndef NDEBUG
      if(!is_well_defined_value(temp) || (temp < T2()) || (!is_well_defined_value(p)))
	ADONIS_ERROR(ValueError,"Bad vals: temperature = "<< temp << "   p = "<< p << ".");
#endif
      //return ChooseRightPrefactor<T2,CHS>::diffprefac()*sqrt(0.5*1./(ChooseRightPrefactor<T2,CHS>::in_gram(BaseClassType::reduced_molar_mass(i,j)))*ntimes<3>(temp))/(p*ntimes<2>(BaseClassType::sigma(i,j))*BaseClassType::Omega_1_1_star_approximation(i,j,BaseClassType::reduced_temperature(i,j,temp)));
      return BaseClassType::template  bin_diff_coeff<CHS>(i,j,p,temp);
    }

    template<class T1, class T2>
    void compute_binary_diffusion_matrix(const T1& p, const T2& temperature){
#ifndef NDEBUG
      if(!is_well_defined_value(temperature) || (temperature < T2()) || (!is_well_defined_value(p)))
	ADONIS_ERROR(ValueError,"Bad vals: temperature = "<< temperature << "   p = "<< p << ".");
#endif
      ComputeTransportProperties<BOOL>::compute_binary_diffusion_matrix(*this,p,temperature);
    }

  private:
    //! I decided to store it as value_type which might be a non-CPPAD type
    value_type BinDiff_[DIM*(DIM+1)/2]; //! store symmetric binary diff. matrix
  };



  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  // · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG · AVG ·
  /**
   * \brief Mix. avg. diffusion. Unlike the multicomponent approximation, mixture-averaged diffusion is based on approximation. Hence it is less accurate and, 
   as stated in [KEE], does not fulfill mass conservation in pure species situations. Yet, it is much cheaper to computate than the mc. description.
   *
   * CAUTION: for single species mixtures (DIM=1) use the self-diffusion coefficient because the mixture-averaged isn't defined in that case!!
  */
  template<class TPTDATATYPE, bool BOOL>
  class MixtureAveragedDiffusion: public GeneralTransport<TPTDATATYPE,MixtureAveragedDiffusion<TPTDATATYPE,BOOL> >{
  public:
    typedef GeneralTransport<TPTDATATYPE,MixtureAveragedDiffusion<TPTDATATYPE,BOOL> > BaseClassType;
    typedef typename BaseClassType::DataType DataType;
    
    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::PointerType PointerType;
    
    //! pure species properties
    typedef PureSpeciesDiffusion<TPTDATATYPE,BOOL> PureDiffusionType;
    

    enum{DIM = SystemDimension<DataType::nspec,DataType::rednspec,BOOL>::Value};
   

    MixtureAveragedDiffusion(PointerType ptr = 0):PureDiff_(ptr),meanmw_(value_type()){
      BaseClassType::ptr_ = ptr;
    }

    MixtureAveragedDiffusion(const MixtureAveragedDiffusion& mad):PureDiff_(mad.PureDiff_),meanmw_(mad.meanmw_){
// #ifndef NDEBUG
//       if(DIM==1)
// 	ADONIS_ERROR(DimensionError,"Mixture-averaged diffusion only defined for mixtures with more than 1 species. \n    Use the self-diffusion coefficient instead!!")
// #endif
      BaseClassType::ptr_ = mad.ptr_;  //prevent slicing
    }

    MixtureAveragedDiffusion& operator=(const MixtureAveragedDiffusion& mad){
      if(this != &mad){
	PureDiff_ = mad.PureDiff_;
	meanmw_ = mad.meanmw_;
	BaseClassType::ptr_ = mad.ptr_;
      }
      return *this;
    }

    void initialize(PointerType ptr){
      BaseClassType::ptr_ = ptr;
      meanmw_ = value_type();
      PureDiff_.initialize(ptr);
    }

    PointerType coefficients() {return &mixDiff_[0];}

     value_type binary_diffusion_coefficients(int i, int j){
      return PureDiff_.coefficients()[SymmetricAccess::offset(i,j,DIM)];
    }

    //! get \f$ k\f$-th mix-average diff. coefficient
    const value_type& operator[](size_t k) const {return mixDiff_[k];}
    value_type& operator[](size_t k) {return mixDiff_[k];}

    //! implements formula 5-45, p. 90 in 
    //! [CHEMKIN-PRO, Reaction Design: San Diego, 2008]
    //! NOTE: This is only well defined when the mixture isn't exaclty 1 species
    template<class T1, class T2, class X>
    void compute_mixture_averaged_diffusion_coefficients(const T1& p, const T2& temp, const X& Xfrac){
#ifndef NDEBUG
      if(!is_well_defined_value(temp) || (temp < T2()) || (!is_well_defined_value(p)))
	ADONIS_ERROR(ValueError,"Bad vals: temperature = "<< temp << "   p = "<< p << ".");
#endif
      //!compute binary diffusion coefficients
      PureDiff_.compute_binary_diffusion_matrix(p,temp);
      
      meanmw_ = value_type(); //reset!!
      ComputeTransportProperties<BOOL>::template mean_molecular_weight_X<DataType>(meanmw_,Xfrac,DataType::molar_masses());
   
      ComputeTransportProperties<BOOL>::template compute_mixture_averaged_diffusion_coefficients<DataType>(mixDiff_,DataType::molar_masses(),meanmw_,Xfrac,PureDiff_.coefficients());
    }

    //! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!A. DEDNER wrote me 
    //! function argument U in conservative form, i.e. \f$ U = [\rho Y_1, \lots, \rho Y_M]\f$
    template <class T1, class X>
    void compute_mixture_averaged_diffusion_coefficients(const T1& temp, const X& U, X& diff){
      T1 frac = 0;
      T1 rho = 0;
      for (unsigned int i=0;i<U.size();++i)
      {
        frac += U[i] / DataType::molar_masses()[i];
	rho += U[i];
      }
      T1 p = rho*PhysicalConstants<T1>::IdealGasConstant*temp;
      X Xfrac;
      for (unsigned int i=0; i < U.size();++i)
      {
        Xfrac[i] = U[i] / (DataType::molar_masses()[i] * frac);
      }
      compute_mixture_averaged_diffusion_coefficients(p,temp,Xfrac);
      for (unsigned int i=0;i < U.size();++i)
      {
        diff[i] = mixDiff_[i];
      }

     }
    

  private:
    PureDiffusionType PureDiff_;
    value_type mixDiff_[DIM];
    value_type meanmw_;

   };
  
} //end namespace

#endif
