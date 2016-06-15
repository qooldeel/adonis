#ifndef THERMAL_DIFFUSION_LAYOUT_HH
#define THERMAL_DIFFUSION_LAYOUT_HH

#include "generaltransport.hh"
#include "../templatemetaprograms/unrollloop.hh"

#include "../common/smartassign.hh"
#include "../common/globalfunctions.hh"
#include "../massactionkinetics/physicalconstants.hh"

#include "fancytransporttmps.hh"
//#include "../templatemetaprograms/functionalities4containersandscalars.hh"

namespace Adonis{


  /** \brief needed when the reduced form (1st approx) is considered */
  template<class TPTDATATYPE, bool BOOL> class TransportVsReducedIndex;
  
   template<class TPTDATATYPE> 
   class TransportVsReducedIndex<TPTDATATYPE,true>{ //reduction
   public:
     typedef typename TPTDATATYPE::index_type SizeType;
     
     static inline SizeType transport_idx(SizeType k){return TPTDATATYPE::rpv_index()[k];}

     static inline SizeType usual_idx(SizeType k){return k;}
   };


  template<class TPTDATATYPE> 
  class TransportVsReducedIndex<TPTDATATYPE,false>{ //No reduction
  public:
    typedef typename TPTDATATYPE::index_type SizeType;
    
    static inline SizeType transport_idx(SizeType k){return k;} //just index
    
    static inline SizeType usual_idx(SizeType k){return k;}
  };


  /**
   * \brief Thermal diffusion ratio [1, p. 93, eq. 5-53] needed for the 
   *         computation of mixture-averaged as well as multicomponent 
   *         thermal diffusion 
   *    
   *  Ref.:
   *  [1] [CHEMKIN-PRO, Reaction Design: San Diego, 2008 -- Theory Manual]
   *  [2] [KEE et al., "Chemically Reacting Flow"]
   *
   * DISTINCTMOLMASS = false is default, i.e.  compute for all species \f$\Theta_k\f$
   *                   true: only \f$\Theta_k\f$ for species below some specific
   *                   molecular mass is computed
   */
  template<class TPTDATATYPE, bool BOOL, bool DISTINCTMOLMASS = false>
  class ThermalDiffusionRatio: public GeneralTransport<TPTDATATYPE,ThermalDiffusionRatio<TPTDATATYPE,BOOL,DISTINCTMOLMASS> >{ 
  public:
    typedef GeneralTransport<TPTDATATYPE,ThermalDiffusionRatio<TPTDATATYPE,BOOL,DISTINCTMOLMASS> > BaseClassType;

    typedef typename BaseClassType::DataType DataType;

    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::PointerType PointerType;

    typedef typename BaseClassType::SizeType SizeType;
    
    typedef TransportVsReducedIndex<TPTDATATYPE,BOOL> RightIndexType;
    
    enum{NSPEC = TPTDATATYPE::nspec, //full system dimension
	 DIM = SystemDimension<DataType::nspec,DataType::rednspec,BOOL>::Value
    };

    ThermalDiffusionRatio(PointerType ptr = 0){
      BaseClassType::ptr_ = ptr;
    }
    
    ThermalDiffusionRatio(const ThermalDiffusionRatio& Pure){
      BaseClassType::ptr_ = Pure.ptr_;
    }
    
    ThermalDiffusionRatio& operator=(const ThermalDiffusionRatio& Pure){
      if(this != &Pure){
	BaseClassType::ptr_ = Pure.ptr_;
      }
      return *this;
    }

    void initialize(PointerType ptr){
      BaseClassType::ptr_ = ptr;
    }
  
    
    //! get \f$ k\f$-th thermal diffusion ratio
    const value_type& operator[](SizeType k) const {
      adonis_assert((0<=k) && (k < DIM));
      return Theta_[k];
    }

    value_type& operator[](SizeType k){
      adonis_assert((0<=k) && (k < DIM));
      return Theta_[k];
    }


    //! [1, eq. 5-54, p. 93] 
    template<class T, class X>
    T theta(SizeType k, SizeType j, const T& temperature, const X& Xfrac) {
      //!transport array is always of full length, so are molar masses
      T Tstar = BaseClassType::reduced_temperature(RightIndexType::transport_idx(k),RightIndexType::transport_idx(j),temperature),
    	Akj = Astar(Tstar);
     
      return ( 15./2.*((2*Akj+5)*(6*Cstar(Tstar)-5))/(Akj*(16*Akj - 12*Bstar(Tstar) + 55))*(BaseClassType::molar_mass(RightIndexType::transport_idx(k))-BaseClassType::molar_mass(RightIndexType::transport_idx(j)))/(BaseClassType::molar_mass(RightIndexType::transport_idx(k)) + BaseClassType::molar_mass(RightIndexType::transport_idx(j)))*perturb_when_zero(Xfrac[k],1.e-16)*perturb_when_zero(Xfrac[j],1.e-16)  //! Xfracs have either full length or reduced length so no special indexing needed for them
);
      
    }

    //another way to access elements of the field Theta_
    PointerType ratios() {return &Theta_[0]; }


    //! [1, eq. 5-53, p. 93] 
    //! Note: \f$\Theta_k\f$ can be <BB>negative</BB>. Its sign describes the 
    //!        direction in which a molecule is driven along the thermal 
    //!        gradient 
    template<class T, class X>
    void compute_thermal_diffusion_ratios(const T& temperature, const X& Xfrac, const value_type& molecmass = 5.e-03){ //!molecmass in SI units!!!!
      ComputeTransportProperties<BOOL>::template compute_thermal_diffusion_ratios<DISTINCTMOLMASS,DataType>(*this,temperature,Xfrac, molecmass);
      
      //by hand -- only full
      // T sum=0;
      // //int fac;
      // for(SizeType k = 0; k < DataType::nspec; ++k){
      // 	sum = 0.;
      // 	for(SizeType j = 0; j < DataType::nspec; ++j){
      // 	  sum += !kronecker_delta(k,j)*theta(k,j,temperature,Xfrac);
      // 	}
      // 	smart_assign(Theta_[k],sum);
      // }


    }

  private:
    value_type Theta_[DIM];
  
    //some functions needed here
    //! calculate \f$A^*_{jk}(T^*_{jk}) = \sum_{n=0}^6 a_n(\ln(T^*_{jk}))^n,  \ldots\$ using Horner's rule, cf. [2, p. 521/22, Tab. 12.2, eqs. 12.141--12.143]
    template<class T>
    T Astar(const T& Tstar) {
      T LnTemp = log(Tstar);
      return ( ((((((1.720853282e-04*LnTemp + -1.313998345e-03)*LnTemp + 7.569367323e-04)*LnTemp + 1.188708609e-02)*LnTemp + -1.671975393e-02)*LnTemp + -7.065517161e-03)*LnTemp) + 1.106910525);
    }

    template<class T>
    T Bstar(const T& Tstar) {
      T LnTemp = log(Tstar);
      return ( ((((((2.492954809e-04*LnTemp + -1.445009039e-03)*LnTemp + -3.030372973e-03)*LnTemp + 2.512965407e-02)*LnTemp + -2.147636665e-03)*LnTemp + -1.140928763e-01)*LnTemp) + 1.199673577);
      
    }

    template<class T>
    T Cstar(const T& Tstar) {
      T LnTemp = log(Tstar);
      return ( ((((((-2.115417788e-04*LnTemp + 1.844922811e-03)*LnTemp + -2.260153363e-03)*LnTemp + -1.625859588e-02)*LnTemp + 3.250097527e-02)*LnTemp + 4.748325276e-02)*LnTemp) + 8.386993788e-01);
     
    }
 
  };

} //end namespace 


#endif
