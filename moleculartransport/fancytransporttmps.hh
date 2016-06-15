#ifndef SOME_FANCY_TRANSPORT_TEMP_META_PROGS_HH
#define SOME_FANCY_TRANSPORT_TEMP_META_PROGS_HH

#include "../templatemetaprograms/unrollloop.hh"
#include "../templatemetaprograms/matrixunroller.hh"
#include "../templatemetaprograms/threeloopsunrolled.hh"

#include "../massactionkinetics/indexinjector.hh"

namespace Adonis{

  /**
   * \tparam N, R full (N) and reduced (R) dimension, respecitvely. CAUTION: it    * is mandatory to pertain this order. Otherwise you'll probably get some 
   * weird results.
   * \tparam B specifies if reduction is applied (true) or full system is 
   *           considered (false = no red)
   */
  template<int N, int R, bool B> class SystemDimension;
  
  template<int N, int R>
  class SystemDimension<N,R,true>{  //!reduced
  public:
    enum{Value = R};
  };

  template<int N, int R>
  class SystemDimension<N,R,false>{ //! full
  public:
    enum{Value = N};
  };


  /**
   * \brief Proper index selection. 
   * 
   * When the reduced \f$\dot{\omega}\f$ is to
   * be constructed, we are only interested in those species and reactions that
   * are RPVs and contain the RPVs, respectively. 
   */
  template<class MECHANISM, bool B> class IndexHandler;
  
  template<class MECHANISM>
  class IndexHandler<MECHANISM,true>{  //reduced ODE model
  public:
    typedef std::size_t SizeType;
    typedef ReducedIndexInjector<const SizeType*,const SizeType*> IndexerType;

    static inline SizeType spec_access(SizeType k){
      return MECHANISM::rpv_index()[k];
    }

    static inline SizeType reac_access(SizeType i){
      return MECHANISM::rpv_reaction_index()[i];
    }

    template<class INDEXER>
    static inline void initialize(INDEXER& ixer){
      ixer.initialize(MECHANISM::rednspec,MECHANISM::rednreac,MECHANISM::rpv_index(),MECHANISM::rpv_reaction_index());
    }
  };

  template<class MECHANISM>
  class IndexHandler<MECHANISM,false>{  //! full ODE model, i.e.
                                        //! full \f$\dot{\omega}\f$
  public:
    typedef std::size_t SizeType;
    typedef CommonIndexInjector IndexerType;

    static inline SizeType spec_access(SizeType k){ //just the usual index
      return k;
    }

    static inline SizeType reac_access(SizeType i){ //just the usual index
      return i;
    }

    template<class INDEXER>
    static inline void initialize(INDEXER& ixer){
      ixer.initialize(MECHANISM::nspec,MECHANISM::nreac);
    }
  };


  /**
   * \brief Choose if you wanna use reduced (true) or normal (false) 
   * description 
   * \tparam B choose (true) or don't choose reduction (false)
   */
  template<bool B> class ComputeTransportProperties;

  template<>
  class ComputeTransportProperties<true>{  //! reduced version
  public:
    template<bool DISTINCTMOLMASS, class TPTDATATYPE, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline void compute_thermal_diffusion_ratios(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac, const typename THERMALDIFFUSIONRATIO::value_type& molec){
      typedef typename THERMALDIFFUSIONRATIO::DataType DataType;
      UnrollLoop<0,DataType::rednspec>::template compute_thermal_diffusion_ratios_reduced<DISTINCTMOLMASS,TPTDATATYPE>(tdr,temp,Xfrac,molec);
    }

    template<class DIFFU, class T1, class T2>
    static inline void compute_binary_diffusion_matrix(DIFFU& diff, const T1& p, const T2& temperature){
#ifndef NDEBUG
      if(!is_well_defined_value(temperature) || (temperature < T2()) || (!is_well_defined_value(p)))
	ADONIS_ERROR(ValueError,"Bad vals: temperature = "<< temperature << "   p = "<< p << ".");
#endif
      typedef typename DIFFU::DataType DataType;
      MatrixUnroller<0,0,DataType::rednspec,DataType::rednspec>::compute_binary_diffusion_matrix_reduced(diff,p,temperature);
    }

    template<class TPTDATATYPE, class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_Y(T& mmw, const IT1& y, const IT2& w){
      return UnrollLoop<0,TPTDATATYPE::rednspec>::template mean_molecular_weight_Y_reduced<TPTDATATYPE>(mmw,y,w);
    }

    template<class TPTDATATYPE, class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_X(T& mmw, const IT1& y, const IT2& w){
      return UnrollLoop<0,TPTDATATYPE::rednspec>::template mean_molecular_weight_X_reduced<TPTDATATYPE>(mmw,y,w);
    }

    template<class TPTDATATYPE, class DMIXIT, class MW, class T, class X, class PITER>
     static inline void compute_mixture_averaged_diffusion_coefficients(DMIXIT dmixit, const MW& mw, const T& meanmw, const X& Xfrac, const PITER& bindiff){
       UnrollLoop<0,TPTDATATYPE::rednspec>::template compute_mixture_averaged_diffusion_coefficients_reduced<TPTDATATYPE>(dmixit,mw,meanmw,Xfrac,bindiff);
  }

    //! NEW: compute \f$ D_k^{\mathrm{mix}}\f$ for a fixed K
    template<int K, class TPTDATATYPE, class MW, class T, class X, class PITER>
    static inline typename X::value_type Dmix(const MW& mw, const T& meanmw, const X& Xfrac, const PITER& bindiff){
      return UnrollLoop<0,TPTDATATYPE::rednspec>::template mix_avg_diff_nominator_reduced<TPTDATATYPE,K>(mw,Xfrac)/(meanmw*UnrollLoop<0,TPTDATATYPE::rednspec>::template mix_avg_diff_denominator_reduced<TPTDATATYPE,K>(Xfrac,bindiff));
    }

    template<class VIS, class T>
    static inline void compute_viscosities(VIS& vis, const T& temp){
#ifndef NDEBUG
      if(!is_well_defined_value(temp) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad vals: temperature = "<< temp <<  ".");
#endif

      typedef typename VIS::DataType DataType; 
      UnrollLoop<0,DataType::rednspec>::compute_viscosities_reduced(vis,temp);
    }

    template<class COND, class T, class RHO, class P>
    static inline void compute_conductivities(COND& cond, const RHO& rho, const P& p, const T& temp){
#ifndef NDEBUG
      if(!is_well_defined_value(temp) || (temp < T()) || (!is_well_defined_value(p)) || (!is_well_defined_value(rho)) || (rho < RHO()))
	ADONIS_ERROR(ValueError,"Bad vals: temperature = "<< temp << "   p = "<< p << "   rho = " << rho <<  ".");
#endif
      typedef typename COND::DataType DataType; 
      UnrollLoop<0,DataType::rednspec>::compute_conductivities_reduced(cond,rho,p,temp);
    }


    template<class TPTDATATYPE, class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void compute_multicomponent_diffusion_matrix(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){
      MatrixUnroller<0,0,TPTDATATYPE::rednspec,TPTDATATYPE::rednspec>::template compute_multicomponent_diffusion_matrix_reduced<TPTDATATYPE>(mcdiff,p,temp,xfrac,mw,meanmw,inv);
    }

    template<class TPTDATATYPE, class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void create_L0000(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){
      ThreeNestedLoopsUnroller<0,0,0,TPTDATATYPE::rednspec,TPTDATATYPE::rednspec,TPTDATATYPE::rednspec>::template create_L0000_reduced<TPTDATATYPE>(Diff,p,temp,x,mw,d);
    }
  };


  //!GENERAL (FULL) SETTING
  template<>
  class ComputeTransportProperties<false>{  //! general setting
  public:
    template<bool DISTINCTMOLMASS, class TPTDATATYPE, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline void compute_thermal_diffusion_ratios(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac, const typename THERMALDIFFUSIONRATIO::value_type& molec){
      UnrollLoop<0,TPTDATATYPE::nspec>::template compute_thermal_diffusion_ratios<DISTINCTMOLMASS>(tdr,temp,Xfrac,molec);
    }

    template<class DIFFU, class T1, class T2>
    static inline void compute_binary_diffusion_matrix(DIFFU& diff, const T1& p, const T2& temperature){
      typedef typename DIFFU::DataType DataType;
      MatrixUnroller<0,0,DataType::nspec,DataType::nspec>::compute_binary_diffusion_matrix(diff,p,temperature);
    }

    template<class TPTDATATYPE, class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_Y(T& mmw, const IT1& y, const IT2& w){
      return UnrollLoop<0,TPTDATATYPE::nspec>::mean_molecular_weight_Y(mmw,y,w);
    }

    template<class TPTDATATYPE, class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_X(T& mmw, const IT1& y, const IT2& w){
      return UnrollLoop<0,TPTDATATYPE::nspec>::mean_molecular_weight_X(mmw,y,w);
    }

    template<class TPTDATATYPE, class DMIXIT, class MW, class T, class X, class PITER>
    static inline void compute_mixture_averaged_diffusion_coefficients(DMIXIT dmixit, const MW& mw, const T& meanmw, const X& Xfrac, const PITER& bindiff){
      UnrollLoop<0,TPTDATATYPE::nspec>::compute_mixture_averaged_diffusion_coefficients(dmixit,mw,meanmw,Xfrac,bindiff);
    }
    
    //! NEW: compute \f$ D_k^{\mathrm{mix}}\f$ for a fixed K
    template<int K, class TPTDATATYPE, class MW, class T, class X, class PITER>
    static inline typename X::value_type Dmix(const MW& mw, const T& meanmw, const X& Xfrac, const PITER& bindiff){
      return UnrollLoop<0,TPTDATATYPE::nspec>::template mix_avg_diff_nominator<TPTDATATYPE,K>(mw,Xfrac)/(meanmw*UnrollLoop<0,TPTDATATYPE::nspec>::template mix_avg_diff_denominator<TPTDATATYPE,K>(Xfrac,bindiff));
    }


    template<class VIS, class T>
    static inline void compute_viscosities(VIS& vis, const T& temp){
      typedef typename VIS::DataType DataType; 
      UnrollLoop<0,DataType::nspec>::compute_viscosities(vis,temp);
    }

    template<class COND, class T, class RHO, class P>
    static inline void compute_conductivities(COND& cond, const RHO& rho, const P& p, const T& temp){
      typedef typename COND::DataType DataType; 
      UnrollLoop<0,DataType::nspec>::compute_conductivities(cond,rho,p,temp);
    }

    template<class TPTDATATYPE, class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void compute_multicomponent_diffusion_matrix(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){
      MatrixUnroller<0,0,TPTDATATYPE::nspec,TPTDATATYPE::nspec>::compute_multicomponent_diffusion_matrix(mcdiff,p,temp,xfrac,mw,meanmw,inv);
    }

    template<class TPTDATATYPE, class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void create_L0000(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){
      ThreeNestedLoopsUnroller<0,0,0,TPTDATATYPE::nspec,TPTDATATYPE::nspec,TPTDATATYPE::nspec>::create_L0000(Diff,p,temp,x,mw,d);
    }
  };


}  //namespace

#endif
