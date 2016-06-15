#ifndef EQUATION_OF_STATE_HH
#define EQUATION_OF_STATE_HH

#include "../common/numericaltraits.hh"
#include "physicalconstants.hh"

namespace Adonis{

  /**
   * provide easy eos generators for, e.g. density, pressure
   */
  template<char C> class EquationOfState;

  //! specializations
  //! ideal gas law 
  template<>
  class EquationOfState<'i'>{   //! ideal gas law 
  public:
    template<class P, class WBAR, class TEMP>
    static inline typename NumericalTypePromotion<P, typename NumericalTypePromotion<WBAR,TEMP>::ReturnType>::ReturnType density(const P& p, const WBAR& Wbar, const TEMP& temp){
      return (p*Wbar/(PhysicalConstants<P>::IdealGasConstant*temp));
    }
    
    template<class RHO, class WBAR, class TEMP>
    static inline typename NumericalTypePromotion<RHO, typename NumericalTypePromotion<WBAR,TEMP>::ReturnType>::ReturnType pressure(const RHO& rho, const WBAR& Wbar, const TEMP& temp){
      return ((rho*PhysicalConstants<RHO>::IdealGasConstant*temp)/Wbar);
    }
  };

  template<>
  class EquationOfState<'I'>{   //! ideal gas law -- same as specialization 'i'
  public:
    template<class P, class WBAR, class TEMP>
    static inline typename NumericalTypePromotion<P, typename NumericalTypePromotion<WBAR,TEMP>::ReturnType>::ReturnType density(const P& p, const WBAR& Wbar, const TEMP& temp){
      return EquationOfState<'i'>::density(p,Wbar,temp);
    }

     template<class RHO, class WBAR, class TEMP>
    static inline typename NumericalTypePromotion<RHO, typename NumericalTypePromotion<WBAR,TEMP>::ReturnType>::ReturnType pressure(const RHO& rho, const WBAR& Wbar, const TEMP& temp){
       return EquationOfState<'i'>::pressure(rho,Wbar,temp);
     }

  };

} //end namespace 

#endif
