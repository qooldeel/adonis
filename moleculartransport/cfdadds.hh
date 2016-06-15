#ifndef NICE_CFD_RELATED_FEATURES_HH
#define NICE_CFD_RELATED_FEATURES_HH

#include <cmath>

#include "../common/numericaltraits.hh"
#include "../common/typeadapter.hh"
#include "../massactionkinetics/data/thermochemicaldata.hh"
#include "../massactionkinetics/thermochemistry.hh"
#include "../massactionkinetics/physicalconstants.hh"

namespace Adonis{

  /**
   * \brief computation of the  Reynolds number [KEE, pdf. 148, Eq. (3.242)]
   *
   * \param rho0 characteristic density (e.g. inlet density)
   * \param v0 characteristic velocity (e.g. inlet velocity)
   * \param x0 characteristic lenght (e.g. pipe diameter)
   * \param mu mixture viscosity (e.g. at inlet)
   */
  template<class T1, class T2>
  typename NumericalTypePromotion<T1,T2>::ReturnType Reynolds_number(const T1& rho0, const T1& v0, const T2& x0, const T1& mu){
    return (rho0*v0*x0)/mu;
  }

  /**
   * \brief Mach number [KEE, pdf. 148, Eq. (3.242)]
   */
  template<int MECH,class T1, class T2, class STATE>
  typename NumericalTypePromotion<T1,T2>::ReturnType Mach_number(const T1& temp0, const T2& v0, const STATE& Y, const T2& v20 = T2()){
    typedef typename TypeAdapter<T1>::BaseType BaseType;
    typedef ThermoData4Mechanism<BaseType,MECH> DataType;
    //calculate speed of sound first
    typedef NASA7CoefficientPolynomial<const BaseType*,const BaseType*> NasaPolyType;
    //only stores 3 ints
     NasaPolyType nasa(DataType::nspec,DataType::thermo(), DataType::temperature_bounds());

     T1 cp(0.), cV(0.), Cp0, wbar(0.);
     
     for(int k = 0; k < DataType::nspec; ++k){
       Cp0 = nasa.C_p(k,temp0);
       cp += Cp0/DataType::molar_masses()[k]*Y[k];
       cV += (Cp0 - PhysicalConstants<T1>::Rgas)/DataType::molar_masses()[k]*Y[k];
       wbar += Y[k]/DataType::molar_masses()[k];
     }
     wbar = 1./wbar;

     return ( sqrt(v0*v0 + v20*v20)/sqrt(cp/cV*(PhysicalConstants<BaseType>::Rgas/wbar)*temp0) );

  }





} //end namespace 

#endif
