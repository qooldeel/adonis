#ifndef PHYSICAL_AND_MATHEMATICAL_CONSTANTS_HH
#define PHYSICAL_AND_MATHEMATICAL_CONSTANTS_HH

#include "../common/universalconstants.hh"

namespace Adonis{
  
  /**
   * \brief Some physical constants. 
   *
   * To avoid scaling problems, use SI units throughout
   *
   * In some cases see
   *
   * References:
   *
   * [1] MOHR, TAYLOR and NEWELL, <I>"CODATA recommended values of the fundamental physical constants"</I>, Rev. Mod. Phys. 80, 633 -- 730 (2005)
   */
  template<class T>
  class PhysicalConstants{
  public:
   
    static const T Pi;

    static const T NBase;

    static const T kBBase;

    // universal gas constant
    static const T IdealGasConstant; //unit: J(K·mol) = (Pa·m³)/(K·mol)
    
    static const T Rgas;  //same as IdealGasConstant

    static const T AvogadrosConstant;
  
    static const T p0;   //! 1bar := 1.e+5 Pa (SI)--is equal everywhere on earth
    
    static const T p0_bar;

    static const T kB; //! Boltzmann constant \f$ k_B = R/N_A\f$ in J/K
    //CGS variant
    static const T kB_cgs; //! Boltzmann constant in erg/K, 1erg = \f$10^{-7}\f$J
    
    static const T speed_of_light_in_vacuum; //cf. [1], p. 637

    static const T gravitational_acceleration;

  };
  
  //!more accurate than M_PI, which is defined in math.h and holds 20 
  //!decimal places, M_PIl holds 34 decimal places for long doubles
  template<class T> const T PhysicalConstants<T>::Pi = UniversalConstants<T>::Pi;   //3.14159265358979323846264338327950288419716939937510;

  //! cf. [1], unit: J/(K·mol) = (Pa·m³)/(K·mol)
  template<class T> const T PhysicalConstants<T>::IdealGasConstant = 8.314472; 
   template<class T> const T PhysicalConstants<T>::Rgas = 8.314472; 

  //floating point without exponent for Avogadro constant
  template<class T> const T PhysicalConstants<T>::NBase = 6.02214179;
  
  //! cf. [2]
  //floating point without exponent for Boltzmann constant
  template<class T> const T PhysicalConstants<T>::kBBase = PhysicalConstants<T>::IdealGasConstant/NBase;

  template<class T> const T PhysicalConstants<T>::AvogadrosConstant = NBase*1.e+23;
  
  template<class T> const T PhysicalConstants<T>::p0 = 1.e+5; 
  template<class T> const T PhysicalConstants<T>::p0_bar = 1.; //1bar = 1e+5Pa

  //! unit: J/K (i.e. SI conform)
  template<class T> const T PhysicalConstants<T>::kB = PhysicalConstants<T>::IdealGasConstant/PhysicalConstants<T>::AvogadrosConstant;
  
  template<class T> const T PhysicalConstants<T>::kB_cgs = 1.e+07*PhysicalConstants<T>::kB;

  template<class T> const T PhysicalConstants<T>::speed_of_light_in_vacuum = 299792458; //in m/s

  template<class T> const T PhysicalConstants<T>::gravitational_acceleration = 9.80665; //in m/s²


  //============================================================================

  template<class T, char C> class PCInRightUnits;


  template<class T>
  class PCInRightUnits<T,'I'>{ //SI 
  public:
    static inline T kB(){
      return PhysicalConstants<T>::kB;
    }
  };
  
  template<class T>
  class PCInRightUnits<T,'C'>{ //CGS system
  public:
    static inline T kB(){
      return PhysicalConstants<T>::kB_cgs;
    }
  };
  

  template<bool C> class Transform2Gram;

  template<>
  class Transform2Gram<true>{
  public:
    template<class T>
    static inline T from_kg(const T& val){return 1.e+3*val;}
  };

  template<>
  class Transform2Gram<false>{
  public:
    template<class T>
    static inline const T& from_kg(const T& val){return val;} //identity
  };

}//end namespace

#endif
