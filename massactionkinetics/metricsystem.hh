#ifndef METRIC_SYSTEM_WHICH_IS_A_SYSTEM_OF_UNITS_OF_MEASUREMENTS_HH
#define METRIC_SYSTEM_WHICH_IS_A_SYSTEM_OF_UNITS_OF_MEASUREMENTS_HH

#include <cmath>
#include "../misc/misctmps.hh"
#include "physicalconstants.hh"

namespace Adonis{

  /**
   * The most widely used metric systems are the SI (syst\`{e}me international d'unit\'es) and the cgs (centimetre-gram-second) system.
   *
   * Here, I list several units and some conversions to other units (only those
   * I am currently working with)
   *  1 C = 1 A·s
   *  1 statC = 1 Fr ~ 3.34e-10 C
   *  1 J = kg·m²/s²
   *  1 dyne = 1 g·cm/s² = 1 statC²/cm²
   *  1 erg = 1 dyne·cm = 1 g·cm²/s² = 1.e-7 J 
   *  1 Debye = 1.e-18 Fr·cm = 1.e-18 statC·cm ~ 3.33564e-30 C·m
   */
  template<class T, char C> class MetricSystem;

  //! convert units into real SI units (m, kg, s, A, K, mol, cd) or 
  //! SI-derived units such as J, C, etc.
  template<class T> 
  class MetricSystem<T,'I'>{  //"I"nternation System of Units (SI)
  public:
    typedef T value_type;
    
    //! debye to C·m
    static inline T debye(const T& val){
      return 1./PhysicalConstants<T>::speed_of_light_in_vacuum*1.e-21*val;  //~3.33564e-30*val;
    }

    //! Angstroem^N to m^N -- only take the power w.r.t. the prefactor
    //!USAGE: \code MetricSystem<T,'I'>::template Angstroem(val) \endcode  
    template<int N>
    static inline T Angstroem(const T& val){
      return ntimes<N>(1.e-10)*val;
    }

   
    //! erg into J 
    static inline T erg(const T& val){
      return 1.e-7*val;
    }
  };


  //! cgs system: convert units into cgs units or cgs-derived units
  template<class T>
  class MetricSystem<T,'C'>{ //"C"GS system
  public:
    typedef T value_type;
    
    //!debye into statC·cm
    static inline T debye(const T& val){
      return 1.e-18*val;
    }

    //! Angstroem^N to cm^N
    //!USAGE: \code MetricSystem<T,'C'>::template Angstroem(val) \endcode 
    template<int N>
    static inline T Angstroem(const T& val){
      return ntimes<N>(1.e-8)*val;
    }
    
    //! This cgs unit has no prefactor like the cgs unit debye
    static inline  const T& erg(const T& val){ return val;} 
      
    
  };

}//end namespace

#endif
