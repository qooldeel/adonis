#ifndef AUXILIARY_STUFF_NEEDED_FOR_MASS_ACTION_KINETICS_HH
#define AUXILIARY_STUFF_NEEDED_FOR_MASS_ACTION_KINETICS_HH

#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"
#include "../common/error.hh"

namespace Adonis{

  template<bool> class TemperatureBoundsCheck;

  template<>
  class  TemperatureBoundsCheck<true>{
  public:
    
    template<class IT,class T>
    static inline void beware_of_bad_temperatures(std::size_t k, const T& temperature, IT boundsit){
#ifndef NDEBUG
      if((temperature < boundsit[3*k]) || (temperature > boundsit[3*k+1]) || (!is_well_defined_value(temperature)))
	ADONIS_ERROR(BoundsError, "Temperature out of bounds: T = "<<temperature << ".");
#endif
    }
  };

  //don't perform check on temperatures
  template<>
  class  TemperatureBoundsCheck<false>{
  public:
    
    template<class IT,class T>
    static inline void beware_of_bad_temperatures(std::size_t k, const T& temperature, IT boundsit){}
  };


  template<class T>
  T cal2joule(const T cal){
    return cal*4.184;
  }

} //end namespace

#endif
