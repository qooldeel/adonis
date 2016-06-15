#ifndef MY_RANDOM_NUMBER_GENERATORS_HH
#define MY_RANDOM_NUMBER_GENERATORS_HH

#include "lcgminstd.hh"  //from Numerical recipes in C
#include "lecuyer32.hh" //a small and good one

namespace Adonis{
  
  /**
   * \brief Choose you favorite random number generator at compile time.
   * \tparam FLPOINTTYPE floating point type, e.g. when values from [0,1) are to be computed
   * \tparam INTPRECTYPE bit type of random number generator (many of the common RNGs still use 32 bit ints)
   * \tparam RNG Specifies which kind of generator is to be chosen
   */
  template<class FLPOINTTYPE, int INTPRECTYPE, char RNG> class RandomNumberGenerator;

  //! specializations
  //! linear feedback shift register generator 
  template<class FLPOINTTYPE>
  class RandomNumberGenerator<FLPOINTTYPE,32,'L'>{ // L'Ecuyer, 32 bit
  public:
    typedef RNGlfsr113<FLPOINTTYPE> Type;
  };

  template<class FLPOINTTYPE>
  class RandomNumberGenerator<FLPOINTTYPE,32,'l'>{ // L'Ecuyer, 32 bit
  public:
    typedef RandomNumberGenerator<FLPOINTTYPE,32,'L'> Type;
  };

  //! linear congruent minimal standard generator
  template<class FLPOINTTYPE>
  class RandomNumberGenerator<FLPOINTTYPE,32,'M'>{ // minimal standard, 32 bit
  public:
    typedef LCGminstd<FLPOINTTYPE> Type;
  };

  template<class FLPOINTTYPE>
  class RandomNumberGenerator<FLPOINTTYPE,32,'m'>{ // L'Ecuyer, 32 bit
  public:
    typedef RandomNumberGenerator<FLPOINTTYPE,32,'M'> Type;
  };


} //end namespace

#endif
