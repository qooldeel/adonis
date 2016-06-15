#ifndef SPECIFY_IN_WHICH_QUANTITY_TO_COMPUTE_HH
#define SPECIFY_IN_WHICH_QUANTITY_TO_COMPUTE_HH

namespace Adonis{

  class Settings4Quantity{
  public:
    //! 'p' partial density, 's' specific moles
    static const char ofcomputation = 's';  
  
    //! for special mechanisms
    static const char h2c6 = 'c';  //already in concentrations
  
    //
    static const int thermoCase = 1; //1 = isochor, 2 = isobar
  };

} //end namespace Adonis 

#endif

