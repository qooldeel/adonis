#ifndef CHOOSE_STUFF_4_2D_MOL_HH
#define CHOOSE_STUFF_4_2D_MOL_HH


//MACROS

#define FULL_MODEL 0

#define SWITCH_ON_REDUCTION 1



#include <string>

namespace Adonis{

  class ChooseFDSetting{
  public:
    static const int MECH = 3;

    static const bool ISREDUCED = false;  //DEPRECATED
   
    static const std::string FILE2READIN;
  };

  const std::string ChooseFDSetting::FILE2READIN = "2Dsettings.dat";
  

} //end namespace

#endif
