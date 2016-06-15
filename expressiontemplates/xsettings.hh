#ifndef DEFINE_COMPILE_TIME_SETTINGS_FOR_XPRESSION_TEMPLATES_HH
#define DEFINE_COMPILE_TIME_SETTINGS_FOR_XPRESSION_TEMPLATES_HH

#include "../common/fancymessages.hh"
#include "../want_2_use.h"

namespace Adonis{
  namespace ExprTmpl{

    class XPCTSettings{
    public:
      
      //!====================================================================
      //!=========== Choose method of floating point summation ==============
      //!====================================================================
      //! 
      //! use either 'k','K' or 'b','B' for Kahan or Babuska-Kahan summation
      
      static const char floatingPointAdditionMethod = 'b';
      
      //=====================================================================
    };







    /**
   * \brief Inform me which addition is applied
   */
    inline void which_addition_is_used(){
      
      if(USE_TUNED_FLOATING_POINT_ARITHMETIC){
	if(XPCTSettings::floatingPointAdditionMethod == 'k' || XPCTSettings::floatingPointAdditionMethod == 'K')
	  FancyMessages().nice_output("\n Kahan algorithm for addition is applied... \n",36,40);
	if(XPCTSettings::floatingPointAdditionMethod == 'b' || XPCTSettings::floatingPointAdditionMethod == 'B')
	  FancyMessages().nice_output("\n Babuska-Kahan algorithm for addition is applied... \n",36,40);
      }
      else 
	FancyMessages().nice_output("\n Conventional addition is applied... \n",34,40);
    }


  } //end namespace Adonis
}   //end namespace ExprTmpl

#endif
