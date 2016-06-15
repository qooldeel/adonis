#ifndef STATE_YOUR_TRANSPORT_SETTINGS_HERE_HH
#define STATE_YOUR_TRANSPORT_SETTINGS_HERE_HH

namespace Adonis{

  class TransportSettings{
  public:
    //!true = calculate in g/mol, false = calculate in kg/mol
    static const bool molarMassInGrams = true; 
    static const bool pressureInBar = true;
    
    //!if 'true' it might be that \f$X_k \approx 0 \longrightarrow det(L) = 0\f$
    //!since one allows pure species situations to occur
    static const bool dontAllowPureSpeciesSituationsToOccur = false; //false;
  };

}//end namespace

#endif
