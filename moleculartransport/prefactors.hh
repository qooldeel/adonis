#ifndef PRACTICAL_PREFACTORS_FOR_TRANSPORT_COEFFICIENTS_HH
#define PRACTICAL_PREFACTORS_FOR_TRANSPORT_COEFFICIENTS_HH

#include <cmath>

#include "../misc/misctmps.hh"
#include "../massactionkinetics/physicalconstants.hh"

namespace Adonis{

  template<class T>
  class PracticalPrefactor{
  public:
    
    static const T prefac_viscoeff_g_per_mol;
    
    static const T prefac_viscoeff_kg_per_mol;

    static const T prefac_diffcoeff_g_per_mol; 

    static const T prefac_diffcoeff_kg_per_mol;

    static const T kg_2_g;

    static const T g_2_kg;
  };

  template<class T> const T PracticalPrefactor<T>::prefac_viscoeff_g_per_mol = 0.3125/M_PI*sqrt(M_PI*PhysicalConstants<T>::kBBase/PhysicalConstants<T>::NBase)/sqrt(1000)*1.e-3; 
  //same as: 5./16.*sqrt(1./1000.)*sqrt(M_PI*PhysicalConstants<double>::Rgas/ntimes<2>(PhysicalConstants<double>::AvogadrosConstant))*1.e+20/M_PI

  template<class T> const T PracticalPrefactor<T>::prefac_viscoeff_kg_per_mol = 0.3125/M_PI*sqrt(M_PI*PhysicalConstants<T>::kBBase/PhysicalConstants<T>::NBase)*1.e-3;
  

  //! note that both below-mentionedprefactors are for the pure diffusion 
  //!coefficient in m²/s, although they scale quite differently!
  //! therefore, it might be more benefical to enter the formulas with molecular
  //! masses in g/mol. Prefac ~ 2.6635e-7
  //! INPUT: molar mass in g/mol, \f$ sigma \matrm{\ in \r{A}, p \matrm{\in bar}}\f$
  template<class T> const T PracticalPrefactor<T>::prefac_diffcoeff_g_per_mol = 0.375/M_PI*sqrt(M_PI*PhysicalConstants<T>::NBase*ntimes<3>(PhysicalConstants<T>::kBBase))*1./sqrt(1000)*1.e-5;
  //same as: 3./8.*sqrt(1000)*1.e+15*sqrt(M_PI*ntimes<3>(PhysicalConstants<double>::Rgas)/ntimes<2>(PhysicalConstants<double>::AvogadrosConstant))/M_PI;

  //! INPUT: as above except for molar mass: Now in kg/mol
  template<class T> const T PracticalPrefactor<T>::prefac_diffcoeff_kg_per_mol = PracticalPrefactor<T>::prefac_diffcoeff_g_per_mol*sqrt(1000.)*1.e-3;
  
  template <class T> const T PracticalPrefactor<T>::kg_2_g = 1e+3;
  template <class T> const T PracticalPrefactor<T>::g_2_kg = 1e-3;



  //! choose the right prefactor at compile time
  template<class T, bool B> class ChooseRightPrefactor;

  template<class T>
  class ChooseRightPrefactor<T,true>{  //in gram(/mol)
  public:
    static inline T visprefac(){return PracticalPrefactor<T>::prefac_viscoeff_g_per_mol;}

    static inline T diffprefac(){return PracticalPrefactor<T>::prefac_diffcoeff_g_per_mol;}

    //! convert e.g. kg/mol (SI) to g/mol needed for diffusion coeff
    static inline T in_gram(const T& val){return PracticalPrefactor<T>::kg_2_g*val;}

    //! since we assume pressure given in Pa = kg/(m·s²), one often intends
    //! to state pressure in bar (1 bar = 100000 Pa)
    static inline T in_bar(const T& p){return 1e-05*p;} // i.e. Pa*1e-05 = bar 
  };
  
  template<class T>         //SI 
  class ChooseRightPrefactor<T,false>{  //in kg(/mol)
  public:
    static inline T visprefac(){return PracticalPrefactor<T>::prefac_viscoeff_kg_per_mol;}

    static inline T diffprefac(){return PracticalPrefactor<T>::prefac_diffcoeff_kg_per_mol;}
  
    //! no prefactor since you intend to calculate in SI 
    static inline T in_gram(const T& val){return static_cast<T>(1);}
  
    //! no special prefactor since you intend to calculate in SI
    static inline T in_bar(const T& p){return static_cast<T>(1);} 
  };
  

} //end namespace

#endif
