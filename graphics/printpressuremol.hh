#ifndef RECORD_PRESSURE_FOR_MOL_WHEN_ITS_NOT_IN_DIFFERENTIAL_FORM_HH
#define RECORD_PRESSURE_FOR_MOL_WHEN_ITS_NOT_IN_DIFFERENTIAL_FORM_HH

#include "../common/universalconstants.hh"
#include "printer.hh"
#include "../massactionkinetics/physicalconstants.hh"

#include "../fdm/gensettings.hh"

namespace Adonis{

  //! forward declaration
  template<bool PRINT> class RecordPressure4MOL;


  template<>
  class RecordPressure4MOL<true>{  //record pressure as p = (rho·R·T)/Wbar
  public:
    
    //!no selectvars necessary since pressure is computed from the WHOLE species
    template<class V>
    static inline void two_D(int i, int j, int nx, PrintSolution& PS, const V& y, int npt, const typename V::value_type* pMolarMasses){
      typedef typename V::value_type value_type;
      if(pMolarMasses == 0){
	ADONIS_ERROR(IncompleteOrMissingDataError,"In order to use this functionality, you MUST provide non-empty molar masses!");
      }
      

      adonis_assert((npt%nx == 0) && ((int)y.size()%npt == 0));
      int numOfSpecs = (int)y.size()/npt - 4;  //4 coz minus rho, v1, v2, T 
      
      value_type Wbar(0);
      for(int k = 0; k < numOfSpecs; ++k){
	Wbar += (
#ifdef NONCONSERVATIVE_FORM		 
		 y[i + j*nx + (4+k)*npt]
#else //conservative form: devision by rho required 
		 y[i + j*nx + (4+k)*npt]/y[i + j*nx]
#endif
		 /pMolarMasses[k] );
      }
      if(is_zero(Wbar))
	Wbar = UniversalConstants<value_type>::smallValue;
      else
	Wbar = 1./Wbar;

      //!calculate ideal gas pressure and write to file that is handled by PS
      PS.get_ofstream() << "   " << std::setprecision(PS.get_precision()) << 
	// rho           R                                   T        
	( (y[i + j*nx]*PhysicalConstants<value_type>::Rgas*
#ifdef NONCONSERVATIVE_FORM
	   y[i + j*nx + 3*npt]
#else  //conservative form: division by rho required
	   y[i + j*nx + 3*npt]/y[i + j*nx]
#endif
	   )/Wbar ) << " ";
    }
  };
  

  template<>
  class RecordPressure4MOL<false>{  //DO NOTHING
  public:
    //!no selectvars necessary since pressure is computed from the WHOLE species
    template<class V>
    static inline void two_D(int i, int j, int nx, PrintSolution& PS, const V& y, int npt, const typename V::value_type* pMolarMasses){}
  };

} //end namespace 

#endif
