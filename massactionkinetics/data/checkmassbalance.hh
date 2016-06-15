#ifndef CHECK_ELEMENTARY_MASS_BALANCE_HH
#define CHECK_ELEMENTARY_MASS_BALANCE_HH

#include "../../common/globalfunctions.hh"
#include "../../common/isclass.hh"
#include "../../common/typeadapter.hh"
#include "../../common/typeselector.hh"

namespace Adonis{

  template<class V, int MECHENCODING> class CheckMassBalance;

  //! H2C6 Mechanism
  template<class V>
  class CheckMassBalance<V,6>{
  public:
    typedef typename ValueTypeSelector<V,IsADuneFieldContainer<V>::Value>::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;

    static inline bool check_mass(const V& z, const value_type& b1 = 2., const value_type& b2 = 1., const BaseType& tol = 1.e-05){
      value_type m1 = 2.*z[0] + 2.*z[4] + z[1] + z[5],
	m2 = 2.*z[2] + z[4] + z[3] + z[5];

      if(Abs(m1-b1) > tol || Abs(m2-b2) > tol){
	//std::cout << "z = " << z  << std::endl;
	std::cout << "MASS CONSERVED?   mass1: "<< m1 << "   mass2: "<< m2 << std::endl; 
	return false;
      }
      else return true;
    }
  };

} //end namespace Adonis

#endif
