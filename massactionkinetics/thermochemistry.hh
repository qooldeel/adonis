#ifndef THERMOCHEMICAL_REFERENCE_DATA_OBJECTS_HH
#define THERMOCHEMICAL_REFERENCE_DATA_OBJECTS_HH

#include <iostream>
#include "auxiliary.hh"
#include "physicalconstants.hh"

namespace Adonis{

  /**
   * \brief 7-coefficient NASA polynomials (we use a more efficient 
   * implementation  using <B>Horner scheme</B>, cf. [2, p. 6 (formulas),
   * p. 17 (example)])
   * \tparam ITER1 stores iterator of thermochemical coefficients in 14- 
   * value blocks (the first 7 correspond to HIGH temperature coeffs the
   * following 7 for LOW temp. coeffs) -- one for each species \f$k.\f$
   * \tparam ITER2 stores iterator of temperature bounds in 3-values blocks 
   * (low, high and switching temperature, respectively) -- one for each 
   * species \f$k\f$ 
   * \tparam BC perform bounds check on temperature (recommended; set true by 
   *            default, i.e. yes.
   *
   * CAUTION: in order to keep storage minimal, I use iterators and don't
   *          supply any range checking utilities ;)
   *
   * References:
   *
   * [1] WARNATZ et al., "Combustion", 4th ed., p. 54
   *
   * [2] BURCAT and RUSCIC, "3rd Millenium Ideal Gas and Condensed Phase Thermodynamical Database for Combustion with Updates from Active Thermodynamical Tables", Tec. Report, 2005, p. 6 and p. 17
   
   */
  template<class ITER1, class ITER2, bool BC = true>
  class NASA7CoefficientPolynomial{
  public:
  
    typedef ITER2 TemperatureBoundsIterType;

    NASA7CoefficientPolynomial(size_t n = 0, const ITER1& i1 = ITER1(), const ITER2& i2 = ITER2()):numberOfSpecies_(n),coeffs_(i1),temps_(i2){}

    void initialize(size_t n, const ITER1& i1, const ITER2& i2){
      numberOfSpecies_ = n;
      coeffs_ = i1;
      temps_ = i2;
    }

    size_t number_of_species() const {return numberOfSpecies_;} 

    ITER2 temperature_bounds() const{ return temps_;}

    //!index of low \f$m\f$-th coefficients of species \f$k\f$ 
    //2nd set of 7 coeffs
    size_t low_index(size_t k, size_t m) const{
      adonis_assert((k < numberOfSpecies_) && (m < 7));
      return k*14 + 7 + m;
    }

    //1st set of 7 coeffs
    size_t high_index(size_t k, size_t m) const {
      adonis_assert((k < numberOfSpecies_) && (m < 7));
      return k*14 + m;
    }

    //! with 'eps' we can enlarge the temperature interval
    template<class T>
    bool check_temperature_bounds(size_t k, const T& temp, const T& eps = T()) const{
      //std::cout << "checked temperature: "<< temp << std::endl;
      return ( ((temps_[3*k]-eps <= temp) && (temp <= temps_[3*k+1]+eps)) ? true : false);
    }

    

    //! heat capacity \f$C_p(T)\f$ of species \f$ k\f$ at temperature \f$temp\f$
    template<class T>
    T heat_capacity(size_t k, const T& temp) const{
      size_t d = 3*k;
      //!check if t_low <= temp <= t_high
      //adonis_assert((temps_[d] <= temp) && (temp <= temps_[d+1]));

      TemperatureBoundsCheck<BC>::beware_of_bad_temperatures(k,temp,temps_);

      T C_p;

      if (temp < temps_[d+2]) //check if temp < switching temperature
	C_p = (((coeffs_[low_index(k,4)]*temp + coeffs_[low_index(k,3)])*temp + coeffs_[low_index(k,2)])*temp + coeffs_[low_index(k,1)])*temp + coeffs_[low_index(k,0)];
      else
	C_p = (((coeffs_[high_index(k,4)]*temp + coeffs_[high_index(k,3)])*temp + coeffs_[high_index(k,2)])*temp + coeffs_[high_index(k,1)])*temp + coeffs_[high_index(k,0)];
      
#ifndef NDEBUG
      if(!is_well_defined_value(C_p))
	ADONIS_ERROR(ValueError,"Bad C_p = "<< C_p <<".");
#endif
      return C_p*PhysicalConstants<T>::IdealGasConstant;
    }

    //! enthalpy \$H_T(T)\$ of species \f$ k\f$ at temperature \f$temp\f$
    template<class T>
    T enthalpy(size_t k, const T& temp) const{
       size_t d = 3*k;
      //!check if t_low <= temp <= t_high
      //adonis_assert((temps_[d] <= temp) && (temp <= temps_[d+1]));

       TemperatureBoundsCheck<BC>::beware_of_bad_temperatures(k,temp,temps_);

      T H_T;
      if (temp < temps_[d+2]) //check if temp < switching temperature
	H_T = ((((coeffs_[low_index(k,4)]/5.*temp + coeffs_[low_index(k,3)]/4.)*temp + coeffs_[low_index(k,2)]/3.)*temp + coeffs_[low_index(k,1)]/2.)*temp + coeffs_[low_index(k,0)])*temp + coeffs_[low_index(k,5)];
      else 
	H_T = ((((coeffs_[high_index(k,4)]/5.*temp + coeffs_[high_index(k,3)]/4.)*temp + coeffs_[high_index(k,2)]/3.)*temp + coeffs_[high_index(k,1)]/2.)*temp + coeffs_[high_index(k,0)])*temp + coeffs_[high_index(k,5)];
    
#ifndef NDEBUG
      if(!is_well_defined_value(H_T))
	ADONIS_ERROR(ValueError,"Bad H_T = "<< H_T <<".");
#endif
      return H_T*PhysicalConstants<T>::IdealGasConstant;
    }

    //! entropy \f$S_T(T)\f$ of species \f$ k\f$ at temperature \f$temp\f$
    template<class T>
    T entropy(size_t k, const T& temp) const{
      size_t d = 3*k;
      //!check if t_low <= temp <= t_high
      //adonis_assert((temps_[d] <= temp) && (temp <= temps_[d+1]));
      TemperatureBoundsCheck<BC>::beware_of_bad_temperatures(k,temp,temps_);
   
      T S_T;
      if (temp < temps_[d+2])
	S_T = (((coeffs_[low_index(k,4)]/4.*temp +  coeffs_[low_index(k,3)]/3.)*temp + coeffs_[low_index(k,2)]/2.)*temp + coeffs_[low_index(k,1)])*temp + coeffs_[low_index(k,0)]*log(temp) + coeffs_[low_index(k,6)];
      else
	S_T = (((coeffs_[high_index(k,4)]/4.*temp +  coeffs_[high_index(k,3)]/3.)*temp + coeffs_[high_index(k,2)]/2.)*temp + coeffs_[high_index(k,1)])*temp + coeffs_[high_index(k,0)]*log(temp) + coeffs_[high_index(k,6)];

#ifndef NDEBUG
      if(!is_well_defined_value(S_T))
	ADONIS_ERROR(ValueError,"Bad S_T = "<< S_T <<".");
#endif
      return S_T*PhysicalConstants<T>::IdealGasConstant;   
    }
    
    //! these abbreviations can also be used instead of the above ones whenever
    //! you feel more comfortable in doing so ;)
    template<class T>
    T C_p(size_t k, const T& temp) const{
      return (*this).heat_capacity(k,temp);
    }

    template<class T>
    T H_T(size_t k, const T& temp) const{
      return (*this).enthalpy(k,temp);
    }

    template<class T>
    T S_T(size_t k, const T& temp) const{
      return (*this).entropy(k,temp);
    }
      
    

  private:
    size_t numberOfSpecies_;
    ITER1 coeffs_;
    ITER2 temps_;
  };


  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //just testing
  template<class T, class ITER>
  inline T heat_capacity(const T& temp, const ITER& a){
    return (a[0] + a[1]*temp + a[2]*temp*temp + a[3]*temp*temp*temp + a[4]*temp*temp*temp*temp)*8.314471;
  }

  template<class T, class ITER>
  inline T enthalpy(const T& temp, const ITER& a){
    return (a[5] + a[0]*temp + a[1]/2.*temp*temp + a[2]/3.*temp*temp*temp + a[3]/4.*temp*temp*temp*temp + a[4]/5.*temp*temp*temp*temp*temp)*8.314471;
  }

  template<class T, class ITER>
  inline T entropy(const T& temp, const ITER& a){
    return (a[6] + a[0]*log(temp) + a[1]*temp + a[2]/2.*temp*temp + a[3]/3.*temp*temp*temp + a[4]/4.*temp*temp*temp*temp)*8.314471;
  }
 
} //end namespace

#endif
