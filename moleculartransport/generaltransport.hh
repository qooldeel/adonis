#ifndef GENERAL_MOLECULAR_TRANSPORT_INTERFACE_HH
#define GENERAL_MOLECULAR_TRANSPORT_INTERFACE_HH

#include <iostream>
#include <cmath>

#include "../common/adonisassert.hh"
#include "../massactionkinetics/physicalconstants.hh"
#include "prefactors.hh"
#include "../misc/misctmps.hh"

#include "../common/globalfunctions.hh"
#include "../templatemetaprograms/matrixunroller.hh"
#include "../templatemetaprograms/threeloopsunrolled.hh"

namespace Adonis{

  /**
   * \brief General base class for transport.
   * \tparam TPTDATATYPE data for transport. Transport data may be found in [1] 
   * and [2]. The data are assumed to be stored as in [2], pp. 44 -- 50; to wit,
   * geometry, \$f \epsilon/k_B, \sigma, \mu, \alpha, Z_{\mathrm{rot}}\f$. 
   * These 6-tuples are stored consecutively in memory (i.e. six entries for 
   * species \f$k \f$, the next 6 entries for species \f$ l\f$ and so forth.
   * \tparam the derived transport class, e.g. for diffusion 
   *
   * References:
   *
   *  [1] KEE,RUPLEY, MILLER, COLTRIN, GRCAR, MEEKS, MOFFAT, LUTZ, DIXON-LEWIS,
   SMOOKE, WARNATZ, EVANS, LARSON, MITCHELL, PETZOLD, REYNOLDS, CARCOTSIOS, 
   STEWART, GLARBORG, WANG and ADIGUN, <I> "Chemkin Collection</I>, Release 3.6,
   Reaction Design, Inc., San Diego, CA (2000)
   
   [2] KEE, DIXON-LEWIS, WARNATZ, COLTRIN and MILLER, <I>"A Fortran Computer Code Package for the Evaluation of Gas-Phase Multicomponent Transport Porperties"</I>, Sandia National Laboratories Report SAND86-8246, (1986), pp. 44--50
  */
  template<class TPTDATATYPE, class IMPL>
  class GeneralTransport{
  public:
    typedef TPTDATATYPE DataType;
    typedef typename TPTDATATYPE::value_type value_type;
    typedef value_type* PointerType;
    typedef GeneralTransport<TPTDATATYPE,IMPL> ThisType;
    typedef std::size_t SizeType;

    //! molar mass of species k (in SI)
    value_type molar_mass (SizeType k){
      return DataType::molar_masses()[k];
    }

    value_type reduced_molar_mass (SizeType j, SizeType k){
      return molar_mass(j)*molar_mass(k)/(molar_mass(j) + molar_mass(k));
    }

    //general functions
    //! geometry (0 = atom, 1 = lin. molecule, 2 = nonlinear molecule) unitless
    SizeType geometry(SizeType k) {return static_cast<SizeType>(ptr_[6*k]);}
    //!Lennard-Jones potential well depth [=] K
    value_type epsilon_kB(SizeType k) {return ptr_[6*k+1];}
    //!Lennard-Jones collision diameter [=] Angström = 1.e-10m
    value_type sigma(SizeType k) {return ptr_[6*k+2];}
    //! dipole moment [=] Debye = 1.e-18statC·cm $\approx$ 3.33564e-30C·m
    value_type mu(SizeType k) {return ptr_[6*k+3];}
    //! polarizability [=] (Angström)³ = 1.e-30m
    value_type alpha(SizeType k) {return ptr_[6*k+4];}
    //! rotational relaxation collision number  uniless
    value_type Z_rot(SizeType k) {return  ptr_[6*k+5];}

    //! cf. [1], eq. (5)
    //! NOTE:  \f$\frac{\epsilon_{jk}}{k_B} = \sqrt(\frac{\epsilon_{j}}{k_B}\frac{\epsilon_{k}}{k_B}) is equivalent to \f$\epsilon_{jk} = \sqrt(\epsilon_{j}\epsilon_{k})\f$ since \f$k_B\f$ is canceled.
    value_type interaction_well_depth(SizeType j, SizeType k){
      return sqrt(epsilon_kB(j)*epsilon_kB(k));
    }

    //! alias which is closer to the CHEMKIN notation
    value_type epsilon(SizeType j, SizeType k){
      return interaction_well_depth(j,k);
    }

    template<class T> //!reduced temperature \f$ T^*\f$
    T reduced_temperature(SizeType k, const T& temp) {
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<".");
#endif
      return 1./epsilon_kB(k)*temp;
    }

    template<class T>
    T reduced_temperature(SizeType j, SizeType k, const T& temp){
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T()))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<".");
#endif
      // std::cout << "epsiolonkB( "<<j<<") = "<< epsilon_kB(j) << ",  epsilonkB("<<k<<") = " <<epsilon_kB(k)<<std::endl;
      return 1./interaction_well_depth(j,k)*temp;
}

    value_type sigma(SizeType j, SizeType k){
      return 0.5*(sigma(j) + sigma(k));
    }
    
   
    //! this is the special case of <TT>delta_tilde_star(i,j)</TT> with <TT>i = j</TT>
     value_type delta_tilde_star(SizeType k){
       return 0.5*ntimes<2>(mu(k))/(epsilon_kB(k)*ntimes<3>(sigma(k)));
     }


    //!UNITLESS term -- SI version, cf. [KEE, p. 496, eq. 12.18]
    value_type delta_tilde_star(SizeType i, SizeType j) {
      return 0.5*mu(i)*mu(j)/(epsilon(i,j)*ntimes<3>(sigma(i,j)));
    }

  
    //! approximation to \f$ \Omega^{(1,1)*}\left(T^*,\tilde \delta^*_{ij}\right), \f$ cf. [KEE, p. 496, eq. (12.19)]. Generalization for the Stockmayer potential.
    template<class T>
    T Omega_1_1_star_approximation(SizeType i, SizeType j, const T& Tstar){
#ifndef NDEBUG
      if((!is_well_defined_value(Tstar)) || (Tstar < T()))
	ADONIS_ERROR(ValueError,"Bad Tstar = "<< Tstar <<".");
#endif
      T d;
      (i != j) ? d = delta_tilde_star(i,j) : d = delta_tilde_star(i);
      return ( (1.0548*pow(Tstar,-0.15504) + pow((Tstar+0.55909),-2.1705))*(1 + (exp(0.093193/Tstar) - exp(-1.5/Tstar))*ntimes<2>(d)/(2 + 2.5*d)) );
    }
    

    //! cf. [KEE, p. 496, eq. 12.20]
    template<class T>
    T Omega_2_2_star_approximation(SizeType i, SizeType j, const T& Tstar){
#ifndef NDEBUG
      if((!is_well_defined_value(Tstar)) || (Tstar < T()))
	ADONIS_ERROR(ValueError,"Bad Tstar = "<< Tstar <<".");
#endif
      T d;
      //! in case of same indices use simpler formula \$f\delta_{k}\f$
      (i != j) ? d = delta_tilde_star(i,j) : d = delta_tilde_star(i);
      return ( (1.0413*pow(Tstar,-0.11930) + pow((Tstar+0.43628),-1.6041))*(1 + (exp(0.095661/Tstar) - exp(-2.0/Tstar))*ntimes<2>(d)/(2 + 2.5*d)) );
    }

    template<bool CHS, class T1, class T2>
    T2 bin_diff_coeff(SizeType i, SizeType j, const T1& p, const T2& temp) {
#ifndef NDEBUG
      if((!is_well_defined_value(temp)) || (temp < T2()) || (!is_well_defined_value(p)))
	ADONIS_ERROR(ValueError,"Bad temp = "<< temp <<" or p = "<< p << ".");
#endif
      //W should enter in g/mol, T in K, p in bar and sigma in Angström
      //DONE: inserted pressure in_bar(p), i.e. if molar masses are inserted in
      // g/mol, i.e. CHS=true, we also assume that pressure is stated in bar. 
      //In other words, to get a practical formula, the T, p, W and sigma 
      //must be input in K, bar, g/mol and Angström, respectively.
      return ChooseRightPrefactor<T2,CHS>::diffprefac()*sqrt(0.5*1./(ChooseRightPrefactor<T2,CHS>::in_gram(ThisType::reduced_molar_mass(i,j)))*ntimes<3>(temp))/(ChooseRightPrefactor<T2,CHS>::in_bar(p)*ntimes<2>(ThisType::sigma(i,j))*ThisType::Omega_1_1_star_approximation(i,j,ThisType::reduced_temperature(i,j,temp)));
    }  //ChooseRightPrefactor<T2,CHS>::in_bar(p)

  protected:
    PointerType ptr_;
   
  private:
    const IMPL& ref2derived() const {return static_cast<const IMPL&>(*this);}
    

    //!delegates to derived class
    inline void initialize(PointerType ptr){
      ref2derived().initialize(ptr);
    }

    inline PointerType result(){
      return ref2derived().coefficients();
    }
  
  };



}//end namespace 

#endif
