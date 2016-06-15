#ifndef THERMAL_ZELDOVICH_MECHANISM_HH
#define THERMAL_ZELDOVICH_MECHANISM_HH

#include <cmath>

#include "../../expressiontemplates/exprvec.hh"
#include "../../massactionkinetics/thermochemistry.hh"
#include "../../massactionkinetics/physicalconstants.hh"
#include "../../templatemetaprograms/unrollloop.hh"
#include "../../io/readinparameters.hh"

namespace Adonis{


  /**
   * \brief some external reaction rates 
   */
  template<class T>
  class ExternalReactionRates{  //zeldovich mechanism
  public:
   
    static inline T k1(const T& temp) {return 1.8e+08*exp(-38370./temp);}
    static inline T km1(const T& temp){return 3.8e+07*exp(-425./temp);}
    static inline T k2(const T& temp){return 1.8e+04*temp*exp(-4680./temp);}
    static inline T km2(const T& temp){return 3.8e+03*temp*exp(-20820./temp);}
    static inline T k3(const T& temp){return 7.1e+07*exp(-450./temp);}
    static inline T km3(const T& temp){return 1.7e+08*exp(-24560./temp);}
  };


  /**
   * Thermal Zeldovich mechanism. The mechanism can be found in [1, p. 556].
   * Temperature-dependent Arrhenius reaction rates on the web [2].
   * Reaction rates can also be found in [3] and I prefer these.

   * Thermal energy is included as last entry in the source term
   *
   * [1] [KEE, "Chemically reacting flows", John Wiley & Sons, 2006]
   * [2] <a href="http://www.jullio.pe.kr/fluent6.1/help/html/ug/node624.htm"> coefficients </a>
   * [3] [BLAUCH et al., "Evaluated Kinetic Data for Combustion Modelling Supplement I", J. Phys. Chem. Ref. Data, 1995]
   */
  template<class T>
  class Zeldovich{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef std::size_t SizeType; 
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<DType> VDType;

    typedef ExternalReactionRates<T> RR;

    Zeldovich(SizeType n =0):rhs_(n),
			     p0_(101325),  //standard atmosphere
			     rho_(0.575), //just default
			     prescale_(1.e-06),  //or 1 for cm^3 /(molecule sec)
			     coef_O_(14),coef_N2_(14),coef_NO_(14),coef_N_(14), coef_O2_(14),coef_OH_(14),coef_H_(14){
      coef_O_ <<= 2.56942078E+00,-8.59741137E-05, 4.19484589E-08,-1.00177799E-11, 1.22833691E-15, 2.92175791E+04, 4.78433864E+00, 3.16826710E+00,-3.27931884E-03, 6.64306396E-06, -6.12806624E-09, 2.11265971E-12, 2.91222592E+04, 2.05193346E+00;
      coef_N2_ <<= 0.02926640E+02, 0.14879768E-02,-0.05684760E-05, 0.10097038E-09,-0.06753351E-13, -0.09227977E+04, 0.05980528E+02, 0.03298677E+02, 0.14082404E-02,-0.03963222E-04, 0.05641515E-07,-0.02444854E-10,-0.10208999E+04, 0.03950372E+02;
      coef_NO_ <<= 0.32606056E+01, 0.11911043E-02,-0.42917048E-06, 0.69457669E-10,-0.40336099E-14, 0.99209746E+04, 0.63693027E+01, 0.42184763E+01,-0.46389760E-02, 0.11041022E-04,  -0.93361354E-08, 0.28035770E-11, 0.98446230E+04, 0.22808464E+01;
      coef_N_ <<= 0.24159429E+01, 0.17489065E-03,-0.11902369E-06, 0.30226245E-10,-0.20360982E-14, 0.56133773E+05, 0.46496096E+01, 0.25000000E+01, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.56104637E+05, 0.41939087E+01;
      coef_O2_ <<=  3.28253784E+00, 1.48308754E-03,-7.57966669E-07, 2.09470555E-10,-2.16717794E-14, -1.08845772E+03, 5.45323129E+00, 3.78245636E+00,-2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12,-1.06394356E+03, 3.65767573E+00;

      coef_OH_ <<=  3.09288767E+00, 5.48429716E-04, 1.26505228E-07,-8.79461556E-11, 1.17412376E-14, 3.85865700E+03, 4.47669610E+00, 3.99201543E+00,-2.40131752E-03, 4.61793841E-06, -3.88113333E-09, 1.36411470E-12, 3.61508056E+03,-1.03925458E-01;
      coef_H_ <<= 2.50000001E+00,-2.30842973E-11, 1.61561948E-14,-4.73515235E-18, 4.98197357E-22, 2.54736599E+04,-4.46682914E-01, 2.50000000E+00, 7.05332819E-13,-1.99591964E-15, 2.30081632E-18,-9.27732332E-22, 2.54736599E+04,-4.46682853E-01; 
      
      //!molar masses of O,N2,NO,N,O2,OH,H
      W_[O] = 0.016;   //kg/mol
      W_[N2] = 0.028;
      W_[NO] = 0.030;
      W_[N] = 0.014;
      W_[O2] = 0.032;
      W_[OH] = 0.017;
      W_[H] = 0.001;

      cp_ = heat_ = wbar_ = T();
      ParameterData PD;
      //!use absolute path here!
      PD.read_from_file("/home/mfein/MARC++/modred/dat/zeldovich.dat"); //from wher
      prescale_ = PD.get_datum<DType>("prescale");
      p0_ = PD.get_datum<DType>("p0");
}


    SizeType dim() const {return 8;}   //6 species + temperature
    SizeType domain_dim() const{return 8;}

    const T* dot_omega() const {return &dotomega_[0];} 
    T* dot_omega() {return &dotomega_[0];} 

    template<class X>
    VType& operator()(const X& vars){
      //!Hard copy
      UnrollLoop<0,8>::assign(&y_[0],vars.begin());
      
      
      //! from [2] already in SI, viz. m³/(gmol s)
      //! (note that a gmol is exactly the same as mol!!)
      
      T k1  = 1.8e+08*exp(-38370./y_[TMP]),
      	km1 = 3.8e+07*exp(-425./y_[TMP]),
      	k2  = 1.8e+04*y_[TMP]*exp(-4680./y_[TMP]),
      	km2 = 3.8e+03*y_[TMP]*exp(-20820./y_[TMP]),
      	k3  = 7.1e+07*exp(-450./y_[TMP]),
      	km3 = 1.7e+08*exp(-24560./y_[TMP]);
     
      //!rates are stated in cm³/(molecule sec), cf. [3], hence for SI we must
      //!premultiply them by $10^{-6}$ (all reactions are bimolecular)
      
      //!from [3]
      // T k1 =  prescale_*3.0e-10*exp(-38300./y_[TMP]),
      // 	km1 = prescale_*7.1e-11*exp(-790./y_[TMP]),
      // 	k2 =  prescale_*1.5e-14*y_[TMP]*exp(-3270./y_[TMP]),
      // 	km2 = prescale_*1.14e-15*pow(y_[TMP],1.13)*exp(-19200./y_[TMP]),
      // 	k3 =  prescale_*4.7e-11,
      // 	km3 = prescale_*3.6e-10*exp(-24910./y_[TMP]);

    



      //grant mass conservation 
      y_[H] = 1. - (y_[O] + y_[N2] + y_[NO] + y_[N] + y_[O2] + y_[OH]);

      //compute in concentrations first
      wbar_ = 0.;  //reset
      for(int k = 0; k < 7; ++k)  
      	wbar_ += y_[k]/W_[k];
      wbar_ = 1./wbar_;
      
      //std::cout << "wbar = "<<wbar_ << std::endl;


      rho_ = (p0_*wbar_)/(PhysicalConstants<DType>::Rgas*y_[TMP]);
      //std::cout << "rho = "<<rho_ << std::endl;

      T invrho = 1./rho_; //used more than once

      for(int k = 0; k < 7; ++k) //compute concentrations
	C_[k] = rho_*y_[k]/W_[k];

      //std::cout << "concentrations = "; print_all(C_,C_+7); 

      //!compute \f$\dot{\omega}\f$ first 
      dotomega_[O] = -k1*C_[O]*C_[N2] + km1*C_[NO]*C_[N] + k2*C_[N]*C_[O2] 
	             - km2*C_[NO]*C_[O];
      dotomega_[N2] = -k1*C_[O]*C_[N2] + km1*C_[NO]*C_[N];
      
      dotomega_[NO] = k1*C_[O]*C_[N2] - km1*C_[NO]*C_[N] + k2*C_[N]*C_[O2] 
	              - km2*C_[NO]*C_[O] + k3*C_[N]*C_[OH] - km3*C_[NO]*C_[H];
      
      dotomega_[N] = k1*C_[O]*C_[N2] - km1*C_[NO]*C_[N] - k2*C_[N]*C_[O2] 
	             + km2*C_[NO]*C_[O] - k3*C_[N]*C_[OH] + km3*C_[NO]*C_[H];

      dotomega_[O2] = -k2*C_[N]*C_[O2] + km2*C_[NO]*C_[O];

      dotomega_[OH] = -k3*C_[N]*C_[OH] + km3*C_[NO]*C_[H];

      dotomega_[H] = k3*C_[N]*C_[OH] - km3*C_[NO]*C_[H];

      //std::cout << "dotomega = "; print_all(dotomega_,dotomega_+7);

      //!SPECIES in mass fractions
      //! compute in mass fractions, i.e. \f$1/\rho\dot{\omega}W_k\f$
      rhs_[O]  = invrho*(dotomega_[O]*W_[O]);
      rhs_[N2] = invrho*(dotomega_[N2]*W_[N2]);
      rhs_[NO] = invrho*(dotomega_[NO]*W_[NO]);
      rhs_[N] = invrho*(dotomega_[N]*W_[N]);
      rhs_[O2]  = invrho*(dotomega_[O2]*W_[O2]);
      rhs_[OH] = invrho*(dotomega_[OH]*W_[OH]);
      rhs_[H] =  invrho*(dotomega_[H]*W_[H]);

      //TEMPERATURE
      //! compute mixture heat capacity -- only consider high temperatures
      //! in mass fractions here!
      cp_ = heat_capacity(y_[TMP],&coef_O_[0])/W_[O]*y_[O] +
	heat_capacity(y_[TMP],&coef_N2_[0])/W_[N2]*y_[N2] +
	heat_capacity(y_[TMP],&coef_NO_[0])/W_[NO]*y_[NO] +
	heat_capacity(y_[TMP],&coef_N_[0])/W_[N]*y_[N]    +
	heat_capacity(y_[TMP],&coef_O2_[0])/W_[O2]*y_[O2] + 
	heat_capacity(y_[TMP],&coef_OH_[0])/W_[OH]*y_[OH] +
	heat_capacity(y_[TMP],&coef_H_[0])/W_[H]*y_[H];
	

      //std::cout << "cp = "<< cp_ << std::endl;

      heat_ = enthalpy(y_[TMP],&coef_O_[0])*dotomega_[O] +
	enthalpy(y_[TMP],&coef_N2_[0])*dotomega_[N2] +
	enthalpy(y_[TMP],&coef_NO_[0])*dotomega_[NO] +
	enthalpy(y_[TMP],&coef_N_[0])*dotomega_[N]   +
	enthalpy(y_[TMP],&coef_O2_[0])*dotomega_[O2]   +
	enthalpy(y_[TMP],&coef_OH_[0])*dotomega_[OH] +
	enthalpy(y_[TMP],&coef_H_[0])*dotomega_[H];
      
      rhs_[TMP] = -invrho*1./cp_*heat_; //mind the minus here!

      return rhs_;
    }

    

    /**
     * Use the same as above but dotomega is used as iterator for later reusage
     */
    template<class X, class ITER>
    VType& operator()(const X& vars, ITER dot){
        //!Hard copy
      UnrollLoop<0,8>::assign(&y_[0],vars.begin());
      
      
      //! from [2] already in SI, viz. m³/(gmol s)
      //! (note that a gmol is exactly the same as mol!!)
      
      T k1  = 1.8e+08*exp(-38370./y_[TMP]),
      	km1 = 3.8e+07*exp(-425./y_[TMP]),
      	k2  = 1.8e+04*y_[TMP]*exp(-4680./y_[TMP]),
      	km2 = 3.8e+03*y_[TMP]*exp(-20820./y_[TMP]),
      	k3  = 7.1e+07*exp(-450./y_[TMP]),
      	km3 = 1.7e+08*exp(-24560./y_[TMP]);
     
      //!rates are stated in cm³/(molecule sec), cf. [3], hence for SI we must
      //!premultiply them by $10^{-6}$ (all reactions are bimolecular)
      
      //!from [3]
      // T k1 =  prescale_*3.0e-10*exp(-38300./y_[TMP]),
      // 	km1 = prescale_*7.1e-11*exp(-790./y_[TMP]),
      // 	k2 =  prescale_*1.5e-14*y_[TMP]*exp(-3270./y_[TMP]),
      // 	km2 = prescale_*1.14e-15*pow(y_[TMP],1.13)*exp(-19200./y_[TMP]),
      // 	k3 =  prescale_*4.7e-11,
      // 	km3 = prescale_*3.6e-10*exp(-24910./y_[TMP]);

    



      //grant mass conservation 
      y_[H] = 1. - (y_[O] + y_[N2] + y_[NO] + y_[N] + y_[O2] + y_[OH]);

      //compute in concentrations first
      wbar_ = 0.;  //reset
      for(int k = 0; k < 7; ++k)  
      	wbar_ += y_[k]/W_[k];
      wbar_ = 1./wbar_;
      
      //std::cout << "wbar = "<<wbar_ << std::endl;


      rho_ = (p0_*wbar_)/(PhysicalConstants<DType>::Rgas*y_[TMP]);
      //std::cout << "rho = "<<rho_ << std::endl;

      T invrho = 1./rho_; //used more than once

      for(int k = 0; k < 7; ++k) //compute concentrations
	C_[k] = rho_*y_[k]/W_[k];

      //std::cout << "concentrations = "; print_all(C_,C_+7); 

      //!compute \f$\dot{\omega}\f$ first 
      dot[O] = -k1*C_[O]*C_[N2] + km1*C_[NO]*C_[N] + k2*C_[N]*C_[O2] 
	             - km2*C_[NO]*C_[O];
      dot[N2] = -k1*C_[O]*C_[N2] + km1*C_[NO]*C_[N];
      
      dot[NO] = k1*C_[O]*C_[N2] - km1*C_[NO]*C_[N] + k2*C_[N]*C_[O2] 
	              - km2*C_[NO]*C_[O] + k3*C_[N]*C_[OH] - km3*C_[NO]*C_[H];
      
      dot[N] = k1*C_[O]*C_[N2] - km1*C_[NO]*C_[N] - k2*C_[N]*C_[O2] 
	             + km2*C_[NO]*C_[O] - k3*C_[N]*C_[OH] + km3*C_[NO]*C_[H];

      dot[O2] = -k2*C_[N]*C_[O2] + km2*C_[NO]*C_[O];

      dot[OH] = -k3*C_[N]*C_[OH] + km3*C_[NO]*C_[H];

      dot[H] = k3*C_[N]*C_[OH] - km3*C_[NO]*C_[H];

      //std::cout << "dotomega = "; print_all(dotomega_,dotomega_+7);

      //!SPECIES in mass fractions
      //! compute in mass fractions, i.e. \f$1/\rho\dot{\omega}W_k\f$
      rhs_[O]  = invrho*(dot[O]*W_[O]);
      rhs_[N2] = invrho*(dot[N2]*W_[N2]);
      rhs_[NO] = invrho*(dot[NO]*W_[NO]);
      rhs_[N] = invrho*(dot[N]*W_[N]);
      rhs_[O2]  = invrho*(dot[O2]*W_[O2]);
      rhs_[OH] = invrho*(dot[OH]*W_[OH]);
      rhs_[H] =  invrho*(dot[H]*W_[H]);

      //TEMPERATURE
      //! compute mixture heat capacity -- only consider high temperatures
      //! in mass fractions here!
      cp_ = heat_capacity(y_[TMP],&coef_O_[0])/W_[O]*y_[O] +
	heat_capacity(y_[TMP],&coef_N2_[0])/W_[N2]*y_[N2] +
	heat_capacity(y_[TMP],&coef_NO_[0])/W_[NO]*y_[NO] +
	heat_capacity(y_[TMP],&coef_N_[0])/W_[N]*y_[N]    +
	heat_capacity(y_[TMP],&coef_O2_[0])/W_[O2]*y_[O2] + 
	heat_capacity(y_[TMP],&coef_OH_[0])/W_[OH]*y_[OH] +
	heat_capacity(y_[TMP],&coef_H_[0])/W_[H]*y_[H];
	

      //std::cout << "cp = "<< cp_ << std::endl;

      heat_ = enthalpy(y_[TMP],&coef_O_[0])*dot[O] +
	enthalpy(y_[TMP],&coef_N2_[0])*dot[N2] +
	enthalpy(y_[TMP],&coef_NO_[0])*dot[NO] +
	enthalpy(y_[TMP],&coef_N_[0])*dot[N]   +
	enthalpy(y_[TMP],&coef_O2_[0])*dot[O2]   +
	enthalpy(y_[TMP],&coef_OH_[0])*dot[OH] +
	enthalpy(y_[TMP],&coef_H_[0])*dot[H];
      
      rhs_[TMP] = -invrho*1./cp_*heat_; //mind the minus here!

      return rhs_;
    }


  private:
    VType rhs_;
    T p0_,rho_;
    DType prescale_;
    VDType coef_O_,coef_N2_,coef_NO_,coef_N_, coef_O2_,coef_OH_,coef_H_;
    DType W_[7];
    T cp_, heat_, wbar_;
    T dotomega_[7]; //only species
    T C_[7];       //only species
    T y_[8]; //for hard copy (species + temp)
   
    enum{O,N2,NO,N,O2,OH,H,TMP};

    

  };

}

#endif
