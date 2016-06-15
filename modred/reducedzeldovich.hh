#ifndef TEST_REDUCED_ZELDOVICH_MECHANISM_BY_HAND_HH
#define TEST_REDUCED_ZELDOVICH_MECHANISM_BY_HAND_HH

#include <cmath>

#include "../common/error.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../linalg/linearsystemsolvers.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../derivatives/jacobianofsource.hh"
#include "../io/readinparameters.hh"
#include "../templatemetaprograms/unrollloop.hh"

#include "../ode/examples/zeldovich.hh" //full model

namespace Adonis{

  template<class T>
  class ReducedZeldovich{
  public:
    
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<DType> VDType;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef ExprTmpl::MyVec<T> VType;
    typedef std::size_t SizeType;
    typedef Norm<'2',BaseType> NormType;

    typedef Zeldovich<DType> FullModelType;  //FULL MODEL
    typedef JacS<DType,Zeldovich,ExprTmpl::MyVec> JacobianType;//ITS JACOBIAN


    //n = #red + 1 (temp) 
    ReducedZeldovich(SizeType n = 3+1):fullode_(8),
				       maxIt_(13),     count_(1),
				       maxNewt_(5),
				       tol_(1.e-10),
				       newtTol_(1.e-08),
				       fulldim_(8),
				       rhs_(n),zM_(8), zM0_(8), C_(7*1), b_(1), residual_(1), z_nu_(8), g_(1), defect_(7), G_(8),Gprime_(8*8),up_(8), low_(8),zMchem_(7) {

      ParameterData PD;
      PD.read_from_file("/home/mfein/MARC++/modred/dat/zeldovich.dat");

      //======== TODO: change setting here====================================
      zM0_ <<= 
      	PD.get_datum<DType>("O"),  
      	PD.get_datum<DType>("N2"),
      	PD.get_datum<DType>("NO"),
      	PD.get_datum<DType>("N"),
      	PD.get_datum<DType>("O2"),
      	PD.get_datum<DType>("OH"),
      	PD.get_datum<DType>("H"),
      	PD.get_datum<DType>("Temperature");

      kn_ = PD.get_datum<DType>("hEstim");//PD.get_datum<DType>("hred"); //

      whichTimeStepper_ = PD.get_datum<char>("timestepper");
      reduced_cp_heat_ = PD.get_datum<bool>("reduced_cp_heat");
      //======================================================================
    
      low_ <<= 0.,0.,0.,0.,0.,0.,0.,PD.get_datum<DType>("Tmin");
      up_ <<=  1.,1.,1.,1.,1.,1.,1.,PD.get_datum<DType>("Tmax");
    
      C_ <<= 1.,1.,1.,1.,1.,1.,1.;   //only dim of chemical species
      b_ <<= 1.;

      zM_ = zM0_;
      
      
      NablaChem_.set_with_init(fulldim_,fulldim_,zM_); 


      // //chemistry stuff
      // coef_N2_.resize(14);
      // coef_NO_.resize(14);
      // coef_O2_.resize(14);
      coef_O_.resize(14);
      coef_N2_.resize(14);
      coef_NO_.resize(14);
      coef_N_.resize(14);
      coef_O2_.resize(14);
      coef_OH_.resize(14);
      coef_H_.resize(14);
      
        coef_O_ <<= 2.56942078E+00,-8.59741137E-05, 4.19484589E-08,-1.00177799E-11, 1.22833691E-15, 2.92175791E+04, 4.78433864E+00, 3.16826710E+00,-3.27931884E-03, 6.64306396E-06, -6.12806624E-09, 2.11265971E-12, 2.91222592E+04, 2.05193346E+00;
      coef_N2_ <<= 0.02926640E+02, 0.14879768E-02,-0.05684760E-05, 0.10097038E-09,-0.06753351E-13, -0.09227977E+04, 0.05980528E+02, 0.03298677E+02, 0.14082404E-02,-0.03963222E-04, 0.05641515E-07,-0.02444854E-10,-0.10208999E+04, 0.03950372E+02;
      coef_NO_ <<= 0.32606056E+01, 0.11911043E-02,-0.42917048E-06, 0.69457669E-10,-0.40336099E-14, 0.99209746E+04, 0.63693027E+01, 0.42184763E+01,-0.46389760E-02, 0.11041022E-04,  -0.93361354E-08, 0.28035770E-11, 0.98446230E+04, 0.22808464E+01;
      coef_N_ <<= 0.24159429E+01, 0.17489065E-03,-0.11902369E-06, 0.30226245E-10,-0.20360982E-14, 0.56133773E+05, 0.46496096E+01, 0.25000000E+01, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.56104637E+05, 0.41939087E+01;
      coef_O2_ <<=  3.28253784E+00, 1.48308754E-03,-7.57966669E-07, 2.09470555E-10,-2.16717794E-14, -1.08845772E+03, 5.45323129E+00, 3.78245636E+00,-2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12,-1.06394356E+03, 3.65767573E+00;

      coef_OH_ <<=  3.09288767E+00, 5.48429716E-04, 1.26505228E-07,-8.79461556E-11, 1.17412376E-14, 3.85865700E+03, 4.47669610E+00, 3.99201543E+00,-2.40131752E-03, 4.61793841E-06, -3.88113333E-09, 1.36411470E-12, 3.61508056E+03,-1.03925458E-01;
      coef_H_ <<= 2.50000001E+00,-2.30842973E-11, 1.61561948E-14,-4.73515235E-18, 4.98197357E-22, 2.54736599E+04,-4.46682914E-01, 2.50000000E+00, 7.05332819E-13,-1.99591964E-15, 2.30081632E-18,-9.27732332E-22, 2.54736599E+04,-4.46682853E-01; 
      // coef_N2_ <<= 0.02926640E+02, 0.14879768E-02,-0.05684760E-05, 0.10097038E-09,-0.06753351E-13, -0.09227977E+04, 0.05980528E+02, 0.03298677E+02, 0.14082404E-02,-0.03963222E-04, 0.05641515E-07,-0.02444854E-10,-0.10208999E+04, 0.03950372E+02;
      // coef_NO_ <<= 0.32606056E+01, 0.11911043E-02,-0.42917048E-06, 0.69457669E-10,-0.40336099E-14, 0.99209746E+04, 0.63693027E+01, 0.42184763E+01,-0.46389760E-02, 0.11041022E-04,  -0.93361354E-08, 0.28035770E-11, 0.98446230E+04, 0.22808464E+01;

      // coef_O2_ <<=  3.28253784E+00, 1.48308754E-03,-7.57966669E-07, 2.09470555E-10,-2.16717794E-14, -1.08845772E+03, 5.45323129E+00, 3.78245636E+00,-2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12,-1.06394356E+03, 3.65767573E+00;

      W_[O] = 0.016;   //kg/mol
      W_[N2] = 0.028;
      W_[NO] = 0.030;
      W_[N] = 0.014;
      W_[O2] = 0.032;
      W_[OH] = 0.017;
      W_[H] = 0.001;

      cp_ = heat_ = wbar_ =  T();
    
      p0_ = PD.get_datum<DType>("p0");
      
      
    }



    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}


    template<class X>  //red = 3 +1 
    VType& operator()(const X& red){
      //[I.] SPECIES RECONSTRUCTION
      zM_[N2] = red[0];
      zM_[NO] = red[1];
      zM_[O2] = red[2];

      zM_[TMP] = red[3]; //temperature

      // 1.) satisfy linear(ized) constraints
      //! assign chemistry part only 
      UnrollLoop<0,7>::assign(zMchem_.begin(),zM_.begin());
      std::cout << "zMchem_ = "<< zMchem_ << std::endl;


      //! assign chemistry to zChem_ -- now only chemistry is considered
      //UnrollLoop<0,7>::assign(zMchem_.begin(),zM_.begin()); 
      
      residual_ = matrix_vector_product(C_,zMchem_) - b_;
      
      bool isperturbed(false);

      if(NormType::norm(residual_) > 1.e-09){
	//!perform "Gauss-Newton" iteration
	//zMchem_nu_= zMchem_;
	int m = (int)b_.size();
	for(SizeType nu = 1; nu <= maxIt_; ++nu){
	  count_ = nu;
	  //g_ = -matrix_vector_product(C_,zMchem_nu_) + b_;
	  g_ = -matrix_vector_product(C_,zMchem_) + b_;
	  defect_ = solve_rank_deficient_ls(C_,m,g_,1); 
	  //zMchem_nu_ += defect_;
	  zMchem_ += defect_;
	  if(NormType::norm(defect_) <= tol_)
	    break;  //leave loop since convergence has been achieved
	}
	if(count_ == maxIt_)
	  ADONIS_INFO(Information, "Max. number of iterations reached (maxit = "<< maxIt_ << ").\n   No convergence with tol = "<< tol_ << ".");
      
	//zMChem_ = zMchem_nu_;
	isperturbed = true;
      }
      //std::cout << "zMchem_ = "<< zMchem_ << std::endl;
      
      //! 2.) Project on bounds
      // assign back to zM_
      if(isperturbed)
	UnrollLoop<0,7>::assign(zM_.begin(),zMchem_.begin());
      zM_ = Max(low_,Min(zM_,up_));

      std::cout << "zM (proj. on bounds) = "<< zM_ << std::endl;

      //! 3.) Perform 1 implicit or expl. Euler time step
      //! note that we use the full model, i.e. the full $\dot{\omega}$ here!!,
      //! that is it is principally available for temperature evaluation of 
      //! reduced temperature
      if(whichTimeStepper_ == 'i' ||  whichTimeStepper_ == 'I'){
      	std::cout << "IMPL. EULER" << std::endl;
      	z_prev_ = zM_;
      	z_nu_ = z_prev_; 
      	count_ = 1; //reset
      	for(SizeType l = 1; l <= maxNewt_; ++l){
      	  G_ = -(z_nu_ - z_prev_ - kn_*fullode_(z_nu_,fullode_.dot_omega()));
      	  Gprime_ = -kn_*NablaChem_.jacobian(z_nu_);
      	  //print_in_matrix_style(Gprime_,fulldim_);
      	  update_diagonal<AddBasicElements>(Gprime_,fulldim_,1.);
      	  good_square_solve(Gprime_,fulldim_,G_,1);
      	  z_nu_ += G_;

      	  if(NormType::norm(G_) <= newtTol_){
      	    //---- ASSIGN -------
      	    zM_ = z_nu_;
      	    //---------------------
      	    break;
      	  }

      	  if(l == maxNewt_){
      	    ADONIS_INFO(Information, "Maximum number of Newton iterations ("<<maxNewt_<<") reached. Decrease kn_....");
      	    l = 0; //set back i and repeat with kn/2 
      	    kn_ *= 0.5;
      	  }

      	  adonis_assert(Abs(kn_) >= 1.e-12); //don't let kn_ become too small

      	  count_++;

      	  if(count_ == 4*maxNewt_){
      	    ADONIS_INFO(Information, "After 4 iterations with decreasing stepsize, there hasn't been encountered any convergence...");
      	    break;  //o.k. 4 times with decreasing step size yielded no result
      	  }
      	} //end NEWTON iter
      }
      
      else if (whichTimeStepper_ == 'e' || whichTimeStepper_ == 'E'){
      	//perform explicit EULER step
      	std::cout << "EX. EULER" << std::endl;
      	zM_ += kn_*fullode_(zM_,fullode_.dot_omega());
      }
      else
      	ADONIS_ERROR(NotImplementedError,"Time stepping method whichTimeStepper_ = "<<whichTimeStepper_ << " undefined.");

      std::cout << "zM (time step) = "<< zM_ << std::endl;

      // [II.] REDUCED SOURCE
      // //! 4.) reduced source terms
     
      T k1  = 1.8e+08*exp(-38370./zM_[TMP]),
      	km1 = 3.8e+07*exp(-425./zM_[TMP]),
      	k2  = 1.8e+04*zM_[TMP]*exp(-4680./zM_[TMP]),
      	km2 = 3.8e+03*zM_[TMP]*exp(-20820./zM_[TMP]),
      	k3  = 7.1e+07*exp(-450./zM_[TMP]),
      	km3 = 1.7e+08*exp(-24560./zM_[TMP]);
      
      
      //zM_[H] = 1. - (zM_[O] + zM_[N2] + zM_[NO] + zM_[N] + zM_[O2] + zM_[OH]);

      //compute in concentrations first
      wbar_ = 0.;  //reset
      
      for(int k = 0; k < 7; ++k) //full concentrations means full wbar!!
	wbar_ += zM_[k]/W_[k];
      wbar_ = 1./wbar_;
      
      rho_ = (p0_*wbar_)/(PhysicalConstants<DType>::Rgas*zM_[TMP]);

      T invrho = 1./rho_; //used more than once
      
      //! full dimension needed here
      for(int k = 0; k < 7; ++k) //compute concentrations
      	Conc_[k] = rho_*zM_[k]/W_[k];

      std::cout << "Conc_ = "; print_all(Conc_,Conc_+7);

      //N2
      dotomega_[0] = -k1*Conc_[O]*Conc_[N2] + km1*Conc_[NO]*Conc_[N]; 
      //NO
      dotomega_[1] = k1*Conc_[O]*Conc_[N2] - km1*Conc_[NO]*Conc_[N] + k2*Conc_[N]*Conc_[O2] 
      	              - km2*Conc_[NO]*Conc_[O] + k3*Conc_[N]*Conc_[OH] - km3*Conc_[NO]*Conc_[H];
      //O2
      dotomega_[2] = -k2*Conc_[N]*Conc_[O2] + km2*Conc_[NO]*Conc_[O];

      //!SPECIES in mass fractions
      //! compute in mass fractions, i.e. \f$1/\rho\dot{\omega}W_k\f$
      rhs_[0] = invrho*(dotomega_[0]*W_[N2]);
      rhs_[1] = invrho*(dotomega_[1]*W_[NO]);
      rhs_[2]  = invrho*(dotomega_[2]*W_[O2]);

      //! Temperature -- do I need full species here??
      if( reduced_cp_heat_){
	cp_ = heat_capacity(zM_[TMP],&coef_N2_[0])/W_[N2]*zM_[N2] +
	  heat_capacity(zM_[TMP],&coef_NO_[0])/W_[NO]*zM_[NO] +
	  heat_capacity(zM_[TMP],&coef_O2_[0])/W_[O2]*zM_[O2];
      
      
	heat_ = enthalpy(zM_[TMP],&coef_N2_[0])*dotomega_[0] +
	  enthalpy(zM_[TMP],&coef_NO_[0])*dotomega_[1] +
	  enthalpy(zM_[TMP],&coef_O2_[0])*dotomega_[2];
      }
      else{
	cp_ = heat_capacity(zM_[TMP],&coef_O_[0])/W_[O]*zM_[O] +
	  heat_capacity(zM_[TMP],&coef_N2_[0])/W_[N2]*zM_[N2] +
	  heat_capacity(zM_[TMP],&coef_NO_[0])/W_[NO]*zM_[NO] +
	  heat_capacity(zM_[TMP],&coef_N_[0])/W_[N]*zM_[N]    +
	  heat_capacity(zM_[TMP],&coef_O2_[0])/W_[O2]*zM_[O2] + 
	  heat_capacity(zM_[TMP],&coef_OH_[0])/W_[OH]*zM_[OH] +
	  heat_capacity(zM_[TMP],&coef_H_[0])/W_[H]*zM_[H];
	

     
	heat_ = enthalpy(zM_[TMP],&coef_O_[0])*fullode_.dot_omega()[O] +
	  enthalpy(zM_[TMP],&coef_N2_[0])*dotomega_[0] + //fullode_.dot_omega()[N2] +
	  enthalpy(zM_[TMP],&coef_NO_[0])*dotomega_[1] + //fullode_.dot_omega()[NO] +
	  enthalpy(zM_[TMP],&coef_N_[0])*fullode_.dot_omega()[N]   +
	  enthalpy(zM_[TMP],&coef_O2_[0])*dotomega_[2] + //fullode_.dot_omega()[O2]   +
	  enthalpy(zM_[TMP],&coef_OH_[0])*fullode_.dot_omega()[OH] +
	  enthalpy(zM_[TMP],&coef_H_[0])*fullode_.dot_omega()[H];
      }


      rhs_[3] = -invrho*1./cp_*heat_; //mind the minus here!

      return rhs_;

    }

  private:
    FullModelType fullode_;
    SizeType maxIt_, count_, maxNewt_;
    BaseType tol_, newtTol_;
    
    SizeType fulldim_;
    VType rhs_, zM_, zM0_, C_, b_,residual_, z_nu_, g_, defect_, G_,Gprime_,
	  up_, low_, z_prev_;
    //only chemisty -- note that these are usually of nchem
    VType zMchem_, zMchem_nu_;
    value_type kn_;
    char whichTimeStepper_;
    bool reduced_cp_heat_;
    JacobianType NablaChem_;

    T p0_,rho_;
    VDType coef_O_,coef_N2_,coef_NO_,coef_N_, coef_O2_,coef_OH_,coef_H_; //coef_N2_,coef_NO_,coef_O2_;
    DType W_[7];
    T cp_, heat_, wbar_;
    T dotomega_[3]; //only species
    T Conc_[7];       //only species


    enum{O,N2,NO,N,O2,OH,H,TMP};
  };

} //end namespace 

#endif
