#ifndef REDUCED_2D_LIGHT_NAVIER_STOKES_HH
#define REDUCED_2D_LIGHT_NAVIER_STOKES_HH

/**
 * This functor can be used by an ODE integrator in the usual way
*/
#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"
#include "../common/typeadapter.hh"
#include "../graphics/printmoldata.hh"
#include "../io/readinparameters.hh"
#include "../templatemetaprograms/unrollloop.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../ode/ode.hh"
#include "../ode/additional/boundarytemperaturedistribution.hh"
#include "../ode/additional/parabolicvelocityprofile.hh"
#include "../containers/staticarray.hh"

//thermochemistry
#include "../massactionkinetics/reactionrates.hh"
#include "../massactionkinetics/stoichiometry.hh"
#include "../massactionkinetics/thermochemistry.hh"
#include "../massactionkinetics/indexinjector.hh"
#include "../massactionkinetics/buildchemicalrhs.hh"
#include "../massactionkinetics/physicalconstants.hh"
#include "../massactionkinetics/eos.hh"
#include "../massactionkinetics/data/thermochemicaldata.hh"
#include "../massactionkinetics/auxiliary.hh"

#include "../moleculartransport/fancytransporttmps.hh"
#include "../moleculartransport/conductivity.hh"
#include "../moleculartransport/diffusion.hh"

#include "choose.hh"
#include "boundary.hh"

#include "reconstructspecies.hh" //===RED
#include "../templatemetaprograms/functionalities4containersandscalars.hh"

namespace Adonis{

  template<class T>
  class ReducedLightNSE{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType; //===RED
    typedef typename TypeAdapter<T>::BaseType BaseType; //for norms, tolerances 
    typedef std::size_t SizeType;
    
     //needed for thermo chemistry
    typedef ThermoData4Mechanism<BaseType,ChooseFDSetting::MECH> DataType;  

    typedef Boundary2D<VType,ChooseFDSetting::MECH> BoundaryType;
    
    //===RED: replace all occurrences of ChooseFDSetting::ISREPLACED
    enum{nchem = SystemDimension<DataType::nspec,DataType::rednspec,true>::Value, 
	 //T + specs
	 nprim = 1+SystemDimension<DataType::nspec,DataType::rednspec,true>::Value,  //1+nchem
	 TMP = 0,          //first entry is temperature
	 SPC = 1          //0 is temperature, from here species start
    };
   
   
    //===RED
    typedef ExprTmpl::MyVec<ExprTmpl::MyVec<T> > TabulationType;
    typedef StaticArray<T,DataType::nspec> FullChemType;
    typedef StaticArray<T,DataType::rednspec> RedChemType;//same as ArrayType
    typedef ExprTmpl::MyVec<DType> VDType;


    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;

    //Other Index
    typedef IndexHandler<DataType,true> AccessType;  //===RED
    typedef typename AccessType::IndexerType IndexerType;

    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,IndexerType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    typedef EquationOfState<'i'> EosType;

    //TRANSPORT COEFFICIENTS
    //NO VISCOSITY NEEDED IN THIS CASE -- +++++FULL+++++ version here anyway
    typedef MixtureAveragedConductivity<DataType,false> ConductivityType;
    typedef MixtureAveragedDiffusion<DataType,false> DiffusionType;

    //====RED: with non-Cppad type since we solve system via LAPACK and integr.
    typedef ReconstructSpecies<DType,DataType::nspec,DataType::rednspec,ChooseFDSetting::MECH> ReconstructorType;
  
    ReducedLightNSE(SizeType dim = 0):rhs_(dim){
      //!construct objects needed for source term
      //! CAUTION: use <TT>static</TT> to maintain each objects' location 
      //!   beyond the call of the constuctor (i.e. beyond the localness of 
      //!   the {...}-block
      //! NOTE: the data are assumed to be stored for the FULL model
      static ThermoType nasa(DataType::nspec,DataType::thermo(),DataType::temperature_bounds());
       
      static StoichType stoichmatrix(DataType::nspec,DataType::nreac,DataType::stoichiometric_matrix());

      static TroeIndex TIndex;
      TIndex.create(DataType::nreac,DataType::troewtb());
	
      //! everything stated as  pointers	
      static FwdType forwardRates(DataType::nspec,DataType::nreac,DataType::arrhenius_values(),DataType::temperature_bounds(), DataType::ntroereac, DataType::troewtb(),DataType::collision_efficiencies(),&TIndex);

      //! everything stated as  pointers	
      static RevType reverseRates(&forwardRates,&stoichmatrix,&nasa);

      //!now either a reduced or the usual index is used for constructing
      //! \f$ \dt{\omega}\f$
      //static CommonIndexInjector Ixer(DataType::nspec,DataType::nreac);
      static IndexerType Ixer;
      AccessType::initialize(Ixer);

      BCS_.initialize(&forwardRates,&reverseRates,&Ixer);

      //READ IN DATA
      ParameterData PD;
      PD.read_from_file(ChooseFDSetting::FILE2READIN);
       //! domain [a,b] x [c,d]
      a_ = PD.get_datum<BaseType>("a");
      b_ = PD.get_datum<BaseType>("b");
      c_ = PD.get_datum<BaseType>("c");
      d_ = PD.get_datum<BaseType>("d");

      nx_ = PD.get_datum<int>("Nx");
      ny_ = PD.get_datum<int>("Ny");
      npt_ = nx_*ny_; //all points, i.e. \f$\Omega_h \cup \partial\Omega_h\f$
      totpoints_ = nprim*npt_;
      hx_ = (b_-a_)/(nx_-1);
      diam_ = d_-c_;
      hy_ = diam_/(ny_-1); 
      halfdiam_ = 0.5*diam_;
      p0_ = PD.get_datum<BaseType>("pconst");
      v1_ = PD.get_datum<BaseType>("v1_in");
      Tinlet_ = PD.get_datum<BaseType>("T_in");
      Tignition_ = PD.get_datum<BaseType>("T_ignition");

      bdytemp_.initialize(PD.get_datum<BaseType>("T_min"),Tignition_,PD.get_datum<BaseType>("hottest_x"));

      //transport coefficients are ready for computation
      //!no viscosity needed
      Mac_.initialize(DataType::transport());
      Mad_.initialize(DataType::transport());  //(i,j)    //5 pt. stencil
      Mad0_.initialize(DataType::transport()); //(i,j-1)
      Mad1_.initialize(DataType::transport()); //(i+1,j)
      Mad2_.initialize(DataType::transport()); //(i,j+1)
      Mad3_.initialize(DataType::transport()); //(i-1,j)

      bdy2d_.initialize(DataType::nspec); //always of full dimension!
      
      veloprofile_.resize(ny_);   //for the velocity profile (along y direction)
      u_.resize(totpoints_);

      cp_ = rho_ = T();  //initialize

      //constant parabolic velocity profile
      BaseType radius(0.);
      for(int j = 0; j < ny_; ++j){
	veloprofile_[j] = zero(parabolic_inlet_velocity(radius,diam_,v1_));
	radius += hy_;
	//std::cout << " "<<j <<".)     v_[0] = "<< veloprofile_[j] << std::endl;
      }
	
      v_[1] = 0.; //v2 is always zero
    
       cp_ = rho_ = T();  //initialize

      //===RED=================================
       Tab_.resize(nx_*ny_); //store zM at each grid point 
       for(int i = 0; i < nx_; ++i){  //intial values for full settings
	 DType x = a_ + i*hx_; 
	 for(int j = 0; j < ny_; ++j){
	   (Tab_[GRID(i,j)]).resize(DataType::nspec+1);
	 
	   Tab_[GRID(i,j)][0] = 0.;
	   Tab_[GRID(i,j)][1] = 0.8;
	   Tab_[GRID(i,j)][2] = 0.2; //= bdy2d_[k];
	   //}
	   //std::cout << std::endl;
	   Tab_[GRID(i,j)][DataType::nspec] = Tinlet_;
	   if(j == ny_-1)
	     Tab_[GRID(i,j)][DataType::nspec] = bdytemp_(x);
	 }
       }
  //! INITIALIZE species Reconstructor
      //!some guess
      VDType guessSpecPlusTemp(DataType::nspec+1);
      //for()
      smart_assign(guessSpecPlusTemp[0],0.);
      smart_assign(guessSpecPlusTemp[1],0.8);
      smart_assign(guessSpecPlusTemp[2],0.2);
      guessSpecPlusTemp[DataType::nspec] = Tignition_; //temperature
      std::cout << "guessSpecPlusTemp = "<<guessSpecPlusTemp<<std::endl;
      Recon_.initialize(guessSpecPlusTemp,PD.get_datum<DType>("hEstim"),PD.get_datum<BaseType>("newtTol"),PD.get_datum<char>("whichSolver"), PD.get_datum<int>("maxit"), PD.get_datum<BaseType>("nmtol"),PD.get_datum<BaseType>("Tol"));

      zM_.resize(DataType::nspec+1);  //Y+T
      Ctilde_.resize(1*DataType::nspec); //only chemistry considered
      btilde_.resize(1);
      for(SizeType s = 0; s < Ctilde_.size(); ++s)
	Ctilde_[s] = DataType::mass_balance_matrix()[s];
      for(SizeType s = 0; s < btilde_.size(); ++s)
	btilde_[s] = DataType::mass_sum()[s];
    
      std::cout << "Ctilde_ = "<< Ctilde_ << std::endl;
      std::cout << "btilde_ = "<< btilde_ << std::endl;
      
      low_.resize(DataType::nspec+1);
      up_.resize(DataType::nspec+1);
      for(int k = 0; k < DataType::nspec; ++k){	
	low_[k] = 0.;
	up_[k]  = 1.;
      }
      low_[DataType::nspec] = DataType::temperature_bounds()[0];
      up_[DataType::nspec] = DataType::temperature_bounds()[1];
      std::cout << "lower bounds: "<< low_ << std::endl << "upper bounds: "<< up_ << std::endl;
       
    }

    int dim() const {return totpoints_;}
    int domain_dim() const {return totpoints_;}
    std::string name() const {return "WITHOUT rho and momentum";}

    //============================
    //evaluate vector field
    //============================
    template<class VEC>
    VType& operator()(const VEC& vars){  
      //HARD COPY
      u_ = vars;
      
      //!BOUNDARY CONDITIONS

      //  !UP
      BaseType tbdy = BaseType();
      BaseType x = BaseType();
      for(int i = 0; i < nx_; ++i){
	//===RED
	reconstruct_full_state(zM_,i,ny_-1,u_);
      

	x = a_ + i*hx_;         //spatial coordinate x
	
	tbdy = bdytemp_(x);   //IGNITION temperature
	u_[OFF(i,ny_-1,TMP)] = tbdy;  //temperature
       
	for(int k = 0; k < nchem; ++k) 	//species 
	  u_[OFF(i,ny_-1,SPC+k)] = bdy2d_[AccessType::spec_access(k)];
      }
      //  !DOWN
      for(int i = 0; i < nx_; ++i){
	//===RED
	reconstruct_full_state(zM_,i,0,u_);
      }
	
      
      //  !RIGHT: NEUMANN
      for(int j = 1; j < ny_-1; ++j){
	//===RED
	reconstruct_full_state(zM_,nx_-1,j,u_);

	u_[OFF(nx_-1,j,TMP)] = u_[OFF(nx_-2,j,TMP)]; //temp
	for(int k = 0; k < nchem; ++k){ //species
	  u_[OFF(nx_-1,j,SPC+k)] = u_[OFF(nx_-2,j,SPC+k)];
	}
      }

      //  !LEFT
      for(int j = 1; j < ny_-1; ++j){
	//===RED
	reconstruct_full_state(zM_,0,j,u_);

	//BaseType tin = bdytemp_(0.);
	u_[OFF(0,j,TMP)] = Tinlet_; //tin;   //temperature

	for(int k = 0; k < nchem; ++k)   //species
	  u_[OFF(0,j,SPC+k)] = bdy2d_[AccessType::spec_access(k)];
      }
	
      //!INTERIOR
      for(int i = 1; i < nx_-1; ++i){
	for(int j = 1; j < ny_-1; ++j){
	  //===RED
	  reconstruct_full_state(zM_,i,j,u_);
	  //std::cout << "zM_ (INTERIOR) = "<< zM_ << std::endl;
	  
	  //=========TODO=====
	  //! ok assign values (5 pt stencil) for temperature from u_
	  //!  to Tab_ and don't take the updated temperature value 
	  take_temperature_from_argument(i,j,u_);
	  //==================
#ifndef NDEBUG  
	  is_temperature_ok(i,j,0,u_);  // 0 = take bounds of first chem. spec
#endif	   


	  v_[0] = veloprofile_[j];
	  v_[1] = 0.;

	  //convective parts v1 always > 0
	  T dx_T =  upwind_downwind_x(i,j,TMP,v_[0],u_),
	      dy_T =  upwind_downwind_y(i,j,TMP,v_[1],u_);
	  
	  for(int k = 0; k < nchem; ++k){
	    dx_spec_[k] = upwind_downwind_x(i,j,SPC+k,v_[0],u_);
	    dy_spec_[k] = upwind_downwind_y(i,j,SPC+k,v_[1],u_);
	  }

	  
	  rho_ = rho(i,j,u_); 
	  cp_ = cp(i,j,u_);
	  adonis_assert(rho_ > T() && cp_ > T());
	  
	  chemistry(dotomega_,i,j,u_); //chemical reactions


	  //!temperature
	  rhs_[OFF(i,j,TMP)] = -( v_[0]*dx_T + v_[1]*dy_T 
			       - 1./(rho_*cp_)*heat_conduction(i,j,u_)
			       - 1./(rho_*cp_)*temperature_diffu_term(i,j,u_)
			       + 1./(rho_*cp_)*heat_production(dotomega_,i,j,u_)
			       );

	  //!species --- +++FULL++++
	  species_diffu_term(DiffFull_,i,j,u_); //now DiffFull_ is filled
	  for(int k = 0; k < nchem; ++k){ //===RED (only)!!!
	    rhs_[OFF(i,j,SPC+k)] = -( v_[0]*dx_spec_[k] + v_[1]*dy_spec_[k] 
				   -1./rho_*DiffFull_[AccessType::spec_access(k)]
				   -1./rho_*dotomega_[k]*DataType::molar_masses()[AccessType::spec_access(k)]
				   );
	  }
	
	} //end interior j
      }   //end interior i
      return rhs_;
    }
 

    //dummy function
    void get_y_prev(typename VType::iterator y_prevIt, const T& h){}
    

  private:
    VType rhs_;
    BuilderType BCS_;

    BaseType a_,b_,c_,d_,hx_,hy_, diam_, halfdiam_;
    int nx_, ny_, npt_, totpoints_;   
   
    BaseType p0_, v1_, Tinlet_, Tignition_; //standard pressure
   
    VType veloprofile_,
      u_;
    T cp_, //for storage
      rho_;  //rho will be evaluated by an EOS
    

    T v_[2]; //velocity component

    T fivepointstencil_[5];
    T star5_[5];

    T dx_spec_[nchem];
    T dy_spec_[nchem];


    BoundaryTemperatureDistribution<BaseType> bdytemp_;
   
    ConductivityType Mac_;
    DiffusionType Mad_,Mad0_,Mad1_,Mad2_,Mad3_;
   
    BoundaryType bdy2d_; 
   
    //===RED
    TabulationType Tab_;
    ReconstructorType Recon_;
    VType zM_;
    RedChemType redspecies_;
    VDType Ctilde_,btilde_;
    VDType low_,up_;
    FullChemType X_, X0_, X1_,X2_,X3_,DiffFull_, Conc_, Cp_; //++++FULL++++
    RedChemType dotomega_;
    
    VType& reconstruct_full_state(VType& zm, int i, int j, const VType& u){
      for(int k = 0; k < DataType::rednspec; ++k)
	redspecies_[k] = u[OFF(i,j,SPC+k)]; //Y_k, reduced composition
      std::cout << "reduced species at ("<< i <<"," << j << "):  " << redspecies_ << std::endl;
      Recon_.assign(redspecies_,Tab_[GRID(i,j)]);
      Recon_.enforce_linearized_constraints_wrt_chemistry(Ctilde_,btilde_);
      Recon_.project_onto_bounds(low_,up_); //Y+T
      zm = Recon_.time_step(Recon_.get_z(),p0_); //Y+T 
       
      Tab_[GRID(i,j)] = zm;  //overwrite with new entry
      return zm;
    }

    //PRIVATE FUNCTIONS -- don't use them from outside

    int OFF(int i, int j, int spec) const{
      return (i + nx_*j + spec*npt_);
    }

    int GRID(int i, int j) const{  //access grid
      return (i + nx_*j);
    }
    //=============================================================
    //====RED: revise functions for model reduction 
    //===      assume reconstruct_full_state(zM_,i,j,u_) has been invoked before
    //===      some full chem. state evaluations are to be performed
    //! CHEMISTRY
    //in the reduced case  only the corresponding fixed quantities
    //are extracted, e.g. molar masses.
    T wbar(int i, int j, const VType& u) const{
      T wb = T();
      for(int k = 0; k < DataType::nspec; ++k){
	wb += Tab_[GRID(i,j)][k]/DataType::molar_masses()[k];
      }
      adonis_assert(wb > 0.);
      
      return 1./wb;
    }
   
    T rho(int i, int j, const VType& u){
      return (p0_*wbar(i,j,u))/(PhysicalConstants<BaseType>::Rgas*Tab_[GRID(i,j)][DataType::nspec]);
    }

    //+++FULL+++
    //MIGHT BE ZERO -- apply perturbations in fcts. mole_fraction
    T mole_fraction(int i, int j, int k, const T& Wbar, const VType& u) const{
      T mfr = (Tab_[GRID(i,j)][k]*Wbar)/DataType::molar_masses()[k]; 
      perturbation(mfr,1.e-16);
      return mfr;
    }

    FullChemType& mole_fraction(FullChemType& x, int i, int j, const VType u){
      for(int k = 0; k < DataType::nspec; ++k){
	x[k] = (Tab_[GRID(i,j)][k]*wbar(i,j,u))/DataType::molar_masses()[k];
	perturbation(x[k],1.e-16);
      }
      return x;
    }

    //if Wbar is reused somewhere
    FullChemType& mole_fraction(FullChemType& x, int i, int j, const T& Wbar, const VType u){
      for(int k = 0; k < DataType::nspec; ++k){
	x[k] = (Tab_[GRID(i,j)][k]*Wbar)/DataType::molar_masses()[k];
	perturbation(x[k],1.e-16);
      }
      return x;
    }

    //! FULL DIMENSION -- needed for chemistry
    FullChemType& concentration(FullChemType& conc, int i, int j, const VType u){
      T Rho =  rho(i,j,u); 
      for(int k = 0; k < DataType::nspec; ++k){
	conc[k] = (Rho*Tab_[GRID(i,j)][k])/DataType::molar_masses()[k];
       }
      return conc;
    }


    RedChemType& chemistry(RedChemType& domeg, int i, int j, const VType& u){
      concentration(Conc_,i,j,u); //C_ is filled now
      for(int k = 0; k < nchem; ++k){ //nchem is either full dim or reddim
	//!=========TODO: in either cases: C_ must be the FULL composition
	//! compare with 
	//! '../ode/examples/automaticallygeneratedsourceterms/redh2c6_attempt.hh'
	domeg[k] = BCS_.net_formation_rate(BCS_.species_index(k),Tab_[GRID(i,j)][DataType::nspec],Conc_); //in the reduced case C_ must be of full dimensino!!!
      }
      return domeg;
    }

    // FULL DIMENSION -- optional
    T heat_production(const RedChemType& dotomega, int i, int j, const VType & u){
      T h = T();
       // for(int k = 0; k < nchem; ++k)
       // 	h += BCS_.H_T(AccessType::spec_access(k),u[OFF(i,j,TMP)])*dotomega[k];
      int nurep = DataType::nspec-DataType::rednspec;
      for(int k = 0; k < DataType::rednspec; ++k)
	h += BCS_.H_T(AccessType::spec_access(k),Tab_[GRID(i,j)][DataType::nspec])*dotomega_[k];
      for(int k = 0; k < nurep; ++k)
	h += BCS_.H_T(Recon_.unrep_index()[k],Tab_[GRID(i,j)][DataType::nspec])*Recon_.dot_omega()[Recon_.unrep_index()[k]];

      return h;
    }

    //!FULL DIMENSION
    T cp(int i, int j, const VType& u){
      T cres = T();
      for(int k = 0; k < DataType::nspec; ++k){ //fill also Cp_
	Cp_[k] = BCS_.C_p(k,Tab_[GRID(i,j)][DataType::nspec]);
	cres += (Cp_[k]/DataType::molar_masses()[k])*Tab_[GRID(i,j)][k];
      }
      //isCpcalculated_ = true;
      return cres;
    }

    //!======== TODO: evaluate in Tab_[i][DataType::nspec] or u ??
    //! +++ FULL +++
    T lambda(int i, int j, const VType& u){
      (*this).mole_fraction(X_,i,j,wbar(i,j,u),u);
      perturbation(X_,1.e-16);
      return Mac_.compute_mixture_averaged_conductivity((*this).rho(i,j,u), p0_, Tab_[GRID(i,j)][DataType::nspec],X_); //X_ overwritten
    }

    T heat_conduction(int i, int j, const VType& u){
      return ( (0.5*(lambda(i+1,j,u)+lambda(i,j,u))*(Tab_[GRID(i+1,j)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(lambda(i,j,u)+lambda(i-1,j,u))*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i-1,j)][DataType::nspec]))/ntimes<2>(hx_) + (0.5*(lambda(i,j+1,u)+lambda(i,j,u))*(Tab_[GRID(i,j+1)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(lambda(i,j,u)+lambda(i,j-1,u))*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i,j-1)][DataType::nspec]))/ntimes<2>(hy_) );
    }

    
    //cp() must be invoked beforehand
     T temperature_diffu_term(int i, int j, const VType& u){
       T td = T(),
	 meanmolmass = wbar(i,j,u),
	 Rho = rho(i,j,u);

       //===RED: perturbation is applied in mole_fraction(...)
       mole_fraction(X_,i,j,meanmolmass,u);
       Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i,j)][DataType::nspec],X_); //diffusion coefficients are computed by now
       
       //isDiffusioncalculated_ = true;
       
       for(int k = 0; k < DataType::nspec; ++k){
	 //std::cout << "++++++ Self-diffusion coefficient = "<< Mad_.binary_diffusion_coefficients(k,k) << std::endl;
	 //!store counterclockwise, starting from the bottom with
	//!(i,j) being the last entry
	fivepointstencil_[0] = mole_fraction(i,j-1,k,wbar(i,j-1,u),u);
	fivepointstencil_[1] = mole_fraction(i+1,j,k,wbar(i+1,j,u),u);
	fivepointstencil_[2] = mole_fraction(i,j+1,k,wbar(i,j+1,u),u);
	fivepointstencil_[3] = mole_fraction(i-1,j,k,wbar(i-1,j,u),u);
	fivepointstencil_[4] = X_[k]; //mole_fraction(i,j,k,meanmolmass,u);

	//! assume Cp_ has been filled previously
	td += ( (Cp_[k]/meanmolmass)*Mad_.binary_diffusion_coefficients(k,k)*Rho*( (0.5*(fivepointstencil_[1]+fivepointstencil_[4])*(Tab_[GRID(i+1,j)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(fivepointstencil_[3]+fivepointstencil_[4])*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i-1,j)][DataType::nspec]))/ntimes<2>(hx_) + (0.5*(fivepointstencil_[2]+fivepointstencil_[4])*(Tab_[GRID(i,j+1)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(fivepointstencil_[4]+fivepointstencil_[0])*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i,j-1)][DataType::nspec]))/ntimes<2>(hy_)) );
      }
      return td;
     }
  
    //+++FULL++++
    FullChemType& species_diffu_term(FullChemType& diffu, int i, int j, const VType& u){
      fivepointstencil_[0] = wbar(i,j-1,u);
      fivepointstencil_[1] = wbar(i+1,j,u);
      fivepointstencil_[2] = wbar(i,j+1,u);
      fivepointstencil_[3] = wbar(i-1,j,u);
      fivepointstencil_[4] = wbar(i,j,u);

      //===RED: already perturbed in mole_fraction(..)
      mole_fraction(X0_,i,j-1,fivepointstencil_[0],u);  //(i,j-1)
      mole_fraction(X1_,i+1,j,fivepointstencil_[1],u);  //(i+1,j)
      mole_fraction(X2_,i,j+1,fivepointstencil_[2],u);  //(i,j+1)
      mole_fraction(X3_,i-1,j,fivepointstencil_[3],u);  //(i-1,j)
      mole_fraction(X_,i,j,fivepointstencil_[4],u);       //(i,j)

      Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i,j-1)][DataType::nspec],X0_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i+1,j)][DataType::nspec],X1_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i,j+1)][DataType::nspec],X2_);
      
      Mad3_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i-1,j)][DataType::nspec],X3_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i,j)][DataType::nspec],X_);
      
      //CAUTION: use self-diffusion coefficient for single gas species
      for(int k = 0; k < DataType::nspec; ++k){
	diffu[k] = (0.5*(rho(i+1,j,u)*DataType::molar_masses()[k]/fivepointstencil_[1]*Mad1_.binary_diffusion_coefficients(k,k) + rho(i,j,u)*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_.binary_diffusion_coefficients(k,k))*(X1_[k]-X_[k]) - 0.5*(rho(i,j,u)*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_.binary_diffusion_coefficients(k,k) + rho(i-1,j,u)*DataType::molar_masses()[k]/fivepointstencil_[3]*Mad3_.binary_diffusion_coefficients(k,k))*(X_[k]-X3_[k]))/ntimes<2>(hx_) + 
	  (0.5*(rho(i,j+1,u)*DataType::molar_masses()[k]/fivepointstencil_[2]*Mad2_.binary_diffusion_coefficients(k,k) + rho(i,j,u)*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_.binary_diffusion_coefficients(k,k))*(X2_[k]-X_[k]) - 0.5*(rho(i,j,u)*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_.binary_diffusion_coefficients(k,k) + rho(i,j-1,u)*DataType::molar_masses()[k]/fivepointstencil_[0]*Mad0_.binary_diffusion_coefficients(k,k))*(X_[k]-X0_[k]))/ntimes<2>(hy_);

      }
      
      return diffu;
    }

    //with wbar
     FullChemType& species_diffu_term_wbar(FullChemType& diffu, int i, int j, const VType& u){
       fivepointstencil_[0] = wbar(i,j-1,u);
       fivepointstencil_[1] = wbar(i+1,j,u);
       fivepointstencil_[2] = wbar(i,j+1,u);
       fivepointstencil_[3] = wbar(i-1,j,u);
       fivepointstencil_[4] = wbar(i,j,u);

       mole_fraction(X0_,i,j-1,fivepointstencil_[0],u);  //(i,j-1)
       mole_fraction(X1_,i+1,j,fivepointstencil_[1],u);  //(i+1,j)
       mole_fraction(X2_,i,j+1,fivepointstencil_[2],u);  //(i,j+1)
       mole_fraction(X3_,i-1,j,fivepointstencil_[3],u);  //(i-1,j)
       mole_fraction(X_,i,j,fivepointstencil_[4],u);       //(i,j)

       Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i,j-1)][DataType::nspec],X0_);
       Mad1_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i+1,j)][DataType::nspec],X1_);
       Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i,j+1)][DataType::nspec],X2_);
      
       Mad3_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i-1,j)][DataType::nspec],X3_);
       Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[GRID(i,j)][DataType::nspec],X_);
       
       star5_[0] = rho(i,j-1,u);
       star5_[1] = rho(i+1,j,u);
       star5_[2] = rho(i,j+1,u);
       star5_[3] = rho(i-1,j,u);
       star5_[4] = rho(i,j,u);

       for(int k = 0; k < DataType::nspec; ++k){
	//!\f$\nabla \cdot (\rho D^{\mathrm{mix}}_k \nabla Y_k)\f$
	diffu[k] = (0.5*(star5_[1]*Mad1_.binary_diffusion_coefficients(k,k) + star5_[4]*Mad_.binary_diffusion_coefficients(k,k))*(Tab_[GRID(i+1,j)][k]-Tab_[GRID(i,j)][k]) - 0.5*(star5_[4]*Mad_.binary_diffusion_coefficients(k,k) + star5_[3]*Mad3_.binary_diffusion_coefficients(k,k))*(Tab_[GRID(i,j)][k]-Tab_[GRID(i-1,j)][k]))/ntimes<2>(hx_) 
	  + (0.5*(star5_[2]*Mad2_.binary_diffusion_coefficients(k,k)+ star5_[4]*Mad_.binary_diffusion_coefficients(k,k))*(Tab_[GRID(i,j+1)][k]-Tab_[GRID(i,j)][k]) - 0.5*(star5_[4]*Mad_.binary_diffusion_coefficients(k,k)+ star5_[0]*Mad0_.binary_diffusion_coefficients(k,k))*(Tab_[GRID(i,j)][k]-Tab_[GRID(i,j-1)][k]))/ntimes<2>(hy_) 
	  //!differential part with mean molecular weight \f$ \overline{W}\f$,
	  //!i.e. \f$ \nabla \cdot (\rho Y_k/\overline{W}D^{\mathrm{mix}}_k \nabla \overline{W})\f$
	  +  (0.5*(star5_[1]*Tab_[GRID(i+1,j)][k]/fivepointstencil_[1]*Mad1_.binary_diffusion_coefficients(k,k) + star5_[4]*Tab_[GRID(i,j)][k]/fivepointstencil_[4]*Mad_.binary_diffusion_coefficients(k,k))*(fivepointstencil_[1]-fivepointstencil_[4]) - 0.5*(star5_[4]*Tab_[GRID(i,j)][k]/fivepointstencil_[4]*Mad_.binary_diffusion_coefficients(k,k) + star5_[3]*Tab_[GRID(i-1,j)][k]/fivepointstencil_[3]*Mad3_.binary_diffusion_coefficients(k,k))*(fivepointstencil_[4]-fivepointstencil_[3]))/ntimes<2>(hx_)
	  +  (0.5*(star5_[2]*Tab_[GRID(i,j+1)][k]/fivepointstencil_[2]*Mad2_.binary_diffusion_coefficients(k,k)+ star5_[4]*Tab_[GRID(i,j)][k]/fivepointstencil_[4]*Mad_.binary_diffusion_coefficients(k,k))*(fivepointstencil_[2]-fivepointstencil_[4]) - 0.5*(star5_[4]*Tab_[GRID(i,j)][k]/fivepointstencil_[4]*Mad_.binary_diffusion_coefficients(k,k) + star5_[0]*Tab_[GRID(i,j-1)][k]/fivepointstencil_[0]*Mad0_.binary_diffusion_coefficients(k,k))*(fivepointstencil_[4]-fivepointstencil_[0]))/ntimes<2>(hy_);
       }
       return diffu;
     }

    //!create upwind or downwind difference for physical quantity t,
    //! depending on the sign of v_l
    T upwind_downwind_x(int i, int j, int t, const T& v_l, const VType& u){
      
      if(v_l <= 0) //first order upwind
	return (u[OFF(i+1,j,t)] - u[OFF(i,j,t)])/hx_;
      else // > 0
	return (u[OFF(i,j,t)] - u[OFF(i-1,j,t)])/hx_;
    }

     T upwind_downwind_y(int i, int j, int t, const T& v_l, const VType& u){
      
      if(v_l <= 0) //first order upwind
	return (u[OFF(i,j+1,t)] - u[OFF(i,j,t)])/hy_;
      else //> 0
	return (u[OFF(i,j,t)] - u[OFF(i,j-1,t)])/hy_;
    }
    
    void is_temperature_ok(int i, int j, int k, const VType& u){
      std::string oper1, oper2;
      for(int l = -1; l < 2; ++l){
	for(int t = -1; t < 2; ++t){
	  int IX = i+l,
	    JX = j+t;
	  
	 
	  if(l>=0)
	    oper1 = "+";
	  if(t>=0)
	    oper2 = "+";
	  
	  if(u[OFF(IX,JX,TMP)] < DataType::temperature_bounds()[3*k] || 
	     u[OFF(IX,JX,TMP)] > DataType::temperature_bounds()[3*k+1]){
	   
	    ADONIS_ERROR(BoundsError,"Temperature bounds violated: T_{i"<<oper1<<l << ",j"<<oper2<<t << "} = "<< u[OFF(IX,JX,TMP)]<<", i.e. T("<<IX<<","<<JX<<").");
	  }
	}
      }
    }


    void take_temperature_from_argument(int i, int j, const VType& u){
      //update 5 point stencil
      Tab_[GRID(i,j)][DataType::nspec] = u[OFF(i,j,TMP)];
      Tab_[GRID(i,j-1)][DataType::nspec] = u[OFF(i,j-1,TMP)];
      Tab_[GRID(i+1,j)][DataType::nspec] = u[OFF(i+1,j,TMP)];
      Tab_[GRID(i,j+1)][DataType::nspec] = u[OFF(i,j+1,TMP)];
      Tab_[GRID(i-1,j)][DataType::nspec] = u[OFF(i-1,j,TMP)];
    }

  };


} //end namespace 

#endif
