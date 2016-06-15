#ifndef REDUCED_NAV_STOKES_VIA_MOL_FOR_PRF_H2_MECH_HH
#define REDUCED_NAV_STOKES_VIA_MOL_FOR_PRF_H2_MECH_HH
/**
 * This functor can be used by an ODE integrator in th usual way
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

//transport quantities
#include "../moleculartransport/viscosity.hh"
#include "../moleculartransport/conductivity.hh"
#include "../moleculartransport/diffusion.hh"
#include "../moleculartransport/fancytransporttmps.hh"

#include "../moleculartransport/thermaldiffusion.hh"

#include "../templatemetaprograms/functionalities4containersandscalars.hh"

#include "gensettings.hh"
#include "usefulextensions.hh"

#ifdef USE_OPENMP
#include "omp.h"
#endif

#include "../common/fancymessages.hh"
#include "reconstructspecies.hh" //===RED

namespace Adonis{

  template<class T>
  class ReducedConaireH2ReactiveNS2D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType; 
    typedef typename TypeAdapter<T>::Type DType; //===RED
    typedef ExprTmpl::MyVec<DType> VDType; 
    typedef typename TypeAdapter<T>::BaseType BaseType; //for norms, tolerances 
    typedef std::size_t SizeType;
    
     //needed for thermo chemistry: primary reference fuel mechanism
    typedef ThermoData4Mechanism<BaseType,10> DataType;  

    
    enum{nchem = DataType::rednspec,   //===RED 
	 //rho,v1,v2,T + specs
	 nprim = 4+DataType::rednspec,  //===RED
	 Encoding = 10
    };
    //===RED
    typedef ExprTmpl::MyVec<ExprTmpl::MyVec<T> > TabulationType;
    typedef StaticArray<T,DataType::nspec> FullChemType;
    typedef StaticArray<T,DataType::rednspec> RedChemType;

    typedef StaticArray<T,2> DimensionType;
    //! element 0 is x-direction, 1 y-direction and so forth
    
    typedef BoundaryTemperatureDistribution<BaseType,FDMSettings::WallTemperatureType> BoundaryWallTemperatureType;

    typedef InletVelocity<BaseType,BaseType,BaseType,FDMSettings::InletVelocityType> InVeloType;


    enum{O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR};

    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;

    //Other Index
    typedef IndexHandler<DataType,true> AccessType; //===RED
    typedef typename AccessType::IndexerType IndexerType;

    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,IndexerType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    typedef EquationOfState<'i'> EosType;

    //TRANSPORT COEFFICIENTS -- +++++FULL+++++ version here anyway
    typedef MixtureAveragedViscosity<DataType,false> ViscosityType;
    typedef MixtureAveragedConductivity<DataType,false> ConductivityType;
    typedef MixtureAveragedDiffusion<DataType,false> DiffusionType;

    typedef ThermalDiffusionRatio<DataType,false,FDMSettings::ThermalDiffusion4LightWeightSpeciesOnly> RatioType;
    
    
    //====RED: with non-Cppad type since we solve system via LAPACK and integr.
     typedef ReconstructSpecies<DType,DataType::nspec,DataType::rednspec,Encoding> ReconstructorType; 

    ReducedConaireH2ReactiveNS2D(SizeType dim = 0):rhs_(dim){
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

      //! include information of possible explicitly given rate coefficients
      BCS_.initialize(&forwardRates,&reverseRates,&Ixer,DataType::explicit_reverse_reaction_coefficients());

      //READ IN DATA
      ParameterData PD;
      PD.read_from_file("data/conaireh2.dat");
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
      t0_ = PD.get_datum<BaseType>("t0");
      tend_ = PD.get_datum<BaseType>("tend");
      halfdiam_ = 0.5*diam_;
      gravi1_ = PD.get_datum<BaseType>("g1");
      gravi2_ = PD.get_datum<BaseType>("g2");
      near_zero_ = Abs(PD.get_datum<BaseType>("near_zero"));
      p0_ = PD.get_datum<BaseType>("pconst");
      v1_in_ = PD.get_datum<BaseType>("v1_in");
      v2_in_ = PD.get_datum<BaseType>("v2_in");
      T_wall_ = PD.get_datum<BaseType>("T_wall");
      T_in_ = PD.get_datum<BaseType>("T_min");
      x_ignition_ = (b_-a_)/PD.get_datum<BaseType>("channel_length_frac");
      
      Tignition_ = PD.get_datum<BaseType>("T_ignition");
      yhottest_ = PD.get_datum<BaseType>("hottest_y");
      bdytemp_.initialize(T_wall_,Tignition_,x_ignition_, PD.get_datum<BaseType>("plateau_width")*hx_); //4th argument is here only meaningful when 'm' profile is used

      useGravi_ = PD.get_datum<bool>("gravForce");
      useMechCompression_ = PD.get_datum<bool>("mechCompress");

      

#ifndef NO_CORRECTION_DIFFUSION_VELOCITY
      std::cout << "COMPUTE CORRECTION VELOCITY:"<< std::endl;
      std::cout << "----------------------------"<<std::endl<<std::endl;
#ifdef CORRECTION_DIFFUSION_VELOCITY_WITH_THERMAL_PART  
      std::cout << "  V_corr = sum(Y_kV_k + Y_k W_k)" << std::endl; 
#else
      std::cout << "  V_corr = sum(Y_kV_k)" << std::endl;
#endif
#else
      std::cout << "NO CORRECTION VELOCITY COMPUTED" << std::endl;
#endif

      std::cout << std::endl;
      
#ifdef LOW_MACH_NUMBER_FORMULATION
      std::cout << "Compute low-Mach-number formulation for mixture density RHO" << std::endl << std::endl;
#else
      std::cout << "Compute d_t rho + D · (rho v) as usual"<< std::endl << std::endl;
#endif
      std::cout <<  "gravitational force: " << yes_no(useGravi_) <<"      mechanical compression: "<< yes_no(useMechCompression_) << std::endl;

#ifdef IGNITION_AT_HOT_WALLS
      std::cout << std::endl << "Ignition at hot walls; temperature profile: "<< FDMSettings::WallTemperatureType << std::endl<< std::endl;
#else
      std::cout << std::endl << "Ignition left side" << std::endl<< std::endl;
#endif

      //transport coefficients are ready for computation
      Mav_.initialize(DataType::transport());
      Mac_.initialize(DataType::transport());
      Mad_.initialize(DataType::transport());  //(i,j)
      Mad0_.initialize(DataType::transport()); //(i,j-1)
      Mad1_.initialize(DataType::transport()); //(i+1,j)
      Mad2_.initialize(DataType::transport()); //(i,j+1)
      Mad3_.initialize(DataType::transport()); //(i-1,j)


      //thermal diffusion ratios
      Theta_.initialize(DataType::transport());
      Theta0_.initialize(DataType::transport());
      Theta1_.initialize(DataType::transport());
      Theta2_.initialize(DataType::transport());
      Theta3_.initialize(DataType::transport());

      Vk_ = T();  


      v1bdy_.resize(ny_);   //for the velocity profile (along y direction)
      u_.resize(npt_);
      cp_ = T();
      tempconc_ = T();
      tempDiff_ = T();
      specDiff_ = T();
      tempHeat_ = T();
      specHeat_ = T();
      tempCond_ = T();

      //O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR  //full vector
      Y_in_[O] = PD.get_datum<BaseType>("O");    
      Y_in_[O2] = PD.get_datum<BaseType>("O2");  
      Y_in_[H] = PD.get_datum<BaseType>("H");   
      Y_in_[OH] = PD.get_datum<BaseType>("OH"); 
      Y_in_[H2] = PD.get_datum<BaseType>("H2");    
      Y_in_[HO2] = PD.get_datum<BaseType>("HO2");   
      Y_in_[H2O2] = PD.get_datum<BaseType>("H2O2");    
      Y_in_[H2O] = PD.get_datum<BaseType>("H2O");      
      Y_in_[N2] = PD.get_datum<BaseType>("N2");
      Y_in_[AR] = PD.get_datum<BaseType>("AR");
      

     
      T sm(0.);
      for(int k = 0; k < DataType::nspec; ++k)
	sm += Y_in_[k];
      if(sm != 1.)
	ADONIS_ERROR(DerivedError,"Mixture data do not sum to 1!");
    

      //fixed parabolic or const inlet profile in x-direction
      // also initial value profile
      veloprofile_.resize(ny_); //along y-axis
      BaseType r = -halfdiam_;  //R = d/2, r = -R
      for(int j = 0; j < ny_; ++j){
        veloprofile_[j] = InVeloType::velocity(r,halfdiam_,v1_in_); //zero(parabolic_inlet_velocity(r,halfdiam_,v1_in_)); 
	r += hy_;
      }

      //fixed smooth temperature profile 
      BaseType xdir(0.);
      walltempprofile_.resize(nx_);  //along x-axis
      for(int i = 0; i < nx_; ++i){
	walltempprofile_[i] = bdytemp_(xdir);
	xdir += hx_;
      }


      //===RED=================================
         Tab_.resize(nx_*ny_); //store zM at each grid point 
       for(int i = 0; i < nx_; ++i){  //intial values for full settings
	 // DType x = a_ + i*hx_; 
	 for(int j = 0; j < ny_; ++j){
	   (Tab_[GRID(i,j)]).resize(DataType::nspec+1);
	 
	   for(int k = 0; k < DataType::nspec; ++k){
	     Tab_[GRID(i,j)][k] = Y_in_[k];
	   }
	   //std::cout << std::endl;
	   Tab_[GRID(i,j)][DataType::nspec] = T_in_;
	   
	   if((j == 0) || (j == ny_-1)){
#ifdef IGNITION_AT_HOT_WALLS
	   Tab_[GRID(i,j)][DataType::nspec] = walltempprofile_[i]; // bdytemp_(x);
#else
	   Tab_[GRID(i,j)][DataType::nspec] = T_wall_;
#endif
	   }
	 }
       }

       //=======================================
        //! INITIALIZE species Reconstructor
      //!some guess
      VDType guessSpecPlusTemp(DataType::nspec+1);
      for(int k = 0; k < DataType::nspec; ++k){
	smart_assign(guessSpecPlusTemp[k],Y_in_[k]);
      }
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
       
      
      FancyMessages().nice_output("\n REDUCTION in use", 35); //35=purple
                                                              //36=cyan

      //! corresponing to the boundary values/initial values (full)
      wbarIn_ = BaseType();
      for(int k = 0; k < DataType::nspec; ++k)
	wbarIn_ += Y_in_[k]/DataType::molar_masses()[AccessType::spec_access(k)];
      
      wbarIn_ = 1./wbarIn_;
      rhoIn_ = (p0_*wbarIn_)/(PhysicalConstants<BaseType>::Rgas*T_in_);

      std::cout << "CFD H2 object created"<< std::endl;
    } //end constructor

    //! \f$\rho, v_1, v_2, T\f$ + nchem species, discretized on \f$\Omega_h\f$
      //! with npt_ nodes
    int dim() const {return totpoints_;}
    int domain_dim() const {return totpoints_;}

    std::string name() const {return "with momentum";}

    //============================
    //evaluate vector field
    //============================
    template<class VEC>
    VType& operator()(const VEC& vars){ 
      //HARD COPY
      u_ = vars;
     
      //!BOUNDARY -- imagine szenario is rotated by 90°
      //!  LEFT: inlet
      DType y(0.);
      for(int j = 1; j < ny_-1; ++j){
	y = j*hy_;
	
	u_[OFF(0,j,0)] = rhoIn_;
	
	u_[OFF(0,j,1)] = veloprofile_[j];
	u_[OFF(0,j,2)] = 0.;

#ifdef IGNITION_AT_HOT_WALLS
	//!====== TODO: (un)comment
	u_[OFF(0,j,3)] = T_in_;	
#else
	//!========== TODO: play around with different ignition regions
	//!small ignition area
	if(is_contained(yhottest_,y-hy_,y+hy_)){
	  u_[OFF(0,j,3)] = Tignition_;
	  
	   // rho now alters with temperature an composition
	  u_[OFF(0,j,0)] = (p0_*wbarIn_)/(PhysicalConstants<BaseType>::Rgas*Tignition_);
	}
	else 
	  u_[OFF(0,j,3)] = T_in_;
#endif


	//!=============================================================
	  
	for(int k = 0; k < nchem; ++k)
	  u_[OFF(0,j,4+k)] = Y_in_[AccessType::spec_access(k)];
      }

   
      for(int i = 0; i < nx_; ++i){	
	//! DOWN · DOWN · DOWN · DOWN · DOWN · DOWN · DOWN: wall

	u_[OFF(i,0,1)] = 0.; //no-slip condition
	u_[OFF(i,0,2)] = 0.;

	//========TODO: maybe other wall temperature
#ifdef IGNITION_AT_HOT_WALLS
	u_[OFF(i,0,3)] = walltempprofile_[i];
	//change with temperature
	u_[OFF(i,0,0)] = (pressure(i,0,u_)*wbar(i,0,u_))/(PhysicalConstants<DType>::Rgas*u_[OFF(i,0,3)]);
	
#else
	u_[OFF(i,0,3)] = T_wall_;  
	u_[OFF(i,0,0)] = (p0_*wbar(i,0,u_))/(PhysicalConstants<BaseType>::Rgas*T_wall_);
#endif
	//!=============================================================
	

	for(int k = 0; k < nchem; ++k)
	  u_[OFF(i,0,4+k)] = u_[OFF(i,1,4+k)]; //Neumann
      

       //!UP · UP · UP  · UP · UP · UP · UP · UP · UP: wall -- not necessarily 
       //!a symmetric scenario
#ifdef IGNITION_AT_HOT_WALLS
	//========TODO: maybe other wall temperature
	u_[OFF(i,ny_-1,3)] = walltempprofile_[i];
	//change with temperature
	u_[OFF(i,ny_-1,0)] = (pressure(i,ny_-1,u_)*wbar(i,ny_-1,u_))/(PhysicalConstants<DType>::Rgas*u_[OFF(i,ny_-1,3)]);
#else
	u_[OFF(i,ny_-1,0)] = (p0_*wbar(i,ny_-1,u_))/(PhysicalConstants<BaseType>::Rgas*T_wall_);
	u_[OFF(i,ny_-1,3)] = T_wall_;
#endif
//!=============================================================
	u_[OFF(i,ny_-1,1)] = 0.;   //no-slip condition
	u_[OFF(i,ny_-1,2)] = 0.;
	
     
	for(int k = 0; k < nchem; ++k)
	  u_[OFF(i,ny_-1,4+k)] = u_[OFF(i,ny_-2,4+k)]; //Neumann
      }


      //! RIGHT: outlet
      for(int j = 1; j < ny_-1; ++j){ 
	//u_[OFF(nx_-1,j,0)] = (pressure(nx_-1,j,u_)*wbar(nx_-1,j,u_))/(PhysicalConstants<DType>::Rgas*u_[OFF(nx_-1,j,3)]);
	//Neumann
	u_[OFF(nx_-1,j,1)] = u_[OFF(nx_-2,j,1)];
	u_[OFF(nx_-1,j,2)] = u_[OFF(nx_-2,j,2)];
	u_[OFF(nx_-1,j,3)] = u_[OFF(nx_-2,j,3)];
	for(int k = 0; k < nchem; ++k)
	  u_[OFF(nx_-1,j,4+k)] = u_[OFF(nx_-2,j,4+k)];
      }

     

      //INTERIOR -- use second order approximations for first, sec. and 
      //            mixed derivatives
#ifdef USE_OPENMP
      //adonis_assert((NUMTHREADS > 0));
#pragma omp parallel for num_threads(NUMTHREADS)
#endif 
      for(int i = 1; i < nx_-1; ++i){ 
	for(int j = 1; j < ny_-1; ++j){
	  //===RED
	  //reconstruct_full_state(zM_,i,j,u_);
	  //FullDotOmegaViaReconstruction_ =  Recon_.dot_omega();

	  //!===RED
	  for(int k =0; k < nchem; ++k)
	    de_negativize(u_[OFF(i,j,4+k)]);

	  //p0_ = (u_[OFF(i,j,0)]*PhysicalConstants<BaseType>::Rgas*u_[OFF(i,j,3)])/wbar(i,j,u_);

	  //! preserve mass Y_nspec = 1 - \sum_k^{nspec-1} Y_k
	  //preserve_mass(N2,i,j,u_);   //N2 is excess species (bath gas)

	  //check direction of velocity
	  T dx_v1 = upwind_downwind_x(i,j,1,u_[OFF(i,j,1)],u_),
	    dy_v1 = upwind_downwind_y(i,j,1,u_[OFF(i,j,2)],u_),
	    dx_v2 = upwind_downwind_x(i,j,2,u_[OFF(i,j,1)],u_),
	    dy_v2 = upwind_downwind_y(i,j,2,u_[OFF(i,j,2)],u_),
	    dx_T =  upwind_downwind_x(i,j,3,u_[OFF(i,j,1)],u_),
	    dy_T =  upwind_downwind_y(i,j,3,u_[OFF(i,j,2)],u_);
	  
#ifdef LOW_MACH_NUMBER_FORMULATION
	  T dx_rho = upwind_downwind_x(i,j,0,u_[OFF(i,j,1)],u_),
	    dy_rho = upwind_downwind_y(i,j,0,u_[OFF(i,j,2)],u_);
#endif
	  
	    for(int k = 0; k < nchem; ++k){
	      dx_spec_[k] = upwind_downwind_x(i,j,4+k,u_[OFF(i,j,1)],u_);
	      dy_spec_[k] = upwind_downwind_y(i,j,4+k,u_[OFF(i,j,2)],u_);
	    }
	    
#ifndef LOW_MACH_NUMBER_FORMULATION //if not defined then use usual computation	
	    //!\f$\rho\f$
	    rhs_[OFF(i,j,0)] = -( (u_[OFF(i+1,j,0)]*u_[OFF(i+1,j,1)] - u_[OFF(i,j,0)]*u_[OFF(i,j,1)])/hx_ + (u_[OFF(i,j+1,0)]*u_[OFF(i,j+1,2)] - u_[OFF(i,j,0)]*u_[OFF(i,j,2)])/hy_ );
#endif

	    //!\f$v_1\f$
	    rhs_[OFF(i,j,1)] = -(u_[OFF(i,j,1)]*dx_v1 + u_[OFF(i,j,2)]*dy_v1 
	    			 - 1./u_[OFF(i,j,0)]*(4./3.*( (0.5*(mu(i+1,j,u_)+mu(i,j,u_))*(u_[OFF(i+1,j,1)]-u_[OFF(i,j,1)]) - 0.5*(mu(i-1,j,u_)+mu(i,j,u_))*(u_[OFF(i,j,1)]-u_[OFF(i-1,j,1)]))/ntimes<2>(hx_)) - 2./3.*(mu(i+1,j,u_)*(u_[OFF(i+1,j+1,2)] - u_[OFF(i+1,j-1,2)]) - mu(i-1,j,u_)*(u_[OFF(i-1,j+1,2)]-u_[OFF(i-1,j-1,2)]))/(4.*hx_*hy_)) - 1./u_[OFF(i,j,0)]*( (0.5*(mu(i,j+1,u_)+mu(i,j,u_))*(u_[OFF(i,j+1,1)] - u_[OFF(i,j,1)]) - 0.5*(mu(i,j-1,u_)+mu(i,j,u_))*(u_[OFF(i,j,1)] - u_[OFF(i,j-1,1)]))/ntimes<2>(hy_) + (mu(i,j+1,u_)*(u_[OFF(i+1,j+1,2)] - u_[OFF(i-1,j+1,2)]) - mu(i,j-1,u_)*(u_[OFF(i+1,j-1,2)] - u_[OFF(i-1,j-1,2)]))/(4.*hx_*hy_)) 
				 - gravitational_force_x(i,j,u_)
				 + mechanical_compression_x(i,j,u_)
				 ); 
	   
	    //! \f$v_2\f$  
	    rhs_[OFF(i,j,2)] = -(u_[OFF(i,j,1)]*dx_v2 + u_[OFF(i,j,2)]*dy_v2 
				 - 1./u_[OFF(i,j,0)]*( (mu(i+1,j,u_)*(u_[OFF(i+1,j+1,1)] - u_[OFF(i+1,j-1,1)]) - mu(i-1,j,u_)*(u_[OFF(i-1,j+1,1)]-u_[OFF(i-1,j-1,1)]))/(4.*hx_*hy_) + (0.5*(mu(i+1,j,u_)+mu(i,j,u_))*(u_[OFF(i+1,j,2)]-u_[OFF(i,j,2)]) - 0.5*(mu(i-1,j,u_)+mu(i,j,u_))*(u_[OFF(i,j,2)]-u_[OFF(i-1,j,2)]))/ntimes<2>(hx_) ) -  1./u_[OFF(i,j,0)]*(4./3.*((0.5*(mu(i,j+1,u_)+mu(i,j,u_))*(u_[OFF(i,j+1,2)] - u_[OFF(i,j,2)]) - 0.5*(mu(i,j-1,u_)+mu(i,j,u_))*(u_[OFF(i,j,2)] - u_[OFF(i,j-1,2)]))/ntimes<2>(hy_)) - 2./3.*((mu(i,j+1,u_)*(u_[OFF(i+1,j+1,1)] - u_[OFF(i-1,j+1,1)]) - mu(i,j-1,u_)*(u_[OFF(i+1,j-1,1)] - u_[OFF(i-1,j-1,1)]))/(4.*hx_*hy_))) 
				 - gravitational_force_y(i,j,u_)
				 + mechanical_compression_y(i,j,u_)
				 ); 
				 
	    //!TEMPERATURE \f$ T \f$
	    
#ifndef NDEBUG  
	    is_temperature_ok(i,j,0,u_);  // 0 = take bounds of first chem. spec
#endif	   
	    //project_temperature_on_its_bounds(i,j,0,u); //brings nothing

	    cp_ = cp(i,j,u_); //Cp_ has also been filled
	   
	    //!CHEMISTRY -- calculate \f$\dot{\omega}\f$
	    chemistry(dotomega_,i,j,u_); //chemical reactions

	    //! note: for low-speed flows, mechanical compression  \f$\frac{Dp}{Dt}\f$ as well viscous dissipation can be neglected [KEE, p. 115 (pdf: 142), bottom], which is realized here
	    tempCond_ = heat_conduction(i,j,u_);
	    tempDiff_ = temperature_diffu_term(i,j,u_);
	    tempHeat_ = heat_production(dotomega_,i,j,u_);
	    rhs_[OFF(i,j,3)] = -(u_[OFF(i,j,1)]*dx_T + u_[OFF(i,j,2)]*dy_T 
				 - 1./(u_[OFF(i,j,0)]*cp_)*tempCond_
				 -               //PLUS/MINUS ??
				 1./(cp_*u_[OFF(i,j,0)])*tempDiff_
				 +                 //PLUS ??
				 1./(cp_*u_[OFF(i,j,0)])*tempHeat_ ); 
	    
	    
	    //! SPECIES  --- +++FULL++++
	    species_diffu_term(diffspec_,i,j,u_); //now diffspec_ is filled
	    //species_diffu_term_wbar(diffspec_,i,j,u_); //now diffspec_ filled
	    for(int k = 0; k < nchem; ++k){//===RED (only)!!!
	      rhs_[OFF(i,j,4+k)] = -(u_[OFF(i,j,1)]*dx_spec_[k] + u_[OFF(i,j,2)]*dy_spec_[k]
                                      -              //PLUS/MINUS???
                                      1./u_[OFF(i,j,0)]*diffspec_[AccessType::spec_access(k)]
                                      -
				     1./u_[OFF(i,j,0)]*dotomega_[k]*DataType::molar_masses()[AccessType::spec_access(k)] );
	    }
			  
#ifdef LOW_MACH_NUMBER_FORMULATION
	    T sm = T();  //==========TODO: DataTyp::nspec and dot from reconstruction????
	    for(int k = 0; k < nchem; ++k){
	      sm +=  (-1./u_[OFF(i,j,0)]*diffspec_[k]
                      -
		      1./u_[OFF(i,j,0)]*dotomega_[k]*DataType::molar_masses()[AccessType::spec_access(k)] );
	    }
	    
	    rhs_[OFF(i,j,0)] = -(u_[OFF(i,j,1)]*dx_rho + u_[OFF(i,j,2)]*dy_rho   -u_[OFF(i,j,0)]*(sm + (- 1./(Tab_[GRID(i,j)][DataType::nspec]*cp_*u_[OFF(i,j,0)])*(tempCond_
																				    -  tempDiff_ +  tempHeat_) ) ) );
	
#endif

	    //====RED: reconstruction of Y, T
	    reconstruct_full_state(zM_,i,j,u_);   //here or at the beginning of for-for-loop???
	    //FullDotOmegaViaReconstruction_ =  Recon_.dot_omega();
	    //======
			  
	} //END INTERIOR y
      }  //END INTERIOR x
	  
      return rhs_;
    }

    //dummy function
    void get_y_prev(typename VType::iterator y_prevIt, const T& h){}
    
  private:
    VType rhs_,
      veloprofile_,
      walltempprofile_;
    BuilderType BCS_;

    BaseType a_,b_,c_,d_,hx_,hy_, t0_, tend_, diam_, halfdiam_,gravi1_,gravi2_, near_zero_;
    T  wbarIn_, rhoIn_;
    int nx_, ny_, npt_, totpoints_;   
   
    BaseType p0_, v1_in_,v2_in_,T_wall_,T_in_, x_ignition_; //standard pressure
   
    //!====RED
    TabulationType Tab_;
    ReconstructorType Recon_;
    VType zM_;
    RedChemType redspecies_;
    VDType Ctilde_,btilde_;
    VDType low_,up_;


    FullChemType X_, X0_, X1_,X2_,X3_,diffspec_, C_, Cp_,
      Y_in_,
      FullDotOmegaViaReconstruction_;  //===RED
    RedChemType  dotomega_;

    VType v1bdy_,
      u_;
    T cp_, tempconc_, tempDiff_, specDiff_, tempHeat_, specHeat_, tempCond_; //for storage
    
    T fivepointstencil_[5];
    T wbarfive_[5];

    T dx_spec_[nchem];
    T dy_spec_[nchem];


    DimensionType Vk_;   // diffusion velocity (2D)
   
    BoundaryWallTemperatureType bdytemp_;
    ViscosityType Mav_;
    ConductivityType Mac_;
    DiffusionType Mad_,Mad0_,Mad1_,Mad2_,Mad3_;
    
    RatioType Theta_, Theta0_, Theta1_, Theta2_, Theta3_;
    
    BaseType yhottest_, Tignition_;

    bool useGravi_,
      useMechCompression_;

    //====RED
      VType& reconstruct_full_state(VType& zm, int i, int j, const VType& u){
      for(int k = 0; k < DataType::rednspec; ++k)
	redspecies_[k] = u[OFF(i,j,4+k)]; //Y_k, reduced composition
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

   
    //===RED
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
    //!PRESERVE MASS
    // void preserve_mass(int i, int j, VType& u){ //u is changed accoringly
    //   T sm1 = T();
      
    //   //! compute \f$ Y_{n} = 1 - \sum_{k = 1}^{n-1} Y_k\f$
    //   for(int k = 0; k < nchem-1; ++k)
    // 	sm1 += de_negativize(u[OFF(i,j,4+k)]);
      
    //   u[OFF(i,j,4+(nchem-1))] = 1. - sm1;
    // }


    // void preserve_mass(int select, int i, int j, VType& u){
    //   T sm = T();
    //   for(int k = 0; k < DataType::nspec; ++k){
    // 	if(k != select){
    // 	  sm += de_negativize(u[OFF(i,j,4+k)]);
    // 	}
    //   }
    //   u[OFF(i,j,4+select)] = 1. - sm;
    // }
    

    void is_temperature_ok(int i, int j, int k, const VType& u){
      char oper1 = '+', oper2 = '+';
      for(int l = -1; l < 2; ++l){
	for(int t = -1; t < 2; ++t){
	  int IX = i+l,
	    JX = j+t;
	  
	  if(l<0)
	    oper1 = '-';
	  else
	    oper1 = '+';
	  if(t<0)
	    oper2 = '-';
	  else
	    oper2 = '+';
	    
	  
	  if(u[OFF(IX,JX,3)] < DataType::temperature_bounds()[3*k] || 
	     u[OFF(IX,JX,3)] > DataType::temperature_bounds()[3*k+1]){
	   
	    ADONIS_ERROR(BoundsError,"Temperature bounds violated: T_{i"<<oper1<<l << ",j"<<oper2<<t << "} = "<< u[OFF(IX,JX,3)]<<", i.e. T("<<IX<<","<<JX<<") = "<< u[OFF(IX,JX,3)] << ".");
	  }
	}
      }
    }


    void project_temperature_on_its_bounds(int i, int j, int k, const VType& u){
      char oper1 = '+', oper2 = '+';
      for(int l = -1; l < 2; ++l){
	for(int t = -1; t < 2; ++t){
	  int IX = i+l,
	    JX = j+t;
	  
	  if(l<0)
	    oper1 = '-';
	  else
	    oper1 = '+';
	  if(t<0)
	    oper2 = '-';
	  else
	    oper2 = '+';
	    

	  if(u[OFF(IX,JX,3)] < DataType::temperature_bounds()[3*k] || 
	     u[OFF(IX,JX,3)] > DataType::temperature_bounds()[3*k+1]){
	     ADONIS_WARNING(Warning,"Temperature bounds violated: T_{i"<<oper1<<l << ",j"<<oper2<<t << "} = "<< u[OFF(IX,JX,3)]<<", i.e. T("<<IX<<","<<JX<<").\n   ==> PROJECT T on bounds and carry on computations...");
	   }
	  

	  Max(DataType::temperature_bounds()[3*k],Min(convert_number(u[OFF(IX,JX,3)]),DataType::temperature_bounds()[3*k+1]));

	}
      }
    }


    //CONVECTION
    //! arguments: v_l velocity component
    //! convective terms are approximated by first-order upwind formulas
    //! if the velocity component is negative
    T upwind_downwind_x(int i, int j, int t, const T& v_l, const VType& u){
      
      if(v_l <= 0) //first order upwind
	return (u[OFF(i+1,j,t)] - u[OFF(i,j,t)])/hx_;
      else 
	return (u[OFF(i,j,t)] - u[OFF(i-1,j,t)])/hx_;
    }

     T upwind_downwind_y(int i, int j, int t, const T& v_l, const VType& u){
      
      if(v_l <= 0) //first order upwind
	return (u[OFF(i,j+1,t)] - u[OFF(i,j,t)])/hy_;
      else 
	return (u[OFF(i,j,t)] - u[OFF(i,j-1,t)])/hy_;
    }
    



    //! CHEMISTRY
    //in the reduced case  only the corresponding fixed quantities
    //are extracted, e.g. molar masses.
    T wbar(int i, int j, const VType& u) const{
      T wb = T();
      for(int k = 0; k < DataType::nspec; ++k){
	wb += Tab_[GRID(i,j)][k]/DataType::molar_masses()[AccessType::spec_access(k)];
      }
      return 1./wb;
    }
    
    //!NOTE: perturbation of X is very important!!!
     T mole_fraction(int i, int j, int k, const T& Wbar, const VType& u) const{
      T mfrac = Tab_[GRID(i,j)][k]*Wbar/DataType::molar_masses()[k];
      perturbation(mfrac,1.e-16);
      return mfrac;
    }
    
    //! needed to evaluate diffusion and viscosity coefficients 
    FullChemType& mole_fraction(FullChemType& x, int i, int j, const T& Wbar, const VType& u){
      for(int k = 0; k < DataType::nspec; ++k){
	x[k] = mole_fraction(i,j,k,Wbar,u);
	perturbation(x[k],1.e-16);
      }
      return x;
    }
    
    //!overload mole fraction vector for <II>any</II> random access container
    template<class RAC>
    RAC& mole_frac(RAC& xfrac, int i, int j, const T& Wbar, const VType& u){
      for(int k = 0; k < DataType::nspec; ++k){
	xfrac[k] = mole_fraction(i,j,k,Wbar,u);
	perturbation(xfrac[k],1.e-16);
      }
      return xfrac;
    }

  

    //! overload concentration 
    template<class RAC>
    RAC& concentration(RAC& conc, int i, int j, VType& u){
      for(int k = 0; k < DataType::nspec; ++k){
	tempconc_ = u[OFF(i,j,0)]*Tab_[GRID(i,j)][k]/DataType::molar_masses()[k];  
	conc[k] = make_zero_when_reasonable(tempconc_);  //also when very slighty negative values are somehow produced
#ifndef NDEBUG
	if((conc[k] < 0.0) && (!is_zero(Abs(conc[k]))))
	  ADONIS_ERROR(ValueError,"[X_"<<k<<"] = "<< conc[k] << " is effectively smaller than ZERO.");
#endif
      }
      return conc;
    }

    //////////////////////////////////////////////////////////////////////////
    //SOME VECTORS ARE EXPLICITLY FILLED DURING CALL ==> NO NEED TO RE-COMP 
    //////////////////////////////////////////////////////////////////////////
     //! compute \f$ \dot{\omega} \f$, unit: mol/(m³s)
   
     RedChemType& chemistry(RedChemType& domeg, int i, int j, VType& u){
      concentration(C_,i,j,u); //C_ is filled now
      for(int k = 0; k < nchem; ++k){ //nchem here is reduced dimension
	//!=========TODO: in either cases: C_ must be the FULL composition
	//! compare with 
	//! '../ode/examples/automaticallygeneratedsourceterms/redh2c6_attempt.hh'
	//! only compute \f$\dot{\omega}_k \f$ when spec. \f$k\f$ participates 
	//! in any reaction; otherwise it is const and therefore  we may take
	//! \f$\dot{\omega}_k = 0.0\f$ 
	domeg[k] = ( (DataType::is_species_reactive()[AccessType::spec_access(k)]==true) ? (BCS_.net_formation_rate(BCS_.species_index(k),Tab_[GRID(i,j)][DataType::nspec],C_)) : 0.0 ); //in the reduced case C_ must be of full dimension!!!
      }
      return domeg;
    }


    
    //! \f$ \dot{\omega} \f$ is computed by now
    T heat_production(const RedChemType& dotomega,int i, int j, const VType & u){
      T h = T();
      // for(int k = 0; k < nchem; ++k)
      // 	h += BCS_.H_T(AccessType::spec_access(k),u[OFF(i,j,3)])*dotomega[k];
      int nurep = DataType::nspec-DataType::rednspec;
      for(int k = 0; k < DataType::rednspec; ++k)
	h += BCS_.H_T(AccessType::spec_access(k),Tab_[GRID(i,j)][DataType::nspec])*dotomega_[k];
      for(int k = 0; k < nurep; ++k)
	h += BCS_.H_T(Recon_.unrep_index()[k],Tab_[GRID(i,j)][DataType::nspec])*Recon_.dot_omega()[Recon_.unrep_index()[k]];
      
      return h;
    }


    T cp(int i, int j, const VType& u){
      T cres = T();
      for(int k = 0; k < DataType::nspec; ++k){
	//! fill Cp_ simultaneously
#ifndef NDEBUG
	
	if(u[OFF(i,j,3)] < DataType::temperature_bounds()[3*AccessType::spec_access(k)] || u[OFF(i,j,3)] > DataType::temperature_bounds()[3*AccessType::spec_access(k) + 1])
	  ADONIS_ERROR(BoundsError, "cp(i,j,u): Temperature out of bound T_{"<<i<<","<<j<< "} = "<< u[OFF(i,j,3)] << ".");
#endif
	Cp_[k] = BCS_.C_p(k,Tab_[GRID(i,j)][DataType::nspec]);
	cres += Cp_[k]/DataType::molar_masses()[k]*Tab_[GRID(i,j)][k];
      }
      //isCpcalculated_ = true;
      return cres;
    }

    //! at constant volume
    T cV(int i, int j, const VType& u){
       T cres = T();
       for(int k = 0; k < DataType::nspec; ++k){
#ifndef NDEBUG
	
	if(u[OFF(i,j,3)] < DataType::temperature_bounds()[3*AccessType::spec_access(k)] || u[OFF(i,j,3)] > DataType::temperature_bounds()[3*AccessType::spec_access(k) + 1])
	  ADONIS_ERROR(BoundsError, "cV(i,j,u): Temperature out of bound T_{"<<i<<","<<j<< "} = "<< u[OFF(i,j,3)] << ".");
#endif
	cres += (BCS_.C_p(k,Tab_[GRID(i,j)][DataType::nspec]) - PhysicalConstants<T>::Rgas)/DataType::molar_masses()[k]*Tab_[GRID(i,j)][k];
      }
      return cres;
    }

    //! EVALUATE TRANSPORT COEFFICIENTS AT point node \f$ (i,j) \f$
    //! compute \f$ \mu(T,X)\f$
    T mu(int i, int j, const VType& u){  
      mole_fraction(X_,i,j,wbar(i,j,u),u);
      return Mav_.compute_mixture_averaged_viscosity(Tab_[GRID(i,j)][DataType::nspec],X_);
    }

    //! \f$ \lambda(\rho, p, T, X) \f$
    T lambda(int i, int j, const VType& u){
      (*this).mole_fraction(X_,i,j,wbar(i,j,u),u);
      perturbation(X_,1.e-16);
      return Mac_.compute_mixture_averaged_conductivity(u[OFF(i,j,0)], pressure(i,j,u_), Tab_[GRID(i,j)][DataType::nspec],X_);
    }

    T heat_conduction(int i, int j, const VType& u){
      return ( (0.5*(lambda(i+1,j,u)+lambda(i,j,u))*(Tab_[GRID(i+1,j)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(lambda(i,j,u)+lambda(i-1,j,u))*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i-1,j)][DataType::nspec]))/ntimes<2>(hx_) + (0.5*(lambda(i,j+1,u)+lambda(i,j,u))*(Tab_[GRID(i,j+1)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(lambda(i,j,u)+lambda(i,j-1,u))*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i,j-1)][DataType::nspec]))/ntimes<2>(hy_) );
    }

    //=======================================================================
    //======TODO:  ADAPT FURTHER with Tab_[·][·] 
    //======       I am too tired by now   :/    ============================
    //=======================================================================

    //! calculate \f$V_{\textrm{corr}} := -V_{\mathrm{c}} = \sum_{k=1}^K Y_kV_k \f$
    //!note: \f$V_{\textrm{corr}}\f$ is a <B>constant</B> correction vector, i.e. independent of species but varying in space and time
    DimensionType& V_corr(DimensionType& v, int i, int j, const VType& u){
      
#ifndef NO_CORRECTION_DIFFUSION_VELOCITY
      v = T(); //reset
      

      fivepointstencil_[0] = wbar(i,j-1,u);
      fivepointstencil_[1] = wbar(i+1,j,u);
      fivepointstencil_[2] = wbar(i,j+1,u);
      fivepointstencil_[3] = wbar(i-1,j,u);
      fivepointstencil_[4] = wbar(i,j,u);

      mole_frac(X0_,i,j-1,fivepointstencil_[0],u);  //(i,j-1)
      mole_frac(X1_,i+1,j,fivepointstencil_[1],u);  //(i+1,j)
      mole_frac(X2_,i,j+1,fivepointstencil_[2],u);  //(i,j+1)
      mole_frac(X3_,i-1,j,fivepointstencil_[3],u);  //(i-1,j)
      mole_frac(X_,i,j,fivepointstencil_[4],u);       //(i,j)

      Mad0_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j-1,u),Tab_[GRID(i,j-1)][DataType::nspec],X0_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(pressure(i+1,j,u),Tab_[GRID(i+1,j)][DataType::nspec],X1_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j+1,u),Tab_[GRID(i,j+1)][DataType::nspec],X2_);
      
      Mad3_.compute_mixture_averaged_diffusion_coefficients(pressure(i-1,j,u),Tab_[GRID(i-1,j)][DataType::nspec],X3_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j,u),Tab_[GRID(i,j)][DataType::nspec],X_);

#ifdef CORRECTION_DIFFUSION_VELOCITY_WITH_THERMAL_PART
      //! thermal diffusion ratios
      Theta0_.compute_thermal_diffusion_ratios(Tab_[GRID(i,j-1)][DataType::nspec],X0_);
      Theta1_.compute_thermal_diffusion_ratios(Tab_[GRID(i+1,j)][DataType::nspec],X1_);
      Theta2_.compute_thermal_diffusion_ratios(Tab_[GRID(i,j+1)][DataType::nspec],X2_);
      Theta3_.compute_thermal_diffusion_ratios(Tab_[GRID(i-1,j)][DataType::nspec],X3_);
      Theta_.compute_thermal_diffusion_ratios(Tab_[GRID(i,j)][DataType::nspec],X_);
#endif

      for(int k = 0; k < DataType::nspec; ++k){
      	//!x-direction
      	v[0] += ( -0.5*(DataType::molar_masses()[k]*Mad1_[k]/fivepointstencil_[1] + DataType::molar_masses()[k]*Mad_[k]/fivepointstencil_[4])*(X1_[k]-X_[0])/hx_  //approximation of \f$ Y_k\mathcal{V}_k\f$ in \f$x\f$-direction
#ifdef CORRECTION_DIFFUSION_VELOCITY_WITH_THERMAL_PART
		  -0.5*(DataType::molar_masses()[k]*Mad1_[k]/fivepointstencil_[1]*Theta1_[k]/Tab_[GRID(i+1,j)][DataType::nspec] + DataType::molar_masses()[k]*Mad_[k]/fivepointstencil_[4]*Theta_[k]/Tab_[GRID(i,j)][DataType::nspec])*(Tab_[GRID(i+1,j)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec])/hx_  //approximation of \f$ Y_k\mathcal{W}_k\f$ in \f$x\f$-direction
#endif
      			  );  

      	//!y-direction
      	v[1] += ( -0.5*(DataType::molar_masses()[k]*Mad2_[k]/fivepointstencil_[2] + DataType::molar_masses()[k]*Mad_[k]/fivepointstencil_[4])*(X2_[k]-X_[0])/hy_  //approximation of \f$ Y_k\mathcal{V}_k\f$ in \f$y\f$-direction
#ifdef CORRECTION_DIFFUSION_VELOCITY_WITH_THERMAL_PART
		  -0.5*(DataType::molar_masses()[k]*Mad2_[k]/fivepointstencil_[2]*Theta2_[k]/Tab_[GRID(i,j+1)][DataType::nspec] + DataType::molar_masses()[k]*Mad_[k]/fivepointstencil_[4]*Theta_[k]/Tab_[GRID(i,j)][DataType::nspec])*(Tab_[GRID(i,j+1)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec])/hy_  //approximation of \f$ Y_k\mathcal{W}_k\f$ in \f$y\f$-direction
#endif
		  );    
      }
#endif //in case of no definition of 'NO_CORRECTION_DIFFUSION_VELOCITY'
      return v;
      
    }
   
   

    //perhaps much better to replace \f$ Y_k/X_k\f$ by \f$W_k/\overline{W}\f$ since \f$X_k\f$ can be zero
    FullChemType& species_diffu_term(FullChemType& diffu, int i, int j, const VType& u){
      fivepointstencil_[0] = wbar(i,j-1,u);
      fivepointstencil_[1] = wbar(i+1,j,u);
      fivepointstencil_[2] = wbar(i,j+1,u);
      fivepointstencil_[3] = wbar(i-1,j,u);
      fivepointstencil_[4] = wbar(i,j,u);

      mole_frac(X0_,i,j-1,fivepointstencil_[0],u);  //(i,j-1)
      mole_frac(X1_,i+1,j,fivepointstencil_[1],u);  //(i+1,j)
      mole_frac(X2_,i,j+1,fivepointstencil_[2],u);  //(i,j+1)
      mole_frac(X3_,i-1,j,fivepointstencil_[3],u);  //(i-1,j)
      mole_frac(X_,i,j,fivepointstencil_[4],u);       //(i,j)

      Mad0_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j-1,u),Tab_[GRID(i,j-1)][DataType::nspec],X0_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(pressure(i+1,j,u),Tab_[GRID(i+1,j)][DataType::nspec],X1_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j+1,u),Tab_[GRID(i,j+1)][DataType::nspec],X2_);
      
      Mad3_.compute_mixture_averaged_diffusion_coefficients(pressure(i-1,j,u),Tab_[GRID(i-1,j)][DataType::nspec],X3_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j,u),Tab_[GRID(i,j)][DataType::nspec],X_);

      
      //! thermal diffusion ratios
      Theta0_.compute_thermal_diffusion_ratios(Tab_[GRID(i,j-1)][DataType::nspec],X0_);
      Theta1_.compute_thermal_diffusion_ratios(Tab_[GRID(i+1,j)][DataType::nspec],X1_);
      Theta2_.compute_thermal_diffusion_ratios(Tab_[GRID(i,j+1)][DataType::nspec],X2_);
      Theta3_.compute_thermal_diffusion_ratios(Tab_[GRID(i-1,j)][DataType::nspec],X3_);
      Theta_.compute_thermal_diffusion_ratios(Tab_[GRID(i,j)][DataType::nspec],X_);

      //!note V_corr is a constant correction vector, i.e. independent of species but varying in space and time
      V_corr(Vk_,i,j,u);       //correction
      
     
      for(int k = 0; k < DataType::nspec; ++k){
	diffu[k] =( (0.5*(u[OFF(i+1,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[1]*Mad1_[k] + u[OFF(i,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_[k])*(X1_[k]-X_[k]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_[k] + u[OFF(i-1,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[3]*Mad3_[k])*(X_[k]-X3_[k]))/ntimes<2>(hx_) + 
	  (0.5*(u[OFF(i,j+1,0)]*DataType::molar_masses()[k]/fivepointstencil_[2]*Mad2_[k] + u[OFF(i,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_[k])*(X2_[k]-X_[k]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_[k] + u[OFF(i,j-1,0)]*DataType::molar_masses()[k]/fivepointstencil_[0]*Mad0_[k])*(X_[k]-X0_[k]))/ntimes<2>(hy_)
	  //now comes the part concerning thermal diffusion...
		    + (0.5*(u[OFF(i+1,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[1]*Mad1_[k]*Theta1_[k]/Tab_[GRID(i+1,j)][DataType::nspec] + u[OFF(i,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_[k]*Theta_[k]/Tab_[GRID(i,j)][DataType::nspec])*(Tab_[GRID(i+1,j)][DataType::nspec] - Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_[k]*Theta_[k]/Tab_[GRID(i,j)][DataType::nspec] + u[OFF(i-1,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[3]*Mad3_[k]*Theta3_[k]/Tab_[GRID(i-1,j)][DataType::nspec])*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i-1,j)][DataType::nspec]))/ntimes<2>(hx_) //x-part
		    + (0.5*(u[OFF(i,j+1,0)]*DataType::molar_masses()[k]/fivepointstencil_[2]*Mad2_[k]*Theta2_[k]/Tab_[GRID(i,j+1)][DataType::nspec] + u[OFF(i,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_[k]*Theta_[k]/Tab_[GRID(i,j)][DataType::nspec])*(Tab_[GRID(i,j+1)][DataType::nspec] - Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[k]/fivepointstencil_[4]*Mad_[k]*Theta_[k]/Tab_[GRID(i,j)][DataType::nspec] + u[OFF(i,j-1,0)]*DataType::molar_masses()[k]/fivepointstencil_[0]*Mad0_[k]*Theta0_[k]/Tab_[GRID(i,j-1)][DataType::nspec])*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i,j-1)][DataType::nspec]))/ntimes<2>(hy_)  //y-part 
#ifndef NO_CORRECTION_DIFFUSION_VELOCITY
	  //...and the correcting term; note: \f$ Y_k = X_k W_k/\overline{W}\f$
	  //... as well as a constant Vk_
	  + 
	  ( Vk_[0]*(u[OFF(i,j,0)]*X_[k]*DataType::molar_masses()[k]/fivepointstencil_[4] - u[OFF(i-1,j,0)]*X3_[k]*DataType::molar_masses()[k]/fivepointstencil_[3])/hx_
	    + Vk_[1]*(u[OFF(i,j,0)]*X_[k]*DataType::molar_masses()[k]/fivepointstencil_[4] - u[OFF(i,j-1,0)]*X0_[k]*DataType::molar_masses()[k]/fivepointstencil_[0])/hy_ )
#endif
		    );

      }
      
      return diffu;
    }


     //! assumes that cp(i,j,u) has been invoked previously
     T temperature_diffu_term(int i, int j, const VType& u){
       T td = T(),
	 meanmolmass = wbar(i,j,u);

       mole_fraction(X_,i,j,meanmolmass,u); //(i,j)
       Mad_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j,u_),Tab_[GRID(i,j)][DataType::nspec],X_); //diffusion coefficients are computed by now
       
        //!store counterclockwise, starting from the bottom with
	//!(i,j) being the last entry
        fivepointstencil_[0] = wbar(i,j-1,u);
	fivepointstencil_[1] = wbar(i+1,j,u);
	fivepointstencil_[2] = wbar(i,j+1,u);
	fivepointstencil_[3] = wbar(i-1,j,u);
	fivepointstencil_[4] = meanmolmass;
	
	mole_frac(X0_,i,j-1,fivepointstencil_[0],u);  //(i,j-1)
	mole_frac(X1_,i+1,j,fivepointstencil_[1],u);  //(i+1,j)
	mole_frac(X2_,i,j+1,fivepointstencil_[2],u);  //(i,j+1)
	mole_frac(X3_,i-1,j,fivepointstencil_[3],u);  //(i-1,j)

       //!note V_corr is a constant correction vector, i.e. independent of species but varying in space and time
	V_corr(Vk_,i,j,u);       //correction


	//! thermal diffusion ratios
	Theta0_.compute_thermal_diffusion_ratios(Tab_[GRID(i,j-1)][DataType::nspec],X0_);


       //isDiffusioncalculated_ = true;
       
	for(int k = 0; k < DataType::nspec; ++k){
   
	//! assume Cp_ has been filled previously
	  td += ( u[OFF(i,j,0)]*Cp_[k]*( Mad_[k]/meanmolmass*( (0.5*(X1_[k]+X_[k])*(Tab_[GRID(i+1,j)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(X3_[k]+X_[k])*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i-1,j)][DataType::nspec]))/ntimes<2>(hx_) + (0.5*(X2_[k]+X_[k])*(Tab_[GRID(i,j+1)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(X_[k]+X0_[k])*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i,j-1)][DataType::nspec]))/ntimes<2>(hy_) 
       //thermal diffusion part...
	+ 						     
							       Theta_[k]/Tab_[GRID(i,j)][DataType::nspec]*( (0.5*(Tab_[GRID(i+1,j)][DataType::nspec]+Tab_[GRID(i,j)][DataType::nspec])*(Tab_[GRID(i+1,j)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(Tab_[GRID(i,j)][DataType::nspec]+Tab_[GRID(i-1,j)][DataType::nspec])*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i-1,j)][DataType::nspec]) )/ntimes<2>(hx_) + (0.5*(Tab_[GRID(i,j+1)][DataType::nspec]+Tab_[GRID(i,j)][DataType::nspec])*(Tab_[GRID(i,j+1)][DataType::nspec]-Tab_[GRID(i,j)][DataType::nspec]) - 0.5*(Tab_[GRID(i,j)][DataType::nspec]+Tab_[GRID(i,j-1)][DataType::nspec])*(Tab_[GRID(i,j)][DataType::nspec]-Tab_[GRID(i,j-1)][DataType::nspec]) )/ntimes<2>(hy_)) ) //!multiplication with \f$D_k^{\mathrm{mix}}/\overline{W}\f$
 //corrective part; note: \f$Y_k/W_k = X_k/\overline{W} \f$
#ifndef NO_CORRECTION_DIFFUSION_VELOCITY
					+ X_[k]/meanmolmass*(Vk_[0]*(Tab_[GRID(i,j)][DataType::nspec] - Tab_[GRID(i-1,j)][DataType::nspec])/hx_ + Vk_[1]*(Tab_[GRID(i,j)][DataType::nspec] - Tab_[GRID(i,j-1)][DataType::nspec])/hy_)						     //multiplication with \f$Y_k/W_k\f$
#endif
					)//multiplic. with \f$\rho C_{pk}^0\f$ 
		 );
      }
      return td;
     }
  

    T& make_zero_when_reasonable(T& val){
      return ( ((val < 0.0) && (is_zero(Abs(val),near_zero_))) ? (val = 0.0) : val );
    }

    typename VType::value_type& make_zero_when_reasonable(int i, int j, VType& u, int primvar){
      return ( ((u[OFF(i,j,primvar)] < 0.0) && (is_zero(Abs(u[OFF(i,j,primvar)]),near_zero_))) ? (u[OFF(i,j,primvar)] = 0.0) : u[OFF(i,j,primvar)] );
    }


    T pressure(int i, int j, const VType& u){
      return ( (u[OFF(i,j,0)]*PhysicalConstants<BaseType>::Rgas*Tab_[GRID(i,j)][DataType::nspec])/wbar(i,j,u) );
    }

    //! see Majda and Braack
    //! \f$ p = \rho^{\gamma}/\exp(-S_0), \f$ where \f$\gamma\f$ is the ratio
    //! of the specific heat capacities
    T hydrodynamical_pressure(int i, int j, const VType& u){
      T gamma = cp(i,j,u)/cV(i,j,u);  //CORRECT ????
      T S0 = T();
      for(int k = 0; k < DataType::nspec; ++k)
	S0 += BCS_.S_T(k,Tab_[GRID(i,j)][DataType::nspec]);   //CORRECT ????
      
      return ( pow(u[OFF(i,j,0)],gamma)/exp(-S0)   );
    }

    //used before: 1/rho*((p(i,j,u)-p(i-1,j,u))/hx) 
    T mechanical_compression_x(int i, int j, const VType& u){
      return ((useMechCompression_==true) ? 1./u[OFF(i,j,0)]*(pressure(i+1,j,u) - pressure(i-1,j,u))/(2*hx_) : 0.);  //via central differences
    }

    //used before: 1/rho*((p(i,j,u)-p(i,j-1,u))/hy) 
    T mechanical_compression_y(int i, int j, const VType& u){
      return ((useMechCompression_==true) ? 1./u[OFF(i,j,0)]*(pressure(i,j+1,u) - pressure(i,j-1,u))/(2*hy_) : 0.); //via central differences
    }

    T gravitational_force_x(int i, int j, const VType& u){
      return ( (useGravi_==true) ? gravi1_  : 0.);
    }

    T gravitational_force_y(int i, int j, const VType& u){
      return ( (useGravi_==true) ? gravi2_  : 0.);
    }
  };

} //end namespace

#endif 
