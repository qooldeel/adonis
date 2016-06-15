#ifndef NAV_STOKES_VIA_MOL_FOR_CONAIRES_PRIMARY_REFERENCE_FUEL_H2_MECH_HH
#define NAV_STOKES_VIA_MOL_FOR_CONAIRES_PRIMARY_REFERENCE_FUEL_H2_MECH_HH
/**
 * This functor can be used by an ODE integrator in th usual way
*/
#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"
#include "../common/typeadapter.hh"
#include "../common/isclass.hh"
#include "../graphics/printmoldata.hh"
#include "../io/readinparameters.hh"
#include "../templatemetaprograms/unrollloop.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../ode/ode.hh"
#include "../ode/additional/boundarytemperaturedistribution.hh"
#include "../ode/additional/parabolicvelocityprofile.hh"
#include "../containers/staticarray.hh"

//ODE
#include "../ode/constants.hh"

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

#include "../ode/examples/automaticallygeneratedsourceterms/mechtype.hh"

#include "gensettings.hh"
#include "usefulextensions.hh"

#include "bdyhandler.hh"

#include "primitivevars.hh"
#include "../ode/repair.hh"

#include "../ode/functorid.hh"

#ifdef USE_OPENMP
#include "omp.h"
#endif

namespace Adonis{

  template<class T>
  class ConaireH2ReactiveNS2D{
  public:
    typedef T value_type;
    typedef ConaireH2ReactiveNS2D<T> ThisType;
    typedef ExprTmpl::MyVec<T> VType; 
    typedef typename VType::iterator IteratorTType;
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<DType> VDType; 
    //typedef typename TypeAdapter<T>::Type DType; //no CppAd
    typedef typename TypeAdapter<T>::BaseType BaseType; //for norms, tolerances 
    typedef std::size_t SizeType;

#ifndef GHOST_POINTS_INCLUDED //not necessary when ghost points are used
    //! mol identifer
    typedef MOLFunctorIdentifier<2,FunctorID::H2conaireFull> mol_type;
#endif
    
    typedef BoundaryFunction2D<FDMSettings::orderOfBdyCond,ThisType> BdyFctType;
    
    //needed for thermo chemistry: primary reference fuel mechanism
    typedef ThermoData4Mechanism<BaseType,10> DataType;  

    enum{
      MECHANISM = DataType::ENCOD
    };

    
    enum{nchem = DataType::nspec, 
	 //rho,v1,v2,T + specs
	 nprim = 4+DataType::nspec};  //4+nchem = number of phys quantities
    
    typedef StaticArray<T,DataType::nspec> ArrayType; 
    typedef StaticArray<T,2> DimensionType;
    //! element 0 is x-direction, 1 y-direction and so forth
    
    typedef BoundaryTemperatureDistribution<BaseType,FDMSettings::WallTemperatureType> BoundaryWallTemperatureType;

    typedef InletVelocity<BaseType,BaseType,BaseType,FDMSettings::InletVelocityType> InVeloType;


    enum{O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR};

    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;

    //Other Index
    typedef IndexHandler<DataType,false> AccessType;
    typedef typename AccessType::IndexerType IndexerType;

    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,IndexerType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    typedef EquationOfState<'i'> EosType;

    //TRANSPORT COEFFICIENTS
    typedef MixtureAveragedViscosity<DataType,false> ViscosityType;
    typedef MixtureAveragedConductivity<DataType,false> ConductivityType;
    typedef MixtureAveragedDiffusion<DataType,false> DiffusionType;

    typedef ThermalDiffusionRatio<DataType,false,FDMSettings::ThermalDiffusion4LightWeightSpeciesOnly> RatioType;
    
    // this may be of some importance
    //typedef ComputeTransportProperties<false> PropType;

    ConaireH2ReactiveNS2D(SizeType dim = 0):rhs_(dim){
#if defined ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS && !defined GHOST_POINTS_INCLUDED
      ADONIS_ERROR(DefinitionError, "You can only use define 'ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS' in conjunction with define 'GHOST_POINTS_INCLUDED'. Check file \"data/conaireh2.dat\" and correct it.");
#endif
      
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
      // can also be '/data/conaireh2.dat' but this is better when invoked from
      // other parts of the code
      PD.read_from_file("../fdm/data/conaireh2.dat");
       //! domain [a,b] x [c,d]
      a_ = PD.get_datum<BaseType>("a");
      b_ = PD.get_datum<BaseType>("b");
      c_ = PD.get_datum<BaseType>("c");
      d_ = PD.get_datum<BaseType>("d");

      nx_ = PD.get_datum<int>("Nx");
      ny_ = PD.get_datum<int>("Ny");
      npt_ = nx_*ny_; //all points, i.e. \f$\Omega_h \cup \partial\Omega_h\f$
      totpoints_ = nprim*npt_;
      length_ = b_-a_;   //length in x-direction of rectangular domain
      hx_ = length_/
#ifndef GHOST_POINTS_INCLUDED
	(nx_-1)
#else
	(nx_-3)
#endif
	;
      diam_ = d_-c_;
      hy_ = diam_/
#ifndef GHOST_POINTS_INCLUDED
	(ny_-1)
#else
	(ny_-3)
#endif
	; 
      t0_ = PD.get_datum<BaseType>("t0");
      tend_ = PD.get_datum<BaseType>("tend");
      halfdiam_ = 0.5*diam_;
      gravi1_ = PD.get_datum<BaseType>("g1");
      gravi2_ = PD.get_datum<BaseType>("g2");
      near_zero_ = Abs(PD.get_datum<BaseType>("near_zero"));
      rhoMin_ = PD.get_datum<BaseType>("rho_min");
      smallVal_ = Abs(PD.get_datum<BaseType>("small_val"));
      p0_ = PD.get_datum<BaseType>("pconst");
      v1_in_ = PD.get_datum<BaseType>("v1_in");
      v2_in_ = PD.get_datum<BaseType>("v2_in");
      v1_wall_ =  PD.get_datum<BaseType>("v1_wall"); //wall velocity
      v2_wall_ = PD.get_datum<BaseType>("v2_wall");  //wall velocity
      T_wall_ = PD.get_datum<BaseType>("T_wall");
      T_in_ = PD.get_datum<BaseType>("T_min");
      x_ignition_ = PD.get_datum<BaseType>("x_ignition");//(b_-a_)/PD.get_datum<BaseType>("channel_length_frac");
      
      adonis_assert(x_ignition_ < length_);

      Tignition_ = PD.get_datum<BaseType>("T_ignition");
      yhottest_ = PD.get_datum<BaseType>("hottest_y");
      bdytemp_.initialize(T_wall_,Tignition_,x_ignition_, (((BoundaryWallTemperatureType::Value == 's') || (BoundaryWallTemperatureType::Value == 'S')) ? PD.get_datum<BaseType>("cubic_spline_from_ignition")*hx_ : PD.get_datum<BaseType>("plateau_width")*hx_),(((BoundaryWallTemperatureType::Value == 'a') || (BoundaryWallTemperatureType::Value == 'A')) ? 2*hx_ : length_)); //4th argument is here only meaningful when 'm' profile is used

      useGravi_ = PD.get_datum<bool>("gravForce");
      useMechCompression_ = PD.get_datum<bool>("mechCompress");
      pressFac_ = PD.get_datum<unsigned int>("FactorpOut");
      approxEqual_ = PD.get_datum<BaseType>("approx_equal");

      WbarBdy_ = pBdy_ = rhoBdy_ = TBdy_ = T(); //just intialize somehow
      td_ = meanmolmass_ = tmpvar_ = T();
      YBdy_ = T();
      CorrectYfrac_.initialize(YBdy_.size());

      //max. allowed values
      v1_max_ = PD.get_datum<BaseType>("v1_max");
      v2_max_ = PD.get_datum<BaseType>("v2_max");
      tooLarge_ =  PD.get_datum<BaseType>("too_large_value");

      //! adapt indices for ghost point treatment and normal domain treatment
      //! accordingly
#ifndef GHOST_POINTS_INCLUDED
      leftbeg_ = 1;
      leftend_ = ny_-2;
      lefti_ = 0;
      downbeg_ = 0;
      downend_ = nx_;
      downj_ = 0;
      upj_ = ny_-1;
      rightbeg_ = 1;
      rightend_ = ny_-2;
      righti_ = nx_-1;
      offs_ = 0;
#else  //ghost points used
      leftbeg_ = 2;
      leftend_ = ny_-3;
      lefti_ = 1;
      downbeg_ = 1;
      downend_ = nx_-1;
      downj_ = 1;
      upj_ = ny_-2;
      rightbeg_ = 2;
      rightend_ = ny_-3;
      righti_ = nx_-2;
      offs_ = 1;
#endif
      

      //! output at construction. To avoid double output when CppAD is used
      //! only print to screen when the constructor of class ODE is
      //! invoked.
      if(are_types_equal<BaseType,T>() == true){
	std::cout << "Numerical types in use: BaseType = "<< typeid(BaseType).name() << "    T = "<< typeid(T).name() << std::endl << std::endl;
	std::cout << "BOUNDARY CONDITION ORDER: " << FDMSettings::orderOfBdyCond << std::endl;
	if((FDMSettings::orderOfBdyCond == 2) && (FDMSettings::adjust2ndOrderNeumann == true))
	  std::cout << "Enforce mass balance of Y at the boundary for the second order zero Neumann condition." << std::endl;


#ifndef NO_CORRECTION_DIFFUSION_VELOCITY
	std::cout << "COMPUTE CORRECTION VELOCITY"<< std::endl;
	std::cout << "----------------------------"<<std::endl<<std::endl;
#else
	std::cout << "NO CORRECTION VELOCITY COMPUTED" << std::endl;
#endif

	std::cout << std::endl;
      
#ifdef LOW_MACH_NUMBER_FORMULATION
	std::cout << "Compute low-Mach-number formulation for mixture density RHO" << std::endl << std::endl;
#else
	std::cout << "Compute d_t rho + D · (rho v) as usual"<< std::endl << std::endl;
#endif
	std::string ghat;
	if(useGravi_)
	  ghat = "   g = ["+Num2str(gravi1_) + ", "+Num2str(gravi2_)+"]";
	std::cout <<  "gravitational force: " << yes_no(useGravi_) << ghat <<"  ||   mechanical compression: "<< yes_no(useMechCompression_) << std::endl;
	if(useMechCompression_==true){
#ifndef V1_CONSTANT 
#ifndef V2_CONSTANT
	  std::cout << "Mechanical compression based on ";
#ifdef HYDRO_PRESSURE
	  std::cout << "HYDRODYNAMIC pressure" << std::endl;
#else
	  std::cout << "THERMODYNAMIC pressure" << std::endl;
#endif
#endif
#endif
	}
	
#ifdef PRESSURE_CONST
	std::cout << "Use constant THERMODYNAMIC pressure p0 = "<< p0_ << " Pa." << std::endl;
#else
	std::cout << "ALGEBRAIC THERMODYNAMIC PRESSURE p = (rho*Rgas*T)/Wbar" << std::endl;
#endif

#ifdef IGNITION_AT_HOT_WALLS
	std::cout << std::endl << "Ignition at hot walls; temperature profile: "<< FDMSettings::WallTemperatureType << std::endl<< std::endl;
#else
	std::cout << std::endl << "Ignition left side" << std::endl<< std::endl;
#endif


#ifdef V1_CONSTANT
	std::cout << "VELOCITY: v1 = const" <<  std::endl << std::endl;
#endif
#ifdef V2_CONSTANT
	std::cout << "VELOCITY: v2 = const" <<  std::endl << std::endl;
#endif

#ifdef RHO_ALGEBRAIC
	ADONIS_ERROR(DerivedError,"Uncomment this 'ADONIS_ERROR' to support algebraic computation of rho instead of differential");
#ifdef PRESSURE_CONST
	std::cout << "ALGEBRAIC RHO: rho = (p0·Wbar)/(Rgas·T)" << std::endl;
#else
	ADONIS_ERROR(CoverageError,"This branch is not covered: Algebraic computation of rho MUST be accompanied by constant pressure!\n    Please uncomment \\include '#define PRESSURE_CONST'");
#endif
#else
	std::cout << "DIFFERENTIAL RHO" << std::endl;
#endif
	std::cout << "Using "<< FDMSettings::orderOfBdyCond << ".-order ZERO NEUMANN condition"<< std::endl;

#ifdef NONCONSERVATIVE_FORM
	std::cout << "Using nonconservative form of the Navier-Stokes equations." << std::endl;
#else
	std::cout << "Using CONSERVATIVE form of the Navier-Stokes equations." << std::endl;  
#endif
#ifdef ORDER_2_CONVECTIVE_TERM
	std::cout << "Use central differences for convective terms (also in the continuity equation)" << std::endl;
#endif
#ifdef PRESSURE_GRADIENT_ORDER_1_FD
	std::cout << "Use 1st order FDs for pressure gradient"  << std::endl;
#endif

	std::cout << "ODE: " << ( (Constant<BaseType>::useSTSZctrl == true) ? ("STEPSIZECTRL") : ("Equidistant stepsize") ) << std::endl;
	std::cout << "ODE: " << ( (Constant<BaseType>::additionConvergenceCheck == true) ? ("Perform ADDITIONAL convergence check") : ("no additional convergence check applied.") ) << std::endl;

	std::cout << "hx = "<< hx_ << "    hy = " << hy_ << std::endl;
	std::cout<< "************************************************"<<std::endl;
	std::cout << "Integration method: ";
#ifndef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
	std::cout << " implicit "<< ((Constant<BaseType>::accuracyOfMethod == 1) ? "EULER." : ((Constant<BaseType>::accuracyOfMethod == 2) ? "TRAPEZOIDAL." : "unknown.")) << std::endl;
	std::cout << ( (Constant<BaseType>::dampMe == true) ? ("---> DAMPED NEWTON based on violation of T and Y_k applied.") : ("---> undamped Newton used") ) << std::endl;
	std::cout << ( (Constant<BaseType>::useSimplifiedNewt == true) ? ("---> SIMPLIFIED NEWTON used.") : ("---> conventional Newton used (calculate Jacobian in each iteration anew)") ) << std::endl;
#else
	FM_.nice_output("2ND ORDER ROSENBROCK METHOD.",36);
#endif
	std::cout<< "************************************************"<<std::endl;

#ifdef SOLVE_FOR_K_MINUS_1_CONSERVATION_EQS
	FM_.nice_output("Solve for K-1 species, Y_excess = 1. - sum_{k != excess}Y_k.", 35);
#endif
	
#ifdef GHOST_POINTS_INCLUDED
	std::cout << "Booooohhhooooohooooohooooohooooo....." <<std::endl;
	std::cout << "++++++ Grid includes GHOST POINTS.... Don't call the Ghostbusters ;)" << std::endl;
#ifdef REFLECTIVE_BDY_CONDITIONS
	std::cout << "   --> Reflective boundary conditions used to describe ghost points" << std::endl;
#else
	std::cout << "   --> Copy boundary conditions to ghost points" << std::endl << std::endl;
#endif
#ifdef ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS
	std::cout << "   Moreover, rhs_ = 0 at bdy where explicit values are given, e.g. at the left bdy for all species or at the upper and lower bdy for v1, v2 and T" << std::endl;
#endif
#endif
	
#ifdef USE_OPENMP
	FM_.nice_output("OpenMP used for PARALLELIZATION.",35);
#endif

      } //are_types_equal -- only output when no CppAD type

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
  


      v1bdy_.resize(ny_);   //for the velocity profile (along y direction)
      u_.resize(npt_);
      cp_ = T();
      tempconc_ = T();
      tempDiff_ = T();
      specDiff_ = T();
      tempHeat_ = T();
      specHeat_ = T();
      tempCond_ = T();
      sp_ = T();
     

      //O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR
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
      veloprofile_.resize(
#ifndef GHOST_POINTS_INCLUDED			  
			  ny_
#else
			  ny_-2
#endif
			  ); //along y-axis
      BaseType r = -halfdiam_;  //R = d/2, r = -R
  
      for(int j = 0; j < (int)veloprofile_.size(); ++j)

	{
        veloprofile_[j] = InVeloType::velocity(r,halfdiam_,v1_in_); //zero(parabolic_inlet_velocity(r,halfdiam_,v1_in_)); 
	r += hy_;
      }

      //fixed smooth temperature profile 
      BaseType xdir(0.);
      walltempprofile_.resize(
#ifndef GHOST_POINTS_INCLUDED
			      nx_
#else
			      nx_-2
#endif
			      );  //along x-axis
    for(int i = 0; i < (int)walltempprofile_.size(); ++i)

	{
	xdir = i*hx_;
	walltempprofile_[i] = bdytemp_(xdir);
      }


      wbarIn_ = BaseType();
      for(int k = 0; k < nchem; ++k)
	wbarIn_ += Y_in_[k]/DataType::molar_masses()[AccessType::spec_access(k)];
      if(is_zero(wbarIn_))
	 wbarIn_ = smallVal_;
      else
	wbarIn_ = 1./wbarIn_;

      rhoIn_ = (p0_*wbarIn_)/(PhysicalConstants<BaseType>::Rgas*T_in_);

      meanVal_ = T();
      pressure_ = p0_;
      useMe_ = 0;       //a temporary variable that can be used throughout
      Veciter_ = rhs_.begin(); //initialization
      // for(SizeType l = 0; l < rhs_.size(); ++l){
      // 	Veciter_[l] = p0_;    //initialization
      // }
      isBeginning_ = true;    

      stepsize_ = T();
      isPrevPressInvoked_=false;
      std::cout << "CFD H2 object created"<< std::endl;
    } //end constructor

    //! \f$\rho, v_1, v_2, T\f$ + nchem species, discretized on \f$\Omega_h\f$
      //! with npt_ nodes
    int dim() const {return totpoints_;}
    int domain_dim() const {return totpoints_;}

    std::string name() const {return "with momentum";}

    template<class X>
    void set_boundary(X& w){
      //!BOUNDARY·BOUNDARY·BOUNDARY·BOUNDARY·BOUNDARY·BOUNDARY·BOUNDARY
      //!BOUNDARY -- imagine szenario is rotated by 90°
      //!BOUNDARY·BOUNDARY·BOUNDARY·BOUNDARY·BOUNDARY·BOUNDARY·BOUNDARY
      
      //! Note: continuity equation only contains 1st order spatial derivatives.
      //!       Hence no boundary conditions are needed, cf. [ANDERSON, §10.3.5
      //!       p. 457]
      //! NOTE: The boundary values for the conservation variables (either given
      //!       in primitive or conservative form) are, unless stated already
      //!       explicitly, solved explicitly for the corresp. boundary value. 
      //!       This value is then copied to the array w.
      //!       Therefore, no complicated ODEs based on specialized FDs are 
      //!       required at the boundaries and we can put there simply 
      //!       rhs_ = 0.0.
      //!  LEFT: inlet
      DType y(0.);
      for(int j = leftbeg_; j <= leftend_; ++j){ //without lower and upper edge!
	y = j*hy_;

#ifdef IGNITION_AT_HOT_WALLS
	rhoBdy_ = rhoIn_;
	TBdy_ = T_in_;  	
#else //ignition at left side
  	
	if((i==lefti_) && (is_contained(yhottest_,y-hy_,y+hy_))){
	  // rho now alters with temperature and composition
	  //here p0 is the atmosperic pressure, see ANDERSON, p. 458, FIG. 10.5]  
	  rhoBdy_ = (p0_*WbarIn_)/(PhysicalConstants<BaseType>::Rgas*Tignition_);
	  TBdy_ = Tignition_;
	}
	else{
	  TBdy_ = T_in_;
	}

#endif //END IGNITION_AT_HOT_WALLS
	//applies to every physical quantity independent of ignition zone
	w[OFF(lefti_,j,0)] = rhoBdy_;  //                times  1 or w[OFF(0,j,0)]
	w[OFF(lefti_,j,1)] = veloprofile_[j-offs_]*FlowVariables<VType>::rho(lefti_,j,w,nx_);
	w[OFF(lefti_,j,2)] = v2_in_*FlowVariables<VType>::rho(lefti_,j,w,nx_);
	w[OFF(lefti_,j,3)] = TBdy_*FlowVariables<VType>::rho(lefti_,j,w,nx_);


	// at bdy there's no dependence on w, hence all derivatives w.r.t.
	// w vanish necessarily. Also prevents badly filled in results
#ifndef  GHOST_POINTS_INCLUDED  //this will be not 0 when ghost points are used
	rhs_[OFF(lefti_,j,0)] = 0.0; //no adaption of i-index needed here
	rhs_[OFF(lefti_,j,1)] = 0.0;
	rhs_[OFF(lefti_,j,2)] = 0.0;
	rhs_[OFF(lefti_,j,3)] = 0.0;
#else
#ifdef ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS
	rhs_[OFF(lefti_,j,0)] = 0.0; //no adaption of i-index needed here
	rhs_[OFF(lefti_,j,1)] = 0.0;
	rhs_[OFF(lefti_,j,2)] = 0.0;
	rhs_[OFF(lefti_,j,3)] = 0.0;
#endif
#endif
	//!=============================================================
	  
	for(int k = 0; k < nchem; ++k){
	  w[OFF(lefti_,j,4+k)] = Y_in_[k]*FlowVariables<VType>::rho(lefti_,j,w,nx_);
#ifndef GHOST_POINTS_INCLUDED
	  rhs_[OFF(lefti_,j,4+k)] = 0.0; //no adaption of i-index needed here
#else
#ifdef ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS
	  rhs_[OFF(lefti_,j,4+k)] = 0.0; //no adaption of i-index needed here  
#endif
#endif
	}
      } //end for 



      for(int i = downbeg_; i < downend_; ++i){	
	//! DOWN · DOWN · DOWN · DOWN · DOWN · DOWN · DOWN: wall
	//#if defined 

	//update species at boundaries first, since w is used a lil' bit later
	//in conjunction with computing rho at the boundary
	for(int k = 0; k < nchem; ++k){
#ifndef GHOST_POINTS_INCLUDED
	  YBdy_[k] = BdyFctType::down(*this,i,w,4+k);
#else
	  YBdy_[k] = w[OFF(i,downj_,4+k)];
#endif
	  
#ifndef GHOST_POINTS_INCLUDED
	  rhs_[OFF(i,downj_,4+k)] = 0.0; //no adaption of j-index needed here
#endif
	}

	
	//! adjust Yfrac at boundary again (resp. ghost points)
	CorrectYfrac_.smart_balance_to_unity(YBdy_,N2,approxEqual_);

	WbarBdy_ = Wbar(YBdy_.begin());

	
#ifndef GHOST_POINTS_INCLUDED 
	//Neumann-BC on pressure, cf. [ANDERSON, p. 458, Eq. (10.17)]
	pBdy_ = BdyFctType::p_bdy_down(*this,i,w);
#else //ghost points used
        pBdy_ = (w[OFF(i,downj_,0)]*PhysicalConstants<DType>::Rgas*walltempprofile_[i-offs_])/WbarBdy_; //use EOL for pressure, p = (rho·R·T)/Wbar
#endif

#ifdef IGNITION_AT_HOT_WALLS
#ifndef GHOST_POINTS_INCLUDED
	rhoBdy_ = (pBdy_*WbarBdy_)/(PhysicalConstants<DType>::Rgas*walltempprofile_[i-offs_]);
#else
	rhoBdy_ = w[OFF(i,downj_,0)];
#endif
	TBdy_ = walltempprofile_[i-offs_];
       
#ifndef NDEBUG  
	if(!is_well_defined_value(rhoBdy_))
	  ADONIS_ERROR(ValueError,"BDY ERR: p = "<< pressure(i,downj_,w) << " wbar = "<< wbar(i,downj_,w) << ".");
#endif

#else //no ignition at hot walls
	rhoBdy_ = (pBdy_*WbarBdy_)/(PhysicalConstants<DType>::Rgas*T_wall_);
	TBdy_ = T_wall_;
#endif

	//! applies to each case (ignition at bdy or not)
#ifndef GHOST_POINTS_INCLUDED
	w[OFF(i,downj_,0)] = rhoBdy_;
#endif
	w[OFF(i,downj_,1)] = v1_wall_*FlowVariables<VType>::rho(i,downj_,w,nx_); //no-slip conditions for velocity components
	w[OFF(i,downj_,2)] = v2_wall_*FlowVariables<VType>::rho(i,downj_,w,nx_);
	w[OFF(i,downj_,3)] = TBdy_*FlowVariables<VType>::rho(i,downj_,w,nx_);
	for(int k = 0; k < nchem; ++k){
	  w[OFF(i,downj_,4+k)] = YBdy_[k]*FlowVariables<VType>::rho(i,downj_,w,nx_);
	}

	
#ifndef GHOST_POINTS_INCLUDED
	//!\f$\partial_t u|_{\partial \Omega} = 0\f$
	rhs_[OFF(i,downj_,0)] = 0.0; //no adaption of j-index needed here
	rhs_[OFF(i,downj_,1)] = 0.0;
	rhs_[OFF(i,downj_,2)] = 0.0;
	rhs_[OFF(i,downj_,3)] = 0.0;
#else //use ghost points
#ifdef ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS
	//!rho will be updated by its PDE
	//rhs_[OFF(i,downj_,0)] = 0.0; //is calculated by EOS on bdy

	rhs_[OFF(i,downj_,1)] = 0.0;  //explicitly given at bdy by const value
	rhs_[OFF(i,downj_,2)] = 0.0;  //explicitly given at bdy by const value
	rhs_[OFF(i,downj_,3)] = 0.0;  //explicitly given at bdy by const value
#endif
#endif
	//!=============================================================
	

	//!UP · UP · UP  · UP · UP · UP · UP · UP · UP: wall -- not necessarily 
	//!a symmetric scenario

	//update species at boundaries first, since w is used a lil' bit later
	//in conjunction with computing rho at the boundary
	for(int k = 0; k < nchem; ++k){
#ifndef GHOST_POINTS_INCLUDED
	  YBdy_[k] = BdyFctType::up(*this,i,w,4+k); //second_order_ZERO_Neumann_dy_up(i,w,4+k);
#else
	  YBdy_[k] = w[OFF(i,upj_,4+k)]; 
#endif

#ifndef GHOST_POINTS_INCLUDED
	  rhs_[OFF(i,upj_,4+k)] = 0.0; //no adaption of i-index needed here
#endif
	}

	//!preserve mass at bdy and try to correct if necessary
	//! adjust Yfrac at boundary again
	CorrectYfrac_.smart_balance_to_unity(YBdy_,N2,approxEqual_);

	WbarBdy_ = Wbar(YBdy_.begin());

	        
#ifndef GHOST_POINTS_INCLUDED
	//Neumann-BC on pressure, cf. [ANDERSON, p. 458, Eq. (10.17)]
	pBdy_ =	BdyFctType::p_bdy_up(*this,i,w);
#else //use ghost points and calculate p via EOS using the bdy rho value which has been computed by the discretization of the PDE
	pBdy_ = (w[OFF(i,upj_,0)]*PhysicalConstants<DType>::Rgas*walltempprofile_[i-offs_])/WbarBdy_;
#endif

#ifdef IGNITION_AT_HOT_WALLS
#ifndef GHOST_POINTS_INCLUDED
	rhoBdy_ = (pBdy_*WbarBdy_)/(PhysicalConstants<DType>::Rgas*walltempprofile_[i-offs_]);
#else //ghost points used. Just value from PDE discretization
	rhoBdy_= w[OFF(i,upj_,0)];
#endif
	TBdy_ = walltempprofile_[i-offs_];

#ifndef NDEBUG  
	if(!is_well_defined_value(rhoBdy_)) 
	  ADONIS_ERROR(ValueError,"BDY ERR: p = "<<  pressure(i,upj_,w) << " wbar = "<< wbar(i,upj_,w) << ".");
#endif

#else //no ignition at hot walls
	rhoBdy_ = (pBdy_*WbarBdy_)/(PhysicalConstants<DType>::Rgas*T_wall_);
	TBdy_ = T_wall_;
#endif

//!applies to every physical quantity regardless ignition zone
#ifndef GHOST_POINTS_INCLUDED
	w[OFF(i,upj_,0)] = rhoBdy_;
#endif
	w[OFF(i,upj_,1)] = v1_wall_*FlowVariables<VType>::rho(i,upj_,w,nx_);   //no-slip condition
	w[OFF(i,upj_,2)] = v2_wall_*FlowVariables<VType>::rho(i,upj_,w,nx_);
	w[OFF(i,upj_,3)] = TBdy_*FlowVariables<VType>::rho(i,upj_,w,nx_);
	for(int k = 0; k < nchem; ++k)
	  w[OFF(i,upj_,4+k)] = YBdy_[k]*FlowVariables<VType>::rho(i,upj_,w,nx_);
	
#ifndef GHOST_POINTS_INCLUDED
	//!\f$\partial_t u|_{\partial \Omega} = 0\f$
	rhs_[OFF(i,upj_,0)] = 0.0; 
	rhs_[OFF(i,upj_,1)] = 0.0; 
	rhs_[OFF(i,upj_,2)] = 0.0;  
	rhs_[OFF(i,upj_,3)] = 0.0;
#else
#ifdef ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS
	//!rho will be updated?
	//rhs_[OFF(i,upj_,0)] = 0.0;  //rho is updated by EOS on bdy

	rhs_[OFF(i,upj_,1)] = 0.0;  //explicitly given at bdy by const value
	rhs_[OFF(i,upj_,2)] = 0.0;  //explicitly given at bdy by const value
	rhs_[OFF(i,upj_,3)] = 0.0;  //explicitly given at bdy by const value
#endif
#endif
      }


      //! RIGHT: outlet
      for(int j = rightbeg_; j <= rightend_; ++j){ //TODO //without lower and upper edge!
	for(int k = 0; k < nchem; ++k){
	  YBdy_[k] = BdyFctType::right(*this,j,w,4+k);
	}

	//preserve mass at boundary and try to adjust
	//! adjust Yfrac at boundary again
	CorrectYfrac_.smart_balance_to_unity(YBdy_,N2,approxEqual_);

	WbarBdy_ = Wbar(YBdy_.begin());
	TBdy_ = BdyFctType::right(*this,j,w,3);
	
#ifndef GHOST_POINTS_INCLUDED
	pBdy_ = BdyFctType::p_bdy_right(*this,j,w);
       	rhoBdy_ = (pBdy_*WbarBdy_)/(PhysicalConstants<DType>::Rgas*TBdy_);
#else
	pBdy_= (w[OFF(righti_,j,0)]*PhysicalConstants<DType>::Rgas*TBdy_)/WbarBdy_;
#endif



#ifndef GHOST_POINTS_INCLUDED
	w[OFF(righti_,j,0)] = rhoBdy_;
#endif
	w[OFF(righti_,j,1)] = BdyFctType::right(*this,j,w,1)*FlowVariables<VType>::rho(righti_,j,w,nx_); 
	w[OFF(righti_,j,2)] = BdyFctType::right(*this,j,w,2)*FlowVariables<VType>::rho(righti_,j,w,nx_); 
	w[OFF(righti_,j,3)] = TBdy_*FlowVariables<VType>::rho(righti_,j,w,nx_);
	
	//repair_out_of_bound_quantity(nx_-1,j,w,3,280.,3500., 300.,true);

#ifndef GHOST_POINTS_INCLUDED  //need not to use 'righti_' here coz this is only relevant when no ghost points are used
	//!\f$\partial_t u|_{\partial \Omega} = 0\f$
	rhs_[OFF(righti_,j,0)] = 0.0;
	rhs_[OFF(righti_,j,1)] = 0.0;
	rhs_[OFF(righti_,j,2)] = 0.0;
	rhs_[OFF(righti_,j,3)] = 0.0;
	//! if ghost points are used, all species are calculated at bdy by their differential equations
#else // use ghost points
	//rhs_[OFF(righti_,j,0)] = 0.0; //updated on bdy by EOS
#endif
	for(int k = 0; k < nchem; ++k){
	  w[OFF(righti_,j,4+k)] = YBdy_[k]*FlowVariables<VType>::rho(righti_,j,w,nx_); //second_order_ZERO_Neumann_dx_right(j,w,4+k);

#ifndef GHOST_POINTS_INCLUDED
	  rhs_[OFF(righti_,j,4+k)] = 0.0;
#endif
	}
      }
    } //end set_boundary


//!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!+++++++++++++++++++++++++++ GHOST POINTS ++++++++++++++++++++++++++++++++
//!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 template<class X>
 void set_ghost_points(X& w){
   //!only do something here when ghost points are used
#ifdef GHOST_POINTS_INCLUDED  
   for(int j = 1; j < ny_-1; ++j){
     for(int l = 0; l < nprim; ++l){  //loop over ALL phys. quantities
#ifdef REFLECTIVE_BDY_CONDITIONS
       w[OFF(0,j,l)] = w[OFF(2,j,l)];               //left ghost points
       w[OFF(nx_-1,j,l)] = w[OFF(nx_-3,j,l)];       //right ghost points
#else //copy from boundary
       w[OFF(0,j,l)] = w[OFF(1,j,l)];               //left ghost points
       w[OFF(nx_-1,j,l)] =  w[OFF(nx_-2,j,l)];      //right ghost points 
#endif
       //rhs_ ghost points values are meaninglesss
       rhs_[OFF(0,j,l)] = 0.0;
       rhs_[OFF(nx_-1,j,l)] = 0.0;
     }
   } //end left and right

   for(int i = 0; i < nx_; ++i){
      for(int l = 0; l < nprim; ++l){ //loop over ALL phys. quantities
#ifdef REFLECTIVE_BDY_CONDITIONS
	w[OFF(i,0,l)] = w[OFF(i,2,l)];              //lower ghost points
	w[OFF(i,ny_-1,l)] = w[OFF(i,ny_-3,l)];      //upper ghost points
#else //copy from boundary
	w[OFF(i,0,l)] = w[OFF(i,1,l)];              //lower ghost points
	w[OFF(i,ny_-1,l)] = w[OFF(i,ny_-2,l)];      //upper ghost points
#endif
	//rhs_ ghost points values are meaninglesss
	rhs_[OFF(i,0,l)] = 0.0;
	rhs_[OFF(i,ny_-1,l)] = 0.0;
      }
   } //end down and up

   //======TODO?: special treatment of corners? 
#endif //GHOST_POINTS_INCLUDED
 }//end set_ghost_points
    
    //! ======================================================================
    //! ===== OPERATOR · OPERATOR · OPERATOR · OPERATOR · OPERATOR · =========
    //! ======================================================================
    //! Implementation of \f$ rhs_ = L_h[U], \f$ where \f$rhs_= \partial_t u\f$ 
    //! \param vars INPUT (will be copied locally to u_)
    //! \return rhs_ OUTPUT (the dynamics)
    template<class VEC>
    VType& operator()(const VEC& vars){ 
      //HARD COPY 
      u_ = vars;

      //! BOUNDARY
      (*this).set_boundary(u_);
      //!GHOST POINTS (if any are used at all)
      (*this).set_ghost_points(u_);
     
      //!INTERIOR · INTERIOR· INTERIOR · INTERIOR· INTERIOR· INTERIOR· INTERIOR
      //!INTERIOR -- use second order approximations for first, sec. and 
      //            mixed derivatives
      //!INTERIOR · INTERIOR· INTERIOR · INTERIOR· INTERIOR· INTERIOR· INTERIOR
#ifdef USE_OPENMP
      //adonis_assert((NUMTHREADS > 0));
#pragma omp parallel for num_threads(NUMTHREADS)
#endif 
      for(int i = 1; i < nx_-1; ++i){ 
	for(int j = 1; j < ny_-1; ++j){
	   

#ifdef REPAIR_ACTION
	  //! assure that values of primitive variables stay within their
	  //! respective physical bounds
	  repair_prim_vars_in_9_points(i,j,u_);

	  //! preserve mass Y_nspec = 1 - \sum_k^{nspec-1} Y_k
	  preserve_mass(N2,i,j,u_);   //N2 is excess species (bath gas)
#endif


	  cp_ = cp(i,j,u_); //Cp_ has also been filled

#ifndef NDEBUG
	  adonis_assert((cp_ > T()) && is_well_defined_value(cp_));
	  adonis_assert(u_[OFF(i,j,0)] > T());  //rho must be positive and bounded away from zero
#endif


#ifndef LOW_MACH_NUMBER_FORMULATION //if not defined then use usual computation	
	    //!\f$\rho\f$  
	  rhs_[OFF(i,j,0)] = 
#ifdef RHO_ALGEBRAIC
#ifdef PRESSURE_CONST
	    (p0_*wbar(i,j,u_))/(PhysicalConstants<T>::Rgas*temperature(i,j,u_));
#endif
#else //no algebraic evaluation
	  -(convection(i,j,0,u_));

#endif
#endif
	  adonis_assert(is_well_defined_value(rhs_[OFF(i,j,0)]));

	  //!\f$v_1\f$
	  rhs_[OFF(i,j,1)] =
#ifdef V1_CONSTANT
	    0.0;
#else
	  -( convection(i,j,1,u_)
	     - gov_fac(i,j,1,u_,cp_)*
	     (4./3.*( (0.5*(mu(i+1,j,u_)+mu(i,j,u_))*(v1(i+1,j,u_)-v1(i,j,u_)) - 0.5*(mu(i-1,j,u_)+mu(i,j,u_))*(v1(i,j,u_)-v1(i-1,j,u_)))/ntimes<2>(hx_)) - 2./3.*(mu(i+1,j,u_)*(v2(i+1,j+1,u_) - v2(i+1,j-1,u_)) - mu(i-1,j,u_)*(v2(i-1,j+1,u_)-v2(i-1,j-1,u_)))/(4.*hx_*hy_)) - 
	     gov_fac(i,j,1,u_,cp_)*
	     ( (0.5*(mu(i,j+1,u_)+mu(i,j,u_))*(v1(i,j+1,u_) - v1(i,j,u_)) - 0.5*(mu(i,j-1,u_)+mu(i,j,u_))*(v1(i,j,u_) - v1(i,j-1,u_)))/ntimes<2>(hy_) + (mu(i,j+1,u_)*(v2(i+1,j+1,u_) - v2(i-1,j+1,u_)) - mu(i,j-1,u_)*(v2(i+1,j-1,u_) - v2(i-1,j-1,u_)))/(4.*hx_*hy_)) 														      - gov_fac(i,j,1,u_,cp_)*
u_[OFF(i,j,0)]*gravitational_force_x(i,j,u_)
																																												  + gov_fac(i,j,1,u_,cp_)*mechanical_compression_x(i,j,u_)
	     ); 
#endif
	   
	  adonis_assert(is_well_defined_value(rhs_[OFF(i,j,1)]));

	    //! \f$v_2\f$  
	    rhs_[OFF(i,j,2)] = 
#ifdef V2_CONSTANT
	      0.0;   //constant velocity field const means zero variations
#else
	    -( convection(i,j,2,u_) 
	       - gov_fac(i,j,2,u_,cp_)*
	       ( (mu(i+1,j,u_)*(v1(i+1,j+1,u_) - v1(i+1,j-1,u_)) - mu(i-1,j,u_)*(v1(i-1,j+1,u_)-v1(i-1,j-1,u_)))/(4.*hx_*hy_) + (0.5*(mu(i+1,j,u_)+mu(i,j,u_))*(v2(i+1,j,u_)-v2(i,j,u_)) - 0.5*(mu(i-1,j,u_)+mu(i,j,u_))*(v2(i,j,u_)-v2(i-1,j,u_)))/ntimes<2>(hx_) ) -  
	       gov_fac(i,j,2,u_,cp_)*
(4./3.*((0.5*(mu(i,j+1,u_)+mu(i,j,u_))*(v2(i,j+1,u_) - v2(i,j,u_)) - 0.5*(mu(i,j-1,u_)+mu(i,j,u_))*(v2(i,j,u_) - v2(i,j-1,u_)))/ntimes<2>(hy_)) - 2./3.*((mu(i,j+1,u_)*(v1(i+1,j+1,u_) - v1(i-1,j+1,u_)) - mu(i,j-1,u_)*(v1(i+1,j-1,u_) - v1(i-1,j-1,u_)))/(4.*hx_*hy_))) 
	       - gov_fac(i,j,2,u_,cp_)*
u_[OFF(i,j,0)]*gravitational_force_y(i,j,u_)
	       + gov_fac(i,j,2,u_,cp_)*
mechanical_compression_y(i,j,u_)
	       ); 
#endif
			
	    adonis_assert(is_well_defined_value(rhs_[OFF(i,j,2)]));
	 
	    //!TEMPERATURE \f$ T \f$
	    
#ifndef NDEBUG  
	    is_temperature_ok(i,j,0,u_);  // 0 = take bounds of first chem. spec
#endif	   

	  
	    
#ifndef NDEBUG
	    //! prevent division by zero or production of NaNs
	    adonis_assert(!is_zero(cp_));
	    adonis_assert(!is_zero(u_[OFF(i,j,0)]));
#endif

	    //!CHEMISTRY -- calculate \f$\dot{\omega}\f$
	    chemistry(dotomega_,i,j,u_); //chemical reactions

	    //! note: for low-speed flows, mechanical compression  \f$\frac{Dp}{Dt}\f$ as well viscous dissipation can be neglected [KEE, p. 115 (pdf: 142), bottom], which is realized here
	    tempCond_ = heat_conduction(i,j,u_);
	    tempDiff_ = temperature_diffu_term(i,j,u_);
	    tempHeat_ = heat_production(dotomega_,i,j,u_);

#ifdef VALUE_PLAUSIBILITY_CHECK
	    
	    check_for_good_values(tempCond_,i,j,-42,"bad term: 'heat_conduction'.");
	    check_for_good_values(tempDiff_,i,j,-42,"bad term: 'temp_diffu_term'.");
	    check_for_good_values(tempHeat_,i,j,-42,"bad term: 'heat_production'");     
#endif
	    rhs_[OFF(i,j,3)] = -( convection(i,j,3,u_)
				  -gov_fac(i,j,3,u_,cp_)*energy_mechanical_expression(i,j,u_)
				  - gov_fac(i,j,3,u_,cp_)*tempCond_
				 +               //PLUS/MINUS ??
				 gov_fac(i,j,3,u_,cp_)*tempDiff_
				 +                 //PLUS ??
				 gov_fac(i,j,3,u_,cp_)*tempHeat_ 
				  //dissipation
				  -
				  gov_fac(i,j,3,u_,cp_)*dissipation(i,j,u_)
				  
				  ); 

	    //std::cout  << "d_t T("<<i<<","<<j<<") ="<< rhs_[OFF(i,j,3)] <<std::endl;

	    if(!is_well_defined_value(rhs_[OFF(i,j,3)])){
	      ADONIS_ERROR(ValueError, "BAD value: rhs_("<<i<<","<<j<<",3) = "<< rhs_[OFF(i,j,3)] << "\n DETAILED: cp = "<< cp_ << "  tempCond_ = "<< tempCond_ <<" \n tempDiff_ = "<< tempDiff_ << "  tempHeat_ = "<< tempHeat_ << "  dissipation = "<< dissipation(i,j,u_)<<".");
	    }
	    // adonis_assert(is_well_defined_value(rhs_[OFF(i,j,3)]));


	    //! SPECIES \f$ Y_k, \qquad k = 1, \ldots, nchem\f$
	     species_diffu_term(diffspec_,i,j,u_); //now diffspec_ is filled
	    // diffspec_= 1.e-05;
	    //species_diffu_term_wbar(diffspec_,i,j,u_); //now diffspec_ filled
	    for(int k = 0; k < nchem; ++k){
#ifdef VALUE_PLAUSIBILITY_CHECK
	      check_for_good_values(diffspec_[k],i,j,-42,"bad term:  'diffspec_[k]'.");
	      check_for_good_values(dotomega_[k],i,j,-42,"bad term: 'dotomega__[]'.");
#endif

	      //species evolution
	      sp_ = -( convection(i,j,4+k,u_)
                                      +              //PLUS/MINUS???
                                      gov_fac(i,j,4+k,u_,cp_)*diffspec_[k]
                                      -
				      gov_fac(i,j,4+k,u_,cp_)*dotomega_[k]*DataType::molar_masses()[AccessType::spec_access(k)] );

	      rhs_[OFF(i,j,4+k)] =
#ifdef SOLVE_FOR_K_MINUS_1_CONSERVATION_EQS
		//! only solve for K-1 species. The excess species will be reconstructed via \f$ Y_{N2} = 1. - \sum_{k\noteq N2} (Y_k) \f$ via repair(N2,y_n) in subsequent iterations
		( (k == N2) ? 0.0 : sp_ ); //if k = excess set rhs to 0
#else
	      sp_;
#endif

	      adonis_assert(is_well_defined_value(rhs_[OFF(i,j,4+k)]));
	    }
	    
#ifdef LOW_MACH_NUMBER_FORMULATION
	    rhs_[OFF(i,j,0)] =
#ifdef RHO_ALGEBRAIC
#ifdef PRESSURE_CONST
	    (p0_*wbar(i,j,u_))/(PhysicalConstants<T>::Rgas*temperature(i,j,u_));
#endif
#else //no algebraic evaluation
	     rho_low_mach(i,j,u_);
	    
#endif
#endif //end LOW_MACH_NUMBER_FORMULATION


#ifdef VALUE_PLAUSIBILITY_CHECK
	    for(int l= 0; l < nprim; ++l){
	      check_for_good_values(rhs_[OFF(i,j,l)],i,j,l);
	    }
#endif

	} //END INTERIOR y
      }  //END INTERIOR x

#ifdef GHOST_POINTS_INCLUDED
#ifdef ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS
//! calculate bdy conditions again to make some rhs_ zero where e.g. explicit
//! values are known (left for all quantities, up and down for v1,v2 and T)
    (*this).set_boundary(u_);
    (*this).set_ghost_points(u_); //make rhs_ zero at ghost points 
#endif
#endif
      return rhs_;
    }

    
    const ArrayType& mixture_in() const{return Y_in_;}
    const BaseType& v1_in() const {return v1_in_;}
    const BaseType& v2_in() const {return v2_in_;}
    const BaseType& T_in() {return T_in_;}
    const T& rho_in() const{
      return rhoIn_;
    }

   //! get previous composition from the numerical scheme
    //! used for calculating the mechanical energy expression
    void get_y_prev(IteratorTType y_prevIt, const T& stepsize){
      Veciter_ = y_prevIt;
      stepsize_ = stepsize;
      isPrevPressInvoked_=true;  //guarantees that it is invoked somewhere
    }

   

  private:
    VType rhs_,
      veloprofile_,
      walltempprofile_;
    BuilderType BCS_;

    BaseType a_,b_,c_,d_,hx_,hy_, t0_, tend_, diam_, halfdiam_,gravi1_,gravi2_, near_zero_, smallVal_, rhoMin_;
    T  wbarIn_, rhoIn_;
    int nx_, ny_, npt_, totpoints_;   
   
    BaseType p0_, v1_in_,v2_in_,v1_wall_,v2_wall_,T_wall_,T_in_, x_ignition_, //standard pressure
      v1_max_,v2_max_, tooLarge_;

    ArrayType X_, X0_, X1_,X2_,X3_,diffspec_, C_, dotomega_, Cp_,
      Y_in_;
    VType v1bdy_,
      u_;
    T cp_, tempconc_, tempDiff_, specDiff_, tempHeat_, specHeat_, tempCond_, sp_; //for storage
   
    T fivepointstencil_[5];
    T wbarfive_[5];


    T jflux_[nchem];
    T meanVal_;
    T pressure_;
    T useMe_;
    IteratorTType Veciter_;
    bool isBeginning_;
    T stepsize_;
    bool isPrevPressInvoked_;
   
    BoundaryWallTemperatureType bdytemp_;
    ViscosityType Mav_;
    ConductivityType Mac_;
    DiffusionType Mad_,Mad0_,Mad1_,Mad2_,Mad3_;
    
    RatioType Theta_, Theta0_, Theta1_, Theta2_, Theta3_;
    
    BaseType yhottest_, Tignition_, length_;

    bool useGravi_,
      useMechCompression_;

    unsigned int pressFac_;
    BaseType approxEqual_;

    T WbarBdy_, pBdy_, rhoBdy_,TBdy_;
    T td_; //temperature diffusion term
    T meanmolmass_; 
    T tmpvar_; //temporary variable
    ArrayType YBdy_;

    int leftbeg_, leftend_, lefti_, downbeg_, downend_, downj_, upj_, rightbeg_, rightend_, righti_, offs_;

    CorrectFractions<ArrayType> CorrectYfrac_;
    
    FancyMessages FM_;

    //friend classes have access to all elements of a host class 
    friend class BoundaryFunction2D<FDMSettings::orderOfBdyCond,ThisType>;
    

    //PRIVATE FUNCTIONS -- don't use them from outside

    int OFF(int i, int j, int spec) const{
      if((i < 0) || (i >= nx_) || (j < 0) || (j >= ny_) || (spec < 0) || (spec >= 4+nchem))
	ADONIS_ERROR(IndexError,"Index violation by trying to access i = "<<i<<", j = "<<j<<", phys. quant = "<< spec << ".");
      return (i + nx_*j + spec*npt_);
    }
    

   T primitive(int i, int j, const VType& u, int quant){
#ifdef NONCONSERVATIVE_FORM
     return u[OFF(i,j,quant)];
#else //retrieve primitive variables from conservative vars
     return ( (quant == 0) ? (u[OFF(i,j,0)]) : (u[OFF(i,j,quant)]/(u[OFF(i,j,0)])) );
#endif
   }


    void assign_flow_variable(int i, int j, VType& u, int quant,const T& val){
      u[OFF(i,j,quant)] = 
#ifdef NONCONSERVATIVE_FORM
	val;
#else //rho*val (and val when quant is rho)
	( (quant == 0) ? (val) : (val*u[OFF(i,j,0)]) );
#endif
    }

    
    T rho(int i, int j, const VType& u){
      return u[OFF(i,j,0)]; 
    }
    

    //!PRESERVE MASS
    void preserve_mass(int i, int j, VType& u){ //u is changed accoringly
      T sm1 = T();
      
      //! compute \f$ Y_{n} = 1 - \sum_{k = 1}^{n-1} Y_k\f$
      for(int k = 0; k < nchem-1; ++k)
	sm1 += de_negativize(Yfrac(i,j,u,k));
      
      u[OFF(i,j,4+(nchem-1))] = 
#ifdef NONCONSERVATIVE_FORM
	1. - sm1;
#else
      (1. - sm1)*u[OFF(i,j,0)];
#endif
    }


    void preserve_mass(int select, int i, int j, VType& u){
      T sm1 = T();
      for(int k = 0; k < DataType::nspec; ++k){
	if(k != select){
	  sm1 += de_negativize(Yfrac(i,j,u,k));
	}
      }
      u[OFF(i,j,4+select)] = 
#ifdef NONCONSERVATIVE_FORM
	1. - sm1;
#else
      (1. - sm1)*u[OFF(i,j,0)];
#endif

    }
    

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
	    
	  
	  if(temperature(IX,JX,u) < DataType::temperature_bounds()[3*k] || 
	     temperature(IX,JX,u) > DataType::temperature_bounds()[3*k+1]){
	   
	    ADONIS_ERROR(BoundsError,"Temperature bounds violated: T_{i"<<oper1<<l << ",j"<<oper2<<t << "} = "<< u[OFF(IX,JX,3)]<<", i.e. T("<<IX<<","<<JX<<") = "<< temperature(IX,JX,u) << ".");
	  }
	}
      }
    }



    //CONVECTION -- usef for conservative and nonconservative form 
    //! arguments: v_l velocity component
    //! convective terms are approximated by first-order upwind formulas
    //! if the velocity component is negative
    T upwind_downwind_x(int i, int j, int t, const T& v_l, const VType& u){
#ifndef ORDER_2_CONVECTIVE_TERM
      if(v_l <= 0) //first order downwind
	return (u[OFF(i+1,j,t)] - u[OFF(i,j,t)])/hx_;
      else        //first order upwind
	return (u[OFF(i,j,t)] - u[OFF(i-1,j,t)])/hx_;
#else
      //you may use 2nd order convective terms,cf. CHEMKIN-PRO Theory, pdf. 207
      return ((u[OFF(i+1,j,t)] - u[OFF(i-1,j,t)])/(2*hx_) );
#endif
    }

     T upwind_downwind_y(int i, int j, int t, const T& v_l, const VType& u){
#ifndef ORDER_2_CONVECTIVE_TERM  
      if(v_l <= 0) 
	return (u[OFF(i,j+1,t)] - u[OFF(i,j,t)])/hy_;
      else 
	return (u[OFF(i,j,t)] - u[OFF(i,j-1,t)])/hy_;
#else
      //you may use 2nd order convective terms,cf. CHEMKIN-PRO Theory, pdf. 207
      return (u[OFF(i,j+1,t)] - u[OFF(i,j-1,t)])/(2*hy_);
#endif
    }
    

    //!general convection of physical quantity t
    T convection(int i, int j, int t, const VType& u){
      adonis_assert((t >= 0) && (t < 4+DataType::nspec));
#ifdef NONCONSERVATIVE_FORM
      //u contains already primitive variables rho, v1, v2, T, Y1,..., Y_K
      return ( (t==0) ? //non-conservative form, [KEE, pdf. 780], [ANDERSON, p. 58, eq. (2.35)]
 //see also [ANDERSON, Computational Fluid Dynamics, Eq. (6.9), p. 219]	
	       //              dx_rho                                           //dy_rho
	       ( ( v1(i,j,u)*upwind_downwind_x(i,j,0,v1(i,j,u),u) + v2(i,j,u)*upwind_downwind_y(i,j,0,v2(i,j,u),u)
		   //                  //dx_v1                                dy_v2 
		   + u[OFF(i,j,0)]*(upwind_downwind_x(i,j,1,v1(i,j,u),u) + upwind_downwind_y(i,j,2,v2(i,j,u),u)) ) ) : 
	       //all phys quantities which are not rho
	       ( v1(i,j,u)*upwind_downwind_x(i,j,t,v1(i,j,u),u) + v2(i,j,u)*upwind_downwind_y(i,j,t,v2(i,j,u),u) ) );
#else //!conservative form. u contains flow variables rho, rho·v1, rho·v2, rho·T, rho·Y1,...,rho·Y_K 
#ifndef ORDER_2_CONVECTIVE_TERM 
      return ( ((v1(i,j,u) <= 0.0)  ? ( (u[OFF(i+1,j,t)]*v1(i+1,j,u) - u[OFF(i,j,t)]*v1(i,j,u))/hx_ ) : ( (u[OFF(i,j,t)]*v1(i,j,u) - u[OFF(i-1,j,t)]*v1(i-1,j,u))/hx_ ))  +
	       ( (v2(i,j,u) <= 0.0) ? ( (u[OFF(i,j+1,t)]*v2(i,j+1,u) - u[OFF(i,j,t)]*v2(i,j,u))/hy_ ) : ( (u[OFF(i,j,t)]*v2(i,j,u) - u[OFF(i,j-1,t)]*v2(i,j-1,u))/hy_ ) ) ); 
#else //2nd order
      return ( (u[OFF(i+j,j,t)]*v1(i+1,j,u) - u[OFF(i-1,j,t)]*v1(i-1,j,u))/(2*hx_) + (u[OFF(i,j+1,t)]*v2(i,j+1,u) - u[OFF(i,j-1,t)]*v2(i,j-1,u))/(2*hy_) );
#endif //2nd order convetive term
#endif
    }


    //!prefactor in front of the terms in the governing equations' terms
    T gov_fac(int i, int j, int quant, const VType& u, const T& cp){
#ifdef NONCONSERVATIVE_FORM
      return ( ((quant == 1) || (quant == 2) || (quant >= 4)) ? (1./u[OFF(i,j,0)]) : ( (quant == 3) ? (1./(cp*u[OFF(i,j,0)])) : 1) );
#else //conservative form
      return ( (quant == 3) ? (1./cp) : 1 );
#endif
    }


    //! viscous dissipation functional. It must be always positive; irreversible work must increase thermal energy in the flow [KEE, pdf. 140] 
    T& dissipation(int i, int j, const VType& u){
     
      useMe_ = 0; //reset if it is still occupied by some other member

      useMe_ = mu(i,j,u)*( //4/3dx(v1 dx v1)
			  4./3.*( 0.5*(v1(i+1,j,u)+v1(i,j,u))*(v1(i+1,j,u)-v1(i,j,u)) - 0.5*(v1(i-1,j,u)+v1(i,j,u))*(v1(i,j,u)-v1(i-1,j,u)) )/ntimes<2>(hx_) 
      //-2/3dy(v2 dx v1)
			  -2./3.*(v2(i,j+1,u)*(v1(i+1,j+1,u)-v1(i-1,j+1,u)) - v2(i,j-1,u)*(v1(i+1,j-1,u)-v1(i-1,j-1,u)) )/(4*hx_*hy_)
	      
			  //+dy(v1 dx v2)
			  + ( v1(i,j+1,u)*(v2(i+1,j+1,u)-v2(i-1,j+1,u)) - v1(i,j-1,u)*(v2(i+1,j-1,u)-v2(i-1,j-1,u)) )/(4*hx_*hy_)
		     
			  //+dx(v2 dx v2)
			  + ( 0.5*(v2(i+1,j,u)+v2(i,j,u))*(v2(i+1,j,u)-v2(i,j,u)) - 0.5*(v2(i-1,j,u)+v2(i,j,u))*(v2(i,j,u)-v2(i-1,j,u)) )/ntimes<2>(hx_)
			  
			  
			  //+dy(v1 dy v1)
			  + (0.5*(v1(i,j+1,u)+v1(i,j,u))*(v1(i,j+1,u)-v1(i,j,u)) - 0.5*(v1(i,j-1,u)+v1(i,j,u))*(v1(i,j,u)-v1(i,j-1,u)) )/ntimes<2>(hy_)

			  
			  //+dx(v2 dy v1)
			  +( v2(i+1,j,u)*(v1(i+1,j+1,u)-v1(i+1,j-1,u)) - v2(i-1,j,u)*(v1(i-1,j+1,u)-v1(i-1,j-1,u)) )/(4*hx_*hy_)

			  //+4/3 dy(v2 dy v2)
			  +4./3.*(0.5*(v2(i,j+1,u)+v2(i,j,u))*(v2(i,j+1,u)-v2(i,j,u)) - 0.5*(v2(i,j-1,u)+v2(i,j,u))*(v2(i,j,u)-v2(i,j-1,u)) )/ntimes<2>(hy_)
			  //-2/3 dx(v1 dy v2)
			  -2./3.*(v1(i+1,j,u)*(v2(i+1,j+1,u)-v2(i+1,j-1,u)) - v1(i-1,j,u)*(v2(i-1,j+1,u)-v2(i-1,j-1,u)))/(4*hx_*hy_)
			  
			   );

#ifdef SHOW_IT_2_ME
      std::cout << "DISSIPATION = "<< useMe_ << std::endl;
#endif

      // if(useMe_ < 0.0){
      // 	FancyMessages().nice_output("Dissipation Phi < 0",35);
      // }

      if(!is_well_defined_value(useMe_)){
	std::cout << "v1(i,j-1,u) = " <<v1(i,j-1,u) <<std::endl;
	std::cout << "v1(i+1,j-1,u) = " <<v1(i+1,j-1,u) <<std::endl;
	std::cout << "v1(i+1,j,u) = " <<v1(i+1,j,u) <<std::endl;
	std::cout << "v1(i+1,j+1,u) = " <<v1(i+1,j+1,u) <<std::endl;
	std::cout << "v1(i,j+1,u) = " <<v1(i,j+1,u) <<std::endl;
	std::cout << "v1(i-1,j+1,u) = " <<v1(i-1,j+1,u) <<std::endl;
	std::cout << "v1(i-1,j,u) = " <<v1(i-1,j,u) <<std::endl;
	std::cout << "v1(i-1,j-1,u) = " <<v1(i-1,j-1,u) <<std::endl;
	std::cout << "v1(i,j,u) = " <<v1(i,j,u) <<std::endl;

	std::cout << "v2(i,j-1,u) = " <<v2(i,j-1,u) <<std::endl;
	std::cout << "v2(i+1,j-1,u) = " <<v2(i+1,j-1,u) <<std::endl;
	std::cout << "v2(i+1,j,u) = " <<v2(i+1,j,u) <<std::endl;
	std::cout << "v2(i+1,j+1,u) = " <<v2(i+1,j+1,u) <<std::endl;
	std::cout << "v2(i,j+1,u) = " <<v2(i,j+1,u) <<std::endl;
	std::cout << "v2(i-1,j+1,u) = " <<v2(i-1,j+1,u) <<std::endl;
	std::cout << "v2(i-1,j,u) = " <<v2(i-1,j,u) <<std::endl;
	std::cout << "v2(i-1,j-1,u) = " <<v2(i-1,j-1,u) <<std::endl;
	std::cout << "v2(i,j,u) = " << v2(i,j,u) <<std::endl;


	ADONIS_ERROR(ValueError,"Viscous dissipation has a bad value: Phi("<<i<<","<<j<<") = "<<useMe_ << ".");
      }

      return useMe_;
    }


    //! CHEMISTRY
    //in the reduced case  only the corresponding fixed quantities
    //are extracted, e.g. molar masses.
    T wbar(int i, int j, const VType& u) const{
      T wb = T();
      for(int k = 0; k < nchem; ++k){
	wb += Yfrac(i,j,u,k)/DataType::molar_masses()[AccessType::spec_access(k)];
      }
      if(is_zero(wb)){
	wb = smallVal_;  //prevent division by zero
      }
      return 1./wb;
    }
    
    //! this can be used if you have calculated some "array" containing mass fractions
    template<class ITER>
    T Wbar(ITER it){
      T wb = T();
      for(int k = 0; k < nchem; ++k){
	wb += it[k]/DataType::molar_masses()[AccessType::spec_access(k)];
      }
      if(is_zero(wb)){
	wb = smallVal_;  //prevent division by zero
      }
      return 1./wb;
    }

    //!NOTE: perturbation of X is very important!!!
     T mole_fraction(int i, int j, int k, const T& Wbar, const VType& u) const{
       T mfrac = T();
       //!This should normally not occur but from the numerical point of view,
       //!it may definitely happen that Y_k is negative (e.g. of order 1.e-24).
       //!The problem with this is, however, that std::sqrt for real values 
       //!expects a <B>positive</B> argument; otherwise a NaN is produced
       if(Yfrac(i,j,u,k) < 0.0){
	 // ADONIS_INFO(Information,"Negative Y_"<<k<<" = "<< u[OFF(i,j,4+k)] << " generated.\n   Considered (but not set to) ZERO.");
	 // //u[OFF(i,j,4+k)] = 0.0;   //already done by repair_Y
	 mfrac = 0.0;   //assume then Y_k is zero
       }
       else{
	 mfrac = Yfrac(i,j,u,k)*Wbar/DataType::molar_masses()[AccessType::spec_access(k)];
       }
       //! only perturb when approx zero.
       return ((is_zero(mfrac,1.e-14)) ? perturbation(mfrac,1.e-16) : mfrac);
    }
    
    //! needed to evaluate diffusion and viscosity coefficients 
    ArrayType& mole_fraction(int i, int j, const T& Wbar,  const VType& u){
      for(int k = 0; k < nchem; ++k){
	X_[k] = mole_fraction(i,j,k,Wbar,u);
	perturbation(X_[k],1.e-16);
      }
      return X_;
    }
    
    //!overload mole fraction vector for <II>any</II> random access container
    template<class RAC>
    RAC& mole_frac(RAC& xfrac, int i, int j, const T& Wbar, const VType& u){
      for(int k = 0; k < nchem; ++k){
	xfrac[k] = mole_fraction(i,j,k,Wbar,u);
	perturbation(xfrac[k],1.e-16);
      }
      return xfrac;
    }

    T concenctration(int i, int j, int k, VType& u, T& rho){
     
#ifndef NDEBUG
      if((rho < 0.0) && (!is_zero(Abs(rho))))
	ADONIS_ERROR(ValueError,"rho("<< i << ","<< j << ") = "<< rho << " is effectively smaller than ZERO.");
      if((Yfrac(i,j,u,k) < 0.0) && (!is_zero(Abs(Yfrac(i,j,u,k)))))
	ADONIS_ERROR(ValueError,"Y_"<< k<< "(" << i << ","<< j << ") = "<< Yfrac(i,j,u,k) << " is effectively smaller than ZERO.");
#endif
      return rho*Yfrac(i,j,u,k)/DataType::molar_masses()[AccessType::spec_access(k)];
    }

    ArrayType& concentration(int i, int j, VType& u){
      for(int k = 0; k < nchem; ++k)
	C_[k] = concenctration(i,j,k,u,u[OFF(i,j,0)]);
      return C_;
    }

    //! overload concentration 
    template<class RAC>
    RAC& concentration(RAC& conc, int i, int j, VType& u){
      for(int k = 0; k < nchem; ++k){
	tempconc_ = concenctration(i,j,k,u,u[OFF(i,j,0)]);
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
   
     ArrayType& chemistry(ArrayType& domeg, int i, int j, VType& u){
      concentration(C_,i,j,u); //C_ is filled now
      for(int k = 0; k < nchem; ++k){ //nchem is either full dim or reddim
	//!=========TODO: in either cases: C_ must be the FULL composition
	//! compare with 
	//! '../ode/examples/automaticallygeneratedsourceterms/redh2c6_attempt.hh'
	//! only compute \f$\dot{\omega}_k \f$ when spec. \f$k\f$ participates 
	//! in any reaction; otherwise it is const and therefore  we may take
	//! \f$\dot{\omega}_k = 0.0\f$ 
	domeg[k] = ( (DataType::is_species_reactive()[k]==true) ? (BCS_.net_formation_rate(k,temperature(i,j,u),C_)) : 0.0 ); //in the reduced case C_ must be of full dimension!!!
      }
      return domeg;
    }


    
    //! \f$ \dot{\omega} \f$ is computed by now
    T heat_production(const ArrayType& dotomega,int i, int j, const VType & u){
      T h = T();
      for(int k = 0; k < nchem; ++k)
	h += BCS_.H_T(AccessType::spec_access(k),temperature(i,j,u))*dotomega[k];
      return h;
    }

     T cp(int i, int j, const VType& u){
      T cres = T();
      for(int k = 0; k < nchem; ++k){
	//! fill Cp_ simultaneously
#ifndef NDEBUG
	
	if(temperature(i,j,u) < DataType::temperature_bounds()[3*AccessType::spec_access(k)] || temperature(i,j,u) > DataType::temperature_bounds()[3*AccessType::spec_access(k) + 1])
	  ADONIS_ERROR(BoundsError, "cp(i,j,u): Temperature out of bound T_{"<<i<<","<<j<< "} = "<< temperature(i,j,u) << ".");
#endif
	Cp_[k] = BCS_.C_p(AccessType::spec_access(k),temperature(i,j,u));
	cres += Cp_[k]/DataType::molar_masses()[AccessType::spec_access(k)]*Yfrac(i,j,u,k);
      }
      //isCpcalculated_ = true;
      return cres;
    }

    
    //!specific heat capacity at fixed pressure of spec k defined as
    //! \f$ c_{pk} = \frac{C_{pk}^0}{W_k}\f$
    T c_pk(int i, int j, const VType& u, int k){
#ifndef NDEBUG
      if(temperature(i,j,u) < DataType::temperature_bounds()[3*AccessType::spec_access(k)] || temperature(i,j,u) > DataType::temperature_bounds()[3*AccessType::spec_access(k) + 1])
	  ADONIS_ERROR(BoundsError, "cp(i,j,u): Temperature out of bound T_{"<<i<<","<<j<< "} = "<< temperature(i,j,u) << ".");
#endif
      return (BCS_.C_p(AccessType::spec_access(k),temperature(i,j,u))/DataType::molar_masses()[AccessType::spec_access(k)] );
    }


    //! at constant volume
    T cV(int i, int j, const VType& u){
       T cres = T();
      for(int k = 0; k < nchem; ++k){
#ifndef NDEBUG
	
	if(temperature(i,j,u) < DataType::temperature_bounds()[3*AccessType::spec_access(k)] || temperature(i,j,u) > DataType::temperature_bounds()[3*AccessType::spec_access(k) + 1])
	  ADONIS_ERROR(BoundsError, "cV(i,j,u): Temperature out of bound T_{"<<i<<","<<j<< "} = "<< temperature(i,j,u) << ".");
#endif
	cres += (BCS_.C_p(AccessType::spec_access(k),temperature(i,j,u)) - PhysicalConstants<T>::Rgas)/DataType::molar_masses()[AccessType::spec_access(k)]*Yfrac(i,j,u,k);
      }
      return cres;
    }

    //! EVALUATE TRANSPORT COEFFICIENTS AT point node \f$ (i,j) \f$
    //! compute \f$ \mu(T,X)\f$
    T mu(int i, int j, const VType& u){  
      return Mav_.compute_mixture_averaged_viscosity(temperature(i,j,u), mole_fraction(i,j,wbar(i,j,u),u));
    }

    //! \f$ \lambda(\rho, p, T, X) \f$
    T lambda(int i, int j, const VType& u){
      return Mac_.compute_mixture_averaged_conductivity(u[OFF(i,j,0)], pressure(i,j,u), temperature(i,j,u),mole_fraction(i,j,wbar(i,j,u),u));
    }

    T heat_conduction(int i, int j, const VType& u){
      return ( (0.5*(lambda(i+1,j,u)+lambda(i,j,u))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(lambda(i,j,u)+lambda(i-1,j,u))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) + (0.5*(lambda(i,j+1,u)+lambda(i,j,u))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(lambda(i,j,u)+lambda(i,j-1,u))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) );
    }


    

    //! Compute \f$ \Nabla \cdot  j_k \f$
    ArrayType& species_diffu_term(ArrayType& diffu, int i, int j, const VType& u){
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

      Mad0_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j-1,u),temperature(i,j-1,u),X0_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(pressure(i+1,j,u),temperature(i+1,j,u),X1_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j+1,u),temperature(i,j+1,u),X2_);
      
      Mad3_.compute_mixture_averaged_diffusion_coefficients(pressure(i-1,j,u),temperature(i-1,j,u),X3_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j,u),temperature(i,j,u),X_);

      
      //! thermal diffusion ratios
      Theta0_.compute_thermal_diffusion_ratios(temperature(i,j-1,u),X0_);
      Theta1_.compute_thermal_diffusion_ratios(temperature(i+1,j,u),X1_);
      Theta2_.compute_thermal_diffusion_ratios(temperature(i,j+1,u),X2_);
      Theta3_.compute_thermal_diffusion_ratios(temperature(i-1,j,u),X3_);
      Theta_.compute_thermal_diffusion_ratios(temperature(i,j,u),X_);

      tmpvar_ = 0; //Reset
      for(int k = 0; k < nchem; ++k){

	diffu[k] =  -( 
		       //+ dx(rho D_k^mix dy Y_k)
		       (0.5*(u[OFF(i+1,j,0)]*Mad1_[k] + u[OFF(i,j,0)]*Mad_[k])*(Yfrac(i+1,j,u,k)-Yfrac(i,j,u,k)) - 0.5*(u[OFF(i-1,j,0)]*Mad3_[k] + u[OFF(i,j,0)]*Mad_[k])*(Yfrac(i,j,u,k)-Yfrac(i-1,j,u,k)) )/ntimes<2>(hx_)  
		       // //+ dy(rho D_k^mix dy Y_k)
		       +  (0.5*(u[OFF(i,j+1,0)]*Mad2_[k] + u[OFF(i,j,0)]*Mad_[k])*(Yfrac(i,j+1,u,k)-Yfrac(i,j,u,k)) - 0.5*(u[OFF(i,j-1,0)]*Mad0_[k] + u[OFF(i,j,0)]*Mad_[k])*(Yfrac(i,j,u,k)-Yfrac(i,j-1,u,k)) )/ntimes<2>(hy_)
		           //+dx((rho D_k^mix Y_k)/Wbar dx Wbar)
		       + (0.5*((u[OFF(i+1,j,0)]*Mad1_[k]*Yfrac(i+1,j,u,k))/fivepointstencil_[1] + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/fivepointstencil_[4])*(fivepointstencil_[1]-fivepointstencil_[4]) - 0.5*((u[OFF(i-1,j,0)]*Mad3_[k]*Yfrac(i-1,j,u,k))/fivepointstencil_[3] + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/fivepointstencil_[4])*(fivepointstencil_[4]-fivepointstencil_[3]) )/ntimes<2>(hx_)
		       //+dy((rho D_k^mix Y_k)/Wbar dy Wbar)
		       + (0.5*((u[OFF(i,j+1,0)]*Mad2_[k]*Yfrac(i,j+1,u,k))/fivepointstencil_[2] + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/fivepointstencil_[4])*(fivepointstencil_[2]-fivepointstencil_[4]) - 0.5*((u[OFF(i,j-1,0)]*Mad0_[k]*Yfrac(i,j-1,u,k))/fivepointstencil_[0] + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/fivepointstencil_[4])*(fivepointstencil_[4]-fivepointstencil_[0]) )/ntimes<2>(hy_)
		       //+ dx((rho D_k^mix Y_k)/p dx p)
		       + (0.5*((u[OFF(i+1,j,0)]*Mad1_[k]*Yfrac(i+1,j,u,k))/pressure(i+1,j,u) + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/pressure(i,j,u))*(pressure(i+1,j,u)-pressure(i,j,u)) - 0.5*((u[OFF(i-1,j,0)]*Mad3_[k]*Yfrac(i-1,j,u,k))/pressure(i-1,j,u) + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/pressure(i,j,u))*(pressure(i,j,u)-pressure(i-1,j,u)) )/ntimes<2>(hx_)
		       //+ dy((rho D_k^mix Y_k)/p dy p)
		       + (0.5*((u[OFF(i,j+1,0)]*Mad2_[k]*Yfrac(i,j+1,u,k))/pressure(i,j+1,u) + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/pressure(i,j,u))*(pressure(i,j+1,u)-pressure(i,j,u)) - 0.5*((u[OFF(i,j-1,0)]*Mad0_[k]*Yfrac(i,j-1,u,k))/pressure(i,j-1,u) + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/pressure(i,j,u))*(pressure(i,j,u)-pressure(i,j-1,u)) )/ntimes<2>(hy_)
		       // - dx ((rho W_k D_k^mix Y_k)/(Wbar p) dx p)
		       - (0.5*((u[OFF(i+1,j,0)]*DataType::molar_masses()[k]*Mad1_[k]*Yfrac(i+1,j,u,k))/(fivepointstencil_[1]*pressure(i+1,j,u)) + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Yfrac(i,j,u,k))/(fivepointstencil_[4]*pressure(i,j,u)))*(pressure(i+1,j,u)-pressure(i,j,u))  - 0.5*((u[OFF(i-1,j,0)]*DataType::molar_masses()[k]*Mad3_[k]*Yfrac(i-1,j,u,k))/(fivepointstencil_[3]*pressure(i-1,j,u)) + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Yfrac(i,j,u,k))/(fivepointstencil_[4]*pressure(i,j,u)))*(pressure(i,j,u)-pressure(i-1,j,u)) )/ntimes<2>(hx_)
		        //- dy ((rho W_k D_k^mix Y_k)/(Wbar p) dy p)
		       - (0.5*((u[OFF(i,j+1,0)]*DataType::molar_masses()[k]*Mad2_[k]*Yfrac(i,j+1,u,k))/(fivepointstencil_[2]*pressure(i,j+1,u)) + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Yfrac(i,j,u,k))/(fivepointstencil_[4]*pressure(i,j,u)))*(pressure(i,j+1,u)-pressure(i,j,u))  - 0.5*((u[OFF(i,j-1,0)]*DataType::molar_masses()[k]*Mad0_[k]*Yfrac(i,j-1,u,k))/(fivepointstencil_[0]*pressure(i,j-1,u)) + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Yfrac(i,j,u,k))/(fivepointstencil_[4]*pressure(i,j,u)))*(pressure(i,j,u)-pressure(i,j-1,u)) )/ntimes<2>(hy_)
		       //+ dx ((rho Wk D_k^mix theta_k)/(Wbar T) dx T)
		       + (0.5*((u[OFF(i+1,j,0)]*DataType::molar_masses()[k]*Mad1_[k]*Theta1_[k])/(fivepointstencil_[1]*temperature(i+1,j,u)) + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Theta_[k])/(fivepointstencil_[4]*temperature(i,j,u)))*(temperature(i+1,j,u)-temperature(i,j,u))  - 0.5*((u[OFF(i-1,j,0)]*DataType::molar_masses()[k]*Mad3_[k]*Theta3_[k])/(fivepointstencil_[3]*temperature(i-1,j,u)) + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Theta_[k])/(fivepointstencil_[4]*temperature(i,j,u)))*(temperature(i,j,u)-temperature(i-1,j,u)) )/ntimes<2>(hx_)
		       //dy ((rho Wk D_k^mix theta_k)/(Wbar T) dy T)
		       + (0.5*((u[OFF(i,j+1,0)]*DataType::molar_masses()[k]*Mad2_[k]*Theta2_[k])/(fivepointstencil_[2]*temperature(i,j+1,u)) + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Theta_[k])/(fivepointstencil_[4]*temperature(i,j,u)))*(temperature(i,j+1,u)-temperature(i,j,u))  - 0.5*((u[OFF(i,j-1,0)]*DataType::molar_masses()[k]*Mad0_[k]*Theta0_[k])/(fivepointstencil_[0]*temperature(i,j-1,u)) + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Theta_[k])/(fivepointstencil_[4]*temperature(i,j,u)))*(temperature(i,j,u)-temperature(i,j-1,u)) )/ntimes<2>(hy_)
		       
		       ); //minus enclosing expression
      
	if(!is_well_defined_value(diffu[k])){
		       ADONIS_ERROR(ValueError,"diffu["<<k<<"] = "<< diffu[k] << "." );
		     }
	 

// CORRECTION · CORRECTION · CORRECTION · CORRECTION · CORRECTION · CORRECTION
#ifndef NO_CORRECTION_DIFFUSION_VELOCITY //USE CORRECTION VELOCITY
	useMe_ = T(); //reset for each k

		     for(int l = 0; l < nchem; ++l){
		       //!the term above has a minus. Hence this one gets a plus
		       useMe_ += ( //! Y_k is fixed now in this iteration
      		 //+ dx(Y_k rho D_l^mix dy Y_l)
				  (0.5*(Yfrac(i+1,j,u,k)*u[OFF(i+1,j,0)]*Mad1_[l] + Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l])*(Yfrac(i+1,j,u,l)-Yfrac(i,j,u,l)) - 0.5*(Yfrac(i-1,j,u,k)*u[OFF(i-1,j,0)]*Mad3_[l] + Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l])*(Yfrac(i,j,u,l)-Yfrac(i-1,j,u,l)) )/ntimes<2>(hx_)  
      		       // //+ dy(Y_k rho D_l^mix dy Y_l)
				  +  (0.5*(Yfrac(i,j+1,u,k)*u[OFF(i,j+1,0)]*Mad2_[l] + Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l])*(Yfrac(i,j+1,u,l)-Yfrac(i,j,u,l)) - 0.5*(Yfrac(i,j-1,u,k)*u[OFF(i,j-1,0)]*Mad0_[l] + Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l])*(Yfrac(i,j,u,l)-Yfrac(i,j-1,u,l)) )/ntimes<2>(hy_)
      		       //+dx((Y_k rho D_l^mix Y_l)/Wbar dx Wbar)
				  + (0.5*((Yfrac(i+1,j,u,k)*u[OFF(i+1,j,0)]*Mad1_[l]*Yfrac(i+1,j,u,l))/fivepointstencil_[1] + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/fivepointstencil_[4])*(fivepointstencil_[1]-fivepointstencil_[4]) - 0.5*((Yfrac(i-1,j,u,k)*u[OFF(i-1,j,0)]*Mad3_[l]*Yfrac(i-1,j,u,l))/fivepointstencil_[3] + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/fivepointstencil_[4])*(fivepointstencil_[4]-fivepointstencil_[3]) )/ntimes<2>(hx_)
      		       //+dy((Y_k rho D_l^mix Y_l)/Wbar dy Wbar)
      		       + (0.5*((Yfrac(i,j+1,u,k)*u[OFF(i,j+1,0)]*Mad2_[l]*Yfrac(i,j+1,u,l))/fivepointstencil_[2] + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/fivepointstencil_[4])*(fivepointstencil_[2]-fivepointstencil_[4]) - 0.5*((Yfrac(i,j-1,u,k)*u[OFF(i,j-1,0)]*Mad0_[l]*Yfrac(i,j-1,u,l))/fivepointstencil_[0] + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/fivepointstencil_[4])*(fivepointstencil_[4]-fivepointstencil_[0]) )/ntimes<2>(hy_)
      		       //+ dx((Y_k rho D_l^mix Y_l)/p dx p)
      		       + (0.5*((Yfrac(i+1,j,u,k)*u[OFF(i+1,j,0)]*Mad1_[l]*Yfrac(i+1,j,u,l))/pressure(i+1,j,u) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/pressure(i,j,u))*(pressure(i+1,j,u)-pressure(i,j,u)) - 0.5*((Yfrac(i-1,j,u,k)*u[OFF(i-1,j,0)]*Mad3_[l]*Yfrac(i-1,j,u,l))/pressure(i-1,j,u) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/pressure(i,j,u))*(pressure(i,j,u)-pressure(i-1,j,u)) )/ntimes<2>(hx_)
      		       //+ dy((Y_k rho D_l^mix Y_l)/p dy p)
      		       + (0.5*((Yfrac(i,j+1,u,k)*u[OFF(i,j+1,0)]*Mad2_[l]*Yfrac(i,j+1,u,l))/pressure(i,j+1,u) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/pressure(i,j,u))*(pressure(i,j+1,u)-pressure(i,j,u)) - 0.5*((Yfrac(i,j-1,u,k)*u[OFF(i,j-1,0)]*Mad0_[l]*Yfrac(i,j-1,u,l))/pressure(i,j-1,u) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/pressure(i,j,u))*(pressure(i,j,u)-pressure(i,j-1,u)) )/ntimes<2>(hy_)
      		       // - dx ((Y_k rho W_l D_l^mix Y_l)/(Wbar p) dx p)
      		       - (0.5*((Yfrac(i+1,j,u,k)*u[OFF(i+1,j,0)]*DataType::molar_masses()[l]*Mad1_[l]*Yfrac(i+1,j,u,l))/(fivepointstencil_[1]*pressure(i+1,j,u)) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Yfrac(i,j,u,l))/(fivepointstencil_[4]*pressure(i,j,u)))*(pressure(i+1,j,u)-pressure(i,j,u))  - 0.5*((Yfrac(i-1,j,u,k)*u[OFF(i-1,j,0)]*DataType::molar_masses()[l]*Mad3_[l]*Yfrac(i-1,j,u,l))/(fivepointstencil_[3]*pressure(i-1,j,u)) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Yfrac(i,j,u,l))/(fivepointstencil_[4]*pressure(i,j,u)))*(pressure(i,j,u)-pressure(i-1,j,u)) )/ntimes<2>(hx_)
      		        //- dy ((Y_k rho W_l D_l^mix Y_l)/(Wbar p) dy p)
				  - (0.5*((Yfrac(i,j+1,u,k)*u[OFF(i,j+1,0)]*DataType::molar_masses()[l]*Mad2_[l]*Yfrac(i,j+1,u,l))/(fivepointstencil_[2]*pressure(i,j+1,u)) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Yfrac(i,j,u,l))/(fivepointstencil_[4]*pressure(i,j,u)))*(pressure(i,j+1,u)-pressure(i,j,u))  - 0.5*((Yfrac(i,j-1,u,k)*u[OFF(i,j-1,0)]*DataType::molar_masses()[l]*Mad0_[l]*Yfrac(i,j-1,u,l))/(fivepointstencil_[0]*pressure(i,j-1,u)) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Yfrac(i,j,u,l))/(fivepointstencil_[4]*pressure(i,j,u)))*(pressure(i,j,u)-pressure(i,j-1,u)) )/ntimes<2>(hy_)
      		       //+ dx ((Y_k rho Wl D_l^mix theta_l)/(Wbar T) dx T)
      		       + (0.5*((Yfrac(i+1,j,u,k)*u[OFF(i+1,j,0)]*DataType::molar_masses()[l]*Mad1_[l]*Theta1_[l])/(fivepointstencil_[1]*temperature(i+1,j,u)) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Theta_[l])/(fivepointstencil_[4]*temperature(i,j,u)))*(temperature(i+1,j,u)-temperature(i,j,u))  - 0.5*((Yfrac(i-1,j,u,k)*u[OFF(i-1,j,0)]*DataType::molar_masses()[l]*Mad3_[l]*Theta3_[l])/(fivepointstencil_[3]*temperature(i-1,j,u)) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Theta_[l])/(fivepointstencil_[4]*temperature(i,j,u)))*(temperature(i,j,u)-temperature(i-1,j,u)) )/ntimes<2>(hx_)
      		       //dy ((rho Wl D_l^mix theta_l)/(Wbar T) dy T)
      		       + (0.5*((Yfrac(i,j+1,u,k)*u[OFF(i,j+1,0)]*DataType::molar_masses()[l]*Mad2_[l]*Theta2_[l])/(fivepointstencil_[2]*temperature(i,j+1,u)) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Theta_[l])/(fivepointstencil_[4]*temperature(i,j,u)))*(temperature(i,j+1,u)-temperature(i,j,u))  - 0.5*((Yfrac(i,j-1,u,k)*u[OFF(i,j-1,0)]*DataType::molar_masses()[l]*Mad0_[l]*Theta0_[l])/(fivepointstencil_[0]*temperature(i,j-1,u)) + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Theta_[l])/(fivepointstencil_[4]*temperature(i,j,u)))*(temperature(i,j,u)-temperature(i,j-1,u)) )/ntimes<2>(hy_)
				   ); //END useMe_
		     }//END FOR l=0,1,,..., nchem
	  
		     //!add correction velocity
		     diffu[k] += useMe_;		     
#endif //END CORRECTION · CORRECTION · CORRECTION · CORRECTION · CORRECTION
      
		     if(!is_well_defined_value(diffu[k])){
		       ADONIS_ERROR(ValueError,"diffu["<<k<<"] = "<< diffu[k] << "    useMe_ = "<< useMe_ << "." );
		     }
	 
		     tmpvar_ += diffu[k];  
      }//END FOR k=0,1,...,nchem

#ifdef  SHOW_IT_2_ME
	std::cout << "SPECIES DIFFU TERM = "<< diffu << std::endl;
#endif	

	// if(!is_zero(tmpvar_,1.e+05)){
	//   for(int k = 0; k < nchem; ++k)
	//     std::cout << "diffu["<<DataType::species_names()[k]<<"] = "<< diffu[k]<< " ";
	//   std::cout << std::endl;
	//   ADONIS_ERROR(ValueError,"Species diffusive fluxes do not sum to zero, i.e. tmpvar_ = "<< tmpvar_ << ".");
	// }

      return diffu;
    }


     //! computes the thermal diffusion term \f$\sum_{k=1}^K c_{pk}j_k\cdot \nabla T \f$
     //! assumes that cp(i,j,u) has been invoked previously
     T& temperature_diffu_term(int i, int j, const VType& u){
       td_ = T(); //reset
       meanmolmass_ = wbar(i,j,u);

       X_ = mole_fraction(i,j,meanmolmass_,u); //(i,j)
       Mad_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j,u),temperature(i,j,u),X_); //diffusion coefficients are computed by now
       
        //!store counterclockwise, starting from the bottom with
	//!(i,j) being the last entry
        fivepointstencil_[0] = wbar(i,j-1,u);
	fivepointstencil_[1] = wbar(i+1,j,u);
	fivepointstencil_[2] = wbar(i,j+1,u);
	fivepointstencil_[3] = wbar(i-1,j,u);
	fivepointstencil_[4] = meanmolmass_;
	
       
	//! thermal diffusion ratios
	Theta_.compute_thermal_diffusion_ratios(temperature(i,j,u),X_);

       
       for(int k = 0; k < nchem; ++k){
	 td_ += c_pk(i,j,u,k)*( -(  
				 //rho D_k^mix(dx Y_k dx T + dy Y_k dy T)
				 u[OFF(i,j,0)]*Mad_[k]*( (0.5*(Yfrac(i+1,j,u,k) + Yfrac(i,j,u,k))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(Yfrac(i-1,j,u,k) + Yfrac(i,j,u,k))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_)
							 + (0.5*(Yfrac(i,j+1,u,k) + Yfrac(i,j,u,k))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(Yfrac(i,j-1,u,k) + Yfrac(i,j,u,k))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) )
				 //+(rho D_k^mix Y_k)/Wbar(dx Wbar dx T + dy Wbar dy T)
				 + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/fivepointstencil_[4]*( (0.5*(fivepointstencil_[1]+ fivepointstencil_[4])*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(fivepointstencil_[3]+fivepointstencil_[4])*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) 
	+ (0.5*(fivepointstencil_[2]+ fivepointstencil_[4])*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(fivepointstencil_[0]+fivepointstencil_[4])*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) )
				 //+(rho D_k^mix Y_k)/p(dx p dx T + dy p dy T)
				 + (u[OFF(i,j,0)]*Mad_[k]*Yfrac(i,j,u,k))/pressure(i,j,u)*( (0.5*(pressure(i+1,j,u) + pressure(i,j,u))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(pressure(i-1,j,u)+pressure(i,j,u))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) 										    + (0.5*(pressure(i,j+1,u) + pressure(i,j,u))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(pressure(i,j-1,u)+pressure(i,j,u))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) ) 

				 //+(rho W_k D_k^mix Y_k)/(Wbar p)(dx p dx T + dy p dy T)

				 + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Yfrac(i,j,u,k))/(fivepointstencil_[4]*pressure(i,j,u))*( (0.5*(pressure(i+1,j,u) + pressure(i,j,u))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(pressure(i-1,j,u)+pressure(i,j,u))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) 
						       + (0.5*(pressure(i,j+1,u) + pressure(i,j,u))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(pressure(i,j-1,u)+pressure(i,j,u))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) ) 
				 //+(rho W_k D_k^mix Theta_k)/(Wbar T) (dx T dx T + dy T dy T)
				 + (u[OFF(i,j,0)]*DataType::molar_masses()[k]*Mad_[k]*Theta_[k])/(fivepointstencil_[4]*temperature(i,j,u))*( (0.5*(temperature(i+1,j,u) + temperature(i,j,u))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(temperature(i-1,j,u)+temperature(i,j,u))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) 
						       + (0.5*(temperature(i,j+1,u) + temperature(i,j,u))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(temperature(i,j-1,u)+temperature(i,j,u))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) ) 
				   )//enclosing MINUS
			       ); 

//CORRECTION · CORRECTION · CORRECTION · CORRECTION · CORRECTION · CORRECTION
#ifndef NO_CORRECTION_DIFFUSION_VELOCITY  //use CORRECTIVE terms
	 useMe_ = T(); //reset

	 for(int l = 0; l < nchem; ++l){
	   useMe_ += (//no minus here
		     //! Y_k is fixed here 
		     //Y_k rho D_k^mix(dx Y_l dx T + dy Y_l dy T)
		     Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*( (0.5*(Yfrac(i+1,j,u,l) + Yfrac(i,j,u,l))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(Yfrac(i-1,j,u,l) + Yfrac(i,j,u,l))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_)
							 + (0.5*(Yfrac(i,j+1,u,l) + Yfrac(i,j,u,l))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(Yfrac(i,j-1,u,l) + Yfrac(i,j,u,l))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) )
				 //+(Y_k rho D_l^mix Y_l)/Wbar(dx Wbar dx T + dy Wbar dy T)
				 + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/fivepointstencil_[4]*( (0.5*(fivepointstencil_[1]+ fivepointstencil_[4])*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(fivepointstencil_[3]+fivepointstencil_[4])*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) 
	+ (0.5*(fivepointstencil_[2]+ fivepointstencil_[4])*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(fivepointstencil_[0]+fivepointstencil_[4])*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) )
				 //+(Y_k rho D_l^mix Y_l)/p(dx p dx T + dy p dy T)
		     + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*Mad_[l]*Yfrac(i,j,u,l))/pressure(i,j,u)*( (0.5*(pressure(i+1,j,u) + pressure(i,j,u))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(pressure(i-1,j,u)+pressure(i,j,u))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) 										    + (0.5*(pressure(i,j+1,u) + pressure(i,j,u))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(pressure(i,j-1,u)+pressure(i,j,u))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) ) 

				 //+(Y_k rho W_l D_l^mix Y_l)/(Wbar p)(dx p dx T + dy p dy T)
				 + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Yfrac(i,j,u,l))/(fivepointstencil_[4]*pressure(i,j,u))*( (0.5*(pressure(i+1,j,u) + pressure(i,j,u))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(pressure(i-1,j,u)+pressure(i,j,u))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) 
						       + (0.5*(pressure(i,j+1,u) + pressure(i,j,u))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(pressure(i,j-1,u)+pressure(i,j,u))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) ) 
				 //+(Y_k rho W_l D_l^mix Theta_l)/(Wbar T) (dx T dx T + dy T dy T)
				 + (Yfrac(i,j,u,k)*u[OFF(i,j,0)]*DataType::molar_masses()[l]*Mad_[l]*Theta_[l])/(fivepointstencil_[4]*temperature(i,j,u))*( (0.5*(temperature(i+1,j,u) + temperature(i,j,u))*(temperature(i+1,j,u)-temperature(i,j,u)) - 0.5*(temperature(i-1,j,u)+temperature(i,j,u))*(temperature(i,j,u)-temperature(i-1,j,u)))/ntimes<2>(hx_) 
						       + (0.5*(temperature(i,j+1,u) + temperature(i,j,u))*(temperature(i,j+1,u)-temperature(i,j,u)) - 0.5*(temperature(i,j-1,u)+temperature(i,j,u))*(temperature(i,j,u)-temperature(i,j-1,u)))/ntimes<2>(hy_) ) 

		      );
	     } //END for l = 1,..., nchem
	 td_ += (c_pk(i,j,u,k)*useMe_);  //add correction
#endif
//END CORRECTION · CORRECTION · CORRECTION · CORRECTION · CORRECTION
					
		 
       }//END for k = 1, ..., nchem
#ifndef NDEBUG
       if(!is_well_defined_value(td_))
	 ADONIS_ERROR(ValueError," Thermal diffusion term has a bad value");
#endif   

#ifdef SHOW_IT_2_ME
       std::cout << "TEMPERATURE DIFFU TERM = "<< std::endl;
#endif

      return td_;
     }
  

    T pressure_previous_in_time(int i, int j, IteratorTType veciter){
      useMe_ = T();  //reset
      for(int k = 0; k < nchem; ++k){ //! calculate Wbar first
	useMe_ += 
#ifdef NONCONSERVATIVE_FORM	  
	  veciter[OFF(i,j,4+k)]
#else   //flow variable is rho Y_k and must be divided by rho 
	  veciter[OFF(i,j,4+k)]/veciter[OFF(i,j,0)]
#endif
/DataType::molar_masses()[AccessType::spec_access(k)];
      }
      if(is_zero(useMe_)){
	useMe_ = smallVal_;  //prevent division by zero
      }
      useMe_ = 1./useMe_;
      adonis_assert((useMe_ > T()) && (is_well_defined_value(useMe_)));

#ifdef NONCONSERVATIVE_FORM
      return ( (veciter[OFF(i,j,0)]*PhysicalConstants<T>::Rgas*veciter[OFF(i,j,3)])/useMe_ );
#else //flow variable is rho T and must be divided by rho
      return ( (veciter[OFF(i,j,0)]*PhysicalConstants<T>::Rgas*(veciter[OFF(i,j,3)]/veciter[OFF(i,j,0)]))/useMe_ );
#endif
    }


    T& make_zero_when_reasonable(T& val){
      return ( ((val < T()) //&& (is_zero(Abs(val),near_zero_))
		) ? (val = T()) : val );
    }

    
  
    void reasonable_Y(int i, int j, VType& u){
      for(int k = 0; k < nchem; ++k){
	make_zero_when_reasonable(primitive(i,j,u,4+k));
      }
    }


    bool is_out_of_bounds(int i, int j, const VType& u, int quantity, const T& low, const T& up){
      bool res = false;
      if((primitive(i,j,u,quantity) < low) || (primitive(i,j,u,quantity) > up) || (!is_well_defined_value(primitive(i,j,u,quantity))))
	res = true;
      return res;
    }


    VType& repair_out_of_bound_quantity(int i, int j, VType& u, int quantity, const BaseType& low, const BaseType& up, const BaseType& small = BaseType(), bool simple = false){
      if(simple==true){ //!project on bounds
	if(primitive(i,j,u,quantity) < low){
	  assign_flow_variable(i,j,u,quantity,low);
	}
	else if (primitive(i,j,u,quantity) > up){
	  assign_flow_variable(i,j,u,quantity,up);
	}
      
	else if(!is_well_defined_value(primitive(i,j,u,quantity))){
	  adonis_assert((low <= small) && (small <= up));
	  assign_flow_variable(i,j,u,quantity,small); //! some small default value 
	}
	else{
	  //do nothing
	}
      }
      else{
	meanVal_ = 0.0;  //reset
	unsigned count(0);
	bool flag = false;  //is quantity out of bounds?...
	if(is_out_of_bounds(i,j,u,quantity,low,up)){
	  flag = true;      //...yes it is out of bounds at point (i,j) 
	  //check surrounding points; calculate mean value and assign it later
	  if(!is_out_of_bounds(i,j-1,u,quantity,low,up)){
	    meanVal_ +=  primitive(i,j-1,u,quantity);
	    count++;
	  }

	  if(!is_out_of_bounds(i+1,j-1,u,quantity,low,up)){
	    meanVal_ +=  primitive(i+1,j-1,u,quantity);
	    count++;
	  }

	  if(!is_out_of_bounds(i+1,j,u,quantity,low,up)){
	    meanVal_ +=  primitive(i+1,j,u,quantity);
	    count++;
	  }

	  if(!is_out_of_bounds(i+1,j+1,u,quantity,low,up)){
	    meanVal_ +=  primitive(i+1,j+1,u,quantity);
	    count++;
	  }

	  if(!is_out_of_bounds(i,j+1,u,quantity,low,up)){
	    meanVal_ +=  primitive(i,j+1,u,quantity);
	    count++;
	  }
	
	  if(!is_out_of_bounds(i-1,j+1,u,quantity,low,up)){
	    meanVal_ +=  primitive(i-1,j+1,u,quantity);
	    count++;
	  }

	  if(!is_out_of_bounds(i-1,j,u,quantity,low,up)){
	    meanVal_ +=  primitive(i-1,j,u,quantity);
	    count++;
	  }

	  if(!is_out_of_bounds(i-1,j-1,u,quantity,low,up)){
	    meanVal_ +=  primitive(i-1,j-1,u,quantity);
	    count++;
	  }
	}

	if(flag == true){ //value out of bounds
	  if((count > 0)){
	    //mean value of surrounding vals
	    assign_flow_variable(i,j,u,quantity,meanVal_/count);
	  }
	
	  if((count == 0) || (!is_well_defined_value(meanVal_)) || (is_zero(meanVal_))){ //all neighboring nodes' quantities are out of bounds or the computed mean value of surrounding nodes turns out to be a NaN, inf or even approx. zero
	    adonis_assert((low <= small) && (small <= up)); //small must be chosen from within bounds
	    assign_flow_variable(i,j,u,quantity,small);
	    FancyMessages().nice_output("ALL SURROUNDING QUANTITIES OUT OF BOUNDS \nSET SMALL VALUE FROM BOUND INTERVAL "+Num2str(primitive(i,j,u,quantity))+".",35);
	  }
	}
      }
      return u;
    }

    T pressure(int i, int j, const VType& u){
#ifdef PRESSURE_CONST
      return p0_;
#else  // no const pressure
#ifndef NDEBUG
      if((!is_well_defined_value(u[OFF(i,j,0)])) || (u[OFF(i,j,0)] < T()) || (!is_well_defined_value(temperature(i,j,u))) || (temperature(i,j,u) < T()) )
	ADONIS_ERROR(ValueError, "Bad values: rho = "<< u[OFF(i,j,0)] << "  temperature = "<< temperature(i,j,u) << " at point ("<<i<<", "<<j<<").");
#endif
      return ( (u[OFF(i,j,0)]*PhysicalConstants<BaseType>::Rgas*temperature(i,j,u))/wbar(i,j,u) );
    
#endif  //PRESSURE_CONST
    }

    //! see Majda and Braack
    //! \f$ p = \rho^{\gamma}/\exp(-S_0), \f$ where \f$\gamma\f$ is the ratio
    //! of the specific heat capacities
    T hydrodynamical_pressure(int i, int j, const VType& u){
#ifdef PRESSURE_CONST
      return p0_;
#else
      T gamma = cp(i,j,u)/cV(i,j,u);  //CORRECT ????
      T S0 = T();
      for(int k = 0; k < nchem; ++k)
	S0 += BCS_.S_T(k,temperature(i,j,u));   //CORRECT ????
      
      return ( pow(u[OFF(i,j,0)],gamma)/exp(-S0)   );
#endif
    }

    //! \f$ p_{\mathrm{hyd}_1} =  \frac{1}{2}\rho v_1^2\f$
    T hydro_pressure_x(int i, int j, const VType& u){
#ifdef PRESSURE_CONST
      return p0_;
#else
      return ( 0.5*u[OFF(i,j,0)]*ntimes<2>(v1(i,j,u)) );
#endif
    }
    //! \f$ p_{\mathrm{hyd}_2} =  \frac{1}{2}\rho v_2^2\f$ 
    T hydro_pressure_y(int i, int j, const VType& u){
#ifdef PRESSURE_CONST
      return p0_;
#else
      return ( 0.5*u[OFF(i,j,0)]*ntimes<2>(v2(i,j,u)) );
#endif
    }

    //! \f$ p_hyd = \frac{1}{2}\rho \|v\|_2^2\f$ unit: kg/(ms²)=:Pa
    T hydro_pressure(int i, int j, const VType& u){
#ifdef PRESSURE_CONST
      return p0_;
#else
      return 0.5*u[OFF(i,j,0)]*(ntimes<2>(v1(i,j,u)) + ntimes<2>(v2(i,j,u)));
#endif
    }

    //! CAUTION: for nearly compressible fluids, you get a checkboard pattern
    //! use here a windward pressure rate
    //! and \f$\partial_x p\f$ blows up [see ANDERSON, p. 252??]
    T mechanical_compression_x(int i, int j, const VType& u){
      //return ((useMechCompression_==true) ? 1./u[OFF(i,j,0)]*(pressure(i+1,j,u) - pressure(i-1,j,u))/(2*hx_) : 0.);  //via central differences
  
      T res = 0.;
      if(useMechCompression_==true){
#ifdef HYDRO_PRESSURE
#ifndef PRESSURE_GRADIENT_ORDER_1_FD
	res = (hydro_pressure(i+1,j,u)-hydro_pressure(i-1,j,u))/(2*hx_);
#else
	if(v1(i,j,u) <= 0)  //v1 is velocity in x-direction
	  res = ( (hydro_pressure(i+1,j,u)-hydro_pressure(i,j,u))/hx_);
	else 
	  res = ( (hydro_pressure(i,j,u)-hydro_pressure(i-1,j,u))/hx_);
#endif //PRESSURE_GRADIENT_ORDER_1_FD
#else
#ifndef PRESSURE_GRADIENT_ORDER_1_FD
	res = (pressure(i+1,j,u)-pressure(i-1,j,u))/(2*hx_);
#else
	if(v1(i,j,u) <= 0)  //v1 is velocity in x-direction
	  res = ( (pressure(i+1,j,u)-pressure(i,j,u))/hx_);
	else 
	  res = ( (pressure(i,j,u)-pressure(i-1,j,u))/hx_);
#endif //PRESSURE_GRADIENT_ORDER_1_FD
#endif
      }
      else 
	res = 0.;

      return res;
	
    }

    //used before: 1/rho*((p(i,j,u)-p(i,j-1,u))/hy) 
    T mechanical_compression_y(int i, int j, const VType& u){
      //return ((useMechCompression_==true) ? 1./u[OFF(i,j,0)]*(pressure(i,j+1,u) - pressure(i,j-1,u))/(2*hy_) : 0.); //via central differences
      T res = 0.;
       if(useMechCompression_==true){
#ifdef HYDRO_PRESSURE
#ifndef PRESSURE_GRADIENT_ORDER_1_FD
	res = (hydro_pressure(i,j+1,u)-hydro_pressure(i,j-1,u))/(2*hy_);
#else
	if(v1(i,j,u) <= 0)  //v1 is velocity in x-direction
	  res = ( (hydro_pressure(i,j+1,u)-hydro_pressure(i,j,u))/hy_);
	else 
	  res = ( (hydro_pressure(i,j,u)-hydro_pressure(i,j-1,u))/hy_);
#endif //PRESSURE_GRADIENT_ORDER_1_FD
#else
#ifndef PRESSURE_GRADIENT_ORDER_1_FD
	res = (pressure(i,j+1,u)-pressure(i,j-1,u))/(2*hy_);
#else
	if(v1(i,j,u) <= 0)  //v1 is velocity in x-direction
	  res = ( (pressure(i,j+1,u)-pressure(i,j,u))/hy_);
	else 
	  res = ( (pressure(i,j,u)-pressure(i,j-1,u))/hy_);
#endif //PRESSURE_GRADIENT_ORDER_1_FD
#endif
      }
      else 
	res = 0.;
      return res;
	
    }

    //! mechanical compression for thermal energy equation
    //! \f$ \partial_t p + \nabla \cdot (pv) \f$
    //! \f$ \partial_t p + \partial_x (pv1) + partial_y (pv2) \f$
    T& energy_mechanical_expression(int i, int j, const VType& u){
      useMe_ = T(); //reset
      //!if ok Veciter_ stores composition and stepsize_ the time stepsize
      //adonis_assert(isPrevPressInvoked_); 
      //adonis_assert(stepsize_ > T());
     
      //! only return \f$\partial_t p  \f$ when it is 
      if((isPrevPressInvoked_) && (stepsize_ > T())){
	useMe_= (pressure(i,j,u)-pressure_previous_in_time(i,j,Veciter_))/stepsize_;
      }
      useMe_ += (
	  
#ifndef PRESSURE_GRADIENT_ORDER_1_FD
#ifdef NONCONSERVATIVE_FORM
	      v1(i,j,u)*((pressure(i+1,j,u)-pressure(i-1,j,u))/(2*hx_)) +
	      v2(i,j,u)*((pressure(i,j+1,u)-pressure(i,j-1,u))/(2*hy_))
#else //conservative form
	      (pressure(i+1,j,u)*v1(i+1,j,u)-pressure(i-1,j,u)*v1(i-1,j,u))/(2*hx_) + (pressure(i,j+1,u)*v2(i,j+1,u)-pressure(i,j-1,u)*v2(i,j-1,u))/(2*hy_)
#endif 

#else //1st order gradient

#ifdef NONCONSERVATIVE_FORM

	      ( ((v1(i,j,u) <= 0.0)  ? ( v1(i,j,u)*(pressure(i+1,j,u) - pressure(i,j,u))/hx_ ) : ( v1(i,j,u)*(pressure(i,j,u) - pressure(i-1,j,u))/hx_ ))  +
		( (v2(i,j,u) <= 0.0) ? ( v2(i,j,u)*(pressure(i,j+1,u) - pressure(i,j,u))/hy_ ) : ( v2(i,j,u)*(pressure(i,j,u) - pressure(i,j-1,u))/hy_ ) ) ) 
    
#else //conservative form for 1st order convection
	      ( ((v1(i,j,u) <= 0.0)  ? ((pressure(i+1,j,u)*v1(i+1,j,u) - pressure(i,j,u)*v1(i,j,u))/hx_ ) : ( (pressure(i,j,u)*v1(i,j,u) - pressure(i-1,j,u)*v1(i-1,j,u))/hx_ ))  +
		( (v2(i,j,u) <= 0.0) ? ( (pressure(i,j+1,u)*v2(i,j+1,u) - pressure(i,j,u)*v2(i,j,u))/hy_ ) : ( (pressure(i,j,u)*v2(i,j,u) - pressure(i,j-1,u)*v2(i,j-1,u))/hy_ ) ) ) 	      

#endif
#endif //PRESSURE_GRADIENT_ORDER_1_FD

		   ); //assignment of useMe_
      //#ifndef NDEBUG
	if(!is_well_defined_value(useMe_))
	  ADONIS_ERROR(ValueError, "mechanical compression is not well defined, i.e. d_t p + nabla·(p v) = "<< useMe_ << ".");
	//#endif

#ifdef SHOW_IT_2_ME
	std::cout << "ENERGY MECH. COMPRESSION = "<< useMe_ << std::endl;
#endif

	//std::cout << "mech_compression("<<i<<","<<j<<") = "<< useMe_ << std::endl;
	return useMe_;
      
    }



    T gravitational_force_x(int i, int j, const VType& u){
      return ( (useGravi_==true) ? gravi1_  : 0.);
    }

    T gravitational_force_y(int i, int j, const VType& u){
      return ( (useGravi_==true) ? gravi2_  : 0.);
    }

//! boundary indeces
int down_j() const{
  return
#ifndef GHOST_POINTS_INCLUDED
    1
#else
    2
#endif
    ;
}

  int down_jj() const{
  return
#ifndef GHOST_POINTS_INCLUDED
    2
#else
    3
#endif
    ;
}

int up_j() const{
  return
#ifndef GHOST_POINTS_INCLUDED
    ny_-2
#else
    ny_-3
#endif
    ;
}

int up_jj() const{
  return
#ifndef GHOST_POINTS_INCLUDED
    ny_-3
#else
    ny_-4
#endif
    ;
}
  
int left_i() const{
  return
#ifndef GHOST_POINTS_INCLUDED
    1
#else
    2
#endif
    ;
}

  int left_ii() const{
  return
#ifndef GHOST_POINTS_INCLUDED
    2
#else
    3
#endif
    ;
}

 int right_i() const{
  return
#ifndef GHOST_POINTS_INCLUDED
    nx_-2
#else
    nx_-3
#endif
    ;
} 

int right_ii() const{
  return
#ifndef GHOST_POINTS_INCLUDED
    nx_-3
#else
    nx_-4
#endif
    ;
} 
  
// ZERO NEUMANN · ZERO NEUMANN · ZERO NEUMANN · ZERO NEUMANN · ZERO NEU
    // ZERO NEUMANN · ZERO NEUMANN · ZERO NEUMANN · ZERO NEUMANN · ZERO NEU
    // ZERO NEUMANN · ZERO NEUMANN · ZERO NEUMANN · ZERO NEUMANN · ZERO NEU
    //! returns in any case the primitve variable
    //!==== first order zero Neumann condition -- PRIMITIVE VARS 
    T first_order_ZERO_Neumann_dy_down(int i, const VType& u, int quantity){
      return 
#ifdef NONCONSERVATIVE_FORM
	u[OFF(i,(*this).down_j(),quantity)];         //u_i,0 = u_i,1
#else  //conservative form
      ((quantity == 0) ? (u[OFF(i,(*this).down_j(),0)]) : (u[OFF(i,(*this).down_j(),quantity)]/u[OFF(i,(*this).down_j(),0)]));
#endif
    }

    T first_order_ZERO_Neumann_dy_up(int i, const VType& u, int quantity){
      return
#ifdef NONCONSERVATIVE_FORM
	u[OFF(i,(*this).up_j(),quantity)];     //u_i,ny-1 = u_i,ny-2
#else //conservative form
      ((quantity == 0) ? (u[OFF(i,(*this).up_j(),0)]) : (u[OFF(i,(*this).up_j(),quantity)]/u[OFF(i,(*this).up_j(),0)]));
#endif
    }

    T first_order_ZERO_Neumann_dx_left(int j, const VType& u, int quantity){
      return 
#ifdef NONCONSERVATIVE_FORM	
	u[OFF((*this).left_i(),j,quantity)];         //u_0,j = u_1,j
#else //conservative form
      ((quantity == 0) ? (u[OFF((*this).left_i(),j,0)]) : (u[OFF((*this).left_i(),j,quantity)]/u[OFF((*this).left_i(),j,0)]) );
#endif
    } 

    T first_order_ZERO_Neumann_dx_right(int j, const VType& u, int quantity){
      return 
#ifdef NONCONSERVATIVE_FORM
	u[OFF((*this).right_i(),j,quantity)];  //u_nx-1,j = u_nx-2,j
#else //conservative form
      ((quantity == 0) ? (u[OFF((*this).right_i(),j,0)]) : (u[OFF((*this).right_i(),j,quantity)]/u[OFF((*this).right_i(),j,0)]) );
#endif
    }

    //!====, see [ANDERSON, "Computational Fluid Dynamics", § 4.2, p. 139]
    //!====, see [ANDERSON, "Computational Fluid Dynamics", § 10.3.5, p.458] 
    T second_order_ZERO_Neumann_dy_down(int i, const VType& u, int quantity){
#ifdef EXTRAPOLATION_BDY
      return 
#ifdef NONCONSERVATIVE_FORM
	( 2*u[OFF(i,(*this).down_j(),quantity)] - u[OFF(i,(*this).down_jj(),quantity)] );
#else
      ((quantity == 0) ? (( 2*u[OFF(i,(*this).down_j(),0)] - u[OFF(i,(*this).down_jj(),0)] )) : (( 2*u[OFF(i,(*this).down_j(),quantity)]/u[OFF(i,(*this).down_j(),0)] - u[OFF(i,(*this).down_jj(),quantity)]/u[OFF(i,(*this).down_jj(),0)] )))
#endif //NONCONSERVATIVE_FORM
#else //no extrapolation
	
      return 
#ifdef NONCONSERVATIVE_FORM
(-1./3.*(u[OFF(i,(*this).down_jj(),quantity)] - 4*u[OFF(i,(*this).down_j(),quantity)]));
#else
      ( (quantity == 0) ? (-1./3.*(u[OFF(i,(*this).down_jj(),0)] - 4*u[OFF(i,(*this).down_j(),0)])) : (-1./3.*(u[OFF(i,(*this).down_jj(),quantity)]/u[OFF(i,(*this).down_jj(),0)] - 4*u[OFF(i,(*this).down_j(),quantity)]/u[OFF(i,(*this).down_j(),0)])) );
#endif //NONCONSERVATIVE_FORM
#endif    
}

    T second_order_ZERO_Neumann_dy_up(int i, const VType& u, int quantity){
#ifdef EXTRAPOLATION_BDY

      return 
#ifdef NONCONSERVATIVE_FORM
( 2*u[OFF(i,(*this).up_j(),quantity)] - u[OFF(i(*this).up_jj(),,quantity)] );
#else //conservative
      ( (quantity == 0) ? (2*u[OFF(i,(*this).up_j(),0)] - u[OFF(i,(*this).up_jj(),0)]) : (2*u[OFF(i,(*this).up_j(),quantity)]/u[OFF(i,(*this).up_j(),0)] - u[OFF(i,(*this).up_jj(),quantity)]/u[OFF(i,(*this).up_jj(),0)]));
#endif //NONCONSERVATIVE_FORM

#else //no extrapolation
      return 
#ifdef NONCONSERVATIVE_FORM
(-1./3.*(u[OFF(i,(*this).up_jj(),quantity)] - 4*u[OFF(i,(*this).up_j(),quantity)]));
#else
      ((quantity == 0) ? (-1./3.*(u[OFF(i,(*this).up_jj(),0)] - 4*u[OFF(i,(*this).up_j(),0)])) : (-1./3.*(u[OFF(i,(*this).up_jj(),quantity)]/u[OFF(i,(*this).up_jj(),0)]- 4*u[OFF(i,(*this).up_j(),quantity)]/u[OFF(i,(*this).up_j(),0)])) );
#endif //NONCONSERVATIVE_FORM
#endif   
 }

    T second_order_ZERO_Neumann_dx_left(int j, const VType& u, int quantity){
#ifdef EXTRAPOLATION_BDY

      return 
#ifdef NONCONSERVATIVE_FORM
	( 2*u[OFF((*this).left_i(),j,quantity)] - u[OFF((*this).left_ii(),j,quantity)]);
#else
      ((quantity == 0) ? (2*u[OFF((*this).left_i(),j,0)] - u[OFF((*this).left_ii(),j,0)]) : (2*u[OFF((*this).left_i(),j,quantity)]/u[OFF((*this).left_i(),j,0)] - u[OFF((*this).left_ii(),j,quantity)]/u[OFF((*this).left_ii(),j,0)]) );  
#endif //NONCONSERVATIVE_FORM
#else  //no extrapolation
      return
#ifdef NONCONSERVATIVE_FORM
 (-1./3.*(u[OFF((*this).left_ii(),j,quantity)] - 4*u[OFF((*this).left_i(),j,quantity)]));
#else //conservative
      ( (quantity == 0) ? ((-1./3.*(u[OFF((*this).left_ii(),j,0)] - 4*u[OFF((*this).left_i(),j,0)]))) : ((-1./3.*(u[OFF((*this).left_ii(),j,quantity)]/u[OFF((*this).left_ii(),j,0)] - 4*u[OFF((*this).left_i(),j,quantity)]/u[OFF((*this).left_i(),j,0)]))) );        
#endif //NONCONSERVATIVE_FORM
#endif    
}
    
    T second_order_ZERO_Neumann_dx_right(int j, const VType& u, int quantity){
#ifdef EXTRAPOLATION_BDY

      return 
#ifdef NONCONSERVATIVE_FORM
( 2*u[OFF((*this).right_i(),j,quantity)] - u[OFF((*this).right_ii(),j,quantity)]);
#else //conservative
      ((quantity == 0) ? (2*u[OFF((*this).right_i(),j,0)] - u[OFF((*this).right_ii(),j,0)]) : (2*u[OFF((*this).right_i(),j,quantity)]/u[OFF((*this).right_i(),j,0)] - u[OFF((*this).right_ii(),j,quantity)]/u[OFF((*this).right_ii(),j,0)]) );
#endif //NONCONSERVATIVE_FORM
#else  //no extrapolation

      return 
#ifdef NONCONSERVATIVE_FORM
(-1./3.*(u[OFF((*this).right_ii(),j,quantity)] - 4*u[OFF((*this).right_i(),j,quantity)]));
#else //conservative
      ((quantity == 0) ? (-1./3.*(u[OFF((*this).right_ii(),j,0)] - 4*u[OFF((*this).right_i(),j,0)])) : (-1./3.*(u[OFF((*this).right_ii(),j,quantity)]/u[OFF((*this).right_ii(),j,0)]- 4*u[OFF((*this).right_i(),j,quantity)]/u[OFF((*this).right_i(),j,0)])) );
#endif //NONCONSERVATIVE_FORM
#endif
    }

  //pressure handling at the boundaries
  //file nseconaire.hh
//lower boundary (0st order) -- for all boundaries the same
  T boundary_pressure_0th_order(int i, const VType& u){  
    return p0_; // atmospheric pressure, function args have no effect
  }


  //lower boundary (1st order)
  T boundary_pressure_1st_order_down(int i, const VType& u){
    return pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      1
#else
      2
#endif
      ,u); // p_i,0 = p_i,1
  }

   //lower boundary (1st order)
  T boundary_pressure_1st_order_up(int i, const VType& u){
    return pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      ny_-2
#else
      ny_-3
#endif   
      ,u); // p_i,ny-1, = p_i,ny-2
  }

  T boundary_pressure_1st_order_left(int j, const VType& u){
    return pressure(
#ifndef GHOST_POINTS_INCLUDED
      1
#else
      2
#endif   
      ,j,u); // p_0,j = p1,j
  }

  T boundary_pressure_1st_order_right(int j, const VType& u){
    return pressure(
#ifndef GHOST_POINTS_INCLUDED
      nx_-2
#else
     nx_-3
#endif   
      ,j,u); // p_nx-1,j = p_nx-2,j
  }


  //lower boundary (2nd order)
  // p_i,0 = (2 p_i,1) - p_i,2, c.f.  [ANDERSON, Computational Fluid Dynamics, p. 458, Eq. (10.17)]
 // u1 = 4/3u2 - 1/3 u3, cf. [ANDERSON, Computational Fluid Dynamics, p. 138, Eq. (4.27)]
  T boundary_pressure_2nd_order_down(int i, const VType& u){
#ifdef EXTRAPOLATION_BDY 
    return ( 2*pressure(i,
 #ifndef GHOST_POINTS_INCLUDED
      1
#else
     2
#endif   
      ,u) - pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      2
#else
     3
#endif   
      ,u) ); 
#else
    return ( -1./3.*(pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      2
#else
     3
#endif   
      ,u) - 4*pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      1
#else
     2
#endif   
      ,u)) );
#endif
      }
  
  T boundary_pressure_2nd_order_up(int i,  const VType& u){
#ifdef EXTRAPOLATION_BDY 
    return ( 2*pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      ny_-2
#else
     ny_-3
#endif   
      ,u) - pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      ny_-3
#else
     ny_-4
#endif   

      ,u) ); 
#else
    return ( -1./3.*(pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      ny_-3
#else
     ny_-4
#endif   
      ,u) - 4*pressure(i,
#ifndef GHOST_POINTS_INCLUDED
      ny_-2
#else
     ny_-3
#endif   

      ,u)) );
#endif
  }

  T boundary_pressure_2nd_order_left(int j, const VType& u){
#ifdef EXTRAPOLATION_BDY 
    return ( 2*pressure(
#ifndef GHOST_POINTS_INCLUDED
      1
#else
     2
#endif   
      ,j,u) - pressure(
#ifndef GHOST_POINTS_INCLUDED
      2
#else
     3
#endif   

      ,j,u) ); 
#else
    return ( -1./3.*(pressure(
#ifndef GHOST_POINTS_INCLUDED
      2
#else
     3
#endif   

      ,j,u) - 4*pressure(
#ifndef GHOST_POINTS_INCLUDED
      1
#else
     2
#endif   
      ,j,u)) );
#endif
  }

  T boundary_pressure_2nd_order_right(int j, const VType& u){
#ifdef EXTRAPOLATION_BDY 
    return ( 2*pressure(
#ifndef GHOST_POINTS_INCLUDED
      nx_-2
#else
     nx_-3
#endif   
	,j,u) - pressure(
#ifndef GHOST_POINTS_INCLUDED
      nx_-3
#else
     nx_-4
#endif   
      ,j,u) ); 
#else
    return ( -1./3.*(pressure(
#ifndef GHOST_POINTS_INCLUDED
      nx_-3
#else
     nx_-4
#endif   

      ,j,u) - 4*pressure(
#ifndef GHOST_POINTS_INCLUDED
      nx_-2
#else
     nx_-3
#endif   
      ,j,u)) );
#endif
  }



    void check_for_good_values(const T& val, int i, int j, int quantity = -1, const std::string& addinfo = std::string(), bool chc_if_zero = false){
      // adonis_assert((0 <= quantity) && (quantity < nprim));
      	std::string quantname;
	if((is_well_defined_value(val) == false)
#ifdef CHECK_4_TOO_LARGE_VALUES
	   || (Abs(val) > tooLarge_)
#endif
	   ){ //val is NaN or Inf
	  if(quantity == 0)
	    quantname = "rho";
	  else if (quantity == 1)
	    quantname = "v_1";
	  else if (quantity == 2)
	    quantname = "v_2";
	  else if (quantity == 3)
	    quantname = "T";
	  else if (quantity >= 4)
	    quantname = "Y_"+DataType::species_names()[4-quantity];
	  else if (quantity < 0) //no specification of species is made
	    quantname = "value";
	  else{
	    //do nothing
	  }
	  std::string tl;
#ifdef CHECK_4_TOO_LARGE_VALUES
	  tl = " -- Value too large";
#endif
	  ADONIS_ERROR(ValueError,"BAD VALUE occurred. "<< ((quantity < 0) ? " ":" d/dt ") << quantname<<" = "<< val << " at grid point ("<<i<<", "<<j<<").\n   "<< (addinfo+tl) << ".");
	}
	if(chc_if_zero){
	  if(is_zero(val))
	    ADONIS_ERROR(ValueError,"Zero value occurs which is prohibitive here"<< ((quantity < 0) ? " ":" in d/dt ") << quantname<<" at grid point ("<<i<<", "<<j<<").\n   "<<addinfo<< ".");
	}
    }    

  
    T& prevent_blow_up(T& val, const BaseType& maxval){
      adonis_assert(maxval >= BaseType());
      return ( (Abs(val) > Abs(maxval)) ? (val = Sgn(val)*maxval) : val);
    }

    //! overload prev. function
    inline T& prevent_blow_up(T& val, const BaseType& minval, const BaseType& maxval){
      adonis_assert(minval <= maxval);
      return ((val < minval) ? (val = minval) : ((val > maxval) ? (val = maxval) : val));
    }


    void repair_prim_vars_in_9_points(int i, int j, VType& u){
       //! temperature o.k.  -- T
	  repair_out_of_bound_quantity(i,j,u,3,280.,3500., 300.,true);
	  repair_out_of_bound_quantity(i,j-1,u,3,280.,3500., 300.,true);
	  repair_out_of_bound_quantity(i+1,j-1,u,3,280.,3500., 300.,true);
	  repair_out_of_bound_quantity(i+1,j,u,3,280.,3500., 300.,true);
	  repair_out_of_bound_quantity(i+1,j+1,u,3,280.,3500., 300.,true);
	  repair_out_of_bound_quantity(i,j+1,u,3,280.,3500., 300.,true);
	  repair_out_of_bound_quantity(i-1,j+1,u,3,280.,3500., 300.,true);
	  repair_out_of_bound_quantity(i-1,j,u,3,280.,3500., 300.,true);
	  repair_out_of_bound_quantity(i-1,j-1,u,3,280.,3500., 300.,true);
	  //safe guard for mass fractions -- Y
	  reasonable_Y(i,j,u);
	  reasonable_Y(i,j-1,u);
	  reasonable_Y(i+1,j-1,u);
	  reasonable_Y(i+1,j,u);
	  reasonable_Y(i+1,j+1,u);
	  reasonable_Y(i,j+1,u);
	  reasonable_Y(i-1,j+1,u);
	  reasonable_Y(i-1,j,u);
	  reasonable_Y(i-1,j-1,u);
	  //! repair when negative density is encountered -- rho
	  //!CAUTION: lower bound zero is desastrous! Take small value instead,
	  //          otherwise you risk production of infs or nans
	  repair_out_of_bound_quantity(i,j,u,0,rhoMin_,1.e+20,smallVal_,true); 
	  repair_out_of_bound_quantity(i,j-1,u,0,rhoMin_,1.e+20, smallVal_,true);
	  repair_out_of_bound_quantity(i+1,j-1,u,0,rhoMin_, 1.e+20, smallVal_,true);
	  repair_out_of_bound_quantity(i+1,j,u,0,rhoMin_, 1.e+20, smallVal_,true);
	  repair_out_of_bound_quantity(i+1,j+1,u,0,rhoMin_, 1.e+20, smallVal_,true);
	  repair_out_of_bound_quantity(i,j+1,u,0,rhoMin_, 1.e+20, smallVal_,true);
	  repair_out_of_bound_quantity(i-1,j+1,u,0,rhoMin_, 1.e+20, smallVal_,true);
	  repair_out_of_bound_quantity(i-1,j,u,0,rhoMin_, 1.e+20, smallVal_,true);
	  repair_out_of_bound_quantity(i-1,j-1,u,0,rhoMin_, 1.e+20, smallVal_,true);

	  //! prevent blow-up of velocity --  v1 and v2
	  prevent_blow_up(v1(i,j,u),v1_max_);
	  prevent_blow_up(v1(i,j-1,u),v1_max_);
	  prevent_blow_up(v1(i+1,j-1,u),v1_max_);
	  prevent_blow_up(v1(i+1,j,u),v1_max_);
	  prevent_blow_up(v1(i+1,j+1,u),v1_max_);
	  prevent_blow_up(v1(i,j+1,u),v1_max_);
	  prevent_blow_up(v1(i-1,j+1,u),v1_max_);
	  prevent_blow_up(v1(i-1,j,u),v1_max_);
	  prevent_blow_up(v1(i-1,j-1,u),v1_max_);
	  //---
	  prevent_blow_up(v2(i,j,u),v2_max_);
	  prevent_blow_up(v2(i,j-1,u),v2_max_);
	  prevent_blow_up(v2(i+1,j-1,u),v2_max_);
	  prevent_blow_up(v2(i+1,j,u),v2_max_);
	  prevent_blow_up(v2(i+1,j+1,u),v2_max_);
	  prevent_blow_up(v2(i,j+1,u),v2_max_);
	  prevent_blow_up(v2(i-1,j+1,u),v2_max_);
	  prevent_blow_up(v2(i-1,j,u),v2_max_);
	  prevent_blow_up(v2(i-1,j-1,u),v2_max_);

    }

    T rho_low_mach(int i, int j, const VType& u){
#ifdef NONCONSERVATIVE_FORM //ONLY PRIMITIVE VARS USED HERE!!!!!!!!
      adonis_assert(!is_zero(cp_));

      T sm = T();
      for(int k = 0; k < nchem; ++k){
	sm += (-diffspec_[k] + dotomega_[k]*DataType::molar_masses()[k]);
      }
      //! NOTE: you need NO minus here when assigning to \f$\partial_t rho\f$
      return ( u[OFF(i,j,1)]*upwind_downwind_x(i,j,0,u) + u[OFF(i,j,2)]*upwind_downwind_y(i,j,0,u) -1./(cp_*u[OFF(i,j,3)])*(tempCond_ - tempDiff_ - tempHeat_) + sm );
      
#else
      ADONIS_ERROR(VariantError,"Low Mach number formulation of rho can only be used in the non-conservative form of the differential equations.");
#endif
    }

    //! get flow variables depending on (non-)conservative form
    T temperature(int i, int j, const VType& u) const{
#ifdef NONCONSERVATIVE_FORM
      return u[OFF(i,j,3)];
#else  //conservative form
      return u[OFF(i,j,3)]/u[OFF(i,j,0)];	
#endif
    }

      
    T v1(int i, int j, const VType& u) const{
#ifdef NONCONSERVATIVE_FORM
      return u[OFF(i,j,1)];
#else  //conservative form
      return u[OFF(i,j,1)]/u[OFF(i,j,0)];	
#endif
    }

    T v2(int i, int j, const VType& u) const{
#ifdef NONCONSERVATIVE_FORM
      return u[OFF(i,j,2)];
#else  //conservative form
      return u[OFF(i,j,2)]/u[OFF(i,j,0)];	
#endif
    }
        
    T Yfrac(int i, int j, const VType& u, int k) const{
#ifdef NONCONSERVATIVE_FORM
      return u[OFF(i,j,4+k)];
#else  //conservative form
      return u[OFF(i,j,4+k)]/u[OFF(i,j,0)];	
#endif
    }   

    //some fluid mechanical numbers
    //! Reynolds number in x direction, cf. [ANDERSON, p. 456, Eqs. (10.14a) and (10.14b)]
    T Reynolds_x(int i, int j, const VType& u) const{
      return ( (u[OFF(i,j,0)]*v1(i,j,u)*hx_)/mu(i,j,u) ); 
    }

    //Reynolds number in y direction
    T Reynolds_y(int i, int j, const VType& u) const{
      return ( (u[OFF(i,j,0)]*v2(i,j,u)*hy_)/mu(i,j,u) ); 
    }

  }; //END OF CLASS

} //end namespace

#endif 
