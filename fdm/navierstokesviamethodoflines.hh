#ifndef NAVIER_STOKES_EQUATIONS_VIA_METHOD_OF_LINES_HH
#define NAVIER_STOKES_EQUATIONS_VIA_METHOD_OF_LINES_HH
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

#include "choose.hh"
#include "boundary.hh"

namespace Adonis{

  template<class T>
  class ReactiveNS2D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType; 
   
    //typedef typename TypeAdapter<T>::Type DType; //no CppAd
    typedef typename TypeAdapter<T>::BaseType BaseType; //for norms, tolerances 
    typedef std::size_t SizeType;
    
     //needed for thermo chemistry
    typedef ThermoData4Mechanism<BaseType,ChooseFDSetting::MECH> DataType;  

    typedef Boundary2D<VType,ChooseFDSetting::MECH> BoundaryType;
    
    enum{nchem = SystemDimension<DataType::nspec,DataType::rednspec,ChooseFDSetting::ISREDUCED>::Value, 
	 //rho,v1,v2,T + specs
	 nprim = 4+SystemDimension<DataType::nspec,DataType::rednspec,ChooseFDSetting::ISREDUCED>::Value};  //4+nchem
    
    typedef StaticArray<T,nchem> ArrayType; 
   
    
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;

    //Other Index
    typedef IndexHandler<DataType,ChooseFDSetting::ISREDUCED> AccessType;
    typedef typename AccessType::IndexerType IndexerType;

    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,IndexerType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    typedef EquationOfState<'i'> EosType;

    //TRANSPORT COEFFICIENTS
    typedef MixtureAveragedViscosity<DataType,ChooseFDSetting::ISREDUCED> ViscosityType;
    typedef MixtureAveragedConductivity<DataType,ChooseFDSetting::ISREDUCED> ConductivityType;
    typedef MixtureAveragedDiffusion<DataType,ChooseFDSetting::ISREDUCED> DiffusionType;

    
    // this may be of some importance
    typedef ComputeTransportProperties<ChooseFDSetting::ISREDUCED> PropType;

    ReactiveNS2D(SizeType dim = 0):rhs_(dim){
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
      t0_ = PD.get_datum<BaseType>("t0");
      tend_ = PD.get_datum<BaseType>("tend");
      halfdiam_ = 0.5*diam_;
      p0_ = PD.get_datum<BaseType>("pconst");
      v1_ = PD.get_datum<BaseType>("v1_in");
      v2_ = PD.get_datum<BaseType>("v2_in");
      Tinlet_ = PD.get_datum<BaseType>("T_in");

      bdytemp_.initialize(PD.get_datum<BaseType>("T_min"),PD.get_datum<BaseType>("T_ignition"),PD.get_datum<BaseType>("hottest_x"));

      //transport coefficients are ready for computation
      Mav_.initialize(DataType::transport());
      Mac_.initialize(DataType::transport());
      Mad_.initialize(DataType::transport());  //(i,j)
      Mad0_.initialize(DataType::transport()); //(i,j-1)
      Mad1_.initialize(DataType::transport()); //(i+1,j)
      Mad2_.initialize(DataType::transport()); //(i,j+1)
      Mad3_.initialize(DataType::transport()); //(i-1,j)

      bdy2d_.initialize(DataType::nspec); //alway in the full dimension!

      // X_.resize(nchem);  //(i,j) //mole fractions vector for the species
      // X0_.resize(nchem); //(i,j-1)
      // X1_.resize(nchem); //(i+1,j)
      // X2_.resize(nchem); //(i,j+1)
      // X3_.resize(nchem); //(i-1,j)
      // diffspec_.resize(nchem);
      
      v1bdy_.resize(ny_);   //for the velocity profile (along y direction)
      u_.resize(npt_);

      // C_.resize(nchem);
      // dotomega_.resize(nchem);
      // Cp_.resize(nchem);
      cp_ = T();

      // isCpcalculated_ = false;
      // isDiffusioncalculated_ = false;
      // isDotOmegacalculated_ = false;
    }

    //! \f$\rho, v_1, v_2, T\f$ + nchem species, discretized on \f$\Omega_h\f$
      //! with npt_ nodes
    int dim() const {return totpoints_;}
    int domain_dim() const {return totpoints_;}

    std::string name() const {return "with momentum";}

    template<class VEC>
    VType& operator()(VEC& u){ //reference here for we wanna change it 
      //HARD COPY
      u_ = u;

      bdy2d_(); //boundary composition activated
      //calculate mean molecular mass
      //++
      //for boundaries
      T wbar;
      PropType::mean_molecular_weight_Y<DataType>(wbar,&bdy2d_[0],DataType::molar_masses());

      //BOUNDARY
      //UP
      BaseType x = BaseType();
      for(int i = 0; i < nx_; ++i){
	x = a_ + i*hx_;         //spatial coordinate x
	BaseType bdytmp = bdytemp_(x);
	//temperature
	u_[OFF(i,ny_-1,3)] = bdytmp; 

	u_[OFF(i,ny_-1,0)] = EosType::density(p0_,wbar,bdytmp); //rho
       
	// //velocity v = (v1,v2)
	u_[OFF(i,ny_-1,1)] = 0.;
	u_[OFF(i,ny_-1,2)] = 0.;	  

	//species 
	for(int k = 0; k < nchem; ++k)
	  u_[OFF(i,ny_-1,4+k)] = bdy2d_[k]; //=====TODO: NEUMANN
      }


      //DOWN: now NEUMANN too   =====TODO
      // for(int i = 0; i < nx_; ++i){
      // 	u_[OFF(i,0,0)] = u_[OFF(i,1,0)];
      	
      // 	u_[OFF(i,0,1)] = v1_; //..but not for veloc
      // 	u_[OFF(i,0,2)] = 0;
      	
      // 	u_[OFF(i,0,3)] = u_[OFF(i,1,3)];
      
      // 	for(int k = 0; k < nchem;++k)
      // 	  u_[OFF(i,0,4+k)] = u_[OFF(i,1,4+k)];

      // }

      //RIGHT: NEUMANN
      	for(int j = 1; j < ny_-1; ++j){
	  u_[OFF(nx_-1,j,0)] = u_[OFF(nx_-2,j,0)]; //rho 
 
	  u_[OFF(nx_-1,j,1)] = u_[OFF(nx_-2,j,1)];  //v1bdy_[j]; 
	  u_[OFF(nx_-1,j,2)] = u_[OFF(nx_-2,j,2)];


	  u_[OFF(nx_-1,j,3)] = u_[OFF(nx_-2,j,3)]; //temp
	  for(int k = 0; k < nchem; ++k){ //species
	    u_[OFF(nx_-1,j,4+k)] = u_[OFF(nx_-2,j,4+k)];
	  }
	}
	
      //LEFT (INLET)
      BaseType radius(0);
      for(int j = 0; j < ny_; ++j){
	v1bdy_[j] = zero(parabolic_inlet_velocity(radius,diam_,v1_));
	radius += hy_;
      }

      for(int j = 1; j < ny_-1; ++j){
	u_[OFF(0,j,1)] = v1bdy_[j];  //v1
	u_[OFF(0,j,2)] = 0;          //v2
	BaseType tin = bdytemp_(0);
	u_[OFF(0,j,3)] =  tin; //Tinlet_; //TODO: T wie ad 1.)
	u_[OFF(0,j,0)] = EosType::density(p0_,wbar,tin); //rho
	
	for(int k = 0; k < nchem; ++k)   //species
	  u_[OFF(0,j,4+k)] = bdy2d_[k];
      }

      
      
      //--------------------------------------------------------
	// //BOUNDARY -- full symmetric scenario
	// // DOWN and UP:
	// BaseType x = BaseType();
	// for(int i = 0; i < nx_; ++i){
	//   x = a_ + i*hx_;         //spatial coordinate x
	//   double bdytmp = bdytemp_(x);
	//   //++ rho
	//   u_[OFF(i,0,0)] = u_[OFF(i,ny_-1,0)] = EosType::density(p0_,wbar,bdytmp);

	//   //velocity v = (v1,v2)
	//   u_[OFF(i,0,1)] = 0.;
	//   u_[OFF(i,ny_-1,1)] = 0.;
	//   u_[OFF(i,0,2)] = 0.;
	//   u_[OFF(i,ny_-1,2)] = 0.;

	  
	//   //temperature
	//   u_[OFF(i,0,3)] = bdytmp;
	//   u_[OFF(i,ny_-1,3)] = bdytmp; 

	//   //species 
	//   for(int k = 0; k < nchem; ++k)
	//     u_[OFF(i,0,4+k)] = u_[OFF(i,ny_-1,4+k)] = bdy2d_[k]; 
	// }

	// //LEFT
	// //! compute parabolic velocity profile
	// BaseType r = -halfdiam_;
	// for(int j = 0; j < ny_; ++j){
	//   v1bdy_[j] = zero(parabolic_inlet_velocity(r,halfdiam_,v1_));
	//   r += hy_;
	// }


	// for(int j = 1; j < ny_-1; ++j){
	//   //velocity
	//   u_[OFF(0,j,1)] = v1bdy_[j];
	//   u_[OFF(0,j,2)] = 0.;

	//   //temperature
	//   u_[OFF(0,j,3)] = Tinlet_; //TODO
	  
	//   for(int k = 0; k < nchem; ++k){ //species
	//     u_[OFF(0,j,4+k)] = bdy2d_[AccessType::spec_access(k)];
	//   }
	// }

	// // RIGHT -- Neuman 
	// for(int j = 1; j < ny_-1; ++j){
	//   u_[OFF(nx_-1,j,0)] = u_[OFF(nx_-2,j,0)]; //rho 
 
	//   u_[OFF(nx_-1,j,1)] = u_[OFF(nx_-2,j,1)];  //v1bdy_[j]; 
	//   u_[OFF(nx_-1,j,2)] = u_[OFF(nx_-2,j,2)];


	//   u_[OFF(nx_-1,j,3)] = u_[OFF(nx_-2,j,3)]; //temp
	//   for(int k = 0; k < nchem; ++k){ //species
	//     u_[OFF(nx_-1,j,4+k)] = u_[OFF(nx_-2,j,4+k)];
	//   }
	// }

	//INTERIOR -- use second order approximations for first, sec. and 
	//            mixed derivatives
	for(int i = 1; i < nx_-1; ++i){
	  x = BaseType();      
	  radius = BaseType();
	  for(int j = 1; j < ny_-1; ++j){
	    
	    //! preserve mass Y_nspec = 1 - \sum_k^{nspec-1} Y_k
	    preserve_mass(i,j,u_);   

	    //check direction of velocity
	    T dx_v1 = upwind_downwind_x(i,j,1,u_[OFF(i,j,1)],u_),
	      dy_v1 = upwind_downwind_y(i,j,1,u_[OFF(i,j,2)],u_),
	      dx_v2 = upwind_downwind_x(i,j,2,u_[OFF(i,j,1)],u_),
	      dy_v2 = upwind_downwind_y(i,j,2,u_[OFF(i,j,2)],u_),
	      dx_T =  upwind_downwind_x(i,j,3,u_[OFF(i,j,1)],u_),
	      dy_T =  upwind_downwind_y(i,j,3,u_[OFF(i,j,2)],u_);

	    for(int k = 0; k < nchem; ++k){
	      dx_spec_[k] = upwind_downwind_x(i,j,4+k,u_[OFF(i,j,1)],u_);
	      dy_spec_[k] = upwind_downwind_y(i,j,4+k,u_[OFF(i,j,2)],u_);
	    }

	    //EQUATION OF STATE
	    u_[OFF(i,j,0)] = EosType::density(p0_,(*this).wbar(i,j,u_),u_[OFF(i,j,3)]);
 
	    //!\f$\rho\f$
	    rhs_[OFF(i,j,0)] = -( (u[OFF(i+1,j,0)]*u[OFF(i+1,j,1)] - u[OFF(i,j,0)]*u[OFF(i,j,1)])/hx_ + (u[OFF(i,j+1,0)]*u_[OFF(i,j+1,2)] - u[OFF(i,j,0)]*u_[OFF(i,j,2)])/hy_ );//-( u_[OFF(i,j,0)]*( (u_[OFF(i+1,j,1)] - u_[OFF(i-1,j,1)])/(2.*hx_) +  (u_[OFF(i,j+1,2)] - u_[OFF(i,j-1,2)])/(2.*hy_) ) );

	    
	    //!\f$v_1\f$
	    // u_[OFF(i,j,1)] = zero(parabolic_inlet_velocity(radius,diam_,v1_));
	    // radius += hy_;
	    

	    rhs_[OFF(i,j,1)] = -( //(ntimes<2>(u_[OFF(i+1,j,1)]) - ntimes<2>(u_[OFF(i-1,j,1)]))/(2.*hx_) + (u_[OFF(i,j+1,1)]*u_[OFF(i,j+1,2)] - u_[OFF(i,j-1,1)]*u_[OFF(i,j-1,2)])/(2.*hy_) 
				 u[OFF(i,j,1)]*dx_v1 + u[OFF(i,j,2)]*dy_v1 
				 - 1./u_[OFF(i,j,0)]*(4./3.*( (0.5*(mu(i+1,j,u_)+mu(i,j,u_))*(u_[OFF(i+1,j,1)]-u_[OFF(i,j,1)]) - 0.5*(mu(i-1,j,u_)+mu(i,j,u_))*(u_[OFF(i,j,1)]-u_[OFF(i-1,j,1)]))/ntimes<2>(hx_)) - 2./3.*(mu(i+1,j,u_)*(u_[OFF(i+1,j+1,2)] - u_[OFF(i+1,j-1,2)]) - mu(i-1,j,u_)*(u_[OFF(i-1,j+1,2)]-u_[OFF(i-1,j-1,2)]))/(4.*hx_*hy_)) - 1./u_[OFF(i,j,0)]*( (0.5*(mu(i,j+1,u_)+mu(i,j,u_))*(u_[OFF(i,j+1,1)] - u_[OFF(i,j,1)]) - 0.5*(mu(i,j-1,u_)+mu(i,j,u_))*(u_[OFF(i,j,1)] - u_[OFF(i,j-1,1)]))/ntimes<2>(hy_) + (mu(i,j+1,u_)*(u_[OFF(i+1,j+1,2)] - u_[OFF(i-1,j+1,2)]) - mu(i,j-1,u_)*(u_[OFF(i+1,j-1,2)] - u_[OFF(i-1,j-1,2)]))/(4.*hx_*hy_)));
	   

	    //! \f$v_2\f$  //========== TODO: = 0 ???????
	    rhs_[OFF(i,j,2)] = -(//(u_[OFF(i+1,j,2)]*u_[OFF(i+1,j,1)] - u_[OFF(i-1,j,2)]*u_[OFF(i-1,j,1)])/(2.*hy_) + (ntimes<2>(u_[OFF(i,j+1,2)]) - ntimes<2>(u_[OFF(i,j-1,2)]))/(2.*hx_) 
	     u[OFF(i,j,1)]*dx_v2 + u[OFF(i,j,2)]*dy_v2 
	     - 1./u_[OFF(i,j,0)]*( (mu(i+1,j,u_)*(u_[OFF(i+1,j+1,1)] - u_[OFF(i+1,j-1,1)]) - mu(i-1,j,u_)*(u_[OFF(i-1,j+1,1)]-u_[OFF(i-1,j-1,1)]))/(4.*hx_*hy_) + (0.5*(mu(i+1,j,u_)+mu(i,j,u_))*(u_[OFF(i+1,j,2)]-u_[OFF(i,j,2)]) - 0.5*(mu(i-1,j,u_)+mu(i,j,u_))*(u_[OFF(i,j,2)]-u_[OFF(i-1,j,2)]))/ntimes<2>(hx_) ) -  1./u_[OFF(i,j,0)]*(4./3.*((0.5*(mu(i,j+1,u_)+mu(i,j,u_))*(u_[OFF(i,j+1,2)] - u_[OFF(i,j,2)]) - 0.5*(mu(i,j-1,u_)+mu(i,j,u_))*(u_[OFF(i,j,2)] - u_[OFF(i,j-1,2)]))/ntimes<2>(hy_)) - 2./3.*((mu(i,j+1,u_)*(u_[OFF(i+1,j+1,1)] - u_[OFF(i-1,j+1,1)]) - mu(i,j-1,u_)*(u_[OFF(i+1,j-1,1)] - u_[OFF(i-1,j-1,1)]))/(4.*hx_*hy_))));
				 
	    //!TEMPERATURE \f$ T \f$
	    
#ifndef NDEBUG  
	    is_temperature_ok(i,j,0,u_);  // 0 = take bounds of first chem. spec
#endif	   
	    //project_temperature_on_its_bounds(i,j,0,u); //brings nothing

	    cp_ = cp(i,j,u_); //Cp_ has also been filled
	   

	    rhs_[OFF(i,j,3)] = -( //(u_[OFF(i+1,j,3)]*u_[OFF(i+1,j,1)] - u_[OFF(i-1,j,3)]*u_[OFF(i-1,j,1)])/(2.*hx_) + (u_[OFF(i,j+1,3)]*u_[OFF(i,j+1,2)] - u_[OFF(i,j-1,3)]*u_[OFF(i,j-1,2)])/(2.*hy_) 
				 u[OFF(i,j,1)]*dx_T + u[OFF(i,j,2)]*dy_T 
- 1./(u_[OFF(i,j,0)]*cp_)*( (0.5*(lambda(i+1,j,u_)+lambda(i,j,u_))*(u_[OFF(i+1,j,3)]-u_[OFF(i,j,3)]) - 0.5*(lambda(i,j,u_)+lambda(i-1,j,u_))*(u_[OFF(i,j,3)]-u_[OFF(i-1,j,3)]))/ntimes<2>(hx_) + (0.5*(lambda(i,j+1,u_)+lambda(i,j,u_))*(u_[OFF(i,j+1,3)] - u_[OFF(i,j,3)]) - 0.5*(lambda(i,j-1,u_)+lambda(i,j,u_))*(u_[OFF(i,j,3)]-u_[OFF(i,j-1,3)]))/ntimes<2>(hy_)) 
				   -               //PLUS/MINUS ??
				 1./(cp_*u[OFF(i,j,0)])*temperature_diffu_term(i,j,u_)
				  +                 //PLUS ??
				  1./(cp_*u_[OFF(i,j,0)])*heat_production(i,j,u_) ); //a call to heat_production computes  dotomega_
	    
	    //! SPECIES \f$ Y_k, \qquad k = 1, \ldots, nchem\f$
	    
	    //species_diffu_term(i,j,u_); //now diffspec_ is filled

	    species_diffu_term(i,j,u_); //now diffspec_ is filled

	    

	    for(int k = 0; k < nchem; ++k){
	      rhs_[OFF(i,j,4+k)] = -( //(u_[OFF(i+1,j,4+k)]*u_[OFF(i+1,j,1)] - u_[OFF(i-1,j,4+k)]*u_[OFF(i-1,j,1)])/(2.*hx_) + (u_[OFF(i,j+1,4+k)]*u_[OFF(i,j+1,2)] - u_[OFF(i,j-1,4+k)]*u_[OFF(i,j-1,2)])/(2.*hy_) //- diffspec_[k] //old
				      
 u[OFF(i,j,1)]*dx_spec_[k] + u[OFF(i,j,2)]*dy_spec_[k]
                                      -              //PLUS/MINUS???
                                      1./u_[OFF(i,j,0)]*diffspec_[k]
                                      -
				      1./u_[OFF(i,j,0)]*dotomega_[k]*DataType::molar_masses()[AccessType::spec_access(k)] );
	    }
	    

	  } //END INTERIOR
	}
	  
      return rhs_;
    }

  private:
    VType rhs_;
    BuilderType BCS_;

    BaseType a_,b_,c_,d_,hx_,hy_, t0_, tend_, diam_, halfdiam_;
    int nx_, ny_, npt_, totpoints_;   
   
    BaseType p0_, v1_,v2_, Tinlet_; //standard pressure
    //VType X_, X0_, X1_,X2_,X3_,diffspec_, C_, dotomega_, Cp_;
    ArrayType X_, X0_, X1_,X2_,X3_,diffspec_, C_, dotomega_, Cp_;
    VType v1bdy_,
      u_;
    T cp_; //for storage
    
    T fivepointstencil_[5];
    T wbarfive_[5];

    T dx_spec_[nchem];
    T dy_spec_[nchem];


    T jflux_[nchem];

    BoundaryTemperatureDistribution<BaseType> bdytemp_;
    ViscosityType Mav_;
    ConductivityType Mac_;
    DiffusionType Mad_,Mad0_,Mad1_,Mad2_,Mad3_;
   
    BoundaryType bdy2d_; 
    //ProjectTemperature<T,ChooseFDSetting::MECH> ProjT_;

    // isCpcalculated_ = false;
    // isDiffusioncalculated_ = false;
    // isDotOmegacalculated_ = false;

    //PRIVATE FUNCTIONS -- don't use them from outside

    int OFF(int i, int j, int spec) const{
      return (i + nx_*j + spec*npt_);
    }

   
    // void species_mass_flux(int i, int j, const VType& u){
      
      
    // }


    //!PRESERVE MASS
    void preserve_mass(int i, int j, VType& u){ //u is changed accoringly
      T sm1 = T();
      
      //! compute \f$ Y_{n} = 1 - \sum_{k = 1}^{n-1} Y_k\f$
      for(int k = 0; k < nchem-1; ++k)
	sm1 += u[OFF(i,j,4+k)];
      
      u[OFF(i,j,4+(nchem-1))] = 1. - sm1;
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
	    
	  
	  if(u[OFF(IX,JX,3)] < DataType::temperature_bounds()[3*k] || 
	     u[OFF(IX,JX,3)] > DataType::temperature_bounds()[3*k+1]){
	   
	    ADONIS_ERROR(BoundsError,"Temperature bounds violated: T_{i"<<oper1<<l << ",j"<<oper2<<t << "} = "<< u[OFF(IX,JX,3)]<<", i.e. T("<<IX<<","<<JX<<").");
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
      for(int k = 0; k < nchem; ++k){
	wb += u[OFF(i,j,4+k)]/DataType::molar_masses()[AccessType::spec_access(k)];
      }
      return 1./wb;
    }
    
     T mole_fraction(int i, int j, int k, const T& Wbar, const VType& u) const{
      return u[OFF(i,j,4+k)]*Wbar/DataType::molar_masses()[AccessType::spec_access(k)];
    }
    
    //! needed to evaluate diffusion and viscosity coefficients 
    ArrayType& mole_fraction(int i, int j, const T& Wbar, const VType& u){
      for(int k = 0; k < nchem; ++k)
	X_[k] = mole_fraction(i,j,k,Wbar,u);
      return X_;
    }
    
    //!overload mole fraction vector for <II>any</II> random access container
    template<class RAC>
    RAC& mole_frac(RAC& xfrac, int i, int j, const T& Wbar, const VType& u){
      for(int k = 0; k < nchem; ++k)
	xfrac[k] = mole_fraction(i,j,k,Wbar,u);

      return xfrac;
    }

    T concenctration(int i, int j, int k, const VType& u, const T& rho) const{
      return rho*u[OFF(i,j,4+k)]/DataType::molar_masses()[AccessType::spec_access(k)];
    }

    ArrayType& concentration(int i, int j, const VType& u){
      for(int k = 0; k < nchem; ++k)
	C_[k] = concenctration(i,j,k,u,u[OFF(i,j,0)]);
      return C_;
    }

    //! overload concentration 
    template<class RAC>
    RAC& concentrations(RAC& conc, int i, int j, const VType& u){
      for(int k = 0; k < nchem; ++k)
	conc[k] = concenctration(i,j,k,u,u[OFF(i,j,0)]);
      return conc;
    }

    //////////////////////////////////////////////////////////////////////////
    //SOME VECTORS ARE EXPLICITLY FILLED DURING CALL ==> NO NEED TO RE-COMP 
    //////////////////////////////////////////////////////////////////////////
     //! compute \f$ \dot{\omega} \f$, unit: mol/(m³s)
    ArrayType& chemical_net_production(int i, int j, const VType& u){
      concentration(i,j,u);  //now C_ contains concentrations
      for(int k = 0; k < nchem; ++k){
	//! ==== TODO: construction of the chem. source term should 
	//! ========== index correctly, right?????????????
	dotomega_[k] = BCS_.net_formation_rate(k,u[OFF(i,j,3)],C_);
      }
      //std::cout << "CHEMISTRY:  dot{omega} = "<< dotomega_ << std::endl;
      return dotomega_;
	
    }
    
     
    //! \f$ \dot{\omega} \f$ is computed by now
    T heat_production(int i, int j, const VType & u){
      T h = T();

      chemical_net_production(i,j,u); //now dotomega_ contains net prod. rate
      //isDotOmegacalculated_ = true;
      
      for(int k = 0; k < nchem; ++k)
	h += BCS_.H_T(AccessType::spec_access(k),u[OFF(i,j,3)])*dotomega_[k];
      return h;
    }

     T cp(int i, int j, const VType& u){
      T cres = T();
      for(int k = 0; k < nchem; ++k){
	//! fill Cp_ simultaneously
#ifndef NDEBUG
	
	if(u[OFF(i,j,3)] < DataType::temperature_bounds()[3*AccessType::spec_access(k)] || u[OFF(i,j,3)] > DataType::temperature_bounds()[3*AccessType::spec_access(k) + 1])
	  ADONIS_ERROR(BoundsError, " Temperature out of bound T_{"<<i<<","<<j<< "} = "<< u[OFF(i,j,3)] << ".");
#endif
	Cp_[k] = BCS_.C_p(AccessType::spec_access(k),u[OFF(i,j,3)]);
	cres += Cp_[k]/DataType::molar_masses()[AccessType::spec_access(k)]*u[OFF(i,j,4+k)];
      }
      //isCpcalculated_ = true;
      return cres;
    }

    //! EVALUATE TRANSPORT COEFFICIENTS AT point node \f$ (i,j) \f$
    //! compute \f$ \mu(T,X)\f$
    T mu(int i, int j, const VType& u){  
      T m = Mav_.compute_mixture_averaged_viscosity(u[OFF(i,j,3)], mole_fraction(i,j,wbar(i,j,u),u));
      return m;
    }

    //! \f$ \lambda(\rho, p, T, X) \f$
    T lambda(int i, int j, const VType& u){
      return Mac_.compute_mixture_averaged_conductivity(u[OFF(i,j,0)], p0_, u[OFF(i,j,3)],mole_fraction(i,j,wbar(i,j,u),u));
    }

   
    //! assumes that cp(i,j,u) has been invoked previously
     T temperature_diffu_term(int i, int j, const VType& u){
       T td = T(),
	 meanmolmass = wbar(i,j,u);

       X_ = mole_fraction(i,j,meanmolmass,u);
       Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,j,3)],X_); //diffusion coefficients are computed by now
       
       //isDiffusioncalculated_ = true;
       
       for(int k = 0; k < nchem; ++k){
	 //!store counterclockwise, starting from the bottom with
	//!(i,j) being the last entry
	fivepointstencil_[0] = mole_fraction(i,j-1,k,wbar(i,j-1,u),u);
	fivepointstencil_[1] = mole_fraction(i+1,j,k,wbar(i+1,j,u),u);
	fivepointstencil_[2] = mole_fraction(i,j+1,k,wbar(i,j+1,u),u);
	fivepointstencil_[3] = mole_fraction(i-1,j,k,wbar(i-1,j,u),u);
	fivepointstencil_[4] = X_[k]; //mole_fraction(i,j,k,meanmolmass,u);

	//! assume Cp_ has been filled previously
	td += ( Cp_[k]/meanmolmass*Mad_[k]*u[OFF(i,j,0)]*( (0.5*(fivepointstencil_[1]+fivepointstencil_[4])*(u[OFF(i+1,j,3)]-u[OFF(i,j,3)]) - 0.5*(fivepointstencil_[3]+fivepointstencil_[4])*(u[OFF(i,j,3)]-u[OFF(i-1,j,3)]))/ntimes<2>(hx_) + (0.5*(fivepointstencil_[2]+fivepointstencil_[4])*(u[OFF(i,j+1,3)]-u[OFF(i,j,3)]) - 0.5*(fivepointstencil_[4]+fivepointstencil_[0])*(u[OFF(i,j,3)]-u[OFF(i,j-1,3)]))/ntimes<2>(hy_)) );
      }
      return td;
     }
  

    //perhaps much better to replace \f$ Y_k/X_k\f$ by \f$W_k/\overline{W}\f$ since \f$X_k\f$ can be zero
    ArrayType& species_diffu_term(int i, int j, const VType& u){
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

      Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,j-1,3)],X0_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i+1,j,3)],X1_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,j+1,3)],X2_);
      
      Mad3_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i-1,j,3)],X3_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,j,3)],X_);
      
      for(int k = 0; k < nchem; ++k){
	diffspec_[k] = (0.5*(u[OFF(i+1,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[1]*Mad1_[k] + u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k])*(X1_[k]-X_[k]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k] + u[OFF(i-1,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[3]*Mad3_[k])*(X_[k]-X3_[k]))/ntimes<2>(hx_) + 
	  (0.5*(u[OFF(i,j+1,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[2]*Mad2_[k] + u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k])*(X2_[k]-X_[k]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k] + u[OFF(i,j-1,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[0]*Mad0_[k])*(X_[k]-X0_[k]))/ntimes<2>(hy_);

	//! OLD Version, division by zero could principally occur
	// //!division by mole fractions can become problematic
	// diffspec_[0] = ( (0.5*(u[OFF(i+1,j,4+k)]/X1_[k]*Mad1_[k] + u[OFF(i,j,4+k)]/X_[k]*Mad_[k])*(X1_[k] - X_[k]) - 0.5*(u[OFF(i-1,j,4+k)]/X3_[k]*Mad3_[k] + u[OFF(i,j,4+k)]/X_[k]*Mad_[k])*(X_[k] - X3_[k]))/ntimes<2>(hx_) 
	// 		 + (0.5*(u[OFF(i,j+1,4+k)]*X3_[k]/Mad2_[k] + u[OFF(i,j,4+k)]/X_[k]*Mad_[k])*(X2_[k] - X_[k]) - 0.5*(u[OFF(i,j-1,4+k)]/X0_[k]*Mad0_[k]+u[OFF(i,j,4+k)]/X_[k]*Mad_[k])*(X_[k] - X0_[k]))/ntimes<2>(hy_));
	

      }
      
      return diffspec_;
    }


    //! version with mean molecular mass
    //! implements 2. and 3. term of [KEE, Eq. (3.128), p. 96]
    //! must be divided by \f$ 1/\rho \f$
    ArrayType& species_diffu_term_Wbar(int i, int j, const VType& u){
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

      Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,j-1,3)],X0_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i+1,j,3)],X1_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,j+1,3)],X2_);
      
      Mad3_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i-1,j,3)],X3_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,j,3)],X_);
      
      for(int k = 0; k < nchem; ++k){
	//!\f$\nabla \cdot (\rho D^{\mathrm{mix}}_k \nabla Y_k)\f$
	diffspec_[k] = (0.5*(u[OFF(i+1,j,0)]*Mad1_[k] + u[OFF(i,j,0)]*Mad_[k])*(u[OFF(i+1,j,4+k)]-u[OFF(i,j,4+k)]) - 0.5*(u[OFF(i,j,0)]*Mad_[k] + u[OFF(i-1,j,0)]*Mad3_[k])*(u[OFF(i,j,4+k)]-u[OFF(i-1,j,4+k)]))/ntimes<2>(hx_) 
	  + (0.5*(u[OFF(i,j+1,0)]*Mad2_[k] + u[OFF(i,j,0)]*Mad_[k])*(u[OFF(i,j+1,4+k)]-u[OFF(i,j,4+k)]) - 0.5*(u[OFF(i,j,0)]*Mad_[k] + u[OFF(i,j-1,0)]*Mad0_[k])*(u[OFF(i,j,4+k)]-u[OFF(i,j-1,4+k)]))/ntimes<2>(hy_) 
	  //!differential part with mean molecular weight \f$ \overline{W}\f$,
	  //!i.e. \f$ \nabla \cdot (\rho Y_k/\overline{W}D^{\mathrm{mix}}_k \nabla \overline{W})\f$
	  +  (0.5*(u[OFF(i+1,j,0)]*u[OFF(i+1,j,4+k)]/fivepointstencil_[1]*Mad1_[k] + u[OFF(i,j,0)]*u[OFF(i,j,4+k)]/fivepointstencil_[4]*Mad_[k])*(fivepointstencil_[1]-fivepointstencil_[4]) - 0.5*(u[OFF(i,j,0)]*u[OFF(i,j,4+k)]/fivepointstencil_[4]*Mad_[k] + u[OFF(i-1,j,0)]*u[OFF(i-1,j,4+k)]/fivepointstencil_[3]*Mad3_[k])*(fivepointstencil_[4]-fivepointstencil_[3]))/ntimes<2>(hx_)
	  +  (0.5*(u[OFF(i,j+1,0)]*u[OFF(i,j+1,4+k)]/fivepointstencil_[2]*Mad2_[k] + u[OFF(i,j,0)]*u[OFF(i,j,4+k)]/fivepointstencil_[4]*Mad_[k])*(fivepointstencil_[2]-fivepointstencil_[4]) - 0.5*(u[OFF(i,j,0)]*u[OFF(i,j,4+k)]/fivepointstencil_[4]*Mad_[k] + u[OFF(i,j-1,0)]*u[OFF(i,j-1,4+k)]/fivepointstencil_[0]*Mad0_[k])*(fivepointstencil_[4]-fivepointstencil_[0]))/ntimes<2>(hy_);
      }
      
      return diffspec_;
    }


  };

} //end namespace

#endif 
