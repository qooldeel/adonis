#ifndef NSE_VIA_MOL_FOR_CONAIRES_PRIMARY_REFERENCE_FUEL_H2_MECH_HH
#define NSE_VIA_MOL_FOR_CONAIRES_PRIMARY_REFERENCE_FUEL_H2_MECH_HH
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

//#include "omp.h"
//#define NUMTHREADS 8

namespace Adonis{

  template<class T>
  class FlameIgnitionConaireH2ReactiveNS2D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType; 
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<DType> VDType; 
    //typedef typename TypeAdapter<T>::Type DType; //no CppAd
    typedef typename TypeAdapter<T>::BaseType BaseType; //for norms, tolerances 
    typedef std::size_t SizeType;
    
     //needed for thermo chemistry: primary reference fuel mechanism
    typedef ThermoData4Mechanism<BaseType,10> DataType;  

    
    enum{nchem = DataType::nspec, 
	 //rho,v1,v2,T + specs
	 nprim = 4+DataType::nspec};  //4+nchem
    
    typedef StaticArray<T,DataType::nspec> ArrayType; 
    typedef StaticArray<T,2> DimensionType;
    //! element 0 is x-direction, 1 y-direction and so forth
    

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

    typedef ThermalDiffusionRatio<DataType,false> RatioType;
    
    // this may be of some importance
    //typedef ComputeTransportProperties<false> PropType;

    FlameIgnitionConaireH2ReactiveNS2D(SizeType dim = 0):rhs_(dim){
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
      p0_ = PD.get_datum<BaseType>("pconst");
      v1_fuel_ = PD.get_datum<BaseType>("v1_fuel");
      v2_fuel_ = PD.get_datum<BaseType>("v2_fuel");
      v1_oxi_ = PD.get_datum<BaseType>("v1_oxi");
      v2_oxi_ = PD.get_datum<BaseType>("v2_oxi");
      T_fuel_ = PD.get_datum<BaseType>("T_fuel");
      T_oxi_ = PD.get_datum<BaseType>("T_oxi");
      T_wall_ = PD.get_datum<BaseType>("T_wall");
      
      std::cout << "v1_fuel = "<< v1_fuel_ << "    v1_oxi = "<< v1_oxi_ << std::endl;

      Tignition_ = PD.get_datum<BaseType>("T_ignition");
      yhottest_ = PD.get_datum<BaseType>("hottest_y");
      bdytemp_.initialize(0.5*(PD.get_datum<BaseType>("T_fuel")+PD.get_datum<BaseType>("T_oxi")),Tignition_,yhottest_);

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

      //O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR
      bdyFuel_[0] = PD.get_datum<BaseType>("O_FUEL");
      bdyFuel_[1]  = PD.get_datum<BaseType>("O2_FUEL");
      bdyFuel_[2]  = PD.get_datum<BaseType>("H_FUEL");
      bdyFuel_[3]  = PD.get_datum<BaseType>("OH_FUEL");
      bdyFuel_[4]  = PD.get_datum<BaseType>("H2_FUEL");
      bdyFuel_[5]  = PD.get_datum<BaseType>("HO2_FUEL");
      bdyFuel_[6]  = PD.get_datum<BaseType>("H2O2_FUEL");
      bdyFuel_[7]  = PD.get_datum<BaseType>("H2O_FUEL");
      bdyFuel_[8]  = PD.get_datum<BaseType>("N2_FUEL");
      bdyFuel_[9]  = PD.get_datum<BaseType>("AR_FUEL");
      
      bdyOxi_[0] = PD.get_datum<BaseType>("O_OXI");
      bdyOxi_[1]  = PD.get_datum<BaseType>("O2_OXI");
      bdyOxi_[2]  = PD.get_datum<BaseType>("H_OXI");
      bdyOxi_[3]  = PD.get_datum<BaseType>("OH_OXI");
      bdyOxi_[4]  = PD.get_datum<BaseType>("H2_OXI");
      bdyOxi_[5]  = PD.get_datum<BaseType>("HO2_OXI");
      bdyOxi_[6]  = PD.get_datum<BaseType>("H2O2_OXI");
      bdyOxi_[7]  = PD.get_datum<BaseType>("H2O_OXI");
      bdyOxi_[8]  = PD.get_datum<BaseType>("N2_OXI");
      bdyOxi_[9]  = PD.get_datum<BaseType>("AR_OXI");
      


     
      
      T air[] = { PD.get_datum<BaseType>("O"),    //0.,
		  PD.get_datum<BaseType>("O2"),   //0.21,   //0.21
		  PD.get_datum<BaseType>("H"),    //0.,
		  PD.get_datum<BaseType>("OH"),   //0.,
		  PD.get_datum<BaseType>("H2"),    //0.,
		  PD.get_datum<BaseType>("HO2"),   //0.
		  PD.get_datum<BaseType>("H2O2"),     //0.,
		  PD.get_datum<BaseType>("H2O"),      //0.01,   //0.01
		  PD.get_datum<BaseType>("N2"),
		  PD.get_datum<BaseType>("AR")
      };

      for(int k = 0; k < DataType::nspec; ++k)
	IV_[k] = air[k];


      inBetween_[O] = PD.get_datum<BaseType>("O_BW");
      inBetween_[O2] = PD.get_datum<BaseType>("O2_BW");
      inBetween_[H] = PD.get_datum<BaseType>("H_BW");
      inBetween_[OH] = PD.get_datum<BaseType>("OH_BW");
      inBetween_[H2] = PD.get_datum<BaseType>("H2_BW");
      inBetween_[HO2] = PD.get_datum<BaseType>("HO2_BW");
      inBetween_[H2O2] = PD.get_datum<BaseType>("H2O2_BW");
      inBetween_[H2O] = PD.get_datum<BaseType>("H2O_BW");
      inBetween_[N2] = PD.get_datum<BaseType>("N2_BW");
      inBetween_[AR] = PD.get_datum<BaseType>("AR_BW");

      BaseType sm(0.);
      for(int k = 0; k < DataType::nspec; ++k)
	sm += inBetween_[k];
      if(sm != 1.)
	ADONIS_ERROR(DerivedError,"Inbetween data do not sum to 1!");
    }

    //! \f$\rho, v_1, v_2, T\f$ + nchem species, discretized on \f$\Omega_h\f$
      //! with npt_ nodes
    int dim() const {return totpoints_;}
    int domain_dim() const {return totpoints_;}

    std::string name() const {return "with momentum";}

    template<class VEC>
    VType& operator()(VEC& vars){ //reference here for we wanna change it 
      //HARD COPY
      u_ = vars;
     
      //!BOUNDARY -- imagine szenario is rotated by 90°
      //!  LEFT: inlet
      DType y(0.);
      for(int j = 1; j < ny_-1; ++j){
	y = j*hy_;

	u_[OFF(0,j,0)] = (101325*wbar(0,j,u_))/(PhysicalConstants<DType>::Rgas*u_[OFF(0,j,3)]);
	
	if(y < yhottest_){
	  u_[OFF(0,j,1)] = v1_oxi_;
	  u_[OFF(0,j,2)] = v2_oxi_;
	  for(int k = 0; k < nchem; ++k)  //TODO=====
	    u_[OFF(0,j,4+k)] = bdyOxi_[k];
	}
	else{
	  u_[OFF(0,j,1)] = v1_fuel_;
	  u_[OFF(0,j,2)] = v2_fuel_;
	  for(int k = 0; k < nchem; ++k)   //TODO====
	    u_[OFF(0,j,4+k)] = bdyFuel_[k];
	}

	//small ignition area
	if(is_contained(yhottest_,y-hy_,y+hy_)){
	  u_[OFF(0,j,3)] = Tignition_;
	  
	  //!==========  TODO: commented the next 2 blocks
	  // u_[OFF(0,j,1)] = 0.5*(v1_oxi_+v1_fuel_); //v1 velocity mean
	  //u_[OFF(0,j,2)] = 0.0;

	   //overwrite now species where mixture is used
	   // for(int k = 0; k < DataType::nspec; ++k)
	   //   u_[OFF(0,j,4+k)] = inBetween_[k];
	
	   // rho now alters with temperature an composition
	   u_[OFF(0,j,0)] = (101325*wbar(0,j,u_))/(PhysicalConstants<DType>::Rgas*Tignition_);
	}
	else 
	  u_[OFF(0,j,3)] = 300.;

	  //}
       

	//u_[OFF(0,j,3)] = bdytemp_(y); //Temperature
      }

      //! DOWN: wall
      for(int i = 0; i < nx_; ++i){
	u_[OFF(i,0,0)] = (pressure(i,0,u_)*wbar(i,0,u_))/(PhysicalConstants<DType>::Rgas*u_[OFF(i,0,3)]);
	u_[OFF(i,0,1)] = 0.;
	u_[OFF(i,0,2)] = 0.;
	u_[OFF(i,0,3)] = T_wall_;
	for(int k = 0; k < nchem; ++k)
	  u_[OFF(i,0,4+k)] = u_[OFF(i,1,4+k)]; //Neumann
      }
      //! RIGHT: outlet
      for(int j = 1; j < ny_-1; ++j){ 
	u_[OFF(nx_-1,j,0)] = (pressure(nx_-1,j,u_)*wbar(nx_-1,j,u_))/(PhysicalConstants<DType>::Rgas*u_[OFF(nx_-1,j,3)]);
	//Neumann
	u_[OFF(nx_-1,j,1)] = u_[OFF(nx_-2,j,1)];
	u_[OFF(nx_-1,j,2)] = u_[OFF(nx_-2,j,2)];
	u_[OFF(nx_-1,j,3)] = u_[OFF(nx_-2,j,3)];
	for(int k = 0; k < nchem; ++k)
	  u_[OFF(nx_-1,j,4+k)] = u_[OFF(nx_-2,j,4+k)];
      }

      //!UP: center
      for(int i = 0; i < nx_; ++i){
	u_[OFF(i,ny_-1,0)] = (pressure(i,ny_-1,u_)*wbar(i,ny_-1,u_))/(PhysicalConstants<DType>::Rgas*u_[OFF(i,ny_-1,3)]);
	u_[OFF(i,ny_-1,1)] = u_[OFF(i,ny_-2,1)]; //Neumann
	u_[OFF(i,ny_-1,2)] = 0;
	u_[OFF(i,ny_-1,3)] = u_[OFF(i,ny_-2,3)]; //Neumann
	for(int k = 0; k < nchem; ++k)
	  u_[OFF(i,ny_-1,4+k)] = u_[OFF(i,ny_-2,4+k)]; //Neumann
      }

      //INTERIOR -- use second order approximations for first, sec. and 
      //            mixed derivatives
      
      //#pragma omp parallel for num_threads(NUMTHREADS)
      for(int i = 1; i < nx_-1; ++i){ 
	for(int j = 1; j < ny_-1; ++j){
	    //p0_ = (u_[OFF(i,j,0)]*PhysicalConstants<BaseType>::Rgas*u_[OFF(i,j,3)])/wbar(i,j,u_);

	    //! preserve mass Y_nspec = 1 - \sum_k^{nspec-1} Y_k
	  preserve_mass(N2,i,j,u_);   //N2 is excess species (bath gas)

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

	  
	    //!\f$\rho\f$
	    rhs_[OFF(i,j,0)] = -( (u_[OFF(i+1,j,0)]*u_[OFF(i+1,j,1)] - u_[OFF(i,j,0)]*u_[OFF(i,j,1)])/hx_ + (u_[OFF(i,j+1,0)]*u_[OFF(i,j+1,2)] - u_[OFF(i,j,0)]*u_[OFF(i,j,2)])/hy_ );

	    //!\f$v_1\f$
	    rhs_[OFF(i,j,1)] = -(u_[OFF(i,j,1)]*dx_v1 + u_[OFF(i,j,2)]*dy_v1 
	    			 - 1./u_[OFF(i,j,0)]*(4./3.*( (0.5*(mu(i+1,j,u_)+mu(i,j,u_))*(u_[OFF(i+1,j,1)]-u_[OFF(i,j,1)]) - 0.5*(mu(i-1,j,u_)+mu(i,j,u_))*(u_[OFF(i,j,1)]-u_[OFF(i-1,j,1)]))/ntimes<2>(hx_)) - 2./3.*(mu(i+1,j,u_)*(u_[OFF(i+1,j+1,2)] - u_[OFF(i+1,j-1,2)]) - mu(i-1,j,u_)*(u_[OFF(i-1,j+1,2)]-u_[OFF(i-1,j-1,2)]))/(4.*hx_*hy_)) - 1./u_[OFF(i,j,0)]*( (0.5*(mu(i,j+1,u_)+mu(i,j,u_))*(u_[OFF(i,j+1,1)] - u_[OFF(i,j,1)]) - 0.5*(mu(i,j-1,u_)+mu(i,j,u_))*(u_[OFF(i,j,1)] - u_[OFF(i,j-1,1)]))/ntimes<2>(hy_) + (mu(i,j+1,u_)*(u_[OFF(i+1,j+1,2)] - u_[OFF(i-1,j+1,2)]) - mu(i,j-1,u_)*(u_[OFF(i+1,j-1,2)] - u_[OFF(i-1,j-1,2)]))/(4.*hx_*hy_)) - PhysicalConstants<BaseType>::gravitational_acceleration + 1./u_[OFF(i,j,0)]*(pressure(i,j,u_) - pressure(i-1,j,u_))/hx_ ); //gravitational force added (acts only in x direction; note that domain is rotated 90° clockwise). Pressure gradient considered here
	   
	    //! \f$v_2\f$  
	    rhs_[OFF(i,j,2)] = -(u_[OFF(i,j,1)]*dx_v2 + u_[OFF(i,j,2)]*dy_v2 
				 - 1./u_[OFF(i,j,0)]*( (mu(i+1,j,u_)*(u_[OFF(i+1,j+1,1)] - u_[OFF(i+1,j-1,1)]) - mu(i-1,j,u_)*(u_[OFF(i-1,j+1,1)]-u_[OFF(i-1,j-1,1)]))/(4.*hx_*hy_) + (0.5*(mu(i+1,j,u_)+mu(i,j,u_))*(u_[OFF(i+1,j,2)]-u_[OFF(i,j,2)]) - 0.5*(mu(i-1,j,u_)+mu(i,j,u_))*(u_[OFF(i,j,2)]-u_[OFF(i-1,j,2)]))/ntimes<2>(hx_) ) -  1./u_[OFF(i,j,0)]*(4./3.*((0.5*(mu(i,j+1,u_)+mu(i,j,u_))*(u_[OFF(i,j+1,2)] - u_[OFF(i,j,2)]) - 0.5*(mu(i,j-1,u_)+mu(i,j,u_))*(u_[OFF(i,j,2)] - u_[OFF(i,j-1,2)]))/ntimes<2>(hy_)) - 2./3.*((mu(i,j+1,u_)*(u_[OFF(i+1,j+1,1)] - u_[OFF(i-1,j+1,1)]) - mu(i,j-1,u_)*(u_[OFF(i+1,j-1,1)] - u_[OFF(i-1,j-1,1)]))/(4.*hx_*hy_))) + 1./u_[OFF(i,j,0)]*(pressure(i,j,u_) - pressure(i,j-1,u_))/hy_ ); //no gravitational acceleration in y-direction
				 
	    //!TEMPERATURE \f$ T \f$
	    
#ifndef NDEBUG  
	    is_temperature_ok(i,j,0,u_);  // 0 = take bounds of first chem. spec
#endif	   
	    //project_temperature_on_its_bounds(i,j,0,u); //brings nothing

	    cp_ = cp(i,j,u_); //Cp_ has also been filled
	   
	    //!CHEMISTRY -- calculate \f$\dot{\omega}\f$
	    chemistry(dotomega_,i,j,u_); //chemical reactions

	    //! note: for low-speed flows, mechanical compression  \f$\frac{Dp}{Dt}\f$ as well viscous dissipation can be neglected [KEE, p. 115 (pdf: 142), bottom], which is realized here
	    rhs_[OFF(i,j,3)] = -(u_[OFF(i,j,1)]*dx_T + u_[OFF(i,j,2)]*dy_T 
				 - 1./(u_[OFF(i,j,0)]*cp_)*heat_conduction(i,j,u_)
				 -               //PLUS/MINUS ??
				 1./(cp_*u_[OFF(i,j,0)])*temperature_diffu_term(i,j,u_)
				 +                 //PLUS ??
				 1./(cp_*u_[OFF(i,j,0)])*heat_production(dotomega_,i,j,u_) ); 
	    
	    //! SPECIES \f$ Y_k, \qquad k = 1, \ldots, nchem\f$
	    species_diffu_term(diffspec_,i,j,u_); //now diffspec_ is filled
	    //species_diffu_term_wbar(diffspec_,i,j,u_); //now diffspec_ filled
	    for(int k = 0; k < nchem; ++k){
	      rhs_[OFF(i,j,4+k)] = -(u_[OFF(i,j,1)]*dx_spec_[k] + u_[OFF(i,j,2)]*dy_spec_[k]
                                      -              //PLUS/MINUS???
                                      1./u_[OFF(i,j,0)]*diffspec_[k]
                                      -
				      1./u_[OFF(i,j,0)]*dotomega_[k]*DataType::molar_masses()[AccessType::spec_access(k)] );
	    }
	    

	  } //END INTERIOR y
	}  //END INTERIOR x
	  
      return rhs_;
    }

    //dummy function
    void get_y_prev(typename VType::iterator y_prevIt){}

  private:
    VType rhs_;
    BuilderType BCS_;

    BaseType a_,b_,c_,d_,hx_,hy_, t0_, tend_, diam_, halfdiam_;
    int nx_, ny_, npt_, totpoints_;   
   
    BaseType p0_, v1_fuel_,v2_fuel_,v1_oxi_,v2_oxi_, T_fuel_,T_oxi_,T_wall_; //standard pressure
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

    DimensionType Vk_;   // diffusion velocity (2D)
   
    BoundaryTemperatureDistribution<BaseType> bdytemp_;
    ViscosityType Mav_;
    ConductivityType Mac_;
    DiffusionType Mad_,Mad0_,Mad1_,Mad2_,Mad3_;
    
    RatioType Theta_, Theta0_, Theta1_, Theta2_, Theta3_;
    
    BaseType yhottest_, Tignition_;

    ArrayType bdyFuel_, bdyOxi_, IV_;

    BaseType inBetween_[DataType::nspec];
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


    void preserve_mass(int select, int i, int j, VType& u){
      T sm = T();
      for(int k = 0; k < DataType::nspec; ++k){
	if(k != select){
	  sm += u[OFF(i,j,4+k)];
	}
      }
      u[OFF(i,j,4+select)] = 1. - sm;
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
    
    //!NOTE: perturbation of X is very important!!!
     T mole_fraction(int i, int j, int k, const T& Wbar, const VType& u) const{
      T mfrac = u[OFF(i,j,4+k)]*Wbar/DataType::molar_masses()[AccessType::spec_access(k)];
      perturbation(mfrac,1.e-16);
      return mfrac;
    }
    
    //! needed to evaluate diffusion and viscosity coefficients 
    ArrayType& mole_fraction(int i, int j, const T& Wbar, const VType& u){
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
    RAC& concentration(RAC& conc, int i, int j, const VType& u){
      for(int k = 0; k < nchem; ++k)
	conc[k] = concenctration(i,j,k,u,u[OFF(i,j,0)]);
      return conc;
    }

    //////////////////////////////////////////////////////////////////////////
    //SOME VECTORS ARE EXPLICITLY FILLED DURING CALL ==> NO NEED TO RE-COMP 
    //////////////////////////////////////////////////////////////////////////
     //! compute \f$ \dot{\omega} \f$, unit: mol/(m³s)
   
     ArrayType& chemistry(ArrayType& domeg, int i, int j, const VType& u){
      concentration(C_,i,j,u); //C_ is filled now
      for(int k = 0; k < nchem; ++k){ //nchem is either full dim or reddim
	//!=========TODO: in either cases: C_ must be the FULL composition
	//! compare with 
	//! '../ode/examples/automaticallygeneratedsourceterms/redh2c6_attempt.hh'
	domeg[k] = BCS_.net_formation_rate(k,u[OFF(i,j,3)],C_); //in the reduced case C_ must be of full dimensino!!!
      }
      return domeg;
    }


    
    //! \f$ \dot{\omega} \f$ is computed by now
    T heat_production(const ArrayType& dotomega,int i, int j, const VType & u){
      T h = T();
      for(int k = 0; k < nchem; ++k)
	h += BCS_.H_T(AccessType::spec_access(k),u[OFF(i,j,3)])*dotomega[k];
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
      return Mav_.compute_mixture_averaged_viscosity(u[OFF(i,j,3)], mole_fraction(i,j,wbar(i,j,u),u));
    }

    //! \f$ \lambda(\rho, p, T, X) \f$
    T lambda(int i, int j, const VType& u){
      return Mac_.compute_mixture_averaged_conductivity(u[OFF(i,j,0)], pressure(i,j,u_), u[OFF(i,j,3)],mole_fraction(i,j,wbar(i,j,u),u));
    }

    T heat_conduction(int i, int j, const VType& u){
      return ( (0.5*(lambda(i+1,j,u)+lambda(i,j,u))*(u[OFF(i+1,j,3)]-u[OFF(i,j,3)]) - 0.5*(lambda(i,j,u)+lambda(i-1,j,u))*(u[OFF(i,j,3)]-u[OFF(i-1,j,3)]))/ntimes<2>(hx_) + (0.5*(lambda(i,j+1,u)+lambda(i,j,u))*(u[OFF(i,j+1,3)]-u[OFF(i,j,3)]) - 0.5*(lambda(i,j,u)+lambda(i,j-1,u))*(u[OFF(i,j,3)]-u[OFF(i,j-1,3)]))/ntimes<2>(hy_) );
    }


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

      Mad0_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j-1,u),u[OFF(i,j-1,3)],X0_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(pressure(i+1,j,u),u[OFF(i+1,j,3)],X1_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j+1,u),u[OFF(i,j+1,3)],X2_);
      
      Mad3_.compute_mixture_averaged_diffusion_coefficients(pressure(i-1,j,u),u[OFF(i-1,j,3)],X3_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j,u),u[OFF(i,j,3)],X_);

#ifdef CORRECTION_DIFFUSION_VELOCITY_WITH_THERMAL_PART
      //! thermal diffusion ratios
      Theta0_.compute_thermal_diffusion_ratios(u[OFF(i,j-1,3)],X0_);
      Theta1_.compute_thermal_diffusion_ratios(u[OFF(i+1,j,3)],X1_);
      Theta2_.compute_thermal_diffusion_ratios(u[OFF(i,j+1,3)],X2_);
      Theta3_.compute_thermal_diffusion_ratios(u[OFF(i-1,j,3)],X3_);
      Theta_.compute_thermal_diffusion_ratios(u[OFF(i,j,3)],X_);
#endif

      for(int k = 0; k < nchem; ++k){
      	//!x-direction
      	v[0] += ( -0.5*(DataType::molar_masses()[k]*Mad1_[k]/fivepointstencil_[1] + DataType::molar_masses()[k]*Mad_[k]/fivepointstencil_[4])*(X1_[k]-X_[0])/hx_  //approximation of \f$ Y_k\mathcal{V}_k\f$ in \f$x\f$-direction
#ifdef CORRECTION_DIFFUSION_VELOCITY_WITH_THERMAL_PART
		  -0.5*(DataType::molar_masses()[k]*Mad1_[k]/fivepointstencil_[1]*Theta1_[k]/u[OFF(i+1,j,3)] + DataType::molar_masses()[k]*Mad_[k]/fivepointstencil_[4]*Theta_[k]/u[OFF(i,j,3)])*(u[OFF(i+1,j,3)]-u[OFF(i,j,3)])/hx_  //approximation of \f$ Y_k\mathcal{W}_k\f$ in \f$x\f$-direction
#endif
      			  );  

      	//!y-direction
      	v[1] += ( -0.5*(DataType::molar_masses()[k]*Mad2_[k]/fivepointstencil_[2] + DataType::molar_masses()[k]*Mad_[k]/fivepointstencil_[4])*(X2_[k]-X_[0])/hy_  //approximation of \f$ Y_k\mathcal{V}_k\f$ in \f$y\f$-direction
#ifdef CORRECTION_DIFFUSION_VELOCITY_WITH_THERMAL_PART
		  -0.5*(DataType::molar_masses()[k]*Mad2_[k]/fivepointstencil_[2]*Theta2_[k]/u[OFF(i,j+1,3)] + DataType::molar_masses()[k]*Mad_[k]/fivepointstencil_[4]*Theta_[k]/u[OFF(i,j,3)])*(u[OFF(i,j+1,3)]-u[OFF(i,j,3)])/hy_  //approximation of \f$ Y_k\mathcal{W}_k\f$ in \f$y\f$-direction
#endif
		  );    
      }
#endif //in case of no definition of 'NO_CORRECTION_DIFFUSION_VELOCITY'
      return v;
      
    }
   
   

    //perhaps much better to replace \f$ Y_k/X_k\f$ by \f$W_k/\overline{W}\f$ since \f$X_k\f$ can be zero
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

      Mad0_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j-1,u),u[OFF(i,j-1,3)],X0_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(pressure(i+1,j,u),u[OFF(i+1,j,3)],X1_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j+1,u),u[OFF(i,j+1,3)],X2_);
      
      Mad3_.compute_mixture_averaged_diffusion_coefficients(pressure(i-1,j,u),u[OFF(i-1,j,3)],X3_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j,u),u[OFF(i,j,3)],X_);

      
      //! thermal diffusion ratios
      Theta0_.compute_thermal_diffusion_ratios(u[OFF(i,j-1,3)],X0_);
      Theta1_.compute_thermal_diffusion_ratios(u[OFF(i+1,j,3)],X1_);
      Theta2_.compute_thermal_diffusion_ratios(u[OFF(i,j+1,3)],X2_);
      Theta3_.compute_thermal_diffusion_ratios(u[OFF(i-1,j,3)],X3_);
      Theta_.compute_thermal_diffusion_ratios(u[OFF(i,j,3)],X_);

      //!note V_corr is a constant correction vector, i.e. independent of species but varying in space and time
      V_corr(Vk_,i,j,u);       //correction
      
     
      for(int k = 0; k < nchem; ++k){
	diffu[k] =( (0.5*(u[OFF(i+1,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[1]*Mad1_[k] + u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k])*(X1_[k]-X_[k]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k] + u[OFF(i-1,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[3]*Mad3_[k])*(X_[k]-X3_[k]))/ntimes<2>(hx_) + 
	  (0.5*(u[OFF(i,j+1,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[2]*Mad2_[k] + u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k])*(X2_[k]-X_[k]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k] + u[OFF(i,j-1,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[0]*Mad0_[k])*(X_[k]-X0_[k]))/ntimes<2>(hy_)
	  //now comes the part concerning thermal diffusion...
	  + (0.5*(u[OFF(i+1,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[1]*Mad1_[k]*Theta1_[k]/u[OFF(i+1,j,3)] + u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k]*Theta_[k]/u[OFF(i,j,3)])*(u[OFF(i+1,j,3)] - u[OFF(i,j,3)]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k]*Theta_[k]/u[OFF(i,j,3)] + u[OFF(i-1,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[3]*Mad3_[k]*Theta3_[k]/u[OFF(i-1,j,3)])*(u[OFF(i,j,3)]-u[OFF(i-1,j,3)]))/ntimes<2>(hx_) //x-part
	  + (0.5*(u[OFF(i,j+1,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[2]*Mad2_[k]*Theta2_[k]/u[OFF(i,j+1,3)] + u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k]*Theta_[k]/u[OFF(i,j,3)])*(u[OFF(i,j+1,3)] - u[OFF(i,j,3)]) - 0.5*(u[OFF(i,j,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4]*Mad_[k]*Theta_[k]/u[OFF(i,j,3)] + u[OFF(i,j-1,0)]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[0]*Mad0_[k]*Theta0_[k]/u[OFF(i,j-1,3)])*(u[OFF(i,j,3)]-u[OFF(i,j-1,3)]))/ntimes<2>(hy_)  //y-part 
#ifndef NO_CORRECTION_DIFFUSION_VELOCITY
	  //...and the correcting term; note: \f$ Y_k = X_k W_k/\overline{W}\f$
	  //... as well as a constant Vk_
	  + 
	  ( Vk_[0]*(u[OFF(i,j,0)]*X_[k]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4] - u[OFF(i-1,j,0)]*X3_[k]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[3])/hx_
	    + Vk_[1]*(u[OFF(i,j,0)]*X_[k]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[4] - u[OFF(i,j-1,0)]*X0_[k]*DataType::molar_masses()[AccessType::spec_access(k)]/fivepointstencil_[0])/hy_ )
#endif
		    );

      }
      
      return diffu;
    }


     //! assumes that cp(i,j,u) has been invoked previously
     T temperature_diffu_term(int i, int j, const VType& u){
       T td = T(),
	 meanmolmass = wbar(i,j,u);

       X_ = mole_fraction(i,j,meanmolmass,u); //(i,j)
       Mad_.compute_mixture_averaged_diffusion_coefficients(pressure(i,j,u_),u[OFF(i,j,3)],X_); //diffusion coefficients are computed by now
       
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
	Theta0_.compute_thermal_diffusion_ratios(u[OFF(i,j-1,3)],X0_);


       //isDiffusioncalculated_ = true;
       
       for(int k = 0; k < nchem; ++k){
   
	//! assume Cp_ has been filled previously
	 td += ( u[OFF(i,j,0)]*Cp_[k]*( Mad_[k]/meanmolmass*( (0.5*(X1_[k]+X_[k])*(u[OFF(i+1,j,3)]-u[OFF(i,j,3)]) - 0.5*(X3_[k]+X_[k])*(u[OFF(i,j,3)]-u[OFF(i-1,j,3)]))/ntimes<2>(hx_) + (0.5*(X2_[k]+X_[k])*(u[OFF(i,j+1,3)]-u[OFF(i,j,3)]) - 0.5*(X_[k]+X0_[k])*(u[OFF(i,j,3)]-u[OFF(i,j-1,3)]))/ntimes<2>(hy_) 
       //thermal diffusion part...
	+ 						     
							      Theta_[k]/u[OFF(i,j,3)]*( (0.5*(u[OFF(i+1,j,3)]+u[OFF(i,j,3)])*(u[OFF(i+1,j,3)]-u[OFF(i,j,3)]) - 0.5*(u[OFF(i,j,3)]+u[OFF(i-1,j,3)])*(u[OFF(i,j,3)]-u[OFF(i-1,j,3)]) )/ntimes<2>(hx_) + (0.5*(u[OFF(i,j+1,3)]+u[OFF(i,j,3)])*(u[OFF(i,j+1,3)]-u[OFF(i,j,3)]) - 0.5*(u[OFF(i,j,3)]+u[OFF(i,j-1,3)])*(u[OFF(i,j,3)]-u[OFF(i,j-1,3)]) )/ntimes<2>(hy_)) ) //!multiplication with \f$D_k^{\mathrm{mix}}/\overline{W}\f$
 //corrective part; note: \f$Y_k/W_k = X_k/\overline{W} \f$
#ifndef NO_CORRECTION_DIFFUSION_VELOCITY
					+ X_[k]/meanmolmass*(Vk_[0]*(u[OFF(i,j,3)] - u[OFF(i-1,j,3)])/hx_ + Vk_[1]*(u[OFF(i,j,3)] - u[OFF(i,j-1,3)])/hy_)						     //multiplication with \f$Y_k/W_k\f$
#endif
					)//multiplic. with \f$\rho C_{pk}^0\f$ 
		 );
      }
      return td;
     }
  


    T pressure(int i, int j, const VType& u){
      return ( (u[OFF(i,j,0)]*PhysicalConstants<BaseType>::Rgas*u[OFF(i,j,3)])/wbar(i,j,u) );
    }

  };

} //end namespace

#endif 
