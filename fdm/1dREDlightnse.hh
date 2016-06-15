#ifndef REDUCED_1D_SPATIAL_VERSION_OF_LIGHT_NSE_HH
#define REDUCED_1D_SPATIAL_VERSION_OF_LIGHT_NSE_HH

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
  class ReducedOneDIMLightNSE{
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
    typedef IndexHandler<DataType,true> AccessType; //===RED
    typedef typename AccessType::IndexerType IndexerType;

    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,IndexerType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    typedef EquationOfState<'i'> EosType;

    //TRANSPORT COEFFICIENTS
    //NO VISCOSITY NEEDED IN THIS CASE -- +++++FULL+++++ version here anyway
    typedef MixtureAveragedConductivity<DataType,false> ConductivityType;
    typedef MixtureAveragedDiffusion<DataType,false> DiffusionType;
    


  //====RED: with non-Cppad type since we solve system via LAPACK and integrate
  typedef ReconstructSpecies<DType,DataType::nspec,DataType::rednspec,ChooseFDSetting::MECH> ReconstructorType;

    ReducedOneDIMLightNSE(SizeType dim = 0):rhs_(dim){
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
    
      nx_ = PD.get_datum<int>("Nx");
      npt_ = nx_; //all points, i.e. \f$\Omega_h \cup \partial\Omega_h\f$
      totpoints_ = nprim*npt_;
      hx_ = (b_-a_)/(nx_-1);
     
      p0_ = PD.get_datum<BaseType>("pconst");
      v1_ = PD.get_datum<BaseType>("v1_in");
     
      Tinlet_ = PD.get_datum<BaseType>("T_in");
      T_burner_ = PD.get_datum<BaseType>("T_ignition");


      //transport coefficients are ready for computation
      //!no viscosity needed
      Mac_.initialize(DataType::transport());
      Mad_.initialize(DataType::transport());  //(i)  //3 pt. stencil
      Mad0_.initialize(DataType::transport()); //(i-1)
     
      Mad2_.initialize(DataType::transport()); //(i+1)
     

      bdy2d_.initialize(DataType::nspec); //alway in the full dimension!
      
     
      u_.resize(totpoints_);

      cp_ = rho_ = T();  //initialize

      //===RED=================================
      Tab_.resize(nx_);

 
      for(int i = 0; i < nx_; ++i){
	(Tab_[i]).resize(DataType::nspec+1);

	//for(int k = 0; k < DataType::nspec; ++k){
	  //std::cout << bdy2d_[k] << ", ";
	Tab_[i][0] = 0.;
	Tab_[i][1] = 0.8;
	Tab_[i][2] = 0.2; //= bdy2d_[k];
	//}
	//std::cout << std::endl;
	Tab_[i][DataType::nspec] = Tinlet_;  //NOTE: here last entry is T!!
	if(i == 0)
	  Tab_[i][DataType::nspec] = T_burner_;
      }

      //! INITIALIZE species Reconstructor
      //!some guess
      VDType guessSpecPlusTemp(DataType::nspec+1);
      //for()
      smart_assign(guessSpecPlusTemp[0],0.);
      smart_assign(guessSpecPlusTemp[1],0.8);
      smart_assign(guessSpecPlusTemp[2],0.2);
      guessSpecPlusTemp[DataType::nspec] = T_burner_; //temperature
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
      
      //!BOUNDARY CONDITIONS -- see Chemkin theory, §12.6, p. 209
      //===RED
      reconstruct_full_state(zM_,0,u_);
      reconstruct_full_state(zM_,nx_-1,u_);
	  
      //······ LEFT ··········
      u_[OFF(0,TMP)] = T_burner_;  
      
     
      for(int k = 0; k < nchem; ++k){
	u_[OFF(0,SPC+k)] = 0.; //bdy2d_[AccessType::spec_access(k)];  //===TODO
      } 
      //··········· RIGHT: Neumann ························
      u_[OFF(nx_-1,TMP)] = u_[OFF(nx_-2,TMP)];
       for(int k = 0; k < nchem; ++k){
	u_[OFF(nx_-1,SPC+k)] = u_[OFF(nx_-2,SPC+k)];  //===TODO
      }
      
// #ifndef NDEBUG       
//        //! check if everything is o.k.
//        check_4_nans(0,u_[OFF(0,TMP)],"T on dOmega_left");
//        check_4_nans(0,u_[OFF(nx_-1,TMP)],"T on dOmega_right");
	 
//        for(int k = 0; k < nchem; ++k){
// 	 check_4_nans(0,u_[OFF(0,SPC+k)],"Spec "+Num2str(SPC+k)+" on dOmega_left");
// 	 check_4_nans(nx_-1,u_[OFF(0,SPC+k)],"Spec "+Num2str(SPC+k)+"on dOmega_right");
//        }
// #endif	 

      //!INTERIOR
      for(int i = 1; i < nx_-1; ++i){
	//===RED
	reconstruct_full_state(zM_,i,u_);
	//std::cout << "zM_ (INTERIOR) = "<< zM_ << std::endl;

#ifndef NDEBUG  
	  is_temperature_ok(i,0,u_);  // 0 = take bounds of first chem. spec
	  
	  check_4_nans(u_);
#endif	   

	  //preserve_mass(i,u_);
	T dx_T = (u_[OFF(i,TMP)]-u_[OFF(i-1,TMP)])/hx_; //always >= 0
	for(int k = 0; k < nchem; ++k){
	  dx_spec_[k] = (u_[OFF(i,SPC+k)]-u_[OFF(i-1,SPC+k)])/hx_;
	}

	rho_ = rho(i,u_); 
	cp_ = cp(i,u_);
	adonis_assert(rho_ > T() && cp_ > T());
	  
	chemistry(dotomega_,i,u_); //chemical reactions

	//temperature
	rhs_[OFF(i,TMP)] = -(v1_*dx_T 
			     - 1./(rho_*cp_)*heat_conduction(i,u_)
			     - 1./(rho_*cp_)*temperature_diffu_term(i,u_)
			     + 1./(rho_*cp_)*heat_production(i,u_)
			     );
      
	//!species --- +++FULL++++
	species_diffu_term(DiffFull_,i,u_); //now DiffFull_ is filled
	for(int k = 0; k < nchem; ++k){ //===RED (only)!!!
	  rhs_[OFF(i,SPC+k)] = -( v1_*dx_spec_[k]  
				  -1./rho_*DiffFull_[AccessType::spec_access(k)]
				  -1./rho_*dotomega_[k]*DataType::molar_masses()[AccessType::spec_access(k)]
				    );
	}
      } //end interior i
      return rhs_;
    }

  private:
    VType rhs_;
    BuilderType BCS_;

    BaseType a_,b_,hx_;
    int nx_, npt_, totpoints_;   
   
  BaseType p0_, v1_, Tinlet_, T_burner_; //standard pressure
   
    VType veloprofile_,
      u_;
    T cp_, //for storage
      rho_;  //rho will be evaluated by an EOS
    

    T threepointstencil_[3];
    T star3_[3];

    T dx_spec_[nchem];
    T dy_spec_[nchem];

   
    ConductivityType Mac_;
    DiffusionType Mad_,Mad0_,Mad2_;
   
    BoundaryType bdy2d_; 
   
  //===RED
    TabulationType Tab_;
    ReconstructorType Recon_;
    VType zM_;
    FullChemType Conc_, Cp_;
    RedChemType redspecies_;
    VDType Ctilde_,btilde_;
    VDType low_,up_;
    RedChemType dotomega_;
    FullChemType X_, X0_,X2_, DiffFull_;   //++++FULL++++

    VType& reconstruct_full_state(VType& zm, int i, const VType& u){
      for(int k = 0; k < DataType::rednspec; ++k)
	redspecies_[k] = u[OFF(i,SPC+k)]; //Y_k, reduced composition
      std::cout << "reduced species at point i = "<< i <<":  " << redspecies_ << std::endl;
      Recon_.assign(redspecies_,Tab_[i]);
      Recon_.enforce_linearized_constraints_wrt_chemistry(Ctilde_,btilde_);
      Recon_.project_onto_bounds(low_,up_); //Y+T
      zm = Recon_.time_step(Recon_.get_z(),p0_); //Y+T 
       
      Tab_[i] = zm;  //overwrite with new entry
      return zm;
    }

    //PRIVATE FUNCTIONS -- don't use them from outside

    int OFF(int i, int spec) const{
      return (i + spec*npt_);
    }

    

    //!PHYSICAL CONSTRAINTS
    // void preserve_mass(int i, VType& u){
    //   T sm1 = T();
    //   for(int k = 0; k < nchem-1; ++k)
    // 	sm1 += u[OFF(i,SPC+k)];
      
    //   T Y_n = 1.-sm1;
    //   u[OFF(i,SPC+(nchem-1))] = Y_n;
    //   adonis_assert(Y_n >= T());  //if not, something went wrong
      
    // }

    //====RED: revise functions for model reduction 
    //===      assume reconstruct_full_state(zM_,i,u_) has been invoked before
    //===      some full chem. state evaluations are to be performed
    //! CHEMISTRY
    //in the reduced case  only the corresponding fixed quantities
    //are extracted, e.g. molar masses.
    //FULL DIMENSION -- in most cases, u is just a dummy variable 
    T wbar(int i, const VType& u) const{ 
      T wb = T();
      for(int k = 0; k < DataType::nspec; ++k){  //FULL DIM
	wb += Tab_[i][k]/DataType::molar_masses()[k];
      }
      adonis_assert(wb > 0.);
      
      return 1./wb;
    }

    T rho(int i, const VType& u){
      return (p0_*wbar(i,u))/(PhysicalConstants<BaseType>::Rgas*Tab_[i][DataType::nspec]); //full
    }

    //+++FULL+++
    //MIGHT BE ZERO
   T mole_fraction(int i, int k, const T& Wbar, const VType& u) const{
     return (Tab_[i][k]*Wbar)/DataType::molar_masses()[k];
    }

   FullChemType& mole_fraction(FullChemType& x, int i, const VType u){
     for(int k = 0; k < DataType::nspec; ++k){
       x[k] = (Tab_[i][k]*wbar(i,u))/DataType::molar_masses()[k];
     }
     return x;
   }
   
   //if Wbar is reused somewhere
    FullChemType& mole_fraction(FullChemType& x, int i,  const T& Wbar, const VType u){
      for(int k = 0; k < DataType::nspec; ++k){
	x[k] = (Tab_[i][k]*Wbar)/DataType::molar_masses()[k];
      }
      return x;
    }

    //! FULL DIMENSION -- needed for chemistry
   FullChemType& concentration(FullChemType& conc, int i, const VType u){
     T Rho =  rho(i,u); 
     for(int k = 0; k < DataType::nspec; ++k){
       conc[k] = (Rho*Tab_[i][k])/DataType::molar_masses()[k];
     }
     return conc;
   }

   RedChemType& chemistry(RedChemType& domeg, int i, const VType& u){
     concentration(Conc_,i,u); //Conc_ is filled now
     for(int k = 0; k < nchem; ++k){  //reduced species
       //! ==== TODO: construction of the chem. source term should 
       //! ========== index correctly, right?????????????
       domeg[k] = BCS_.net_formation_rate(BCS_.species_index(k),Tab_[i][DataType::nspec],Conc_);
     }
     return domeg;
    }

    // FULL DIMENSION -- optional
   T heat_production(int i, const VType & u){
      T h = T();  
      // for(int k = 0; k < DataType::nspec; ++k)
      // 	h += BCS_.H_T(k,Tab_[i][DataType::nspec])*Recon_.dot_omega()[k];  //dotomega[k];
     
      int nurep = DataType::nspec-DataType::rednspec;
      for(int k = 0; k < DataType::rednspec; ++k)
	h += BCS_.H_T(AccessType::spec_access(k),Tab_[i][DataType::nspec])*dotomega_[k];
      for(int k = 0; k < nurep; ++k)
	h += BCS_.H_T(Recon_.unrep_index()[k],Tab_[i][DataType::nspec])*Recon_.dot_omega()[Recon_.unrep_index()[k]];

      return h;
    }

    //!FULL DIMENSION
    T cp(int i,  const VType& u){
      T cres = T();
      for(int k = 0; k < DataType::nspec; ++k){ //fill also Cp_
	Cp_[k] = BCS_.C_p(k,Tab_[i][DataType::nspec]);
	cres += (Cp_[k]/DataType::molar_masses()[k])*Tab_[i][k];
      }
      //isCpcalculated_ = true;
      return cres;
    }


    //!======== TODO: evaluate in Tab_[i][DataType::nspec] or u ??
    //! +++ FULL +++
    T lambda(int i, const VType& u){
      (*this).mole_fraction(X_,i,wbar(i,u),u);
      //perturb(X_,1.e-16);
      perturbation(X_,1.e-16);
      return Mac_.compute_mixture_averaged_conductivity((*this).rho(i,u), p0_, Tab_[i][DataType::nspec],X_); //X_ overwritten
    }
  
   T heat_conduction(int i, const VType& u){
    
     return ( (0.5*(lambda(i+1,u)+lambda(i,u))*(Tab_[i+1][DataType::nspec]-Tab_[i][DataType::nspec]) - 0.5*(lambda(i,u)+lambda(i-1,u))*(Tab_[i][DataType::nspec]-Tab_[i-1][DataType::nspec]))/ntimes<2>(hx_)  
	      );
    }
   
   //cp() must be invoked beforehand
     T temperature_diffu_term(int i, const VType& u){
       T td = T(),
	 meanmolmass = wbar(i,u),
	 Rho = rho(i,u);

       mole_fraction(X0_,i-1,wbar(i-1,u),u);
       mole_fraction(X_,i,meanmolmass,u);
       mole_fraction(X2_,i+1,wbar(i+1,u),u);

       //===RED
       perturbation(X0_,1.e-16);
       perturbation(X_,1.e-16);
       perturbation(X2_,1.e-16);

       Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[i][DataType::nspec],X_); //diffusion coefficients are computed by now
     

       //isDiffusioncalculated_ = true;
       
       //++++ 
       for(int k = 0; k < DataType::nspec; ++k){
	 //std::cout << "++++++ Self-diffusion coefficient = "<< Mad_.binary_diffusion_coefficients(k,k) << std::endl;
	 //!store counterclockwise, starting from the bottom with
	//!(i) being the last entry
	 threepointstencil_[0] = X0_[k];//mole_fraction(i-1,k,wbar(i-1,u),u);
	 threepointstencil_[1] = X_[k];   //mole_fraction(i,k,wbar(i,u),u);
	 threepointstencil_[2] = X2_[k];   //mole_fraction(i+1,k,wbar(i+1,u),u);

	 //!===RED
	//! assume Cp_ has been filled previously
	 td += ( (Cp_[k]/meanmolmass)*Mad_.binary_diffusion_coefficients(k,k)*Rho*( (0.5*(threepointstencil_[2]+threepointstencil_[1])*(Tab_[i+1][DataType::nspec]-Tab_[i][DataType::nspec]) - 0.5*(threepointstencil_[1]+threepointstencil_[0])*(Tab_[i][DataType::nspec]-Tab_[i-1][DataType::nspec]))/ntimes<2>(hx_) ) );
	
	 //std::cout << "td = "<< td << std::endl;
      }
      return td;
     }
  
    //+++FULL++++
   FullChemType& species_diffu_term(FullChemType& diffu, int i, const VType& u){
      threepointstencil_[0] = wbar(i-1,u);
      threepointstencil_[1] = wbar(i,u);
      threepointstencil_[2] = wbar(i+1,u);
      
      mole_fraction(X0_,i-1,threepointstencil_[0],u);  //(i-1)
      mole_fraction(X_, i,  threepointstencil_[1],u);  //(i)
      mole_fraction(X2_,i+1,threepointstencil_[2],u);  //(i+1)
      
      //===RED
      perturbation(X0_,1.e-16);
      perturbation(X_,1.e-16);
      perturbation(X2_,1.e-16);


      Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[i-1][DataType::nspec],X0_);
    
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[i][DataType::nspec],X_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,Tab_[i+1][DataType::nspec],X2_);
      

      //CAUTION: use self-diffusion coefficient for single gas species
      for(int k = 0; k < DataType::nspec; ++k){ //+++FULL+++
	diffu[k] = (0.5*(rho(i+1,u)*DataType::molar_masses()[k]/threepointstencil_[2]*Mad2_.binary_diffusion_coefficients(k,k) + rho(i,u)*DataType::molar_masses()[k]/threepointstencil_[1]*Mad_.binary_diffusion_coefficients(k,k))*(X2_[k]-X_[k]) - 0.5*(rho(i,u)*DataType::molar_masses()[k]/threepointstencil_[1]*Mad_.binary_diffusion_coefficients(k,k) + rho(i-1,u)*DataType::molar_masses()[k]/threepointstencil_[0]*Mad0_.binary_diffusion_coefficients(k,k))*(X_[k]-X0_[k]))/ntimes<2>(hx_);

      }
      
      return diffu;
    }

   
    //with wbar
     // ArrayType& species_diffu_term_wbar(ArrayType& diffu, int i, const VType& u){
       
     //   threepointstencil_[0] = wbar(i-1,u);
     //   threepointstencil_[1] = wbar(i,u);
     //   threepointstencil_[2] = wbar(i+1,u);
      
     //   mole_fraction(X0_,i-1,threepointstencil_[0],u);  //(i-1)
     //   mole_fraction(X_, i,  threepointstencil_[1],u);  //(i)
     //   mole_fraction(X2_,i+1,threepointstencil_[2],u);  //(i+1)
      
     //   Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i-1,TMP)],X0_);
    
     //   Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,TMP)],X_);
     //   Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i+1,TMP)],X2_);
      
       
     //   star3_[0] = rho(i-1,u);
     //   star3_[1] = rho(i,u);
     //   star3_[2] = rho(i+1,u);
      
     //   for(int k = 0; k < nchem; ++k){
     // 	//!\f$\nabla \cdot (\rho D^{\mathrm{mix}}_k \nabla Y_k)\f$
     // 	 diffu[k] = ( (0.5*(star3_[2]*Mad2_[k] + star3_[1]*Mad_[k])*(u[OFF(i+1,SPC+k)]-u[OFF(i,SPC+k)]) - 0.5*(star3_[1]*Mad_[k] + star3_[0]*Mad0_[k])*(u[OFF(i,SPC+k)]-u[OFF(i-1,SPC+k)]))/ntimes<2>(hx_) 
	
     // 	  //!differential part with mean molecular weight \f$ \overline{W}\f$,
     // 	  //!i.e. \f$ \nabla \cdot (\rho Y_k/\overline{W}D^{\mathrm{mix}}_k \nabla \overline{W})\f$
     // 		      +  (0.5*(star3_[2]*u[OFF(i+1,SPC+k)]/threepointstencil_[2]*Mad2_[k] + star3_[1]*u[OFF(i,SPC+k)]/threepointstencil_[1]*Mad_[k])*(threepointstencil_[2]-threepointstencil_[1]) - 0.5*(star3_[1]*u[OFF(i,SPC+k)]/threepointstencil_[1]*Mad_[k] + star3_[0]*u[OFF(i-1,SPC+k)]/threepointstencil_[0]*Mad0_[k])*(threepointstencil_[1]-threepointstencil_[0]))/ntimes<2>(hx_) );
	
     //   }
     //   return diffu;
     // }

    void is_temperature_ok(int i,int k, const VType& u){
      std::string oper1;
      for(int l = -1; l < 2; ++l){
	int IX = i+l;
	if(l>=0)
	  oper1 = "+";
	
	if(u[OFF(IX,TMP)] < DataType::temperature_bounds()[3*k] || 
	     u[OFF(IX,TMP)] > DataType::temperature_bounds()[3*k+1]){
	   
	  ADONIS_ERROR(BoundsError,"Temperature bounds violated: T_{i"<<oper1<<l  << "} = "<< u[OFF(IX,TMP)]<<", i.e. T("<<IX<<").");
	  }
      }
    }

    void check_4_nans(int i, const T& val, const std::string& s){
      if(std::isnan(convert_fancy_num_2_num(val)))
	ADONIS_ERROR(ComputationError,"NaN occurred at point i = "<< i << s << ".");
	
    }

    void check_4_nans(const VType& u){
      for(int i = 0; i < nx_; ++i){
	if(std::isnan(convert_fancy_num_2_num(u[OFF(i,TMP)])))
	  ADONIS_ERROR(ComputationError,"NaN occurred at point i = "<< i << " (TEMPERATURE).");
	for(int k = 0; k < nchem; ++k)
	  if(std::isnan(convert_fancy_num_2_num(u[OFF(i,SPC+k)])))
	    ADONIS_ERROR(ComputationError,"NaN occurred at point i = "<< i << " (red. SPECIES).");
      }

    }
  };


} //end namespace 


#endif
