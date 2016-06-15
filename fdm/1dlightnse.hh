#ifndef ONE_DIMENSIONAL_VERSION_OF_LIGHT_NSE_HH
#define ONE_DIMENSIONAL_VERSION_OF_LIGHT_NSE_HH

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

namespace Adonis{

template<class T>
  class OneDIMLightNSE{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;

    typedef typename TypeAdapter<T>::BaseType BaseType; //for norms, tolerances 
    typedef std::size_t SizeType;
    
     //needed for thermo chemistry
    typedef ThermoData4Mechanism<BaseType,ChooseFDSetting::MECH> DataType;  

    typedef Boundary2D<VType,ChooseFDSetting::MECH> BoundaryType;
    
    enum{nchem = SystemDimension<DataType::nspec,DataType::rednspec,ChooseFDSetting::ISREDUCED>::Value, 
	 //T + specs
	 nprim = 1+SystemDimension<DataType::nspec,DataType::rednspec,ChooseFDSetting::ISREDUCED>::Value,  //1+nchem
	 TMP = 0,          //first entry is temperature
	 SPC = 1          //0 is temperature, from here species start
    };
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
    //NO VISCOSITY NEEDED IN THIS CASE
    typedef MixtureAveragedConductivity<DataType,ChooseFDSetting::ISREDUCED> ConductivityType;
    typedef MixtureAveragedDiffusion<DataType,ChooseFDSetting::ISREDUCED> DiffusionType;

    // this may be of some importance
    typedef ComputeTransportProperties<ChooseFDSetting::ISREDUCED> PropType;


    OneDIMLightNSE(SizeType dim = 0):rhs_(dim){
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
      Mad_.initialize(DataType::transport());  //(i)
      Mad0_.initialize(DataType::transport()); //(i-1)
     
      Mad2_.initialize(DataType::transport()); //(i+1)
     

      bdy2d_.initialize(DataType::nspec); //alway in the full dimension!
      
     
      u_.resize(totpoints_);

      cp_ = rho_ = T();  //initialize

     
    }

    int dim() const {return totpoints_;}
    int domain_dim() const {return totpoints_;}
    std::string name() const {return "WITHOUT rho and momentum";}

    //============================
    //evaluate vector field
    //============================
    template<class VEC>
    VType& operator()(VEC& vars){ //reference here for we wanna change it 
      //HARD COPY
      u_ = vars;
      
      //!BOUNDARY CONDITIONS -- see Chemkin theory, §12.6, p. 209
      //······ LEFT ··········
      u_[OFF(0,TMP)] = T_burner_;  
      
      mole_fraction(X0_,0,u_);
      mole_fraction(X_,1,u_);
      Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,u_[OFF(0,TMP)],X0_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,u_[OFF(1,TMP)],X_);
      for(int k = 0; k < nchem; ++k){
	u_[OFF(0,SPC+k)] = bdy2d_[AccessType::spec_access(k)];  //===TODO
	  //	- (0.5*(rho(0,u_)*DataType::molar_masses()[AccessType::spec_access(k)]/wbar(0,u_)*Mad0_[k]/(rho(0,u_)*v1_) + rho(1,u_)*DataType::molar_masses()[AccessType::spec_access(k)]/wbar(1,u_)*Mad0_[k]/(rho(1,u_)*v1_))  );
      } 
      //··········· RIGHT: Neumann ························
      u_[OFF(nx_-1,TMP)] = u_[OFF(nx_-2,TMP)];
       for(int k = 0; k < nchem; ++k){
	u_[OFF(nx_-1,SPC+k)] = u_[OFF(nx_-2,SPC+k)];  //===TODO
      }
      

	  

      //!INTERIOR
      for(int i = 1; i < nx_-1; ++i){

#ifndef NDEBUG  
	  is_temperature_ok(i,0,u_);  // 0 = take bounds of first chem. spec
#endif	   

	preserve_mass(i,u_);
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
			     + 1./(rho_*cp_)*heat_production(dotomega_,i,u_)
			     );
      
	//!species
	species_diffu_term(diffspec_,i,u_); //now diffspec_ is filled
	//species_diffu_term_wbar(diffspec_,i,u_); //now diffspec_ is filled
	for(int k = 0; k < nchem; ++k){
	  rhs_[OFF(i,SPC+k)] = -( v1_*dx_spec_[k]  
				    -1./rho_*diffspec_[k]
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
    ArrayType X_, X0_,X2_,diffspec_, C_, dotomega_, Cp_;
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
   
    //PRIVATE FUNCTIONS -- don't use them from outside

    int OFF(int i, int spec) const{
      return (i + spec*npt_);
    }

    //!PHYSICAL CONSTRAINTS
    void preserve_mass(int i, VType& u){
      T sm1 = T();
      for(int k = 0; k < nchem-1; ++k)
	sm1 += u[OFF(i,SPC+k)];
      
      T Y_n = 1.-sm1;
      u[OFF(i,SPC+(nchem-1))] = Y_n;
      adonis_assert(Y_n >= T());  //if not, something went wrong
      
    }

   //! CHEMISTRY
    //in the reduced case  only the corresponding fixed quantities
    //are extracted, e.g. molar masses.
    T wbar(int i, const VType& u) const{
      T wb = T();
      for(int k = 0; k < nchem; ++k){
	wb += u[OFF(i,SPC+k)]/DataType::molar_masses()[AccessType::spec_access(k)];
      }
      adonis_assert(wb > 0.);
      
      return 1./wb;
    }

    T rho(int i, const VType& u){
      return (p0_*wbar(i,u))/(PhysicalConstants<BaseType>::Rgas*u[OFF(i,TMP)]);
    }


   T mole_fraction(int i, int k, const T& Wbar, const VType& u) const{
      return (u[OFF(i,SPC+k)]*Wbar)/DataType::molar_masses()[AccessType::spec_access(k)];
    }

   ArrayType& mole_fraction(ArrayType& x, int i, const VType u){
     for(int k = 0; k < nchem; ++k){
       x[k] = (u[OFF(i,SPC+k)]*wbar(i,u))/DataType::molar_masses()[AccessType::spec_access(k)];
     }
     return x;
   }
   
   //if Wbar is reused somewhere
    ArrayType& mole_fraction(ArrayType& x, int i,  const T& Wbar, const VType u){
      for(int k = 0; k < nchem; ++k){
	x[k] = (u[OFF(i,SPC+k)]*Wbar)/DataType::molar_masses()[AccessType::spec_access(k)];
      }
      return x;
    }

   ArrayType& concentration(ArrayType& conc, int i, const VType u){
     T Rho =  rho(i,u); 
     for(int k = 0; k < nchem; ++k){
       conc[k] = (Rho*u[OFF(i,SPC+k)])/DataType::molar_masses()[AccessType::spec_access(k)];
     }
     return conc;
   }

   ArrayType& chemistry(ArrayType& domeg, int i, const VType& u){
     concentration(C_,i,u); //C_ is filled now
      for(int k = 0; k < nchem; ++k){
	//! ==== TODO: construction of the chem. source term should 
	//! ========== index correctly, right?????????????
	domeg[k] = BCS_.net_formation_rate(k,u[OFF(i,TMP)],C_);
      }
      return domeg;
    }

   T heat_production(const ArrayType& dotomega, int i, const VType & u){
      T h = T();
       for(int k = 0; k < nchem; ++k)
	h += BCS_.H_T(AccessType::spec_access(k),u[OFF(i,TMP)])*dotomega[k];
      return h;
    }

    T cp(int i,  const VType& u){
      T cres = T();
      for(int k = 0; k < nchem; ++k){ //fill also Cp_
	Cp_[k] = BCS_.C_p(AccessType::spec_access(k),u[OFF(i,TMP)]);
	cres += (Cp_[k]/DataType::molar_masses()[AccessType::spec_access(k)])*u[OFF(i,SPC+k)];
      }
      //isCpcalculated_ = true;
      return cres;
    }

     //! EVALUATE TRANSPORT COEFFICIENTS AT point node \f$ (i) \f$
    T lambda(int i, const VType& u){
      return Mac_.compute_mixture_averaged_conductivity((*this).rho(i,u), p0_, u[OFF(i,TMP)],(*this).mole_fraction(X_,i,wbar(i,u),u)); //X_ overwritten
    }
  
   T heat_conduction(int i, const VType& u){
     return ( (0.5*(lambda(i+1,u)+lambda(i,u))*(u[OFF(i+1,TMP)]-u[OFF(i,TMP)]) - 0.5*(lambda(i,u)+lambda(i-1,u))*(u[OFF(i,TMP)]-u[OFF(i-1,TMP)]))/ntimes<2>(hx_)  );
    }
   
   //cp() must be invoked beforehand
     T temperature_diffu_term(int i, const VType& u){
       T td = T(),
	 meanmolmass = wbar(i,u),
	 Rho = rho(i,u);

       mole_fraction(X_,i,meanmolmass,u);
       Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,TMP)],X_); //diffusion coefficients are computed by now
       
       //isDiffusioncalculated_ = true;
       
       for(int k = 0; k < nchem; ++k){
	 //!store counterclockwise, starting from the bottom with
	//!(i) being the last entry
	threepointstencil_[0] = mole_fraction(i-1,k,wbar(i-1,u),u);
	threepointstencil_[1] = X_[k];   //mole_fraction(i,k,wbar(i,u),u);
	threepointstencil_[2] = mole_fraction(i+1,k,wbar(i+1,u),u);

	//! assume Cp_ has been filled previously
	td += ( (Cp_[k]/meanmolmass)*Mad_[k]*Rho*( (0.5*(threepointstencil_[2]+threepointstencil_[1])*(u[OFF(i+1,TMP)]-u[OFF(i,TMP)]) - 0.5*(threepointstencil_[1]+threepointstencil_[0])*(u[OFF(i,TMP)]-u[OFF(i-1,TMP)]))/ntimes<2>(hx_) ) );
      }
      return td;
     }
  
   ArrayType& species_diffu_term(ArrayType& diffu, int i, const VType& u){
      threepointstencil_[0] = wbar(i-1,u);
      threepointstencil_[1] = wbar(i,u);
      threepointstencil_[2] = wbar(i+1,u);
      
      mole_fraction(X0_,i-1,threepointstencil_[0],u);  //(i-1)
      mole_fraction(X_, i,  threepointstencil_[1],u);  //(i)
      mole_fraction(X2_,i+1,threepointstencil_[2],u);  //(i+1)
      
      Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i-1,TMP)],X0_);
    
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,TMP)],X_);
      Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i+1,TMP)],X2_);
      

      
      for(int k = 0; k < nchem; ++k){
	diffu[k] = (0.5*(rho(i+1,u)*DataType::molar_masses()[AccessType::spec_access(k)]/threepointstencil_[2]*Mad2_[k] + rho(i,u)*DataType::molar_masses()[AccessType::spec_access(k)]/threepointstencil_[1]*Mad_[k])*(X2_[k]-X_[k]) - 0.5*(rho(i,u)*DataType::molar_masses()[AccessType::spec_access(k)]/threepointstencil_[1]*Mad_[k] + rho(i-1,u)*DataType::molar_masses()[AccessType::spec_access(k)]/threepointstencil_[0]*Mad0_[k])*(X_[k]-X0_[k]))/ntimes<2>(hx_);

      }
      
      return diffu;
    }

   
    //with wbar
     ArrayType& species_diffu_term_wbar(ArrayType& diffu, int i, const VType& u){
       
       threepointstencil_[0] = wbar(i-1,u);
       threepointstencil_[1] = wbar(i,u);
       threepointstencil_[2] = wbar(i+1,u);
      
       mole_fraction(X0_,i-1,threepointstencil_[0],u);  //(i-1)
       mole_fraction(X_, i,  threepointstencil_[1],u);  //(i)
       mole_fraction(X2_,i+1,threepointstencil_[2],u);  //(i+1)
      
       Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i-1,TMP)],X0_);
    
       Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i,TMP)],X_);
       Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,u[OFF(i+1,TMP)],X2_);
      
       
       star3_[0] = rho(i-1,u);
       star3_[1] = rho(i,u);
       star3_[2] = rho(i+1,u);
      
       for(int k = 0; k < nchem; ++k){
	//!\f$\nabla \cdot (\rho D^{\mathrm{mix}}_k \nabla Y_k)\f$
	 diffu[k] = ( (0.5*(star3_[2]*Mad2_[k] + star3_[1]*Mad_[k])*(u[OFF(i+1,SPC+k)]-u[OFF(i,SPC+k)]) - 0.5*(star3_[1]*Mad_[k] + star3_[0]*Mad0_[k])*(u[OFF(i,SPC+k)]-u[OFF(i-1,SPC+k)]))/ntimes<2>(hx_) 
	
	  //!differential part with mean molecular weight \f$ \overline{W}\f$,
	  //!i.e. \f$ \nabla \cdot (\rho Y_k/\overline{W}D^{\mathrm{mix}}_k \nabla \overline{W})\f$
		      +  (0.5*(star3_[2]*u[OFF(i+1,SPC+k)]/threepointstencil_[2]*Mad2_[k] + star3_[1]*u[OFF(i,SPC+k)]/threepointstencil_[1]*Mad_[k])*(threepointstencil_[2]-threepointstencil_[1]) - 0.5*(star3_[1]*u[OFF(i,SPC+k)]/threepointstencil_[1]*Mad_[k] + star3_[0]*u[OFF(i-1,SPC+k)]/threepointstencil_[0]*Mad0_[k])*(threepointstencil_[1]-threepointstencil_[0]))/ntimes<2>(hx_) );
	
       }
       return diffu;
     }

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

  };


} //end namespace 


#endif
