#ifndef REDUCED_OZONE_O3_MECHANISM_HH
#define REDUCED_OZONE_O3_MECHANISM_HH

//! COMPUTE IN PARTIAL DENSITIES ONLY
#include "../../../massactionkinetics/reactionrates.hh"
#include "../../../massactionkinetics/stoichiometry.hh"
#include "../../../massactionkinetics/thermochemistry.hh"
#include "../../../massactionkinetics/indexinjector.hh"
#include "../../../massactionkinetics/buildchemicalrhs.hh"
#include "../../../massactionkinetics/physicalconstants.hh"
#include "../../../massactionkinetics/data/thermochemicaldata.hh"

#include "../../../expressiontemplates/exprvec.hh"
#include "../../../templatemetaprograms/unrollloop.hh"
#include "../../../common/typeadapter.hh"
#include "../../../common/smartassign.hh"
#include "../../../common/adonisassert.hh"
#include "../../../common/typeadapter.hh"

#include "../../../containers/staticarray.hh"

//RED
#include "../../../fdm/reconstructspecies.hh"
#include "../io/readinparameters.hh"


namespace Adonis{

 template<class T> 
  class ReducedO3Decomposition{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;

    enum{o3decomp = 3};  //! mechanism signature

    typedef ThermoData4Mechanism<DType,o3decomp> DataType;  
   
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;
    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
   //RED 
   typedef IndexHandler<DataType,true> AccessType;
   typedef typename AccessType::IndexerType IndexerType;
   typedef BuildChemicalSourceTerm<FwdType,RevType,IndexerType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;


   //RED
   typedef std::size_t SizeType;
   typedef typename TypeAdapter<T>::BaseType BaseType;
   typedef StaticArray<T,DataType::rednspec> RedChemType;
   typedef StaticArray<T,DataType::nspec> FullChemType;
   typedef ReconstructSpecies<DType,DataType::nspec,DataType::rednspec,o3decomp> ReconstructorType;
   typedef ExprTmpl::MyVec<DType> VDType;

   ReducedO3Decomposition(size_t dim = 0, const T& rho = 1.):rhs_(dim),y_(DataType::rednspec+1),rho_(rho),qheatReac_(T()),cp_(T()),wbar_(T()){
      if(dim > 0){
		//!construct objects needed for source term
	//! CAUTION: use <TT>static</TT> to maintain each objects' location 
	//!   beyond the call of the constuctor (i.e. beyond the localness of 
	//!   the {...}-block
	static ThermoType nasa(DataType::nspec,DataType::thermo(),DataType::temperature_bounds());
       
	static StoichType stoichmatrix(DataType::nspec,DataType::nreac,DataType::stoichiometric_matrix());

	static TroeIndex TIndex;
	TIndex.create(DataType::nreac,DataType::troewtb());
	
	//! everything stated as  pointers	
	static FwdType forwardRates(DataType::nspec,DataType::nreac,DataType::arrhenius_values(),DataType::temperature_bounds(), DataType::ntroereac, DataType::troewtb(),DataType::collision_efficiencies(),&TIndex);

	//! everything stated as  pointers	
	static RevType reverseRates(&forwardRates,&stoichmatrix,&nasa);

	//this is just the normal index
	static IndexerType Ixer; 
	AccessType::initialize(Ixer);

	BCS_.initialize(&forwardRates,&reverseRates,&Ixer);
	
	//===RED
	ParameterData PD;
	PD.read_from_file("/home/mfein/MARC++/fdm/2Dsettings.dat");
	VType yt0(DataType::nspec+1);
	
	yt0 <<= 0., 0.8, 0.2, PD.get_datum<T>("T_ignition"); //no tabularization needed here
	  Recon_.initialize(yt0,PD.get_datum<DType>("hEstim"),PD.get_datum<BaseType>("newtTol"),PD.get_datum<char>("whichSolver"), PD.get_datum<int>("maxit"), PD.get_datum<BaseType>("nmtol"),PD.get_datum<BaseType>("Tol"));
      
	Ctilde_.resize(1*3);
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

	zM_.resize(DataType::nspec+1);
	zM_ = yt0;  //initial guess

	p0_ = PD.get_datum<T>("pconst");
	//=========
      }
    }  //end constructor

    //4 (3 species (O2, O3+ temperature). Then O2 = 1 - (O2 + O3)
    size_t dim() const {return 2;}
    size_t domain_dim() const {return 2;}

    T qheatReac() const {return qheatReac_;}
    
    template<class E> 
    VType& operator()(const E& u){
      //std::cout << "Invoke O3 operator "<< std::endl;
      y_ = u;   //copy, y_ is the reduced type here
     
      //=====RED
      reconstruct_full_state(zM_,y_);
      // cp_ = Recon_.cp();
      // rho_ = Recon_.rho();
      //========

      qheatReac_ = cp_ = rho_ = wbar_ = T();        //reset
     
      T temperature = zM_[DataType::nspec];//y_[DataType::rednspec];//  //y_[DataType::rednspec];
           
      rho_ = rho();
      
      for(int k = 0; k < DataType::nspec; ++k){
	Conc_[k] = rho_*zM_[k]/DataType::molar_masses()[k];
	cp_ += BCS_.C_p(k,temperature)/DataType::molar_masses()[k]*zM_[k];
	//qheatReac_ += BCS_.H_T(k,zM_[DataType::nspec])*Recon_.dot_omega()[k];
      }
      //cp_ = Recon_.cp();
      
      for(size_t k = 0; k < DataType::rednspec; ++k){	
	//! compute \f$ \dot{\omega}_k \ \forall k\f$ 
	dotomega_[k] = BCS_.net_formation_rate(BCS_.species_index(k),temperature,Conc_);
  
	rhs_[k] = 1./rho_*(dotomega_[k]*DataType::molar_masses()[AccessType::spec_access(k)]);

	//rhs_[k] = 1./rho_*Recon_.dot_omega()[AccessType::spec_access(k)]*DataType::molar_masses()[AccessType::spec_access(k)];

      }
      //take reduced dot omega and full dotomega from species reconstruction
      qheatReac_ =  BCS_.H_T(0,zM_[DataType::nspec])*dotomega_[0]+ BCS_.H_T(1,zM_[DataType::nspec])*Recon_.dot_omega()[1] + BCS_.H_T(2,zM_[DataType::nspec])*Recon_.dot_omega()[2];

	

      //qheatReac_ = Recon_.heat();
      //! heat source MIND THE MINUS HERE!!!!!!!!!!
      rhs_[DataType::rednspec] = -qheatReac_/(rho_*cp_);
      //std::cout << "TEMPERATURE = "<< rhs_[3] << std::endl;
     
      return rhs_;
    }

    
    //! fake time dependent operator -- needed for my EXPLICIT methods
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
    }



  private:
   VType rhs_, 
     y_, zM_;
   T rho_, qheatReac_,cp_, p0_, wbar_;
    BuilderType BCS_;
    
   //===RED
   RedChemType dotomega_, redspecies_;
   FullChemType Conc_;
   ReconstructorType Recon_;
   VDType Ctilde_, btilde_, low_,up_;
  

   VType& reconstruct_full_state(VType& zm, const VType& vars){
     for(int k=0; k < DataType::rednspec; ++k)
       redspecies_[k] = vars[k];   //first entries are variables

     Recon_.assign(redspecies_,zM_); //zM_ overwritten with start assignment
     Recon_.enforce_linearized_constraints_wrt_chemistry(Ctilde_,btilde_);
     Recon_.project_onto_bounds(low_,up_);
     zm = Recon_.time_step(Recon_.get_z(),p0_);
   
     return zm;
   }

   T& wbar(){
     wbar_ = T();
     for(int k = 0; k < DataType::nspec; ++k){ //full concentrations means full wbar!!
       wbar_ += zM_[k]/DataType::molar_masses()[k];
     }
     wbar_ = 1./wbar_;
     return wbar_;
   }

   T rho(){
     return (p0_*wbar())/(PhysicalConstants<DType>::Rgas*zM_[DataType::nspec]);
   }

   //===

};

} //end namespace 

#endif
