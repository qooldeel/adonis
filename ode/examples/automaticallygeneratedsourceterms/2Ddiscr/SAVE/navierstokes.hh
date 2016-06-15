#ifndef REACTIVE_NAVIER_STOKES_EQUATIONS_HH
#define REACTIVE_NAVIER_STOKES_EQUATIONS_HH

#include "../../../../massactionkinetics/reactionrates.hh"
#include "../../../../massactionkinetics/stoichiometry.hh"
#include "../../../../massactionkinetics/thermochemistry.hh"
#include "../../../../massactionkinetics/indexinjector.hh"
#include "../../../../massactionkinetics/buildchemicalrhs.hh"
#include "../../../../massactionkinetics/physicalconstants.hh"
#include "../../../../massactionkinetics/data/thermochemicaldata.hh"

#include "../../../../expressiontemplates/exprvec.hh"
#include "../../../../templatemetaprograms/unrollloop.hh"
#include "../../../../common/typeadapter.hh"
#include "../../../../common/smartassign.hh"
#include "../../../../common/adonisassert.hh"
#include "../../../../io/readinparameters.hh"

//transport quantities
#include "../../../../moleculartransport/viscosity.hh"

namespace Adonis{

  /**
   * \brief Reactive Navier Stokes on a rectangular domain [a,b] x [c,d],
   * discretized with finite differences
   */
  template<class T>
  class ReactiveNavierStokes2D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;
    typedef std::size_t SizeType;

    //================== TODO: set mechanism here ============================
    enum{thermoMech = 3};   //3 = O3, 9 = H2Gri
    //========================================================================

    typedef ThermoData4Mechanism<DType,thermoMech> DataType;  

    //needed for chemistry
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;
    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,CommonIndexInjector> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;


    ReactiveNavierStokes(SizeType n = 0):rhs_(n),C_(DataType::nspec),Nx_(0),Ny_(0),npt_(0),a_(DType()),hx_(DType()),hy_(DType()), ignitionx_(DType()), Tignition_(DType()),cp_(T()){
      ParameterData PD;
      PD.read_from_file("2Dsettings.dat");
      T a = PD.get_datum<DType>("a"),
	b = PD.get_datum<DType>("b"),
	c = PD.get_datum<DType>("c"),
	d = PD.get_datum<DType>("d");
      
      a_ = a;
      Nx_ = PD.get_datum<SizeType>("Nx");
      Ny_ = PD.get_datum<SizeType>("Ny");

      npt_ = Nx_*Ny_;
      hx_ = (b-a)/(Nx_-1);
      hy_ = (d-c)/(Ny_-1);
      
      ignitionx_ = PD.get_datum<DType>("ignition_x");
      Tignition_ = PD.get_datum<DType>("T_ignition");

      //! these are the conservative variables
      rho_.resize(npt);
      mom1_.resize(npt);
      mom2_.resize(npt);
      thermEy_.resize(npt);
    
      partDens_.resize(DataType::nspec*npt_); //species
     

      if(n > 0){
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
	static CommonIndexInjector Cii(DataType::nspec,DataType::nreac);
	BCS_.initialize(&forwardRates,&reverseRates,&Cii);
      }

      //!transport quantities
      Mav_.initialize(DataType::transport());

    } //end constructor
    

    //! OPERATOR -- input primite variables [rho,v1,v2,T,Y1,...,YK]
    template<class E>
    VType& operator()(const E& prim){
      //!reset
      cp_ = T();
      
      //! flow -- conservation variables
      for(SizeType i = 0; i < Nx_; ++i){
	for(SizeType j = 0; j < Ny_; ++j){
	  rho_[(i,j)] = prim[(i,j)];
	  mom1_[(i,j)] = rho_*prim[(i,j,1)];
	  mom2_[(i,j)] = rho_*prim[(i,j,2)];
	 
	  for(SizeType k = 0; k < DataType; ++k){
	    cp_ += BCS_.C_p(k,prim[(i,j,3)])/DataType::molar_masses()[k]*prim[(i,j,4+k)];
	    //! concentrations
	    C_[k] = rho_[(i,j)]*prim[(i,j,4+k)]/DataType::molar_masses()[k];
	  }
	  
	  thermEy_[(i,j)] = rho_*cp_*prim[(i,j,3)];
	}
      }

      //! species -- conservation variables
      for(SizeType k = 0; k < DataType::nspec; ++k){
	for(SizeType i = 0; i < Nx_; ++i){
	  for(SizeType j = 0; j < Ny_; ++j){
	    partDens_[(i,j,k)] = rho_*prim[(i,j,4+k)];
	  }
	}
      }
      
      //! BOUNDARY CONDITIONS
      for(SizeType i = 0; i < Nx_; ++i){
	//down
	//(i,0)
	//up
	//(i,Ny_-1)
      }
      for(SizeType j = 0; j < Ny_; ++ij){
	//left
	//(0,j)
	//right
	//(Nx_-1,j)
      }
      

      return rhs_;
    }
    

  private:
    VType rhs_, C_;
    SizeType Nx_, Ny_, npt_;
    DType a_,hx_, hy_,ignitionx_,Tignition_;
    T cp_;
    VType rho_,    //! \rho
      mom1_,       //! \rho v_1
      mom2_,       //! \rho v_2
      thermEy_,    //! \rho c_p T
      partDens_;   //! \rho Y_k
      

    BuilderType BCS_;  //calculate chemistry
   
    MixtureAveragedViscosity<DataType,false> Mav_;

    SizeType operator()(SizeType i, SizeType j){
      return (i + Nx_*j);
    }
    
    SizeType operator()(SizeType i, SizeType j, SizeType spec){
      return ((i,j) + spec*npt_);
    }
  };

}

#endif
