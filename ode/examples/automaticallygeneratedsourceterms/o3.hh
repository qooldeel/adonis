#ifndef CHEMICAL_SOURCE_TERM_FOR_O3_DECOMPOSITION_HH
#define CHEMICAL_SOURCE_TERM_FOR_O3_DECOMPOSITION_HH

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

#include "../../../containers/staticarray.hh"

namespace Adonis{

  //! this is only the source term of th species conservation, i.e.
  //! \f[ \partial_t(\rho Y_k) + \nabla\cdot(\rho Y_kv + j_K) = \dot{\omega}_kW_k, \quad k = 1,\ldots, K.\f]
  //! note that we \f$ \rho \f$ is const. in time and that we assume isobaric 
  //! settings --> use <B> PRIMITIVE </B> variables for species and thermal energy
  template<class T> 
  class O3Decomposition{
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
    typedef BuildChemicalSourceTerm<FwdType,RevType,CommonIndexInjector> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;


    O3Decomposition(size_t dim = 0, const T& rho = 1.):rhs_(dim),conc_(3),y_(4),rho_(rho),qheatReac_(T()),cp_(T()){
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
	static CommonIndexInjector Cii; //(DataType::nspec,DataType::nreac)
	Cii.initialize(DataType::nspec,DataType::nreac);

	BCS_.initialize(&forwardRates,&reverseRates,&Cii);
	
      }
    }  //end constructor

    //4 (3 species (O2, O3+ temperature). Then O2 = 1 - (O2 + O3)
    size_t dim() const {return 4;}
    size_t domain_dim() const {return 4;}

    T qheatReac() const {return qheatReac_;}
    
    template<class E> 
    VType& operator()(const E& u){
      //std::cout << "Invoke O3 operator "<< std::endl;
      y_ = u;   //copy
      qheatReac_ = cp_ = T();        //reset

      T temperature = y_[3];
      
      ////!calculate new rho
      T mmw = T();
      UnrollLoop<0,3>::mean_molecular_weight_Y(mmw,&y_[0],DataType::molar_masses());
      rho_ = 101325*mmw/(PhysicalConstants<T>::IdealGasConstant*temperature);
	
      y_[2] = 1. - (y_[0]+y_[1]);   //obeys mass conservation

      //adonis_assert(y_[0] > T());

      for(size_t k = 0; k < 3; ++k){
	//equals Jochen's sum2
	cp_ += ( BCS_.C_p(k,temperature)*y_[k]/DataType::molar_masses()[k] );
	//! compute concentrations
	conc_[k] = rho_*y_[k]/DataType::molar_masses()[k];  
      }

      for(size_t k = 0; k < 3; ++k){	
	
	//! compute \f$ \dot{\omega}_k \ \forall k\f$ 
	dotomega_[k] = BCS_.net_formation_rate(k,temperature,conc_);
	//std::cout << "rhs_[ "<<k<<"] = "<< rhs_[k] << std::endl;
	//rhs_[k] *= DataType::molar_masses()[k]/rho_;  //! Y_k
	rhs_[k] = (DataType::molar_masses()[k]/rho_)*dotomega_[k]; //back from [X_k] to Y_k
	qheatReac_ += BCS_.H_T(k,temperature)*dotomega_[k]; //sum1
      }
      
      //! heat source MIND THE MINUS HERE!!!!!!!!!!
      rhs_[3] = -qheatReac_/(rho_*cp_);
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
    VType rhs_, conc_;
    VType y_;
    T rho_, qheatReac_,cp_;
    BuilderType BCS_;
    StaticArray<T,3> dotomega_;
  };
}

#endif
