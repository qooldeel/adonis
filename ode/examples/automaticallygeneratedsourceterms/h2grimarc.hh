#ifndef EXCERPT_OF_H2_COMBUSTION_OF_GRI_MECH_3_0_HH
#define EXCERPT_OF_H2_COMBUSTION_OF_GRI_MECH_3_0_HH

#include<iostream>

#include "../../../massactionkinetics/quantity.hh"
#include "../../../massactionkinetics/qsettings.hh"
#include "../../../massactionkinetics/tsettings.hh"
#include "../../../massactionkinetics/reactionrates.hh"
#include "../../../massactionkinetics/stoichiometry.hh"
#include "../../../massactionkinetics/thermochemistry.hh"
#include "../../../massactionkinetics/indexinjector.hh"
#include "../../../massactionkinetics/buildchemicalrhs.hh"
#include "../../../massactionkinetics/physicalconstants.hh"
#include "../../../massactionkinetics/data/thermochemicaldata.hh"

#include "../../../expressiontemplates/exprvec.hh"
#include "../../../common/typeadapter.hh"
#include "../../../common/smartassign.hh"

#include "../../../templatemetaprograms/unrollloop.hh"

namespace Adonis{


  /**
   * \brief Implementation of the "reactive" part of a isobaric, quasi-1D flame
   * propagation, cf. [1]
   *
   * References:
   *   [1] [Chemkin-Pro, Reaction Design: San Diego 2008, "Theory Manual", § 12, pp. 201]
   */
  template<class T>
  class H2GriMech3{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;

    enum{h2gri = 9};  //! mechanism signature
    
    typedef ThermoData4Mechanism<DType,h2gri> DataType;
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;
    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,CommonIndexInjector> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    H2GriMech3(size_t dim = 0):p_(T()),temp_(T()),rho_(T()),A_(T()),cpk_(T()),rhs_(dim){ //the plus one denotes the energy equation 
      if(dim > 0){
	conc_.resize(DataType::nspec);
	spec_.resize(DataType::nspec);
     
	//============== p = const. ===========================================
	p_ = PhysicalConstants<T>::p0;

	//! mean piston speed 
	//v_ = 16.;    //m/sec for medium speed petrol, (high speed: 20-25m/sec

	//area in m²
	A_ = 0.01;

	rho_ = 0.3;   // at the beginning -- some meaningful default value
	//====================================================================

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
	
      } //end if dim > 0
    }//end constructor

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}
   
    size_t nspec() const{return BCS_.number_of_species();}
  
    
    template<class E>
    VType& operator()(const E& u){
      // std::cout << "rhs_.size() = "<< rhs_.size() << std::endl;
      //=======================================================================
      T mmw = T(),
	energysource = T(), 
	cp = T();

      temp_ = u[DataType::nspec]; //here we store the thermal energy
      
      //!mean molecular weight
      for(size_t k = 0; k < DataType::nspec; ++k){
	mmw += u[k]/DataType::molar_masses()[k];
      }
      adonis_assert(Abs(mmw) > 0.);
      
      mmw = 1./mmw;
      //std::cout << "mean molecular weight = "<< mmw << std::endl;

      //! ++++++++++  EOS +++++++++++++++++++++++++++++++++++++++++++++++
      //rho_ = p_*mmw/(PhysicalConstants<T>::IdealGasConstant*temp_);
      //std::cout << " DENSITY rho = "<< rho_ << std::endl; 
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      //!compute concentration from mass fraction
      for(size_t k = 0; k < DataType::nspec; ++k){
	conc_[k] = rho_*u[k]/DataType::molar_masses()[k];
      }
      //std::cout << "conc = "<< conc_ << std::endl;
 
      chemical_source(conc_,temp_);   //!compute \f$ \dot w_k\f$ 
      
      for(size_t k = 0; k < DataType::nspec; ++k){
	//! \f$ \frac{d Y_k}{dt}\f$
	spec_[k] = rhs_[k]*DataType::molar_masses()[k]/rho_; //unit: 1/sec

	cpk_ =  BCS_.C_p(k,temp_)/DataType::molar_masses()[k];
	energysource += rhs_[k]*BCS_.H_T(k,temp_);
	cp += cpk_*u[k]; //u[k]= Y[k] 
      }

      //! \f$ \frac{dT}{dt}\f$ 
      energysource *= -1./(rho_*cp);  

      //o.k, spec is smaller by one.
      UnrollLoop<0,DataType::nspec>::assign_rac2rac(rhs_,spec_);
      //std::cout << "species = "<< spec_ << std::endl; 
      rhs_[DataType::nspec] = energysource;


      //No transformation back to mass fractions, since we already compute in
      //mass fractions Y_k :)
      
      //========================================================================
      
      return rhs_;
    }


  private:
    T p_,temp_,rho_,A_, cpk_;
    VType rhs_, conc_, spec_;
    BuilderType BCS_;
    
    template<class CONC, class THETA>
    VType& chemical_source(const CONC& conc, const THETA& temperature){
      
      //! loop only over species (first nspec elements)
      for(size_t k = 0; k < nspec(); ++k){
	rhs_[k] = BCS_.net_formation_rate(BCS_.species_index(k),temperature,conc);
      }
      return rhs_;
    }
  };


} //end namespace 

#endif
