#ifndef AUTOMATICALLY_GENERATED_SOURCE_TERM_4_FICTITIOUS_H2_COMBUSTION_HH
#define AUTOMATICALLY_GENERATED_SOURCE_TERM_4_FICTITIOUS_H2_COMBUSTION_HH

#include<iostream>

#include "../../../massactionkinetics/quantity.hh"
#include "../../../massactionkinetics/qsettings.hh"
#include "../../../massactionkinetics/tsettings.hh"
#include "../../../massactionkinetics/reactionrates.hh"
#include "../../../massactionkinetics/stoichiometry.hh"
#include "../../../massactionkinetics/thermochemistry.hh"
#include "../../../massactionkinetics/indexinjector.hh"
#include "../../../massactionkinetics/buildchemicalrhs.hh"

#include "../../../expressiontemplates/exprvec.hh"
#include "../../../massactionkinetics/data/thermochemicaldata.hh"
#include "../../../common/typeadapter.hh"

#include "../../../common/smartassign.hh"
#include "../../../common/adonisassert.hh"
#include "../../../templatemetaprograms/unrollloop.hh"


namespace Adonis{
  
  /**
   * \brief This is an isothermal combustion mechanism
   */
  template<class T>
  class H2MechIn6Species{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;


    enum{h2c6 = 6};  //! mechanism signature or encoding

    typedef ThermoData4Mechanism<DType,h2c6> DataType;
    //THIS MECHANISM IS ALREADY BASED ON CONCENTRATIONS
    typedef Quantity<'c'> QType;
   
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;

    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;

    //typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    //NOTE: This mechanism is isothermal!!!!!
    typedef IsoThermalReverseReactionRates<VType,StoichType> RevType;


    typedef BuildChemicalSourceTerm<FwdType,RevType,CommonIndexInjector> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;
  
  
    H2MechIn6Species(size_t dim = 0):rhs_(dim){
    
      if(dim >1){
	QType::resize(dim,conc_);
      }
      if(dim >0){
	
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
	//static RevType reverseRates(&forwardRates,&stoichmatrix,&nasa);
	//!NOTE: mechanism is isothermal!!!
	static T rr[] = {216., 337.5, 1400., 10800., 33750., 0.7714}; //!these are fixed reverse rates for this mechanism
	static RevType reverseRates(rr,rr+DataType::nspec,&stoichmatrix);

	static CommonIndexInjector Cii(DataType::nspec,DataType::nreac);

	BCS_.initialize(&forwardRates,&reverseRates,&Cii);
	
      } //end if dim > 0
    }//end constructor

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    size_t nspec() const{return BCS_.number_of_species();}

   
    template<class E>
    VType& operator()(const E& u){
      
      //convert to right quantity
      QType::concentration(nspec(),conc_,u,DataType::molar_masses(),DataType::density_of_mixture());   //transform into concentrations. When u = [X] do nothing
      
     
      
      for(size_t k = 0; k < nspec(); ++k){
	//! note that the 2nd argument can mean anything since we consider only
	//! an <B> isothermal </B> mechanism here
	//std::cout << "index = "<<BCS_.species_index(k) << std::endl;
	//! Caution: type <TT> T(1000.) </TT>.Otherwise 1000. is considered a 
	//! double which it is not when using integration ;)
	rhs_[k] = BCS_.net_formation_rate(BCS_.species_index(k),T(1000.),QType::choose_right_quantity(u,conc_));  //! if u is concentration, choose take u else conc_
      }

      
      QType::back_to_input_quantity(nspec(),rhs_,DataType::molar_masses(),DataType::density_of_mixture());

      return rhs_;
     }

    //! fake time dependent operator
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
    }

    //! subtraction of equations cancels OH
    //! 2[H2] + [H] - 2[O2] - [O] + [H2O] = b1-b2
    //!
    template<class V, class W>
    typename V::value_type mass_balance_for_02(const V& b, const W& u) const{
      typedef typename V::value_type value_type;
      value_type O2;
      //                          H2     H      O      H2O  (OH was canceled)
      O2 = (-0.5*((b[0]-b[1]) - 2*u[0] - u[1] + u[3] - u[4]));
      
      std::cout << "[O2] = "<< O2 << std::endl;

      if(O2 < value_type())
	ADONIS_ERROR(PositivityError, "Damn it! [O2] < 0.");
      return O2;
    }

    template<class V, class W>
    typename V::value_type mass_balance_for_H2O(const V& b, const W& u) const{
      typedef typename V::value_type value_type;
      value_type H2O;
      H2O = (b[0]-b[1] - 2*u[0]- u[1] + 2*u[2] + u[3]);
      
      std::cout << "[H2O] = "<< H2O << std::endl;

      if(H2O < value_type())
	ADONIS_ERROR(PositivityError, "Damn it! [H2O] < 0.");
      return H2O;
    }

  private:
    VType rhs_, conc_;
    BuilderType BCS_; 

  };

} //end namespace

#endif
