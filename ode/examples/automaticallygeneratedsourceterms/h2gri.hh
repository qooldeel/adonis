#ifndef CHEMICAL_SOURCE_TERM_GENERATED_WITH_THERMOCHEMICAL_DATA_HH
#define CHEMICAL_SOURCE_TERM_GENERATED_WITH_THERMOCHEMICAL_DATA_HH

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

#include "../../../massactionkinetics/temperaturesource.hh"

namespace Adonis{

  //test me 
  template<class T>
  class H2Gri{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;


    enum{h2gri = 9};  //! mechanism signature

    typedef ThermoData4Mechanism<DType,h2gri> DataType;
    typedef Quantity<Settings4Quantity::ofcomputation> QType;
   
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;

    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;
    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;

    typedef BuildChemicalSourceTerm<FwdType,RevType,CommonIndexInjector> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    typedef TemperatureSourceTerm<Settings4Quantity::thermoCase,T,BuilderType> TempSourceType;

    inline H2Gri(size_t dim = 0):rhs_(dim){
      

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
	static RevType reverseRates(&forwardRates,&stoichmatrix,&nasa);

	//this is just the normal index
	static CommonIndexInjector Cii(DataType::nspec,DataType::nreac);

	BCS_.initialize(&forwardRates,&reverseRates,&Cii);
	
	TempSource_.initialize(T(0.),&BCS_);

      } //end if dim > 0
    }//end constructor

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    
    RevPointerType rev_address() const {return BCS_.address_reverse();}
    
    size_t nspec() const{return BCS_.number_of_species();}
  

    template<class E>
    VType& operator()(const E& u){
      /*T phi_sum = T(),
	sum1 = T(),
	sum2 = T(),
	temp;
      */

      //convert to right quantity
      QType::concentration(nspec(),conc_,u,DataType::molar_masses(),DataType::density_of_mixture());   //transform into concentrations. When u = [X] do nothing
      
      //!THEMPERATURE
      // temp = u[nspec()];
      TempSource_.reset();
      TempSource_.set_temperature(u[nspec()]);
				  
      conc_[nspec()] = u[nspec()]; //the last entry will be temperature again

     
      for(size_t k = 0; k < nspec(); ++k){
	//                                                       temp
	rhs_[k] = BCS_.net_formation_rate(BCS_.species_index(k),TempSource_.get_temperature(),QType::choose_right_quantity(u,conc_));  //!if u is already concentration return u, in other cases return conc_ 
	

	//! TEMPERATURE
	//	phi_sum +=  u[k];
	//sum2 += (u[k]*BCS_.C_p(k,temp));
	TempSource_.compute_phi_sum(u[k]);
	TempSource_.compute_sum2(k,u[k]);

	TempSource_.compute_sum1(k,rhs_[k]);
	/*#if TEMP_ISO_CASE == 1
	  
	  sum1 += ((BCS_.H_T(k,temp) - PhysicalConstants<T>::IdealGasConstant*temp)*rhs_[k]);
	  #endif
	  
	  #if TEMP_ISO_CASE == 2
	  sum1 += (BCS_.H_T(k,temp)*rhs_[k]);
	  #endif
	*/
	
      } //END FOR over species
				  
      
      TempSource_.temperature(rhs_[nspec()],DataType::density_of_mixture(),PhysicalConstants<T>::p0);
				  /*
				    #if TEMP_ISO_CASE == 1
      rhs_[nspec()] = -1./DataType::density_of_mixture()*(sum1/(sum2 - PhysicalConstants<T>::IdealGasConstant*phi_sum));
#endif

#if TEMP_ISO_CASE == 2
      rhs_[nspec()] = -PhysicalConstants<T>::IdealGasConstant*temp/PhysicalConstants<T>::p0*((sum1*phi_sum)/sum2);

#endif
*/

      //! CONVERT BACK 2 INPUT QUANTITY (rhs_ back to specific moles, etc. If
      //! input species is already concentration, do nothing
      QType::back_to_input_quantity(nspec(),rhs_,DataType::molar_masses(),DataType::density_of_mixture());
      
      return rhs_;
    }


    //! fake time dependent operator
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
    }


  private:
    VType rhs_, conc_;
    BuilderType BCS_;
    
    TempSourceType TempSource_;
  };

} //end namespace 

#endif
