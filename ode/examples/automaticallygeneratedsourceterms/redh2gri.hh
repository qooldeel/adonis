#ifndef REDUCED_H2GRI_SOURCE_HH
#define REDUCED_H2GRI_SOURCE_HH


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
#include "../../../common/globalfunctions.hh"
#include "../../../io/readinparameters.hh"

namespace Adonis{
  
  
  template<class T>
  class ReducedH2Gri{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;

    enum{h2gri = 9};  //! mechanism signature

    typedef ThermoData4Mechanism<DType,h2gri> DataType;
    //the mechanism will be stated in, e.g., specific moles 
    typedef Quantity<Settings4Quantity::ofcomputation> QType;

     enum{
       reducedDimension = DataType::rednspec,
       fullDimension = DataType::nspec
    };


    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;

    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;
    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;

    typedef ReducedIndexInjector<const size_t*, const size_t*> InjectorType;

    typedef BuildChemicalSourceTerm<FwdType,RevType,InjectorType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;
    
    //! temperature source term
    //typedef TemperatureSourceTerm<Settings4Quantity::thermoCase,T,BuilderType> TempSourceType;
    
    //== DIFFERENCE Manifold stuff -- thermally dependent mechanism ===========
    //============= Thermally dependent mechanism =============================
    typedef ComputeManifold<DType,DataType::nspec,reducedDimension,false> ManifoldType;
    //==========================================================================
    typedef typename ManifoldType::CompositionType CompositionType;
    typedef typename ManifoldType::VecD VecD; //std::vector<double>
    typedef typename ManifoldType::VecS VecS; //std::vector<std::string> 
    typedef typename ManifoldType::CoordinateType CoordinateType; //Dune::FieldVector<int,·>
    
    typedef Dune::FieldVector<T,DataType::nspec> FullType;
    

    //************** TODO: change here *****************************************
    enum{iniTrans = 1,
	 coarseEv = 1
    };
    //**************************************************************************


    ReducedH2Gri(size_t dim = 0, 
		 size_t iBd = iniTrans,
		 size_t cE = coarseEv):rhs_(dim),countInvocations_(0),initialLayerBound_(iBd), coarseEvals_(cE){
     

       if(dim > 0) {
	 reduced_.resize(reducedDimension);

	 
	 Param_.read_from_file("dat/sim_start_H2_GRI.dat");

	 //initialize Jochen's more::MoRo
	 VecD pvvalues(reducedDimension);
	 
	 VecS names(reducedDimension); //names of of rpv
	 CoordinateType rpvindex; 
	 
	 for(size_t i = 0; i < reducedDimension; ++i){
	   names[i] = DataType::species_names()[DataType::rpv_index()[i]];
	   pvvalues[i] = Param_.get_datum<typename VecD::value_type>(names[i]);
	   rpvindex[i] = DataType::rpv_index()[i];
	 }
	 
	 std::cout << "names = "; print_all(names); 
	 std::cout << "pvvalues = "; print_all(pvvalues);
	 std::cout << "rpvindex = "<< rpvindex << std::endl;
	 
     
	 //std::cout << "names = "; print_all(names); std::cout << std::endl << " pvvalues = "; print_all(pvvalues); std::cout << std::endl << " rpvindex = "<< rpvindex <<std::endl; 


	 std::cout << "H2GRI: invoke JOCHEN's first:"<< std::endl;

	 Mani_.activate(names,pvvalues,rpvindex); //first(...) invoked
	 
	 std::cout << "exit Jochen's first with z_0 = "<< Mani_.get_z() << std::endl;

      //build rhs 
     
	
	 QType::resize(dim,conc_);  //resize it or not
	 
	 //note that the data operator on the full dimensions
	 static ThermoType nasa(DataType::nspec,DataType::thermo(),DataType::temperature_bounds());
	
	 
	 static StoichType stoichmatrix(DataType::nspec,DataType::nreac,DataType::stoichiometric_matrix());


	 static TroeIndex TIndex;
	 TIndex.create(DataType::nreac,DataType::troewtb());
	 
	 //! everything stated as  pointers	
	 static FwdType forwardRates(DataType::nspec,DataType::nreac,DataType::arrhenius_values(),DataType::temperature_bounds(), DataType::ntroereac, DataType::troewtb(),DataType::collision_efficiencies(),&TIndex);

	 static RevType reverseRates(&forwardRates,&stoichmatrix,&nasa);

	 //reduced injector
	 static InjectorType Inj(reducedDimension,DataType::rednreac,DataType::rpv_index(), DataType::rpv_reaction_index());
	
	 BCS_.initialize(&forwardRates,&reverseRates,&Inj);

	 //TempSource_.initialize(Param_.get_datum<T>("T"),&BCS_);
       }
    } //end constructor


    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    size_t nspec() const{return BCS_.number_of_species();}

    //you might change invocations 
    size_t count_invocations() {return countInvocations_;}


    //! ============== EVALUATION ==================
    //! the parenthesis-operator
    template<class RTYPE>
    VType& operator()(const RTYPE& u){ //! note that u is the reduced type now!
      countInvocations_++; //increase counter 
     
      std::cout << "transform to proper quantities"<< std::endl;

      //convert to right quantity or not
      QType::concentration(nspec(),conc_,u,DataType::molar_masses(),DataType::density_of_mixture()); 

    
      if(countInvocations_ <= initialLayerBound_ || countInvocations_%coarseEvals_ == 0){
	UnrollLoop<0,reducedDimension>::assign_rac2rac(reduced_,QType::choose_right_quantity(u,conc_)); 
      
	std::cout << "reduced_ = "; print_all(reduced_);

	std::cout << "invoke JOCHEN's .warm(…)"<<std::endl;
	Mani_.evaluate(reduced_);  //compute zM
	std::cout << "Invocation .warm(…) DONE!"<<std::endl;
	std::cout << "CAN YOU READ THIS?????????"<<std::endl;

	UnrollLoop<0,DataType::nspec>::assign_rac2rac(zM_,Mani_.get_z());
	//in order to prevent numerical inaccuracies, assign back the identity
	//values (values of the rpvs)
	UnrollLoop<0,reducedDimension>::assign_smaller_2_greater_rac(zM_,u,DataType::rpv_index());
	std::cout << "zM_ = " << zM_ << std::endl;

	//construct reduced RHS
	for(size_t k = 0; k < reducedDimension; ++k){

	  //! note that \f$T\f$ does not serve as rpv and we do not need to
	  //! compute temperature explicitly here. Hence, zM_[nspec] is T
	  rhs_[k] = BCS_.net_formation_rate(BCS_.species_index(k),zM_[DataType::nspec],zM_);
	  
	}
      }

      QType::back_to_input_quantity(reducedDimension,rhs_,DataType::molar_masses(),DataType::density_of_mixture(),DataType::rpv_index());

      return rhs_;
    }

    //! fake time dependent operator
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
    }

  private:
    VType rhs_,
      conc_;
    BuilderType BCS_; 
    
    //TempSourceType TempSource_;

    ManifoldType Mani_;
    VecD reduced_;
    FullType zM_; 
     
    size_t countInvocations_,
      initialLayerBound_, coarseEvals_;
    
    ParameterData Param_;
  };


}

#endif
