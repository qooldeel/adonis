#ifndef REDH2C6_AN_ATTEMPT_BY_ME_HH
#define REDH2C6_AN_ATTEMPT_BY_ME_HH

#include<iostream>

#include "../../../massactionkinetics/quantity.hh"
#include "../../../massactionkinetics/qsettings.hh"
#include "../../../massactionkinetics/tsettings.hh"
#include "../../../massactionkinetics/reactionrates.hh"
#include "../../../massactionkinetics/stoichiometry.hh"
#include "../../../massactionkinetics/thermochemistry.hh"
#include "../../../massactionkinetics/indexinjector.hh"
#include "../../../massactionkinetics/buildchemicalrhs.hh"
#include "../../../massactionkinetics/data/checkmassbalance.hh"

#include "../../../expressiontemplates/exprvec.hh"
#include "../../../massactionkinetics/data/thermochemicaldata.hh"
#include "../../../common/typeadapter.hh"

#include "../../../common/smartassign.hh"
#include "../../../templatemetaprograms/unrollloop.hh"

//for reduced stuff
#include "../dunestuff/fvector.hh"
#include "../../../modred/speciesreconstruction.hh"
#include "../../../modred/verysimple.hh"

//#include "../../../modred/lebiedzapproach_test.hh"



//that's the full source term needed anyway...
#include "../../../ode/examples/automaticallygeneratedsourceterms/h2c6.hh"
#include "../../../io/readinparameters.hh"

#include "../../../modred/simstartingpoint.hh"

namespace Adonis{

   //==========================================================================
  //===================== REDUCED SOURCE TERM ================================
  //==========================================================================

  /**
   * \brief Implementation of the reduced function of the H2C6 mechanism
   *
   */
  template<class T>
  class Red6{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<DType> VDType;
    
    enum{h2c6 = 6};  //! mechanism signature or encoding

    typedef ThermoData4Mechanism<DType,h2c6> DataType;
   

     //THIS MECHANISM IS ALREADY BASED ON CONCENTRATIONS
    typedef Quantity<'c'> QType;
   
    enum{
      reducedDimension = DataType::rednspec,
      fullDimension = DataType::nspec
    };


    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;

    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;

    //typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    //NOTE: This mechanism is isothermal!!!!!
    typedef IsoThermalReverseReactionRates<VType,StoichType> RevType;

    //=========================== DIFFERENCE ==================================
    typedef ReducedIndexInjector<const size_t*, const size_t*> InjectorType;
    //=========================================================================

    typedef BuildChemicalSourceTerm<FwdType,RevType,InjectorType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;
  

    //=========================== DIFFERENCE ==================================
   
#if RED_METH_TO_APPLY == 'M'
    //                             full dim       //red dim         //isothermal
    typedef SpeciesReconstruction<DType,DataType::nspec,reducedDimension,H2MechIn6Species,h2c6> ManifoldType;
#endif
    
#if RED_METH_TO_APPLY == 'L'
    // typedef LebdiedzOptimizationApproach<DataType::nspec,reducedDimension,H2MechIn6Species,h2c6> ManifoldType;
    ADONIS_ERROR(DerivedError, "NOT BEEN USED ANY MORE");
#endif

#if RED_METH_TO_APPLY == 'F'
    typedef VerySimpleSpeciesReconstruction<DType,DataType::nspec,reducedDimension,h2c6> ManifoldType;
#endif 

    typedef Dune::FieldVector<T,DataType::nspec> FullType;
    //=========================================================================

    //************** TODO: change here *****************************************
    enum{iniTrans = 1,
	 coarseEv = 1
    };
    //**************************************************************************



    typedef typename ManifoldType::BTType BTType;
    typedef typename ManifoldType::UType UType;
    
    const BTType& B_T() const {return Mani_.B_T();}
    const UType& U() const {return Mani_.U();}


    //============================= DIFFERENCE ===============================
    Red6(size_t dim = 0, const DType& rpv1 = 0.3, const DType& rpv2 = 0.6,
			    size_t iBd = iniTrans,
	 size_t cE =coarseEv):iszMfilled_(false),rhs_(dim),rpvh2_(rpv1), rpvh2o_(rpv2),countInvocations_(0),initialLayerBound_(iBd), coarseEvals_(cE){
    

      //======================== DIFFERENCE ==================================
      //! get default starting point for SIM
      SIMStartingPoint<DType,h2c6> simstart(DataType::nspec);
      VDType Xpast = simstart.get_point();
      
#if RED_METH_TO_APPLY == 'M'
      VDType Cmass(DataType::mass_balance_matrix(), DataType::mass_balance_matrix()+DataType::rednspec*DataType::nspec),
	bmass(DataType::mass_sum(),DataType::mass_sum()+DataType::rednspec);    
      
      ParameterData PD;
      PD.read_from_file("/home/mfein/MARC++/modred/dat/kn.dat");
      
      DType kn = PD.get_datum<DType>("Delta_t_inner"); //1.25e-04;  //1.e-13; //1.25e-5; //1.e-8; //1.25e-4 


      Mani_.initialize(Xpast,kn,1.e-08,DataType::rpv_index(),Cmass,bmass,'I');
#endif

#if RED_METH_TO_APPLY == 'F'
      VDType Cmass(DataType::mass_balance_matrix(), DataType::mass_balance_matrix()+DataType::rednspec*DataType::nspec),
	bmass(DataType::mass_sum(),DataType::mass_sum()+DataType::rednspec);    
      
      Mani_.initialize(Xpast,DataType::nelem,Cmass,bmass,DataType::rpv_index());
#endif

      eval_.resize(DataType::nspec);
     
      //======================================================================


      //initialize automatization for rhs
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

	//=========================== DIFFERENCE ==============================
	static InjectorType Inj(reducedDimension,DataType::rednreac,DataType::rpv_index(), DataType::rpv_reaction_index());
	//=====================================================================
	BCS_.initialize(&forwardRates,&reverseRates,&Inj);
	

      } //end if dim > 0
    }//end constructor

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    size_t nspec() const{return BCS_.number_of_species();}

    //you might change invocations 
    size_t count_invocations() {return countInvocations_;}


    // const FullType& get_zM() const {
    //   adonis_assert(iszMfilled_); //beware of returning uninitialized data
    //   return zM_;
    // }

    const VDType& simpoint() const {return Mani_.get_z();} 

    //! ============== EVALUATION ==================
    //! the parenthesis-operator
    // template<class RTYPE, class TME>
    // VType& operator()(const TME& kn, const RTYPE& r){ //! r reduced type 
    template<class RTYPE>
    VType& operator()(const RTYPE& r){ //! r reduced type 
      
    //countInvocations_++; //increase counter 
      //convert to right quantity
      //QType::concentration(nspec(),conc_,r,DataType::molar_masses(),DataType::density_of_mixture());   //transform into concentrations. When u = [X] do nothing
      
      
      //========================== DIFFERENCE ================================
      //calculate SIM point -- note that Jochen's tool needs as input a 
      //std::vector. Hence, we copy the stuff (less efficient)
      
      //!only compute manifoldpoint every n times (or so) and also at the 
      //!beginning to overcome possible initial transients
      //if(countInvocations_ <= initialLayerBound_ || countInvocations_%coarseEvals_ == 0){
       
	//! perhaps better choices are possible 
	eval_ = Mani_.get_z();  //zM

	//! compute new SIM point zM here
	//Mani_.evaluate(QType::choose_right_quantity(u,conc_), eval_);
	Mani_.evaluate(QType::choose_right_quantity(r,conc_), eval_);
	
	////get computed full composition from reduced 
	//zM_ = Mani_.get_z();   //that's the SIM point in right concentrations
	UnrollLoop<0,DataType::nspec>::assign_rac2rac(zM_,Mani_.get_z());
	iszMfilled_ = true;
	 //=================== DIFFERENCE ====================
	 //in order to prevent numerical inaccuracies, assign back the identity
	 //values (values of the rpvs)
	 //UnrollLoop<0,reducedDimension>::assign_smaller_2_greater_rac(zM_,r,DataType::rpv_index());
      
	 std::cout << "zM_ = " << zM_ << std::endl;
	 CheckMassBalance<FullType,6>::check_mass(zM_);
	 
	 //====================================================================
      

      //======================== DIFFERENCE ==================================
      
	 for(size_t k = 0; k < reducedDimension; ++k){
	   //! note that the 2nd argument can mean anything since we consider only
	   //! an <B> isothermal </B> mechanism here
	   //! 
	   
	   rhs_[k] = BCS_.net_formation_rate(BCS_.species_index(k),T(1000.),zM_);  //! if u is concentration, choose take u else conc_
	 }
	 
	 //std::cout << "###### Can you read this line, pal ?  #####"<< std::endl;
	 
    //} //end of selected invocations
      
    //QType::back_to_input_quantity(reducedDimension,rhs_,DataType::molar_masses(),DataType::density_of_mixture(),DataType::rpv_index());

      return rhs_;
    }

    /*//! fake time dependent operator
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
      }*/
    //! fake autonomous operator here
    // template<class FV>
    // VType& operator()(const FV& z){   
    //   rhs_ = (*this).operator()(1.0,z); //time is just a dummy value ;)
    //   return rhs_;
    // }


  private:
    bool iszMfilled_;
    VType rhs_, conc_;
    BuilderType BCS_; 

    //=========================== DIFFERENCE ================================
    T rpvh2_, rpvh2o_;
    ManifoldType Mani_;
    FullType zM_; //thats the full composition having 6 species
    //=======================================================================
  
    size_t countInvocations_;
    
    size_t initialLayerBound_, coarseEvals_;
    VType eval_;
  };

  

}
#endif
