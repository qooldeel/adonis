#ifndef REDUCED_FICTICIOUS_H2C6_MECH_HH
#define REDUCED_FICTICIOUS_H2C6_MECH_HH

#include<iostream>
#include<fstream>

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
#include "../../../templatemetaprograms/unrollloop.hh"

//for reduced stuff
#include "../../../modred/manifold.hh"



namespace Adonis{

   //==========================================================================
  //===================== REDUCED SOURCE TERM ================================
  //==========================================================================

  /**
   * \brief Implementation of the reduced function of the H2C6 mechanism
   *
   */
  template<class T>
  class ReducedH2MechIn6Species{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;

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
   
    //                             full dim       //red dim         //isothermal
    typedef ComputeManifold<DType,DataType::nspec,reducedDimension,true> ManifoldType;
    typedef typename ManifoldType::CompositionType CompositionType;
    typedef typename ManifoldType::VecD VecD; //std::vector<double>
    typedef typename ManifoldType::VecS VecS; //std::vector<std::string> 
    typedef typename ManifoldType::CoordinateType CoordinateType; //Dune::FieldVector<int,·>
    
    typedef Dune::FieldVector<T,DataType::nspec> FullType;
    //=========================================================================

    //************** TODO: change here *****************************************
    enum{iniTrans = 1,//1,//10,
	 coarseEv = 1// 1 //35
    };
    //**************************************************************************

    //============================= DIFFERENCE ===============================
    ReducedH2MechIn6Species(size_t dim = 0, const DType& rpv1 = 0.3, const DType& rpv2 = 0.6,
			    size_t iBd = iniTrans,
			    size_t cE =coarseEv):rhs_(dim),rpvh2_(rpv1), rpvh2o_(rpv2),reduced_(reducedDimension),countInvocations_(0),initialLayerBound_(iBd), coarseEvals_(cE){
    

      //======================== DIFFERENCE ==================================
      VType pvvalues(reducedDimension);
      pvvalues[0] = rpv1;
      pvvalues[1] = rpv2;
      
      VecS names(reducedDimension);
      names[0] = "H2";  
      names[1] = "H2O";
      
      CoordinateType rpvindex;     
      rpvindex[0] = 0; rpvindex[1] = 4;
      
      Mani_.activate(names,pvvalues,rpvindex); //first(...) invoked
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
      rpvgrid_.open("h2c6_rpvgrid.dat");
    }//end constructor

    
    ~ReducedH2MechIn6Species(){rpvgrid_.close();}

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    size_t nspec() const{return BCS_.number_of_species();}

    //you might change invocations 
    size_t count_invocations() {return countInvocations_;}


    //! ============== EVALUATION ==================
    //! the parenthesis-operator
    template<class RTYPE, class TME>
    VType& operator()(const TME& kn, const RTYPE& u){ //! note that u is the reduced type now!
      countInvocations_++; //increase counter 

     

	

      //convert to right quantity
      QType::concentration(nspec(),conc_,u,DataType::molar_masses(),DataType::density_of_mixture());   //transform into concentrations. When u = [X] do nothing
      
      
      //========================== DIFFERENCE ================================
      //calculate SIM point -- note that Jochen's tool needs as input a 
      //std::vector. Hence, we copy the stuff (less efficient)
      
      //!only compute manifoldpoint every n times (or so) and also at the 
      //!beginning to overcome possible initial transients
      if(countInvocations_ <= initialLayerBound_ || countInvocations_%coarseEvals_ == 0){
	//! Jochen's tool is currently based on std::vectors. Hence I have to 
	//! copy the stuff to the reduced_ std::vector
         UnrollLoop<0,reducedDimension>::assign_rac2rac(reduced_,QType::choose_right_quantity(u,conc_));  //smart assign ;)
	 
	 std::cout << "INPUT: "; print_all(reduced_.begin(),reduced_.end());
	 Mani_.evaluate(reduced_);
	 
	 ////get computed full composition from reduced 
	 //zM_ = Mani_.get_z();   //that's the SIM point in right concentrations
	 UnrollLoop<0,DataType::nspec>::assign_rac2rac(zM_,Mani_.get_z());
	 
	 //=================== DIFFERENCE ====================
	 //in order to prevent numerical inaccuracies, assign back the identity
	 //values (values of the rpvs)
	 UnrollLoop<0,reducedDimension>::assign_smaller_2_greater_rac(zM_,u,DataType::rpv_index());
      
	 std::cout << "zM_ = " << zM_ << std::endl;
	 
	 //!output
	 for(int r = 0; r < DataType::rednspec; ++r){
	   rpvgrid_ << zM_[DataType::rpv_index()[r]] << " "; 
	 }
	 rpvgrid_ << std::endl;

	 //! check that components don't get too big
	 for(size_t i = 0; i < reducedDimension;++i){
	   if(zM_[i] >= 1.e+12){
	     ADONIS_ERROR(DerivedError, "entries of zM_ are becoming TOO BIG!\n   I doubt this will result in something meaningful");
	   }
	 }
	   
	 //====================================================================
      

      //======================== DIFFERENCE ==================================
      
	 for(size_t k = 0; k < reducedDimension; ++k){
	   //! note that the 2nd argument can mean anything since we consider only
	   //! an <B> isothermal </B> mechanism here
	   //! 
	   
	   rhs_[k] = BCS_.net_formation_rate(BCS_.species_index(k),T(1000.),zM_);  //! if u is concentration, choose take u else conc_
	 }
	 
	 //std::cout << "###### Can you read this line, pal ?  #####"<< std::endl;
	 
      } //end of selected invocations
      
      QType::back_to_input_quantity(reducedDimension,rhs_,DataType::molar_masses(),DataType::density_of_mixture(),DataType::rpv_index());

      return rhs_;
     }

    /*//! fake time dependent operator
    template<class TIME, class FV>
    VType& operator()(const TIME& time, const FV& z){   
      rhs_ = (*this).operator()(z);
      return rhs_;
      }*/
    //! fake autonomous operator here
    template<class FV>
    VType& operator()(const FV& z){   
      rhs_ = (*this).operator()(1.0,z); //time is just a dummy value ;)
      return rhs_;
    }


    const FullType& get_zM() const {return zM_;}

  private:
    VType rhs_, conc_;
    BuilderType BCS_; 

    //=========================== DIFFERENCE ================================
    T rpvh2_, rpvh2o_;
    ManifoldType Mani_;
    VecD reduced_;
    FullType zM_; //thats the full composition having 6 species
    //=======================================================================
  
    size_t countInvocations_;
    
    size_t initialLayerBound_, coarseEvals_;
		
    std::ofstream rpvgrid_;
  };

  

}
#endif
