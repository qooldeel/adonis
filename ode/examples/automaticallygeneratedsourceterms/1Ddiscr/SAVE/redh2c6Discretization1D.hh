#ifndef REDUCED_1_D_SPATIAL_DISCR_OF_FICTICIOUS_H2C6_MECH_HH
#define REDUCED_1_D_SPATIAL_DISCR_OF_FICTICIOUS_H2C6_MECH_HH

#include <string>

#include "../redh2c6_attempt.hh"
//#include "../redh2c6.hh"   //include reduced source term


#include "../../../../common/typeadapter.hh"
#include "../../../../moleculartransport/diffusion.hh"
#include "../../../../massactionkinetics/physicalconstants.hh"
#include "../../../../io/readinparameters.hh"

#include "../../../../common/universalconstants.hh"
#include "../../../../common/globalfunctions.hh"

namespace Adonis{

  template<class T>
  class ReducedH2MechIn6SpeciesSpatial1D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;
    //typedef ReducedH2MechIn6Species<T> RedFunType; //with Jochen's tool
    typedef Red6<T>  RedFunType;

    enum{h2c6 = 6};  //! mechanism signature or encoding

    typedef ThermoData4Mechanism<DType,h2c6> DataType;
    //! true = reduction
    typedef MixtureAveragedDiffusion<DataType,true> DiffusionType;

    //quantity -- here already in concentrations!
    //typedef Quantity<'c'> QType;  //not needed for H2C6 example at all

    enum{
      space = Constant<DType>::spacePoints,
      
      reducedDimension = DataType::rednspec,
      fullDimension = DataType::nspec
    };

    //====================== TODO: Set evaluation bounds here =================
    enum{
      initTrans = 1,
      coarseEval = 1
    };
    //=========================================================================
    
    ReducedH2MechIn6SpeciesSpatial1D(size_t dim = 0, const DType& rpv1 = 0.3, const DType& rpv2 = 0.6, const std::string& fname = "statesH2C6.dat"):rhs_(dim),redh2c6_((dim != 0) ? reducedDimension : 0, rpv1,rpv2,initTrans,coarseEval){  //invokde reduced function here via initTrans and coarseEval
            
      species_.resize(reducedDimension);
      species_prev_.resize(reducedDimension);
      species_next_.resize(reducedDimension);
      
      Diffusion_.resize(reducedDimension);
      
      delta_x_ = Constant<DType>::rightbdy1D/(space-1);
      
      //!Diffusion coeffs
      diffcoefs_.resize(reducedDimension);
      diffcoefs_ = T(1.); //all diffusion coefficients have the same value
	
      mixavgdiff_.initialize(DataType::transport());


      para_.read_from_file(fname);
      rho_ = para_.get_datum<DType>("rho");
      temperature_ = para_.get_datum<DType>("temperature");
      diffprefac_ = para_.get_datum<DType>("diffprefac");
    }

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    template<class E>
    VType& operator()(const E& u){
      //!rho and T are constant throughout
      T p = rho_*PhysicalConstants<T>::IdealGasConstant*temperature_;   //evaluate pressure 

      size_t x;
      //BOUNDARIES
      // rhs_[0] = 0.;            //LEFT BDY
      // rhs_[space-1] = 0.;      //RIGHT BDY
      
      rhs_[0] = rhs_[space-1];    //PERIODIC BDY
      //rhs_[0] = rhs_[space-1] = ?????

      //T Wbar = T();
      
      //! note that we already calculate in concentrations, hence the mole fraction is calculated via \f[ X_k = \frac{[X_k]}{\sum_{i=1}^K [X_j]}\f]
      T sumconc = T();

      VType Xfrac(reducedDimension);
      
      //INTERIOR
      for(size_t i = 1; i < space-1; ++i){ 
	
	//extract species at discretisation point i
	sumconc = T(); //reset
	for(size_t k = 0; k< reducedDimension; ++k){
	  x = k*space + i;
	  species_[k] = u[x];
	  species_prev_[k] = u[x-1];
	  species_next_[k] = u[x+1];
	  
	  
	  Xfrac[k] = u[x];
	  //!compute mean molecular weight given specific moles:
	  //Wbar += u[x]; //note that u[x] = Yk/Wk
	  sumconc += u[x];  
	}
	
	//Wbar = 1./Wbar;
	//Xfrac *= Wbar;
	adonis_assert(Abs(sumconc) > UniversalConstants<T>::aboutZero);
	Xfrac /= sumconc;  

	std::cout << "Xfrac = "<< Xfrac << std::endl;
	
	mixavgdiff_.compute_mixture_averaged_diffusion_coefficients(p,temperature_,Xfrac);
	
	
	diffcoefs_[0] = diffprefac_* mixavgdiff_.coefficients()[0];
	diffcoefs_[1] = diffprefac_* mixavgdiff_.coefficients()[1];

	//diffcoefs_[0] = diffcoefs_[1] = 1.; //TEST ONLY

	 std::cout << "DIFF COEFFS = " << diffcoefs_ << std::endl;

	 adonis_assert(is_well_defined_value(mixavgdiff_.coefficients()[0]) &&
		       is_well_defined_value(mixavgdiff_.coefficients()[1]));


	 //exit(1); //FOR TESTING ONLY !!!!!!!!!!!!!!!!!!

	Diffusion_[0] = diffcoefs_[0]*(species_next_[0] - 2.*species_[0] + species_prev_[0])/ntimes<2>(delta_x_);
	Diffusion_[1] = diffcoefs_[1]*(species_next_[1] - 2.*species_[1] + species_prev_[1])/ntimes<2>(delta_x_);

	

	//add diffusion and source term
	Diffusion_ += redh2c6_(species_); 

	rhs_[i] = Diffusion_[0];
	rhs_[space+i] = Diffusion_[1];

      }//end interior

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
      diffcoefs_,
      species_,
      species_prev_,
      species_next_,
      Diffusion_;

    DType delta_x_;

    RedFunType redh2c6_;
    DiffusionType mixavgdiff_;

    T rho_, temperature_,diffprefac_;
    ParameterData para_;
  };

}

#endif
