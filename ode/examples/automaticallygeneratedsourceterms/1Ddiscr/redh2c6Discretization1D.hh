#ifndef REDUCED_1_D_SPATIAL_DISCR_OF_FICTICIOUS_H2C6_MECH_HH
#define REDUCED_1_D_SPATIAL_DISCR_OF_FICTICIOUS_H2C6_MECH_HH

#include <string>

#include "../redh2c6_attempt.hh"



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
      space = 51,
      
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
      
      delta_x_ = 0.01/(space-1);
      
      //!Diffusion coeffs
      diffcoefs_.resize(reducedDimension);
      diffcoefs_ = T(1.); //all diffusion coefficients have the same value
	
      Mad_.initialize(DataType::transport());
    
      Xfrac_.resize(reducedDimension);

      para_.read_from_file(fname);
      rho_ = para_.get_datum<DType>("rho");
      p0_ = 101325.; //default
      temperature_ = para_.get_datum<DType>("temperature");
      diffprefac_ = para_.get_datum<DType>("diffprefac");
    }

    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    template<class E>
    VType& operator()(const E& u){
      //!rho and T are constant throughout
      // p0_ = rho_*PhysicalConstants<T>::IdealGasConstant*temperature_;   //evaluate pressure 
      p0_ =  (PhysicalConstants<T>::IdealGasConstant*temperature_)*total_concentration(0,u);   

      //PERIODIC BOUNDARIES
      for(int k = 0; k < reducedDimension; ++k)
	species_[k] = u[OFF(0,k)];
      
      mole_fraction(Xfrac_,0,u);
      std::cout << "Xfrac_ = "<< Xfrac_ << std::endl;
      
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,temperature_,Xfrac_);
      for(int k=0; k < reducedDimension; ++k) //assume u_{-1} = u_0
	Diffusion_[k] = diffcoefs_[k]*(u[OFF(1,k)] - 2.*u[OFF(0,k)] + u[OFF(0,k)] )/ntimes<2>(delta_x_);

      //add diffusion and source term
      Diffusion_ += redh2c6_(species_); 
      
      for(int k = 0; k < reducedDimension; ++k){
	rhs_[OFF(0,k)] = Diffusion_[k];
	rhs_[OFF(space-1,k)] = rhs_[OFF(0,k)]; //PERIODIC
      }

      //rhs_[0] = rhs_[space-1];    //OLD:PERIODIC BDY
     
      
      //INTERIOR
      for(size_t i = 1; i < space-1; ++i){ 
	 p0_ =  (PhysicalConstants<T>::IdealGasConstant*temperature_)*total_concentration(i,u); 
	//extract species at discretisation point i
	for(size_t k = 0; k< reducedDimension; ++k){
	  //x = k*space + i;
	  species_[k] = u[OFF(i,k)];
	}

	mole_fraction(Xfrac_,i,u);
	std::cout << "Xfrac_ = "<< Xfrac_ << std::endl;
	
	Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,temperature_,Xfrac_);
	
	
	diffcoefs_[0] = Mad_[0];
	diffcoefs_[1] = Mad_[1];

	 std::cout << "DIFF COEFFS = " << diffcoefs_ << std::endl;

	 adonis_assert(is_well_defined_value(Mad_[0]) &&
		       is_well_defined_value(Mad_[1]));


	 //exit(1); //FOR TESTING ONLY !!!!!!!!!!!!!!!!!!

	 for(int k=0; k < reducedDimension; ++k)
	   Diffusion_[k] = diffcoefs_[k]*(u[OFF(i+1,k)] - 2.*u[OFF(i,k)] + u[OFF(i-1,k)] )/ntimes<2>(delta_x_); //overwrite Diffusion again
	

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

    VType Xfrac_;
    DType delta_x_;

    RedFunType redh2c6_;
    DiffusionType Mad_;

    T rho_, p0_, temperature_,diffprefac_;
    ParameterData para_;

     int OFF(int i, int k){
       return (k*space+i);
    }

     T total_concentration(int i , const VType& u){
      T totconc = T();
      for(int k = 0; k < reducedDimension; ++k)
	totconc += u[OFF(i,k)];
      if(is_zero(totconc))
	ADONIS_ERROR(ZeroDivision, "Total concentration is ~ 0 at point "<<i<<".");
       
      return totconc;
    }

     VType& mole_fraction(VType& xfrac, int i, const VType& u){
      T totconc = total_concentration(i,u);
      for(int k = 0; k < reducedDimension; ++k)
	xfrac[k] = u[OFF(i,k)]/totconc;
      return xfrac;
    }

  };

}

#endif
