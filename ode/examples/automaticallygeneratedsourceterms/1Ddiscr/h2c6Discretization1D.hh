#ifndef ONE_D_SPATIAL_DISCR_OF_FICTICIOUS_H2C6_MECH_HH
#define ONE_D_SPATIAL_DISCR_OF_FICTICIOUS_H2C6_MECH_HH

#include <iostream>
#include <string>

#include "../../../../massactionkinetics/quantity.hh"
#include "../../../../massactionkinetics/qsettings.hh"
#include "../../../../massactionkinetics/tsettings.hh"
#include "../../../../massactionkinetics/reactionrates.hh"
#include "../../../../massactionkinetics/stoichiometry.hh"
#include "../../../../massactionkinetics/thermochemistry.hh"
#include "../../../../massactionkinetics/indexinjector.hh"
#include "../../../../massactionkinetics/buildchemicalrhs.hh"

#include "../../../../expressiontemplates/exprvec.hh"
#include "../../../../massactionkinetics/data/thermochemicaldata.hh"
#include "../../../../common/typeadapter.hh"

#include "../../../../common/smartassign.hh"
#include "../../../../common/universalconstants.hh"
#include "../../../../common/globalfunctions.hh"
#include "../../../../templatemetaprograms/unrollloop.hh"

#include "../h2c6.hh" //source term

#include "../../../../moleculartransport/diffusion.hh"
#include "../../../../massactionkinetics/physicalconstants.hh"
#include "../../../../io/readinparameters.hh"

namespace Adonis{

  /***
   * Note that no <TT> BuildChemicalSourceTerm </TT> is needed here, since we
   * use the corresponding function for the source term directly
   */
  template<class T>
  class H2MechIn6SpeciesSpatial1D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;

    typedef H2MechIn6Species<T> FunType;  

    enum{h2c6 = 6};  //! mechanism signature or encoding

    typedef ThermoData4Mechanism<DType,h2c6> DataType;
    //! false = full diffusion coefficients
    typedef MixtureAveragedDiffusion<DataType,false> DiffusionType;

    //quantity -- here already in concentrations!
    //typedef Quantity<'c'> QType;
   
    
    
    enum{H2 = DataType::H2,
	 H = DataType::H,       
	 O2 = DataType::O2,      
	 O  = DataType::O,       
	 H2O = DataType::H2O,     
	 OH = DataType::OH
    };

    enum{
      space = 51
    };

    //!dim is is for 1D nspec*discpoints
    H2MechIn6SpeciesSpatial1D(size_t dim = 0, const std::string& fname = "statesH2C6.dat"):rhs_(dim),h2c6_((dim != 0) ? DataType::nspec : 0){
      if(dim >=1){
	//! these are in the normal dimensions!
	// QType::resize(DataType::nspec,conc_);
	// QType::resize(DataType::nspec,conc_prev_);
	// QType::resize(DataType::nspec,conc_next_);
	
	species_.resize(DataType::nspec);
	species_prev_.resize(DataType::nspec);
	species_next_.resize(DataType::nspec);

	Diffusion_.resize(DataType::nspec);
      
	delta_x_ = 0.01/(space-1); // 
      	
	//!Diffusion coeffs
	diffcoefs_.resize(DataType::nspec);
	diffcoefs_ = T(1.); //all diffusion coefficients have the same value

	Mad_.initialize(DataType::transport());

	para_.read_from_file(fname);
	rho_ = para_.get_datum<DType>("rho");
	temperature_ = para_.get_datum<DType>("temperature");
	diffprefac_ = para_.get_datum<DType>("diffprefac");
	
	adonis_assert(Abs(diffprefac_) > DType());

	Xfrac_.resize(6);
	p0_ = T();
      }
    }
    
    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const {return dim();}

    size_t nspec() const{return DataType::nspec;}



    template<class E>
    VType& operator()(const E& u){
      //!rho and T are constant throughout
      //T p = rho_*PhysicalConstants<T>::IdealGasConstant*temperature_;   //evaluate pressure 
      p0_ = (PhysicalConstants<T>::IdealGasConstant*temperature_)*total_concentration(0,u); 
      
     
      //PERIODIC BOUNDARIES
      //rhs_[0] = 0.;            //LEFT BDY
      //rhs_[space-1] = 0.;      //RIGHT BDY
      
      for(int k = 0; k < 6; ++k)
	species_[k] = u[OFF(0,k)];
      
      mole_fraction(Xfrac_,0,u);
      //std::cout << "Xfrac_ = "<< Xfrac_ << std::endl;
      
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,temperature_,Xfrac_);
      for(int k=0; k < 6; ++k) //assume u_{-1} = u_0
	Diffusion_[k] = diffcoefs_[k]*(u[OFF(1,k)] - 2.*u[OFF(0,k)] + u[OFF(0,k)] )/ntimes<2>(delta_x_);

      //add diffusion and source term
      Diffusion_ += h2c6_(species_); 
      
      for(int k = 0; k < 6; ++k){
	rhs_[OFF(0,k)] = Diffusion_[k];
	rhs_[OFF(space-1,k)] = rhs_[OFF(0,k)]; //PERIODIC
      }

      //========================== TODO: fix 'em ? =============================
      //rhs_[0] = rhs_[space-1];    //PERIODIC BDY
      //========================================================================

      //T Wbar = T();
      //! note that we already calculate in concentrations, hence the mole fraction is calculated via \f[ X_k = \frac{[X_k]}{\sum_{i=1}^K [X_j]}\f]
      
      //INTERIOR
      for(size_t i = 1; i < space-1; ++i){ 
	p0_ =  (PhysicalConstants<T>::IdealGasConstant*temperature_)*total_concentration(i,u); 
	//extract species at discretisation point i
	for(size_t k = 0; k< DataType::nspec; ++k){
	  //x = k*space + i;
	  species_[k] = u[OFF(i,k)];
	}
	
	mole_fraction(Xfrac_,i,u);

	Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,temperature_,Xfrac_);

	for(size_t k = 0; k< DataType::nspec; ++k){
	  diffcoefs_[k] = Mad_.coefficients()[k];
	}
      
	//std::cout << " DIFF. COEFFS = "<< diffcoefs_ << std::endl;
	//exit(1) just for testing when something seems to be wrong

	//only transform the species when it's not already concentration
	// QType::concentration(nspec(),conc_,species_,DataType::molar_masses(),DataType::density_of_mixture()); 
	// QType::concentration(nspec(),conc_prev_,species_prev_,DataType::molar_masses(),DataType::density_of_mixture()); 
	// QType::concentration(nspec(),conc_next_,species_next_,DataType::molar_masses(),DataType::density_of_mixture()); 

	


	for(size_t k = 0; k < DataType::nspec; ++k){
	  Diffusion_[k] = diffcoefs_[k]*(u[OFF(i+1,k)] - 2.*u[OFF(i,k)] + u[OFF(i-1,k)] )/ntimes<2>(delta_x_); 
	}
	
	

	//add diffusion and source term
	Diffusion_ += h2c6_(species_);
	

	rhs_[i] = Diffusion_[0];
	rhs_[space + i] = Diffusion_[1];
	rhs_[2*space+i] = Diffusion_[2];
	rhs_[3*space+i] = Diffusion_[3];
	rhs_[4*space+i] = Diffusion_[4];
	rhs_[5*space+i] = Diffusion_[5];
	

	//=============== TODO: only in case where u isn't concentration! =====
	//transform back to input quantity
	//for(size_t k = 0; k < DataType::nspec; ++k)
	//  rhs_[k*space + i] = QType::back_to_input_quantity(rhs_[k*space + i],DataType::molar_masses()[k],DataType::density_of_mixture());
	
      }

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
      conc_,
      conc_prev_,
      conc_next_,
      diffcoefs_,
      species_,
      species_prev_,
      species_next_,
      Diffusion_,
      Xfrac_;
    

    DType delta_x_;

    FunType h2c6_; 

    T p0_;

    DiffusionType Mad_;
    T rho_, temperature_,diffprefac_;
    ParameterData para_;
    
    int OFF(int i, int k){
       return (k*space+i);
    }

     T total_concentration(int i , const VType& u){
      T totconc = T();
      for(int k = 0; k < 6; ++k)
	totconc += u[OFF(i,k)];
      if(is_zero(totconc))
	ADONIS_ERROR(ZeroDivision, "Total concentration is ~ 0 at point "<<i<<".");
       
      return totconc;
    }

     VType& mole_fraction(VType& xfrac, int i, const VType& u){
      T totconc = total_concentration(i,u);
      for(int k = 0; k < 6; ++k)
	xfrac[k] = u[OFF(i,k)]/totconc;
      return xfrac;
    }


  };


}

#endif
