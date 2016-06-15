#ifndef SECOND_APPROX_REDUCED_1_D_SPATIAL_DISCR_OF_FICTICIOUS_H2C6_MECH_HH
#define SECOND_APPROX_REDUCED_1_D_SPATIAL_DISCR_OF_FICTICIOUS_H2C6_MECH_HH

#include <string>

#include "../redh2c6_attempt.hh"
#include "../h2c6.hh"
#include "../../../../derivatives/jacobianofsource.hh"
//#include "../redh2c6.hh"   //include reduced source term


#include "../../../../common/typeadapter.hh"
#include "../../../../common/smartassign.hh"
#include "../../../../moleculartransport/diffusion.hh"
#include "../../../../massactionkinetics/physicalconstants.hh"
#include "../../../../io/readinparameters.hh"

#include "../../../../common/universalconstants.hh"
#include "../../../../common/globalfunctions.hh"

#include "../../../../dunexternal/tangentspace.hh" 
#include "../../../../dunexternal/normalspace.hh"
#include "../../../../templatemetaprograms/unrollloop.hh"
#include "../../../../misc/operations4randomaccesscontainers.hh"

namespace Adonis{

  template<class T>
  class SecondReducedH2MechIn6SpeciesSpatial1D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;
     typedef ExprTmpl::MyVec<DType> VDType;
    //typedef ReducedH2MechIn6Species<T> RedFunType; //with Jochen's tool
    typedef Red6<T>  RedFunType;

    typedef H2MechIn6Species<DType> FullFunType;
    typedef JacS<DType,H2MechIn6Species> JacobianFullType; 
    typedef Dune::FieldMatrix<DType,4,4> LSMtxType;
    typedef Dune::FieldVector<DType,4> DeltaMtxType;
    typedef Dune::FieldMatrix<DType,4,6> NTType;
    typedef Dune::FieldMatrix<DType,2,4> OtherMtxType;
    typedef Dune::FieldMatrix<DType,6,6> JacType;
    typedef TangentSpaceRobust<DType,2,6> TangentSpaceType;
    typedef typename TangentSpaceType::BType BType;

    enum{h2c6 = 6};  //! mechanism signature or encoding

    typedef ThermoData4Mechanism<DType,h2c6> DataType;
   
    //! false = full evaluation of diffusion coeffs
    typedef MixtureAveragedDiffusion<DataType,false> DiffusionFullType;

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
    
    SecondReducedH2MechIn6SpeciesSpatial1D(size_t dim = 0, const DType& rpv1 = 0.3, const DType& rpv2 = 0.6, const std::string& fname = "statesH2C6.dat"):rhs_(dim),redh2c6_((dim != 0) ? reducedDimension : 0, rpv1,rpv2,initTrans,coarseEval),h2c6_((dim != 0) ? 6 : 0){  //invokde reduced function here via initTrans and coarseEval
            
      species_.resize(reducedDimension);
      species_prev_.resize(reducedDimension);
      species_next_.resize(reducedDimension);
      
      Diffusion_.resize(reducedDimension);
      
      delta_x_ = Constant<DType>::rightbdy1D/(space-1);
      
      //!Diffusion coeffs
      diffcoefs_.resize(6);
      diffcoefs_ = T(1.); //all diffusion coefficients have the same value
	
      mad_.initialize(DataType::transport());


      para_.read_from_file(fname);
      rho_ = para_.get_datum<DType>("rho");
      temperature_ = para_.get_datum<DType>("temperature");
      diffprefac_ = para_.get_datum<DType>("diffprefac");
    
      //! some default values to set up Jacobian , e.g. equilibrium point
      ExprTmpl::MyVec<DType> xeq(6);
      xeq <<= 0.27, 0.05, 0.135, 0.02, 0.7, 0.01;
      jac_.set_with_init(6,6,xeq);
    
      Diff_.resize(6);
      z_.resize(6);
      red_.resize(2);
      corr_.resize(2);
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
      DType sumconc = 0.;

      VDType Xfrac(6);
      
      //INTERIOR
      for(size_t i = 1; i < space-1; ++i){ 
	
	//extract species at discretisation point i
	for(size_t k = 0; k< reducedDimension; ++k){
	  x = k*space + i;
	  species_[k] = u[x];
	  species_prev_[k] = u[x-1];
	  species_next_[k] = u[x+1];
	  
	  smart_assign(red_[k],species_[k]);
	  
	}
	
	//!compute SIM POINT
	smart_assign(z_,redh2c6_.simpoint()); //full concentrations
	sumconc = UnrollLoop<0,6>::sum<1>(z_);
	Xfrac = z_/sumconc;  

	
	
	//! full Xfrac required!!
	//compute all mixture averaged diffusion coeffs
	mad_.compute_mixture_averaged_diffusion_coefficients(p,temperature_,Xfrac);
	
	//assign diffusion coeffs
	for(int l = 0; l < 6; ++l){
	  diffcoefs_[l] = diffprefac_* mad_.coefficients()[l];
	  smart_assign(Diff_[l],diffcoefs_[l]); 
	}

	
        std::cout << "FULL DIFF COEFFS = " << diffcoefs_ << std::endl;

	adonis_assert(is_well_defined_value(mad_.coefficients()[0]) &&
		      is_well_defined_value(mad_.coefficients()[4]));
	
	
	//exit(1); //FOR TESTING ONLY !!!!!!!!!!!!!!!!!!

	Diffusion_[0] = diffcoefs_[0]*(species_next_[0] - 2.*species_[0] + species_prev_[0])/ntimes<2>(delta_x_);
	Diffusion_[1] = diffcoefs_[4]*(species_next_[1] - 2.*species_[1] + species_prev_[1])/ntimes<2>(delta_x_);

	

	//add diffusion and source term
	Diffusion_ += redh2c6_(species_); //reduced composition

	//! ===================== compute correction =========================
	//! problem: correction may be big an mass conservation is violated
	std::cout << "Start>>===================== CPA part ========================"<<std::endl;
	fill_from_rac(Jac_,jac_.jacobian(z_));
	TangentSpaceType Tspace(transpose(redh2c6_.B_T())); //reference only
        NTType NT = normalspace_transposed<2>(Tspace.compute(red_,z_));
	NtJU_ = NT*Jac_*redh2c6_.U();
	std::cout << "N^TJU = "<< std::endl << NtJU_ << std::endl;
	NT *= -1.; //-N^T
	Diff_ += h2c6_(z_); //! D{z} + S(z)
	matrix_vector_product(delta_,NT,Diff_); //delta stores rhs now
	delta_ =  solve_rank_deficient_ls(NtJU_,delta_);
	BtJU_ = redh2c6_.B_T()*Jac_*redh2c6_.U();
	matrix_vector_product(corr_,BtJU_,delta_);
	
	normalize(corr_); //NORMALIZE corr
	std::cout << "corr_ = "<< std::endl<< corr_ << std::endl;

	std::cout << "End<<======================================================"<< std::endl;
	//!===================================================================


	rhs_[i] = Diffusion_[0] + corr_[0];
	rhs_[space+i] = Diffusion_[1] + corr_[1];

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
    FullFunType h2c6_;
    DiffusionFullType mad_;

    T rho_, temperature_,diffprefac_;
    ParameterData para_;
    JacobianFullType jac_;

    LSMtxType NtJU_;
    DeltaMtxType delta_;
    VDType Diff_, z_, red_;
    VDType corr_;
    OtherMtxType BtJU_;
    JacType Jac_;
  };

}

#endif
