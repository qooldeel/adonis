#ifndef MULTICOMPONENT_TRANSPORT_HH
#define  MULTICOMPONENT_TRANSPORT_HH

#include "generaltransport.hh"
#include "diffusion.hh"
#include "fancytransporttmps.hh"
#include "../templatemetaprograms/threeloopsunrolled.hh"

namespace Adonis{

  // · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC
  // · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC
  // · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC · MC
  //!MULTICOMPONEN TRANSPORT ideally with mc. diffusion, thermal conductivity
  //! and thermal diffusion. Until now, only the mc. diffusion is calculated
  template<class TPTDATATYPE, bool BOOL>
  class MultiComponentTransport: public GeneralTransport<TPTDATATYPE,MultiComponentTransport<TPTDATATYPE,BOOL> >{
  public:
    typedef GeneralTransport<TPTDATATYPE,MultiComponentTransport<TPTDATATYPE,BOOL> > BaseClassType;
    typedef typename BaseClassType::DataType DataType;

    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::PointerType PointerType;
    
    //! pure species properties
    typedef PureSpeciesDiffusion<TPTDATATYPE,BOOL> PureDiffType;
   
    enum{DIM = SystemDimension<DataType::nspec,DataType::rednspec,BOOL>::Value};


    MultiComponentTransport(PointerType ptr = 0):PureDiff_(ptr),meanMolecularWeight_(value_type()){
      BaseClassType::ptr_ = ptr;
    }

    MultiComponentTransport(const MultiComponentTransport& MCT):PureDiff_(MCT.PureDiff_),meanMolecularWeight_(MCT.meanMolecularWeight_){
      BaseClassType::ptr_ = MCT.ptr_;
    }
    
    MultiComponentTransport& operator=(const MultiComponentTransport& MCT){
      if(this != &MCT){
	PureDiff_ = MCT.PureDiff_;
	meanMolecularWeight_ = MCT.meanMolecularWeight_;
	BaseClassType::ptr_ = MCT.ptr_;
      }
      return *this;
    }

    void initialize(PointerType ptr){
      BaseClassType::ptr_ = ptr;
      PureDiff_.initialize(ptr);
      meanMolecularWeight_ = value_type();
    }
    
   

    //! solves the multicomponent L matrix system, and computes the 
    //! multicomponent transport coefficients -- NOTE: 
    //!up to now, only multicomp. diffusion is considered.
    template<class T1, class T2, class X>
    void properties(const T1& p, const T2& temp, X& Xfrac){
      //! avoid singularities, due to [KEE] a value of \f$ approx 10^{-12}\f$ works well in practice
      AvoidSingularities<TransportSettings::dontAllowPureSpeciesSituationsToOccur,DIM,value_type>::in_computation(Xfrac.begin());

      //! Alternative: This uses different \f$delta\f$ for each \f$ X_k\f$
      //AvoidSingularities<TransportSettings::dontAllowPureSpeciesSituationsToOccur,N,value_type>::in_computation(Xfrac.begin(),&perturbations_[0]);

      //! compute binary diffusion coefficient matrix, needed for all \f$L\f$
      PureDiff_.compute_binary_diffusion_matrix(p,temp);
      //! compute \f$ L^{00,00} \f$
    
      ComputeTransportProperties<BOOL>::template create_L0000<DataType>(&L_0000_[0],p,temp,Xfrac,DataType::molar_masses(),PureDiff_.get_matrix_ptr());
      
      //! compute inverse of \f$ L^{00,00} \f$ 
      //!copy L0000, since it is needed for the computation of thermal 
      //! conductivities and thermal diffusion later on
      UnrollLoop<0,DIM*DIM>::assign(&inv_[0],&L_0000_[0]); 
      
      int n = DIM,
	lda = n,
	info;

      int ipiv[DIM];
      //LU decomposition of L0000
      F77Type<TypeTraits<value_type>::Value>::GETRF(&n,&n,&inv_[0],&lda,&ipiv[0],&info);
      adonis_assert(info == 0);

      int lwork = n;
      value_type work[DIM];

      //calulate inverse of L0000 based on the previous LU decomposition
      F77Type<TypeTraits<value_type>::Value>::GETRI(&n,&inv_[0],&lda,&ipiv[0],&work[0],&lwork,&info);

      if(info < 0)
	ADONIS_ERROR(LapackError,"The " << info <<"th argument has an illegal value.");
	
      if(info > 0)
	ADONIS_ERROR(LapackError,"U(" << info <<", "<<info<<") = 0 (in C++: U["<<info-1<<"]["<<info-1<<"] = 0). Matrix is singular to working precision ==>\n   No Inverse exists!");

      // //JUST TEST
      // std::cout << "INVERSE OF L0000 = "<< std::endl;
      // for(int i = 0; i < DIM; ++i){
      // 	for(int j = 0; j < DIM; ++j){
      // 	  std::cout << inv_[i*DIM+j] << " ";
      // 	}
      // 	std::cout << std::endl;
      // 	}


      
      //! calculate mean molecular weight for w.r.t. mole fractions <TT> Xfrac</TT>
      meanMolecularWeight_ = value_type(); //reset to zero
      ComputeTransportProperties<BOOL>::template mean_molecular_weight_X<DataType>(meanMolecularWeight_,Xfrac.begin(),DataType::molar_masses());

      //std::cout << "MEAN MOLECULAR WEIGHT = "<< meanMolecularWeight_ << std::endl;

      //! calculate multicomponent diffusion coefficients now
      ComputeTransportProperties<BOOL>::template compute_multicomponent_diffusion_matrix<DataType>(&mcdiff_[0],p,temp,Xfrac,DataType::molar_masses(),meanMolecularWeight_,&inv_[0]);

    }
    

    //! give back some stuff
    PointerType  get_matrix_ptr() {return &PureDiff_[0];}
    PointerType get_L0000_ptr(){return &L_0000_[0];}
    //PointerType get_mc_diffusion_matrix_ptr(){return &mcdiff_[0];}
    PointerType coefficients(){return &mcdiff_[0];}
    

  private:
    PureDiffType PureDiff_;
    value_type mcdiff_[DIM*DIM]; //! m.c.diffusion is a non-symmetric matrix
    value_type L_0000_[DIM*DIM]; //! L0000 is a non-symmetric matrix
    value_type inv_[DIM*DIM];    //! inverse of L0000, non-symm., too
    value_type meanMolecularWeight_;
    value_type perturbations_[DIM];
  };

}

#endif
