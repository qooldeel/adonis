#ifndef ONE_DIMENSIONAL_REDUCED_H2C6_MECHANISM_WITH_DIRECT_SOURCE_TERM_HH
#define ONE_DIMENSIONAL_REDUCED_H2C6_MECHANISM_WITH_DIRECT_SOURCE_TERM_HH

#include "../expressiontemplates/exprvec.hh"
#include "../io/readinparameters.hh"
#include "../massactionkinetics/data/thermochemicaldata.hh"
#include "../moleculartransport/diffusion.hh"
#include "../ode/sourceterms.hh"
#include "../linalg/linearsystemsolvers.hh"
#include "../derivatives/jacobianofsource.hh"
#include "../templatemetaprograms/unrollloop.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../ode/examples/automaticallygeneratedsourceterms/h2c6.hh"


#include "speciesreconstruction.hh"   //RED
#include "../ode/sourceterms.hh"      //RED

//#include "1dfullh2c6.hh"   //RED

namespace Adonis{

  template<class T>
  class ReducedMOLH2C61D{    //RED
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<DType> VDType;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef std::size_t SizeType;

    typedef ThermoData4Mechanism<DType,6> DataType;
    
    typedef ExprTmpl::MyVec<ExprTmpl::MyVec<T> > TableType;
    
    typedef MixtureAveragedDiffusion<DataType,true> DiffusionType; //RED
    typedef SpeciesReconstruction<DType,6,2,H2Combustion6Spex,6> Reconstructor; //RED

    typedef Norm<'2',BaseType> NormType;

   
    enum{nchem=2}; //RED

    ReducedMOLH2C61D(SizeType n = 0){
      ParameterData PD;
      PD.read_from_file("/home/mfein/MARC++/modred/dat/h2c6.dat");
      nx_ = PD.get_datum<int>("Nx");
      hx_ = (PD.get_datum<DType>("b")-PD.get_datum<DType>("a"))/(nx_-1);
      p0_ = PD.get_datum<DType>("p0");
      rho_ = PD.get_datum<DType>("rho"); //default value
      T0_ = PD.get_datum<DType>("temperature");

      dim_ = nchem*nx_;   //RED species, each discretized with nx_ points
      rhs_.resize(dim_);
      u_.resize(dim_);
     
      Conc_.resize(nchem);
      Xfrac0_.resize(nchem);  
      Xfrac_.resize(nchem);  
      Xfrac1_.resize(nchem);  
      Mad0_.initialize(DataType::transport());
      Mad_.initialize(DataType::transport());
      Mad1_.initialize(DataType::transport());
      
     
      wbar_ = 0; 
      
      dotomega_.resize(nchem);
      zM_.resize(6);           //RED
      VDType Xpast(6);
      
      Xpast <<=                                  //RED
      	PD.get_datum<DType>("H2"),             //RED
      	PD.get_datum<DType>("H"),              //RED  
      	PD.get_datum<DType>("O2"),             //RED
      	PD.get_datum<DType>("O"),              //RED
      	PD.get_datum<DType>("H2O"),            //RED
      	PD.get_datum<DType>("OH");             //RED

      
      Tabularization_.resize(nx_); 
      VType IV(6*nx_);      //some initial values for tabularization
      ExprTmpl::MyVec<T> u0(6), u1(6);
    
      u0 <<= 0.455, 0.779, 0.237, 0.363, 0.148, 0.015;
      u1 <<= 0.2, 0.95, 0.31, 0.03, 0.3, 0.05;
      
      DType x(0);
      for(int i = 0; i < nx_ ;++i){
	(Tabularization_[i]).resize(6);  //store full state
	x = i*hx_;  //CAUTION: previously: x += i*hx
	if(0.002 <= x && x <= 0.008){
	  for(int k = 0; k < 6; ++k)
	    IV[k*nx_+i] = u0[k];
	}
	if(x < 0.002 || x > 0.008){
	  for(int k = 0; k < 6; ++k)
	    IV[k*nx_+i] = u1[k];
	}
      
	for(int k = 0; k < 6; ++k){
	  Tabularization_[i][k] = IV[k*nx_+i];
	}
      }


      // SIMStartingPoint<DType,6> simstart(DataType::nspec);
      // VDType Xpast = simstart.get_point();
      // zM_ = zM0_;   //RED 

      VDType Cmass(DataType::mass_balance_matrix(), DataType::mass_balance_matrix()+DataType::rednspec*DataType::nspec),
	bmass(DataType::mass_sum(),DataType::mass_sum()+DataType::rednspec); 

      Recon_.initialize(Xpast,PD.get_datum<DType>("kn"),          //RED
			PD.get_datum<BaseType>("newtTol"),DataType::rpv_index(),
			Cmass, bmass,PD.get_datum<char>("whichTimestepper"), 
			PD.get_datum<int>("maxIt"),PD.get_datum<BaseType>("nmtol"), PD.get_datum<BaseType>("tol"));
		       

      //fixed reaction rates
      pk[0] = 2.0;
      pk[1] = 1.0;
      pk[2] = 1.0;
      pk[3] = 1000.0;
      pk[4] = 1000.0;
      pk[5] = 100.0;
      pk[6] = 216.0;
      pk[7] = 337.5;
      pk[8] = 1400.0;
      pk[9] = 10800.0;
      pk[10] = 33750.0;
      pk[11] = 0.7714;

      species_.resize(2);
    }

    template<class X>
    VType& operator()(const X& redconc){
      //hardcopy
      u_ = redconc;

      p0_ = (PhysicalConstants<T>::IdealGasConstant*T0_)*total_concentration(0,u_);            //total_concentration(zM_); 
      
     
     
     
      //!PERIODIC BOUNDARY CONDITIONS
      for(int k = 0; k < nchem; ++k){
       	u_[OFF(0,k)] = 	u_[OFF(nx_-1,k)]; 
      }
      // ! important here for smooth boundary transitions
      // !use u_{-1} = u_0
      mole_fraction(Xfrac0_,0,u_); //same as Xfrac0
      mole_fraction(Xfrac_,0,u_);
      mole_fraction(Xfrac1_,1,u_);
      
      Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,T0_,Xfrac0_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,T0_,Xfrac_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(p0_,T0_,Xfrac1_);
      
      reconstruct_full_state(zM_,0,u_);          //RED
      reduced_chemistry(&dotomega_[0],0,zM_);   //RED
      std::cout << "zM = "<< zM_ << std::endl; 

      for(int k = 0; k < nchem; ++k){

	rhs_[OFF(0,k)] = (0.5*(Mad1_[k]+Mad_[k])*(u_[OFF(1,k)]- u_[OFF(0,k)]) - 0.5*(Mad_[k]+Mad0_[k])*(u_[OFF(0,k)]-u_[OFF(0,k)]))/ntimes<2>(hx_)//Mad_[k]*(u_[OFF(i+1,k)] - 2.*u_[OFF(i,k)] + u_[OFF(i-1,k)])/ntimes<2>(hx_) 
	  + dotomega_[k];
	rhs_[OFF(nx_-1,k)] = rhs_[OFF(0,k)]; //periodicity
      }
      

      //INTERIOR
      for(int i = 1; i < nx_-1; ++i){
	
	//! calculate full state and assign it to zM_ 
	reconstruct_full_state(zM_,i,u_);          //RED
	// for(size_t k = 0; k< nchem; ++k){
	//   species_[k] = u_[OFF(i,k)];
	// }

	// zM_ = Recon_.get_z();       //RED
	// Recon_.evaluate(species_,zM_);   //RED
	// UnrollLoop<0,6>::assign_rac2rac(zM_,Recon_.get_z());   //RED
	std::cout << "zM = "<< zM_ << std::endl;           //RED
      
	p0_ = PhysicalConstants<T>::IdealGasConstant*T0_*total_concentration(i,u_);             //total_concentration(zM_); //


	
      	mole_fraction(Xfrac0_,i-1,u_);
      	mole_fraction(Xfrac_,i,u_);
      	mole_fraction(Xfrac1_,i+1,u_);
	
      	Mad0_.compute_mixture_averaged_diffusion_coefficients(p0_,T0_,Xfrac0_);
      	Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,T0_,Xfrac_);
      	Mad1_.compute_mixture_averaged_diffusion_coefficients(p0_,T0_,Xfrac1_);


      	reduced_chemistry(&dotomega_[0],i,zM_);   //RED

      	//std::cout << "dotomega_ = "; print_all(dotomega_,dotomega_+nchem);

      	for(int k = 0; k < nchem; ++k){
	  std::cout << "D(i+1)_"<<k<<" = "<< Mad1_[k] << "     D(i)_"<<k<<" = " << Mad_[k] << "    D(i-1)_"<<k<< Mad0_[k] << std::endl;
      	  rhs_[OFF(i,k)] = (0.5*(Mad1_[k]+Mad_[k])*(u_[OFF(i+1,k)]- u_[OFF(i,k)]) - 0.5*(Mad_[k]+Mad0_[k])*(u_[OFF(i,k)]-u_[OFF(i-1,k)]))/ntimes<2>(hx_)//Mad_[k]*(u_[OFF(i+1,k)] - 2.*u_[OFF(i,k)] + u_[OFF(i-1,k)])/ntimes<2>(hx_) 
      	    + dotomega_[k];
      	}
      }

      return rhs_;
    }
 

  private:
    int nx_, dim_;
    DType hx_;
    VType rhs_, Xfrac_, Xfrac0_, Xfrac1_,u_, Conc_;
    T p0_, T0_, rho_, wbar_;
    DiffusionType Mad_, Mad0_, Mad1_;
    VType species_;

    VType dotomega_;
    VType zM_, zM0_, z_;              //RED
    T pk[12]; //reaction rates 
    T rV[12];      //reaction velocities

    Reconstructor Recon_;
    
    TableType Tabularization_;

    //! private members
    int OFF(int i, int spec){
      return (i + spec*nx_);
    }

    T rho(int i, const VType& u){
      T r = T();
      for(int k = 0; k < nchem; ++k)
	r += u[OFF(i,k)]*DataType::molar_masses()[k];
      if(is_zero(r))
       	ADONIS_ERROR(ZeroDivision, "rho("<<i<<",u) ~ 0");
      return r;
    }
    
    
    T total_concentration(int i , const VType& u){
      T totconc = T();
      for(int k = 0; k < nchem; ++k)
	totconc += u[OFF(i,k)];
      if(is_zero(totconc))
	ADONIS_ERROR(ZeroDivision, "Total concentration is ~ 0 at point "<<i<<".");
       
      return totconc;
    }

    //full composition
    T total_concentration(const VType& zm){
      T tc = T();
      for(SizeType k = 0; k < zm.size(); ++k)
	tc += zm[k];
       if(is_zero(tc))
	ADONIS_ERROR(ZeroDivision, "Total concentration is ~ 0");
       return tc;
    }

    T wbar(int i, const VType& u){
      return rho(i,u)/total_concentration(i,u);
    }

    VType& mole_fraction(VType& xfrac, int i, const VType& u){
      T totconc = total_concentration(i,u);
      for(int k = 0; k < nchem; ++k)
	xfrac[k] = u[OFF(i,k)]/totconc;
      return xfrac;
    }
      
    VType& reconstruct_full_state(VType& zm, int i, const VType& u){
      
      for(size_t k = 0; k< nchem; ++k){
	species_[k] = u[OFF(i,k)];
      }

      // //! symmetric problem -- take mirror values
      // int half = (nx_-1)/2,
      // 	nix(0);
      // if(i >= half){
      // 	nix = 2*half-i; // same as: half - dist;   //= 2*half -i
      // 	for(size_t k = 0; k< nchem; ++k){
      // 	  species_[k] = u[OFF(nix,k)];
      // 	}
      // }

      zm = Recon_.evaluate(species_,Tabularization_[i]);   //RED
      Tabularization_[i] = zm; //overwrite with new value
      //zm = Recon_.get_z();       //RED
      //Recon_.evaluate(species_,zm);
      //UnrollLoop<0,6>::assign_rac2rac(zm,Recon_.get_z());   //RED
     

      return zm;
    }
    

    template<class IT, class ZM>
    void reduced_chemistry(IT doto, int i, const ZM& z){   //RED
     
      rV[0]  =  pk[0]*z[0];
      rV[1]  =  pk[6]*z[1]*z[1];
      
      rV[4]  =  pk[2]*z[4];
      rV[5]  =  pk[8]*z[1]*z[5];
      rV[6]  =  pk[3]*z[0]*z[3];
      rV[7]  =  pk[9]*z[1]*z[5];
     
      rV[10] =  pk[5]*z[0]*z[3];
      rV[11] =  pk[11]*z[4];


      
      //the 6 species involved in the reaction:
      doto[0] = - rV[0]  + rV[1]              //H2
      - rV[6]  + rV[7]
      - rV[10] + rV[11];

    
      doto[1] = - rV[4]  + rV[5]               //H2O
      + rV[10] - rV[11];
    }
    

  };

} //end namespace 

#endif
