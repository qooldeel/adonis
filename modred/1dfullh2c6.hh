#ifndef ONE_DIMENSIONAL_FULL_H2C6_MECHANISM_WITH_DIRECT_SOURCE_TERM_HH
#define ONE_DIMENSIONAL_FULL_H2C6_MECHANISM_WITH_DIRECT_SOURCE_TERM_HH

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

namespace Adonis{

  template<class T>
  class MOLH2C61D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type DType;
    typedef ExprTmpl::MyVec<DType> VDType;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef std::size_t SizeType;

    typedef ThermoData4Mechanism<DType,6> DataType;
    //full diffusion
    typedef MixtureAveragedDiffusion<DataType,false> DiffusionType; //cge
    typedef Norm<'2',BaseType> NormType;

    typedef H2MechIn6Species<T> ModelType;

    enum {H2,H,O2,O,H2O,OH};

    MOLH2C61D(SizeType n = 0):model_(6){
      ParameterData PD;
      PD.read_from_file("/home/mfein/MARC++/modred/dat/h2c6.dat");
      nx_ = PD.get_datum<int>("Nx");
      hx_ = (PD.get_datum<DType>("b")-PD.get_datum<DType>("a"))/(nx_-1);
      p0_ = PD.get_datum<DType>("p0");
      rho_ = PD.get_datum<DType>("rho"); //default value
      T0_ = PD.get_datum<DType>("temperature");

      dim_ = 6*nx_;   //6 species, each discretized with nx_ points
      rhs_.resize(dim_);
      u_.resize(dim_);
     
      Conc_.resize(6);
      Xfrac0_.resize(6);  
      Xfrac_.resize(6);  
      Xfrac1_.resize(6);  
      Yfrac0_.resize(6);  
      Yfrac_.resize(6);  
      Yfrac1_.resize(6);  
      Mad0_.initialize(DataType::transport());
      Mad_.initialize(DataType::transport());
      Mad1_.initialize(DataType::transport());
      
     
      wbar_ = 0; 
      
      dotomega_.resize(6);

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

      p1_ = p2_ = T();
      D1_ = D2_ = DType();
    }

    template<class X>
    VType& operator()(const X& concfull){
      //hardcopy
      u_ = concfull;

      //calculate p0
      //p0_ = rho_*PhysicalConstants<T>::IdealGasConstant*T0_;
      p0_ = PhysicalConstants<T>::IdealGasConstant*T0_*total_concentration(0,u_);
      p1_ = PhysicalConstants<T>::IdealGasConstant*T0_*total_concentration(0,u_);
      p2_ = PhysicalConstants<T>::IdealGasConstant*T0_*total_concentration(1,u_);

      //PERIODIC BOUNDARY CONDITIONS
      for(int k = 0; k < 6; ++k){
       	u_[OFF(0,k)] = 	u_[OFF(nx_-1,k)]; 
      }
      //! important here for smooth boundary transitions
      //use u_{-1} = u_0
      mole_fraction(Xfrac0_,0,u_);
      mole_fraction(Xfrac_,0,u_);
      mole_fraction(Xfrac1_,1,u_);
	
      Mad0_.compute_mixture_averaged_diffusion_coefficients(p1_,T0_,Xfrac0_);
      Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,T0_,Xfrac_);
      Mad1_.compute_mixture_averaged_diffusion_coefficients(p2_,T0_,Xfrac1_);

      mass_fraction(Yfrac0_,0,u_);
      mass_fraction(Yfrac_,0,u_);
      mass_fraction(Yfrac1_,1,u_);
	
		
      chemistry(&dotomega_[0],0,u_);   //dotomega_ filled

      //std::cout << "dotomega_ = "; print_all(dotomega_,dotomega_+6);

      rh_[0] = rho(0,u_);
      rh_[1] = rho(0,u_);
      rh_[2] = rho(1,u_);
	
      wb_[0] = wbar(0,u_);
      wb_[1] = wbar(0,u_);
      wb_[2] = wbar(1,u_);

      for(int k = 0; k < 6; ++k){
	  
	rhs_[OFF(0,k)] = 1./DataType::molar_masses()[k]*(0.5*(rh_[2]*DataType::molar_masses()[k]/wb_[2]*Mad1_[k] + rh_[1]*DataType::molar_masses()[k]/wb_[1]*Mad_[k])*(Xfrac1_[k]-Xfrac_[k])  - 0.5*(rh_[1]*DataType::molar_masses()[k]/wb_[1]*Mad_[k] + rh_[0]*DataType::molar_masses()[k]/wb_[0]*Mad0_[k])*(Xfrac_[k]-Xfrac0_[k]) 
							 +  0.5*(rh_[2]*DataType::molar_masses()[k]/wb_[2]*Mad1_[k]*(Xfrac1_[k]-Yfrac1_[k])/p2_ + rh_[1]*DataType::molar_masses()[k]/wb_[1]*Mad_[k]*(Xfrac_[k]-Yfrac_[k])/p0_)*(p2_-p0_) - 0.5*(rh_[1]*DataType::molar_masses()[k]/wb_[1]*Mad0_[k]*(Xfrac_[k]-Yfrac_[k])/p0_ + rh_[0]*DataType::molar_masses()[k]/wb_[0]*Mad0_[k]*(Xfrac0_[k]-Yfrac0_[k])/p1_)*(p0_-p1_))/ntimes<2>(hx_)
	    + dotomega_[k];
	
	
	rhs_[OFF(nx_-1,k)] = rhs_[OFF(0,k)]; //periodicity
      }

      //==============================================================0
      //INTERIOR
      for(int i = 1; i < nx_-1; ++i){
	p0_ = PhysicalConstants<T>::IdealGasConstant*T0_*total_concentration(i,u_);
	p1_ = PhysicalConstants<T>::IdealGasConstant*T0_*total_concentration(i-1,u_);
	p2_ = PhysicalConstants<T>::IdealGasConstant*T0_*total_concentration(i+1,u_);

	mole_fraction(Xfrac0_,i-1,u_);
	mole_fraction(Xfrac_,i,u_);
	mole_fraction(Xfrac1_,i+1,u_);
	
	Mad0_.compute_mixture_averaged_diffusion_coefficients(p1_,T0_,Xfrac0_);
	Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,T0_,Xfrac_);
	Mad1_.compute_mixture_averaged_diffusion_coefficients(p2_,T0_,Xfrac1_);

	mass_fraction(Yfrac0_,i-1,u_);
	mass_fraction(Yfrac_,i,u_);
	mass_fraction(Yfrac1_,i+1,u_);
	
	
	// for(int k = 0; k < 6; ++k) //!the same but auto-generated
	//   Conc_[k] = u_[OFF(i,k)];
	// dotomega_ = model_(Conc_);
	
	chemistry(&dotomega_[0],i,u_);   //dotomega_ filled

	//std::cout << "dotomega_ = "; print_all(dotomega_,dotomega_+6);

	rh_[0] = rho(i-1,u_);
	rh_[1] = rho(i,u_);
	rh_[2] = rho(i+1,u_);
	
	wb_[0] = wbar(i-1,u_);
	wb_[1] = wbar(i,u_);
	wb_[2] = wbar(i+1,u_);

	for(int k = 0; k < 6; ++k){
	  
	  rhs_[OFF(i,k)] = 1./DataType::molar_masses()[k]*(0.5*(rh_[2]*DataType::molar_masses()[k]/wb_[2]*Mad1_[k] + rh_[1]*DataType::molar_masses()[k]/wb_[1]*Mad_[k])*(Xfrac1_[k]-Xfrac_[k])  - 0.5*(rh_[1]*DataType::molar_masses()[k]/wb_[1]*Mad_[k] + rh_[0]*DataType::molar_masses()[k]/wb_[0]*Mad0_[k])*(Xfrac_[k]-Xfrac0_[k]) 
			    +  0.5*(rh_[2]*DataType::molar_masses()[k]/wb_[2]*Mad1_[k]*(Xfrac1_[k]-Yfrac1_[k])/p2_ + rh_[1]*DataType::molar_masses()[k]/wb_[1]*Mad_[k]*(Xfrac_[k]-Yfrac_[k])/p0_)*(p2_-p0_) - 0.5*(rh_[1]*DataType::molar_masses()[k]/wb_[1]*Mad0_[k]*(Xfrac_[k]-Yfrac_[k])/p0_ + rh_[0]*DataType::molar_masses()[k]/wb_[0]*Mad0_[k]*(Xfrac0_[k]-Yfrac0_[k])/p1_)*(p0_-p1_))/ntimes<2>(hx_)
	    + dotomega_[k];
	}
      }//end interior

      return rhs_;
    }

  private:
    ModelType model_;
    int nx_, dim_;
    DType hx_;
    VType rhs_, Xfrac_, Xfrac0_, Xfrac1_,u_, Conc_, Yfrac_,Yfrac0_,Yfrac1_;
    T p0_,p1_,p2_, T0_, rho_, wbar_;
    DiffusionType Mad_, Mad0_, Mad1_;
    
    VType dotomega_;

    DType D1_,D2_;
    
    T three_[3];
    T rh_[3];
    T wb_[3];

    T pk[12]; //reaction rates 
    T rV[12];      //reaction velocities

   

    //! private members
    int OFF(int i, int spec){
      return (i + spec*nx_);
    }

    T rho(int i, const VType& u){
      T r = T();
      for(int k = 0; k < 6; ++k)
	r += u[OFF(i,k)]*DataType::molar_masses()[k];
      if(is_zero(r))
       	ADONIS_ERROR(ZeroDivision, "rho("<<i<<",u) ~ 0");
      return r;
    }
    
    
    T total_concentration(int i , const VType& u){
      T totconc = T();
      for(int k = 0; k < 6; ++k)
	totconc += u[OFF(i,k)];
      if(is_zero(totconc)){
	ADONIS_ERROR(ZeroDivision, "Total concentration is ~ 0 at point "<<i<<".");
      }
      return totconc;
    }

    T wbar(int i, const VType& u){
      return rho(i,u)/total_concentration(i,u);
    }

    VType& mole_fraction(VType& xfrac, int i, const VType& u){
      T totconc = total_concentration(i,u);
      for(int k = 0; k < 6; ++k)
	xfrac[k] = u[OFF(i,k)]/totconc;
      return xfrac;
    }
      
    VType& mass_fraction(VType& yfrac, int i, const VType& u){
      T Rho = rho(i,u);
      for(int k = 0; k < 6; ++k)
	 yfrac[k] = (u[OFF(i,k)]*DataType::molar_masses()[k])/Rho;

       return yfrac;
    }

    template<class IT>
    void chemistry(IT doto, int i, const VType& u){
      // rV[0]  =  pk[0]*z[0];
      // rV[1]  =  pk[6]*z[1]*z[1];
      // rV[2]  =  pk[1]*z[2];
      // rV[3]  =  pk[7]*z[3]*z[3];
      // rV[4]  =  pk[2]*z[4];
      // rV[5]  =  pk[8]*z[1]*z[5];
      // rV[6]  =  pk[3]*z[0]*z[3];
      // rV[7]  =  pk[9]*z[1]*z[5];
      // rV[8]  =  pk[4]*z[2]*z[1];
      // rV[9]  =  pk[10]*z[3]*z[5];
      // rV[10] =  pk[5]*z[0]*z[3];
      // rV[11] =  pk[11]*z[4];

      rV[0]  =  pk[0]*u[OFF(i,0)];
      rV[1]  =  pk[6]*u[OFF(i,1)]*u[OFF(i,1)];
      rV[2]  =  pk[1]*u[OFF(i,2)];
      rV[3]  =  pk[7]*u[OFF(i,3)]*u[OFF(i,3)];
      rV[4]  =  pk[2]*u[OFF(i,4)];
      rV[5]  =  pk[8]*u[OFF(i,1)]*u[OFF(i,5)];
      rV[6]  =  pk[3]*u[OFF(i,0)]*u[OFF(i,3)];
      rV[7]  =  pk[9]*u[OFF(i,1)]*u[OFF(i,5)];
      rV[8]  =  pk[4]*u[OFF(i,2)]*u[OFF(i,1)];
      rV[9]  =  pk[10]*u[OFF(i,3)]*u[OFF(i,5)];
      rV[10] =  pk[5]*u[OFF(i,0)]*u[OFF(i,3)];
      rV[11] =  pk[11]*u[OFF(i,4)];


      
      //the 6 species involved in the reaction:
      doto[0] = - rV[0]  + rV[1]              //H2
      - rV[6]  + rV[7]
      - rV[10] + rV[11];

      doto[1] = 2.0*rV[0] - 2.0*rV[1]          //H
      + rV[4] - rV[5]
      + rV[6] - rV[7]
      - rV[8] + rV[9];
    
      doto[2] = - rV[2] + rV[3]                //O2                 
      - rV[8] + rV[9];
    
      doto[3] = + 2.0*rV[2] - 2.0*rV[3]        //O
      - rV[6]     + rV[7]
      + rV[8]     - rV[9]
      - rV[10]    + rV[11];
    
      doto[4] = - rV[4]  + rV[5]               //H2O
      + rV[10] - rV[11];

      doto[5] = rV[4] - rV[5]                  //OH
      + rV[6] - rV[7]
      + rV[8] - rV[9];
    }
    

  };

} //end namespace 

#endif
