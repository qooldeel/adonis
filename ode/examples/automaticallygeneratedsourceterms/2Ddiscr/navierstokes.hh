#ifndef REACTIVE_NAVIER_STOKES_EQUATIONS_HH
#define REACTIVE_NAVIER_STOKES_EQUATIONS_HH

#include "../../../../massactionkinetics/reactionrates.hh"
#include "../../../../massactionkinetics/stoichiometry.hh"
#include "../../../../massactionkinetics/thermochemistry.hh"
#include "../../../../massactionkinetics/indexinjector.hh"
#include "../../../../massactionkinetics/buildchemicalrhs.hh"
#include "../../../../massactionkinetics/physicalconstants.hh"
#include "../../../../massactionkinetics/eos.hh"
#include "../../../../massactionkinetics/data/thermochemicaldata.hh"


#include "../../../../expressiontemplates/exprvec.hh"
#include "../../../../templatemetaprograms/unrollloop.hh"
#include "../../../../common/typeadapter.hh"
#include "../../../../common/smartassign.hh"
#include "../../../../common/adonisassert.hh"
#include "../../../../io/readinparameters.hh"
#include "../../../../misc/misctmps.hh"

//transport quantities
#include "../../../../moleculartransport/viscosity.hh"
#include "../../../../moleculartransport/conductivity.hh"

#include "../../../additional/parabolicvelocityprofile.hh"
#include "../../../additional/boundarytemperaturedistribution.hh"

//reaction mechanisms
//#include "../o3.hh"

namespace Adonis{

  /**
   * \brief Reactive Navier Stokes on a rectangular domain [a,b] x [c,d],
   * discretized with finite differences
   */
  template<class T>
  class ReactiveNavierStokes2D{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;    
    typedef typename TypeAdapter<T>::Type DType;
    //typedef ExprTmpl::MyVec<DType> VDouType;
    typedef std::size_t SizeType;

    //================== TODO: set mechanism here ============================
    enum{thermoMech = 3};   //3 = O3, 9 = H2Gri
    
    //========================================================================

    typedef ThermoData4Mechanism<DType,thermoMech> DataType;  

    //needed for chemistry
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;
    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,CommonIndexInjector> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    typedef EquationOfState<'i'> EosType;

    ReactiveNavierStokes2D(SizeType n = 0):rhs_(n),Nx_(0),Ny_(0),npt_(0),a_(DType()),hx_(DType()),hy_(DType()),Tignition_(DType()), Tmin_(DType()),halfdiam_(DType()),cp_(T()),cp1_(T()),cp2_(T()),p0_(T()),density_(T()),vx_(T()),vy_(T()),paravelo_(T()),Tin_(T()),qheat_(T()),qh_(T()),mb_(T()),mmw_(T()),rho_(T()), qheatReac_(T()){
      if(n > 0){
	C_.resize(DataType::nspec);
	chem_.resize(DataType::nspec+1); //last entry holds temperature

	ParameterData PD;
	//! CAUTION: from file you intend to invoke 'Reactive..2D' !!!!!
	PD.read_from_file("examples/automaticallygeneratedsourceterms/2Ddiscr/2Dsettings.dat");
	DType a = PD.get_datum<DType>("a"),
	  b = PD.get_datum<DType>("b"),
	  c = PD.get_datum<DType>("c"),
	  d = PD.get_datum<DType>("d");
	
	a_ = a;
	Nx_ = PD.get_datum<SizeType>("Nx");
	Ny_ = PD.get_datum<SizeType>("Ny");
	
	npt_ = Nx_*Ny_;
	hx_ = (b-a)/(Nx_-1);
	DType diam = d-c;
	hy_ = diam/(Ny_-1);
	halfdiam_ = diam/2.;

	

	Tignition_ = PD.get_datum<DType>("T_ignition");
	Tmin_ = PD.get_datum<DType>("T_min");
	

	p0_ = PD.get_datum<T>("pconst");
	vx_ = PD.get_datum<T>("v1_in");
	vy_ = PD.get_datum<T>("v2_in");
	Tin_ = PD.get_datum<T>("T_in");
	
	bdytemp_.initialize(Tmin_,Tignition_,PD.get_datum<DType>("hottest_x"));
	v1bdy_.resize(Ny_);

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
	static RevType reverseRates(&forwardRates,&stoichmatrix,&nasa);
	//this is just the normal index
	static CommonIndexInjector Cii(DataType::nspec,DataType::nreac);
	BCS_.initialize(&forwardRates,&reverseRates,&Cii);
      

	//!transport quantities
	Mav_.initialize(DataType::transport());
	Mac_.initialize(DataType::transport());
	Mad_.initialize(DataType::transport());
	Mad2_.initialize(DataType::transport());
	Mad3_.initialize(DataType::transport());
	Mad4_.initialize(DataType::transport());
	Mad5_.initialize(DataType::transport());

	X_.resize(DataType::nspec);
	X1_.resize(DataType::nspec);
	X2_.resize(DataType::nspec);
	X3_.resize(DataType::nspec);
	X4_.resize(DataType::nspec);
	X5_.resize(DataType::nspec);

	Yfrac_.resize(DataType::nspec);
	Yfrac2_.resize(DataType::nspec);
	Yfrac3_.resize(DataType::nspec);
	Yfrac4_.resize(DataType::nspec);
	Yfrac5_.resize(DataType::nspec);

	Y_.resize(DataType::nspec);
	Y1_.resize(DataType::nspec);
      }    
    } //end constructor
    

    SizeType dim() const {return rhs_.size();}
    SizeType domain_dim() const {return dim();}

    //! OPERATOR -- input primite variables [rho,v1,v2,T,Y1,...,YK]
    template<class E>
    VType& operator()(const E& prim){
      
      
      //! ++++++++++++++  BOUNDARY CONDITIONS -- in primitive variables ++++++
      //! down: (i,0), top: (i,Ny_-1), left: (0,j), right: (Nx_-1,j)
      
      // UP/DOWN · UP/DOWN · UP/DOWN · UP/DOWN · UP/DOWN · UP/DOWN · UP/DOWN
      // UP/DOWN · UP/DOWN · UP/DOWN · UP/DOWN · UP/DOWN · UP/DOWN · UP/DOWN
      DType x = DType(); 
      T Wbar, bdyTemp = T();
      

      //DOWN and UP BDY -- Dirichlet (not necessarily zero-rho and species)
      for(SizeType i = 0; i < Nx_; ++i){
	x = a_ + i*hx_;
        
	//symmetric domain
	//!species
	for(SizeType k = 0; k < DataType::nspec-1; ++k){
	  Y_[k] = prim[offset(i,0, 4+k)];
	  Y1_[k] = prim[offset(i,Ny_-1, 4+k)];
	}
	//!preserve mass
	Y_[DataType::nspec-1] = 1. - UnrollLoop<0,DataType::nspec-1>::template sum<1>(Y_);
	Y1_[DataType::nspec-1] = 1. - UnrollLoop<0,DataType::nspec-1>::template sum<1>(Y1_);

	bdyTemp = bdytemp_(x);

	//! calculate chemical source, including temperature
	csource1_ = chemical_source_term(Y_,bdyTemp);
	csource2_ = chemical_source_term(Y1_,bdyTemp);
	
	//species
	for(SizeType k = 0; k < DataType::nspec; ++k){
	  rhs_[offset(i,0, 4+k)] = csource1_[k];
	  rhs_[offset(i,Ny_-1, 4+k)] = csource2_[k];
	}

	Wbar = T(); //reset
	mmw_ = T();
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar,Y_.begin(),DataType::molar_masses());
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar,Y1_.begin(),DataType::molar_masses());
	

        //! velocity: assume zero velocity at boundary
	rhs_[offset(i,0, 1)] = rhs_[offset(i,Ny_-1, 1)] = 0.; //! v1
	rhs_[offset(i,0, 2)] = rhs_[offset(i,Ny_-1, 2)] = 0.; //! v2
     
	//! temperature 
	rhs_[offset(i,0,3)] = csource1_[DataType::nspec];
	rhs_[offset(i,Ny_-1,3)] = csource2_[DataType::nspec];
	//!rho -- since v1 = v2 = 0 here
	rhs_[offset(i,0, 0)] = rhs_[offset(i,Ny_-1, 0)] = 0.;

      }
      
      //LEFT · LEFT · LEFT · LEFT · LEFT · LEFT · LEFT · LEFT · LEFT · LEFT
      //LEFT · LEFT · LEFT · LEFT · LEFT · LEFT · LEFT · LEFT · LEFT · LEFT
      //L·L·L·L·L   LEFT -- INLET (cold bdy). Dirichlet  L·L·L·L·L 
      //! velocity profile at left bdy
      DType r = -halfdiam_;
      for(SizeType j = 0; j < Ny_; ++j){
	v1bdy_[j] = zero(parabolic_inlet_velocity(r,halfdiam_,vx_));
	r += hy_;
      }

      T Wbar2 = T(), Wbar3 = T(), Wbar4 = T(), Wbar5 = T();
      
      for(SizeType j = 1; j < Ny_-1; ++j){  //0 and Ny_-1 already covered by $\Gamma_{\mathrm{wall}$
	
      
	Wbar = Wbar2 = Wbar3 = Wbar4 = Wbar5 = T();
      	for(SizeType k = 0; k < DataType::nspec; ++k){
	  Yfrac_[k] = DataType::default_values()[k];
	  Yfrac2_[k] = prim[offset(1,j, 4+k)];
	  Yfrac3_[k] = DataType::default_values()[k];
	  Yfrac4_[k] = DataType::default_values()[k];
	  Yfrac5_[k] =  prim[offset(2,j, 4+k)];
	}
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar,&Yfrac_[0],DataType::molar_masses());
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar2,&Yfrac2_[0],DataType::molar_masses());
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar3,&Yfrac3_[0],DataType::molar_masses());
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar4,&Yfrac4_[0],DataType::molar_masses());
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar5,&Yfrac5_[0],DataType::molar_masses());
       

	// compute mole fractions
	for(SizeType k = 0; k < DataType::nspec; ++k){
	  X_[k] = Yfrac_[k]*Wbar/DataType::molar_masses()[k];
	  X2_[k] = Yfrac2_[k]*Wbar2/DataType::molar_masses()[k];
	  X3_[k] = Yfrac3_[k]*Wbar3/DataType::molar_masses()[k];
	  X4_[k] = Yfrac4_[k]*Wbar4/DataType::molar_masses()[k];
	  X5_[k] = Yfrac5_[k]*Wbar5/DataType::molar_masses()[k];
	}

	//calculate rho at boundary from EOS for ideal gas
	rho_ = p0_*Wbar/(PhysicalConstants<T>::IdealGasConstant*Tin_);

	//rho -- v2 = 0, hence second rhs term vanishes
	rhs_[offset(0,j, 0)] = -((prim[offset(1,j, 0)]*prim[offset(1,j, 1)] - rho_*v1bdy_[j])/hx_);   



	//! chemistry with diffusion
	Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,Tin_,X_);
	Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(1,j, 3)],X2_);

	csource1_ = chemical_source_term(Yfrac_,Tin_);
	
	cp_ = cp1_ = qheat_ = qh_ = T();
	for(SizeType k = 0; k < DataType::nspec; ++k){
	  

	  cp_ += (BCS_.C_p(k,Tin_)/DataType::molar_masses()[k]*Yfrac_[k]);
	  cp1_ += (BCS_.C_p(k,prim[offset(1,j, 3)])/DataType::molar_masses()[k]*Yfrac2_[k]);

	  qheat_ += (BCS_.H_T(k,Tin_)*csource1_[k]);

	  //! note that at the left bdy X_k is constant along the y-axis, assume
	  //! T_i-1,j = T_i,j. This leads to vanishing terms
	  qh_ += (BCS_.C_p(k,Tin_)/DataType::molar_masses()[k]*(-1.)*(rho_*DataType::molar_masses()[k]/Wbar*Mad_[k])*(0.5*(X2_[k]+X_[k])*(prim[offset(1,j, 3)]-Tin_))/ntimes<2>(hx_) 
);
	}
	

	
	//! v1
	// how to discr. {mu_i+1/2,j(v1_i+1,j - v1_i,j) - mu_i-1/2,j(v1_i,j - v1_i-1,j)}/hx^2 at bdy ?? for the second term, take mu_i,j and v1 like in 1st
	rhs_[offset(0,j, 1)] = -1./rho_*(ntimes<2>(prim[offset(1,j, 1)]) - rho_*(ntimes<2>(v1bdy_[j]))/hx_ -(2.*2./3.*0.5*(Mav_.compute_mixture_averaged_viscosity(Tin_,X_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(1,j, 3)],X2_))*(prim[offset(1,j, 1)] - v1bdy_[j]) - (Mav_.compute_mixture_averaged_viscosity(Tin_,X_))*(prim[offset(1,j, 1)] - v1bdy_[j])/ntimes<2>(hx_) + (Mav_.compute_mixture_averaged_viscosity(Tin_,X3_)*(v1bdy_[j+1]-v1bdy_[j]) - Mav_.compute_mixture_averaged_viscosity(Tin_,X3_)*(v1bdy_[j]-v1bdy_[j-1]))/ntimes<2>(hy_))
);
	//!always 0 velocity in y-direction
	rhs_[offset(0,j, 2)] = 0.; 


	//!species	
	for(SizeType k = 0; k < DataType::nspec; ++k){
	  //! some terms just vanish due to similar values
	  rhs_[offset(0,j, 4+k)] = -1./rho_*( (prim[offset(1,j, 0)]*prim[offset(1,j, 4+k)]*prim[offset(1,j, 1)] - rho_*Yfrac_[k]*v1bdy_[j])/hx_  -(0.5*(prim[offset(1,j, 0)]*DataType::molar_masses()[k]/Wbar2*Mad2_[k] + rho_*DataType::molar_masses()[k]/Wbar*Mad_[k])*(X2_[k]-X_[k]))/ntimes<2>(hx_) - csource1_[k]*DataType::molar_masses()[k]);								      
	}

	//!temperature -- the dy(lambda dy T) vanishes since Tin is constant 
	//! in y-direction. The same applies to X_k
	//assume T_i,j = T_i-1,j. This leads to a cancellation of terms
	rhs_[offset(0,j, 3)] = -1./(rho_*cp_)*( (prim[offset(1,j, 0)]*cp1_*prim[offset(1,j, 3)]*prim[offset(1,j, 1)] - rho_*cp_*Tin_*v1bdy_[j])/hx_ -(0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(1,j, 0)],p0_,prim[offset(1,j, 3)],X2_) + Mac_.compute_mixture_averaged_conductivity(rho_,p0_,Tin_,X_))*(prim[offset(1,j, 3)] - Tin_))/ntimes<2>(hx_) - (qh_ + qheat_) );

      }


      //R·R·R·R·R· RIGHT -- OUT zero NEUMANN for all species, temperature, velocities R·R·R·R·R
      //!--nearly the same as input, except we use the 0-Neumann condition to
      //! get data: u_Nx,j = u_Nx-2,j (2nd order FD gradient), u_Nx,j is ghost node
      Wbar2 = Wbar3 = Wbar4 = Wbar5 = T();
      for(SizeType j = 1; j < Ny_-1; ++j){  //0 and Ny_-1 already covered by $\Gamma_{\mathrm{wall}$
	
	//calculate rho at boundary from EOS for ideal gas
	rho_ = p0_*Wbar/(PhysicalConstants<T>::IdealGasConstant*prim[offset(Nx_-1,j, 3)]);

	//rho, always: v2 = 0 
	rhs_[offset(0,j, 0)] = -((prim[offset(Nx_-2,j, 0)]*prim[offset(Nx_-2,j, 1)] - rho_*prim[offset(Nx_-1,j, 1)])/(2.*hx_));   

	

	Wbar = Wbar2 = Wbar3 = Wbar4 = Wbar5 = T();
      	for(SizeType k = 0; k < DataType::nspec; ++k){
	  Yfrac_[k] = prim[offset(Nx_-1,j, 4+k)];
	  Yfrac2_[k] = prim[offset(Nx_-2,j, 4+k)]; //ghost point from Neumann
	  Yfrac3_[k] = prim[offset(Nx_-1,j-1, 4+k)];
	  Yfrac4_[k] = prim[offset(Nx_-1,j+1, 4+k)];
	}
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar,&Yfrac_[0],DataType::molar_masses());
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar2,&Yfrac2_[0],DataType::molar_masses());
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar3,&Yfrac3_[0],DataType::molar_masses());
	UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(Wbar4,&Yfrac4_[0],DataType::molar_masses());

       
	
	csource1_ = chemical_source_term(Yfrac_,prim[offset(Nx_-1,j, 3)]);
	
	cp_ = cp1_ = qheat_ = qh_ = T();
	for(SizeType k = 0; k < DataType::nspec; ++k){
	  X_[k] = Yfrac_[k]*Wbar/DataType::molar_masses()[k];
	  X2_[k] = Yfrac2_[k]*Wbar2/DataType::molar_masses()[k];
	  X3_[k] = Yfrac3_[k]*Wbar3/DataType::molar_masses()[k];
	  X4_[k] = Yfrac4_[k]*Wbar4/DataType::molar_masses()[k];
	 

	  cp_ += (BCS_.C_p(k,prim[offset(Nx_-1,j, 3)])/DataType::molar_masses()[k]*Yfrac_[k]);
	  cp1_ += (BCS_.C_p(k,prim[offset(Nx_-2,j, 3)])/DataType::molar_masses()[k]*Yfrac2_[k]);

	  qheat_ += (BCS_.H_T(k,prim[offset(Nx_-1,j, 3)])*csource1_[k]);
	}
	
	//! chemistry with diffusion
	Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(Nx_-1,j, 3)],X_);
	Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(Nx_-2,j, 3)],X2_);
	Mad3_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(Nx_-1,j-1, 3)],X3_);
	Mad4_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(Nx_-1,j+1, 3)],X4_);

	//compute penultimate term of thermal energy:
	for(SizeType k = 0; k < DataType::nspec; ++k){
	  qh_ += (BCS_.C_p(k,prim[offset(Nx_-1,j, 3)])/DataType::molar_masses()[k]*(-1.)*(rho_*DataType::molar_masses()[k]/Wbar*Mad_[k])*( (0.5*(X2_[k]+X_[k])*(prim[offset(Nx_-2,j, 3)]-prim[offset(Nx_-1,j, 3)]) - 0.5*(X_[k]+X2_[k])*(prim[offset(Nx_-1,j, 3)] - prim[offset(Nx_-1,j, 3)]))/ntimes<2>(hx_) + (0.5*(X4_[k]+X_[k])*(prim[offset(Nx_-1,j+1, 3)] - prim[offset(Nx_-1,j, 3)]) - 0.5*(X_[k]+X3_[k])*(prim[offset(Nx_-1,j, 3)] - prim[offset(Nx_-1,j-1, 3)]))/ntimes<2>(hy_)) 
		  );
	}


	
	//! v1, note that we exploit v2 = 0 everywhere. Thus terms involving
	//! v2 can be neglected
	rhs_[offset(Nx_-1,j, 1)] = -1./rho_*( rho_*(ntimes<2>(prim[offset(Nx_-2,j, 1)]) - ntimes<2>(prim[offset(Nx_-1,j, 1)]))/(2.*hx_) -( (2.*2./3.*0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(Nx_-2,j, 3)],X2_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(Nx_-1,j, 3)],X_))*(prim[offset(Nx_-2,j, 1)] - prim[offset(Nx_-1,j, 1)]) - 0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(Nx_-1,j, 3)],X_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(Nx_-2,j, 3)],X2_))*(prim[offset(Nx_-1,j, 1)] - prim[offset(Nx_-2,j, 1)]))/ntimes<2>(hx_) 
																     + (0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(Nx_-1,j+1, 3)],X4_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(Nx_-1,j, 3)],X_))*(prim[offset(Nx_-1,j+1, 1)]-prim[offset(Nx_-1,j, 1)]) - 0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(Nx_-1,j, 3)],X_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(Nx_-1,j-1, 3)],X3_))*(prim[offset(Nx_-1,j, 1)]-prim[offset(Nx_-1,j-1, 1)]))/ntimes<2>(hy_))
);
	//!always 0 velocity in y-direction
	rhs_[offset(0,j, 2)] = 0.; 


	//!species	
	for(SizeType k = 0; k < DataType::nspec; ++k){
	  //! some terms just vanish due to similar values
	  rhs_[offset(0,j, 4+k)] = -1./rho_*( (prim[offset(Nx_-2,j, 0)]*prim[offset(Nx_-2,j, 4+k)]*prim[offset(Nx_-2,j, 1)] - rho_*Yfrac_[k]*prim[offset(Nx_-1,j, 1)])/(2.*hx_)  
					      -( (0.5*(prim[offset(Nx_-2,j, 0)]*DataType::molar_masses()[k]/Wbar2*Mad2_[k] + rho_*DataType::molar_masses()[k]/Wbar*Mad_[k])*(X2_[k]-X_[k]) + 0.5*(rho_*DataType::molar_masses()[k]/Wbar*Mad_[k] + prim[offset(Nx_-2,j, 0)]*DataType::molar_masses()[k]/Wbar2*Mad2_[k])*(X_[k] - X2_[k]))/ntimes<2>(hx_) 
						 + (0.5*(prim[offset(Nx_-1,j+1, 0)]*DataType::molar_masses()[k]/Wbar4*Mad4_[k] + rho_*DataType::molar_masses()[k]/Wbar*Mad_[k])*(X4_[k]-X_[k]) + 0.5*(rho_*DataType::molar_masses()[k]/Wbar*Mad_[k] + prim[offset(Nx_-1,j-1, 0)]*DataType::molar_masses()[k]/Wbar3*Mad3_[k])*(X_[k]-X3_[k]))/ntimes<2>(hy_))- csource1_[k]*DataType::molar_masses()[k]);								      
	}

	//!temperature  -- note v2 = 0 i.e. some terms involving partial 
	//! derivatives of v2 vanish
	rhs_[offset(0,j, 3)] = -1./(rho_*cp_)*( (prim[offset(Nx_-2,j, 0)]*cp2_*prim[offset(Nx_-2,j, 3)]*prim[offset(Nx_-2,j, 1)] - rho_*cp_*prim[offset(Nx_-1,j, 3)]*prim[offset(Nx_-1,j, 1)])/(2.*hx_) -( (0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(Nx_-2,j, 0)],p0_,prim[offset(Nx_-2,j, 3)],X2_) + Mac_.compute_mixture_averaged_conductivity(rho_,p0_,prim[offset(Nx_-1,j, 3)],X_))*(prim[offset(Nx_-2,j, 3)] -prim[offset(Nx_-1,j, 3)])
 - 0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(Nx_-1,j, 0)],p0_,prim[offset(Nx_-1,j, 3)],X_) + Mac_.compute_mixture_averaged_conductivity(prim[offset(Nx_-2,j, 0)],p0_,prim[offset(Nx_-2,j, 3)],X2_))*(prim[offset(Nx_-1,j, 3)] - prim[offset(Nx_-2,j, 3)]))/ntimes<2>(hx_) 
																									   + ( 0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(Nx_-1,j+1, 0)],p0_,prim[offset(Nx_-1,j+1, 3)],X4_) + Mac_.compute_mixture_averaged_conductivity(rho_,p0_,prim[offset(Nx_-1,j, 3)],X_))*(prim[offset(Nx_-1,j+1, 3)] - prim[offset(Nx_-1,j, 3)]) 
																									       - 0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(Nx_-1,j, 0)],p0_,prim[offset(Nx_-1,j, 3)],X_) + Mac_.compute_mixture_averaged_conductivity(rho_,p0_,prim[offset(Nx_-1,j-1, 3)],X3_))*(prim[offset(Nx_-1,j, 3)] - prim[offset(Nx_-1,j-1, 3)]) )/ntimes<2>(hy_) ) - (qh_ + qheat_) );
	
      }



      //! ------------------- END BOUNDARY ----------------------------------
      

      
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //! ++++++++++++++ INTERIOR ++++++++++++++++++++++++++++++++++++++++++++
      //!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //! note: discretization in conservative variables
      for(SizeType i = 1; i < Nx_-1; ++i){
	for(SizeType j = 1; j < Ny_-1; ++j){
	  //!density
	  rhs_[offset(i,j, 0)] = -( (prim[offset(i+1,j, 0)]*prim[offset(i+1,j, 1)] - (prim[offset(i,j, 0)]*prim[offset(i,j, 1)]))/hx_ + (prim[offset(i,j+1, 0)]*prim[offset(i,j+1, 2)] - prim[offset(i,j, 0)]*prim[offset(i,j, 2)])/hy_ );
	 
	  

	  //! Note since we calculate in mass fractions, we have to transform
	  //! them into mole fractions, i.e. Xk = Yk·Wbar/Wk
	  Wbar = T();
	  Wbar2 = Wbar3 = Wbar4 = Wbar5 = T();
	  cp_ = cp1_ = cp2_ = T();
	  qheat_ = T();
	  qh_ = T();

	  for(SizeType k = 0; k < DataType::nspec; ++k){
	    Wbar += prim[offset(i,j, 4+k)]/DataType::molar_masses()[k];
	    Wbar2 += prim[offset(i-1,j, 4+k)]/DataType::molar_masses()[k];
	    Wbar3 += prim[offset(i+1,j, 4+k)]/DataType::molar_masses()[k];
	    Wbar4 += prim[offset(i,j-1, 4+k)]/DataType::molar_masses()[k];
	    Wbar5 += prim[offset(i,j+1, 4+k)]/DataType::molar_masses()[k];
	  
	    Y_[k] = prim[offset(i,j, 4+k)]; 
	  }
	  Wbar = 1./Wbar;
	  Wbar2 = 1./Wbar2;
	  Wbar3 = 1./Wbar3;
	  Wbar4 = 1./Wbar4;
	  Wbar5 = 1./Wbar5;

	  //! ENFORCE MASS BALANCE
	  mb_ = UnrollLoop<0,DataType::nspec-1>::template sum<1>(Y_);
	  Y_[DataType::nspec-1] = 1. - mb_;


	  for(SizeType k = 0; k < DataType::nspec; ++k){
	    //each Xl contains the species at the specific discr. point 
	    X1_[k] = prim[offset(i,j, 4+k)]*Wbar/DataType::molar_masses()[k];
	    //! x-direction
	    X2_[k] = prim[offset(i-1,j, 4+k)]*Wbar2/DataType::molar_masses()[k];
	    X3_[k] = prim[offset(i+1,j, 4+k)]*Wbar3/DataType::molar_masses()[k];
	    //! y-direction
	    X4_[k] = prim[offset(i,j-1, 4+k)]*Wbar4/DataType::molar_masses()[k];
	    X5_[k] = prim[offset(i,j+1, 4+k)]*Wbar5/DataType::molar_masses()[k];
	  
	    //! concentrations, rho·Yk/Wk  Yk = prim[offset(i,j, 4+k)]
	    C_[k] = prim[offset(i,j, 0)]*Y_[k]/DataType::molar_masses()[k];
	  
	    //! compute dot omega
	    dotOmega_[k] = BCS_.net_formation_rate(BCS_.species_index(k),prim[offset(i,j, 3)],C_);
	    
	    //! this is needed in the computation of the thermal energy term
	    cp_ += (BCS_.C_p(k,prim[offset(i,j, 3)])/DataType::molar_masses()[k]*prim[offset(i,j, 4+k)]);
	    cp1_ += (BCS_.C_p(k,prim[offset(i+1,j, 3)])/DataType::molar_masses()[k]*prim[offset(i,j, 4+k)]);

	    cp2_ += (BCS_.C_p(k,prim[offset(i,j+1, 3)])/DataType::molar_masses()[k]*prim[offset(i,j, 4+k)]);

	    qheat_ += (BCS_.H_T(k,prim[offset(i,j, 3)])*dotOmega_[k]);
	  
	    
	  }
	    
	 
	  
	  Mad_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(i,j, 3)],X1_);
	  Mad2_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(i-1,j, 3)],X2_);
	  Mad3_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(i+1,j, 3)],X3_);
	  Mad4_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(i,j-1, 3)],X4_);
	  Mad5_.compute_mixture_averaged_diffusion_coefficients(p0_,prim[offset(i,j+1, 3)],X5_);
	  


	  //! compute \f$ \sum_{k=1}^K c_{pk}j_k\cdot \nabla T\f$ 
	  for(SizeType k = 0; k < DataType::nspec; ++k){
	    qh_ += ( BCS_.C_p(k,prim[offset(i,j, 3)])/DataType::molar_masses()[k]*prim[offset(i,j, 0)]*((-1.)*DataType::molar_masses()[k]/Wbar*Mad_[k])*
		     //! dx Xk dx T
		     ( (0.5*(X1_[k]+X3_[k])*(prim[offset(i+1,j, 3)] - prim[offset(i,j, 3)]) - 0.5*(X1_[k]+X2_[k])*(prim[offset(i,j, 3)] - prim[offset(i-1,j, 3)]))/ntimes<2>(hx_) 
		     //dy Xk dy T
		       + (0.5*(X1_[k]+X5_[k])*(prim[offset(i,j+1, 3)] - prim[offset(i,j, 3)]) - 0.5*(X4_[k]+X1_[k])*(prim[offset(i,j, 3)] - prim[offset(i,j-1, 3)]))/ntimes<2>(hy_) )
		     );
	  }
	  	
	  //rhs_[offset(i,j, 0)] = EosType::density(p0_,Wbar,prim[offset(i,j, 3)]);
  

	  //! momentum rho·v1
	  rhs_[offset(i,j, 1)] = -( (prim[offset(i+1,j, 0)]*ntimes<2>(prim[offset(i+1,j, 1)]) - prim[offset(i,j, 0)]*ntimes<2>(prim[offset(i,j, 1)]))/hx_ + (prim[offset(i,j+1, 0)]*prim[offset(i,j+1, 1)]*prim[offset(i,j+1, 2)] - prim[offset(i,j, 0)]*prim[offset(i,j, 1)]*prim[offset(i,j, 2)])/hy_  
			      //!dx tau_11:
			      //! 2·dy mu dx v1
			      - 2./3.*(2*(0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(i+1,j, 3)],X3_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j, 3)],X1_))*(prim[offset(i+1,j, 1)] - prim[offset(i,j, 1)]) -  0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(i-1,j, 3)],X2_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j, 3)],X1_))*(prim[offset(i,j, 1)] - prim[offset(i-1,j, 1)]))/ntimes<2>(hx_)
			      //! mixed derivatives involving v2 !
				       - (Mav_.compute_mixture_averaged_viscosity(prim[offset(i+1,j, 3)],X3_)*(prim[offset(i+1,j+1, 2)] - prim[offset(i+1,j-1, 2)]) - Mav_.compute_mixture_averaged_viscosity(prim[offset(i-1,j, 3)],X2_)*(prim[offset(i-1,j+1, 2)] - prim[offset(i-1,j-1, 2)]))/(4.*hx_*hy_) )
			      //dy tau_12:
			      //2·dy mu dy v1 
			      -((0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j+1, 3)],X5_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j, 3)],X1_))*(prim[offset(i,j+1, 1)] - prim[offset(i,j, 1)]) - 0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j, 3)],X1_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j-1, 3)],X4_))*(prim[offset(i,j, 1)] - prim[offset(i,j-1, 1)]))/ntimes<2>(hy_)  //dy mu dy v1
			      //dy mu dx v2
				+  (Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j+1, 3)],X5_)*(prim[offset(i+1,j+1, 2)] - prim[offset(i-1,j+1, 2)]) - Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j-1, 3)],X4_)*(prim[offset(i+1,j-1, 3)] - prim[offset(i-1,j-1, 2)]))/(4.*hx_*hy_))
			      ) ;
	  
	  //! momentum rho·v2
	  rhs_[offset(i,j, 2)] = -( (prim[offset(i+1,j, 0)]*prim[offset(i+1,j, 2)]*prim[offset(i+1,j, 1)] - prim[offset(i,j, 0)]*prim[offset(i,j, 2)]*prim[offset(i,j, 1)])/hx_  + (prim[offset(i,j+1, 0)]*ntimes<2>(prim[offset(i,j+1, 2)]) - prim[offset(i,j, 0)]*ntimes<2>(prim[offset(i,j, 2)]))/hy_
			      //!dx tau_21:
			      //!dx mu dy v1
			      -( (Mav_.compute_mixture_averaged_viscosity(prim[offset(i+1,j, 3)],X3_)*(prim[offset(i+1,j+1, 1)] - prim[offset(i+1,j-1, 1)]) - Mav_.compute_mixture_averaged_viscosity(prim[offset(i-1,j, 3)],X2_)*(prim[offset(i-1,j+1, 1)] - prim[offset(i-1,j-1, 1)]))/(4.*hx_*hy_) 
				 //!dx mu dx v2
				 +  (0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(i+1,j, 3)],X3_)+ Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j, 3)],X1_))*(prim[offset(i+1,j, 2)] - prim[offset(i,j, 2)]) - 0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j, 3)],X1_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(i-1,j, 3)],X2_))*(prim[offset(i,j, 2)] - prim[offset(i-1,j, 2)]))/ntimes<2>(hx_) )
			      //! dy tau_22:
			      //! 2·dy mu dy v2
			      -2./3.*(2.*(0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j+1, 3)],X5_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j, 3)],X1_))*(prim[offset(i,j+1, 2)] - prim[offset(i,j, 2)]) - 0.5*(Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j-1, 3)],X4_) + Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j, 3)],X1_))*(prim[offset(i,j, 2)] - prim[offset(i,j-1, 2)]))/ntimes<2>(hy_) 
				      //! dy mu dx v1
				      - (Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j+1, 3)],X5_)*(prim[offset(i+1,j+1, 1)] -prim[offset(i-1,j+1, 1)]) - Mav_.compute_mixture_averaged_viscosity(prim[offset(i,j-1, 3)],X4_)*(prim[offset(i+1,j-1, 1)] - prim[offset(i-1,j-1, 1)]))/(4.*hx_*hy_) )
			      ); 

	  //coefficient between nodes
	  cp1_ = 0.5*(cp1_ + cp_);
	  cp2_ = 0.5*(cp2_ + cp_);
	  
	  //! thermal energy
	  rhs_[offset(i,j, 3)] = -( (cp1_*(prim[offset(i+1,j, 0)]*prim[offset(i+1,j, 3)]*prim[offset(i+1,j, 1)] - prim[offset(i,j, 0)]*prim[offset(i,j, 3)]*prim[offset(i,j, 1)]))/hx_  //! dx(rho·cp·T·v1)
			      + (cp2_*(prim[offset(i,j+1, 0)]*prim[offset(i,j+1, 3)]*prim[offset(i,j+1, 2)] - prim[offset(i,j, 0)]*prim[offset(i,j, 3)]*prim[offset(i,j, 2)]))/hy_ //! dy(rho·cp·T·v2) 
			      //! dx lambda dx T
			      -(0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(i+1,j, 0)],p0_,prim[offset(i+1,j, 3)],X3_) + Mac_.compute_mixture_averaged_conductivity(prim[offset(i,j, 0)],p0_,prim[offset(i,j, 3)],X1_))*(prim[offset(i+1,j, 3)] - prim[offset(i,j, 3)]) - 0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(i,j, 0)],p0_,prim[offset(i,j, 3)],X1_) + Mac_.compute_mixture_averaged_conductivity(prim[offset(i-1,j, 0)],p0_,prim[offset(i-1,j, 3)],X2_))*(prim[offset(i,j, 3)] - prim[offset(i-1,j, 3)]))/ntimes<2>(hx_)
			      //! dy lambda dx T
			      -(0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(i,j+1, 0)],p0_,prim[offset(i,j+1, 3)],X5_) + Mac_.compute_mixture_averaged_conductivity(prim[offset(i,j, 0)],p0_,prim[offset(i,j, 3)],X1_))*(prim[offset(i,j+1, 3)] - prim[offset(i,j, 3)]) - 0.5*(Mac_.compute_mixture_averaged_conductivity(prim[offset(i,j, 0)],p0_,prim[offset(i,j, 3)],X1_) + Mac_.compute_mixture_averaged_conductivity(prim[offset(i,j-1, 0)],p0_,prim[offset(i,j-1, 3)],X4_))*(prim[offset(i,j, 3)] - prim[offset(i,j-1, 3)]))/ntimes<2>(hy_)
				    - (qh_ + qheat_) //CAUTION, cf. [KEE, p. 115, eq. 3.213]: MINUS or PLUS
				    );

	  //SPECIES -- in conservative form
	  for(SizeType k = 0; k < DataType::nspec; ++k){
	    rhs_[offset(i,j, 4+k)] = -( (prim[offset(i+1,j, 0)]*prim[offset(i+1,j,4+k)]*prim[offset(i+1,j, 1)] - prim[offset(i,j, 0)]*prim[offset(i,j,4+k)]*prim[offset(i,j, 1)])/hx_ //dx (rho Yk v1)
				  +  (prim[offset(i,j+1, 0)]*prim[offset(i,j+1,4+k)]*prim[offset(i,j+1, 2)] - prim[offset(i,j, 0)]*prim[offset(i,j,4+k)]*prim[offset(i,j, 2)])/hy_ //dy (rho Yk v2)
				  -(1./DataType::molar_masses()[k]*(0.5*(prim[offset(i+1,j, 0)] + prim[offset(i,j, 0)])/(0.5*(Wbar3+Wbar))*0.5*(Mad_[k]+Mad3_[k])*(X3_[k]-X1_[k]) - 0.5*(prim[offset(i,j, 0)] + prim[offset(i-1,j, 0)])/(0.5*(Wbar+Wbar2))*0.5*(Mad_[k]+Mad2_[k])*(X1_[k]-X2_[k]))/ntimes<2>(hx_)    //dx(rho Wk/Wbar Dmix dx Xk)
				    + 1./DataType::molar_masses()[k]*(0.5*(prim[offset(i,j+1, 0)] + prim[offset(i,j, 0)])/(0.5*(Wbar5+Wbar))*0.5*(Mad5_[k]+Mad_[k])*(X5_[k]-X1_[k]) - 0.5*(prim[offset(i,j, 0)] + prim[offset(i,j-1, 0)])/(0.5*(Wbar+Wbar4))*0.5*(Mad_[k]+Mad4_[k])*(X1_[k] - X4_[k]))/ntimes<2>(hy_) )   //dy(rho Wk/Wbar Dmix dy Xk)
				  - dotOmega_[k]*DataType::molar_masses()[k]
				  );
 	  }
	  
	  
	  //!BACK to primitive form
	  rhs_[offset(i,j, 1)] /= prim[offset(i,j, 0)];
	  rhs_[offset(i,j, 2)] /= prim[offset(i,j, 0)];
	  rhs_[offset(i,j, 3)] /= (prim[offset(i,j, 0)]*cp_);
	  for(SizeType k = 0; k < DataType::nspec; ++k)
	    rhs_[offset(i,j, 4+k)] /= prim[offset(i,j, 0)];

	}
      }
      
      //! ------------- END INTERIOR ----------------------------------------
      
      

      return rhs_;
    }
    
    //dummy function
    void get_y_prev(typename VType::iterator y_prevIt, const T& h){}
  

  private:
    VType rhs_, C_, chem_, csource1_, csource2_, v1bdy_;
    SizeType Nx_, Ny_, npt_;
      DType a_,c_,hx_, hy_, Tignition_, Tmin_, halfdiam_;
    T cp_, cp1_,cp2_, p0_, density_,
      vx_, vy_, paravelo_, Tin_, qheat_, qh_,
      mb_, mmw_, rho_, qheatReac_;
    
    VType X_, X1_, X2_, X3_,X4_,X5_,
      Y_, Y1_;
    
    BuilderType BCS_;  //calculate chemistry
   
    MixtureAveragedViscosity<DataType,false> Mav_;
    MixtureAveragedConductivity<DataType,false> Mac_;
    MixtureAveragedDiffusion<DataType,false> Mad_, Mad2_, Mad3_, Mad4_, Mad5_;

    T dotOmega_[DataType::nspec];
    VType Yfrac_, Yfrac2_, Yfrac3_, Yfrac4_, Yfrac5_;

    BoundaryTemperatureDistribution<DType> bdytemp_;

    SizeType offset(SizeType i, SizeType j){
      return (i + Nx_*j);
    }
    
    SizeType offset(SizeType i, SizeType j, SizeType spec){
      return (offset(i,j) + spec*npt_);
    }
 
    template<class YFRAC, class TEMP>
    VType& chemical_source_term(const YFRAC& y, const TEMP& temperature){
      //reset
      qheatReac_ = T(); 
      cp_ = T();
      mmw_ = T();

      UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(mmw_,&y[0],DataType::molar_masses());
      
      rho_ = p0_*mmw_/(PhysicalConstants<T>::IdealGasConstant*temperature);
      
      //!mass conservation
      Yfrac_[DataType::nspec-1] = 1. - UnrollLoop<0,DataType::nspec-1>::template sum<1>(y);
      
      for(size_t k = 0; k < DataType::nspec; ++k){
	cp_ += (BCS_.C_p(k,temperature)*y[k]/DataType::molar_masses()[k]);
	//! compute concentrations
	C_[k] = rho_*y[k]/DataType::molar_masses()[k];  
      }
      
      for(size_t k = 0; k < DataType::nspec; ++k){		
	//! compute \f$ \dot{\omega}_k \ \forall k\f$ 
	chem_[k] = BCS_.net_formation_rate(k,temperature,C_);
	qheatReac_ += BCS_.H_T(k,temperature)*chem_[k]; 
	chem_[k] *= (DataType::molar_masses()[k]/rho_); //back from [X_k] to Y_k
      } 
      
      //!Temperature source term: MIND THE MINUS HERE!!!!!!!
      chem_[DataType::nspec] = -qheatReac_/(rho_*cp_);

      return chem_;
    }
   

  }; //end class 

}

#endif
