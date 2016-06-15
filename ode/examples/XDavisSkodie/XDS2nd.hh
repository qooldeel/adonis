#ifndef XDS_2_APPROX_AKA_CPA_NUMERICAL_TREATMENT_HH
#define XDS_2_APPROX_AKA_CPA_NUMERICAL_TREATMENT_HH

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>


#include "../../../expressiontemplates/exprvec.hh"
#include "../../constants.hh"
#include "../../../io/readinparameters.hh"
#include "../../../modred/manifold.hh"

#include "../../../common/typeadapter.hh"

//==========
#include "../extendedDavisSkodje.hh"   //contains the original rhs
#include "../weirdtransport.hh"          //D2(x) isn't constant 
//==========

#include "../../../derivatives/jacobianofsource.hh"
#include "../../../modred/rvp.hh"
#include "../../../modred/indexmanipulation.hh" 
#include "../../../dunexternal/dunextensions.hh"
#include "../../../misc/useful.hh"
#include "../../../common/smartassign.hh"

#include "../../../misc/operations4randomaccesscontainers.hh"
#include "../../../misc/misctmps.hh"

#include "../../../dunexternal/scalingandregularisation.hh"
#include "../../../noderivativemethods/finitedifferences.hh"


namespace Adonis{

  template<class T>
  class XDS2ndApproxAkaCPA{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;

    typedef Constant<double> ConstantType;

    typedef typename TypeAdapter<T>::Type DType;

    enum{fudim = 2, reddim = 1, undim = fudim - reddim};

    typedef ComputeManifold<DType,fudim,reddim,true> ManifoldType;
    typedef typename ManifoldType::CompositionType CompositionType;
    typedef typename ManifoldType::VecD VecD; //std::vector<double>
    typedef typename ManifoldType::VecS VecS; //std::vector<std::string> 
    typedef typename ManifoldType::CoordinateType CoordinateType;
  

    //types needed for 2nd approximation
    typedef Dune::FieldMatrix<DType,reddim,fudim> BTType;
    typedef Dune::FieldMatrix<DType,fudim,undim> UType;
    typedef Dune::FieldMatrix<DType,undim,fudim> NTType;
    typedef Dune::FieldMatrix<DType,undim,undim> PseudoInverseType;
    typedef Dune::FieldMatrix<DType,reddim,undim> BTJUType; 
    typedef Dune::FieldVector<DType,undim> DeltaUType;
    typedef Dune::FieldVector<DType,fudim> FullType;
    typedef Dune::FieldVector<DType,reddim> ReducedType;

    typedef Dune::FieldMatrix<DType,fudim,fudim> JacobianMatrixType;

    //originial right hand side
    typedef RPDavisSkodje<DType> OriginalFunctionType;

    typedef FiniteDifference<OriginalFunctionType> NonDerivativeType;
    
    typedef JacS<DType,RPDavisSkodje,ExprTmpl::MyVec> ADJacType;
  
  private:
    DType D1_,   //the diffusion coefficient, note that D2 is most likely nonlinear!
      a_,
      b_,
      c_,
      d_,
      e_,
      eps_,
      rpv_,  //rpv of the reduced problem
      h_;   //spatial stepsize

    ManifoldType Mani_;
    std::ofstream of_;
    CompositionType zM_;
    VecD reduced_;

    OriginalFunctionType fun_;  //original function 

    BTType B_T_;
    UType U_;

    VType non_;   //the fudim-reddim unrepresented species
  
    //=============
    VType rhs_;
    ADJacType Jf_;
    //=============

    
    // TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · · TEST
    // TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · · TEST 
    // TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · · TEST
    //!this is just for checking
    template<class Scal>
    NTType check_exact_N_T(const Scal& z1) const {
      NTType Nt;
      
      Scal fprime = 1./(ntimes<2>(1+a_*z1)),
	length = 1./(std::sqrt(1 + ntimes<2>(fprime)));

      Nt[0][0] = -length*fprime; 
      Nt[0][1] = length;

      return Nt;
    }

    //!this is just for checking
    template<class Scal>
    FullType check_exact_zM(const Scal& z1) const {
      FullType zMani;
      zMani[0] = z1;
      zMani[1] = z1/(1+a_*z1);
      
      return zMani;
    }


    //! just for checking -- exact CPA
    template<class Z, class D>
    typename Z::value_type analytical_cpa_value(const Z& z, int i, const value_type& x, D& D2) const {
      typedef typename Z::value_type val_t;
      
      val_t fprime =  1./(ntimes<2>(1 + a_*z[i])),
	diffuse1 = D1_*(z[i+1]- 2.*z[i] + z[i-1])/(ntimes<2>(h_));

      
      val_t cpa = -d_*z[i] +  diffuse1 +     //1st approx
	c_/(c_*fprime + 1)*( d_*z[i]*fprime - z[i]/(ntimes<2>(1+b_*z[i])) - fprime*diffuse1 + ( D2(x+0.5*h_)*(z[i+1]/(1.+a_*z[i+1]) - z[i]/(1.+a_*z[i])) - D2(x-0.5*h_)*(z[i]/(1.+a_*z[i]) - z[i-1]/(1.+a_*z[i-1])) )/(ntimes<2>(h_)) ); //correction

      return cpa;
    }

    template<class Z, class D>
    typename Z::value_type analytical_correction(const Z& z, int i, const value_type& x, D& D2) const{
      typedef typename Z::value_type val_t;
      
      val_t fprime =  1./(ntimes<2>(1 + a_*z[i])),
	diffuse1 = D1_*(z[i+1]- 2.*z[i] + z[i-1])/(ntimes<2>(h_));

      return ( c_/(c_*fprime + 1)*( d_*z[i]*fprime - z[i]/(ntimes<2>(1+b_*z[i])) - fprime*diffuse1 + ( D2(x+0.5*h_)*(z[i+1]/(1.+a_*z[i+1]) - z[i]/(1.+a_*z[i])) - D2(x-0.5*h_)*(z[i]/(1.+a_*z[i]) - z[i-1]/(1.+a_*z[i-1])) )/(ntimes<2>(h_)) ) );
    }
    // TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · · TEST
    // TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · · TEST 
    // TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · · TEST


     //! approximation of \f$ u_{xx} \approx \frac{u_{i+1} -2u_i + u_{i-1}}{(\Delta x)^2}\f$
    template<class U>
    T central_fd_2nd(const U& u, int i, const T& h) const{ //const diffusion D1
      return ( (u[i+1] - 2.*u[i] + u[i-1])/(ntimes<2>(h)) );
    }


    template<class NLD, class U, class X, class H>
    typename NLD::value_type strange_nonlin_diff(NLD& D2, const U& u, int i, const X& x, const H& h) const{
      return ( (D2(x+0.5*h)*(u[i+1] - u[i]) - D2(x-0.5*h)*(u[i] - u[i-1]))/(ntimes<2>(h)) );
    }

   
  public:
    enum{
      numOfSpecies = 1,          //REDUCED
      space = ConstantType::spPts,            //number of space points
      domainDim = numOfSpecies*space
    };

    XDS2ndApproxAkaCPA(std::size_t dim = 0, const DType& rpv = 1.0):rhs_(dim),  reduced_(numOfSpecies), rpv_(rpv),fun_(fudim),h_(1./(space-1)){
      Parameter<DType,7> P("Parameter.dat");
      a_ = P[0]; b_ = P[1]; c_ = P[2]; d_ = P[3]; e_ = P[4]; 
      D1_ = P[5]; eps_ = P[6];

      Jf_.set(fudim,fudim); //set up Jacobian

      std::vector<std::string> names(numOfSpecies);
      names[0] = "z_1"; // see RenDS/START/grid.dat as clue

      std::vector<DType> values(numOfSpecies);
      values[0] = rpv_; 

      CoordinateType rpvindex;
      rpvindex[0] = 0;


      Mani_.activate(names,values,rpvindex);

      B_T_ = Mani_.B_T();
      U_ = Mani_.U();

      non_.resize(space); //unrepresented species
      
      //the initial concentration for the unrepresented species are taken to be the same as in the original system -- in that case we take the initial concenctration as given in the paper, p. 251, bottom 1st text column
      for(unsigned i = 0; i < space; ++i){
	T step = i*h_;
	non_[i] = step/(1 + a_*step);   //in this case non_ corresponds to z2
      }
    }

    ~XDS2ndApproxAkaCPA(){}
    
    unsigned dim() const {return rhs_.size();}
    unsigned domain_dim() const{return domainDim;}


    
    //================== PARENTHESIS-OPERATOR ==================================

    /////////////////////////////////////////////////////
    //operator()(·) -- here, warmstart is invoked
    /////////////////////////////////////////////////////
    template<class X>
    inline VType& operator()(const X& y){
      std::cout <<"##### INVOKE PARENTHESIS OPERATOR" << std::endl;

      //!NOTE: all boundary values are constant, hence \f$\frac{d}{dt} y_{\operatorname{boundary}} = 0.\f$
      //left bdy                                  //right bdy
      rhs_[0] = 0.;            rhs_[space-1] = 0.; 
      

      //BC's for unrepresented species, note that the BCs are constant, i.e. 
      // z2(t,x=0) = 0. and z2(t,x=1) = 1./(1+a_). Hence, \f$\frac{d}{dt} y_{\operatorname{boundary}} = 0.\f$ 
      non_[0] = 0.;            non_[space-1] = 0.;

      
      ExprTmpl::MyVec<DType> jac;
      ExprTmpl::MyVec<DType> eval(fudim);

      
      UType JU;
      PseudoInverseType NTJU;  //rhs of lin. system for solving for delta u
      DeltaUType rhs,  //!DON'T CONFUSE rhs with rhs_ !!
	deltau;
      BTJUType BTJU;

      FullType Diffusion, 
	SPlusD;            //source term plus diffusion
      ReducedType Correction;

      JacobianMatrixType Jacobi;

      //++++++++++++++++++++++ Jacobian +++++++++++++++++++++++++++++++++++++++
      //NonDerivativeType FiDi(fun_);  //the Jacobian J via central finite diff.
      

      
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      //======= the strange DIFFUSION for D2 ============================
      DiffusiveTransport<T> D2(D1_,e_); //fields consist of const refs!
     
      //interior
      T z1_z2 = T();
      T x = T();


      //INTERIOR
      for(int i = 1; i < space-1; ++i){

	//std::cout << "B^T = "<< B_T_ << std::endl<< "U = " << U_ <<std::endl; 

	//============ retrieve point ==========================
	smart_assign(reduced_[0],y[i]);
      
	std::cout << std::endl <<">>>> RETRIEVE POINT: \n reduced_ = "; print_all(reduced_,11);
	Mani_.evaluate(reduced_);  //compute new manifold point
	
	zM_ = Mani_.get_z();         // zM = [z1, z2]^T
	//======================================================

	//copy it to eval such that you circumvent .size() which is needed
	//in the .jacobian(·)-call but which is absent in Dune::FieldVector
	ManipulateDynamicArray<DType,int,0,fudim>::fill(&eval[0],zM_.begin(),zM_.end());
	
	Jf_.jacobian(eval); //Evaluate Jacobian of original source term
	
	
	transform_2_matrix(Jf_.get_jacobian(),Jacobi,'r'); //row-wise
	
	
	//transform_2_matrix(FiDi.central(eval),Jacobi,'r');

	std::cout << "J(zM) = " << std::endl << Jacobi << std::endl;

       
	JU =  Jacobi*U_;
	NTJU = Mani_.get_N_T() * JU;
	std::cout << "N^TJU = "<< NTJU << std::endl;

	non_[i] = zM_[1];   //note z2 is the unrepresented species

	//! Diffusion of the full chemistry, i.e. FullType -- diff for z1
	smart_assign(Diffusion[0],  D1_ * central_fd_2nd(y,i,h_));
	
	x = i*h_; //!actual space point
	
	

	//! note: has to be evaluated at the <I> unrepresented </I> species z2
	smart_assign(Diffusion[1], strange_nonlin_diff(D2,non_,i,x,h_));
	
	std::cout << "Diffusion = " << Diffusion << std::endl;
	//unit_vector(Diffusion); std::cout << "norm. Diffusion = " << Diffusion << std::endl;
	
	std::cout << "f(zM) = " << fun_(zM_) << std::endl;

	SPlusD = vector_vector_addition(fun_(zM_),Diffusion);
	matrix_vector_product(rhs,Mani_.get_N_T(),SPlusD);
	
	//-rhs (field vector multiplication with scalar)
	rhs *= -1.;
	
	std::cout << "rhs = "<< rhs << std::endl;

	//!NOTE: both rhs, NTJU  and deltau are scalars
	//deltau[0] = rhs[0]/NTJU[0][0];
	deltau = solve_rank_deficient_ls(NTJU,rhs); //! this is more automatic
	std::cout << "solution delta_ u = "<< deltau << std::endl;

	//==================================================================
	//========= A-POSTERIORI SCALING (normalization) ===================
	//==================================================================

	//! ATTENTION: don't norm a scalar -- this isn't sooo good ;-)
	//unit_vector(deltau); //normalized u
	//std::cout << "normalized delta_ u = "<< deltau << std::endl;

	BTJU = B_T_ * JU;  //!is always the same 
	std::cout << "B^TJU = "<< std::endl << BTJU << std::endl;

	//! ATTENTION: \f$ B^TJU \f$ is a 1 x 1 matrix (a scalar) !! 
	row_column_equilibration(BTJU); //column_row_equilibration(BTJU)
	std::cout << "row-column scaled B^TJU = "<< std::endl<< BTJU << std::endl; 
	//==================== END A-POST. SCALING =========================

	//! B^TJU·\delta u, the correction term in eq. (37) -- all are scalars
	//Correction[0] = BTJU[0]*deltau[0];
	matrix_vector_product(Correction,BTJU,deltau);

	

	//unit_vector(Correction);
	//std::cout << "normalized CORRECTION = "<< Correction << std::endl;


	//! NOW COMES THE ACTUAL SOURCE 
	std::cout << "==========================================="<<std::endl;

	//!TEST ONLY
	std::cout << "zM_ = "<< zM_ << std::endl;
	std::cout << "exact zM = "<<check_exact_zM(zM_[0]) << std::endl;
	std::cout << "Normal space basis N^T = " << std::endl <<
	  Mani_.get_N_T(); 
	std::cout << "EXACT normal sp. basis = "<< std::endl << 
	  check_exact_N_T(zM_[0]) << std::endl;


	z1_z2 = zM_[1] - y[i]/(1. + a_*y[i]);      //z2 - z1/(1+az1)

	

	rhs_[i] = c_/eps_*z1_z2 - d_*y[i] + D1_*central_fd_2nd(y,i,h_) //1st 
	          + Correction[0];   //2nd -- this is used for the CPA 
	
	
	
	//TEST ONLY
	//std::cout << "Numerical CPA val.  = " << rhs_[i] << std::endl;
	//std::cout << "Analytical CPA val. = "<< analytical_cpa_value(y,i,x,D2) << std::endl;

	std::cout << "Num. CORRECTION = "<< Correction << std::endl;
	std::cout << "Ana. Correction = "<< analytical_correction(y,i,x,D2) << std::endl;

	std::cout << "==========================================="<<std::endl;


      }//end spatial discretization, i.e. loop over domain \Omega
 
      return rhs_;
    }


    //! fake time dependent operator
    template<class V>
    VType& operator()(const T& time, const V& y){
      return (rhs_ = (*this).operator()(y));
    }
  }; //end class


} //end namespace 

#endif
