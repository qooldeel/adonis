#ifndef CORRECTION_AKA_2ND_APPROX_H2C6_1DIM_WITH_TRANSPORT_HH
#define CORRECTION_AKA_2ND_APPROX_H2C6_1DIM_WITH_TRANSPORT_HH

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

#include "../../sourceterms.hh"   //contains the original rhs
#include "../../../derivatives/jacobianofsource.hh"
#include "../../../modred/rvp.hh"
#include "../../../modred/indexmanipulation.hh" 
#include "../../../dunexternal/dunextensions.hh"
#include "../../../misc/useful.hh"
#include "../../../common/smartassign.hh"

#include "../../../misc/operations4randomaccesscontainers.hh"


#include "../../../dunexternal/scalingandregularisation.hh"
#include "../../../noderivativemethods/finitedifferences.hh"


namespace Adonis{

  template<class T>
  class H2MechWithDiffusion1dCORRECTION{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;
    
    typedef typename TypeAdapter<T>::Type DType;
    
    typedef ComputeManifold<DType,6,2,true> ManifoldType;
    typedef typename ManifoldType::CompositionType CompositionType;
    typedef typename ManifoldType::VecD VecD; //std::vector<double>
    typedef typename ManifoldType::VecS VecS; //std::vector<std::string> 
    typedef typename ManifoldType::CoordinateType CoordinateType;

    //types needed for 2nd approximation
    
    typedef Dune::FieldMatrix<DType,2,6> BTType;
    typedef Dune::FieldMatrix<DType,6,4> UType;
    typedef Dune::FieldMatrix<DType,4,6> NTType;
    typedef Dune::FieldMatrix<DType,4,4> PseudoInverseType;
    typedef Dune::FieldMatrix<DType,2,4> BTJUType; 
    typedef Dune::FieldVector<DType,4> DeltaUType;
    typedef Dune::FieldVector<DType,6> FullType;
    typedef Dune::FieldVector<DType,2> ReducedType;

    typedef Dune::FieldMatrix<DType,6,6> JacobianMatrixType;
        
    typedef H2Combustion6Spex<DType> OriginalFunctionType;

    //typedef JacS<DType,H2Combustion6Spex,ExprTmpl::MyVec> JacobianType;
    typedef FiniteDifference<OriginalFunctionType> NonDerivativeType;
    
    typedef JacS<DType,H2Combustion6Spex,ExprTmpl::MyVec> ADJacType;

  private:
    T D1_, D2_, D3_, D4_, D5_, D6_,  //the diffusion coefficients
      rpvh2o_, rpvh2_;

    ManifoldType Mani_;
    std::ofstream of_;
    //JacobianType Jf_;

    OriginalFunctionType fun_;  //original function 

    BTType B_T_;
    UType U_;
    
    VType non1_, //the unrepresented species
      non2_,
      non3_,
      non4_;
    

    //! approximation of \f$ u_{xx} \approx \frac{u_{i+1} -2u_i + u_{i-1}}{(\Delta x)^2}\f$
    template<class U>
    T fd1D2nd(const U& u, int i, const T& h){ //const diffusion
      return ( (u[i+1] - 2.*u[i] + u[i-1])/(h*h) );
    }
    
   
    
  public:
    
    enum{ 
      space = Constant<double>::spacePoints,
      numOfSpecies = 2,                      //only 2 now !!
      domainDim = numOfSpecies*space,
      fulldim = 6
    };

    
    H2MechWithDiffusion1dCORRECTION(unsigned dim = 0, const DType& rpvh2 = 0.3, const DType& rpvh2o = 0.6):rhs_(dim), rpvh2_(rpvh2), rpvh2o_(rpvh2o), fun_(6){
      
      Parameter<double,6> Dc("dat/Diffcoefs.dat");
      D1_ = Dc[0]; D2_ = Dc[1];  D3_ = Dc[2]; D4_ = Dc[3]; D5_ = Dc[4];
      D6_ = Dc[5];
     
      Jf_.set(fulldim,fulldim);

      VType pvvalues(numOfSpecies);
      pvvalues[0] = rpvh2;
      pvvalues[1] = rpvh2o;
      
      VecS names(numOfSpecies);
      names[0] = "H2";  
      names[1] = "H2O";

      //std::cout << " INVOKE JOCHEN's routine \"first(…) \"" << std::endl;
      
      CoordinateType rpvindex;
      rpvindex[0] = 0; rpvindex[1] = 4;


      Mani_.activate(names,pvvalues,rpvindex); //first(...) invoked
      
      B_T_ = Mani_.B_T();
      U_ = Mani_.U();

      //the unrepresented species
      non1_.resize(space);
      non2_.resize(space);
      non3_.resize(space);
      non4_.resize(space);

      //the initial concentration are taken to be the same as in the original system:
      DType onmani[] = {0.4549999999992859, 0.7780585738931507, 0.2366143850825262, 0.3628298037265891, 0.1479999999999196, 0.01594142610843904};
      
      for(unsigned i = 0; i < space; ++i){
	non1_[i] = onmani[1];
	non2_[i] = onmani[2];
	non3_[i] = onmani[3];
	non4_[i] = onmani[5];
      }
    }


    ~H2MechWithDiffusion1dCORRECTION(){
      //of_.close();
    }

    unsigned dim() const {return rhs_.size();}
    unsigned domain_dim() const{return domainDim;}

    /////////////////////////////////////////////////////
    //operator()(·) -- here, warmstart is invoked
    /////////////////////////////////////////////////////
    template<class X>
    VType& operator()(const X& z){
     
      std::vector<T> rV(12);

      T pk[12]; //reaction rates 
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
      
      const DType h = 1./(space-1);

      //bdy constant? Sim characterized via rpv (no time- and space-dpcy)
      //no flux across bdy
      for(int i = 0; i < numOfSpecies; ++i){
	rhs_[i*space] = 0.;                           //left
	rhs_[(i+1)*(space-1)] = 0.;                   //right
      }
     
      //BC's for unrepresented species:
      non1_[0] = 0.; non1_[space-1] = 0.;  //left and right respectively
      non2_[0] = 0.; non2_[space-1] = 0.;
      non3_[0] = 0.; non3_[space-1] = 0.;
      non4_[0] = 0.; non4_[space-1] = 0.;

      std::vector<T> reduced(numOfSpecies);
      CompositionType zM; //thats the full composition having 6 species

      ExprTmpl::MyVec<DType> jac;
      ExprTmpl::MyVec<DType> eval(6);

      UType JU;
      PseudoInverseType NTJU;  //rhs of lin. system for solving for delta u
      DeltaUType rhs, deltau;
      BTJUType BTJU;

      FullType Diffusion, 
	SPlusD;            //source term plus diffusion
      ReducedType Correction;
      
      //ExprTmpl::MyVec<T> Corr(2);

      JacobianMatrixType Jacobi;

      //++++++++++++++++++++++ Jacobian +++++++++++++++++++++++++++++++++++++++
      //NonDerivativeType FiDi(fun_);  //the Jacobian J via central finite diff.
      

      
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      //INTERIOR
      for(int i = 1; i < space-1; ++i){
	
	//fill reduced again -- maybe before assigning stuff to rhs_
	//WORST CASE: evaluate at every point
	smart_assign(reduced[0],z[i]);        
	smart_assign(reduced[1],z[space+i]);    
      

	Mani_.evaluate(reduced);   //evaluate new point here

	zM = Mani_.get_z();  //very first entry is from intialisation
	
	//std::cout << "z_sim = "<< zM << std::endl;
	
	ManipulateDynamicArray<DType,int,0,6>::fill(&eval[0],zM.begin(),zM.end());
	std::cout << "z_sim = "<< eval << std::endl;

	Jf_.jacobian(eval); //Evaluate Jacobian of original source term
	transform_2_matrix(Jf_.get_jacobian(),Jacobi,'r'); //row-wise
	
	//transform_2_matrix(FiDi.central(eval),Jacobi,'r');

	std::cout << "J(zM) = " << std::endl << Jacobi << std::endl;
	
	//scale Jacobian
	//column_row_equilibration(Jacobi);
	//std::cout << "row-column scaled Jacobian "<< std::endl<< Jacobi << std::endl; 

	JU =  Jacobi*U_; //matrix_matrix_multiplication<RowMajor,6>(Jf_.get_jacobian(),U_);
	
	NTJU = Mani_.get_N_T() * JU;
       
	//assign zM to current space point
	non1_[i] = zM[1];
	non2_[i] = zM[2];
	non3_[i] = zM[3];
	non4_[i] = zM[5];

	smart_assign( Diffusion[0], D1_ * fd1D2nd(z,i,h) );  //rpv1
	
	smart_assign( Diffusion[1], D2_ * fd1D2nd(non1_,i,h) );
	smart_assign( Diffusion[2], D3_ * fd1D2nd(non2_,i,h) );
	smart_assign( Diffusion[3], D4_ * fd1D2nd(non3_,i,h) );

	smart_assign( Diffusion[4], D5_ * fd1D2nd(z,space+i,h) );//rpv2
	smart_assign( Diffusion[5], D6_ * fd1D2nd(non4_,i,h) );


	//std::cout << "D{zM} = " << Diffusion << std::endl;
	//std::cout << "S(zM) = " << fun_(zM) << std::endl;

	//SPlusD = S(zM) + D(zM) 
	SPlusD = vector_vector_addition(fun_(zM),Diffusion);
	
	//rhs = N^T · (S(zM) + D(zM) )
	matrix_vector_product(rhs,Mani_.get_N_T(),SPlusD);
	
	//-rhs (field vector multiplication with scalar)
	rhs *= -1.;

	deltau = solve_rank_deficient_ls(NTJU,rhs);
	std::cout << "solution delta_ u = "<< deltau << std::endl;
	
	//A-POSTERIORI SCALING (normalization)
	unit_vector(deltau); //normalized u
	std::cout << "normalized delta_ u = "<< deltau << std::endl;

	BTJU = B_T_ * JU;
	std::cout << "B^TJU = "<< std::endl << BTJU << std::endl;
	row_column_equilibration(BTJU);
	//column_row_equilibration(BTJU);
	std::cout << "row-column scaled B^TJU = "<< std::endl<< BTJU << std::endl; 
	


	//B^TJU·\delta u, the correction term in eq. (37)
	matrix_vector_product(Correction,BTJU,deltau);
	std::cout << " CORRECTION for diffusive drawn-off effects = "<< Correction << std::endl;
	
	//unit_vector(Correction);
	//std::cout << "normalized CORRECTION = "<< Correction << std::endl;
       
	

	//note that here are species in dependence of others, only z_1 and z_5 occur (represented by corresponding forward and backward reactions)
	rV[0] =   pk[0] * z[i];             //zM[0];          
	rV[1]  =  pk[6] * zM[1]*zM[1];

	rV[4]  =  pk[2] * z[space+i];       //zM[4];             
	rV[5]  =  pk[8] * zM[1]*zM[5];
	rV[6]  =  pk[3] * z[i]*zM[3];       // zM[0]
	rV[7]  =  pk[9] * zM[1]*zM[5];

	rV[10] =  pk[5] * z[i]*zM[3];       // zM[0] 
	rV[11] =  pk[11] * z[space+i];      // zM[4]; 


	T diff1 =  D1_*fd1D2nd(z,i,h), diff2 = D5_*fd1D2nd(z,space+i,h);

	//with diffusion -- the reaction progress variables
	rhs_[i] = - rV[0]  + rV[1]              //H2
	  - rV[6]  + rV[7]
	  - rV[10] + rV[11]   + diff1    //up to here: 1st  approximation
	  //+ Diffusion[0] +       
	  + Correction[0];                      //B^TJU·delta_u
	
	//std::cout <<"  D1·ð²_x H_2 = " << diff1 << std::endl;
	//std::cout << "Diffusion[0] = "<<Diffusion[0]<<std::endl;
	
	rhs_[space+i] = - rV[4]  + rV[5]               //H2O
	  + rV[10] - rV[11]     + diff2  //up to here: 1st  approximation
	  //+ Diffusion[4] +    
	  + Correction[1];                      //B^TJU·delta_u
      
	//std::cout <<"  D5·ð²_x( H_2O) = " << diff2 << std::endl;
	 //std::cout << "Diffusion[4] = "<<Diffusion[4]<<std::endl;
	 
      } //end spatial discretisation

      return rhs_;
    }
    

     template<class TIME, class FV>
     VType& operator()(const TIME& time, const FV& z){   
       rhs_ = (*this).operator()(z);
       return rhs_;
     }

  private:
    VType rhs_;
     ADJacType Jf_;
  };

}

#endif 
