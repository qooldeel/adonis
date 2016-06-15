#ifndef VERY_SIMPLE_LOCAL_SPECIES_RECONSTRUCTION_HH
#define VERY_SIMPLE_LOCAL_SPECIES_RECONSTRUCTION_HH

#include <iostream>

#include "../common/globalfunctions.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../dunestuff/fmatrix.hh"
#include "../dunexternal/dunextensions.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../containers/staticarray.hh"
#include "rvp.hh"
#include "indexmanipulation.hh"
#include "bounds.hh"

namespace Adonis{

  template<class T, int N, int R, int ENCOD>
  class VerySimpleSpeciesReconstruction{
  public:
    typedef ExprTmpl::MyVec<T> VType;
    typedef T value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef std::size_t SizeType;
    
    typedef Dune::FieldMatrix<value_type,R,N> BTType;
    typedef Dune::FieldMatrix<value_type,N,N-R> UType;
    typedef Dune::FieldMatrix<value_type,N,N> DiagType; 

    typedef StaticArray<value_type,N> FullType;
    
    typedef Norm<'2',BaseType> NormType;

    typedef ConcatenateMatrix<R,N,VType> MatrixConcatenatorType;
    typedef ConcatenateMatrix<N,N,VType> SecondMatrixConcatenatorType;
    
    typedef MechanismBounds<BaseType,ENCOD> MechBdsType;

    VerySimpleSpeciesReconstruction():nelem_(0),maxIt_(0),count_(0),tol_(0){}

    //! when Ct and b are constant
    template<class V1, class IXIT>
    void initialize(const V1& v, int nelem, const VType& Ct, const VType& b, const IXIT& ixptr,
		    SizeType maxIt = 7, const BaseType& tol = 1.e-12){
      adonis_assert((int)b.size() == nelem);
      zM_.resize(N);
      for(int i = 0; i < N; ++i)
      	smart_assign(zM_[i],v[i]);
      
      int rhsdim = nelem+R+N;
      btilde_.resize(rhsdim);
      
      nelem_ = nelem;
      b_ = b;

      IndexManipulator<N,R,ExprTmpl::MyVec> IM;
      IM.create_unrepresented_index(ixptr);
      ReactionProgressVariables<value_type,N,R> RPV(BT_,U_); //references 
      RPV.create_B_T(ixptr);     //now B^T is filled 
      RPV.create_U(IM.get_index()); //now U is filled

      // std::cout << "B^T = "<< BT_ << std::endl;
      // std::cout << "U = " << U_ << std::endl;

      MatrixConcatenatorType Join;
      Ctilde_ = Join.columnwise(nelem,N,Ct,BT_);

      // std::cout << "Ctilde_ = "<< Ctilde_ << std::endl;
      // print_in_matrix_style(Ctilde_,N);

      maxIt_ = maxIt;
      count_ = 0;
      tol_ = tol;

      perturbation_.resize(N);
     
      residual_.resize(rhsdim);
      identity(Diag_); //identity matrix 
       //! enlarge Ctilde_
      Ctilde_ = concat_.columnwise(nelem_+R,N,Ctilde_,Diag_);

      g_.resize(rhsdim);
      mechbds_.initialize(N);
    }

    VType& get_z() {return zM_;}
    const VType& get_z() const{return zM_;}


    template<class V1, class W>
    VType& evaluate(const V1& red, const W& u){
      //! assigns also from CppAD
      UnrollLoop<0,N>::assign(zM_.begin(),u.begin()); 

      //! overwrie btilde with b, followed by red
      std::cout << "btilde_ (before) = "<< btilde_ << std::endl; 
      concatenate(btilde_,b_,red);
      
      concatenate(btilde_,Max(mechbds_.low(),Min(zM_,mechbds_.up())));
      std::cout << "btilde_ (after) = "<< btilde_ << std::endl;
      
     
      print_in_matrix_style(Ctilde_,N,"Ctilde_ (enlarged)",12,true);

      residual_ = matrix_vector_product(Ctilde_,zM_) - btilde_;
      
      BaseType nm = NormType::norm(residual_);
      std::cout << " redisual = "<< residual_ << std::endl << " ||residual_|| = "<< nm << std::endl;
      if(nm > tol_){
	std::cout << "PERTURBATION IS CALCULATED" << std::endl;
      	perturbation_ =  solve_linear_constraints(zM_);
      }
      else{ //o.k. everythings o.k. then perturbation_ = u
	UnrollLoop<0,N>::assign(perturbation_.begin(),u.begin());
      }

      return perturbation_;      
    }

    private:
    int nelem_;
    SizeType maxIt_, count_;
    BaseType tol_;
   
    VType Ctilde_, btilde_, b_, perturbation_, zM_;
    BTType BT_;
    UType U_;
    DiagType Diag_;
   
    
    VType g_, uk_, delta_, residual_; 
    MechBdsType mechbds_;

    SecondMatrixConcatenatorType concat_;

    //!actually this is exactly what is performed in massbalanceinvariant.hh
    //! only with Ctilde_ and btilde_
    template<class VEC>
    VType& solve_linear_constraints(const VEC& u){ //argument: full state
      int m = static_cast<int>(btilde_.size()),
	nrhs = 1;
    
      //g(u) := Ctilde_u -btilde_
      //g'(u) = Ctilde_
      //solve g'(u) du = -g(u) iteratively 
      uk_ = u;
      for(SizeType k = 1; k <= maxIt_; ++k){
	count_ = k;
	std::cout << k << ".)   ";
	//!negative residuals
	g_ = -matrix_vector_product(Ctilde_,uk_) + btilde_;
	delta_ = solve_rank_deficient_ls(Ctilde_,m,g_,nrhs); //Ctilde_ is copied
	std::cout << "delta = "<< delta_ << std::endl;
	
	uk_ += delta_;
	
	if(NormType::norm(delta_) <= tol_){
	  break;  //leave loop since convergence has been achieved
	}
      }
      
      if(count_ == maxIt_){
	//ADONIS_INFO(Information, "Max. number of iterations reached (maxit = "<< maxit << ").\n   No convergence with tol = "<< tol << ".");
	ADONIS_ERROR(IterationError, "Max. number of iterations reached (#max. iers = "<< maxIt_ << ").\n   No convergence detected with tol = "<< tol_ << ".");
	
      }
      
      return uk_;   //this is a perturbation of u that satisfies linear
      //mass balances
    }
    
  };

}//end namespace

#endif
