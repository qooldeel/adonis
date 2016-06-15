#ifndef MY_SPECIES_RECONSTRUCTION_XTENDED_HH
#define MY_SPECIES_RECONSTRUCTION_XTENDED_HH

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../common/elementaryoperations.hh"
#include "../common/typeadapter.hh"
#include "../common/smartassign.hh"

#include "../expressiontemplates/exprvec.hh"
#include "../templatemetaprograms/unrollloop.hh"
#include "../linalg/linearsystemsolvers.hh"
#include "../misc/operations4randomaccesscontainers.hh"

#include "../modred/indexmanipulation.hh"
#include "../modred/rvp.hh"
#include "../dunestuff/fmatrix.hh"
#include "../dunexternal/dunextensions.hh"


#include "fulldotomegawithtemperature.hh"
namespace Adonis{


  /**
   * \brief This is only an alternative attempt aiming at replacing the
   * NLP formulation
   *
   * \tparam T value type
   * \tparam N full dimension
   * \tparam R reduced dimension
   * \tparam ENCO integer encoding for mechanism needed e.g. for bounds on vars 
   */
  template<class T, int N, int R, int ENCOD, char NORM = '2'> 
  class ReconstructSpecies{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType; 
    typedef typename TypeAdapter<T>::Type DType;
     typedef ExprTmpl::MyVec<DType> VDType; 
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef Norm<NORM,typename TypeAdapter<T>::BaseType> NormType;
    typedef std::size_t SizeType;
    
    typedef ExprTmpl::MyVec<size_t> IndexVType; 
    typedef Dune::FieldMatrix<T,R,N> BTType;
    typedef Dune::FieldMatrix<T,N,N-R> UType;

    typedef DotOmegaWithTemperatureFull<T,ENCOD> DotOmegaTempType;
    typedef typename DotOmegaTempType::DataType DataType;
    typedef typename DotOmegaTempType::ArrayType FullType;
    typedef JacobianOfDotOmegaWithTemperatureFull<T,ENCOD> JacobianType; 
     
    enum{NP1=N+1};

    //========= maxnumber of newton iters ============
    enum{maxNewt = 5};
    //================================================

    
   
    ReconstructSpecies():fulldim_(0),count_(0),isConverged_(false),isPerturbed_(false),kn_(T()),kEstim_(T()),newtTol_(T()),indexptr_(0),whichTimeStepper_('e'),maxit_(0),nmtol_(0.),Tol_(0.){ adonis_assert(N>R);}

   
    size_t size() const {return zM_.size();}

    //! guessSpeciesTemp contains full species + temperature (last entry)
    template<class V>
    void initialize(const V& guessSpeciesTemp, const T& kn, const BaseType& newtTol, char wts = 'e', int maxit= 12, const BaseType& nmtol = 1.e-09, const BaseType& Tol = 1.e-09 ){
     
      adonis_assert((int)fun_.dim() == N+1); 
      adonis_assert((int) guessSpeciesTemp.size() == N+1);
      fulldim_ = N+1;  //species + temperature

      isConverged_ = false;
      isPerturbed_ = false;
      zM_.resize(fulldim_);
      
      G_.resize(fulldim_); 
      Gprime_.resize(fulldim_*fulldim_);   //Jacobian
     
      zMchem_.resize(N); //chemical part only
     
      
      for(int i = 0; i < fulldim_; ++i)
      	smart_assign(zM_[i],guessSpeciesTemp[i]); //full composition with zM initialized
     
          
      NablaChem_.initialize(zM_);
      
      kn_ = kn;
      kEstim_ = kn;
      newtTol_ = newtTol;
      whichTimeStepper_ = wts;
      maxit_ = maxit;
      nmtol_ = nmtol;
      Tol_ = Tol;

      //! just some indexing stuff
      rpvindex_.resize(R); 
      UnrollLoop<0,R>::assign(rpvindex_.begin(),DataType::rpv_index());
      IndexManipulator<N,R,ExprTmpl::MyVec> IM;
      IM.create_unrepresented_index(rpvindex_);
      ReactionProgressVariables<T,N,R> RPV(BT_,U_); //references 
      RPV.create_B_T(rpvindex_);     //now B^T is filled
      unrepindex_ = IM.get_index();
      RPV.create_U(unrepindex_); //now U is filled
 
    }
    
    const IndexVType& rpv_index() const {return rpvindex_;}
    const IndexVType& unrep_index() const {return unrepindex_;}

    const BTType& B_T() const {return BT_;}
    const UType& U() const {return U_;}

    const VType& get_z() const {return zM_;}
    VType& get_z() {return zM_;}

    const char which_time_stepper() const {return whichTimeStepper_;} 

    template<class FULLSTATE, class REDSTATE>
    VType& assign(const REDSTATE& redstate, const FULLSTATE& fullstate){
      UnrollLoop<0,NP1>::assign_rac2rac(zM_,fullstate);
      std::cout << "zM_ = "<< zM_ << std::endl;
      // std::cout << "redstate = "<< redstate << std::endl;
      UnrollLoop<0,R>::assign_smaller_2_greater_rac(zM_,redstate,DataType::rpv_index());
      //smart_assign(zM_[N], redstate[R]);  //temperature
     
      std::cout << "B^T·zM_ = r: "<< zM_ << std::endl;

      return zM_;
    }

    template<class CTILDE, class BTILDE>
    VType& enforce_linearized_constraints_wrt_chemistry(const CTILDE& Ctilde, const BTILDE& btilde){
      //!assign chemistry first
      UnrollLoop<0,N>::assign(zMchem_.begin(),zM_.begin()); 
      std::cout << "zMchem_ = "<< zMchem_ << std::endl;
      (residual_.size() != btilde.size()) ? residual_.resize(btilde.size()) : do_nothing();
      (mvp_.size() != btilde.size()) ? mvp_.resize(btilde.size()) : do_nothing();
      residual_ = matrix_vector_product(Ctilde,zMchem_) - btilde;
     
      std::cout << "Ctilde·zMchem-btilde = " << residual_ << std::endl;
      isPerturbed_ = false; //reset
      BaseType nm = NormType::norm(residual_);
      if(nm  > nmtol_){ //o.k. calculate perturbation that satisfies constraints
	std::cout << "||residual_|| = "<< nm << std::endl;
	int m = (int)btilde.size(); //number of rows of matrix Ctilde
	for(SizeType nu = 1; nu <= maxit_+1; ++nu){
	  count_ = nu;
	   mvp_ = -matrix_vector_product(Ctilde,zMchem_) + btilde;
	   defect_ = solve_rank_deficient_ls(Ctilde,m,mvp_,1); 
	   zMchem_ += defect_;
	  if(NormType::norm(defect_) <= Tol_)
	    break;  //leave loop since convergence has been achieved
	}
	if(count_ == maxit_)
	  ADONIS_ERROR(IterationError, "Max. number of iterations reached (maxit = "<< maxit_ << ").\n   No convergence with tol = "<< Tol_ << ".");
      
	isPerturbed_ = true;
      }
      
      if(isPerturbed_){ //only assign back in case of perturbation
	UnrollLoop<0,N>::assign(zM_.begin(),zMchem_.begin());
	std::cout << "zMchem_ = "<< zMchem_ << std::endl;
      }
      return zM_;
    }

    //both low and up are of size N+1, the last entry standing for temperature
    template<class LOW, class UP>
    VType& project_onto_bounds(const LOW& low, const UP& up){
      zM_ = Max(low,Min(zM_,up));
      std::cout << "zM (proj. on bounds) = "<< zM_ << std::endl;
      return zM_;
    }

    //! perform explicit or implicit Euler step with step size kn_
    //! \param X mass fractions, last entry temperature
    //! \param p pressure
    template<class X, class P>
    X& time_step(X& y, const P& press){
      adonis_assert(y.size() == N+1);
      
      fun_.set_pressure(press);

      if(whichTimeStepper_ == 'i' ||  whichTimeStepper_ == 'I'){
	std::cout << "IMPL. EULER:" << std::endl;
	y_prev_ = y;
	y_nu_ = y_prev_;
	count_ = 1; //reset
	BaseType nm = BaseType();
      	for(int l = 1; l <= maxNewt; ++l){ 
	  //y_nu_[N] = y[N]; //fixed temperature
	  count_ = l;
	  std::cout << " " << l << ".)  "<< y_nu_ << std::endl; 
	  G_ = -(y_nu_ - y_prev_ - kn_*fun_(y_nu_));
	  Gprime_ = -kn_*NablaChem_.jacobian(y_nu_);
      	  //print_in_matrix_style(Gprime_,fulldim_);
      	  update_diagonal<AddBasicElements>(Gprime_,fulldim_,1.);
      	  good_square_solve(Gprime_,fulldim_,G_,1);
      	  y_nu_ += G_;

	  nm = NormType::norm(G_);
	  std::cout << "||delta y|| = "<< nm << std::endl;
	  if(nm <= newtTol_){
      	    break;
      	  }
      	} //end NEWTON iter
	if(count_ == maxNewt)
	  ADONIS_INFO(Information,"Newton iteration did not converge withing "<< maxNewt << " steps and tolerance = "<< newtTol_ << ".");
	y = y_nu_;
	//y_prev_ = y;
      }
      else if (whichTimeStepper_ == 'e' || whichTimeStepper_ == 'E'){
	std::cout << "EX. EULER" << std::endl;
      	y += kn_*fun_(y);
      }
      else{
	ADONIS_ERROR(NotImplementedError,"Time stepping method whichTimeStepper_ = "<<whichTimeStepper_ << " undefined.");
      }
      std::cout << "concentrations + temperature: "<< y << std::endl;
      return y;
    }

    //! quantities that are already available when fun_(·) is called
    //!access full \$\dot{\omega}\$ from <TT>DotOmegaWithTemperatureFull</TT>
    const FullType& dot_omega() const {return fun_.dot_omega();}
   
    const FullType& concentrations() const {return fun_.concentrations();}
    const T& wbar() const {return fun_.wbar();}
    const T& rho() const {return fun_.rho();}
    const T& heat() const {return fun_.heat();}
    const T& cp() const {return fun_.cp();}

  private:
    int fulldim_, count_;
    mutable bool isConverged_, isPerturbed_;
    VType Gprime_, G_;
    //VType dotomega_;
    VType zM_, zMchem_, z_nu_,z_prev_,residual_,defect_, mvp_;
    VType y_prev_, y_nu_;
    T kn_,
      kEstim_,
      newtTol_;
    size_t* indexptr_;
    char whichTimeStepper_;
    BaseType maxit_,nmtol_,Tol_;
   
    JacobianType NablaChem_;
    BTType BT_;
    UType U_;
   
    DotOmegaTempType fun_; //everything for thermo-chemistry is initialized
    IndexVType rpvindex_,unrepindex_;
  };



}

#endif
