#ifndef XDS_1_APPROXIMATION_NUMERICAL_TREATMENT_HH
#define XDS_1_APPROXIMATION_NUMERICAL_TREATMENT_HH

#include <vector>
#include <string>

#include "../../../expressiontemplates/exprvec.hh"
#include "../../../misc/misctmps.hh"
#include "../../../common/globalfunctions.hh"

#include "../weirdtransport.hh"

#include "../../../io/readinparameters.hh"

#include "../../constants.hh"
#include "../../../modred/manifold.hh"

#include "../../../common/typeadapter.hh"

namespace Adonis{

  template<class T>
  class XDS1stApprox{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;

    typedef Constant<double> ConstantType;

    typedef typename TypeAdapter<T>::Type DType;

    typedef ComputeManifold<DType,2,1,true> ManifoldType;
    typedef typename ManifoldType::CompositionType CompositionType;
    typedef typename ManifoldType::VecD VecD; //std::vector<double>
    typedef typename ManifoldType::VecS VecS; //std::vector<std::string> 
    typedef typename ManifoldType::CoordinateType CoordinateType;

    private:
    T a_,
      c_,
      d_,
      D1_,
      eps_;
    
    VType rhs_;
    ManifoldType mani_;
    CompositionType zM_;
    VecD reduced_;

    template<class U>
    T central_fd_2nd(const U& u, int i, const T& h){
      return ( (u[i+1] - 2.*u[i] + u[i-1])/(ntimes<2>(h)) );
    }

  public:
    enum{
      numOfSpecies = 1,          //REDUCED
      space = ConstantType::spPts,            //number of space points
      domainDim = numOfSpecies*space
    };

    
    XDS1stApprox(std::size_t dim = 0):rhs_(dim), reduced_(numOfSpecies){
      Parameter<double,7> P("Parameter.dat");
      a_ = P[0]; c_ = P[2]; d_ = P[3]; D1_ = P[5]; eps_ = P[6];
     

      // invoke first(…)
      std::vector<std::string> names(numOfSpecies);
      names[0] = "z_1"; // see RenDS/START/grid.dat as clue

      std::vector<DType> values(numOfSpecies);
      values[0] = 1.; 

      CoordinateType rpvindex;
      rpvindex[0] = 0;


      mani_.activate(names,values,rpvindex);

    }
    
    size_t dim() const {return rhs_.size();}
    size_t domain_dim() const{return domainDim;}

     /**
     * \brief the discretised system in space 
     */
    template<class X>
    VType& operator()(const X& y){
      // std::cout << "rhs_.size() = "<<rhs_.size() << "  y.size() = "<<y.size() << std::endl; 
     
      const T h = 1./(space-1);  //interval [0,1]
      
      //!NOTE: all boundary values are constant, hence \f$\frac{d}{dt} y_{\operatorname{boundary}} = 0.\f$
      //left bdy                                  //right bdy
      rhs_[0] = 0.;             rhs_[space-1] = 0.; 
      
      //interior
      T z1_z2 = T();
      
       //==================== SPACE DISCRETISATION ============================
      
      for(int i = 1; i < space-1; ++i){
	
	//============ retrieve point ==========================
	smart_assign(reduced_[0],y[i]);
	std::cout << "reduced = "; print_all(reduced_,11);
	mani_.evaluate(reduced_);  //compute new manifold point
	
	zM_ = mani_.get_z();         // zM = [z1, z2]^T
	//======================================================

	std::cout << "==========================================="<<std::endl;
	std::cout << "zM_ = "<< zM_ << std::endl;
	
	//for a point on *the specific* Sim [z1,z1/(1+az1)], z1_z2 = 0
	
	z1_z2 = zM_[1] - y[i]/(1. + a_*y[i]);      //z2 - z1/(1+az1)

#ifndef NDEBUG
	DType dt;
	smart_assign(dt,zM_[1]);
	if(!is_well_defined_value(dt))
	  ADONIS_ERROR(DerivedError,"zM_[1] isn't well defined. zM_[1] = "<<zM_[1]<<".");
	
	smart_assign(dt,z1_z2);
	if(!is_well_defined_value(dt))
	   ADONIS_ERROR(DerivedError,"z1_z2 isn't well defined. z1_z2 = "<<z1_z2<<".");
#endif
	std::cout << "z_2 - z_1/(1+a·z_1) = " << z1_z2 << std::endl;
	std::cout << "==========================================="<<std::endl;

	//std::cout <<"c_/eps_ = "<< c_/eps_ << std::endl;  //o.k.

	rhs_[i] = c_/eps_*z1_z2 - d_*y[i] + D1_*central_fd_2nd(y,i,h); 
      }
      //======================================================================
      
      return rhs_;
    }
    
    //! fake time dependent operator
    template<class V>
    VType& operator()(const T& time, const V& y){
      return (rhs_ = (*this).operator()(y));
    }
    
  };

}

#endif
