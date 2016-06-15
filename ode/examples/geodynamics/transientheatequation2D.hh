#ifndef TRANSIENT_HEAT_EQUATION_IN_2D_FOR_THE_EARTH_S_LITHOSPHERE_HH
#define TRANSIENT_HEAT_EQUATION_IN_2D_FOR_THE_EARTH_S_LITHOSPHERE_HH

#include "../../../common/typeadapter.hh"
#include "../../../expressiontemplates/exprvec.hh"
#include "../../../io/readinparameters.hh"
#include "../../../misc/misctmps.hh"

namespace Adonis{

  /**
   * Solve the following differential equation numerically:
   * \f[ \partial_t T = \kappa(\partial_{xx} T + \partial_{yy} T) + \frac{Q}{\rho c_p}\f]
   */
  template<class T>
  class GeoTransientHeatEquation{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;   
    typedef typename TypeAdapter<T>::Type DType;
    typedef std::size_t SizeType;

    GeoTransientHeatEquation(SizeType n = 0):rhs_(n), kappa_(T()),Tsurf_(T()),Qsurf_(T()),rec_(T()),npt_(0),Nx_(0),Ny_(0),hx_(DType()),hy_(DType()),depth_(DType()){
      ParameterData PD;
      PD.read_from_file("examples/geodynamics/data2readin.dat");
      kappa_ = PD.get_datum<T>("kappa");
      Tsurf_ = PD.get_datum<T>("Tsurf");
      Qsurf_ = PD.get_datum<T>("Qsurf");
      rec_ = Reciprocal(PD.get_datum<T>("rho")*PD.get_datum<T>("cp"));
      Nx_ = PD.get_datum<SizeType>("Nx");
      Ny_ = PD.get_datum<SizeType>("Ny");
      hx_ = (PD.get_datum<DType>("b")-PD.get_datum<DType>("a"))/Nx_;
      depth_ = PD.get_datum<DType>("d")-PD.get_datum<DType>("c");
      hy_ = depth_/Ny_;
      npt_ = Nx_*Ny_;
      
      T_.resize(npt_);
    }

    const SizeType domain_dim() const {return rhs_.size();}
    const SizeType dim() const {return rhs_.size();}

    template<class U>
    VType& operator()(const U& u){
      T_ = u;   //hard copy: current temperature profile
      
      //BOUNDARY CONDITIONS ---------------------------------------------------
      // UP is constant, DOWN alters like u
      for(SizeType i = 0; i  < Nx_; ++i){
	//T_[offset(i,0)] = u[offset(i,0)]; //o.k. T_ is u anyway
	T_[offset(i,Ny_-1)] = Tsurf_; 
	//rhs_[offset(i,Ny_-1)] = 0;   //constant d/dt·
      }
     
      //! Zero boundary Neumann condition, 2nd order, substitute ghost nodes
      //! into discretization scheme
      for(SizeType j = 1; j < Ny_-1; ++j){
      	//!LEFT  
	rhs_[offset(0,j)] = kappa_*( (2*T_[offset(1,j)] - 2*T_[offset(0,j)])/ntimes<2>(hx_) + (T_[offset(0,j+1)] - 2*T_[offset(0,j)] + T_[offset(0,j-1)])/ntimes<2>(hy_)); 

	
	//!first order NEUMANN -- higher disrc. error at boundary
	// T_[offset(0,j)] = T_[offset(1,j)]; //1st order
	// rhs_[offset(0,j)] = 0;
 
      	//!RIGHT
      	rhs_[offset(Nx_-1,j)] = kappa_*( (2*T_[offset(Nx_-2,j)] - 2*T_[offset(Nx_-1,j)])/ntimes<2>(hx_) + (T_[offset(Nx_-1,j+1)] - 2*T_[offset(Nx_-1,j)] + T_[offset(Nx_-1,j-1)])/ntimes<2>(hy_) );
      
	//!first order NEUMANN -- higher disrc. error at boundary
	// T_[offset(Nx_-1,j)] = T_[offset(Nx_-2,j)]; //1st order
	// rhs_[offset(Nx_-1,j)] = 0;
      }

      //!INTERIOR -------------------------------------------------------------
      for(SizeType i = 1; i < Nx_-1; ++i){
      	for(SizeType j = 1; j < Ny_-1; ++j){
      	  rhs_[offset(i,j)] = kappa_*( (T_[offset(i+1,j)] - 2*T_[offset(i,j)] + T_[offset(i-1,j)])/ntimes<2>(hx_) + (T_[offset(i,j+1)] - 2*T_[offset(i,j)] + T_[offset(i,j-1)])/ntimes<2>(hy_));
      	}
      }


      return rhs_;
    }

    //dummy function
    void get_y_prev(typename VType::iterator y_prevIt, const T& h){}

  private:
    VType rhs_;
    T kappa_, Tsurf_, Qsurf_, rec_;
    SizeType npt_, Nx_, Ny_;
    DType hx_, hy_,depth_;
    VType T_;
  
    SizeType offset(SizeType i, SizeType j){
      return (i + Nx_*j);
    }

    T radiogenic_heat_production(const T& y){
      return ( Qsurf_*exp(-y/depth_) );
    }

  };


} //end namespace

#endif
