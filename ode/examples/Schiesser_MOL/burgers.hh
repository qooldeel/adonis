#ifndef BURGERS_EQUATION_SCHIESSER_1D_HH
#define BURGERS_EQUATION_SCHIESSER_1D_HH

#include "../../../common/error.hh"
#include "../../../expressiontemplates/exprvec.hh"
#include "../../../containers/staticarray.hh"
#include "../../../misc/misctmps.hh"
#include "../../../io/readinparameters.hh"
#include "../../functorid.hh"

namespace Adonis{

  /**
   * \brief Analytical solution of Burger's equation, see [SCHIESSER and GRIFFITHS, "A Compendium of PDEs Models - Method of Lines Analysis with Matlab", chap. 5]
   */
  template<class T>
  inline T analytical_burger_s_equation(const T& x, const T& t, const T& vis){
    T a = (0.05/vis)*(x-0.5+4.95*t),
      b = (0.25/vis)*(x-0.5+0.75*t),
      c = ( 0.5/vis)*(x-0.375),
      ea = exp(-a),
      eb = exp(-b),
      ec = exp(-c);
      return ( (0.1*ea+0.5*eb+ec)/(ea+eb+ec) );
  }

  /**
   * Taken from [SCHIESSER and GRIFFITHS, "A Compendium of PDEs Models - Method of Lines Analysis with Matlab", chap. 5]
   * Solve <i>Burger's equation</i> equation \f[ \partial_t u + u \partial_x u - \nu \partial_{xx}^2 u = 0, \qquad \nu > 0\f] with corresponding BCs and ICs. The value \f$ \nu \f$ is the (constant) kinematic viscosity.
   * Burger's equation is a good example of a nonlinear hyperbolic (conservation) equation for testing numerical methods.
   */
  template<class T>
  class OriginalBurgersEquationSchiesser{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;

    OriginalBurgersEquationSchiesser(int n = 0):rhs_(n),ux_(n), uxx_(n){
      ParameterData PM;
      //! use absolute path here!
      PM.read_from_file("/home/mfein/MARC++/ode/datafiles/burgers.dat");
      Nx_ = PM.get_datum<int>("Nx");
      xl_ = PM.get_datum<double>("a");
      xu_ = PM.get_datum<double>("b");
      vis_ = PM.get_datum<double>("vis");
      delta_x_ = (xu_-xl_)/(Nx_-1);
    }

    std::size_t dim() const {return rhs_.size();}
    std::size_t domain_dim() const {return rhs_.size();}

    //! this follows closely the implementation of Schiesser. Note however that
    //! it takes much more time since firstly \f$\partial_x u\f$ and then 
    //! \f$\partial_{xx}^2\f$ are computed over the whole domain!
    template<class X>
    VType& operator()(const X& u){
      //calculate ux_
      (*this).dss002(ux_, xl_,xu_,Nx_,u);
      //BC at x = 0
      ux_[0] = 0.0;
      //BC at x = 1
      ux_[Nx_-1] = 0.0;
      //calculate uxx_
      (*this).dss002(uxx_, xl_,xu_,Nx_,ux_);
      //!PDE -- now over the whole domain because of the boundary treatment
      //!       in private function dss002 
      for(int i = 0; i < Nx_; ++i){
	rhs_[i] = vis_*uxx_[i] - u[i]*ux_[i];
      }
	
	
      return rhs_;
    }

 
  private:
    VType rhs_, ux_, uxx_;
    int Nx_;
    T xl_, xu_, vis_, delta_x_;
    
    //private function
    VType& dss002(VType& ux, const T& xl, const T& xu, int n, const VType& u){
      T dx = (xu-xl)/(n-1),
      r2fdx = 1./(2.*dx);
      int nm1 = n-1;
      
      ux[0] = r2fdx*(-3.   *u[0]     +4.   *u[1]     -1.   *u[2]);
      
      for(int i = 1; i < nm1; ++i){
	ux[i] = r2fdx* (-1.   *u[i-1]     +0.   *u[i]     +1.   *u[i+1]);
      }

      ux[n-1] = r2fdx*( 1.   *u[n-3]     -4.   *u[n-2]     +3.   *u[n-1]);

      return ux;
    }
    
    
  };





  /**
   * \brief My implementation of Burger's equation using the solution \f$u\f$
   *        explicitly for the boundary value. In this case no ODEs are 
   *        required at the boundary, .i.e. rhs_|bdy = 0.
   */
  template<class T>
  class BurgersEquationSchiesser{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;

    //! handle set_boundary in integrator
    typedef MOLFunctorIdentifier<1,FunctorID::burgerMarc> mol_type;

    BurgersEquationSchiesser(int n = 0):rhs_(n),u_(n){
      ParameterData PM;
      //! use absolute path here!
      PM.read_from_file("/home/mfein/MARC++/ode/datafiles/burgers.dat");
      Nx_ = PM.get_datum<int>("Nx");
      xl_ = PM.get_datum<double>("a");
      xu_ = PM.get_datum<double>("b");
      vis_ = PM.get_datum<double>("vis");
      delta_x_ = (xu_-xl_)/(Nx_-1);
    }

    std::size_t dim() const {return rhs_.size();}
    std::size_t domain_dim() const {return rhs_.size();}


    template<class X>
    void set_boundary(X& w){
      //BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC
      //BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC
      //BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC · BC
      
      // explicit solution w.r.t. u_left and u_right, respectively, requires
      // no ODE for the boundary points --> rhs|bdy = 0.0
      //  LEFT, Zero Neumann
      w[0] = 1./3.*(4.*w[1] - w[2]);
      rhs_[0] = 0.0;
      //  RIGHT, Zero Neumann
      w[Nx_-1] =  1./3.*(-w[Nx_-3] + 4.*w[Nx_-2]);
      rhs_[Nx_-1] = 0.0;
    }
   
    template<class X>
    VType& operator()(const X& x){
      //hard copy
      u_ = x;

      (*this).set_boundary(u_);

      //INTERIOR: \f$ \partial_t u = \underbrace{\frac{\mu}{\rho}}_{:= \mathrm{vis}} \ cdot \partial_{xx}^2 u - u\partial_x u  \f$
      for(int i = 1; i < Nx_-1; ++i){
	rhs_[i] = vis_*(u_[i+1] - 2*u_[i] + u_[i-1])/ntimes<2>(delta_x_) - u_[i]*(u_[i+1] - u_[i-1])/(2.*delta_x_);
      }

      return rhs_;
    }


 
  private:
    VType rhs_, u_;
    int Nx_;
    T xl_, xu_, vis_, delta_x_;
    
  };

} //end namespace 


#endif
