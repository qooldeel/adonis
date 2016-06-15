#ifndef PRIMITIVE_VARIABLES_IN_A_FD_MOL_FLOW_HH
#define PRIMITIVE_VARIABLES_IN_A_FD_MOL_FLOW_HH

#include "gensettings.hh"

namespace Adonis{

  /**
   * \brief Calculate primitive variable
   *
   * Note that \f[ \rho\frac{D w}{Dt} = \partial_t(\rho w) + \nabla \cdot (\rho w v).\f] In particular \f[ \rho c_p\frac{D T}{Dt} = c_p (\partial_t(\rho T) + \nabla \cdot (\rho T v) ),\f] yielding \f$ \rho T\f$ rather than \f$ \rho c_p T\f$ as dependent variable for the thermal energy. The reason why I use \f$c_p\f$ outside the derivatives is simply the fact that \f$c_p\f$ is a polynomial in \f$T.\f$ Physically, this is o.k. as can be seen from the definition of the substantial derivative.
   */
  template<class U> 
  class FlowVariables{
  public:
    typedef typename U::value_type value_type;
    typedef FlowVariables<U> ThisType;

    //! note that 1D, 2D, 3D can be easily overloaded via rho(i,u), rho(i,j,u), etc.
    static inline value_type primitive(int i, int j, const U& u, int quant, int Nx, int Ny){
#ifdef NONCONSERVATIVE_FORM
      return u[i + Nx*j + quant*Nx*Ny];
#else  //conservative form; all vars are devided by rho except rho itself (quant = 0)
      //adonis_assert(u[i + Nx*j] > 0.);
      return ( (quant == 0) ? (u[i + Nx*j]) : (u[i + Nx*j + quant*Nx*Ny]/u[i + Nx*j]) );
#endif
    }
    
    //! gives back rho or 1, depending on whether nonconservative or conservative form of equations is to be considered
    static inline const value_type rho(int i, int j, const U& u, int Nx){
#ifdef NONCONSERVATIVE_FORM
      return value_type(1);
#else
      //adonis_assert(u[i + Nx*j] > 0.);
      return u[i + Nx*j];
#endif
    }

    //should be only invoked after function 'primitive'. val represents the raw value
    static inline void assign(int i, int j, U& u, int quant, int Nx, int Ny, const value_type& val){
      // adonis_assert(ThisType::rho(i,j,u,Nx) > 0.);
      u[i + Nx*j + quant*Nx*Ny] =
#ifdef NONCONSERVATIVE_FORM
      val;
#else  //rho*val (and val when quant is rho)
      ( (quant == 0) ? (val) : (val*ThisType::rho(i,j,u,Nx)) );
#endif   
    }
    

  };

} //end namespace 

#endif //include guard
