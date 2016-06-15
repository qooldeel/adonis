#ifndef BOUNDARY_CONIDTION_HANDLER_FDM_2D_HH
#define BOUNDARY_CONIDTION_HANDLER_FDM_2D_HH

#include "../ode/repair.hh"
#include "../fdm/gensettings.hh"
namespace Adonis{
  

  //! ================= ZERO NEUMANN BOUNDARY CONDITION ====================
  //!just take inlet values everywhere
  template<class FCTR>
  inline typename FCTR::value_type mech_quantity_value(FCTR& fct, int quant){
    typedef typename FCTR::value_type value_type;
    value_type ret;
    if(quant==0)
	ret = fct.rho_in();
    else if(quant==1)
	ret = 
#ifdef NONCONSERVATIVE_FORM
	  fct.v1_in();
#else
    fct.v1_in()*fct.rho_in();
#endif
    else if(quant==2)
      ret =
#ifdef NONCONSERVATIVE_FORM
	fct.v2_in();
#else
    fct.v2_in()*fct.rho_in();
#endif
    else if(quant==3)
      ret = 
#ifdef NONCONSERVATIVE_FORM
	fct.T_in();
#else
    fct.T_in()*fct.rho_in();
#endif
    else if(quant>=4)
      ret = 
#ifdef NONCONSERVATIVE_FORM	
	fct.mixture_in()[quant];
#else
    fct.mixture_in()[quant]*fct.rho_in();
#endif
    else{}
  
    return ret;
  }

/**
 * \brief Equal treatment of primitive and conservative variables (no destinction required)
*/
  template<int N, class FCTR> class BoundaryFunction2D;

  template<class FCTR>
  class BoundaryFunction2D<0,FCTR>{ //specialization for Dirichlet
  public:
    typedef typename FCTR::value_type value_type;
    typedef typename FCTR::VType VType;

    static inline value_type down(FCTR& fct, int i, const VType& u, int quant){
      return mech_quantity_value<FCTR>(fct,quant);
    }
    
    static inline value_type up(FCTR& fct, int i, const VType& u, int quant){
      return mech_quantity_value<FCTR>(fct,quant);
    }

    static inline value_type left(FCTR& fct, int j, const VType& u, int quant){
      return mech_quantity_value<FCTR>(fct,quant);
    }

    static inline value_type right(FCTR& fct, int j, const VType& u, int quant){
      return mech_quantity_value<FCTR>(fct,quant);
    }

    //boundary pressure handling
    static inline value_type p_bdy_down(FCTR& fct, int i, const VType& u){
      return fct.boundary_pressure_0th_order(i,u);
    }
    
    static inline value_type p_bdy_up(FCTR& fct, int i, const VType& u){
      return fct.boundary_pressure_0th_order(i,u);
    }

    static inline value_type p_bdy_left(FCTR& fct, int i, const VType& u){
      return fct.boundary_pressure_0th_order(i,u);
    }

    static inline value_type p_bdy_right(FCTR& fct, int i, const VType& u){
      return fct.boundary_pressure_0th_order(i,u);
    }
    
  };

  template<class FCTR>
  class BoundaryFunction2D<1,FCTR>{ //1st order specialization
  public:
    typedef typename FCTR::value_type value_type;
    typedef typename FCTR::VType VType;
    
    static inline value_type down(FCTR& fct, int i, const VType& u, int quant){
      return fct.first_order_ZERO_Neumann_dy_down(i,u,quant);
    }
    
    static inline value_type up(FCTR& fct, int i, const VType& u, int quant){
      return fct.first_order_ZERO_Neumann_dy_up(i,u,quant);
    }

    static inline value_type left(FCTR& fct, int j, const VType& u, int quant){
      return fct.first_order_ZERO_Neumann_dx_left(j,u,quant);
    }

    static inline value_type right(FCTR& fct, int j, const VType& u, int quant){
      return fct.first_order_ZERO_Neumann_dx_right(j,u,quant);
    }
    
     //boundary pressure handling
    static inline value_type p_bdy_down(FCTR& fct, int i, const VType& u){
      return fct.boundary_pressure_1st_order_down(i,u);
    }
    
    static inline value_type p_bdy_up(FCTR& fct, int i, const VType& u){
      return fct.boundary_pressure_1st_order_up(i,u);
    }

    static inline value_type p_bdy_left(FCTR& fct, int j, const VType& u){
      return fct.boundary_pressure_1st_order_left(j,u);
    }

    static inline value_type p_bdy_right(FCTR& fct, int j, const VType& u){
      return fct.boundary_pressure_1st_order_right(j,u);
    }
  };


  template<class FCTR>
  class BoundaryFunction2D<2,FCTR>{ //2nd order specialization
  public:
    typedef typename FCTR::value_type value_type;
    typedef typename FCTR::VType VType;
    
    static inline value_type down(FCTR& fct, int i, const VType& u, int quant){
      return fct.second_order_ZERO_Neumann_dy_down(i,u,quant);
    }
    
    static inline value_type up(FCTR& fct, int i, const VType& u, int quant){
      return fct.second_order_ZERO_Neumann_dy_up(i,u,quant);
    }

    static inline value_type left(FCTR& fct, int j, const VType& u, int quant){
      return fct.second_order_ZERO_Neumann_dx_left(j,u,quant);
    }

    static inline value_type right(FCTR& fct, int j, const VType& u, int quant){
      return fct.second_order_ZERO_Neumann_dx_right(j,u,quant);
    }

     //boundary pressure handling
    static inline value_type p_bdy_down(FCTR& fct, int i, const VType& u){
      return fct.boundary_pressure_2nd_order_down(i,u);
    }
    
    static inline value_type p_bdy_up(FCTR& fct, int i, const VType& u){
      return fct.boundary_pressure_2nd_order_up(i,u);
    }

    static inline value_type p_bdy_left(FCTR& fct, int j, const VType& u){
      return fct.boundary_pressure_2nd_order_left(j,u);
    }

    static inline value_type p_bdy_right(FCTR& fct, int j, const VType& u){
      return fct.boundary_pressure_2nd_order_right(j,u);
    }
    
  };
  

} //end namespace

#endif
