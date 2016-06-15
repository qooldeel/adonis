#ifndef PHYSICAL_QUANTITY_TO_COMPUTE_WITH_HH
#define PHYSICAL_QUANTITY_TO_COMPUTE_WITH_HH

#include "../common/adonisassert.hh"
#include "../common/typeadapter.hh"

namespace Adonis{
  
  /**
   * \brief Relations between physical quantities
   *
   * chemcial source terms (the RHS of the equation \f$ \dot U(t) = F(U(t)\f$) 
   * must incorporate the species concentrations \f$ [X_k].\f$
   * However, in much applications, e.g. when using the Euler equations, one
   * uses partial densities \f$ \rho_i := \rho w_i\f$ in others specific moles
   * of species, \f$ \zeta_k.\f$
   * In those cases it is necessary to transform the quantities into 
   * concentrations and back.
   * The formula to do so is stated as \f[ [X_k] := \frac{n_k}{V} = \frac{n_k\rho}{\sum_{k=1}^K M_kn_k} = \frac{n_k\rho}{\sum_{k=1}^K m_k} = \frac{\rho w_k}{M_k} = \rho \zeta_k.\f]
   */
  template<char C>    
  class Quantity{   //! default case -- already in concentrations -- do nothing
  public:
    
    template<class T1, class T2>
    static inline const T1& concentration(const T1& uk, const T2& Mk, const T2& rho){
      return uk;
    }

    template<class T1, class T2>
    static inline const T1& back_to_input_quantity(const T1& Xk, const T2& Mk, const T2& rho){
      return Xk;
    }

    template<class V>
    static inline void resize(size_t n, V& v){ } //do nothing
  
    template<class U, class X, class MW, class T>
    static inline void concentration(size_t n, U& u, const X& x, const MW& mw, const T& rho){}  //do nothing
  
    template<class U, class MW, class T>
    static inline void back_to_input_quantity(size_t n, U& x, const MW& mw, const T& rho){} //do nothing

    template<class U, class MW, class T, class IX>
    static inline void back_to_input_quantity(size_t n, U& x, const MW& mw, const T& rho, const IX& index){} //do nothing


    //! x is already the concentration -- take it 
    template<class U, class W>
    static inline const U& choose_right_quantity(const U& input, const W& z){
      return input;
    }
  };

  //! specialisations
  //! <B>p</B>artial density (of species \f$k\f$) as input quantity
  template<>
  class Quantity<'p'>{ 
  public:
    
    typedef Quantity<'p'> ThisType;

    template<class T1, class T2>
    static inline T1 concentration(const T1& uk, const T2& Mk, const T2& rho){
      return uk/Mk; 
    }

    template<class T1, class T2>
    static inline T1 back_to_input_quantity(const T1& Xk, const T2& Mk, const T2& rho){
      return Xk*Mk;
    }
  
     template<class V>
     static inline void resize(size_t n, V& v){ 
       v.resize(n);
     } 
  
    template<class U, class X, class MW, class T>
    static inline void concentration(size_t n, U& u, const X& x, const MW& mw, const T& rho){
      for(size_t k = 0; k < n; ++k){
	u[k] = ThisType::concentration(x[k],mw[k],rho);
      }
    } 
  
    template<class U, class MW, class T>
    static inline void back_to_input_quantity(size_t n, U& x, const MW& mw, const T& rho){
      for(size_t k = 0; k < n; ++k){
	x[k] = ThisType::back_to_input_quantity(x[k],mw[k],rho);
      }
    }

    
    template<class U, class MW, class T, class IX>
    static inline void back_to_input_quantity(size_t n, U& x, const MW& mw, const T& rho, const IX& index){
      for(size_t k = 0; k < n; ++k){   //!n is index size now !
	x[k] = ThisType::back_to_input_quantity(x[k],mw[index[k]],rho);
      }
    }
    
    //! choose z which is the transformed quantity
    template<class U, class W>
    static inline const W& choose_right_quantity(const U& input, const W& z){
      return z;
    }

  };

  
  //!<B>s</B>pecific mole (of species \f$k\f$) as input quantity
  template<>
  class Quantity<'s'>{       
  public:

    typedef Quantity<'s'> ThisType;
    
    template<class T1, class T2>
    static inline T1 concentration(const T1& uk, const T2& Mk, const T2& rho){
      return rho*uk;
    }
    
    template<class T1, class T2>
    static inline T1 back_to_input_quantity(const T1& Xk, const T2& Mk, const T2& rho){
      return Xk/rho;
    }
  
    template<class V>
    static inline void resize(size_t n, V& v){ 
      v.resize(n);
    } 
  
    template<class U, class X, class MW, class T>
    static inline void concentration(size_t n, U& u, const X& x, const MW& mw, const T& rho){
      for(size_t k = 0; k < n; ++k){
	u[k] = ThisType::concentration(x[k],mw[k],rho);
      }
    } 
    
    template<class U, class MW, class T>
    static inline void back_to_input_quantity(size_t n, U& x, const MW& mw, const T& rho){
      for(size_t k = 0; k < n; ++k){
	x[k] = ThisType::back_to_input_quantity(x[k],mw[k],rho);
      }
    }

    template<class U, class MW, class T, class IX>
    static inline void back_to_input_quantity(size_t n, U& x, const MW& mw, const T& rho, const IX& index){
      for(size_t k = 0; k < n; ++k){   //!n is index size now !
	x[k] = ThisType::back_to_input_quantity(x[k],mw[index[k]],rho);
      }
    }

    //! choose z which is the transformed quantity
    template<class U, class W>
    static inline const W& choose_right_quantity(const U& input, const W& z){
      return z;
    }

  };


  
} //end namespace 

#endif
