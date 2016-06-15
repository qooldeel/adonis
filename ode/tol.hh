#ifndef USER_DEFINED_TOLERANCE_HH
#define USER_DEFINED_TOLERANCE_HH

#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"

namespace Adonis{

  /**
   * \brief User defined tolerance
   * \tparam NORM the norm object (1-norm, 2-norm)
   */
  template<class NORM>
  class UserDefinedTolerance{
  public:
    typedef typename NORM::return_type return_type;
    
  private:
    return_type tol_;
    
  public:
    /**
     *\brief set the tolerance. I choose a default value of 1.e-12. 
     *In general: tol_ should be chosen well above the machine round-off error, 
     * e.g. 
     *     1.e04*std::numeric_limits<return_type>::epsilon().x 
     */
    UserDefinedTolerance(const return_type& tol = return_type(1.e-12)):tol_(tol){}
    
    //! get tolerance 
    return_type& tol() {return tol_;}
    const return_type& tol() const{return tol_;}


    /**
     * \brief calculates norm
     * \tparam VType random access container (rac's)
     * NOTE: relies on definition of subtraction of two racs!!
     */
    template<class VType>
    inline return_type norm(const VType& prev, const VType& now) const{
      adonis_assert(prev.size() == now.size());

      VType res(prev.size());  
      res = now - prev;
      return NORM::norm(res);
    }
    

    /**
     *\brief calculates norm of a given random access container
     */
    template<class VType>
    inline return_type norm(const VType& v) const{
      return NORM::norm(v);
    }
    

    /**
     * decides whether an iteration is to be terminated via \f[ \|y_{n+1}-y_{n}\| \leq NTOL, \f] where \f$ \| \cdot \|.\f$
     * \tparam VType random access type (rac)
     * \param prev the previous iteration 
     * \param now the actual iteration
     */
    template<class VType>
    inline bool terminate_iteration(const VType& prev, const VType& now) const{
      return ( (norm(prev,now) <= tol_) ? true : false );
    }

  
    /**
     * \brief another criterion when to terminate an iteration by checking the norms of previous random access containers
     */
    inline bool terminate_iteration(const return_type& normPrev, const return_type& normNow) const{
      return ( (Abs(normNow-normPrev) <= tol_) ? true : false);
    }

    /**
     *\brief Can be used when, e.g., the solution of the Newton system \f$\delta := y_n^{\nu+1} - y_n^{\nu} \f$ is below a user-specified threshold.  
     */
    template<class VType>
    inline bool terminate_iteration(const VType& v) const{
      // std::cout << "NORM NORM NORM NORM NORM ||delta|| = " << norm(v) << std::endl;
      return  ( (norm(v) <= tol_) ? true : false);
    }

 
    //! overloaded function -- assignment of a result 
    template<class VType>
    inline bool terminate_iteration(const VType& residual, VType& u, const VType& ul) const{
      u = ul;
      return  ( (norm(residual) <= tol_) ? true : false);
    }

  };
}

#endif
