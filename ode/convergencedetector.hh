#ifndef DETECT_CONVERGENCE_WRT_TRANSIENT_DIFFERENTIAL_EQUATION_HH
#define DETECT_CONVERGENCE_WRT_TRANSIENT_DIFFERENTIAL_EQUATION_HH

#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../common/adonisassert.hh"

namespace Adonis{


  /**
   * \brief checks convergence of timestepping method by comparing the difference of entries of the current and last timestep solution
   */
  template<class V, bool CHECK, int DIM> 
  class ConvergenceDetector{
  public:
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef std::size_t SizeType;
    
    static inline bool is_converged(const V& y_n, const V& y_prev, int numOfQuants, const BaseType& eps){ return false;}  //do nothing -- always "not converged"
  
    static inline void info() {std::cout << "No additional convergence check applied."<<std::endl;}
  };


  //!partial specialization -- true general case
  template<class V, int DIM> 
   class ConvergenceDetector<V,true,DIM>{
   public:
     typedef typename V::value_type value_type;
     typedef typename TypeAdapter<value_type>::BaseType BaseType;
     typedef std::size_t SizeType;

     static inline bool is_converged(const V& y_n, const V& y_prev, int numOfQuants, const BaseType& eps){
       adonis_assert((y_n.size()%numOfQuants == 0) && (y_prev.size()%numOfQuants == 0) && (y_n.size() == y_prev.size()));
       adonis_assert(eps > BaseType());

       bool res(true);
       //!check if all quantities from one timestep to the next differ by less than a prescribed tolerance
       for(SizeType l = 0; l < y_n.size(); ++l){
	 if(Abs(y_n[l]-y_prev[l]) < eps){
	   res &= true;
	 }
	 else{
	   res &= false;   //once false the output is false "not converged"
	   break;          //now you can leave the loop 
	 }
       }
       return res;

     }

    static inline void info() {std::cout << "ADDITIONAL convergence check applied for MOL problems of dim = "<< DIM <<"."<<std::endl;}
  };

  //!partial specialization -- 2D case
   template<class V> 
   class ConvergenceDetector<V,true,2>{
   public:
     typedef typename V::value_type value_type;
     typedef typename TypeAdapter<value_type>::BaseType BaseType;
     typedef std::size_t SizeType;

     static inline bool is_converged(const V& y_n, const V& y_prev, int numOfQuants, const BaseType& eps){
       adonis_assert((y_n.size()%numOfQuants == 0) && (y_prev.size()%numOfQuants == 0) && (y_n.size() == y_prev.size()));
       adonis_assert(eps > BaseType());

       bool res(true);
       //!assumption: rho stored at the very beginning!!! 
       //check if rho does not change at each point in space between timesteps
       //cf. [ANDERSON, p. 461]
       SizeType dim = y_n.size()/((SizeType)numOfQuants);
       for(SizeType l = 0; l < dim; ++l){
	 if(Abs(y_n[l]-y_prev[l]) < eps){
	   res &= true;
	 }
	 else{
	   res &= false;   //once false the output is "not converged"
	   break;   //now you can leave the loop
	 }
       }
       return res;

     }

     static inline void info() {std::cout << "ADDITIONAL convergence check applied for MOL problems of dim = 2."<<std::endl;}
     
  };

  

}  //end namespace

#endif
