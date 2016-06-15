#ifndef AVOID_SINGULARITIES_HH
#define AVOID_SINGULARITIES_HH

#include <iostream>
#include <limits>

#include "../templatemetaprograms/unrollloop.hh"

namespace Adonis{
  
  /**
   * /brief Defines a small yet numerically insignificant number \f$ \epsilon > 0 \f$
   */
  template<class T>
  class SmallNumericallyInsignificant{
  public:
    static const T Number;
  };

  template<class T> const T SmallNumericallyInsignificant<T>::Number = 1.e+04;


  /**
   * \brief Performs or doesn't perform \f$ X_k \leftarrow X_k + \epsilon, \quad \epsilon > 0\f$
   */
  template<bool B, int N, class BTYPE> class AvoidSingularities;

  template<int N, class BTYPE>
  class  AvoidSingularities<true,N,BTYPE>{
  public:
    //! each  \f$ X_k\f$ gets its own \f$delta_k\f$ with randomly chosen base
    template<class ITER1, class ITER2>
    static inline void in_computation(ITER1 it, ITER2 delta){
      UnrollLoop<0,N>::random(delta,1.e+04*std::numeric_limits<BTYPE>::epsilon());
      UnrollLoop<0,N>::calibrate_against_singularity_multiple_perturbations(it,delta);
    }

    //! one delta for all \f$ X_k\f$
    template<class ITER>
    static inline void in_computation(ITER it){
      UnrollLoop<0,N>::calibrate_against_singularity(it,SmallNumericallyInsignificant<BTYPE>::Number*std::numeric_limits<BTYPE>::epsilon());
    }
  };

  template<int N, class BTYPE>
  class  AvoidSingularities<false,N,BTYPE>{
  public:
    template<class ITER1, class ITER2>
    static inline void in_computation(ITER1 it, ITER2 perturb){} //do nothing

    template<class ITER>
    static inline void in_computation(ITER it){} //do nothing
  };
  

}

#endif
