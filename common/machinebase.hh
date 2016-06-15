#ifndef MACHINE_BASE_HH
#define MACHINE_BASE_HH

#include "adonisassert.hh"
#include "globalfunctions.hh"

namespace Adonis{

  /**
   * \brief Machine Base. 
   * MachineBase<2,double>::base() and MachineBase<2,double>::rev_base() returns the Machine Base and the Reverse Machine Base respectively.
   */
  template<size_t B = 2, class T = double>
  class MachineBase{
  public:
    typedef size_t IndexType;
    typedef T value_type;
    
    //no const's here of any kind for cv-qualifier aren't supported here
    static inline size_t base() {return B;}
    static inline T rev_base()  {return 1./static_cast<T>(B);}

  };

  


  template<class T, size_t B>
  inline T pow_of_machine_base(int n){
    adonis_assert(MachineBase<B>::base() != 0);
    
    if(n == 0)
      return 1;    //since 2^0 and 0.5^0 are one, respectively 

    T mbase, res;
    
    if(n < 0){
      n = std::abs(n);
      mbase = MachineBase<B>::rev_base();
      res = mbase;
    }
    else{
      mbase = MachineBase<B>::base();
      res = mbase;
    }
    for(int i = 1; i < n; ++i){
      res *= mbase;
    }
    return res;
  }



}

#endif
