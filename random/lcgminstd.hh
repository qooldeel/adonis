#ifndef MINIMAL_32_BIT_LINEAR_CONGRUENTIAL_GENERATOR_MINIMAL_STANDARD_HH
#define MINIMAL_32_BIT_LINEAR_CONGRUENTIAL_GENERATOR_MINIMAL_STANDARD_HH

#include "generalrng.hh"

namespace Adonis{

  /**
   * \brief Minimal standard random number 32 bit generator [1,2]
   *
   *  [1] W. H. Press et al., "Numerical Recipes in C", Chap. 7, Cambridge University Press, 1992, p. 279
   *  [2] D. E. Knuth, "The Art of Computer Programming", Vol. 2, Addison-Wesley, 1997
   */
  template<class T = double>
  class LCGminstd: public GeneralRNGBase<LCGminstd<T>,32>{
  public:
    typedef T value_type;
    typedef GeneralRNGBase<LCGminstd<T>,32> BaseType;
    typedef typename BaseType::IntType IntType;
    
    LCGminstd(IntType seed = std::time(0)):k_(0){
      set_seed(seed);
    }

    void set_seed(IntType seed){
      BaseType::seed_ = seed;
    }


    IntType draw_number(){
      (*this).prepare();
      return BaseType::seed_;
    }

    T draw_number_between_0_and_1(){
      (*this).prepare();
      T res = static_cast<T>(BaseType::seed_)/2147483647U; //convert to fl.point
      BaseType::seed_ ^= 123459876U;	
      return res;
    }

    T draw_number_from_range(const T& low, const T& up){
      
      return (low + (*this).draw_number_between_0_and_1()*(up-low));
    }
    
  private:
    IntType k_;

    void prepare(){
      BaseType::seed_ ^= 123459876U;
      k_ = BaseType::seed_/127773U;
      BaseType::seed_ = 16807U*(BaseType::seed_-k_*127773U)-2836U*k_;
      if(BaseType::seed_ < 0) BaseType::seed_ += 2147483647U;
    }
  };

}

#endif
