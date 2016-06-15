#ifndef L_ECUYER_32_AND_64_BIT_RANDOM_NUMBER_GENERATOR_HH
#define L_ECUYER_32_AND_64_BIT_RANDOM_NUMBER_GENERATOR_HH

#include "generalrng.hh"


namespace Adonis{

  /* LINEAR FEEDBACK SHIFT REGISTER (LFSR) RANDOM NUMBER GENERATORS  */
 
  /**
   * \brief A fast and reliable random number generator (RNG) due to L'Ecuyer [1, Figure 1] -- 32 bit version
   * Here we use an adapted version (in its original form, the algorithm's 
   * output is multiplied with 2.3283064365386963e-10
   *
   * In my version I take mod \f$ 2^31\f$ in order to get values between 0 and
   * \f$ 2^31-1\f$
   *
   * Reference:
   * [1] L'Ecuyer, "Tables of Maximally-Equidistributed Combined Lfsr Generators", Mathematics of Computation 68, no. 225, 1998
   */
  template<class T = double>
  class RNGlfsr113: public GeneralRNGBase<RNGlfsr113<T>,32>{
  public:
    typedef T value_type;
    typedef GeneralRNGBase<RNGlfsr113<T>,32> BaseType;
    typedef typename BaseType::IntType IntType;
  
    RNGlfsr113(IntType seed = std::time(0)){//seeded with system time
      set_seed(seed);
    }

    void set_seed(IntType seed){
      BaseType::seed_ = seed;
      //Note: initial seeds z1_, z2_, z3_, z4_  MUST be larger than 1, 7, 15, and 127, respectively. In case of violation, just take some (arbitrary) numbers which are slighly greater than the threshold values 
      z1_ = (seed < 1) ? 3 : seed;
      z2_ = (seed < 7) ? 11 : seed;
      z3_ = (seed < 15) ? 19 : seed;
      z4_ = (seed < 127) ? 129 : seed;
    }


    IntType draw_number(){
      (*this).shifter();
      return ((z1_ ^ z2_ ^ z3_ ^ z4_) % 2147483648); // values on [0,2^31-1]
    }


    //! gives values on [0,1)
    T draw_number_between_0_and_1(){
      (*this).shifter();
      return (z1_ ^ z2_ ^ z3_ ^ z4_) * 2.3283064365386963e-10;
    }

    T draw_number_from_range(const T& low, const T& up){
   
      return (low + (*this).draw_number_between_0_and_1()*(up-low));
    }

  private:
    IntType z1_, z2_, z3_, z4_; 

    void shifter(){
      IntType b;
      b  = ((z1_ << 6) ^ z1_) >> 13;
      z1_ = ((z1_ & 4294967294U) << 18) ^ b;
      b  = ((z2_ << 2) ^ z2_) >> 27;
      z2_ = ((z2_ & 4294967288U) << 2) ^ b;
      b  = ((z3_ << 13) ^ z3_) >> 21;
      z3_ = ((z3_ & 4294967280U) << 7) ^ b;
      b  = ((z4_ << 3) ^ z4_) >> 12;
      z4_ = ((z4_ & 4294967168U) << 13) ^ b;
    }
  };


  //! THIS GIVES OVERFLOWS --> 128 bit ints required
  /**
   * \brief 64 bit random number generator  as given in [1, Figure 2]
   *
   * [1] [Pierre L'Ecuyer, "Tables of maximally equidistributed combined LFSR generators"]
   */
  /*  //DOES NOT WORK PROPERLY SINCE TO LARGE CONSTS FOR UI64!!!!!!! (valgrind errs occur perhaps because operations in shifter() generate numbers beyond 64 bit unsigned ints)
  template<class T = double>
  class RNGlfsr258: public GeneralRNGBase<RNGlfsr258<T>,64>{
  public:
    typedef T value_type;
    typedef GeneralRNGBase<RNGlfsr258<T>,64> BaseType;
    typedef typename BaseType::IntType IntType;
  
    RNGlfsr258(IntType seed){
      BaseType::seed_ = seed;
      //Note: initial seeds z1_, z2_, z3_, z4_  MUST be larger than 1, 7, 15, and 127, respectively. In case of violation, just take some (arbitrary) numbers which are slighly greater than the threshold values 
      z1_ = (seed < 1) ? 3 : seed;  //is this needed here?
	z2_ = (seed < 7) ? 11 : seed;
	z3_ = (seed < 15) ? 19 : seed;
	z4_ = (seed < 127) ? 129 : seed;
    }



    IntType draw_number(){
      (*this).shifter();
      //on [0, 2^63-1]
      return ((z1_ ^ z2_ ^ z3_ ^ z4_ ^ z5_)%((IntType)9.223372037e+18) );
    }

    //! gives values on [0,1)
    T draw_number_between_0_and_1(){
      (*this).shifter();
      return ((z1_ ^ z2_ ^ z3_ ^ z4_ ^ z5_) * 5.4210108624275221e-20);
    }

  private:
    IntType z1_, z2_, z3_, z4_, z5_; 

    void shifter(){
      IntType b;
      b = (((z1_ << 1) ^ z1_) >> 53);
      z1_ = (((z1_ & 18446744073709551614LU) << 10) ^ b);
      b = (((z2_ << 24) ^ z2_) >> 50);
      z2_ = (((z2_ & 18446744073709551104LU) << 5) ^ b);
      b = (((z3_ << 3) ^ z3_) >> 23);
      z3_ = (((z3_ & 18446744073709547520LU) << 29) ^ b);
      b = (((z4_ << 5) ^ z4_) >> 24);
      z4_ = (((z4_ & 18446744073709420544LU) << 23) ^ b);
      b = (((z5_ << 3) ^ z5_) >> 33);
      z5_ = (((z5_ & 18446744073701163008LU) << 8) ^ b);
    }
  };
  */
} //end namespace 

#endif
