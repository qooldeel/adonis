#ifndef GENERAL_RANDOM_NUMBER_GENERATOR_BASE_CLASSES_HH
#define GENERAL_RANDOM_NUMBER_GENERATOR_BASE_CLASSES_HH

#include <ctime>
#include "../common/typeselector.hh"

namespace Adonis{

  /**
   * \brief andom number generation basis class (for 32 or 64 bit ints) using 
   * the Barton-Nackman trick. Seeming to be a little bit outdated, it still 
   * has a certain appeal of a fancy programming technique ;)
   */
  template<class OBJ, int INTTYPE>
  class GeneralRNGBase{
  public:
    //common functionalities
    typedef typename IntegerTypeSelection<INTTYPE>::IntegerType IntType;
    //typedef typename OBJ::value_type value_type;

    IntType get_seed(){return seed_;}
    const IntType get_seed() const{return seed_;}
    void reset_seed(IntType seed){seed_ = seed;}

  private:
    //delegate to derived class
    OBJ& ref_2_derived(){
      return static_cast<OBJ&>(*this); //give back derived class reference
    }

    //use properties of derived implementation
    IntType draw_number(){ 
      return ref_2_derived().draw_number(); //invoke derived class's member
    }                                       //has to be implemented in every
                                            //derived class 


                                   
  protected:
   IntType seed_;  //this is common to all derived classes
  
  };

} //end namespace

#endif
