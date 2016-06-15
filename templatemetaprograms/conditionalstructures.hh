#ifndef CONDITIONAL_CONTROL_STRUCTURES_THE_META_WAY_HH
#define CONDITIONAL_CONTROL_STRUCTURES_THE_META_WAY_HH


#include "../expressiontemplates/exprvec.hh"


//see how it works
namespace Adonis{

  //forward declaration
  template<bool C> class IfElse{};

  template<>
  class IfElse<true>{            //TRUE
  public:
    static const bool Value = 1;   //true value 

    static inline void statement(){
      //do something
      std::cout << "That's TRUE, hunk!" <<std::endl;
    }
  
    //overload fct 'statement'
    static inline void statement(const std::string& s){
      ADONIS_ERROR(DerivedError,s);
    }

  };



  template<>
  class IfElse<false>{          //FALSE
  public:
    static const bool Value = 0;  //false value

    static inline void statement(){
      //do something else 
      std::cout << " Well, not really..." <<std::endl;
    }
  
    static inline void statement(const std::string& s){ //do nothing
    }
  };


  //replacement for if/else statement via 
  // IfElse<condition>::statement();



  
  //#################### Switch statements ###################################

  template<int I>
  class Switch{
  public:
    static inline void s(){     //default: ....
      std::cout << " This is the default statement (when no template specialisation is supplied)." <<std::endl;
    }
  };
  
  template<>
  class Switch<0>{           //case 0: ....
  public:
    static inline void s(){
      std::cout << " Here you define Statement_0 "<< std::endl;
    }
  };


  template<>          
  class Switch<1>{           //case 1: ....
  public:
    static inline void s(){
      std::cout << " Here reads Statement_1. "<< std::endl;
    }
  };


  template<>
  class Switch<2>{           //case 2: .... 
  public:
    static inline void s(){
      std::cout << " To conclude, a 2nd statement is implemented as well. \n Actually this can arbitrarily be extended due to your flavour..."<< std::endl;
    }
  };

  //...et cetera

  //replacement for switch(i):
  //Switch<I>::s();    //I = 0, 1, 2, ...
}

#endif
