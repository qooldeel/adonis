#ifndef IF_THEN_ELSE_HH
#define IF_THEN_ELSE_HH

//cf. [Vandevoorde/Josuttis, C++ Templates, §17, pp. 308] for this dashing trick
/*
  Replaces  ternary expressions in C/C++, for instance 

  if(x > 0) return x; else return -x;

  or equivalently:
  
  (x > 0) ? x : -x;

  Thus three arguments are involved ("if" <- "then", "else")
  The '<-' denotes the dependencies of the last two on the first argument according to the outcome of the conditional statement

  QESTION: What's all this Mumbo Jumbo about ????
  
  ANWSER: When used cleverly at the right places, it can help to increase performance tremendously! 
 */

namespace Adonis{
  
  //basic template. T1 and T2 depend on B
  template<bool B, typename T1, typename T2>
  class IfThenElse;

  //partial specialisation (the "if TRUE <- then" case)
  template<typename T1, typename T2>
  class IfThenElse<true, T1, T2>{
  public:
    typedef T1 ResultT;    //yields 2nd argument
  };

  //partial specialisation (the "if FALSE <- else" case)
  template<typename T1, typename T2>
  class IfThenElse<false, T1, T2>{
  public:
    typedef T2 ResultT;    //yields 3rd argument
  };


}//end namespace


#endif
