#ifndef COMPILE_TIME_ARRAY_HH
#define COMPILE_TIME_ARRAY_HH

#include <iostream>

namespace Adonis{
 
  //! do nothing 
  class Self{
  public:
    typedef Self Type;
    enum{value};
  };
  
  //! Array will be, e.g., <TT> Element<N,Element<M,Element<P> > >  </TT>
  /**
   * \brief Compile-time array
   * USAGE:
   *\code 
      typedef Element<0, Element<3, Element<6, Element<7> > > > ArrayType;

      // element access
      cout << AccessElement<Array, 0>::result << endl;
      cout << AccessElement<Array, 1>::result << endl;
      cout << AccessElement<Array, 2>::result << endl;
      cout << AccessElement<Array, 3>::result << endl;
    
   *\endcode 
   *
   * Adapted code from 
   * <a href="http://daniel-albuschat.blogspot.com/2009/06/compile-time-array-in-c-using-templates.html">Compile-time array in C++ using templates</a>!
   * 
   */
  template<int I, class T = Self>
  class Element{
  public:
    typedef T Type;
    enum { value = I };
  };


  /**
   * \tparam A underlying <BB> integral-type </BB> array 
   * \tparam INDEX index in array
   * enums force a purce compile-time execution, while <TT> static const int resutl =  AccessElement<typename Array::Type,INDEX-1>::result; </TT> forces the compiler to instantiate and allocate the definition for the static member, yielding a bad "compile-time" execution by virtue of inherent lvalues associated with static constants.
   */
  template<class A, int INDEX>
  class AccessElement{
  public:
    enum{ result = AccessElement<typename A::Type,INDEX-1>::result };
  };
  
  template<class A>
  class AccessElement<A, 0>{
  public:
    enum{ result  = A::value };
  };



} //end namespace 

#endif
