#ifndef MY_EXPRESSION_VECTOR_TRAITS_HH
#define MY_EXPRESSION_VECTOR_TRAITS_HH

namespace Adonis{

  template<class T>
  class ExprVecTraits{
  public:
    typedef ExprTmpl::MyVec<T> ExprVecType;
  };


} //end namespace 

#endif
