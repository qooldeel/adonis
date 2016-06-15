#ifndef A_RECURSIVE_RANDOM_ACCESS_CONTAINER_HH
#define A_RECURSIVE_RANDOM_ACCESS_CONTAINER_HH


/**
 * \brief Recursive STL-compliant random access container
 * \code
    // create a 3-dim matrix object 
    typedef RecursiveRandomAccessContainer<double,3>::VType Matrix3DType;
    * \endcode  
 */
namespace Adonis{

  template<class T, int I, template<class S, class A = std::allocator<S> > 
	   class V>                //  = ExprTmpl::MyVec>
  class RecursiveRandomAccessContainer{
  public:
    typedef typename RecursiveRandomAccessContainer<T,I-1,V>::VType RecType;
    typedef V<RecType,std::allocator<RecType> > VType;
  };


  //!Partial specialisation: end of recursion. VType = V<double>
  template<class T, template<class S, class A = std::allocator<S> > 
	   class V>
  class RecursiveRandomAccessContainer<T,1,V>{
  public:
    typedef V<T> VType;
  };


  //!Partial specialisation: Special case 0 just means that the grid is a scalar of type T
  template<class T, template<class S, class A = std::allocator<S> > 
	   class V>
  class RecursiveRandomAccessContainer<T,0,V>{
  public:
    typedef T VType;
  };


}  //end namespace 

#endif 
