#ifndef MISCELLANEOUS_STUFF_FOR_MY_ODE_SOLVERS_HH
#define MISCELLANEOUS_STUFF_FOR_MY_ODE_SOLVERS_HH

namespace Adonis{
  
  
  template<int DIM, class ITER>
  class ReadOutCompositionFromODEScheme{
  public:
    typedef ITER IteratorType;

    template<class FUN, class T>
    static inline void get_previous(FUN& fun, const ITER y_prevIt, const T& stepsize){}
  };

  //! partial specialization for 2D
  template<class ITER>  
  class ReadOutCompositionFromODEScheme<2,ITER>{
  public:
     typedef ITER IteratorType;

    template<class FUN, class T>
    static inline void get_previous(FUN& fun, const ITER y_prevIt, const T& stepsize){
      fun.get_y_prev(y_prevIt,stepsize);
    }
  };


  /**
   * \brief Boundary treatment
   */
  template<class FUN, bool B> class BoundaryFromFunctor;

  template<class FUN>
  class BoundaryFromFunctor<FUN,true>{
  public:
    template<class V>
    static inline void set_boundary(V& u, FUN& fun){
      fun.set_boundary(u);
    }
  };

  //! do nothing
  template<class FUN>
  class BoundaryFromFunctor<FUN,false>{
  public:
    template<class V>
    static inline void set_boundary(V& u, FUN& fun){}
  };
  
  
} //end namespace 

#endif //include guard
