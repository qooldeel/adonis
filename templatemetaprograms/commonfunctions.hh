#ifndef FREQUENTLY_USED_FUNCTIONS_DEFINED_AS_META_PROGRAMS_HH
#define FREQUENTLY_USED_FUNCTIONS_DEFINED_AS_META_PROGRAMS_HH

#include <iostream>

#include "conditionalstructures.hh"

namespace Adonis{

  
  /**
   * \brief compile-time recursion to compute \f$ N!\f$ 
   */
  template<int N>
  class Factorial{
    public:
    enum{result = N*Factorial<N-1>::result};
    
  };

  
  //specialisation -- stop recursion
  template<>
  class Factorial<0>{
  public:
    enum{result = 1};
  };
  

  /**
   *\ brief Compile time computation of Fibonacci numbers
   */
  template<int N>
  class Fibonacci{
  public:
    enum{fib = Fibonacci<N-1>::fib + Fibonacci<N-2>::fib};
  };
  
  //specialisations for n= 1... 
  template<>
  class Fibonacci<1>{
  public:
    enum{fib = 1};
  };

  //... and n = 0
  template<>
  class Fibonacci<0>{
  public:
    enum{fib = 0};
  };


 
  /**
   * Compile time Kronecker delta
   */
  template<int I, int J> 
  class Kronecker{
  public:                        //or: IfElse<((I==J)?true:false)>::Value;
    static inline bool delta() {return IfElse<I==J>::Value;}
  };

  
  /**
   * \brief Dense symmatrix 'store by row' access (store upper triangular part).
   * 'store by column' (store lower triangular part) is accomplished by just
   *  interchanging <TT> I </TT> and <TT> J </TT>, i.e. <TT>
   * SymmetricDenseMatrixAccess<J,I,N>::offset() </TT>, but this should be of
   * no consequence since we're dealing with symmetric matrices
   */
  template<size_t I, size_t J, size_t N> //!N is the order of the sym. matrix
  class SymmetricDenseMatrixAccess{
  public:
    static inline int offset(){
      return ((J>=I) ? I*N - (I+1)*I/2 + J : SymmetricDenseMatrixAccess<J,I,N>::offset());
    }
  };


  //! Compile time access on elements of a dense matrix when it is stored in 
  //! linear memory
  template<size_t I, size_t J, size_t ROWS, size_t COLS, char SF> class MatrixStorageFormat4LinearMemory;

  //! row-major order, select via 'r' or 'R'
  template<size_t I, size_t J, size_t ROWS, size_t COLS>
  class MatrixStorageFormat4LinearMemory<I,J,ROWS,COLS,'r'>{
  public:
    static inline size_t offset(){return I*COLS+J;}
  };

  template<size_t I, size_t J, size_t ROWS, size_t COLS>
  class MatrixStorageFormat4LinearMemory<I,J,ROWS,COLS,'R'>{ 
  public:
    static inline size_t offset(){return MatrixStorageFormat4LinearMemory<I,J,ROWS,COLS,'r'>::offset();}
  };

  
  //! column-major order, select via 'c' or 'C'
  template<size_t I, size_t J, size_t ROWS, size_t COLS>
  class MatrixStorageFormat4LinearMemory<I,J,ROWS,COLS,'c'>{
  public:
    static inline size_t offset(){return I + J*ROWS;}
  };

  template<size_t I, size_t J, size_t ROWS, size_t COLS>
  class MatrixStorageFormat4LinearMemory<I,J,ROWS,COLS,'C'>{
  public:
    static inline size_t offset(){return MatrixStorageFormat4LinearMemory<I,J,ROWS,COLS,'c'>::offset();}
  };

}//end namespace 



#endif
