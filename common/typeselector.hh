#ifndef SELECT_PROPRIATE_TYPE_HH
#define SELECT_PROPRIATE_TYPE_HH

#include "../dunestuff/fmatrix.hh"
#include "../containers/fdiagonal.hh"

#include "types.hh"

namespace Adonis{
  
  //! select appropriate integer length 
  template<int NUMOFBYTES> class IntegerTypeSelection;
  //! partial specializations 
  template<>
  class IntegerTypeSelection<8>{
  public:
    typedef ui8 IntegerType;
  };

  template<>
  class IntegerTypeSelection<16>{
  public:
    typedef ui16 IntegerType;
  };

  template<>
  class IntegerTypeSelection<32>{
  public:
    typedef ui32 IntegerType;
  };
  
  template<>
  class IntegerTypeSelection<64>{
  public:
    typedef ui64 IntegerType;
  };

  //signed ints
  template<>
  class IntegerTypeSelection<-8>{
  public:
    typedef si8 IntegerType;
  };

  template<>
  class IntegerTypeSelection<-16>{
  public:
    typedef si16 IntegerType;
  };

  template<>
  class IntegerTypeSelection<-32>{
  public:
    typedef si32 IntegerType;
  };
  
  template<>
  class IntegerTypeSelection<-64>{
  public:
    typedef si64 IntegerType;
  };


  template<class OBJ, bool B> class ValueTypeSelector;

  template<class OBJ>
  class ValueTypeSelector<OBJ,true>{
  public:
    typedef typename OBJ::field_type value_type;
  };


  template<class OBJ>
  class ValueTypeSelector<OBJ,false>{
  public:
    typedef typename OBJ::value_type value_type;
  };



  /**
   * Compile time determination of matrix type via Template Meta Programming
   * depending on the char value of the forth argument.
   */
  template<class T, int M, int N, char C> class MatrixTypeSelector;

  //! partial specialisation(s) 
  template<class T, int M, int N>
  class MatrixTypeSelector<T,M,N,'f'>{         //dense fieldmatrix    
  public:
    typedef Dune::FieldMatrix<T,M,N> MatrixType;
  
    template<class V>
    static inline MatrixType& fill_1st_row(V& fzM, MatrixType& A){
      int sz = (int)std::distance(fzM.begin(),fzM.end());
      adonis_assert(sz == A.cols);

      for(int j = 0; j < A.cols; ++j)
	A[0][j] = fzM[j];
    
      return A;
    }
 
    template<class V>
    static inline MatrixType& fill_1st_column(V& fzM, MatrixType& A){
      int sz = (int)std::distance(fzM.begin(),fzM.end());
      adonis_assert(sz == A.rows);

      for(int i = 0; i < A.rows; ++i)
	A[i][0] = fzM[i];
    
      return A;
    }
  };

  template<class T, int M, int N>
  class MatrixTypeSelector<T,M,N,'d'>{         //dense fielddiagonal 
  public:
    typedef FieldDiagonal<T,M> MatrixType;
  
    //!fill diagonal 
    template<class V>
    static inline MatrixType& fill_1st_row(V& fzM, MatrixType& D){
      size_t sz = (int)std::distance(fzM.begin(),fzM.end());
      adonis_assert(sz == D.size());

      for(size_t j = 0; j < D.size(); ++j)
	D[j] = fzM[j];
    
      return D;
    }
 
    template<class V>
    static inline MatrixType& fill_1st_column(V& fzM, MatrixType& D){
      return fill_1st_row(fzM,D);   //actually  the same as fill_1st_column
    }
  };


} //end namespace

#endif

