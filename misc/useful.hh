#ifndef USEFUL_TOOLS_CAN_NEVER_HARM_HH
#define USEFUL_TOOLS_CAN_NEVER_HARM_HH


#include "../common/adonisassert.hh"
#include "../dunestuff/fmatrix.hh"

namespace Adonis{
  

  /**
   *\brief wraps random access container that hasn't got a []-operator.
   *
   * const version
   */
  template<class BV>
  class ConstWrapVector2HaveBracketOperator{ 
  public:
    typedef BV ContainerType;
    typedef typename BV::value_type value_type;

     ConstWrapVector2HaveBracketOperator(const BV& v):v_(v){}

    template<class INT>
    const value_type& operator[](INT i) const {
      return v_(static_cast<unsigned>(i));
    }
    
    const BV& get_container() const {return v_;}

  private:
    const BV& v_;
    
  };

  //!mutable version
  template<class BV>
  class WrapVector2HaveBracketOperator{ 
  public:
    typedef BV ContainerType;
    typedef typename BV::value_type value_type;
    
    WrapVector2HaveBracketOperator(BV& v):v_(v){}

    template<class INT>
    value_type& operator[](INT i){
      return v_(static_cast<unsigned>(i));
    }
    
    BV& get_container() {return v_;}

  private:
    BV& v_;
    
  };
  

  /**
   * \brief Fill FieldMatrix from any linear STL-compliant random access container
   *\tparam SO Storage organisation (either <TT>RowMajor</TT> or <TT>ColumnMajor</TT> 
   */
  template<class SO>
  class FillMatrix{
  public:
    template<int M, int N, class V>
    static inline void now(Dune::FieldMatrix<typename V::value_type,M,N>& mtx, const V& v){
      adonis_assert((int)v.size() == M*N); 
      
      for(int i = 0; i < M; ++i){
	for(int j = 0; j < N; ++j){
	  mtx[i][j] = v[SO::offset(i,j,((SO::row_major) ? N : M))];
	}
      }
    }
  
    template<class T, int M, int N>
    static inline void now(Dune::FieldMatrix<T,M,N>& mtx, T* ptr){
     
      for(int i = 0; i < M; ++i){
	for(int j = 0; j < N; ++j){
	  mtx[i][j] = ptr[SO::offset(i,j,((SO::row_major) ? N : M))];
	}
      }
    }
    
    //! fill from given iterators
    template<class T, int M, int N, class ITER>
    static inline void now(Dune::FieldMatrix<T,M,N>& mtx, ITER i1, ITER i2){
      adonis_assert((int)std::distance(i1,i2) == M*N);

      for(int i = 0; i < M; ++i){
	for(int j = 0; j < N; ++j){
	  mtx[i][j] = *(i1 + SO::offset(i,j,((SO::row_major) ? N : M)));
	}
      }
    }
    

  };




  /**
   * \brief Storage organization of matrix like containers. This and the next class can be used when <B>dense</B> matrices are stored in a <B>linear</B> random access container such as a vector. 
   *
   * Consider, e.g. the matrix \f{eqnarray*}{ 
   A &:=& \begin{bmatrix} 1 & 2 & 3 & 4\\ 5 & 6 & 7 & 8 \\ 9 & 10 & 11 & 12 \end{bmatrix}, 
   \f}
   which is thought to be stored contiguously in memory. Then there are 2 possibilities in doing so, namely by C-like <B>row major</B> arrangement, i.e. \f[ 1 2 3 4 5 6 7 8 9 10 11 12\f] or by F77-like <B>column major</B> arrangement, i.e. \f[1 5 9 2 6 10 3 7 11 4 8 12.\f]
   *
   * NOTE: Depending on the storage organization the indexing of array elements is different (with row major order, knowledge of the number of <B>columns</B>, with column order the number of <B>rows</B> is required).   
   *
   * static, and no constructor need to be invoked
   */
  class RowMajor{
  public:
    template<class IX, class JX, class DX>
    static inline IX offset(const IX& i, const JX& j, const DX& cols){  //offset index 
      return i*cols + j;
    }
    
    template<class IX, class JX>
    static inline const JX& proper_dim(const IX& rows, const JX& cols){
      return cols; //! returns the orientation dimension i 'offset'.Here: cols
    }
    
    enum{row_major = true};
  };

  /**
   * \brief a column major ordered array is a transposed row major ordered array
   */
  class ColumnMajor{
  public:
    template<class IX, class JX, class DX>
    static inline IX offset(const IX& i, const JX& j, const DX& rows){  //offset index 
      return i + j*rows;
    }
    
    template<class IX, class JX>
    static inline const IX& proper_dim(const IX& rows, const JX& cols){
      return rows; //! returns the orientation dimension i 'offset'.Here: rows
    }

    enum{row_major = false};
  };
  
  
  class SymmetricAccess{
  public:
    template<class IX, class DX>
    static inline IX offset(const IX& i, const IX& j, const DX& n){
      return ( (i>=j) ? j*n - (j+1)*j/2 + i : SymmetricAccess::offset(j,i,n) );
    }
  };

}

#endif //include guard
