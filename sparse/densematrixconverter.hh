#ifndef DENSE_MATRIX_CONVERTER_HH
#define DENSE_MATRIX_CONVERTER_HH

#include <vector>
#include "../common/adonis_assert.hh"

namespace Adonis{

  /**
   * Convert dense matrices into sparse format
   * Assumption: matrices are stored row-wise (C-style)
   */
  template<class T, char C> class DenseMatrixConverter;

  template<class T> 
  class DenseMatrixConverter<T,'c'>{ //!column compressed form
  public:
    typedef T value_type;
    typedef std::vector<T> VType;
    typedef std::vector<INT> IxVType;

    template<class SPARSE, class INT, class DENSE>
    static inline void fill(SPARSE& sparseMtx, INT rows, INT cols, const DENSE& denseMtx, bool createPatternAnew = true){
      adonis_assert(rows != 0 || cols != 0); //makes no sense to do nothing here
      //! o.k. clearly 
      if(static_cast<INT>(sparseMtx.rows()) != rows || static_cast<INT>(sparseMtx.cols()) != cols || createPatternAnew == true){
	

      }
    }

  };
  

} //end namespace

#endif
