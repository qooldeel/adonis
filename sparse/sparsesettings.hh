#ifndef SETTINGS_FOR_SPARSE_MATRICES_HH
#define SETTINGS_FOR_SPARSE_MATRICES_HH

namespace Adonis{
  
  //#define SHOW_BAD_DIAGONAL_ENTRIES

  class SparseSettings{
  public:
    //!'U','u' = use Umfpack for solution
    static const char SparseSolver = 'U';  

    //!'C','c' = column compressed, 'R','r' = row compressed
    static const char TypeOfCompression = 'C';  

  };

} //end namespace 

#endif
