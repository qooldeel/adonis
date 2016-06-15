#ifndef SOME_SPARSE_STUFF_THAT_MIGHT_ALSO_BE_BENEFICIAL_HH
#define SOME_SPARSE_STUFF_THAT_MIGHT_ALSO_BE_BENEFICIAL_HH

#include <bitset>

#include "sparsematrix.hh"

namespace Adonis{

 /**
   * \brief Select between popular sparse matrix formats
   */
  template<char CHAR, class T, class INTEGER> 
  class SparseMatrixSelector{  //default setting, do nothing
  public:
    typedef std::size_t SizeType;

    template<class IT, class U>
    static inline void copy(U& u, IT valit){}

    template<class U>
    static inline void resize(U& u, SizeType dim){}

    template<class VALIT, class W, class MAPPER>
    static inline void reorder_values(VALIT valit, const W& w, const MAPPER& mapper){}

    template<class INT, class MAPPER>
    static inline INT val_index_order(INT k, const MAPPER& mapperIdx){return 0;}

    template<class IT>
    static inline IT row_iter(IT index, IT position){return 0;}
    template<class IT>
    static inline IT col_iter(IT index, IT position){return 0;}

    template<class IT, class MAPIT, class PATT>
    static inline void create_format_data(SizeType m, SizeType n, IT rowIt, IT colPtrIt, MAPIT mapperIt, const PATT& p){}

    template<class IT, class VALIT>
    static inline void fill_in_values(SizeType nnz, VALIT val, VALIT a, IT map){
    }  

     //! select proper sizes for the storage arrays
    static inline SizeType dim_row_array(SizeType m, SizeType n, SizeType nnz){
      return 0;
    }
    
    static inline SizeType dim_col_array(SizeType m, SizeType n, SizeType nnz){
      return 0;
    }
  };

  template<class T, class INTEGER>
  class SparseMatrixSelector<'C',T, INTEGER>{ // compressed column storage
  public:
    typedef std::size_t SizeType;
    typedef T value_type;
    typedef CompressedSparseStorage<T,'C',INTEGER> SparseMatrixType;
    static const char Value = 'C';

    template<class IT, class U>
    static inline void copy(U& u, IT valit){
      for(std::size_t k = 0; k < u.size(); ++k){
	smart_assign(u[k], *(valit++));
      }
    }

    template<class U>
    static inline void resize(U& u, SizeType dim){
      u.resize(dim);
    }

    template<class VALIT, class W, class MAPPER>
    static inline void reorder_values(VALIT valit, const W& w, const MAPPER& mapper){
      adonis_assert(w.size() == mapper.size());

      for(std::size_t k = 0; k < w.size(); ++k){
	smart_assign(*(valit++), w[mapper[k]]);
      }
    }


    template<class INT, class MAPPER>
    static inline INT val_index_order(INT k, const MAPPER& mapperIdx){
      return (INT)mapperIdx[k];
    }

    //!the functions row_iter and col_iter are primarily used in conjunction 
    //! with create_format_data
    template<class IT>
    static inline IT row_iter(IT index, IT position){
      return index;
    }

    template<class IT>
    static inline IT col_iter(IT index, IT position){
      return position;
    }

    //! this need only be invoked once!
    //! Fill rowindex (nnz), colPtrIt (n_+1) and mapperIt (nnz). After invocation, the last iterator contains the mapping from which the order of the sparse matrix values in compressed row style can be easily transferred to compressed column form later on. 
    //!CAUTION: Mind the proper dimensions of the row and col iterator when 
    //!         switching between 'C' and 'R'. There is no size check!
    //! Complexity: \f$O(m \cdot n^2)\f$ when matrix is nearly dense
    //!             \f$O(m \cdot \mathrm{nnz})\f$ when matrix is suffic. sparse.
    template<class IT, class MAPIT, class PATT>
    static inline void create_format_data(SizeType m, SizeType n, IT rowIt, IT colPtrIt, MAPIT mapperIt, const PATT& p){
  
      SizeType ct = 0,
	dist,
	entriesPerCol = 0,
	nnz = 0,
	numOfEmptyRows = 0,
	numOfEmptyCols = 0;

      bool jfound;
      jfound = 0;   

      typedef typename PATT::value_type SetType;
      typedef typename SetType::const_iterator SetIteratorType;
     
      bool foundColEntryj = false;
    
      for(SizeType j = 0; j < n; ++j){     //cols
	ct = 0;
	*colPtrIt = entriesPerCol; 
	foundColEntryj = false;
	jfound = 0; //reset
	for(SizeType i = 0; i < m; ++i){   //rows
	  for(SetIteratorType sit = p[i].begin(); sit != p[i].end(); ++sit){
                  
	    if((*sit == j)){
	      dist = (int)std::distance(p[i].begin(),sit);           
	      *rowIt = i;
	      rowIt++;
            
	      *mapperIt = ct + dist; 
	      mapperIt++;
	      entriesPerCol++;
	      foundColEntryj = true;
	      jfound |= 1; //logical OR: once 1 is set, it won't be 0 again
	    }
	    else{
	      jfound |= 0; //no entry in this row of sparsity pattern
	    }
                  
	  }
	  ct += p[i].size();
	  if(j==0){  //only compute once
	    nnz += p[i].size();
	     
	    if(p[i].size() == 0)
	      numOfEmptyRows++;  //empty row detected
	  }
	}
	if(foundColEntryj) //otherwise a column is completely empty
	  colPtrIt++;

     
	if(jfound==0){ //column j is not contained in pattern
	  numOfEmptyCols++;
	}
      
      
      }
      
      *(colPtrIt++)  = nnz; //number of nonzeros assigned

      if(numOfEmptyRows != 0)
        std::cout << "There were "<< numOfEmptyRows << " empty rows detected in sparsity pattern." << std::endl;
      if(numOfEmptyCols != 0)
        std::cout << "There were "<< numOfEmptyCols << " empty cols detected in sparsity pattern." << std::endl;
    }


    template<class IT, class VALIT>
    static inline void fill_in_values(SizeType nnz, VALIT val, VALIT a, IT map){
      for(SizeType k = 0; k < nnz; ++k){
	val[k] = a[map[k]];
      }
    }  

     //! select proper sizes for the storage arrays
    static inline SizeType dim_row_array(SizeType m, SizeType n, SizeType nnz){
      return nnz;
    }
    
    static inline SizeType dim_col_array(SizeType m, SizeType n, SizeType nnz){
      return n+1;
    }

  };

  
  template<class T, class INTEGER>
  class SparseMatrixSelector<'c',T,INTEGER>{ // compressed column storage
  public:
    typedef std::size_t SizeType;
    typedef T value_type;
    typedef CompressedSparseStorage<T,'C',INTEGER> SparseMatrixType;
    static const char Value = 'c';

    template<class IT, class U>
    static inline void copy(U& u, IT valit){
      SparseMatrixSelector<'C',T,INTEGER>::copy(u,valit);
    }

    template<class U>
    static inline void resize(U& u, SizeType dim){
      SparseMatrixSelector<'C',T,INTEGER>::resize(u,dim);
    }

    template<class VALIT, class W, class MAPPER>
    static inline void reorder_values(VALIT valit, const W& w, const MAPPER& mapper){
       SparseMatrixSelector<'C',T,INTEGER>::reorder_values(valit,w,mapper);
     }

    template<class INT, class MAPPER>
    static inline INT val_index_order(INT k, const MAPPER& mapperIdx){
      return SparseMatrixSelector<'C',T,INTEGER>::val_index_order(k,mapperIdx);
    }

    template<class IT>
    static inline IT row_iter(IT index, IT position){
      return SparseMatrixSelector<'C',T,INTEGER>::row_iter(index,position);
    }

     template<class IT>
    static inline IT col_iter(IT index, IT position){
      return SparseMatrixSelector<'C',T,INTEGER>::col_iter(index,position);
    }

    template<class IT, class MAPIT, class PATT>
    static inline void create_format_data(SizeType m, SizeType n, IT rowIt, IT colPtrIt, MAPIT mapperIt, const PATT& p){
      SparseMatrixSelector<'C',T,INTEGER>::create_format_data(m,n,rowIt,colPtrIt,mapperIt,p);
      }

    template<class IT, class VALIT>
    static inline void fill_in_values(SizeType nnz, VALIT val, VALIT a, IT map){
      SparseMatrixSelector<'C',T,INTEGER>::fill_in_values(nnz,val,a,map);
    }  

    static inline SizeType dim_row_array(SizeType m, SizeType n, SizeType nnz){
       return  SparseMatrixSelector<'C',T,INTEGER>::dim_row_array(m,n,nnz);
    }
    
    static inline SizeType dim_col_array(SizeType m, SizeType n, SizeType nnz){
      return  SparseMatrixSelector<'C',T,INTEGER>::dim_col_array(m,n,nnz);
    }
  };



  //! this format is the C/C++ "standard" format
  template<class T, class INTEGER>
  class SparseMatrixSelector<'R',T,INTEGER>{ // compressed row storage
  public:
    typedef std::size_t SizeType;
    typedef T value_type;
    typedef CompressedSparseStorage<T,'R',INTEGER> SparseMatrixType;
    static const char Value = 'R';

    template<class IT, class U>
    static inline void copy(U& u, IT valit){} //do nothing

    template<class U>
    static inline void resize(U& u, SizeType dim){} //do nothing
	
    template<class VALIT, class W, class MAPPER>
    static inline void reorder_values(VALIT valit, const W& w, const MAPPER& mapper){} //do nothing
    
    template<class INT, class MAPPER>
    static inline INT val_index_order(INT k, const MAPPER& mapperIdx){
      return k;  //just return index
    }
      
    template<class IT>
    static inline IT row_iter(IT index, IT position){
      return position;
    }

    template<class IT>
    static inline IT col_iter(IT index, IT position){
      return index;
    }

    //!invoked only once!
    //! Complexity: \f$O(\mathrm{nnz})\f$. We only need to go through the
    //!             sparsity pattern and that's it.
    template<class IT, class MAPIT, class PATT>
    static inline void create_format_data(SizeType m, SizeType n, IT rowIt, IT colIt, MAPIT mapperIt, const PATT& p){
      SizeType k = 0,
	sz = 0,
	numOfEmptyRows = 0;
      typedef typename PATT::value_type SetType;
      typedef typename SetType::const_iterator SetIteratorType;

      for(SizeType i = 0; i < p.size(); ++i){
	*rowIt = sz;
	for(SetIteratorType it = p[i].begin(); it != p[i].end(); ++it){
	  *mapperIt = k;  //identity
	  *colIt = *it;
	  colIt++;
	  mapperIt++;
	  k++;
	  
	}
	if(p[i].size() == 0)
	  numOfEmptyRows++;

	sz += p[i].size();
	if(p[i].size() != 0) //otherwise, index is copied
	  rowIt++;
      }
      *(rowIt++)=k;

      if(numOfEmptyRows != 0)
        std::cout << "There were "<< numOfEmptyRows << " empty rows detected in sparsity pattern." << std::endl;
    }


    template<class IT, class VALIT>
    static inline void fill_in_values(SizeType nnz, VALIT val, VALIT a, IT map){
      //do nothing; values already correcly ordered
    }  

    //! select proper sizes for the storage arrays
    static inline SizeType dim_row_array(SizeType m, SizeType n, SizeType nnz){
      return m+1;
    }
    
    static inline SizeType dim_col_array(SizeType m, SizeType n, SizeType nnz){
      return nnz;
    }

  };

 
  template<class T, class INTEGER>
  class SparseMatrixSelector<'r',T,INTEGER>{ // compressed row storage
  public:
    typedef std::size_t SizeType;
    typedef T value_type;
    typedef CompressedSparseStorage<T,'R',INTEGER> SparseMatrixType;
    static const char Value = 'r';

    template<class IT, class U>
    static inline void copy(U& u, IT valit){
      SparseMatrixSelector<'R',T,INTEGER>::copy(valit,u);
    }

    template<class U>
    static inline void resize(U& u, SizeType dim){
      SparseMatrixSelector<'R',T,INTEGER>::resize(u,dim);
    }

    template<class VALIT, class W, class MAPPER>
    static inline void reorder_values(VALIT valit, const W& w, const MAPPER& mapper){
      SparseMatrixSelector<'R',T,INTEGER>::reorder_values(valit,w,mapper);
    }
    
    template<class INT, class MAPPER>
    static inline INT val_index_order(INT k, const MAPPER& mapperIdx){
      return SparseMatrixSelector<'R',T,INTEGER>::val_index_order(k,mapperIdx);
    }

    template<class IT>
    static inline IT row_iter(IT index, IT position){
      return SparseMatrixSelector<'R',T,INTEGER>::row_iter(index,position);
    }

    template<class IT>
    static inline IT col_iter(IT index, IT position){
      return SparseMatrixSelector<'R',T,INTEGER>::col_iter(index,position);
    }

    //invoke only once!
    template<class IT, class MAPIT, class PATT>
    static inline void create_format_data(SizeType m, SizeType n, IT rowIt, IT colIt, MAPIT mapperIt, const PATT& p){
      SparseMatrixSelector<'R',T,INTEGER>::create_format_data(m,n,rowIt,colIt,mapperIt,p);
    }

    template<class IT, class VALIT>
    static inline void fill_in_values(SizeType nnz, VALIT val, VALIT a, IT map){
      SparseMatrixSelector<'R',T,INTEGER>::fill_in_values(nnz,val,a,map);
    }  
    
    static inline SizeType dim_row_array(SizeType m, SizeType n, SizeType nnz){
       return  SparseMatrixSelector<'R',T,INTEGER>::dim_row_array(m,n,nnz);
    }
    
    static inline SizeType dim_col_array(SizeType m, SizeType n, SizeType nnz){
      return  SparseMatrixSelector<'R',T,INTEGER>::dim_col_array(m,n,nnz);
    }

  };


  /**
   * \brief Matrix wrapper
   */
  template<char CHAR, class T> 
  class MatrixTypeDetector{ //general case
  public:
    template<class MTX>
    static inline T* value_ptr(MTX& mtx){return mtx.values();}
  };

  //! Dense cases
  template<class T>
  class MatrixTypeDetector<'d',T>{
  public:
    template<class MTX>
    static inline T* value_ptr(MTX& mtx){return &mtx[0];}
  };
 
  template<class T>
  class MatrixTypeDetector<'D',T>{
  public:
    template<class MTX>
    static inline T* value_ptr(MTX& mtx){return &mtx[0];}
  }; 

} //end namespace

#endif
