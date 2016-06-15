#ifndef SPARSE_MATRIX_CONTAINER_HH
#define SPARSE_MATRIX_CONTAINER_HH

#include <iostream>
#include <vector>
#include "umftypetraits.hh"
#include "sparseutilities.hh"

#include "../common/globalfunctions.hh"
#include "../common/error.hh"
#include "../common/typeadapter.hh"
#include "../common/adonisassert.hh"
#include "../misc/useful.hh"

namespace Adonis{

  /**
   * \brief For rectangular sparse matrices, the position pointer has a 
   * different size.
   * 
   */
  template<char C> class WhatCompression;

  //! specialization
  template<>
  class WhatCompression<'r'>{ //!row compressed
  public:
    template<class INT>
    static inline INT dim(INT rows, INT cols){
      return rows;
    }
  };
  
  template<>
  class WhatCompression<'R'>{ //!row compressed--  just invoke 'r' version
  public:
    template<class INT>
    static inline INT dim(INT rows, INT cols){
      return WhatCompression<'r'>::dim(rows,cols);
    }
  };
  
  template<>
  class WhatCompression<'c'>{ //!column compressed
  public:
    template<class INT>
    static inline INT dim(INT rows, INT cols){
      return cols;
    }
  };
  
  template<>
  class WhatCompression<'C'>{ //!row compressed -- just invoke 'c' version
  public:
    template<class INT>
    static inline INT dim(INT rows, INT cols){
      return WhatCompression<'c'>::dim(rows,cols);
    }
  };

  


  /**
   * \brief Compressed Sparse Row or Column format. The latter is used by 
   * UMFPACK.
   *
   * References:
   *  [1]  [Y. SAAD, "Iterative Methods for Sparse Linear Systems", 2nd ed., SIAM, chap. 3.4, pp. 84]
   */
  template<class T, char RC, class INT = int>
  class CompressedSparseStorage{
  public:
    typedef INT IndexType;
    typedef T value_type;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef std::size_t SizeType;
    typedef value_type* PointerType;
    typedef IndexType* IndexPointerType;

    typedef std::vector<typename UmfpackTypeTraits<BaseType>::Type> VecUmfBaseType;

    static const char compressionType = RC;

    //alias
    static const char Value = RC;
 
    //! the char-int value is transformed
    const char value() const {return Value_;}
    
    //!default constructor
    CompressedSparseStorage():rows_(0),cols_(0),nz_(0), posdim_(0){
      index_ = new IndexType[nz_];
      //! the last position is the number of nonzeros, anyway, here 0 too
       //posdim_ = WhatCompression<RC>::dim(rows_,cols_)+1;
      position_ = new IndexType[posdim_];
      a_ = new T[nz_];

      for(IndexType i = 0; i < nz_; ++i){
	index_[i] = IndexType();
	a_[i] = T();
      }

    
      for(IndexType t = 0; t < posdim_;++t)
	position_[t] = IndexType();
    }


    //! constructor: pos can just contain the positions or the positions+1
    //!In any case, the last entry of<TT> position_ </TT>contains <TT>nz_</TT> 
    template<class IV1, class IV2, class V>
    CompressedSparseStorage(IndexType m, IndexType n, IndexType nz, const IV1& ix, const IV2& pos, const V& val):rows_(m),cols_(n),nz_(nz), posdim_(WhatCompression<RC>::dim(rows_,cols_)+1){
      SizeType dimix = static_cast<SizeType>(std::distance(ix.begin(),ix.end())),
	dimval = static_cast<SizeType>(std::distance(val.begin(),val.end()));

      adonis_assert(dimix == dimval && dimix == static_cast<SizeType>(nz));

      index_ = new IndexType[nz];
      
      IndexType pdimm1 = WhatCompression<RC>::dim(m,n);
      position_ = new IndexType[posdim_];
      a_ = new T[nz];

      for(IndexType i = 0; i < nz; ++i){
	index_[i] = ix[i];
	a_[i] = val[i];
      }

      for(IndexType t = 0; t < pdimm1; ++t)
	position_[t] = pos[t];
      position_[pdimm1] = nz;
    }

    //! Read from dense matrix that is stored row-wise in a linear container
    template<class VEC>
    void fill_from_row_wise_dense(INT rows, INT cols, const VEC& Afull){
      adonis_assert(static_cast<INT>(rows*cols) == static_cast<INT>(std::distance(Afull.begin(),Afull.end())));
      if(RC == 'c' || RC == 'C'){ //column compressed form
	//typedef T value_type;
	typedef std::vector<T> VType;
	typedef std::vector<INT> IxVType;
	
	INT nz = 0;
	VType av;
	IxVType ix, pos;
	
	av.reserve(rows*cols);
	ix.reserve(rows*cols);
	pos.reserve(cols+1);
	
	nz = 0;
	bool detected = false;

	pos.push_back(0);
	for(INT j = 0; j < cols; ++j){  //! just swap for-loops
	  detected = false;
	  for(INT i = 0; i < rows; ++i){
	    if(!is_zero(Afull[RowMajor::offset(i,j,cols)])){
	      ix.push_back(i);
	      av.push_back(Afull[RowMajor::offset(i,j,cols)]);
	      nz++;
	      detected = true;   //o.k. found some non-negative entry in column j
	    }
	  }
	  if(!detected){
	    ADONIS_ERROR(LinearAlgebraError,"Full matrix contains a (nearly) ZERO column.\n   This may cause severe problems!");
	    //ADONIS_INFO(Information,"Full matrix contains a (nearly) ZERO column");
	    //pos.push_back(pos.back()); //just duplicate last entry of pos
	  }
	  pos.push_back(nz);
	}
	//! o.k. now we've filled the vectors. now we assign them to the array
	delete[] index_;
	delete[] a_;
	delete[] position_; 

	index_ = new INT[nz];
	a_ = new T[nz];
	position_ = new INT[cols+1];

	//now assign values 
	for(INT i = 0; i < nz; ++i){
	  index_[i] = ix[i];
	  a_[i] = av[i];
	}
	for(INT t = 0; t < cols; ++t)
	  position_[t] = pos[t];
	position_[cols] = nz;
	//! assign 
	rows_= rows;
	cols_ = cols;
	nz_ = nz;
	posdim_ = cols+1;
      }
      if(RC == 'r' || RC == 'R'){ //column compressed form
	//typedef T value_type;
	typedef std::vector<T> VType;
	typedef std::vector<INT> IxVType;
	
	INT nz = 0;
	VType av;
	IxVType ix, pos;
	
	av.reserve(rows*cols);
	ix.reserve(rows*cols);
	pos.reserve(cols+1);
	
	nz = 0;
	bool detected = false;
	
	pos.push_back(0);
	for(INT i = 0; i < rows; ++i){  //DIFFERENCE to 'c','C'
	  detected = false;
	  for(INT j = 0; j < cols; ++j){
	    if(!is_zero(Afull[RowMajor::offset(i,j,cols)])){
	      ix.push_back(j);         //DIFFERENCE to 'c','C': store column
	      av.push_back(Afull[RowMajor::offset(i,j,cols)]);
	      nz++;
	      detected = true;   //o.k. found some non-negative entry in column j
	    }
	  }
	  if(!detected){
	    ADONIS_ERROR(LinearAlgebraError,"Full matrix contains a (nearly) ZERO column.\n   This may cause severe problems!");
	    //ADONIS_INFO(Information,"Full matrix contains a (nearly) ZERO column");
	    //pos.push_back(pos.back()); //just duplicate last entry of pos
	  }
	  pos.push_back(nz);
	}
	//! o.k. now we've filled the vectors. now we assign them to the array
	delete[] index_;
	delete[] a_;
	delete[] position_; 
	
	index_ = new INT[nz];
	a_ = new T[nz];
	position_ = new INT[rows+1];         //DIFFERENCE to 'c','C'

	for(INT i = 0; i < nz; ++i){
	  index_[i] = ix[i];
	  a_[i] = av[i];
	}
	for(INT t = 0; t < rows; ++t)   //DIFFERENCE to 'c','C'
	  position_[t] = pos[t];
	position_[rows] = nz;           //DIFFERENCE to 'c','C'
	//! assign 
	rows_= rows;
	cols_ = cols;
	nz_ = nz;
	posdim_ = rows+1;
      }
    }
    
  

    //! 2nd constructor -- construct from iterators
    template<class InputIter1, class InputIter2, class InputIter3>
    CompressedSparseStorage(IndexType m, IndexType n, IndexType nz, InputIter1 f1, InputIter1 l1, InputIter2 f2, InputIter2 l2,InputIter3 f3, InputIter3 l3):rows_(m),cols_(n),nz_(nz),posdim_(WhatCompression<RC>::dim(m,n)+1){
      SizeType dimix = static_cast<SizeType>(std::distance(f1,l1)),
	dimval = static_cast<SizeType>(std::distance(f3,l3));

      adonis_assert(dimix == dimval && dimix == static_cast<SizeType>(nz));

      index_ = new IndexType[nz];
      
      IndexType pdimm1 = posdim_-1;// WhatCompression<RC>::dim(m,n);
      position_ = new IndexType[posdim_];
      a_ = new T[nz];
      
      InputIter1 it1 = f1;
      InputIter3 it3 = f3;
      for(IndexType i = 0; i < nz; ++i){
	index_[i] = *it1;
	a_[i] = *it3;
	++it1;
	++it3;
      }

      InputIter2 it2 = f2;
      for(IndexType t = 0; t < pdimm1; ++t){
	position_[t] = *it2;
	++it2;
      }
      position_[pdimm1] = nz;
    }

    //! copy constructor
    CompressedSparseStorage(const CompressedSparseStorage& spmtx){
      rows_ = spmtx.rows_;
      cols_ = spmtx.cols_;
      nz_ = spmtx.nz_;
      posdim_ = spmtx.posdim_;

      index_ = new IndexType[nz_];
      position_ = new IndexType[posdim_];
      a_ = new T[nz_];

      for(IndexType i = 0; i< nz_; ++i){
	index_[i] = spmtx.index_[i];
	a_[i] = spmtx.a_[i];
      }
      //! o.k. just copy from spmtx
      for(IndexType t = 0; t < posdim_; ++t)
	position_[t] = spmtx.position_[t];
    }


    //! copy assignment
    CompressedSparseStorage& operator=(const CompressedSparseStorage& spmtx){
      if(this != &spmtx){  //no self-assignment
	//adonis_assert(rows_ == spmtx.rows_ && cols_ == spmtx.cols_);
	//adonis_assert(nz_ == spmtx.nz_);
	
	if(nz_ != spmtx.nz_){
	  //std::cout << "===== OPERATOR=: CREATE NEW ARRAYS..."<< std::endl;
	  delete[] index_;
	  delete[] position_;
	  delete[] a_;

	  index_ = new IndexType[spmtx.nz_];
	  position_ = new IndexType[spmtx.posdim_];
	  a_ = new T[spmtx.nz_];
	}

	rows_ = spmtx.rows_;
	cols_ = spmtx.cols_;
	nz_ = spmtx.nz_;
	posdim_ = spmtx.posdim_;
	
	for(IndexType i = 0; i< nz_; ++i){
	  index_[i] = spmtx.index_[i];
	  a_[i] = spmtx.a_[i];
	}
	
	//! o.k. just copy from spmtx
	for(IndexType t = 0; t < posdim_; ++t){
	  position_[t] = spmtx.position_[t];
	}
      }
      
      return *this;
    }

    //! destructor
    ~CompressedSparseStorage(){
      delete[] a_;
      delete[] position_;
      delete[] index_;
    }


    void resize(SizeType m, SizeType n, SizeType nnz){
      rows_ = m;
      cols_ = n;
      nz_ = nnz;
      posdim_ = WhatCompression<RC>::dim(rows_,cols_)+1;

      if((index_ != 0) && (position_ != 0) && (a_ != 0)){
            delete[] index_;
            delete[] position_;
            delete[] a_;
      }
      
      index_ = new IndexType[nz_];
      position_ = new IndexType[posdim_];
      a_ = new T[nz_];

      // initialize arrays  --- better it is to prevent spurious results ;)
      for(IndexType i = 0; i < nz_; ++i){
	index_[i] = 0;
	a_[i] = T(0);
      }
      for(IndexType i = 0; i < posdim_; ++i)
	position_[i] = 0;
    }

    //! CAUTION: no size check, it's up to you to provide properly sized
    //!          input data!
    template<class IT>
    void fill_index(IT dataIt){
      for(IndexType i = 0; i < nz_; ++i)
	index_[i] = *(dataIt++);
    }

    template<class IT>
    void fill_position(IT posIt){
      for(IndexType i = 0; i < posdim_; ++i)
	position_[i] = *(posIt++);
    }

    template<class VALIT>
    void fill_values(VALIT valIt){
      for(IndexType i = 0; i < nz_; ++i)
	a_[i] = *(valIt++);
    }

    //! return pointer to beginning of arrays for compressed storage
    PointerType values(){return &a_[0];}
    IndexPointerType index(){return &index_[0];}
    IndexPointerType position(){return &position_[0];}

    const PointerType values() const {return &a_[0];}
    const IndexPointerType index() const {return &index_[0];}
    const IndexPointerType position() const{return &position_[0];}

    //! number of nonzeros in matrix
    const IndexType nonz() const {return nz_;}
    const IndexType rows() const {return rows_;}
    const IndexType cols() const {return cols_;}
    //! alias 
    const IndexType size1() const {return rows_;}
    const IndexType size2() const {return cols_;}

    const value_type& operator[](SizeType i) const {
      adonis_assert(static_cast<INT>(i) < nz_);
      return a_[i];
    }

    value_type& operator[](SizeType i){
      adonis_assert(static_cast<INT>(i) < nz_);
      return a_[i];
    }

    const IndexType dimension_of_position_array() const{return posdim_;}

    void restructure(SizeType m, SizeType n, SizeType nz){
      delete[] index_;
      delete[] position_;
      delete[] a_;

      index_ = new IndexType [nz];
      IndexType posdim = WhatCompression<RC>::dim(m,n)+1;
      position_ = new IndexType [posdim];
      a_ = new T[nz];

      //fill with defaults for proper initialization
      for(SizeType i = 0; i < nz; ++i){
	index_[i] = 0;
	a_[i] = T();
      }
      SizeType pm1 = posdim-1;
      for(SizeType p = 0; p < pm1; ++p)
	position_[p] = 0;
      position_[pm1] = nz;
    }

    // void reserve(SizeType rows, SizeType cols){
    //   delete[] index_;
    //   delete[] position_;
    //   delete[] a_;

    // }


    friend std::ostream& operator<<(std::ostream& os, const CompressedSparseStorage& spmtx){
      os << std::endl<< "Sparse matrix: rows: "<< spmtx.rows_ << "  cols: "<< spmtx.cols_ << "  #(nonzeros): "<< spmtx.nz_ << std::endl; 
      os << "Storage organisation: Compressed Sparse " << ((RC == 'c' || RC == 'C')? "Column " : "Row ") << "Scheme:"<< std::endl;
     
      os << std::endl<< "nonzeros:"<<std::endl;
      for(IndexType i = 0; i < spmtx.nz_; ++i)
	os << spmtx.a_[i] << "  ";
      os << std::endl;
     
      os << "index: "<< std::endl;
      for(IndexType i = 0; i < spmtx.nz_; ++i)
	os << spmtx.index_[i] << "  ";
      os << std::endl;

      os << "position:"<<std::endl;
      IndexType dm = spmtx.posdim_; //WhatCompression<RC>::dim(spmtx.rows_,spmtx.cols_)+1;
      for(IndexType t = 0; t < dm; ++t)
	os << spmtx.position_[t] << "  ";
      os << std::endl;
    
      return os;
    }


    //! solve system \f$ A\cdot x = b, \f$ where \f$A\f$ is sparse 
    //! DOES THIS ONLY WORK FOR n x n MATRICES ????????
    template<class V> 
    V& solve(V& x, const V& b){
#if USE_UMFPACK
      if(static_cast<INT>(x.size() == 0))
	x.resize(cols_);

      //typedef typename UmfpackTypeTraits<INT>::Type IndexType;
      typedef typename UmfpackTypeTraits<T>::Type Type;
      typedef typename Type::BaseType BaseType;

      int status, sys;
      BaseType Control [UMFPACK_CONTROL], Info[UMFPACK_INFO];
#else
      ADONIS_ERROR(NotImplementedError,"Sorry, I was too lazy to come up with a substitute for UMFPACK.");
#endif
    }

  private:
    IndexType rows_, cols_, 
	    nz_,   //number of nonzeros 
	    posdim_; //size of position array
    IndexType* index_;      //the index of size nz
    IndexType* position_;   //position marking beginnings of rows/columns in a_ 
    T* a_;
   
    enum{Value_ = RC};

    //! this is needed for the solution of a system \f$ A\cdot x = b\f$,
    //! especially when using Umfpack
    VecUmfBaseType are_, aim_,
	    bre_,bim_,
	    xre_,xim_;

  };
    
} //end namespace
    
#endif
