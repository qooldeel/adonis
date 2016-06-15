#ifndef COPY_FIELDCONTAINER_2_DYNAMIC_ARRAY_WHEN_USING_LAPACK_HH
#define COPY_FIELDCONTAINER_2_DYNAMIC_ARRAY_WHEN_USING_LAPACK_HH

#include <iostream>

#include "../dunestuff/fmatrix.hh"
#include "../common/typeadapter.hh"
#include "../common/smartassign.hh"

namespace Adonis{
/**
   * \brief Tiny class which copies a FieldMatrix to a dynamic array. This is needed for usage within Lapack since <I>without</I> copying to a dynamic array <B>segmentation faults</BB> may occur. 
   * User can specify if he/she wants to store the matrix row-by-row or  column-by-column (default. Recall that FORTRAN is column-based).
   * Use it in case of complications when using raw Dune::FieldMatrix within LAPACK.
   *\code
          
   * \endcode
   */
  template<class T, int M = 0, int N = 0,
	   class NUMT = typename TypeAdapter<T>::Type> //before: TypeTraits<T>
  class Copy2DynamicArray{
  public:
    typedef NUMT* PtrType;
    typedef PtrType iterator;
    typedef const PtrType const_iterator; 
    
    typedef Dune::FieldMatrix<T,M,N> FieldType;


    Copy2DynamicArray(const FieldType& A = FieldType(), char rcwise = 'c'):A_(A),rc_(rcwise),dim_(M*N),a_(0){
      try{
	
	a_ = new NUMT [dim_];

      }
      catch(std::bad_alloc){
	ADONIS_ERROR(DerivedError, "Not enough memory available on disk to launch 'Copy2DynamicArray(const Dune::FieldMatrix<T,M,N>&, char rcwise)'.");
      }
      
      size_t x = 0;

      switch(rcwise){
      case 'r':                           //store A row-wise (default)
	for(int i = 0; i < M; ++i){   
	  for(int j = 0; j < N; ++j){
	    smart_assign(a_[x++],A[i][j]);
	  }
	}
	break; 
      case 'c':                           //store A column-wise 
	
	for(int j = 0; j < N; ++j){   
	  for(int i = 0; i < M; ++i){
	     smart_assign(a_[x++],A[i][j]);
	  }
	}
	break;
      default:
	ADONIS_ERROR(DerivedError,"rcwise = "<< rcwise << " is illegal! \n   Possible values: 'r'(ow-wise) or 'c'(olumn-wise) storage.");
      }
    }

    ~Copy2DynamicArray(){delete[] a_;}  //proper clean up of dynamic vector
    
    
    //copy construction
    Copy2DynamicArray(const Copy2DynamicArray& C2DA):A_(C2DA.A_),rc_(C2DA.rc_),dim_(C2DA.dim_),a_(0){
      
      try{
	a_ = new NUMT [dim_];
      }catch(std::bad_alloc){
	ADONIS_ERROR(DerivedError, "Not enough memory available on disk to launch 'Copy2DynamicArray(const Copy2DynamicArray&)'.");
      }
      for(int i = 0; i < dim_; ++i)
	a_[i] = C2DA.a_[i];
      
    }
    
    //no copy assignent because of 'const FieldType&'
    

    size_t size() const {return static_cast<size_t>(dim_);}
    
    inline void clear(){
      delete[] a_;
      dim_ = 0;
      a_ = 0;    //pointer is null pointer now
    } 
    

    /**
     * \brief Restructure the array to hold solution later on
     * \param n columns of lhs coefficient matrix = rows of solution matrix)  
     * Depending on rc_ the columns and rows, respectively, are stored consecutively with zeros between them. 
     *
     * This is primarily needed for <B>over/underdetermined</B> linear systems.
     */
    inline void restructure(int n){
      
      int m = M,   //A_.rows,
	nrhs = N,   //A_.cols,
	ldb = (m > n) ? m : n,  //max of lhs coefficient matrix
	diff = ldb - m;
	
      
      //std::cout<<"m = "<<m<<std::endl<<"nrhs = "<<nrhs<<std::endl<<"n = "<<n<<std::endl<<"ldb = "<<ldb<<std::endl<<std::endl<<diff<<std::endl<<"rc_ = "<<rc_<<std::endl;

      
      
      (*this).clear();                     //call destructor; dim_ = 0         
      
      try{
	a_ = new NUMT [dim_ = ldb*nrhs];      //new dimension
      }
      catch(std::bad_alloc){
	ADONIS_ERROR(DerivedError, "Not enough memory available on disk 4 resizing within Copy2DynamivVector.");
      }

      for(int i = 0; i < dim_; ++i)
	a_[i] = NUMT();                       //fill with T()'s


      if(rc_ == 'c'){                      //column-wise (FORTRAN style)
	int k = 0,
	  idx = 0;
	  
	for(int j = 0; j < nrhs; ++j){
	  for(int i = 0; i < m; ++i){
	    if(k == m){
	      k = 0;  //reset
	      idx += diff;
	      //std::cout << "idx = "<<idx << std::endl;
	    }
	    
	     smart_assign(a_[idx++],A_[i][j]);
	    ++k;
	  }
	}
      }
     
      else if (rc_ == 'r'){              //row-wise filling
	int k = 0,
	  idx = 0;
	
	for(int i = 0; i < m; ++i){
	  for(int j = 0; j < nrhs; ++j){
	    if(k == nrhs){
	      k = 0;  //reset
	      idx += diff;
	      //std::cout << "idx = "<<idx << std::endl;
	    }
	    
	    smart_assign(a_[idx++],A_[i][j]);
	    ++k;
	  }
	}
      }
      else{
	ADONIS_ERROR(DerivedError,"rc_ = "<< rc_ << " is illegal because it hasn't been defined so far!");
      }

    }
    

    inline PtrType get_ptr(){ return &a_[0];}  //give back address

    iterator begin() {return &a_[0];}
    iterator end() {return &a_[dim_];}
    const_iterator begin() const{return &a_[0];}
    const_iterator end() const{return &a_[dim_];}
    


    friend std::ostream& operator<<(std::ostream& os, const Copy2DynamicArray& C2DV){
      for(iterator i = C2DV.begin(); i != C2DV.end(); ++i) 
	os << *i << " ";
      os << std::endl;
      return os;
    }

  private:
    const FieldType& A_;
    char rc_;
    int dim_;
    NUMT* a_; 
  };


}

#endif 
