#ifndef BAND_MATRIX_STORAGE_FORMAT_HH
#define BAND_MATRIX_STORAGE_FORMAT_HH

#include <cstdlib>
#include <cmath>
#include <iostream> 
#include <complex>
#include <algorithm>
#include <typeinfo>

#include "../common/error.hh"
#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../misc/useful.hh"
#include "../marcvecmatrix/MAVEC.h"

namespace Adonis{

  /**
   * \brief Band matrix storage (<b>square matrices</b> only
   *
   *  I use a storage format that I've slightly altered from [1]
   *
   * References:
   *
   *  [1]  [YANG, "C++ and Object-Oriented Numeric Computing", Springer, 2001,Chap. 11]
   */
  template<class T> 
  class BandMtx{                // banded matrix
  private:
    int dim_,                                      //matrix order
      leftbw_,                                    // left bandwidth 
      rightbw_;                                     // right bandwidth 
    T** bdmx_;                                      // entries within the band 
   
  public: 
    typedef T value_type;
    typedef typename TypeAdapter<T>::BaseType BaseType;

    //default constructor
    BandMtx():dim_(0),leftbw_(0),rightbw_(0){
      bdmx_ = new T* [dim_]; 
      for (int i = 0; i < dim_; ++i) {
	bdmx_[i] = new T [leftbw_ + rightbw_ + 1];
	bdmx_[i] += leftbw_;                              
      }
      for (int i = 0; i < dim_; ++i)  //connection between full and banded matrix
	for (int j = -leftbw_; j <= rightbw_; ++j) 
	  bdmx_[i][j] = T();
    }

    BandMtx(int n, int p, int r, T** t):dim_(n),leftbw_(p),rightbw_(r){          // n: # of rows = # of columns
    // p: left bandwidth, r: right bandwidth, t: entries within the band 
      bdmx_ = new T* [n]; 
      for (int i = 0; i < n; ++i) {
	bdmx_[i] = new T [leftbw_ + rightbw_ + 1];
	bdmx_[i] += leftbw_;                              
      }
      for (int i = 0; i < n; ++i)  //connection between full and banded matrix
	for (int j = -leftbw_; j <= rightbw_; ++j) 
	  bdmx_[i][j] = t[i][j];
    }

    BandMtx(int n, int p, int r, const T& t = T()):dim_(n),leftbw_(p),rightbw_(r){         // initialize all entries to t
      bdmx_ = new T* [n]; 
      for (int i = 0; i < n; ++i) {
	bdmx_[i] = new T [leftbw_ + rightbw_ + 1];
	bdmx_[i] += leftbw_;                              
      }
      for (int i = 0; i < n; ++i) 
	for (int j = -leftbw_; j <= rightbw_; ++j) 
	  bdmx_[i][j] = t;
    }

    BandMtx(int n, int p, int r, const matrix<T>& mtx):dim_(n),leftbw_(p),rightbw_(r){    // m: full matrix format
       bdmx_ = new T* [n]; 
       for (int i = 0; i < n; ++i) {
	 bdmx_[i] = new T [leftbw_ + rightbw_ + 1];
	 bdmx_[i] += leftbw_;                              
       }
       for (int i = 0; i<n; ++i) {
	 for (int j = -leftbw_; j <= rightbw_; ++j) 
	   bdmx_[i][j] = 0;
       }
       for (int i = 0; i<n; ++i) {                   //mtx may be non-symmetric
	 int ip = std::max(i-leftbw_, 0); 
	 int ir = std::min(i+rightbw_, n-1); 
	 for (int j = ip; j <= ir; ++j) 
	   bdmx_[i][j-i] = mtx[i][j];
       }
    }

    template<class VEC>  //from full matrix stored as rac
    BandMtx(int n, int p, int r, const VEC& mtx):dim_(n),leftbw_(p),rightbw_(r){
      if(typeid(typename VEC::value_type) != typeid(T))
	ADONIS_ERROR(TypeError, "Not the same types used!");
      bdmx_ = new T* [n]; 
       for (int i = 0; i < n; ++i) {
	 bdmx_[i] = new T [leftbw_ + rightbw_ + 1];
	 bdmx_[i] += leftbw_;                              
       }
       for (int i = 0; i<n; ++i) {
	 for (int j = -leftbw_; j <= rightbw_; ++j) 
	   bdmx_[i][j] = 0;
       }
       for (int i = 0; i<n; ++i) {                   //mtx may be non-symmetric
	 int ip = std::max(i-leftbw_, 0); 
	 int ir = std::min(i+rightbw_, n-1); 
	 for (int j = ip; j <= ir; ++j) 
	   bdmx_[i][j-i] = mtx[RowMajor::offset(i,j,n)];
       }
    }
     
    BandMtx(const BandMtx& bd):dim_(bd.dim_),leftbw_(bd.leftbw_),rightbw_(bd.rightbw_){                      // copy constructor
      bdmx_ = new T* [dim_]; 
      for (int i = 0; i < dim_; ++i) {
	bdmx_[i] = new T [leftbw_ + rightbw_ + 1];
	bdmx_[i] += leftbw_;                              
      }
      for (int i = 0; i < dim_; ++i) 
	for (int j = -leftbw_; j <= rightbw_; ++j) 
	  bdmx_[i][j] = bd[i][j];
    }
    

    BandMtx& operator=(const BandMtx& bd){            // overload =
      if(this != &bd){
	//delete actual matrix
	for (int i = 0; i < dim_; ++i) 
	  delete[] (bdmx_[i] -= leftbw_);
	delete[] bdmx_; 
	
	//assign new dimensions
	dim_ = bd.dim_;
	leftbw_ = bd.leftbw_;
	rightbw_ = bd.rightbw_;
	bdmx_ = new T* [dim_]; 
	 
	for (int i = 0; i < dim_; ++i) {
	  bdmx_[i] = new T [leftbw_ + rightbw_ + 1];
	  bdmx_[i] += leftbw_;                              
	}
	for (int i = 0; i < dim_; ++i) 
	  for (int j = -leftbw_; j <= rightbw_; ++j) 
	    bdmx_[i][j] = bd[i][j];
      }
      return *this;
    }
    
    ~BandMtx(){  //destructor
      for (int i = 0; i < dim_; ++i) 
	delete[] (bdmx_[i] -= leftbw_);
      delete[] bdmx_; 
    }
    
    T* operator[](int i) const { return bdmx_[i]; } // i-th row of bdmx_
   
    template<class VEC>
    VEC operator*(const VEC& vec) const{        //! matrix vector multiplication
      if (dim_ != static_cast<int>(vec.size()))
	ADONIS_ERROR(DimensionError, "Order of banded matrix and vector do not match.");
      if(typeid(typename VEC::value_type) != typeid(T))
	ADONIS_ERROR(TypeError, "Not the same types used!");

      VEC tm(dim_);
      for (int i = 0; i < dim_; ++i) {
	int ip = std::max(i-leftbw_, 0); 
	int ir = std::min(i+rightbw_, dim_- 1);
	for (int j = ip; j <= ir; ++j) 
	  tm[i] += bdmx_[i][j-i]*vec[j]; 
      }
      return tm;
    }

    //! just some members that return values of the fields
    const int left_band_width() const {return leftbw_;}
    const int right_band_width() const {return rightbw_;}
    const int rows() const {return dim_;}
    const int cols() const {return dim_;}
    const int dim() const {return dim_;}

    //! launches band matrix to screen in its peculiar storage format
    friend std::ostream& operator<<(std::ostream& os, const BandMtx& bd){
      os << std::endl << "Band matrix of order "<< bd.dim_<< std::endl;
      for(int i = 0; i < bd.dim_; ++i){  //rows
	for(int j = -bd.leftbw_; j <= bd.rightbw_; ++j){ //columns
	  os << bd[i][j] << "  ";
	}
	os << std::endl;
      }
      return os;
    }

    //! preconditioning for band matrix
    template<class VEC>                      // solve Pz = r, return z
    VEC preconditioning(const VEC & r, int precn) const {
      if(typeid(typename VEC::value_type) != typeid(T))
	ADONIS_ERROR(TypeError, "Not the same types used!");

      if (precn == 0) {                    // no preconditioning
	return r;
      } else if (precn == 1) {             // diagonal preconditioning
	VEC z(dim_);
	for (int i = 0; i < dim_; ++i) z[i] = r[i]/bdmx_[i][0];
	return z;
      } else if (precn == 2) {             // symmetric SOR preconditioning
	const T omega = 1;                 // SSOR parameter for preconditioning
	VEC z(dim_);
	for (int i = 0; i < dim_; ++i) {
	  T sum = T();
	  int ip = std::max(i-leftbw_, 0);
	  for (int j = ip; j < i; ++j)  sum += omega*bdmx_[i][j-i]*z[j];
	  z[i] = (r[i] - sum)/bdmx_[i][0];
	}
	for (int i = dim_-1; i >= 0; --i) {
	  T sum = T();
	  int ir = std::min(i+rightbw_, dim_-1); 
	  for (int j = i+1; j <= ir; ++j) sum += omega*bdmx_[i][j-i]*z[j];
	  z[i] -= sum/bdmx_[i][0];
	}
	return z;
      } else {
	ADONIS_ERROR(Error,"Specified preconditioner for BandMtx not implemented");
      }
    }

    //!banded Gauss elimination with partial pivoting\f[A\cdot = b\f] 
    //!where $A$ is a band matrix and $b$ is the right hand side vector 
    //!\tparam VEC an STL-compliant random access container
    //!\param bb rhs vector will be overwritten during solution process
    template<class VEC>
    VEC& GaussElimPP(VEC& bb) const{      
     if(typeid(typename VEC::value_type) != typeid(T))
	ADONIS_ERROR(TypeError, "Not the same types used!");

     if (dim_ != static_cast<int>(bb.size())) 
       ADONIS_ERROR(DimensionError,"matrix-vector sizes do not match");

     BandMtx<T> tx(dim_,leftbw_,
		   std::min(dim_-1, leftbw_+rightbw_));
     for (int i = 0; i < dim_; ++i) 
       for (int j = -leftbw_; j <= rightbw_; ++j) 
	 tx[i][j] = bdmx_[i][j];

      VEC pvt(dim_);   // store pivoting info

      // LU decomposition with column partial pivoting
      const int nrowsmone = tx.dim() - 1;
      for (int k = 0; k < nrowsmone; ++k) {
	int kbrow = std::min(nrowsmone - k, tx.left_band_width()); 
	int kbcol = std::min(nrowsmone - k, tx.right_band_width()); 

	// find the pivot in the k-th column
	int pc = k;
	BaseType aet = Abs(tx[k][0]);
	for (int i = 1; i <= kbrow; ++i) {
	  if (Abs(tx[k+i][-i]) > aet) {
	    aet = Abs(tx[k+i][-i]);
	    pc = k + i;
	  }
	}
	if (!aet) 
	  ADONIS_ERROR(ZeroDivision,"Pivot is zero in banded Gauss elimination.");
	pvt[k] = pc;

	// interchange pivot row and k-th row in U, not in L
	if (pc != k) {
	  for (int j = 0; j <= kbcol; ++j) 
	    std::swap(tx[pc][k+j-pc], tx[k][j]);
	}

	// now eliminate column entries 
	for (int i = 1; i <= kbrow; ++i) { 
	  int kpi = k + i;
	  if (!is_zero(tx[kpi][-i])) {
	    T dmul = tx[kpi][-i]/tx[k][0];
	    tx[kpi][-i] = dmul;
	    for (int j = 1; j <= kbcol; ++j) 
	      tx[kpi][j-i] -= dmul*tx[k][j];
	  } 
	}
      }
      pvt[nrowsmone] = nrowsmone;

      // Forward substitution LY = b
      for (int k = 0; k < nrowsmone; ++k) {
	int kbrow = std::min(nrowsmone - k, tx.left_band_width()); 
	int pvtk = pvt[k];
	T sb = bb[pvtk];
	if (k != pvtk) std::swap(bb[k], bb[pvtk]);
	for (int j = 1; j <= kbrow; ++j) 
	  bb[k+j] -= tx[k+j][-j]*sb;
      }

      // Backward substitution U x = y
      for (int k = nrowsmone; k>= 0; --k) {
	int kb = std::min(nrowsmone -k, tx.right_band_width()); 
	for (int j = 1; j <= kb; ++j) bb[k] -= tx[k][j]*bb[k+j];
	bb[k] /= tx[k][0];
      }
      
      return bb;
    }

    //! just an alias
    template<class VEC>
    VEC& solve(VEC& rhs) const{
      return (*this).GaussElimPP(rhs);
    }
};

} //end namespace

#endif
