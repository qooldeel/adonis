#ifndef SPARSE_CPPAD_DRIVER_HH
#define SPARSE_CPPAD_DRIVER_HH

#if USE_CPPAD
#include <cppad/cppad.hpp> //include the CPPAD stuff
#endif

#include "../common/globalfunctions.hh"

namespace Adonis{

#if USE_CPPAD

  /**
   * \brief Driver for efficient sparse Jacobian generation via CppAD.
   *  Prevents excessive temporary dynamic allocation and storage in 
   *  dense matrix formats for Jacobians. 
   *
   *  Original code can be found on the web [1].
   *
   *  References:
   *
   *  [1]  <a href="http://www.coin-or.org/CppAD/Doc/doxydoc/html/sparse__jacobian_8hpp_source.html#l01064"> sparse_jacobian.hpp </a>
   *
   * CAUTION: Currently you have to add to CppAD's local\ad_fun.hpp a public
   * members size_tsparse_jac_for(x,p,jac,work){return SparseJacobianFor(x,p,jac,work)
   * and size_t sparse_jac_rev(x,p,jac,work){return SparseJacobianRev(x,p,jac,work)
   * to make a swift Jacobian computation accessible.
   */
  template<class T, class PATT, class V>
  class SparseCppADDriver{
  public:
    typedef CppAD::ADFun<T> ADFunType;
    typedef T value_type;
    typedef typename PATT::value_type SetType;//either bool or set<size_t>
    typedef typename CppAD::internal_sparsity<SetType>::pattern_type PatternType;
    typedef std::size_t SizeType;
    
    SparseCppADDriver(ADFunType& f, const PATT& p):f_(f),p_(p){
      CppAD::CheckSimpleVector<T, V>();

      SizeType m = f_.Range(),
	n = f_.Domain(),
	i, j, k;
      
  
      if( n <= m ){  
	CppAD::vector<size_t>& row(work_.user_row);
	CppAD::vector<size_t>& col(work_.user_col);
	CppAD::vector<size_t>& sort_col(work_.sort_col);
 
	// forward mode, columns are sorted
	adonis_assert( row.size() == 0 );
	adonis_assert( col.size() == 0 );
	adonis_assert( sort_col.size() == 0 );

	// need an internal copy of sparsity pattern
	bool transpose = true;
	CppAD::sparsity_user2internal(inPatt_, p_, m, n, transpose);
 
	k = 0;
	for(j = 0; j < n; j++)
	  {    inPatt_.begin(j);
	    i = inPatt_.next_element();
	    while( i != inPatt_.end() )
	      {    row.push_back(i);
		col.push_back(j);
		sort_col.push_back(k);
		k++;
		i = inPatt_.next_element();
	      }
	  } 
        K_ = k;
	J_.resize(K_);
	row.push_back(m);
	col.push_back(n);
	sort_col.push_back(K_);

	// // now we have folded this into the following case
	// SparseJacobianFor(x, inPatt_, J_, work_);

	// // now set the non-zero return values
	// for(k = 0; k < K_; k++)
	//   jac[row[sort_col[k]] * n + col[sort_col[k]]] = J_[k];
      }
      else{
	CppAD::vector<size_t>& row(work_.user_row);
	CppAD::vector<size_t>& col(work_.user_col);
	CppAD::vector<size_t>& sort_row(work_.sort_row);
 
	// reverse mode, rows are sorted
	adonis_assert( row.size() == 0 );
	adonis_assert( col.size() == 0 );
	adonis_assert( sort_row.size() == 0 );
 
	// need an internal copy of sparsity pattern
	bool transpose = false;
	CppAD::sparsity_user2internal(inPatt_, p_, m, n, transpose);
 
	k = 0;
	for(i = 0; i < m; i++)
	  {    inPatt_.begin(i);
	    j = inPatt_.next_element();
	    while( j != inPatt_.end() )
	      {    row.push_back(i);
		col.push_back(j);
		sort_row.push_back(k);
		k++;
		j = inPatt_.next_element();
	      }
	  } 
	K_ = k;
	J_.resize(K_);
	row.push_back(m);
	col.push_back(n);
	sort_row.push_back(K_);
 
	// // now we have folded this into the following case
	// SparseJacobianRev(x, inPatt_, J_, work_);
 
	// // now set the non-zero return values
	// for(k = 0; k < K; k++)
	//   jac[row[sort_row[k]] * n + col[sort_row[k]]] = J_[k];
      }

      //return jac
    }

   

    V& sparse_jacobian(const V& x){
      adonis_assert(x.size() == (SizeType)(f_.Domain()));
      (f_.Domain() <= f_.Range()) ? (f_.sparse_jacobian_for(x,inPatt_,J_,work_)) : (f_.sparse_jacobian_rev(x,inPatt_,J_,work_));
      return J_;	
    }
    
  private:
    ADFunType& f_;
    const PATT& p_;
    V J_;
    CppAD::sparse_jacobian_work work_;
    PatternType inPatt_;
    SizeType K_;
  };


#endif

} //end namespace 

#endif
