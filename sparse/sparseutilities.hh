#ifndef SPARSE_MATRIX_UTILITIES_HH
#define SPARSE_MATRIX_UTILITIES_HH

#include <iostream>
#include <vector>

#include "../common/globalfunctions.hh"
#include "../common/numerictypechecker.hh"
#include "../common/adonisassert.hh"
#include "../common/error.hh"
#include "../common/smartassign.hh"
#include "../ode/newton.hh"
#include "../expressiontemplates/exprvec.hh"
#include "sparsesettings.hh"
#include "umftypetraits.hh"

namespace Adonis{

  /*
   * \brief creates compressed column storage from compressed row storage
   * Note that when you use UMFPACK, the column compressed scheme is 
   * required. Otherwise this routine performs expensive copying into the 
   * desired sparsity pattern
   *
   * Complexity: \f$ nz + cols\cdot nz. \f$ 
   *
   * \param rows, cols dimension of a matrix \f$A\f$ 
   * \param nz number of non-zeros of \f$ A\f$
   * \param  ixit,pit,valit compressed sparse representation of \f$ A\f$ (given
   *         as <I> iterators </I> )
   * \param indexNew, pNew, a new containers for column compressed 
   *                  representation (given as <I> iterators </I> )
   */
  template<class INT, class ITI, class ITP, class ITVAL, class ITINEW, class ITPNEW, class ITVALNEW>
  inline void row2column(INT rows, INT cols, INT nz, ITI ixit, ITP pit, ITVAL valit,  ITINEW indexNew, ITPNEW pNew, ITVALNEW a){
    //! note that 'a' and 'indexNew' hold nz entries and 'pNew' cols+1


    //! test version
    typedef std::vector<double>  VType;
   
    VType incr(nz);  //! just a temporary

    //! complexity: nz
    INT count = 0;
    for(INT k = 0; k < rows; ++k){
      while(count < pit[k+1]){
	incr[count++] = k;
      }
    } 
    //print_all(incr);

    //! complexity cols x nz
    count = 0;
    INT occur = 0, ety = 0;
    for(INT j = 0; j < cols; ++j){
      occur = 0;
      for(INT k = 0; k < nz; ++k){
	if(ixit[k] == j){
	  a[count] = valit[k];
	  indexNew[count] = incr[k];
	  count++;
	  occur++;
	  pNew[j] = ety;
	}
      }
      ety += occur;
    }
    pNew[cols] = nz;  //! Remember: pNew must have cols+1 entries!
    
    // std::cout << "Values (column compressed):"<<std::endl; 
    // print_all(a,a+nz);
    // std::cout << "Index (column compressed):"<<std::endl; 
    // print_all(indexNew,indexNew+nz);
    // std::cout << "Position pointers (column compressed):"<<std::endl; 
    // print_all(pNew,pNew+cols+1);
  }



  /**
   * \brief Very similar to <TT>row2column</TT>, except that all occurrences of
   * 'rows' and 'cols' are interchanged, and pNew is assumed to be an iterator
   * of a linear container having 'rows+1' elements. 
   */
  template<class INT, class ITI, class ITP, class ITVAL, class ITINEW, class ITPNEW, class ITVALNEW>
  inline void column2row(INT rows, INT cols, INT nz, ITI ixit, ITP pit, ITVAL valit,  ITINEW indexNew, ITPNEW pNew, ITVALNEW a){
    //! note that 'a' and 'indexNew' hold nz entries and 'pNew' rows+1
    
    
    //! test version
    typedef std::vector<double>  VType;
   
    VType incr(nz);  //! just a temporary

    //! complexity: nz
    INT count = 0;
    for(INT k = 0; k < cols; ++k){
      while(count < pit[k+1]){
	incr[count++] = k;
      }
    } 
    //print_all(incr);

    //! complexity cols x nz
    count = 0;
    INT occur = 0, ety = 0;
    for(INT j = 0; j < rows; ++j){
      occur = 0;
      for(INT k = 0; k < nz; ++k){
	if(ixit[k] == j){
	  a[count] = valit[k];
	  indexNew[count] = incr[k];
	  count++;
	  occur++;
	  pNew[j] = ety;
	}
      }
      ety += occur;
    }
    pNew[rows] = nz;  //! Remember: pNew must have rows+1 entries!
    
    ////! check output
    // std::cout << "Values (column compressed):"<<std::endl; 
    // print_all(a,a+nz);
    // std::cout << "Index (column compressed):"<<std::endl; 
    // print_all(indexNew,indexNew+nz);
    // std::cout << "Position pointers (column compressed):"<<std::endl; 
    // print_all(pNew,pNew+rows+1);
   }


  /**
   * \brief diagonal matrix times sparse matrix; only requires O(nnz) operations
   * res = D1*a, 
   * i.e. a is on the right hand side of the *-operator
   */
  template<class V, class DIAG, class PATT>
  inline V& right_diagonal_sparse_matrix_multiplication(V& res, const V& a, const DIAG& D1, const PATT& patt){
    typedef typename PATT::value_type SetType;
    typedef typename SetType::const_iterator SetIterType;
    adonis_assert(res.size() == a.size());

    std::size_t k = 0;
    for(std::size_t i = 0; i < patt.size(); ++i){
      for(SetIterType it = patt[i].begin(); it != patt[i].end(); ++it){
	res[k] = a[k]*D1[i];
	k++;
      }
    }
    return res;
  }
    

 /**
   * \brief sparse matrix times diagonal matrix; only requires O(nnz) operations
   * res = a*D2, 
   * i.e. a is on the left hand side of the *-operator
   */
  template<class V, class DIAG, class PATT>
  inline V& left_diagonal_sparse_matrix_multiplication(V& res, const V& a, const DIAG& D2, const PATT& patt){
    typedef typename PATT::value_type SetType;
    typedef typename SetType::const_iterator SetIterType;
    adonis_assert(res.size() == a.size());

    std::size_t k = 0;
    for(std::size_t i = 0; i < patt.size(); ++i){
      for(SetIterType it = patt[i].begin(); it != patt[i].end(); ++it){
	res[k] = a[k]*D2[*it];
	k++;
      }
    }
    return res;
  }


  

  //! forward declaration
  template<bool VIAAD> class DerivativeCalculationMethodSelection;

  //! no we assume that the sparsity pattern is computed via AD
  template<>
  class DerivativeCalculationMethodSelection<true>{
  public:

    enum{
      Value = 1
    };
    
    template<class JACOBJ, class EVAL> 
    static inline typename JACOBJ::ReturnType& compute(JACOBJ& jac, const EVAL& y){
      return jac.sparse_jacobian(y);
    }
    
    template<class ADJAC, class FDJAC>
    static inline ADJAC& right_jacobian_object(ADJAC& o1, FDJAC& o2){
      return o1;
    }
  };
 

  //! ... and via FD
  template<>
  class DerivativeCalculationMethodSelection<false>{
  public:
    
    enum{
      Value = 0
    };

    template<class JACOBJ, class EVAL> 
    static inline typename JACOBJ::ReturnType& compute(JACOBJ& jac, const EVAL& y){
      return jac.fwd(y);  //CAUTION: central(y), neg. values can be produced
                          //due to f(x-eps), which is fatal for computing, e.g.
                          //temperature in the Navier-Stokes equations
    }             

    template<class ADJAC, class FDJAC>
    static inline FDJAC& right_jacobian_object(ADJAC& o1, FDJAC& o2){
      return o2;
    }
  };
 
  
} //end namespace

#endif
