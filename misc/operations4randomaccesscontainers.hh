#ifndef OPERATIONS_4_STL_COMPLIANT_VECTORS_HH
#define OPERATIONS_4_STL_COMPLIANT_VECTORS_HH

#include <iostream>
#include <cmath>
#include <algorithm>
#include <typeinfo>

#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"
#include "../common/numerictypechecker.hh"
#include "../common/smartassign.hh"
#include "../misc/useful.hh"
#include "../misc/misctmps.hh"
#include "../common/typeadapter.hh"
#include "../common/isclass.hh"
#include "../common/typeselector.hh" //for field_types
#include "../common/universalconstants.hh"
#include "../common/elementaryoperations.hh"

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
#include "../accuracy/floatingpointarithmetic.hh"
#include "../expressiontemplates/xsettings.hh"
#endif

namespace Adonis{
  
  template<class U, class V, class W>
  inline U& concatenate(U& u, const V& v, const W& w){
    std::size_t dim = v.size() + w.size();
    //adonis_assert(u.size() == dim);
    (u.size() != dim) ? u.resize(dim) : do_nothing(); //resize when necessary

    for(std::size_t i = 0; i < v.size(); ++i)
      u[i] = v[i];
    for(std::size_t i = v.size(); i < dim; ++i)
      smart_assign(u[i],w[i-v.size()]);

    return u;
  }

  //! expand u by v
  template<class U, class V>
  inline U& concatenate(U& u, const V& v){
    std::size_t dim = u.size() + v.size(),
      olddim = u.size();
    
    (u.size() != dim) ? u.resize(dim) : do_nothing(); //resize when necessary
    
    for(std::size_t i = olddim; i < dim; ++i)
      smart_assign(u[i],v[i-olddim]);

    return u;
  }

  /**
   * \brief elementwise function evaluation
   * USAGE:
   * \code
      std::vector<double> v(3);
      v[0] = 1; v[1] = 2; v[2] = 3;
      entrywise_function_evaluation<double,std::exp>(v);
      for(size_t i = 0; i < 3; ++i)
       std::cout << v[i] << " ";
      std::cout << std::endl;
   * \endcode
   */
  template<class T, T FCT(T), class VCT>
  inline VCT& entrywise_function_evaluation(VCT& v){
    typedef typename VCT::iterator iterator;
    for(iterator it = v.begin(); it != v.end(); ++it)
      *it = FCT(*it); 
    return v;
  }

  /**
   * \brief Create randomly computed random access container
   */
  template<class T, template<class S, class A = std::allocator<S> > class V>
  inline void random_container(V<T>& c){
    typedef typename V<T>::iterator iterator;
    for(iterator it = c.begin(); it != c.end(); ++it){
      *it = static_cast<T>(rand())/RAND_MAX;
    }
  }

  /**
   * \brief Depending on type (scalars as well as containers), these template meta programs decide whether the norm of a STL-compliant container is used or the norm of a scalar (real or complex). Note that for scalars, the corresponding norm is just the (real or complex) absolute value ;)
   */
  template<class X, class RT,  bool B> class l1Metric;
  template<class X, class RT,  bool B> class l2Metric;
  template<class X, class RT,  bool B> class linfMetric;

  template<class X, class RT>
  class l1Metric<X,RT,true>{  //is container
  public:
    static inline RT metric(const X& v){
      RT norm = RT(); 
      
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
      RT s = RT(),
	c = RT();
#endif

      for(typename X::const_iterator it = v.begin(); it != v.end(); ++it){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	AdditionType::add(Abs(*it),s,c);
#else
	norm += Abs(*it);
#endif
      }

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<RT,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
      norm = s; 
#endif

      return norm;
    }
  };

  template<class X, class RT>
  class l1Metric<X,RT,false>{  //no container
  public:
    static inline RT metric(const X& x){
      return Abs(x);
    }
  };

 
  template<class X, class RT>
  class l2Metric<X,RT,true>{  //is container
  public:
    static inline RT metric(const X& v){
      RT norm = RT(); 
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
      RT s = RT(),   //sum
	c = RT();    //rounding correction
#endif

      for(typename X::const_iterator it = v.begin(); it != v.end(); ++it){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	AdditionType::add(ntimes<2>(Abs(*it)),s,c);
#else
	norm += ntimes<2>(Abs(*it));
#endif     
      }

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<RT,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
      norm = s;   //now assign sum
#endif

      return sqrt(norm);
    }
  };

  template<class X, class RT>
  class l2Metric<X,RT,false>{  //no container
  public:
    static inline RT metric(const X& x){
      return Abs(x);
    }
  };


  template<class X, class RT>
  class linfMetric<X,RT,true>{  //is container
  public:
    static inline RT metric(const X& v){
      RT norm = RT();
      for(typename X::const_iterator it = v.begin(); it != v.end(); ++it){
	norm = std::max(norm,Abs(*it));
      }
      return norm;
    }
  };

  template<class X, class RT>
  class linfMetric<X,RT,false>{  //no container
  public:
    static inline RT metric(const X& x){
      return Abs(x);
    }
  };
  

  /**
   * \brief A less constrained usage of norms than further above
   * \tparam V random access container 
   * \tparam RT <I> real </I> return type. Even if complex objects are employed, the return type is always real!
   */
  template<class V, class RT>
  inline RT norm_1(const V& v){
    return l1Metric<V,RT,IsContainer<V>::Value>::metric(v);
  }
  
  template<class V, class RT>
  inline RT norm_2(const V& v){
   return l2Metric<V,RT,IsContainer<V>::Value>::metric(v);
  }
  
  template<class V, class RT>
  inline RT norm_inf(const V& v){
     return linfMetric<V,RT,IsContainer<V>::Value>::metric(v);
  }

  /**
   * \brief Select Norm in an easy way
   * \tparam C '1', '2' or 'i'('I'), whether you want to use the \f$ l_1, l_2 \f$ or \f$ l_{\infty}\f$ norm, respectively
   * \tparam RT <I> real </I> return type
   */
  template<char C, class RT = double> class Norm;

  //!partial specialisations
  template<class RT>
  class Norm<'1',RT>{
  public:
    typedef RT return_type;
    static const char Value = '1';
    
    template<class V>
    static inline RT norm(const V& v){ return norm_1<V,RT>(v);}
  };

  template<class RT>
  class Norm<'2',RT>{
  public:
    typedef RT return_type;
    static const char Value = '2';

    template<class V>
    static inline RT norm(const V& v){ return norm_2<V,RT>(v); }
  };

  
  template<class RT>
  class Norm<'i',RT>{
  public:
    typedef RT return_type;
    static const char Value = 'i';

    template<class V>
    static inline RT norm(const V& v){ return norm_inf<V,RT>(v); }
  };

  template<class RT>
  class Norm<'I',RT>{
  public:
    typedef RT return_type;
    static const char Value = 'I';

    template<class V>
    static inline RT norm(const V& v){ return norm_inf<V,RT>(v); }
  };


  //===========================================================================
  
  template<class X, class RT, bool B> class lpMetric;

  template<class X, class RT>
  class lpMetric<X,RT,true>{ //! is container
  public:
    static inline RT metric(const X& v, const RT& p){
      RT norm = RT();

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;
      RT s = RT(),
	c = RT();
#endif
      for(typename X::const_iterator it = v.begin(); it != v.end(); ++it){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	AdditionType::add(std::pow(Abs(*it),p),s,c);	
#else
	norm += std::pow(Abs(*it),p);
#endif
      }
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<RT,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
      norm = s;   //now assign sum
#endif   
      return std::pow(norm,1./p);
    }
  };

  //! the scalar case -- then the norm just coincides with the Abs-value
  template<class X, class RT>
  class lpMetric<X,RT,false>{
  public:
    static inline RT metric(const X& x, const RT& p){
       return Abs(x);
     }
  };
  
  template<class X, class RT>
  inline RT norm_p(const X& x, const RT& p){
    return lpMetric<X,RT,IsContainer<X>::Value>::metric(x,p);
  }


  /**
   * \brief the \f$ p\f$-norm, i.e. \f[ \left(\|x\|_p := \sum_{i = 1}^n |x_i|^p\right^{\frac{1}{p}}, \f] where \f$ p \geq 1 \f$ is a <I>real</I> number
   *
   * NOTE: in cases where you only use the 1, 2, or \f$\infty\f$-norm, I 
   *       suggest to use <TT> Norm<'1',·>, Norm<'2',·> </TT> and 
   *       <TT> Norm<'i',·> </TT> or <TT> Norm<'I',·> </TT>, respectively, 
   *       since these do not employ the 'std::pow' function.
   */
  template<class RT>
  class PNorm{
  public:
    typedef RT return_type;
    
    template<class X>
    static inline RT norm(const X& x, const RT& p){
      //! we don't have a norm for \f$ 0 < p < 1\f$
#ifndef NDEBUG
      if(p < 1)
	ADONIS_ERROR(DerivedError, "p = "<< p << "\n   We only have a norm for 1 <= p < inf."); 
#endif
      return norm_p<X,RT>(x,p);
    }
  };


  //! classical dot product (also called inner product)
  template<class V, class W>
  typename V::value_type dot(const V& v, const W& w){
    size_t dim = std::distance(v.begin(),v.end());
    adonis_assert(dim == static_cast<size_t>(std::distance(w.begin(),w.end())));
    adonis_assert(typeid(typename V::value_type) == typeid(typename W::value_type));
    typedef typename V::value_type value_type;

    value_type d = value_type();

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType;

    value_type s = value_type(),
      c = value_type();
#endif

    for(size_t k = 0; k < dim; ++k){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      d = AdditionType::add(v[k]*Conj(w[k]),s,c);
#else
      d += v[k]*Conj(w[k]);
#endif
    }
    //! if Babuska-Kahan summation has been chosen, update s and assign it to d
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    CorrectRoundingErrorAfterwards<value_type,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
    CorrectRoundingErrorAfterwards<value_type,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(d,s);  // in any case, addition has been performed
#endif
    
    return d;
  }

  


  /**
   * \brief Normalised vector \f$ v \leftarrow \frac{1}{\|v \|_2}\cdot v. \f$
   */
   template<class T, template<class D, class A = std::allocator<D> > 
	    class V>
   inline void normalize(V<T>& v){
     typedef typename TypeAdapter<T>::BaseType BaseType;
     BaseType nm = Norm<'2',BaseType>::norm(v); 
     //adonis_assert(nm > T());    //no division by zero
     BaseType secure_nm = (Abs(nm) <= UniversalConstants<BaseType>::aboutZero) ? 1 : nm; //if nm == 0 then devide by 1

     //normalisation of vector v
     for(size_t i = 0; i < v.size(); ++i){
       v[i] /= secure_nm;
     }
   }


  /**
   * \brief Calculate unit vector of a given random access container, i.e. \f$     \hat{u} = \frac{u}{\|u\|}.\f$
   */
  template<class V>
  inline void unit_vector(V& v){
    
    //Can be used with FieldContainers or STL-compliant stuff
    typedef typename TypeAdapter<typename ValueTypeSelector<V,IsADuneFieldContainer<V>::Value>::value_type>::BaseType BaseType;

    typedef Norm<'2',BaseType> NormType;

    BaseType nm = NormType::norm(v),
      secure_nm = (Abs(nm) <= UniversalConstants<BaseType>::aboutZero) ? 1 : nm; //if nm == 0 then devide by 1

    for(typename V::iterator it = v.begin(); it != v.end(); ++it)
      *it /= secure_nm;
 
  }
  
/**
   * \brief prints array v in matrix style
   * \tparam V random access container 
   * \param v matrix stored in vector
   * \param numOfCols number of columns of matrix v
   */
  template<class V>
  inline void print_matrix(const V& v, int numOfcols){
    adonis_assert((int)v.size()%numOfcols == 0);
 
    int idx = 0;
    for(int i = 0; i < (int)v.size(); ++i){
      std::cout << v[i] <<  "   ";
      if(idx++ == numOfcols-1){
	idx = 0;  //reset
	std::cout << std::endl;
      }
    }
  }



  /**
   *\brief transposes matrix stored in vector of type VType
   * \param  A matrix to be transposed
   * \param cols number of columns
   * The rows are not stored 
   */
  template<class VType,class INT>
  inline VType transpose(const VType& A, INT cols){
    INT vdim = static_cast<INT>(std::distance(A.begin(), A.end()));
    adonis_assert((vdim%cols == 0) && (vdim >= cols));

    VType Trans;
    Trans.reserve(vdim);
    

    INT rows = vdim/cols;

    for(INT j = 0; j < cols; ++j){
      for(INT i = 0; i < rows; ++i){
	Trans.push_back(Conj(A[j + i*cols]));
      }
    }

    return Trans;
  }


  //! NO complex conjugate considered here 
  template<class VType,class INT>
  inline VType transpose_array(const VType& A, INT cols){
    INT vdim = static_cast<INT>(std::distance(A.begin(), A.end()));
    adonis_assert((vdim%cols == 0) && (vdim >= cols));

    VType Trans;
    Trans.reserve(vdim);
    

    INT rows = vdim/cols;

    for(INT j = 0; j < cols; ++j){
      for(INT i = 0; i < rows; ++i){
	Trans.push_back(A[j + i*cols]);
      }
    }

    return Trans;
  }


  template<class VType,class INT>
  inline VType& transpose_array(VType& Anew, const VType& A, INT cols){
    INT vdim = static_cast<INT>(std::distance(A.begin(), A.end()));
    adonis_assert((vdim%cols == 0) && (vdim >= cols));
    adonis_assert(vdim == static_cast<INT>(std::distance(Anew.begin(), Anew.end())) );
    

    INT rows = vdim/cols;

    for(INT j = 0; j < cols; ++j){
      for(INT i = 0; i < rows; ++i){
	Anew[j*cols+i] = A[i*cols+j];
      }
    }

    return Anew;
  }

  /**
    \brief Consecutively joins two vectors \f$w\f$ and \f$u\f$ into a new 
    *      one \f$v\f$
   */
  template<class V, class W, class U>
  inline V& join(V& v, const W& w, const U& u){
    size_t dim1 = static_cast<size_t>(std::distance(w.begin(),w.end()));
    size_t dim2 = static_cast<size_t>(std::distance(u.begin(),u.end()));
    size_t dim =  dim1 + dim2; 
    //! only resize when necessary
    (v.size() != dim) ? v.resize(dim) : do_nothing();
    for(size_t i = 0; i < dim1; ++i)  //first w...
      smart_assign(v[i],w[i]);
    for(size_t i = 0; i < dim2; ++i)  //...then u
      smart_assign(v[w.size()+i],u[i]);
    return v;
  }
  

  /**
 * \brief Invert a random access container
 */
  template<class V>
  inline V& invert(V& inv, const V& x){
    
    typedef typename V::iterator IterType;
    typedef typename V::const_iterator ConstIterType;
    //typedef typename ValueTypeSelector<V,IsADuneFieldContainer<V>::Value>::value_type value_type;
    
    ConstIterType xit = x.begin();
    
    for(IterType it = inv.begin(); it != inv.end(); ++it){
      adonis_assert(!is_zero(*xit)); //avoid division by zero
    //if(is_zero(*xit)) *it = *xit;
      
      *it = 1./(*xit); 
      ++xit;
    }
    
    return inv;
  }


  //:-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/
  //FROM HERE ON I WAS TOO LAZY TO EXTEND THE STUFF TO USING COMENSATED ADDITION
  //:-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/ :-/

  /**
   * \brief The one_norm for matrices which are stored row by row, i.e.
   \f[  A = \begin{pmatrix} 1 & 2 & 3 & 4 \\ 5 & 6 & 7 & 8 \\ 9 & 10 & 11 & 12 \end{pmatrix} \f]
   is stored in an array a[] = {1,2,3,4,5,6,7,8,9,10,11,12};  
   */
  template<class V>
  inline typename TypeAdapter<typename V::value_type>::BaseType one_norm(const V& v, int cols){
    
    typedef typename TypeAdapter<typename V::value_type>::BaseType BaseType;

    adonis_assert(v.size() > 0);

    if(v.size() == 1)
      return Abs(v[0]);
    
    BaseType sum = BaseType(),
      s_succ = BaseType(),
      res = BaseType();

  
    for(int i = 0; i < cols;++i){
      for(int j = i; j < (int)v.size(); j+=(cols)){
	sum += Abs(v[j]);
      }
      s_succ = sum;
      
      //std::cout<<" s_succ = "<< s_succ << std::endl; //control

      res = std::max(res, s_succ);
   
      sum = BaseType(); //reset 
    }
  
    return res;
  }

  


/**
 * \brief Product of dense matrix stored as vector with a vector
 */
  template<class V, class M, class INT>
inline V matrix_vector_product(const M& mtx, INT rows, const V& v){
  INT cols = static_cast<INT>(mtx.size())/rows;
  adonis_assert(cols == static_cast<INT>(v.size()));

  V u(rows);

  for(INT i = 0; i < rows; ++i)
    for(INT j = 0; j < cols; ++j)
      u[i] += mtx[RowMajor::offset(i,j,cols)]*v[j]; 

  return u;
}


  /**
   * \brief the dimensions of the multiplication can be recovered directly from the dimensions of m and v, respectively.
   */  
  template<class V, class M>
  inline V matrix_vector_product(const M& mtx, const V& v){
    adonis_assert(mtx.size()%v.size() == 0 && v.size() > 0 && mtx.size() > 0);
    
    size_t rows = mtx.size()/v.size();
    //std::cout << "matrix_vector_product: rows = "<< rows << std::endl;
    
    V u(rows);
    
    for(size_t i = 0; i < rows; ++i)
      for(size_t j = 0; j < v.size(); ++j)
	u[i] += mtx[RowMajor::offset(i,j,v.size())]*v[j]; 
    
    return u;
  }

  /**
   * \brief u = mtx·v
   */
   template<class V, class M>
   inline V& matrix_vector_multiplication(V& u, const M& mtx, const V& v){
     adonis_assert(mtx.size()%v.size() == 0 && v.size() > 0 && mtx.size() > 0);
    
    size_t rows = mtx.size()/v.size();
    
    (u.size() != rows) ? u.resize(rows) : do_nothing();
  
    
    for(size_t i = 0; i < rows; ++i)
      for(size_t j = 0; j < v.size(); ++j)
	u[i] += mtx[RowMajor::offset(i,j,v.size())]*v[j]; 
    
    return u;
  }



  template<class V, class M>
  inline V vector_matrix_product(const V& v, const M& m){
    adonis_assert(m.size()%v.size() == 0 && v.size() > 0 && m.size() > 0);

    size_t cols =  m.size()/v.size();

    V u(cols);

    for(size_t i = 0; i < v.size(); ++i)
      for(size_t j = 0; j < cols; ++j)
	u[j] += v[i]*m[RowMajor::offset(i,j,cols)];

    return u;

  }

 /**
   * \brief DUNE-matrix-vector-multiplication -- using a temp. object
   */
  template<class M, class V>
  inline V matrix_vector_multiplication(const M& m, const V& v){
    adonis_assert(m.cols == int(v.size()));
    
    V res(m.rows);
    
    for(int i = 0; i < m.rows; ++i)
      for(int j = 0; j < m.cols; ++j)
	res[i] += m[i][j]*v[j];
    
    return res;
  }

  /**
   * \brief Can be thought, e.g., of multiplying a diagonal matrix and a vector or vice versa (since diagonal matrix multiplication is commutative)
   * 
   * NOTE: a more elegant way can be followed by using expression templates here, e.g. ExprTmpl::MyVec
   */
  template<class V, class D>
  inline V vector_vector_product(const D& d, const V& w){
    adonis_assert(d.size() == w.size());
    V res(d.size());
    
    for(size_t i = 0; i < d.size(); ++i)
      res[i] = d[i]*w[i];
    
    return res;
  }

  /**
   * \brief create diagonal matrix (stored densely!!). This is only worth, when
   * some sort of Lapack routine is to be invoked
   */
  template<class V, class W>
  inline V& diagonal_matrix(V& diag, const W& w){
    typedef typename V::value_type value_type;
    
    std::size_t order = w.size();
    adonis_assert(diag.size() == ntimes<2>(order));

    for(std::size_t i = 0; i < order; ++i){
      for(std::size_t j = 0; j < order; ++j){
	if(i == j)
	  diag[RowMajor::offset(i,j,order)] = w[i];
	else
	  diag[RowMajor::offset(i,j,order)] = value_type();
      }
    }
    return diag;
  }

  //diag*A
  template<class V, class W, class U>
  inline V& left_diagonal_matrix_multiplication(V& res, const W& diag, const U& A){
    std::size_t m = diag.size(),
      n = A.size()/m;
    
    adonis_assert(A.size()%m == 0);
    
    (res.size() != A.size()) ? res.resize(A.size()) : do_nothing();

    for(std::size_t i = 0; i < m; ++i)
      for(std::size_t j = 0; j < n; ++j)
	res[RowMajor::offset(i,j,n)] = diag[i]*A[RowMajor::offset(i,j,n)];
    return res;
  }

  //A*diag
  template<class V, class W, class U>
  inline V& right_diagonal_matrix_multiplication(V& res, const W& diag, const U& A){
    std::size_t n = diag.size(),
      m = A.size()/n;
    
    adonis_assert(A.size()%n == 0);
    
    (res.size() != A.size()) ? res.resize(A.size()) : do_nothing();

    for(std::size_t i = 0; i < m; ++i)
      for(std::size_t j = 0; j < n; ++j)
	res[RowMajor::offset(i,j,n)] = diag[j]*A[RowMajor::offset(i,j,n)];
    return res;
  }



  template<class U, class V, class W>
  inline void vector_vector_addition(U& u, const V& v, const W& w){
    unsigned dim = container_size(v);
    adonis_assert(container_size(u) == dim && dim == container_size(w));
    
    for(unsigned i = 0; i < dim; ++i)
      u[i] = v[i] + w[i];
  }

  template<class U, class V, class W>
  inline void vector_vector_subtraction(U& u, const V& v, const W& w){
    unsigned dim = container_size(v);
    adonis_assert(container_size(u) == dim && dim == container_size(w));
    
    for(unsigned i = 0; i < dim; ++i)
      u[i] = v[i] - w[i];
  }

  //! this function may replace some of the preceding routines
  //! USAGE:
  //! \code
  //!   double a[] = {2, 0.5, -1.5, 0.75};
  //!   double b[] = {-0.65, 0.45, 2.24, -3.5};
  //!   std::vector<double> u(4), v(a,a+4), w(b,b+4);
  //!  //perform u = v + w
  //!  print_all(vector_operation<AddBasicElements>(u,v,w));
  //! \endcode
  template<template<class T> class OP, class U, class V, class W>
  inline U& vector_operation(U& u, const V& v, const W& w){
    std::size_t dim = container_size(v);
    adonis_assert(container_size(u) == dim && dim == container_size(w));
    
    for(std::size_t i = 0; i < dim; ++i)
      u[i] = OP<typename U::value_type>::apply(v[i],w[i]);
    return u;
  }



  /**
   * \brief see [GOLUB/VANLOAN, "Matrix Computations", 3rd ed., algo 1.2.3, p. 21] 
   * only lower triang part of A is stored (column-wise) resp. upper triag (row-wise)
   * \f$ y \leftarrow y + A\cdot x\f$
   */
  template<class V>
  inline V symm_matrix_vector_multiplication(const V& A, const V& x){
    adonis_assert(static_cast<size_t>(midnight_formula(1.,1.,-2.*A.size()).first) == x.size());
    
    V y(x.size());
    
    for(size_t j = 0; j < x.size(); ++j){
      for(size_t i = 0; i < j; ++i){
	y[i] += A[i*x.size() - (i+1)*i/2 + j]*x[j];      
      }
      for(size_t i = j; i < x.size(); ++i){
	y[i] += A[j*x.size() - (j+1)*j/2 + i]*x[j];
      }
    }
    return y;
  }

  //! empty y before by setting, e.g. y = 0
  template<class V>
  inline V& symm_matrix_vector_multiplication(V& y, const V& A, const V& x){
    adonis_assert(static_cast<size_t>(midnight_formula(1.,1.,-2.*A.size()).first) == x.size());

    (y.size() == 0) ? y.resize(x.size()) : do_nothing();
    
    for(size_t j = 0; j < x.size(); ++j){
      for(size_t i = 0; i < j; ++i){
	y[i] += A[i*x.size() - (i+1)*i/2 + j]*x[j];      
      }
      for(size_t i = j; i < x.size(); ++i){
	y[i] += A[j*x.size() - (j+1)*j/2 + i]*x[j];
      }
    }
    return y;
  }
  

  /**
   * identity of symmetric matrix (only store lower triangle)
   * NOTE: n is the order of I, not the size of the underlying 
   *       storage array!
   */
  template<class V, class INT>
  inline V& symm_identity(V& id, INT n){
    (id.size() == 0) ? id.resize(gauss_sum(n)) : do_nothing();
    for(INT i = 0; i < n; ++i)
      id[SymmetricAccess::offset(i,i,n)] = 1;
    return id;
  }

  /**
   * \brief This might be beneficial when, e.g. a dense symmetric linear system
   * is to be solved via LAPACK, which needs an \f$n \times n\f$ array as input
   */
  template<class SYM, class V>
  inline V& symmetric_matrix_2_full_matrix(V& v, const SYM& Avec){
    size_t n = static_cast<size_t>((midnight_formula(1.,1.,-2.*Avec.size())).first);
    //std::cout << "n = "<< n << std::endl;
    
    (v.size() != ntimes<2>(n)) ? v.resize(ntimes<2>(n)) : do_nothing();

    for(size_t i = 0; i < n; ++i){
      for(size_t j = 0; j < n; ++j){
	v[RowMajor::offset(i,j,n)] = Avec[SymmetricAccess::offset(i,j,n)];
      }
    }
    
    return v;
  }

			       
  template<class SYM, class V>
  inline V& full_matrix_2_symmetric_matrix(SYM& Avec, const V& v){
    size_t n = static_cast<size_t>(std::sqrt(v.size()));
    adonis_assert(v.size() == n*n);
    (Avec.size() == 0) ? Avec.resize(gauss_sum(n)) : do_nothing();
    for(size_t i = 0; i < n; ++i){
      for(size_t j = 0; j <= i; ++j){
	Avec[SymmetricAccess::offset(i,j,n)] = v[RowMajor::offset(i,j,n)];
      }
    }
    return Avec;
  }


  /**
   * \brief a more general method for adding, subtracting, multiplying and deviding  a vector, i.e.
  
   * \code
      MyVec<double> u(3), v(3), w(3);
      
      for(unsigned i = 0; i < 3; ++i){
        v[i] = 0.5*(i+1);
        w[i] = std::pow(0.5*(i+1),2.);
      }
      
   
      vector_vector_operation<AddBasicElements<double> >(u,v,w);

      std::cout << "u = " << u << std::endl;
   * \endcode 
   */
  template<class OP, class U, class V, class W>
  inline void vector_vector_operation(U& u, const V& v, const W& w){
    unsigned dim = container_size(v);
    adonis_assert(container_size(u) == dim && dim == container_size(w));
    
    for(unsigned i = 0; i < dim; ++i)
      u[i] = OP::apply(v[i],w[i]);
  }



  /***/
  template<class V, class W>
  inline void build_vector(V& v, const W& w, unsigned dim1, unsigned i){
    adonis_assert(v.size()%dim1 == 0);

    unsigned wdim = std::distance(w.begin(),w.end());

    for(unsigned j = 0; j < wdim; ++j){
      v[ColumnMajor::offset(i,j,dim1)] = w[j];
    }
  }

  template<class V>
  inline void absify(V& v){
    //typedef typename V::value_type value_type;
    for(size_t i = 0; i < v.size(); ++i){
      //(v[i] >= value_type(0.)) ? v[i] : -v[i]; //apply abs value
      v[i] = Abs(v[i]);
    }
  }


  /**
   * \brief computes \f$ \sqrt{\|v\|_A} := \sqrt{(v,Av)}\f$, where (.,.) 
   *  is the usual scalar product in \f$ R^n\f$ 
   * Note that A is a full matrix since A·v requires 2n² flops like when 
   * A is stored in some symmetric format (we neglect storage savings here) 
   */
  template<class M, class V>
  inline typename TypeAdapter<typename M::value_type>::BaseType scaled_norm(const V& v, const M& A){
    adonis_assert(std::distance(A.begin(),A.end())%std::distance(v.begin(),v.end()) == 0);
    //typedef typename TypeAdapter<typename M::value_type>::BaseType BaseType;
    return std::sqrt(dot(v,matrix_vector_product(A,v)));
  }

  //! version for symmetric matrices
  template<class V>
  inline typename TypeAdapter<typename V::value_type>::BaseType symm_scaled_norm(const V& v, const V& A){
    adonis_assert(static_cast<size_t>(midnight_formula(1.,1.,-2.*A.size()).first) == v.size());
    return std::sqrt(dot(v,symm_matrix_vector_multiplication(A,v)));
  }


  /**
   * \brief a norm wrapper to switch between scaled norm
   *  and usual one. In general, the scaled norm is defined via \f[ \| p\|_A := \sqrt{\langle p, Ap\rangle},\f] where \f$ A\f$ is a s.p.d. scaling  matrix   
   */
  template<char S, char NT> class NormWrapper;

  template<char NT>
  class NormWrapper<'s',NT>{ //!'s'caled norm; second argument is redundant here
  public:
    template<class V>
    static inline typename TypeAdapter<typename V::value_type>::BaseType norm(const V& v, const V& mtx, const typename TypeAdapter<typename V::value_type>::BaseType& p = 2){
      return symm_scaled_norm(v,mtx);
    }
  };

  template<char NT>
  class NormWrapper<'S',NT>{ //!second argument is redundant here
  public:
    template<class V>
    static inline typename TypeAdapter<typename V::value_type>::BaseType norm(const V& v, const V& mtx, const typename TypeAdapter<typename V::value_type>::BaseType& p = 2){
      return NormWrapper<'s',NT>::norm(v,mtx);
    }
  };

   template<char NT>
   class NormWrapper<'o',NT>{ //!'o'rdinary norm
   public:
     template<class V>
     static inline typename TypeAdapter<typename V::value_type>::BaseType norm(const V& v, const V& mtx = V(), const typename TypeAdapter<typename V::value_type>::BaseType& p = 2){   //!second fct argument is dispensable here
       typedef  typename TypeAdapter<typename V::value_type>::BaseType BaseType;
       return Norm<NT,BaseType>::norm(v);
     }
   };
  
   template<char NT>
   class NormWrapper<'O',NT>{ //!usual norm
   public:
     template<class V>
     static inline typename TypeAdapter<typename V::value_type>::BaseType norm(const V& v, const V& mtx = V(), const typename TypeAdapter<typename V::value_type>::BaseType& p = 2){   //!second fct argument is dispensable here
       //typedef  typename TypeAdapter<typename V::value_type>::BaseType BaseType;
       return NormWrapper<'o',NT>::norm(v,mtx);
     }
   };
  
 
  //! the p-norm
  template<char NT>
  class NormWrapper<'p',NT>{
  public:
    template<class V>
    static inline typename TypeAdapter<typename V::value_type>::BaseType norm(const V& v, const V& mtx, const typename TypeAdapter<typename V::value_type>::BaseType& p){
      typedef  typename TypeAdapter<typename V::value_type>::BaseType BaseType;
      return PNorm<BaseType>::norm(v,p);
    }
  };

  template<char NT>
  class NormWrapper<'P',NT>{
  public:
    template<class V>
    static inline typename TypeAdapter<typename V::value_type>::BaseType norm(const V& v, const V& mtx, const typename TypeAdapter<typename V::value_type>::BaseType& p){
      return NormWrapper<'p',NT>::norm(v,mtx,p);
    }
  };

  //! the dyadic or outer product
  template<class V>
  inline V dyad(const V& v, const V& w){
    V dy(v.size()*w.size());
    for(size_t i = 0; i < v.size(); ++i){
      for(size_t j = 0; j < w.size(); ++j){
	dy[RowMajor::offset(i,j,w.size())] = v[i]*w[j];//row major by default
      }
    }
    return dy;
  }

  template<class V>
  inline V& dyad(V& dy, const V& v, const V& w){
    (dy.size() != v.size()*w.size()) ? dy.resize(v.size()*w.size()) : do_nothing();
    for(size_t i = 0; i < v.size(); ++i){
      for(size_t j = 0; j < w.size(); ++j){
	dy[RowMajor::offset(i,j,w.size())] = v[i]*w[j];//row major by default
      }
    }
    return dy;
  }

  //!symmetric dyad -- store only lower triangular part (column-wise) resp.
  //!upper triangular part row-wise
  template<class V>
   inline V dyad(const V& v){
    V dy(gauss_sum(v.size()));
    for(size_t i = 0; i < v.size(); ++i){
      for(size_t j = 0; j <= i; ++j){
	//dy[j*v.size() - (j+1)*j/2 + i] 
	dy[SymmetricAccess::offset(i,j,v.size())]= v[i]*v[j];//row major by default
      }
    }
    return dy;
  }


  template<class V, class INT>
  inline V& matrix_matrix_multiplication(V& C, const V& A, INT arows, const V& B, INT bcols){
    (static_cast<INT>(C.size()) != arows*bcols) ? C.resize(arows*bcols) : do_nothing();

    adonis_assert(static_cast<INT>(A.size())%arows == 0 &&  static_cast<INT>(B.size())%bcols == 0 && arows != 0 && bcols != 0);
    INT n = A.size()/arows;
    adonis_assert(static_cast<INT>(B.size()/n) == bcols);
    
    for(INT i = 0; i < arows; ++i){
      for(INT j = 0; j < bcols; ++j){
	for(INT k = 0; k < n; ++k){
	  C[RowMajor::offset(i,j,bcols)] += A[RowMajor::offset(i,k,n)]*B[RowMajor::offset(k,j,bcols)];
	}
      }
    }
    return C;
  }

  
  template<class V>
  inline bool nearly_zero(const V& v){
    bool b = true;
    for(typename V::const_iterator cit = v.begin(); cit != v.end(); ++cit){
      if(!is_zero(*cit)){
	b = false;
	break;
      }
    }
    return b;
  }

  /**
   * \brief angle \f$\theta\f$ between two vectors \f$u\f$ and \f$v\f$
   *
   * \return V::value_type value type in radians \f$[0,\pi]\f$
   *
   * NOTE: std::acos is only defined for various real (<B>not</B> for complex
   *       numbers!
   */
  template<class V>
  inline typename V::value_type angle(const V& u, const V& v){
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
     BaseType nm = Norm<'2',BaseType>::norm(u);
     BaseType mn = Norm<'2',BaseType>::norm(v);
     BaseType secure_nm = (Abs(nm) <= UniversalConstants<BaseType>::aboutZero) ? 1 : nm; 
     BaseType secure_mn = (Abs(mn) <= UniversalConstants<BaseType>::aboutZero) ? 1 : mn; 
    
     return acos(dot(u,v)/(secure_nm*secure_mn));
  }

  /**
   * \brief Scalar projection: calculate scalar projection of one vector 
   *  \f$v\f$ in the direction of another vector \$fu\f$f
   */
  template<class V>  
  inline typename V::value_type scalar_projection(const V& v, const V& u){
    adonis_assert(v.size() == u.size());
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
     BaseType nm = Norm<'2',BaseType>::norm(u); 
     BaseType secure_nm = (Abs(nm) <= UniversalConstants<BaseType>::aboutZero) ? 1 : nm; //if nm == 0 then devide by 1
     return (dot(v,u)/secure_nm);
  }

  /**
   * \brief vector projection of \f$v\f$ in the direction of \f$u\f$
   */
  template<class V>  
  inline V vector_projection(const V& v, const V& u){
    adonis_assert(v.size() == u.size());
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    V res(v.size());
    BaseType nm = Norm<'2',BaseType>::norm(u); 
    BaseType secure_nm = (Abs(nm) <= UniversalConstants<BaseType>::aboutZero) ? 1 : nm; //if nm == 0 then devide by 1
    res = dot(v,u)/(ntimes<2>(secure_nm))*u;
    return res;
  }

  /**
   * \brief If \f$v \in R^{n}\f$ then \f$P = vv^T/vTv\f$ is the orthogonal 
   * projection onto \f$ span\{v\}, \f$ cf. [1, p. 75, § 2.6.1]
   *
   * References: 
   *
   *  [1]   [GOLUB and VAN LOAN, "Matrix Computations", 3rd ed., The John Hopkins University Press, 1996] 
   *
   * \tparam V storage type (linear rac, STL-compliant w.r.t. [] and size())
   * \param v vector whose ortho. projection is to be computed
   * \return V dense <I> symmetric </I> matrix type (stores only upper/lower half row-/columnwise)
   */
  template<class V>  
  inline V orthogonal_projection_onto_span_of_vector(const V& v){
    typedef typename V::value_type value_type;
    value_type sp = dot(v,v);
    value_type sec_sp = (Abs(sp) <= UniversalConstants<value_type>::aboutZero) ? 1 : sp; //applies also for complex numbers via their operator=(const T&)
    V P = dyad(v)/sec_sp;
    return P;
  }

  
}//end of namespace

#endif
