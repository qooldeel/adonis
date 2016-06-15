#ifndef SPARSITY_PATTERN_HH
#define SPARSITY_PATTERN_HH

#include <iostream>
#include <vector>
#include <set>
#include <cmath>

#include "../common/globalfunctions.hh"
#include "../misc/useful.hh"
#include "../common/adonisassert.hh"

namespace Adonis{

  /** 
   * \brief sparsity pattern of a dense matrix (currently in vector form
   * stored row-wise)
   *
   * \tparam T <TT>bool</TT>:  brute-force style, i.e. store zeros and nonzeros
   *                           via false and true values (O(m*n) storage)
   *           <TT>std::set<std::size_t> </TT> more effective storage: store 
   *                           only nonzeros of cols entries j per row i 
  **/
  template<class T> class VectorSet;

  template<>
  class VectorSet<std::set<std::size_t> >{
  public:
    typedef std::set<std::size_t> Type;
    typedef VectorSet<std::set<std::size_t> > ThisType;

    static inline const int dim(int rows, int cols) {return rows;}

    template<class SVEC,class MTX>
    static inline const void scan_col(int i, int col, SVEC& svec, const MTX& mtx, int& countNonz){
      std::set<std::size_t> s;
      for(int j = 0; j < col; ++j){
	if(!is_zero(mtx[RowMajor::offset(i,j,col)])){
	  s.insert(j);
	  countNonz++;
	}
      }
      svec[i] = s;
    }
 
    //!sparsity pattern
    //!If transpose = true, the sparsity pattern of \f$ F'(x)^T\f$ is computed
    template<class SVEC, class ADFUN> //m = Rangedim, n = Domaindim (=#vars)
    static inline SVEC& sparsity(SVEC& s, ADFUN& f, bool transpose = false){
      std::size_t m = f.Range(),
	n = f.Domain();

      if(m >= n){ //fwd mode
	SVEC r_set(n);
	for(std::size_t j = 0; j < n; ++j)
	  r_set[j].insert(j);
	s = f.ForSparseJac(n,r_set,transpose);
      }
      else{ //rev mode
	SVEC s_s(m);
	for(std::size_t i = 0; i < m; ++i)
	  s_s[i].insert(i);
	s  = f.RevSparseJac(m,s_s,transpose);
      }

      return s;
    }

    template<class SVEC>
    static inline int num_of_nonzeros(const SVEC& s){
      int nnz = 0;
      for(std::size_t i = 0; i < s.size(); ++i){
	nnz += static_cast<int>(s[i].size());
      }
      return nnz;
    }

     template<class SVEC> 
     static inline int number_of_nonempty_elements(const SVEC& s){
       int nempt = 0;
       for(std::size_t i = 0; i < s.size(); ++i){
	 if(s[i].size() != 0)
	   nempt++;
       }
       return nempt;
     }


    //! this may be beneficial if, e.g., the Newton matrix \f$ G'(u) = I - hS'(u)\f$ is to be considered, i.e. on the diagonal, there are non-zero entries assumed 
//! this member must be called after fct <TT>sparsity</TT>
    template<class SVEC>
    static inline SVEC& diag_nz(int m, int n, SVEC& s){
#ifndef NDEBUG
      if(m != n)
	ADONIS_ERROR(DimensionError,"This operation is only allowed for n = m.");
#endif
      for(int i =0; i < n; ++i)
	s[i].insert(i);
      return s;
    }

    static inline void info() {
      std::cout << "================================================="<<std::endl;
      std::cout << "Use 'std::set<std::size_t>' as template argument." << std::endl;
      std::cout << "================================================="<<std::endl;
}
    
    template<class SVEC>
    static inline void print_2_file(const SVEC& s, const std::string& outfile, std::size_t col){  //third argument is dummy
      typedef typename SVEC::value_type SimpleType;
      typedef typename SimpleType::const_iterator SetIterType;
      SetIterType setIt;
      std::ofstream ofs(outfile.c_str(),std::ios_base::out);
      for(std::size_t i = 0; i < s.size(); ++i){
      	setIt = (s[i]).begin();
      	for(size_t j = 0; j < s.size(); ++j){
      	  if((j == (*setIt)) && (s[i].size() != 0)){
      	    ofs << 1 << " ";
      	    setIt++;
      	  }
      	  else{
      	    ofs << 0 << " ";
      	  }
      	}
      	ofs << std::endl;
      }
      ofs.close();
      ThisType::info();
    }

    template<class RT>
    static inline const void print(std::ostream& os, int i, int cols, const RT& v){
      typedef typename RT::value_type value_type;
      typedef typename value_type::const_iterator const_iterator;
      
      os << i << ".)  ";
      for(const_iterator cit = (v[i]).begin(); cit != (v[i]).end(); ++cit){
	os << *cit << "  ";
      }
      os << std::endl;
    }
  };



  template<>
  class VectorSet<bool>{
  public:
    typedef bool Type;
    typedef VectorSet<bool> ThisType;

    static inline const int dim(int rows, int cols) {return rows*cols;}
 
    template<class SVEC,class MTX>
    static inline const void scan_col(int i, int col, SVEC& svec, const MTX& mtx, int& countNonz){
      for(int j = 0; j < col; ++j){
	if(!is_zero(mtx[RowMajor::offset(i,j,col)])){
	  svec[RowMajor::offset(i,j,col)] = true;
	  countNonz++;
	}
	else 
	  svec[RowMajor::offset(i,j,col)] = false;
      }
    }

     //!sparsity pattern
    template<class SVEC, class ADFUN> //m = Rangedim, n = Domaindim (=#vars)
    static inline SVEC& sparsity(SVEC& s, ADFUN& f, bool transpose = false){
      std::size_t m = f.Range(),
	n = f.Domain();

      if(m >= n){ //fwd mode
	SVEC r_bool(n * n);
	std::size_t i, j;
	for(i = 0; i < n; ++i){	
	  for(j = 0; j < n; ++j)
	    r_bool[ i * n + j] = false;
	  r_bool[ i * n + i] = true;
	}
	s = f.ForSparseJac(n, r_bool,transpose);
      }
      else{ //rev mode
	SVEC s_b(m * m);
	std::size_t i,ell;
	for(i = 0; i < m; ++i){	
	  for(ell = 0; ell < m; ++ell)
	    s_b[i * m + ell] = false;
	  s_b[i * m + i] = true;
	}
	s  = f.RevSparseJac(m, s_b, transpose);
      }
      return s;
    }


    template<class SVEC>
    static inline int num_of_nonzeros(const SVEC& s){
      int nnz = 0;
      for(std::size_t i = 0; i < s.size(); ++i){
	if(s[i] == true)
	  nnz++;
      }
      return nnz;
    }

     template<class SVEC> 
     static inline int number_of_nonempty_elements(const SVEC& s){
       return (int)s.size(); //here just num of nonzeros
     }

    //! this may be beneficial if, e.g., the Newton matrix \f$ G'(u) = I - hS'(u)\f$ is to be considered, i.e. on the diagonal, there are non-zero entries assumed
    //! this member must be called after fct <TT>sparsity</TT>
    template<class SVEC>
    static inline SVEC& diag_nz(int m, int n, SVEC& s){
#ifndef NDEBUG
      if(m != n)
	ADONIS_ERROR(DimensionError,"This operation is only allowed for n = m.");
#endif
      for(int i = 0; i < n; ++i){	
	for(int j = 0; j < n; ++j){
	  s[i*m + i] = true;
	}
      }
      return s;
    }

    static inline void info() {
      std::cout << "================================"<<std::endl;
      std::cout << "Use 'bool' as template argument." << std::endl;
      std::cout << "================================"<<std::endl;
    }

    template<class SVEC>
    static inline void print_2_file(const SVEC& s, const std::string& outfile, std::size_t col){
      print_matrix_pattern_2_file(s,col,outfile);
      ThisType::info();
    }

    template<class VEC>
    static inline const void print(std::ostream& os, int i, int cols, const VEC& v){
      for(int j = 0; j < cols; ++j) 
	os  << v[RowMajor::offset(i,j,cols)] << "  ";
      os << std::endl;
    }
  };


  /**
   * \brief sparsity pattern of a dense \f$ m \times n \f$ matrix
   * \tparam PATTERNTYPE either a std::vector<std::set<std::size_t> > or 
   *                     bool (brute-force pattern storage)
   * 
   * Note: The std::set<std::size> version is much more economical, especially
   *       when large sparse matrices are considered
   */
  template<class PATTERNTYPE>
  class SparsityPattern{
  public:
    typedef PATTERNTYPE SimpleType;
    typedef SimpleType value_type;
    typedef std::vector<SimpleType> SVecType;
    typedef SVecType VType;
    
    

  private:
    typedef VectorSet<SimpleType> VSetType; //don't make accessible
    int m_, n_, countNonz_;
    SVecType p_;
  
  public:
    typedef typename VSetType::Type ElementType; //either bool or set<size_t>

    SparsityPattern(int rows = 0, int cols = 0):m_(rows),n_(cols),countNonz_(0),p_(VSetType::dim(rows,cols)){}

    void initialize(int rows, int cols){
      m_ = rows;
      n_ = cols;
      p_.resize(VSetType::dim(rows,cols));
    }


    //create from full matrix
    template<class FULLMTX> 
    SVecType& create(const FULLMTX& mtx){
      countNonz_ = 0; //reset counter
      for(int i = 0; i < m_; ++i){
      	VSetType::scan_col(i,n_,p_,mtx,countNonz_);
      }
      return p_;
    }

    //! Compute sparsity pattern of Jacobian without evaluating it. This is
    //! independent of the evaluation point.
    //! There is no dependence on any evaluation point. Just detect the 
    //! dependence on the independent variables. Check out
    //! \param f CppAD::ADFun object
    //! \param transpose if true, compute sparsity pattern of transposed Jac
    //! <a href="http://www.coin-or.org/CppAD/Doc/sparse_jacobian.cpp.htm"> Sparse Jacobian via CppAD, part 1 </a>
    //! <a href="http://www.coin-or.org/CppAD/Doc/cppad_sparse_jacobian.cpp.htm"> Sparse Jacobian via CppAD, part 2 </a>
    template<class ADFUN>
    SVecType& calc_sparsity(ADFUN& f,bool transpose = false){
      return VSetType::sparsity(p_,f,transpose);
    }

    const int number_of_nonzeros() const{
      return VSetType::num_of_nonzeros(p_);
    }
    //const int number_of_nonzeros() const {return countNonz_;}


    const int number_of_nonempty_elements() const{
      return VSetType::number_of_nonempty_elements(p_);
    }

    const ElementType& operator[](std::size_t i) const {
      adonis_assert(i < p_.size());
      return p_[i];
    }

    ElementType& operator[](std::size_t i) {
      adonis_assert(i < p_.size());
      return p_[i];
    }


    SVecType& diagonal_never_zero(){
      return VSetType::diag_nz(m_,n_,p_);
    }

    //can be used for an external sparsity pattern of the same type
    SVecType& diagonal_never_zero(SVecType& p_s) const{
      return VSetType::diag_nz(m_,n_,p_s);
    }

    const SVecType& get_p() const {return p_;}
    SVecType& get_p() {return p_;}
    
    const std::size_t size() const {return p_.size();} //container size

    const int rows() const {return m_;}
    const int cols() const {return n_;}
    const int size1() const {return m_;}
    const int size2() const {return n_;}
    bool is_initialized() const {return ( (n_ == 0 || m_ == 0) ? false : true);}

    
    //! very strict density measurement
    bool is_dense() const{return ( ((*this).number_of_nonzeros() == m_*n_)?true:false);}
    //! a more rigorous density measurement -- maybe replace log by sqrt
    bool is_strictly_dense() const{ return ( ((*this).number_of_nonzeros() <= log2(m_*n_))?false:true);}
    
    void resize(int m, int n, const SimpleType& sty = SimpleType()){
      m_ = m;
      n_ = n;
      countNonz_ = 0;
      p_.resize(VSetType::dim(m,n),sty);
    }

    SimpleType& operator[](int i) {
      adonis_assert(i < VSetType::dim(m_,n_));
      return p_[i];
    }

    const SimpleType& operator[](int i) const {
      adonis_assert(i < VSetType::dim(m_,n_));
      return p_[i];
    }

    const void print_sparsity_pattern_2_file(const std::string& outfile){
      VSetType::print_2_file(p_,outfile,n_);
    }


    friend std::ostream& operator<<(std::ostream& os, const SparsityPattern& sp){
      if(sp.number_of_nonzeros() != 0){
	for(int i = 0; i < sp.m_; ++i){
	  VSetType::print(os,i,sp.n_,sp.p_);
	}
      }
      else
	os << "There were NO non-zeros found, pal."<<std::endl;
      return os;
    }
    
  };

} //end namespace

#endif
