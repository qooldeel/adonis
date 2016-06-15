#ifndef VERSATILE_DERIVATIVE_CALCULATOR_API_HH
#define VERSATILE_DERIVATIVE_CALCULATOR_API_HH

#if USE_CPPAD
#include <cppad/cppad.hpp> //include the CPPAD stuff
#endif

#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../common/fancymessages.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../sparse/sparsitypattern.hh"
#include "../sparse/sparsetmp.hh"

namespace Adonis{

  template<class T, template<typename D1> 
	    class SOURCE, template<typename D2, class A = std::allocator<D2> > 
	   class VEC, class OBJ>
  class AbstractDerivativeCalculatorGeneral{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type BaseType;
    typedef VEC<BaseType> VecBaseType;
    typedef std::size_t SizeType;
    typedef VEC<SizeType> VecIndexType;
    typedef SOURCE<BaseType> FunType;


#if USE_CPPAD
    //common functionalities
    typedef CppAD::AD<T> ADType;
    typedef VEC<ADType> VecADType;
    typedef CppAD::ADFun<BaseType> ADFunType;
#else
   
#endif   

    template<class W>
    void set_with_init(SizeType domdim, SizeType randim, const W& init){
#if USE_CPPAD
      SOURCE<ADType> chemistry(randim); //source term
      VecADType X(domdim),  //domain space vector
	Y(randim);                    //range space vector
      
      //! for some functors it might be beneficial to initialize X first, since
      //! otherwise [0,0,....,0] might be an illegal value to start with 
      if(init.size() != 0){
	for(size_t i = 0; i < domdim; ++i)
	  X[i] = init[i];
      }
      else{
	ADONIS_INFO(Information, "No initialization of X. This may cause errors  in the source term if \n   physical meaningful relations are violated through this.");
      }

      CppAD::Independent(X);            //indenpendent variables, recording
                                        //must be declared here!! 
      Y = chemistry(X);
      adseq_.Dependent(X,Y);  //store sequence and stop recording 
      adseq_.optimize();
#else
      
#endif
    }


  protected:
    mutable bool isInitialized_;
#if USE_CPPAD
    ADFunType adseq_;
#else
    
#endif

  private:

    //delegate to derived class
    OBJ& ref_2_derived(){
      return static_cast<OBJ&>(*this); //give back derived class reference
    }

    //use properties of derived implementation
    //invoke derived class's member has to be implemented in every derived class
    template<class X>
    VecBaseType evaluate_jacobian(const X& x){
      return ref_2_derived().evaluate_jacobian(x);
    }

  };

  /* ==================================================================== */
  
  //forward declaration
  template<bool SPARSE, class T, template<typename D1> 
	    class SOURCE, template<typename D2, class A = std::allocator<D2> > 
	   class VEC, char CHAR>
  class JacobianMatrixCalculator: public AbstractDerivativeCalculatorGeneral<T,SOURCE,VEC,JacobianMatrixCalculator<SPARSE,T,SOURCE,VEC,CHAR> >{}; 


  //specializations -- use SPARSE matrix stuff
  /**
   * \brief Driver for efficient sparse Jacobian generation. When 
   *  compiler flag <TT>USE_CPPAD</TT> is set, CppAD is used as AD tool:
   * 
   *  It prevents excessive temporary dynamic allocation and storage in 
   *  dense matrix formats for Jacobians. 
   *
   *  Original code can be found on the web [1].
   *
   *  References:
   *
   *  [1]  <a href="http://www.coin-or.org/CppAD/Doc/doxydoc/html/sparse__jacobian_8hpp_source.html#l01064"> sparse_jacobian.hpp </a>
   *
   * CAUTION: Currently you have to add to CppAD's local/ad_fun.hpp the public
   * members size_t sparse_jac_for(x,p,jac,work){return SparseJacobianFor(x,p,jac,work)
   * and size_t sparse_jac_rev(x,p,jac,work){return SparseJacobianRev(x,p,jac,work)
   * to make a swift Jacobian computation accessible via CppAd.
   *
   * NOTE: up to now an error is thrown, if you attempt to use sparse FDs, 
   *       simply because I was too lazy to implement it. Feel free to complete.
   */
  template<class T, template<typename D1> 
	    class SOURCE, template<typename D2, class A = std::allocator<D2> > 
	   class VEC, char CHAR>
  class JacobianMatrixCalculator<true,T,SOURCE,VEC,CHAR>: public AbstractDerivativeCalculatorGeneral<T,SOURCE,VEC,JacobianMatrixCalculator<true,T,SOURCE,VEC,CHAR> >{
  private:
    typedef SparseMatrixSelector<CHAR,T> MtxSelectorType;

  public:
    typedef AbstractDerivativeCalculatorGeneral<T,SOURCE,VEC,JacobianMatrixCalculator<true,T,SOURCE,VEC,CHAR> > BaseClassType;
    typedef typename BaseClassType::SizeType SizeType;
    typedef typename BaseClassType::VecBaseType VecBaseType;
    typedef typename BaseClassType::FunType FunType;
    typedef typename BaseClassType::VecIndexType VecIndexType;
    typedef typename BaseClassType::BaseType BaseType;
     
    typedef typename MtxSelectorType::SparseMatrixType MatrixType;

    typedef std::set<SizeType> SetType;
    typedef SparsityPattern<SetType> SparsePatternType;
    typedef typename SetType::const_iterator ConstSetIterType;

    static const char Value = CHAR;

    JacobianMatrixCalculator():K_(0)
{
#if USE_CPPAD
      FancyMessages().nice_output("\n AD Jacobian obtained via AD CppAD in use.",33);
#else

      ADONIS_ERROR(NotImplementedError,"Sorry pal, but the evaluation of sparse FD 1st order derivatives is not supported yet ;)");
#endif     
      BaseClassType::isInitialized_ = false;
    }

    template<class W>
    void initialize(SizeType ddim, SizeType rdim, const W& init){
      if(!BaseClassType::isInitialized_){
	//std::cout << "" << std::endl;

	BaseClassType::set_with_init(ddim,rdim,init);
#if USE_CPPAD
    
	SizeType m = BaseClassType::adseq_.Range(),
	  n = BaseClassType::adseq_.Domain(),
	  i, j, k;

	// set up sparsity pattern of Jacobian
	SP_.initialize(m,n);
	SP_.calc_sparsity(BaseClassType::adseq_); //sparsity pattern of jacobian

	if(m==n){
	  pattNewt_ = SP_;   //same size...
	  pattNewt_.diagonal_never_zero(); //...but maybe more nonzeros
	  assert(pattNewt_.size() == SP_.size());
	}

	//do the computationally expensive suff from CppAD::SparseJacobian here
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
	  CppAD::sparsity_user2internal(inPatt_, SP_.get_p(), m, n, transpose);
 
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
	  CppAD::sparsity_user2internal(inPatt_, SP_.get_p(), m, n, transpose);
 
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

	//============================================================
	//create sparsity data for iteration Matrix
	SizeType nnzNewt = (SizeType)pattNewt_.number_of_nonzeros();
	iterationMatrix_.resize(rdim,ddim,nnzNewt);
	
	mapperIdx_.resize(nnzNewt);
	MtxSelectorType::create_format_data(rdim,ddim,MtxSelectorType::row_iter(iterationMatrix_.index(),iterationMatrix_.position()),MtxSelectorType::col_iter(iterationMatrix_.index(),iterationMatrix_.position()),mapperIdx_.begin(),pattNewt_.get_p() );

	//info 
	std::cout << "-----------------------------------------------------"<<std::endl;
	std::cout << "INFO about G' (CHAR = '"<< CHAR << "'):" << std::endl << "  nnz = "<< nnzNewt << std::endl <<           "  INDEX =    ";
	print_all(iterationMatrix_.index(),iterationMatrix_.index()+nnzNewt);
	std::cout << "  POSITION = ";
	print_all(iterationMatrix_.position(),iterationMatrix_.position()+iterationMatrix_.dimension_of_position_array());
	std::cout << "  MAPPER =   ";
	print_all(mapperIdx_);
	std::cout << "-----------------------------------------------------"<<std::endl;
	
#else
	
#endif
	BaseClassType::isInitialized_ = true; //o.k., initialization done once
      } //end !isInitialized
    }


    const SizeType storage_dimension() const {return K_;} //nnz

    template<class X>
    VecBaseType& evaluate_jacobian(const X& x){
#if USE_CPPAD
      adonis_assert(x.size() == (SizeType)(BaseClassType::adseq_.Domain()));  
      (BaseClassType::adseq_.Domain() <= BaseClassType::adseq_.Range()) ? (BaseClassType::adseq_.sparse_jacobian_for(x,inPatt_,J_,work_)) : (BaseClassType::adseq_.sparse_jacobian_rev(x,inPatt_,J_,work_));
#else
     
#endif
      return J_;
    }

    //! compute \f$ G' := aI - hbS',\f$ where the Jacobian \f$S'\f$ is 
    //! supposed to have already been computed
    template<class COEF>
    MatrixType& compute_G_prime(const T& h, const COEF& a, const COEF& b){
      //only defined for square sparse matrices
      adonis_assert(iterationMatrix_.rows() == iterationMatrix_.cols());

      bool sameSize = false;
      T bh = b*h;
      SizeType k = 0; //positioin in G'
      SizeType diffsz = 0;
    
      ConstSetIterType pattJacIt,
	setIt;        
      SizeType idx;

      for(SizeType i = 0; i < pattNewt_.size(); ++i){
	pattJacIt = SP_[i].begin();
	(pattNewt_[i].size() != SP_[i].size()) ? sameSize = false : sameSize = true;
        
     
        
	for(setIt = pattNewt_[i].begin(); setIt != pattNewt_[i].end(); ++setIt){

	  if(i == *setIt){
	    diffsz += pattNewt_[i].size() - SP_[i].size(); //either zero or one
	    if(!sameSize){
	      //[] fills value array of compressed sparse matrix
	      iterationMatrix_[k] = a*1;
	    }          
	    else{
	      iterationMatrix_[k] = a*1 - bh*J_[k-diffsz];
	    }
	  }
	  else{
	    iterationMatrix_[k] = -bh*J_[k-diffsz];    
	  }
          
	  k++;
	}
    
      }
      return iterationMatrix_;    
    }

    const SparsePatternType& get_pattern() const {return pattNewt_;}
    const MatrixType& get_matrix() const {return iterationMatrix_;}

  private:
    
    SizeType K_;        //number of nonzeros
#if USE_CPPAD
    typedef typename CppAD::internal_sparsity<SetType>::pattern_type PatternType;
    CppAD::sparse_jacobian_work work_;
    PatternType inPatt_;
#else

#endif
    VecBaseType J_;  //will be overwritten by the sparse Jacobian values
    // the diagonal is never empty
    SparsePatternType pattNewt_;
    //outside compiler switch
    SparsePatternType SP_;
    VecIndexType mapperIdx_;  //can store the ordering of the values in a 
                              //sparse matrix format in case a conversion to
                              //another format is required
    MatrixType iterationMatrix_;
  };

} //end namespace

#endif
