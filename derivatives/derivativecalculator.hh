#ifndef VERSATILE_DERIVATIVE_CALCULATOR_API_HH
#define VERSATILE_DERIVATIVE_CALCULATOR_API_HH

#if USE_CPPAD
#include <cppad/cppad.hpp> //include the CPPAD stuff
#else
#include "../noderivativemethods/finitedifferences.hh"
#endif

#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../common/tempusfugit.hh"
#include "../common/fancymessages.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../sparse/sparsitypattern.hh"
#include "../sparse/sparsetmp.hh"


namespace Adonis{

  template<class T, class INTEGER, template<typename D1> 
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
    typedef INTEGER IndexType;


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
	adonis_assert(init.size() == domdim);
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
    // template<class X>
    // VecBaseType evaluate_jacobian(const X& x){
    //   return ref_2_derived().evaluate_jacobian(x);
    // }

  };

  /* ==================================================================== */
  
  //forward declaration
  /**
   * \brief Computation of derivatives of a vector field 
   *
   * \tparam SPARSE boolean CT switch (true: use sparse evaluation and 
   *         corresponding container, false: use dense analogon)
   * \tparam T value type 
   * \tparam SOURCE vector field whose 1st derivatives are to be computed
   * \tparam VEC STL-compliant vector type for storage
   * \tparam CHAR select matrix format (only relevant when SPARSE is true. 
   *              At disposal: 'c','C' and 'r','R' for compressed sparse 
   *                            column and row format, respectively.
   * \tparam INTEGER index type of (sparse) matrix container
   *
   */
  template<bool SPARSE, class T, template<typename D1> 
	    class SOURCE, template<typename D2, class A = std::allocator<D2> > 
	   class VEC, char CHAR, class INTEGER = int>
  class JacobianMatrixCalculator: public AbstractDerivativeCalculatorGeneral<T,INTEGER,SOURCE,VEC,JacobianMatrixCalculator<SPARSE,T,SOURCE,VEC,CHAR,INTEGER> >{}; 


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
   * NOTE 1: Up to now an error is thrown, if you attempt to use sparse FDs, 
   *       simply because I was too lazy to implement it. Feel free to complete.
   *
   * NOTE 2: We assume that each operation is executed in a C/C++-style manner,
   *         i.e. row-wise style. If CHAR is 'C' or 'c', one has to reorder the
   *         values of the values array for compressed sparse matrices. 
   *
   * NOTE 3: needs about 8 times nnz of space. Thus from nnz >= 9, the memory
   *         savings of using the sparse format are increasing.
   */
  template<class T, template<typename D1>  //! SPARSE evaluation
	    class SOURCE, template<typename D2, class A = std::allocator<D2> > 
	   class VEC, char CHAR, class INTEGER>
  class JacobianMatrixCalculator<true,T,SOURCE,VEC,CHAR,INTEGER>: public AbstractDerivativeCalculatorGeneral<T,INTEGER,SOURCE,VEC,JacobianMatrixCalculator<true,T,SOURCE,VEC,CHAR,INTEGER> >{
  public:
    typedef AbstractDerivativeCalculatorGeneral<T,INTEGER,SOURCE,VEC,JacobianMatrixCalculator<true,T,SOURCE,VEC,CHAR,INTEGER> > BaseClassType;
    typedef typename BaseClassType::SizeType SizeType;
    typedef typename BaseClassType::VecBaseType VecBaseType;
    typedef typename BaseClassType::FunType FunType;
    typedef typename BaseClassType::VecIndexType VecIndexType;
    typedef typename BaseClassType::BaseType BaseType;
         
  private:
    typedef SparseMatrixSelector<CHAR,T,INTEGER> MtxSelectorType;
    typedef std::set<SizeType> SetType;
    typedef SparsityPattern<SetType> SparsePatternType;
    typedef typename SetType::const_iterator ConstSetIterType;

  public:
    typedef typename MtxSelectorType::SparseMatrixType MatrixType;
    typedef SparsePatternType PatternType;
    static const char Value = CHAR;

    //the 2nd argument is a dummy here
    JacobianMatrixCalculator(bool determineGprime = true, const FunType& fun = FunType()):status_(0), determineGprime_(determineGprime), firstGprimeComputation_(false)
{
   FancyMessages().nice_output("SPARSE derivatives in use: ",33);
#if USE_CPPAD
      FancyMessages().nice_output("\n AD Jacobian obtained via AD CppAD.");
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
	std::cout <<std::endl << "SETTING UP SPARSITY PATTERN... \""<< std::endl; 
	CPUTime<BaseType> stepCPUtm;
	BaseType stepElapsedCpuTime(0);
	stepCPUtm.start(false); //print nothing to screen

	SizeType m = BaseClassType::adseq_.Range(),
	  n = BaseClassType::adseq_.Domain();
	std::cout << "1.) Compute sparsity pattern of Jacobian"<<std::endl;
	// set up sparsity pattern of Jacobian, this is always done
	SP_.initialize(m,n);
	SP_.calc_sparsity(BaseClassType::adseq_); //sparsity pattern of jacobian

	std::cout << "2.) Create triplet info for CppAD::SparseJacobianForward/Reverse drivers" <<std::endl;
	// create arguments for SparseJacobianForward/Reverse driver in 
	// triplet format. This is solely needed for evaluating the Jacobian!
	SizeType nnzJac = (SizeType)SP_.number_of_nonzeros();
	row_.resize(nnzJac);
	col_.resize(nnzJac);
	J_.resize(nnzJac);

	// create triplet row and columns vector respectively
	SizeType k = 0;
	bool nzstsz;
	for(SizeType i = 0; i < SP_.size(); ++i){
	  if(SP_[i].size() == 0)
	    nzstsz = false;
	  else
	    nzstsz = true;
	  for(ConstSetIterType cit = SP_[i].begin(); cit != SP_[i].end(); ++cit){
	    if(nzstsz){  //row contains elements
	      row_[k] = i;
	    }
	    col_[k] = *cit;
	    k++;
	  }
	}

#ifdef PRINT_2_SCREEN //really too much output for large-scale scenarios
	std::cout << "Jacobian pattern: "<< std::endl;
	std::cout << SP_ << std::endl;

	//this is needed for the CppAD members SparseJacobianForward/Reverse
	std::cout << "TRIPLET FORMAT index arrays:" << std::endl;
	std::cout << "row_ = "; 
	for(SizeType l = 0; l < row_.size(); ++l)
	  std::cout << row_[l] << "  ";
	std::cout << std::endl;
	std::cout << "col_ = "; 
	for(SizeType l = 0; l < col_.size(); ++l)
	  std::cout << col_[l] << "  ";
	std::cout << std::endl;
#endif

	//! fill index arrays of compressed storage matrix. If compressed
	//! column storage is selected, one has to permute the value array
	//! via the mapperIdx_
	//============================================================
	//! create sparsity data for iteration matrix \f$G'\f$ -- SQUARE only
	//! and only when it is desired.
	if( (determineGprime_==true) && (m==n)){
	  std::cout << "3.) Create sparsity pattern for G'" << std::endl;
	  pattNewt_ = SP_;   //same size...
	  pattNewt_.diagonal_never_zero(); //...but maybe more nonzeros
	  assert(pattNewt_.size() == SP_.size());

	  SizeType nnzNewt = (SizeType)pattNewt_.number_of_nonzeros();
	  iterationMatrix_.resize(rdim,ddim,nnzNewt);
	  
	  MtxSelectorType::resize(vals_,nnzNewt); //resize or not, dep. on CHAR

	  std::cout << "4.) Create indices for CS"<<CHAR<< " format, including mapper for values. This may take a while..." << std::endl;
	  mapperIdx_.resize(nnzNewt);
	  //! fill indices and mapper from given sparsity pattern. Can be either
	  //! CSR or CSC format, but here we mainly consider the CSC scheme.
	  MtxSelectorType::create_format_data(rdim,ddim,MtxSelectorType::row_iter(iterationMatrix_.index(),iterationMatrix_.position()),MtxSelectorType::col_iter(iterationMatrix_.index(),iterationMatrix_.position()),mapperIdx_.begin(),pattNewt_.get_p() );

	  //info 
	  std::cout << std::endl<< "-----------------------------------------------------"<<std::endl;
	  std::cout << "INFO about G' (CHAR = '"<< CHAR << "'):" << std::endl;
#ifdef PRINT_2_SCREEN //this is really much for large-scale matrices
	 
	  std::cout  << "  nnz = "<< nnzNewt << std::endl <<           "  INDEX =    ";
	  print_all(iterationMatrix_.index(),iterationMatrix_.index()+nnzNewt);
	  std::cout << "  POSITION = ";
	  print_all(iterationMatrix_.position(),iterationMatrix_.position()+iterationMatrix_.dimension_of_position_array());
	  std::cout << "  MAPPER =   ";
	  print_all(mapperIdx_);
#endif
     
	  std::cout << "____________________________"<<std::endl;
	  std::cout << "Sparsity Pattern info: size: "<< pattNewt_.size() << "  #(nonz): "<< pattNewt_.number_of_nonzeros() << " (out of "<< pattNewt_.size() << " x "<< pattNewt_.size() << " entries)"<<std::endl; 
	  std::cout << "For comparison: Jacobian has "<< SP_.number_of_nonzeros() << " nonzeros." << std::endl;
      std::cout << std::endl<< "-----------------------------------------------------"<<std::endl;
     
	  
	}
	else{ //!consider only JACOBIAN here -- default case
	  if((determineGprime_ == true) && (m!=n)){
	    FancyMessages Fmsg;
	    Fmsg.nice_output("**** NOTE: No SQUARE matrix considered. You cannot compute G' in the following!",35);
	  }
	  std::cout << "3.) Jacobian considered."<<std::endl;
	  iterationMatrix_.resize(rdim,ddim,nnzJac);
	
	  MtxSelectorType::resize(vals_,nnzJac); //resize or not, dep. on CHAR
	  std::cout << "4.) Create indices for CS"<<CHAR<< " format, including mapper for values. This may take a while..." << std::endl;
	  //! create indices for either compressed row or column format
	  mapperIdx_.resize(nnzJac);
	  MtxSelectorType::create_format_data(rdim,ddim,MtxSelectorType::row_iter(iterationMatrix_.index(),iterationMatrix_.position()),MtxSelectorType::col_iter(iterationMatrix_.index(),iterationMatrix_.position()),mapperIdx_.begin(),SP_.get_p() );
#ifdef PRINT_2_SCREEN //may be really large
	  //nonzeros not yet filled with meaningful values
	  std::cout << "Jacobian<'"<<CHAR<<"'> = "<< iterationMatrix_ << std::endl;
	  std::cout << "mapperIdx:  ";
	  print_all(mapperIdx_);
#endif
	  std::cout << "____________________________"<<std::endl;
	  std::cout << "Sparsity Pattern info: size: "<< SP_.size() << "  #(nonz): "<< SP_.number_of_nonzeros() << " (out of "<< SP_.size() << " x "<< SP_.size() << " entries)"<<std::endl; 
	}
#else
	
#endif
	BaseClassType::isInitialized_ = true; //o.k., initialization done once   
	stepElapsedCpuTime = stepCPUtm.stop(false);
	 std::cout << "Spent time for pattern generation of sparse matrix: "<< time_in_human_readable_form(stepElapsedCpuTime) << std::endl;
	 stepElapsedCpuTime = BaseType(); //reset
	std::cout << "Sparse initialization DONE!"<<std::endl;
      } //end !isInitialized
    }


    const char format_info() const {return CHAR;}

    const BaseType* value_pointer() const {return iterationMatrix_.values();}
    BaseType* value_pointer()  {return iterationMatrix_.values();}
    
    const SizeType dimension_of_values() const {return (SizeType)iterationMatrix_.nonz();}
    
    template<class X>
    VecBaseType& evaluate_jacobian(const X& x, bool reassign = true){
#if USE_CPPAD
      //std::cout << "x.size() = "<< x.size() << "   BaseClassType::adseq_.Domain() = "<< BaseClassType::adseq_.Domain() << std::endl;
      adonis_assert(x.size() == (SizeType)(BaseClassType::adseq_.Domain()));  
      (BaseClassType::adseq_.Domain() <= BaseClassType::adseq_.Range()) ? ( status_ = BaseClassType::adseq_.SparseJacobianForward(x,SP_.get_p(),row_,col_,J_,work_)) : (status_ = BaseClassType::adseq_.SparseJacobianReverse(x,SP_.get_p(),row_,col_,J_,work_));
      
      //! more effort, but currently CppAD does not support iterators for 
      //! the SparseJacobianForward/Reverse drivers, so I have to assign it to
      //! the compressed jacobian structure in addition...
      if(reassign){
	for(SizeType k = 0; k < J_.size(); ++k){
	  smart_assign(iterationMatrix_[k],J_[k]);
	}
      }
#else
     
#endif
      return J_;
    }


    //! compute \f$ G' := aI - hbS',\f$ where the Jacobian \f$S'\f$ is 
    //! supposed to have already been computed.
    //! CAUTION: the values are supposed to be in CSR format since the 
    //!          sparsity pattern only applies to this format.
    template<class COEF>
    MatrixType& compute_G_prime(const T& h, const COEF& a, const COEF& b){
      //only defined for square sparse matrices
      adonis_assert(iterationMatrix_.rows() == iterationMatrix_.cols());
      if(firstGprimeComputation_ == false){
	std::cout << "1st evaluation of G'" << std::endl;
      //now always set to true to prevent the above message from printing to screen again
	firstGprimeComputation_ = true;    
      }
      if(determineGprime_){
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

	    if(i == *setIt){ //! diagonal elements
	      diffsz += pattNewt_[i].size() - SP_[i].size(); //either zero or one
	      if(!sameSize){
		//[] fills value array of compressed sparse matrix
		iterationMatrix_[k] = a*1;
	      }          
	      else{
		iterationMatrix_[k] = a*1 - bh*J_[k-diffsz];
	      }
	    }
	    else{ //! off-diagonal elements
	      iterationMatrix_[k] = -bh*J_[k-diffsz];    
	    }
          
	    k++;
	  }
    
	}
      }
      return iterationMatrix_;    
    }


    //! in the column compressed case ('c','C'), reorder value array and
    //! reassign to value field of compressed matrix container.
    //! This must be invoked shortly before you wanna start an UMFPACK routine.
    void reorder_values(){
      MtxSelectorType::copy(vals_,iterationMatrix_.values()); //hardcopy to vals
      MtxSelectorType::reorder_values(iterationMatrix_.values(),vals_,mapperIdx_);
    }

    //! output pattern either for G' or S', depending on your choice
    const SparsePatternType& get_pattern() const {
      //!                            underlying pattern must be  for square mtx
      return ( ((determineGprime_==true) && (SP_.rows() == SP_.cols())) ? pattNewt_ : SP_ );
    }
    const MatrixType& get_matrix() const {return iterationMatrix_;}
    MatrixType& get_matrix() {return iterationMatrix_;}

    const VecBaseType& get_jacobian() const {return J_;}
    VecBaseType& get_jacobian() {return J_;}

    const VecIndexType& get_mapper_idx() const {return mapperIdx_;} 

  private:
    
    SizeType status_;        //number of nonzeros
    bool determineGprime_, firstGprimeComputation_;
#if USE_CPPAD
    CppAD::sparse_jacobian_work work_;
    CppAD::vector<SizeType> row_, col_;

#else

#endif
    VecBaseType J_,  //will be overwritten by the sparse Jacobian values
      vals_;         //may be used as intermediate variable to hold values
    // the diagonal is never empty
    SparsePatternType pattNewt_;
    SparsePatternType SP_;
    VecIndexType mapperIdx_;  //can store the ordering of the values in a 
                              //sparse matrix format in case a conversion to
                              //another format is required
    MatrixType iterationMatrix_;
  };


//***************************************************************************//
//***************************************************************************//
//***************************************************************************//
 template<class T, template<typename D1>  //! DENSE evaluation
	    class SOURCE, template<typename D2, class A = std::allocator<D2> > 
	  class VEC, char CHAR, class INTEGER> //! 4th arg devoid of meaning for dense matrices
 class JacobianMatrixCalculator<false,T,SOURCE,VEC,CHAR,INTEGER>: public AbstractDerivativeCalculatorGeneral<T,INTEGER,SOURCE,VEC,JacobianMatrixCalculator<false,T,SOURCE,VEC,CHAR,INTEGER> >{
  public:
   typedef AbstractDerivativeCalculatorGeneral<T,INTEGER,SOURCE,VEC,JacobianMatrixCalculator<false,T,SOURCE,VEC,CHAR,INTEGER> > BaseClassType;
    typedef typename BaseClassType::SizeType SizeType;
    typedef typename BaseClassType::VecBaseType VecBaseType;
    typedef typename BaseClassType::FunType FunType;
    typedef typename BaseClassType::VecIndexType VecIndexType;
    typedef typename BaseClassType::BaseType BaseType;


   //dense matrix, row-wise storage
    typedef VecBaseType MatrixType;
    //!there is no sparsity pattern information. To stay conform, just define it
    //!to be a T-type
    typedef T PatternType;  

  public:
   static const char Value = 'D'; //dense
  
   JacobianMatrixCalculator(bool determineGprime, FunType& fun):determineGprime_(determineGprime),dim_(0),firstGprimeComputation_(false)
 #if USE_CPPAD									 
 #else
								,fun_(fun),eps_(std::pow(std::numeric_limits<BaseType>::epsilon(),1./3))
 #endif
{
  BaseClassType::isInitialized_ = false;
  FancyMessages().nice_output("DENSE derivatives in use: ",33);
#if USE_CPPAD
      FancyMessages().nice_output("\n AD Jacobian obtained via AD CppAD.");
#else

      FancyMessages().nice_output("\n FD APPROXIMATION used.");
#endif     
      BaseClassType::isInitialized_ = false;
}

    template<class W>
    void initialize(SizeType ddim, SizeType rdim, const W& init){
      if(!BaseClassType::isInitialized_){
	if((determineGprime_==true) && (ddim==rdim)){
	  dim_ = ddim;  //this is only needed for the calculation of G'
	}
	std::cout << "dim_ = "<< dim_ << std::endl;
	
	iterationMatrix_.resize(rdim*ddim); //m x n matrix storage

	BaseClassType::set_with_init(ddim,rdim,init);
#if USE_CPPAD	

#else

#endif
      
	BaseClassType::isInitialized_ = true;
	std::cout << "Dense initialization done!"<<std::endl;
      }
    }
 
    //! this is just a dummy function since tmpl arg CHAR is devoid of meaning here
    //! 'D' like dense
   const char format_info() const {return 'D';} //{return CHAR;} 

   const BaseType* value_pointer() const {return &iterationMatrix_[0];}
   BaseType* value_pointer() {return &iterationMatrix_[0];}

   const SizeType dimension_of_values() const {return iterationMatrix_.size();}

   //! 2nd argument does not has any effect here
   template<class X> 
   VecBaseType& evaluate_jacobian(const X& x, bool reassign = true){
#if USE_CPPAD
     adonis_assert(BaseClassType::adseq_.Domain() == x.size());
     iterationMatrix_ = BaseClassType::adseq_.Jacobian(x);
#else
      //central FD for the sake of accuracy
     ArgumentTypeDependency<IsContainer<X>::Value,VecBaseType>::central_fd_jacobian(iterationMatrix_,x,fun_.domain_dim(),fun_.dim(),eps_,fun_);
#endif
      return iterationMatrix_;
    }

    //! compute \f$ G' := aI - hbS',\f$ where the Jacobian \f$S'\f$ is 
    //! supposed to have already been computed
    template<class COEF>
    MatrixType& compute_G_prime(const T& h, const COEF& a, const COEF& b){
      if(firstGprimeComputation_ == false){
	std::cout << "1st evaluation of G'" << std::endl;
      //now always set to true to prevent the above message from printing to screen again
	firstGprimeComputation_ = true;    
      }
      if(determineGprime_){
	adonis_assert(dim_!=0); //dim only filled when G' is to be computed
	iterationMatrix_ *= -(b*h);
	update_diagonal<AddBasicElements>(iterationMatrix_,dim_,a);
      }
      return iterationMatrix_;
    }


    //always empty here, since this is only needed in conjunction with sparsity
    void reorder_values(){} 

    //fake function here, just return a T
   const PatternType get_pattern() const {return 42;}

   const MatrixType& get_matrix() const {return iterationMatrix_;}
   MatrixType& get_matrix() {return iterationMatrix_;}

   const VecBaseType& get_jacobian() const {return iterationMatrix_;}
   VecBaseType& get_jacobian() {return iterationMatrix_;}

   //! dummy function
   const SizeType get_mapper_idx() const {return 0;} 

 private:
   bool determineGprime_;
   SizeType dim_;
   bool firstGprimeComputation_;
#if USE_CPPAD

#else
   FunType& fun_;
   BaseType eps_;
#endif	
    MatrixType iterationMatrix_;
}; 

} //end namespace

#endif
