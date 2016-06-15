

////////////////////////////////////////////////////////////////////////
/**
* \brief The identity of the the compressed matrix.
*  Usage: some applications require the specification of the compressed matrix type explicitly. As 
*         'C','c' as well as 'R','r' both are legal and reasonable identifiers for the compressed 
*         sparse column and compressed sparse row format, respectively, we make the upper case the *         default. This becomes especially beneficial when partial specializations become larger.
* 
*  \code
    const char Chr = 'r';  //lower case, but the specialization is only laid out for the upper case 'R'

	class DoSomething<CompressedSparseStorage<Chr>::Value>.... //CT conversion to upper case
	*  \endcode 
*/
template<char CHAR> 
class CompressedSparseStorageIdentity{
public:
	static const char Value = 'D';    //default is dense, i.e. if CHAR does not hit any specification,
                                      //density is assumed          
	};	

//! column compressed storage
template<>
class CompressedSparseStorageIdentity<'C'>{ //by default, an upper case is considered standard
public:
	static const char Value = 'C';
};    

template<>
class CompressedSparseStorageIdentity<'c'>{  //lower case 
public:
	static const char Value = 'C';           //upper case
};    

//! row compressed storage (= Yale format)
template<>
class CompressedSparseStorageIdentity<'R'>{ //by default, an upper case is considered standard
public:
	static const char Value = 'R';
};    

template<>
class CompressedSparseStorageIdentity<'r'>{  //lower case 
public:
	static const char Value = 'R';           //upper case
};    

template<>
class CompressedSparseStorageIdentity<'Y'>{ //Yale format (is the same as 'R', 'r')
public:
	static const char Value = 'R';
};    

template<>
class CompressedSparseStorageIdentity<'y'>{  //lower case 
public:
	static const char Value = 'R';           //upper case
};    


//////////////////////////////////////////////////////////////////////////
/**
* \brief Select special functionalities at compile time. 
* NOTE: this template meta program only checks the FLAG, not the type of the matrix under consideration *       (e.g. dense, row or column-compressed form). Therefore, don't use it as stand-alone!
*
* Usage: for example in a partial class specialization, we can select whether the sparsity pattern within
*        UMFPACK is computed once, or in each iteration.
*/
template<bool FLAG> class LSSFunctionalityHandler;

template<>
class LSSFunctionalityHandler<true>{ //compute pattern ONCE
public:
	
	//only comute pattern for UMFPACK once
	template<class MTX, class V>
	static inline int pattern_once(const MTX& mtx, V& ax, V& az, void* Symbolic, const double Ctrl[], const double Info[]){
		#if USE_UMFPACK
		//!we assume that Ax is not present (a (double*) NULL pointer). This means that any entry in A is assumed to be "large"
		//! no splitting used for complex system matrix when only computed once 
			//return UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::symbolic_analysis(mtx.rows(),mtx.cols(),mtx.position(),mtx.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::proper_iterator(mtx.values(),&ax_[0]),&az_[0],&Symbolic_,Ctrl_,Info_);
			return UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::symbolic_analysis(mtx.rows(),mtx.cols(),mtx.position(),mtx.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::proper_iterator((typename MTX::value_type*) NULL,(double*) NULL, (double*) NULL,&Symbolic_,Ctrl_,Info_);
		#else
		ADONIS_ERROR(UmfpackError,"This is meant only for usage with some function in conjunction with UMFPACK.");
		#endif
	}

	template<class MTX, class V> //do nothing in iteration
	static inline int pattern_in_each_iteration(const MTX& mtx, V& ax, V& az, void* Symbolic, const double Ctrl[], const double Info[]){
	 return 42;
	}
	
	//do not delete 
	void free_symbolic(void* Symbolic){}
	
};

template<>
class LSSFunctionalityHandler<false>{ //in EACH iteration
public:

    template<class MTX, class V> 
	static inline int  pattern_once(onst MTX& mtx, V& ax, V& az, void* Symbolic, const double Ctrl[], const double Info[]){
	 return 42;
	}  //do nothing
	
	//!compute pattern in each iteration anew, assume that the matrix has already been split into ax_ and az_
	template<class MTX, class V> 
	static inline int pattern_in_each_iteration(const MTX& mtx, V& ax, V& az, void* Symbolic, const double Ctrl[], const double Info[]){
	#if USE_UMFPACK
         
 		//return umfpack_di_symbolic(mtx.rows(), mtx.cols(), mtx.position(), mtx.index(), mtx.values(), &Symbolic, Ctrl,Info);
			
			return UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::symbolic_analysis(mtx.rows(),mtx.cols(),mtx.position(),mtx.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::proper_iterator(mtx.values(),&ax[0]),&az[0],&Symbolic,Ctrl,Info);
		
		
		#else
		ADONIS_ERROR(UmfpackError,"This is meant only for usage with some function in conjunction with UMFPACK.");
		#endif
	}
	
	void free_symbolic(void* Symbolic){
	UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::free_symbolic(&Symbolic);
	}
};



//////////////////////////////////////////////////////////////////////////
/**
* \brief Class for solving <B>square</B> linear systems \f$A\cdot x = b, A \in R^{n \times n} \ \mathrm{or}\ A \in C^{n \times n}.\f$
* GENERAL case (default): dense systems
* 
* \tparam MTXTYPE matrix type \f$A\f$
* \tparam RHSTYPE right hand side type \f$b\f$
*
* Usage:
* \code
 int dim = 4
  Matrix<double> A(dim,dim);
  //...fill A
  MyVec<double> b(dim);;
  //...fill b
  
  //the last template argument is (currently) a dummy here
  SquareLinearSystem<Matrix<double>,MyVec<double>,'D',true> SQLS(A,b,1);  //construct object at beginning
   
   //... some code that updates A and b e.g. during an iterative process
  std::cout << "solution = "<< SQLS.solve() << std::endl;                 //solve system and output sol.
   
* \endcode
*/
template<CLASS MTXTYPE, class VTYPE, char CHAR, bool FLAG>
class SquareLinearSystem{
public:
	typedef MTXTYPE MatrixType;
		typedef typename MTX::value_type value_type;	
        typedef typename TypeAdapter<value_type>::BaseType BaseType;	
		typedef typename VTYPE::value_type Rhs_value_type;
		typedef std::size_t SizeType;
		typedef VTYPE VecType;
		
	SquareLinearSystem(MTXTYPE& mtx, VTYPE& b, int nrhs):A_(mtx),b_(b),status_(0),nrhs_(nrhs),dim_(int(b.size())), ipiv_(dim_)
	#if USE_LAPACK
   ,lda_(std::max(1,dim_)),ldb_(std::max(1,dim_)),trans_('N')
   #else
   , xx_(dim_)
	#endif
	{
	 adonis_assert(typeid(value_type) == typeid(Rhs_value_type)); //same value_types mandatory
	 adonis_assert(mtx.size() == ntimes<2>(b.size()));
			FancyMessages().nice_output("Solve square linear system in DENSE format");
			#if USE_LAPACK
			FancyMessages().nice_output("LAPACK is used.");
			 if(!F77Type<TypeTraits<value_type>::Value>::IsComplex){ //real
	if(useTransposedForm)
	  trans_ = 'T';
      }
      else{                                                   //complex
	if(useTransposedForm)
	  trans_ = 'C';
      }
			#else
			ADONIS_ERROR(ImplementationError,"No alternative to solve dense systems has been specified so far ;) ");
			#endif
	}

	
	VTYPE& solve(){
	#if USE_LAPACK
	//!perform LU-decomposition
	F77Type<TypeTraits<value_type>::Value>::GETRF(&dim_,&dim_,&mtx_[0],&lda_,&ipiv_[0],&status_);
      if(status_ != 0)
	    ADONIS_ERROR(LapackError,"An error occured within Lapack's "<<typeid(value_type).name() << "xtrf");
	  //!solve square system
	  F77Type<TypeTraits<value_type>::Value>::GETRS(&trans_,&dim_,&nrhs_,&mtx_[0],&lda_,&ipiv_[0],&b_[0],&ldb_,&status_);
#ifndef NDEBUG
      if(status_ != 0)
	ADONIS_ERROR(LapackError,"An error occured within Lapack's "<<typeid(value_type).name() << "xtrs");
#endif
	#else
	//! like in LAPACK, the mtx_ and b_ are overwritten by the PLU decomposition and the solution vector, respectively.
	//GAUSSIAN elimination with partial pivoting, taken from Daoqi Yang's book. It should work with 
	// complex types as well
     //MTXTYPE tmpx = mtx_; //local copy

  // ipiv_ contains the pivot elements now  
  for (int k = 0; k < dim_; ++k){ 
     ipiv_[k] = k;  //overwrite ipiv_
	 xx_[k] = value_type(); //reset before each call 
	 }

  int nrowsmone = nrows - 1;
  for (int k = 0; k < nrowsmone; ++k) {  // main loop

    // find the pivot in column k in rows pvt[k], 
    // pvt[k+1], ..., pvt[n-1]
    int pc  = k; 
    BaseType aet = Abs(*(mtx + k*dim_ + k) ;//(tmpx[pvt[k]][k]);
    for (int i = k + 1; i < dim_; ++i) {
      if (Abs(*(mtx_ + ipiv_[i]*dim_ + k) > aet){ //[pvt[i]][k]) > aet) {
        aet = Abs(*(mtx_ + ipiv_[i]*dim_+k)  );//tmpx[pvt[i]][k]); 
        pc = i;
      }
    }
    if (is_zero(aet)) 
      ADONIS_ERROR(DivisionByZeroError, "Pivot is zero in Gaussian elimination with partial pivoting.\n   ==> matrix is SINGULAR."); 
    if (pc != k) std::swap(ipiv_[k], ipiv_[pc]);  //defined in algorithm
    int pvtk = ipiv_[k];                  // pivot row
    value_type pivot = *(mtx_ + pvtk*dim_ + k); //tmpx[pvtk][k];            // pivot

    // now eliminate column entries logically 
    // below tmpx[pvt[k]][k]
	value_type mult;
    for (int i = k + 1; i < dim_; ++i) {
      int pvti = ipiv_[i];
      if ( *(mtx_ + pvti*dim_ + k) != value_type()){ //tmpx[pvti][k] != 0) {
         mult = *(mtx_ + pvti*dim_ +k)/pivot;  //tmpx[pvti][k]/pivot;
        *(mtx_ + pvti*dim_ + k) = mult;  //tmpx[pvti][k] = mult;
        for (int j = k + 1; j < dim_; ++j) 
            *(mtx_ + pvti*dim_+j) -= mult*(*(mtx_+pvtk*dim_+j));//tmpx[pvti][j] -= mult*tmpx[pvtk][j];
      }
    }
  }

  // forwad substitution for L y = Pb.
  for (int i = 1; i < dim_; ++i)  
    for (int j = 0; j < i; ++j) 
      b_[ipiv_[i]] -= (*(mtx_ + ipiv_[i]*dim_+j))*b_[ipiv_[j]];//tmpx[pvt[i]][j]*bb[pvt[j]];

  // back substitution for Ux = y
  for (int i = nrowsmone; i >= 0; --i) {
    for (int j = i+1; j < dim_; ++j) 
      b_[ipiv_[i]] -=  (*(mtx_ + ipiv_[i]*dim_ + j)*xx_[j]); //tmpx[pvt[i]][j]*xx[j];
    xx_[i] = b_[ipiv_[i]] / ( *(mtx_ + ipiv_[i]*dim_ + i) );//tmpx[pvt[i]][i];
  }

  b_ = xx_;             // put solution
	#endif
    return b_;
	}
	
	
private:
MTXTYPE& A_;  //use references here since we want to use LAPACK here in the first place
VTYPE& b_;
int status_, 
nrhs_, dim_;
IndexVecType ipiv_;

#if USE_LAPACK
int lda_,ldb_;
char trans_;     
    //Lapack is column-oriented, hence we use the transposed form (C-style)
    enum{useTransposedForm = true}; 
	#else 
	//! use Gaussian elimination with partial pivoting
	VTYPE xx_; //stores solution in correct otder 
#endif

//! make copy stuff private
	SquareLinearSystem(const LinearSolver&);
    SquareLinearSystem& operator=(const LinearSolver&);
};


//!=========================================================================================
//!============== partial specializations for SPARSE square systems ========================
//!=========================================================================================
/**
* \brief Solves <B>sparse</B> square linear system \f$A\cdot x = b,\f$ where \f$A\f$ is given in 
* compressed sparse column format (the UMFPACK-style). 
* \tparam MTXTYPE compressed sparse column object
* \tparam VTYPE right hand side vector type
* \tparam FLAG if true, the pattern is computed only once, otherwise in each call to the 'solve' routine
*/
template<CLASS MTXTYPE, class VTYPE, bool FLAG>
class SquareLinearSystem<MTXTYPE,VTYPE,'C',FLAG>{
public:
		typedef MTXTYPE MatrixType;
		typedef typename MTX::value_type value_type;	
        typedef typename TypeAdapter<value_type>::BaseType BaseType;	
		typedef typename VTYPE::value_type Rhs_value_type;
		typedef std::size_t SizeType;
		typedef VTYPE VecType;

		//when Complex linear algebra is applied, we need two additional vectors storing the real and imaginary parts of the complex coefficients
	    typedef ExprTmpl::MyVec<SizeType> IndexVecType;
    typedef ExprTmpl::MyVec<BaseType> BaseVecType;

	
		SquareLinearSystem(MTXTYPE& mtx, VTYPE& b, int nrhs = 1):A_(mtx),b_(b),status_(0),nrhs_(nrhs)
		{
			adonis_assert(typeid(value_type) == typeid(Rhs_value_type)); //same value_types are mandatory
			FancyMessages().nice_output("Solve square linear system in compressed sparse column (CSC) format");
			#if USE_UMFPACK
			FancyMessages().nice_output("UMFPACK is used.");
			adonis_assert(nrhs == 1); //only this seems to be possible within Umfpack at present
			if(mtx.nonz()! = 0){
            //only resize when size is not equal to nonz and if complex systems are to be solved!
            UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::resize(mtx.nonz(),ax_);
      UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::resize(mtx.nonz(),az_);

      UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::resize(b.size(),xz_);
      UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::resize(b.size(),bz_);
			}
			//! see <a> http://www.math.umbc.edu/~rouben/2003-09-math625/umfpack-ex3.c</a>
				//*Symbolic_ = NULL;
				//*Numceric_ = NULL;
				status_ =  LSSFunctionalityHandler<FLAG>::pattern_once(A_,ax_,az_,&Symbolic_,Control_,Info_);
				status_umfpack(status_,0); //check status -- throw warning or even error when UMFPACK_OK status isn't met
				#else
				ADONIS_ERROR(ImplementationError,"No alternative to solve sparse systems based on the CSC format has been specified so far ;) ");
			#endif	
		} //end construction
		
		//destruction
		~SquareLinearSystem(){
		#if USE_UMFPACK
		//clear me- if not delete yet
		if(&Symbolic_ != NULL)
			UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::free_symbolic(&Symbolic_);
			if(&Numeric_ != NULL)
			  UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::free_numeric(&Numeric_);  //??????????
		#endif
		}
	 
		
		VTYPE& solve(){
		#if USE_UMFPACK
			//! split when complex numbers are used
			UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::split_complex_container_into_real_values(mtx.nonz(),mtx.values(),&ax_[0],&az_[0]);
			
			//! perform symbolic operation in each call if not computed at the very beginning
			status_ = LSSFunctionalityHandler<FLAG>::pattern_in_each_iteration(A_,ax_,az_,&Symbolic_,Ctrl_,Info_);			
			status_umfpack(status_,0); //everything o.k.?
			
			status_ =  UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::lu_decomposition(A_.position(),A_.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::proper_iterator(A_.values(),&ax_[0]),&az_[0],Symbolic_,&Numeric_,Ctrl_,Info_);  //umfpack_di_numeric(A_.position(),A_.index(), A_.values(),Symbolic_,Numeric_,Control_,Info_);
			status_umfpack(status_,1); //everything o.k.?
			UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::solve(status_,A_.position(),A_.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::proper_iterator(A_.values(),&ax_[0]),&az_[0],&xx_[0],&xz_[0],UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::proper_iterator(&b_[0],&bx_[0]),&bz_[0],Numeric_,Ctrl_,Info_);
		    status_umfpack(status_,2);//umfpack_di_solve(UMFPACK_A,A_.position(),A_.index(),A_.values(),solution.begin(),Numeric_,Control_,Info_);
			UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::join((int)b_.size()*nrhs_,&b_[0],&xx_[0],&xz_[0]);
		
		//================TODO:  free Symbolic if FLAG not true and Numeric here ?????????
	       LSSFunctionalityHandler<FLAG>::free_symbolic(&Symbolic_); //only delete it when FLAG = false
		   UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::free_numeric(&Numeric_);  //delete LU decomposition object	
		#else
		
		#endif
		return b_;
	}
	
private:
MTXTYPE& A_;  //references since A and b are updated during iteration. They remain unaffected by UMFPACK
VTYPE& b_;
int status_, //status of the solver, depends of the specific implementation
nrhs_;       //number of RHSs; UMFPACK can currently not handle multiple right hand sides

#if USE_UMFPACK
//!That's Umfpack's C-style: a void pointer can point to <I>any</I> object.
//! In C++, the use of templates would be more convenient
void* Symbolic_;
void* Numeric_;

//! both for double an complex routines, these arrays are of type double
double Control_[UMFPACK_CONTROL];
double Info_[UMFPACK_INFO];

//additional storage for complex square systems
BaseVecType ax_, az_,
      xx_, xz_,
      bx_, bz_;

//make this private and only usable when UMFPACK usage is switched on
void status_umfpack(const int status, int kind){
    if(status != UMFPACK_OK){  //UMFPACK reports an error. Give notice about kind of error
	std::string str = "umfpack_*_";
	if(kind == 0)
		str = "symbolic";
    else if (kind == 1)
	    str = "numeric";
		else if (kind == 2)
		  str = "solve";
		  else
		   str = "???";
	switch(status){
	case 42:
	//do nothing, that's my notification that a TMP is used not used.
	break;
    case UMFPACK_ERROR_argument_missing:
      ADONIS_ERROR(UmfpackError,"One or more required arguments is missing in '"<<str<<"'.");  
      break;
    case UMFPACK_ERROR_out_of_memory:
      ADONIS_ERROR(UmfpackError, "Insufficient memory to perform the analysis in '"<<str<<"'.");
    break;
    //! symbolic part
    case UMFPACK_ERROR_n_nonpositive:
        ADONIS_ERROR(UmfpackError,"n is less than or equal to zero in '"<< str << "'.");
        break;
        case UMFPACK_ERROR_invalid_matrix:
        ADONIS_ERROR(UmfpackError,"Invalid matrix encountered in '" << str<< "'.\n   For details, please check the UMFPACK manual\n   http://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK/Doc/UserGuide.pdf");
        break:
        case UMFPACK_ERROR_internal_error:
        ADONIS_ERROR(UmfpackError,"Something very serious went wrong. This is a BUG.\n    Please contact DrTimothyAldenDavis@gmail.com");
        break;
        //! numeric part
        case UMFPACK_WARNING_singular_matrix
        FancyMessages().nice_output("Umfpack WARNING: Numeric factorization was successful, but the matrix is singular in '" << str<< "'.\n    You'll get a devide by zero in 'umpfack_di_solve' and your solution will contain Inf's and/or NaN's.",36);
        break;
        case UMFPACK_ERROR_invalid_Symbolic_object:
        ADONIS_ERROR(UmfpackError,"Symbolic object provided as input is invalid in '"<< str<< "'.");
        break;
        case UMFPACK_ERROR_different_pattern:
        ADONIS_ERROR(UmfpackError,"The pattern (Ap and/or Ai) has changed since the call to 'umfpack_di_symbolic' which produced the symbolic object. \n   This message comes from '"<<str << "'.");
        break;
        //! solution part
        case UMFPACK_ERROR_invalid_system:
        ADONIS_ERROR(UmfpackError,"The sys argument is not valid, or the matrix A is NOT SQUARE in '"<<str<<"'.");
        break;
        case UMFPACK_ERROR_invalid_Numeric_object:
        ADONIS_ERROR(UmfpackError,"The Numeric object is not valid in '"<< str << "'.");
        break;
        default:
        ADONIS_ERROR(UmfpackError,"UNKNOWN error encountered in '"<<str << "'.");
        
     }   
	}
}
#endif
//! make copy stuff private
	SquareLinearSystem(const LinearSolver&);
    SquareLinearSystem& operator=(const LinearSolver&);
	
};
