#ifndef SOLVE_SQUARE_DENSE_OR_SPARSE_LINEAR_SYSTEMS_HH
#define SOLVE_SQUARE_DENSE_OR_SPARSE_LINEAR_SYSTEMS_HH

#include "../want_2_use.h"
#include "../common/error.hh"
#include "../common/typeadapter.hh"
#include "../common/numerictypechecker.hh"
#include "../common/error.hh"
#include "../common/universalconstants.hh"
#include "../expressiontemplates/exprvec.hh"

#include "sparsematrix.hh"
#include "sparsesettings.hh"
#include "umfdriver.hh"
#include "umftypetraits.hh"

#include "../random/randomnumbergenerators.hh"
#include "../fdm/gensettings.hh" //error handling when singularity encountered

#if USE_LAPACK
#ifdef LAPACK_USAGE
#include "../lapackwrapper/lapackdrivertraits.hh"
#include "../lapackwrapper/wrappertraits.hh"
#endif
#endif



//! you shouln't include this file but I want to try something...
#ifdef DESINGULARIZE_MATRIX
#if USE_UMFPACK
#include "../../SOFTWARE/UMFPACK/Source/umf_internal.h"

#ifndef COMPLEX
template<class T> //use double stuff
class DeSingValueSelection{
public:
  template<class NUMT, class INT>
  static inline void assign(INT i, NUMT* numeric, const T& val){
    numeric->D[i] = val;
  }
};

template<class T>  //!rule out possibility of using a complex when not complex
class DeSingValueSelection<std::complex<T> >{
public:
  template<class NUMT, class INT> //do nothing
  static inline void assign(INT i, NUMT* numeric, const std::complex<T>& val){}
};


#else //use complex stuff
template<class T> //!rule out possibility of using a double when complex
class DeSingValueSelection{
public:
  template<class NUMT, class INT> //do nothing
  static inline void assign(INT i, NUMT* numeric, const T& val){}
};

template<class T>
class DeSingValueSelection<std::complex<T> >{
public:
  template<class NUMT, class INT>
  static inline void assign(INT i, NUMT* numeric, const std::complex<T>& val){
    REAL_COMPONENT((numeric->D[i]))  = val.real();
    IMAG_COMPONENT((numeric->D[i]))  = val.imag();
  }
};
#endif

//! overload abs for Umfpack's DoubleComplex struct 
#ifdef COMPLEX
 inline double Abs(const DoubleComplex& z){
   return (sqrt(Adonis::ntimes<2>(z.component[0]) + Adonis::ntimes<2>(z.component[1])));
 }
#endif

#endif
#endif


//! error treatment in case of a singular n x n matrix
template<bool> class ErrorTreatmentOfSingularMatrix;

//!throw error when singular matrix is detected
template<>
class ErrorTreatmentOfSingularMatrix<true>{ 
public:

  enum{Value=true};


  static inline void info(const std::string& str1 = std::string(), const std::string& str2 = std::string(), int n1 = 0, int n2 = 0){
    ADONIS_ERROR(Adonis::LinearAlgebraError, "Matrix is singular up to numeric precision. Although, Umfpack's numeric factorization was successful, the matrix is singular in '" << str1 << "' (There are exact ZEROS on the diagonal of U). \n    You'll get a devide-by-zero in 'umpfack_"<< str2 << "_solve' and your solution will contain Inf's and/or NaN's. \n   Therefore, I don't see any reason why I should continue futher computations....sorry :( ");
  }

  static inline bool value() {return true;}
};

//!throw warning only  when singular matrix is detected
template<>
class ErrorTreatmentOfSingularMatrix<false>{ 
public:

  enum{Value=false};

  static inline void info(const std::string& str1 = std::string(), const std::string& str2 = std::string(), int n1 = 0, int n2 = 0){
    ADONIS_WARNING(Adonis::Warning, "Umfpack: Numeric factorization was successful, but the matrix is singular in '" << str1 << "' (There are exact ZEROS on the diagonal of U). \n    You'll get a devide-by-zero in 'umpfack_"<< str2 << "_solve' and your solution will contain Inf's and/or NaN's.");
  }

  static inline bool value() {return false;}
};



//!Use it here only, since 'Element' is also a typedef within UMFPACK
namespace Adonis{ 

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
  template<bool FLAG, class INT> class LSSFunctionalityHandler;

  template<class INT>
  class LSSFunctionalityHandler<true,INT>{ //compute pattern ONCE
  public:

#if USE_UMFPACK	
    //only compute pattern for UMFPACK once
    template<class MTX>
    static inline INT pattern_once(const MTX& mtx, const double* ax, const double* az, void** Symbolic, const double* Ctrl,  double* Info){
      //!we assume that Ax is not present (a (double*) NULL pointer). This means that any entry in A is assumed to be "large"
      //! no splitting used for complex system matrix when only computed once 
      //return UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::symbolic_analysis(mtx.rows(),mtx.cols(),mtx.position(),mtx.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::proper_iterator(mtx.values(),&ax_[0]),&az_[0],&Symbolic_,Ctrl_,Info_);
      return UmfpackDriver<typename MTX::value_type,INT>::symbolic_analysis(mtx.rows(),mtx.cols(),mtx.position(),mtx.index(),(double*) NULL,(double*) NULL,Symbolic,Ctrl,Info);
    }

	template<class MTX> //do nothing in iteration
	static inline INT pattern_in_each_iteration(const MTX& mtx, const double* ax, const double* az, void** Symbolic, const double* Ctrl, double* Info){
	return 42;
      }
	
      //do not delete 
    template<class MTX> //need the template argument here
      static inline void free_symbolic(void** Symbolic){}
	
    
    template<class MTX>
    static inline void free_symbolic_once(void** Symbolic){
      UmfpackDriver<typename MTX::value_type,INT>::free_symbolic(Symbolic);
    }

#endif //end umfpack specific stuff
      
  };

    template<class INT>
    class LSSFunctionalityHandler<false,INT>{ //in EACH iteration
    public:

#if USE_UMFPACK
      template<class MTX> 
      static inline INT pattern_once(const MTX& mtx, const double* ax, const double* az, void** Symbolic, const double* Ctrl, double* Info){
	return 42;
      }  //do nothing
	
      //!compute pattern in each iteration anew, assume that the matrix has already been split into ax_ and az_
      template<class MTX> 
      static inline INT pattern_in_each_iteration(const MTX& mtx, const double* ax, const double* az, void** Symbolic, const double* Ctrl, double* Info){
		
	return UmfpackDriver<typename MTX::value_type,INT>::symbolic_analysis(mtx.rows(),mtx.cols(),mtx.position(),mtx.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<typename MTX::value_type>::IsComplex>::proper_iterator(mtx.values(),&ax[0]),&az[0],Symbolic,Ctrl,Info);
      }
	
      template<class MTX>
      static inline void free_symbolic(void** Symbolic){
	UmfpackDriver<typename MTX::value_type,INT>::free_symbolic(Symbolic);
      }

      template<class MTX>
      static inline void free_symbolic_once(void** Symbolic){} //do nothing
#endif //end umfpack specific stuff

    };



    //////////////////////////////////////////////////////////////////////////
    /**
     * \brief Class for solving <B>square</B> linear systems \f$A\cdot x = b, A \in R^{n \times n} \ \mathrm{or}\ A \in C^{n \times n}.\f$
     * GENERAL case (default): dense systems
     * 
     * \tparam MTXTYPE matrix type \f$A\f$ on input. Is overwritten by LU 
     *                 decomposition on output.
     * \tparam RHSTYPE right hand side type \f$b\f$ on input. Is overwritten 
     *                 by solution on output.
     * \tparam CHAR matrix type. Default: 'D'ense
     * \tparam FLAG You may set something at compile time (default: false).
     *         Depends on specific application. Dummy argument for most applic.
     * \tparam INT index type of underlying matrix type (default: int)
     * Usage:
     * \code
     int dim = 4
     Matrix<double> A(dim,dim);
     //...fill A
     MyVec<double> b(dim);
     //...fill b
  
     //the last template argument is (currently) a dummy here
     SquareLinearSystem<Matrix<double>,MyVec<double>,'D',true> SQLS(A,b,1);  //construct object at beginning
   
     //... some code that updates A and b e.g. during an iterative process
     std::cout << "solution = "<< SQLS.solve() << std::endl;                 //solve system and output sol.
   
     * \endcode
     */
  template<class MTXTYPE, class VTYPE, char CHAR='D', bool FLAG=false, class INT=int>
  class SquareLinearSystem{
  public:
    typedef MTXTYPE MatrixType;
    typedef typename MTXTYPE::value_type value_type;	
    typedef typename TypeAdapter<value_type>::BaseType BaseType;	
    typedef typename VTYPE::value_type Rhs_value_type;
    typedef VTYPE VecType;
    typedef typename RandomNumberGenerator<value_type,32,'L'>::Type RandomNumberGeneratorType;
    typedef ExprTmpl::MyVec<INT> IndexVecType;
		
    static const char Value = 'D';  //the matrix ID

    SquareLinearSystem(MTXTYPE& mtx, VTYPE& b, int nrhs, INT systype = 0):A_(mtx),b_(b),status_(0),nrhs_(nrhs),dim_(int(b.size())/nrhs), ipiv_(dim_)
#if USE_LAPACK
									 ,lda_(std::max(1,dim_)),ldb_(std::max(1,dim_)),trans_('N')
#else
									 , xx_(dim_)
#endif
    {
      adonis_assert(nrhs > 0);
      adonis_assert(typeid(value_type) == typeid(Rhs_value_type)); //same value_types mandatory
      if(mtx.size() != ntimes<2>(b.size()))
	ADONIS_ERROR(DimensionError,"Dense matrix is NOT square.");
      FancyMessages().nice_output("Solve square linear system in DENSE format");
#if USE_LAPACK
      FancyMessages().nice_output("LAPACK is used.");
      if(!F77Type<TypeTraits<value_type>::Value>::IsComplex){ //real
	if(useTransposedForm)
	  trans_ = 'T';
      }
      else{                                                   //complex
	//bring array into right form for COMPLEX matrices
	  AinLapackOrder_.resize(mtx.size()); 
	  
	  if(useTransposedForm)
	    trans_ = 'N';
	}
#else
	FancyMessages().nice_output("Use Daoqi Yang's Gaussian elimination with partial pivoting (slightly adapted).");
#endif
      }

    template<class VAL>
    void ls_solver_settings(std::size_t i, const VAL& val){} //do nothing so far

    void system_to_solve(INT sys){
#if USE_LAPACK
      trans_ = ascii(sys);
#endif
    } 
    
    
    std::string system_to_solve(){
      std::string str = "System type: ";
#if USE_LAPACK
      str += trans_;
      str += " (LAPACK).";
#else
      str += "A·x = b (GaussPP).";
#endif
      return str;
    }
     
    //! do nothing right now. May be updated by you ;
    void de_singularize_matrix(){
#ifdef DESINGULARIZE_MATRIX
      for(int i = 0; i < dim_; ++i){
	A_[i*dim_ +i] = (value_type)randomgenerator_.draw_number_from_range(UniversalConstants<BaseType>::smallLow,UniversalConstants<BaseType>::smallUp); //update diagonal with nonzero random values  
      }
#endif 
    }
	
      VTYPE& solve(){
#if USE_LAPACK
	if(NumericDataTypeChecker<typename MTXTYPE::value_type>::IsComplex == true){
	  //! ONLY Necessary for COMPLEX numbers
	  //! reorder matrix array = transpose array (NO conjugate transpose)
	  transpose_array(AinLapackOrder_,A_,dim_);
	  A_ = AinLapackOrder_;
	  
	}

	//!perform LU-decomposition
	F77Type<TypeTraits<value_type>::Value>::GETRF(&dim_,&dim_,&A_[0],&lda_,&ipiv_[0],&status_);
	if(status_ != 0)
	  ADONIS_ERROR(LapackError,"An error occured within Lapack's "<<typeid(value_type).name() << "xtrf");
	//!solve square system
	F77Type<TypeTraits<value_type>::Value>::GETRS(&trans_,&dim_,&nrhs_,&A_[0],&lda_,&ipiv_[0],&b_[0],&ldb_,&status_);
#ifndef NDEBUG
	if(status_ != 0)
	  ADONIS_ERROR(LapackError,"An error occured within Lapack's "<<typeid(value_type).name() << "xtrs");
#endif
#else
	//! like in LAPACK, the A_ and b_ are overwritten by the PLU decomposition and the solution vector, respectively.
	//GAUSSIAN elimination with partial pivoting, taken from Daoqi Yang's book. It should work with 
	// complex types as well
	//MTXTYPE tmpx = A_; //local copy

	// ipiv_ contains the pivot elements now  
	for (int k = 0; k < dim_; ++k){ 
	  ipiv_[k] = k;  //overwrite ipiv_
	  xx_[k] = value_type(); //reset before each call 
	}

	int nrowsmone = dim_ - 1;
	for (int k = 0; k < nrowsmone; ++k) {  // main loop

	  // find the pivot in column k in rows pvt[k], 
	  // pvt[k+1], ..., pvt[n-1]
	  int pc  = k; 
	  BaseType aet = Abs(A_[k*dim_ + k]);//(tmpx[pvt[k]][k]);
	  for (int i = k + 1; i < dim_; ++i) {
	    if (Abs(A_[ipiv_[i]*dim_ + k]) > aet){ //[pvt[i]][k]) > aet) {
	      aet = Abs(A_[ipiv_[i]*dim_+k]);//tmpx[pvt[i]][k]); 
	      pc = i;
	    }
	  }
	  if (is_zero(aet)){ 
	    ADONIS_WARNING(ZeroDivision, "Pivot is zero in Gaussian elimination with partial pivoting.\n   ==> matrix is SINGULAR or badly scaled.");
	    (*this).de_singularize_matrix();
	  } 
	  if (pc != k) std::swap(ipiv_[k], ipiv_[pc]);  //defined in algorithm
	  int pvtk = ipiv_[k];                  // pivot row
	  value_type pivot = A_[pvtk*dim_ + k]; //tmpx[pvtk][k];            // pivot

	  // now eliminate column entries logically 
	  // below tmpx[pvt[k]][k]
	  value_type mult;
	  for (int i = k + 1; i < dim_; ++i) {
	    int pvti = ipiv_[i];
	    if (!is_zero(A_[pvti*dim_ + k])){ //tmpx[pvti][k] != 0) {
	      mult = A_[pvti*dim_ +k]/pivot;  //tmpx[pvti][k]/pivot;
	      A_[pvti*dim_ + k] = mult;  //tmpx[pvti][k] = mult;
	      for (int j = k + 1; j < dim_; ++j) 
		A_[pvti*dim_+j] -= mult*A_[pvtk*dim_+j];//tmpx[pvti][j] -= mult*tmpx[pvtk][j];
	    }
	  }
	}

	//!sovle system after LU decomposition
	// forwad substitution for L y = Pb.
	for (int i = 1; i < dim_; ++i)  
	  for (int j = 0; j < i; ++j) 
	    b_[ipiv_[i]] -= (A_[ipiv_[i]*dim_+j])*b_[ipiv_[j]];//tmpx[pvt[i]][j]*bb[pvt[j]];

	// back substitution for Ux = y
	for (int i = nrowsmone; i >= 0; --i) {
	  for (int j = i+1; j < dim_; ++j) 
	    b_[ipiv_[i]] -=  (A_[ipiv_[i]*dim_ + j]*xx_[j]); //tmpx[pvt[i]][j]*xx[j];
	  xx_[i] = b_[ipiv_[i]] / (A_[ipiv_[i]*dim_ + i]);//tmpx[pvt[i]][i];
	}

	b_ = xx_;             // put solution
#endif
	return b_;
      }
	
	
    bool is_singular() {return false;}

	private:
	    MTXTYPE& A_;  //use references here since we want to use LAPACK here in the first place
	  VTYPE& b_;
	  int status_, 
	    nrhs_, dim_;
	  IndexVecType ipiv_;

#ifdef DESINGULARIZE_MATRIX
    RandomNumberGeneratorType randomgenerator_; //seeded with time(0)
#endif

#if USE_LAPACK
	  int lda_,ldb_;
	  char trans_;     
	  //Lapack is column-oriented, hence we use the transposed form (C-style)
    MTXTYPE AinLapackOrder_;
	  enum{useTransposedForm = true}; 
#else 
	  //! use Gaussian elimination with partial pivoting
	  VTYPE xx_; //stores solution in correct otder 
#endif

	  //! make copy stuff private
	  SquareLinearSystem(const SquareLinearSystem&);
	  SquareLinearSystem& operator=(const SquareLinearSystem&);
      };


  //!==========================================================================
  //!============== partial specializations for SPARSE square systems =========
  //!==========================================================================
	/**
	 * \brief Solves <B>sparse</B> square linear system \f$A\cdot x = b,\f$ where \f$A\f$ is given in 
	 * compressed sparse column format (the UMFPACK-style). 
	 * \tparam MTXTYPE compressed sparse column object
	 * \tparam VTYPE right hand side vector type
	 * \tparam FLAG if true, the pattern is computed only once, otherwise in each call to the 'solve' routine
	 */
	template<class MTXTYPE, class VTYPE, bool FLAG, class INT>
	  class SquareLinearSystem<MTXTYPE,VTYPE,'C',FLAG,INT>{
#if USE_UMFPACK
	  typedef UmfpackDriver<typename MTXTYPE::value_type,INT> DriverType;
#endif
	public:
	  typedef MTXTYPE MatrixType;
	  typedef typename MTXTYPE::value_type value_type;	
	  typedef typename TypeAdapter<value_type>::BaseType BaseType;	
	  typedef typename VTYPE::value_type Rhs_value_type;
	  typedef VTYPE VecType;

	  typedef typename RandomNumberGenerator<value_type,32,'L'>::Type RandomNumberGeneratorType;
	  
	  //when Complex linear algebra is applied, we need two additional vectors storing the real and imaginary parts of the complex coefficients
	  typedef ExprTmpl::MyVec<INT> IndexVecType;
	  typedef ExprTmpl::MyVec<BaseType> BaseVecType;

	
	  static const char Value = 'C';  //the matrix ID

	  SquareLinearSystem(MTXTYPE& mtx, VTYPE& b, int nrhs = 1, INT systype = 0):A_(mtx),b_(b),status_(0),nrhs_(nrhs), isSing_(false),numOfSings_(0)
#if USE_UMFPACK
										   ,sys_(check_system_type(systype)) // Ax = b , UMFPACK_At // A'x = b (transpose/cplx conj)
								  ,ax_(1),az_(1),  //do set to 1 in any case -- just two scalars
								   xx_(b.size()),xz_(1),bx_(b.size()),bz_(1)
								  
#endif
	  {
	    adonis_assert(nrhs > 0);
	    if(mtx.rows() != mtx.cols())
	      ADONIS_ERROR(DimensionError,"Sparse matrix is NOT square.");
	    adonis_assert(typeid(typename MTXTYPE::IndexType) == typeid(INT));
	    adonis_assert((MTXTYPE::Value == 'C') || (MTXTYPE::Value == 'c'));
	    adonis_assert(typeid(value_type) == typeid(Rhs_value_type)); //same value_types are mandatory
	    std::string str = ((FLAG==true) ? "ONCE." : "in EACH call to 'solve'.");
	    FancyMessages().nice_output("Solve square linear system in compressed SPARSE column (CSC) format.\nPattern computed " + str);
#if USE_UMFPACK
	    FancyMessages().nice_output("UMFPACK is used.");
	    FancyMessages().nice_output(" System type: "+(*this).system_info()+".");
	    
	    
	      

	    adonis_assert(nrhs == 1); //only this seems to be possible within Umfpack at present
	    if(mtx.rows() != (INT)b.size())
	      ADONIS_ERROR(DimensionError,"Sparse matrix order and RHS vector size do not match.");
	    //! default settings MUST be loaded into Control_ array
	    DriverType::defaults(Control_); 

	    Wi_.resize(mtx.rows());
	    W_.resize(DriverType::w_dim(mtx.rows(),sys_,Control_,Info_));

	    if(mtx.nonz() != 0){
	      //only resize when size is not equal to nonz and if complex systems are to be solved!
	      UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::resize(mtx.nonz(),ax_);
	      UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::resize(mtx.nonz(),az_);

	      UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::resize(b.size(),xz_);
	      UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::resize(b.size(),bz_);
	    }
	    //! see <a> http://www.math.umbc.edu/~rouben/2003-09-math625/umfpack-ex3.c</a>
	    //*Symbolic_ = NULL;
	    //*Numceric_ = NULL;
	    //====TODO
	    status_ =  LSSFunctionalityHandler<FLAG,INT>::pattern_once(A_,(double*) NULL,(double*) NULL,&Symbolic_,Control_,Info_);
	    status_umfpack(status_,0); //check status -- throw warning or even error when UMFPACK_OK status isn't met
#else
	    ADONIS_ERROR(ImplementationError,"Sorry, NO alternative to solve sparse systems based on the CSC format\n   has been specified so far! Feel free to do so ;) ");
#endif	
	  } //end construction
		
	  //destruction
	  ~SquareLinearSystem(){
#if USE_UMFPACK
	    //! clear me- if not delete yet
	    LSSFunctionalityHandler<FLAG,INT>::template free_symbolic_once<MTXTYPE>(&Symbolic_);
	    // if(&Symbolic_ != NULL)
	    //   DriverType::free_symbolic(&Symbolic_);
	    // if(&Numeric_ != NULL)
	    //   DriverType::free_numeric(&Numeric_);  //??????????
#endif
	  }
	 
	  //! you may set some solver specific features with this member
	  template<class VAL>
	  void ls_solver_settings(std::size_t i, const VAL& val){
#if USE_UMFPACK
	    //! cf. UMFPACK UserGuide, p. 68: individual settings can be 
	    //! modified by changing specific entries in the Control_ array. 
	    //! If Control_ = 0, use default control settings (these are 
	    //! suitable for all matrices, ranging from those with symmetric 
	    //! to highly unsymmetric patterns. See p. 22 for available 
	    //! parameters and their default settings. For instance, 
	    //! (0,2) sets umfpack printlevel to 2 (only useful with 
	    //! umfpack_*_report_* routines, cf. pp. 115. 
	    adonis_assert(i < (std::size_t)UMFPACK_CONTROL);
	    //adonis_assert(typeid(VAL) == typeid(double));
	    Control_[i] = (double)val;    
#endif
	  }

	  void system_to_solve(INT sys){
#if USE_UMFPACK
	    //! see Umfpack documentation, p. 59: UMFPACK_A = 0, UMFPACK_At = 1
	    //! UMFPACK_Aat = 2, etc.
	    sys_ = sys;
#endif
	  }

	  std::string system_to_solve(){
	    std::string str;
#if USE_UMFPACK
	    str = "System type (UMFPACK): "+ system_info() + ".";
#else
	    str = "A·x = b.";  //default
#endif
	    return str;
	  }

	  void de_singularize_matrix(){
#ifdef DESINGULARIZE_MATRIX
#if USE_UMFPACK
	    //! only when matrix is detected to be singular
	    //! use this then between umfpack_*_numeric and umfpack_*_wsolve
	    if(status_ == UMFPACK_WARNING_singular_matrix){
	      //do not allow too many desingularizations...
	      numOfSings_++;  //increment
	      if(numOfSings_ > 10) //if we have shitted too much, throw error
	      ADONIS_ERROR(DerivedError, "Too many de-singularizations performed. I think the result will be spoiled too much \n   Try to repose your problem in order to avoid singular linear systems!");
	      
	      std::cout << std::endl;
	      FancyMessages().nice_output("***** Trying to de-singularize matrix. This may result in nothing better...",35);

	      NumericType* NumT = (NumericType*)Numeric_;
	      //value_type mean = value_type();
	      BaseType min_d(0), // = Abs(NumT->D[0]), 
		max_d(0); // = Abs(NumT->D[0]); //min and max abs values of diagonal of U
	      INT count(0);
	      for(INT i = 0; i < A_.rows(); ++i){
		//std::cout << "NumT->D["<<i<<"] = "<< (NumT->D[i]) << std::endl;
		//! consider only good values on diagonal of U here
		if( (is_well_defined_value(NumT->D[i])) && (!is_zero(NumT->D[i])) ){
		  if(count == 0){ //min/max value is first proper value on diag
		    min_d = Abs(NumT->D[i]);
		    max_d = Abs(NumT->D[i]);
		  }

		  min_d = Min(min_d,Abs(NumT->D[i]));
		  max_d = Max(max_d,Abs(NumT->D[i]));
		  //std::cout << count << ".)   min_d = "<< min_d << "   max_d = "<< max_d << std::endl;
		  
		  //mean += (NumT->D[i]); //can be a complex number
		  count++;
		}
		
	      }
	      if(count==0)
		ADONIS_ERROR(FatalError,"Columns of matrix are ALL linearly DEPENDENT. This is really serious");
	      else{
		//mean /= count;
	      }
	      // std::cout << "mean = "<< mean << std::endl;
	      //!assign value to "bad" values on diagonal of U
	      for(INT i = 0; i < A_.rows(); ++i){
		if( (is_well_defined_value(NumT->D[i])==false) || (is_zero(NumT->D[i])) ){
		 
		  //put small random values on diagonal
		  DeSingValueSelection<value_type>::assign(i,NumT,(BaseType)(Sgn((NumT->D[i])))*((NumT->D[i])+randomgenerator_.draw_number_from_range(UniversalConstants<BaseType>::smallLow,UniversalConstants<BaseType>::smallUp))); //numbers from range not considered zero
     
		}
	      }
	      NumT->nnzpiv = A_.rows();
	      NumT->min_udiag = min_d;
	      NumT->max_udiag = max_d;
	     
	      //NumT->Upattern[A_.rows()-1] = 1;
	      //std::cout << std::setprecision(5) << "min_udiag = "<< min_d << "   max_udiag = "<< max_d << std::endl;
	      
	      NumT->rcond = min_d/max_d;
	      printf("%c[%d;%d;%dm", 0x1B, 1,  35, 40);
	      //std::cout << "***** bad value substitution = "<< mean << std::endl;
	      std::cout << "***** min_diag = "<< min_d << ",  max_diag = "<<max_d << std::endl;
	      std::cout << "***** RCOND_new = "<< std::setprecision(15) << NumT->rcond  << ".  Estimated rank = "<< count  << " (full rank: "<< A_.rows() << ")." << std::endl;
	      std::cout << "***** Resume further calculations..."<<std::endl;
	      printf("%c[%dm", 0x1B, 0);

	      if(is_zero(NumT->rcond)){
		//ADONIS_ERROR(BadConditioningError,"RCOND_new is still approx. ZERO.");
		 printf("%c[%d;%d;%dm", 0x1B, 1,  33, 40);
		 std::cout << "WARNING:  RCOND_new is still approx. ZERO, making the LHS matrix close to singular." << std::endl;
		 printf("%c[%dm", 0x1B, 0);
	      }

	      Numeric_ = (void*)NumT; //assign back (necessary?)

	      status_ = UMFPACK_OK; //assume nonsingular matrix has been created
	    }
#endif	 
#endif
	  }

	  VTYPE& solve(){
	    isSing_ = false; // by default system is considered regular
#if USE_UMFPACK
	    //! split system matrix when complex numbers are used
	    UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::split_complex_container_into_real_values(A_.nonz(),A_.values(),&ax_[0],&az_[0]);
	    // std::cout << "ax_ = "<< ax_ << std::endl;
	    //std::cout << "az_ = "<< az_ << std::endl;

	    //! split rhs when complex numbers are used
	    UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::split_complex_container_into_real_values(b_.size(),&b_[0],&bx_[0],&bz_[0]);		
	    //std::cout << "bx_ = "<< bx_ << std::endl;
	    //std::cout << "bz_ = "<< bz_ << std::endl;
	    
	    //! perform symbolic operation in each call if not computed at the very beginning
	    status_ = LSSFunctionalityHandler<FLAG,INT>::pattern_in_each_iteration(A_,&ax_[0],&az_[0],&Symbolic_,Control_,Info_);			
	    //status_ = DriverType::symbolic_analysis(A_.rows(),A_.cols(),A_.position(),A_.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::proper_iterator(A_.values(),&ax_[0]),&az_[0],&Symbolic_,Control_,Info_);
	    status_umfpack(status_,0); //everything o.k.?
			
	    status_ =  DriverType::lu_decomposition(A_.position(),A_.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::proper_iterator(A_.values(),&ax_[0]),&az_[0],Symbolic_,&Numeric_,Control_,Info_);  //umfpack_di_numeric(A_.position(),A_.index(), A_.values(),Symbolic_,Numeric_,Control_,Info_);
	    
	    status_umfpack(status_,1); //everything o.k.?
	    
	    (*this).de_singularize_matrix();

	    //! wsolve is faster for iterative systems since space allocations 
	    //! are done once and not in each call
	    //! Numeric_ is not modified in the solution process
	    status_ = DriverType::wsolve(sys_,A_.position(),A_.index(),UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::proper_iterator(A_.values(),&ax_[0]),&az_[0],&xx_[0],&xz_[0],UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::proper_iterator(&b_[0],&bx_[0]),&bz_[0],Numeric_,Control_,Info_,&Wi_[0],&W_[0]);

	    // NumericType* NumT = (NumericType*)Numeric_; 
	    // std::cout << "------------------------------------" << std::endl;
	    // std::cout << "nnzpiv = "<< NumT->nnzpiv<< std::endl;
	    // std::cout << "rcond  = "<< NumT->rcond << std::endl;
	    // std::cout << "------------------------------------" << std::endl;

	    status_umfpack(status_,2);//umfpack_di_solve(UMFPACK_A,A_.position(),A_.index(),A_.values(),solution.begin(),Numeric_,Control_,Info_);
	    
	    UmfpackUsesComplexNumbers<NumericDataTypeChecker<value_type>::IsComplex>::join((int)b_.size()*nrhs_,&b_[0],&xx_[0],&xz_[0]);
		
      
	    LSSFunctionalityHandler<FLAG,INT>::template free_symbolic<MTXTYPE>(&Symbolic_); //only delete it when FLAG = false
	    DriverType::free_numeric(&Numeric_);  //delete LU decomposition object	
#else
		
#endif
	    return b_;
	  }
	

	  bool is_singular() {return isSing_;}

	private:
	  MTXTYPE& A_;  //references since A and b are updated during iteration. They remain unaffected by UMFPACK
	  VTYPE& b_;
	  INT status_; //status of the solver, depends of the specific implementation
	  int  nrhs_;       //number of RHSs; UMFPACK can currently not handle multiple right hand sides
	  bool isSing_;        //matrix singular
	  int numOfSings_;     //number of singularities encountered

#if USE_UMFPACK
	  INT sys_;  //!tells UMFPACK which system is going to be solved 
	  //!That's Umfpack's C-style: a void pointer can point to <I>any</I> object.
	 
	  //! both for double an complex routines, these arrays are of type double

	  //additional storage for complex square systems
	  BaseVecType ax_, az_,
	    xx_, xz_,
	    bx_, bz_;
	  
	  IndexVecType Wi_;
	  BaseVecType W_;

	  //double* Control_; 
	  double Control_[UMFPACK_CONTROL]; //if this is set, you must fill the array with values before proceding; cf. UMFPACK documentation
          double Info_[UMFPACK_INFO];
	 

	  //! In C++, the use of templates would be more convenient
	  void* Symbolic_;
	  void* Numeric_;

#ifdef DESINGULARIZE_MATRIX
	  RandomNumberGeneratorType randomgenerator_;
#endif

	  typedef ErrorTreatmentOfSingularMatrix<FDMSettings::throwErrorWhenSingularMatrixDetected> ErrorSingularityTreatmentType;

	  //! private functions
	  //make this private and only usable when UMFPACK usage is switched on
	  void status_umfpack(const INT status, int kind){
	    if(status != UMFPACK_OK){  //UMFPACK reports an error. Give notice about kind of error
	      std::string str = "umfpack_"+DriverType::routine_type()+"_";
	      if(kind == 0)
		str += "symbolic";
	      else if (kind == 1)
		str += "numeric";
	      else if (kind == 2)
		str += "solve";
	      else
		str += "???";
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
		break;
	      case UMFPACK_ERROR_internal_error:
		  ADONIS_ERROR(UmfpackError,"Something very serious went wrong. This is a BUG.\n    Please contact DrTimothyAldenDavis@gmail.com");
		break;
		//! numeric part
	      case UMFPACK_WARNING_singular_matrix:
		isSing_ = true;
		FancyMessages().nice_output("********* UMFPACK: SINGULAR MATRIX ENCOUNTERED...",33);

		if(kind>=1){  //only output for Numeric
		  //! copy information from the Numeric object
		  INT lnz, unz, n_row, n_col,nz_udiag;
		  DriverType::get_lunz(&lnz,&unz,&n_row,&n_col,&nz_udiag,Numeric_);
		  //select which part of the Numeric object you wanna to see
                  ExprTmpl::MyVec<double> Dx(n_row), Dz;
		  if(NumericDataTypeChecker<typename MTXTYPE::value_type>::IsComplex == true){
		    Dz.resize(n_row);
		    
		  }
		  
		  

		  DriverType::get_numeric((INT*)NULL,(INT*)NULL,(double*)NULL,(double*)NULL,(INT*)NULL,(INT*)NULL,(double*)NULL,(double*)NULL,(INT*)NULL,(INT*)NULL,&Dx[0],((NumericDataTypeChecker<typename MTXTYPE::value_type>::IsComplex == true) ? &Dz[0] : (double*)NULL),(INT*)NULL,(double*)NULL,Numeric_);
		  printf("%c[%d;%d;%dm", 0x1B, 1,  33, 40);
		  std::cout << " *****CAUTION: Bad diagonal entries of U detected: " << std::endl; //contains zeros
#ifdef SHOW_BAD_DIAGONAL_ENTRIES //defined in sparsesettings.hh
		  std::cout << " ----------------------------------- "<<std::endl;
#endif
		  BaseType min_d,
		    max_d;
		  std::size_t ct = 0;
		   if(NumericDataTypeChecker<typename MTXTYPE::value_type>::IsComplex == false){ //real case
		      min_d = Abs(Dx[0]);
		      max_d = Abs(Dx[0]);
		   }
		   else{//complex case
		     min_d = sqrt(ntimes<2>(Dx[0])+ntimes<2>(Dz[0]));
		     max_d = sqrt(ntimes<2>(Dx[0])+ntimes<2>(Dz[0]));
		   }

		  for(INT i = 0; i < n_row; ++i){
		    if(NumericDataTypeChecker<typename MTXTYPE::value_type>::IsComplex == false){ //real case
		      min_d = Min(min_d,Abs(Dx[i]));
		      max_d = Max(max_d,Abs(Dx[i]));

		      if(!is_well_defined_value(Dx[i]) || is_zero(Dx[i])){
#ifdef SHOW_BAD_DIAGONAL_ENTRIES
			std::cout << "D["<<i<<"] = "<<Dx[i] << "  ";
#endif
			ct++;
		      }

		    }
		    else{ //complex
		      min_d = Min(min_d,sqrt(ntimes<2>(Dx[i])+ntimes<2>(Dz[i])));
		      max_d = Max(max_d,sqrt(ntimes<2>(Dx[i])+ntimes<2>(Dz[i])));

		      if((!is_well_defined_value(Dx[i])) || (is_zero(Dx[i]) && (!is_well_defined_value(Dz[i]) || is_zero(Dz[i])))){
#ifdef SHOW_BAD_DIAGONAL_ENTRIES
			std::cout << "D["<<i<<"] = ("<<Dx[i] << ",  "<<Dz[i] << ")  ";
#endif
			ct++;
		      }
		    }
		    
		  }
#ifdef SHOW_BAD_DIAGONAL_ENTRIES
		  std::cout << std::endl;
		  std::cout << " ----------------------------------- "<<std::endl;
#endif
		  //rcond is zero/NaN if there is any zero or NaN on the diagonal
		  std::cout << "  There are " << ct << " bad values on D (out of "<< A_.rows() << " in total)."<<std::endl;
		  std::cout << "  RCOND = "<< min_d/max_d << std::endl;
		  printf("%c[%dm", 0x1B, 0);
		}
		//either warning or error as it was choosen in fdm/gensettings.hh
		ErrorSingularityTreatmentType::info(str,DriverType::routine_type());

		
	      break;
	      case UMFPACK_ERROR_invalid_Symbolic_object:
		ADONIS_ERROR(UmfpackError,"Symbolic object provided as input is invalid in '"<< str<< "'.");
		break;
	      case UMFPACK_ERROR_different_pattern:
		ADONIS_ERROR(UmfpackError,"The pattern (Ap and/or Ai) has changed since the call to 'umfpack_"+ DriverType::routine_type() +"_symbolic' which produced the symbolic object. \n   This message comes from '"<<str << "'.");
		break;
		//! solution part
	      case UMFPACK_ERROR_invalid_system:
		ADONIS_ERROR(UmfpackError,"The sys argument is not valid, or the matrix A is NOT SQUARE in '"<<str<<"'.");
		break;
	      case UMFPACK_ERROR_invalid_Numeric_object:
		ADONIS_ERROR(UmfpackError,"The Numeric object is not valid in '"<< str << "'.");
		break;
	      default:
		ADONIS_ERROR(UmfpackError,"UNKNOWN error encountered in '"<<str << "', probably because some data were unintialized.");
        
	      }   
	    }
	  }

	  const INT check_system_type(const INT sys) const{
	    if((sys < 0) || (sys > 14))
	      ADONIS_ERROR(ValueError, "System type "<< sys << " is not defined  within UMFPACK.\n   Check \"umfpack.h\" for allowed encodings of types of systems to be solved." );
	    return sys;
	  }

	  std::string system_info() const{
	    std::string systeminfo;
	    switch(sys_){
	    case UMFPACK_A:             //0
	      systeminfo = "Ax=b";
	      break;
	    case UMFPACK_At:            //1
	      systeminfo = "A'x=b";
	      break;
	    case UMFPACK_Aat:           //2
	      systeminfo = "A.'x=b, where (.') is the array transpose";
	      break;
	    case UMFPACK_Pt_L:          //3
	      systeminfo = "P'Lx=b";
	      break;
	    case UMFPACK_L:             //4
	      systeminfo = "Lx=b";
	      break;
	    case UMFPACK_Lt_P:          //5
	      systeminfo = "L'Px=b";
	      break;
	    case UMFPACK_Lat_P:         //6
	      systeminfo = "L.'Px=b, where (.') is the array transpose";
	      break;
	    case UMFPACK_Lt:            //7
	      systeminfo = "L'x=b";
	      break; 
	    case UMFPACK_Lat:           //8
	      systeminfo = "L.'x=b, where (.') is the array transpose";
	      break;
	    case UMFPACK_U_Qt:          //9
	      systeminfo = "UQ'x=b";
	      break;
	    case UMFPACK_U:             //10
	      systeminfo = "Ux=b";
	      break;
	    case UMFPACK_Q_Ut:          //11
	      systeminfo = "QU'x=b";
	      break;
	    case UMFPACK_Q_Uat:         //12
	      systeminfo = "QU.'x=b, where (.') is the array transpose";
	      break;
	    case UMFPACK_Ut:            //13
	      systeminfo = "U'x=b";
	      break;
	    case UMFPACK_Uat:           //14
	      systeminfo = "U.'x=b, where (.') is the array transpose";
	      break;
	    default:
	      systeminfo = "Type of linear system undefined so far";

	    }
	    return systeminfo;
	  }
#endif //end UMFPACK 
	  //! make copy stuff private
	  SquareLinearSystem(const SquareLinearSystem&);
	  SquareLinearSystem& operator=(const SquareLinearSystem&);
	
	};

	

} //end namespace 

#endif
