#ifndef ERROR_CLASS_USING_MACROS_HH
#define ERROR_CLASS_USING_MACROS_HH

#include <iostream>

#include<string>
#include<fstream>
#include<sstream>
#include<cstdlib>  //for function 'exit' 

#include "fancymessages.hh"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define ATTHISPLACELETMEKNOW __FILE__ ":"  TOSTRING(__LINE__) ": "


namespace Adonis{

//you may apply the following...
//USAGE: error(ATTHISPLACELETMEKNOW, "some text describing  the error");
  inline void error(const std::string& location, const std::string& msg){
    std::cout<< "Error at "<< location  <<  msg<< std::endl;
  }


  //...or you can also use my error class just below

  
  class Error{
  private:
    std::string s_;
  
  public:
    //also default destructor
    Error(const std::string& s = std::string()):s_(s){} 

    void notes(const std::string& s){s_ = s; }

    const std::string& notice() const {return s_;}
    
    void gimme() const;
  
    void info() const;

    void leave_programme() const;
  
    void warning() const;

    void caution() const;

    void create_error_log_file() const {
      FancyMessages FM;
      FM.error_log(s_);
    }
  };

  /**
   *  \brief This macro wraps the following code
   *  \code
   *  FancyMessages FM;
   *  FM.nice_error(s_);
   * \endcode
   */
#define WRAP_FANCY_MESSAGE(Obj,message) do{Obj O; \
    O.nice_error(message);		  \
  }while(0)

/**
   *  \brief This macro wraps the following code
   *  \code
      FancyMessages FM;
      FM.nice_info(s_);
   * \endcode
   */
#define WRAP_FANCY_INFO(Obj,message) do{Obj O; \
    O.nice_info(message);		  \
  }while(0)


  /**
   *  \brief This macro wraps the following code
   *  \code
      FancyMessages FM;
      FM.nice_caution(s_);
   * \endcode
   */
#define WRAP_FANCY_CAUTION(Obj,message) do{Obj O; \
    O.nice_caution(message);		  \
  }while(0)


  /**
   *  \brief This macro wraps the following code
   *  \code
      FancyMessages FM;
      FM.nice_warning(s_);
   * \endcode
   */
#define WRAP_FANCY_WARNING(Obj,message) do{Obj O; \
    O.nice_warning(message);		  \
  }while(0)

  class FancyNancy: public FancyMessages{};


  inline void Error::gimme() const{
    
    WRAP_FANCY_MESSAGE(FancyNancy,s_);
  }


  inline void Error::info() const{
    WRAP_FANCY_INFO(FancyNancy,s_);
  }

  inline void Error::caution() const{
    WRAP_FANCY_CAUTION(FancyNancy,s_);
  }

  inline void Error::warning() const{
    WRAP_FANCY_WARNING(FancyNancy,s_);
  }

  inline void Error::leave_programme() const{
    (*this).gimme();
    (*this).create_error_log_file(); //write error file
    abort();   //terminate process
  }


  //non-member function
  inline std::ostream& operator<<(std::ostream &os, const Error& e){
    return os << e.notice();
  }




  /*  
   *    @file
   *   Learn how to use parameterized Macros on
   *  
   *    <a href="http://gcc.gnu.org/onlinedocs/gcc-2.95.3/cpp_1.html#SEC16"> Click here</a>
   *
   *  USAGE, e.g.: 
   * \code  
     try{
          .   
	  .
	 }
	 catch(Error& e){
            std::cerr << "Your error specification."<<std::endl;
	    return 1;
	 }
	 catch(...){std::cerr << "Unspecified error occurred."<<std::endl}
   *	 \endcode
  */
  #define SPECIFY_ERROR_LOCATION "In file:'" << __FILE__<<"':\n   | \n   In function: `" << __func__ << "':\n       | \n       At line: " << __LINE__ << " \n \n Description: \n ------------ \n   " 


//another Macro which stops all subsequent calculation on an instant 
#define ADONIS_ERROR(Err, say) do{Err e; std::ostringstream os;		\
    os << "Type: " << #Err << std::endl<< SPECIFY_ERROR_LOCATION << say; e.notes(os.str());\
 e.leave_programme();  \
  } while(0)   


#define ADONIS_WARNING(Err, say) do{Err e; std::ostringstream os;		\
    os << "Type: " << #Err << std::endl<< SPECIFY_ERROR_LOCATION << say; e.notes(os.str());\
 e.warning();  \
  } while(0)   

#define ADONIS_INFO(Err, say) do{Err e; std::ostringstream os;		\
    os << "Type: " << #Err << std::endl<< SPECIFY_ERROR_LOCATION << say; e.notes(os.str());\
 e.info();  \
  } while(0)   


#define ADONIS_CAUTION(Err, say) do{Err e; std::ostringstream os;		\
    os << "Type: " << #Err << std::endl<< SPECIFY_ERROR_LOCATION << say; e.notes(os.str());\
 e.caution();  \
  } while(0)   


  //derived default classes defining specialised errors 
  class BadConditioningError: public Error{};

  class BoundsError: public Error{};

  class CoverageError: public Error{};
  
  class DerivedError: public Error{};

  class ClearError: public Error{};
  
  class ComputationError: public Error{};

  class ConfigurationError: public Error{};

  class DefinitionError:  public Error{};
  
  class DimensionError: public Error{};

  class DerivativeError: public Error{};

  class IndexError: public Error{};
  
  class NotSupportedAnyMoreError: public Error{};

  class CancellationError: public Error{};

  class ZeroDivision: public ComputationError{};

  class LapackError: public Error{};

  class OutOfRange: public Error{};

  class MassBalanceError: public Error{};

  class MemoryError: public Error{};

  class IndexRangeViolated: public Error{};

  class IllegalParameters: public Error{};

  class IllegalInitialization: public Error{};

  class ImplementationError: public Error{};

  class IterationError: public Error{};

  class IncompleteOrMissingDataError: public Error{};

  class NothingHasBeenDeclaredYet: public Error{};

  class IOError: public Error{};

  class IODimensionError: public IOError{};
  
  class FatalError: public Error{};

  class FileError: public IOError{};

  class GraphicsError: public Error{};

  class GnuplotError: public GraphicsError{};
  
  class LinearAlgebraError: public Error{};

  class MainProgramError: public Error{};
  
  class MaxIterationError: public Error{};

  class NotImplementedError: public Error{};

  class NotYetDefinedError: public Error{};

  class PositivityError: public Error{};

  class BewareOfTooMuchEffortsError: public Error{};

  class RangeError: public Error{};

  class UninitializedValuesError: public Error{};
  
  class UmfpackError: public Error{};

  class UndefinedCaseError: public Error{};
  
  class UndefinedFunctionError: public Error{};

  class InfeasibleError: public Error{};

  class SegmentationFaultError: public Error{};
  
  class SettingsError: public Error{};

  class TypeError: public Error{};

  class ValueError: public Error{};

  class VariantError: public Error{}; 
  

  //!a warning
  class Warning: public Error{}; 

  //!can be used for general information
  class Information: public Error{};

  //! caution
  class Caution: public Warning{};


  /**
   * \brief TMP conditional whether you throw an error or, e.g. a warning instead
   * \tparam B decide whether error is thrown (true) or not (false)
   * \tparem ERR error class from above, e.g. ZeroDivision, Warning, Caution, etc.
   */
  template<bool B, class ERRCLSTYPE> class ThrowErrorOrNot;

  //! throw error assuming that any further computation leads to spurious or unreasonable results
  template<class ERRCLSTYPE>
  class ThrowErrorOrNot<true,ERRCLSTYPE>{
  public:
    enum{Value=true};
    
    static inline void info(const std::string& str1 = std::string(), const std::string& str2 = std::string(), const std::string& str3 = std::string() ){
      ADONIS_ERROR(ERRCLSTYPE, "*****start error: \n "<< str1 << "\n " << str2 << "\n "<< str3 << "\n end error *****");
    }

    static inline bool value() {return Value;}
  };

  // throw warning only. 
  template<class ERRCLSTYPE>
  class ThrowErrorOrNot<false,ERRCLSTYPE>{
  public:
    enum{Value=false};
    
    static inline void info(const std::string& str1 = std::string(), const std::string& str2 = std::string(), const std::string& str3 = std::string() ){
      ADONIS_WARNING(ERRCLSTYPE, "*****start warning: \n "<< str1 << "\n " << str2 << "\n "<< str3 << "\n end warning *****");
    }

    static inline bool value() {return Value;}
  };
  
} //end of namespace



#endif
