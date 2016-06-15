#ifndef ALLOWED_NUMERIC_TYPES_HH
#define ALLOWED_NUMERIC_TYPES_HH

#include<complex>
#include<string>
#include<typeinfo>
#include "error.hh"

namespace Adonis{

  /**
   * \brief Virtually any numeric method is based on floating point arithmetic. Hence we can check at compile time whether a type is 'non-numeric', such as int, bool, long int, etc. to avoid stupid results. 
   */
  template<class T> class NumericDataTypeChecker{ //! default
  public:
    static const bool IsComplex = false;
    static const bool IsIntegralType = false;
    static inline void certify(){} //do nothing -- accepted types
  };

  template<class T> 
  class NumericDataTypeChecker<std::complex<T> >{
  public:
    static const bool IsComplex = true;
    static const bool IsIntegralType = false;
    static inline void certify(){} //do nothing -- accepted types
  };

  //! more specialisations
  template<>
  class NumericDataTypeChecker<short>{ //the same as 'short int'
  public:
    static const bool IsComplex = false;
    static const bool IsIntegralType = true;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(short int).name() <<"\" is not a numeric data type");
    }
  };

  
  template<>
  class NumericDataTypeChecker<int>{ //same a as signed
  public:
    static const bool IsComplex = false;
    static const bool IsIntegralType = true;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(int).name() <<"\" is not a numeric data type");
    }
  };

  template<>
  class NumericDataTypeChecker<unsigned>{  //same as 'unsigned int' and 'size_t'
  public:
    static const bool IsComplex = false;
    static const bool IsIntegralType = true;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(unsigned).name() <<"\" is not a numeric data type");
    }
  };
  

  template<>
  class NumericDataTypeChecker<long>{ //same as 'long int', on 64 bit: 'size_t'
  public:
    static const bool IsComplex = false;
    static const bool IsIntegralType = true;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(long int).name() <<"\" is not a numeric data type");
    }
  };


  template<>
  class NumericDataTypeChecker<bool>{
  public:
    static const bool IsComplex = false;
    static const bool IsIntegralType = true;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(bool).name() <<"\" is not a numeric data type");
    }
  };
  
  template<>
  class NumericDataTypeChecker<char>{
  public:
    static const bool IsComplex = false;
    static const bool IsIntegralType = true;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(char).name() <<"\" is not a numeric data type");
    }
  };
  
  template<>
  class NumericDataTypeChecker<wchar_t>{ //wide char
  public:
    static const bool IsComplex = false;
    static const bool IsIntegralType = true;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(wchar_t).name() <<"\" is not a numeric data type");
    }
  };

  
  //the same holds for COMPLEX numbers seeded with integral types
  template<>
  class NumericDataTypeChecker<std::complex<short> >{ //the same as 'short int'
  public:
    static const bool IsComplex = true;
    static const bool IsIntegralType = false;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(std::complex<short>).name() <<"\" is not a numeric data type");
    }
  };

  
  template<>
  class NumericDataTypeChecker<std::complex<int> >{ //same a as signed
  public:
    static const bool IsComplex = true;
    static const bool IsIntegralType = false;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(std::complex<int>).name() <<"\" is not a numeric data type");
    }
  };

  template<>
  class NumericDataTypeChecker<std::complex<unsigned> >{  //same as 'unsigned int' and 'size_t'
  public:
    static const bool IsComplex = true;
    static const bool IsIntegralType = false;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(std::complex<unsigned>).name() <<"\" is not a numeric data type");
    }
  };
  

  template<>
  class NumericDataTypeChecker<std::complex<long> >{ //same as 'long int'
  public:
    static const bool IsComplex = true;
    static const bool IsIntegralType = false;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(std::complex<long>).name() <<"\" is not a numeric data type");
    }
  };


  template<>
  class NumericDataTypeChecker<std::complex<bool> >{
  public:
    static const bool IsComplex = true;
    static const bool IsIntegralType = false;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(std::complex<bool>).name() <<"\" is not a numeric data type");
    }
  };
  
  template<>
  class NumericDataTypeChecker<std::complex<char> >{
  public:
    static const bool IsComplex = true;
    static const bool IsIntegralType = false;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(std::complex<char>).name() <<"\" is not a numeric data type");
    }
  };
  
  template<>
  class NumericDataTypeChecker<std::complex<wchar_t> >{ //wide char
  public:
    static const bool IsComplex = true;
    static const bool IsIntegralType = false;
    static inline void certify(){
      ADONIS_ERROR(TypeError, "Type \""<<typeid(std::complex<wchar_t>).name() <<"\" is not a numeric data type");
    }
  };



  /**
   * \brief Check if 2 objects are of the same type
   *
   * EXAMPLE:
   *\code
     std::cout << "Same types? "<< AreSameTypes<double,MyVec<double> >::used() << std::endl;
     std:: cout << "Same types? "<< AreSameTypes<double,MyVec<double>::value_type >::used() << std::endl;
   * \endcode
   */
  template<class T1, class T2>
  class AreSameTypes{
  public:
    static inline bool used(){
      return ( (typeid(T1) == typeid(T2)) ? true : false );
    }
  };


  //! a bool does not have references, it's just the value itself
  template<class T>
  class WhenBool{
  public:
    typedef T& Reference;
    typedef const T& ConstReference;
  };

  template<>
  class WhenBool<bool>{
  public:
    typedef bool Reference;
    typedef bool ConstReference;
  };

} //end namespace

#endif
