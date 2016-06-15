#ifndef CHECK_WHETHER_A_TEMPLATE_PARAMETER_IS_A_CLASS_TYPE_HH
#define CHECK_WHETHER_A_TEMPLATE_PARAMETER_IS_A_CLASS_TYPE_HH

/** \brief This cool code snippet tests whether T is a class or a non-class type
 *  It uses the SFINAE concept which can be found in 
      @Book{vanjos2003,
        author = {David Vandevoorde and Nicolai M. Josuttis},
	title = {C++ Templates -- The Complete Guide},
	publisher = {Addison-Wesley},
	year = {2003},
	isbn = {0-201-73484-2}
      }
 *  herein see §8.3.1 p.106/7 and §15.2.2 pp.266      
 *
 *  USAGE:
 *   class MyClass{};           //test with some classes, structs, functions,..
 *   struct MyStruct{};
 *   void myfunc(){}
 *   enum MyEnum{e1} enM;
 *
 *    std::cout << is_class<MyClass>() <<std::endl;   //returns 1
 *    std::cout << is_class(enM) <<std::endl;         //returns 0
 *    MyStruct s;
 *    is_class(s);                                    //returns 1
 *    std::cout << is_class<int>() <<std::endl;       //returns 0
 *    std::cout << is_class(myfunc) <<std::endl;      //returns 0
 */

namespace Adonis{

  template<class T>
  class IsClass{
  private:
    typedef char one_;
    typedef struct{char a[2];} two_;
    template<class C> static one_ test(int C::*); //this is valid only for 
                                                  //classes and fails for
                                                  //types, such as 'int'
    template<class C> static two_ test(...); //substitution-failure-is-not-an-
                                             //error (SFINAE) principle. This
                                             //has no problems with e.g. 'int'
  public:
    enum{Yes = sizeof(test<T>(0)) == 1};
    enum{No = !Yes};
  
    //!alternative 
    static const bool Value = sizeof(test<T>(0)) == sizeof(one_);
  };

  //!non-compile time access
  template<class T>
  bool is_class(){   //check by passing type as template argument
    if(IsClass<T>::Yes){
      return true;
    }
    else{
      return false;
    }
  }

  template<class T>
  bool is_class(T){  //check by passing type as function call argument
    return is_class<T>();
  }
  

  /**
   * \brief Check if a type is a container. The most prominent fact about a well-implemented container is that it's endowed with at least an <I> iterator </I>. 
   * Ad.: Maybe the classname is misleading; it should be renamed as HasIterator since an expression template does have an iterator  
   */
  template<class T>
  class IsContainer{
  private:
    typedef char one_;
    typedef struct{char a[2];} two_;
    //!Dont forget the wildcard '*' at the end
    //!NOTE: SFINAE is only applicaple w.r.t. checking for typedefs 
    template<class C> static one_& test(typename C::iterator*); //this is valid only for
                                                  //containers and fails for
                                                  //types, such as 'int', 'std::complex' etc.
    template<class C> static two_& test(...); //substitution-failure-is-not-an-
                                             //error (SFINAE) principle. This
                                             //has no problems with e.g. 'int'
  public:
    enum{Yes = sizeof(test<T>(0)) == 1};
    enum{No = !Yes};
    
    //!alternative 
    static const bool Value = sizeof(test<T>(0)) == sizeof(one_);

  };

  
  /**
   * \brief check if a given type is a DUNE field type. Those currently have a 'field_type' instead of a 'value_type' which is often a nuisance 
   */
  template<class T>
  class IsADuneFieldContainer{
  private:
    typedef char one_;
    typedef struct{char a[2];} two_;
    //!Dont forget the wildcard '*' at the end
    template<class C> static one_& test(typename C::field_type*); //this is valid only for DUNE's FieldVector and FieldMatrix
    template<class C> static two_& test(...); //substitution-failure-is-not-an-
                                             //error (SFINAE) principle. This
                                             //has no problems with e.g. 'int'
  public:
    enum{Yes = sizeof(test<T>(0)) == 1};
    enum{No = !Yes};
    
    //!alternative 
    static const bool Value = sizeof(test<T>(0)) == sizeof(one_);
  };


  /**
   * \brief checks if a given functor has been specified as an implementation of some method of lines.
   */
  template<class T>
  class IsMOLFunctor{
  private:
    typedef char one_;
    typedef struct{char a[2];} two_;
    //!Dont forget the wildcard '*' at the end
    template<class C> static one_& test(typename C::mol_type*); //this is valid only for functors that have been explicitly marked as MOL functors
    template<class C> static two_& test(...); //substitution-failure-is-not-an-
                                             //error (SFINAE) principle. This
                                             //has no problems with e.g. 'int'
  public:
    enum{Yes = sizeof(test<T>(0)) == 1};
    enum{No = !Yes};
    
    //!alternative 
    static const bool Value = sizeof(test<T>(0)) == sizeof(one_);
  };

}//end namespace


#endif
