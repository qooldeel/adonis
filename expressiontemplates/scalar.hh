#ifndef SCALAR_CLASS_FOR_EXPRESSIONS_HH
#define SCALAR_CLASS_FOR_EXPRESSIONS_HH

#include <iostream>

namespace Adonis{
  namespace ExprTmpl{

//An object representing a scalar, i.e. a constant or literal
    template<typename T>
    class AScalar{
    private:
      const T& scal_;   //value of scalar
      
    public:
      
      typedef T value_type;

      //constructor initialises value
      AScalar(const T& s = 0. ):scal_(s){}    //also default constructor
  

      //define dereferencing operator which is NEEDED to hit the requirements of an iterator
      T operator*() const{
	return scal_;
      }

      //increment operator needed for iteration
      void operator++(){}   //do nothing for a scalar
  

      const T& get_scalar() const{
	return scal_;
      }
  
      //define an index operation similar to vectors so that you can use it like a vector ;-)
      inline const T& operator[](size_t) const{
	return scal_;
      }

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      const T& operator()(int, T& s, T& c) const{
	return scal_;
      }
#endif

      //size of scalar is zero by definition
      size_t size() const{return 0;}  //0 
    

      //+++++++++++++ define iterator +++++++++++++++
      //typedef T* iterator;
      typedef const T* const_iterator;
      
      //iterator begin() {return scal_;}  //nothing
      //iterator end() {return scal_;}

      const_iterator begin() const{return &scal_;}
      const_iterator end() const{return &scal_+1;} //plus one !!!
      
      //+++++++++++++++++++++++++++++++++++++++++++++

    };

    //end namespace(s)
  }
}

#endif
