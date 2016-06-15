#ifndef SIMPLE_ALLOCATION_HH
#define SIMPLE_ALLOCATION_HH

#include <iostream>

namespace Adonis{

  template<class T> 
  class MyAllocator{
  public:
    typedef T* PointerType;
    typedef std::size_t Size_Type;

    MyAllocator(PointerType pstart = 0, PointerType pend = 0):memStart_(pstart),memEnd_(pend){}

    PointerType allocate(SizeType n){
      return  memStart_ + n;
    }

  protected:
    PointerType memStart_,
      memEnd_;
  };

}//namespace 

#endif
