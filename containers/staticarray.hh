#ifndef STATIC_ARRAY_HH
#define STATIC_ARRAY_HH

#include <iostream>

#include "../templatemetaprograms/unrollloop.hh"
#include "../common/adonisassert.hh"

namespace Adonis{
  
  template<class T, int N>
  class StaticArray{
  public:
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;

    enum{dimension = N};

    const std::size_t size() const{
      return static_cast<std::size_t>(N);
    }

    StaticArray(){}

    StaticArray(const T& val){
      UnrollLoop<0,N>::assign_value(&a_[0],val);
    }

    template<class IT>
    StaticArray(IT start, IT end){
      adonis_assert(static_cast<int>(std::distance(start,end) == N));
      
      UnrollLoop<0,N>::assign(&a_[0],start);
    }

    StaticArray& operator=(const T& val){
      UnrollLoop<0,N>::assign_value(&a_[0],val);
      return *this;
    }
    
    StaticArray& operator=(const StaticArray& b){
      if(this != &b)
	UnrollLoop<0,N>::assign(&a_[0],&b[0]);
      return *this;
    }

    T& operator[](int i){
      adonis_assert( i >= 0 && i < N);
      return a_[i];
    }

    const T& operator[](int i) const{
      adonis_assert( i >= 0 && i < N);
      return a_[i];
    }

    T& operator()(int i){
      adonis_assert( i >= 0 && i < N);
      return a_[i];
    }

    const T& operator()(int i) const{
      adonis_assert( i >= 0 && i < N);
      return a_[i];
    }


    iterator begin(){
      return &a_[0];
    }

    iterator end(){
      return &a_[N];
    }

    const_iterator begin() const{
      return &a_[0];
    }

    const_iterator end() const{
      return &a_[N];
    }


    
    friend std::ostream& operator<<(std::ostream& os, const StaticArray& a){
      os << "[  ";
      //for(int i = 0; i < N; ++i)
      for(const_iterator cit = a.begin(); cit != a.end(); ++cit)	
	os << *cit << "  ";
      os << "]" << std::endl;

      return os;
    }

   
  private:
    T a_[N > 0 ? N : 1];
  };

}

#endif
