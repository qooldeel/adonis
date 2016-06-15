#ifndef SIMPLE_DOT_PRODUCT_TEST_ONLY_HH
#define SIMPLE_DOT_PRODUCT_TEST_ONLY_HH


#include <iostream>
#include "../common/globalfunctions.hh"
#include "../common/adonisassert.hh"

namespace Adonis{

  /*  A simple programme to calculate the dot product in a dashing way ;-)  

      DRAWBACK: vector size N has to be known at compile time!
  
      IN GENERAL: vector sizes are often known a priori.

      ADVANTAGE: speeds up computations tremendously!
    
      --------------------------------------------------------------------
      --------------------------------------------------------------------
      
      USAGE:
      
        int main() {
      
          MetaVec<3, double> a, b;
 
          a[0] = 1; a[1] = 2; a[2] = 3;
          b[0] = 5; b[1] = 6; b[2] = 7;

          std::cout << " Dot product <a,b> = " << dot(a,b) <<std::endl; 
        }


      OUTPUT: Dot product <a,b> = 38

      LOOKING BEHIND THE SCENES:

        The inlined function 'dot(a,b)' is expanded during compilation time as
        follows:

	dot(a,b) 
	= MetaDot<3>::result(u,v)
	= u[2]*u[2] + MetaDot<2>::result(u,v)
	= u[2]*u[2] + u[1]*u[1] + MetaDot<1>::result(u,v)
	= u[2]*u[2] + u[1]*u[1] + u[0]*u[0]


      REFERENCES:
        [VANDEVOORDE/JOSUTTIS, C++ Templates, §17.7, pp. 314]
        [YANG, C++ and OO Numeric Computing, §7.7.4, pp. 271]
  */

  //Very tiny NON-DYNAMIC vector class (i.e. you cannot use it for vectors whose size isn't known at compilation time!
  template<int N, typename T>
  class MetaVec{          
  private:
    T entry[N];
  public:
    T& operator[](int i) {
      adonis_assert(i>=0 && i < N);
      return entry[i];
    }
    
    inline int size() const{ return N;}
    

    friend inline std::ostream& operator<<(std::ostream& os, MetaVec& mv){
      os << "  MetaVec = { ";
      for(int i = 0; i < N; ++i){
	os << mv[i] << "  "; 
      }
      os << "}"<<std::endl;
      return os;
    }
  };

 

  //####################   TEMPLATE METAPROGRAMME #############################
  
  
  //WORKS !
  template<int DIM>
  class MetaDot{
  public:
     template<int N, typename T>
     static T result(MetaVec<N,T>& u, MetaVec<N,T>& v){
       return u[DIM-1]*v[DIM-1] + MetaDot<DIM-1>::result(u,v); //recurrence
    }
  };


  //partial specialisation -- end criterion of recurrence 
  template<>
  class MetaDot<1>{
  public:
     template<int N, typename T>
     static T result(MetaVec<N,T>& u, MetaVec<N,T>& v){
       return u[0]*v[0];
    }
  };
  
 
  //##########################################################################



  //now we define a CONVENIENT function
  template<int N, typename T>
  inline T dot(MetaVec<N,T>& a, MetaVec<N,T>& b){
    adonis_assert(a.size() == b.size());
    return MetaDot<N>::result(a,b);
  }
  

 
} //end namespace 



#endif
