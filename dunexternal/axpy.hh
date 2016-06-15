#ifndef SCALAR_A_TIMES_X_PLUS_Y_HH
#define SCALAR_A_TIMES_X_PLUS_Y_HH

/** \brief The (S)AXPY ((scalar) a times x plus y) operation.
 *  The usual way to achieve it would be to calculate z = a*x + y.
 *  This works well but performs very poorly since we need a loop to calculate
 *  a*x, a temporary to store the result, say z1, another loop to compute 
 *  z1+y and a temporay z2 for storage and yet another loop to assign z2 to z 
 *  via z = z2.
 *
 *  Much more efficiency is attained by deferring the "freed-of-temporaries-" 
 *  evaluation of all those operations mentioned above and later recombine and 
 *  assign them in one "ideal" loop (instead of three)! 
 *  This is actually the idea of which Expression Templates are made ! 
 *
 *  The only thing which needs to be done is to define a copy constructor and
 *  copy assignment in class Essential, where the evaluation of the saxpy is 
 *  finally done (Note: In order to avoid ambiguity, it might be necessary not
 *  to define an a*x operator in Essential again).
 *
 *   VARIATIONS of SAXPY that are possible:
 *       a*x + y
 *       x*a + y
 *       x + a*y
 *       x + y*a
 *
 *      No "extended" saxpy is possible (e.g. NO x+y, a*x + b*y etc.)
 *
 *    USAGE (efficient usage of saxpy operation on 'Essential' objects):
 *   \code
        std::cout << "Copy assignments on Axpy:" << std::endl;
        Essential<double,6,2> saxpy_copy =  1.5*E3 + E4;
        std::cout << saxpy_copy << std::endl;

        Essential<double,6,2> saxpy(0.); 
        saxpy =  1.5*E3 + E4;
        std::cout << saxpy << std::endl;

	std::cout << "Copy construction on Axpy:" << std::endl;
	Essential<double, 6,2> R1(0.5*E1 + E2);
	std::cout << R1 << std::endl;
*/


#include "axpytraits.hh"  //contains forward declarations and typedefs


namespace Adonis{

  //AX: calculate a*x -- no temporaries in use ;-) 
  template<class K, int M, int N,bool B>  
  class Ax{                      
  public:
    typedef typename AxpyTraits<K,M,N,B>::Type Type;

    Ax(const K& a, const Type& x):a_(a), x_(x){}
    
    //access fields of Ax
    const K& a() const{return a_;}
    const Type& x() const{return x_;}

  private:
    const K& a_;       //reference to avoid copying
    const Type& x_;
  };

  
  //operators on Ax
  template<class K, int M, int N,bool B> 
  inline Ax<K,M,N,B> operator*(const K& a, const Essential<K,M,N,B>& x){   //a*x
    return Ax<K,M,N,B>(a,x);
  }

  template<class K, int M, int N,bool B> 
  inline Ax<K,M,N,B> operator*(const Essential<K,M,N,B>& x, const K& a){  //x*a  
    return Ax<K,M,N,B>(a,x);
  }



  //AXPY: calculate a*x+y by using class Ax  
  template<class K, int M, int N,bool B>
  class Axpy{
  public:
    typedef typename AxpyTraits<K,M,N,B>::Type Type; 
    typedef typename AxpyTraits<K,M,N,B>::AxType AxType;
    
    //access fields of Axpy
    const K& a() const{return a_;}
    const Type& x() const{return x_;}
    const Type& y() const{return y_;}

    Axpy(const AxType& ax, const Type& w):a_(ax.a()), x_(ax.x()), y_(w){}    
  
  private:
    const K& a_;           //the scalar  
    const Type& x_;        //references to Essential<K,M,N> objects
    const Type& y_;
  };
  
  
  //operators on Axpy
  template<class K, int M, int N,bool B>
  inline Axpy<K,M,N,B> operator+(const Ax<K,M,N,B>& ax, const Essential<K,M,N,B>& x){
    return Axpy<K,M,N,B>(ax,x);   // a*x+y
  }

  template<class K, int M, int N,bool B>
  inline Axpy<K,M,N,B> operator+(const Essential<K,M,N,B>& x, const Ax<K,M,N,B>& ax){
    return Axpy<K,M,N,B>(ax,x);   // x+a*y
  }

  


  /** \brief A LESS EFFICIENT way of defining the saxpy operation
   * \code
      Axpy(const K& a, const Type& x, const Type& y):a_(a),x_(x),y_(y){} 
    
      Type axpy() const{   //does the job with temporaries
        Type Result = x_;             
        //Result = K();   //necessary !
      
        //Result += x_;  //contains x_ now
        Result *= a_;  //multiply a_ with x_
        Result += y_;  //finally add y_
    
        return Result;
      }
  */
  
}//end of namespace

#endif
