#ifndef EXPRESSION_TEMPLATES_4_MATRIX_VEC_AND_VEC_MATRIX_OPERATIONS_HH
#define EXPRESSION_TEMPLATES_4_MATRIX_VEC_AND_VEC_MATRIX_OPERATIONS_HH

#include "../common/adonisassert.hh"

namespace Adonis{

  namespace ExprTmpl{

    /**
     * \brief Multiply a Dune::FieldMatrix-compliant object (i.e. element access via [][] and matrix dimensions via .rows and .cols and typedef field_type) with an STL-compliant random access container (i.e. possesses at least .size() member and typedef value_type)
     *
     * \tparam M Matrix object
     * \tparam O Operation object
     * \tparam V Vector object
     *
     * Note: it is not suggested to use it as stand-alone-class. Instead use it in conjunction with an expression vector.
     *
     * TODO: in every expression container class you must overload the operator= for the matrix-vector and vector-matrix operations, respectively, e.g. 
     * \code
      template<class M, class O, class V>
      MyVec& operator=(MOV<M,O,V> ex){
      for(int i = 0; i < ex.rows(); ++i){
       for(int j = 0; j < ex.cols(); ++j){
        v_[i] += ex(i,j);
       }
     } 
     return *this;
     } 

  //v*M
   template<class M, class O, class V>
   MyVec& operator=(VOM<M,O,V> ex){
    for(int i = 0; i < ex.rows(); ++i){
      for(int j = 0; j < ex.cols(); ++j){
	v_[j] += ex(i,j);
      }
    }
    return *this;
  }
     *  \endcode
     *
     * USAGE: Suppose you have an expression vector class named MyVec.
     * \code 
     double smurf[] = {1,3,9,6};
     ExprTmpl::MyVec<double> v(smurf, smurf+4);
  
     Dune::FieldMatrix<double,3,4> M;
     M[0][0] = 2; M[0][1] = 5;  M[0][2] = -3; M[0][3] = 1;
     M[1][0] = 0; M[1][1] = 2;  M[1][2] = 3;  M[1][3] = 2;
     M[2][0] = 1; M[2][1] = -1; M[2][2] = -3; M[2][3] = 6;
  
     ExprTmpl::MyVec<double> vtimesM(3);  //necessary to initialise!
  
     vtimesM = M*v;

     std::cout << "M*v = "<< std::endl;
     std::cout << vtimesM << std::endl;

     * \endcode
     */
    template<class M, class O, class V>
    class MOV{
    private:
      const M& m_;        //stores constant references
      const V& v_;

    public:
      typedef typename V::value_type value_type;
      
      int rows() const {return m_.rows;}
      int cols() const {return m_.cols;}
      
      MOV(const M& m = M(), const V& v = V()):m_(m), v_(v){}
      
      value_type operator()(int i, int j){
	return O::apply(m_[i][j],v_[j]);
      }
    };

    
    /**
     \brief Multiplication of a STL-compliant random access container with a Dune::FieldMatrix-compliant object, a.k.a. left-multiplication with DUNE-matrix
     */
    template<class M, class O, class V>
    class VOM{
    private:
      const M& m_;        //stores constant references
      const V& v_;

    public:
      typedef typename V::value_type value_type;

      int rows() const {return m_.rows;}
      int cols() const {return m_.cols;}
      
      VOM(const M& m = M(), const V& v = V()):m_(m), v_(v){}
      
      value_type operator()(int i, int j){
	return O::apply(v_[i],m_[i][j]);
      }
    };
 
    /**
     \brief More general multiplication.
     \tparam RT Return type (either T1 or T2)
     \tparam T1 Object1 
     \tparam T2 Object2 
     */
    template<class RT, class T1, class T2>
    class Multiplication{
    public:
      Multiplication(){}
      
      static inline RT apply(const T1& a, const T2& b){
	return a*b;
      }
    };
 

    //Finally overload operator* for M*v
    template<class M, class V>
    inline MOV<M,Multiplication<typename V::value_type,typename M::field_type,typename V::value_type>,V> operator*(const M& m, const V& v){
      adonis_assert(m.cols == (int)v.size());

      typedef MOV<M,Multiplication<typename V::value_type,typename M::field_type,typename V::value_type>,V > ExpressionType;

      return ExpressionType(m,v);
    }


    
    //as well as for v*M
    template<class M, class V>
    inline VOM<M,Multiplication<typename V::value_type,typename M::field_type,typename V::value_type>,V> operator*(const V& v, const M& m){
      adonis_assert((int)v.size() == m.rows);
  
      typedef VOM<M,Multiplication<typename V::value_type,typename M::field_type,typename V::value_type>,V> ExpressionType;

      return ExpressionType(m,v);
    }

   
 
  }//end namespaces
    
}

#endif
