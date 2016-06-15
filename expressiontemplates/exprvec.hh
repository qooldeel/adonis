#ifndef MY_EXPRESSION_VECTOR_CLASS_HH
#define MY_EXPRESSION_VECTOR_CLASS_HH


#include <iostream>
#include <vector>
#include <cassert>
#include <stdlib.h>  //for 'exit' function
#include <limits> 
#include <cmath>
#include <string>

#include "scalar.hh"
#include "../common/error.hh"
#include "../common/globalfunctions.hh"

#include "../common/adonisassert.hh"

#include "matrixvecexpression.hh" //expression templates for M*v and v*M
#include "../common/numerictypechecker.hh"

#include "xsettings.hh"

#ifdef USE_TUNED_FLOATING_POINT_ARITHMETIC
#include "../accuracy/floatingpointarithmetic.hh"
#endif

#include "../misc/commaoverloading.hh"

namespace Adonis{
  
  class ExpressVecError: public Error{};

  namespace ExprTmpl{


    template<class T> class Add;  //forward declaration

    /** \brief A WORD OF CAUTION: YOU CANNOT CHECK FOR INVALID ITERATORS. ONE HAS TO AVOID USING THEM! (EVEN ASSIGNMENT AND COPYING IS NOT ALLOWED).
THE ONLY THING YOU CAN DO IS, FOR INSTANCE, TO CHECK AGAINST THE PAST-THE-END ITERATOR OF THE OBJECT UNDER CONSIDERATION (E.G. if(iterator == container.end()) 
*
*REMEDY: when "Invalid read of ..." appears, carefully check the dimension of all expression vectors involved in the expression. ALWAYS make sure they have the same size!!
*
*  
* <TT> VExpr </TT> wraps <TT> LOR </TT> which in turn administrates operations between <TT> MyVec, VExpr</TT> and <TT>  AScalar</TT>, respectively (and between a mixture of those aforementioned objects, of course) 
*
* For a detailed discussion on expression templates, check out the following
*
* References:
*
*  [1] VELDHUIZEN, "Techniques for Scientific C++", Tech. Report, 1999, p. 24-29 
*
*  [2] VELDHUIZEN, "Expression Templates", C++ Report, Vol. 7 No. 5 June 1995, pp. 26-31
*
*  The last reference can also be found in 
*  [3] LIPPMAN, "C++ Gems: Programming Pearls from the C++ Report", Cambridge University Press, 1997 
*
* URL:
*  <a href="http://www10.informatik.uni-erlangen.de/~pflaum/pflaum/ProSeminar/exprtmpl.html"> Expression Templates </a>
*/
    //############# Expressions involving Vector iterating objects.########## 
    //Actually this is a WRAPPER class for more interesting types               
    template<typename T, class A>
    class VExpr{
    private:
      A iter_;                    //stores iterator
  
    public:

      typedef T value_type;
      typedef typename A::ReturnType ReturnType;

      VExpr(const A& a = A()):iter_(a) {} //constructor


#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      ReturnType operator()(int i, T& s, T& c) const{
	return iter_(i,s,c);
      }
#endif

      //!the following 2 members are needed for my functors
      //DEREFERENCING operator: return value of iterator
      T operator*() const{
	return *iter_;
      }


      //bracket operator
      T operator[](int i) const{
	return iter_[i];
      }

      //increment iterator
      void operator++(){
	++iter_;
      }

      //postfix -- should be avoided cause less efficient
      void operator++(int){
	iter_++;
      }
      
       void operator--(){ --iter_;}
      
      //suffix increment
      void operator--(int){ iter_--;}
    };
 
    
    //########################### My VECTOR Class #################################
    /**
     * \brief Expression template vector class that is just derived from 
     * the STL vector. Hence it shares all functionalities of the latter and
     * there is no need to redefine any of such functions, such as copy 
     * construction, copy assignment and so forth.
     *
     * \tparam T value type
     * \tparam A allocator type. Dummy 
    */
    template<typename T, class A = std::allocator<T> >
    class MyVec: public std::vector<T,A>{
    public:
      //! some abbreviations for commonly used types 
      typedef A allocator_type;
      typedef T value_type;
      typedef std::size_t IndexType;
      typedef std::size_t SizeType;
      typedef std::vector<T,A> VectorType;
      typedef typename VectorType::size_type size_type;
      
      typedef MyVec<T> ThisType;
      typedef ThisType MyOwnVectorType;

      typedef typename VectorType::reference reference;
      typedef typename VectorType::const_reference const_reference;
      
      //Simple ITERATORs
      typedef typename VectorType::iterator iterator;
      typedef typename VectorType::const_iterator const_iterator;
      

      //1st constructor -- also default construction
      MyVec(IndexType n = 0, const T& val = T(), const allocator_type& alloc = allocator_type()):VectorType(n,val,alloc){}
      
      //2nd constructor
      template<class InputIterator>
      MyVec(InputIterator first, InputIterator last, const allocator_type& alloc = allocator_type()):VectorType(first,last,alloc){}
  


      //for arbitrary container assignment
      template<class OBJ>
      MyVec& operator=(const OBJ& obj){
	//std::cout << "ExprTmpl::MyVec: ARBITRARY CONTAINER ASSIGNMENT" << std::endl;
	SizeType length = static_cast<size_t>(std::distance(obj.begin(),obj.end()));//IMPORTANT, don't forget it !!!!!!!
	if((*this).size() == 0){
	  (*this).reserve(length);
	  for(typename OBJ::const_iterator it = obj.begin(); it != obj.end(); ++it)
	    (*this).push_back(*it);
	}
	else{
	  adonis_assert((*this).size() == length);
	  unsigned index = 0;
	  for(typename OBJ::const_iterator it = obj.begin(); it != obj.end(); ++it){
	    this->operator[](index++) = *it;
	    
	  }
	}
	return *this;
      }
      

      //! usage: MyVec<double> v(5); 
      //!        v = 1.45;  // now v = [1.45,1.45,..., 1.45]
      MyVec& operator=(const T& d){
        for(iterator it = (*this).begin(); it!= (*this).end(); ++it)
          *it = d;
        
        return *this;
      }

      //!special assignment. Usage:
      //! MyVec<T> v(3);
      //! v <<= 1., 2., 3.;
      //! I decided to overload this operator in order to avoid 
      //! disambiguity with conventional assignment (=) operators.  
      CommaOverloading<T,iterator> operator<<=(const T& val){
      SizeType count = 0;
      (*this)[0] = val;
      return CommaOverloading<T,iterator>((*this).begin()+1,(*this).size(),count);
    }

      //!compare operators
      bool operator==(const MyVec& v) const{
	if(this == &v)     //it it's the identity return true
	  return true;
	if((*this).size() != v.size())     //if vec's don't have the same size
	  return false;

	bool cool = true;
	for(IndexType i = 0; i < (*this).size(); ++i){
	  if((*this)[i] != v[i]) {       //if entries differ
	    cool = false;
	    break;
	  }
	}  
	return cool;
      }

      bool operator!=(const MyVec& v) const{
	return !operator==(v);
      }

      //approximately equal -- I overload bitwise xor for that purpuse, due to
      //its similarity to the \f$ \approx\f$-relation
      bool operator^=(const MyVec& v){
	if(this == &v)     //it it's the identity return true
	  return true;
	if((*this).size() != v.size())     //if vec's don't have the same size
	  return false;
	
	bool cool = true;
	for(IndexType i = 0; i < (*this).size(); ++i){
	  if(!is_equal(this->operator[](i),v[i])) {       //if entries differ approximately
	    cool = false;
	    break;
	  }
	}  
	return cool;

      }

      //!introduce paranthesis operator that wraps [] and can be used in the
      //!same matter
       reference operator()(size_type i){
	 adonis_assert(i < (*this).size());
	 return this->operator[](i);
       }

      //read only access
      const_reference operator()(size_type i) const{
	assert(i < (*this).size());
	return (*this)[i];
      }
      
     
      /**
       * \brief fill vector from a matrix object which supports member functions 'rows()' and 'cols()' to determine the matrix sizes.
       * NOTE: row-wise storage, as usual
       */
      template<class MatrixType>
      void fill_from_matrix(const MatrixType& M){
	if((*this).size() == 0)
	  (*this).resize(M.rows()*M.cols());

        else
	  adonis_assert(M.rows()*M.cols() == (*this).size());

	int idx = 0;
	
	for(int i = 0; i < (int)M.rows(); ++i){
	  for(int j = 0; j < (int)M.cols(); ++j){
	    this->operator[](idx++) = M[i][j];
	  }
	}
      }


      //output 
      std::ostream& print(std::ostream& os) const{
	std::string nice = ", ";   //make output a feast for the eyes ;-) 
	os<< " [ ";
	for(IndexType i = 0; i < (*this).size(); ++i){
	  if(i == ((*this).size())-1)             //if last element is read,
	    nice = " ";                  //output only a space
	  os << (*this)[i] << nice;  
	}
	os <<"]"<< std::endl;
      
	return os;
      }


      //==================== assign Expression of type VExpr to vector ========
      template<class E>
      MyVec& operator=(VExpr<T,E> result){ //const  VExpr<T,E>& result when no iterators are used, since these are non-const references ;)
	//This can be circumvented by using operator[i], cf. (#)
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	T s = T(),
	  c = T();
	
	for(IndexType i = 0; i < (*this).size(); ++i){
	  this->operator[](i) = result(i,s,c); 
	  //!the following 2 lines will only be meaningful, if the 2nd template argument is 'b', 'B' 
	  CorrectRoundingErrorAfterwards<T,XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
	  
	  if(c != T()){  //! if c == 0 then do not assign s since 
	                 //! it might be 0 because no addition was involved!
	    CorrectRoundingErrorAfterwards<T,XPCTSettings::floatingPointAdditionMethod>::assign(this->operator[](i),s);
	  }

	  s = c = T();   //reset afterwards!!
	}                  

#else

	////!Here you can use operator=(const VExpr<T,E>& expr) 
	// for(IndexType i = 0; i < (*this).size(); ++i){
	//   this->operator[](i) = result[i];                     // (#)
	// }

//!NOTE: if iterators (which must be non-const. here) are used
//!      then you have to rewrite operator=(const VExpr<T,E>& expr) 
//!      as operator=(VExpr<T,E> expr) 
//!only one "IDEAL" loop instead of a vast amount of loops & temp. objs!
	for(iterator it = (*this).begin(); it != (*this).end(); ++it){

	  //expression can be built without valgrind-error!
	  //it is in the assignment that we can check it properly
	  //bounds_check(*result);
	  
	  *it = *result; //assign expression to vector where it is finally stored
	  ++result;      //increment result iterator as well 
    
	}
#endif

	return *this;
      }




      //====================  UNARY OPERATORS ================================  
      template<class E>
      MyVec& operator+=(VExpr<T,E> result){     //+= Expression
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	T s = T(),
	  c = T();
	
	for(IndexType i = 0; i < (*this).size(); ++i){
	  //s = 0; // since all information are contained in result(i,s,c)
	  //! result(i,s,c) and this->operator[](i) must be in exactly this order!!
	  this->operator[](i) = Add<T>::apply(result(i,s,c),this->operator[](i),s,c);

	  //!only meaningful for Babuska-Kahan
	  CorrectRoundingErrorAfterwards<T,XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
	  if(c != T()){  //! if c == 0 then do not assign s since 
	                 //! it might be 0 because no addition was involved!
	    CorrectRoundingErrorAfterwards<T,XPCTSettings::floatingPointAdditionMethod>::assign(this->operator[](i),s); 
	  }
	  s = c = T();  //reset for next index i
	}

#else
	//!standard +=
	for(iterator it = (*this).begin(); it != (*this).end(); ++it, ++result){
	  //increment result as well
	  *it += *result; 
	}

#endif
	return *this;
      }
  

      //! += MyVec corresponds to the additive assignement of 2 numbers per 
      //! index <TT> i </TT> (unlike when an expression is assigned!) 
      MyVec& operator+=(const MyVec& v){        
    
	//assert(size() == v.size()); //OR:
	if((*this).size() != v.size()){
	  ADONIS_ERROR(ExpressVecError, "The two MyVec sizes don't match: left size = "<< (*this).size() << "  right size = "<<  v.size() << "."); 
	}

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	
	typedef IntelligentArithmeticOperation<XPCTSettings::floatingPointAdditionMethod> PlusEqual;

	for(IndexType i = 0; i < (*this).size(); ++i){
	  PlusEqual::pluseq(this->operator[](i),v[i]); //!only addition of 2 numbers at a time
	}
	
#else
	const_iterator jt = v.begin();
	for(iterator it = (*this).begin(); it != (*this).end(), jt != v.end(); ++it, ++jt){
	  *it += *jt; //apply '+='
	}

#endif
	return *this;
      }

      //the same for minus
      template<class E>
      MyVec& operator-=(VExpr<T,E> result){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	T s = T(),
	  c = T();
	
	T val;
	for(IndexType i = 0; i < (*this).size(); ++i){
	  //s = 0; // since all information are contained in result(i,s,c)
	  //! result(i,s,c) and this->operator[](i) must be in exactly this order!!
	  val = result(i,s,c); //evaluate expression containing possible '+' 
	  
	  //!only meaningful for Babuska-Kahan
	  CorrectRoundingErrorAfterwards<T,XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
	  if(c != T()){  //! if c == 0 then do not assign s since 
	                 //! it might be 0 because no addition was involved!
	    CorrectRoundingErrorAfterwards<T,XPCTSettings::floatingPointAdditionMethod>::assign(val,s);
	  }
	  
	  this->operator[](i) -= val;
	  
	  s = c = T();  //reset for next index i
	}

#else

	for(iterator it = (*this).begin(); it != (*this).end(); ++it){
     
	  //bounds_check(*result);

	  *it -= *result; //apply '-='
	  ++result;      //increment result as well 
	}
#endif
	return *this;
      }
  
      MyVec& operator-=(const MyVec& v){
   
	//assert(size() == v.size());
	if((*this).size() != v.size()){
	  ADONIS_ERROR(ExpressVecError, "Expression and MyVec size don't match."); 
	}
	#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	
	typedef IntelligentArithmeticOperation<XPCTSettings::floatingPointAdditionMethod> MinusEqual;

	for(IndexType i = 0; i < (*this).size(); ++i){
	  MinusEqual::minuseq(this->operator[](i),v[i]); //!only addition of 2 numbers at a time
	}
	
#else
       
	const_iterator jt = v.begin();
	for(iterator it = (*this).begin(); it != (*this).end(), jt != v.end(); ++it, ++jt){
	  *it -= *jt; //apply '-='
	}
#endif
	return *this;
      }
  

      //...and for multiplication
      template<class E>
      MyVec& operator*=(VExpr<T,E> result){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	T s = T(),
	  c = T();
	
	T val;
	for(IndexType i = 0; i < (*this).size(); ++i){
	  //s = 0; // since all information are contained in result(i,s,c)
	  //! result(i,s,c) and this->operator[](i) must be in exactly this order!!
	  val = result(i,s,c); //evaluate expression containing possible '+' 
	  
	  //!only meaningful for Babuska-Kahan
	  CorrectRoundingErrorAfterwards<T,XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);

	  if(c != T()){  //! if c == 0 then do not assign s since 
	                 //! it might be 0 because no addition was involved!
	    CorrectRoundingErrorAfterwards<T,XPCTSettings::floatingPointAdditionMethod>::assign(val,s);
	  }
	  
	  this->operator[](i) *= val;  //MULTIPLY=
	  
	  s = c = T();  //reset for next index i
	}

#else

	for(iterator it = (*this).begin(); it != (*this).end(); ++it){

      
	  *it *= *result; //apply '*='
	  ++result;      //increment result as well 
    }
#endif
	return *this;
      }
  
      MyVec& operator*=(const MyVec& v){
   
	//assert(size() == v.size());
	if((*this).size() != v.size()){
	  ADONIS_ERROR(ExpressVecError, "Expression and MyVec size don't match."); 
	}
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	
	typedef IntelligentArithmeticOperation<XPCTSettings::floatingPointAdditionMethod> TimesEqual;

	for(IndexType i = 0; i < (*this).size(); ++i){
	  TimesEqual::timeseq(this->operator[](i),v[i]); //!only addition of 2 numbers at a time
	}
	
#else	
	const_iterator jt = v.begin();
	for(iterator it = (*this).begin(); it != (*this).end(), jt != v.end(); ++it, ++jt){
	  *it *= *jt; //apply '*='
	}
#endif
	return *this;
      }
  

      //!ONLY MyVec -- Scalar operations. Same for tuned and normal!
      //SPECIAL with multiplication we have also multiplication with a scalar.
      MyVec& operator*=(const T& scal){
	for(iterator it = (*this).begin(); it != (*this).end(); ++it){
	  *it *= scal; //apply '*= scal'
	}
	return *this;
      }

       MyVec& operator*(const T& scal){ //e.g. x*1.75 (similar to x *= 1.75)
	 return ( *this *= scal);
       }

      

      //DIVISION by non-zero scalar
      MyVec& operator/=(const T& scal){
	if(is_zero(scal))  //(scal == T(0.))
	  ADONIS_ERROR(ZeroDivision, "Division by zero occurred, hunk.");

	for(iterator it = (*this).begin(); it != (*this).end(); ++it){
	  *it /= scal; //apply '*= scal'
	}
	return *this;
      }
      

      MyVec& operator/(const T& scal){ //e.g. x/1.75 (similar to x /= 1.75)
	if(scal == T(0.))
	   ADONIS_ERROR(ZeroDivision, "Division by zero occurred, hunk.");
	
	return ( *this /= scal);
      }


      //! these functions are only possible for size-1 vectors!!
      MyVec& operator-=(const T& scal){
	if((*this).size() != 1)
	  ADONIS_ERROR(DimensionError, "Only possible if vector size is 1.");
	this->operator[](0) -= scal; 
	return *this;
      }

      MyVec& operator+=(const T& scal){
	if((*this).size() != 1)
	  ADONIS_ERROR(DimensionError, "Only possible if vector size is 1.");
	this->operator[](0) += scal;
	return *this;
      }

      //unary minus: -x  (i.e. x *= -1)
      MyVec& operator-(){
	return ((*this) *= static_cast<T>(-1));
      }

      //======= MATRIX-VECTOR / VECTOR-MATRIX multiplication assignment =======
      //!the crucial stuff is defined in 'matrixvecexpression.hh'
      //! note that the assignment is for the resulting vector of the matrix-
      //! vector resp. vector-matrix multiplication
      template<class M, class O, class V>
      MyVec& operator=(MOV<M,O,V> ex){
	adonis_assert(typeid(typename MOV<M,O,V>::value_type) == typeid(value_type));

#if  USE_TUNED_FLOATING_POINT_ARITHMETIC
	value_type s = value_type(),
	  c = value_type();

        typedef IntelligentArithmeticOperation<XPCTSettings::floatingPointAdditionMethod> AdditionType;

	for(int i = 0; i < ex.rows(); ++i){
	  for(int j = 0; j < ex.cols(); ++j){
	    AdditionType::add(ex(i,j),s,c);  //! s is the sum
	  }
	  //!only meaningful for Babuska-Kahan -- nothing to be done for Kahan
	  CorrectRoundingErrorAfterwards<value_type,XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
	  this->operator[](i) = s;             //! finally assign sum
	                         //! note summation is always applied!
	  s = c = value_type();  //! reset sum and correction
	}
	
#else

	for(int i = 0; i < ex.rows(); ++i){
	  for(int j = 0; j < ex.cols(); ++j){
	    this->operator[](i) += ex(i,j);
	  }
	}
#endif
	return *this;
      }


      //v*M
      template<class M, class O, class V>
      MyVec& operator=(VOM<M,O,V> ex){
	  adonis_assert(typeid(typename VOM<M,O,V>::value_type) == typeid(value_type));
	  
#if  USE_TUNED_FLOATING_POINT_ARITHMETIC
	typedef IntelligentArithmeticOperation<XPCTSettings::floatingPointAdditionMethod> AdditionType;
	
	NumericDataTypeChecker<value_type>::certify();
	
	
	VectorType s((*this).size()), c((*this).size());

	for(int i = 0; i < ex.rows(); ++i){
	  for(int j = 0; j < ex.cols(); ++j){
	    this->operator[](j) = AdditionType::add(ex(i,j),s[j],c[j]);
	  }
	}
	

	//! only meaningful for Babuska-Kahan summation
	CorrectRoundingErrorAfterwards<value_type,XPCTSettings::floatingPointAdditionMethod>::matvecmult_special_x((*this).begin(),s,c);

#else
	for(int i = 0; i < ex.rows(); ++i){
	  for(int j = 0; j < ex.cols(); ++j){
	    this->operator[](j) += ex(i,j);
	  }
	}
#endif
	return *this;
      }

      //============================================================
    };//end MyVec class


    
    //##########################################################################
    //#### APPLICATIVE TEMPLATE CLASSES provide "atomic" operations(+,*,..) ###
    //#########################################################################
    template<typename T>
    class Add{              //perform addition between two expression fields []
    public:
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      typedef T& ReturnType;
#else
      typedef T ReturnType;
#endif

      // Add(){}

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      static inline ReturnType apply(const T& a, const T& b, T& s, T& c){
	//! elementary addition using either Kahan or K.-Babuska summation
	return BalancingAlgorithm<T,XPCTSettings::floatingPointAdditionMethod>::addition(a,b,s,c);
      }
#endif
      
      //!since e.g. operator() cannot be 'static' one has to circumvent this drawback by defining an 'apply' function:
      static inline T apply(const T& a, const T& b){
	//	std::cout << "element-standard addition applied" << std::endl;
	return a+b;     //PLUS
      }
    
    };

    template<typename T>
    class Subtract{       //perform subtraction between two expression fields []
    public:
      typedef T ReturnType;

      //Subtract(){}
 
      //! fake operator -- leave s and c unchanged    
#if USE_TUNED_FLOATING_POINT_ARITHMETIC 
      static inline T apply(const T& a, const T& b, T& s, T& c){
	return a-b;
      }
#endif
 
      static inline T apply(const T& a, const T& b){
	return a-b;    //MINUS
      }
      
};

    template<typename T>
    class Multiply{    //perform mulitplication between two expression fields []
    public:
      typedef T ReturnType;

      //Multiply(){}

      //! fake operator -- leave s and c unchanged  
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      static inline T apply(const T& a, const T& b, T& s, T& c){
	return a*b;
      }
#endif
      
      static inline T apply(const T& a, const T& b){
	return a*b;   //TIMES
      }
    

    };
    
    template<typename T>
    class Divide{    //perform mulitplication between two expression fields []
    public:
      typedef T ReturnType;
      
      //Divide(){}

      //! fake operator -- leave s and c unchanged  
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      static inline T apply(const T& a, const T& b, T& s, T& c){
	adonis_assert(b != T());
	return a/b;
      }
#endif

      static inline T apply(const T& a, const T& b){
	if(b == T(0.))
	  ADONIS_ERROR(ZeroDivision, "Division by zero occurred, hunk.");
	
	return a/b;   //DIVISION
      }    

    };
    
    


    //####### BINARY OPERATION (Op) ON 2 EXPRESSIONS (Left, Right)#############
    template<typename T, class Left, class Op, class Right>
    class LOR{
    private:
      Left iter1_;
      Right iter2_;

    public:

      typedef T value_type;
      typedef typename Op::ReturnType ReturnType;

      LOR(const Left& i1 = Left(), const Right& i2 = Right()):iter1_(i1), iter2_(i2){}

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      ReturnType operator()(int i, T& s, T& c) const{
	return Op::apply(iter1_(i,s,c), iter2_(i,s,c), s, c); 
      }
#endif
      
      //these are analogues to the definitions in VExpr<.,.>. see above
      //dereferencing operator
      T operator*() const {   //has to be const 
	return Op::apply(*iter1_, *iter2_);
      }
  

      //bracket operator
      T operator[](int i) const{
	return Op::apply(iter1_[i], iter2_[i]);
      }

      //increment iterator
      void operator++(){
	++iter1_;
	++iter2_;
      }
       
      //should be avoided
       void operator++(int){ 
	iter1_++;
	iter2_++;
      }

      //decrement iterator
      void operator--(){ 
	--iter1_;
	--iter2_;
      }

      //should be avoided
      void operator--(int){ 
	iter1_--;
	iter2_--;
      }
    };



    /**
     * \brief Wraps a STL-compliant iterator
     */
    template<class T, class ITER>
    class XWrapIter{
    private:
      ITER iter_;

    public:
      typedef T value_type;

      XWrapIter(const ITER& iter = ITER()):iter_(iter){}

      
      void operator++(){
	++iter_;
      }
    
      //suffix plus 1
      void operator++(int){
	iter_++;
      }

      //!NOTE: only *iter and iter_[i] available for standard iters!     
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      T operator()(int i, T& s, T& c) const{
	return iter_[i];
      }

#endif
     
      T operator*() const{
	return *iter_;
      }

      T operator[](int i) const{
	return iter_[i];
      }	

    };

    
     //UUUUUUUUUUUUUU UNARY operation -- Just for fun UUUUUUUUUUUUUUUUUUUUUUU
    template<typename T, class Op, class U>
    class UnaryOp{
    private:
      U iter_;
    
    public:
      typedef T value_type;
      typedef typename Op::ReturnType ReturnType;
      
      UnaryOp(const U& it):iter_(it){} //constructor


#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      ReturnType operator()(int i, T& s, T& c) const{
	return Op::apply(iter_(i,s,c), s, c); 
      }
#endif
      


      T operator*() const {
	return Op::apply(*iter_);
      }

      T operator[](int i) const{
	return Op::apply(iter_[i]);
      }

      void operator++(){ ++iter_;}

      void operator++(int){ iter_++;}

      void operator--(){ --iter_;}

      void operator--(int){ iter_--;}   
    };

    //APPLICATIVE TEMPLATE CLASSES for unary stuff
    template<typename T>
    class AbsApp{
    public:
      AbsApp(){}
      
      static inline T apply(const T& a){
	return std::abs(a);      //invoke abs on every expression entry
      }
    };


    
    template<typename T>
    class NegativeNumber{
    public:
      typedef T ReturnType;

#if USE_TUNED_FLOATING_POINT_ARITHMETIC 
      static inline T apply(const T& a, const T& s, const T& c){
	return -a;
      }
#else  
      static inline T apply(const T& a){
	return -a;
      }
#endif
    };


    
    template<typename T, class E>
    inline  VExpr<T, UnaryOp<T, NegativeNumber<T>, VExpr<T,E> > > operator-(const VExpr<T,E>& expr){
      typedef UnaryOp<T, NegativeNumber<T>, VExpr<T,E> > UO;

      return VExpr<T,UO>(UO(expr));
    }
    
    
    //! note that this is already covered by <TT> operator-() </TT>, see above
    //! but for the sake of completeness, I added it here...
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    template<typename T>
    inline  VExpr<T, UnaryOp<T, NegativeNumber<T>, XWrapIter<T,typename MyVec<T>::const_iterator> > > operator-(const MyVec<T>& v){
      typedef UnaryOp<T, NegativeNumber<T>, XWrapIter<T,typename MyVec<T>::const_iterator> > UO;
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      
      return VExpr<T,UO>(UO(WrapIterType(v.begin())));
    }
#else

    template<typename T>
    inline  VExpr<T, UnaryOp<T, NegativeNumber<T>, typename MyVec<T>::const_iterator> > operator-(const MyVec<T>& v){
      typedef UnaryOp<T, NegativeNumber<T>, typename MyVec<T>::const_iterator> UO;
      return VExpr<T,UO>(UO(v.begin()));
    }
#endif



    //UUUUUUUU define UNARY OPERATORS/FUNCTIONS here UUUUUUUUUUUUU 
    template<typename T, class A>
    inline VExpr<T, UnaryOp<T, AbsApp<T>, VExpr<T,A> > > abs(const VExpr<T,A>& e ){
 
      typedef UnaryOp<T, AbsApp<T>, VExpr<T,A> > Uop;
      return VExpr<T,Uop>(Uop(e));
    } 



    //UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU



    //!######### ACTUAL OPERATIONS ON EXPRESSIONS DISGUISING MyVec operations ##
    //!NOTE1: MyVec<T>::iterator, std::vector<T>::iterator, etc. do not work as        template parameter as long as you don't add 'typename'.
    //!NOTE2: it is mandatory to use the CONST_ITERATOR as template type.

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    //ADDITION
    //1.) "MyVec + MyVec"
    template<typename T>
    inline VExpr<T, LOR<T, XWrapIter<T,typename MyVec<T>::const_iterator>, Add<T>, XWrapIter<T,typename MyVec<T>::const_iterator > > >
    operator+(const MyVec<T>& v1, const MyVec<T>& v2){
  
      //assert(v1.size() == v2.size());
      if(v1.size() != v2.size()){ 
	ADONIS_ERROR(ExpressVecError, "MyVec sizes do not agree.");
      }

      //abbreviate the Binary Expression
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T, WrapIterType, Add<T>, WrapIterType > ET;
      //std::cout << "MyVec + MyVec" << std::endl;
      return VExpr<T,ET>(ET(WrapIterType(v1.begin()),WrapIterType(v2.begin())));
    }

    //2.) "Expression + MyVec"
    template<typename T, class L>
    inline VExpr<T, LOR<T, VExpr<T,L>, Add<T>,  XWrapIter<T,typename MyVec<T>::const_iterator> > >
    operator+(const VExpr<T,L>& expr, const MyVec<T>& v2){
  
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T, VExpr<T,L>, Add<T>, WrapIterType > ET;
      
      //std::cout << "Expression + MyVec: wi2 = " <<std::endl;

      return VExpr<T,ET>(ET(expr, WrapIterType(v2.begin())));
    }

    //3.) "MyVec + Expression"
    template<typename T, class R>
    inline VExpr<T, LOR<T,XWrapIter<T,typename MyVec<T>::const_iterator>, Add<T>, VExpr<T,R> > >
    operator+(const MyVec<T>& v1,const VExpr<T,R>& expr){
  
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T,WrapIterType, Add<T>, VExpr<T,R> > ET;
      
      //std::cout << "MyVec + Expression" <<std::endl;
      return VExpr<T,ET>(ET(WrapIterType(v1.begin()),expr));
    }


    //SUBTRACTION
    //1.) "MyVec - MyVec"
    template<typename T>
    inline VExpr<T, LOR<T, XWrapIter<T,typename MyVec<T>::const_iterator>, Subtract<T>, XWrapIter<T,typename MyVec<T>::const_iterator> > >
    operator-(const MyVec<T>& v1, const MyVec<T>& v2){
      
      if(v1.size() != v2.size()){
	ADONIS_ERROR(ExpressVecError, "MyVec sizes do not agree.");
      }

      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      //abbreviate the Binary Expression
      typedef LOR<T, WrapIterType, Subtract<T>, WrapIterType >  ET;
      return VExpr<T,ET>(ET(WrapIterType(v1.begin()), WrapIterType(v2.begin())));
    }

    //2.) "Expression - MyVec"
    template<typename T, class L>
    inline VExpr<T, LOR<T, VExpr<T,L>, Subtract<T>, XWrapIter<T,typename MyVec<T>::const_iterator> > >
    operator-(const VExpr<T,L>& expr, const MyVec<T>& v2){
  
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T, VExpr<T,L>, Subtract<T>, WrapIterType> ET;
      return VExpr<T,ET>(ET(expr, WrapIterType(v2.begin())));
    }

    //3.) "MyVec - Expression"
    template<typename T, class R>
    inline VExpr<T, LOR<T,XWrapIter<T,typename MyVec<T>::const_iterator>, Subtract<T>, VExpr<T,R> > >
    operator-(const MyVec<T>& v1,const VExpr<T,R>& expr){
  
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T,WrapIterType, Subtract<T>, VExpr<T,R> > ET;
      return VExpr<T,ET>(ET(WrapIterType(v1.begin()),expr));
    }

    //MIND that it might be that the AScalar only lives within the operator
    //AScalar * MyVec
    template<typename T>
    inline VExpr<T, LOR<T, AScalar<T>, Multiply<T>, XWrapIter<T,typename MyVec<T>::const_iterator> > >
    operator*(const T& s, const MyVec<T>& v){
      //std::cout << "invoke scal*MyVec (use tuned flt.pt. arithmetic)"<<std::endl;

      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T, AScalar<T>, Multiply<T>, WrapIterType> ET;
      return VExpr<T,ET>(ET(AScalar<T>(s), WrapIterType(v.begin())));
    }
    
    //MyVec*AScalar
    template<typename T>
    inline VExpr<T, LOR<T, XWrapIter<T,typename MyVec<T>::const_iterator> , Multiply<T>,AScalar<T> > >
    operator*(const MyVec<T>& v, const T& s){
    
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T, WrapIterType, Multiply<T>, AScalar<T> > ET;
      return VExpr<T,ET>(ET(WrapIterType(v.begin()), AScalar<T>(s)));
    }


     //MyVec / AScalar
    template<typename T>
    inline VExpr<T, LOR<T, XWrapIter<T,typename MyVec<T>::const_iterator> , Divide<T>,AScalar<T> > >
    operator/(const MyVec<T>& v, const T& s){
      if( s == T(0.))
	ADONIS_ERROR(ZeroDivision, "Denominator = " << s <<".");
      
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T, WrapIterType, Divide<T>, AScalar<T> > ET;
      return VExpr<T,ET>(ET(WrapIterType(v.begin()), AScalar<T>(s)));
    }

     //MULTIPLICATION BETWEEN OBJECTS, BOTH NOT BEING A SCALAR
    //1.) "MyVec * MyVec"
    template<typename T>
    inline VExpr<T, LOR<T, XWrapIter<T,typename MyVec<T>::const_iterator>, Multiply<T>, XWrapIter<T,typename MyVec<T>::const_iterator> > >
    operator*(const MyVec<T>& v1, const MyVec<T>& v2){
  
      if(v1.size() != v2.size()){
	ADONIS_ERROR(ExpressVecError, "MyVec sizes do not agree.");
      }
      
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      //abbreviate the Binary Expression
      typedef LOR<T, WrapIterType, Multiply<T>, WrapIterType >  ET;
      return VExpr<T,ET>(ET(WrapIterType(v1.begin()), WrapIterType(v2.begin())));
    }

    //2.) "Expression * MyVec"
    template<typename T, class L>
    inline VExpr<T, LOR<T, VExpr<T,L>, Multiply<T>, XWrapIter<T,typename MyVec<T>::const_iterator> > >
    operator*(const VExpr<T,L>& expr, const MyVec<T>& v2){
  
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T, VExpr<T,L>, Multiply<T>, WrapIterType> ET;
      return VExpr<T,ET>(ET(expr, WrapIterType(v2.begin())));
    }

    //3.) "MyVec * Expression"
    template<typename T, class R>
    inline VExpr<T, LOR<T,XWrapIter<T,typename MyVec<T>::const_iterator>, Multiply<T>, VExpr<T,R> > >
    operator*(const MyVec<T>& v1,const VExpr<T,R>& expr){
  
      typedef XWrapIter<T,typename MyVec<T>::const_iterator> WrapIterType;
      typedef LOR<T,WrapIterType, Multiply<T>, VExpr<T,R> > ET;
      return VExpr<T,ET>(ET(WrapIterType(v1.begin()),expr));
    }


    
//STANDARD OPERATIONS, I.E. WITHOUT TUNED SUMMATION
#else  
     //1.) "MyVec + MyVec"
    template<typename T>
    inline VExpr<T, LOR<T, typename MyVec<T>::const_iterator, Add<T>, typename MyVec<T>::const_iterator > >
    operator+(const MyVec<T>& v1, const MyVec<T>& v2){
  
      //assert(v1.size() == v2.size());
      if(v1.size() != v2.size()){ 
	ADONIS_ERROR(ExpressVecError, "MyVec sizes do not agree.");
      }

      //abbreviate the Binary Expression
      typedef LOR<T, typename MyVec<T>::const_iterator, Add<T>, typename MyVec<T>::const_iterator >  ET;
      return VExpr<T,ET>(ET(v1.begin(), v2.begin()));
    }

    //2.) "Expression + MyVec"
    template<typename T, class L>
    inline VExpr<T, LOR<T, VExpr<T,L>, Add<T>, typename MyVec<T>::const_iterator > >
    operator+(const VExpr<T,L>& expr, const MyVec<T>& v2){
  
      typedef LOR<T, VExpr<T,L>, Add<T>, typename MyVec<T>::const_iterator > ET;
      return VExpr<T,ET>(ET(expr, v2.begin()));
    }

    //3.) "MyVec + Expression"
    template<typename T, class R>
    inline VExpr<T, LOR<T,typename MyVec<T>::const_iterator, Add<T>, VExpr<T,R> > >
    operator+(const MyVec<T>& v1,const VExpr<T,R>& expr){
  
      //std::cout<<"SIZE-VEC = "<<v1.size() <<", EXPR-SIZE = "<<expr.size()<<std::endl;
      
      typedef LOR<T,typename MyVec<T>::const_iterator, Add<T>, VExpr<T,R> > ET;
      return VExpr<T,ET>(ET(v1.begin(),expr));
    }

    //SUBTRACTION
     //1.) "MyVec - MyVec"
    template<typename T>
    inline VExpr<T, LOR<T, typename MyVec<T>::const_iterator, Subtract<T>, typename MyVec<T>::const_iterator > >
    operator-(const MyVec<T>& v1, const MyVec<T>& v2){
      
      if(v1.size() != v2.size()){
	ADONIS_ERROR(ExpressVecError, "MyVec sizes do not agree.");
      }

      //abbreviate the Binary Expression
      typedef LOR<T, typename MyVec<T>::const_iterator, Subtract<T>, typename MyVec<T>::const_iterator >  ET;
      return VExpr<T,ET>(ET(v1.begin(), v2.begin()));
    }

    //2.) "Expression - MyVec"
    template<typename T, class L>
    inline VExpr<T, LOR<T, VExpr<T,L>, Subtract<T>, typename MyVec<T>::const_iterator > >
    operator-(const VExpr<T,L>& expr, const MyVec<T>& v2){
  
      typedef LOR<T, VExpr<T,L>, Subtract<T>, typename MyVec<T>::const_iterator > ET;
      return VExpr<T,ET>(ET(expr, v2.begin()));
    }

    //3.) "MyVec - Expression"
    template<typename T, class R>
    inline VExpr<T, LOR<T,typename MyVec<T>::const_iterator, Subtract<T>, VExpr<T,R> > >
    operator-(const MyVec<T>& v1,const VExpr<T,R>& expr){
  
      typedef LOR<T,typename MyVec<T>::const_iterator, Subtract<T>, VExpr<T,R> > ET;
      return VExpr<T,ET>(ET(v1.begin(),expr));
    }

    //MULTIPLICATION
    //MIND that it might be that the AScalar only lives within the operator
    //AScalar * MyVec
    template<typename T>
    inline VExpr<T, LOR<T, AScalar<T>, Multiply<T>, typename MyVec<T>::const_iterator > >
    operator*(const T& s, const MyVec<T>& v){
  
      typedef LOR<T, AScalar<T>, Multiply<T>,  typename MyVec<T>::const_iterator> ET;
      return VExpr<T,ET>(ET(AScalar<T>(s), v.begin()));
    }
    
    //MyVec*AScalar
    template<typename T>
    inline VExpr<T, LOR<T, typename MyVec<T>::const_iterator , Multiply<T>,AScalar<T> > >
    operator*(const MyVec<T>& v, const T& s){
    
      typedef LOR<T, typename MyVec<T>::const_iterator, Multiply<T>, AScalar<T> > ET;
      return VExpr<T,ET>(ET(v.begin(), AScalar<T>(s)));
    }

    //MyVec / AScalar
    template<typename T>
    inline VExpr<T, LOR<T, typename MyVec<T>::const_iterator , Divide<T>,AScalar<T> > >
    operator/(const MyVec<T>& v, const T& s){
      if( s == T(0.))
	ADONIS_ERROR(ZeroDivision, "Denominator = " << s <<".");
      
      typedef LOR<T, typename MyVec<T>::const_iterator, Divide<T>, AScalar<T> > ET;
      return VExpr<T,ET>(ET(v.begin(), AScalar<T>(s)));
    }


    //MULTIPLICATION BETWEEN OBJECTS, BOTH NOT BEING A SCALAR
    //1.) "MyVec * MyVec"
    template<typename T>
    inline VExpr<T, LOR<T, typename MyVec<T>::const_iterator, Multiply<T>, typename MyVec<T>::const_iterator > >
    operator*(const MyVec<T>& v1, const MyVec<T>& v2){
  
      if(v1.size() != v2.size()){
	ADONIS_ERROR(ExpressVecError, "MyVec sizes do not agree.");
      }
     
      //abbreviate the Binary Expression
      typedef LOR<T, typename MyVec<T>::const_iterator, Multiply<T>, typename MyVec<T>::const_iterator >  ET;
      return VExpr<T,ET>(ET(v1.begin(), v2.begin()));
    }

    //2.) "Expression * MyVec"
    template<typename T, class L>
    inline VExpr<T, LOR<T, VExpr<T,L>, Multiply<T>, typename MyVec<T>::const_iterator > >
    operator*(const VExpr<T,L>& expr, const MyVec<T>& v2){
  
      typedef LOR<T, VExpr<T,L>, Multiply<T>, typename MyVec<T>::const_iterator > ET;
      return VExpr<T,ET>(ET(expr, v2.begin()));
    }

    //3.) "MyVec * Expression"
    template<typename T, class R>
    inline VExpr<T, LOR<T,typename MyVec<T>::const_iterator, Multiply<T>, VExpr<T,R> > >
    operator*(const MyVec<T>& v1,const VExpr<T,R>& expr){
  
      typedef LOR<T,typename MyVec<T>::const_iterator, Multiply<T>, VExpr<T,R> > ET;
      return VExpr<T,ET>(ET(v1.begin(),expr));
    }

    
#endif
    


    //!EXPRESSION-EXPRESSION OPERATIONS REMAIN UNAFFECTED in normal addition as well as compensated addition (and, of course, all other operations defined so far)
    //ADDITION
    //4.) "Expression + Expression"
    template<typename T, class L, class R>
    inline VExpr<T, LOR<T, VExpr<T,L>, Add<T>, VExpr<T,R> > >
    operator+(const VExpr<T,L>& expr1,const VExpr<T,R>& expr2){
  
      typedef LOR<T, VExpr<T,L>, Add<T>, VExpr<T,R> > ET;
      //std::cout << "Expression + Expression" <<std::endl;
      return VExpr<T,ET>(ET(expr1,expr2));
    }

    
    //SUBTRACTION
    //4.) "Expression - Expression"
    template<typename T, class L, class R>
    inline VExpr<T, LOR<T, VExpr<T,L>, Subtract<T>, VExpr<T,R> > >
    operator-(const VExpr<T,L>& expr1,const VExpr<T,R>& expr2){
  
      typedef LOR<T, VExpr<T,L>, Subtract<T>, VExpr<T,R> > ET;
      return VExpr<T,ET>(ET(expr1,expr2));
    }


    //MULTIPLICATION    
    //"AScalar*Expression"
    template<typename T, class R>
    inline VExpr<T, LOR<T, AScalar<T>, Multiply<T>, VExpr<T,R> > >
    operator*(const T& s, const VExpr<T,R>& expr){
      //std::cout << "scal*expr (all versions) "<<std::endl;

      typedef LOR<T, AScalar<T>, Multiply<T>,VExpr<T,R> > ET;
      return VExpr<T,ET>(ET(AScalar<T>(s), expr));
    }

    //"Expression*AScalar"
    template<typename T, class R>
    inline VExpr<T, LOR<T,  VExpr<T,R>, Multiply<T>,AScalar<T> > >
    operator*(const VExpr<T,R>& expr, const T& s){
  
      typedef LOR<T, VExpr<T,R>, Multiply<T>, AScalar<T> > ET;
      return VExpr<T,ET>(ET(expr, AScalar<T>(s)));
    }


    //DIVISION:
    //"Expression / AScalar"
    template<typename T, class R>
    inline VExpr<T, LOR<T,  VExpr<T,R>, Divide<T>,AScalar<T> > >
    operator/(const VExpr<T,R>& expr, const T& s){
      if( s == T(0.))
	ADONIS_ERROR(ZeroDivision, "Denominator = " << s <<".");
      
      typedef LOR<T, VExpr<T,R>, Divide<T>, AScalar<T> > ET;
      return VExpr<T,ET>(ET(expr, AScalar<T>(s)));
    }
      



    //Multiplication operators between "higher" objects than scalars
    //4.) "Expression * Expression"
    template<typename T, class L, class R>
    inline VExpr<T, LOR<T, VExpr<T,L>, Multiply<T>, VExpr<T,R> > >
    operator*(const VExpr<T,L>& expr1,const VExpr<T,R>& expr2){
      
      typedef LOR<T, VExpr<T,L>, Multiply<T>, VExpr<T,R> > ET;
      return VExpr<T,ET>(ET(expr1,expr2));
    }

  


    //----------------- OUTPUT using '<<' ------------------------------------- 
    //defined outside class MyVec but usin MyVec's 'print(...)'
    template<typename T>
    inline std::ostream& operator<<(std::ostream& os, const MyVec<T>& m){
      return m.print(os);
    }
    


    //end namespaces
  }
}


#endif
