#ifndef GLOBAL_FUNCTIONS_AND_ALGORITHMS_ETC_HH
#define GLOBAL_FUNCTIONS_AND_ALGORITHMS_ETC_HH

#include <iostream>
#include <cmath>   
#include <math.h>
#include <cstdlib>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>
//#include<stdlib.h>
#include <sstream>   //for string streams  
#include <fstream>
#include <iomanip>
#include <curses.h> //for printing in color  
#include <typeinfo>

#include "adonisassert.hh"
#include "isclass.hh"
#include "typeadapter.hh"
#include "numerictypechecker.hh"
#include "universalconstants.hh"
#include "fancymessages.hh"
#include "error.hh"
#include "numerictypechecker.hh"

#include "../marcvecmatrix/myfunctions.h"

/*
  DON'T USE 'using namespace' AND DON'T FORGET SPECIFYING THE NAMESPACE, E.G. IN 'std:pair, std::ostream' ETC.. THIS HELPS PREVENTING COMPLICATIONS ;-)
 */


#define LENGTH(a)(sizeof(a)/sizeof(*a)) //don't use it in case of a[]


namespace Adonis{

  template<int N> class ChooseFltRadix;

  template<>
  class ChooseFltRadix<2>{  //that's the usual binary case
  public:
    static const int Value = 2;
  };

  
  template<>
  class ChooseFltRadix<8>{  //qubit case
  public:
    static const int Value = 8;
  };

  /**
   * \brief standard detection of a NaN value; equivalent to IEEE's isnan(val)
   */  
  template<class T>
  inline bool is_NaN(const T& val){
    return (val != val);  
  }


  /**
   * \brief return absolute value of given number (either real number or cmplx)
   *
   *  Should also work for <I> multiprecision types </I>
   * Advantage: my Abs class also works with multiprec. types 
   */
  template<class T>
  inline T Abs(const T& x){
    return ( (x >= T()) ? x : -x);
  }

  //! partial specialisation
  template<class T>
  inline T Abs(const std::complex<T>& x){
    return std::abs(x); 
  }

#if USE_CPPAD
  template<class T>
  inline CppAD::AD<T> Abs(const CppAD::AD<T>& ad){
    return CppAD::abs(ad); 
  }
#endif

 template<class T>
  inline T Sigmoid(const T& t){
    return 1./(1+std::exp(-t));
  }

  //! produces a S-shaped curve continuously reaching from 0 to +1
#if USE_CPPAD
  template<class T>
  inline CppAD::AD<T> Sigmoid(const CppAD::AD<T>& at){
    return 1./(1+CppAD::exp(-at));
  }
#endif

  //! square a number (no matter which type, as long as it possesses a
  //! multiplication operator
  template<class T>
  inline T Sqr(const T& x){
    return x*x;
  }

  template<class T>
  inline bool is_negative(const T& val){
    return ((val < T()) ? true : false);
  }

  template<class T>
  inline bool is_negative(const std::complex<T>& z){
    return false;  //well, a complex number does not obey an ordering like reals
  }

 
  template<class T>
  inline bool IsNan(const T& val){
    return ((val != val) ? true : false);
  }

  template<class T>
  inline bool IsInf(const T& val){
    return ((Abs(val) >= std::numeric_limits<T>::infinity()) ? true : false);
  }

  template<class T>
  inline bool IsInf(const std::complex<T>& val){
    return ((Abs(val) >= std::numeric_limits<T>::infinity()) ? true : false);
  }

#if USE_CPPAD
  template<class T>
  inline bool IsInf(const CppAD::AD<T>& val){
    return ((Abs(val) >= std::numeric_limits<typename TypeAdapter<T>::BaseType>::infinity()) ? true : false);
  }
#endif

  /**
   * \brief If a value is NaN or Inf, then false is returned otherwise true
   */


  template<class T>
  inline bool is_well_defined_value(const T& val){
    if(IsNan(val) || IsInf(val))
      return false;
    else
      return true;
  }


  //! check if RAC or scalar is a well defined value, i.e. neither NaN nor Inf
  template<class C, bool B> class IsWellDefinedValue;

  template<class C>
  class IsWellDefinedValue<C,true>{  //! is container
  public:
    typedef C ReturnType;
    typedef typename C::value_type value_type;

    
    static inline void check(const C& v){
      for(size_t i = 1; i < v.size(); ++i){
	if(!is_well_defined_value(v[i])){
	  ADONIS_ERROR(ValueError,"Corrupt RAC value of type '"<<typeid(typename C::value_type).name() << "' : v["<<i<<"] = "<< v[i] << ".");
	}
      }
    }
    
    static inline bool is_plausible_val(const C& v){
      bool b = true;
      for(size_t i = 1; i < v.size(); ++i){
	if(!is_well_defined_value(v[i])){
	  b = false;
	  break;  //leave loop earlier 
	}
      }
      return b;
    }
  };
  
  template<class C>
  class IsWellDefinedValue<C,false>{ //! is conventional number
  public:
    typedef C ReturnType;
    typedef C value_type;   //! C is already value_type

  
    static inline void check(const C& a){
      if(!is_well_defined_value(a)){
	ADONIS_ERROR(ValueError,"Corrupt SCALAR value of type '"<< typeid(C).name()<< "': a = "<< a << ".");
      }
    }

    static inline bool is_plausible_val(const C& a){
      return ( (!is_well_defined_value(a)) ? false : true);
    }	
  };



  template<class W>
  inline void value_plausibility_check(const W& x){
#ifndef NDEBUG
    IsWellDefinedValue<W,IsContainer<W>::Value>::check(x);
#endif
  }
  
  template<class W>
  inline bool is_plausible_value(const W& x){
    return IsWellDefinedValue<W,IsContainer<W>::Value>::is_plausible_val(x);
  }

  
  //!check \f$ |\mathrm{val}_1| > fac\cdot |\mathrm {val2}| \f$
  template<class T1, class T2>
  inline bool much_larger_than_in_abs(const T1& val1, const T2& val2, const typename TypeAdapter<T1>::BaseType& fac){
    adonis_assert(fac > 1);
    return ((Abs(val1) > fac*Abs(val2)) ? true : false);
  }


  /**
   * \brief convert in to char 
   *
   * \param i representation of ASCII symbol: 0-31     control comms
   *                                          32-127   representable symbols
   *                                          128-255  extended symbols
   *
   * ASCII: 
   */
  inline char ascii(int i){
    return (char)i;     
  }

  inline std::string yes_no(bool b, bool logicalAnswer = false){
    return ( ((logicalAnswer == false) ? ((b==true) ? "yes" : "no") : ((b==true) ? "true" : "false")) );
  }

  template<class T>
  inline T switching_function(const T& x, const typename TypeAdapter<T>::BaseType& a = 1, const typename TypeAdapter<T>::BaseType& b = 0.05){
    //adonis_assert(a > 0.);
    return 1./(1+exp(-a*x));    //b is only a dummy here
    // return 0.5+0.4*tanh((x-a)/b);
  }
  
  template<class VIX>
  inline VIX index_union(VIX& index1, VIX& index2){
    VIX v(index1.size()+index2.size());
    typedef typename VIX::iterator iterator;
    std::sort(index1.begin(),index1.end());
    std::sort(index2.begin(),index2.end());

    iterator it;
    it = std::set_union(index1.begin(),index1.end(),index2.begin(),index2.end(),v.begin());
    v.resize(it-v.begin());
    return v;
  }


  inline bool is_even(int k){
    if(k%2==0)
      return true;
    else
      return false;
  }
  
  /**
   * Converts hex address to int such that a human can easily decipher it ;)
   * \return long int because addresses can be quite large integers
   * USAGE:
   * \code
   double d = 0.75;
   double* ptr = &d;
   std::cout << "integer address = "<< hex2int(ptr) << std::endl;
   //further examples
   std::cout << "3E8 ^= "<< hex2int("3E8") << std::endl;
   std::string sixthirteen = "0x265";
   std::cout << sixthirteen << " ^= "<< hex2int(sixthirteen) << std::endl;
   * \endcode
   */
  template<class ST>
  inline long int hex2int(const ST& hx){
    long int i(0);
    std::stringstream ss;
    ss << hx;
    ss >> std::hex >> i;
    return i;
  }

  inline std::string& proper_file_suffix(std::string& s, const std::string& suffix = ".dat"){
    return (s+=suffix);
  }

  //! efficient binomial coefficient calculation.
  //! source:
  //! <a href="http://de.wikipedia.org/wiki/Binomialkoeffizient">Wiki entry (german)</a>
  size_t binomial_coefficient(size_t n, size_t k){
    if((n == k) || (n == 0) || (k == 0)) return 1;
    if (n < k) return 0;

    size_t res(0);
    
    if(2*k > n){
      res = binomial_coefficient(n,n-k);
    }
    else{
      res = n-k+1;
      for(size_t i = 2; i <=k; ++i){
	res *= (n-k+i);
	res /= i;
      }
    }
    return res; 
  }

  std::string nice_bool_out(const bool george){
    return ((george) ? "TRUE" : "FALSE"); 
  }

  //a small value 
  template<typename T>
  inline T epsilon(const T& value){
    return value*std::numeric_limits<T>::epsilon();
  }

  //returns the machine precision e.g. eps_mech<double>()
  template<typename T>
  inline T eps_mach(){
    return std::numeric_limits<T>::epsilon();
  }

  
  template<class INT>
  inline INT gauss_sum(const INT& n){
    return n*(n+1)/2;
  }

  inline void do_nothing(){}

  template<class ARG>
  inline void do_nothing(const ARG& a){}

  //! resize a STL-compliant container, when it is empty, i.e. its size() is 0
  template<class STLV>
  inline void resize_me_when_empty(size_t n, STLV& v, const typename STLV::value_type& d = typename STLV::value_type()){
    (std::distance(v.begin(),v.end()) == 0) ? v.resize(n,d) : do_nothing();
  }

  /**
   *  \brief Converts a given string to a number of given precision.
   */
  template<class T>
  inline T& convert_str_2_num(T& number, const std::string& s){
    std::istringstream is(s);
    is >> number;
    return number;
  }

  /**
   * \brief Converts a given string to a number of given precision.
   * USAGE: (you MUST explicitly state the type T!)
   * \code
   std::cout << convert_str_2_num("3.14159265358979323846264338327950288419716939937510") << std::endl;
   * \endcode
   */
  template<class T>
  inline T convert_str_2_num(const std::string& s){
    T number;
    std::istringstream is(s);
    is >> number;
    return number;
  }
  

  /**
   * \brief convert string to number
   */
  // template<class T>
  // inline std::string& convert_num_2_str(std::string& s, const T& number){
  //   std::ostringstream os;
  //   os << number
  //   return (s = os.str());
  // }

  
  template<class T>
  inline std::string convert_num_2_str(const T& number){
    //std::string s;
    std::ostringstream oss;
    oss << number;
    //s = oss.str();
    return oss.str(); //s;
  }
    
  /**
   * \brief convert some fancy type to a numeric type 
   * USAGE:
   * \code 
   CppAD::AD<double> ad = 3.75;
   double dc = convert_fancy_num_2_num(ad);
     
   std::complex<double> z(3.,4.);
   CppAD::AD<std::complex<double> > adcplx = z; 
   std::complex<double> zc = convert_fancy_num_2_num(adcplx);
   * \endcode
   *
   * RECOMMENDATION:
   *   Don't use it in that way but rather use function <TT> convert_number </TT>
   */
  template<class CT>
  inline typename TypeAdapter<CT>::BaseType convert_fancy_num_2_num(const CT& ad){
    typedef typename TypeAdapter<CT>::BaseType BaseType;
    return ( (IsNan(ad)) ? std::numeric_limits<BaseType>::signaling_NaN() : ( (IsInf(ad)) ? std::numeric_limits<BaseType>::infinity() : convert_str_2_num<BaseType>(convert_num_2_str(ad))) );
  }


  template<class T>
  inline const T& convert_number(const T& a){
    return a; //just return value
  }
  
  //! note: nan and inf are converted to zero when read in as string, so this little trick creates proper values
#if USE_CPPAD
  template<class T>
  inline T convert_number(const CppAD::AD<T>& a){
    return ( (IsNan(a)) ? std::numeric_limits<T>::signaling_NaN() : ((IsInf(a)) ? std::numeric_limits<T>::infinity() : convert_fancy_num_2_num(a)) );
  }
#endif


  template<class T>
  inline T Floor(const T& d){
    return std::floor(d);
  }
  
  //! applies separately to real and imaginary part
  template<class T>
  inline std::complex<T> Floor(const std::complex<T>& c){
    return std::complex<T>(std::floor(c.real()),std::floor(c.imag()));
  }

  //! T suffices as return type here
#if USE_CPPAD
  template<class T>
  inline T Floor(const CppAD::AD<T>& ad){
    return std::floor(convert_fancy_num_2_num(ad));
  }
#endif
 

  //! the same applies for the ceiling function
  template<class T>
  inline T Ceil(const T& d){
    return std::ceil(d);
  }
  
  //! applies separately to real and imaginary part
  template<class T>
  inline std::complex<T> Ceil(const std::complex<T>& c){
    return std::complex<T>(std::ceil(c.real()),std::ceil(c.imag()));
  }

  //! return type T suffices here
#if USE_CPPAD
  template<class T>
  inline T Ceil(const CppAD::AD<T>& ad){
    return std::ceil(convert_fancy_num_2_num(ad));
  }
#endif
 


  template<class T>
  inline T fix(const T& x){
    return ( (x < T()) ? Ceil(x) : Floor(x) );
  }

  //! dummy. No relation defined for complex numbers and also Abs-fkt not
  //! benefitial here
  template<class T>
  inline std::complex<T> fix(const std::complex<T>& z){
    return Floor(z);
  }
  
  /** \brief Generates linearly spaced vector container. 
   *   \tparam V random access container 
   *
   *   \param d1 left point of container 
   *   \param d2 right point of container
   *   \param n  number of points in the container
   *
   *    \return random access container of type V of size n.
   */
  template<class V>
  inline V linspace(const typename V::value_type& d1, const typename V::value_type& d2,int n){
    
    if(n < 2){
      n = 1;     //set n = 1
      std::cout << n << std::endl;
    }
    
    V v(n);
    
    if(n >=2){
      for(int i = 0; i < n-1; ++i){
	v[i] = d1+i*(d2-d1)/(n-1); 
      }
    }
    v[n-1] = d2;   //include last point

    return v;
  }

  inline char found_meaningful_char_in_str(const std::string& s){
    char c = '#';
    for(unsigned i = 0; i < s.size(); ++i){
      if((!my_function_collection::delimiter(s[i])) && 
	 (!my_function_collection::is_comment(s[i]))){
	c = s[i];
	break;
      }
    }
    return c;
  }

  inline bool is_radix_point(char c){
    return ((c == '.') ? true : false);
  }

  /**
   *  \brief show whitespaces in a string str which will be made visible by a non-whitespace char c (default: '^')
   */
  template<class ST>
  inline void show_whitespaces(const ST& str, const char& c = '^'){
    adonis_assert(!my_function_collection::delimiter(c)); //otherwise you won't see anything 
    for(unsigned i = 0; i < str.size(); ++i){
      if(my_function_collection::delimiter(str[i]))
	std::cout << c;
      else 
	std::cout << str[i]; 
    }
    std::cout<<std::endl;
  }

  /**
   * \brief Whitespaces in a string may be troublesome. Therefore, this function replaces all whitespaces by a designated non-whitespace character (default:'_')
   */
  template<class ST>
  inline void replace_whitespaces(ST& str, const char& c = '_'){
    adonis_assert(!my_function_collection::delimiter(c));
    for(unsigned i = 0; i < str.size(); ++i){
      if(my_function_collection::delimiter(str[i]))
	str[i] = c;
    }
  }
  


  //=========== MATRIX stored as random access container providing [] ==========
  
  /**
   *\brief return the \f$ i \f$th row of a matrix stored as rac. 
   *
   * NOTE: Row storage assumed 
   */
  template<class V, class W>
  inline V& row(unsigned i, V& v, unsigned ncols, const W& m){
    adonis_assert(v.size() == m.size()/ncols);
    unsigned nrows = m.size()/ncols,
      start = i*nrows;
    
    //adonis_assert( i < nrows);
    unsigned c = 0;
    for(unsigned l = start; l < (i+1)*nrows; ++l)
      v[c++] = m[l];

    return v;
  }

  /**
   *\brief return the \f$ i \f$th column of a matrix stored as rac. 
   *
   * NOTE: Row storage assumed 
   */
  template<class V, class W>
  inline V& column(unsigned j, V& v, unsigned ncols, const W& m){
    adonis_assert(v.size() == ncols);
    unsigned nrows = m.size()/ncols;
     
    //adonis_assert( j < ncols);
    unsigned c = 0;
    for(unsigned k = 0; k < ncols; ++k)
      v[c++] = m[k*nrows + j];
      
    return v;
  }



   
  template<class T>
  std::vector<T> 
  inline extract_from_entry(const std::vector<T>& v, int i){
    adonis_assert(i >= 0 && i <= int(v.size()));
 
    std::vector<T> store(int(v.size())-i);
	
    for(int dx = i; dx < int(v.size()); dx++)
      store[dx-i] = v[dx];
	
    return store;
  }


  template<class T, class Container>
  inline Container fill_me_hunk(const T* x, int length){ //for const T*
    Container tofill(length);
    int idx=0;
    for(typename Container::iterator i = tofill.begin(); i != tofill.end(); ++i)
      *i = x[idx++];
    
    return tofill;
  }

  template<class T, class Container>
  inline Container fill_me_hunk(T* x, int length){   //non const argument
    Container tofill(length);
    int idx=0;
    for(typename Container::iterator i = tofill.begin(); i != tofill.end(); ++i)
      *i = x[idx++];
    
    return tofill;
  }

 
  


  //uses merge sort to sort an array of type 'Type' and of size 'length'
  template<typename Type>
  inline std::vector<Type> sort_sequence(Type sordid[], int length){
    int mid = length/2;
    int rest = length - mid;

	

    std::vector<Type> fst(mid), snd(rest), Con(length);
	
    for(int i = 0; i < length; i++){
      if(i < mid)
	fst[i] = sordid[i];
      else snd[i-mid] = sordid[i];
    }

    stable_sort(fst.begin(), fst.end());
    stable_sort(snd.begin(), snd.end());
	
    merge(fst.begin(), fst.end(), snd.begin(), snd.end(), Con.begin());
	
    int idx = 0;
    for(typename std::vector<Type>::iterator it = Con.begin(); it != Con.end(); ++it){
      Con[idx++] = *it;   
	
    }
	
    return Con;
  }



  //another implementation which can be used with different containers providing a const and normal iterator. The contents of these iterators are read into a std::vector and then returned after the sequence has been sorted.
  template<typename Type, class Container>
  inline std::vector<Type> sort_sequence(Container& Co){//use reference as parameter NOT const ref !!
    //void sort_sequence(Container& Co){
    int length = int(Co.size());
    int mid = length/2;
    int rest = length - mid;

	

    std::vector<Type> fst(mid), snd(rest), v(length);

    //works with iterator and const_iterator
    int i= 0;
    for(typename Container::const_iterator it = Co.begin(); it != Co.end(); ++it){
      if(i < mid)
	fst[i] = *it;
      else snd[i-mid] = *it;
      

      i++;
    }
	
    stable_sort(fst.begin(), fst.end());
    stable_sort(snd.begin(), snd.end());
	
    merge(fst.begin(), fst.end(), snd.begin(), snd.end(), v.begin());
       
	
	

    //v should contain the sorted sequence now
    return v;
  }
	

  /*
    template<class Iter, class T>
    Iter fast_find(Iter begin, Iter end, T val)
    {
    // Finds the lower bound in at most log(last - first) + 1 comparisons
    Iter i = std::lower_bound(begin, end, val);

    if (i != end && *i == val)
    return i; // found
    else
    return end; // not found
    }
  */


  //a slightly changed function for STL-containers. 
  //NOTE: if an item is found but is several times contained in the container, only the position of the very first occurence is stored as second argument of the return type 'pair'! This is not really dramatic since the Container has to be sorted first before applying the 'lowerBound' function and therefore if one is interested in the actual number of the found items, one simply calculates the difference between its position and the position of the following entry therein.
  template <class Iter, class T>
  inline std::pair<Iter,int> lowerBound (Iter first, Iter last, const T& value )
  //int lowerBound (Iter first, Iter last, const T& value ) //o.k. test
  {
    Iter it;
  
  
    //instead I use ints
    int count, step; //size_t count,step;
  
  
    count = int(distance(first,last));
  
  
    int alpha = 0, omega = count, mid = 0, mstore = -42;

    /* //original loop  
       while (count>0)
       {
       it = first; step=count/2; advance (it,step);
       if (*it<value)                   // or: if (comp(*it,value)), for the comp version
       { first=++it; count-=step+1;  }
       else count=step;
       }
    */

    while (count>0)
      {
	it = first; 
	step = count/2;  
	advance(it,step);       //same as ' it += step; '
     
	//  mid point: do NOT compute via 'mid =  (alpha + omega)/2;' here!!
	mid = alpha + ((omega-alpha)/2);                
     
	if (*it<value){  //repeat search in top half
	  first = ++it;        
	  count -= step+1;  

	  alpha = mid +1;
	}
	else{    //repeat search in bottom half (for '>')
	  count=step; 
	  
	  omega = mid; //omega = mid -1; //if you use 'mid =  (alpha + omega)/2;'
	  if(*it == value) mstore = mid;
	      
	}
      }
 
    //return first; 
 
    std::pair<Iter,int> withposition(first,mstore);          //return first 
    return withposition;
    //return mstore;
  }

    
  template<typename T>
  inline bool less_than_or_equal(const T& c1, const T& c2){
    if( c1 <= c2) return true;
    else return false;
  }

  template<typename T>
  inline bool greater_than_or_equal(const T& c1, const T& c2){
    if( c1 >= c2) return true;
    else return false;
  }


  //return sqare of a number 
  template<typename T>
  inline T sqr(T x){
    return x*x;
  }

  


  //! if a value is smaller than zero than set it to zero, otherwise return 
  //! argument value reference
  template<class T>
  inline T& de_negativize(T& val){
    return ((val < T()) ? (val = T(0)): val);
  }

  template<class T>
  inline const std::complex<T>& de_negativize(const std::complex<T>& z){
    return z;   //do nothing -- there is no order relation defined on C
  }

  //!check if a value is approximately zero and use Abs function here so
  //! that we can use it for complex numbers as well
  template<class T> 
  inline bool is_zero(const T& check, const typename TypeAdapter<T>::BaseType& fac = 1.e+03){
    typedef typename TypeAdapter<T>::BaseType BaseType;
    adonis_assert(fac > Abs(BaseType()));  //assert positivity of fac
    //something above machine epsilon
    //! Note that numeric_limits<complex> is (0,0), Hence use BaseType
    if((check == T()) || Abs(check) <= fac*Abs(std::numeric_limits<BaseType>::epsilon()))
      return true;
    else return false;
  }

  //! same as is_zero, but might be more intuitive
  template<class T>
  inline bool approximately_zero(const T& num, const T& fac = 1.e+03){
    return is_zero(num,fac);
  }

  //! make entry zero whenever a value (or abs. value) is small
  template<class T>
  inline T zero(const T& z, const typename TypeAdapter<T>::BaseType& fac = 1.e+03){
    return (is_zero(z,fac) ? T() : z);
  }


  template<class T , class D>
  inline bool is_approximately_equal(const T& val, const D& desired, const typename TypeAdapter<T>::BaseType& eps = 1.e-11){
    adonis_assert((eps > 0) && (eps <= 1.e-06));  //only allow sufficiently small value here
    return ( (Abs(val-desired) <= eps) ? true : false);
  }

 

  template<class T>
  inline T Coth(const T& x){
    T sh = sinh(x);
    adonis_assert(!is_zero(sh));
    return cosh(x)/sh;
  }

  template<class T> 
  inline T NaturalPow(const T& x, int e){
    int p = Abs(e);
    T pow(1);

    for(int i = 0; i < p; ++i)
      pow *= x;

    return ((e >= 0) ? pow : T(1)/pow);
  }



  //! not every C++-Standard has a log2-function, let alone, if it exists, 
  //! is defined for complex argument values 
  //! Computation of base b logarithm from k base logarithm: \f[ \log_b(x) = \frac{\log_k(x)}{\log_k(b)}, \f] where, usually, \f$ k = 10\f$ or \f$k = e.\f$
  //! We choose \f$k = e\f$ for convenience.
  template<class T>
  inline T Log2(const T& x){  //! x can also be complex
    adonis_assert(!is_zero(x));
    return std::log(x)/std::log(2.);
  }

#if USE_CPPAD
  template<class T> //use CppAD functions to avoid work in that context
  inline CppAD::AD<T> Log2(const CppAD::AD<T>& x){
    adonis_assert(!is_zero(x));
    return CppAD::log(x)/CppAD::log(2.);
  }
#endif

  template<int N, class T>
  inline T Logb(const T& x){
    adonis_assert(!is_zero(x));
    return ((N==2) ? log2(Abs(x)) : std::log(Abs(x)));
  }

#if USE_CPPAD
  template<int N, class T>
  inline CppAD::AD<T> Logb(const CppAD::AD<T>& x){
    adonis_assert(!is_zero(x));
    return ((N==2) ? Log2(Abs(x)) : CppAD::log(Abs(x)));
  }
#endif

  //!signum function for real numbers 
  template<class T>
  inline int Sgn(const T& x){
    return ((x > T()) ? 1 : ((x < T()) ? -1 : 0));  //(x > T()) - (x < T())
  }

  //!generalization for complex numbers
  template<class T>
  inline std::complex<T> Sgn(const std::complex<T>& z){
    return (Abs(z) > T() ?  z/Abs(z) : std::complex<T>());
  }

  template<class T>
  inline T Reciprocal(const T& x){
    adonis_assert(!is_zero(x));
    return 1./x;
  }

  template<class T>
  inline std::complex<T> Reciprocal(const std::complex<T>& z){
    adonis_assert(!is_zero(Abs(z)));
    return (conj(z)/(z.real()*z.real() + z.imag()*z.imag()));
  }


  template<class T>
  inline T Acosh(const T& x){return std::log(x + std::sqrt(x*x-1.));}
  
  //! cf. <a href="http://mathworld.wolfram.com/InverseHyperbolicCosine.html">Inverse Hyperbolic Cosine </a>
  template<class T>
  inline std::complex<T> Acosh(const std::complex<T>& z){
    return std::log(z+std::sqrt(z+1.)*std::sqrt(z-1.)); //!log's defined for C
  }

  template<class T>
  inline T Asin(const T& x) {return std::asin(x);}

  //! complex asin, cf. german wiki
  template<class T>
  inline std::complex<T> Asin(const std::complex<T>& z){
    T a2 = z.real()*z.real(),
      b2 = z.imag()*z.imag(),
      a2plusb2 = a2+b2,
      sqrtrdx = std::sqrt((a2plusb2 - 1)*(a2 + b2 - 1) + 4*b2);
    return std::complex<T>(Sgn(z.real())/2.*std::acos(sqrtrdx - a2plusb2), Sgn(z.imag())/2.*Acosh(sqrtrdx + a2plusb2));
  }

  template<class T>
  inline T Acos(const T& x) {return std::acos(x);}

  template<class T>
  inline std::complex<T> Acos(const std::complex<T>& z){
    return ( UniversalConstants<T>::Pi/2. - Asin(z) );
  }

  //! checks whether \f$ x \in [a,b] \f$ 
  template<class T>
  inline bool is_contained(const T& x, const T& a, const T& b){
    adonis_assert(a <= b);
    return ((a <= x) && (x <= b));
  }

  template<class T>
  inline bool is_contained(std::complex<T>& z, const T& a, const T& b){
    adonis_assert(a <= b);
    return ((a <= Abs(z)) && (Abs(z) <= b));
  }

  //! checks whether \f$ x \in (a,b) \f$ 
  template<class T>
  inline bool is_properly_contained(const T& x, const T& a, const T& b){
    adonis_assert(a <= b);
    return ((a < x) && (x < b));
  }

  template<class T>
  inline bool is_properly_contained(std::complex<T>& z, const T& a, const T& b){
    adonis_assert(a <= b);
    return ((a < Abs(z)) && (Abs(z) < b));
  }


  //!complex conjugate, for non-complex numbers just return the number
  template<class T> //do nothing, since conj(a) = conj((a,0)) = a
  inline const T& Conj(const T& a){return a;}

  template<class T>
  inline std::complex<T> Conj(const std::complex<T>& z){return std::conj(z);}

  //! to my knowledge, this is a characteristic function of the reals only
  template<class T>
  inline size_t Heaviside(const T& val){
    return (val < T() ? T() : 1);
  }  

  template<class T>
  inline T natural_power(const T& base, int exponent){
    if(exponent == 0)
      return 1;
    
    int e = Abs(exponent);

    T b = 1;   //!for complex number, this is the unit element \f$ 1 \equiv (1,0)\f$
    for(int i = 0; i < e; ++i)
      b *= base;
    
    //! note that the reciprocal of a complex number is given by 1./z as well
    return ( (exponent < 0) ? 1./b : b );
  }

 


  
  /**
   * \brief \f$ \mathrm{sinc}(x):= \frac{\sin(x)}{x} \f$ 
   * Note that sinc(x) has a <I>removable singularity</I>, i.e. if \f$ x= 0\f$ then sinc(0) = 0. This follows from the series expansion of the sine-fct.
   *
   * See [1] FREITAG and BUSAM, "Funktionentheorie I", Springer, p. 130 
   */
  template<class T>  //! real version 
  inline T Sinc(const T& x){
    if(x == T())
      return 1;
    else 
      return sin(x)/x;
  }

  template<class T>  //! complex version 
  inline std::complex<T> Sinc(const std::complex<T>& z){
    if(z == std::complex<T>())
      return std::complex<T>(1,0);
    else 
      return sin(z)/z;
  }



  //! not that reals are a totally ordered set...
  template<class T>
  const T& order_relation(const T& d){return d;}

  //! whereas the complex numbers are not. A preorder is attained by employing
  //! the absolute value 
  template<class T>
  T order_relation(const std::complex<T>& z){
    return Abs(z);
  }

  //! take maximum element of container or scalar
  template<class C, bool B> class MaxElement;

  template<class C>
  class MaxElement<C,true>{  //! is container
  public:
    typedef C ReturnType;
    typedef typename C::value_type value_type;

    static inline C get_maximum(const C& v, const C& w){
      adonis_assert(v.size() == w.size());
      C u;
      u.reserve(u.size()); //reserve some space to avoid reallocation
      for(size_t i = 0; i < v.size(); ++i){
	u.push_back((order_relation(v[i]) > order_relation(w[i])) ? v[i] : w[i]);
      }
      return u;
    }
  
    //! one element
    static inline value_type get_maximum(const C& v){
      value_type mx = v[0];
      for(size_t i = 1; i < v.size(); ++i)
	mx = ((order_relation(mx)>order_relation(v[i])) ? mx : v[i]);
      return mx;
    }
  };
  
  template<class C>
  class MaxElement<C,false>{ //! is conventional number
  public:
    typedef C ReturnType;
    typedef C value_type;   //! C is already value_type

    static inline C get_maximum(const C& a, const C& b){
      return ((order_relation(a) > order_relation(b)) ? a : b);
    }
  
    static inline const C& get_maximum(const C& a){
      return a;
    }
  };


  template<class W>
  inline typename MaxElement<W,IsContainer<W>::Value>::ReturnType Max(const W& x, const W& y){
    return MaxElement<W,IsContainer<W>::Value>::get_maximum(x,y);
  }

  //! overload only one element is returned!
  template<class W>
  inline typename MaxElement<W,IsContainer<W>::Value>::value_type Max(const W& x){
    return MaxElement<W,IsContainer<W>::Value>::get_maximum(x);
  }
  

  // find minimum element
  //! take maximum element of container or scalar
  template<class C, bool B> class MinElement;

  template<class C>
  class MinElement<C,true>{  //! is container
  public:
    typedef C ReturnType;
    typedef typename C::value_type value_type;

    static inline C get_minimum(const C& v, const C& w){
      adonis_assert(v.size() == w.size());
      C u;
      u.reserve(u.size()); // reserve space to avoid reallocation
      for(size_t i = 0; i < v.size(); ++i){
	u.push_back((order_relation(v[i]) < order_relation(w[i])) ? v[i] : w[i]);
      }
      return u;
    }
  
    //! one element
    static inline value_type get_minimum(const C& v){
      value_type mn = v[0];
      for(size_t i = 1; i < v.size(); ++i)
	mn = ((order_relation(mn)<order_relation(v[i])) ? mn : v[i]);
      return mn;
    }
  };
  
  template<class C>
  class MinElement<C,false>{ //! is conventional number
  public:
    typedef C ReturnType;
    typedef C value_type;   //! C is already value_type

    static inline C get_minimum(const C& a, const C& b){
      return ((order_relation(a) < order_relation(b)) ? a : b);
    }
  
    static inline const C& get_minimum(const C& a){
      return a;
    }
  };


  template<class W>
  inline typename MinElement<W,IsContainer<W>::Value>::ReturnType Min(const W& x, const W& y){
    return MinElement<W,IsContainer<W>::Value>::get_minimum(x,y);
  }

  //! overload only one element is returned!
  template<class W>
  inline typename MinElement<W,IsContainer<W>::Value>::value_type Min(const W& x){
    return MinElement<W,IsContainer<W>::Value>::get_minimum(x);
  }
  

  /**
   * \brief Equidistant spacing for a designated interval.
   * \param nodes Number of nodes of interval i.e. nodes-1 subintervals
   * \param x,y left and right end of interval respectively.
   */
  template<class T, class INT>
  inline T interspace(INT nodes, const T& x, const T& y){
    return Adonis::Abs(y-x)/(nodes-1);
  }


  /**
   * \brief "Invert" container by inverting its single elements
   */
  template<class C>
  inline void invert(C& c){
    for(typename C::iterator it = c.begin(); it != c.end(); ++it){
      adonis_assert(!is_zero(*it));
      *it = 1./(*it);
    }
  }


  template<class T>
  inline bool is_different(const T& a, const T& b){
    return ((a != b) ? true : false);
  }

  //specialisation
  template<class T>
  inline bool is_different(const std::complex<T>& z1, const std::complex<T>& z2){
    return (((z1.real() == z2.real()) && (z1.imag() == z2.imag())) ? true : false);//((Abs(z1-z2) != 0) ? true : false);
  }


  //!for reals and complex numbers
  template<class T>
  inline bool is_equal(const T& a, const T& b, const T& fac = 1.e+03){
    if(Abs(a-b) < fac*std::numeric_limits<typename TypeAdapter<T>::BaseType>::epsilon())  //take a value well above round-off
      return true;
    else
      return false;
  }

  //! for random access containers having an operator[]
  template<class V>
  inline bool is_equal(const V& cont){
    unsigned dim = std::distance(cont.begin(),cont.end());
    bool ise = true;
    for(unsigned i = 0; i < dim-1; ++i)
      if(!is_equal(cont[i],cont[i+1])){
	//std::cout << "entry "<< i <<" and" << i+1 << " are not equal. Leave loop." <<std::endl;
	ise = false;
	break;
      }
    return ise;
  }

  //!not usable for complex numbers since they don't possess in the usual sense
  template<class T>
  inline bool is_enclosed(const T& d, const T& a, const T& b){
    return (((is_different(a,b)) && (a <= d) && (d <= b)) ? true : false );
  }


  //! if a value is exactly smaller than another
  template<class T>
  inline bool is_smaller(const T& a, const T& b){
    return ((a < b) ? true : false);
  }


  //!return index in random access container
  template<class T, class RAC>
  inline int find(const T& t, const RAC& r){
    int dim = (int)std::distance(r.begin(),r.end());
    for(int i=0; i < dim; ++i){
      if(r[i] == t){
	return i;
	break;
      }
    }
    return -1;  //failure: no index found
  }

  inline void skull(const std::string& softwarename){

    std::cout<<"********************************************************************************"<<std::endl;
    std::cout<<" *** "<< softwarename <<" FAILURE *** "<< softwarename <<" FAILURE *** "<< softwarename <<" FAILURE *** "<<std::endl;
    //std::cout<<setw(80)<<"*"<<std::endl;
    std::cout<<"********************************************************************************"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"      #######"<<std::endl;
    std::cout<<"    #         #          "<< softwarename <<" HAS ABORTED UNTIMELY..."<<std::endl;
    std::cout<<"  #             # "<<std::endl;
    std::cout<<" #   ##     ##   #       ... POSSIBLE REASON(S):"<<std::endl;
    std::cout<<" #   ##     ##   #       "<<std::endl;
    std::cout<<" #       #       #       (1) SEVERE ERRORS DURING COMPUTATION."<<std::endl;
    std::cout<<"   #    ###    #         (2) POORLY SET OPTIONS."<<std::endl;
    std::cout<<"     #       #           "<<std::endl;
    std::cout<<"      # # # #            "<<std::endl;
    std::cout<<std::endl;      
    std::cout<<"      # # # #            (I resume all subsequent calculations, Marc...)"<<std::endl;
    std::cout<<"       #####"<<std::endl;
    std::cout<<std::endl;         
    //std::cout<<setw(80)<<"*"<<std::endl;
    std::cout<<"********************************************************************************"<<std::endl;
    std::cout<<" *** "<< softwarename <<"FAILURE *** "<< softwarename <<"FAILURE *** "<< softwarename <<"FAILURE *** "<<std::endl; 
    std::cout<<"********************************************************************************"<<std::endl;
    std::cout<<std::endl;
  }

  inline void print_placeholder_2_screen(std::ostream& os, unsigned n, char ph = '*', const std::string& inbetween = "  "){
    for(unsigned i = 0; i < n; ++i)
      os << ph << inbetween;
  }


  //check if a char is a delimiter ('\t', ' ' or ',')
  inline bool is_dlm(const char& c){
    bool gb = false;
    switch(c){
    case ' ':       //space
      gb = true;
      break;
    case '\t':      //tabulator
      gb = true;
      break;
    case ',':      //comma
      gb = true;
      break; 
    default: gb = false; 
    }
    return gb;
  }



  /**
   * \brief check if file is there by using the find command
   */
  template<class S>
  inline int check_4_file(const S& s){
    int y; 
    y = system(("find " + s).c_str());
    return y;  //success in finding: 0. 
  }

  /**
   * concatenate 2 files whose names are given by strings
   * See cat's man page
   */
  inline void cat(const std::string& s1, const std::string& s2, const std::string& combined){
    system(("cat " + s1 + " " +s2 + " > " + combined).c_str());
  }


  /**
   * \brief simple read function to read in from files having only one column.
   * NOTE: if a file has several columns then the second, third, etc. are ignored
   *\tparam given STL-compliant container
   */
  template<class V>
  inline void read_n_assign(const std::string& s, V& w){
    std::ifstream infile(s.c_str(),std::ios_base::in);
    typedef typename V::value_type value_type;
    
    adonis_assert(infile); //file has either been mispelled or is not existant
    adonis_assert(infile.is_open());
    
    std::string line;
    
    while(!infile.eof()){
      getline(infile,line,'\n');
      
      if((!my_function_collection::is_whiteline(line)) && (!my_function_collection::is_comment(line))){
	
	w.push_back(static_cast<value_type>(atof(line.c_str()))); //assign to
	//some other container at the same time
      }  
    }
    
    infile.close();
  }
  

  /**
   * \brief the kronecker delta, i.e. \f[ \delta_{ij} = 1 \Rightleftarrow i = j, \delta_{ij} = 0 \Rightleftarrow i \not= j.\f]
   * Example:
   * \code
   for(int i = 0; i < 5; ++i){
   for(int j = 0; j < 5; ++j){
   std::cout << kronecker_delta<double>(i,j)<<" ";
   }
   std::cout << std::endl;
   }
   * \endcode
   */
  template<class IDX>
  inline IDX kronecker_delta(IDX i, IDX j){
    return (i == j ? 1 : 0);
  }


  /**
   * \brief Solves <I> real </I> quadratic equation \f$ ax^2 + bx + c = 0.\f$
   * 
   * The classical midnight formula, \f$ x_{1,2} = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}\f$ suffers massive cancellation. Therefore, I've implemented a formula which is more efficient both w.r.t. accuracy and cancellation, see. [1, p. 10]
   *
   * References:
   *
   *   [1] [HIGHAM, "Accuracy and Stability of Numerical Algorithms", 2nd ed.,SIAM, 2002]
   *
   * NOTE: Cancellation can, however, still occur in calculating \f$ b^2 - 4ac\f$ for \f$ b^2 \approx 4ac\f$ (the case of nearly equal roots) and due to [1] the 
   only workaround is to use extended precision or some trick tantamount to the use of such a precision.
  */
  template<class T>
  inline std::pair<T,T> midnight_formula(const T& a, const T& b, const T& c){
    T det = b*b  - 4*a*c;
    if(det < T())
      ADONIS_ERROR(InfeasibleError,"b² - 4ac < 0  ==> no real solution");
  
    //! robust version to avoid catastrophic cancellation
    std::pair<T,T> p;
    T q = -0.5*(b + Sgn(b)*std::sqrt(det));
    p.first = c/q;   //the plus-part
    p.second = q/a;  //the minus-part

    // T quot = 2*a;
    // adonis_assert(Abs(quot) > T());

    // 
    // p.first = (-b + std::sqrt(det))/quot;
    // p.second = (-b - std::sqrt(det))/quot;
  
    return p;
  }


  /**
   * \brief transform hexadecimal address of object to decimal integer
   */
  template<class O>
  inline int human_readable_address(const O& address){
    std::string s = my_function_collection::num2str(&address);
    return static_cast<int>(strtol(s.c_str(),0, 16));
  }

  //! don't use my_function_collection:: namespace
  template<class T>
  inline std::string Num2str(const T& number){
    std::stringstream giveback;
    giveback << std::setprecision (PREC)  << number; 
    return giveback.str(); 
  }

  /**
   * \brief returns fracional part of a number.
   * frac(223.5670) = 0.5678
   * NOTE: frac(-223.5678) = 0.4322 since floor(-223.5678) = -224
   */
  template<class T>
  inline T frac(const T& x){
    return (x - std::floor(x));
  }
  
  /**
   * \brief returns 'actual' fractional part of number despite its sign, i.e.
   * actual_frac(-223.5678) = 0.5678
   */
  template<class T>
  inline T actual_frac(const T& x){
    if(Abs(x) < 1)
      return Abs(x);
    else
      return (Abs(x - static_cast<int>(x)));
  }

  /**
   * \brief Rounds to the <B>nearest</B> integer, e.g.
   * round(223.5678) = 224  and round(-223.5678= = -224
   * round(4.367) = 4 and round(-4.367) = -4
   * round(0.73) = 1 and round(-0.73) = -1
   * round(0.45) = 0 and round(-0.45) = 0.
   */
  template<class T>
  inline int round(const T& x){
    int r = static_cast<int>(x);
    if(actual_frac(x) >= 0.5){
      if(r < 0)	--r;
      else ++r;
    }
    return r;
  }

  //! complex version
  template<class T>
  inline std::complex<int> round(const std::complex<T>& z){
    return ( std::complex<int>(round(z.real()),round(z.imag())) );
  }
 
  /**
   * \brief similar to 'round' except for the case that the fraction is equal to 0.5. Then we have, for instance,
   *  round_2_nearest_plus_round_2_even(-17.5) = -18
   *  round_2_nearest_plus_round_2_even(-16.5) = -16
   *  round_2_nearest_plus_round_2_even(14.5) = 14
   *  round_2_nearest_plus_round_2_even(21.5) = 21
   */
  template<class T>
  inline int round_2_nearest_plus_round_2_even(const T& x){
    int r = static_cast<int>(x), 
      hrk;
    T af = actual_frac(x);
    if(af > 0.5){      //round to nearest
      r = round(x);
      if((x > -1) && (x < -0.5))
	r = -1;
    }
    
    if(af == 0.5){    //round to even 
      if(r < 0){
	hrk = r-1;
	(hrk%2 == 0) ? r = hrk : r;
      }
      else{
	hrk = r+1;
	(hrk%2 == 0) ? r = hrk : r;
      }
    }
    return r;
  }

  //! complex version
  template<class T>
  inline std::complex<int> round_2_nearest_plus_round_2_even(const std::complex<T>& z){
    return ( std::complex<int>(round_2_nearest_plus_round_2_even(z.real()),round_2_nearest_plus_round_2_even(z.imag())) );
  }
  
  /**
   * \brief Get dimension of container which does not posses as 'size' function,
   * as in the case of a FieldVector
   */
  template<class C>
  inline size_t container_size(const C& c){
    return std::distance(c.begin(),c.end());
  }


  /**
   * Calculate the sum of a STL-compliant container's entries
   */
  template<class C>
  inline typename C::value_type sum(const C& cont){
    typedef typename C::value_type value_type;
    typedef typename C::const_iterator const_iterator;
    value_type tm = value_type();
    for(const_iterator cit = cont.begin(); cit != cont.end(); ++cit){
      adonis_assert(*cit >= value_type());
      tm += *cit;
    }
    return tm;
  }



  /**
   *  \brief Print STL-compliant containers with information about their size. 
   *  This might be benficial especially if you are dealing with pure STL containers which, up to now, do not possess output operators
   *  \tparam C Container type 
   * \param c Container of type C to be printed
   * \param prec precision of output
   * \param FurtherOutput boolean type; outputs a line with some more or less informative text, including the size of the container (no such output by default)  
   * \param s String denoting the container name (no name by default)
   */
  template<class C>  
  inline void print_all(const C& c, int prec = 9, bool FurtherOutput = false, const std::string& s = std::string()){
    if(FurtherOutput)
      std::cout<< std::endl << "# Printed values of container \""<<s<<"\" of size "<<c.size()<<":" <<std::endl;
  
    for(typename C::const_iterator i = c.begin(); i != c.end(); ++i){
      std::cout << std::setprecision(prec) << *i <<"  ";
    }
    std::cout << std::endl;
  }


  template<class IT>
  inline void print_all(IT s, IT e, int prec = 9, bool FurtherOutput = false, const std::string& str = std::string()){
    if(FurtherOutput)
      std::cout<< std::endl << "# Printed values of container \""<<str<<"\" of size "<< std::distance(s,e) <<":" <<std::endl;
  
  
    for(IT it = s; it != e; ++it)
      std::cout << *it << "  ";
    std::cout << std::endl;
  }

  /**
   * \brief print container to file 
   */
  template<class C>
  inline void print_2_file(const std::string& s, const C& c, size_t prec = 9){
    typedef typename C::const_iterator const_iterator;
    
    std::ofstream off(s.c_str(),std::ios_base::out);

    for(const_iterator it = c.begin(); it != c.end(); ++it){
      off << std::setprecision(prec)<<  *it << '\t';
    }
    off << std::endl;
    
    off.close();
  }

  /**
   * \brief print random access container. It is assumed that the matrix is stored row-wise. Hence the second argument specifies the number of columns by which the rows can be extracted.
   */
  template<class C, class IX>
  inline void print_in_matrix_style(const C& c, const IX& cols, const std::string& s = std::string(), size_t prec = 9, bool matlabstyle = false){
    adonis_assert(std::distance(c.begin(),c.end())%cols == 0);
  
    std::cout<< std::endl << "# Printed values of container \""<<s<<"\" of size "<<c.size()<<" in MATRIX style:" <<std::endl;

    std::string delim("  "), rowbreak, mtxstart, mtxend;
    if(matlabstyle){
      delim = ", ";
      rowbreak = ";";
      mtxstart = "[";
      mtxend = "]";
      if(s.size() == 0)
	std::cout << " A = ";
      else 
	std::cout << s << " = ";
      std::cout << mtxstart;
    }
    
    IX count = 0;
    for(typename C::const_iterator it = c.begin(); it != c.end(); ++it){
      std::cout << std::setprecision(prec)<< *it; //'\t';
      count++;
      if (count == cols){
	if(it == c.end()-1)
	  std::cout << mtxend;
	else
	  std::cout << rowbreak << std::endl;
	count = 0; //reset
      }
      else{
	std::cout << delim;
      }
    }
    std::cout << std::endl;
  }


  /**
   * \brief Output sparsity pattern of a given dense matrix in vector form, stored
   *        in row major order
   */
  template<class V>
  inline void print_matrix_pattern(const V& mtx, std::size_t col, const std::string& zeroentryrepres = "   "){ //last argument: "0  "
    adonis_assert(mtx.size()%col == 0);

    std::size_t row = mtx.size()/col;
    for(std::size_t i = 0; i < row; ++i){
      for(std::size_t j = 0; j < col; ++j){
	std::cout << ((is_zero(mtx[i*col+j])) ? zeroentryrepres : "x  " );
      }
      std::cout << std::endl;
    }
  }

  template<class V>
  inline void print_matrix_pattern_2_file(const V& mtx, std::size_t col, const std::string& outfile, bool transpose = false, int prec = 9){
    adonis_assert(mtx.size()%col == 0);
  
    std::size_t row = mtx.size()/col;
    std::ofstream ofs(outfile.c_str(),std::ios_base::out);
    for(std::size_t i = 0; i < row; ++i){
      for(std::size_t j = 0; j < col; ++j){
	ofs << std::setprecision(prec) << ((transpose==false) ?  (mtx[i*col+j]) : (mtx[i + j*row]) )<< " ";
      }
      ofs << std::endl;
    }
    ofs.close();
  }

  /**
   * \brief Compare sparsity patterns of 2 given dense matrices in vector form, 
   * stored in row major order
   */
  template<class V>
  inline bool compare_matrix_pattern(const V& mtx1, const V& mtx2, std::size_t col){
    adonis_assert(mtx1.size() == mtx2.size());
    adonis_assert(mtx1.size()%col == 0);
  
    bool flag = true;
    std::size_t row = mtx1.size()/col;
    for(std::size_t i = 0; i < row; ++i){
      for(std::size_t j = 0; j < col; ++j){
	if(is_zero(mtx1[i*col+j]) == is_zero(mtx2[i*col+j])) //row major
	  flag = true;
	else{
	  flag = false;
	  break;  //leave loop
	}
      }
    }
    return flag;
  }

  //! print dense matrix which is stored in linear memory
  template<class ITER, class IX>
  inline void print_in_matrix_style(ITER start, ITER end, const IX& cols, const std::string& s = std::string(), size_t prec = 9){
    adonis_assert(std::distance(start,end)%cols == 0);
  
    std::cout<< std::endl << "# Printed values of container \""<<s<<"\" of size "<<std::distance(start,end)<<" in MATRIX style:" <<std::endl;

    IX count = 0;
    for(ITER it = start; it != end; ++it){
      std::cout << std::setprecision(prec)<< std::showpoint << *it << "  "; //'\t';
      count++;
      if (count == cols){
	std::cout << std::endl;
	count = 0; //reset
      }
    }
    std::cout << std::endl;

  }

  /**
   *\brief fill container with default types 
   */
  template<class C>
  inline void set_back(C& c, const typename C::value_type& d = typename C::value_type()){
    for(typename C::iterator i = c.begin(); i != c.end(); ++i)
      *i = d;
  }



  //! transform radians to degrees
  template<class T>
  inline T rad2deg(const T& r){
    adonis_assert(!NumericDataTypeChecker<T>::IsComplex);
    adonis_assert(!NumericDataTypeChecker<T>::IsIntegralType);//T != int, short,..
    return r*(180./UniversalConstants<T>::Pi);
  }

  //! transform degrees to radians
  template<class T>
  inline T deg2rad(const T& d){
    adonis_assert(!NumericDataTypeChecker<T>::IsComplex);
    adonis_assert(!NumericDataTypeChecker<T>::IsIntegralType);
    return d*UniversalConstants<T>::Pi/180.;
  }

  /**
   * \brief resets container back to hold default values. This is similar to the above routine 'set_back' but might be more expensive since both clear and resize are linear on the number of elements
   */
  /*template<class C>
    inline void reset(C& c, const typename C::value_type& d = const typename C::value_type()){
    unsigned dim = c.size();
    c.clear();  //all elements of the container are droppend, leaving c with size 0
    c.resize(dim,d);  //resize it again with old size but filled with default vals
    }*/

  /**
   * \brief Canonical base vector \f$ e_i := [0,0,...,0,1,0,...0]^T \f$ whose \f$i \f$-th entry is one
   */
  template<class V, class INT>
  inline V& canonical_base_vector(V& v, const INT& i, bool filledWithZeros = false ){
    adonis_assert(i >= 0 && i < static_cast<INT>(v.size()));

    if(!filledWithZeros)
      set_back(v); //fill with zeros, unless it contains already zeros 
    v[i] = 1;

    return v;
  }

  template<class V>
  inline V canonical_base_vector(size_t i, size_t n){
    adonis_assert(i < n);
    V v(n);
    v[i] = 1.;
    return v;
  }

  /**
   *\brief Update diagonal of a square matrix stored as STL-compliant vector
   *
   *\code
   //updates diagonal by an additive value of 3.5
   update_diagonal<AddBasicElements<T> >(mtx,n,3.5);
   *\endcode
   */
  template<template<class T> class OP, class V, class INT>
  inline void update_diagonal(V& v, INT order, const typename V::value_type& d){
    // std::cout << "v.size() = "<< v.size() << "   order² = "<< order*order << std::endl;
    adonis_assert(static_cast<INT>(std::distance(v.begin(),v.end())) == order*order);
 
  
    for(INT i = 0; i < order; ++i)
      OP<typename V::value_type>::apply(v[i*(order+1)],d);
  }


  template<class T>
  inline T& away_from_zero(T& d, const typename TypeAdapter<T>::BaseType& eps = 1.e-13){
    adonis_assert(eps <= 1.e-12 && eps > typename TypeAdapter<T>::BaseType()); //sufficiently small
    return (d += eps);
  }

  /**
   * \brief Perturb the elements of a container
   */
  template<class C>
  inline C& perturb(C& c, const typename C::value_type& eps){
    for(typename C::iterator it = c.begin(); it != c.end(); ++it)
      *it += eps;
    return c;
  }

  //! with scalar eps
  template<template<class D> class OP, class C>
  inline C apply_perturbation(const C& c, const typename C::value_type& eps){
    typedef typename C::value_type value_type;
    C ptb(c.size());
    for(size_t i = 0; i < c.size(); ++i)
      ptb[i] = OP<value_type>::apply(c[i],eps);
    return ptb;
  }

  //!with vector-valued eps
  template<template<class D> class OP, class C>
  inline C apply_perturbation(const C& c, const C& epsvec){
    typedef typename C::value_type value_type;
    C ptb(c.size());
    for(size_t i = 0; i < c.size(); ++i)
      ptb[i] = OP<value_type>::apply(c[i],epsvec[i]);
    return ptb;
  }

  template<class T>
  inline T perturb_when_zero(const T& x, const typename TypeAdapter<T>::BaseType eps = 1.e-16){
    return ( (is_zero(Abs(x))) ? x + eps : x );
  }

  template<class T>
  inline T perturb_when_zero(const std::complex<T>& z, const typename TypeAdapter<T>::BaseType eps = 1.e-16){
    return ( (is_zero(Abs(z))) ?  std::complex<T>(z.real() + eps, z.imag() + eps) : z);
  }

  /**
   * \brief extract reduced vector from given random access container, assumed the index is given
   */
  template<class V, class W, class IX>
  inline void extract_red(V& v, const W& w, const IX& index){
    adonis_assert(v.size() != 0 && v.size() == static_cast<size_t>(std::distance(index.begin(),index.end())) && v.size() <= w.size());
    unsigned rdim = static_cast<unsigned>(std::distance(index.begin(),index.end()));
    for(unsigned i = 0; i < rdim; ++i)
      v[i] = w[index[i]];
  }

  /**
   * \brief fill zM' entries with reduced composition due to the rpv indices 'idx'
   */
  template<class V, class W>
  inline void replace_with_red(V& zM, const V& r, const W& idx){
    unsigned rdim = static_cast<unsigned>(std::distance(idx.begin(),idx.end()));
    adonis_assert(static_cast<unsigned>(std::distance(zM.begin(),zM.end())) >= static_cast<unsigned>(std::distance(r.begin(),r.end())) && rdim == static_cast<unsigned>(std::distance(r.begin(),r.end())));

    for(unsigned i = 0; i < rdim; ++i){
      zM[idx[i]] = r[i];
    }
  }

  /**
   * \brief Build zM form represented and unrepresented values
   * \param zM manifold point zM
   * \param fullDim the full dimension of the system
   * \param index the index of the rpvs
   * \param unrep the index of the unrepresented species (a.k.a. non-rpvs)
   * \param r rpvs
   * \param u non-rpvs
   */
  template<class IR, class UR, class V>
  inline void build_zM(V& zM, const int& fullDim, const IR& index, const UR& unrep, const V& r, const V& u){
  
    int rdim = static_cast<int>(std::distance(index.begin(),index.end())),
      udim = static_cast<int>(std::distance(unrep.begin(),unrep.end()));
  
    adonis_assert((int)rdim+udim == fullDim);
  
    if((int)std::distance(zM.begin(),zM.end()) != fullDim)
      zM.resize(fullDim);
  
    replace_with_red(zM,r,index);       

    for(int i = 0; i < udim; ++i)
      zM[unrep[i]] = u[i];
  }


  /**
   *\brief extract a solution from a given matrix stored in a single vector
   */
  template<class V>
  inline V extract_solution(const V& B, int m, int n, int nrhs){
    adonis_assert(m >= nrhs && nrhs != 0 && m%nrhs == 0);
    //int n = m-nrhs; //number of columns
    
    V solution(n*nrhs);
  
    int k = 0;
    for(int i = 0; i < nrhs; ++i){
      for(int j = 0; j < n; ++j){
	solution[k++] = B[j + i*m];
      }
    }
    return solution;
  }




  /**
   * \brief \f$\ell_2\f$ error of full and reduced solution with index vector 
   */
  template<class V, class IV>
  inline typename TypeAdapter<typename V::value_type>::BaseType l_2_error(const V& full, const V& red, const IV& index){
    //adonis_assert((int)std::distance(index.begin(),index.end()) == (int)red.size());
    typedef typename TypeAdapter<typename V::value_type>::BaseType BaseType;
    
    BaseType val = BaseType();
    for(unsigned i = 0; i < red.size(); ++i){
      val += ((Abs(full[index[i]] - red[i]))*(Abs(full[index[i]] - red[i])));
    }
    adonis_assert(val > BaseType());
    return std::sqrt(val);
  }

  /*
   * \brief Check whether first argument is smaller,equal or greater than second
   * and lauch some colored output to screen.
   */
  template<class T>
  void is_faster(const T& redtime, const T& fulltime){
    FancyMessages FM;
    std::cout << "RUNNING TIME: ";
    if(redtime < fulltime){
      FM.nice_output("FASTER :)", 32);
      FM.nice_output(" speed-up: "+Num2str(fulltime/redtime), 32);
    }
    if(is_zero(Abs(redtime-fulltime)))
      FM.nice_output("approx. EQUAL :/",33);
    if(redtime > fulltime)
      FM.nice_output("SLOWER :(", 35);
    std::cout << std::endl;

  }


  template<class T>
  inline T slow_down_rate(const T& rtime, const T& ftime){
    if(is_zero(ftime)){
      ADONIS_ERROR(Information, "'ftime' is approx. zero sec. \n   No ratio can be calculated.");
    }

    T ratio = rtime/ftime;

    if(ratio > 1){
      printf("%c[%d;%d;%dm", 0x1B, 1,35, 40);  //violet
      std::cout << "SLOW DOWN :-("<<std::endl;
      printf("%c[%dm", 0x1B, 0);
    }
    else{
      printf("%c[%d;%d;%dm", 0x1B, 1,32, 40); //green 
      std::cout << "SPEED UP  :-)  "<< 1./ratio<<" times faster than before."<<  std::endl;
      printf("%c[%dm", 0x1B, 0);
    }
    return ratio;
  }

  template<class T>
  inline T time_ratio(const T& rtime, const T& ftime){
    if(ftime == 0.)
      ADONIS_ERROR(DerivedError, "Time isn't measurable, i.e. = 0 sec");
    std::cout << "time1 = "<<rtime << "  time2 = "<< ftime << std::endl;

    return rtime/ftime;
  }


  template<class IT> 
  inline size_t reduced_index(size_t i, const IT& it){
    return it[i];
  }

  /**
   * \brief Approximations to Dirac \f$ \delta\f$ distribution which is defined by
   * \f$ \delta(x) = +\infty, x = 0\f$ and \f$ \delta(x) = 0, x\not= 0.\f$ 
   *
   * Recall that \f$ \delta(x) = \lim_{\sigma \rightarrow 0} = \frac{1}{\sqrt{2\pi}\sigma}\exp(-x^2/(2\sigma^2)) \f$
   * Another approximation is the Lorentzian, i.e. \f$ \delta(x) = \lim_{\epsilon \rightarrow 0} = \frac{1}{\epsilon} \frac{\epsilon}{x^2 + \epsilon^2}\f$
   *
   * NOTE: specializations 'g' and 's' might yield bad result, since exp and sin functions are applied
   */
  template<char C> class DiracDeltaDistribution;

  template<>
  class DiracDeltaDistribution<'g'>{ //! Gaussian approx.
  public:
    template<class X>
    static inline X approximation(const X& x, const typename TypeAdapter<X>::BaseType& eps) {
    
      return 1./(sqrt(2.*M_PI)*eps)*exp(-(x*x/(2.*eps*eps)));
    }
  };

  template<>
  class DiracDeltaDistribution<'l'>{   //! Lorentzian approx.
  public:
    template<class X>
    static inline X approximation(const X& x, const typename TypeAdapter<X>::BaseType& eps) {
      return 1./(M_PI)*(eps/(x*x + eps*eps));
    }
  };

  template<>
  class DiracDeltaDistribution<'s'>{   //! Sinc approx.
  public:
    template<class X>
    static inline X approximation(const X& x, const typename TypeAdapter<X>::BaseType& eps) {
      return 1./eps*Sinc(x/eps);
    }
  };


  //! show me some information concerning derivative calculations
  inline void which_derivative_is_used(){
#if USE_CPPAD
    FancyMessages().nice_output("\n Derivatives obtained via CPPAD.\n",36,40);
#else
    FancyMessages().nice_output("\n Derivatives obtained via FINITE DIFFERENCES.\n",34,40);
#endif
  }


  //! compare two types for similarity
  template<class T1, class T2>
  inline bool are_types_equal() {
    return (typeid(T1)==typeid(T2));
  }

  //! checks if a file exists. The file is automatically closed at the end of
  //! the function scope
  //! source:
  //! <a href="http://www.cplusplus.com/forum/general/1796/">Original function</a>
  inline bool does_file_exist(const std::string& fname){
    std::ifstream ifile(fname.c_str(), std::ios_base::in);
    return ifile;
  }

  template<class INT>
  inline std::string ordinal_number(INT i){
    //adonis_assert(NumericDataTypeChecker<INT>::IsIntegralType == true);
    std::string str;
    if(is_equal(Abs(i),(INT)1))
      str="st ";
    else if(is_equal(Abs(i),(INT)2))
      str="nd ";
    else if(is_equal(Abs(i),(INT)3))
      str="rd ";
    else
      str="th ";
    return str;
  }

}//end namespace

#endif
