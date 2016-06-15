/*CAUTION: header file (.h) of TEMPLATE includes member definitions as well. Otherwise there arise problems: Templates are invoked at compilation time and thus  the compiler needs to know the complete classes including their template type(s) which isn't guaranteed during compilation. This is a known drawback of most compilers and can only be avoided as done below [i.e. not seperating template class declaration and member definitions]!!!*/

#ifndef MAVEC_H
#define MAVEC_H

#include "param.h"  //also o.k. !!
#include "myfunctions.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>  //for output precisions
#include <stdlib.h> 
#include <time.h>
#include <limits>  //numeric_limits<double>::epsilon(), etc. 
#include <list>     //for storing object of different size, e.g. matrices


using namespace my_function_collection;

template<class T> class vec {
protected: int length;                    //how long's the vector ?
           T* entry;                     //store the vector's entries  

public: vec(int,  T* );               //constructor 
	 vec(int = 0, T = T());         //default constructor (before 'T = 0')
        vec(const vec&);            //copy constructor
        ~vec(){delete[] entry;}    //destructor
  
  int size() const {return length;}    
  vec& operator=(const vec&);        //copy assignment:overload =
  T& operator[](int ) const;//{return entry[i];}
  vec& operator+=(const vec&);       //overload +=
  vec& operator-=(const vec&);       //overload -=
  vec& operator*=(const vec&);       //overload *=
  vec& operator*=(const T&);       //overload *=
  vec& operator/=(const vec&);       //overload /=
  vec& operator/=(const T&);       //overload /=
 

  bool operator!=(const vec&);       //overload !=
  bool operator==(const vec&);       //overload ==
  
  void resize(int, T val = T());
  
  T euknorm() const;
  T onenorm() const;
  T maxnorm() const;

  T maxelement() const; //finds largest element in a vector
  T minelement() const; //finds smallest element in a vector

  void output() const;                   //outputs vector
  std::ostream& print(std::ostream&) const;  //output-aternative
  
  void print_in_matlab_style(const std::string name) const{
    std::cout << name << " = [";
    std::string delim = ", ";
    for(int i = 0; i < length; ++i){
      if(i == length)
	delim = "";
      std::cout << entry[i] << delim;
    }
    std::cout << "]" << std::endl;
  }
  vec<T> univec(int pos, int n);                //ith-unity vector of size n
  
  // two definitions of an inner product 
  template<class S> friend S inner(const vec<S>&, const vec<S>&); 
  template<class S> friend S operator*(const vec<S>&, const vec<S>&); 
  
  template<class S> friend vec<S> operator+(const vec<S>&, const vec<S>&); 
  template<class S> friend vec<S> operator-(const vec<S>&, const vec<S>&);
  template<class S> friend vec<S> operator-(const vec<S>&);
  template<class S> friend vec<S> operator*(S, const vec<S>&);
  template<class S> friend vec<S> operator*(const vec<S>&, S);
  template<class S> friend vec<S> operator/(const vec<S>&, S);
  template<class S> friend vec<S> cross(const vec<S>&, const vec<S>&); //cross product - only defined in R^3 

  //template<class S> friend vec<S> house(vec<S>&);
 
  T meanvalue() const;
  T variance () const;  //un-biased variance
  T covar(const vec<T>&) const;
  T standard_deviation() const;

  T prodvalue() const;
  T sumup() const;

  void set2zero() const;
 

  //sort with mergesort (O(nlog(n)))
  void merge(int, int, int) const;
  void stablemerge(int, int, int) const;
  
  void stablemerge_with_index(int, int, int) const;

  void mergesort(int, int) const;
  void mergesort() const;
  void shellsort() const;
  void randomvector(int);

  vec<int> position(T);
  
  vec<vec<T> > find_different_entries();
  // vec<T> find_different_entries();

  void write_2_file(const std::string&) const;
  vec<T> reverse();

  int find(T) const;

  int search(T,int, int) const;
  int search(T) const;
  bool different() const;
 
  vec<int> permutation();
  
  vec<T> canonical_basis(int);
  vec<T> subvec(int, int);
}; 

//several constructors for several different constructions of object "vec"
template<class T>                            //constructor
vec<T>::vec(int n, T* a) {
    try{
	entry = new T [length =n]; 
    }
    catch(std::bad_alloc){
	error("Unfortunately not enough memory available to create object >>vec<T>::vec(in, T*)<<.");
    }
    for (int i = 0; i < length; i++)  entry[i]= *(a +i); //same as a[i];
}


template<class T>                            //default constructor
vec<T>::vec(int n, T a) {
  try{
      entry = new T [length =n]; 
  }  
  catch(std::bad_alloc){
      error("I'am upset since I have not enough memory available to create object >>vec<T>::vec(in, T)<<.");
  }
  for (int i = 0; i < length; i++)  entry[i]= a;
}


template<class T>                          //copy constructor: later vec<T> v; vec<T> x = v; id est: piecewise copy
vec<T>::vec(const vec & ve) {
    try{
	entry = new T [length = ve.length]; 
    }
    catch(std::bad_alloc){
	error("Damn it ! I have not enough memory available to apply copy-constructor>>vec<T>::vec(const vec&)<<.");
    }
    for (int i = 0; i < length; i++)  entry[i] = ve[i]; 
}

template<class T> 
vec<T>& vec<T>::operator=(const vec& ve) { //copy assignment
  if (this != &ve) {                       //no self assignment (eg ve = ve)
    /*if (length  != ve.length ) error("bad vector sizes when using EQUALITY operator >> = <<.  ");
   for (int i = 0; i < length; i++) entry[i] = ve[i];
   
    */
   
//........that's only some text that can be neglected.......
    //if(length > 0 && length != ve.length) std::cout<<"CAUTION: Assignment (>>vec<T>::operator= <<) of NON-DIMENSIONLESS vectors with different sizes may lead to unexected - even ERRONEOUS - results ! "<<std::endl;  
  
    //adjust length of current vector to length of argument vector ve
    try{
	if(length > ve.length) delete[] (entry + ve.length);
	if(length < ve.length){
	    delete[] entry;
	    entry = new T[ve.length];
	}
    }
    catch(std::bad_alloc){
	error("Damn it ! I have not enough memory available to apply assignment>>vec<T>::operator=(const vec&)<<.");
    }    
    
    for(int i = 0; i < ve.length; i++) entry[i] = ve.entry[i];
    length = ve.length;

  }

  return *this; //points to the object for which a member fct is invoked, id est " *this " is the object. (returns reference to the object). Tus "this" represents represents a pointer whose member fct is being executed (It checks if a parameter passed to the member fct is the object itself.
}  

template<class T>
T& vec<T>::operator[](int i) const{
    if( i < 0 || i >= length){
	error("Index i = '"+num2str(i)+"' out of range (since vec.size() = "+num2str(length)+") when using >>vec<T>::operator[]<<. To prevent erroneous results further computations are stopped !");
    }
   return entry[i];
}


template<class T> 
vec<T> & vec<T>::operator+=(const vec& ve) {
  if (length != ve.length ) error("bad vector sizes while applying >> += <<" );
  for (int i = 0; i < length; i++) entry[i] += ve[i];
  return *this;
}
 
template<class T> 
vec<T> & vec<T>::operator-=(const vec& ve) {
  if (length != ve.length ) error("bad vector sizes");
  for (int i = 0; i < length; i++) entry[i] -= ve[i];
  return *this;
}

template<class T> 
vec<T> & vec<T>::operator*=(const vec& ve) {    //Vectormultiplication: WEIRD
  if (length != ve.length ) error("Bad vector sizes while applying >> *= <<");
  for (int i = 0; i < length; i++) entry[i] *= ve[i];
  return *this;
}

template<class T> 
vec<T> & vec<T>::operator*=(const T& a) {    //Vector-scalar multiplication
  for (int i = 0; i < length; i++) entry[i] *= a;
  return *this;
}


template<class T> 
vec<T> & vec<T>::operator/=(const T& a) {//Vector-scalar division
  if (a == 0.0) error("Division by 0 occured in >>vec<T>::operator/=(const T&)<<. ");
  for (int i = 0; i < length; i++) entry[i] /= a;
  return *this;
}

// some overloadings of boolean operators: Test if vectors are equal (entries,  rows number, column number ). CAUTION: When outputting an expression of the form 'v1 == v2', never forget to use parentheses, since '<<' has higher priority than '==', i.e.: std::cout<< (v1 == v2) <<std::endl;   (NEVER: std::cout<< v1 == v2 <<std::endl; !!!)

template<class T> 
bool vec<T>::operator==(const vec& logic) {
  //the next line is optional ! 
 if(this == &logic){return true;}
  //the following is needed ! otherwise it outputs false at the end !! 
 if(length != logic.length) return false;
  for(int i =0; i < length; i++){
    if(entry[i] != logic[i]) return false;
  }  
  return true;
}


template<class T> 
inline bool vec<T>::operator!=(const vec& logic) {
  return ! operator==(logic);
  //if(this != &logic) {return false;}
}

template<class T>
void vec<T>::resize(int nouveau, T val){
  
  if(nouveau != length){
      vec<T> aux = *this;
      int alen = length;     

//from here uncomment the next 2 lines while uncommenting the 3rd one
      //maybe look at 'operator=' how itcould work
      delete[] entry;
      try{
	  entry = new T[length = nouveau];
      }catch(std::bad_alloc){
	  error("Damn it ! I have not enough memory available to apply assignment >>vec<T>::resize(int, T=T())<<.");
      }
     
      if(nouveau > alen){
	  for(int i = 0; i < nouveau; i++){
	      if(i < alen) entry[i] = aux[i];
	      else entry[i] = val; 
	  }
      }
      if(nouveau < alen){
	  std::cout<<"ATTENTION: You will probably LOOSE significant data if you use >>vec<T>::resize<< with a smaller length as before..... "<<std::endl;
	  
	  for(int i = 0; i < length; i++){
	      entry[i] = aux[i];
	  }
      }
  }
}

  
  


template<class T> T vec<T>::euknorm() const{        //2-Norm
  T no = entry[0]*entry[0];
  for(int i=1; i < length; i++) no += entry[i]*entry[i];
  return sqrt(no);
}

template<class T> T vec<T>::onenorm() const{        //1-Norm
  T no = abs(entry[0]);
  for(int i=1; i < length; i++) no += abs(entry[i]);
  return no;
}

template<class T> T vec<T>::maxnorm() const{        //max(infinity)-Norm
  T no = abs(entry[0]);
  for(int i=1; i < length; i++) no = std::max(no, abs(entry[i]));
  return no;
}


template<class T> T vec<T>::maxelement() const{
  T res = entry[0];
  for(int i = 1; i < length; i++) res = std::max(res, entry[i]);
  return res;
}

template<class T> T vec<T>::minelement() const{
  T res = entry[0];
  for(int i = 1; i < length; i++) res = std::min(res, entry[i]);
  return res;
}


template<class T>                 //a friend
  T inner(const vec<T>& u, const vec<T>& v){
  if(u.length != v.length) error(" Ooops ! Check your vector sizes, Marc ! ");
  T prodsum = u[0]*v[0];
  for(int i=1; i< u.length; i++) prodsum += u[i]*v[i];
  return prodsum;
}

template<class T>  //inner product with operator*
T operator*(const vec<T>& u, const vec<T>& v){
  if(u.length != v.length) error(" Ooops ! Check your vector sizes, Marc ! ");
  T prodsum = u[0]*v[0];
  for(int i=1; i< u.length; i++) prodsum += u[i]*v[i];
  return prodsum;
}

template<class T>
 vec<T> cross(const vec<T>& x, const vec<T>& y){
  if((x.length != y.length) || (x.length != 3 && y.length != 3)) error("Cross Product is not defined !");
  vec<T> cro(3);
  cro[0] = x[1]*y[2]-x[2]*y[1];
  cro[1] = x[2]*y[0]-x[0]*y[2];
  cro[2] = x[0]*y[1]-x[1]*y[0];
  return cro;
}

template<class T>                      
inline void vec<T>::output() const {
  for(int i = 0; i < length; i++)
    std::cout<< entry[i] << std::endl;vec<T>::length;
}

template<class T>  //a friend
  vec<T> operator+(const vec<T>& u, const vec<T>& v){
  if(u.length != v.length) error(" Ooops ! Check your vector sizes, Marc ! ");
  vec<T> sum(u.length);
  for(int i = 0; i <u.length; i++) sum[i]=u[i]+v[i];
  return sum;
  
//alternativ (possible because of copy constructor):
  /*vec<T> sum = u;
  sum += v;
  return sum;*/
}


template<class T>  //a friend
  vec<T> operator-(const vec<T>& u, const vec<T>& v){
  if(u.length != v.length) error(" Ooops ! Check your vector sizes, Marc ! ");
  vec<T> diff(u.length);
  for(int i = 0; i <u.length; i++) diff[i]=u[i]-v[i];
  return diff;
}
   


template<class T>  //a friend
  inline vec<T> operator-(const vec<T>& u){
  for(int i = 0; i <u.length; i++) u[i] = -u[i];
  return u;
}

template<class T>  // a friend: scalar-vector multiplication
vec<T> operator*(const vec<T>& u, T scal){
  vec<T> scalvec(u.length);
    for(int i = 0; i <u.length; i++) scalvec[i] = scal*u[i];
  return scalvec;
  }

template<class T>  // a friend: vector-scalar  multiplication
vec<T> operator*(T scal, const vec<T>& u){
  vec<T> scalvec(u.length);
    for(int i = 0; i <u.length; i++) scalvec[i] = u[i]*scal;
  return scalvec;
  } 

template<class T>  // a friend: vector-scalar  division
vec<T> operator/(const vec<T>& u, T scal){
  if(scal == 0.0) error("Division by ZERO has occured while executing vec<T> operator/(const vec<T>&, T)."); 
  vec<T> scalvec(u.length);
    for(int i = 0; i <u.length; i++) scalvec[i] = u[i]/scal;
  return scalvec;
  } 

template<class T>
vec<T> vec<T>::univec(int pos, int n) {
  for(int k=0; k< n; k++) {
    if(k == pos-1) entry[k] = 1.;         //pos-1 coz one always counts beginnig with 0
    else  entry[k]=0.;
       }
   }

template<class T> 
std::ostream& vec<T>::print(std::ostream &ou) const {
  for(int i=0; i<length; i++){ 
    ou <<"   "<<std::left<<std::setw(BROAD)<<std::setprecision(VAL)<< entry[i] << std::endl;
 //ou << std::std::endl;
  }
 return ou;
} 

template<class T>  //defined OUTSIDE the class !!
std::ostream& operator<<(std::ostream &ou, const vec<T>& v){
  return v.print(ou);
}


template<class T> //computes the empirical expectation 
T vec<T>::meanvalue() const{
 if(length == 0) return 0; 
 T mean = entry[0];
  for(int i = 1; i < length; i++){
    mean += entry[i];
  }
    
  return (mean/(T(length))); //T(...) shall erase problems when deviding by int  
}

template<class T> //computes the variance 
T vec<T>::variance() const{
  if(length == 0) return 0; 
  if(length == 1) return 1;  //normally it is not defined for one dim
  T expect = (*this).meanvalue();
  T var = (entry[0] - expect)*(entry[0] - expect);
  for(int i=1; i < length; i++){
      var += ((entry[i] -expect)*(entry[i] -expect));
  }
  return (var/(T(length-1)));
}

template<class T>
T vec<T>::standard_deviation() const{
    return (sqrt((*this).variance()));
}

template<class T>
T vec<T>::covar(const vec<T>& y) const{
    //vec<T> x = *this;
    T mu =  (*this).meanvalue();
  T nu = y.meanvalue();
  if(length != y.length) error("The sizes of the covariance vectors have to be equal !. "); 
 if(length == 0 || length == 1 ) return 0; 
  T co =  (entry[0] - mu)*(y[0]-nu);
  for(int i=1; i < length; i++){
    co += (entry[i]-mu)*(y[i]-nu);
  }
  return (co/(T(length-1)));   //for x.length<=20 the nominator can be replaced 
}                            // by x.length (but then change if statement !)



//multiplies all entries
template<class T> 
T vec<T>::prodvalue() const{
 if(length == 0) return 0; 
 T pr = entry[0];
  for(int i = 1; i < length; i++){
    pr *= entry[i];
  }
  return pr;  
}


//adds all entries
template<class T> 
T vec<T>::sumup() const{
 if(length == 0) return 0; 
 T sum = entry[0];
  for(int i = 1; i < length; i++){
    sum += entry[i];
  }
  return sum;  
}


template<class T>  //if an entry is below a certain bound it is regarded as 0
void  vec<T>::set2zero() const {
  for(int i= 0; i < length; i++){
      if(abs(entry[i]) <= SMALL) entry[i] = 0.0;
  }
}



//bitonic merging, see [R.SEDGEWICK, Algo. in C++, Algo 8.2]. l = lower index, r = upper index. NOT STABLE 
template<class T>
void vec<T>::merge(int l, int m, int r) const{
  int n = length;
  int i,j;
  vec<T> aux(n);  //auxilary vec
  for(i = m+1; i > l; i--) aux[i-1] = entry[i-1];
  for(j = m; j < r; j++) aux[r+m-j] = entry[j+1];
  for(int k = l; k <= r; k++)
    if(aux[j] < aux[i])
      entry[k] = aux[j--];
    else entry[k] = aux[i++]; 
}

//STABLE merge version, i.e. data with same keys remain in same order:
template<class T>
void vec<T>::stablemerge(int l, int m, int r) const{
  int i,j,k;
  i = 0; j = l;
  int n = length;
  vec<T> aux((n+1)/2); //auxilary array
  while( j <= m) 
    aux[i++] = entry[j++]; //copy first part into aux
  
  i=0; k = l;
  while(k < j && j <= r) //copy second biggest element back
      //if(aux[i] <= entry[j]) entry[k++] = aux[i++];
      if(lessthan(aux[i], entry[j])) entry[k++] = aux[i++];
      else entry[k++] = entry[j++]; 
 
  while(k < j)  //copy back remaining part of aux
      entry[k++] = aux[i++];
}

/*
//stablemerge which stores index (permutation) 
template<class T>
void vec<T>::stablemerge_with_index(int l, int m, int r) const{
  int i,j,k;
  i = 0; j = l;
  int n = length;
  vec<int> index(n), idx((n+1)/2);
  for(int inc=0; inc<n; inc++)
      index[inc] = inc;

  vec<T> aux((n+1)/2); //auxilary array
  while( j <= m){ 
      aux[i] = entry[j]; //copy first part into aux
      idx[i] = index[j];
      i++, j++;
  }
  i=0; k = l;
  while(k < j && j <= r) //copy second biggest element back
      //if(aux[i] <= entry[j]) entry[k++] = aux[i++];
      if(lessthan(aux[i], entry[j])){
	  entry[k] = aux[i];
	  index[k] = idx[i];
	  k++, i++;
      }
      else{
	  entry[k] = entry[j]; 
	  index[k] = index[j];
	  k++, j++;
      }
  while(k < j){  //copy back remaining part of aux
      entry[k] = aux[i];
      index[k] = idx[i];
      k++, i++;
  }
  std::cout<<std::endl;
  std::cout<<"INDEX"<<std::endl;
  std::cout<<index<<std::endl;
  std::cout<<std::endl;
}
*/




//recursive implementation
template<class T>
void vec<T>::mergesort(int l,int r) const{
  if(r <= l) return;
  if(r > length-1) error("Bad move in >>vec<T>::mergesort<<. Indexing of l and r, i.e. m = r+l exceeds size of input vector !! Adjust this and call me again, hunk ;-)");
  int m = (r+l)/2;
  mergesort(l,m);
  mergesort(m+1,r);
  stablemerge(l,m,r);//merge(l,m,r);
}

//this sorts and merges the WHOLE array form l=0 to r=n-1 
template<class T>
void vec<T>::mergesort() const{
  mergesort(0,length-1);
}

//returns permutation after a vector has been sorted
template<class T>
vec<int> vec<T>::permutation(){
    int n = length, nmone = n-1;
    vec<int> index(n);
    for(int i = 0; i < n; i++)
	index[i] = i;
    
    
    //use BUBBLESORT to sort (disadvantage: O(n^2) !!):
    /*for(int i = nmone; i >= 0; i--)
	for(int j = 0; j<i; j=j+1)
	    if(entry[j+1]<entry[j]){
		permute(entry[j+1],entry[j]);
		permute(index[j+1],index[j]); //permute index accordingly
	    }
    */
    //use SELECTIONSORT (disadvantage: O(n^2), advantage: only n swaps compared to (n^2)/2 swaps in bubblesort !!):
    for(int i = 0; i < nmone; i=i+1){
	int mini = i;
	for(int j=i+1; j<n; j=j+1)
	    if(entry[j]<entry[mini]) mini=j;
	permute(entry[i],entry[mini]);
	permute(index[i],index[mini]);
    }

    return index;
}




template<class T> 
void vec<T>::shellsort() const{
  int n = length; 
  int gap, i, j;
  for(gap = n/2; gap > 0; gap/=2)
    for(i = gap; i <n; i++)
      for(j = i-gap; j >= 0; j -= gap)
	if(entry[j+gap] < entry[j]) std::swap(entry[j],entry[j+gap]); 
}

template<class T>
void vec<T>::randomvector(int eger){
  srand(time(0));
  for(int i =0; i< length; i++){
    entry[i] =  T(double(rand())/RAND_MAX);           
 }
}

// Returns the position(s) of an entry within an vector. If this entry is not part of the vector then an error message is launched. 
template<class T>
vec<int> vec<T>::position(T whereami){
    vec<int> pos;
    int count = 0;
    for(int i = 0; i < length; i++){
	if(whereami == entry[i]){
	    count++;
            pos.resize(count, i);  //enlarge vector by one entry with value i
        }
        
        if(i == length -1 && pos.size() == 0){
	    std::cout<< "Specified entry '"<<whereami<<"' could not be found while >>vec<T>::position<< was invoked. I stop all further computations to prevent severe errors !"<<std::endl;
            exit(1); 

        }
    }

    return pos;
}


//searches for different entries in an array, orders them and determines the dimension of each entry block
template<class T>
vec<vec<T> > vec<T>::find_different_entries(){
//vec<T> vec<T>::find_different_entries(){
    (*this).mergesort();   //fast sorting first  
    
 
    vec<T> simpleones;       //store entries without multiplicity
    vec< vec<T> > twovec(2);
    int lone = length-1;
    int count = 1, subdim = 1, newdim = 1;     
    vec<T> dims;  //store multiplicity of entry seperately

    for(int i = 0; i < lone; i++){
	if(entry[i] == entry[i+1]){
	    subdim++;
        }
        else{  
            //count++;
            dims.resize(count, subdim);
            simpleones.resize(count,entry[i]);
        
            subdim = 1; //reset
	    newdim++;
	    count++;
       }
        
 
  //compare last entries of vectores
        if(i == length-2){
            dims.resize(count, subdim);
            simpleones.resize(count, entry[lone]); 
         }
        //if(i == lone){// && entry[lone] != simpleones[simpleones.size()-1]){
	//   simpleones.resize(count++,entry[lone]);
	    //  dims.resize(count++, subdim);
            //break;
        }
      
	twovec[0] = simpleones;
	twovec[1] = dims;
	


	 
	return twovec;  //simpleones;//twovec;
	
 }


template<class T>
void vec<T>::write_2_file(const std::string& filename) const{
  std::ofstream outputTheStuff(filename.c_str(), std::ios_base::out);
    outputTheStuff << (*this) << '\n';
}


template<class T>
vec<T> vec<T>::reverse(){
    vec<T> aux(length);
    int lmone = length-1;
    for(int i = lmone; i >= 0; i--)
	aux[lmone-i] = entry[i];
    return aux;
}

//returns the first position where "seek" is encountered otherwise returns -1 (= not found at all). An O(n) search procedure. A faster one is 'search' that implements binary search. However, one should stick to member 'find' if the array is small, let's say less than 1000 entries (see comment for routine 'search') 
template<class T>
int vec<T>::find(T seek) const{
  
    for(int i = 0; i < length; i++){
	if(entry[i] == seek){
	    return i;
	    break;
	}
    }
    return -1; //not found
    
    //return pos;
}


//binary search is worthwile when seeking for entries in very HUGE arrays. In small ones it rather fails !!! (see [SEDGEWICK, Algorithms in C++, 3rd ed., pp 56]
template<class T>
int  vec<T>::search(T tosearch, int l, int r) const{
   int m;
  while(r >= l){
    m = (l+r)/2;
    if(tosearch == entry[m]){
      return m;
      //break;
    }
    if(tosearch < entry[m]) 
      r = m-1; 
    else l = m+1;
  }
  return -1; //searched for 'tosearch' did not yield any result 

}

//searching the whole array
template<class T>
int  vec<T>::search(T tosearch) const{
  int left = 0;
  int right = length-1;
  return ((*this).search(tosearch,left,right));

}


//O(nlog(n) + n) routine for determining if all entries are different or not
template<class T>
bool vec<T>::different() const{
  vec<T> refer = *this;  //if here is a reference, *this will also be sorted !
   refer.mergesort(); //sort entries
  int b = refer.size();
  for(int i = 0; i < b-1; i++){
    if(refer[i] == refer[i+1]){
      return false;
      break;
    }
  }
  return true;
}

template<class T>
vec<T> vec<T>::canonical_basis(int j){
    if(j < 0 || j >= length)
	error("In >>vec<T>::canonical_basis(int)<<: Argument exceeds vector range. \n");
    vec<T> basis(length); 
    basis[j] = T(1);
    return basis;
}

template<class T>
vec<T> vec<T>::subvec(int low, int up){
    if(low > up)
      std::swap(low,up);
    if((low < 0 || low >= length) || up < 0 || up >= length)
	error("In >>vec<T>::subvec(int, int)<<: Bound(s) out of range!\n");
    
    int subsz = up-low+1;
    vec<T> sub(subsz);
    for(int i = 0; i < subsz; i++)
	sub[i] = entry[low+i];

    return sub;
}



/*************************************************************************************************************************************************************************   CLASS FOR (REAL) GENERAL, UNSTRUCTURED & DENSE MATRICES    ********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************/

/*class AbstractMatrix{
 protected:
  int row;
 public:
  virtual ~AbstractMatrix(){}

  virtual rock() const = 0;
};
*/


template<class T> class matrix{
 private: int row; 
         int column; 
         T** entry;
                             
public: matrix(int , int, T**);             //constructor
	 matrix(int n = 0, int m = 0, T d = T());   //constructor:all entries = 0 
	 template<class VecContainer>
	 matrix(int,const VecContainer&); //construct from vector (either your own vector class or the STL vector, fore instance), given the number of of rows
 
 
         matrix(const matrix&);
         ~matrix();
  
	 
	 //different methods to get the number of rows/columns
	 int rock() const {return row;}         
	 int roll() const {return column;} 
	 
	 int rows() const {return row;}
	 int cols() const {return column;}
	 
	 int size1() const {return row;}
	 int size2() const {return column;}

  //clear matrix, i.e. empty it and assign 0 dimensions
  inline void clear(){
    for(int i = 0; i < row; ++i) delete[] entry[i];
    delete[] entry;
    row = 0; column = 0;
  }

  /* void init(int n, int m, const T& d = T()){ */
  /*   assert((row == 0) && (col == 0)); */
    
  /*    try{ */
  /*     row = n; */
  /*     column = m; */
  /* 	entry = new T* [row];  */
  /* 	for (int i =  0; i < row; i++) { */
  /* 	    entry[i] = new T [column];           */
  /* 	    for (int j = 0; j < column; j++) entry[i][j] = d; */
  /* 	} */
  /*   } */
  /*   catch(std::bad_alloc){ */
  /* 	error("Damn it Marc, I haven't enough memory to construct a matrix object via >>matrix<T>::matrix(int, int, T**)<<. "); */
  /*   } */
  /* } */


  //! v stores entries row-wise
  template<class V>
    void read_from_vector(const V& v){
    assert(row*column == (int)(v.size()));
    for(int i = 0; i < row; ++i)
      for(int j = 0; j < column; ++j)
	entry[i][j] = v[i*column + j];
		  
  }

  vec<int> dimensions();

  matrix& operator=(const matrix&);     // overload = (copy) (remember: &
  matrix& operator+=(const matrix&);    // overload +=, ...
  matrix& operator-=(const matrix&);   //(matrix division not defined)
  matrix& operator*=(const matrix&); 
  matrix& operator*=(const T&);  
  matrix& operator/=(const T&); 

  bool operator==(const matrix&);
  bool operator!=(const matrix&);

  vec<T> operator*(const vec<T>&) const;  // matrix-vector multiplication
  T frobenius() const;                              //Frobenius norm
  T onenorm() const;                               //2-Norm   
  T maxnorm() const;
  T off() const;                    //"norm" of the off-diagonal elements 

  matrix<T> transpose() const;                    //transpose of matrix
  matrix<T> identity();
  template<class S> friend matrix<S> ident(int); //alternative identity 

  void rearrangedata(int, int);
  void print_matrix();

  void print_in_matlab_style(const std::string name, int prec = 15) const{
    std::cout << name << " = [" << std::setprecision(prec);
    std::string delim = ", ";
    for(int i = 0; i < row; ++i){
      for(int j = 0; j < column; ++j){
	delim = ", ";
	if(j == column-1){
	  delim = ";";
	  if(i == row-1)
	    delim = "]";
	}
	
	std::cout << entry[i][j] <<delim;
      }
      std::cout << std::endl;
    }
  }

  // void inmatrix();
  T* operator[](int i) const { return entry[i]; }  // subscripting, row i     
    matrix<T> operator+(const matrix<T>&); 
    matrix<T> operator-(const matrix<T>&);
  
 template<class S> friend matrix<S> operator+(const matrix<S>&, const matrix<S>&); //matrix addition 
  template<class S> friend matrix<S> operator*(const matrix<S>&, const matrix<S>&); //matrix multiplication
  template<class S> friend vec<S> operator*(const vec<S>&, const matrix<S>&); //vector-matrix multiplication

template<class S> friend matrix<S> operator*(S, const matrix<S>&);
template<class S> friend matrix<S> operator*(const matrix<S>&, S);
//scalar-matrix, matrix-scalar multiplication, respectively and matrix-scalar division
template<class S> friend matrix<S> operator/(const matrix<S>&, S);

template<class S> friend matrix<S> operator-(const matrix<S>&);
template<class S> friend matrix<S> dyad(const vec<S>&, const vec<S>&); 
std::ostream& print(std::ostream&) const;  //output-aternative
 std::istream& read(std::istream&) const;
  
  
  /******** LU-decomposition, case: row = column ******************/
template<class S> friend matrix<S> ludecomp(matrix<S>&);
template<class S> friend matrix<S> lu(matrix<S>&); //lu with part. pivoting
template<class S> friend matrix<S> lupp(matrix<S>&); //lu with part. pivoting
template<class S> friend matrix<S> lucp(matrix<S>&);
  matrix<T> lufac(vec<int>&,  T& ) const;
  vec<T> gausselimination(vec<T>&) const;
  vec<T> forbackwardsub(vec<int>&, vec<T>&) const;

 /******** LU-decomposition, case: row > column ******************/
template<class S> friend matrix<S> ludemgn(matrix<S>&);  


 /******** LU-decomposition, case: column > row ******************/
template<class S> friend matrix<S> ludengm(matrix<S>&);  


template<class S> friend S det(matrix<S>&);
  bool sdd() const;  //strict digonal dominance
  bool posdef() const; // positive definiteness
  vec<T> getcolumn(int) const; //extract a specific column...
  vec<T> getrow(int) const;   //...row of a given matrix
  
vec<int> partition(int, int);
vec<int> partmult(int, int);     //partition for matrixmultiplication

template<class S> friend  matrix<S>  blockadd(const matrix<S>&, const matrix<S>& );
 template<class S> friend  matrix<S>  blockmult(const matrix<S>&, const matrix<S>& );
   
 matrix<T> blockplus(const matrix<T>&);
  //matrix<T> blockmult(const matrix<T>&)



matrix<T> francisqrstep() const;
matrix<T> wilkinsonqrstep() const; //for symmetric matrices
void set2zero() const;
matrix<T> Sigma() const;

  //matrix<T> eigenqr() const;
   vec<T> eigenqr() const;
  //T eigenqr() const;

  //vec<T> symmetricqr() const;
  matrix<T> symmetricqr() const;
  //int symmetricqr() const;
  //matrix<T> implicitql() const;   //see Wilkinson, Reinsch
  vec<T> implicitql() const;
matrix<T> doubleimplicitshiftfrancis() const;


  //calculation of the sensitivity matrix in case that the eigenvalues have widely varying magnitudes. Balancing does not harm ! Mind that symmetric matrices need not to be balanced at all !!
  matrix<T> balancing() const;

 //GAUSSian elimination with partial pivoting (for n times n matrices) to solve the system Ax = b
  vec<T> gausselimpartpivot(vec<T>& ) const;


template<class V>
  V& solve(V& b) const{
  if(row != column) error("I only consider the n times n case here in >>matrix<T>::solve<<. ");
 
  if(row != (int)b.size()) error("Vector size and coefficient matrix dimension do not match in >>matrix<T>::solve<<. ");
  
  matrix<T> A = *this;  //reference
  int n = A.row; int m = A.column;

  //check for strict diagonal dominance
   matrix<T> B = A.transpose();
   bool decide = B.sdd();
  
  if(decide == true){//strict diagonal dominance
   for(int k = 0; k< n-1; k++){
    if(A[k][k] == 0) error("Pivot equals 0 in >>matrix<T>::solve(V&)<< . "); 
    for(int i = k+1; i<n; i++){       
      if(A[i][k]!=0){  
	 A[i][k] /= A[k][k];
	 for(int j = k+1; j<m; j++){   
         A[i][j] -= A[i][k]*A[k][j];
         }
      }
    }
   }   
   //forward substitution for Ly =b. (y stored in b)
   for(int i = 1; i < n; i++)
     for(int j = 0; j < i; j++) 
       b[i] -= A[i][j]*b[j];

   //back substitution for Ux = y. (x stored in b)
   for(int i = n-1; i >= 0; i--){ 
     for(int j = i+1; j < n; j++)b[i] -= A[i][j]*b[j];
     b[i] /= A[i][i];  
   }   
   return b;
   }//end if:   
  

  else if(decide == false){

   vec<int> piv(n);      //pivot(row) vector (contains row numbers)
  for(int k = 0; k < n; k++){   //..and stores swapping information
    piv[k] = k;                 // piv=[0,1,2,......,n-1]   
  }                                 
  
  for(int k = 0; k < n-1; k++){      //here the loop starts
    int seek = k;                    //find entry in column k, row pvt[k]
    T mag = abs(A[piv[k]][k]);       // with largest magnitude mag & take it as
    
    for(int i = k+1; i < n; i++){   // pivot element = max.|column entry|
      if(abs(A[piv[i]][k]) > mag){  //found such element. 
	mag = abs(A[piv[i]][k]);    //pivot search
      
      seek = i;  
      }
        
    }
    if(!mag) error("Pivot is ZERO in LU factorization.\n    Matrix is singular up to machine precision.");
    if(seek != k) {
      //ONLY HERE: the rows of A need NOT to be interchanged in order to reduce the operations. The permutation information is solely stored  in vector "piv" to which is later referred in the backward and forward substititution !
      permute(piv[k], piv[seek]);
    }    
     for(int i = k+1; i < n; i++){       //entry elimination below A[pvt[k]][k]
      if(A[piv[i]][k] != 0){
        A[piv[i]][k] /=  A[piv[k]][k]; //"zero out" entries in the column
        for(int j = k+1; j < m; j++)
	  A[piv[i]][j] -= A[piv[i]][k]*A[piv[k]][j];
      }
    }
   }
  
  //check if determinant is zero (here just verify if product of diagonal entries is zero regardless the sign). If yes there exist no solution for the matrix is defective.
  T dia = 1.0;
  for(int i = 0; i < n; i++) dia *= A[i][i];
  if(std::abs(dia) < 1.e-14) error("Matrix is singular -- at least within my precision. Hence there is no solution when applying  >>matrix<T>::solve<< !"); 

  //forward substitution Ly = Pb
  for(int i = 1; i < n; i++){
    for(int j = 0; j < i; j++){
      b[piv[i]] -= A[piv[i]][j]*b[piv[j]];
    }
  }
  //back substitution Ux = y
  V x(n); //stores solution in correct order 
  for(int i = n-1; i >= 0; i--){
    for(int j = i+1; j < m; j++) b[piv[i]] -= A[piv[i]][j]*x[j];
    //if the pivot element is zero then we assign a tiny number to the pivot element which is desirerable for some applications. 
    if(A[piv[i]][i] == 0.0){
       std::cout<<"WARNING: I have to cope with a singular matrix in >>matrix<T>::solve(V&)<<. Hence I assigned a tiny number as pivot element. "<<std::endl;
      A[piv[i]][i] = TINY;;


      // error("Matrix is SINGULAR in >>matrix<T>::gausselimpartpivot(vec<T>&)<<. ");
    }
    x[i] = b[piv[i]]/A[piv[i]][i];
  }
  b = x; 

 
  } //end else
 return b;
 }




  //rank of a matrix: Note that rank(A) = rank(A^T).
  template<class S> friend int rank(matrix<S>&);
  int rank() const;
  bool symm() const;
  template<class S> friend S tr(matrix<S>&);  //trace of a n*n matrix
 
  //cut off from the column specified as argument (inclusive that column)
  matrix<T> cut_off_from_column(int) const;
  matrix<T> cut_off_from_row(int) const;
  matrix<T> delcolumn(int) const;
  matrix<T> delrow(int) const;
 
  //skip from specific row/column
  matrix<T> skip_from_row(int) const;
  matrix<T> skip_from_column(int) const;

  //turns a given matrix into .txt file
  void mat2txt(const std::string&) const;

  //matrix<T> wolpertinger(int) const;
  std::list<matrix<T> > wolpertinger(int);

  vec<T> order_rows_by_column_entries(int);

  std::list<matrix<T> > order_stimuli(int); //store differently sized objects
  void write_2_file(const std::string&, const std::string& comment = std::string()) const;
  
  matrix<T> reverse_by_rows();
 
  template<class S> friend int size_of_pointer(matrix<S>*);

  //extract columns/rows
  matrix<T> choice_by_columns(const vec<int>&); 
  matrix<T> choice_by_rows(const vec<int>&); 

  void refill_column(vec<T>, int);
  void refill_column(T, int);
  void refill_row(vec<T>, int);
  void refill_row(T, int);

  matrix<T> reduce_to_3_rows(int, int, int, int, int, int);

  matrix<T> extract_via_same_values(int, T);
};


//Member initialization
 
template<class T> 
matrix<T>::matrix(int n,  int m, T** e){ // construct from T**
    try{
      row = n;
      column = m;
	entry = new T* [row]; 
	for (int i =  0; i < row; i++) {
	    entry[i] = new T [column];          
	    for (int j = 0; j < column; j++) entry[i][j] = e[i][j];
	}
    }
    catch(std::bad_alloc){
	error("Damn it Marc, I haven't enough memory to construct a matrix object via >>matrix<T>::matrix(int, int, T**)<<. ");
    }
}


template<class T> 
matrix<T>::matrix(int n, int m, T a){              // construct from a type T
    try{
	entry = new T* [row=n];         //remember: T** d; create new space: d=new T* [n] 
	column = m;
	for (int i =  0; i< row; i++) {
	    entry[i] = new T [column];
	    for (int j = 0; j < column; j++) entry[i][j] = a;
	}
    }
    catch(std::bad_alloc){
	error("I am really sorry, Marc but I haven't got enough memory to construct a matrix object via >>matrix<T>::matrix(int, int, T)<<. ");
    }
}


//a third constructor which constructs from a vector of size row*column containing the entries ROW by ROW (IPOPT C-style) of a general unstructured dense matrix. This means that special structures such as symmetry, banded matrices etc.are neither checked nor explicitly supported.  
template<class T> 
template<class VecContainer>
matrix<T>::matrix(int numofrows, const VecContainer& v):row(numofrows){
  int vsi = int(v.size());
  
  if(vsi%numofrows != 0 || numofrows > vsi)
    error("In constructor >>matrix::matrix(int, const vec<T>&)<<:\n     Dimension of first argument is illegal. Change it, hunk\n");

  //assign to fields of class matrix
  column = vsi/numofrows;
  entry = new T* [row];
  
  int k = 0;       //reaches from 0 to (row*colum)-1
  for(int i = 0; i < row; ++i){
    entry[i] = new T [column];
    for(int  j = 0; j < column; ++j){
     entry[i][j] = v[k++];
    }
  }
}




template<class T>                                 //destructor
inline  matrix<T>::~matrix(){
          for(int i=0; i< row; i++) delete[] entry[i]; 
         delete[] entry;
       }
     

template<class T> 
matrix<T>::matrix(const matrix& M) {                      // copy constructor
    try{
	entry = new T* [row = M.row]; 
	column = M.column;
	for (int i =  0; i< row; i++) {
	    entry[i] = new T [column];
	    for (int j = 0; j < column; j++) entry[i][j] = M[i][j]; 
	}
    }
    catch(std::bad_alloc){
	error("I am really sorry, Marc but I haven't got enough memory to construct a matrix object via copy-constructor >>matrix<T>::matrix(const matrix&)<<. ");
    }
}

template<class T> 
matrix<T>& matrix<T>::operator=(const matrix& M) {           // copy assignment(for making objects equal)
  if (this != &M) { 
   // if (row != M.row || column != M.column) error("Matrix sizes do not match, Marc !");
    //for (int i = 0; i < row; i++) 
    // for (int j = 0; j < column; j++) entry[i][j]  = M[i][j];


    //this is just information...
    //if((row > 0 && column > 0) && (row != M.row || column != M.column)) std::cout<<"CAUTION: Assignment (>>matrix<T>::operator= <<) of  NON-DIMENSIONLESS matrices of different sizes may lead to unexpected - even ERRONEOUS - results ! "<<std::endl;
/******************** Make dimensions equal  **************/
     
    //Delete the matrix if dimensions do not match
    try{
	if(row != M.row || column != M.column){
	    for(int i = 0; i < row; i++) delete[] entry[i];
	    delete[] entry;
    
	    //Create it new with values of M
	    entry = new T* [M.row];
	    for(int i = 0; i < M.row; i++) 
		entry[i] = new T [M.column];
	}
    }
    catch(std::bad_alloc){
	error("I am really sorry, Marc but I haven't got enough memory to construct a matrix object via assignment >>matrix<T>::operator=(const matrix&)<<. ");
    }

/********************************/

   //After adjustment: assign values of M to adjusted matrix
     for (int i = 0; i < M.row; i++) 
       for (int j = 0; j < M.column; j++) entry[i][j]  = M.entry[i][j];
  
     row = M.row;
     column = M.column; 
  }

  return *this;
}

template<class T> 
matrix<T>& matrix<T>::operator+=(const matrix& M) {         // add-assign
  if (row != M.row || column != M.column) error("Matrix sizes do not match, Marc !  ");
  for (int i = 0; i < row; i++)
    for (int j = 0; j < column; j++) entry[i][j] += M[i][j];
  return *this;
}


template<class T> 
matrix<T>& matrix<T>::operator-=(const matrix& M) {         // diff-assign
  if (row != M.row || column != M.column) error("Matrix sizes do not match, Marc !  ");
  for (int i = 0; i < row; i++)
    for (int j = 0; j < column; j++) entry[i][j] -= M[i][j];
  return *this;
}

template<class T> 
matrix<T>& matrix<T>::operator*=(const matrix& M) {         // mult-assign
  if (column != M.row) error("Matrix sizes do not match, Marc !  ");
  //Note: "*=" is ascribed to "*"
  *this = *this * M;  
   return *this;
}

template<class T> 
matrix<T>& matrix<T>::operator*=(const T& a) {       
//Note: "*=" is ascribed to "*"
  *this = *this * a;  
   return *this;
}

template<class T> 
matrix<T>& matrix<T>::operator/=(const T& a) {       
  if(a == 0.0) error(" Hey lad, division by 0 w.r.t. >>matrix<T>::operator/=(const T&)<< is not really cool !");
  for(int i = 0; i < row; i++) 
    for(int j = 0; j< column; j++)
      entry[i][j] /= a;
   
  return *this;
}


//DEF.: A matrix M is equal to another N if they have equal row/column numbers and equal entries. Some overloadings of boolean operators: Test if matrices are equal w.r.t entries (while checking rows numbers, column numbers ). 
//CAUTION: When outputting an expression of the form 'M1 == M2', never forget to use parentheses, since '<<' has higher priority than '==', i.e.: std::cout<< (M1 == M2) <<std::endl;   (NEVER: std::cout<< M1 == M2 <<std::endl; !!!)
template<class T> 
bool matrix<T>::operator==(const matrix& logic) {
  //the next line is optional ! 
 if(this == &logic){return true;}
  //the following is needed ! otherwise it outputs false at the end !! 
 if(row != logic.row || column != logic.column) return false;
  for(int i =0; i < row; i++){
    for(int j = 0; j < column; j++){
     if(entry[i][j] != logic[i][j]) return false;
    }  
  } 
 return true;
}


template<class T> 
inline bool matrix<T>::operator!=(const matrix& logic) {
  //since the operator  is a bool we can just negate it via "!"  
 return ! operator==(logic);
}

template<class T> 
matrix<T> matrix<T>::operator+(const matrix& M) {  // usage: m = m1 + m2
 matrix<T> sum = *this;           // user-defined copy constructorits friends
 sum += M;                                    // is important here
 return sum;                                    // otherwise m1 would be changed
}

template<class T>  //No friend
matrix<T> matrix<T>::operator-(const matrix& M) {
  if(row != M.row || column != M.column) error("matrix sizes do not match. Hence I cannot compute any matix difference");
  matrix<T> diff = *this;
  diff -= M;
  return diff;
} 


template<class T> //a friend again  (A = -A)
inline matrix<T> operator-(const matrix<T>& A) {
  for(int i = 0; i< A.row; i++){
    for(int j = 0; j< A.column; j++){
      A[i][j]= -A[i][j];
    }
  }
  return A;
} 


template<class T>  // a friend: scalar-matrix  multiplication
matrix<T> operator*(T scal, const matrix<T>& M){
  matrix<T> sm(M.row, M.column);
  for(int i = 0; i <M.row; i++){ 
    for(int j=0; j< M.column; j++){
        sm[i][j] = scal*M[i][j];
    }  
  }
  return sm;
} 


template<class T>  // a friend: matrix-scalar  multiplication
matrix<T> operator*( const matrix<T>& M, T scal){
  matrix<T> sm(M.row, M.column);
  for(int i = 0; i <M.row; i++){ 
    for(int j=0; j< M.column; j++){
	sm[i][j] = M[i][j]*scal;     
    }  
  }
  return sm;
} 

template<class T>  // a friend: matrix-scalar division
matrix<T> operator/( const matrix<T>& M, T scal){
  if( scal == 0.0) error("Divions by 0.0 in >>matrix<T> operator/(const matrix<T>&, T)<<. ");
  matrix<T> sm(M.row, M.column);
  for(int i = 0; i <M.row; i++){ 
    for(int j=0; j< M.column; j++){
      sm[i][j] = M[i][j]/scal;     
    }  
  }
  return sm;
} 





template<class T>   //friend member:   s = v^T*M
vec<T> operator*(const vec<T>& w, const matrix<T>& M){
  if(w.size() != M.row) error("Vector sizes do not match when using 'vec*matrix'");
  vec<T> s(M.column);
  for(int i = 0; i<M.row; i++){
    for(int j = 0; j<M.column; j++){
      s[j]+= w[i]*M[i][j];
    }
  }
  return s;
}

template<class T>
vec<T> matrix<T>::operator*(const vec<T>& v) const {  //NO friend: u =M*v
 if(column != v.size()) error("Bad vector matrix dimensions when using >>matrix<T>::operator*(const vec<T>&)<<. Having slept during linear algebra courses ?!"); 
  vec<T> u(row);
  for(int i=0; i < row; i++){
    for(int j=0; j<column; j++){
       u[i]+= entry[i][j]*v[j];
    }
  }
   return u;
  } 


template<class T>         //dyadic product
matrix<T> dyad(const vec<T>& v, const vec<T>& u){
  matrix<T> D(v.size(), u.size());
  for(int i = 0; i < v.size(); i++){
    for(int j = 0; j < u.size(); j++){
      D[i][j] += v[i]*u[j];   //update A = A + vu^T
    }
  }
 return D;
}


template<class T> //not a friend: remember the scope resolution(::)
matrix<T> matrix<T>::transpose() const{ 
  matrix<T> M(column, row);      //transpose matrix

   for(int i=0; i< row; i++){
    for(int j=0; j< column; j++) M[j][i]=entry[i][j]; // M[i][j] =entry[j][i];
    }
  return M;
}

template<class T>
matrix<T> matrix<T>::identity(){
if(row != column) error("Identity matrix only defined for square matrices. "); 
  int n= row; int m = column; 
 for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      entry[i][j] = 0.;
      entry[i][i] = 1.;
    }
  }
  return *this;
 }


template<class T>         //alternative definition of identity matrix 
matrix<T> ident(int n){    //not for output !!
  matrix<T> Id(n, n);
  for(int i = 0; i< n; i++) {
    for(int j = 0; j< n; i++) {
      Id[i][j] = 0;
      Id[i][i] = 1;
    }
  }   
  return Id;
} 

//Suppose a real data matrix is given whose columns represent experimental observations, then the following routine can rearrange the columns by swapping till the desired structure is achieved.
template<class T>
void  matrix<T>::rearrangedata(int cola, int colb){
  if(cola >= column || colb >= column) error("Column index exceeded in >>matrix<T>::rearrangedata()<<. ");
  for(int i = 0; i < row; i++)
    std::swap(entry[i][cola],entry[i][colb]);   
}


template<class T>  //matrix multiplication -- a friend: no scope resolution
matrix<T> operator*(const matrix<T>& A, const matrix<T>& B){
 if(A.column != B.row) error("Matrix dimensions do not match, Marc ! ");
 matrix<T> C(A.row, B.column);
  for(int i=0; i<A.row; i++){
    for(int j=0; j<B.column; j++){ 
      //T sum = 0; 
     for(int k=0; k<A.column; k++){
	C[i][j] +=A[i][k]*B[k][j];
       
     }
     //C[i][j]=sum;
    }
  }  
  return C;
}

/* Matrix Norms; for formulas see [Golub/VanLoan, section 2.3.2.] */
template<class T>
T matrix<T>::frobenius() const{
  T frobi =0.0; 
  for(int i=0; i< row; i++) {
    for(int j=0; j< column; j++) {
      frobi += abs(entry[i][j])*abs(entry[i][j]);
    }
  }
return sqrt(frobi);
}

template<class T>
T matrix<T>::onenorm() const{
  T oger = 0.0;
  for(int j = 0; j< column; j++){
    T tmp = 0.0;
      for(int i =0; i< row; i++) tmp += abs(entry[i][j]);
      oger = std::max(oger, tmp);
  }
  return oger;
}

template<class T>
T matrix<T>::maxnorm() const{
  T infi = 0.0;
  for(int i = 0; i< row; i++){
    T tmp = 0.0;
    for(int j =0; j< column; j++) tmp += abs(entry[i][j]);
    infi = std::max(infi, tmp);
    
  }
  return infi;
}


template<class T>
T matrix<T>::off() const{
  T offdiag = 0.;
 
  //if the matrix is symmetric, then one needs to regard the upper 3ang. part
  matrix<T> Dia = *this;
  if(row == column && Dia.symm()){
   for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      if( j >= i+1){
	if(j!= i) offdiag += entry[i][j]*entry[i][j];
      }
    }
   }
   offdiag *= 2.;
  }
  else{
  for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      if(j != i) offdiag += entry[i][j]*entry[i][j];

	}
  }
  }  
 return offdiag;
}


/*Matrix output routines*/
template<class T> 
 void matrix<T>::print_matrix() {
 
  for(int i=0; i<row; i++){
    for(int j=0; j<column; j++) std::cout.precision(7) << entry[i][j] << '\t'; 
    std::cout << '\n';
  }
}

template<class T> 
std::ostream& matrix<T>::print(std::ostream &ou) const {
  //guarantees that very small values are set to ZERO
  /*for(int i=0; i<row; i++){
   for(int j=0; j<column; j++){
     if(abs(entry[i][j]) <= 1e-13) entry[i][j] = 0;

   }
   }*/

/* "std::setw(n)" = column for next output object is n spaces broad. "left" ("right") calibrate column objects left (right). "right" = default value. */
  int IL = int(BROAD);           //since 'BROAD' is declared const in 'param.h'
  if(int(sizeof("ou")) >= BROAD){
    IL += int(sizeof("ou"));    
  }
  for(int i=0; i<row; i++){
    for(int j=0; j<column; j++) 
    
      ou <<"   "<<std::left<<std::setw(IL)<<std::setprecision(VAL)<< entry[i][j] << '\t';
     ou << '\n';//<< std::std::endl;
    
 }
 return ou;
} 

template<class T>  //defined OUTSIDE the class
std::ostream& operator<<(std::ostream &ou, const matrix<T>& M){
  return M.print(ou);
}


template<class T>
std::istream& matrix<T>::read(std::istream &inside) const {

  for(int i=0; i<row; i++){
    for(int j=0; j<column; j++) 
      inside>>entry[i][j];
  }
  /*char c = 0;
  
  int rabbies = 0, cholera = 0;
  

  inside>>c;
  while(!inside.eof()){
    if(c == '#'){
      inside.putback(c);
    }
    while(c != '\n'){
       inside >> entry[rabbies][cholera]; 
        if(c == '\t'){ 
	  inside.putback(c);
         cholera++;
        }         
    }  
    rabbies++;
    } */
  return inside;
}



template<class T>  //ATTENTION: non-constant reference
std::istream& operator>>(std::istream& inside, matrix<T>& M){ 
 return M.read(inside);
}



/*Strict diagonal dominance of a square matrix. The following routine calculates the sum of all non-diagonal entries for each i = 1,...,n and returns a vector, containing the sums w.r.t. indices i.*/

/*template<class T>
vec<T> strictdiagdom(matrix<T>& A){
  if(A.row != A.column) error("A has to be SQUARE");
  vec<T> v(A.row);
  int i,j; 
   for(i=0; i<A.row; i++){
   
    for(j =0; j < A.column; j++){
      if(j!=i){ 
      v[i] += abs(A[i][j]);
      }
     }
   }
   return v;
   }*/


//Check if a given n times n matrix is strictly diagonally dominant 
template<class T>
bool matrix<T>::sdd() const{
  if(row != column) error("Gi'me a square matrix in >>matrix<T>::sdd()<<. ");
  matrix<T> A = *this;
 
 vec<T> v(A.row);
  int i,j; 
   for(i=0; i < A.row; i++){
   
    for(j =0; j < A.column; j++){
      if(j!=i){ 
      v[i] += abs(A[i][j]);
      }
     }
   }
  
  bool decide = true; //by default
  for(int i = 0; i < A.row; i++){
    if(abs(A[i][i]) <= v[i]){//if |a_ii| is not > sum(|a_ij|) then A         
      decide = false;        //can never be s. d. dominant !!
       break;
    }
  } 
  return decide;
}

//ATTENTION: The following routine only obeys the definition of positive definiteness and is correct only in the case of SYMMETRIC matrices. It does not detect an check if a given square matrix A is positive definite, i.e. it holds x^TAx > 0 for ANY vector x, that is not the zero vector   
template<class T>
bool matrix<T>::posdef() const{
  matrix<T> A = *this;
  if(A.row != A.column) error("Expect equal row and column sizes in >>matrix<T>::posdef()<< since positive definiteness is only defined for square matrices. ");
  //generate a random vector
  vec<T> x(A.row);
  vec<T> null(A.row);
  srand(time(0));   //start pseudorandom generator with system ticks
 for(int i = 0; i< A.row; i++){ //fill v with random integers 
   x[i] = double(rand())/RAND_MAX;
 }
 if(x == null){ //if it happens that v's entries are all zero... 
   for(int i = 0; i < A.row; i++) //..then take the unity vector as default vec
     x[i] = 1.0;
 }

 if(x*A*x > 0.0) return true;   //or inner(x, A*x) > 0.0
 else return false;
}

//extract the jth column of a given matrix
template<class T>
vec<T> matrix<T>::getcolumn(int j) const{
  if(j >= column || j < 0) error("Your index exceeds column dimension while using >>matrix<T>::getcolumn(int)<< !"); 
 vec<T> v(row);
  for(int i = 0; i < row; i++)
    v[i] = entry[i][j];
  return v;
}

//extract the ith row of a given matrix
template<class T>
vec<T> matrix<T>::getrow(int i) const{
  if(i >= row || i < 0) error("Your index exceeds row dimension while using >>matrix<T>::getrow(int)<< !"); 
 vec<T> v(column);
  for(int j = 0; j < column; j++)
    v[j] = entry[i][j];
  return v;
}






template<class T> 
matrix<T> ludecomp(matrix<T>& A){    //Golub-Van Loan (LU without pivoting)
 int n,m;                             
 n = A.row; m = A.column; 
   if(n!=m) error("Matrix has to be SQUARE.");
  
   for(int k = 0; k< n-1; k++){
    if(A[k][k] == 0) error("Pivot equals 0"); //mind the "==" here !!!
    for(int i = k+1; i<n; i++){       
      if(A[i][k]!=0){  
	 A[i][k] /= A[k][k];
	 for(int j = k+1; j<m; j++){   
        A[i][j] -= A[i][k]*A[k][j];
      }
    }
  }
 }  
 return A;   
}

//case row m > column n
template<class T> 
matrix<T> ludemgn(matrix<T>& A){  
  int m = A.row;
  int n = A.column;

  if(m<=n) error("I want to see a m times n, m > n, matrix in >>ludemgn<<. ");
 
  for(int k = 0; k < n; k++){
    if(A[k][k] == 0.0) error("Pivot is 0 in >>ludemgn<<. ");
    for(int i = k+1; i < m; i++){
      if(A[i][k] != 0.0){
         A[i][k] /= A[k][k];
         if(k < n){
           for(int j = k+1; j < n; j++){
	     if(A[k][j] != 0.0){ 
               A[i][j] -= A[i][k]*A[k][j];
             }
           }
         }
      } 
    } 
  }
  return A;
}


//case row m < column n
template<class T> 
matrix<T> ludengm(matrix<T>& A){  
  int m = A.row;
  int n = A.column;

  if(n<=m) error("I want to see a m times n, m > n, matrix in >>ludemgn<<. ");
 
  for(int k = 0; k < m; k++){
    if(A[k][k] == 0.0) error("Pivot is 0 in >>ludemgn<<. ");
    for(int j = k+1; j < n; j++){
      if(A[k][j] != 0.0){
         A[k][j] /= A[k][k];
         if(k < m){
           for(int i = k+1; i < m; i++){
	     if(A[i][k] != 0.0){ 
               A[i][j] -= A[i][k]*A[k][j];
             }
           }
         }
      } 
    } 
  }
  return A;
}



template<class T>            //LU decomposition with (partial) pivoting
matrix<T> lupp(matrix<T>& A){   
  int n,m; 
  n= A.row; m= A.column;
  if(n != m) error("Matrix is not SQUARE.");
   
  vec<int> piv(n);      //pivot(row) vector (contains row numbers)
  for(int k = 0; k < n; k++){   //..and stores swapping information
    piv[k] = k;                 // piv=[0,1,2,......,n-1]   
  }                                 
  
  for(int k = 0; k < n-1; k++){      //here the loop starts
    int seek = k;                    //find entry in column k, row pvt[k]
    T mag = abs(A[piv[k]][k]);       // with largest magnitude mag & take it as
    
    for(int i = k+1; i < n; i++){   // pivot element = max.|column entry|
      if(abs(A[piv[i]][k]) > mag){  //found such element. 
	mag = abs(A[piv[i]][k]);    //pivot search
      
      seek = i;  
      }
        
    }
    if(!mag) error("pivot is ZERO in lu(matrix<T>&).");
    if(seek != k) {
      for(int l = 0;l<m;l++) 
        permute(A[piv[seek]][l], A[piv[k]][l]); //Swap row contents
    }    
     for(int i = k+1; i < n; i++){       //entry elimination below A[pvt[k]][k]
      if(A[piv[i]][k] != 0){
        A[piv[i]][k] /=  A[piv[k]][k]; //"zero out" entries in the column
        for(int j = k+1; j < m; j++)
	  A[piv[i]][j] -= A[piv[i]][k]*A[piv[k]][j];
      }
    }
   }
  return A;
}


/*LU decomposition with partial pivoting. HOWEVER, partial pivoting can be avoided as soon as the transposed (square) matrix is strictly diagonally dominant [see Golub/VanLoan, Thm. 3.4.3].    */
template<class T>      
matrix<T> lu(matrix<T>& A){   
  int n,m;                   
  n= A.row; m= A.column;
  if(n != m) error("Matrix is not SQUARE.");
 
  vec<int> pivot(n);
  T d;

  matrix<T> B = A.transpose();
  bool crucial = B.sdd();  //Is transpose of A stricly diagonally dominant ?
  if(crucial == true)  return ludecomp(A); //no pivoting is needed !!
  else return  A.lufac(pivot, d);//lupp(A);
}

/*********************************************************************************************************************************************************************************************************************************************/

template<class T> 
matrix<T> matrix<T>::lufac(vec<int>& piv, T& d) const{
  if (row != column) error("No square matrix in >>matrix<T>::lufac<<. ");
  matrix<T> A = *this;
  
  d=1.;              //for determinant calculation("sign track keeping")   
  int n,m; 
  n= A.row; m= A.column;
  if(n != m) error("Matrix must always be SQUARE for determinant calculation.");
 
 
  for(int k = 0; k < n; k++){   
    piv[k] = k;                  
  }                                 
  
  for(int k = 0; k < n-1; k++){   //main loop  
    int seek = k;                   
    T mag = abs(A[k][k]);       
    
    for(int i = k+1; i < n; i++){   
      if(abs(A[i][k]) > mag){  
	mag = abs(A[i][k]);    
      
      seek = i;  
      }
        
    }
    //if(!mag) error("pivot is ZERO in lu(matrix<T>&).");
    if(seek != k) {
      for(int l = 0; l < m;l++) {
        permute(A[seek][l], A[k][l]); //permute rows 
      } 
      permute(piv[seek], piv[k]);   //update pivot vector (needed later on)
      d = -d;                  //for determinant(remind permutation number).
    }                         //even: sign = +; odd: sign = -
     for(int i = k+1; i < n; i++){       
      if(A[i][k] != 0){
        A[i][k] /=  A[k][k]; 
        for(int j = k+1; j < m; j++)
	  A[i][j] -= A[i][k]*A[k][j];
      }
    }
  }
  A.set2zero();   //if some entries are too small then set 'em to 0 
 return A;
}

//this solves the square system Ax = b after LU decomposition was invoked on A.
template<class T> 
vec<T> matrix<T>::forbackwardsub(vec<int>& piv, vec<T>& b) const{
  matrix<T> A = *this;
  int n = row;
  int m = column;
  if( n != b.size() ) error("Vec-matrix dimensions do not match in >>matrix<T>::forbackwardsub<<. ");

//solve Ax = b. First compute y = M_n-1E_n-1....M_1E_1b. See [GOLUB/VANLOAN, § 3.4.3., p. 112]

 vec<T> d = b;
   //forward substitution [GOLUB/VANLOAN, p.113]. First permute the right hand side vector b according to the permutation pivot done with the LU factorization. Then update b
 
 vec<int> null(n);
 if(piv  != null){   //if pivot vector is zero then no permutation was applied so far and we need no swapping 
   for(int k = 0; k < n; k++){
     //std::swap(b[k], b[piv[k]]);
     d[k] = b[piv[k]];
     }   
     b = d;  
  }

   //update y = M_n-1E_n-1.....M_1E_1b
   for(int k = 0; k < n; k++){
      for(int j = k+1; j  < n; j++)
        b[j] -= b[k]*A[j][k];
   } 

 
    //backward substitution Ux = y
  vec<T> x(row);
  for(int i = n-1; i >= 0; i--){
    for(int j = i+1; j < m; j++)
      b[i] -= A[i][j]*x[j];
    //in many cases such as determining eigenvectors it is desirable to perform a Gaussian elimination even if the matrix is singular 
    if(A[i][i] == 0.0){
      std::cout<<"WARNING: Singular Matrix in >>matrix<T>::gausselimination<<. "<<std::endl;
      A[i][i] = TINY;
    }
      x[i] = b[i]/A[i][i];
  }
  b = x;  
 

 return b;
  
}


//GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING (n times n case). The ride hand side vector b is later overwritten with the solution vector x of the system A*x = b.
template<class T>
vec<T> matrix<T>::gausselimination(vec<T>& b) const{
  //if(row != column) error("Guarantee me a square matrix in >>matrix<T>::gausselimination(vec<T>&)<<. ");
  if(column != b.size()) error("RHS vector size and coefficient matrix dimension do not match !"); 
 int n = row;
  matrix<T> A = *this;
  
  //LU decomposition; check if pivoting can be avoided
  vec<int> pivot(n);
  T d;
   matrix<T> B = A.transpose();
  bool crucial = B.sdd();  //Is transpose of A stricly diagonally dominant ?
  if(crucial == true)  A =  ludecomp(A); //no pivoting is needed !!
  else A = A.lufac(pivot, d);//...or with pivoting
   
  vec<T> sol = A.forbackwardsub(pivot, b);
  

  return sol;
}

template<class T>
vec<T> matrix<T>::gausselimpartpivot(vec<T>& b) const{
  if(row != column) error("I only consider the n times n case here in >>matrix<T>::gausselimpartpivot<<. ");
 
  if(row != b.size()) error("Vector size and coefficient matrix dimension do not match in >>matrix<T>::gausselimpartpivot<<. ");
  
  matrix<T> A = *this;
  int n = A.row; int m = A.column;

  //check for strict diagonal dominance
   matrix<T> B = A.transpose();
   bool decide = B.sdd();
  
  if(decide == true){//strict diagonal dominance
   for(int k = 0; k< n-1; k++){
    if(A[k][k] == 0) error("Pivot equals 0 in >>matrix<T>::gausselimpartpivot(vec<T>&)<< . "); 
    for(int i = k+1; i<n; i++){       
      if(A[i][k]!=0){  
	 A[i][k] /= A[k][k];
	 for(int j = k+1; j<m; j++){   
         A[i][j] -= A[i][k]*A[k][j];
         }
      }
    }
   }   
   //forward substitution for Ly =b. (y stored in b)
   for(int i = 1; i < n; i++)
     for(int j = 0; j < i; j++) 
       b[i] -= A[i][j]*b[j];

   //back substitution for Ux = y. (x stored in b)
   for(int i = n-1; i >= 0; i--){ 
     for(int j = i+1; j < n; j++)b[i] -= A[i][j]*b[j];
     b[i] /= A[i][i];  
   }   
   return b;
   }//end if:   
  

  else if(decide == false){

   vec<int> piv(n);      //pivot(row) vector (contains row numbers)
  for(int k = 0; k < n; k++){   //..and stores swapping information
    piv[k] = k;                 // piv=[0,1,2,......,n-1]   
  }                                 
  
  for(int k = 0; k < n-1; k++){      //here the loop starts
    int seek = k;                    //find entry in column k, row pvt[k]
    T mag = abs(A[piv[k]][k]);       // with largest magnitude mag & take it as
    
    for(int i = k+1; i < n; i++){   // pivot element = max.|column entry|
      if(abs(A[piv[i]][k]) > mag){  //found such element. 
	mag = abs(A[piv[i]][k]);    //pivot search
      
      seek = i;  
      }
        
    }
    if(!mag) error("pivot is ZERO in lu(matrix<T>&).");
    if(seek != k) {
      //ONLY HERE: the rows of A need NOT to be interchanged in order to reduce the operations. The permutation information is solely stored  in vector "piv" to which is later referred in the backward and forward substititution !
      permute(piv[k], piv[seek]);
    }    
     for(int i = k+1; i < n; i++){       //entry elimination below A[pvt[k]][k]
      if(A[piv[i]][k] != 0){
        A[piv[i]][k] /=  A[piv[k]][k]; //"zero out" entries in the column
        for(int j = k+1; j < m; j++)
	  A[piv[i]][j] -= A[piv[i]][k]*A[piv[k]][j];
      }
    }
   }
  
  //check if determinant is zero (here just verify if product of diagonal entries is zero regardless the sign). If yes there exist no solution for the matrix is defective.
  T dia = 1.0;
  for(int i = 0; i < n; i++) dia *= A[i][i];
  if(dia == 0.0) error("Matrix is singular -- at least within my precision. Hence there is no solution when applying  >>matrix<T>::gausselimpartpivot<< !"); 

  //forward substitution Ly = Pb
  for(int i = 1; i < n; i++){
    for(int j = 0; j < i; j++){
      b[piv[i]] -= A[piv[i]][j]*b[piv[j]];
    }
  }
  //back substitution Ux = y
  vec<T> x(n); //stores solution in correct order 
  for(int i = n-1; i >= 0; i--){
    for(int j = i+1; j < m; j++) b[piv[i]] -= A[piv[i]][j]*x[j];
    //if the pivot element is zero then we assign a tiny number to the pivot element which is desirerable for some applications. 
    if(A[piv[i]][i] == 0.0){
       std::cout<<"WARNING: I have to cope with a singular matrix in >>matrix<T>::gausselimpartpivot(vec<T>&)<<. Hence I assigned a tiny number as pivot element. "<<std::endl;
      A[piv[i]][i] = TINY;;


      // error("Matrix is SINGULAR in >>matrix<T>::gausselimpartpivot(vec<T>&)<<. ");
    }
    x[i] = b[piv[i]]/A[piv[i]][i];
  }
  b = x; 

  return b;
 } //end else

}

/*Complete pivoting is more robust than partial and yields better error bounds. However, when used with large systems its perfo  T duenn;
  if(n == 1) return A[0][0];
  if(n == 2){
    duenn = A[0][0]*A[1][1];
    return duenn;
  }rmance is quite inacceptable thus favouring partial pivoting to complete [--> Golub/VanLoan].  */

/*template<class T>            //LU decomposition with complete pivoting
matrix<T>lucp(matrix<T>& A){  //void matrix<T>::lu(matrix<T>& A){
  int n,m; 
  n= A.row; m= A.column;
  if(n != m) error("Matrix is not SQUARE.");
 
  vec<int> px(n);
  vec<int> qy(n);
  for(int k = 0; k < n; k++) px[k] = qy[k] = k;template<class T>  //if an entry is below a certain bound it is regarded as 0
void  matrix<T>::set2zero() const {
  for(int i= 0; i < row; i++){
    for(int j =0; j < column; j++){
      if(abs(entry[i][j]) <= SMALL) entry[i][j] = 0.0;
    }
  }
}

     
  for(int k = 0; k < n-1; k++){
    int pct = k, qdt = k;
    double aet = 0;
    for(int i = k; i<n;i++){
      for(int j=k; j < n;j++){
	double tmp = abs(A[px[i]][qy[j]]); 
        if(tmp > aet) {aet = tmp; pct =i; qdt = j;}
      }
    }
    if(!aet) error("Pivot = 0.");;
    std::swap(px[k], px[pct]);
    std::swap(qy[k], qy[qdt]);

    for(int i = k+1; i< n;i++){
      if (A[px[i]][qy[k]] !=0){
	T mult = A[px[i]][qy[k]]/A[px[k]][qy[k]];
        A[px[i]][qy[k]] = mult;
        for(int j = k+1; j<m; j++)
	  A[px[i]][qy[j]] -= mult*A[px[k]][qy[j]];
      }
    }
  }
  return A;
  }*/
  


template<class T>               //Determinant calculation after having invoked
T det(matrix<T>& A){             
  int n,m; 
  n= A.row; m= A.column;
  if(n != m) error("Matrix must always be SQUARE for determinant calculation.");  //if(n == 1) return A[0][0]; 
  
 vec<int> pivot(n);
  T d;
 
  matrix<T> B = A.transpose();
  bool crucial = B.sdd();  //Is transpose of A stricly diagonally dominant ?
  if(crucial == true){
          A = ludecomp(A); //no pivoting is needed !!
          d = 1.0;
  }
  else A = A.lufac(pivot, d);//lupp(A);
  
  //actual determinant calculation
 for(int k = 0; k < A.row; k++){ 
    d *= A[k][k];           //ATTENTION: det(A) = +/- A[0][0]*...*A[n-1][n-1]
    } 
  return d;
 }


//Rank of a matrix, calculated via "HOUSEHOLDER QR WITH COLUMN PIVOTING". Given A in R^(m*n), m>=n, the following algorithm computes r = rank(A). See [GOLUB/VANLOAN, § 5.4.1, algo 5.4.1, pp 249]
template<class T> 
int matrix<T>::rank() const{
   int m = row; 
   int n = column;
   int r = 0; //this is the rank
   if(m < n) error("I expect row >= column in >>matrix<T>::rank()<<. Just transpose the matrix since Rank(A) = Rank(A^T) !"); 
  
  matrix<T> A = *this;

  vec<T> c(n);
  vec<T> tm(m);  
  for(int j = 0; j < n; j++){
    for(int i = 0; i < m ; i++){
      tm[i] = A[i][j];   
    }
    c[j] = tm*tm;
  } 
  
  int k = 0; //find the smallest 0 <= k < n such that c[k] = max{c[0],.,c[n-1]}
  T tau = c[0];
  for(int i = 1; i < n; i++){
    if(c[i] > tau){
      tau = c[i];
      k = i;
    }
  }
 
  
  int dim = int(std::min(m,n));
   vec<int> piv(dim);
   
  while(tau > SMALL){   //SMALL = 1.e-13, i.e. "nearly zero"
     
     //for(r = 0; r < dim; r++){  //alternative loop
     
       /* if(tau <= SMALL){   //Marc's stopping criterion
       //r = k;  //maybe not that necessary
       break;
       }*/
       // r++;           //put this at the bottom of this loop !!
        piv[r] = k;
      for(int i = 0; i < m; i++) permute(A[i][r], A[i][k]);
    permute(c[r], c[k]);

    
    //Householder QR
     vec<T> xf(m-r);
    for(int i = 0; i < m-r; i++){ 
      xf[i] = A[i+r][r];
    }

       vec<T> ho(m-r+1);
    ho = house(xf);
    T beta = ho[m-r];
    
 
     vec<T> v(m-r);    //Householder vector
     for(int i = 0; i < m-r; i++){
       v[i] = ho[i];
     }
    
     vec<T> essv(m-r-1); //essential part of Householder vector      
     for(int i=0; i < m-r-1; i++){
       essv[i] = v[i+1]; 
     }


     
     matrix<T> B(m-r, m-r);     //Dyadic product of Householder vector
     B =  beta*dyad(v,v);
     
     matrix<T> C(m-r, n-r);      //matrix A(j:m,j:n)
       for(int i = 0; i < m-r; i++){        
        for(int l = 0; l < n-r; l++){
	  C[i][l] = A[i+r][l+r];    
        }
       }

       C -= B*C;   //update
     
     for(int i = r; i < m; i++){       //update is assigned back to A  
        for(int l = r; l < n; l++){
	 A[i][l] = C[i-r][l-r];    
        }
     }


     for(int i = r+1; i < m; i++)
        A[i][r] = essv[i-(r+1)];

     ///////////////////
  
     for(int i = r+1; i < n; i++)
       c[i] -= pow(A[r][i],2);
  
     if(r < n){
      tau = c[r+1];
      k = r+1;
      for(int i = r+2; i < n; i++){
        if(c[i] > tau){
        tau = c[i];
        k = i;
       }
      }
     }
     else tau = 0;  

     r++;    //increment rank if necessary
 }//end while-loop  

  return r;
}

//rank determination in case that m < n occurs
template<class T> 
int rank(matrix<T>& A){
  int m = A.row;
  int n = A.column;
  int ra;
  if(m >= n){
     ra = A.rank();
  }
  else if (m < n){
    matrix<T> B = A.transpose();
    ra = B.rank();
  }
  return ra;
}


template<class T>
bool matrix<T>::symm() const{
  if(row != column) error("Symmetry is only defined for square matrices. ");
  bool george = true;
  for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      if(i != j){
	if(entry[i][j] != entry[j][i]){
	  george = false;
          break;
        } 
      }   
    }
  }
  return george;
}

template<class T>
T tr(matrix<T>& A){
  if(A.row != A.column) error("Trace is only defined for square matrices. ");
  T trace = A[0][0];
  for(int i = 1; i < A.row; i++) trace += A[i][i];
  return trace;
}


/*Partition of a matrix into q times r blocks of more or less equal size. Returns a vector containing integer data of number of equally sized blocks w.r.t. block row dimension (N1),  number of such row blocks (mul) and the remainig row dimension (sigma). (Anologuously for column which constitutes entry 2 to 5 of the returned vector.*/

template<class T>
vec<int> matrix<T>::partition(int q, int r){ //q times r processes
  int m = row;
  int n = column;
  if(q>m || r > n) error("Cannot partition into greater blocks than dimensions of the matrix are !"); 

  T vd = m/T(q); //Problem: quotient of 2 ints is often not treated as double
  int vf = (int)vd;
  if(abs(vd-vf) >= 0.5) {vd = ceil(vd);}  //rounding
  else vd = floor(vd);

  int N1 = (int)vd;           //m = N1*mul + sigma
  
  int mul; 
  int sigma;
  if(N1 == 1){
    mul = q-1;
    sigma = m - mul*N1;
  } 
  else{
  sigma = m%N1;
   mul = (m-sigma)/N1;
  }

  //analoguously for column:

  T cd = n/T(r);           
  double cf = (int)cd;
  if(abs(cd-cf) >= 0.5) cd = ceil(cd);
  else cd = floor(cd);

  int N2 = (int)cd;          //n = N2*kul + tau
  
  int tau;
  int kul;

  if(N2 == 1){
    kul = r-1;
    tau = n - kul*N2;
  } 
  else{
  tau = n%N2;
   kul = (n-tau)/N2;
  }
  //int tau = n%N2;
  //int kul = (n-tau)/N2;

  vec<int> sp(6);   //stores N1, mul, sigma, N2, kul & tau
  sp[0] = N1; sp[1] = mul; sp[2] = sigma;
  sp[3] = N2; sp[4] = kul; sp[5] = tau;

  return sp;
}




template<class T> 
 matrix<T> blockadd(const matrix<T>& A, const matrix<T>& B){
  if(A.row != B.row || A.column != B.column) error("Check dimensions in blockadd( const matrix<T>&, const matrix<T>&) !! ");
  
 int m = A.row;
 int n = A.column;
 matrix<T> C(m,n);
 //matrix<T> C = *this;  //NON-FRIEND VERSION 

  int q = q_part;   //defined in 'param.h'
  int r = r_part; 
  matrix<T> D = A; 
  matrix<T> E = B; 
  
  vec<int> ju = D.partition(q,r);
  vec<int> store = E.partition(q,r);
  if(ju != store) error("Not equally partitioned.");
  
  /* vector ju =[dim of rowblock, number of that block, dimension of row rest; dim of columnblock, number of that block, dimension of column rest] */
  
  int N1 = ju[0];  //row dimension of block 
  int lam1 = ju[1];   //number of that block w.r.t. row
  int rest1 = ju[2]; // remaining rows

   int N2 = ju[3];  //column dimension of block 
  int lam2 = ju[4];   //number of that block w.r.t. column
  int rest2 = ju[5]; // remaining rows

  if(m != lam1*N1+rest1) error("Hmmm. something's weird with the partition!");
  if(n != lam2*N2+rest2) error("Hmmm. something's weird with the partition!");
  //matrix<T> C(lam1*N1+rest1, lam2*N2+rest2);
   
  for(int alpha = 1; alpha <= lam1; alpha++){
    for(int i = (alpha-1)*N1 ; i < alpha*N1 + rest1; i++){
      //the same for columns
       for(int beta = 1; beta <=lam2; beta++){
         for(int j = (beta-1)*N2; j < beta*N2 + rest2; j++){ 
      
             C[i][j] = A[i][j] + B[i][j];
         }   
       }
    }
   }
  return C;
}

//as a member fucntion
template<class T>
matrix<T> matrix<T>::blockplus(const matrix<T>& A){
  if(row != A.row || column != A.column) error("Adjust dimensions .");
  int m = row;
  int n = column;  

matrix<T> C = *this;
  int q = q_part;  // defined in param.h
  int r = r_part; 
  matrix<T> D = C; 
  matrix<T> E = A; 
  
  vec<int> ju = D.partition(q,r);
  vec<int> store = E.partition(q,r);
  if(ju != store) error("Not equally partitioned.");
  
  // vector ju =[dim of rowblock, number of that block, dimension of row rest; dim of columnblock, number of that block, dimension of column rest] 
  
  int N1 = ju[0];  //row dimension of block 
  int lam1 = ju[1];   //number of that block w.r.t. row
  int rest1 = ju[2]; // remaining rows

   int N2 = ju[3];  //column dimension of block 
  int lam2 = ju[4];   //number of that block w.r.t. column
  int rest2 = ju[5]; // remaining rows

  if(m != lam1*N1+rest1) error("Hmmm. something's weird with the partition!");
  if(n != lam2*N2+rest2) error("Hmmm. something's weird with the partition!");
 
   
 //store dimensions and "remainders" in 2 seperate vectors (for clarity only)
 //Only vectors of dim 2 since matrix addition is per definitionem only for totally equally blocked matrices !!
 
 vec<int> nve(2);
  nve[0] = N1; nve[1] = N2; 
  vec<int> resve(2);
  resve[0] = rest1;  resve[1] = rest2;  
  vec<int> lave(2);
   lave[0] = lam1; lave[1] = lam2;
   
   vec<int> addi(2); //stores how many new rows and columns to be added to A, B  
/*GENERAL modification like V. STRASSEN proposed: if the remainders are less or greater (depends on partion parameters) then add zero rows and column such that the "expanded" rests are of the same dimensions as the corresponding N_i.   */

  //Ajust "rest" to the corresponding N_i:
   for(int i = 0; i< 2; i++) {
     if(resve[i] >= nve[i]){
      if(resve[i]%nve[i] == 0){  //if rest is a multiple of blocksize
        if(resve[i]/(double(nve[i])) == 1){
	  resve[i] = nve[i];  
        addi[i] = 0;  //equals block size -- no changes
	} 
        else if( resve[i]/(double(nve[i])) > 1){
	  int add = resve[i]/nve[i]; //i.e.:int 7/3 = 2, 10/5 = 2
	   resve[i] = add*nve[i]; 
           lave[i] += add;     //block number increases
           addi[i] = 0; //no additional 0-row/cols 
        }
      }  
      else{ //if modulo != 0, i.e. no multiple of blocksize
        //resve[i] += (nve[i]-resve[i]%nve[i]);} 
       int zu = (int)resve[i]/nve[i];
       zu++;    //zu+1
       lave[i] += zu;
       addi[i] = zu*nve[i] - resve[i]; //resve[i]%nve[i];
       resve[i] = zu*nve[i];
       }
      } 
   
   if(resve[i] < nve[i]){
     if(resve[i] == 0){
     addi[i] = 0;
     //lave[i] = lave[i];
     } 
     else{ lave[i]++;      //increment number of blocks
       addi[i] = nve[i]-resve[i];
       //resve[i] = 0;
     } 
    }
   }
  
  //expanded matrices(may conrtain zero rows and columns towards the end):
  matrix<T> Cex(row + addi[0], column + addi[1]);
  matrix<T> Aex(A.row + addi[0], A.column + addi[1]);
 
  //fill expanded matrices:
  for(int i = 0; i < A.row; i++){
    for(int j = 0; j < A.column; j++){
       Aex[i][j] = A[i][j];
     }
  }
  

 for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      Cex[i][j] = C[i][j];
    }
 }

// Reassign values:
 lam1 = lave[0]; lam2 = lave[1];
 

   
  for(int alpha = 1; alpha <= lam1; alpha++){
    for(int i = (alpha-1)*N1 ; i < alpha*N1; i++){
      //the same for columns
       for(int beta = 1; beta <=lam2; beta++){
         for(int j = (beta-1)*N2; j < beta*N2; j++){ 
      
             Cex[i][j] += Aex[i][j];
         }   
       }
    }
   }
  
  
 //reassignment to matrix without zeros rows/columns
  matrix<T> Sum(row, column);  
 
 for(int i = 0; i < row; i++){
    for(int j = 0; j < row; j++){
      Sum[i][j] = Cex[i][j];
    }
  }
return Sum;
 
}


template<class T>
matrix<T> blockmult(const matrix<T>& A, const matrix<T>& B){
  if(A.column != B.row) error("The matrices to multiply should have the same number of columns resp. rows !");
  matrix<T> C(A.row, B.column);

  //partition parameters are defined in "param.h" (included in "shit.cc")
  int q = qq_part;
  int r = r_part;
  int s = s_part;

  matrix<T> D = A; matrix<T> E = B;
  //columns of A and rows of B have to be equally partitioned !!!!!!!
     vec<int> au = D.partition(q,r); 
     vec<int> bu = E.partition(r,s);

  //for matrix A: 
  int N1 = au[0];  //row dimension of block 
  int lam1 = au[1];   //number of that block w.r.t. row
  int rest1 = au[2]; // remaining rows

   int N2 = au[3];  //column dimension of block 
  int lam2 = au[4];   //number of that block w.r.t. column
  int rest2 = au[5]; // remaining rows
 
  // the same for matrix B:
  int N3 = bu[0];  
  int lam3 = bu[1];   
  int rest3 = bu[2]; 

  int N4 = bu[3];  
  int lam4 =  bu[4];   
  int rest4 = bu[5]; 


  if(au[3] != bu[0] || au[4] != bu[1] || au[5] != bu[2]) error("Something's wrong with your partition. ");
  
  //store dimensions and "remainders" in 2 seperate vectors (for clarity only)
 vec<int> nve(4);
  nve[0] = N1; nve[1] = N2; nve[2] = N3; nve[3] = N4;
  vec<int> resve(4);
  resve[0] = rest1;  resve[1] = rest2;  resve[2] = rest3;  resve[3] = rest4; 
  vec<int> lave(4);
   lave[0] = lam1; lave[1] = lam2; lave[2] = lam3; lave[3] = lam4;
   
   vec<int> addi(4); //stores how many new rows and columns to be added to A, B  
/*GENERAL modification as V. STRASSEN proposed: if the remainders are less or greater (depends on partion parameters) than the row/column blocks then add zero rows and columns such that the "expanded" rests are of the same dimensions as the corresponding N_i.   */

  //Ajust "rest" to the corresponding N_i:
   for(int i = 0; i< 4; i++) {
     //if(resve[i] == 0) addi[i] = 0; 
    if(resve[i] >= nve[i]){
      if(resve[i]%nve[i] == 0){  //if rest is a multiple of blocksize
        if(resve[i]/(double(nve[i])) == 1){
	  resve[i] = nve[i];  
        addi[i] = 0;  //equals block size -- no changes
	} 
        else if( resve[i]/(double(nve[i])) > 1){
	  int add = resve[i]/nve[i]; //i.e.:int 7/3 = 2, 10/5 = 2
	   resve[i] = add*nve[i]; 
           lave[i] += add;     //block number increases
           addi[i] = 0; //no additional 0-row/cols 
        }
      }  
      else{ //if modulo != 0, i.e. no multiple of blocksize
        //resve[i] += (nve[i]-resve[i]%nve[i]);} 
       int zu = (int)resve[i]/nve[i];
       zu++;    //zu+1
       lave[i] += zu;
       addi[i] = zu*nve[i] - resve[i]; //resve[i]%nve[i];
       resve[i] = zu*nve[i];
        }
      } 
   
   if(resve[i] < nve[i]){
      if(resve[i] == 0){
      addi[i] = 0;
      //lave[i]=lave[i];
     } 
     else{ 
     lave[i]++;      //increment number of blocks
       addi[i] = nve[i]-resve[i];
       resve[i] = 0;
     }
   }
 }
  //expanded matrices(may conrtain zero rows and columns towards the end):
  matrix<T> Aex(A.row + addi[0], A.column + addi[1]);
  matrix<T> Bex(B.row + addi[2], B.column + addi[3]);
 
  //fill expanded matrices:
  for(int i = 0; i < A.row; i++){
    for(int j = 0; j < A.column; j++){
       Aex[i][j] = A[i][j];
     }
  }
  

 for(int i = 0; i < B.row; i++){
    for(int j = 0; j < B.column; j++){
      Bex[i][j] = B[i][j];
    }
 }



 //blockmultiplication of Aex and Bex with result Cex
 matrix<T> Cex(A.row + addi[0], B.column + addi[3]);

 // Reassign values:
 lam1 = lave[0]; lam2 = lave[1]; lam3 = lave[2]; lam4 = lave[3];
 

 /*for(int alpha = 1; alpha <= N1; alpha++){
    for(int i = (alpha-1)*lam1 ; i < alpha*lam1; i++){ //rowblocks of A
      for(int beta = 1; beta <= N4; beta++){
        for(int j = (beta-1)*lam4 ; j < beta*lam4; j++){//col.blx. of B
          for(int gamma = 1; gamma <=N2; gamma++){
	    for(int k = (gamma-1)*lam2; k < gamma*lam2; k++){//col.blx  A

	      Cex[i][j] += Aex[i][k]*Bex[k][j];
     
            }
          }
        }
      }
    }
    }*/
   

 for(int alpha = 1; alpha <= lam1; alpha++){
    for(int i = (alpha-1)*N1 ; i < alpha*N1; i++){ //rowblocks of A
      for(int beta = 1; beta <= lam4; beta++){
        for(int j = (beta-1)*N4 ; j < beta*N4; j++){//col.blx. of B
          for(int gamma = 1; gamma <=lam2; gamma++){
	    for(int k = (gamma-1)*N2; k < gamma*N2; k++){//col.blx  A

	      Cex[i][j] += Aex[i][k]*Bex[k][j];
     
            }
          }
        }
      }
    }
   }
   
 //matrix without added zeros -- use this matrix for further computations
  matrix<T> Result(A.row, B.column); 
  for(int i = 0; i < A.row; i++){
   for(int j = 0; j < B.column; j++){
     C[i][j] = Cex[i][j];
   }
  } 

//return Aex;

  return C;
}


template<class T>  //if an entry is below a certain bound it is regarded as 0
void  matrix<T>::set2zero() const {
  for(int i= 0; i < row; i++){
    for(int j =0; j < column; j++){
      if(abs(entry[i][j]) <= SMALL) entry[i][j] = 0.0;
    }
  }
}

/* Let Data be a real data matrix from whose column vectors (for example concentration over time [=rows]) the covariance matrix Sigma is to be calculated. Remark: Sigma is a symmetric column times column matrix ! */

template<class T> 
matrix<T> matrix<T>::Sigma() const{ 
 vec<T>** d = new vec<T>* [column];   //how many vectors
  for(int i = 0; i < column; i++){
    d[i] = new vec<T>(row);
  }
 
  vec<T> x(row);
  for(int j = 0; j < column; j++){
    for(int i = 0; i < row; i++){
      x[i] = entry[i][j];   //columns of the data matrix   
   }
    *d[j] = x;        //columnvecs of which covariance is to be computed
  }
  matrix<T> S(column, column);
  vec<T> ersatz1(row), ersatz2(row);

  for(int i = 0; i < column; i++){
    for(int j = 0; j < column; j++){
      ersatz1 = *d[i];
      ersatz2 = *d[j];
      S[i][j] = ersatz1.covar(ersatz2);   //calculate and assign covariances
    }
  }
 
  //deallocate space after use:
  for(int i = 0; i < column; i++) delete d[i];
  delete[] d;
  
  return S;
}

//cutting off (EXCLUSIVE column/row 'scissors')
template<class T>
matrix<T> matrix<T>::cut_off_from_column(int scissors) const{
    if(scissors >= column) error("Your specified column exceeds matrix dimension in >>matrix<T>::cut_off_from_column<< ! ");
    matrix<T> Cut(row,scissors);
    for(int i = 0; i<row; i++)
	for(int j = 0; j <scissors; j++)
            Cut[i][j] = entry[i][j];

    return Cut;
}


template<class T>
matrix<T> matrix<T>::cut_off_from_row(int scissors) const{
    if(scissors >= row) error("Your specified column exceeds matrix dimension in >>matrix<T>::cut_off_from_row<< ! ");
    matrix<T> Cut(scissors,column);
    for(int i = 0; i<scissors; i++)
	for(int j = 0; j < column; j++)
            Cut[i][j] = entry[i][j];

    return Cut;
}

template<class T>
  matrix<T> matrix<T>::delcolumn(int c) const{
    if(c >= column) error("Your specified column exceeds matrix dimension in >>matrix<T>::delcolumn<< ! ");
    int ncol = column-1; 
    matrix Red(row,ncol);
    for(int i = 0; i< row; i++){
	for(int j = 0; j < ncol; j++){
	    if(j >= c && c != ncol){   
               Red[i][j] = entry[i][j+1];
            }
            else Red[i][j] = entry[i][j];
        }
    }
    return Red;
}
template<class T>
  matrix<T> matrix<T>::delrow(int c) const{
    if(c >= row) error("Your specified column exceeds matrix dimension in >>matrix<T>::delrow<< ! ");
    int nro = row-1; 
    matrix Red(nro,column);
    for(int i = 0; i< nro; i++){
	for(int j = 0; j < column; j++){
	    if(i >= c && c != nro){   
               Red[i][j] = entry[i+1][j];
            }
	    else  Red[i][j] = entry[i][j];
        }
    }
    return Red;
}


//skip stuff from specific row/column (inclusive !), i.e. entries of a matrix whill be ignored before row/column 'c'
template<class T>
matrix<T> matrix<T>::skip_from_row(int c) const{
    if(c >= row) error("Row index exceeds matrix dimensions in >>matrix<T>::skip_from_row<<");
   if(c == 0) return *this; 
   else{
       matrix<T> N(row-c,column);
       for(int i = c; i < row; i++){
	   for(int j = 0; j < column; j++){
	       N[i-c][j] = entry[i][j];
	   }
       }
       return N;
   }
}

//the same for column
template<class T>
matrix<T> matrix<T>::skip_from_column(int c) const{
    if(c >= column) error("Column index exceeds matrix dimensions in >>matrix<T>::skip_from_column<<");
    if(c == 0) return *this;
    else{
	matrix<T> N(row,column-c);
	for(int i = 0; i < row; i++){
	    for(int j = c; j < column; j++){
		N[i][j-c] = entry[i][j];
	    }
	}
	return N;
    }
}





template<class T>
void matrix<T>::mat2txt(const std::string& word) const{
  std::ofstream outfile(word.c_str(), std::ios_base::out);
    outfile<<(*this);
    outfile.close();
}


//this routine calculates the sums (or mean values etc.)  of the matrix row blocks (a block is a part of column c (specified as argument), that contains the same values, e.g. several (indep.) measurements for the same time point t_i (= one block).   
//All the sums of these blocks along all colums, exept c, are computed in
// block style !!
//Note: the member function wolpertinger, e.g.,  sums up all multiple measurements done at the same time point.
//the list contains to matrices: the first one the experimental measurements, the second the calculated sample variances

template<class T>
//martix<T> matrix<T>::wolpertinger(int c) const{
std::list<matrix<T> > matrix<T>::wolpertinger(int c){
    if(row == 0 || column == 0) error("Nothing to do in >>matrix<T>::wolpertinger. "); 
  
//ecvalute different time steps and the dimension of each time block
//stored in simpletime and dims (in order) 
  int subdim = 1, newdim = 1;
   int index = 0;


   vec<int> dims(row);
   vec<T> simpletime(row);   //everything exept maybe int

    for(int i = 0; i < row-1; i++){
      if(entry[i][c] == entry[i+1][c]){
	  subdim++;   
      }
      else{
	  dims[index] = subdim;
          simpletime[index] = entry[i][c]; 
	  //since of the first if statement we cannot incremtent fully
          
          
          subdim = 1; //reset dimension for next block
           index++;
          newdim++; 
      }
       if( i == row-2){
          simpletime[index] = entry[row-1][c];         
          dims[index] = subdim;  
          break;
       }
      
    }

    ///////above 99.99999 % correct	  
    
   
    int colmone = column -1;
   
    vec<T>* dragonfly = new vec<T> [newdim]; //products (means) of each block
    vec<T>* firefly = new vec<T> [newdim]; //newdim = new row dim

    vec<T> prodsto(colmone), varsto(colmone);

      int plusdim = 0;   

     //main iteration
      for(int dx = 0; dx < newdim; dx++){ //row = newdim !
       vec<T> temp(dims[dx]);
       for(int k = 0; k < colmone; k++){  //colum-1 (first col = time) !

	  //fill temp
          for(int l = 0; l < dims[dx]; l++){
            temp[l] = entry[l+plusdim][k+1]; //col from 1 to column-1       
         
          }
      
	  prodsto[k] = temp.meanvalue();    //temp.sumup();
                      //temp.prodvalue();   //first block column multiplied
       
	  varsto[k] = temp.variance();    //variances as well 
       }
       dragonfly[dx] = prodsto;
       firefly[dx] = varsto;

       plusdim += dims[dx]; //next row index in matrix (*this)
     }

   
    //products calculated............Let's hope
    //fill somehow new matrix

 matrix<T> Newmat(newdim, column); 
 matrix<T> Varmat(newdim, colmone);
//matrix<T> HALLO(3, 2);

 for(int i = 0; i < newdim; i++){ 
     vec<T> mandrake = dragonfly[i];  
     vec<T> vars = firefly[i];
     
     for(int j = 0; j < column; j++){ 
       if(j == 0){
           Newmat[i][j] = simpletime[i];     //time in first col
       }
       else{
	   Newmat[i][j] = mandrake[j-1];      //the rest the added etc. 
	   Varmat[i][j-1] = vars[j-1]; 
       }
     }
 }
    
 delete[] dragonfly;
 delete[] firefly;

 std::list<matrix<T> > Stoi;   //empty list
 Stoi.push_back(Newmat);  //add matrix at the end of list
 Stoi.push_back(Varmat); 
//Stoi.push_back(HALLO); //add Varmat at end of list: now: Newmat, Varmat

 return Stoi;
 
// return  Newmat; //dragonfly[1]; //o.k.
  
 
}



//orders w.r.t. stimuli (perturbation = other experiment). 

template<class T>
std::list<matrix<T> > matrix<T>::order_stimuli(int parcol){
    if(parcol >= column) error("Your specified column violates matrix range in >>matrix<T>::order_stimuli<<. ");
    
    std::list<matrix<T> > Lstore; //empty list 

    vec<T> col = (*this).getcolumn(parcol), slut = col;
    vec<vec<T> > moo = col.find_different_entries();
    vec<T> so = moo[0], di = moo[1];  //different entry vector, multiplicity of entry
    int sosi = so.size(),  colmone = column-1;   //disi = di.size()
    
    //matrix<T>* MS = new matrix<T> [sosi];
    
    matrix<T> Red = (*this).delcolumn(parcol);  //delete stimuli column    

    //main loop
    for(int i = 0; i < sosi; i++){  //number of different entries 
	vec<int> pos = slut.position(so[i]); //position(s)of 1st entry (corr. size of di)
        int dis = int(di[i]);
	matrix<T> S(dis, colmone);  //temporary matrix without stimulus column
	for(int l = 0; l < dis; l++){     //multiplicity of entry x
	    for(int j = 0; j < colmone; j++){
		S[l][j] = Red[pos[l]][j]; 
   
	    }
	}
	    //} 
	//MS[i] = S; 
	Lstore.push_back(S);
    }
    
    //delete[] MS;
    
    return Lstore;
    
}



template<class T>
void matrix<T>::write_2_file(const std::string& filename, const std::string& comment) const{
  std::ofstream outputTheStuff(filename.c_str(), std::ios_base::out);
  outputTheStuff << "#"<<comment << std::endl;  
  outputTheStuff << (*this) << '\n';
}


template<class T>
matrix<T> matrix<T>::reverse_by_rows(){
    int rmone = row-1;
    matrix<T> Aux(row, column);
    for(int i = rmone; i >= 0; i--){
	for(int j = 0; j < column; j++)
	    Aux[rmone-i][j] = entry[i][j];
    } 
    return Aux;
}


//a normal array: sizeof(array)/sizeof(array[0]) but not with matrix pointer sinze otherwise column index or whatever is called...
template<class T>
int size_of_pointer(matrix<T>* Mpoint){
    return sizeof(Mpoint);
}

//you specify an integer vector 'cidx' denoting the columns which shall be extracted from a given matrix thus yielding a new n times (cidx.size()) matrix... 
template<class T> 
matrix<T> matrix<T>::choice_by_columns(const vec<int>& cidx){
    int vlen = cidx.size();
    if(vlen > column)
	error("In >>matrix<T>::choice_by_columns(const vec<int>&)<<: \nargument vector contains MORE columns than corresponding matrix object!\n");

    matrix MNew(row, vlen);
    for(int i = 0; i < row; i++){
	for(int j = 0; j < vlen; j++){
	    if(cidx[j] >= column){
	       error("In >>matrix<T>::choice_by_columns(const vec<int>&)<<: \nargument vector contains column(s) which exceed(s) column dimension \nof matrix object!\n");
	    }
	    MNew[i][j] = entry[i][cidx[j]];
	}
    }
    return MNew;
}

//you specify an integer vector 'ridx' denoting the rows which shall be extracted from a given matrix thus yielding a new (ridx.size()) times n  matrix... 
template<class T> 
matrix<T> matrix<T>::choice_by_rows(const vec<int>& ridx){
    int vlen = ridx.size();
    if(vlen > row)
	error("In >>matrix<T>::choice_by_rows(const vec<int>&)<<: \nargument vector contains MORE rows than corresponding matrix object!\n");

    matrix MNew(vlen, column);
    for(int i = 0; i < vlen; i++){
	if(ridx[i] >= row){
	       error("In >>matrix<T>::choice_by_rows(const vec<int>&)<<: \nargument vector contains row(s) which exceed(s) row dimension \nof matrix object!\n");
	    }
	for(int j = 0; j < column; j++){
	    
	    MNew[i][j] = entry[ridx[i]][j];
	}
    }
    return MNew;
}



template<class T>
void matrix<T>::refill_column(vec<T> v, int co){
if(v.size() != row) error("Vector size and matrix row size do not match when using >>matrix<T>::refill_column(vec<T>, int)<<.  ");  
if(co >= column) error("Your specified column index actually exceeds column size in >>matrix<T>::refill_column(vec<T>, int)<<. ");
for(int i = 0; i < row; i++)
    entry[i][co] = v[i];
}

template<class T>
void matrix<T>::refill_column(T d, int co){  
if(co >= column) error("Your specified column index actually exceeds column size in >>matrix<T>::refill_column(vec<T>, int)<<. ");
for(int i = 0; i < row; i++)
    entry[i][co] = d;
}

template<class T>
void matrix<T>::refill_row(vec<T> v, int ro){
if(v.size() != column) error("Vector size and matrix column size do not match when using >>matrix<T>::refill_row(vec<T>, int)<<.  ");  
if(ro >= row) error("Your specified row index actually exceeds column size in >>matrix<T>::refill_row(vec<T>, int)<<. ");
for(int j = 0; j < column; j++)
    entry[ro][j] = v[j];
}

template<class T>
void matrix<T>::refill_row(T d, int ro){  
if(ro >= row) error("Your specified row index actually exceeds row size in >>matrix<T>::refill_row(vec<T>, int)<<. ");
for(int j = 0; j < row; j++)
    entry[ro][j] = d;
}

//return all matrix dimensions simultaneously
template<class T>
vec<int> matrix<T>::dimensions(){
    vec<int> two(2);
    two[0] = row;
    two[1] = column;
    return two;
}

//The following algorithm first creates 3 int vectors out of given bounds. Hereafter it calculates the mean values of all indices which are contained in these 3 vecs and assigns the output to an 3 times column matrix. This procedure is especeially fruitful when dealing with spatio-temporal models (NUC, CYTOPLASM, PLASMAMEMBRANE). Beyond that I think it will hardly find any application ;-)
//lower and upper bounds as argument
template<class T>
matrix<T> matrix<T>::reduce_to_3_rows(int l1, int u1, int l2, int u2,int l3, int u3){
    if(u1 == l2 || u2 == l3 || u3 == l1)
	error("In >>matrix<T>::reduce_to_3_rows(...)<<: Some argument bounds coincide. \n Instead the bounds should read from left to right with the left bounds SMALLER than the right ones. Otherwise, unpredictable results my occur !\n Bounds for the same subvec may be equal but not lower and upper bounds of different subvecs !\n ");
     
   if((l1 >= row ||  u1 >= row || l2 >= row ||  u2 >= row || l3 >= row ||  u3 >= row)||(l1 < 0 ||  u1 < 0 || l2 < 0 ||  u2 < 0 || l3 < 0 ||  u3 < 0))
	error("In >>matrix<T>::reduce_to_3_rows(...)<<: Hey hunk, all argument bounds should stay within the ROW range of the matrix!\n");	
  
    vec<int> index(row);
    for(int i = 0; i < row; i++)
	index[i] = i;

    //create subvectors from bounds
    vec<int> nuc = index.subvec(l1, u1), cyt = index.subvec(l2, u2), plm = index.subvec(l3, u3);
    
    matrix<T> Three(3,column);
    //sizes of subvecs
    int ns = nuc.size(), cs = cyt.size(), ps = plm.size();
   
    //needed to perform addition in main loop below
    T nent,cent,pent;
    
    //first add up corresponding entries with index from 3 subvecs
    //running time = O(column*(ns+cs+ps +3)) = O(column*(row+3))
    for(int j = 0; j < column; j++){ 
            
	nent = 0, cent = 0, pent = 0; //don't forget to reset for each column!
	
	for(int k = 0; k < ns; k++)
	    nent += entry[nuc[k]][j];
	for(int k = 0; k < cs; k++)
	    cent += entry[cyt[k]][j];
        for(int k = 0; k < ps; k++)
	    pent += entry[plm[k]][j];
	
    //then assign it (devided by size of subvec = meanvalue) to outputmatrix
	Three[0][j] = nent/ns;
	Three[1][j] = cent/cs;
	Three[2][j] = pent/ps;
    }
    return Three;
}

template<class T>
matrix<T> matrix<T>::extract_via_same_values(int co, T val){
    if(co >= column)
	error("In >>matrix<T>::extract_via_same_values(int, T)<<: Argument exceeds matrix dimension.\n");
    vec<T> v =(*this).getcolumn(co);
    vec<int> pos(row);
    int numb = 0;
    for(int i = 0; i < v.size(); i++){
	if(v[i] == val){
	    pos[numb++] = i;
	}
    }
    
    if(numb == 0)
	error("In >>matrix<T>::extract_via_same_values(int, T)<<: There is no value " + num2str(val) + " contained in column " + num2str(co) + "!\n");
    
    
    vec<int> s(numb);
    for(int i = 0; i < numb; i++)
	s[i] = pos[i]; 
    
    matrix<T> M(numb, column);
    
    for(int i = 0; i < numb; i++){
	for(int j = 0; j < column; j++){
	    M[i][j] = entry[s[i]][j];
	}
    }

    //std::cout<<numb<<std::endl;   //just for control
    return M;
}






/*************************************************************************************************************** CLASS DBLOCK ******************************************** GENERAL USAGE: A MATRIX OBJECT HAS TO BE NAMED IN ADDITION OR IDENTIFIED IN ANOTHER WAY.
MORE SPECIFIC USAGE: IF A FILE CONTAINS SEVERAL DIFFERENT DATA SERIES THEN THIS CLASS IS MEANT TO BE A DEVICE FOR STORING DATA AND A NAME ASSOCIATED WITH THE DATA ************************************************************************************************************************************************************/

template<class A_type, class B_type>
class dblock{
private:
    A_type signifier_;        //name of data, for instance, or whatever
    matrix<B_type> data_;     //contains values
public:
    dblock(A_type, matrix<B_type>);
    dblock(A_type = A_type(), B_type = B_type());  //default constructor

    dblock(const dblock&);
    dblock& operator=(const dblock&);
    
    ~dblock(){};
    
   
    std::ostream& print(std::ostream&) const; 

    A_type get_name();
    matrix<B_type> get_data();  //returns data

    B_type& operator()(int, int) const;
    void impose_name(std::string);
  
};

//constructor
template<class A_type, class B_type>
    dblock<A_type, B_type>::dblock(A_type s, matrix<B_type> M):signifier_(s),data_(M){}

//default constructor
template<class A_type, class B_type>
    dblock <A_type, B_type>::dblock(A_type s, B_type b):signifier_(s){
	for(int i = 0; i < data_.rock(); i ++)
	    for(int j = 0; j < data_.roll(); j++)
		data_[i][j] = b;
	

}


//copy constructor
template<class A_type, class B_type>
    dblock<A_type, B_type>::dblock(const dblock& beau):signifier_(beau.signifier_), data_(beau.data_){}

//copy assignment
template<class A_type, class B_type>
    dblock<A_type, B_type>& dblock<A_type, B_type>::operator=(const dblock& beau){
    if(this != &beau){
	signifier_ = beau.signifier_;
	data_ = beau.data_;
    }
    return *this;
}


//output stuff
template<class A_type, class B_type> 
    std::ostream& dblock<A_type, B_type>::print(std::ostream &ou) const {
    ou << "# Data of \"" << signifier_ <<"\":"<< std::endl;
    ou << data_;
    
 return ou;
} 

template<class A_type, class B_type>  //defined OUTSIDE the class !!
std::ostream& operator<<(std::ostream &ou, const dblock<A_type, B_type>& d){
  return d.print(ou);
}

//end ouputstuff


template<class A_type, class B_type>
    A_type dblock<A_type, B_type>::get_name(){
    return signifier_;
}
    
template<class A_type, class B_type>
    matrix<B_type> dblock<A_type, B_type>::get_data(){
    return data_;
}

template<class A_type, class B_type>
    B_type& dblock<A_type, B_type>::operator()(int i, int j) const{
    if(i >= data_.rock() || j >= data_.roll()){
	error("Index out of range when using >>dblock<A_type, B_type>::operator()(int, int)<<.");
    }
    return data_[i][j];
}

template<class A_type, class B_type>
    void dblock<A_type, B_type>::impose_name(std::string nn){
    signifier_ = nn;
}


#endif


