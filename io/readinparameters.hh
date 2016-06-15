#ifndef READ_IN_PARAMETERS_FROM_A_DISTINCT_FILE_HH
#define READ_IN_PARAMETERS_FROM_A_DISTINCT_FILE_HH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../marcvecmatrix/myfunctions.h"


namespace Adonis{

  /**
   * \brief Extracts values which are stored consecutively in some file 
   * and stores them in an array.
   *
   * NOTE: This is meant for a <BB> small and clearly arranged</BB> set of data
   *       only. Otherwise use class below.
   *
   * A typical file reads as follows:
   * \code 
      0.75
      1.55
      3.25
      0.001
      0.345
   * \endcode
   * 
   */
  template<class T, int N>
  class Parameter{
  private:
    T para_[N];
    
  public:
    typedef T value_type;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;

    const int size() const{return N;}

    Parameter(const T& d = T()){ //default constructor
      for(int i = 0; i < N; ++i)
	para_[i] = d;
    }

    template<class Iter>
    Parameter(Iter first, Iter last){
      adonis_assert((int)std::distance(first,last) == N);
      size_t idx = 0; 
      for(Iter it = first; it != last; ++it){
	para_[idx++] = *it;
      }
    }

    //construct from file
    Parameter(const std::string& file){
      std::ifstream infile(file.c_str(), std::ios_base::in);
      
      //adonis_assert(infile.is_open());
      if(!infile.is_open())
	ADONIS_ERROR(FileError,"File \""<<file<<"\" not found");


      std::string line;
      int idx = 0;

      while(!infile.eof()){
	getline(infile,line,'\n');

	if((!my_function_collection::is_whiteline(line)) && (!my_function_collection::is_comment(line))){
	  
	  para_[idx++] = static_cast<T>(atof(line.c_str()));
	 
	  }  
      }
      adonis_assert(idx == N);
    }

  
    const T& operator[](int i) const{
      adonis_assert((i >= 0) && (i < N));
      return para_[i];
    }

    T& operator[](int i){
      adonis_assert((i >= 0) && (i < N));
      return para_[i];
    }
    
    iterator begin(){
      return &para_[0];
    }
    
    iterator end(){
      return &para_[N];
    }
    
    const_iterator begin() const{
      return &para_[0];
    }
    
    const_iterator end() const{
      return &para_[N];
    }
    


    friend std::ostream& operator<<(std::ostream& os, const Parameter& P){
      os << " Parameter(s): ";
      for(const_iterator it = P.begin(); it != P.end(); ++it)
	os << *it << "  ";
    
      os << std::endl;
      return os;
  
    }
  
  };

  
  /**
   * \brief Detects and transforms the actual parameter
   */
  template<class S>
  inline void determine_read_value(S& ds, const S& line, unsigned iPlusOne){
    //save guards
    std::size_t countRadixPoints = 0;

    for(unsigned i = iPlusOne; i < line.size(); ++i){
      if(is_radix_point(line[i]))
	countRadixPoints++;
      if(countRadixPoints > 1)
	ADONIS_ERROR(IOError, "Too many radix points '.'");

      
      if(!my_function_collection::delimiter(line[i]) && !my_function_collection::is_comment(line[i])){
	ds.push_back(line[i]);

	//!o.k the next (i+1)st char turns out to be a 'whiteline' char while
	//!the current isn't. Then we stop further reading 
	if(i != line.size()-1 && (my_function_collection::delimiter(line[i+1]) || my_function_collection::is_comment(line[i+1]))){
	  //std::cout << "i = "<< i << std::endl;
	  break;
	}
	
	if((i  > 0) && (!my_function_collection::delimiter(line[i-1])) && (my_function_collection::delimiter(line[i])))
	  //i = line.size()-1;
	  //std::cout << "i = "<< i << std::endl;
	  break;  //ok leave loop for we've read in the value
      }
      
      //if comment encountered, then leave loop for rest of line isn't a param.
      if(my_function_collection::is_comment(line[i])){
	break;
      }
    }  
  }
  

 /**
   * \brief Reads in parameters from a distinct file format. 
   *A (weird) looking file might look as follows:
   * This is a file named <I>parameter.txt</I>
   * \code 
   #THIS IS A TEST PARAMETER FILE VERS 0.1 alpha
# a '#' denotes a comment and the line containing it will be rejected
	      
 
   c1:  0.245   #lower bracket value	

   		
   
   c2	:		17.2
   #tolerance
   tol: 1.e-3

   SIM- file : ~/MARC++/chemistrybasedmanifold/io/data/SIMDAT_441_E.dat # here I stored the files
   # shape parameter can be negative as well 4 it enters quadratically
   shap  ec:-0.75	 # gives also string "shapec"    
   beta_imq:                                                                     0.98
   *\endcode
   *
   * USAGE:
   * \code 
     ParameterData Data;
     Data.read_from_file("parameter.txt");

     std::cout << Data.get_datum<double>("c1") << std::endl; 
     std::cout << Data.get_datum<std::string>("SIM-file") << std::endl; 
   * \endcode
   *
   * Access time: since a binary tree is used to store the elements, we have
   *  O(\f$\log_2(n)\f$) for input size \f$n\f$
   */ 
class ParameterData{
public:
  static const unsigned PREC = 16;
  typedef std::string KeyType;
  typedef std::map<KeyType,KeyType> ContainerType;
  typedef ContainerType::iterator iterator;
  typedef ContainerType::const_iterator const_iterator;
  
  ParameterData():readIn_(false){}
  
    
  iterator begin(){
    return  Store_.begin();
  }
  
  iterator end(){
    return Store_.end();
  }
  
  const_iterator begin() const{
    return Store_.begin();
  }
  
  const_iterator end() const{
    return Store_.end();
  }
  
  unsigned size() const {return Store_.size();}
    
  bool has_been_read_in() const {return readIn_;}
    
  //! If a KeyType isn't present then the value is a zero string and you'll end up with an error thrown by  OutputTypeAdapter<X>::convert!
  //! \tparam X Type into which KeyType shall be transformed
  template<class X>
  X get_datum(const KeyType& kt) {
    //std::cout << "Show me whitespaces (if any): "; show_whitespaces(Store_[kt]);
#ifndef NDEBUG
    if(!(*this).found(kt)) //! if not found
      ADONIS_ERROR(DerivedError,"Parameter \""<<kt << "\" not found. Misspelled it!? Check it out, pal ;)");
#endif
    return OutputTypeAdapter<X>::convert(Store_[kt]);
  }
  
   
  //! overload previous function to take a default value
  template<class X>
  X get_datum(const KeyType& kt, const X& dflt){
    if(!(*this).found(kt)){ //! if not found
      ADONIS_INFO(Information,"Parameter \""<<kt << "\" not found.... Taking default value instead");
      return dflt;
    }
    else
      return OutputTypeAdapter<X>::convert(Store_[kt]);
  }

  //! read in parameters and store them with keyword in a map (binary tree)
  void read_from_file(const KeyType& file){
    if(!readIn_){  //only read in if it hasn't been read in so far
      KeyType line,
	key,
	ds;
      
      unsigned index = 0;
      
      std::ifstream readme(file.c_str(), std::ios_base::in);
      //adonis_assert(readme.is_open());
      if(!readme.is_open())
	ADONIS_ERROR(FileError,"File \""<<file<<"\" not found");

      while(!readme.eof()){
	getline(readme,line,'\n');
	
	key.clear(); //clear strings before they will get resized anew
	ds.clear();

	//! 'header' of each line, i.e. the parameter name including colon (:)
	for(unsigned i = 0; i < line.size(); ++i){
	  if(my_function_collection::is_comment(line[i])){//(rest of ) line is a comment 
	    break; 
	  }
	  if(!my_function_collection::delimiter(line[i]) && 
	     !my_function_collection::is_colon(line[i])){
	    key.push_back(line[i]);
	  }
	  
	  if(my_function_collection::is_colon(line[i])){
	    index = i;
	    break;
	  }
	}
	
	//now search and assemble parameter value in rest of 'line'
	determine_read_value(ds,line,index+1); 
	
	//! values are stored <I> alphabetically </I> w.r.t. key 'key'
	if(!my_function_collection::is_whiteline(key) && !my_function_collection::is_comment(key)){
	  Store_[key] = ds;    //!map need not to be resized ;)
	}
      }
      
      readme.close();
      
      readIn_ = true;   //read-in completed 
    }   
  }
        
  friend std::ostream& operator<<(std::ostream& os, const ParameterData& data){
      os << "------------------------------------------------------------------------------"<<std::endl;
      os << std::setprecision(PREC) << "# Key"<< '\t' <<"Value" <<std::endl <<"------------------------------------------------------------------------------"<<std::endl;
      for(const_iterator it = data.begin(); it != data.end(); ++it){
	os << (*it).first << '\t' << (*it).second << std::endl;
      }
      //os <<  std::setprecision(6); //"reset" output precison
      os << "------------------------------------------------------------------------------" <<std::endl;
      return os;
    }
    
  private:
    ContainerType Store_;
    mutable bool readIn_;

  //! not that map's find(key) delivers iterators. Unfortunately, you'll get a
  //! SEG. FAULT when using it for a non-existant entry in the binary tree.
  //! I circumvent this drawback by the following routine.
  bool found(const KeyType& key){
    bool fd = false;
    for(const_iterator cit = Store_.begin(); cit != Store_.end(); ++cit){
      if((*cit).first == key){ //! constant access via map's iterators
	//std::cout << (*cit).first << " ==> "<< (*cit).second << std::endl;
	fd = true;
	break;
      }
    }
    return fd;
  }
};
  
} //end of namespace 

#endif
