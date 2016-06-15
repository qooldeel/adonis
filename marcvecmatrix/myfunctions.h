/*CAUTION: 1) NEVER devide a template into .h and .cc files because there is the following problem when doing so: template stuff is invoked at compilation time and for that the complete class has to be known, i.e. especially the type T, etc. has to be known at compilation time. 
 2)Inline functions can (and should) be declared directly in .h file !!*/

#ifndef MYFUNC_H
#define MYFUNC_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip> //for 'setprecision' function

#include <stdio.h>    //to use function 'void exit(int status)'
#include <stdlib.h>  

#include "param.h"




namespace my_function_collection{ //define my own namespace

/************************** INLINE FCTS **********************************/
inline void error(const std::string& c){    //error(char* c){
    std::cout<<"\n"<<std::endl;
  std::cout<<"********  ERROR **** ERROR **** ERROR **** ERROR **** ERROR *******"<<std::endl;
  std::cout<<"\n"<<std::endl; 
  std::cout<<"    "<< c << " I'll stop my computations, Marc ! \n";
  std::cout<<"\n"<<std::endl;
  std::cout<<"*******  ERROR **** ERROR **** ERROR **** ERROR **** ERROR *******"<<std::endl;
  std::cout<<"\n"<<std::endl; 
  exit(1);
}

//a warning does not terminate computations like an error-invocation 
inline void warning(const std::string& c){  
    std::cout<<"\n"<<std::endl;
  std::cout<<"------- WARNING ----- WARNING ----- WARNING ----- WARNING -------"<<std::endl;
  std::cout<<"\n"<<std::endl; 
  std::cout<<"    "<< c << " FURTHER COMPUTATIONS MIGHT GET SPOILED ! \n";
  std::cout<<"\n"<<std::endl;
   std::cout<<"------- WARNING ----- WARNING ----- WARNING ----- WARNING -------"<<std::endl;
  std::cout<<"\n"<<std::endl; 
}

//counts the number of lines in a file, that are separated by '\n'
 inline int lines(std::string filename){
   std::ifstream countlines(filename.c_str(), std::ios_base::in);
   if(!countlines){
     std::cerr<<"Could not open inputfile "<<"'"<<filename<<"'"<<std::endl;
     exit(1);
   }
   std::string s;
   int count = 0;
   while(getline(countlines,s,'\n')){
     count++;
   }
   countlines.close();
   return count;
 }


inline bool find(std::string sg, std::string subsg){
    // if(subsg.size() > sg.size()) error("\"Sub\"std::string is longer than it's superstring in >>bool find(string, string)<< as defined in >>myfunctions.h<<.");
    size_t encountered;     //integral type of string
    encountered = sg.find(subsg);
    
    if(encountered != std::string::npos)
	return true;
    else return false;
}

inline int substrpos(std::string sg, std::string subsg){
     if(subsg.size() > sg.size()) error("\"Sub\"string is longer than it's superstring in >>bool find(string, string)<< as defined in >>myfunctions.h<<.");
    size_t encountered;     //integral type of string
    encountered = sg.find(subsg);
    
    if(encountered != std::string::npos)
      return int(encountered);
    else return -1;   //-1 for no string position !!
}

//replaces all occurances of substring 'sub' in whole string 'insert' with string 'rep' 
inline std::string substitute(std::string& str, const std::string& sub, const std::string& rep){
    size_t found;
    //std::string result;
    
    //sanity check
    //if(substrpos(str,sub) == -1) error("No occurence of substring at all in global function >>substitute<< as implemented in 'myfunctions.h'.");

    while((found = str.find(sub)) != std::string::npos){
	//result = str.replace(found,sub.size(), rep);
	str.replace(found,sub.size(), rep);
    }

    // return result;
    return str;


    /*
    int found = substrpos(insert,sub);
    if(found == -1) error("Hey buddy, the substring you wanted to search is not included in the whole string in global function >>substitute<< (is implemented in 'myfunctions.h'.");
       return insert.replace(insert.find(sub), sub.size(), rep);
    */

} 

//returns present working directory without basename
inline std::string get_present_directory(){
    std::string store;
    system("basename $(pwd) > XSTORE.txt"); //remind the '$'
    std::ifstream ilsa("XSTORE.txt", std::ios_base::in);
    getline(ilsa,store);
    ilsa.close();
    system("rm XSTORE.txt");
    return store;
}
//...and with basename
inline std::string p_w_d(){
    std::string store;
    system("pwd > YSTORE.txt");       //no '$' !!
    std::ifstream ilsa("YSTORE.txt", std::ios_base::in);
    getline(ilsa,store);
    ilsa.close();
    system("rm YSTORE.txt");
    return store;

}
 
//returns 'count' (position of last slash) as well as if a '/''s been found at all
 inline std::pair<int,bool> slash(const std::string& filename){
  int count = 0, dim = int(filename.length());
     bool recognize = false;
     
     //find out if there are any delimiters '/' in the filespecifier...
     for(int i = 0; i < dim; i++){
	if(filename[i] == '/'){
	   recognize = true;
	   break;   //leave loop
	}
    }

    //...if '/' found then extract everything after the last '/' 
    if(recognize == true){
	for(int i = dim-1; i >= 0; i--){
	    count++;
	    if(filename[i] == '/'){
		count--; //since '/' is encountered
		break;
	    }
	}
    }
    
    std::pair<int,bool> couple(count,recognize);
    return couple;
}




//this inline procedure returns the filename out of the total path, e.g. if the whole path is named "~/FOO/Yeah/sucks.txt" then 'suck.txt' will be extracted. Conversely, if no delimiter '/' is used then 'filename' is returned since this is already the filename
inline std::string get_filename(const std::string& filename){
   
    std::pair<int,bool> two = slash(filename);   
    int dim = int(filename.length());
    
	if(two.second == true){
	    std::string onlyName(two.first,' ');
	    for(int i = 0; i < two.first; i++){
		onlyName[i] = filename[dim-two.first+i];
	    }	
	    return onlyName;
	}
	else return filename;
}

//returns an empty string if no '/' is encountered (useful for init-function of class inputdata...DO NOT FORGET 'inline' here ;-)
inline std::string extract_filename(const std::string& filename){
    std::pair<int,bool> two = slash(filename);   
    int dim = int(filename.length());
    
	    std::string onlyName(two.first,' ');
	    for(int i = 0; i < two.first; i++){
		onlyName[i] = filename[dim-two.first+i];
	    }	
	    return onlyName;
}








/************************** TEMPLATE FCTS *******************************/
 template<class T>  T absolut( T a ){ // see introduction remark
  if(a<0) return -a;
  else return a; 
}

template<class T> double squareroot(T& s){     //Babylonian iteration of the square root
  T* x = new T [10];              //number of iterations: quad. convergence
  x[0]=1.0; 
  if (s < 0) error(" Negative reals have no real roots.");
  for(int i = 0; i<=9; i++){
    x[i+1]=0.5*(x[i]+s/(x[i]));
  }
  return x[9];
  delete[] x;
}
 
template<class T> T maxi(T* v, int n){     
for(int i=0; i<=n-1; i++){
    if (v[i]>=v[i+1])  v[i+1] = v[i]; 
    else v[i+1] = v[i+1];
    
 }
  return v[n-1];
}    

 
template<class T> T mini(T* v, int n){     //alternative: T v[] 
for(int i=0; i<=n-1; i++){
  if (v[i]<= v[i+1])  v[i+1]=v[i];   
  else v[i+1]=v[i+1];
 }
  return v[n-1];
}    

template<class T> int sign(const T& si){
  if(si > 0) return +1;
  else if (si == 0) return 0;//0.0;
  else  return -1;
} 

template<class T> const T power(const T& base, int exponent){
  T res = 1;
  T pof = base;
  while(exponent){
    if(exponent % 2) res *= pof; 
      pof *= pof;
      exponent /= 2;
  }
  return res;
}   


template<class T> void permute(T& a, T& b){
  T tmp = a;
  a = b;
  b = tmp;
}

template<class T> bool lessthan(T a, T b){ return  a<b;}

//convert s.th. to a string
 template<class T> 
   inline std::string num2str(const T& number){
   std::stringstream giveback;
   giveback << std::setprecision (PREC)  << number; 
   return giveback.str();
 }



//convert string to s.th. (int, long int, double): invoke it via str2num<double>(string). For non-rounded output use, e.g. 'std::cout.precision(17);' etc..
 template<class T> 
   inline T str2num(const std::string& s){
    return T(atof(s.c_str()));
} 


/*************************** NON-TEMPLATE FCTS  **************************/

 inline int gcd(int a, int b){
   if(b==0) return a;
   else return gcd(b,a%b);
 }
 
 inline bool whitecomma(char c){
   bool result(false); 
  switch(c){ //the following 6 cases equal the 'isspace' function of <fstream>
      case ' ':          //space
	  result = true;
	  break;
      case '\t':         //horizontal tab
	  result = true;
	break;
      case '\n':        //newline
	  result = true;
	  break;
      case '\v':        //vertical tab 
	  result = true;
	  break;
      case '\f':        //form feed
	  result =  true;
	  break;
      case '\r':        //carriage return
	  result = true;
	  break;
      case ',':        //if a comma is detected
	  result =  true;
          break;
      default: 
          result = false;
  }
  return result;
 }
 
 inline bool delimiter(char c){
   if(isspace(c)) return true; // c=whitespace (' ','\t','\n','\f','\v','\r')?
   bool result(false);
  switch(c){
    case ',':                //possible field delimiter (like whitespace)
      result = true;  
      break;
    case '"':               //text delimiter
      result =  true;
      break;
      // case '#':
      //result =  true;
      //break;
    default: 
     result = false;  
  }
  return result;
 }


 //! alias 
 inline bool is_delimiter(char c){ 
   return delimiter(c);
}

 inline bool is_colon(char c){
   if(c == ':')
     return true;
   else 
     return false;
 }
 
 
 inline bool is_whiteline(const std::string& s){
   bool george = true;
   for(size_t i = 0; i < s.size(); ++i){
     if(!delimiter(s[i])){
       george = false;
       break;            //no nead for further checking
     }
   }
   return george;
 }

 inline bool is_comment(const std::string& s){
   bool george = false;
   for(size_t i = 0; i < s.size(); ++i){
     if(s[i] == '#'){          //o.k. found comment 
       george = true;
       break;
     }
   }
   return george;
 }


 //indicate comment by '#', akin to shell scripts, perl, gnuplot, etc.
 inline bool is_comment(const char& c){
   return ( (c == '#') ? true : false );
 }

 

 inline bool is_comment_or_block(const char& c){
   if(c == '#' || c == '$')
     return true;
   else return false;
 }

 
 inline bool is_hash_dollar_or_at(const char& c){
   if(is_comment_or_block(c) || c =='@')
     return true;
   else return false;
 }


 inline std::string ignore_whitespaces(std::string s){
   int bou = int(s.size());
   int i, count = 0;
   std::string aux(bou,' ');
   for(i = 0; i < bou; i++){
      if(!whitecomma(s[i])){
	aux[count] = s[i]; 
	count++;
      }
   }
   return aux.substr(0,count);
 }

 inline bool fancychar(char c){
   if(delimiter(c) || c == '$') return true;
   else return false;
 }

inline  bool no_whitespaces(std::string s){
   int b = int(s.size());
    // bool cool = true;
   for(int i = 0; i < b; i++){
     if(whitecomma(s[i])){
	    //cool = false;
       return false;
       break;   //leave the loop
     }
   }
    return true;
 }

inline  bool only_whitespaces(std::string s){
   int b = int(s.size());
   for(int i = 0; i < b; i++){
     if(!delimiter(s[i])){
       return false;
       break;
     }
   }
   return true;
 }


inline  bool separator(char c){
   switch(c){ //without '\n'
    case ' ':
      return true;
      break;
    case '\t':
      return true;
      break;
    case '\v':
      return true;
      break;
    case '\f':
      return true;
      break;
    case '\r':
      return true;
      break;
    case ',':        //if a comma is detected
      return true;
    default: return false;
    }
 }

 inline bool specialsymb(char c){
   if(delimiter(c)) return true;
   switch(c){ 
   case '#':       // if '#' is found
     return true;
     break;   
   case '$':       //and dollar found
     return true;
     break;
   default: return false;
   }
}
 
inline bool mysigns(std::string& s){
  int b = int(s.size());
  for(int i = 0; i < b; i++){
    if((s[i]) == '$' || s[i] == '&'){
      return true;
      break;
    }
  }
  return false;
}

inline bool contains_hash_dollar(const std::string& s){
  int b = int(s.size());
  for(int i = 0; i < b; i++){
    if((s[i]) == '#' || s[i] == '$'){
      return true;
      break;
    }
  }
  return false;
}


inline bool contains_dollar(const std::string& s){
  size_t b = s.size();
  for(size_t i = 0; i < b; i++){
    if((s[i]) == '$'){
      return true;
      break;
    }
  }
  return false;
}



inline bool contains_comment(const std::string& s){
  size_t b = s.size();
  for(size_t i = 0; i < b; i++){
    if((s[i]) == '#'){
      return true;
      break;
    }
  }
  return false;
}

 
 inline bool lines_to_skip(const std::string& s){
   if(only_whitespaces(s))
     return true;
   else{
     size_t b = s.size();
     for(size_t i = 0; i < b; i++){
       if(s[i] == '#' || s[i] == '$'){
	 return true;
	 break;
       }
     }
   }
   return false;
 }


 inline std::string skip_Comment(std::string& s){
   int b = int(s.size()),  count = 0;
   for(int i = 0; i < b; i++){
     if(s[i] == '/' && s[i+1] == '/' && i != b-1){
	  break;//o.k. comment found; leave loop earlier
     }
     count++;
   }
   return s.substr(0,count);
    //return count;

 }



}//end namespace

#endif
