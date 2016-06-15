#ifndef OUTPUT_FEATURES_WHEN_WRITING_TO_FILES_HH
#define OUTPUT_FEATURES_WHEN_WRITING_TO_FILES_HH

#include <iostream>

namespace Adonis{

  template<bool B> class WhiteLine;

  template<> class WhiteLine<true>{
  public:
    typedef std::ofstream OfStreamType;
    
    static inline void write(OfStreamType& ofs){
      ofs << std::endl;
    }

    static inline void write(){
      std::cout << std::endl;
    }
  };

  template<> class WhiteLine<false>{
  public:
    typedef std::ofstream OfStreamType;
    
    static inline void write(OfStreamType& ofs){}  //don't do anything
 
    static inline void write(){} //do nothing
  };
  

}  //namespace Adonis

#endif
