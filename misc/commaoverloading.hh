#ifndef COMMA_OVERLOADING_IN_C_PLUS_PLUS_HH
#define COMMA_OVERLOADING_IN_C_PLUS_PLUS_HH

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "../common/adonisassert.hh"
#include "../common/error.hh"

namespace Adonis{

  /**
   * \brief Overload commas; code adapted from [1]
   * 
   * Example:
   * \code 
         MyVec<double> v(4);
         v <<= 0.5, -2.25, 1.5, 0.75; 
   * \endcode
   *
   * Reference:
   *
   * [1] [T. Veldhuizen, "Techniques for Scientific C++", Tech. Report, Indiana University Computer Science, 2000]
   *
   */
  template<class T, class ITER>
  class CommaOverloading{
  public:
    typedef std::size_t SizeType;
    typedef T value_type;
    typedef ITER iter_type;
    
    //! \param iter iterator of random access container
    //! \param n length of container (for size check)
    //! \param count counter that counts commata (for size check)  
    CommaOverloading(ITER iter, const SizeType n, SizeType count):it_(iter),length_(n),count_(count){}

    CommaOverloading operator,(const T& tval){ //! overload comma here
      count_++;
      
      //!check if there are more entries than lenght_ are assigned (less are 
      //! notest problematic, since following entries are zero)
#ifndef NDEBUG    
      if(count_ >= length_) {
	std::cout << "count_ = "<<count_ << "   length_ = "<< length_ << std::endl;
	ADONIS_ERROR(IterationError,"Length of rac and number of assigned types disagree:\n   You wanna assign "<<count_+1<<" elements to a countainer of size "<< length_<<".");
      }
#endif
      
      *it_ = tval;
      return CommaOverloading(it_+1,length_,count_);
    }

  private:
    ITER it_;
    SizeType length_;   //! length of random access container
    SizeType count_;          //! count entries to the right of assignment 
                               //! operator (assume that count_ = 0)

    CommaOverloading(){}  //no default constructor
  };

} //end namespace 

#endif
