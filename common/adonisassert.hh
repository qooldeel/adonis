#ifndef ADONIS_ASSERT_HH

# undef	ADONIS_ASSERT_HH
# undef	adonis_assert
# undef __ASSERT_VOID_CAST

#define ADONIS_ASSERT_HH      1

#include <cstdlib>          //! needed for abort() etc
#include <features.h>       //!this includes <sys/cdefs.h> from /usr/include/sys
//======================
#include<sstream>
#include "fancymessages.hh"  //!cool error representation ;-)

//========================


#if defined __cplusplus && __GNUC_PREREQ (2,95)
# define __ASSERT_VOID_CAST static_cast<void>
#else
# define __ASSERT_VOID_CAST (void)
#endif

/** \brief Acts in a similar way like the <cassert>  
 *  void adonis_assert (int expression);
 *
 *  If an assertion fails, a message will be printed and the program will be aborted.
 *  When the macro NDEBUG is defined then all assertions will be ignored.
 *  Thus assertions can be freely used during program development and later discarded for the sake of run-time efficiency by just defining NDEBUG !
 *  Code without NDEBUG will produce lower performance since the assert condition must be evaluated.
 *
 *  If NDEBUG is defined, do nothing.
 *  If not, and EXPRESSION is zero, print an error message and abort.  
 *
 * check: '/usr/include/assert.h' to get some insight into the matter...
*/

#ifdef	NDEBUG

# define adonis_assert(expr)   (__ASSERT_VOID_CAST (0))

#else // Not NDEBUG.

//!the do-while trick where Obj directly corresponds to Adonis::FancyMessages
#define NICE_ERROR(Obj,assertion) do{Obj O; std::ostringstream oss; \
    oss << "Description: \n------------\n "<<file<<":" <<line<<": "<<function<<": ASSERTION `"<< assertion <<"' FAILED, HUNK."<<std::endl; O.nice_error(oss.str());O.error_log(oss.str()); \
  }while(0)



/**
 *  \brief Print message and abort 
 *  
 *  \param assertion Test if this assertion is fulfilled
 *  \param file Invokes macro__FILE__
 *  \param line Invokes macro __FILE__
 *  \param function If you're system supports the ISO C99 standard, macro __func__  or __PRETTY_FUNCTION__ will be executed (otherwise no fct representation is avialable).
 */

class FunnyAss: public Adonis::FancyMessages{};

inline void print_message(const char* assertion, const char* file,
			   unsigned int line, const char* function){

  NICE_ERROR(FunnyAss, assertion);

  abort();   //terminate processes
}

//! instead of 'ASSERT_FUNCTION' you can write '__func__', assumed you have the C99 standard
# define adonis_assert(expr)							\
  ((expr)								\
   ? __ASSERT_VOID_CAST (0)						\
   : print_message(__STRING(expr), __FILE__, __LINE__, ASSERT_FUNCTION))


//! Checks which function description is available (either your compiler finds
//! some nice representations or the function representation is just empty.
# if defined __cplusplus ? __GNUC_PREREQ (2, 6) : __GNUC_PREREQ (2, 4)
#   define ASSERT_FUNCTION	__PRETTY_FUNCTION__
# else
# if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
# define ASSERT_FUNCTION	__func__
# else
# define ASSERT_FUNCTION  ((__const char*) 0)
#endif
#endif

#endif // NDEBUG. 

#endif  // end of ifndef ADONIS_ASSERT_HH
