#ifndef DEFINE_BASIC_MACHINE_TYPES_HERE_HH
#define DEFINE_BASIC_MACHINE_TYPES_HERE_HH


namespace Adonis{
                                //# of bytes(= 8 bits):
  typedef unsigned char ui8;    //sizeof(unsigned char) = 1
  typedef unsigned short ui16;  //sizeof(unsigned short) = 2
  typedef unsigned ui32;        //sizeof(unsigned) = 4
  typedef unsigned long ui64;   //sizeof(unsigned long) = 8

  typedef signed char si8;      //sizeof(signed char) = 1
  typedef signed short si16;    //sizeof(signed short) = 2
  typedef signed si32;          //sizeof(signed) = 4
  typedef signed long si64;     //sizeof(signed long) = 8

  typedef float fl32;           //sizeof(float) =  4
  typedef double fl64;          //sizeof(double) = 8
  typedef long double fl128;    //sizeof(long double) = 16
  
}

#endif
