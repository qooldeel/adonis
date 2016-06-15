#ifndef TRAITS_4_CLASS_AX_AND_AXPY_HH
#define TRAITS_4_CLASS_AX_AND_AXPY_HH

namespace Adonis{

  //forward declarations
  template<class K, int NFULL, int NRED, bool B> class Essential;
  template<class K, int NFULL, int NRED, bool B> class Ax;
  template<class K, int NFULL, int NRED, bool B> class Axpy;

  template<class K, int NFULL, int NRED, bool B>
  class AxpyTraits{
  public:
    typedef Essential<K,NFULL,NRED,B> Type;
    typedef Ax<K,NFULL,NRED,B> AxType;
    typedef Axpy<K,NFULL,NRED,B> AxpyType;
  };

}

#endif
