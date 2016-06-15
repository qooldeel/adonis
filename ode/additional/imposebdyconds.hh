#ifndef IMPOSE_PHYSICAL_BOUNDARY_CONDITIONS_HH
#define IMPOSE_PHYSICAL_BOUNDARY_CONDITIONS_HH

#include "../../common/globalfunctions.hh"
#include "../../fdm/gensettings.hh"
#include "../expressiontemplates/exprvec.hh"

namespace Adonis{


  template<class V, int DIM, bool ISRECTANGULAR>
  class BoundaryConditions{
  public:
    
    void init(const std::string& filename){}

    void impose(V& u) {}
  };

  //! partial specialization for rectangular 2D boundary
  template<class V>
  class BoundaryConditions<V,2,true>{
  public:
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef ExprTmpl::MyVec<T> VType; 
    typedef ExprTmpl::MyVec<DType> VDType; 
    
  private:
    

  };

}

#endif
