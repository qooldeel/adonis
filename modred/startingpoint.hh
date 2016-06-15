#ifndef DEFINE_STARTING_POINT_FOR_ODE_SOLUTION_HH
#define DEFINE_STARTING_POINT_FOR_ODE_SOLUTION_HH

#include <iostream>

namespace Adonis{

  template<class V>
  class StartingPoint{
  public:
    typedef typename V::value_type value_type;

    StartingPoint(size_t n = 0):n_(n),start_(n){}

    V& get_point(){
      //===== TODO: set your initial point here ===============================
       adonis_assert(n_ == 6);
      value_type onmani[] = {0.455, 0.779, 0.237, 0.363, 0.148, 0.015};
      for(size_t i = 0; i< n_ ;++i)
      	start_[i] = onmani[i];
     
      
      // start_[0] = 0.4549999999961185; // H2
      // start_[1] = 0.7780667626813738;
      // start_[2] = 0.2365973949160861;
      // start_[3] = 0.3628719728410923;
      // start_[4] = 0.1479999999996536; // H2O
      // start_[5] = 0.01593323732708175;
      

      //=======================================================================
      
      return start_;
    } 

    friend std::ostream& operator<<(std::ostream& os, const StartingPoint& sp){
      for(size_t i = 0; i < sp.n_; ++i)
	os << sp.start_[i] << " ";
      os << std::endl;
      return os;
    }

  private:
	size_t n_;
    V start_;

  };


} //end namespace

#endif
