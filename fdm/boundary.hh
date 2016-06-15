#ifndef SPECIES_BOUNDARY_VALUES_MOL_2D_HH
#define SPECIES_BOUNDARY_VALUES_MOL_2D_HH

#include "../common/adonisassert.hh"

namespace Adonis{

  //! \tparam V vector type
  //! \tparam MECH encoding for mechansim
  template<class V, int MECH>
  class Boundary2D{
  public:
    typedef typename V::value_type value_type;
    
    Boundary2D(int n = 0):v_(n){}

    void initialize(int n){
      adonis_assert(n == 3);
      if(v_.size() != 3)
	v_.resize(3);
      //makes more sense here, since we only assigne values once and allow
      //acces via []
      v_[0] = 0.;     //O
      v_[1] = 0.8;    //O2
      v_[2] = 0.2;    //O3

    }

    const value_type& operator[](int i) const {
      adonis_assert(i >= 0 && i < (int)v_.size());
      return v_[i];
    }

    value_type& operator[](int i){
#ifndef NDEBUG
      //adonis_assert(i >= 0 && i < (int)v_.size());
      if(i<0 || i >= (int)v_.size())
	ADONIS_ERROR(IndexError,"Index i = "<<i<<" out of bounds [0,"<<v_.size()-1<<"].");
#endif
      return v_[i];
    }

    const V& get_vector() const {return v_;}

    V& operator()(){
      
      v_[0] = 0.;     //O
      v_[1] = 0.8;    //O2
      v_[2] = 0.2;    //O3

      return v_;
    }

  private:
    V v_;

  };

}

#endif
