#ifndef PHYSICAL_TRANSPORT_HH
#define PHYSICAL_TRANSPORT_HH

namespace Adonis{

  /**
   * \brief the diffusion \f$ D_2(x) := D_1 + e\cdot x, \ D_1 \in \mathcal{R}_+, x \in \Omega\f$ as it can be found in [REN/POPE, Comb. Flame, vol. 147, 2006]
   */
  template<class T>
  class DiffusiveTransport{
  public:
    typedef T value_type;
    
    DiffusiveTransport(const T& D1 = T(), const T& e = T()):D1_(D1),e_(e){}
    
    T operator()(const T& x){
      return D1_ + e_*x;
    }
    
  private:
    const T& D1_;       //has no copy constructor…
    const T&  e_;       //…only const members
  };

}

#endif
