#ifndef PARAMETERS_FOR_THE_REN_POPE_TOY_EXAMPLE_HH
#define PARAMETERS_FOR_THE_REN_POPE_TOY_EXAMPLE_HH

#include <iostream>
#include <vector>
#include "../common/adonisassert.hh"

namespace Adonis{

  /**
   * \brief The fixed parameters \f$ a, b, c, d, e, D_1, \epsilon\f$ are stored in the described order
   */
  template<class T>
  class SpecificParameters4RP{
  public:
    static inline T* case1(){
      static T a[7] = {0,1,1,2,0,0,0.01};
      return a;
    }

    static inline T* case2(){
      static T a[7] = {1,1,1,1,0,0,0.01};
      return a;
    }

    static inline T* case3(){
      static T a[7] = {0,0,1,1,3,1,0.01};
      return a;
    }
    
    static inline T* case4(){
      static T a[7] = {1,1,1,2,3,1,0.01};
      return a;
    }

  };

  
  /**
   * \brief TMP selector for the parameter set
   */
  template<class T, unsigned B>
  class ChooseFromParameterCollection{
  public:
    static inline void param(){
      ADONIS_ERROR(UndefinedCaseError,"Nothing specified for B = "<<B<<".");
    }
  };

  //!partial specialisation, the 2nd tmpl argument denotes the case
  template<class T>
  class ChooseFromParameterCollection<T,1>{
  public:
    static inline T* param(){return SpecificParameters4RP<T>::case1();}
  };

  template<class T>
  class ChooseFromParameterCollection<T,2>{
  public:
    static inline T* param(){return SpecificParameters4RP<T>::case2();}
  };

  template<class T>
  class ChooseFromParameterCollection<T,3>{
  public:
    static inline T* param(){return SpecificParameters4RP<T>::case3();}
  };
  
  
  template<class T>
  class ChooseFromParameterCollection<T,4>{
  public:
    static inline T* param(){return SpecificParameters4RP<T>::case4();}
  };

  
  /**
   *\brief Choose parameter set 
   */
  template<unsigned N, class V>
  class FixedRPParameter{
  public:
    typedef typename V::value_type value_type;

    FixedRPParameter():p_(7){
      for(unsigned i = 0; i < 7; ++i)
	p_[i] = ChooseFromParameterCollection<value_type,N>::param()[i];
    }

    const value_type& a() const {return p_[0];}
    const value_type& b() const {return p_[1];}
    const value_type& c() const {return p_[2];}
    const value_type& d() const {return p_[3];}
    const value_type& e() const {return p_[4];}
    const value_type& D1() const {return p_[5];}
    const value_type& epsilon() const {return p_[6];}

    const size_t size() const {return p_.size();}

    bool is_sim_linear() const{
      if(p_[0] == 0)
	return true;
      else 
	return false;
    }

    template<class IDX>
    const value_type& operator[](IDX i) const{
      adonis_assert(i >= 0 && i < 7);
      return p_[i];
    }
    
    friend std::ostream& operator<<(std::ostream& os, const FixedRPParameter& p){
      os << "Parameter set case " << N <<" (a,b,c,d,e,D1,eps):"<<std::endl;
      for(size_t i = 0; i < p.size(); ++i)
	os<<p[i]<< "  ";
      os<<std::endl;
      return os;
    }
    
  private:
    V p_;
  };

}

#endif
