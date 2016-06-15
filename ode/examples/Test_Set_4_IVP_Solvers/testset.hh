#ifndef TEST_SUITE_FOR_INITIAL_VALUE_PROBLEM_SOLVERS_HH
#define TEST_SUITE_FOR_INITIAL_VALUE_PROBLEM_SOLVERS_HH

#include "../../../expressiontemplates/exprvec.hh"
#include "../../../misc/misctmps.hh"

namespace Adonis{

  /**
   * \brief stiff system
   *
   * This example can be found in 
   * 
   * [1] ERIKSSON, JOHNSON and LOGG, "Explicit time-stepping for stiff ODEs", SIAM, J. Sci. Comput., 25, 2003
   *
   * [2] LIOEN and DE SWART, "Test Set for Initial Value Problem Solvers", Tech. report MAS-R9832, CWI, Amsterdam, 1998
   */
  template<class T>
  class AkzoNobelProblem{
  public:
    typedef T value_type;
    typedef ExprTmpl::MyVec<T> VType;

    size_t dim() const {return v_.size();}
    size_t domain_dim() const {return 6;}
    

    AkzoNobelProblem(size_t dim = 0):v_(dim){}

    template<class X>
    VType& operator()(const X& u) {
      
      //for(int i = 0; i < 6; ++i) std::cout << "Input expression: u["<<i<<"] = "<< u[i] << std::endl;

      //!the F-term
      T Fterm = 3.3*(0.9/737. - u[1]),
	goodsqrtu1 = sqrt(Abs(u[1]));  //beware of negative approximation values when dealing with sqrt-roots ;) ;) ;) -- Otherwise a NaN is produced

      //!reaction rates
      T r1 = 18.7*ntimes<4>(u[0])*goodsqrtu1,
	r2 = 0.58*u[2]*u[3],
	r3 = 0.58/34.4*u[0]*u[4],
	r4 = 0.09*u[0]*ntimes<2>(u[3]),
	r5 = 0.42*ntimes<2>(u[5])*goodsqrtu1;

      //std::cout << "r1 = "<< r1 <<  std::cout << "r2 = "<< r2 << std::endl << "r3 = "<< r3 << std::endl << "r4 = "<< r4 << std::endl << "r5 = "<< r5 << std::endl;

      v_[0] = -2.*r1 + r2 - r3 -r4;
      v_[1] = -0.5*r1 - r4 - 0.5*r5 + Fterm;
      v_[2] = r1 - r2 + r3;
      v_[3] = -r2 + r3 - 2.*r4;
      v_[4] = r2 - r3 + r5;
      v_[5] = -r5;


      //std::cout << "v_ = "<< v_ << std::endl;

      return v_;

    }

    //! fake time dependent operator
    template<class V>
    VType& operator()(const T& time, const V& y){
      return (v_ = (*this).operator()(y));
    }
    


  private:
    VType v_;

  };

}

#endif
