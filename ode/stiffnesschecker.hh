#ifndef DETERMINE_WHETHER_A_NONLINEAR_PROBLEM_IS_STIFF_HH
#define DETERMINE_WHETHER_A_NONLINEAR_PROBLEM_IS_STIFF_HH

#include <limits>

#include "../common/globalfunctions.hh"
#include "../expressiontemplates/expressionvectortraits.hh"
#if USE_CPPAD
#include "../derivatives/jacobianofsource.hh"
#else
#include "../noderivativemethods/finitedifferences.hh"
#endif
#include "../linalg/linearsystemsolvers.hh"

#include "sourceterms.hh"

namespace Adonis{
  


  /**
   * \brief Formulae for computing the stiffness rate
   */
  template<class EIGVAL, class BT, unsigned I> 
  class StiffnessRate{
  public:
    static inline BT  computation(const EIGVAL& eig, const BT& a, const BT& b,const BT& prefac = 1.e+03){
      ADONIS_ERROR(UndefinedCaseError, "Case '"<<I<<"' hasn't been defined so far, pal.");
      return 0;
    }

  };

  /**
   * \brief The larger this ratio becomes (compared to 1) the stiffer the ode system gets. 
   * 
   * Implements the quotient formula \f[ \sigma = \frac{\max\{|Re(\lambda_i)|\}}{\min\{|Re(\lambda_i)|\}} \]
   * 
   *  If \f$ \sigma > 1.e-08\f$ one can say, to some extend, that the system is stiff (at the given evaluation point of the Jacobian), i.e. the larger \f$ \sigma\f$ gets, the stiffer the underlying system (at least for some \f$ Re(\lambda_i)< 0.\f$  
   */
  template<class EIGVAL, class BT>
  class StiffnessRate<EIGVAL,BT,1>{
  public:
    static inline BT computation(const EIGVAL& eig, const BT& a, const BT& b, const BT& prefac = 1.e+03){
      BT re_min = Abs(eig[0].real()), re_max = Abs(eig[0].real());
      for(unsigned i = 1; i < eig.size(); ++i){
	//if(eig[i].real() < 0){
	  re_min = std::min(re_min,Abs(eig[i].real()));
	  re_max = std::max(re_max,Abs(eig[i].real()));
	  //}
      }
      std::cout << "max{|Re(lambda_i)|} = "<< re_max << "  min{|Re(lambda_i)|} = "<< re_min << std::endl; 
     
      return ( re_max/((re_min < prefac*std::numeric_limits<BT>::epsilon()) ? 10.*std::numeric_limits<BT>::epsilon() : re_min) ) ;
    }
  };
  
  /**
   * \brief For this formula, cf. [ASCHER/PETZOLD, p. 51]: the smaller the output gets compared to -1, the stiffer the ODE system gets.
   *
   * This formula is probably less reliable than the quotient formula above.
   */
  template<class EIGVAL, class BT>
  class StiffnessRate<EIGVAL,BT,2>{
  public:
    static inline BT computation(const EIGVAL& eig, const BT& a, const BT& b,const BT& prefac = 1.e+03){
      BT re_min = eig[0].real();
      bool fac = true;
      for(unsigned i = 1; i < eig.size(); ++i){
		re_min = std::min(re_min,eig[i].real());
		//!if it turns out that all \f$Re(\lambda_i) \f$ are approximately equal, then give back -1 (no stiffness detected)
		fac *= is_equal(eig[i-1].real(),eig[i].real());
      }
      std::cout << "min{Re(lambda_i)} = " << re_min << std::endl;
      std::cout << "are all Re(lambda_i) approx. equal: " << fac << std::endl;
      
      return  (fac ? -1 : Abs(b-a)*re_min);
    }
  };
  
  

  /**
   * \brief Construction of Jacobian with or without time-dependency
   */
  template<class T, class JO, bool TD> class WithTimeDependency;

  template<class T, class JO>
  class WithTimeDependency<T,JO,true>{
  public:
    static inline void set_up(JO& Df, unsigned dd, unsigned rd, const T& time){ 
      Df.set(time,dd,rd);
    }
  };
  
  
  template<class T, class JO>
  class WithTimeDependency<T,JO,false>{
  public:
    static inline void set_up(JO& Df, unsigned dd, unsigned rd, const T& time = T()){ 
      Df.set(dd,rd);
    }
  };

  

  /**
   * \brief <B> Definition (Stiffness -- due to ref. [1]): </B> System \f[ y'(t) = f(t,y(t)), \ t \in [a,b], \quad y(a) = y_0 \f] is <B> stiff </B> iff the eigenvalues of the local Jacobian \f$ \partial_yf(t,y(t))\f$ obey \f[ \sigma(t) := \f[ |b-a|\cdot \min_j \operatorname{Re}(\lambda_j) \ll -1. \f]
   *
   * NOTE: Since the computation of eigenvalues is a rather tedious matter, you shouldn't check it to often. Normally, if reaction rates, for instance, are fixed (over a longer period of time) it suffices to measure stiffness at the beginning.
   * REFERENCES:
   * [1] [ASHER,PETZOLD, <I> Computer methods for O.D.E.s and D.A.E.s</I>, ch. 3, p. 51]
   
   *\code
    H2Combustion6Spex<double> rhs(6);
    StiffnessDetector<double,H2Combustion6Spex> IPS(rhs);
    double sta[] = {0.2325000000003616, 0.1956183381186958, 0.1206907567545383, 0.08473682461055236, 0.6655000000002098, 0.008381661880161094};
    MyVec<double> Eval(sta,sta+6);
    std::cout << "stiff-rate = "<<IPS.check(0., 5., Eval) << std::endl;
   *\endcode
   
   * How to determine stiffness ?
   * After the above code snippet has been applied, we define a system \f$ y'(t) = f(t,y(t))\f$ to  be nonstiff iff \f$ \sigma(t) \geq -1\f$, mildly stiff iff \f$ -1000 \leq \sigma(t) < -1  \f$ and heavily stiff \f$\forall \sigma(t) < -1000. \f$
   */
  template<class T, template<class D> class FUN, bool WTIME = false, unsigned FORMULATYPE = 2>
  class StiffnessDetector{
  public:
    typedef typename ExprVecTraits<T>::ExprVecType ExprVecType;
    typedef FUN<T> FunType;
#if USE_CPPAD
    typedef JacS<T,FUN,ExprTmpl::MyVec> JacobianType;
#else
    typedef FiniteDifference<FunType> FDType;
#endif
    typedef typename TypeTraits<T>::BaseType BaseType;
    typedef typename ExprVecTraits<std::complex<T> >::ExprVecType ComplexExprVecType;

    StiffnessDetector(FunType& f):f_(f){}
    
    
    
    template<class X>
    typename X::value_type check(const BaseType& a, const BaseType& b, const X& y, const typename X::value_type& time = typename X::value_type(), const BaseType& prefac = 1.e+03){
     
      ExprVecType jac;
#if USE_CPPAD
      JacobianType Df;
      //Df.set(f_.domain_dim(),f_.dim());
      WithTimeDependency<T,JacobianType,WTIME>::set_up(Df,f_.domain_dim(),f_.dim(),time);
      //1.) compute Jacobian Df(y) at y
      Df.jacobian(y);
      jac = Df.get_jacobian();
#else
      FDType FDiff(f_);
      jac= FDiff.jacobian(y);
#endif
      

      //2.) compute eigenvalues of Df(y)
      int n = static_cast<int>(f_.dim());
      ComplexExprVecType EVals(n);
      eig_general(jac,n,EVals);
      
      std::cout << "Eigen values = "<< EVals << std::endl;
      
      //3.) find eigenvalue with smallest real part and use formula of your choice
      return StiffnessRate<ComplexExprVecType,BaseType,FORMULATYPE>::computation(EVals,a,b,prefac);
     
    }

  
  private:
    FunType& f_;  //reference
  };

}

#endif
