#ifndef SIMPLY_BOUNDED_AND_EQUALITY_CONSTRAINED_PROBLEM_INTERFACE_HH
#define SIMPLY_BOUNDED_AND_EQUALITY_CONSTRAINED_PROBLEM_INTERFACE_HH

#include "../expressiontemplates/exprvec.hh"
#include "../common/typeadapter.hh"
#include "../common/adonisassert.hh"

namespace Adonis{

  /**
   * \brief General interface using the somehow deprecated (due to C++11)
   *  but nevertheless astonishing Barton-Nackman trick ;) 
   *
   * This interface is solely defined for equality-constrained, simply bounded
   * problems, i.e. optimization problems of the form \f[
   *   \min_{x \in R^n}\{f(x)| c(x) = 0 \wedge l \leq x \leq u\},
   *   \f]
   * where the expression after the '|' are to be understood componentwise.
   */
  template<class T, class C>
    class GeneralSimplyConstrainedProblem{
  protected:
    typedef ExprTmpl::MyVec<T> VType;
    typedef typename TypeAdapter<T>::Type Type;
    typedef ExprTmpl::MyVec<Type> OrdinaryTypeVType; 
   
    size_t n_,   //!number of variables in the problem
      m_;        //!number of equality constraints
    VType eq_,       //!holds the m equality constraints
      low_,      //!holds the n lower bounds for variable x
      up_;       //!holds the n upper bounds for variable x  
    VType red_;

    void initialize(){
      if(m_ != 0) 
	eq_.resize(m_);
      if(n_ != 0){
	low_.resize(n_);
	up_.resize(n_);
      }
    }

  private:
    C& ref2derived(){return static_cast<C&>(*this);}
    
    //! add some implementations here that differ in each derived class
    template<class E> 
    VType& eq(const E& x){
      return ref2derived().equalityconstraints(x);
    }

    VType& lowerbds(){
      return ref2derived().low();
    }
    
     VType& upperbds(){
      return ref2derived().up();
    }
    
    template<class E>
    T cost(const E& x){
      return ref2derived().objective(x);
    }

    //!common functionality
  public:
    size_t number_of_variables() const {return n_;}
    size_t number_of_equalities() const {return m_;}

    template<class RV>
    void set_additional_bounds(const RV& red){
      if(red.size() != 0){
	//red_ = red;
	 red_.resize(red.size());
	 for(size_t i = 0; i < red.size(); ++i)
	   red_[i] = red[i];
      }
    }

    VType& get_red() {return red_;}
  };

} //end namespace 

#endif
