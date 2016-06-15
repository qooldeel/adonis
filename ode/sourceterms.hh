#ifndef SOURCE_TERMS_FOR_CHEMISTRY_HH
#define SOURCE_TERMS_FOR_CHEMISTRY_HH

#include <cassert>
#include <cmath>

#include "constants.hh"

#include "../common/typeadapter.hh"
#include "../misc/misctmps.hh"
#include "rpparameter.hh" //parameter for the Ren Pope source term

#include "../linalg/linearsystemsolvers.hh"
#include "../containers/densesymmetricmatrix.hh"
#include "../common/typeadapter.hh"

namespace Adonis{

  /**
   * \brief Base class for chemical source terms.
   *NOTE: Sometimes a solver for class ODE is time dependent. In case of autonomous systems just overload <TT> operator()()</TT> with a first dummy argument representing time 
   */
  template<typename T> 
  class SourceTraits{
  public:
    typedef ExprTmpl::MyVec<T> VType;
  };


  //No constructor at all in base class that uses the Barton-Nackman trick
  template<class CHILD, class T, template<class D, class A = std::allocator<D> >
	   class V>
  class BaseSource{
  protected:
    typedef T field_type;
    typedef T value_type;
#if USE_CPPAD
    typedef CppAD::AD<T> ADType;
    typedef CppAD::ADFun<T> ADFunType;
    typedef V<ADType> VecADType; 
#endif 
    typedef V<T> VecType;
  
    typedef CHILD LeafType;
    
    typedef typename TypeAdapter<T>::BaseType time_type;

    VecType res_;
    int n_;
    


  private:
    /*inline const CHILD& Ref2Derived() const{
      return static_cast<const CHILD&>(*this);
      }*/

    inline CHILD& Ref2Derived(){           //read 'n' write permission
      return static_cast<CHILD&>(*this);
    }
    
    //use feature of derived class, namely operator()(·)
    template<class X> // typename X::value_type
    T parenthesis1D(const X& x){
      return  Ref2Derived()(x);   //invoke derived class' operator()(·)
    }

    
    
  public: 
    T& operator[](int i){
      assert((i >= 0) && (i < n_));
      return res_[i];
    }
    
    const T& operator[](int i) const{
      assert((i >= 0) && (i < n_));
      return res_[i];
    }
  
   
  
    void resize(size_t sz, const T& c = T()){
      res_.resize(sz,c);
    }

    const VecType& get_container() const{return res_;}
    VecType& get_container() {return res_;}

    const VecType& state() const{return res_;}
    VecType& state() {return res_;}

#if USE_CPPAD
    //!Common functionality: only <B>one</B> version to set up the AD without virtual functions
    template<class ADF>
    void set_up_4_cppad(ADF& q, unsigned domainDim, unsigned rangeDim){
      adonis_assert(rangeDim == 1); //only for 1D
      VecADType X(domainDim), Y(rangeDim);
      CppAD::Independent(X);
      
      // = X[0]*X[0] //.assign(X[0]*X[0])
      //parenthesis1D(X);
      //Y.assign(parenthesis1D(X));
      Y.assign(X[0]*X[0]);//does NOT work: assign(parenthesis1D(X));
      q.Dependent(X,Y);
      q.optimize();
    }
#endif

    //works
    template<class X>
    T do_something(const X& x){
      T s =0;
      s += parenthesis1D(x);
      return s;
    }
  };

  

  //============================================================================
  //========================= DERIVED CLASSES ==================================
  //============================================================================
  
  //************************** TEST SPARSE JACOBIANS **************************/
  /**
   * \brief Simple source term producing  a sparse matrix
   */
 template<class T>
  class MechWithSparseJac: public BaseSource<MechWithSparseJac<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<MechWithSparseJac<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
   MechWithSparseJac(int n = 0){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 3;}
    size_t domain_dim() const{return 3;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

   template<class X>
   VecType& operator()(const X& x){
      BaseType::res_[0] = x[0]*x[0]*x[1];
      BaseType::res_[1] = -2.5*x[2];
      BaseType::res_[2] = x[1] + x[2]*x[2];

      return BaseType::res_;
   }

 };

  /**
   * \brief Simple source term producing  another  sparse matrix, taken
   *        from CppAD homepage
   */
 template<class T>
  class MechWithSparseJac_1: public BaseSource<MechWithSparseJac_1<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<MechWithSparseJac_1<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
   MechWithSparseJac_1(int n = 0){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 3;}
    size_t domain_dim() const{return 4;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

   template<class X>
   VecType& operator()(const X& x){
      BaseType::res_[0] = x[0] + x[1];
      BaseType::res_[1] = x[2] + x[3];
      BaseType::res_[2] = x[0] + x[1] + x[2] + x[3] * x[3] / 2.;

      return BaseType::res_;
   }

 };

 template<class T>
  class MechWithSparseJac_2: public BaseSource<MechWithSparseJac_2<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<MechWithSparseJac_2<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
   MechWithSparseJac_2(int n = 0){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 3;}
    size_t domain_dim() const{return 2;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

   template<class X>
   VecType& operator()(const X& x){
      BaseType::res_[0] = x[0]*x[0];
      BaseType::res_[1] = 1.75;
      BaseType::res_[2] = x[0]*x[1];

      return BaseType::res_;
   }

 };

  template<class T>
  class MechWithSparseJac_3: public BaseSource<MechWithSparseJac_3<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<MechWithSparseJac_3<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
   MechWithSparseJac_3(int n = 0){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 5;}
    size_t domain_dim() const{return 4;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

   template<class X>
   VecType& operator()(const X& x){
      BaseType::res_[0] = x[0]*x[2] + x[3];
      BaseType::res_[1] = 1.75;
      BaseType::res_[2] = x[0] + 0.75*x[3];
      BaseType::res_[3] = x[3]*x[3];
      BaseType::res_[4] = x[2];

      return BaseType::res_;
   }

 };

  template<class T>
  class MechWithSparseJac_4: public BaseSource<MechWithSparseJac_4<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<MechWithSparseJac_4<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
   MechWithSparseJac_4(int n = 0){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 4;}
    size_t domain_dim() const{return 5;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

   template<class X>
   VecType& operator()(const X& x){
     BaseType::res_[0] = x[1]*x[2] -0.5*x[4];
     BaseType::res_[1] = 1.75;
     BaseType::res_[2] = x[0]*x[0] + 2.5*x[2]*x[2];
     BaseType::res_[3] = -(x[1]*x[1]*x[1]);
     
     
      return BaseType::res_;
   }
    
 };

  template<class T>
  class MechWithSparseJac_5: public BaseSource<MechWithSparseJac_5<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<MechWithSparseJac_5<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
   MechWithSparseJac_5(int n = 0){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 5;}
    size_t domain_dim() const{return 5;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

   template<class X>
   VecType& operator()(const X& x){
     BaseType::res_[0] = 2.45;
     BaseType::res_[1] = x[1]*x[2]*x[2];
     BaseType::res_[2] = 0;
     BaseType::res_[3] = x[1] + x[3]*x[4];
     BaseType::res_[4] = 0.;
     
     
      return BaseType::res_;
   }
    
 };


  /***************************************************************************/
  /**************************** RHS of ODEs *********************************/
  /**************************************************************************/
 /**
   * \brief implements a simple nonlinear test mechanism with fixed reaction rates \f[ A + A \leftrightarrow B \leftrightarrow C\f]
   */
 template<class T>
  class ToyMechanism: public BaseSource<ToyMechanism<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<ToyMechanism<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
   ToyMechanism(int n = 0, const value_type& k1 = 1., const value_type& k2 = 0.01, const value_type& km1 = 1.e-05, const value_type& km2 = 1.e-05):k1_(k1),k2_(k2),km1_(km1),km2_(km2){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 3;}
    size_t domain_dim() const{return 3;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

   template<class X>
   VecType& operator()(const X& x){
      BaseType::res_[0] = -k1_*Adonis::ntimes<2>(x[0]) + km1_*x[1];
      BaseType::res_[1] = k1_*Adonis::ntimes<2>(x[0]) - km1_*x[1] -k2_*x[1] + km2_*x[2];
      BaseType::res_[2] = k2_*x[1] - km2_*x[2];

      return BaseType::res_;
   }

   template<class X>
   value_type check_mass(const X& x) const{
     return x[0] + x[1] + x[2];    //should be 1.0
   }

 private:
   value_type k1_, k2_, km1_, km2_;
 };

  /**
   * \brief This is just for testing sparse Jacobians 
   */
  template<class T>
  class FunWithSparseJac5: public BaseSource<FunWithSparseJac5<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<FunWithSparseJac5<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    FunWithSparseJac5(int n = 0){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 5;}
    size_t domain_dim() const{return 5;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

   template<class X>
   VecType& operator()(const X& y){
     BaseType::res_[0] = y[1] + Adonis::ntimes<2>(y[2]);
     BaseType::res_[1] = 4*Adonis::ntimes<3>(y[1]) - y[2]*y[4];
     BaseType::res_[2] = -y[0]*y[3];
     BaseType::res_[3] = 2*y[1];
     BaseType::res_[4] = Adonis::ntimes<2>(y[1])*y[3];

      return BaseType::res_;
   }


 private:
  
  };



  /*
   * /brief 6-species H2 combustion mechanism, see, e.g.
   * V. Reinhardt et al., Proceedings of the European Combustion Meeting, 2007
   */
  template<class T>
  class H2Combustion6Spex: public BaseSource<H2Combustion6Spex<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<H2Combustion6Spex<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    H2Combustion6Spex(int n = 0){
      adonis_assert(n == 0 || n == 6); //only zero or correct dim construction
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 6;}
    size_t domain_dim() const{return 6;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }


    //! a different random access container may be used here
    template<class FV>
    VecType& operator()(const FV& z){
      //=======================================================================
      
      const size_t nor = 12;   //number of reactions involved
      
      T pk[nor]; //reaction rates 
      T	rV[nor];      //reaction velocities

      pk[0] = 2.0;
      pk[1] = 1.0;
      pk[2] = 1.0;
      pk[3] = 1000.0;
      pk[4] = 1000.0;
      pk[5] = 100.0;
      pk[6] = 216.0;
      pk[7] = 337.5;
      pk[8] = 1400.0;
      pk[9] = 10800.0;
      pk[10] = 33750.0;
      pk[11] = 0.7714;


      rV[0]  =  pk[0]*z[0];
      rV[1]  =  pk[6]*z[1]*z[1];
      rV[2]  =  pk[1]*z[2];
      rV[3]  =  pk[7]*z[3]*z[3];
      rV[4]  =  pk[2]*z[4];
      rV[5]  =  pk[8]*z[1]*z[5];
      rV[6]  =  pk[3]*z[0]*z[3];
      rV[7]  =  pk[9]*z[1]*z[5];
      rV[8]  =  pk[4]*z[2]*z[1];
      rV[9]  =  pk[10]*z[3]*z[5];
      rV[10] =  pk[5]*z[0]*z[3];
      rV[11] =  pk[11]*z[4];

      
      //the 6 species involved in the reaction:
      BaseType::res_[0] = - rV[0]  + rV[1]              //H2
      - rV[6]  + rV[7]
      - rV[10] + rV[11];

      BaseType::res_[1] = 2.0*rV[0] - 2.0*rV[1]          //H
      + rV[4] - rV[5]
      + rV[6] - rV[7]
      - rV[8] + rV[9];
    
      BaseType::res_[2] = - rV[2] + rV[3]                //O2                 
      - rV[8] + rV[9];
    
      BaseType::res_[3] = + 2.0*rV[2] - 2.0*rV[3]        //O
      - rV[6]     + rV[7]
      + rV[8]     - rV[9]
      - rV[10]    + rV[11];
    
      BaseType::res_[4] = - rV[4]  + rV[5]               //H2O
      + rV[10] - rV[11];

      BaseType::res_[5] = rV[4] - rV[5]                  //OH
      + rV[6] - rV[7]
      + rV[8] - rV[9];
    

      //=======================================================================
    
      return BaseType::res_; 
    }
  
    //! endow with a fake time operator to use it with time-dependent solvers
    //! use template TIME to make CPPAD distinguish between independent and dependent variables !!!
    template<class TIME, class FV>
     VecType& operator()(const TIME& time, const FV& z){
      //BasicDataType<value_type>::check();
      
      BaseType::res_ = (*this).operator()(z);
       //std::cout << "orig source: "<< BaseType::res_ << std::endl;
       return BaseType::res_;
     }


    T mass_conservation_for_H(){
      return 2*BaseType::res_[0]+2*BaseType::res_[4]+BaseType::res_[1]+BaseType::res_[5]; 
    }

    T  mass_conservation_for_O(){
      return 2*BaseType::res_[2]+BaseType::res_[4]+BaseType::res_[3]+BaseType::res_[5];
    }
    

    template<class W>
    VecType mass_conservation_in_one_equation(const W& r, bool st = false){
      adonis_assert(r.size() == 2);
      
      VecType A(4), B(4);
      A[0] = 1; A[1] = 2; A[2] = 1; A[3] = 2;
      B[0] = 3.0 -2*r[0] - 3*r[1];
      //std::cout<< "b = "<<B[0] << std::endl;

      solve(A,1,B,1,st);
      return B;
    }
  };

  
  /**
   *\brief The chemical source term from the Ren Pope paper.
   *  Z. Ren and S. Pope, Comb. Flame 147, 2006
   */
  template<class T>
  class RenPopeDavisSkodjie: public BaseSource<RenPopeDavisSkodjie<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<RenPopeDavisSkodjie<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    
    typedef typename TypeAdapter<value_type>::Type BasicType; 
    typedef Constant<BasicType> ConstParamType;

    //============= CHOOSE your cases here ================================
    typedef FixedRPParameter<ConstParamType::renpopeparameter,ExprTmpl::MyVec<T> > ParamType;
    //=====================================================================

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    RenPopeDavisSkodjie(int n = 0):param_(),sim_(2){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 2;}
    size_t domain_dim() const{return 2;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }

    inline void info() const {std::cout << "RenPope-Parameter case: "<< ConstParamType::renpopeparameter << std::endl;} 

    inline void parameter() const {std::cout<< param_ << std::endl;}

    //! a different random access container may be used here
    template<class FV>
    VecType& operator()(const FV& z){
      //=======================================================================
      T inv = 1./param_.epsilon(),
	f = z[0]/(1+param_.a()*z[0]),
	z1_z2 = z[1]-f;
      
      BaseType::res_[0] = param_.c()*inv*(z1_z2) - param_.d()*z[0];
      
      BaseType::res_[1] = -inv*(z1_z2) - z[0]/(ntimes<2>(1+param_.b()*z[0]));
      
      //=======================================================================
      return BaseType::res_;
    }

    //also the same type
    template<class TIME, class FV>
    VecType& operator()(const TIME& t,const FV& z){
      return (BaseType::res_ = (*this).operator()(z));
    }

    template<class Z>
    VecType& sim(const Z& z){
      sim_[0] = z[0];
      sim_[1] = z[0]/(1 + param_.a()*z[0]);
      return sim_;
    }

    template<class Z>
    T fprime(const Z& z) const {return 1/(ntimes<2>(1+param_.a()*z[0]));}

    template<class Z>
    VecType& tangent(const Z& z){
      sim_[0] = 1.;
      sim_[1] = (*this).fprime(z);
      return sim_;
    }

    template<class W> //massconservations unknown
    VecType mass_conservation_in_one_equation(const W& r, bool st = false){
      ADONIS_ERROR(DerivedError, "Mass conservations unknown");
      return r;
    }

  private:
    ParamType param_;
    VecType sim_;
};



/**
   *\brief The chemical source term from the Ren Pope paper.
   *  Z. Ren and S. Pope, Comb. Flame 147, 2006
   */
  template<class T>
  class FirstApproxRenPopeDavisSkodjie: public BaseSource<FirstApproxRenPopeDavisSkodjie<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<FirstApproxRenPopeDavisSkodjie<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    
    typedef typename TypeAdapter<value_type>::Type BasicType; 
    typedef Constant<BasicType> ConstParamType;

    //============= CHOOSE your cases here ================================
    typedef FixedRPParameter<ConstParamType::renpopeparameter,ExprTmpl::MyVec<T> > ParamType;
    //=====================================================================

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    FirstApproxRenPopeDavisSkodjie(int n = 0):param_(){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 1;}
    size_t domain_dim() const{return 1;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }


    inline void parameter() const {std::cout<< param_ << std::endl;}

    //! a different random access container may be used here
    template<class FV>
    VecType& operator()(const FV& z){
      //=======================================================================
      // The SIM is given by [z1,z1/(1+az1)]^T. Neglecting the second equation and noting that z2 = z1/(1+az1), we obtain simply -d·z1 for the 1st approximation
      
      BaseType::res_[0] =  - param_.d()*z[0];
      
      //=======================================================================
      return BaseType::res_;
    }

    template<class TIME, class FV>
    VecType& operator()(const TIME& t,const FV& z){
      return (BaseType::res_ = (*this).operator()(z));
    }

    template<class W> //massconservations unknown
    VecType mass_conservation_in_one_equation(const W& r, bool st = false){
      ADONIS_ERROR(DerivedError, "Mass conservations unknown");
      return r;
    }

  private:
    ParamType param_;
};

 /**
   * \brief 2nd order approximation
   */
 template<class T>
  class SecondApproxRenPopeDavisSkodjie: public BaseSource<SecondApproxRenPopeDavisSkodjie<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<SecondApproxRenPopeDavisSkodjie<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    

    typedef typename TypeAdapter<value_type>::Type BasicType; 
    typedef Constant<BasicType> ConstParamType;

    //============= CHOOSE your cases here ================================
   typedef FixedRPParameter<ConstParamType::renpopeparameter,ExprTmpl::MyVec<T> > ParamType;
    //=====================================================================

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    SecondApproxRenPopeDavisSkodjie(int n = 0):param_(){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 1;}
    size_t domain_dim() const{return 1;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }


    inline void parameter() const {std::cout<< param_ << std::endl;}

    //! a different random access container may be used here
    template<class FV>
    VecType& operator()(const FV& z){
      //=======================================================================
      // The SIM is given by [z1,z1/(1+az1)]^T. Neglecting the second equation and noting that z2 = z1/(1+az1), we obtain simply -d·z1 for the 1st approximation
      T fprime = 1./(ntimes<2>(1+param_.a()*z[0]));
      BaseType::res_[0] =  - param_.d()*z[0] + param_.c()/(param_.c()*fprime + 1)*(param_.d()*z[0]*fprime - z[0]/(ntimes<2>(1+param_.b()*z[0])));
      
      //=======================================================================
      return BaseType::res_;
    }

   template<class TIME, class FV>
   VecType& operator()(const TIME& t,const FV& z){
     return (BaseType::res_ = (*this).operator()(z));
   }

    template<class W> //massconservations unknown
    VecType mass_conservation_in_one_equation(const W& r, bool st = false){
      ADONIS_ERROR(DerivedError, "Mass conservations unknown");
      return r;
    }

  private:
    ParamType param_;
};




 /**
   *\brief Nonlinear chemical reaction, taken from
   *  
   * [HULL,ENRIGHT,FELLEN,SEDGEWICK, "Comparing numerical methods for ODEs", SIAM, J. Numer. Anal. vol. 9, no. 4, 1972]
   *
   * \f$ y(t = 0) = [1,0,0]^T. \f$
   */
  template<class T>
  class NonlinearChemicalReaction: public BaseSource<NonlinearChemicalReaction<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<NonlinearChemicalReaction<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    NonlinearChemicalReaction(int n = 0){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 3;}
    size_t domain_dim() const{return 3;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }



    //! a different random access container may be used here
    template<class X>
    VecType& operator()(const value_type& t, const X& y){
      //=======================================================================
     
      BaseType::res_[0] = -y[0];
      
      BaseType::res_[1] = y[0] - ntimes<2>(y[2]);

      BaseType::res_[2] = ntimes<2>(y[2]);
      //=======================================================================
      return BaseType::res_;
    }

    //!overload operator to work with autonomous solvers as well
    template<class X>
     VecType& operator()(const X& y){
      return (BaseType::res_ = (*this).operator()(static_cast<value_type>(1),y));
    }
};


  /**
   * \brief Models the orbit of a satellite in the gravity field induced by two huge bodies (here: earth and moon).
   *
   * ADVANTAGE: also works with expressions (see operator()() ;)
   *
   * Ref.: [WALTER, "Gewoehnliche Differentialgleichungen", 7. Aufl., p. 5, eq. (10)]
   *
   * NOTE: Eq. (10) models the movement of a satellite, viz. the curve in 2D space \f$ t \maps (x(t),y(t))\f$, around two fixed bodies (e.g. earth (E) and moon (M)). The coordinate system itself is <B>moving</B> as well (mind that when plotting; therefore E is centered at \f$ [0,0]^T \f$ and M at \f$ [1,0]^T \f$.)!    
   */
 template<class T>
  class SatelliteOrbit: public BaseSource<SatelliteOrbit<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<SatelliteOrbit<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
   SatelliteOrbit(int n = 0):mu_(0.01213), muprime_(1-mu_){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 4;}
    size_t domain_dim() const{return 4;}

    
    //!returns the rhs component of index i evaluated at t and c 
    template<class C>
    T& operator()(int i, const T& t, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const T& t, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }



    //! a different random access container may be used here
    template<class X>
    VecType& operator()(const value_type& t, const X& y){
      //=======================================================================
      value_type d1 = ntimes<3>(sqrt(ntimes<2>(y[0]+mu_)+ ntimes<2>(y[2]))),
	d2 = ntimes<3>(sqrt(ntimes<2>(y[0]-muprime_)+ ntimes<2>(y[2])));
      
      BaseType::res_[0] = y[1];
      
      BaseType::res_[1] = y[0] + 2*y[3] - muprime_*(y[0] + mu_)/d1 - mu_*(y[0] - muprime_)/d2;

      BaseType::res_[2] = y[3];

      BaseType::res_[3] = y[2] - 2*y[1] - muprime_*y[2]/d1 - mu_*y[2]/d2;

      //=======================================================================
      return BaseType::res_;
    }

    
   template<class X>
   VecType& operator()(const X& y){
     return (BaseType::res_ = (*this).operator()(static_cast<value_type>(1),y));
   }

 private:
   value_type mu_,  //! \f$ \mu \approx \frac{m_{\textrm{moon}}}{m_{\textrm{earth}}}\f$
     muprime_;      //!\f$ \mu' := 1 - \mu\f$
};




  /**
   *\brief Five body problem: the motion of the 5 outer planets about the sun. Each component has 3 spacial coordinates, where the 15 coordinates satisfy a second order ODE, cf. [HULL,ENRIGHT,FELLEN,SEDGEWICK, "Comparing numerical methods for ODEs", SIAM, J. Numer. Anal. vol. 9, no. 4, pp. 619--620, 1972]. 
   *
   * Thus, after transforming the system to first order system, we end up with 30 components.
   */
  template<class T>
  class FiveBodyProblem: public BaseSource<FiveBodyProblem<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<FiveBodyProblem<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    
    

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    FiveBodyProblem(int n = 0):k2_(2.95912208286), //gravitational constant
			       m0_(1.00000597682) //mass of sun and the 4 inner planets
    {                                   //masses of ...
      m_[0] = .000954786104043;        //Jupiter
      m_[1] = .000285583733151;        //Saturn
      m_[2] = .0000437273164546;       //Uranus
      m_[3] = .0000517759138449;       //Neptune
      m_[4] = .00000277777777778;      //Pluto

      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 30;}
    size_t domain_dim() const{return 30;}

    
    //!returns the rhs component of index i evaluated at t and c 
    template<class C>
    T& operator()(int i, const T& t, const C& c){
      return operator()(t,c)[i];  
    }

    template<class C>
    const T& operator()(int i, const T& t, const C& c) const{
      return asLeaf()(t, c)[i];  
    }



    //! a different random access container may be used here
    template<class X>
    VecType& operator()(const value_type& t, const X& y){
      //=======================================================================
     
      unsigned idx = 0;
      for(unsigned j = 0; j < 5; ++j){          //# of species 
	for(unsigned i = 0; i < 3; ++i){        //each species has 3 coordinates
	  
	  BaseType::res_[idx] = y[odd_index(i,j)];  //y'_ij
	  BaseType::res_[idx+1] = k2_*(first_term(i,j,y) + second_term(i,j,y));  
	  //std::cout << "idx = "<<idx << std::endl;
	  idx += 2;
	}
      
      }

      //=======================================================================
      return BaseType::res_;
    }

    //!overload parenthesis operator for autonomous solvers
    //!Introduce a dummy argument in the function body 
    template<class X>
    VecType& operator()(const X& y){
      return (BaseType::res_ = (*this).operator()(static_cast<value_type>(1),y));
    }

  private:
    value_type k2_,
      m0_;
    value_type m_[5];

    //these are functions implementing Problem class C5, p.619 of the above reference [HULL et al.]

    //! returns index for species \f$ y_{ij} \f$ within <TT>BaseType::res_</TT>
    //! i = coordinate index, j = species index
    template<class I>
    inline I index(I i, I j) const{
      return i*2+6*j;  
    }

    //!returns the index of the \f$ y'_{ij}\f$
    template<class I>  
    inline I odd_index(I i, I j){
      return index(i,j)+1;      
    }

    template<class I, class X>
    inline value_type r_j(I j, const X& y){  //j is fixed 
      value_type r = 0.;
      for(I i = 0; i < 3; ++i){
	r += ntimes<2>(y[index(i,j)]);
      }
      return sqrt(r);
    }

    template<class I, class X>
    inline value_type d_kj(I k, I j, const X& y){ //k and j are fixed
      value_type d = 0.;
      for(I i = 0; i < 3; ++i){
	d += ntimes<2>(y[index(i,k)] - y[index(i,j)]);
      }
      return sqrt(d);
    }

    template<class I, class X>
    inline value_type first_term(I i, I j, const X& y){
      return -(m0_ + m_[j])*y[index(i,j)]/(ntimes<3>(r_j(j,y)));
    }

    template<class I, class X>
    inline value_type second_term(I i, I j,const X& y){
      value_type st = 0.;
      for(I k = 0; k < 5; ++k){
	if(k!=j){
	  st += ( m_[k]*( (y[index(i,k)] - y[index(i,j)]))/(ntimes<3>(d_kj(j,k,y))) - y[index(i,k)]/(ntimes<3>(r_j(k,y))) );
	}
      }
      return st;
    }

};


  /**
   * \brief Example taken from [KINCAID/CHENEY, "Numerical Analysis", § 8]
   * Intial values are given by \f$ [x(1),y(1)]^T = [3, 1]^T. \f$ 
   * Integration is done on the interval \f$ -2 \leq t \leq 1\f$ with stepsize \f$ |h| = 0.1, \f$ i.e. start with \f$ h = -0.1\f$ as in the algo of this example and use \f$ t_0 = 1, t_f = -2.\f$
    */
  template<class T>
  class KincaidCheney_Ex3: public BaseSource<KincaidCheney_Ex3<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<KincaidCheney_Ex3<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
   
    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    KincaidCheney_Ex3(int n = 0){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 2;}
    size_t domain_dim() const{return 2;}

    
    //!returns the rhs component of index i evaluated at t and c 
    template<class C>
    T& operator()(int i, const T& t, const C& c){
      return operator()(t,c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const T& t, const C& c) const{
      return asLeaf()(t,c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }


    //! a different random access container may be used here
    template<class X>
    VecType& operator()(const T& t, const X& x){
      //=======================================================================
     
      
      BaseType::res_[0] = x[0] + ntimes<2>(x[1]) - ntimes<3>(t);
      
      BaseType::res_[1] = x[1] + ntimes<3>(x[0]) + cos(t);
      
      //=======================================================================
      return BaseType::res_;
    }
  };


template<class T>
  class MotivatingStiffExample: public BaseSource<MotivatingStiffExample<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<MotivatingStiffExample<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
   
    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    MotivatingStiffExample(int n = 0){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 1;}
    size_t domain_dim() const{return 1;}

    
   
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }


    //! a different random access container may be used here
    template<class X>
    VecType& operator()(const X& x){
      //=======================================================================
     
      
      BaseType::res_[0] = -15*x[0];
      
      //=======================================================================
      return BaseType::res_;
    }

  //exact solution known
  T exact(const T& time){
    return std::exp(-15*time);
  }
  
  
};

  /**
   * \brief Cost function \f$ f: \mathbb{R}^n \longrightarrow \mathbb{R}\f$
   */
   template<class T>
  class DennisSchnabelExa6_3_1: public BaseSource<DennisSchnabelExa6_3_1<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<DennisSchnabelExa6_3_1<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    
     typedef DenseSymmetricMatrix<T> SymmType;

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
     DennisSchnabelExa6_3_1(int n = 0):grad_(2), hess_(2){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 1;}
    size_t domain_dim() const{return 2;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }



    //! a different random access container may be used here
    template<class X>
    typename VecType::value_type operator()(const X& y){
      //=======================================================================
     
      BaseType::res_[0] = ntimes<4>(y[0]) +  ntimes<2>(y[0]) +  ntimes<2>(y[1]);
      //=======================================================================
      return BaseType::res_[0];
    }

     //hand-coded gradient
     template<class X>
     VecType& gradient(const X& x){
       //=======================================================================
       grad_[0] = 4*ntimes<3>(x[0]) + 2*x[0];
       grad_[1] = 2*x[1];
       //=======================================================================
       return grad_;
     }
     
     template<class X>
     SymmType& hessian(const X& x){
       //=======================================================================
       hess_[0] = 12*ntimes<2>(x[0]) + 2;
       hess_[1] = 2;
       hess_[2] = 0;
       //=======================================================================
       return hess_;
     }

   private:
     VecType grad_;
     SymmType hess_;
};

 /**
   * \brief The famous Rosenbrock's function which has its <B> global </B> minimum at \f$ (x_1,x_2)^T = (1,1)^T\f$
   *
   * Ref.: [FLETCHER, "Practical Methods of Optimization", 2nd ed., § 1, p. 7/8]
   */
  template<class T>
  class Rosenbrock: public BaseSource<Rosenbrock<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<Rosenbrock<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    typedef DenseSymmetricMatrix<T> SymmType;
   

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    Rosenbrock(int n = 0):grad_(2), completeHess_(4),hess_(3){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
      
      //well not what I wanted
      //BaseType::set_up_4_cppad(q_,(*this).domain_dim(),(*this).dim());
     }

    size_t dim() const {return 1;}
    size_t domain_dim() const{return 2;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    


    //! a different random access container may be used here
    template<class X> //typename VecType::value_type
     T operator()(const X& x){ //typename X::value_type 
      //=======================================================================
      BaseType::res_[0] = 100*ntimes<2>(x[1] - ntimes<2>(x[0])) + ntimes<2>(1- x[0]);
      //=======================================================================
      return BaseType::res_[0];
    }

     //hand-coded gradient
     template<class X>
     VecType& gradient(const X& x){
       //=======================================================================
       grad_[0] = -400*x[0]*(x[1] - ntimes<2>(x[0])) - 2*(1 - x[0]);
       grad_[1] = 200*(x[1] - ntimes<2>(x[0]));
       
       //UseAD::evaluate(diffOper_,x);
       // UseAD::D1(diffOper_,grad_);   //get operator

       //=======================================================================
       return grad_;
     }
     
    template<class X>
    SymmType& hessian(const X& x){
      //=======================================================================
      hess_[0] = 1200*ntimes<2>(x[0]) - 400*x[1] + 2;
      hess_[1] = 200;
      hess_[2] = -400*x[0];
      
      //UseAD::evaluate(diffOper_,x,0);
      //UseAD::D2(diffOper_,hess_);    //get operator
      
      //=======================================================================
      return hess_;
    }

    //store complete Hessian -- store row-wise
    template<class X>
    VecType& complete_hessian(const X& x){
      completeHess_[0] = 1200*ntimes<2>(x[0]) - 400*x[1] + 2;
      completeHess_[1] = -400*x[0];
      completeHess_[2] = completeHess_[1];
      completeHess_[3] = 200;

      return completeHess_;
    }
    
    
  private:
    VecType grad_,
      completeHess_;
    SymmType hess_;
    //ADFunType q_;
  };
 

  /**
   * \brief Exercise 2.2., p. 40, taken from Fletcher's book.
   * 
   * A local minimizer is \f$ x^* ~ (0.6959, -1.3479)^T \f$
   *
   * CAUTION: the Hessian is in general <I>not</I> positive definite!
   */
  template<class T>
  class FEx2_2: public BaseSource<FEx2_2<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<FEx2_2<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    
     typedef DenseSymmetricMatrix<T> SymmType;

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
     FEx2_2(int n = 0):grad_(2), hess_(2){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 1;}
    size_t domain_dim() const{return 2;}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }



    //! a different random access container may be used here
    template<class X>
    typename VecType::value_type operator()(const X& x){
      //=======================================================================
     
      BaseType::res_[0] = ntimes<4>(x[0]) + x[0]*x[1] + ntimes<2>(1 + x[1]);
      //=======================================================================
      return BaseType::res_[0];
    }

     //hand-coded gradient
     template<class X>
     VecType& gradient(const X& x){
       //=======================================================================
       grad_[0] = 4*ntimes<3>(x[0]) + x[1];
       grad_[1] = x[0] + 2*(1 + x[1]);
       //=======================================================================
       return grad_;
     }
     
     template<class X>
     SymmType& hessian(const X& x){
       //=======================================================================
       hess_[0] = 12*x[0];
       hess_[1] = 2;
       hess_[2] = 1;
       //=======================================================================
       return hess_;
     }

   private:
     VecType grad_;
     SymmType hess_;
};
 

  /**
   * \brief This pathological example can be found in [2, exercise 4.3, p. 98].
   * Apparently, this seems to be a multidimensional variant of the Rosenbrock 
   * function, see also [2, exercise 7.1, p. 191]. 
   * The solution is \f$ [1,1,1,...,1]^T \f$
   *
   * Troublesome starting value: \f$ [-1.2,1,-1.2,1,...,-1.2,1]^T\f$, cf. [1,p.176]
   * NOTE: This is only defined for even \f$ n\f$ !!
   *
   *  References:
   * [1] [DIXON & MILLS, "Effect of rounding errors on the variable metric method", J. Optimization Theory and Applications, Vol. 80(1), 1994, p. 175-179]
   *
   * [2] NOCEDAL & WRIGHT, "Numerical Optimization",2nd ed., Springer, 2006]
   */
  template<class T>
  class ExtendedRosenbrock: public BaseSource<ExtendedRosenbrock<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<ExtendedRosenbrock<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
 

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}
    
    //==========================================================
    enum{
      dodim = 10 //10, 50        //TODO: change domain dimension here
    };
    //==========================================================

    ExtendedRosenbrock(size_t n = 0):val_(T()){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 1;}
    size_t domain_dim() const{return static_cast<size_t>(dodim);}

     //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

    //! a different random access container may be used here
    template<class X> //typename VecType::value_type
    T operator()(const X& x, const T& alpha = 10.){//e.g. alpha = 1, 10, 100
      val_ = T();  //reset!!
      //std::cout << "n = "<< (*this).domain_dim() << std::endl;
      size_t up = (*this).domain_dim()/2;  //! sum up till n/2 is reached
      // std::cout << "up = "<< up << std::endl;
      
      for(size_t i = 0; i < up; ++i){
	val_ += ( ntimes<2>((1 - x[2*i])) + alpha*ntimes<2>((x[2*i+1] - ntimes<2>(x[2*i]))) ) ;
      }
      //=======================================================================
      BaseType::res_[0] = val_;
      //=======================================================================
      return BaseType::res_[0];
    }

  private:
    T val_;
  };



  /**
   * \brief Hodgkin-Huxley model for nerve conduction of <I> Loligo vulgaris</I>.
   * It consists of 4 ODEs with highly nonlinear terms
   *
   * The model has been taken from [1, p.322/23]
   *
   * Firing can be found, e.g., in [2, pp. 239]  
   *
   * References:
   *
   *   [1] [EDELSTEIN-KESHET, "Mathematical Models in Biology", SIAM Classics in Applied Mathematics, 2005]
   *
   *   [2] [MURRAY, "Mathematical Biology", Vol I, 3rd ed., Springer 2002]
   */
  template<class T>
  class HodgkinHuxley: public BaseSource<HodgkinHuxley<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<HodgkinHuxley<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
   
    typedef typename BaseType::time_type time_type;

    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    HodgkinHuxley(int n = 0, const T& fire = 10., const T& period = 100):C_(1.),gNa_(120.),gK_(36.),gL_(0.3),vNa_(-115.),vK_(12.),vL_(-10.5989),firing_(fire),period_(period){
      //parameter();
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return 4;}
    size_t domain_dim() const{return 4;}

    
    //!returns the rhs component of index i evaluated at t and c 
    template<class C>
    T& operator()(int i, const T& t, const C& c){
      return operator()(t,c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const T& t, const C& c) const{
      return asLeaf()(t,c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }


    //! a different random access container may be used here
    template<class X>
    VecType& operator()(const T& t, const X& x){
      //=======================================================================
      //!d/dt v
      BaseType::res_[0] = -1./C_*(gNa_*ntimes<3>(x[2])*x[3]*(x[0]-vNa_) + gK_*ntimes<4>(x[1])*(x[0] - vK_) + gL_*(x[0]-vL_) + input_current(t) ) ;
      //! d/dt n
      BaseType::res_[1] = alpha_n(x[0])*(1. - x[1]) - beta_n(x[0])*x[1];
      //! d/dt m
      BaseType::res_[2] = alpha_m(x[0])*(1. - x[2]) - beta_m(x[0])*x[2];
      //! d/dt h
      BaseType::res_[3] = alpha_h(x[0])*(1. - x[3]) - beta_h(x[0])*x[3];
      
      //=======================================================================
      return BaseType::res_;
    }

    //!autonomous case
    template<class X>
    VecType& operator()(const X& y){
      return (BaseType::res_ = (*this).operator()(static_cast<value_type>(1),y));
    }

  private:
    //!the constants
    T C_, gNa_, gK_, gL_, vNa_, vK_, vL_,
      firing_, period_; 
    
    T alpha_m(const T& v){ return 0.1*(v+25)*1./(exp((v+25.)/10.) - 1.); }
    T alpha_h(const T& v){ return 0.07*exp(v/20.); }
    T alpha_n(const T& v){ return 0.01*(v+10.)*1./(exp((v+10.)/10.) - 1.); }
    T beta_m(const T& v) { return 4.*exp(v/18.); }
    T beta_h(const T& v) { return 1./(exp((v+30.)/10.) + 1.); }
    T beta_n(const T& v) { return 0.125*exp(v/80.); } 

    T input_current(const T& time){ 
      return ((time < period_) ? firing_ : 0.);
    }
  };
  


  // template<class T>
  // class HodgkinHuxley: public BaseSource<HodgkinHuxley<T>,T,ExprTmpl::MyVec>{   
  // public:
  //   typedef BaseSource<HodgkinHuxley<T>,T,ExprTmpl::MyVec> BaseType;
  //   typedef typename BaseType::field_type field_type;
  //   typedef typename BaseType::value_type value_type;
  //   typedef typename BaseType::VecType VecType;
  //   typedef typename BaseType::LeafType LeafType;
   
  //   typedef typename BaseType::time_type time_type;

  //   //!need to be defined here since the operator depend on a template parameter
  //   LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
  //   const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


  //   //constructor 
  //   HodgkinHuxley(int n = 0):Cm(1.0),GNa(120.0), GK(36.0), Gm(0.3),
  // 			     ENa(115.0), EK(-12.0),Em(10.613){
  //     BaseType::res_.resize(n);
  //     BaseType::n_ = n;
  //   }

  //   size_t dim() const {return 4;}
  //   size_t domain_dim() const{return 4;}
  

  //      template<class C>
  //   const T& operator()(int i, const T& t, const C& c) const{
  //     return asLeaf()(t,c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
  //   }


  //   //! a different random access container may be used here
  //   template<class X>
  //   VecType& operator()(const T& t, const X& x){
  //     //=======================================================================
  //     T V=x[0],
  // 	m=x[1],
  // 	h=x[2],
  // 	n=x[3];

  //     //!d/dt v
  //     BaseType::res_[0] = (Isource(t) + GNa*m*m*m*h*(ENa-V) + GK*n*n*n*n*(EK-V) + Gm*(Em-V))/Cm;
  //     //! d/dt n
  //     BaseType::res_[1] = alpham(V)*(1-m)-betam(V)*m;
  //     //! d/dt m
  //     BaseType::res_[2] = alphah(V)*(1-h)-betah(V)*h;
  //     //! d/dt h
  //     BaseType::res_[3] = alphan(V)*(1-n)-betan(V)*n;
      
  //     //=======================================================================
  //     return BaseType::res_;
  //   }

  //   //!autonomous case
  //   template<class X>
  //   VecType& operator()(const X& y){
  //     return (BaseType::res_ = (*this).operator()(static_cast<value_type>(1),y));
  //   }


  // private:
  //   T Cm;
  //   T GNa, GK, Gm;
  //   T ENa, EK, Em;
    
  //   T alphan (const T& V) const
  //   {
  //     return (10-V)/(100.0*(exp((10-V)/10)-1));
  //   }
    
  //   T betan (const T& V) const
  //   {
  //     return 0.125*exp(-V/80);
  //   }
    
  //   T alpham (const T& V) const
  //   {
  //     return (25-V)/(10.0*(exp((25-V)/10)-1));
  //   }
    
  //   T betam (const T& V) const
  //   {
  //     return 4*exp(-V/18);
  //   }
    
  //   T alphah (const T& V) const
  //   {
  //     return 0.07*exp(-V/20);
  //   }

  //   T betah (const T& V) const
  //   {
  //     return 1.0/(exp((30-V)/10)+1);
  //   }
    
  //   T Isource (const T& t) const
  //   {
  //     if (t<100)
  // 	return 10.0;
  //     else
  // 	return 0.0;
  //   }
    
  // };


  /**
   * \brief Implementation of \f$u_t = D(u_{xx}+u_{yy})\f$
   */
  template<class T>
  class TwoDParabolic: public BaseSource<TwoDParabolic<T>,T,ExprTmpl::MyVec>{   
  public:
    typedef BaseSource<TwoDParabolic<T>,T,ExprTmpl::MyVec> BaseType;
    typedef typename BaseType::field_type field_type;
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::VecType VecType;
    typedef typename BaseType::LeafType LeafType;
    
    enum{
      NX = 15,
      NY = 10
    };


    //!need to be defined here since the operator depend on a template parameter
    LeafType& asLeaf(){return static_cast<LeafType&>(*this);}
    const LeafType& asLeaf() const{return static_cast<const LeafType&>(*this);}


    //constructor 
    TwoDParabolic(int n = 0):npt_(NX*NY), hx_(1./(NX-1)), hy_(0.1/(NY-1)){
      BaseType::res_.resize(n);
      BaseType::n_ = n;
    }

    size_t dim() const {return BaseType::res_.size();}
    size_t domain_dim() const{return BaseType::res_.size();}

    
    //!returns the rhs component of index i evaluated at c 
    template<class C>
    T& operator()(int i, const C& c){
      return operator()(c)[i];  //or: return (*this)(c)[i]
    }

    template<class C>
    const T& operator()(int i, const C& c) const{
      return asLeaf()(c)[i];  //or: return (*this)(c)[i], or:operator()(c)[i]
    }
    

    template<class X>
    VecType& operator()(const X& var){
      u_ = var; //hard copy

      //bdy treatment
      for(int i = 0; i < NX; ++i){
	u_[OFF(i,0)] = Dirichlet_bdy(i,0);       //down
	u_[OFF(i,NY-1)] = Dirichlet_bdy(i,NY-1); //up
      }
      for(int j = 0; j < NY-1; ++j){
	u_[OFF(0,j)] = Dirichlet_bdy(0,j);       //left
	u_[OFF(NX-1,j)] = Dirichlet_bdy(NX-1,j); //right
      }
      
      //interior 
      for(int i = 1; i < NX-1; ++i){
	for(int j = 1; j < NY-1; ++j){
	  BaseType::res_[OFF(i,j)] = -(D_coeff(i,j)*( (u_[OFF(i+1,j)] - 2*u_[OFF(i,j)] + u_[OFF(i-1,j)])/(hx_*hx_) + (u_[OFF(i,j+1)] - 2*u_[OFF(i,j)] + u_[OFF(i,j-1)])/(hy_*hy_) ) );
	}
      }
      
      return BaseType::res_;
    }

   
  private:
    int npt_;
    T hx_, hy_;
    VecType u_;

    int OFF(int i, int j) const{
      return (i + NX*j);
    }

    T Dirichlet_bdy(int i, int j){
      return 0;
    }
  
    T D_coeff(int i, int j){
      return static_cast<T>((i+1.5*j)/3.75);
    }


  };






} //end namespace 

#endif 
