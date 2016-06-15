#ifndef COMPUTE_JACOBIAN_OF_SOURCE_VIA_CPPAD_HH
#define COMPUTE_JACOBIAN_OF_SOURCE_VIA_CPPAD_HH

#if USE_CPPAD
#include <cppad/cppad.hpp> //include the CPPAD stuff
#endif

#include "../common/adonisassert.hh"
#include "../common/error.hh"
#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../common/fancymessages.hh"
#include "../misc/useful.hh" 
#include "../misc/misctmps.hh"
#include "../expressiontemplates/exprvec.hh"



namespace Adonis{

#if USE_CPPAD
  /**
   * \brief Get Jacobian of a source term of the form 
   * <TT> template<typename T> class Fun </TT>, i.e.
   * we expect the functor to be templatized only by the precision type and 
   * nothing else
   *
   *Sometimes it might be beneficial when the dimension of the domain and range space, respectively, need not to be known at compile time.
   *
   * General remark: The following does not work with expression (e.g. VExpr, LOR, ...) since these are solely built on iterators. That means that even when endowed with a size function you cannot resize them, as needed by CPPAD. This would leed to an errror like  usr/local/include/cppad/check_simple_vector.hpp:130: error: ‘class Adonis::ExprTmpl::VExpr<CppAD::AD<double>, Dune::FieldIterator<const Dune::FieldVector<CppAD::AD<double>, 4>, const CppAD::AD<double> > >’ has no member named ‘resize’
   * This holds for Dune's Field types as well (which are constant (compile time) containers and cannot be resized. 
   * In this case those dimensions are specified via the function 'set()', e.g.
   * \code
     JacS<double,Oregonator,ExprTmpl::MyVec> Deriv;
     Deriv.set(3,5);
     
   * \endcode
   * This is meant for non-xtended functions 
   */
  template<typename T, template<typename D1> 
	    class SOURCE, template<typename D2, class A = std::allocator<D2> > 
	    class VEC = ExprTmpl::MyVec> 
   class JacS{
   public:
    typedef T value_type;
     typedef CppAD::ADFun<T> ADFunType;
     typedef VEC<T> VecType;
     typedef CppAD::AD<T> ADType;
     typedef VEC<ADType> VecADType;
    //!in the case of CppAD this will be the Base of CppAD::AD<T>, to wit T. 
    typedef typename TypeAdapter<T>::Type BaseType; 
    typedef VEC<BaseType> BaseVecType;           
    typedef VecType ReturnType; //type of calculated Jacobian
    
   private:
    CppAD::ADFun<T> adseq_;        //AD sequence 
    VecType jac_;                   //store calculated Jacobian
  
    
  public:
    
    //!default constructor, uninitialized -- do nothing at all
    JacS(){
      //      BasicDataType<T>::check();
    }  


    //!use 'set()', 'set(..)' to initialise  the default constructed object 
    void set(size_t domdim, size_t randim){
       SOURCE<ADType> chemistry(randim); //source term
       VecADType X(domdim),  //domain space vector
	Y(randim);                    //range space vector

      
      CppAD::Independent(X);            //indenpendent variables, recording
                                        //must be declared here!! 


      Y = chemistry(X);
      
      adseq_.Dependent(X,Y);  //store sequence and stop recording 
    
      adseq_.optimize();
    }


    template<class W>
    void set_with_init(size_t domdim, size_t randim, const W& init){
      SOURCE<ADType> chemistry(randim); //source term
      VecADType X(domdim),  //domain space vector
	Y(randim);                    //range space vector
      
      //! for some functors it might be beneficial to initialize X first, since
      //! otherwise [0,0,....,0] is an illegal value to start with 
      if(init.size() != 0){
	for(size_t i = 0; i < domdim; ++i)
	  X[i] = init[i];
      }
      
      CppAD::Independent(X);            //indenpendent variables, recording
                                        //must be declared here!! 

      

      Y = chemistry(X);
      
      adseq_.Dependent(X,Y);  //store sequence and stop recording 
    
      adseq_.optimize();
    }




     
     //time dependency
    void set(const T& time, size_t domdim, size_t randim){
       //std::cout << "TYPE is "<<typeid(T).name() << std::endl;
       
      //BasicDataType<T>::check();
      SOURCE<ADType> chemistry(randim); //source term
      
      VecADType X(domdim),  //domain space vector
	Y(randim);                    //range space vector

      
      CppAD::Independent(X);            //indenpendent variables, recording
                                        //must be declared here!! 

      
     
      Y = chemistry(time,X);            //Note the time dependency
     
      adseq_.Dependent(X,Y);  //store sequence and stop recording 
    
      adseq_.optimize();
    }


    //time dependency
    template<class W>
    void set_with_init(const T& time, size_t domdim, size_t randim, const W& init){
       //std::cout << "TYPE is "<<typeid(T).name() << std::endl;
       
       //BasicDataType<T>::check();
      SOURCE<ADType> chemistry(randim); //source term
      
      VecADType X(domdim),  //domain space vector
	Y(randim);                    //range space vector

      if(init.size() != 0){
	for(size_t i = 0; i < domdim; ++i)
	  X[i] = init[i];
      }


      
      CppAD::Independent(X);            //indenpendent variables, recording
                                        //must be declared here!! 

      
     
      Y = chemistry(time,X);            //Note the time dependency
     
      adseq_.Dependent(X,Y);  //store sequence and stop recording 
    
      adseq_.optimize();
    }

    //!with parameter
    void set(size_t domdim, size_t randim, const VEC<T>& p){
            
      VecADType X(domdim),  //domain space vector
	Y(randim);                    //range space vector

      
      
      //ATTENTION: Declare Independet Variable before defining range space Y !!
      CppAD::Independent(X);            //indenpendent variables, recording
      
      SOURCE<ADType> chemistry(randim); //source term constructor

      Y = chemistry(X, p);    //PARAMETERIZED 

      adseq_.Dependent(X,Y);  //store sequence and stop recording 
   
      adseq_.optimize();
    
    }
 
    template<class W>
    void set_with_init(size_t domdim, size_t randim, const VEC<T>& p, const W& init){
            
      VecADType X(domdim),  //domain space vector
	Y(randim);                    //range space vector

      if(init.size() != 0){
	for(size_t i = 0; i < domdim; ++i)
	  X[i] = init[i];
      }
      
      //ATTENTION: Declare Independet Variable before defining range space Y !!
      CppAD::Independent(X);            //indenpendent variables, recording
      
      SOURCE<ADType> chemistry(randim); //source term constructor

      Y = chemistry(X, p);    //PARAMETERIZED 

      adseq_.Dependent(X,Y);  //store sequence and stop recording 
   
      adseq_.optimize();
    
    }


    ADFunType& get_AD_sequence() {return adseq_;}
    const ADFunType& get_AD_sequence() const {return adseq_;}
  
     inline VecType& jacobian(const VecType& v){
       
#ifndef NDEBUG
       if(adseq_.Domain() == 0)
	ADONIS_ERROR(DerivedError, "CppAD::ADFun object is still uninitialised (# of independent variables in operation sequence is ZERO and does not match size of function argument vector).\n Hence I cannot create Jacobian, master Marc.\n You should invoke member 'NablaS::set(...)' to revise that ;-)");

     
       if(adseq_.Domain() != v.size())
	  ADONIS_ERROR(DerivedError,"Size of evaluation point and those of AD sequence do NOT match");
#endif  

      jac_ = adseq_.Jacobian(v);
      return jac_;
     }




      inline VecType& sparse_jacobian(const VecType& v){

#ifndef NDEBUG
      if(adseq_.Domain() == 0)
	ADONIS_ERROR(DerivedError, "CppAD::ADFun object is still uninitialised (# of independent variables in operation sequence is ZERO and does not match size of function argument vector).\n Hence I cannot create sparse Jacobian, master Marc.\n You should invoke member 'NablaS::set(...)' to revise that ;-)");


      if(adseq_.Domain() != v.size())
	  ADONIS_ERROR(DerivedError,"Size of evaluation point and those of AD sequence do NOT match");
#endif
	 
      jac_ = adseq_.SparseJacobian(v); 
      return jac_;
      }
   

    //! overload sparse_jacobian to work with a predefined sparsity pattern
    //! cf. <a href="http://www.coin-or.org/CppAD/Doc/sparse_jacobian.xml">CppAD Online Documentation</a>
    template<class SPATTERN>
    inline VecType& sparse_jacobian(const VecType& v, const SPATTERN& p){
      #ifndef NDEBUG
      if(adseq_.Domain() == 0)
	ADONIS_ERROR(DerivedError, "CppAD::ADFun object is still uninitialised (# of independent variables in operation sequence is ZERO and does not match size of function argument vector).\n Hence I cannot create sparse Jacobian, master Marc.\n You should invoke member 'NablaS::set(...)' to revise that ;-)");


      if(adseq_.Domain() != v.size())
	  ADONIS_ERROR(DerivedError,"Size of evaluation point and those of AD sequence do NOT match");
#endif
      jac_ = adseq_.SparseJacobian(v,p); 
      return jac_;
    }
 
    inline const VecType& get_jacobian() const{
      return jac_;
    }

    inline VecType& get_jacobian() {
      return jac_;
    }

 
    friend std::ostream&  operator<<(std::ostream& os, const JacS& Jac){
      os << Jac.get_jacobian() << std::endl;
      return os;
    }

    /**
     * Assume \f$ f: R^m \rightarrow R^n, f:= [f_1,\ldots, f_n], f_i(x_1, \ldots, x_m) \f$. Then this member can be used to take the corresponging derivatives due to 'index'. For instance, if index = [0,4], then it calculates the derivatives of \f$ f_1\f$ and \f$ f_5\f$  with respect to \$ x_1\f$ and \f$ x_5\f$, which is equivalent to \f$ B^T\cdot Df(z).\f$ 
     */
    template<class V, class W>
    inline void cut_out_jacobian(const V& x, const W& index){
      cut_out_jacobian_from_reduced_function(adseq_,jac_,x,index);
    }

 
    //gives back some info about the the cppad stuff; in the case of errors
    //you might want to get some information about the ADFun sequence.
    inline void info() const{
      std::cout<<" Following numbers arose in operation sequence 'CppAD::ADFun': "<<std::endl<<  " ------------------------------------------------------------- "<<std::endl;

      std::cout<<std::endl<< "  # of variables (total):                  "<<adseq_.size_var()<<std::endl;
      std::cout<<            "  # of INDEPENDENT variables:              "<<adseq_.Domain()<<std::endl;

      std::cout<<            "  # of DEPENDENT variables:                "<<adseq_.Range()<<std::endl;
      
      std::cout<<            "  # of parameters:                         "<<adseq_.size_par()<<std::endl;
      std::cout<<            "  # of VecAD indices:                      "<< adseq_.size_VecAD()<<std::endl; 
     
      std::cout<<            "  # of TAYLOR coefficients (per variable): "<<adseq_.size_taylor()<<std::endl<<std::endl; 
    }


  };

#endif // USE_CPPAD

}//end namespace 

#endif
