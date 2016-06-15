#ifndef XTENDED_SCALING_FUNCTIONALITIES_FOR_BOTH_SPARSE_AND_DENSE_MTXS_HH
#define XTENDED_SCALING_FUNCTIONALITIES_FOR_BOTH_SPARSE_AND_DENSE_MTXS_HH

#include <math.h>
#include "../common/globalfunctions.hh"
#include "../common/fancymessages.hh"
#include "../misc/useful.hh"

namespace Adonis{

  /**
   * \brief select at compile time if some feature gonna to be executed
   */
  template<bool SCALE, class V> class CTscalingFeatures;

  template<class V>
  class CTscalingFeatures<true,V>{ //scaling switched on
  public:
    typedef std::size_t SizeType;
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    static const bool SwitchValue = true;
    
    static inline void resize(SizeType& dim, SizeType n, V& D1, V& D2, const BaseType& tol){
      dim = n;
      D1.resize(n);
      D2.resize(n,tol);
    }

    template<class Y>
    static inline void compute_D2(V& D2, const Y& y, const BaseType& tol1, const BaseType& tol2){
      //std::cout << "D2_ = ";
      for(SizeType i = 0; i < D2.size(); ++i){
	//scalbn uses integer power
	D2[i] = scalbn(1.,(int)rint(Logb<2>(tol1*Abs(y[i]) + tol2))); //overwrite D2_
	//std::cout << D2_[i] << " ";
      }
      //std::cout << std::endl;
    }
    
    template<class X>
    static inline void retrieve_solution(X& xtilde, const V& D2){
	for(SizeType i = 0; i < D2.size(); ++i){
	  xtilde[i] *= D2[i];
	}
    }

  };
  

  template<class V>
  class CTscalingFeatures<false,V>{ //scaling switched OFF
  public:
    typedef std::size_t SizeType;
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    static const bool SwitchValue = false;

    static inline void resize(SizeType& dim, SizeType n, V& D1, V& D2, const BaseType& tol){}  //do nothing

    template<class Y>
    static inline void compute_D2(V& D2, const Y& y, const BaseType& tol1, const BaseType& tol2){}  //do nothing
    
    template<class X> //do nothing
    static inline void retrieve_solution(X& xtilde, const V& D2){}

  };
  



  /**
   * \brief Scale SQUARE system -- General handler class
   */
  template<class V, bool SCALE, class OBJ>  
  class BasicScaling{
  public:
    typedef std::size_t SizeType;
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;

    //! use features of derived class
  private:
    OBJ& ref_2_derived(){return static_cast<OBJ&>(*this);}
    

  protected:
    SizeType dim_;      
    V D1_, D2_;  //diagonal scaling mtxs, stored as RACs
   
  public: //common functions that can be invoked by derived classes
    

    const V& get_D1() const {return D1_;}
    const V& get_D2() const {return D2_;}

  
    void resize(SizeType n, const BaseType& tol){
      CTscalingFeatures<SCALE,V>::resize(dim_,n,D1_,D2_,tol);
    }

    //! there seems to be a BUG with logb(x). Anyway, logb(x) = log2(|x|) on
    //! those platforms for which FLT_RADIX = 2 (these will be most ;) )
    template<class Y>
    void compute_D2(const Y& y, const BaseType& tol1, const BaseType& tol2){
      CTscalingFeatures<SCALE,V>::compute_D2(D2_,y,tol1,tol2);
    }

    template<class X>
    void retrieve_solution(X& xtilde){
      CTscalingFeatures<SCALE,V>::retrieve_solution(xtilde,D2_);
    }
    
  };


  /***************************************************************************/
  /**************************** DERIVED CLASSES ******************************/
  /***************************************************************************/
  //default class, no scaling at all
  /**
   * \brief Scale (sparse) square linear system or not
   *
   * \tparam SCALE    switch on/off scaling (true/false)
   * \tparam SPARSE   use sparse or dense matrix (true=sparse, false=dense)
   * \tparam V        rhs side vector type
   * \tparam PATT     sparsity pattern type (only meaningful for sparse mtxs)
   */
  template<bool SCALE, bool SPARSE, class V, class PATT>
  class ScaleSquareLinearSystem: public BasicScaling<V,SCALE,ScaleSquareLinearSystem<SCALE,SPARSE,V,PATT> >{
  public:
    typedef BasicScaling<V,SCALE,ScaleSquareLinearSystem<SCALE,SPARSE,V,PATT> > BaseClassType; 
    typedef typename BaseClassType::SizeType SizeType;
    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::BaseType BaseType;
    
    //dummy arguments
    ScaleSquareLinearSystem(const PATT& patt, const BaseType& tol, SizeType n){
      BaseClassType::dim_ = n; //proper initialization, though it won't be used
    }
   
    template<class MTXVALIT, class RHSIT, class Y>
    void scale(MTXVALIT mtxvalit, RHSIT b, const Y& y, const BaseType& tol1, const BaseType& tol2){} //do nothing
    
  private:
  };

  //! partial specializations -- Apply scaling
  //! SPARSE case
  //! CAUTION: We expect the <B>values</B> of the sparse matrix to be in 
  //!          compressed sparse row storage (CSR) form!
  template<class V, class PATT>//scale                                  scale
  class ScaleSquareLinearSystem<true,true,V,PATT>: public BasicScaling<V,true,ScaleSquareLinearSystem<true,true,V,PATT> >{
    typedef typename PATT::SimpleType SetType;  //private aliases
    typedef typename SetType::const_iterator SetIterType;
  public://                scale                       scale
    typedef BasicScaling<V,true,ScaleSquareLinearSystem<true,true,V,PATT> >BaseClassType;
    typedef typename BaseClassType::SizeType SizeType;
    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::BaseType BaseType;
    
    ScaleSquareLinearSystem(const PATT& patt, const BaseType& tol, SizeType n):pattern_(patt){
      adonis_assert(patt.size() == n);
      adonis_assert(patt.rows() == patt.cols());
      BaseClassType::resize(patt.size(),tol);
      FancyMessages().nice_output("==============================================\n Apply SPARSE SCALING of sparse linear system \n==============================================\n",36);
    }

    //! scalbn(x,n) = x*FLT_RADIX^n
    template<class MTXVALIT, class RHSIT, class Y>
    void scale(MTXVALIT mtxValIt, RHSIT bIt, const Y& y, const BaseType& tol1, const BaseType& tol2){
      MTXVALIT startIt = mtxValIt; //start iterator position

      BaseClassType::compute_D2(y,tol1,tol2);
      //! A = A*D2, calculate D1 and scale rhs vector b in one sweep
      BaseType sum = BaseType();
      for(SizeType i = 0; i < pattern_.size(); ++i){
	sum = BaseType();  //reset
	for(SetIterType it = pattern_[i].begin(); it != pattern_[i].end(); ++it){
	  *(mtxValIt) *= BaseClassType::D2_[*it];
	  sum += Abs(*mtxValIt);

	  mtxValIt++;
	}
	BaseClassType::D1_[i] = scalbn(1.,-(int)ceil(Logb<2>(sum))); //overwrite D1
	*(bIt++) *= BaseClassType::D1_[i]; //overwrite  rhs
      }
      

      mtxValIt = startIt;  //reset to starting iterator position!!
      //A = D1_*A, where A = A*D2
      for(SizeType i = 0; i < pattern_.size(); ++i){
	for(SetIterType it = pattern_[i].begin(); it != pattern_[i].end(); ++it){
	  *(mtxValIt) *= BaseClassType::D1_[i];
	  mtxValIt++;
	}
      }
    }
    
  private:
    const PATT& pattern_;
  };


  //! partial specialization
  //! DENSE case
  template<class V, class PATT> //<SCALE,SPARSE> = <true,false>           scale
  class ScaleSquareLinearSystem<true,false,V,PATT>: public BasicScaling<V,true,ScaleSquareLinearSystem<true,false,V,PATT> >{
  public://                scale                        scale     
    typedef BasicScaling<V,true,ScaleSquareLinearSystem<true,false,V,PATT> >BaseClassType;
    typedef typename BaseClassType::SizeType SizeType;
    typedef typename BaseClassType::value_type value_type;
    typedef typename BaseClassType::BaseType BaseType;
    
    //! the first 2 args are mere dummies
    ScaleSquareLinearSystem(const PATT& patt, const BaseType& tol, SizeType n){
      BaseClassType::resize(n,tol);
      FancyMessages().nice_output("==============================================\n Apply DENSE  SCALING of sparse linear system \n==============================================\n",35);
    }

    template<class MTXVALIT, class RHSIT, class Y>
    void scale(MTXVALIT mtxValIt, RHSIT bIt, const Y& y, const BaseType& tol1, const BaseType& tol2){
      //! note the *(a+i) is equivalent to a[i], but more general in case of 
      //! a lacking []-operator for a given RAC
      BaseClassType::compute_D2(y,tol1,tol2);
      //dense looping...
      //A = A*D2, calculate D1 and scale rhs vector b in one sweep
      BaseType sum = BaseType();
      for(SizeType i = 0; i < y.size(); ++i){
	sum = BaseType(); //reset
	for(std::size_t j = 0; j < y.size(); ++j){
	  *(mtxValIt + RowMajor::offset(i,j,y.size())) *= BaseClassType::D2_[j];
	  sum += Abs(*(mtxValIt + RowMajor::offset(i,j,y.size())));
	  //mtxValIt++;
	}
	BaseClassType::D1_[i] = scalbn(1.,-(int)ceil(Logb<2>(sum))); //overwrite D1
	*(bIt++) *= BaseClassType::D1_[i]; //overwrite  rhs
      }
       
      //A = D1_*A, where A = A*D2
      for(SizeType i = 0; i < y.size(); ++i){
      	for(SizeType j = 0; j < y.size(); ++j){
      	  *(mtxValIt + RowMajor::offset(i,j,y.size())) *= BaseClassType::D1_[i];
      	}
      }
    }

  };

} //end namespace

#endif
