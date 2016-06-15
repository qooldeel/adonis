#ifndef COMPUTE_MANIFOLD_VIA_OPTIMIZATION_HH
#define COMPUTE_MANIFOLD_VIA_OPTIMIZATION_HH

#include <iostream>
#include <typeinfo>
#include <string>


#include "../dunexternal/dunextensions.hh"
#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../misc/useful.hh"
#include "../misc/misctmps.hh"

#include "../dunexternal/whatweneed.hh"


#include "../dunexternal/normalspace.hh"

#include "rvp.hh"
#include "indexmanipulation.hh"

//! Jochen's stuff -- since included in Makefile via…
//! … <TT> -I${HOME}/MODELREDUCTION/MoRe </TT> 
//#include <MoRe.hpp>   //with BDF integrator and IPOPT
#include <MoRe_core.hpp>  //only gGN


namespace Adonis{

  /**
   * \brief In 'mtv' the last entries in each 'column' represent time. You can get rid of the temperature in <I> isothermal</I> mechanisms by invoking this routine
   */
  template<class W, class V>
  inline void discard_temperature(W& v, const V& mtv, int reddim){
    assert(mtv.size()%reddim == 0);
    unsigned dim = (mtv.size()/reddim) - 1; //mtv.size()-1;
    
    //std::cout << "dim w.o. temperature = " << dim << std::endl<<std::endl;
    //necessary because of push_back:
    v.clear();

    unsigned k = 0;
    unsigned c = 1;
    for(unsigned i = 0; i < mtv.size(); i += c){
      k++;
      
      if(k <= dim){
	v.push_back(mtv[i]);
	c = 1;
      }
      if(k == dim){
	c = 2;
	k = 0;
      }
   
    }
  }


  /**
   * \brief needed for isothermal setting in conjunction with Jochen's tool
   */
  template<int N, bool B> class DiscardTemperature;

  template<int N>
  class DiscardTemperature<N,true>{ //isothermal (\f$T = \textrm{const}\f$)
  public:
    template<class W, class V>
    static inline W& tangentspace_wo_temp(W& v, const V& mtv,int reddim){
      adonis_assert(mtv.size()%(N+1) == 0);
      discard_temperature(v,mtv,reddim);
      return v;
    }
  };

  template<int N>
  class DiscardTemperature<N,false>{ //temperature-dependent
  public:
    template<class W, class V>   //do nothing
    static inline const V& tangentspace_wo_temp(W& v, const V& mtv,int reddim){
      adonis_assert(mtv.size()%N == 0); //mtv is always ns+1
      return mtv;
    }
  };


  /**
   * \brief Compute \f$ z^{M}(r), T(r), N(r)\f$ via variational principle.
   * [LEBIEDZ, SIEHR: "A variational principle for computing slow invariant manifolds in dissipative dynamical systems", SIAM J. SCI. COMPUT., Vol 33, 2011]
   *
   * \tparam T precision of computation
   * \tparam N dimension of full composition
   * \tparam R dimension of reduced composition (R < N)
   * \tparam ISOTHERMAL (default: true, i.e. \f$ T = \textrm{const} \f$, false refers to temperature-dependency via Arrhenius' law)
   * This is important when preprocessor <TT> REDUCTION_TO_BE_USED </TT> is 2 (2nd, i.e. CPA, approximation)
   *
   * NOTE: everything is computed with temperature as the \f$(n_s +1)\f$-st species.
   *
   * \code
      //'temperature-dependent' mechanism:
      ComputeManifold<double,7,2,false> Mani; (then N = n_s + 1)
      //same mechanism without temperature: 
       ComputeManifold<double,6,2,true> Mani; (then N = n_s)
   * \endcode
   */
  template<class T, int N, int R, bool ISOTHERMAL = true>
  class ComputeManifold{
  public:
    typedef Essential<T,N,R,false> Type;   //! false: \f$ N^T\f$ is considered
    typedef typename Type::MatrixType4NormalSpace MatrixType4NormalSpace;

    typedef typename Type::CompositionType CompositionType;
    typedef typename Type::TangentSpaceType TangentSpaceType;
    typedef typename Type::TransposedNormalSpaceType TransposedNormalSpaceType;

    typedef more::MoRe<T> VariationalPrinciple4SIM;
    typedef js::VecS VecS;
    typedef js::VecD VecD;
    
    typedef std::string StrType; 

    typedef Dune::FieldVector<int,R> CoordinateType;
    typedef Dune::FieldVector<int,N-R> UnrepIndexType;
    typedef Dune::FieldMatrix<T,R,N> BTType;
    typedef Dune::FieldMatrix<T,N,N-R> UType;


    //default constructor
    ComputeManifold():isActivated_(false),isFilled_(false),isDone_(false){}

    template<class C, class S>
    ComputeManifold(const S& pvname, const C& vd, const CoordinateType& rpvindex):pvname_(pvname.begin(),pvname.end()),isFilled_(false),isDone_(false){
      activate(pvname,vd,rpvindex);
    }

    template<class ITER, class C>
    ComputeManifold(ITER i1, ITER i2,  const C& vd, const CoordinateType& rpvindex):pvname_(i1,i2),isFilled_(false),isDone_(false){
      activate(pvname_,vd,rpvindex);
    }
    
    template<class C,class S>
    void activate(const S& sv, const C& vd, const CoordinateType& rpvindex){
      adonis_assert(((int)sv.size() == (int)(std::distance(vd.begin(),vd.end()))) && ((int)sv.size() == R) && ((int)std::distance(rpvindex.begin(),rpvindex.end()) == R ));
      pvname_.resize(R);
      
      //std::cout << "ComputeManifold.activate(…) -- FIRST --vd (input) = ";print_all(vd,16);

      ManipulateDynamicArray<StrType,int,0,R>::fill(&pvname_[0],sv.begin(),sv.end());
      
      //std::cout<< "Names = "; print_all(pvname_);
      pvvalue_.resize(R);
      ManipulateDynamicArray<T,int,0,R>::fill(&pvvalue_[0],vd.begin(),vd.end());
      
      //std::cout << "ComputeManifold.activate(…) --FIRST-- assigned (output) = "; print_all(pvvalue_,16);

      more_.first(z_,mtv_,pvname_,pvvalue_);
      // std::cout << "FIRST CALL MORE: z_ = "; print_all(z_);
      // abort();

      ManipulateDynamicArray<CompositionType,int,0,N>::fill(essex_.get_composition(),z_);

      // std::cout <<"more_first z_ to essex_ = "<< essex_.get_composition()<<std::endl;

#if REDUCTION_TO_BE_USED == 2
      // print_all(DiscardTemperature<N-1,ISOTHERMAL>::tangentspace_wo_temp(mtv_without_temperature_,mtv_,R));

      FillMatrix<ColumnMajor>::now(essex_.get_tangentspace(),DiscardTemperature<N,ISOTHERMAL>::tangentspace_wo_temp(mtv_without_temperature_,mtv_,R));

      essex_.get_normalspace() = normalspace_transposed<R>(essex_.get_tangentspace());

#else
      //!all other reduction methods apart from CPA, have useless T, and N^T,
      //!coz we don't need 'em
      essex_.get_tangentspace() = T();
      essex_.get_normalspace() = T();
#endif

      //form index, B^T and U 
      (*this).form_data_4_rpvs(rpvindex);


      isActivated_ = true;
   
      isFilled_ = true;
 
      // std::cout << "zM_ (FIRST) = "<< (*this).get_z() << std::endl;

    }

    

    template<class C>
    Type& evaluate(const C& reduced){
      adonis_assert((int)std::distance(reduced.begin(),reduced.end()) == R);
#ifndef NDEBUG
      if(!isActivated_){
	ADONIS_ERROR(IllegalInitialization, "You must first invoke 'more.first(....)' with a (default) value for the reduced composition.");
      }
#endif
      adonis_assert((int)pvname_.size() == R);
      adonis_assert(R == (int)container_size(reduced));
      
      //std::cout << "ComputeManifold.evaluate(…) -- reduced (input) = ";print_all(reduced,16);

      ManipulateDynamicArray<T,int,0,R>::fill(&pvvalue_[0],reduced.begin(),reduced.end());

      //std::cout << "ComputeManifold.evaluate(…) -- pvvalue (output) = ";print_all(pvvalue_,16);
      
      more_.warm(z_,mtv_,pvname_,pvvalue_); //compute stuff

      //assign stuff to proper containers
      ManipulateDynamicArray<CompositionType,int,0,N>::fill(essex_.get_composition(), z_);   //! \f$ z(r) \f$
      
      //std::cout << "mtv_ = "; print_all(mtv_);
      
      //mtv_without_temperature_.clear();
#if REDUCTION_TO_BE_USED == 2
      FillMatrix<ColumnMajor>::now(essex_.get_tangentspace(),DiscardTemperature<N,ISOTHERMAL>::tangentspace_wo_temp(mtv_without_temperature_,mtv_,R));
				   //mtv_); //!\f$ T(r)\f$
      
      essex_.get_normalspace() = normalspace_transposed<R>(essex_.get_tangentspace());  //!\f$ N^T(r)\f$
#endif


      isFilled_ = true;

      return essex_;
    }

    const CompositionType& get_z() const {adonis_assert(isFilled_); return essex_.get_composition();}
    CompositionType& get_z() {adonis_assert(isFilled_); return essex_.get_composition();}
    
    const TangentSpaceType& get_T() const {adonis_assert(isFilled_); return essex_.get_tangentspace();}
    TangentSpaceType& get_T() {adonis_assert(isFilled_); return essex_.get_tangentspace();}
    
    const TransposedNormalSpaceType& get_N_T() const {adonis_assert(isFilled_); return essex_.get_normalspace();}
    TransposedNormalSpaceType& get_N_T() {adonis_assert(isFilled_); return essex_.get_normalspace();}
    

    const VecS& get_names() const {return pvname_;}
    VecS& get_names() {return pvname_;}

    void print_names() const {print_all(pvname_);}


    void form_data_4_rpvs(const CoordinateType& idx){
      if(!isDone_){
	rpvindex_ = idx;
	
	IndexManipulator<N,R,std::vector> IM;
	IM.create_unrepresented_index(rpvindex_);
	
	ReactionProgressVariables<T,N,R> RPV(B_T_,U_);
	RPV.create_B_T(rpvindex_); //B^T is formed
	RPV.create_U(IM.get_index()); // U is formed
      
	isDone_ = true;
	//std::cout << "rpvindex_ = "<< rpvindex_ << std::endl;
	//std::cout << "B^T = "<< std::endl<< B_T_ << std::endl;
	//std::cout << "U = "<< std::endl<< U_ << std::endl;
      }
    }

    const BTType& B_T() const {
      adonis_assert(isDone_);
      return B_T_;
    }
    BTType& B_T() {
      adonis_assert(isDone_);
      return B_T_;
    }

    const UType& U() const {
      adonis_assert(isDone_);
      return U_;
    }
    UType& U() {
      adonis_assert(isDone_);
      return U_;
    }
    
    const CoordinateType& index() const {
      adonis_assert(isDone_);
      return rpvindex_;
    }
    CoordinateType& index() {
      adonis_assert(isDone_);
      return rpvindex_;
    }

  private:
    VariationalPrinciple4SIM more_;
    Type essex_;
    VecS pvname_;
    VecD z_,                   //full value                   
		  mtv_,        //column vectors spanning tangent space 
 		  pvvalue_;   //reduced value
    VecD mtv_without_temperature_;

    mutable bool isActivated_;
    mutable bool isFilled_;
    mutable bool isDone_;

    CoordinateType rpvindex_;
    BTType B_T_;
    UType U_;
  };

}

#endif 
