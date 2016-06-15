#ifndef MATRIX_UNROLLER_HH
#define MATRIX_UNROLLER_HH

#include "commonfunctions.hh"
#include "compiletimearray.hh"
#include "../common/globalfunctions.hh"
#include "../moleculartransport/transportsettings.hh"
#include "../moleculartransport/prefactors.hh"
#include "../moleculartransport/transportutilities.hh"
#include "../common/smartassign.hh"
#include "../common/typeadapter.hh"
#include "../common/elementaryoperations.hh"

namespace Adonis{


  template<int I, int J, int N, int M>
  class MatrixUnroller{
  private:
    enum{i_ = (I+1) != N,
	 j_ = (J+1) != M};

  public:
    typedef MatrixUnroller<I,j_ ? (J+1):M,N,M> InnerLoopType;
    typedef MatrixUnroller<i_ ? (I+1):N,J,N,M> OuterLoopType;
   
    template<class D, class X, class MW>
    static inline void inner_loop(D& d, const X& x, const MW& mw){
      // std::cout << "  Inner Loop: "<< J << std::endl; //TEST
      d[I][J] = x[J]*mw[I];
      InnerLoopType::inner_loop(d,x,mw);
    }

    
    template<class D, class X, class MW>
    static inline void assign(D& d, const X& x, const MW& mw){
      //std::cout  << "Outer Loop: "<< I << std::endl; //TEST
      MatrixUnroller<I,0,N,M>::inner_loop(d,x,mw); //from J = 0 to M !!

      //UNCOMMENTED IT FOR THE SAKE OF BETTER TESTING
      OuterLoopType::assign(d,x,mw);
      
    }

    //============================ d - B*y ====================================
    template<class OP, class VIT, class MXIT>
    static inline typename OP::value_type inner_loop_vec_op_matrix_vec(const MXIT& B, const VIT& y){
      return  B[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()]*y[J] + InnerLoopType::template inner_loop_vec_op_matrix_vec<OP>(B,y);
    }

     template<class OP, class VIT, class MXIT>
     static inline void vec_op_matrix_vec(VIT d, const MXIT& B, const VIT& y){
       
       OP::apply(d[I],MatrixUnroller<I,0,N,M>::template inner_loop_vec_op_matrix_vec<OP>(B,y));
       OuterLoopType::template vec_op_matrix_vec<OP>(d,B,y);
     }

    //====================

    //!note that for complex arrays the transposition it the <BB> conjugate</BB> transposition
    template<class IT1, class IT2>
    static inline void inner_loop_transpose(IT1 tr, const IT2& a){
      //! tr has reverse rows an columns, and also swap N and M for iter tr
      smart_assign(tr[MatrixStorageFormat4LinearMemory<J,I,M,N,'r'>::offset()],Conj(a[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()]));
      InnerLoopType::inner_loop_transpose(tr,a);
    }

     template<class IT1, class IT2>
     static inline void transpose(IT1 tr, const IT2& a){
       MatrixUnroller<I,0,N,M>::inner_loop_transpose(tr,a);
       OuterLoopType::transpose(tr,a);
     }

    template<class V1, class V2, class V3>
    static inline V1& inner_loop_multiply_vector_entry_with_vector(V1& v, const V2& w, const V3& mtx){
      v[J] += w[I]*mtx[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()];
      return InnerLoopType::inner_loop_multiply_vector_entry_with_vector(v,w,mtx);
    }

    template<class V1, class V2, class V3>
    static inline V1& multiply_vector_entry_with_vector(V1& v, const V2& w, const V3& mtx){
      MatrixUnroller<I,0,N,M>::inner_loop_multiply_vector_entry_with_vector(v,w,mtx);
      return OuterLoopType:: multiply_vector_entry_with_vector(v,w,mtx);
    }

    //with Sdiag
    template<class V1, class V2, class V3,class V4>
    static inline V1& inner_loop_multiply_vector_entry_with_vector(V1& v, const V2& w, const V3& mtx, const V4& Sdiag){
      v[J] += w[I]*Sdiag[I]*mtx[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()];
      return InnerLoopType::inner_loop_multiply_vector_entry_with_vector(v,w,mtx,Sdiag);
    }

    template<class V1, class V2, class V3, class V4>
    static inline V1& multiply_vector_entry_with_vector(V1& v, const V2& w, const V3& mtx,const V4& Sdiag){
      MatrixUnroller<I,0,N,M>::inner_loop_multiply_vector_entry_with_vector(v,w,mtx,Sdiag);
      return OuterLoopType::multiply_vector_entry_with_vector(v,w,mtx,Sdiag);
    }

    //left multiplication with diagonal matrix -- NO size checking performed!!
    template<class DIAG, class MTX>
    static inline MTX& inner_loop_diag_matrix_multiplication(const DIAG& diag, MTX& mtx){
      //std::cout << "diag["<<I<<"] = "<< diag[I] << std::endl << "mtx["<<J<<"] = "<< mtx[J] << std::endl;
      mtx[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()] *= diag[I];
      return InnerLoopType::inner_loop_diag_matrix_multiplication(diag,mtx);
    }

    template<class DIAG, class MTX>
    static inline MTX& diag_matrix_multiplication(const DIAG& diag, MTX& mtx){
      MatrixUnroller<I,0,N,M>::inner_loop_diag_matrix_multiplication(diag,mtx);
      return OuterLoopType::diag_matrix_multiplication(diag,mtx);
    }

    
    template<class DYAD, class V1, class V2>
    static inline DYAD& inner_loop_dyad(DYAD& dy, const V1& v, const V2& w){
      dy[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()] = v[I]*w[J];
      return InnerLoopType::inner_loop_dyad(dy,v,w);
    }

    template<class DYAD, class V1, class V2>
    static inline DYAD& dyad(DYAD& dy, const V1& v, const V2& w){
      MatrixUnroller<I,0,N,M>::inner_loop_dyad(dy,v,w);
      return OuterLoopType::dyad(dy,v,w);
    }

    //======= symmetric matrix style===============
    //! NOTE: storage by row, cf. [1]
    //! References:
    //![1] [GOLUB/VAN LOAN, <I>"Matrix Computations"</I>, 3rd ed., §1.2, pp.20] 
    template<class DIFF, class T1, class T2>
    static inline void inner_loop_bin_diffusion_reduced(DIFF& diff, const T1& p, const T2& temp){
      typedef typename DIFF::DataType DataType;
      typedef typename DataType::RpvIndexType RpvIndexType;
      //std::cout << "  inner: "<< J << std::endl;
      //! extract the right index out of a given one
      smart_assign( diff.get_binary_diffusion_coefficient(SymmetricDenseMatrixAccess<I,J,DataType::rednspec>::offset()), diff.template binary_diffusion_coefficient<TransportSettings::molarMassInGrams>(AccessElement<RpvIndexType,I>::result,AccessElement<RpvIndexType,J>::result,p,temp) );
      InnerLoopType::inner_loop_bin_diffusion_reduced(diff,p,temp);
    }

    template<class DIFF, class T1, class T2>
    static inline void inner_loop_bin_diffusion(DIFF& diff, const T1& p, const T2& temp){
      //std::cout << "  inner: "<< J << std::endl;
      smart_assign( diff.get_binary_diffusion_coefficient(SymmetricDenseMatrixAccess<I,J,DIFF::NSPEC>::offset()), diff.template binary_diffusion_coefficient<TransportSettings::molarMassInGrams>(I,J,p,temp) );
      InnerLoopType::inner_loop_bin_diffusion(diff,p,temp);
    }

    //outer loop
    template<class DIFF, class T1, class T2>
    static inline void compute_binary_diffusion_matrix_reduced(DIFF& diff, const T1& p, const T2& temp){
	//std::cout << "OUTER: "<< I << std::endl;
      //! from J = I to M!!
      MatrixUnroller<I,I,N,M>::inner_loop_bin_diffusion_reduced(diff,p,temp);
      OuterLoopType::compute_binary_diffusion_matrix_reduced(diff,p,temp);
    }

     //outer loop
    template<class DIFF, class T1, class T2>
    static inline void compute_binary_diffusion_matrix(DIFF& diff, const T1& p, const T2& temp){
	//std::cout << "OUTER: "<< I << std::endl;
      //! from J = I to M!!
      MatrixUnroller<I,I,N,M>::inner_loop_bin_diffusion(diff,p,temp);
      OuterLoopType::compute_binary_diffusion_matrix(diff,p,temp);
    }

    //==============================================
  

    //--------- compute mc diffusion ---------------
    template<class TPTDATATYPE, class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void inner_loop_mc_diffusion_reduced(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){
      typedef typename TPTDATATYPE::RpvIndexType RpvIndexType;
      
      smart_assign(mcdiff[MatrixStorageFormat4LinearMemory<I,J,N,M,'R'>::offset()], xfrac[I]*16./25.*(temp/(ChooseRightPrefactor<T1,TransportSettings::pressureInBar>::in_bar(p)))*meanmw/mw[AccessElement<RpvIndexType,J>::result]*(inv[MatrixStorageFormat4LinearMemory<I,J,N,M,'R'>::offset()] - inv[MatrixStorageFormat4LinearMemory<I,I,N,M,'R'>::offset()]) );
    
      InnerLoopType::template inner_loop_mc_diffusion_reduced<TPTDATATYPE>(mcdiff,p,temp,xfrac,mw,meanmw,inv);
    }

    
    template<class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void inner_loop_mc_diffusion(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){
      
      smart_assign(mcdiff[MatrixStorageFormat4LinearMemory<I,J,N,M,'R'>::offset()], xfrac[I]*16./25.*(temp/(ChooseRightPrefactor<T1,TransportSettings::pressureInBar>::in_bar(p)))*meanmw/mw[J]*(inv[MatrixStorageFormat4LinearMemory<I,J,N,M,'R'>::offset()] - inv[MatrixStorageFormat4LinearMemory<I,I,N,M,'R'>::offset()]) );
      
      InnerLoopType::inner_loop_mc_diffusion(mcdiff,p,temp,xfrac,mw,meanmw,inv);
    }

    template<class TPTDATATYPE, class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void compute_multicomponent_diffusion_matrix_reduced(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){
      MatrixUnroller<I,0,N,M>::template inner_loop_mc_diffusion_reduced<TPTDATATYPE>(mcdiff,p,temp,xfrac,mw,meanmw,inv); //! from J = 0 to M since this matrix isn't symmetric
      OuterLoopType::template compute_multicomponent_diffusion_matrix_reduced<TPTDATATYPE>(mcdiff,p,temp,xfrac,mw,meanmw,inv);
    }
    

    template<class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void compute_multicomponent_diffusion_matrix(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){
      MatrixUnroller<I,0,N,M>::inner_loop_mc_diffusion(mcdiff,p,temp,xfrac,mw,meanmw,inv); //! from J = 0 to M since this matrix isn't symmetric
      OuterLoopType::compute_multicomponent_diffusion_matrix(mcdiff,p,temp,xfrac,mw,meanmw,inv);
    }

    //----------------------------------------------
  
     //············ compute mixture-averaged viscosity ······················
    template<class TPTDATATYPE, bool BOOL, class X, class VIS>
    static inline typename TypeAdapter<typename X::value_type>::Type inner_loop_eta_denom(const X& Xfrac, const VIS& vis){
      return convert_number(Xfrac[J])*Phi4MixVis<I,J,TPTDATATYPE,BOOL>::calculate(vis) +  InnerLoopType::template inner_loop_eta_denom<TPTDATATYPE,BOOL>(Xfrac,vis);
    }

    template<class TPTDATATYPE, bool BOOL, class X, class VIS>
    static inline typename TypeAdapter<typename X::value_type>::Type compute_mixture_averaged_viscosity(const X& Xfrac, const VIS& vis){
      return convert_number(Xfrac[I])*vis[I]/(MatrixUnroller<I,0,N,M>::template inner_loop_eta_denom<TPTDATATYPE,BOOL>(Xfrac,vis)) + OuterLoopType::template compute_mixture_averaged_viscosity<TPTDATATYPE,BOOL>(Xfrac,vis);
    }

    //······································································

  };


  //· END OF RECURSION · END OF RECURSION · END OF RECURSION · END OF RECURSION 
  //end of INNER Recursion
  template<int I, int M, int N>
  class MatrixUnroller<I,M,N,M>{
  public:
    template<class D, class X, class MW>
    static inline void inner_loop(D& d, const X& x, const MW& mw){ } //do nothing

    template<class OP, class VIT, class MXIT>
    static inline typename OP::value_type inner_loop_vec_op_matrix_vec(const MXIT& B, const VIT& y){
      typedef typename OP::value_type value_type;
      return value_type();
    }


     template<class IT1, class IT2>
     static inline void inner_loop_transpose(IT1 tr, const IT2& a){}

     template<class V1, class V2, class V3>
     static inline V1& inner_loop_multiply_vector_entry_with_vector(V1& v, const V2& w, const V3& mtx){ return v;}

    //with Sdiag
    template<class V1, class V2, class V3, class V4>
    static inline V1& inner_loop_multiply_vector_entry_with_vector(V1& v, const V2& w, const V3& mtx, const V4& Sdiag){ return v;}


    template<class DIAG, class MTX>
    static inline MTX& inner_loop_diag_matrix_multiplication(const DIAG& diag, MTX& mtx){ return mtx;}

    template<class DYAD, class V1, class V2>
    static inline DYAD& inner_loop_dyad(DYAD& dy, const V1& v, const V2& w){
      return dy;
    } 

     //======= symmetric matrix style ===============
    template<class DIFF, class T1, class T2>
    static inline void inner_loop_bin_diffusion_reduced(DIFF& diff, const T1& p, const T2& temp){}

    template<class DIFF, class T1, class T2>
    static inline void inner_loop_bin_diffusion(DIFF& diff, const T1& p, const T2& temp){}
    //===============================================

    
    //--------- compute mc diffusion ---------------
    template<class TPTDATATYPE, class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
  static inline void inner_loop_mc_diffusion_reduced(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){}
    
    template<class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
  static inline void inner_loop_mc_diffusion(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){}
    
    //-----------------------------------------------

    //············ compute mixture-averaged viscosity ······················
    template<class TPTDATATYPE, bool BOOL, class X, class VIS>
    static inline typename TypeAdapter<typename X::value_type>::Type inner_loop_eta_denom(const X& Xfrac, const VIS& vis){
      typedef typename TypeAdapter<typename X::value_type>::Type Type;
      return Type();
    }
    //········································································

  };

  //end of OUTER Recursion
  template<int J, int N, int M>
  class MatrixUnroller<N,J,N,M>{
  public:
    template<class D, class X, class MW>
    static inline void assign(D& d, const X& x, const MW& mw){}


    template<class OP, class VIT, class MXIT>
    static inline void vec_op_matrix_vec(VIT d, const MXIT& B, const VIT& y){}

    template<class IT1, class IT2>
    static inline void transpose(IT1 tr, const IT2& a){}
    
    template<class V1, class V2, class V3,class V4>
    static inline V1& multiply_vector_entry_with_vector(V1& v, const V2& w, const V3& mtx,const V4& Sdiag){return v;}
    
    //with Sdiag
    template<class V1, class V2, class V3>
    static inline V1& multiply_vector_entry_with_vector(V1& v, const V2& w, const V3& mtx){return v;}

    template<class DIAG, class MTX>
    static inline MTX& diag_matrix_multiplication(const DIAG& diag, MTX& mtx){
      return mtx;
    }

    template<class DYAD, class V1, class V2>
    static inline DYAD& dyad(DYAD& dy, const V1& v, const V2& w){return dy;}

    //======= symmetric matrix style ===============
     template<class DIFF, class T1, class T2>
    static inline void compute_binary_diffusion_matrix_reduced(DIFF& diff, const T1& p, const T2& temp){}

    template<class DIFF, class T1, class T2>
    static inline void compute_binary_diffusion_matrix(DIFF& diff, const T1& p, const T2& temp){}
    //==============================================
  

    //--------- compute mc diffusion ---------------
    template<class TPTDATATYPE, class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void compute_multicomponent_diffusion_matrix_reduced(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){}

    template<class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void compute_multicomponent_diffusion_matrix(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){}

    //-----------------------------------------------


    //······················ mix. avg. viscosity ·······························
    template<class TPTDATATYPE, bool BOOL, class X, class VIS>
    static inline typename TypeAdapter<typename X::value_type>::Type compute_mixture_averaged_viscosity(const X& Xfrac, const VIS& vis){
      typedef typename TypeAdapter<typename X::value_type>::Type Type;
      return Type();
    }
    //··········································································
   
  };



  //end of recursion -- COMPLETE END -- seems not to be mandatory at all
  template<int N, int M>
  class MatrixUnroller<N,M,N,M>{
  public:
    template<class D, class X, class MW>
    static inline void inner_loop(D& d, const X& x, const MW& mw){} //do nothing
    

    template<class OP, class VIT, class MXIT>
    static inline typename OP::value_type inner_loop_vec_op_matrix_vec(const MXIT& B, const VIT& y){
      typedef typename OP::value_type value_type;
      return value_type();
    }

    template<class IT1, class IT2>
    static inline void inner_loop_transpose(IT1 tr, const IT2& a){}

    template<class IT1, class IT2>
    static inline void transpose(IT1 tr, const IT2& a){}

    template<class OP, class VIT, class MXIT>
    static inline void vec_op_matrix_vec(VIT d, const MXIT& B, const VIT& y){}

    template<class D, class X, class MW>
    static inline void  assign(D& d, const X& x, const MW& mw){} //do nothing
  
    
    //======= symmetric matrix style ===============
    template<class DIFF, class T1, class T2>
    static inline void inner_loop_bin_diffusion_reduced(DIFF& diff, const T1& p, const T2& temp){}
    
    template<class DIFF, class T1, class T2>
    static inline void inner_loop_bin_diffusion(DIFF& diff, const T1& p, const T2& temp){}

    template<class DIFF, class T1, class T2>
    static inline void compute_binary_diffusion_matrix_reduced(DIFF& diff, const T1& p, const T2& temp){}

    template<class DIFF, class T1, class T2>
    static inline void compute_binary_diffusion_matrix(DIFF& diff, const T1& p, const T2& temp){}
    //===============================================
  
    //--------- compute mc diffusion ---------------
    template<class TPTDATATYPE, class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
  static inline void inner_loop_mc_diffusion_reduced(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){}

    template<class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
  static inline void inner_loop_mc_diffusion(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){}



    template<class TPTDATATYPE, class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void compute_multicomponent_diffusion_matrix_reduced(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){}

    template<class DITER, class X, class MW, class T1, class T2, class T3, class INVITER>
    static inline void compute_multicomponent_diffusion_matrix(DITER mcdiff, const T1& p, const T2& temp, const X& xfrac, const MW& mw, const T3& meanmw, const INVITER& inv){}

    //-----------------------------------------------


     //············ compute mixture-averaged viscosity ······················
    template<class TPTDATATYPE, bool BOOL, class X, class VIS>
    static inline typename TypeAdapter<typename X::value_type>::Type inner_loop_eta_denom(const X& Xfrac, const VIS& vis){
      typedef typename TypeAdapter<typename X::value_type>::Type Type;
      return Type();
    }


    template<class TPTDATATYPE, bool BOOL, class X, class VIS>
    static inline typename TypeAdapter<typename X::value_type>::Type compute_mixture_averaged_viscosity(const X& Xfrac, const VIS& vis){
      typedef typename TypeAdapter<typename X::value_type>::Type Type;
      return Type();
    }
    //········································································
  };


} //end namespace

#endif
