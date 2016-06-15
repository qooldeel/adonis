#ifndef THREE_NESTED_LOOPS_UNROLLED_HH
#define THREE_NESTED_LOOPS_UNROLLED_HH

#include "commonfunctions.hh"
#include "compiletimearray.hh"
#include "../moleculartransport/prefactors.hh"
#include "../moleculartransport/transportsettings.hh"
#include "../common/smartassign.hh"

namespace Adonis{
  
  template<int I, int J, int K, int N, int M, int R>
  class ThreeNestedLoopsUnroller{
  private:
    enum{i_ = (I+1) != N,
	 j_ = (J+1) != M,
	 k_ = (K+1) != R};

  public:
    typedef ThreeNestedLoopsUnroller<i_?(I+1):N,J,K,N,M,R> OuterLoopType;
    typedef ThreeNestedLoopsUnroller<I,j_?(J+1):M,K,N,M,R> MediumLoopType;
    typedef ThreeNestedLoopsUnroller<I,J,k_?(K+1):R,N,M,R> InnerLoopType;

    
    //! at first, I thought that the \f$\delta_{i,j}\f$ of the \f$L^{00,00}\f$ matrix were actually the reduced dipole moment for the collision, \f$ \delta^*_{i,j}\f$. In the meantime, I've been convinced more and more that this symbol is the Kronecker delta!
    //HERE, YOU CAN USE TUNED FLT PT ARITHMETIC -- NOT DONE SO FAR:)
    template<class TPTDATATYPE, class X, class MW, class PITER>
    static inline typename X::value_type inner_loop_L0000_reduced(const X& x, const MW& mw, PITER d){
      typedef typename TPTDATATYPE::RpvIndexType RpvIndexType;
 
      //std::cout << "D[i,k] = "<<d[SymmetricDenseMatrixAccess<I,K,N>::offset()] <<std::endl;
      //std::cout << "MOLAR MASS["<<AccessElement<RpvIndexType,I>::result <<"] = " << mw[AccessElement<RpvIndexType,I>::result] <<std::endl;
      //!in [KEE et al.]                                                    x[I]
      return ( x[K]/(mw[AccessElement<RpvIndexType,I>::result]*d[SymmetricDenseMatrixAccess<I,K,N>::offset()])*(mw[AccessElement<RpvIndexType,J>::result]*x[J]*(1.-Kronecker<I,K>::delta()) - mw[AccessElement<RpvIndexType,I>::result]*x[J]*(Kronecker<I,J>::delta() - Kronecker<J,K>::delta())) + InnerLoopType::template inner_loop_L0000_reduced<TPTDATATYPE>(x,mw,d) );
    }

    
    template<class X, class MW, class PITER>
    static inline typename X::value_type inner_loop_L0000(const X& x, const MW& mw, PITER d){
      return ( x[K]/(mw[I]*d[SymmetricDenseMatrixAccess<I,K,N>::offset()])*(mw[J]*x[J]*(1.-Kronecker<I,K>::delta()) - mw[I]*x[J]*(Kronecker<I,J>::delta() - Kronecker<J,K>::delta())) + InnerLoopType::inner_loop_L0000(x,mw,d) ); 
    }

    template<class TPTDATATYPE, class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void medium_loop_L0000_reduced(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){
      //! \f$ L^{00,00}\f$ is assumed to be stored row-wisely
      smart_assign(Diff[MatrixStorageFormat4LinearMemory<I,J,N,M,'R'>::offset()],16./25.*(temp/(ChooseRightPrefactor<PRESS,TransportSettings::pressureInBar>::in_bar(p)))*ThreeNestedLoopsUnroller<I,J,0,M,N,R>::template inner_loop_L0000_reduced<TPTDATATYPE>(x,mw,d)); //K = 0 to R
      MediumLoopType::template medium_loop_L0000_reduced<TPTDATATYPE>(Diff,p,temp,x,mw,d);
    }

    
     template<class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void medium_loop_L0000(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){
      //! \f$ L^{00,00}\f$ is assumed to be stored row-wisely
      smart_assign(Diff[MatrixStorageFormat4LinearMemory<I,J,N,M,'R'>::offset()],16./25.*(temp/(ChooseRightPrefactor<PRESS,TransportSettings::pressureInBar>::in_bar(p)))*ThreeNestedLoopsUnroller<I,J,0,M,N,R>::inner_loop_L0000(x,mw,d)); //K = 0 to R
      MediumLoopType::medium_loop_L0000(Diff,p,temp,x,mw,d);
    }
    
    template<class TPTDATATYPE, class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void create_L0000_reduced(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){
      ThreeNestedLoopsUnroller<I,0,K,N,M,R>::template medium_loop_L0000_reduced<TPTDATATYPE>(Diff,p,temp,x,mw,d);  //!from J = 0 to M. Unlike the binary diffusion matrix, L0000 is non-symmetric! 
      OuterLoopType::template create_L0000_reduced<TPTDATATYPE>(Diff,p,temp,x,mw,d);
    }

    template<class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void create_L0000(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){
      ThreeNestedLoopsUnroller<I,0,K,N,M,R>::medium_loop_L0000(Diff,p,temp,x,mw,d);  //!from J = 0 to M. Unlike the binary diffusion matrix, L0000 is non-symmetric! 
      OuterLoopType::create_L0000(Diff,p,temp,x,mw,d);
    }

    //LLLLLLLLLLLLLLLLLLLL LINEAR ALGEBRA STUFF LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    //················ MATRIX-MATRIX PRODUCT: c = a*b ·························
    //matrix - matrix multiplication
    //TUNED FLPT ARITHMETIC CAN USED HERE - NOT DONE SO FAR ;)
    template<class A, class B, class C>
    static inline void inner_loop_matrix_matrix_prod(C& c, const A& a, const B& b){
     
      c[I][J] += a[I][K]*b[K][J];
      InnerLoopType::inner_loop_matrix_matrix_prod(c,a,b);
    }
    
    template<class A, class B, class C>
    static inline void medium_loop_matrix_matrix_prod(C& c, const A& a, const B& b){
      //!perform one complete inner loop
      ThreeNestedLoopsUnroller<I,J,0,N,M,R>::inner_loop_matrix_matrix_prod(c,a,b);
      MediumLoopType::medium_loop_matrix_matrix_prod(c,a,b);
    }

    //outer loop 
    /**
     * Matrix-matrix multiplication \f$c = a\cdot b\f$.
     * NOTE: as in the conventional for loop version, the second triple of the
     * template parameters reads as <·,·,·,a.rows(),b.columns(),a.columns()>.
     * More precise, if \f$a \in R^{m \times n}, b \in \in R^{n \times p}\f$ and
     * \f$ c \in R^{m \times p}\f$, then the multiplication can be started via
     * \code
       TinyMatrix<double,m,n> a;
       TinyMatrix<double,n,p> b;
       TinyMatrix<double,m,p> c(0.);  //MUST be initialized since we use +=
       ThreeNestedLoopsUnroller<0,0,0,m,p,n>::matrix_matrix_prod(c,a,b);
     *\endcode
     */
    template<class A, class B, class C>
    static inline C& matrix_matrix_prod(C& c, const A& a, const B& b){
      //! complete medium loop
      ThreeNestedLoopsUnroller<I,0,K,N,M,R>::medium_loop_matrix_matrix_prod(c,a,b);
      OuterLoopType::matrix_matrix_prod(c,a,b);
      
      return c;
    }


    /////
    template<class OP, class AIT, class BIT>
    static inline typename OP::value_type inner_loop_mat_op_mat_mat(const AIT& a, const BIT& b){
      return a[MatrixStorageFormat4LinearMemory<I,K,N,R,'r'>::offset()]*b[MatrixStorageFormat4LinearMemory<K,J,R,M,'r'>::offset()] + InnerLoopType::template inner_loop_mat_op_mat_mat<OP>(a,b);
    }
    

    template <class OP, class AIT, class BIT, class CIT>
    static inline void medium_loop_mat_op_mat_mat(CIT c, const AIT& a, const BIT& b){
      //! test matrix-matrix product
      //c[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()] = ThreeNestedLoopsUnroller<I,J,0,N,M,R>::template inner_loop_mat_op_mat_mat<OP>(a,b);
      //std::cout << c[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()] << " ";
     
      OP::apply(c[MatrixStorageFormat4LinearMemory<I,J,N,M,'r'>::offset()], ThreeNestedLoopsUnroller<I,J,0,N,M,R>::template inner_loop_mat_op_mat_mat<OP>(a,b));
       
      MediumLoopType::template medium_loop_mat_op_mat_mat<OP>(c,a,b);
    }


    template <class OP, class AIT, class BIT, class CIT>
    static inline void mat_op_mat_mat(CIT c, const AIT& a, const BIT& b){
      ThreeNestedLoopsUnroller<I,0,K,N,M,R>:: template medium_loop_mat_op_mat_mat<OP>(c,a,b);
      OuterLoopType::template mat_op_mat_mat<OP>(c,a,b);
    }

    //··········································································
  };


  //· END OF RECURSION · END OF RECURSION · END OF RECURSION · END OF RECURSION
  //· END OF RECURSION · END OF RECURSION · END OF RECURSION · END OF RECURSION 
  //· END OF RECURSION · END OF RECURSION · END OF RECURSION · END OF RECURSION

  //end INNER loop
  template<int I, int J, int N, int M, int R>
  class ThreeNestedLoopsUnroller<I,J,R,N,M,R>{
  public:
    template<class TPTDATATYPE, class X, class MW, class PITER>
    static inline typename X::value_type inner_loop_L0000_reduced(const X& x, const MW& mw, PITER d){
      return typename X::value_type();
    }
  
    
    template<class X, class MW, class PITER>
    static inline typename X::value_type inner_loop_L0000(const X& x, const MW& mw, PITER d){
      return typename X::value_type();
    }

     //················ MATRIX-MATRIX PRODUCT ··································
    template<class A, class B, class C>
    static inline void inner_loop_matrix_matrix_prod(C& c, const A& a, const B& b){} //do nothing
    //··········································································
  
    template<class OP, class AIT, class BIT>
    static inline typename OP::value_type inner_loop_mat_op_mat_mat(const AIT& a, const BIT& b){
      typedef typename OP::value_type value_type;
      return value_type();
    }

    
  };
  
  
 

  //end MEDIUM loop
  template<int I, int K, int N, int M, int R>
  class ThreeNestedLoopsUnroller<I,M,K,N,M,R>{
  public:
    template<class TPTDATATYPE, class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void medium_loop_L0000_reduced(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){} //do nothing

    template<class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void medium_loop_L0000(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){} //do nothing

    //················ MATRIX-MATRIX PRODUCT ··································
    template<class A, class B, class C>
    static inline void medium_loop_matrix_matrix_prod(C& c, const A& a, const B& b){} //do nothing
    //··········································································
  

    template <class OP, class AIT, class BIT, class CIT>
    static inline void medium_loop_mat_op_mat_mat(CIT c, const AIT& a, const BIT& b){}
  };

  //end OUTER loop
  template<int J, int K, int N, int M, int R>
  class ThreeNestedLoopsUnroller<N,J,K,N,M,R>{
  public:
    template<class TPTDATATYPE, class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void create_L0000_reduced(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){} //do nothing
    
    template<class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void create_L0000(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){} //do nothing
    
    

    //················ MATRIX-MATRIX PRODUCT ··································
    template<class A, class B, class C>  //do nothing
    static inline C& matrix_matrix_prod(C& c, const A& a, const B& b){
      return c;
    }
    //··········································································
  
    template <class OP, class AIT, class BIT, class CIT>
    static inline void mat_op_mat_mat(CIT c, const AIT& a, const BIT& b){}
  };

  //OVERALL STOP
  template<int N, int M, int R>
  class ThreeNestedLoopsUnroller<N,M,R,N,M,R>{
  public:
    template<class TPTDATATYPE, class X, class MW, class PITER>
    static inline typename X::value_type inner_loop_L0000_reduced(const X& x, const MW& mw, PITER d){
      return typename X::value_type();
    }

     template<class X, class MW, class PITER>
    static inline typename X::value_type inner_loop_L0000(const X& x, const MW& mw, PITER d){
      return typename X::value_type();
    }

    template<class TPTDATATYPE, class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void medium_loop_L0000_reduced(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){} //do nothing
    
    template<class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void medium_loop_L0000(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){} //do nothing
    
    
    template<class TPTDATATYPE, class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void create_L0000_reduced(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){} //do nothing

     template<class MCITER, class PRESS, class TEMP, class X, class MW, class PITER>
    static inline void create_L0000(MCITER Diff, const PRESS& p, const TEMP& temp, const X& x, const MW& mw, PITER d){} //do nothing



    //················ MATRIX-MATRIX PRODUCT ··································
    template<class A, class B, class C>
    static inline void inner_loop_matrix_matrix_prod(C& c, const A& a, const B& b){} //do nothing

    template<class A, class B, class C>
    static inline void medium_loop_matrix_matrix_prod(C& c, const A& a, const B& b){} //do nothing

    template<class A, class B, class C>  //do nothing
    static inline C& matrix_matrix_prod(C& c, const A& a, const B& b){
      return c;
    }
    //··········································································
    
    template<class OP, class AIT, class BIT>
    static inline typename OP::value_type inner_loop_mat_op_mat_mat(const AIT& a, const BIT& b){
      typedef typename OP::value_type value_type;
      return value_type();
    }
    
    template <class OP, class AIT, class BIT, class CIT>
    static inline void medium_loop_mat_op_mat_mat(CIT c, const AIT& a, const BIT& b){}

    template <class OP, class AIT, class BIT, class CIT>
    static inline void mat_op_mat_mat(CIT c, const AIT& a, const BIT& b){}

  };

} //end namespace

#endif
