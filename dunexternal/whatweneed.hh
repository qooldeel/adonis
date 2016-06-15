#ifndef THINGS_WE_NEED_FOR_CALCULATING_THE_REDUCED_EVOLUTION_EQUATIONS_HH
#define THINGS_WE_NEED_FOR_CALCULATING_THE_REDUCED_EVOLUTION_EQUATIONS_HH


/** @brief DESCRIPTION:
 *
 *  In the following I refer to Ren/Pope, "Combustion and Flame 147, 2006"
 *
 * The evolution equation of the reduced composition r := B^T·z is given via 
 * eq. (44).
 *
 * Important quantities to describe eq. (44):
 * 
 *     z_M(r) := z = \phi(r) \in \mathbb{R}^{n_s} describes the full composition
 *            w.r.t. the reduced variable vector r \in \mathbb{R}^{n_r}.
 *
 *     T(r) \in \mathbb{R}^{n_s \times n_r} denotes the tangent space.
 *
 *     N(r) \in \mathbb{R}^{n_s \times n_u} denotes the normal space.
 *
 *     
 *  \tparam D  object to typify FieldMatrices and FieldVector.
 *  \tparam NFULL n_s, i.e. the dimension of the full system.
 *  \tparam NRED n_r, i.e. # of rpv's and dimension of reduced system, respectively
 *            (n_u := n_s - n_r is created implicitly). 
 *  \tparam B select either usual or transposed representation of normal space
 */


#include "../dunestuff/fmatrix.hh"
#include "../templatemetaprograms/conditionalstructures.hh"

#include "axpy.hh"
#include<cassert>

namespace Adonis{

  /**
   *  \brief This fancy template meta programme allows you, e.g. to overload a function with different return types due to the boolean of the COMPILE TIME expression 'B'
   * \tparam T1 return type 
   * \tparam T2 another return type
   * \tparam B boolean which if 'true' takes T1 as return type else T2
   */
  template<class T1, class T2, bool B>
  class FieldMatrixTypeSelector{};

  template<class T1, class T2>
  class FieldMatrixTypeSelector<T1,T2,true>{
  public:
    typedef T1 MatrixType4NormalSpace;
  };
  
  template<class T1, class T2>
  class FieldMatrixTypeSelector<T1,T2,false>{
  public:
    typedef T2 MatrixType4NormalSpace;
  };


  /**
   * \brief Base class for essential using the Barton Nackman trick, cf.
   *   [BARTON/NACKMAN, "Scientific and Engineering C++, 1994"]
   */
  template<class K, int N, int R, class ESS>
  class BaseEssentialInheritance{
  private:
    inline const ESS& R2D() const{       //refer to derived class 
      return static_cast<const ESS&>(*this);
    }

    
  protected:
    typedef Dune::FieldVector<K,N> CompositionType;
    typedef Dune::FieldMatrix<K,N,R> TangentSpaceType;
    
    CompositionType zM_;
    TangentSpaceType T_;

  public:
    inline int full_dim() const{return N;}      //common functionalities
    inline int red_dim() const{return R;}

    const CompositionType& get_composition() const {return zM_;}
    CompositionType& get_composition() {return zM_;}
    
    const TangentSpaceType& get_tangentspace() const {return T_;}
    TangentSpaceType& get_tangentspace(){return T_;}

  };



  //! if <TT> B = false </TT> then \f$ N^T \f$ will be used 
  template<class D, int NFULL, int NRED, bool B = true>
  class Essential: public BaseEssentialInheritance<D,NFULL,NRED,Essential<D,NFULL,NRED,B> >{
  public:

    typedef BaseEssentialInheritance<D,NFULL,NRED,Essential<D,NFULL,NRED,B> > Basis;

    typedef typename Basis::CompositionType CompositionType;
    typedef typename Basis::TangentSpaceType TangentSpaceType;
     
    typedef Dune::FieldMatrix<D, NFULL, NFULL-NRED> NormalSpaceType;
    typedef Dune::FieldMatrix<D, NFULL-NRED, NFULL> TransposedNormalSpaceType;
    typedef Dune::FieldVector<D, NRED>  ReducedType;
  
    //! choose matrix type of normal space according to boolean B 
    typedef typename FieldMatrixTypeSelector<NormalSpaceType,TransposedNormalSpaceType,B>::MatrixType4NormalSpace MatrixType4NormalSpace; 

  private:
    
    MatrixType4NormalSpace N_;

    //NormalSpaceType N_;     //!Normal Space in C-style representation
   

  public:
    Essential(){ 
      IfElse<(NRED > NFULL)>::statement("# of reduced variables is greater or equal than # of full variables. \n   No reduction possible ;-)"); 
     
    }   //default construction

    //construct from the same value
    Essential(const D& d): N_(d){
     
      IfElse<(NRED > NFULL)>::statement("# of reduced variables is greater or equal than # of full variables. \n   No reduction possible ;-)"); 
      
      /**
       * \brief Initialise the Dune::FieldObjects here.
       *  NOTE: you cannot use them in the initialiser list, i.e.... :Basis::zM_(d)!. Nor can you use a constructor of the base class which initialises zM_ and T_!
       */
      Basis::zM_ = d; 
      Basis::T_ = d;
    }
    
    //initialise each field with a different value 
    Essential(const D& d1, const D& d2, const D& d3): N_(d3){
      IfElse<(NRED >= NFULL)>::statement("# of reduced variables is greater or equal than # of full variables. \n   No reduction possible ;-)"); 
    
      Basis::zM_ = d1; 
      Basis::T_ = d2;
    }
    
    //construct from CompositionType, Dune::FieldMatrix and NormalSpaceType
    Essential(const CompositionType& z, const Dune::FieldMatrix<D, NFULL, NRED>& tspace, const NormalSpaceType& nspace): N_(nspace){
      IfElse<(NRED >= NFULL)>::statement("# of reduced variables is greater or equal than # of full variables. \n   No reduction possible ;-)"); 
     
      Basis::zM_ = z;
      Basis::T_ = tspace;
    }

    //just contruct from fiven value
    Essential& operator=(const D& d){
      Basis::zM_ = d;
      Basis::T_ = d;
      N_ = d;
      
      return *this;
    }

    const MatrixType4NormalSpace& get_normalspace() const {return N_;}
    MatrixType4NormalSpace& get_normalspace(){return N_;}
    
    //const NormalSpaceType& get_normalspace() const {return N_;}
    //NormalSpaceType& get_normalspace(){return N_;}
   


    // const typename FieldMatrixTypeSelector<NormalSpaceType,TransposedNormalSpaceType,true>::MatrixType4NormalSpace& get_normalspace() const{ return N_;}
    


    //const TransposedNormalSpaceType& get_normalspace() const {return TransN_;}


    int full_dim() const {return NFULL;}

    int reduced_dim() const {return  NRED;}

    int unrepresented_dim() const {return NFULL-NRED;}
    


    //======= Copy constructor and copy assignment from Axpy to Essential ======

    Essential(const Axpy<D,NFULL,NRED,B>& A){
      
      for(int i = 0; i < NFULL; ++i){
	Basis::get_composition()[i] =  A.a()*A.x().Basis::get_composition()[i] + A.y().Basis::get_composition()[i];
      
	/*for(int j = 0; j < NRED; j++)
	  Basis::get_tangentspace()[i][j] = A.x().Basis::get_tangentspace()[i][j];

	for(int k = 0; k < NFULL-NRED; ++k)
	  get_normalspace()[i][k] = A.x().get_normalspace()[i][k];
	*/  
      }
      //does the same job, except those FieldMatrices need some update
       Basis::get_tangentspace() =  A.x().Basis::get_tangentspace();
       get_normalspace() = A.x().get_normalspace();
    }

    //basically the same as the copy constructor except that we have a return
    //value
    Essential& operator=(const Axpy<D,NFULL,NRED,B>& A){
     
      for(int i = 0; i < NFULL; ++i){
	Basis::get_composition()[i] =  A.a()*A.x().Basis::get_composition()[i] + A.y().Basis::get_composition()[i];
      
	/*for(int j = 0; j < NRED; j++)
	  Basis::get_tangentspace()[i][j] = A.x().Basis::get_tangentspace()[i][j];

	for(int k = 0; k < NFULL-NRED; ++k)
	  get_normalspace()[i][k] = A.x().get_normalspace()[i][k];
	*/  
      }	
      //does the same job, except those FieldMatrices need some update
      Basis::get_tangentspace() =  A.x().Basis::get_tangentspace();
      get_normalspace() = A.x().get_normalspace();

      return *this;
    }


    //==========================================================================



    friend std::ostream& operator<<(std::ostream& os, const Essential<D,NFULL,NRED,B>& ess){
      os << " # Manifold point z_M, tangent space (n_s x n_r) and normal space (n_s x n_u,   n_u := n_s - n_r), respectively :" << std::endl << std::endl;
      os << ess.Basis::zM_ << std::endl<< std::endl; 
      os << "   " << ess.Basis::T_ << std::endl<< std::endl;
      os << "   " << ess.N_ << std::endl<< std::endl;
      return os;
    }
  

    //define some operators on essential such +=,-=, *=...
    Essential<D,NFULL,NRED,B>& operator+=(const Essential<D,NFULL,NRED,B>& x){
      Basis::zM_ += x.Basis::zM_;
      //Basis::T_ += x.Basis::T_;         //is this and the next line desireable ?
      //N_ += x.N_;
      return *this;
    }
    
    Essential<D,NFULL,NRED,B>& operator-=(const Essential<D,NFULL,NRED,B>& x){
      Basis::zM_ -= x.Basis::zM_;
      //Basis::T_ -= x.Basis::T_;         //is this and the next line desireable ?
      //N_ -= x.N_;
      return *this;
    }
    
    
    Essential<D,NFULL,NRED,B>& operator*=(const D& s){
      Basis::zM_ *= s;
      //Basis::T_ *= s;         //is this and the next line desireable ?
      //N_ *= s;
      return *this;
    }
    
    
    Essential<D,NFULL,NRED,B>& operator/=(const D& s){
      if(s == D())
	 ADONIS_ERROR(ZeroDivision, "Divison by Zero (s ~ "<<s<<").");
	 
       Basis::zM_ /= s;
      //Basis::T_ /= s;         //is this and the next line desireable ?
      //N_ /= s;
      return *this;
    }

    
    /*Essential<D,NFULL,NRED,B>& axpy(const D& a, const Essential<D,NFULL,NRED,B>& x){
      Basis::get_composition() += x.Basis::get_composition();
      
      return *this;
    }
    */

 };




}//end namespace 


#endif
