#ifndef THE_MATRICES_B_T_and_U_NEEDED_4_COMPUTING_H_T_HH
#define THE_MATRICES_B_T_and_U_NEEDED_4_COMPUTING_H_T_HH

#include "../dunestuff/fmatrix.hh"

namespace Adonis{
  
  /**
   * \brief Create matrices \f$ B^T \f$ and \f$ U\f$ representing the representing the reaction progress variables and the inrepresented ones, respectively.
   *
   * \code
   //dimenson of full PDE system: 5, reduced system: 3
   
   std::vector<int> ix(3);  //Dune::FieldMatrix<int,3> ix; //also working
   ix[0] = 1; ix[1] = 3; ix[2] = 4;
   
   std::vector<int> unrep(2); 
   unrep[0] = 0; unrep[1] = 2;

   Dune::FieldMatrix<double,3,5> B_T;
   Dune::FieldMatrix<double,5,2> U;

   ReactionProgressVariables<double,5,3> RVP(B_T,U);
   RVP.create_B_T(ix);
   RVP.create_U(unrep);
   
   std::cout << B_T << std::endl;
   std::cout << U << std::endl;
   
   * \endcode
   */
  template<class T, int FULL, int RED>
  class ReactionProgressVariables{
  public:
    typedef Dune::FieldMatrix<T,RED,FULL> BTransposedType;
    typedef Dune::FieldMatrix<T,FULL,FULL-RED> UType;
    typedef Dune::FieldVector<int,RED> IndexBType;
    typedef Dune::FieldVector<int,FULL-RED> IndexUType;

    ReactionProgressVariables(BTransposedType& B_T, UType& U):B_T_(B_T), U_(U){}

    /**
     * \brief Create \f$ B^T \f$ via index vector.
     *
     * \param indices The indices of the reaction progress variables within the full composition vector \f$ z = z^{\mathcal{M}} \in \mathbb{R}^{n_{\operatorname{s}}}\f$.
     *
     * Example: \f$ z^{\mathcal{M}} = [z_1,z_2,z_3,z_4,z_5]^T\f$. Then the integer vector \f$ [1,3,4] \f$ suggests that the components \f$ z_1, z_3\f$ and \f$ z_4\f$ represent the variables in which the reduced system is to be described. 
     *
     * NOTE: make sure to arrange the indices in an appropriate order. Otherwise \f$ B^T\f$ swaps the equations while \f$ U \f$ remains always unchanged. Therefore it is recommended to store the rvp indices in increasing order. Example: 1,3,4 (not 3,1,4. This would lead to a swapped representation of the original equations while \f$ U\f$ keeps its "normal" appearance). 
     */
    inline void create_B_T(const IndexBType& indices){ //works also for iters
      for(int i = 0; i < RED; ++i){
	for(int j = 0; j < FULL; ++j){
	  if(j == indices[i])
	    B_T_[i][j] = 1.;
	  else 
	    B_T_[i][j] = 0.;
	}
      }
    }


    /**
     * \brief for STL-compliant random access containers
     */
    template<class VEC>
    inline void create_B_T(const VEC& indices){
      for(int i = 0; i < RED; ++i){
	for(int j = 0; j < FULL; ++j){
	  if(j == static_cast<int>(indices[i]))
	    B_T_[i][j] = 1.;
	  else 
	    B_T_[i][j] = 0.;
	}
      }
    }
    

    /**
     * \brief Here you can define \f$ B^T\f$ by hand.
     */
    inline void create_B_T(){
      B_T_ = 0.;
      //B_T_[0][0] = 1.; //define the represented species 
      //B_T_[0][1] = 1.;
      //....
    }


    /**
     * \brief Create \f$ U \f$ via index vector.
     *
     * \param indices The indices of the reaction progress variables within the full composition vector \f$ z = z^{\mathcal{M}} \in \mathbb{R}^{n_{\operatorname{s}}}\f$
     *
     * NOTE: make sure to arrange the indices in an appropriate order. Otherwise \f$ U f$ swaps the equations. See above.  
     */
    inline void create_U(const IndexUType& indices){
      const int nu = FULL-RED;
      for(int i = 0; i < FULL; ++i){
	for(int j = 0; j < nu; ++j){
	  if(i == indices[j])
	    U_[i][j] = 1.;
	  else 
	    U_[i][j] = 0.;
	}
      }
    }

    /**
     * \brief for STL-compliant random access containers
     */
    template<class VEC>
    inline void create_U(const VEC& indices){
      adonis_assert((int)indices.size() == FULL-RED);
      const int nu = FULL-RED;
      for(int i = 0; i < FULL; ++i){
	for(int j = 0; j < nu; ++j){
	  if(i == static_cast<int>(indices[j]))
	    U_[i][j] = 1.;
	  else 
	    U_[i][j] = 0.;
	}
      }
    }


    /**
     * \brief Here you can define \f$ U \f$ by hand.
     */
    inline void create_U(){
      //U_[0][0] = ; 
      //U_[0][1] = ;
      //....
    }
    
    
  private:
    BTransposedType& B_T_;   //references
    UType& U_;

  };


}

#endif 
