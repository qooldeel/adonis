#ifndef STOICHIOMETRIC_COEFFICIENTS_SUPPLIED_FROM_SOMEWHERE_HH
#define STOICHIOMETRIC_COEFFICIENTS_SUPPLIED_FROM_SOMEWHERE_HH

#include <iostream>
#include "../common/adonisassert.hh"

namespace Adonis{

  /**
   * \brief The stoichiometric matrix only contains forward stoichiometric
   * coefficients. Each <--> reaction is considered as two ones, namely --> and <--. When a reaction consists only of a --> reaction, then the LHS of the  <-- reaction consists of zero stoichiometric coeffs only (since no products are formed), thus leading to a (# of species) x (2· # of --> reactions) matrix
   *
   * Note that the stoichiometric matrix is symmetric. Hence it is obvious to store only, say, the \$\nu'\f$, since the \f$\nu''\$ can be easily queried.
    */
  template<class IT1>
  class StoichiometricMatrix{
  public:
    
    StoichiometricMatrix(size_t nspec = 0, size_t nreac = 0, const IT1& nu = IT1()):numberOfSpecies_(nspec),numberOfReactions_(nreac),nu_(nu){}

    size_t number_of_species() const {return numberOfSpecies_;}

    //! number of forward reactions
    size_t number_of_forward_reactions() const {return numberOfReactions_;}
    

    //! species k and reaction index ( index = 1, ..., 2·I)
    size_t nu_forward(size_t k, size_t index) const{
      adonis_assert(k<numberOfSpecies_);
      adonis_assert(index < 2*numberOfReactions_);
      return static_cast<size_t>(nu_[k + index*numberOfSpecies_]);
    }
    
    //! if index \f$i\f$ is even, then \f$ \nu''_{k,i} = \nu'_{k,i+1}.\f$ 
    //! If index \f$i \f$ is odd, then \f$ \nu''_{k,i} = \nu'_{k,i-1}.\f$
    size_t nu_reverse(size_t k, size_t index) const{
      adonis_assert(k<numberOfSpecies_);
      adonis_assert(index < 2*numberOfReactions_);
      return ( (index%2 == 0) ? nu_forward(k,index+1) : nu_forward(k,index-1) );
    }
    
    //! note that the net stoichiometric coefficent can be negative as well;
    //! Hence we must use a signed integer here 
    int nu_net(size_t k, size_t index) const{
      adonis_assert(k<numberOfSpecies_);
      adonis_assert(index < 2*numberOfReactions_);
      
      return ( static_cast<int>(nu_reverse(k,index)) - static_cast<int>(nu_forward(k,index)) );
    }

    //! calculate molecularity (without 3rd bodies) of reaction \f$i\f$
    size_t molecularity_of_reaction(size_t i) const{
      adonis_assert(i < 2*numberOfReactions_);

      size_t ix =  i*numberOfSpecies_,
	molec = 0;
      for(size_t k = 0; k <  numberOfSpecies_; ++k){
	if(nu_[ix] != 0)
	  molec += nu_[ix];
	ix++;
	if(molec == 3) // o.k. more than trimolecular reactions don't exist
	  break;       // and we can leave the loop early
      }
      return molec;
    }

    friend std::ostream& operator<<(std::ostream& os, const StoichiometricMatrix& SM){
      os << "Stoichiometric matrix nu', nu'' {rows: reactions (even row index: -->, odd: <--), columns: species }"<< std::endl;
     
      for(size_t i = 0; i < 2*SM.numberOfReactions_; ++i){ //!2 times !!
	for(size_t k = 0; k < SM.numberOfSpecies_; ++k){
	  os << SM.nu_forward(k,i) << "  "<< SM.nu_reverse(k,i)<< "  ";
       }
	os << std::endl;
     }
      return os;
    }

  private:
    size_t numberOfSpecies_, 
    numberOfReactions_;  //! consider only forward (-->) reactions
    IT1 nu_;                 //! store only \f$ \nu'_{ki}, \ i = 1,\ldots, 2*I\f$ 
  };

}

#endif
