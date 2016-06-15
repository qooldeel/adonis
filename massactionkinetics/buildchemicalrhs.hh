#ifndef BUILD_CHEMICAL_COMBUSTION_SOURCE_TERM_HH
#define BUILD_CHEMICAL_COMBUSTION_SOURCE_TERM_HH

#include <iostream>
#include <cmath>

#include "../common/adonisassert.hh"
#include "../accuracy/floatingpointarithmetic.hh"
#include "../common/globalfunctions.hh"
#include "../common/error.hh"
#include "indexinjector.hh"

namespace Adonis{

  /**
   * \brief Creates the chemical source term.
   *
   * Note: allowing dimensions known at compile time could make the code much faster...
   */
  template<class FORWRATE, class REVRATE, class INDEXER>
  class BuildChemicalSourceTerm{
  public:
    typedef FORWRATE* FwdPointerType;
    typedef REVRATE* RevPointerType; 
    
    typedef INDEXER* IndexerPointerType;

    typedef typename REVRATE::TemperatureBoundsIterType TemperatureBoundsIterType;
    typedef bool* BoolPointerType;

  private:
    
    FwdPointerType fwdRatesPtr_;
    RevPointerType revRatesPtr_;
    IndexerPointerType ixPtr_;
    BoolPointerType expRevCoeffsPtr_;

    //! if we do not care about given reverse rate coeffs then just return false,
    //! and everything is computed the usual way. Otherwise, use data given by the 
    //! specified bool pointer
    bool exp_rev_coeff_access(size_t i) const{
      if(expRevCoeffsPtr_ == 0)
	return false;
      else
	return expRevCoeffsPtr_[i];
    }

  public:
    typedef typename REVRATE::StoichPointerType StoichPointerType;

#if USE_TUNED_FLOATING_POINT_ARITHMETIC 
    typedef typename REVRATE::AdditionType AdditionType;
#endif

    //! also default
    BuildChemicalSourceTerm(FwdPointerType fw = 0, RevPointerType rv = 0,IndexerPointerType ip = 0, BoolPointerType exrev = 0):fwdRatesPtr_(fw),revRatesPtr_(rv),ixPtr_(ip),expRevCoeffsPtr_(exrev){}  
    
    //! if default constructor is used then use this member to initialize object
    void initialize(FwdPointerType fwd, RevPointerType rev, IndexerPointerType ip, BoolPointerType exrev = 0){
      if((fwdRatesPtr_ == 0) && (revRatesPtr_ == 0)){
	fwdRatesPtr_ = fwd;
	revRatesPtr_ = rev;
	ixPtr_ = ip;
	expRevCoeffsPtr_ = exrev;
      }
    }

    
    size_t number_of_species() const{return ixPtr_->number_of_species();}

    size_t species_index(size_t k) const {return ixPtr_->s_ix(k);}

    size_t reaction_index(size_t i) const {return ixPtr_->r_ix(i);}
    
    template<class T> 
    T C_p(size_t k, const T& temp) const{
#ifndef NDEBUG     
      if(temp < revRatesPtr_->bounds_on_T()[3*k] || temp > revRatesPtr_->bounds_on_T()[3*k+1])
	ADONIS_ERROR(BoundsError, "T = "<<temp << " left its bounds.");
#endif
      return revRatesPtr_->C_p(k,temp);
    }
   
    template<class T> 
    T H_T(size_t k, const T& temp) const{
      return revRatesPtr_->H_T(k,temp);
    }

    template<class T> 
    T S_T(size_t k, const T& temp) const{
      return revRatesPtr_->S_T(k,temp);
    }

    //! define other thermodynamic quantites 
    template<class ITY, class ITW, class TEMP>
    typename ITY::value_type cp(size_t nspec, const TEMP& temp, ITY y, ITW w){
      typedef typename ITY::value_type value_type;
      value_type cp = value_type();
      for(size_t k = 0; k < nspec; ++k){
	cp += (*this).C_p(k,temp)/w[k]*y[k];
      }
#ifndef NDEBUG
      if(!is_well_defined_value(cp))
	ADONIS_ERROR(ValueError,"Bad cp = "<< cp <<".");
#endif
      return cp;
    }


    /** \brief Prefactor for rate-of-progress expression \f$ q_i\f$ (see below) 
     *	for reaction \f$i\f$ at concentration \f$u\f$
     *   in case of 3rd body occurance: \f$ \sum_{k=1}^K\alpha_{k,i}[X_k]\f$
     *   else: 1.
     */
    template<class U>  //! U is concentration here
    typename U::value_type effective_3rd_body_total_mixture_concentration(size_t i, const U& u) const{
      
      size_t nfreac = fwdRatesPtr_->number_of_forward_reactions(),
	nos = fwdRatesPtr_->number_of_species();
      
      adonis_assert(i < nfreac);

      typedef typename U::value_type value_type;
      
      if(fwdRatesPtr_->get_troe_data()[i+nfreac] != -1){  //found 3rd body in reaction i
	value_type tb = value_type();
	
#if USE_TUNED_FLOATING_POINT_ARITHMETIC 
	value_type s = value_type(),
	  c = value_type();
#endif	

	for(size_t k = 0; k < nos; ++k){
	  
#if USE_TUNED_FLOATING_POINT_ARITHMETIC 
	  tb = AdditionType::add(fwdRatesPtr_->get_collision_efficiencies()[static_cast<int>(fwdRatesPtr_->get_troe_data()[i+nfreac])*nos + k]*u[k],s,c);
	  

#else
	  tb += fwdRatesPtr_->get_collision_efficiencies()[static_cast<int>(fwdRatesPtr_->get_troe_data()[i+nfreac])*nos + k]*u[k];	
#endif
	}
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	CorrectRoundingErrorAfterwards<value_type,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
	if(c != value_type()){  //! if c == 0 then do not assign s since 
	                 //! it might be 0 because no addition was involved!
	  CorrectRoundingErrorAfterwards<value_type,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(tb,s);
	}
#endif
	//std::cout << "3rd body contribution [M] = "<< tb << std::endl;
	//adonis_assert(tb != value_type());
#ifndef NDEBUG
      if(!is_well_defined_value(tb))
	ADONIS_ERROR(ValueError,"Bad tb = "<< tb <<".");
#endif
	return tb;
      }
      else 
	return 1;
      
    }

    /** The rate of progress for reaction \f$i\f$ is given by  \f[ q_i = \sum_{k=1}^K \alpha_{k,i}[X_k]\left(k_{f,i}\prod_{k=1}^K[X_k]^{\nu'_{k,i}} - k_{r,i}\prod_{k=1}^K[X_k]^{\nu''_{k,i}} \right),\f] where the first term is 1 as soon as <I> no </I> 3rd body occurs in reaction \f$i.\f$
     *
     * REMARK: in an unchemical way, one can say that \f$q_i = \frac{d[X_{ki}]}{dt}\f$  
     */
    template<class U, class TEMP> //! U is concentration here
    TEMP rate_of_progress(size_t i, const TEMP& temp, const U& u, bool rev = false) const{
      //typedef typename U::value_type value_type;
      
      size_t nos = revRatesPtr_->number_of_species();

      //! note that the nos-th entry corresponds to temperature
      TEMP k_fi = fwdRatesPtr_->k_forward(i,temp,u),
	//if rev is true, then calculate reverse rates by explicit coeffs
	//! note: the forward rates need not to be argumented by 'rev'
	k_ri = revRatesPtr_->k_reverse(i,temp,u,rev); 

      //! test me
      //std::cout << "k_fi = "<< k_fi << "   k_ri = "<< k_ri << std::endl;

      TEMP P1 = 1, P2 = 1;

      size_t nuF, nuR;  //nu' and nu'', respectively

      
      for(size_t k = 0; k < nos; ++k){
	//! ATTENTION: you need the --> stoichiometric coeffs here, i.e. 2*i
	nuF = revRatesPtr_->nu_forward(k,2*i);
	nuR = revRatesPtr_->nu_reverse(k,2*i);

	if(nuF != 0)
	  P1 *= natural_power(u[k],nuF); //pow(u[k],static_cast<TEMP>(nuF));
	if(nuR != 0)
	  P2 *= natural_power(u[k],nuR); //pow(u[k],static_cast<TEMP>(nuR));
      }
#ifndef NDEBUG
      if((!is_well_defined_value(P1)) || (!is_well_defined_value(P2)))
	ADONIS_ERROR(ValueError,"Bad: P1 = "<< P1 <<",  P2 = "<< P2 << ".");
      else if ((!is_well_defined_value(k_fi)) || (!is_well_defined_value(k_ri)))
	ADONIS_ERROR(ValueError,"Bad k_fi = "<< k_fi << " or k_ri = "<<k_ri << ".");
      else {}
#endif
      return ( (*this).effective_3rd_body_total_mixture_concentration(i,u)*(k_fi*P1 - k_ri*P2) );

    }
    

    /**
     * \brief Calculation of the net production or destruction rate (in mol/m\f$^3 \cdot\f$s) of species \f$k\f$ which is given as \f$\dot{\omega}_k = \dot{[X_k]} = \sum_{i=1}^I \nu_{ki}q_i, \f$ where \f$q_i\f$ is the rate of progress (cf. above) and \f$ \nu_{ki} = \nu''_{ki} - \nu'_{ki}\f$ is the net stoichiometric coefficient. 
     */
    template<class U, class TEMP> //! U is concentration here!!
    TEMP net_formation_rate(size_t k, const TEMP& temp, const U& u) const{
      //typedef typename U::value_type value_type;

      TEMP wkprime = TEMP();
      
      size_t nor = ixPtr_->number_of_forward_reactions();//revRatesPtr_->number_of_forward_reactions();
      
      int nu_ki;

#if USE_TUNED_FLOATING_POINT_ARITHMETIC 
	TEMP s = TEMP(),
	  c = TEMP();
#endif	


      for(size_t i = 0; i < nor; ++i){
	
	//! ATTENTION: you need the --> stoichiometric coeffs here, i.e. 2*i
	//! NOTE: the reaction index only must be supplied for bidirectional 
	//!       reactions, i.e. as they occur in the mechanism file
	nu_ki = revRatesPtr_->nu_net(k,2*ixPtr_->r_ix(i));	
	
	if(nu_ki != 0){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC 
	  wkprime = AdditionType::add(nu_ki*(*this).rate_of_progress(ixPtr_->r_ix(i),temp,u,(*this).exp_rev_coeff_access(i)),s,c);
#else
	  wkprime += nu_ki*(*this).rate_of_progress(ixPtr_->r_ix(i),temp,u,(*this).exp_rev_coeff_access(i));
#endif
	}	
      }

      
#if USE_TUNED_FLOATING_POINT_ARITHMETIC      
      CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
      if(c != TEMP()){  //! if c == 0 then do not assign s since 
	                 //! it might be 0 because no addition was involved!
	CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(wkprime,s);
      }
#endif

#ifndef NDEBUG
      if(!is_well_defined_value(wkprime))
	ADONIS_ERROR(ValueError,"Bad wkprime = "<< wkprime <<"; temperature = "<< temp << ".");
#endif
    
      return wkprime;
    }

};


} //end namespace

#endif
