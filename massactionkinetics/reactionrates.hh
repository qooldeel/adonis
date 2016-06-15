#ifndef REACTION_RATES_FOR_NON_ISOTHERMAL_REACTIONS_HH
#define REACTION_RATES_FOR_NON_ISOTHERMAL_REACTIONS_HH

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>


#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../common/typeselector.hh"
#include "../common/isclass.hh"
#include "../common/universalconstants.hh"
#include "../accuracy/floatingpointarithmetic.hh"
#include "../expressiontemplates/xsettings.hh"
#include "../misc/misctmps.hh"
#include "physicalconstants.hh"

#include "metricsystem.hh"
#include "auxiliary.hh"

namespace Adonis{
  

  /**
   * \brief Arrhenius law
   * Note: temperature is considered as part of the differential equation and
   * may therefore have a different type. The formula is given by 
   * \f[ k_{f,i}(T) := AT^{\beta}\exp(-E_a/(RT)).\]
   */
  template<class T1, class T2>
  inline T1 arrhenius_formula(const T1& t, const T2& A, const T2& beta, const T2& Ea){
#ifndef NDEBUG
    if((!is_well_defined_value(t)) || (t < T1()))
      ADONIS_ERROR(ValueError,"Temperature is corrupted: T = "<< t << ".");   
#endif
    //adonis_assert((t != T1()) && (Ea != T2())); //!avoid <TT> exp(Nan) </TT> 
    return A*pow(t,beta)*exp(-Ea/(PhysicalConstants<T2>::IdealGasConstant*t));
  }


  //============================================================================
  //====================== CONVERT TO SI UNITS =================================
  //============================================================================
  /**
   * \brief Convert Arrhenius data (\f$ A, E_a\f$) to corresponding units, e.g.
   * SI units (kg, m, sec, K, etc.)
   */
  template<bool A, bool E> class ConvertArrheniusDataToConsistentUnits;

  template<>
  class ConvertArrheniusDataToConsistentUnits<true,false>{ //!adjust A, but not Ea
  public:
    
    template<class T>
    static inline T prefactor_A(const T& A, size_t molecularity, bool istroe){
      T Anew = A;

      if(molecularity == 2){   //bimolecularity 
	Anew = A*1.e-6;
	if(istroe)            //for <B> high pressure </B> troe reactions
	  Anew = A;           //the 3rd bodies don't count. This corresponds
      }                       //to a molecularity-1  

      if(molecularity == 3){   //trimolecularity 
	Anew = A*1.e-12;
	if(istroe)             //see above
	  Anew = A*1.e-6;
      }
    
    return Anew;
    }
  
    template<class T>
    static inline const T& E_activation(const T& Ea) {return Ea;}
  

    //!for LOW pressure troe reactions the third body counts!
    template<class T>
    static inline T prefactor_A_low_troe(const T& A, size_t molecularity){
      T Anew = A;

      if(molecularity == 2){   //bimolecularity (3rd body is counted)
	Anew = A*1.e-6;
      }

      if(molecularity == 3){   //trimolecularity (3rd body is counted)
	Anew = A*1.e-12;
      }
    
    return Anew;
    }
  };


  template<>
  class ConvertArrheniusDataToConsistentUnits<true,true>{ //!adjust A and Ea
  public:
    
    template<class T>
    static inline T prefactor_A(const T& A, size_t molecularity, bool istroe){
      return ConvertArrheniusDataToConsistentUnits<true,false>::prefactor_A(A,molecularity,istroe);
    }
   
    template<class T>
    static inline T E_activation(const T& Ea){
      return Ea*1.e+3;
    } 

    //!for LOW pressure troe reactions the third body counts!
    template<class T>
    static inline T prefactor_A_low_troe(const T& A, size_t molecularity){
       return ConvertArrheniusDataToConsistentUnits<true,false>::prefactor_A_low_troe(A,molecularity);
     }
  
  };

  template<>
  class ConvertArrheniusDataToConsistentUnits<false,true>{ //!not A, but adjust Ea
  public:
    
    template<class T>
    static inline const T& prefactor_A(const T& A, size_t molecularity, bool istroe){
      return A;
    }

    template<class T>
    static inline T E_activation(const T& Ea){
      return ConvertArrheniusDataToConsistentUnits<true,true>::E_activation(Ea);
    }
  
    template<class T>
    static inline const T& prefactor_A_low_troe(const T& A, size_t molecularity){
      return A;
    }
  };

  
  template<>
  class ConvertArrheniusDataToConsistentUnits<false,false>{ //!not A, not Ea
  public:
    
    template<class T>
    static inline const T& prefactor_A(const T& A, size_t molecularity, bool istroe){
      return A;
    }

    template<class T>
    static inline const T& E_activation(const T& Ea){
      return Ea;
    }
  
    template<class T>
    static inline const T& prefactor_A_low_troe(const T& A, size_t molecularity){
      return A;
    }
    
  };
  

  /**
   * \brief 
   */
  template<char C> class ArrheniusDataIndexing;

  template<>
  class ArrheniusDataIndexing<'a'>{  //!'a'lternating, i.e. \f$...,A_i, \beta_i, E_i, A_{i+1}, \beta_{i+1}, E_{i+1},...\f$   
  public:
    typedef size_t SizeType;
    //! second argument is just dummy here
    static inline const size_t A_offset(size_t index, size_t nfreac){
      return 3*index;
    }

    static inline const size_t beta_offset(size_t index, size_t nfreac){
      return 3*index+1;
    }

    static inline const size_t Ea_offset(size_t index, size_t nfreac){
      return 3*index+2;
    }
  };

 
  template<>
  class ArrheniusDataIndexing<'n'>{  //!'n'on-alternating, i.e. \f$A_1,...,A_n,\beta_1,...,\beta_n, E_1,..., E_n.\f$   
  public:
    typedef size_t SizeType;

    static inline const size_t A_offset(size_t index, size_t nfreac){
      return index;
    }
    
    static inline const size_t beta_offset(size_t index, size_t nfreac){
      return index+nfreac;
    }
    
    static inline const size_t Ea_offset(size_t index, size_t nfreac){
      return index+2*nfreac;
    }
  };


  /**
   * \brief Convert the Arrhenius Data using a singleton (which guarantees that
   * a class has only one instance, thus avoiding multiple invocation)
   *
   * \tparam T precision type of new Arrhenius data
   * \tparam IT1 iterator for old arrhenius data
   * \tparam IT2 iterator for troe and third body information (IMPORTANT: it is 
   *             assumed that first the information if a reaction \f$i\f$ (out 
   *             of the <I>nfreac</I> reactions) is a Troe-reaction (1 =yes, 0
   *             = no) and then the information of the 3rd bodies of size <I>
   *             nreac</I> is given (namely if a reaction contains a 3rd body,
   *             and if so which one)
   * \tparam SM the stoichiometric matrix 
   * \tparam A convert prefactor \f$A\f$ (true = yes, false = no)
   * \tparam E convert activation energy \f$E_a\f$ (true = yes, false = no)
   * \tparam C 'a' or 'n' depending on what storage format for the Arrhenius
   *           data is desired
   *
   * EXAMPLE:
   *\code
     TransformArrheniusDataIntoRightUnits<double,const double*,const double*,StoichType,Aalter,Ealter,storageForm>::instance(nor,AbetaEa,troewtb,SM).convert();
   *\endcode
   *
   * References: 
   *
   *  [1] A. ALEXANDRESCU, "Modern C++ Design", §6, pp. 113, esp. §6.3, §6.4 and §6.5 are worthwhile reading

   */
  template<class T, class IT1, class IT2, class SM, class IT3, bool A, bool E, char C = 'a'>
  class TransformArrheniusDataIntoRightUnits{
  public:
    typedef T value_type;
    typedef std::vector<T> ContainerType;
     
    typedef ArrheniusDataIndexing<C> IndexingType;
    typedef ConvertArrheniusDataToConsistentUnits<A,E> ConvertType;


  private:
    size_t nfreac_;
    const IT1& arrh_;
    const IT2& troewtb_;  
    const SM& stoichmatrix_;
    size_t ntroe_;
    const IT3& troecoeffs_;
    ContainerType newdata_;

    //! make all constructors (and, if existant, all copy 
    //! constructors and copy assignments private)
    TransformArrheniusDataIntoRightUnits(size_t nf = 0, const IT1& arrh = IT1(), const IT2& troewtb = IT2(), const SM& sm = SM(), size_t nt = 0, const IT3& tcoeff = IT3()):nfreac_(nf),arrh_(arrh),troewtb_(troewtb),stoichmatrix_(sm),ntroe_(nt),troecoeffs_(tcoeff),newdata_(3*nf + nt*7){}

    TransformArrheniusDataIntoRightUnits(const TransformArrheniusDataIntoRightUnits&);

    TransformArrheniusDataIntoRightUnits& operator=(const TransformArrheniusDataIntoRightUnits&);

    //! make destructor private
    ~TransformArrheniusDataIntoRightUnits(){}

  public:
   
    //! the Meyer way of creating only one instance -- return local instance
    //! cf. Ref. [1]
    static TransformArrheniusDataIntoRightUnits& instance(size_t nf, const IT1& arrh, const IT2& troewtb, const SM& sm, size_t nt, const IT3& troecoeffs){
      
      static TransformArrheniusDataIntoRightUnits inst(nf,arrh,troewtb,sm,nt,troecoeffs);
      return inst;
    }

    //! transform A, Ea, A_low and Ea_low to SI units 
    ContainerType& convert(const std::string& separator = ",   ",  int prec = 5){
      size_t molec,
	d;

      //if a reaction is troe, store its molecularity(with 3rd body)
      std::vector<size_t> troemolecularity;

      //HIGH PRESSURE coefficients
      for(size_t i = 0; i < nfreac_; ++i){  
	d = 3*i;
	//! the forward reactions in the stoichiometric matrix are the even ones
	molec = stoichmatrix_.molecularity_of_reaction(2*i);

	if(troewtb_[i+nfreac_] != -1){
	  molec += 1;
	  adonis_assert(molec <= 3); 
	}

        newdata_[IndexingType::A_offset(i,nfreac_)] = ConvertType::prefactor_A(arrh_[d],molec,troewtb_[i]);
	
	//remains unchanged
	newdata_[IndexingType::beta_offset(i,nfreac_)] = arrh_[d+1];

	newdata_[IndexingType::Ea_offset(i,nfreac_)] = ConvertType::E_activation(arrh_[d+2]);

	
	if(troewtb_[i]){
	  std::cout << "TROE reaction at index "<< i << " with molecularity "<< molec<< std::endl;
	  troemolecularity.push_back(molec); //store consecutively
	}

	//!test only
	std::cout << std::setprecision(prec) << std::scientific << newdata_[IndexingType::A_offset(i,nfreac_)] << separator << newdata_[IndexingType::beta_offset(i,nfreac_)]<< separator <<newdata_[IndexingType::Ea_offset(i,nfreac_)] << ((i == nfreac_-1) ? " " : separator)<< std::endl;
      }
      
      adonis_assert((int)troemolecularity.size() == (int)(ntroe_));

      //! LOW PRESSURE DATA -- only if ntroe_ > 0
      for(size_t t = 0; t < ntroe_; ++t){
	size_t tix = 7*t, 
	  ix = 3*nfreac_ + tix; //index from which troe coeffs are stored
	                        //store it in any case after the Arrhenius
	                        //coefficients in the order 
	                        //LOW_1|TROE_2|LOW2|TROE2|...

	  newdata_[ix] = ConvertType::prefactor_A_low_troe(troecoeffs_[tix],troemolecularity[t]);
	  newdata_[ix+1] = troecoeffs_[tix+1]; //dimensionless number
	  
	  newdata_[ix+2] = ConvertType::E_activation(troecoeffs_[tix+2]);
      
	  //!test only
	  std::cout << std::endl<< "Alow               betalow           Ealow"<<std::endl <<  std::scientific <<  newdata_[ix] << separator <<newdata_[ix+1] << separator <<  newdata_[ix+2]  <<  ((t == ntroe_-1) ? " " : separator) << std::endl << std::endl;
      }

      return newdata_;
    }
    
  };

  /**
   * \brief convert transport data (given in the Sandia format) into SI units
   *
   * Implemented as singleton. For the question <I>how to use singletons </I>, see aforementioned class.
   */
  template<class T, class ITER, char METRICSYS = 'I'>
  class TransformTransportDataIntoRightUnits{
  public:
    typedef T value_type;
    typedef std::vector<T> ContainerType;

    typedef MetricSystem<T,METRICSYS> UnitsType;
    
  private:
    size_t nspec_;
    const ITER& transportData_;
    ContainerType newdata_;
  
    //! make all constructors (and, if existant, all copy 
    //! constructors and copy assignments private)
    TransformTransportDataIntoRightUnits(size_t n = 0, const ITER& tpd = ITER()):nspec_(n),transportData_(tpd),newdata_(6*n){}

    TransformTransportDataIntoRightUnits(const TransformTransportDataIntoRightUnits&);

    TransformTransportDataIntoRightUnits& operator=(const TransformTransportDataIntoRightUnits&);

    ~TransformTransportDataIntoRightUnits(){}

  public:
    static TransformTransportDataIntoRightUnits& instance(size_t n, const ITER& it){
      static TransformTransportDataIntoRightUnits inst(n,it);

      return inst;
    }
  
    ContainerType& convert(const std::string& separator = ",   ", int prec = 5){
      size_t ix;
      std::string s;
      for(int i = 0; i< prec; ++i)
	s.push_back(' ');  //a whitespace added via push_back(char)
      

      std::cout <<"#Geometry"<< s<< "epsilon/kappa_B"<< s << "sigma"<<s<<"mu"<<s<<"alpha"<< s <<"Z_rot"<< std::endl;
      
      for(size_t t = 0; t < nspec_; ++t){
	ix = 6*t;
	newdata_[ix] = transportData_[ix];
	newdata_[ix+1] = transportData_[ix+1];
	
	//transform the next 3 to SI units
	//newdata_[ix+2] = 1.e-10*transportData_[ix+2]; //Angström to meter
	newdata_[ix+2] = UnitsType::template Angstroem<1>(transportData_[ix+2]);
	//newdata_[ix+3] = 1.e-28*transportData_[ix+3]; //debye in SI
	//newdata_[ix+3] = 1.e-49*transportData_[ix+3];  //debye to kg·m²/s²
	//newdata_[ix+3] = 3.33564e-30*transportData_[ix+3]; //debye to C·m
	newdata_[ix+3] = UnitsType::debye(transportData_[ix+3]);
	//newdata_[ix+4] = 1.e-30*transportData_[ix+4]; //Angström³ to m³ 
	newdata_[ix+4] = UnitsType::template Angstroem<3>(transportData_[ix+4]);

	newdata_[ix+5] = transportData_[ix+5]; 
      
	//! test -- output 
	std::cout<< std::setprecision(prec) << std::scientific <<  newdata_[ix] << separator << newdata_[ix+1] << separator << newdata_[ix+2] << separator << newdata_[ix+3] << separator << newdata_[ix+4] << separator <<  newdata_[ix+5] << ((t == nspec_-1) ? " " : separator) << std::endl;

      }
      return newdata_;
    }
  
  };



  //==========================================================================
  //============ NEEDED FOR CONSTRUCTION OF CHEM. SOURCE TERM ================
  //==========================================================================

  /**
   *\brief Create index binary tree for TROE reactions.
   *
   * Assumption: container troewtb conists of (# of forward reactions)*2
   * entries. The first (# of forward reactions) entries contain 1 or 0, depending on whether reaction \f$i\f$ is a troe reaction or not. The second (# of forward reactions) entries correspond to 3rd body detection (-1 = no 3rd body, 0, 1, 2, ... are M(0), M(1), M(2),...
   * 
   * Here we store, in advance, the reaction where the TROE reaction occurs (as key value), with the troe index (as value). The troe index is used to retrieve the troe coefficients which are stored for M(0), M(1), M(2), ... (in this order)  
   * 
   * NOTE: a call to operator[](size_t) has a complexity of  \f$ \log_2(n)\f$, where \f$n\f$ denotes the input length
   
   */
  class TroeIndex{
  public:
    typedef size_t key_type;
    typedef size_t value_type;
    typedef std::map<key_type,value_type> MapType;

    TroeIndex():numberOfTroeReactions_(0),alreadyCreated_(false){}

    template<class IT>
    MapType& create(size_t nfreac, const IT& troewtb){
      if(!alreadyCreated_){ //only doe this if nothing is created so far
	numberOfTroeReactions_ = 0;
	for(size_t i = 0; i < nfreac; ++i){
	  adonis_assert(troewtb[i] == 0 || troewtb[i] == 1);
	  if(troewtb[i] == 1){
	    troeindextree_[i] = numberOfTroeReactions_++;
	  
	    //test only
	    //std::cout << "map["<<i<<"] = "<< troeindextree_[i] << std::endl;
	  }
	
	}
	alreadyCreated_ = true;
      }
      return troeindextree_;
    }

    size_t number_of_TROE_reactions() const {return numberOfTroeReactions_;}

    const MapType& troe_index() const {return troeindextree_;}

    //! caution: std::map has only a value_type& operator[](i) !!
    value_type& operator[](const key_type& i){
      adonis_assert(alreadyCreated_);
      
      return troeindextree_[i];
      //alternatively, you can commet the above line and uncomment the following
      //return troeindextree_.find(i)-> second;
    }

    //! Kid map's operator[] via map's const_iterator find(const key_type&) const. This works well. Also logarithmic in input size like operator[]
    const value_type& operator[](const key_type& i) const{
       adonis_assert(alreadyCreated_);
      
       return troeindextree_.find(i)->second;
    }

    friend std::ostream& operator<<(std::ostream& os, const TroeIndex& ti){
      os << std::endl << "#Key(*):                 value(**):"<<std::endl;
      for(MapType::const_iterator it = ti.troeindextree_.begin(); it != ti.troeindextree_.end(); ++it)
	os << (*it).first << "                      "<<(*it).second << std::endl;
      os << std::endl << "#nota bene: (*)  refers to index w.r.t. chem. reactions" << std::endl << "#           (**) refers to index w.r.t. # of troe reacs."<<std::endl;
      return os;
    }

  private:
    MapType troeindextree_;
    size_t numberOfTroeReactions_;
    mutable bool alreadyCreated_;
  };

  

  /**
   * \brief the (forward) reaction rate for reaction \f$ i\f$ can be computed 
   * via the three-parameter Arrhenius formula (where the 3 values are 
   * tabulated with the mechanism):
   * \f[ k_{f,i} = AT^{\beta}\exp(-E_a/RT)k\]
   * These values are stored in the following order: prefactor \f$ A\$, exponent \f$ \beta\f$ and activation energy \f$E_a. \f$
   * ATTENTION: Only \f$ \beta \f$ is dimensionless. Usually, it is adviceable
   * to transform units into SI units.
   * Thus for \f$ A_i\f$ we have: bimolecular: \f$ \cdot 1.e-6\f$ and 
   * trimolecular: 1.e-12.
   * Note: for TROE reactions we have: bimolecular: \cdot 1., trimolec.: 1.e-6
   *
   * This can be done via the class <TT>TransformArrheniusDataIntoRightUnits</TT>, which should be invoked once in advance.
   * Caution: sometimes \f$ E_a\f$ is given in kJ/mol (in lieu of J/mol): to
   * get rid of this, the value for \f$ E_a\f$ has to be multiplied by 1e3.
   *
   * \tparam ARRHIT iterator for Arrhenius parameters (for HIGH and LOW pressure
   *    reactions as well as the corresponding troe fit coefficients)
   * \tparam TBDIT iterator for temperature bounds 
   * \tparam TROEIT is a reaction a troe reac, succeeded by 3rd body information
   * \tparam COLLISIONIT iterator for collision efficiencies
   * \tparam C storage style
   */
  template<class ARRHIT, class TBDIT, class TROEIT, class COLLISIONIT, char C = 'a'>
  class ForwardReactionRates{
  public:
    typedef ArrheniusDataIndexing<C> IndexingType;
    typedef typename IndexingType::SizeType SizeType;
    typedef TroeIndex TroeIndexType;
    typedef TroeIndexType* TroeIndexTypePointer;
    typedef TroeIndex::value_type troe_index_value_type;

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType; //! floating point addition type
#endif

    ForwardReactionRates(size_t nspec = 0, size_t nreac = 0,const ARRHIT& a = ARRHIT(), const TBDIT& tbd = TBDIT(), size_t ntroe = 0, const TROEIT& troewtb = TROEIT(), const COLLISIONIT& coll = COLLISIONIT(), TroeIndexTypePointer ti = 0):nspec_(nspec),nfreac_(nreac),arrh_(a),tbd_(tbd),numberOfTroeReacs_(ntroe),troewtb_(troewtb),collision_(coll),troeIndexPtr_(ti){}

    size_t number_of_species() const {return nspec_;}
    size_t number_of_forward_reactions() const {return nfreac_;}
    
    size_t number_of_troe_reactions() const {return numberOfTroeReacs_;}


    SizeType rev_coeffs_offset() const {return 3*nfreac_ + 7*numberOfTroeReacs_;}

    //! forward coefficient \f$ k_{f,i}(temp)\f$ of forward reaction i
    //! this is simply the evaluation of the Arrhenius formula
    template<class TEMP>
    TEMP k_forward(size_t i, const TEMP& temp, bool rev = false) const {
      //! note that we do not check, whether <TT> temp </TT> is suitably chosen
      //! since we \f$k_{f,i}\f$ doesn't depend on species index \f$k\f$
      adonis_assert(i < nfreac_);
      
      //!test me by uncommenting the next line
      //std::cout << std::scientific << "A = "<< arrh_[IndexingType::A_offset(i,nfreac_)]<< "   beta = "<< arrh_[IndexingType::beta_offset(i,nfreac_)] << "   Ea = "<< arrh_[IndexingType::Ea_offset(i,nfreac_)]<<std::endl;
   
      if(rev==false)
	return arrhenius_formula(temp,arrh_[IndexingType::A_offset(i,nfreac_)],arrh_[IndexingType::beta_offset(i,nfreac_)],arrh_[IndexingType::Ea_offset(i,nfreac_)]);
      else{
	SizeType beyondFwdCoeffs = (*this).rev_coeffs_offset();

	return arrhenius_formula(temp,arrh_[beyondFwdCoeffs+IndexingType::A_offset(i,nfreac_)],arrh_[beyondFwdCoeffs+IndexingType::beta_offset(i,nfreac_)],arrh_[beyondFwdCoeffs+IndexingType::Ea_offset(i,nfreac_)]);
      }
	  
    }
    
    
    const ARRHIT& get_thermodata() const {return arrh_;}
    const TROEIT& get_troe_data() const {return troewtb_;}
    const COLLISIONIT& get_collision_efficiencies() const {return collision_;}

    //!the troe coefficients are stored in 7-entries blocks: 0,1,2 correspond 
    //!to A_low,beta_low,Ea_low and 3,4,5,6 to alpha,T***,T* and T**
    //!offset index for troe stuff from newdata_ in previous class 
    //! <TT> t </TT> is the troe-index {0,1,2,..,# of troe reactions} and 
    //! not the index of the troe reaction within the nfreac reactions!
    //! if the second argument is <TT>true</TT> then one computes the reverse
    //! reaction rates for given reverse rate coefficients (in CHEMKIN, the 
    //! the keyword '\rev' is used. 
    size_t troe_index(size_t t, bool reverse = false) const {
      if(reverse==false)
	return 3*nfreac_ + 7*t;  //!they are stored <I> after</I>  Arrh. coeffs. 
      else
	return (*this).rev_coeffs_offset() + 3*nfreac_ + 7*t;
    }

    size_t index_A_low(size_t t, bool reverse = false) const {
	return troe_index(t,reverse);
    }
    
    size_t index_beta_low(size_t t, bool reverse = false) const {
	return troe_index(t,reverse) +1;
    }

    size_t index_Ea_low(size_t t, bool reverse = false) const {
      return troe_index(t,reverse) +2;
    }

    //!these are the fitting parameters for eq. (9.117), p. 390: \$\alpha\$
    size_t index_troe_alpha(size_t t, bool reverse = false) const {
      return troe_index(t,reverse) + 3;
    }

    //! \f$T^{***}\f$
    size_t index_troe_t3star(size_t t, bool reverse = false) const {
      return troe_index(t,reverse) + 4;
    }
    
    //! \f$T^{*}\f$
    size_t index_troe_t1star(size_t t, bool reverse = false) const {
      return troe_index(t,reverse) + 5;
    }

    //! \f$T^{**}\f$
    size_t index_troe_t2star(size_t t, bool reverse = false) const {
      return troe_index(t,reverse) + 6;
    }

    
     /** 
     *  \brief This is needed for the calculation of \f$k_f\f$ for TROE 
     *  reactions which use \f$ [M] := \sum_{k=1}^K [X_k]\f$ as well.
     * See e.g. 
     *         [CHEMKIN-PRO, Theory Manual, p. 41]
     *
     *  NOTE: depending on entities to be used one has to tranform them into
     *  the right quantities first. For instance, U can be a container/
     *  expression whose atomic entries are <I> specific moles </I> 
     *  \f$ \zeta_k\f$. The following formula should be kept in mind: 
     *  \f[ [X_k] = \frac{n_k}{V} = \frac{\rho w_k}{M_k} = \rho \zeta_k.\f]
     * 
     * In any case we need the species' concentration, \f$ [[X_k]]_{k=1}^K\f$ 
     * here for u. Hence one has to determine the quantities to calculate with, 
     * e.g. \f$ \rho w_k\f$, \f$ [X_k]\f$ or 
     */
    template<class TEMP, class U> //! U is considered the concentration here
    TEMP k_forward(size_t i, const TEMP& temp, const U& u, bool rev = false) const {
     
     
      if(troewtb_[i] == 1){  //!reaction is a troe reaction which has a forward
                             //!reaction rate that must be computed differently
	                     //!Precisely, species concentrations are involved
	
	//std::cout<<"compute k_f for TROE reaction (i = "<<i<<")."<< std::endl;
	
	troe_index_value_type tix = troeIndexPtr_->operator[](i);

	//! previously with <TT>counttroe_</TT>  
	TEMP k0 = arrhenius_formula(temp,arrh_[(*this).index_A_low(tix,rev)], arrh_[(*this).index_beta_low(tix,rev)], arrh_[(*this).index_Ea_low(tix,rev)]),
	  
	  kinfty = (*this).k_forward(i,temp);

	//std::cout << std::scientific << "k0 = "<< k0 << "   kinfty = "<< kinfty << std::endl;

	TEMP totalConcentrationOfMixture = TEMP();

      //!third body detected at reaction \f$i\f$
      //!troewtb_[i+nfreac_] = 0,1,... specifies the third body M(0),M(1),...
	if(troewtb_[i+nfreac_] != -1){  //this is always(?) the case
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	  TEMP s = TEMP(), 
	    c = TEMP();
#endif

	  int tindex = static_cast<int>(troewtb_[i+nfreac_]);
	
	  //std::cout << "3rd body M("<<tindex << ")" << std::endl;

	  //! note that for third bodies we have \f$ [M] = \sum_{k=1}^K\alpha_{i,k}[X_k]\f$
	  for(size_t k = 0; k < nspec_; ++k){
	    adonis_assert((tbd_[3*k] <= temp) && (temp <= tbd_[3*k+1]));
	  
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	    totalConcentrationOfMixture = AdditionType::add(collision_[tindex*nspec_ + k]*u[k],s,c);
#else
	    totalConcentrationOfMixture += collision_[tindex*nspec_ + k]*u[k];  
#endif	    
	  }
      
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	  CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s,c);
	  if(c != TEMP()){  //! if c == 0 then do not assign s since 
	    //! it might be 0 because no addition was involved!
	    CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(totalConcentrationOfMixture,s);
	  }
#endif 
	}
      	
	//std::cout << "[M] = "<< totalConcentrationOfMixture << std::endl;
     
	TEMP pr = k0/kinfty*totalConcentrationOfMixture,  // reduced pressure
	  Fcent = (1-arrh_[(*this).index_troe_alpha(tix,rev)])*exp(-temp/(arrh_[(*this).index_troe_t3star(tix,rev)])) + arrh_[(*this).index_troe_alpha(tix,rev)]*exp(-temp/arrh_[(*this).index_troe_t1star(tix,rev)]) + exp(-arrh_[(*this).index_troe_t2star(tix,rev)]/temp);   //formula (9.117)
	
	//std::cout << "pr = "<< pr << "   Fcent = "<< Fcent << std::endl;

	//! for real numbers, the logarithm isn't defined for arguments < 0
#ifndef NDEBUG
	//adonis_assert((Fcent > TEMP()) && (pr > TEMP()));
	if(Fcent <= TEMP()){
	  ADONIS_ERROR(ValueError,"Fcent = "<< Fcent << "\n   (must be > 0 since it acts as argument to log).");
	}
	if(pr <= TEMP()){
	  //prior to some error, give some information
	  ADONIS_ERROR(ValueError,"pr := k_0·[M]/k_infty = "<< pr << "\n   (must be > 0 since it acts as argument to log10)\n   [M] = "<< totalConcentrationOfMixture << "\n  [X] = "<< u << "\n   k_0 = "<< k0 << "\n   k_infty = "<< kinfty << ".");
	}
#endif						     
	TEMP l10Fc = log10(Fcent), //formulas (9.114) -- (9.116)
	
	  cp = -0.4 - 0.67*l10Fc,
	  np = 0.75 - 1.27*l10Fc,
	  dp = 0.14,
	
	  term = log10(pr) + cp,
	  den = np - dp*term;
      
	//std::cout << "c = "<<cp << "   n = "<< np << "   d = "<< dp << "   den = "<<den <<std::endl; 

	adonis_assert(den != TEMP());
	
	TEMP F = l10Fc/(1.0+(ntimes<2>(term/den)));
	
	//std::cout << "F = "<<F << " (before exp)"<<std::endl;
 
	F = exp(F*log(10.));
	
	//std::cout << "F = "<< F << std::endl;
#ifndef NDEBUG
	adonis_assert((is_well_defined_value(F)) && (is_well_defined_value(pr)) && (is_well_defined_value(kinfty)));
#endif	

	return kinfty*pr/(1. + pr)*F;

      }
      else  //! ok no Troe -- use standard computation of \f$ k_{r,i}\f$
	return (*this).k_forward(i,temp,rev);
    }




  private:
    size_t nspec_;
    size_t nfreac_; //!number of forward direction reactions
    ARRHIT arrh_;
    TBDIT tbd_;
    size_t numberOfTroeReacs_;
    TROEIT troewtb_;
    COLLISIONIT collision_;
    TroeIndexTypePointer troeIndexPtr_;
  };
  



  
  
  /**
   * \brief Compute the reverse rate constants \f$ k_{r,i}\f$ of reaction \f$i\f$
   * Remark: everything is constructed from iterators and/or pointers. 
   *
   * Reference:
   *   [1] KEE, COLTRIN and Glarborg, <I> Chemically reacting flow</I>, Wiley-Interscience, 2003, p. 386
   */
  template<class FORWARD, class SM, class NASAPOLY>
  class ReverseReactionRates{
  public:
    //o.k. -- we store pointers that are initialized to 0 (NULL)
    typedef FORWARD* FwdPointerType;
    typedef SM* StoichPointerType;
    typedef NASAPOLY* NasaPointerType;
    
    typedef typename NASAPOLY::TemperatureBoundsIterType TemperatureBoundsIterType;
    
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType; //! floating point addition type
#endif

    ReverseReactionRates(FwdPointerType f = 0,  //default: point to no object 
			 StoichPointerType s = 0,  
			 NasaPointerType nasa = 0):fwdRatesPointer_(f),stoichMatrixPointer_(s),thermoNasaPointer_(nasa){}


   
    //!this isn't working if you are working with (const) references instead
    //! of raw pointers
    ReverseReactionRates& operator=(const ReverseReactionRates& rrr){
      if( this!= &rrr){
	fwdRatesPointer_ = rrr.fwdRatesPointer_;
	stoichMatrixPointer_ = rrr.stoichMatrixPointer_;
	thermoNasaPointer_ = rrr.thermoNasaPointer_;
      }
      return *this;
    }

   
    size_t number_of_species() const {return stoichMatrixPointer_->number_of_species();}  

    size_t number_of_forward_reactions() const {return fwdRatesPointer_->number_of_forward_reactions();}

    TemperatureBoundsIterType bounds_on_T() const {return thermoNasaPointer_->temperature_bounds();}

    template<class T> 
    T C_p(size_t k, const T& temp) const{
#ifndef NDEBUG   
      if((temp < thermoNasaPointer_->temperature_bounds()[3*k]) || (temp > thermoNasaPointer_->temperature_bounds()[3*k+1]) || (!is_well_defined_value(temp)))
	ADONIS_ERROR(BoundsError, "T = "<<temp << " left its bounds.");
#endif
      return thermoNasaPointer_->heat_capacity(k,temp);
    }
   
    template<class T> 
    T H_T(size_t k, const T& temp) const{
      return thermoNasaPointer_->enthalpy(k,temp);
    }

    template<class T> 
    T S_T(size_t k, const T& temp) const{
      return thermoNasaPointer_->entropy(k,temp);
    }


    size_t nu_forward(size_t k, size_t i) const{
      return stoichMatrixPointer_->nu_forward(k,i);
    }

    size_t nu_reverse(size_t k, size_t i) const{
      return stoichMatrixPointer_->nu_reverse(k,i);
    }

    int nu_net(size_t k, size_t i) const{
      return stoichMatrixPointer_->nu_net(k,i);
    }

    //! calculate the reverse coefficient, compare [1, eq. (9.95), p.386], of
    //!  reaction \f$i\f$
    template<class TEMP>
    TEMP k_reverse(size_t i, const TEMP& temp, bool rev = false) const {
      TEMP DSri = TEMP(), //the change in entropy and enthalpy -- cf. formulas
	DHri = TEMP();    // (9.94), p. 386

      int nu_ki,
	nu_i = 0;

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      TEMP s1 = TEMP(), s2 = TEMP(),
	c1 = TEMP(), c2 = TEMP();
#endif

      for(size_t k = 0; k < stoichMatrixPointer_->number_of_species(); ++k){
	adonis_assert(thermoNasaPointer_->check_temperature_bounds(k,temp,UniversalConstants<TEMP>::enlargeRange));

	//! Recall that the stoichiometric matrix stores \nu from both -->
	//! and <-- reactions. Here, only forward reactions will be considered. 
	//! This means that \f$ i \f$ is even (--> reactions).
	nu_ki = stoichMatrixPointer_->nu_net(k,2*i);

	if(nu_ki != 0){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	  DSri = AdditionType::add(nu_ki*thermoNasaPointer_->entropy(k,temp),s1,c1);
	  DHri = AdditionType::add(nu_ki*thermoNasaPointer_->enthalpy(k,temp),s2,c2);
#else
	  DSri += nu_ki*thermoNasaPointer_->entropy(k,temp);
	  DHri += nu_ki*thermoNasaPointer_->enthalpy(k,temp);
#endif  
	//eq. (9.92), p. 386 -- only integers (no floating point addition)
	  nu_i += nu_ki;    
	}
      }
      
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s1,c1);
      if(c1 != TEMP()){  //! if c == 0 then do not assign s since 
	                 //! it might be 0 because no addition was involved!
	//! Anyway, here only addition is applied
	CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(DSri,s1);
      }
      
      CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s2,c2);
      if(c2 != TEMP()){
	CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(DHri,s2);
      }
#endif      

      //not that even when temp = 0 hence the exponent is -inf, exp returns
      // the right result, namely zero :)
      TEMP K_pi = exp(DSri/PhysicalConstants<TEMP>::IdealGasConstant - DHri/(PhysicalConstants<TEMP>::IdealGasConstant*temp)),
	
	K_ci = K_pi*natural_power(PhysicalConstants<TEMP>::p0/(PhysicalConstants<TEMP>::IdealGasConstant*temp),nu_i);   //eq. (9.91), p. 386
      
      adonis_assert(K_ci != TEMP()); //prohibit division by zero
      
      return fwdRatesPointer_->k_forward(i,temp,rev)/K_ci;  //eq. (9.95), p. 386
    }
    

    //!should be prefered 
    template<class TEMP,class U>
    TEMP k_reverse(size_t i, const TEMP& temp, const U& u, bool rev = false) const {
      TEMP DSri = TEMP(), //the change in entropy and enthalpy -- cf. formulas
	DHri = TEMP();    // (9.94), p. 386

      int nu_ki,
	nu_i = 0;

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      TEMP s1 = TEMP(), s2 = TEMP(),
	c1 = TEMP(), c2 = TEMP();
#endif

      for(size_t k = 0; k < stoichMatrixPointer_->number_of_species(); ++k){
#ifndef NDEBUG
	if(!thermoNasaPointer_->check_temperature_bounds(k,temp,UniversalConstants<TEMP>::enlargeRange))
	  ADONIS_ERROR(RangeError,"Temperature out of bounds: T = "<< temp << ".");
#endif	

	//adonis_assert(thermoNasaPointer_->check_temperature_bounds(k,temp,UniversalConstants<TEMP>::enlargeRange));

	//! Recall that the stoichiometric matrix stores \nu from both -->
	//! and <-- reactions. Here, only forward reactions will be considered. 
	//! This means that \f$ i \f$ is even (--> reactions).
	nu_ki = stoichMatrixPointer_->nu_net(k,2*i);

	if(nu_ki != 0){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
	  DSri = AdditionType::add(nu_ki*thermoNasaPointer_->entropy(k,temp),s1,c1);
	  DHri = AdditionType::add(nu_ki*thermoNasaPointer_->enthalpy(k,temp),s2,c2);
#else
	  DSri += nu_ki*thermoNasaPointer_->entropy(k,temp);
	  DHri += nu_ki*thermoNasaPointer_->enthalpy(k,temp);
#endif  
	//eq. (9.92), p. 386 -- only integers (no floating point addition)
	  nu_i += nu_ki;    
	}
      } 
      
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s1,c1);
      if(c1 != TEMP()){  //! if c == 0 then do not assign s since 
	                 //! it might be 0 because no addition was involved!
	//! Anyway, here only addition is applied
	CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(DSri,s1);
      }
      
      CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s2,c2);
      if(c2 != TEMP()){
	CorrectRoundingErrorAfterwards<TEMP,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(DHri,s2);
      }
#endif      

      //not that even when temp = 0 hence the exponent is -inf, exp returns
      // the right result, namely zero :)
      TEMP K_pi = exp(DSri/PhysicalConstants<TEMP>::IdealGasConstant - DHri/(PhysicalConstants<TEMP>::IdealGasConstant*temp)),
	
	K_ci = K_pi*natural_power(PhysicalConstants<TEMP>::p0/(PhysicalConstants<TEMP>::IdealGasConstant*temp),nu_i);   //eq. (9.91), p. 386
      
      adonis_assert(K_ci != TEMP()); //prohibit division by zero
      
      if(fwdRatesPointer_->get_troe_data()[i] == 1)
	return fwdRatesPointer_->k_forward(i,temp,u,rev)/K_ci;  //eq. (9.95), p. 386
      else{  //no Troe reaction, hence no special forward to compute
	if(rev == false) //that's the usual way of computation
	  return fwdRatesPointer_->k_forward(i,temp,rev)/K_ci;
	else
	  return fwdRatesPointer_->k_forward(i,temp,rev); //no equilib. constant
      }
    }
    
    
    StoichPointerType get_stoichiometic_matrix() const {return &stoichMatrixPointer_;}

  private:
    FwdPointerType fwdRatesPointer_;          
    StoichPointerType stoichMatrixPointer_;     
    NasaPointerType thermoNasaPointer_;  
  };
  

  
  /**
   * \brief Isothermal ReactionRates -- here we don't need thermochemical data
   *
   * For isothermal mechanisms, the stuff presented in Kee et. al., pp. 386, concerning computing reverse reaction rates, isn't needed any more...
   */
  template<class V, class STOICH>
  class IsoThermalReverseReactionRates{
  public:
    typedef STOICH* StoichPointerType;
    typedef typename ValueTypeSelector<V,IsADuneFieldContainer<V>::Value >::value_type value_type;
    
    typedef const value_type* TemperatureBoundsIterType; //pro forma only coz
    // reaction rate has it too

    //! this is not needed here but included for reasons of compliance 
    #if USE_TUNED_FLOATING_POINT_ARITHMETIC
    typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> AdditionType; //! floating point addition type
#endif

    IsoThermalReverseReactionRates(const V& kf = V(), StoichPointerType sm = 0):kf_(kf),sm_(sm){}

    template<class InputIterator>
    IsoThermalReverseReactionRates(InputIterator s, InputIterator e, StoichPointerType sm = 0):kf_(static_cast<size_t>(std::distance(s,e))),sm_(sm){
     
      size_t dim = static_cast<size_t>(std::distance(s,e));
      InputIterator it = s;
      for(size_t i = 0; i < dim; ++i, ++it){
	kf_[i] = *it;
      }
    }
    
    size_t number_of_species() const {return sm_->number_of_species();}  

    size_t number_of_forward_reactions() const {return sm_->number_of_forward_reactions();}

    size_t nu_forward(size_t k, size_t i) const{
      return sm_->nu_forward(k,i);
    }

    size_t nu_reverse(size_t k, size_t i) const{
      return sm_->nu_reverse(k,i);
    }

    int nu_net(size_t k, size_t i) const{
      return sm_->nu_net(k,i);
    }


    //! 4th argument is just a dummy
    value_type k_reverse(size_t i, const value_type& temp, bool rev = false) const {
      return kf_[i];
    }
    
    //! 4th argument is just dummy here
    template<class U>
    value_type k_reverse(size_t i, const value_type& temp, const U& u, bool rev  = false) const {
      return kf_[i];
    }
    
  private:
    V kf_;
    StoichPointerType sm_;
  };


}//end namespace

#endif
