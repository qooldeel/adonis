#ifndef CREATE_INJECTIVE_MAPPING_FOR_REDUCED_SPECIES_AND_REAC_INDEX_HH
#define CREATE_INJECTIVE_MAPPING_FOR_REDUCED_SPECIES_AND_REAC_INDEX_HH

#include "../common/adonisassert.hh"

namespace Adonis{

  /**
   * \brief The Base class 
   *
   *  Access the right index. Using <TT> CommonIndexInjector</TT> is equivalent 
   *  to using \f$ i \in I:=\{ 0,1,\ldots,n\}.\f$. When using ReducedIndexInjector an 
   *  injective mapping \f$ \iota: I_1 \subset I \Rightarrow I \$ from a
   *  reduced index \f$I_1\f$ to the normal index \f$I\f$ is created. 
   */
  template<class K>
  class IndexInjector{
  private:
    //reference is established to derived class via <I> Barton-Nackman </I>
    K& ref2derived(){return static_cast<K&>(*this);}

    //this handles implementations in derived classes, i.e.
    // <TT> spec_index, reac_index</TT> must defined in each derived class
    // independently
    size_t species(size_t k){ return ref2derived().spec_index(k);}
    
    size_t reaction(size_t i){return ref2derived().reac_index(i);}

  protected:
    size_t nspec_,
      nreac_;

    //define common functionality which uses individual implementations of
    //derived classes
  public:
    size_t number_of_species() {return nspec_;}
    size_t number_of_forward_reactions() {return nreac_;}

    size_t s_ix(size_t k) {return species(k);} //delegate stuff to derived class
    size_t r_ix(size_t i) {return reaction(i);} //via ref2derived().

  };


  /**
   * \brief Index mapper for usual indexing. See description above.
   */
  class CommonIndexInjector: public IndexInjector<CommonIndexInjector>{
  public:
    typedef IndexInjector<CommonIndexInjector> BaseClassType;

    CommonIndexInjector(size_t nspec = 0, size_t nreac = 0){
      //! note that you cannot write these in the initializer list ;)
      BaseClassType::nspec_ = nspec;
      BaseClassType::nreac_ = nreac;
    }

    CommonIndexInjector(const CommonIndexInjector& cii){
       BaseClassType::nspec_ = cii.nspec_;
       BaseClassType::nreac_ = cii.nreac_;
    }

    CommonIndexInjector& operator=(const CommonIndexInjector& cii){
      if(this != &cii){
	BaseClassType::nspec_ = cii.nspec_;
	BaseClassType::nreac_ = cii.nreac_;
      }
      return *this;
    }

    void initialize(size_t nspec, size_t nreac){
      BaseClassType::nspec_ = nspec;
      BaseClassType::nreac_ = nreac;
    }

    size_t spec_index(size_t k) {return k;} //just give back k
    size_t reac_index(size_t i) {return i;} //just give back i
    
  };


  /**
   * \brief Create reduced index mapping. See above. 
   *
   * \tparam IT1 stores index of reduced species
   * \tparam IT2 stores index of reduced reactions
   *
   * To be more precise, <TT> s_ </TT> stores the domain of the mapping
   * \f$ S: s_ \Rightarrow K \f$ and <TT> r_</TT> the domain of the mapping
   * \f$ R: r_ \Rightarrow I \f$, where \f$K\f$ and \f$I\f$ are the index sets
   * of variables (e.g. chemical species) and relations on them (e.g. reactions)
   * respectively.
   */
  template<class IT1, class IT2>
  class ReducedIndexInjector: public IndexInjector<ReducedIndexInjector<IT1,IT2> >{
  public:
    typedef IndexInjector<ReducedIndexInjector<IT1,IT2> > BaseClassType;

    ReducedIndexInjector(size_t nspec = 0, size_t nreac = 0, const IT1& i1 = IT1(), const IT2& i2 = IT2()):s_(i1),r_(i2){
      BaseClassType::nspec_ = nspec;
      BaseClassType::nreac_ = nreac;
    }

    ReducedIndexInjector(const ReducedIndexInjector& rii):s_(rii.s_),r_(rii.r_){
      BaseClassType::nspec_ = rii.nspec_;
      BaseClassType::nreac_ = rii.nreac_;
    }
    

    ReducedIndexInjector& operator=(const ReducedIndexInjector& rii){
      if(this != &rii){
	s_ = rii.s_;
	r_ = rii.r_;
	BaseClassType::nspec_ = rii.nspec_;
	BaseClassType::nreac_ = rii.nreac_;
      }
      return *this;
    }

    void initialize(size_t nspec, size_t nreac, const IT1& i1, const IT2& i2){
       BaseClassType::nspec_ = nspec;
       BaseClassType::nreac_ = nreac;
       s_ = i1;
       r_ = i2;
    }


    //!give back appropriate reduced indices (subset of indices)
    size_t spec_index(size_t k) {
      adonis_assert(k < BaseClassType::nspec_);
      return s_[k];
    } 
    size_t reac_index(size_t i) {
      adonis_assert(i < BaseClassType::nreac_);
      return r_[i];
    } 

  private:   //store iterators only 
    IT1 s_;  //subset index for species (e.g. only rpv)
    IT2 r_;  //subset index for reactions (e.g. only reactions with rpv)
 
  };

} //end namespace 

#endif
