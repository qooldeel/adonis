#ifndef MOL_IDENTIFIER_TO_CHECK_WHETHER_FUNCTOR_DESCRIBES_MOL_HH
#define  MOL_IDENTIFIER_TO_CHECK_WHETHER_FUNCTOR_DESCRIBES_MOL_HH

namespace Adonis{


  
  /**
   * \brief Type to be used to mark a functor as one describing a method of line
   *        related discretization
   */
  template<int DIM, int ENCOD>
  class MOLFunctorIdentifier{
  public:
    enum{
      Dimension=DIM,
      Id=ENCOD
    };
  };



  /**
   * \brief The actual MOL functor encodings
   */
  class FunctorID{
  public:
    enum{
      burgerSchiesser = 0,
      burgerMarc,
      burger,
      cylindricalCDRSchiesser,
      cylindricalCDRMarc,
      H2conaireFull,
      
      USED_IDs //! must be the last entry
    };
  };

} //end namespace


#endif

