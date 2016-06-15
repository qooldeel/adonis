#ifndef MECHANISM_TYPE_NOTIFIER_FOR_AUTOMATICALLY_GENERATED_CHEMISTRY_HH
#define MECHANISM_TYPE_NOTIFIER_FOR_AUTOMATICALLY_GENERATED_CHEMISTRY_HH

#include <string>

namespace Adonis{


  /**
   * \brief Type for identifying a automatically generated mechanism
   */
  class ChemMechTypePresent{
  public:
    typedef std::string StringType;

    ChemMechTypePresent(const StringType& name = "No name"):nameOfMechanism_(name){}
  private:
    StringType nameOfMechanism_;
  };

} // end namespace

#endif
