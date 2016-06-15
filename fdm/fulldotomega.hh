#ifndef FULL_DOT_OMEGA__ONLY_NEEDED_FOR_REACTIVE_FLOWS_HH
#define FULL_DOT_OMEGA__ONLY_NEEDED_FOR_REACTIVE_FLOWS_HH

#include "../expressiontemplates/exprvec.hh"
#include "../massactionkinetics/data/thermochemicaldata.hh"

namespace Adonis{

  /**
   * \brief Compute \f$ \dot{\omega} \equiv \dot{\omega}(T,[X])\f$, i.e. 
   * the chemical net production rate which is needed for the reaction part
   * of reactive flows.
   */
  template<class T, int ENCOD> 
  class FullDotOmega{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type DType;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef ExprTmpl::MyVec<T> VType;

    typedef ThermoData4Mechanism<DType,ENCOD> DataType;

    typedef StaticArray<T,DataType::nspec> ArrayType; 

    //To build thermo-chemistry
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;
    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;

    //Other Index
    typedef IndexHandler<DataType,false> AccessType;  //always full here
    typedef typename AccessType::IndexerType IndexerType;

    typedef ReverseReactionRates<FwdType,StoichType,ThermoType> RevType;
    typedef BuildChemicalSourceTerm<FwdType,RevType,IndexerType> BuilderType;
    typedef typename BuilderType::RevPointerType RevPointerType;

    typedef EquationOfState<'i'> EosType;  //ideal gas law


    FullDotOmega(){
       //!construct objects needed for source term
      //! CAUTION: use <TT>static</TT> to maintain each objects' location 
      //!   beyond the call of the constuctor (i.e. beyond the localness of 
      //!   the {...}-block
      //! NOTE: the data are assumed to be stored for the FULL model
      static ThermoType nasa(DataType::nspec,DataType::thermo(),DataType::temperature_bounds());
       
      static StoichType stoichmatrix(DataType::nspec,DataType::nreac,DataType::stoichiometric_matrix());

      static TroeIndex TIndex;
      TIndex.create(DataType::nreac,DataType::troewtb());
	
      //! everything stated as  pointers	
      static FwdType forwardRates(DataType::nspec,DataType::nreac,DataType::arrhenius_values(),DataType::temperature_bounds(), DataType::ntroereac, DataType::troewtb(),DataType::collision_efficiencies(),&TIndex);

      //! everything stated as  pointers	
      static RevType reverseRates(&forwardRates,&stoichmatrix,&nasa);

      //!now either a reduced or the usual index is used for constructing
      //! \f$ \dt{\omega}\f$
      //static CommonIndexInjector Ixer(DataType::nspec,DataType::nreac);
      static IndexerType Ixer;
      AccessType::initialize(Ixer);

      BCS_.initialize(&forwardRates,&reverseRates,&Ixer);
    }

    const int size() const {return DataType::nspec;}

    template<class X>
    ArrayType& operator()(const typename X::value_type& temperature, const X& concentrations){
      for(int k = 0; k < DataType::nspec; ++k){
	//dotomega_[k] = BCS_.net_formation_rate(k,temperature,concentrations);
	smart_assign(dotomega_[k],BCS_.net_formation_rate(k,temperature,concentrations));
      }
      return dotomega_;
    }


  private:
    ArrayType dotomega_;
    BuilderType BCS_;
  };

} //end namespace

#endif
