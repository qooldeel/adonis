#ifndef CONC_FULL_DOT_OMEGA_WITH_TEMPERATURE_ONLY_NEEDED_FOR_REACTIVE_FLOWS_HH
#define CONC_FULL_DOT_OMEGA_WITH_TEMPERATURE_ONLY_NEEDED_FOR_REACTIVE_FLOWS_HH

#include "../massactionkinetics/physicalconstants.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../massactionkinetics/data/thermochemicaldata.hh"
#include "../templatemetaprograms/unrollloop.hh"

namespace Adonis{

  /**
   * \brief Compute \f$ \dot{\omega} \equiv \dot{\omega}(T,[X])\f$, i.e. 
   * the chemical net production rate which is needed for the reaction part
   * of reactive flows.
   * This class implements: \f[ [\dot{X}_k] = \dot{\omega}_k(T,[X]), \quad k = 1, \ldots, n \\ \dot{T} = \sum_{k=1}^n H_k^0(T)\dot{\omega}_k(T,[X]).\f]
   * Note that the <BB>last</BB> entry always corresponds to temperature
   */
  template<class T, int ENCOD> 
  class DotOmegaWithTemperatureFull{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type DType;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef ExprTmpl::MyVec<T> VType;
    typedef std::size_t SizeType;

    typedef ThermoData4Mechanism<DType,ENCOD> DataType;

    enum{TMP=DataType::nspec,  //index where temperature is to be found
	 DIM=DataType::nspec+1
    };

    typedef StaticArray<T,DataType::nspec> ArrayType; 
    typedef StaticArray<T,DIM> DomainType;

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

    //!ASSUMPTION: first nspec entries of rhs_ contain concentrations,
    //!            last one temperature
    DotOmegaWithTemperatureFull():heat_(0.),rho_(0.),cp_(0),p_(101325.),wbar_(0.),
				  rhs_(DataType::nspec+1){ //concentrations + temperature
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

    const SizeType size() const {return rhs_.size();}
    const SizeType dim() const {return rhs_.size();}
    const SizeType domain_dim() const {return rhs_.size();}

    const ArrayType& dot_omega() const {return dotomega_;}

    const ArrayType& concentrations() const {return concentrations_;}

    const T& wbar() const {return wbar_;}
    const T& rho() const {return rho_;}
    const T& heat() const {return heat_;}
    const T& cp() const {return cp_;}


    //use this before calling operator() to cope with variable pressure
    template<class T1>  
    void set_pressure(const T1& press){
      smart_assign(p_,press);
    }

    //! rhs_ contains concentrations, where the last entry is temperature
    //! MASS fraction version, i.e. eval contains mass fractions and temperature
    template<class X>
    VType& operator()(const X& eval){
      heat_ = rho_ = cp_ = T(); //reset
      UnrollLoop<0,DIM>::assign(Y_.begin(),eval.begin()); //copy
     
      //grant mass conservation
      Y_[DIM-2] = 1. - UnrollLoop<0,DIM-2>::template sum<1>(Y_);
   
      //! calculate concentrations
      wbar_ = T(); //reset
      UnrollLoop<0,DataType::nspec>::mean_molecular_weight_Y(wbar_,Y_.begin(),DataType::molar_masses());
      
      rho_ = (p_*wbar_)/(PhysicalConstants<DType>::Rgas*Y_[TMP]);
      
      T invrho = 1./rho_;
      //concentrations and calculate mixture cp, rho is needed here
      for(int k = 0; k < DataType::nspec; ++k){
	cp_ += (BCS_.C_p(k, Y_[TMP])/DataType::molar_masses()[k]*Y_[k]);
	concentrations_[k] = rho_*Y_[k]/DataType::molar_masses()[k];
      }

      for(int k = 0; k < DataType::nspec; ++k){
	smart_assign(dotomega_[k],BCS_.net_formation_rate(k,Y_[TMP],concentrations_)); //compute dot omega
	//TODO: outside this loop?
	rhs_[k] = invrho*(dotomega_[k]*DataType::molar_masses()[k]);
	heat_ += BCS_.H_T(k,Y_[TMP])*dotomega_[k];
      }
      
      rhs_[DataType::nspec] = -invrho*1./cp_*heat_; //mind the minus here!

      return rhs_;
    }


    // //! xt contains concentrations, where the last entry is temperature
    // //! CONCENTRATION version, i.e. eval containes concentrations and temp.
    // template<class X>
    // VType& operator()(const X& eval){
    //   heat_ = rho_ = cp_ = T(); //reset
    //   UnrollLoop<0,DataType::nspec>::assign(concentrations_.begin(),&eval[0]);
      
    //   for(int k = 0; k < DataType::nspec; ++k){
    // 	smart_assign(dotomega_[k],BCS_.net_formation_rate(k,eval[TMP],concentrations_)); //compute dot omega
    // 	xt_[k] = dotomega_[k];
    // 	heat_ += BCS_.H_T(k,eval[TMP])*dotomega_[k];
	
    // 	rho_ += eval[k]*DataType::molar_masses()[k];  
    //   }

    //   //calculate mixture cp, rho is needed here
    //   for(int k = 0; k < DataType::nspec; ++k){
    // 	cp_ += (BCS_.C_p(k, eval[TMP])*eval[k])/rho_;
    //   }

    //   xt_[DataType::nspec] = -1./(rho_*cp_)*heat_; //mind the minus here!

    //   return xt_;
    // }



  private:
    T heat_, rho_, cp_, p_, wbar_;
    VType rhs_;
    ArrayType dotomega_, concentrations_;
    DomainType Y_;
    BuilderType BCS_;
  };





  //! DERIVATIVE 
  #if USE_CPPAD
  template<class T, int ENCOD, template<typename D2, class A = std::allocator<D2> >  class VEC = ExprTmpl::MyVec>
  class JacobianOfDotOmegaWithTemperatureFull{
  public:
    typedef T value_type;
    typedef typename TypeAdapter<T>::Type DType;
    typedef CppAD::ADFun<T> ADFunType;
    typedef VEC<T> VecType;
    typedef CppAD::AD<T> ADType;
    typedef VEC<ADType> VecADType;
    //!in the case of CppAD this will be the Base of CppAD::AD<T>, to wit T. 
    typedef typename TypeAdapter<T>::Type BaseType; 
    typedef VEC<BaseType> BaseVecType;  

    typedef typename DotOmegaWithTemperatureFull<DType,ENCOD>::DataType DataType;


  private:
    CppAD::ADFun<T> adseq_;        //AD sequence 
    VecType jac_;   

  public:
    JacobianOfDotOmegaWithTemperatureFull(){}

    template<class W>
    void initialize(const W& init){
      int dim = DataType::nspec+1;  //concentrations + temperature
      DotOmegaWithTemperatureFull<ADType,ENCOD> chemistry;
      VecADType X(dim),Y(dim);
      for(int i = 0; i < dim; ++i)
	X[i] = init[i];
      CppAD::Independent(X);
      Y = chemistry(X);
      adseq_.Dependent(X,Y);
      adseq_.optimize();
    }
  
     inline VecType& jacobian(const VecType& v){
       
#ifndef NDEBUG
       if(adseq_.Domain() == 0)
	 ADONIS_ERROR(DerivedError, "CppAD::ADFun object is still uninitialised (# of independent variables in operation sequence is ZERO and does not match size of function argument vector).\n Hence I cannot create Jacobian, master Marc.\n You should invoke member 'NablaS::set(...)' to revise that ;-)");

     
       if(adseq_.Domain() != v.size())
	 ADONIS_ERROR(DerivedError,"Size of evaluation point and those of AD sequence do NOT match");
#endif  

       jac_ = adseq_.Jacobian(v);
       return jac_;
     }

    
    inline VecType& sparse_jacobian(const VecType& v){

#ifndef NDEBUG
      if(adseq_.Domain() == 0)
	ADONIS_ERROR(DerivedError, "CppAD::ADFun object is still uninitialised (# of independent variables in operation sequence is ZERO and does not match size of function argument vector).\n Hence I cannot create sparse Jacobian, master Marc.\n You should invoke member 'NablaS::set(...)' to revise that ;-)");


      if(adseq_.Domain() != v.size())
	  ADONIS_ERROR(DerivedError,"Size of evaluation point and those of AD sequence do NOT match");
#endif
	 
      jac_ = adseq_.SparseJacobian(v); 
      return jac_;
      }


  };
  
#endif






} //end namespace

#endif
