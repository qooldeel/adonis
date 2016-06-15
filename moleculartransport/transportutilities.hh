#ifndef STUFF_NEEDED_IN_THE_COURSE_OF_CALCULATING_TRANSPORT_PROPERTIES_HH
#define STUFF_NEEDED_IN_THE_COURSE_OF_CALCULATING_TRANSPORT_PROPERTIES_HH

namespace Adonis{

  //! compute \f$ \phi_{jk}\f$ needed for the calulation of the mixture-averaged
  //! viscosity \$f \mu. \f$
  template<int K, int J, class TPTDATATYPE, bool B> class Phi4MixVis;

  
  //!use reduction (true)
  template<int K, int J, class TPTDATATYPE>
  class Phi4MixVis<K,J,TPTDATATYPE,true>{  
  public:
    typedef TPTDATATYPE DataType;
    typedef typename DataType::value_type value_type;
    
    typedef typename DataType::RpvIndexType RpvIndexType;


    template<class IT>
    static inline value_type calculate(const IT& vis){
      return (1./(sqrt(8.*(1 + DataType::molar_masses()[AccessElement<RpvIndexType,K>::result]/(DataType::molar_masses()[AccessElement<RpvIndexType,J>::result]))))*ntimes<2>((1 + sqrt(vis[K]/vis[J])*sqrt(sqrt(DataType::molar_masses()[AccessElement<RpvIndexType,J>::result]/(DataType::molar_masses()[AccessElement<RpvIndexType,K>::result]))))));
    }
  };

  
  //!use full representation (false)
  template<int K, int J, class TPTDATATYPE>
  class Phi4MixVis<K,J,TPTDATATYPE,false>{  
  public:
    typedef TPTDATATYPE DataType;
    
    typedef typename DataType::value_type value_type;
   
    template<class IT>
    static inline value_type calculate(const IT& vis){
      return (1./(sqrt(8.*(1 + DataType::molar_masses()[K]/(DataType::molar_masses()[J]))))*ntimes<2>((1 + sqrt(vis[K]/vis[J])*sqrt(sqrt(DataType::molar_masses()[J]/(DataType::molar_masses()[K]))))));
    }
  };

} //end namespace

#endif
