#ifndef TEMPERATURE_SOURCE_TERM_HH
#define TEMPERATURE_SOURCE_TERM_HH

#include "thermochemistry.hh"
#include "../expressiontemplates/xsettings.hh"
#include "../accuracy/floatingpointarithmetic.hh"

namespace Adonis{

  /**
   * \brief The following class snippets are meant to compute the 
   * temperature source term for isochoric (I = 1) and isobaric (I = 2) 
   * realistic thermochemical scenarios
   */
  template<int I, class T> class ComputeTemperatureTerms;

  template<class T> 
  class ComputeTemperatureTerms<1,T>{
  public:
    template<class BCSPTR>
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
     static inline void calculate_sum1(T& sm1, size_t k, BCSPTR bcsPtr, const T& temp, const T& rhsk, T& s, T& c){
      sm1 = IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::add((bcsPtr->H_T(k,temp) - PhysicalConstants<T>::IdealGasConstant*temp)*rhsk,s,c); 
    }
#else
    static inline void calculate_sum1(T& sm1, size_t k, BCSPTR bcsPtr, const T& temp, const T& rhsk){
       sm1 += (bcsPtr->H_T(k,temp) - PhysicalConstants<T>::IdealGasConstant*temp)*rhsk;
    }
#endif
    

    /*  static inline void calculate_sum1(T& sm1, size_t k, BCSPTR bcsPtr, const T& temp, const T& rhsk, T& s, T& c){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      sm1 = IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::add((bcsPtr->H_T(k,temp) - PhysicalConstants<T>::IdealGasConstant*temp)*rhsk,s,c); 
#else
      sm1 += (bcsPtr->H_T(k,temp) - PhysicalConstants<T>::IdealGasConstant*temp)*rhsk;
#endif
    }
    */  

    template<class T2, class T3> 
    static inline void T_t(T& rhsk, const T2& rho, const T3& p, const T& sum1, const T& sum2, const T& phi_sum, const T& temp){
      rhsk = -1./rho*(sum1/(sum2 - PhysicalConstants<T>::IdealGasConstant*phi_sum));
    }

  };

  template<class T> 
  class ComputeTemperatureTerms<2,T>{
  public:
    template<class BCSPTR>
    static inline void calculate_sum1(T& sm1, size_t k, BCSPTR bcsPtr, const T& temp, const T& rhsk, T& s, T& c){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      sm1 = IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::add(bcsPtr->H_T(k,temp)*rhsk,s,c); 
#else
      sm1 += (bcsPtr->H_T(k,temp)*rhsk);
#endif
    }
  
    template<class T2, class T3> 
    static inline void T_t(T& rhsk, const T2& rho, const T3& p, const T& sum1, const T& sum2, const T& phi_sum, const T& temp){
      rhsk =  -PhysicalConstants<T>::IdealGasConstant*temp/PhysicalConstants<T>::p0*((sum1*phi_sum)/sum2);
      
    }
  };



  template<int I, class T, class BUILDCHEM>
  class TemperatureSourceTerm{
  public:
    typedef BUILDCHEM* BcsPointerType;

  private:
    T phi_sum_,
      sum1_,
      sum2_,
      temp_;
    
    BcsPointerType bcsPtr_;

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    T s_,
      c_,
      s1_,
      c1_,
      s2_,
      c2_;
#endif
    
  public:

#if USE_TUNED_FLOATING_POINT_ARITHMETIC
    typedef IntelligentArithmeticOperation<ExprTmpl::XPCTSettings::floatingPointAdditionMethod> OperationType;
#endif
    

    TemperatureSourceTerm(const T& t = T(), BcsPointerType bcsptr = 0): phi_sum_(T()),sum1_(T()),sum2_(T()),temp_(t),bcsPtr_(bcsptr)
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
								      ,s_(T()),c_(T()),s1_(T()),c1_(T()),s2_(T()),c2_(T())									
#endif
{} 

    
    void initialize(const T& t, BcsPointerType bcsptr){
      temp_ = t;
      bcsPtr_ =  bcsptr;
    }

    
    void set_temperature(const T& t){
      temp_ = t;
    }


    void reset(){
      phi_sum_ = sum1_ = sum2_ = temp_ = T();
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      s_ = c_ = s1_ = c1_ =  s2_ = c2_ = T();
#endif
      
    }



    void compute_phi_sum(const T& uk){
#if  USE_TUNED_FLOATING_POINT_ARITHMETIC
      phi_sum_ = OperationType::add(uk,s_,c_);
#else
      phi_sum_ += uk;
#endif
    }
 

    void compute_sum2(size_t k, const T& uk){
#if  USE_TUNED_FLOATING_POINT_ARITHMETIC
      sum2_ = OperationType::add(uk*bcsPtr_->C_p(k,temp_),s1_,c1_);
#else 
      sum2_ += uk*bcsPtr_->C_p(k,temp_);
#endif
    }

    
    void compute_sum1(size_t k, const T& rhsk){
#if USE_TUNED_FLOATING_POINT_ARITHMETIC
      ComputeTemperatureTerms<I,T>::calculate_sum1(sum1_,k,bcsPtr_,temp_,rhsk,s2_,c2_);
#else
      ComputeTemperatureTerms<I,T>::calculate_sum1(sum1_,k,bcsPtr_,temp_,rhsk);
#endif
    }

   

    template<class T2, class T3> 
    void temperature(T& rhsk, const T2& rho, const T3& p){
      #if  USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<T,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s_,c_);
      if(c_ != T())
	CorrectRoundingErrorAfterwards<T,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(phi_sum_,s_);
#endif

      #if  USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<T,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s1_,c1_);
      if(c1_ != T())
	CorrectRoundingErrorAfterwards<T,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(sum2_,s1_);
#endif

      #if  USE_TUNED_FLOATING_POINT_ARITHMETIC
      CorrectRoundingErrorAfterwards<T,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::update_s(s2_,c2_);
      if(c2_ != T())
	CorrectRoundingErrorAfterwards<T,ExprTmpl::XPCTSettings::floatingPointAdditionMethod>::assign(sum1_,s2_);
#endif
      

      ComputeTemperatureTerms<I,T>::T_t(rhsk,rho,p,sum1_,sum2_,phi_sum_,temp_);
    }
    

    const T& get_temperature() const {return temp_;}
    T& get_temperature() {return temp_;}
  };
}

#endif
