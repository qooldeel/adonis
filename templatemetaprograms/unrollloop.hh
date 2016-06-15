#ifndef LOOP_TO_UNROLL_HH
#define LOOP_TO_UNROLL_HH

#include<iostream>

#include "compiletimearray.hh"
#include "commonfunctions.hh"
#include "../expressiontemplates/exprvec.hh"
#include "../accuracy/floatingpointarithmetic.hh"
#include "../common/globalfunctions.hh"
#include "../common/smartassign.hh"
#include "../moleculartransport/transportsettings.hh"
#include "../misc/misctmps.hh"

namespace Adonis{
  /**
    DESCRIPTION: A template to unroll a for-loop, e.g.

       for(int i = 0; i < 10; ++i){
       
        //statement stand here                        (1)
       }
          
       
    ADVANTAGE: Increases performance 
    
    DISADVANTAGE: Complicated code
    
    USAGE:  for-loop in (1) is unwrapped 

         UnrollLoop<0,10>::loop();

     

    REF.: [LIPPMAN, C++ Report, § "Using C++ Template Metaprograms",p.468/69 ]


    * here, only "standard" loops are implemented, i.e. copying a random access
    * container to another, etc. 
   */
  template<int I, int N>
  class UnrollLoop{
  private:
    enum{unwind_ = (I+1) != N};

  public:

    //!simplify notation ;)
    typedef UnrollLoop<unwind_ ? (I+1) : N, N> NextIterType; 



    template<class IT>
    static inline void print(IT it){
      std::cout << *it << " ";
      ++it;
      NextIterType::print(it);
    }
   
    template<class C>
    static inline void print_container(const C& c){
      std::cout << c[I] << " ";
      NextIterType::print_container(c);
    }

    //! for iterators
    template<class IT1,class IT2> 
    static inline void assign(IT1 it1, IT2 it2)  {
      //++it1; ++it2;  //increase iterators
      smart_assign(*it1,*it2);
      ++it1; ++it2;
      UnrollLoop<unwind_ ? (I+1) : N, N>::assign(it1,it2);
    }


    template<class IT, class T>
    static inline void assign_value(IT it, const T& val){
      smart_assign(*it,val);
      ++it;
      NextIterType::assign_value(it,val);
    }


    //! assign random access container x2 2 random access container x1
    template<class X1, class X2>
    static inline void assign_rac2rac(X1& x1, const X2& x2){
      smart_assign(x1[I],x2[I]);   //x1[I] = x2[I]
      UnrollLoop<unwind_ ? (I+1) : N, N>::assign_rac2rac(x1,x2);
    }
    
    //!assign rac 2 rac via, but only indices given in an index vector
    //! this is only for equally sized containers !
    template<class IX, class X1, class X2>
    static inline void assign_rac2rac(X1& x1, const X2& x2, const IX& index){
      smart_assign(x1[index[I]],x2[index[I]]);
      UnrollLoop<unwind_ ? (I+1) : N, N>::assign_rac2rac(x1,x2,index);
    }

    template<class IX, class X1, class X2>
    static inline void assign_greater_2_smaller_rac(X1& x1, const X2& x2, const IX& index){
      smart_assign(x1[I],x2[index[I]]);
      UnrollLoop<unwind_ ? (I+1) : N, N>::assign_greater_2_smaller_rac(x1,x2,index);
    }

    template<class IX, class X1, class X2>
    static inline void assign_smaller_2_greater_rac(X1& x1, const X2& x2, const IX& index){
      smart_assign(x1[index[I]],x2[I]);
      UnrollLoop<unwind_ ? (I+1) : N, N>::assign_smaller_2_greater_rac(x1,x2,index);
    }

    static inline void print_cool_stuff(){
      std::cout << I << " ";
      UnrollLoop<unwind_ ? (I+1) : N, N>::print_cool_stuff();
    }

    /* //THIS IS NOT WORKING; j CANNOT APPEAR IN A CONSTANT-EXPRESSION!!!
      static inline void attemp(const int j){
       UnrollLoop<0,j>::print_cool_stuff();
       NextIterType::attemp(j);
       }*/


    template<int M>
    static inline void matrix_style(){
      UnrollLoop<0,M>::print_cool_stuff(); //loop in loop
      std::cout << std::endl;           //!don't forget the template here
      UnrollLoop<unwind_ ? (I+1) : N, N>::template matrix_style<M>();
    }

    //!conventional addition of linear object of size N having operator-[]
    template<class A>
    static inline typename A::value_type conventional_addition(const A& a){
      return a[I] + NextIterType::conventional_addition(a);
    }

    template<int P, class A>
    static inline typename A::value_type sum(const A& a){
      return ntimes<P>(a[I]) + NextIterType::template sum<P>(a);
    }
    
    template<int P, class A, class D>
    static inline typename A::value_type sum(const D& diag, const A& a){
      return diag[I]*ntimes<P>(a[I]) + NextIterType::template sum<P>(diag,a);
    }
    
    template<class X>
    static inline typename X::value_type product(const X& x){
      return x[I]*NextIterType::product(x);
    }

    template<class IT1, class IT2>
    static inline typename IT1::value_type sum_product(IT1 x, IT2 w){
      return x[I]*w[I] + NextIterType::sum_product(x,w);
    }

    template<class V, class IT>
    static inline typename V::value_type sum_prod(const V& x, IT w){
      return x[I]*w[I] + NextIterType::sum_prod(x,w);
    }

    //! possibiliy of raise entries to a given integer power
    template<int M1, int M2, class IT1, class IT2>
    static inline typename IT1::value_type sum_product_xtended(IT1 x, IT2 w){
      return ntimes<M1>(x[I])*ntimes<M2>(w[I]) + NextIterType::template sum_product_xtended<M1,M2>(x,w);
    }

    template<int M1, int M2, class V, class IT>
    static inline typename V::value_type sum_product_x(const V& x, IT w){
      return ntimes<M1>(x[I])*ntimes<M2>(w[I]) + NextIterType::template sum_product_xtended<M1,M2>(x,w);
    }
    

    //!tuned floating point addition of a linear object of size N 
    //! having operator-[]
    template<char C, class A>
    static inline typename A::value_type& tuned_addition(const A& a, typename A::value_type& s, typename A::value_type& c){
      TunedFloatingPointAddition<C,A>::addition(I,a,s,c);
      return UnrollLoop<unwind_ ? (I+1) : N, N>::template tuned_addition<C>(a,s,c);
    }

    //Hessian stuff
    template<class V, class X, class L, class SEQ>
    static inline V& hessian_terms_added(V& hess, V& hessI, const X& x, const L& lambda, SEQ& adseq){
      hessI = adseq.Hessian(x,I);
      hess += lambda[I]*hessI;
      return NextIterType::hessian_terms_added(hess,hessI,x,lambda,adseq);
    }

    template<class V, class X, class L, class SEQ, class S>
    static inline V& hessian_terms_added(V& hess, V& hessI, const X& x, const L& lambda, SEQ& adseq, const S& Sdiag){
      hessI = adseq.Hessian(x,I);
      hess += lambda[I]*Sdiag[I]*hessI;
      return NextIterType::hessian_terms_added(hess,hessI,x,lambda,adseq,Sdiag);
    }

    //! K-th row, N serves as number of cols
    //! assume A is stored row-wise (C-style)
    template<int K, class V, class A>  
    static inline V& get_row(V& v, const A& a){
      v[I] = a[K*N+I];
      return NextIterType::template get_row<K>(v,a);
    }
    

    template<int K, class TPTDATATYPE, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline T calculate_Theta_reduced(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac){
      return !Kronecker<K,I>::delta()*tdr.theta(K,I,//AccessElement<RpvIndexType,K>::result,AccessElement<RpvIndexType,I>::result, // // would yield error
temp,Xfrac) + NextIterType::template calculate_Theta_reduced<K,TPTDATATYPE>(tdr,temp,Xfrac);
    }

    template<int K, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline T calculate_Theta(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac){
      return !Kronecker<K,I>::delta()*tdr.theta(K,I,temp,Xfrac) + NextIterType::template calculate_Theta<K>(tdr,temp,Xfrac);
    }

    
    template<bool DISTINCTMOLMASS, class TPTDATATYPE, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline void compute_thermal_diffusion_ratios_reduced(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac, const typename THERMALDIFFUSIONRATIO::value_type& molecmass){
      typedef typename THERMALDIFFUSIONRATIO::DataType DataType;
      (DISTINCTMOLMASS==false) ? smart_assign(tdr[I],UnrollLoop<0,DataType::rednspec>::template calculate_Theta_reduced<I,TPTDATATYPE>(tdr,temp,Xfrac)) :
	smart_assign(tdr[I], ((DataType::molar_masses()[I] <= molecmass) ? UnrollLoop<0,DataType::rednspec>::template calculate_Theta_reduced<I,TPTDATATYPE>(tdr,temp,Xfrac) : 0));
      NextIterType::template compute_thermal_diffusion_ratios_reduced<DISTINCTMOLMASS,TPTDATATYPE>(tdr,temp,Xfrac,molecmass);
    }


    template<bool DISTINCTMOLMASS, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline void compute_thermal_diffusion_ratios(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac, const typename THERMALDIFFUSIONRATIO::value_type& molecmass){
      typedef typename THERMALDIFFUSIONRATIO::DataType DataType;
      (DISTINCTMOLMASS==false) ? smart_assign(tdr[I],UnrollLoop<0,DataType::nspec>::template calculate_Theta<I>(tdr,temp,Xfrac)) : 
	smart_assign(tdr[I], ((DataType::molar_masses()[I] <= molecmass) ? UnrollLoop<0,DataType::nspec>::template calculate_Theta<I>(tdr,temp,Xfrac) : 0));
      NextIterType::template compute_thermal_diffusion_ratios<DISTINCTMOLMASS>(tdr,temp,Xfrac,molecmass);
    }

    
    //! mean molecular weight for <B> mass fracion </B> 
    template<class TPTDATATYPE, class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_Y_reduced(T& mmw, const IT1& y, const IT2& w){
      typedef typename TPTDATATYPE::RpvIndexType RpvIndexType;
      mmw += convert_number(y[I]/w[AccessElement<RpvIndexType,I>::result]); 
      return NextIterType::template mean_molecular_weight_Y_reduced<TPTDATATYPE>(mmw,y,w);
    }

    template<class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_Y(T& mmw, const IT1& y, const IT2& w){
      mmw += convert_number(y[I]/w[I]);  
      return NextIterType::mean_molecular_weight_Y(mmw,y,w);
    }

     //! mean molecular weight for <B> mole fracion </B> 
    template<class TPTDATATYPE, class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_X_reduced(T& mmw, const IT1& x, const IT2& w){
      typedef typename TPTDATATYPE::RpvIndexType RpvIndexType;
      mmw += convert_number(x[I]*w[AccessElement<RpvIndexType,I>::result]); 
      return NextIterType::template mean_molecular_weight_X_reduced<TPTDATATYPE>(mmw,x,w);
    }

    template<class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_X(T& mmw, const IT1& x, const IT2& w){
      mmw += convert_number(x[I]*w[I]);
      return NextIterType::mean_molecular_weight_X(mmw,x,w);
    }
      

    template<class ITER, class T>
    static inline void calibrate_against_singularity(ITER it, const T& delta){
      it[I] += delta;
      NextIterType::calibrate_against_singularity(it,delta);
    }

    template<class ITER, class JITER>
    static inline void calibrate_against_singularity_multiple_perturbations(ITER it, const JITER& delta){
      it[I] += convert_number(delta[I]);
      NextIterType::calibrate_against_singularity_multiple_perturbations(it,delta);
    }

   
    template<class VIS, class T>
    static inline void compute_viscosities_reduced(VIS& vis, const T& temp){
      typedef typename VIS::DataType DataType;
      typedef typename DataType::RpvIndexType RpvIndexType;

      smart_assign(vis.get_viscosity(I),vis.template viscosity_coefficient<TransportSettings::molarMassInGrams>(AccessElement<RpvIndexType,I>::result,temp));
      NextIterType::compute_viscosities_reduced(vis,temp);
    }

    template<class VIS, class T>
    static inline void compute_viscosities(VIS& vis, const T& temp){
      smart_assign(vis.get_viscosity(I),vis.template viscosity_coefficient<TransportSettings::molarMassInGrams>(I,temp));
      NextIterType::compute_viscosities(vis,temp);
    }

    //!only select indices from the rpv index when computing this quantity
    template<class COND, class T, class RHO, class P>
    static inline void compute_conductivities_reduced(COND& cond, const RHO& rho, const P& p, const T& temp){
      typedef typename COND::DataType DataType;
      typedef typename DataType::RpvIndexType RpvIndexType;
      
      smart_assign(cond.get_conductivity(I),cond.template conductivity_coefficient<TransportSettings::molarMassInGrams>(AccessElement<RpvIndexType,I>::result,rho,p,temp));
      NextIterType::compute_conductivities_reduced(cond,rho,p,temp);
    }

    template<class COND, class T, class RHO, class P>
    static inline void compute_conductivities(COND& cond, const RHO& rho, const P& p, const T& temp){
      smart_assign(cond.get_conductivity(I),cond.template conductivity_coefficient<TransportSettings::molarMassInGrams>(I,rho,p,temp));
      NextIterType::compute_conductivities(cond,rho,p,temp);
    }

    template<class TPTDATATYPE,int K, class MW, class X>
    static inline typename X::value_type mix_avg_diff_nominator_reduced(const MW& mw, const X& Xfrac){
      typedef typename TPTDATATYPE::RpvIndexType RpvIndexType;
      //std::cout << " W_"<< AccessElement<RpvIndexType,I>::result << " = "<<mw[AccessElement<RpvIndexType,I>::result] <<std::endl;
      return (!Kronecker<K,I>::delta())*Xfrac[I]*mw[AccessElement<RpvIndexType,I>::result] + NextIterType::template mix_avg_diff_nominator_reduced<TPTDATATYPE,K>(mw,Xfrac);
    }

    template<int K, class MW, class X>
    static inline typename X::value_type mix_avg_diff_nominator(const MW& mw, const X& Xfrac){
      return (!Kronecker<K,I>::delta())*Xfrac[I]*mw[I] + NextIterType::template mix_avg_diff_nominator<K>(mw,Xfrac);
    }

    template<class TPTDATATYPE,int K, class X, class PITER>
    static inline typename X::value_type mix_avg_diff_denominator_reduced(const X& Xfrac, const PITER& bindiff){
      //typedef typename TPTDATATYPE::RpvIndexType RpvIndexType;
      return (!Kronecker<K,I>::delta())*Xfrac[I]/(bindiff[SymmetricDenseMatrixAccess<I,K,N>::offset()]) + NextIterType::template mix_avg_diff_denominator_reduced<TPTDATATYPE,K>(Xfrac,bindiff);
    }

    template<int K, class X, class PITER>
    static inline typename X::value_type mix_avg_diff_denominator(const X& Xfrac, const PITER& bindiff){
      return (!Kronecker<K,I>::delta())*Xfrac[I]/(bindiff[SymmetricDenseMatrixAccess<I,K,N>::offset()]) + NextIterType::template mix_avg_diff_denominator<K>(Xfrac,bindiff);
    }

    //! implements formula 5-45, p. 90 in 
    //! [CHEMKIN-PRO, Reaction Design: San Diego, 2008]
    //! NOTE: This is only well defined when the mixture isn't exaclty 1 species
    template<class TPTDATATYPE, class DMIXIT, class MW, class T, class X, class PITER>
    static inline void compute_mixture_averaged_diffusion_coefficients_reduced(DMIXIT dmixit, const MW& mw, const T& meanmw, const X& Xfrac, const PITER& bindiff){
      smart_assign(dmixit[I],UnrollLoop<0,N>::template mix_avg_diff_nominator_reduced<TPTDATATYPE,I>(mw,Xfrac)/(meanmw*UnrollLoop<0,N>::template mix_avg_diff_denominator_reduced<TPTDATATYPE,I>(Xfrac,bindiff)));
      
      NextIterType::template compute_mixture_averaged_diffusion_coefficients_reduced<TPTDATATYPE>(dmixit,mw,meanmw,Xfrac,bindiff);
    }
    
    
    template<class DMIXIT, class MW, class T, class X, class PITER>
    static inline void compute_mixture_averaged_diffusion_coefficients(DMIXIT dmixit, const MW& mw, const T& meanmw, const X& Xfrac, const PITER& bindiff){
      smart_assign(dmixit[I],UnrollLoop<0,N>::template mix_avg_diff_nominator<I>(mw,Xfrac)/(meanmw*UnrollLoop<0,N>::template mix_avg_diff_denominator<I>(Xfrac,bindiff)));
      
      NextIterType::compute_mixture_averaged_diffusion_coefficients(dmixit,mw,meanmw,Xfrac,bindiff);
    }


   

    template<class ITER, class T>
    static inline void random(ITER it, const T& eps){
      it[I] =  eps*static_cast<T>(rand())/RAND_MAX;
      NextIterType::random(it,eps);
    }

    //just test ;)
    static inline void loop(){
      
      //======== concrete loop statement stands here ================
     
      std::cout << I<<".) two times I = "<<2.*I<<std::endl; 
    
      //=============================================================

      UnrollLoop<unwind_ ? (I+1) : N, N>::loop();
    }
  
  };

  
  //Specialisation for base case of recursion (END RECURSION)
  template<int N >
  class UnrollLoop<N, N>{
  public:

    template<class IT>
    static inline void print(IT it){
      std::cout << std::endl; //do nothing any more
    }

    template<class C>
    static inline void print_container(const C& c){
      std::cout << std::endl; //do nothing any more
    }

    template<class IT1, class IT2> 
    static inline void assign(IT1 it1, IT2 it2){} //do nothing
    
    template<class IT, class T>   //do nothing
    static inline void assign_value(IT it, const T& val){}

    template<class X1, class X2>  //do nothing
    static inline void assign_rac2rac(X1& x1, const X2& x2){}
    
    template<class IX, class X1, class X2>
    static inline void assign_greater_2_smaller_rac(X1& x1, const X2& x2, const IX& index){} //do nothing

    template<class IX, class X1, class X2>
    static inline void assign_smaller_2_greater_rac(X1& x1, const X2& x2, const IX& index){}

    template<class IX, class X1, class X2> //do nothing
    static inline void assign_rac2rac(X1& x1, const X2& x2, const IX& index){}

    static inline void print_cool_stuff(){std::cout << std::endl;}

    template<int M>
    static inline void matrix_style(){}


    template<class A>
    static inline typename A::value_type conventional_addition(const A& a){
       return typename A::value_type(0);
     }

    template<int P, class A>
    static inline typename A::value_type sum(const A& a){
      return typename A::value_type(0);
    }

    template<int P, class A, class D>
    static inline typename A::value_type sum(const D& diag, const A& a){
      return typename A::value_type(0);
    }
    
    template<class X>
    static inline typename X::value_type product(const X& x){
      return typename X::value_type(1);
    }

    template<class IT1, class IT2>
    static inline typename IT1::value_type sum_product(IT1 x, IT2 w){
      return typename IT1::value_type(0);
    }
    
    template<class V, class IT>
    static inline typename V::value_type sum_prod(const V& x, IT w){
      return typename V::value_type(0);
    }


    template<int M1, int M2, class IT1, class IT2>
    static inline typename IT1::value_type sum_product_xtended(IT1 x, IT2 w){
      return typename IT1::value_type(0);
    }

    template<int M1, int M2, class V, class IT>
    static inline typename V::value_type sum_product_x(const V& x, IT w){
      return V::value_type(0);
    }


    template<char C, class A>
    static inline typename A::value_type& tuned_addition(const A& a, typename A::value_type& s, typename A::value_type& c){
      CorrectRoundingErrorAfterwards<typename A::value_type,C>::update_s(s,c);
      //!in any case s contains the sum
      return s;
    }
    
    template<class V, class X, class L, class SEQ>
    static inline V& hessian_terms_added(V& hess, V& hessI, const X& x, const L& lambda, SEQ& adseq){ return hess;}

    template<class V, class X, class L, class SEQ, class S>
    static inline V& hessian_terms_added(V& hess, V& hessI, const X& x, const L& lambda, SEQ& adseq, const S& Sdiag){ return hess;}


    template<int K, class V, class A>  
    static inline V& get_row(V& v, const A& a){return v;}


    template<int K, class TPTDATATYPE, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline T calculate_Theta_reduced(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac){
      return T();
    }
      
    
    template<int K, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline T calculate_Theta(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac){
      return T();
    }

    template<bool DISTINCTMOLMASS, class TPTDATATYPE, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline void compute_thermal_diffusion_ratios_reduced(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac, const typename THERMALDIFFUSIONRATIO::value_type& molecmass){}


    template<bool DISTINCTMOLMASS, class THERMALDIFFUSIONRATIO, class T, class X>
    static inline void compute_thermal_diffusion_ratios(THERMALDIFFUSIONRATIO& tdr, const T& temp, const X& Xfrac, const typename THERMALDIFFUSIONRATIO::value_type& molecmass){}


    template<class TPTDATATYPE, class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_Y_reduced(T& mmw, const IT1& y, const IT2& w){
      return (mmw = 1./mmw);
    }
    
    
    template<class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_Y(T& mmw, const IT1& y, const IT2& w){
      return (mmw = 1./mmw);
    }


    template<class TPTDATATYPE, class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_X_reduced(T& mmw, const IT1& x, const IT2& w){ 
      return mmw;
    }

     template<class T, class IT1, class IT2>
    static inline T& mean_molecular_weight_X(T& mmw, const IT1& x, const IT2& w){ 
      return mmw;
    }


    template<class ITER, class T> //do nothing
    static inline void calibrate_against_singularity(ITER it, const T& delta){}

    
    template<class ITER, class JITER>
    static inline void calibrate_against_singularity_multiple_perturbations(ITER it, const JITER& delta){} //do nothing


    template<class VIS, class T>  //do nothing
    static inline void compute_viscosities_reduced(VIS& vis, const T& temp){}
    
  
    template<class VIS, class T>  //do nothing
    static inline void compute_viscosities(VIS& vis, const T& temp){}

    template<class COND, class T, class RHO, class P>
    static inline void compute_conductivities_reduced(COND& cond, const RHO& rho, const P& p, const T& temp){} //do nothing
    
    template<class COND, class T, class RHO, class P>
    static inline void compute_conductivities(COND& cond, const RHO& rho, const P& p, const T& temp){} //do nothing


    template<class TPTDATATYPE,int K, class MW, class X>
    static inline typename X::value_type mix_avg_diff_nominator_reduced(const MW& mw, const X& Xfrac){
      typedef typename X::value_type value_type;
      return value_type();
    }
    
    template<int K, class MW, class X>
    static inline typename X::value_type mix_avg_diff_nominator(const MW& mw, const X& Xfrac){
      typedef typename X::value_type value_type;
      return value_type();
    }
      
    template<class TPTDATATYPE,int K, class X, class PITER>
    static inline typename X::value_type mix_avg_diff_denominator_reduced(const X& Xfrac, const PITER& bindiff){
       typedef typename X::value_type value_type;
       return value_type();
    }

    template<int K, class X, class PITER>
    static inline typename X::value_type mix_avg_diff_denominator(const X& Xfrac, const PITER& bindiff){
       typedef typename X::value_type value_type;
       return value_type();
    }


    template<class TPTDATATYPE, class DMIXIT, class MW, class T, class X, class PITER> //!do nothing
    static inline void compute_mixture_averaged_diffusion_coefficients_reduced(DMIXIT dmixit, const MW& mw, const T& meanmw, const X& Xfrac, const PITER& bindiff){}

    template<class DMIXIT, class MW, class T, class X, class PITER> //!do noth.
    static inline void compute_mixture_averaged_diffusion_coefficients(DMIXIT dmixit, const MW& mw, const T& meanmw, const X& Xfrac, const PITER& bindiff){
    }



    template<class V, class T>
    static inline void random(V& v, const T& eps){} //do nothing

    //!just test ;)
    static inline void loop(){}   //do nothing
  };



}//end namespace 


#endif
