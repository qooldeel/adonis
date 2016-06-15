#ifndef 1D_MOL_WRAPPER_HH
#define 1D_MOL_WRAPPER_HH

#include "../common/globalfunctions.hh"
#include "../common/fancymessages.hh"
#include "../common/tempusfugit.hh"
#include "../io/readinparameters.hh"
#include "choose.hh"
#include "../ode/ode.hh"

namespace Adonis{

  /**
   * \brief Just wrap implicit ODE solver with space-discretized function
   */
  template<class T, template<class D> class DISCRFUN, bool SPARSE, bool SCALE>
  class MOL1D{
  public:

    typedef typename TypeAdapter<T>::BaseType BaseType;

    MOL1D(){}

    template<class V>
    void solve(const V& u0){
      std::cout << " SCALING?: ";
      (SCALE) ? FancyMessages().nice_output("yes",32) : FancyMessages().display_in_bold("no");
      
       ParameterData PD;
       PD.read_from_file(ChooseFDSetting::FILE2READIN);
       
       unsigned maxsteps = PD.get_datum<unsigned>("maxsteps"),
	 maxNewt = PD.get_datum<unsigned>("maxNewt");

       BaseType t0 =  PD.get_datum<BaseType>("t0"),
	 tend =  PD.get_datum<BaseType>("tend"),
	 
	 ntol = PD.get_datum<BaseType>("ntol"),
	 atol = PD.get_datum<BaseType>("atol"),
	 hEstim = PD.get_datum<BaseType>("hEstim"),
	 Cscal = PD.get_datum<BaseType>("Cscal"),
	 hmin = PD.get_datum<BaseType>("hmin"),
	 hmax = PD.get_datum<BaseType>("hmax");
   
       int printprec = PD.get_datum<int>("printprec");
       std::string filename = PD.get_datum<std::string>("filename");

       typedef ODE<T,DISCRFUN,1,SPARSE,SCALE> ModelType; 
       typedef typename ModelType::IndexVecType IndexVecType;
       IndexVecType selectVariables; //use when only some variables are to be
       // written to file for printing 

       ModelType ivp(t0,tend,u0);
       CPUTime<BaseType> cpu;
       cpu.start();
       ivp.implicit_euler(maxsteps,maxNewt,ntol,atol,filename,hEstim,Cscal,hmin,hmax,selectVariables,ChooseFDSetting::FILE2READIN,printprec);
       BaseType elapsedTime = cpu.stop();
       time_in_human_readable_form(elapsedTime);
    }

  };

}


#endif
