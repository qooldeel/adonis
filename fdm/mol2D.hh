#ifndef METHOD_OF_LINES_2D_WRAPPER_HH
#define METHOD_OF_LINES_2D_WRAPPER_HH

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
  template<int DIM,class T, template<class D> class DISCRFUN, bool SPARSE, bool SCALE, bool PATTERNCREATIONVIAAD=false,bool REPAIR=false,char REPAIRTYPE='n',bool PRINTPRESSURE=false>
  class MOL{
  public:

    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef ExprTmpl::MyVec<T> VType;

    MOL():elapsedTime_(BaseType()){}

    //! \param u0 initial condition for mol
    //! \param fname specification file for implicit euler method and spatial
    //! settings, etc.
    //! \param outname if this string is empty, take <TT>filename</TT> from read-in file ChooseFDSetting::FILE2READIN
    //! \param timefile Write elapsed time, date, etc. into a file after MOL
    //!        computations have been performed (if string is not empty)
    //! \param pMolarMasses These may be used when the pressure is also to be
    //!        recorded alongside the other fluid variables.
    //! \param fromLastStepOn Print from the last recorded slide (its number)
    //! \param excessSpecIndex Provide index of excess species, e.g. N2
    template<class V>
    VType& solve(const V& u0, const std::string& fname = ChooseFDSetting::FILE2READIN, const std::string& outname = std::string(), const std::string& timefile = std::string(), const typename V::value_type* pMolarMasses = 0, int fromLastStepOn = 0, int excessSpecIndex = 0){
      if(fromLastStepOn > 0)
	std::cout << "PRINT SLIDES BEGINNING WITH NO. "<< fromLastStepOn << std::endl;
      std::cout << " SCALING?: ";
      (SCALE) ? FancyMessages().nice_output("yes",32) : FancyMessages().display_in_bold("no");
      
      if((PRINTPRESSURE==true)){
	std::cout << std::endl << "+++++ You've chosen 'PRINT PRESSURE TO FILE' for a ";
	switch(DIM){
	case 0: 
	  std::cout << "0-D mechanism." << std::endl;
	  break;
	case 1:
	  std::cout << "1-D mechanism." << std::endl;
	  break;
	case 2:
	  std::cout << "2-D mechanism." << std::endl;
	  break;
	default:
	  ADONIS_ERROR(DimensionError," Dimension "<< DIM << " not (yet) supported!");
	};
      }
      else{
	std::cout << "No pressure will be printed to file." << std::endl;
      }

       ParameterData PD;
       PD.read_from_file(fname);
       
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
       std::string filename;
       if(outname.size() == 0)
	 filename = PD.get_datum<std::string>("filename");
       else
	 filename = outname;

       //reassign t0 and u0 if you wanna compute from the last slide and time
       //step, respectively
       if(fromLastStepOn > 0){
	 t0 = PD.get_datum<BaseType>("t0_last");
	 std::cout << std::endl << "NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE"<<std::endl;
	 std::cout << "STARTING TIME IS NOW t0_last = "<< t0 << std::endl;
	 std::cout << "LAST READ IN SLIDE FILE: "<< PD.get_datum<std::string>("lastSlide") << std::endl << ", i.e. next slide to be printed carries no. "<< fromLastStepOn << std::endl;
	 std::cout << "NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE · NOTE"<<std::endl << std::endl;
       }

       typedef ODE<T,DISCRFUN,DIM,SPARSE,SCALE,PATTERNCREATIONVIAAD,REPAIR,REPAIRTYPE,PRINTPRESSURE> ModelType; 
       typedef typename ModelType::IndexVecType IndexVecType;
       IndexVecType selectVariables; //use when only some variables are to be
       // written to file for printing 

       ModelType ivp(t0,tend,u0);
       CPUTime<BaseType> cpu;    //also measure CPU time
       cpu.start();
       solution_ = ivp.implicit_euler(maxsteps,maxNewt,ntol,atol,filename,hEstim,Cscal,hmin,hmax,selectVariables,fname,printprec,pMolarMasses,fromLastStepOn,excessSpecIndex);
       elapsedTime_ = cpu.stop(true,timefile);
       time_in_human_readable_form(elapsedTime_);

       return solution_;
    }


    const BaseType& elapsed_time() const {return elapsedTime_;}

  private:
    BaseType elapsedTime_;
    VType solution_;
  };

}

#endif
