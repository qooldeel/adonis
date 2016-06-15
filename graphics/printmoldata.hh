#ifndef PRINT_MOL_DATA_HH
#define PRINT_MOL_DATA_HH

#include <iostream> 
#include <iomanip>
#include <string>

#include "outputfeatures.hh"
#include "printsettings.hh"

#include "../io/readinparameters.hh"
#include "../common/adonisassert.hh"
#include "../common/error.hh"
#include "../common/fancymessages.hh"
#include "../common/globalfunctions.hh"
#include "../common/tempusfugit.hh"
#include "../common/typeadapter.hh"
#include "printer.hh"
#include "printpressuremol.hh"

#include "../fdm/gensettings.hh"

namespace Adonis{
  
  /**
   * CAUTION: the 'get_file' member assumes a file, that provides fields with
   *          Nx, a, b, ....
   *
   *  \tparam D dimension of method of lines (0 = no MOL, 1 = 1D, 2 = 2D) 
   *  \tparam  T type precision
   *
   * Note that the class <TT> Printer </TT> is already endowed with a precision
   *  output via <TT> std::setprecision </TT> through the prec argument
   */
  template<int D, class T, bool PRINTPRESSURE = false> class PrintMOLData;

  //! specialisations
  template<class T, bool PRINTPRESSURE>
  class PrintMOLData<0,T,PRINTPRESSURE>{   //NO method of lines 
  public:
    PrintMOLData(){}
    
    void init(const std::string& molSettingFile, const std::string& outfilename, PrintSolution& ps, int prec){
      ps.initialize(outfilename,prec);
    }

    
    template<class SOL, class SELECTVARS>
    inline void write_2_file(const T& h, unsigned timestep, const T& time, const SOL& y, const std::string& outfilename, PrintSolution& ps, const SELECTVARS& selectVars, const typename SOL::value_type* molarmasses = 0, const typename TypeAdapter<T>::BaseType& addCPUTime = typename TypeAdapter<T>::BaseType(), int fromLastStepOn = 0, const T& kthStep = T()){ // several dummy arguments are used for this specialization
      ps.write_2_file(time,y,selectVars);
    }

    template<class SOL>
    std::size_t num_of_quants(const SOL& y) const {return y.size();}
  };
  
  //1D method of lines
  template<class T, bool PRINTPRESSURE>
  class PrintMOLData<1,T,PRINTPRESSURE>{
  public:
    PrintMOLData():Nx_(0),a_(T()), b_(T()), hx_(T()),nthstep_(0),savestep_(0),saveinterval_(T()),gotFile_(false),isDone_(false){}
    
    void init(const std::string& molSettingFile, const std::string& outfilename, PrintSolution& ps, int prec){ //the 2 mid arguments are meaningless here
#ifndef	NDEBUG
      if(molSettingFile.size() == 0)
	ADONIS_ERROR(IOError,"Settings file "<< molSettingFile << " is empty.\n   Illegal coz u have chosen 1D METHOD OF LINES option.");
#endif
      
      ps.change_writing_precision(prec);

      if(!gotFile_){
	ParameterData Param;
	Param.read_from_file(molSettingFile);
	Nx_ = Param.get_datum<int>("Nx");
	a_ = Param.get_datum<T>("a");
	b_ = Param.get_datum<T>("b");
	hx_ = (b_-a_)/
#ifndef GHOST_POINTS_INCLUDED
	  (Nx_-1)
#else
	  (Nx_-3)
#endif
	  ;
	
	nthstep_ = Param.get_datum<int>("nthstep");
	adonis_assert(nthstep_ >= 1); //you can print at least every step
	//!primacy of 'nthstep_' over saveinterval stuff
	if(nthstep_ == 0){
	  savestep_ = Param.get_datum<int>("savestep");
	  saveinterval_ = Param.get_datum<T>("saveinterval");
	}
      }
      gotFile_ = true;
    }
    
    
     //!Assumptions: 
    //!  1.) species are stored consecutively, each being discretized in space
    //!  2.) space is \f$ \overline{\Omega} = [0,b_{\mathrm{right}}]\f$
    template<class SOL, class SELECTVARS>
    inline void write_2_file(const T& h, unsigned timestep, const T& time, const SOL& y, const std::string& outfilename, PrintSolution& ps, const SELECTVARS& selectVars, const typename SOL::value_type* molarmasses = 0, const typename TypeAdapter<T>::BaseType& addCPUTime = typename TypeAdapter<T>::BaseType(), int fromLastStepOn = 0,const T& kthStep = T()){
      if((timestep == 0) || (isDone_ == false)){  //only assign fromLastStepOn at first timestep
	if(fromLastStepOn > 0)
	  savestep_ = fromLastStepOn;
	isDone_ = true;
      }
      //! primacy of 'nthstep' over everything else
      //if(nthstep_ != 0){
      if(!is_zero(Abs(kthStep))) //reset nthstep_ when kthStep is used
	nthstep_ = 1;      
      if(timestep%nthstep_ == 0){
	//std::cout << "timestep: " << timestep << "  print 1D solution" << std::endl;
	str_ = "timestep: "+ my_function_collection::num2str(timestep)+"  print 1D solution";
	FancyMessages().nice_output(str_,33);
	create_output_format(h,timestep,time,y,outfilename,ps,selectVars,molarmasses,addCPUTime);
      }
      //}
      // else{ //o.k. you want output based on savestep*saveinterval
      // 	if(time > savestep_*saveinterval_)
      // 	  create_output_format(time,y,outfilename,ps,selectVars);
      // } 
    }

    template<class SOL>
    std::size_t num_of_quants(const SOL& y) const {
      adonis_assert(Nx_ > 0);
      return y.size()/Nx_;
    }
    
  private:
    size_t Nx_;
    T a_, b_,
      hx_;
    int nthstep_, savestep_;
    T saveinterval_;
    mutable bool gotFile_;
    bool isDone_;
    std::string str_;


    template<class SOL, class SELECTVARS>
    void create_output_format(const T& h, unsigned timestep, const T& time, const SOL& y, const std::string& outfilename, PrintSolution& ps, const SELECTVARS& selectVars, const typename SOL::value_type* molarmasses = 0, const typename TypeAdapter<T>::BaseType& addCPUTime = typename TypeAdapter<T>::BaseType()){
      ps.open_file(outfilename+my_function_collection::num2str(savestep_)+".gnu",ps.get_precision());
      ps.get_ofstream() << "# TIME: "<< std::setprecision(ps.get_precision()) << time <<
#ifdef GHOST_POINTS_INCLUDED
	"                #GHOST POINTS USED (not printed!)" << 
#endif
	std::endl;
      ps.get_ofstream() << "# time step: "<< timestep << std::endl;
      ps.get_ofstream() << "# last stepsize: "<< h << std::endl;
      if(!is_zero(addCPUTime)){
	ps.get_ofstream() << std::setprecision(ps.get_precision()) << "# Elapsed CPU time so far: "<< time_in_human_readable_form(addCPUTime) << std::endl;
	ps.get_ofstream() << std::setprecision(ps.get_precision()) << "#   (i.e. "<< addCPUTime << " secs)" << std::endl;
      }
	size_t  numberOfVariables = y.size()/Nx_;
	T x = T();
	// create a header (only comment)
	ps.get_ofstream() << "#  x     Variables:..." << std::endl;

	ps.get_ofstream() << std::endl;
#ifndef GHOST_POINTS_INCLUDED
	for(size_t i = 0; i < Nx_; ++i)
#else //output without ghost points
	for(size_t i = 1; i < Nx_-1; ++i) 
#endif
	  {
	  x =
#ifndef GHOST_POINTS_INCLUDED
	    i
#else
	    (i-1)
#endif
	    *hx_;
	  ps.get_ofstream() << std::setprecision(ps.get_precision()) << " " << x << "   ";
	  if(selectVars.size() != 0){
	    for(size_t r = 0; r < selectVars.size(); ++r){
	      ps.get_ofstream() << std::setprecision(ps.get_precision()) << y[selectVars[r]*Nx_+i] << "   ";
	    }
	  }
	  else{
	    for(size_t l = 0; l < numberOfVariables; ++l){
	      ps.get_ofstream() << std::setprecision(ps.get_precision()) << y[l*Nx_+i] << "  ";
	    }
	  }
	  ps.get_ofstream() << std::endl;
	} //end for
		
	ps.close_file();
	++savestep_;        //!increment savestep   
    }

  };
  

  //! 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D
  //! 2D method of lines
  //! 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D
  template<class T, bool PRINTPRESSURE>
  class PrintMOLData<2,T,PRINTPRESSURE>{
  public:
    PrintMOLData():Nx_(0),Ny_(0),npt_(0),a_(T()), b_(T()), c_(T()), d_(T()), hx_(T()), hy_(T()),nthstep_(0),savestep_(0),saveinterval_(T()),gotFile_(false),isDone_(false),numOfSlide_(0){}
    
    void init(const std::string& molSettingFile, const std::string& outfilename, PrintSolution& ps, int prec){ //the 2 middle arguments are meaningless here
#ifndef	NDEBUG
      if(molSettingFile.size() == 0)
	ADONIS_ERROR(IOError,"Settings file "<< molSettingFile << " is empty.\n   Illegal coz u have chosen 2-DIMENSIONAL METHOD OF LINES option.");
#endif
      ps.change_writing_precision(prec);
      //std::cout << "PRECISION: "<< ps.get_precision() << std::endl;

      if(!gotFile_){
	ParameterData Param;
	Param.read_from_file(molSettingFile);
	Nx_ = Param.get_datum<int>("Nx");
	Ny_ = Param.get_datum<int>("Ny");
	npt_ = Nx_*Ny_;
	a_ = Param.get_datum<T>("a");
	b_ = Param.get_datum<T>("b");
	c_ = Param.get_datum<T>("c");
	d_ = Param.get_datum<T>("d");
	hx_ = (b_-a_)/
#ifndef GHOST_POINTS_INCLUDED
	  (Nx_-1)
#else
	  (Nx_-3)
#endif
	  ;
	hy_ = (d_-c_)/
#ifndef GHOST_POINTS_INCLUDED
	  (Ny_-1)
#else
	  (Ny_-3)
#endif
	  ;
	nthstep_ = Param.get_datum<int>("nthstep");
	adonis_assert(nthstep_ >= 1); //you can print at least every step
	//!primacy of 'nthstep_' over saveinterval stuff
	if(nthstep_ == 0){
	  savestep_ = Param.get_datum<int>("savestep");
	  saveinterval_ = Param.get_datum<T>("saveinterval");
	}
      }
      gotFile_ = true;
      offs_ =
#ifndef GHOST_POINTS_INCLUDED
	0
#else
	1
#endif
	;
    }


    //!Assumptions: 
    //!  1.) species are stored consecutively, each being discretized in space
    //!  2.) space is \f$ \overline{\Omega} = [0,b_{\mathrm{right}}]\f$
    template<class SOL, class SELECTVARS>
    inline void write_2_file(const T& h, unsigned timestep, const T& time, const SOL& y, const std::string& outfilename, PrintSolution& ps, const SELECTVARS& selectVars, const typename SOL::value_type* molarmasses = 0, const typename TypeAdapter<T>::BaseType& addCPUTime = typename TypeAdapter<T>::BaseType(), int fromLastStepOn = 0, const T& kthStep = T()){
      if((timestep == 0) || (isDone_ == false)){  //only assign fromLastStepOn at first timestep or if it hasn't assigned so far
	if(fromLastStepOn > 0){
	  savestep_ = fromLastStepOn;
	  numOfSlide_ = fromLastStepOn;
	}
	isDone_ = true;
      }
      //! primacy of 'nthstep' over everything else
      //if(nthstep_ != 0){
      if(!is_zero(Abs(kthStep))) //reset nthstep_ when kthStep is used
	nthstep_ = 1; 
      if(timestep%nthstep_ == 0){
	//std::cout << "timestep: " << timestep << "  print 2D solution" << std::endl;
	numOfSlide_++; //increase slide counter

	str_ = "timestep: "+my_function_collection::num2str(timestep)+"  print 2D solution. Creating slide no. "+Num2str(numOfSlide_);
	FancyMessages().nice_output(str_,33);
	create_output_format(h,timestep,time,y,outfilename,ps,selectVars,molarmasses,addCPUTime);
     
      }
      //}
      // else{ //o.k. you want output based on savestep*saveinterval
      // 	if(time > savestep_*saveinterval_)
      // 	  create_output_format(time,y,outfilename,ps,selectVars);
      // } 
    }

    template<class SOL>
    std::size_t num_of_quants(const SOL& y) const {
      adonis_assert((Nx_ > 0) && (Ny_ > 0));
      return y.size()/(Nx_*Ny_);
    }

  private:
    size_t Nx_, Ny_, npt_;
    T a_, b_,
     c_, d_,
      hx_, hy_;
    int nthstep_, //nth time step 
         savestep_;
    T saveinterval_;
    mutable bool gotFile_;
    bool isDone_;
    unsigned numOfSlide_;
    std::string str_;
    int offs_;

    template<class SOL, class SELECTVARS>
    void create_output_format(const T& h, unsigned timestep, const T& time, const SOL& y, const std::string& outfilename, PrintSolution& ps, const SELECTVARS& selectVars, const typename SOL::value_type* molarmasses = 0, const typename TypeAdapter<T>::BaseType& addCPUTime = typename TypeAdapter<T>::BaseType()){
      ps.open_file(outfilename+my_function_collection::num2str(savestep_)+".gnu",ps.get_precision());
      //std::cout << "Print file "<< outfilename << "..."<<std::endl;
      ps.get_ofstream() << "# TIME: "<< std::setprecision(ps.get_precision()) << time <<
#ifdef GHOST_POINTS_INCLUDED
	"                #GHOST POINTS USED (not printed!)" << 
#endif
	std::endl;
      ps.get_ofstream() << "# time step: " << timestep << std::endl;
      ps.get_ofstream() << "# last step size: "<< h << std::endl;
      if(!is_zero(addCPUTime)){
	ps.get_ofstream() << std::setprecision(ps.get_precision()) << "# Elapsed CPU time so far: "<< time_in_human_readable_form(addCPUTime) << std::endl;
	ps.get_ofstream() << std::setprecision(ps.get_precision()) << "#   (i.e. "<< addCPUTime << " secs)" << std::endl;
      }
      

	size_t  numberOfVariables = y.size()/npt_;
	T x1 = T(),
	  x2 = T();
	ps.get_ofstream() << "#  x     y      Variables:...  " << ((PRINTPRESSURE==true)? "pressure" : "")  << std::endl;
#ifndef GHOST_POINTS_INCLUDED
	for(size_t i = 0; i < Nx_; ++i)
#else  //output without ghost points
	for(size_t i = 1; i < Nx_-1; ++i)  
#endif
	  {
	    x1 = (i-offs_)*hx_;
#ifndef GHOST_POINTS_INCLUDED
	  for(size_t j = 0; j < Ny_; ++j)
#else //output without ghost points
	  for(size_t j = 1; j < Ny_-1; ++j)  
#endif
	    {
	      x2 = (j-offs_)*hy_;
	    
	    ps.get_ofstream() << std::setprecision(ps.get_precision()) << " " << x1 << "   " << x2 << "   ";
	    if(selectVars.size() != 0){  //! full mode; select only spec. vars
	      for(size_t r = 0; r < selectVars.size(); ++r){
		ps.get_ofstream() << std::setprecision(ps.get_precision()) << 
#ifdef NONCONSERVATIVE_FORM
		  y[selectVars[r]*npt_+ i +j*Nx_] 
#else //conservation; devision by rho. ASSUMPTION: selectVars doesn't contain rho
		  ( (selectVars[r] == 0) ? (y[i +j*Nx_]) : ( (y[selectVars[r]*npt_+ i +j*Nx_]/y[i +j*Nx_])) )
#endif
<< "   ";
	      }
	    }
	    else{
	      for(size_t l = 0; l < numberOfVariables; ++l){
		ps.get_ofstream() << std::setprecision(ps.get_precision()) << 
#ifdef NONCONSERVATIVE_FORM
		  y[l*npt_ + i + j*Nx_] 
#else  //conservation; devision by rho (except for l = 0 (this is rho itself)
		  ((l == 0) ? (y[l*npt_ + i + j*Nx_] ) : ((y[l*npt_ + i + j*Nx_] )/(y[i + j*Nx_] )))		  
#endif
<< "   ";
	      }
	    }
	    //!print pressure in addition or not
	    RecordPressure4MOL<PRINTPRESSURE>::two_D(i,j,Nx_,ps,y,npt_,molarmasses);

	    ps.get_ofstream() << std::endl;
	  }  //end for j
	  
	  //! add a white line or not  
	  WhiteLine<PrintSettings::printWhiteLine>::write(ps.get_ofstream());	 	}  //end for i

	ps.close_file();
	++savestep_;        //!increment savestep  
    }

  };
  

} //end namespace

#endif
