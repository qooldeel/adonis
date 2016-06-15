#ifndef MEASURE_ELAPSED_CPU_TIME_HH
#define MEASURE_ELAPSED_CPU_TIME_HH

#include <iostream>
#include <ctime>       //standard time measurement


#ifdef USE_NATIVE_LINUX_TIME_MEASUREMENT  //if defined 
#include <sys/resource.h> //for cpu time, headers for getrusage
#include <sys/types.h>
#include <sys/time.h>     //for wall time
#endif

#include <cmath>
#include <string>
#include <fstream>

#include "numerictypechecker.hh"
#include "adonisassert.hh"

#include "globalfunctions.hh"

namespace Adonis{

  /**
   * \brief When time measurements are longer it is rather beneficial to print 
   *  it in a more human readable form ;-)
   *
   * NOTE: all the following code snippets implement time measuring methods <B> up to seconds</B>. 
   *
   * Example :
   * \code
   //write somewhere in header file:
   #define USE_NATIVE_LINUX_TIME_MEASUREMENT 1

        CPUTime<double> cpu;   
	cpu.start();           //start time measurement
	
	//..... your computation is done here

	cpu.stop();            //o.k. stop time measurement

	//when computation time is very long, then you get a nice output via:

	time_in_human_readable_form(cpu.get_time());

	// you can reuse the object be resetting it again:
	cpu.reset();
	//...etc.
   * \endcode 
   */
  template<class T>
  inline std::string time_in_human_readable_form(const T& tempus, bool additional2screen = false){
    std::string str;

    NumericDataTypeChecker<T>::certify();
    
    bool inMins = false, inHrs = false, inDays = false;
    
    T days = T(), hours = T(), mins = T(), secs = T(),
      tm = tempus; 

    //the case when computations took minutes rather than secs or hours…
    if(tempus >= 60. && tempus < 3600.){
      tm /= 60.;
      mins = std::floor(tm);
      secs = (tm - mins)*60.;
      inMins = true;
    }
    //… and when computations took hours rather than minutes 
    else if(tempus >= 3600. && tempus < 86400.){
      tm /= 3600.;
      hours = std::floor(tm);
      T mn = (tm - hours)*60.;
      mins = std::floor(mn);
      secs = (mn - mins)*60.;
      inHrs = true;
    }
    //… and when computations take days
    else if(tempus >= 86400.){
      tm /= 86400.;
      days = std::floor(tm);
      T hr = (tm - days)*24.;
      hours = std::floor(hr);
      T mn = (hr-hours)*60.;
      mins = std::floor(mn);
      secs = (mn -mins)*60.;
      inDays = true;
    }
  
    if(inDays){
      if(additional2screen){
	std::cout <<days<<" d ";
	if(hours > 0.)
	  std::cout << hours << " h ";
	if(mins > 0.)
	  std::cout << mins << " min ";
	if(secs > 0.) 
	  std::cout << secs << " sec";
      }
      
      str = Num2str(days) + " d " + Num2str(hours) + " h " + Num2str(mins) + " min " +  Num2str(secs) + " secs.";
    }
    else if(inHrs){
      if(additional2screen){
	std::cout << hours << " h ";
	if(mins > 0.)
	  std::cout << mins << " min ";
	if(secs > 0.)
	  std::cout << secs << " sec ";
      }

      str = Num2str(hours) + " h " + Num2str(mins) + " min " +  Num2str(secs) + " secs.";
    }
    else if (inMins){
      if(additional2screen)
	std::cout << mins << " min "<< secs << " sec "; 
      
      str = Num2str(mins) + " min " +  Num2str(secs) + " secs.";
    }
    else{
      if(additional2screen)
	std::cout << tempus << " sec "; 

      str = Num2str(tempus) + " secs.";
    }
    if(additional2screen)
      std::cout << std::endl;
    
    return str;
  }


  /**
   * \brief Measures elapsed CPU time
   *
   *  \tparam T precision in which to output time. NOTE: if set to int and the time is below 1 sec than the output may be zero (unless it is in the millisecond 
   range anyway)
   *\code
     CPUTime<double> cpu;
     cpu.start();
     //...some computations to be done
     cpu.stop();
   *\endcode
   * NOTE (in case of <TT> USE_NATIVE_LINUX_TIME_MEASUREMENT 0 </TT> only): 
   *       Watch out for overflow. If clock_t is a 32-bit signed integer type
   *       and CLOCKS_PER_SEC is 1000000 (1e6), then it will overflow in less
   *       than 36 minutes; assuming signed integer overflow wraps around, you
   *       can get repeated results in less than 72 minutes.
   */
  template<class T>
  class CPUTime{
  private:
#if USE_NATIVE_LINUX_TIME_MEASUREMENT
    struct rusage ruse_;
#else
    //!'clock_t' is typedef-synonym of 'long int' 
    clock_t c1_,     //time stamp before computations
      c2_;           //time stamp after computations
#endif

    T time_;
    
    bool hasStarted_,
      hasStopped_;

  public:
    typedef T value_type;

    CPUTime():time_(T()),hasStarted_(false),hasStopped_(false){
      NumericDataTypeChecker<T>::certify();
      
#if USE_NATIVE_LINUX_TIME_MEASUREMENT 
      //do nothing then
#else
      //std::cout << "use clocks for standard time measurement" << std::endl;
      c1_ = 0;
      c2_ = 0;
#endif
}
    

    inline void reset(){  //!reset time stamps
#if USE_NATIVE_LINUX_TIME_MEASUREMENT      
      //do nothing then
#else
      c1_ = 0;
      c2_ = 0;
#endif
      time_ = T();
      hasStarted_ = false;
      hasStopped_ = false;
    }

    inline void start(bool printtoscreen = true){   //! set time stamp BEFORE computations
      if(!hasStarted_){
#if USE_NATIVE_LINUX_TIME_MEASUREMENT 
	getrusage(RUSAGE_SELF,&ruse_);
	time_ = static_cast<T>(ruse_.ru_utime.tv_sec) + static_cast<T>(ruse_.ru_utime.tv_usec)/1.e+06;
#else
	c1_ = clock();
#endif
	hasStarted_ = true;
	if(printtoscreen)
	  std::cout << " Starting measuring CPU time..." << std::endl;
      }
    }

    inline T& stop(bool printtoscreen = true, const std::string& fnameout = std::string()){   //! set time stamp AFTER computations and calculate
      if(!hasStopped_ && hasStarted_){
#if USE_NATIVE_LINUX_TIME_MEASUREMENT 
	if(printtoscreen)
	  std::cout << " I use native LINUX time measurement via 'getrusage', hunk." << std::endl;
	getrusage(RUSAGE_SELF,&ruse_);
	time_ *= -1;  //! now: time_ = -starttime
	time_ += static_cast<T>(ruse_.ru_utime.tv_sec) + static_cast<T>(ruse_.ru_utime.tv_usec)/1.e+06; //now: time = endtime - starttime 
#else  
	if(printtoscreen)
	  std::cout << " Standard time measurement via <ctime> applied, pal. "<<std::endl;
	c2_ = clock();
	time_ = static_cast<T>(c2_-c1_)/CLOCKS_PER_SEC;
#endif
	if(printtoscreen)
	  std::cout << std::endl <<" Elapsed CPU time = " << time_ << " secs."<<std::endl;
	if(fnameout.size() != 0){
	  std::ofstream of(fnameout.c_str());
	  of << "TIME ELAPSED: " << time_in_human_readable_form(time_) << std::endl;
	  of.close();
	}

	hasStopped_ = true;
      }

      return time_;
    }
  
    inline const T& get_time() const{
      return time_;
    }

  };




  /**
   * \brief Wall time calculation
   */
  template<class T>
  class WallTime{
  private:
#if USE_NATIVE_LINUX_TIME_MEASUREMENT
    struct timeval tp_;
#else
    time_t t1_,
      t2_;
#endif
    T wall_;
    
    bool hasStarted_,
    hasStopped_;
      
  public:
    typedef T value_type;

    WallTime():wall_(T()),hasStarted_(false),hasStopped_(false){
      NumericDataTypeChecker<T>::certify();
      
#if USE_NATIVE_LINUX_TIME_MEASUREMENT
      //do nothing
#else
      t1_ = 0;
      t2_ = 0;
#endif
    }

    inline void reset(){  //!reset time stamps
#if USE_NATIVE_LINUX_TIME_MEASUREMENT      
      //do nothing then
#else
      t1_ = 0;
      t2_ = 0;
#endif
      wall_ = T();
      hasStarted_ = false;
      hasStopped_ = false;
    }


    inline void start(bool printtoscreen = true){   //! starts setting the clock timer
      if(!hasStarted_){
#if USE_NATIVE_LINUX_TIME_MEASUREMENT
	gettimeofday(&tp_,NULL);
	wall_ = static_cast<T>(tp_.tv_sec) + static_cast<T>(tp_.tv_usec)/1.e+06;
#else
	time(&t1_);
#endif  
	if(printtoscreen)
	  std::cout << " Starting measuring Wall time..." << std::endl;
	hasStarted_ = true;
      }
    }
  
    inline T& stop(bool printtoscreen = true, const std::string& fnameout = std::string()){
      if(!hasStopped_ && hasStarted_){
#if USE_NATIVE_LINUX_TIME_MEASUREMENT
	gettimeofday(&tp_,NULL);
	wall_ *= -1.;  
	wall_ += static_cast<T>(tp_.tv_sec) + static_cast<T>(tp_.tv_usec)/1.e+06;
#else
	time(&t2_);
	wall_ = difftime(t2_,t1_);
#endif
	if(printtoscreen)
	  std::cout << std::endl <<" Elapsed Wall time = "<< wall_ << " secs."<<std::endl;

	if(fnameout.size() != 0){
	  std::ofstream of(fnameout.c_str());
	  of << "TIME ELAPSED: " << time_in_human_readable_form(wall_) << std::endl;
	  of.close();
	}

	hasStopped_ = true;
      }
      return wall_;
    }
    
    inline const T& get_time() const{
      return wall_;
    }

  };

}

#endif
