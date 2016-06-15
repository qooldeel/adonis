#ifndef THIS_CAN_BE_USED_2_PRINT_OUT_THE_VALUES_OF_AN_ODE_SOLVER_HH
#define THIS_CAN_BE_USED_2_PRINT_OUT_THE_VALUES_OF_AN_ODE_SOLVER_HH

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

#include "../common/adonisassert.hh"
#include "../common/isclass.hh"
#include "../marcvecmatrix/myfunctions.h"

#include "gnuplot.hh"

namespace Adonis{

  /**
   * \brief TMP to decide whether the solution is a vector container or a type
   */
  template<class X, class OFF, bool B> class IsSolutionAContainer;

  template<class X, class OFF>
  class IsSolutionAContainer<X,OFF,true>{
  public:
    static inline void write(OFF& of, const X& yn, int prec){
      for(typename X::const_iterator it = yn.begin(); it != yn.end(); ++it)
	of << std::setprecision(prec) << *it << "  ";
      of << std::endl;
    } 
  };

  //! is a type
  template<class X, class OFF>
  class IsSolutionAContainer<X,OFF,false>{
  public:
    static inline void write(OFF& of, const X& yn, int prec){
      of << std::setprecision(prec)  << yn << std::endl;
    }
  };
  

  /**
   * \brief print solution of ode solver.
   *
   * USAGE:
   * \code
      double time = 0.234;
      std::vector<double> v(4);
      for(unsigned i = 0; i < 4; ++i) v[i] = (i+1)*0.5;
      
      PrintSolution PS("outputfilename.dat");
    
      for(int i = 0; i < 20; ++i){
        time += std::sin(3.1415*time + 0.5);
        for(int j = 0; j < 4; ++j)
	  v[j] += std::cos(2*(1.5*j)*0.45*3.1415);
        PS.write_2_file(time,v);
      }
      PS.plot(2);  //plot time against species 1 (line 2)
   * \endcode
   */
  class PrintSolution{
  public:
    typedef std::ofstream OfStreamType;
    typedef std::string StringType;
   
    enum{defprec = 7};   //! default precision for output, when now precision
                         //! is provided

    //default constructor
    PrintSolution():prec_(defprec),hasWritten2File_(false){}

    PrintSolution(const StringType& name,int prec = defprec):prec_(prec),name_(name),hasWritten2File_(false){
      of_.open(name_.c_str(), std::ios_base::out);
    }
    
    //! you can leave a comment at the very first line, e.g. "time     species"
    PrintSolution(const StringType& name, const StringType& comment, int prec = defprec):prec_(prec),name_(name),hasWritten2File_(false){
      of_.open(name_.c_str(), std::ios_base::out);
      of_ << "# "<<comment<< std::endl;
    }

    ~PrintSolution(){of_.close();}
  

    void initialize(const StringType& fname, int prec = defprec){
      prec_ = prec;
      of_.open(fname.c_str(), std::ios_base::out);
      name_ = fname;
    }

    void change_writing_precision(int prec){
      prec_ = prec;
    }

    const int get_precision() const {return prec_;} 

    //! alternative description of an intialization
    void open_file(const StringType& fname, int prec){
      prec_ = prec;
      of_.open(fname.c_str(), std::ios_base::out);
      name_ = fname;
    }
    void close_file(){
      of_.close();
    }

    const StringType& name() const{return name_;}
    StringType& name() {return name_;}

    //! o.k., gimme back the ofstream
    OfStreamType& get_ofstream() {return of_;}

    const bool has_written_2_file() const {return hasWritten2File_;}

    template<class V, class T>
    inline void write_2_file(const T& time, const V& yn){
      of_ << " " << time <<"    ";
      
      IsSolutionAContainer<V,OfStreamType,IsContainer<V>::Value>::write(of_,yn,prec_);
     
      hasWritten2File_ = true; //has been written to file
    }

    
    //!overload function -- only print selected values.
    //!NOTE: Only possible for containers with operator[]
    template<class V, class W, class T> 
    inline void write_2_file(const T& time, const V& yn, const W& index){
      adonis_assert(IsContainer<V>::Value); 
      adonis_assert(std::distance(index.begin(),index.end()) <= std::distance(yn.begin(),yn.end()));

      int iDim = static_cast<int>(std::distance(index.begin(),index.end()));
      
      if(iDim != 0){
	of_ << " " << std::setprecision(prec_) << time <<"    ";
     
	for(int i = 0; i < iDim; ++i){
	  adonis_assert(static_cast<int>(index[i]) < static_cast<int>(yn.size()));
	  of_ << std::setprecision(prec_) << yn[index[i]] << "  ";
	}
	of_ << std::endl;
      }
      else
	(*this).write_2_file(time,yn);
      
      hasWritten2File_ = true;
    }

    template<class T, class INT>
    void write_times(INT i, const T& time, const T& kn){
      of_ << "Timestep "<< i << ".)     t = "<< std::setprecision(prec_) << time << "   dt = " << kn << std::endl;
      hasWritten2File_ = true;
    }

    
    //!plot two species against each other, not necessarily time 
    template<class INT> 
    inline void plot(const INT& i, const INT j, const StringType& add2title = StringType(), const StringType& precomms = StringType(), const StringType& suffixcomms = StringType()) const{
      adonis_assert(hasWritten2File_);
      adonis_assert(i != j);  //do not plot the same a column agains itself ;)
      
      //set labels for x - and y-axes, respectively
      StringType xlab = "y_" + my_function_collection::num2str(i),
	ylab = "y_" + my_function_collection::num2str(j);
      if(i == 1){
	xlab = "time t";
	ylab = "y_" + my_function_collection::num2str(j-1);
      }

      std::cout<< "\nSOLUTION OF ODE SOLVER WILL BE PLOTTED AUTOMATICALLY, PAL..."<<std::endl;
      GNUplot GP;
      GP(precomms+"set key box left; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \"" + xlab + "\"; set ylabel \""+ ylab +"\"; plot \""+name_+"\" using "+my_function_collection::num2str(i)+":"+ my_function_collection::num2str(j)+ " with lines lw 3 title \" ""  "+add2title+ "\""+suffixcomms);
    }



    //! print all species in one plot -- This is only reasonable for a smaller number of species :) 
    template<class INT>
    void plot_all(INT dimension, const StringType& add2title = StringType(), const StringType& precomms = StringType(), const StringType& suffixcomms = StringType()) const{
      
      StringType ylab = "U(t)",
	xlab = "time t";

      if(dimension == 1){
	ylab = "U_1(t)";
      }
    
      std::cout<< "\nSOLUTION OF ODE SOLVER WILL BE PLOTTED AUTOMATICALLY, PAL..."<<std::endl;
      GNUplot GP;
      //! no species will be titled inside the plot for the sake of better reading
      GP(precomms+"set title \""+add2title+"\"; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \"" + xlab + "\"; set ylabel \""+ ylab +"\"; plot for[i=2:" + my_function_collection::num2str(dimension+1)+ "] \""+name_+"\" using 1:i with lines  lw 3 notitle" +suffixcomms);

    }

  private:
    int prec_;
    OfStreamType of_;
    StringType name_;
    
    mutable bool hasWritten2File_;

    //PrintSolution(){} //no default constructor
  };

  



  /**
   * \brief Template metaprogram that can be used if you want to invoke 'PrintSolution''s member function or not, thus saving time.
   */
  template<bool B> class PrintMeOrNot; //!forward declaration

  template<>
  class PrintMeOrNot<true>{
  public:
    typedef std::string StringType;

    template<class T, class Y>
    static inline void w_2_f(PrintSolution& ps, const T& time, const Y& y){
      ps.write_2_file(time,y);
    }
    
    //! MOL variant
    template<class INT, class T, class Y, class PRINTMOL, class IXVEC>
    static inline void w_2_f(PRINTMOL& pmd, const T& h, INT i, const T& time, const Y& y, const std::string& outfilename, PrintSolution& ps, const IXVEC& slcvars){
      pmd.write_2_file(h,i,time,y,outfilename,ps,slcvars);
    }

    //! plot i against j
    template<class INT>
    static inline void graphics(const PrintSolution& ps, INT i, INT j, const StringType& addSomething = StringType(), const StringType& precomms = StringType(), const StringType& suffixcomms = StringType()){
      ps.plot(i,j,addSomething,precomms,suffixcomms);
    }
  
    
    template<class INT>
    static inline void all_graphics(const PrintSolution& ps, INT dimension,const StringType& add2title = StringType(), const StringType& precomms = StringType(), const StringType& suffixcomms = StringType()){
      ps.plot_all(dimension,add2title,precomms,suffixcomms);
    }
  };
  

  //! if you don't intend to print and plot something, then nothing will be done
  template<>
  class PrintMeOrNot<false>{
  public:
    typedef std::string StringType;

    template<class T, class Y>
    static inline void w_2_f(PrintSolution& ps, const T& time, const Y& y){}
    
    //! MOL variant - do nothing
    template<class INT, class T, class Y, class PRINTMOL, class IXVEC>
    static inline void w_2_f(PRINTMOL& pmd, const T& h, INT i, const T& time, const Y& y, const std::string& outfilename, PrintSolution& ps, const IXVEC& slcvars){}

    //! plot i against j
    template<class INT>
    static inline void graphics(const PrintSolution& ps, INT i, INT j, const StringType& addSomething = StringType(), const StringType& precomms = StringType(), const StringType& suffixcomms = StringType()){
      std::cout << "NOTHING'S BEEN WRITTEN TO FILE, NOTHING WILL BE PLOTTED." << std::endl;
    }
 
    template<class INT>
    static inline void all_graphics(const PrintSolution& ps, INT dimension,const StringType& add2title = StringType(), const StringType& precomms = StringType(), const StringType& suffixcomms = StringType()){
      std::cout << "NOTHING'S BEEN WRITTEN TO FILE, NOTHING WILL BE PLOTTED." << std::endl;
    }
  };


} //end namespace 

#endif
