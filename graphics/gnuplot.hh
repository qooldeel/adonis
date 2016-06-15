#ifndef GNU_PLOT
#define GNU_PLOT

#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "../common/adonisassert.hh"

namespace Adonis{
  
  
  class Gnuplotfile{
  private:
    std::string myfile_;   //a file you wanna use to execute gnuplot

  public:
    Gnuplotfile(const std::string& s = std::string()): myfile_(s){}
    
    std::string file() const;
  };


  inline std::string Gnuplotfile::file() const{
    return myfile_;
  }



  /** 
      \brief A very simple file, say 'marc.dat', for customising some commands may be
   ___________________________________________________ 
     #My file 4 customizing frequently used commands
     # a hash sign denotes commands in gnuplot

     set xrange[0 : 7]
     set yrange[-2 : 2.5]

     plot sin(2.*x)/x
   ___________________________________________________
  
   class 'Gnuplotfile' describes an object meant for storing the name of such a file.
   class 'GNUplot' is meant for actual plotting -- either by specifying directly the command or by naming the file that contains the plot commands.

   Example:
   ---------

   Gnuplotfile Afile("marc.dat");

   GNUplot GP;
   
   //either
   GP("plot sin(2.*x)/x");           //use command and function operator
   
   //or 
   GP.plot_from_file(Afile);         //plot from file (arg.: object Gnuplotfile)

   //or
   GP.plot_from_file("marc.dat");    //plot from file (arg.: name of file)
   
  */
  class GNUplot{
  private:
    FILE *gnuplotpipe;  //FILE variable pointer
  public:
    GNUplot() throw(std::string);
    ~GNUplot();
    void operator()(const std::string& command); //send any command to gnuplot
    void plot_from_file(const Gnuplotfile& gnuf);
    void plot_from_file(const std::string& f);
  };



  //constructor
  GNUplot::GNUplot() throw(std::string){ 
    gnuplotpipe=popen("gnuplot -persist", "w");
    if(!gnuplotpipe){ 
      ADONIS_ERROR(GnuplotError,"I couldn't find gnuplot. \n Check if the file is in any directory I have access to.");
    }  
  }

  //destructor
  GNUplot::~GNUplot(){ 
    fprintf(gnuplotpipe,"exit\n");
    pclose(gnuplotpipe);
  }
  
  //invoke some gnuplot command
  inline void GNUplot::operator()(const std::string& command){
    fprintf(gnuplotpipe,"%s\n",command.c_str());
    fflush(gnuplotpipe); //necessary for plotting
  }
  
  //plot file based on object 'Gnuplotfile'...
  inline void GNUplot::plot_from_file(const Gnuplotfile& gnuf){
    fprintf(gnuplotpipe,"%s\n", ( "load '" + gnuf.file() + "'").c_str());
    fflush(gnuplotpipe);
  }
 
  //..or directly from filename 
  inline void GNUplot::plot_from_file(const std::string& f){
    fprintf(gnuplotpipe,"%s\n", ( "load '" + f + "'").c_str());
    fflush(gnuplotpipe);
  }
 

} //end namespace

#endif

