#include "../common/globalfunctions.hh"

//! test stuff from this directory
#include "readinparameters.hh"
#include "readinmatrix.hh"

//! test FDM specific stuff for io here
#include "../fdm/gensettings.hh"

using namespace std;
using namespace Adonis;

//! main program
int main(){
  cout << "Test io interfaces for good:"<<endl;


   cout << "Test function 'show_2D_MOL_boundaries'"<< endl;
     ParameterData Param;
     Param.read_from_file("../fdm/data/conaireh2.dat");

     int Nx = Param.get_datum<int>("Nx"),
       Ny = Param.get_datum<int>("Ny");
     
     double a = Param.get_datum<double>("a"),
       b = Param.get_datum<double>("b"),
       c = Param.get_datum<double>("c"),
       d = Param.get_datum<double>("d");

     double hx = (b-a)/
#ifndef GHOST_POINTS_INCLUDED
	(Nx-1)
#else
	(Ny-3)
#endif
       ;
     double hy = (d-c)/
#ifndef GHOST_POINTS_INCLUDED
	(Ny-1)
#else
	(Ny-3)
#endif
       ;



     const int quantity = 3;

     //================ TODO =================================================
     //================ alter number and pressure if necessary ===============
     //=======================================================================
     int numberOfQuantities = 4 + 10 + 1;  //rho,v1,v1,T, + chem specs + pressure
     string testdir = "testfiles/",
       datafile = "CONAIRE_H2_2D_FULL190.gnu"; //5//10//20//50 //190 //280 //338
     string slide =  testdir+datafile;
     //=======================================================================
     cout << "Nx (from file) = "<< Nx << "   Ny (from file) = "<< Ny <<endl;
     //correct hx, hy for ghostpoints
     cout << "hx = "<< hx << "   hy = "<< hy << endl;

     cout << endl << " Extract slide stored in file '"<<slide <<"'"<<endl;
     ReadInMatrix<double> RIM;
     RIM.read_from_file(slide);
     cout << "rows = " << RIM.rows() << "   cols = "<< RIM.cols() << endl;

     RIM.write_2_file("Test_read_in_file.dat");
     
     RIM.show_2D_MOL_boundaries(Nx,Ny,hx,hy,quantity,numberOfQuantities);

     cout << "Temperature:"<<endl;
     RIM.show_min_max_value_2D_MOL(Nx,Ny,quantity,numberOfQuantities);
     cout<< "rho:"<<endl;
     RIM.show_min_max_value_2D_MOL(Nx,Ny,0,numberOfQuantities);
     
} //end main
  
  

