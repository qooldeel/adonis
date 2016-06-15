

//!!!!!!!!!! ===== CPPAD USAGE ================================================
#include "../want_2_use.h"  //Should always be put at the very beginning
//==============================================================================

/// · REDUCTION · REDUCTION · REDUCTION · REDUCTION · REDUCTION · REDUCTION ·
//======================== WHICH REDUCTION U WANNA USE ? ====================

#include "which_reduction_to_be_used.h"

//==============================================================================
// · REDUCTION · REDUCTION · REDUCTION · REDUCTION · REDUCTION · REDUCTION ·

#include <iostream>
#include <typeinfo>
#include <string>

#include "../expressiontemplates/exprvec.hh"


//#include "manifold.hh"

#include "../ode/ode.hh"   //for solving odes

//extended Davis Skodie -- NON ANALYTICAL form
// #include "../ode/examples/XDavisSkodie/XDSfull.hh"
// #include "../ode/examples/XDavisSkodie/XDS1st.hh"
// #include "../ode/examples/XDavisSkodie/XDS2nd.hh"

//H2C6 example
#include "../ode/examples/h2c6/h2c61Dim.hh"  //original system
// #include "../ode/examples/h2c6/h2c61Dim_1st.hh"  //1st approx
// #include "../ode/examples/h2c6/h2c61Dim_withcorr.hh"
#include "../ode/examples/zeldovich.hh"


#include "../common/tempusfugit.hh"

#include "../marcvecmatrix/myfunctions.h"
#include "../misc/operations4randomaccesscontainers.hh"

#include "../ode/stiffnesschecker.hh"

//////////// only the sources ////////////////
//#include "../ode/examples/h2c6/h2c61Dim_1st_pure_source.hh"

//////////////////////////////////////////////

//H2C6 example revisited with automatically computed rhs
//!original source term
#include "../ode/examples/automaticallygeneratedsourceterms/h2c6.hh"
//!reduced source term
//#include "../ode/examples/automaticallygeneratedsourceterms/redh2c6.hh"

//// ATTEMPT BY ME /////
#include "../ode/examples/automaticallygeneratedsourceterms/redh2c6_attempt.hh"
//////////////

//H2GRI 
#include "../ode/examples/automaticallygeneratedsourceterms/h2gri.hh"
//#include "../ode/examples/automaticallygeneratedsourceterms/redh2gri.hh"

//!spatiotemporal stuff
//1D case
#include "../ode/examples/automaticallygeneratedsourceterms/1Ddiscr/h2c6Discretization1D.hh"   

#include "../ode/examples/automaticallygeneratedsourceterms/1Ddiscr/redh2c6Discretization1D.hh"
#include "../ode/examples/automaticallygeneratedsourceterms/1Ddiscr/second_appx_redh2c6Discretization1D.hh"

#include "../templatemetaprograms/unrollloop.hh"
#include "startingpoint.hh"

#include "../massactionkinetics/data/thermochemicaldata.hh"

#include "../io/readinparameters.hh"

#include "verysimple.hh"

#include "handreduced.hh"

#include "reducedzeldovich.hh"
#include "2speciesreducedzeldovich.hh"

#include "1dredh2c6.hh"
#include "1dfullh2c6.hh"

using namespace std;
using namespace Adonis;
using namespace ExprTmpl;
using namespace my_function_collection; 





// =========================================
// ========== MAIN PROGRAM =================
//==========================================
int main(int argc, char** argv)
{
  if(argc != 2) //argv[0] = Program name, argv[1] = number of executable
    ADONIS_ERROR(MainProgramError, "\""<< argv[0] <<"\" takes exactely 1 argument, pal.");

  int executable = atoi(argv[1]);
 
  const int ZEROD = 0;
  const int ONED = 1;
  const int TWOD = 2;


  //======================================================================
  //======================== PARAMTERERS =================================
  //======================================================================
  ParameterData PD;
  PD.read_from_file("dat/solversettings.dat");
  double ODE_t0 = PD.get_datum<double>("ODE_t0"),
    ODE_tend = PD.get_datum<double>("ODE_tend");

  size_t ODE_IE_maxsteps = PD.get_datum<size_t>("ODE_IE_maxsteps"),
    ODE_IE_maxNewtoniters = PD.get_datum<size_t>("ODE_IE_maxNewtoniters"),
    ODE_E4S_maxFixIter = PD.get_datum<size_t>("ODE_E4S_maxFixIter");

  double ODE_IE_ntol = PD.get_datum<double>("ODE_IE_ntol"),
    ODE_IE_rtol = PD.get_datum<double>("ODE_IE_rtol"),
    ODE_IE_hEstim = PD.get_datum<double>("ODE_IE_hEstim"),
    ODE_IE_Cscal =  PD.get_datum<double>("ODE_IE_Cscal"),
    ODE_IE_hmin = PD.get_datum<double>("ODE_IE_hmin"),
    ODE_IE_hmax = PD.get_datum<double>("ODE_IE_hmax"),
    
    ODE_E4S_errtol = PD.get_datum<double>("ODE_E4S_errtol"),
    ODE_E4S_tolFixPt = PD.get_datum<double>("ODE_E4S_tolFixPt"),
    ODE_E4S_hEst = PD.get_datum<double>("ODE_E4S_hEst"),
    ODE_E4S_St = PD.get_datum<double>("ODE_E4S_St"),
    ODE_E4S_cs = PD.get_datum<double>("ODE_E4S_cs"),
    ODE_E4S_kmin = PD.get_datum<double>("ODE_E4S_kmin"),
    ODE_E4S_kmax = PD.get_datum<double>("ODE_E4S_kmax");
    
  double  ODE_IE_REDUC_ntol = PD.get_datum<double>("ODE_IE_REDUC_ntol"),
    ODE_IE_REDUC_rtol = PD.get_datum<double>("ODE_IE_REDUC_rtol"),
    ODE_IE_REDUC_hEstim = PD.get_datum<double>("ODE_IE_REDUC_hEstim"),
    ODE_E4S_REDUC_errtol = PD.get_datum<double>("ODE_E4S_REDUC_errtol"),
    ODE_E4S_REDUC_tolFixPt = PD.get_datum<double>("ODE_E4S_REDUC_tolFixPt"),
    ODE_E4S_REDUC_hEst = PD.get_datum<double>("ODE_E4S_REDUC_hEst"),
    ODE_E4S_REDUC_cs = PD.get_datum<double>("ODE_E4S_REDUC_cs");
    
  bool totex = PD.get_datum<bool>("latex");
  string plotsize = " size "+PD.get_datum<string>("width")+PD.get_datum<string>("unit")+","+PD.get_datum<string>("height")+PD.get_datum<string>("unit")+" ";
   
  
    
  //======================================================================

  
  if(executable == 0){ 
    cout<< endl;    cout << "==============================================================="<<endl;
    cout << "======================= H2C6 EXAMPLE =========================="<<endl;
    cout << "==============================================================="<<endl;

   

    enum{hydrogen = 0, water = 4};

    ParameterData AddData;
    AddData.read_from_file("dat/Molsettings.dat");

    const unsigned nos = 6,                 //# of orig. species
      red = 2,                              //# of reduced species   
      //unrep = nos - red,                    //# of unrep. species
      spacepts = AddData.get_datum<unsigned>("Nx"),
      fulldim = nos*spacepts,
      reddim = red * spacepts;
    
    // const unsigned i1 = 0,   //rpv indizes denoting species in full system
    //  i2 = 4;
    
    double t0 = 0.,
      tend =  ODE_tend; //2.;  //0.2
    
    //1D domain
    double x1 = AddData.get_datum<double>("b");

    double hx = x1/(spacepts-1);

    //Note that the SIM is neither dependent on time t nor on space x...
   
    
    MyVec<double> ICf(fulldim);

    
    
    // solution of U'(t) = S(U(t)), t --> inf
    //double ustar[] = {0.269999, 0.0499998, 0.135, 0.0199999, 0.700001, 0.00999997 };


    //EQUILIBRIUM
    MyVec<double> ustar(6);
    ustar[0] = 0.27;      // H2
    ustar[1] = 0.05;
    ustar[2] = 0.135;
    ustar[3] = 0.02;
    ustar[4] = 0.7;       // H2O
    ustar[5] = 0.01;
    //! start value
    StartingPoint<MyVec<double> > startpoint(6);
    
    MyVec<double> U0 = startpoint.get_point();
    
    MyVec<double> u0(6), u1(6);
    
    u0 <<= 0.455, 0.779, 0.237, 0.363, 0.148, 0.015;

    u1 <<= 0.2, 0.95, 0.31, 0.03, 0.3, 0.05;

    double x = 0.;
    //initial conditions -- constant at each grid point x_i
   
    MyVec<double> zM0(6);

    for(unsigned i = 0; i < spacepts;++i){
      x = i*hx;  //CAUTION: previously: x += i*hx

      if(0.002 <= x && x <= 0.008){
	for(int k = 0; k < 6; ++k)
	  ICf[k*spacepts+i] = u0[k];
      }
      if(x < 0.002 || x > 0.008){
	for(int k = 0; k < 6; ++k)
	  ICf[k*spacepts+i] = u1[k];
      }
   
      if(i == 0){
	for(int k = 0; k < 6; ++k)
	  zM0[k] = ICf[k*spacepts];
      }
	

      ////! older version
      // if( (x <= 0.002) || (x > 0.008) ){
      // 	ICf[i] = ustar[0];
      // 	ICf[spacepts+i] = ustar[1];
      // 	ICf[2*spacepts+i] = ustar[2];
      // 	ICf[3*spacepts+i] = ustar[3];
      // 	ICf[4*spacepts+i] = ustar[4];
      // 	ICf[5*spacepts+i] = ustar[5];
      // }
      // else{
      // 	ICf[i] = U0[0];
      // 	ICf[spacepts+i] = U0[1];
      // 	ICf[2*spacepts+i] = U0[2];
      // 	ICf[3*spacepts+i] = U0[3];
      // 	ICf[4*spacepts+i] = U0[4];
      // 	ICf[5*spacepts+i] = U0[5];
      // }

      //! old ones
      // ICf[i] = onmani[0];
      // ICf[spacepts+i] = onmani[1];
      // ICf[2*spacepts+i] = onmani[2];
      // ICf[3*spacepts+i] = onmani[3];
      // ICf[4*spacepts+i] = onmani[4];
      // ICf[5*spacepts+i] = onmani[5];
    }
    
    unsigned maxit = 20000,
      newtmaxit = 7; //1000;  //represents infinity
  
    double timefull = 0., time1st = 0.;  // time2nd = 0.;
    
    CPUTime<double> cpu;   //measure elapsed CPU time
  

    string sol1dfullout = "Sol1D_FULL_",
      sol1dredout = "Sol1D_RED_1st_",
      sol1dredCPAout = "Sol1D_RED_CPA_";
   
    
    //============================ original ====================================
    typedef MyVec<int> IndexVecType;
     MyVec<double> Sol_full;

     size_t niterFull = 0,
       niter1st = 0;

#if FULL_MODEL 
    cout << endl << " === ORIGINAL SYSTEM ===" << endl;
    //original spatio-temporal 1D system 
    
    //!AUTO
    typedef ODE<double,H2MechIn6SpeciesSpatial1D,ONED> DiffEqType;
    IndexVecType rpvix(2);
    rpvix[0] = 0; rpvix[1] = 4;

    //!HAND-CODED usually much faster through direct implementation
    //typedef ODE<double, H2MechWithDiffusion1D> DiffEqType;
    
    DiffEqType ivp_full(t0,tend,ICf);
    
    cpu.start();
    
   

#if SOLVER_TO_BE_USED == 'i'
     
      //++++++ TODO: change inputs
      double ntol = 1.e-03, //1.e-04,//1.e-03,
	rtol = 1.e-02;  //1.e-03; //1.e-02;
      
      //=============== estimate of step size =======================
      double hEstim = ODE_IE_hEstim; //0.01; //for equidist. //0.0025;
      //++++++++++++++++++++++++++++

      Sol_full = ivp_full.implicit_euler(maxit,newtmaxit,ntol,rtol,sol1dfullout,hEstim,1.,1.e-12,1.e+03,rpvix,"dat/Molsettings.dat");

      niterFull = ivp_full.number_of_iterations();
#endif


#if SOLVER_TO_BE_USED == 'j'
      //++++++ TODO: change inputs
      double TOL = 1.e-06,
	tol = 1.e-06,
	k1 = 0.0000125,
	S_T = 0.75,
	c_fix = 0.99;
      
      size_t mxit = 7;
      
      Sol_full = ivp_full.explicit_4_stiff(TOL,tol,k1,S_T,c_fix,mxit,"full_EJL.dat");

      niterFull = ivp_full.number_of_iterations();
#endif

    cpu.stop(); //_without_message();
    
    timefull = cpu.get_time();

#endif //end full model     


#if APPLY_ANY_REDUCTION
    //========================= REDUCTIONS =====================================
    //========================= REDUCTIONS =====================================
    //========================= REDUCTIONS =====================================
    
    printf("%c[%d;%d;%dm", 0x1B, 1,33, 40);
    cout << "APPLY REDUCTION "<< REDUCTION_TO_BE_USED << endl;
    printf("%c[%dm", 0x1B, 0);
    
    
    MyVec<int> index(reddim);  //index for spatial 1d case
    for(unsigned i = 0; i < spacepts; ++i){
      index[i] = hydrogen*spacepts + i;
      index[spacepts+i] = water*spacepts+i;
    }

    
    //ATTENTION: with constant values, diffusion might be zero when discretised
    //IC for reduced discription
    MyVec<double> ICr(reddim);
    for(unsigned i = 0; i < spacepts;++i){
      ICr[i] =  ICf[i];                    //rpv index 0
      ICr[spacepts+i] = ICf[4*spacepts+i];      //rpv index 4 of original system
    }
  

    IndexVecType redrpvix;  //zero vector

    //TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    //TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    //TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
 #if SOLVER_TO_BE_USED == 'i'
    //o.k.: spacepts = 5, tol1 = 0.5 (1.e-01), tol2 = 1.e-01 (0.3), hest = 0.005
    
    double tol1 =  1.e-06,//1.e-03, //1.e-04,//1.e-03,     //newton iter
      tol2 = 1.e-06; //1.e-02; //1.e-03;//1.e-02;             //local error
    
    double hest =  ODE_IE_REDUC_hEstim; //0.00012;//0.0025; //0.00175;     
    double Cstszctrl = 1.;  // 1.1111111111; 
    //TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    //TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    //TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#endif
    
#if (REDUCTION_TO_BE_USED == 0 || REDUCTION_TO_BE_USED > 2)
      ADONIS_ERROR(DerivedError,"Reduction method '" << whichReduction <<"' has not been defined so far");
#endif

#if REDUCTION_TO_BE_USED == 1
      cout << endl << " ===== 1ST APPROXIMATION =================" << endl;
   
      //1st approximation of spatio-temporal system. Invoke Jochen's first
      //AUTO
      ODE<double,ReducedH2MechIn6SpeciesSpatial1D,ONED> ivp_1st(t0,tend,ICr);
      //HANDCODED -- usually (much) faster
      //ODE<double,H2MechWithDiffusion1d1STAPPROX>  ivp_1st(t0,tend,ICr); 
      cpu.reset();  //RESET time
      cpu.start();

      MyVec<double> Sol_1st;
      
#if SOLVER_TO_BE_USED == 'i'
      Sol_1st = ivp_1st.implicit_euler(maxit,20,tol1,tol2,sol1dredout,hest,Cstszctrl,1.e-12,0.1,redrpvix,"dat/Molsettings.dat");
      //ivp_1st.adaptive_rkf(maxit, hest, 1.e-04, "rkfspace1st.dat");
      niter1st = ivp_1st.number_of_iterations();
#endif
      
#if SOLVER_TO_BE_USED == 'j'
	double TOL1 = 1.e-05,
	tol = 1.e-05,
	k11 = 0.0000125,
	  S_T1 = 1,
	c_fix1 = 0.99956;
      
      size_t mxit1 = 4;
      
      Sol_1st = ivp_1st.explicit_4_stiff(TOL1,tol,k11,S_T1,c_fix1,mxit1,"1st_EJL.dat");
      niter1st = ivp_1st.number_of_iterations();
#endif
      

      cpu.stop();
    
      time1st = cpu.get_time();
#endif    

  
    
    
    //============================== 2nd =======================================
#if REDUCTION_TO_BE_USED == 2
      cout << endl << " === 2nd APPROXIMATION ===" << endl;
      // ODE<double,H2MechWithDiffusion1dCORRECTION,ONED>  ivp_2nd(t0,tend,ICr);  //errors
      ODE<double,SecondReducedH2MechIn6SpeciesSpatial1D,ONED> ivp_2nd(t0,tend,ICr);
      cpu.reset();  //RESET time
      cpu.start();
      MyVec<double> Sol_2nd;

#if SOLVER_TO_BE_USED == 'i'
      Sol_2nd = ivp_2nd.implicit_euler(maxit,20,tol1,tol2,sol1dredCPAout,hest,Cstszctrl,1.e-12,1.e+03,redrpvix,"dat/Molsettings.dat");

#endif

#if SOLVER_TO_BE_USED == 'j'
	Sol_2nd = ivp_2nd.explicit_4_stiff(TOL1,tol1,k11,S_T1,c_fix1,mxit1,"2nd_EJL.dat");
#endif

      cpu.stop();

      time2nd = cpu.get_time();
#endif

#endif //APPLY_ANY_REDUCTION
    
      //++++++++++++++++++++++++++++ PLOTTING ++++++++++++++++++++++++++++
    //PLOT STUFF HERE
    //take last slide, i.e. t = tend with x \in [0,1]
    std::ofstream myo("Solution.dat",std::ios_base::out);

#if APPLY_ANY_REDUCTION
#if REDUCTION_TO_BE_USED == 1
    std::ofstream first("1stapprox.dat",std::ios_base::out);
#endif
#if REDUCTION_TO_BE_USED == 2
    std::ofstream second("2ndapprox.dat",std::ios_base::out);
#endif
#endif


    string heading = "# space x           H2             H2O";

    myo << heading <<endl;

#if APPLY_ANY_REDUCTION
#if REDUCTION_TO_BE_USED == 1
    first << heading << endl;
#endif
#if REDUCTION_TO_BE_USED == 2
    second << heading << endl;
#endif
#endif

    x = 0;
    for(unsigned i = 0; i < spacepts; ++i){
      x = i*hx;
 #if FULL_MODEL     
      myo << " "<< x << "    "<<Sol_full[i] << "    "<<Sol_full[4*spacepts+i] <<endl;
#endif

#if APPLY_ANY_REDUCTION
#if REDUCTION_TO_BE_USED == 1
      first << " "<< x << "    "<<Sol_1st[i] << "   "<<Sol_1st[spacepts+i] <<endl; 
#endif
#if REDUCTION_TO_BE_USED == 2
      second << " "<< x << "    "<<Sol_2nd[i] << "   "<<Sol_2nd[spacepts+i] <<endl; 
#endif
#endif 
    
    }

  
  
    cout << endl<< "System type    ||     elapsed CPU time [secs]    ||     #iters"<< endl
	 << "--------------------------------------------------------------"<< endl;
#if FULL_MODEL
	 cout << " Original      ||        "<<  timefull << "     ||     "<<niterFull <<endl; 
#endif

#if APPLY_ANY_REDUCTION
#if REDUCTION_TO_BE_USED == 1
      cout << " 1st approx.   ||        "; 
      time_in_human_readable_form(time1st);
      cout << "     ||     "<< niter1st << endl;
#endif
#if REDUCTION_TO_BE_USED == 2
      cout << " 2nd approx.   ||        ";   
      time_in_human_readable_form(time2nd); 
#endif

#if FULL_MODEL //a slow down etc. can only be computed when the full model  
#if REDUCTION_TO_BE_USED == 1            //has been computed previously
      cout<< endl<< "l_2 error (1st) = "<< l_2_error(Sol_full,Sol_1st,index) << endl;
      cout <<"slowdown (1st) = "<< slow_down_rate(time1st,timefull) <<endl;
#endif
#if REDUCTION_TO_BE_USED == 2
      cout<< endl<< "l_2 error (2nd) = "<< l_2_error(Sol_full,Sol_2nd,index) << endl;
      cout <<"slowdown (2nd) = "<< slow_down_rate(time2nd,timefull) <<endl;
#endif
#endif
#endif //full model
    // ============== TODO: name for output files ============================= 
    bool createEpsFile = false;
    
    string myepsfilename = "VGN_DiffusionConst1st.eps";
    //=========================================================================

    string toeps;
      
    
    if(createEpsFile) {
      toeps = "set terminal postscript eps color enhanced; set output";
      toeps +=  ("'" + myepsfilename + "';");
    }
      
    
    cout << "~~~~~ PLOT STUFF"<<endl;
    //plot stuff
    GNUplot wildebeest;
 
#if APPLY_ANY_REDUCTION
#if PLOT_BOTH_REDUCTIONS_SIMULTANEOUSLY 
    wildebeest(toeps +
	       "set key box left; set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"x\"; set ylabel \"z(x,t="+my_function_collection::num2str(tend)+")\"; set title \"H2 and H2O with " + my_function_collection::num2str(spacepts) + " spatial points\";plot \"Solution.dat\" using 1:2 with lines lc rgb \"red\" lw 2 title \"H_2 full\", \"Solution.dat\" using 1:3 with lines lc rgb \"green\" lw 2 title \"H_2O full\", \"1stapprox.dat\" using 1:2 with lines lc rgb \"violet\" lw 2 title \"H_2 1st\", \"1stapprox.dat\" using 1:3 with lines lc rgb \"blue\" lw 2 title \"H_2O 1st\", \"2ndapprox.dat\" using 1:2 with lines lc rgb \"yellow\" lw 2 title \"H_2 2nd\", \"2ndapprox.dat\" using 1:3 with lines lc rgb \"orange\" lw 2 title \"H_2O 2nd\"");
#else
    cout << "Plot only full vs. 1st reduction"<<endl;
 wildebeest(toeps +
	       "set key box left; set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"x\"; set ylabel \"z(x,t="+my_function_collection::num2str(tend)+")\"; set title \"H2 and H2O with " + my_function_collection::num2str(spacepts) + " spatial points\";plot \"Solution.dat\" using 1:2 with lines lc rgb \"red\" lw 2 title \"H_2 full\", \"Solution.dat\" using 1:3 with lines lc rgb \"green\" lw 2 title \"H_2O full\", \"1stapprox.dat\" using 1:2 with lines lc rgb \"violet\" lw 2 title \"H_2 1st\", \"1stapprox.dat\" using 1:3 with lines lc rgb \"blue\" lw 2 title \"H_2O 1st\"");

#endif
#else
    
      wildebeest(toeps +
		 "set key box left; set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"x\"; set ylabel \"z(x,t="+my_function_collection::num2str(tend)+")\"; set title \"H2 and H2O with " + my_function_collection::num2str(spacepts) + " spatial points\";plot \"Solution.dat\" using 1:2 with lines lc rgb \"red\" lw 2 title \"H_2 full\", \"Solution.dat\" using 1:3 with lines lc rgb \"green\" lw 2 title \"H_2O full\"");
#endif //APPLY_ANY_REDUCTION

    myo.close();

#if APPLY_ANY_REDUCTION
#if REDUCTION_TO_BE_USED == 1
    first.close();
#endif
#if REDUCTION_TO_BE_USED == 2
    second.close();
#endif
#endif
    


    cout << "zM0 = "<< zM0 << endl;
  } //end executable








  //==========================================================================
  //================ EXTENDED DAVIS SKODIE ===================================
  //==========================================================================

  else if(executable == 1){
    /* std::cout << "COMPLETELY NUMERICALLY COMPUTED extended Davis Skodie"<< std::endl;
    
    typedef Constant<double> ConstantType;


    const unsigned nos = 2,                 //# of orig. species
      red = 1,                              //# of reduced species   
      unrep = nos - red,                    //# of unrep. species
      spacepts = ConstantType::spPts,
      fulldim = nos*spacepts,
      reddim = red * spacepts;
    
    
    MyVec<int> index(reddim);  //index for spatial 1d case
    for(unsigned i = 0; i < spacepts; ++i){
      index[i] =  i;
    }


    double t0 = 0.,
      tend = 1.;

    Parameter<double,7> Para("Parameter.dat");
    const double a = Para[0],
      hx = 1./(spacepts-1);

    
    MyVec<double> ICFull(fulldim);

    //take the ICs as stated in the paper 
    for(unsigned i = 0; i < spacepts; ++i){
      ICFull[i] = i*hx;
      ICFull[spacepts+i] = ICFull[i]/(1+a*ICFull[i]); //lies on the manifold 
    }
      
    

    double timefull = 0., time1st = 0., time2nd = 0.;
    
    CPUTime<double> cpu;   //measure elapsed CPU time
  
    unsigned maxit = 20000,
      newtmaxit = 7;
    
    
    //==================== ORIGIANL PROBLEM ===================================
    cout << endl << " === ORIGINAL SYSTEM ===" << endl;
    
    double ntol = 1.e-04,
      rtol = 1.e-03;
    
    double hEstim = 0.0025;
    //original spatio-temporal 1D system 
    ODE<double,XDSFull> ivp_full(t0,tend,ICFull);
    
    cpu.start();
    
    MyVec<double> Sol_full = ivp_full.implicit_euler(maxit,newtmaxit,ntol,rtol,"OrigRP.dat",hEstim);
    cpu.stop();
    
    timefull = cpu.get_time();
    
    

    
    //=========================================================================
    //=============================== REDUCTIONS ==============================
    //=========================================================================
   
    MyVec<double> ICRed(reddim);
    
    for(int i = 0; i < reddim; ++i)
      ICRed[i] = ICFull[i];
  
    double tol1 =  1.e-04, //0.5     //newton iter
      tol2 = 1.e-03;//1.e-01;  //0.3      //local error
    double hest = 0.0075; 
    
    
 
#if REDUCTION_TO_BE_USED == 1
    //========================= 1ST APPROX ====================================
    
    cout << "Applying 1st approximation to REN-POPE..." <<endl;
    
    ODE<double,XDS1stApprox>  ivp_1st(t0,tend,ICRed); 
    
    cpu.reset();  //RESET time
    cpu.start();
   
    MyVec<double> Sol_1st = ivp_1st.implicit_euler(20000,7,tol1,tol2,"1stRP.dat",hest);
    // ivp_1st.adaptive_rkf(maxit, hest, 1.e-04, "rkfspace1st.dat");
    cpu.stop();
    
    time1st = cpu.get_time();
#endif


    //SECOND APPROX · SECOND APPROX · SECOND APPROX · SECOND APPROX ·
#if REDUCTION_TO_BE_USED == 2
    cout << endl << " Applying 2nd approximation to REN-POPE..." << endl;
    ODE<double,XDS2ndApproxAkaCPA>  ivp_2nd(t0,tend,ICRed); 
      
    cpu.reset();  //RESET time
    cpu.start();
     
    MyVec<double> Sol_2nd = ivp_2nd.implicit_euler(20000,7,tol1,tol2,"2ndRP.dat",hest);
    
    cpu.stop();
    
    time2nd = cpu.get_time();
#endif





    //PLOT STUFF HERE
    //take last slide, i.e. t = tend with x \in [0,1]
    std::ofstream myo("RPFULL.dat",std::ios_base::out);
#if REDUCTION_TO_BE_USED == 1
    std::ofstream first("RP1ST.dat",std::ios_base::out);
#endif
#if REDUCTION_TO_BE_USED == 2
    std::ofstream second("RP2ND.dat",std::ios_base::out);
#endif
    string heading = "# space x           H2             H2O";

    myo << heading <<endl;
#if REDUCTION_TO_BE_USED == 1
    first << heading << endl;
#endif
#if REDUCTION_TO_BE_USED == 2
    second << heading << endl;
#endif
    
    double x = 0;
    for(unsigned i = 0; i < spacepts; ++i){
      x = i*hx;
      
      myo << " "<< x << "    "<<Sol_full[i] << "    "<<Sol_full[spacepts+i] <<endl;
#if REDUCTION_TO_BE_USED == 1
       first << " "<< x << "    "<<Sol_1st[i]  <<endl; 
#endif
#if REDUCTION_TO_BE_USED == 2
       second << " "<< x << "    "<<Sol_2nd[i]  <<endl; 
#endif
    
    }
  
  
    cout << "System type    ||     elapsed CPU time [secs]"<< endl
	 << "-------------------------------------------"<< endl
	 << " Original      ||        "<<  timefull << endl; 
#if REDUCTION_TO_BE_USED == 1
      cout << " 1st approx.   ||        ";
    time_in_human_readable_form(time1st);
#endif
#if REDUCTION_TO_BE_USED == 2
    cout << " 2nd approx.   ||        ";
    time_in_human_readable_form(time2nd);
#endif
  
#if REDUCTION_TO_BE_USED == 1
      cout<< endl<< "l_2 error (1st) = "<< l_2_error(Sol_full,Sol_1st,index) << endl;
      cout <<"slowdown (1st) = "<< slow_down_rate(time1st,timefull) <<endl;
#endif

#if REDUCTION_TO_BE_USED == 2
      cout<< endl<< "l_2 error (2nd) = "<< l_2_error(Sol_full,Sol_2nd,index) << endl;
      cout <<"slowdown (2nd) = "<< slow_down_rate(time2nd,timefull) <<endl;
      //cout << "No l2 error measurement available at the moment..."<<endl;
#endif

    //============================= 2 EPS ===================================
    bool createEpsFile = false;
    string wc = "case2";
    string myepsfilename = "RP"  + wc + ".eps";
    //=======================================================================

    string toeps;
      
   
    if(createEpsFile) {
      toeps = "set terminal postscript eps color enhanced; set output";
      toeps +=  ("'" + myepsfilename + "';");
    } 
    
    //plot stuff
    GNUplot wildebeest;
  
    wildebeest("set key box left; set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"x\"; set ylabel \"z(x,t=1)\"; set title \"Extended Davis Skodie " + my_function_collection::num2str(spacepts) + " spatial points\";plot \"RPFULL.dat\" using 1:2 with lines lc rgb \"red\" lw 2 title \"z_1 full\", \"RPFULL.dat\" using 1:3 with lines lc rgb \"green\" lw 2 title \"z_2 full\", \"RP1ST.dat\" using 1:2 with lines lc rgb \"violet\" lw 2 title \"z_1 first\", \"RP2ND.dat\" using 1:2 with lines lc rgb \"yellow\" lw 2 title \"z_1 second\"");


    myo.close();
#if REDUCTION_TO_BE_USED == 1
    first.close();
#endif
#if REDUCTION_TO_BE_USED == 2
    second.close();
#endif
    */






  }

  else if (executable == 2){
    
  }

  else if (executable == 3){
    typedef double DType;
    typedef ThermoData4Mechanism<DType,6> DataType;
    const int reddim = DataType::rednspec;
    const int fulldim = DataType::nspec;
    cout<< "fulldim = "<< fulldim << "   reddim = "<<reddim<<endl;
   
    //two initial conditions satisfying the mass balance relations
    MyVec<DType> iv(6);
    iv <<= 0.4549999999992859, 0.7780585738931507, 0.2366143850825262, 0.3628298037265891, 0.1479999999999196, 0.01594142610843904;
 
    // //! not working well
    // random_container(iv);
    // iv[2] = -0.5*((2-1) - 2*iv[0] - iv[1] + iv[3] - iv[4]); //mass balance
 
    //iv <<= 3.0004541410314495e-01, 1.8944296646551151e-01, 1.4354979042038563e-01, 1.0238618155887168e-01, 5.9995196772784154e-01,1.0562269872515486e-02;  //near equilibrium

    

    MyVec<DType> IV(iv);

    cout << "IV satisfies linear mass balance?: "<< CheckMassBalance<MyVec<DType>,6>::check_mass(IV) << endl;
    MyVec<double> C(DataType::mass_balance_matrix(),DataType::mass_balance_matrix()+12),
      b(DataType::mass_sum(),DataType::mass_sum()+2);
    MyVec<double> r(2);
    r = matrix_vector_product(C,IV)-b;
    cout << "C*u-b = "<< r << endl;

    //! compute random composition and try again
    //MyVec<DType> IV(fulldim);  //o.k. l2-err = 0.0386 and speedup of 3.428
    //random_container(IV);

    cout << "Test automatically REDUCED source term:"<<endl;
    MyVec<int> index(DataType::rednspec);
    index[0] = 0; index[1] = 4;

    MyVec<DType> RV(reddim);

    UnrollLoop<0,reddim>::assign_greater_2_smaller_rac(RV,IV,index);

    cout << "RV = "<< RV<<endl;

    cout << setprecision(9) <<endl;

    CPUTime<DType> cpu;
    DType tm,
      tmred;

    // //just some mumbo jumbo
    // cpu.start();
    // ReducedH2MechIn6Species<DType> redsource(reddim);    //same stuff
    // cout << "evaluation at point = "<< redsource(1.5,RV) << endl;
    // tm = cpu.stop();
  
    // cout << "cpu-time [AUTOMAT]: " << tm<<endl;
    
    // cpu.reset();
    
    // cout<<endl<<"cross-check with hand-coded function:"<<endl;
    
    // cpu.start();
    // H2MechPURESOURCE1d1STAPPROX<DType> crosscheckfunction(reddim);
    // cout<<"evaluation point = "<< crosscheckfunction(RV) <<endl;
    // tm = cpu.stop();
    // cout << "cpu-time [BY HAND]: " << tm<<endl<<endl;
    



    //=================== ODE system =========================================
    //H2MechIn6Species  //automatically built source
    //H2Combustion6Spex //original hand coded full stuff
    
    typedef ODE<DType,H2MechIn6Species> FullDiffEqType; //automatic
    //typedef ODE<DType, H2Combustion6Spex> FullDiffEqType; //by hand
    
    string comm,
      fnamefull,
      fname1st;

    MyVec<DType> Sol;
    cpu.reset();  //reset time 

    cpu.start();
    
    //time horizon
    DType t0 = ODE_t0,
      tend = ODE_tend;  //23.66;  //long time inteval
    

    cout<<"#################### ORIGINAL (FULL) SYSTEM ####################### "<<endl;
    FullDiffEqType ivpFull(t0,tend,IV);
    
    fnamefull = "h2c6_AUTO_source.dat";


#if SOLVER_TO_BE_USED == 'i'
    comm = "Implicit Euler";
    
    size_t maxsteps = ODE_IE_maxsteps,
      maxNewtoniters = ODE_IE_maxNewtoniters;
    DType ntol = ODE_IE_ntol,
      rtol = ODE_IE_rtol;

    DType hestim = ODE_IE_hEstim,
      Cscal = ODE_IE_Cscal,
      hmin = ODE_IE_hmin,
      hmax = ODE_IE_hmax;
    
    Sol = ivpFull.implicit_euler(maxsteps,maxNewtoniters,ntol,rtol,fnamefull,hestim,Cscal,hmin,hmax);
#endif

#if SOLVER_TO_BE_USED == 'j'
     //======== explicit method for stiff odes ============================
    comm = "Explicit method 4 stiff ODEs";
    
    DType errtol = ODE_E4S_errtol,
      tolFixPt = ODE_E4S_tolFixPt,
      hEst = ODE_E4S_hEst,
      St = ODE_E4S_St,
      cs = ODE_E4S_cs;
    
    size_t maxFixIter = ODE_E4S_maxFixIter;

    DType kmin = ODE_E4S_kmin,
      kmax = ODE_E4S_kmax;
      
    Sol = ivpFull.explicit_4_stiff(errtol,tolFixPt,hEst,St,cs,maxFixIter,fnamefull,kmin,kmax);
#endif 
    
   
    // ============== TODO: name for output files ============================= 
    bool createEpsFile = false;
    
    string ofile1 = "fullSource.eps",
      ofile2 = "redSource.eps",
      ofile3 = "simt.eps";
    //=========================================================================

    string toeps,
      toeps1,
      toeps2,
      toeps3;
      
    
    if(createEpsFile) {
      toeps = "set terminal postscript eps color enhanced; set output";
      toeps1 = toeps + ("'" + ofile1 + "';");
      toeps2 = toeps + ("'" + ofile2 + "';");
      toeps3 = toeps + ("'" + ofile3 + "';");
    }
    
    

    GNUplot plotMePal;
    plotMePal(toeps1 + "set title \"" + comm + " -- Order of method: " + num2str(ivpFull.order_of_method()) + "\" ; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \" time t \"; set ylabel \" U(t=T)\"; plot for[i=2:" + num2str(Sol.size()+1) + "] \"" + fnamefull + "\" using 1:i with lines lw 1 notitle");
     tm = cpu.stop();

     cout << "Elapsed time (with ODE construction and plotting: ";
     time_in_human_readable_form(tm);


     
     cout << endl<<"##########################################"<<endl;
     cout <<       "########## REDUCED ODE SYSTEM ############"<<endl;
     cout <<       "##########################################"<<endl;
     
#if RED_METH_TO_APPLY == 'M' || RED_METH_TO_APPLY == 'F'
     /////////// ATTEMPT BY ME ////////////////
     typedef ODE<DType,Red6> ReducedDiffEqType;
     //////////////////////////////////////////
#endif
#if RED_METH_TO_APPLY == 'L'
      typedef ODE<DType,ReducedH2MechIn6Species> ReducedDiffEqType; //autom.
     // typedef ODE<DType,H2MechPURESOURCE1d1STAPPROX> ReducedDiffEqType; //hand
#endif

     MyVec<DType> SolRed;

     cpu.reset();
     cpu.start();
     ReducedDiffEqType ivpRed(t0,tend,RV);


     
     
     fname1st = "Red_h2c6_AUTO_source_1ST.dat";
     
     //! explicit Euler
     SolRed = RV;
     Red6<double> RedFun(2);

     double time(t0),
       kn(1.e-04);
     ofstream of(fname1st.c_str());
     while(time < tend){
       of << time << "  " << SolRed[0] << "  " << SolRed[1] << endl;
       SolRed += kn*RedFun(SolRed);
       time += kn;
     }
     of.close();

// #if SOLVER_TO_BE_USED == 'i'
//      comm = "REDUCTION. use Implicit Euler";

//      maxsteps = ODE_IE_maxsteps;
//      maxNewtoniters = ODE_IE_maxNewtoniters;
     
//      ntol = ODE_IE_REDUC_ntol;
//      rtol = ODE_IE_REDUC_rtol;
//      hestim = ODE_IE_REDUC_hEstim;
     
//      Cscal = ODE_IE_Cscal;
//      hmin = ODE_IE_hmin;
//      hmax = ODE_IE_hmax;
    
  
//      SolRed = ivpRed.implicit_euler(maxsteps,maxNewtoniters,ntol,rtol,fname1st,hestim,Cscal,hmin,hmax);
// #endif


// #if SOLVER_TO_BE_USED == 'j'   
//      comm = "REDUCTION. use explicit 4 stiff";
     
//      //arguments for reduction
//      errtol = ODE_E4S_REDUC_errtol;
//      tolFixPt = ODE_E4S_REDUC_tolFixPt;
//      hEst = ODE_E4S_REDUC_hEst;
//      St = ODE_E4S_St;
//      cs = ODE_E4S_REDUC_cs;
    
//      maxFixIter = ODE_E4S_maxFixIter;

//      kmin = ODE_E4S_kmin;
//      kmax = ODE_E4S_kmax;

//      SolRed = ivpRed.explicit_4_stiff(errtol,tolFixPt,hEst,St,cs,maxFixIter,fname1st,kmin,kmax);
// #endif
     
     cout<<endl<<"Compare full and reduced solutions with each other:"<<endl;
     cout << "U*(t =  "<< tend<< ") = "<< Sol << endl; 
     cout << "U*_{red}(t =  "<< tend<< ") = "<< SolRed << endl; 

     cout << endl << "Plot compound solution... "<< endl;
    
     

     string outtex,
       texfile1 = PD.get_datum<string>("outfilename");
     if(totex){
       outtex = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile1+".tex'; set format xy \"$%g$\";";  
     }

     string stit;  //"set title \"" + comm + " -- Order of method: " + num2str(ivpRed.order_of_method()) + "\" ; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics;";
     GNUplot redplot;
     redplot(outtex+ "set xlabel \"$t$\"; plot\"" + fname1st + "\" using 1:2 with lines lw 2 title \"$[\\\\ce{H2}]$ red\", \""+fname1st+"\" using 1:3 with lines lw 2 title \"$[\\\\ce{H2O}]$ red\", \""+fnamefull+ "\" using 1:2 with lines lw 2 title \"$[\\\\ce{H2}]$ full\", \""+fnamefull+"\" using 1:6 with lines lw 2 title \"$[\\\\ce{H2O}]$ full\"");

     tmred = cpu.stop();
     //cout << "1ST: Elapsed time (with ODE construction and plotting): ";
     //time_in_human_readable_form(tmred);
     

     //........................................................................
     
     cout << endl<< "------------------------------------------"<<endl;
     cout <<        "|  full time:     ";
     time_in_human_readable_form(tm);
     cout<<endl;
     cout <<        "|  1st time:      ";
     time_in_human_readable_form(tmred);
     cout<<endl;
     cout<<         "------------------------------------------"<<endl;
     cout<< "l_2 error = "<< l_2_error(Sol,SolRed,index) << endl;
     cout << slow_down_rate(tmred,tm)<<endl;
   
#if RED_METH_TO_APPLY == 'L'
     cout << endl << "MoRe: LEBIEDZ's local method"<<endl;
#endif
#if RED_METH_TO_APPLY == 'M'
       cout << endl << "MoRe: MARC's method" <<endl;
#endif


     // cout <<endl<< "plot SIM(t)"<<endl;
     // GNUplot simt;
     // simt(toeps3 + "set view 60,75; set xlabel \"time t\"; set ylabel \"zM_1\"; set zlabel \"zM_2\"; set title \"SIM(t)\"; splot \""+ fname1st +"\"");
 

  }







  /////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// H2 GRI ////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  else if (executable == 4){
    cout << "----------------------------------------------------------"<<endl;
    cout << "------------- H2 excerpt from GRI Mech 3.0 ---------------"<<endl;
    cout << "----------------------------------------------------------"<<endl;
    
    ParameterData Para;
    Para.read_from_file("dat/sim_start_H2_GRI.dat");

    typedef double DType;
    typedef ThermoData4Mechanism<DType,9> DataType;
    for(int k = 0; k < DataType::rednspec; ++k){
      cout << DataType::species_names()[DataType::rpv_index()[k]] << ":   "<< Para.get_datum<DType>(DataType::species_names()[DataType::rpv_index()[k]]) << endl;
    }

    MyVec<DType> Init(Para.size());

    for(size_t i = 0; i < Para.size()-1; ++i){
      Init[i] = Para.get_datum<DType>(DataType::species_names()[i]); 
    }
    //temperature
    Init[DataType::nspec] = Para.get_datum<DType>("T"); 

    cout << "Init = "<< Init << endl;
    

    // //············· general stuff ··················
    // DType t0 = 0.,
    //   tend = 3.0e-4;

    // CPUTime<DType> cpu;

    // string ieFile,
    //   comm;

    // DType rtm;
    // //··············································


    //==================== REDUCTION =======================================
    // MyVec<DType> RV(DataType::rednspec);
    
    // for(size_t i = 0; i < DataType::rednspec; ++i){
    //   RV[i] = Init[DataType::rpv_index()[i]];
    // }

    // cout << "RV = "<< RV <<endl;

    // typedef ODE<DType,ReducedH2Gri> ReducedDiffEqType; 
    // ReducedDiffEqType ivpRed(t0,tend,RV);

    //  MyVec<DType> SolRed;

    //  cpu.reset();
    //  cpu.start();
     


     
     
    //  ieFile = "Red_h2_GRI_AUTO_source_1ST.dat";
     

    //  comm = "REDUCTION. use Implicit Euler";

    //  size_t maxsteps = 50000,
    //    maxNewtoniters = 7;
    //  DType ntol = 1.e-04,
    //    rtol = 1.e-03,
    //    hestim = 0.001,
    //    Cscal = 1.,
    //    hmin = 1.15e-09,
    //    hmax = 1.35e-02;
    
    //  cpu.reset();
    //  cpu.start();

  
    //  SolRed = ivpRed.implicit_euler(maxsteps,maxNewtoniters,ntol,rtol,ieFile,hestim,Cscal,hmin,hmax);

    //  GNUplot rplotter;
    //  rplotter("set title \"" + comm + " -- Order of method: " + num2str(ivpRed.order_of_method()) + "\" ; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \" time t \"; set ylabel \" U(t=T)\"; plot for[i=2:" + num2str(SolRed.size()+1) + "] \"" + ieFile + "\" using 1:i with lines lw 1 notitle");

    //  rtm = cpu.stop();

  }


  // TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST · TEST ·
  else if (executable == 5){
    

    //works
   // MechanismBounds<double,6> mechbds;
    // mechbds.initialize(6);

    // cout << "low = "<< mechbds.low() << endl << "up = "<< mechbds.up() <<endl;
  }
  else if (executable == 6){
    cout << "solve U'(t) = S(U(t)) for t --> inf:"<<endl;
    typedef ExprTmpl::MyVec<double> VType;
    typedef VType::value_type value_type;
    
    const int dim = 6; 
    StartingPoint<VType> start(dim);
    cout << "Startpoint = "<< start.get_point()<<endl;

     typedef ODE<value_type,H2MechIn6Species> DiffEqType; //automatic
    //typedef ODE<value_type, H2Combustion6Spex> FullDiffEqType; //by hand
    
    MyVec<value_type> Sol;
    
    CPUTime<value_type> cpu;
    cpu.start();  //reset time 
    
    //time horizon
    value_type t0 = 0.,
      tend = 2.25;  //0.2;  //23.66;  //long time inteval
    
    DiffEqType ivp(t0,tend,start.get_point());
    
    unsigned maxsteps = 50000,
      maxNewtoniters = 7;
    value_type ntol = 1.e-04,
      rtol = 1.e-03;

    value_type hestim = 0.01,
      Cscal = 1.,
      hmin = 1.15e-09,
      hmax = 2.;
    
    string ieFile = "SourceH2C6.dat";

    Sol = ivp.implicit_euler(maxsteps,maxNewtoniters,ntol,rtol,ieFile,hestim,Cscal,hmin,hmax);
    
    cpu.stop();
    
    cout << "U* = " << Sol << endl;
    
    GNUplot plotMePal;
    plotMePal("set title \" H2C6 chemical source -- Order of method: " + num2str(ivp.order_of_method()) + "\" ; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \" time t \"; set ylabel \" U(t) in mass fractions\"; plot for[i=2:7] \"" + ieFile + "\" using 1:i with lines lw 2 notitle");

  }
  else if (executable == 7){
    cout << "----------- TEST FACILITY ---------------------"<<endl;
     cout << "Simple nonlinear test mechanism" << endl;
    
     cout<<          "=============================================="<<endl;
    cout << endl << "================ FULL MODEL ===============" <<endl;
    cout<<          "=============================================="<<endl;

    double t0 = 0.,
      tend = 300;
    double hEstim = 0.75; //0.015;
     // double Cscal = 1.111,
     //  hmin = 0.5,
     //  hmax = 1.;
    
    unsigned maxNewt = 7;

    MyVec<double> x0(3);
    x0[0] = 0.35;
    x0[1] = 0.25;    //x[0] + x[1] + x[2] = 1.0
    x0[2] = 0.4;

    ToyMechanism<double> fullfun(3);
   
    CPUTime<double> cpu;
    cpu.start();
    JacS<double,ToyMechanism> Jac;
    Jac.set(3,3);

    ofstream ie("ToyMechanism.dat");
    MyVec<double> U(x0), Uprev(x0), G(3),Gp(3*3);
    double time(t0);
    unsigned c(1);

    //USE IMPLICIT EULER
    while(time < tend){
      cout << c << ".)   time: "<< time << "   dt = "<< hEstim << endl;
      ie << setprecision(12) << time << "  "<< U[0] << "   "<< U[1] << "   "<<U[2] <<endl;
      U = Uprev;
      unsigned count(1);
      for(unsigned i = 1; i <= maxNewt; ++i){
	count = i;
	G = -(U - Uprev - hEstim*fullfun(U));
	Gp *= -hEstim*Jac.jacobian(U);
	update_diagonal<AddBasicElements>(Gp,3,1.);
	good_square_solve(Gp,3,G,1);
	U += G;
	if(Norm<'2',double>::norm(G) <= 1.e-08)
      	  break;
      } //end Newton
      
      if(count == maxNewt){
	ADONIS_INFO(Information, "Newton iteration did not converge in "<< maxNewt << " iterations.");
      }
      Uprev = U;

      time += hEstim;
      ++c; //iteration counter
    }
    ie.close();
   
    double futi = cpu.stop();

    cout<<          "=============================================="<<endl;
    cout << endl << "================ REDUCED MODEL ===============" <<endl;
    cout<<          "=============================================="<<endl;
   
    cpu.reset(); //reset timer

    MyVec<double> rpv(1);
    rpv <<= x0[2];
    ReducedToyMechanism<double> redfun(1);
  
    time = t0;  //reuse it
  
    ofstream of("RedToy.dat");
    of << "# time     y2"<<endl;
    cpu.start();
    unsigned cred(1);
    //USE EXPLICIT EULER
    while(time < tend){
      of << " " <<time << "     " << rpv[0] << endl;
      cout << cred << ".)     t = "<<time << "  dt = "<< hEstim << endl;
      rpv += hEstim*redfun(rpv);
      time += hEstim;
      ++cred; //counter
    }
    of.close();
    double reti = cpu.stop();

    cout << "Full time (after "<<c<<" iters): "; time_in_human_readable_form(futi);
    cout << "Red. time (after "<<cred<<" iters): ";   time_in_human_readable_form(reti);
    //slow_down_rate(reti,futi);
    
    cout << "Solution (FULL) = "<< U << endl;
    cout << "Solution (RED)  = "<<rpv[0] <<endl;

    cout << "l_2-error(tend = "<<tend <<") = "<<setprecision(12)<< Abs(U[2]-rpv[0]) <<endl;

    is_faster(reti,futi); //check if first argument is smaller,equal,greather than second


    
    string outtex,
      texfile = PD.get_datum<string>("outfilename");
    if(totex){  //change it in dat/solversettings.dat
       outtex = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile+".tex'; set format xy \"$%g$\";";  
     }
    string y3lab;
    bool setylabel = PD.get_datum<bool>("setylabel");
    if(setylabel)
      y3lab = "$y_3$"; 


    GNUplot show;
    show(outtex+"set key box left; set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"$t$\"; set ylabel \"$"+y3lab+"$\"; plot 'ToyMechanism.dat' using 1:4 with lines lw 2 lc rgbcolor \"#DC143C\" title 'full', 'RedToy.dat' using 1:2 with lines lw 2 lc rgbcolor \"#32CD32\" title 'reduced'");
    
    //show("plot 'ToyMechanism.dat' using 1:2 with lines lw 2 title 'y1', 'ToyMechanism.dat' using 1:3 with lines lw 2 title 'y2', 'ToyMechanism.dat' using 1:4 with lines lw 2 title 'y3'");
    


    
    

    // typedef double Type;
    // const int nelem = 2;
    // const int nred = 2;
    // const int nspec = 6;

    // const int mech = 6;

    // typedef ThermoData4Mechanism<double,mech> DataType; //H2C6 Mech

    // MyVec<double> Ct(nelem*6), b(nelem);
    // Ct <<= 2,1,0,0,2,1,
    //   0,0,2,1,1,1;

    // b <<= 2,1;

    // // Dune::FieldMatrix<double,nred,nspec> Bt;
    // // BT = 0;
    // // BT[0][0] = 1.;
    // // BT[1][4] = 1.;

    
      

    // typedef MyVec<Type> VType;
    // typedef VerySimpleSpeciesReconstruction<double,6,2,mech> SpecConstType;
    // SpecConstType VSR;
    
    // VType zM0(6);
    // random_container(zM0);
    // cout << "zM0 = "<<zM0 << endl;

    // VSR.initialize(zM0,2,Ct,b,DataType::rpv_index());
    
    // MyVec<double> r(2), full(6);
    // r <<= 0.2, 0.8;
    // full <<= 0.1, 0.3, 0.01, 0.6, 0.9, 0.7;

    // VType z = VSR.evaluate(r,full);
    // cout << "x_tilde = "<< z << endl; 

    // // //!works as expected
    // // VType diag(6*6);
    // // diagonal_matrix(diag,full);  //create diagonal matrix
    // // cout << "diagonal matrix = "<<endl;
    // // print_in_matrix_style(diag,6);
    
  }
  else if(executable == 8 ){
    cout << "ZELDOVICH MECHANISM" <<endl;
    cout<<          "=============================================="<<endl;
    cout << endl << "================ FULL MODEL ===============" <<endl;
    cout<<          "=============================================="<<endl;

    ParameterData ZeldDat;
    ZeldDat.read_from_file("dat/zeldovich.dat");
    double t0 = ZeldDat.get_datum<double>("t0"),
      tend = ZeldDat.get_datum<double>("tend");
    double hEstim = ZeldDat.get_datum<double>("hEstim"); 
    
    adonis_assert(hEstim <= tend);

    unsigned maxNewt = 5;

    const int dim = 8;

    MyVec<double> x0(dim);
    x0 <<= 
      ZeldDat.get_datum<double>("O"),  
      ZeldDat.get_datum<double>("N2"),
      ZeldDat.get_datum<double>("NO"),
      ZeldDat.get_datum<double>("N"),
      ZeldDat.get_datum<double>("O2"),
      ZeldDat.get_datum<double>("OH"),
      ZeldDat.get_datum<double>("H"),
      ZeldDat.get_datum<double>("Temperature");
      
    adonis_assert((1. - UnrollLoop<0,7>::sum<1>(x0)) == 0.);

    cout << "x0 = "<< x0 << endl;
    
    cout << "used prefactor for k_i: "<<  ZeldDat.get_datum<double>("prescale") << endl << "p0: "<<  ZeldDat.get_datum<double>("p0") << endl;

    Zeldovich<double> fullfun(dim);
    cout << "Zeldovich evaluated at x0: "<< fullfun(x0) <<endl;

    CPUTime<double> cpu;
    cpu.start();
    JacS<double,Zeldovich> Jac;
    Jac.set(dim,dim);

    ofstream ie("Zeldovich.dat");
    MyVec<double> U(x0), Uprev(x0), G(dim),Gp(dim*dim);
    double time(t0);
    unsigned c(1);

    //USE IMPLICIT EULER
    while(time < tend){
      // U += hEstim*fullfun(U);
      
      cout << c << ".)   time: "<< time << "   dt = "<< hEstim << endl;
      ie << setprecision(12) << time << "  ";
      for(int i = 0; i < dim; ++i)
      	ie << U[i] << "  ";
       ie << endl;
      
      U = Uprev;
      unsigned count(1);
      for(unsigned i = 1; i <= maxNewt; ++i){
      	count = i;
      	G = -(U - Uprev - hEstim*fullfun(U));
      	Gp *= -hEstim*Jac.jacobian(U);
      	update_diagonal<AddBasicElements>(Gp,dim,1.);
      	good_square_solve(Gp,dim,G,1);
      	U += G;
      	if(Norm<'2',double>::norm(G) <= 1.e-08)
      	  break;
      } //end Newton
      
      if(count == maxNewt){
      	ADONIS_INFO(Information, "Newton iteration did not converge in "<< maxNewt << " iterations.");
      }
      Uprev = U;

      time += hEstim;
      ++c; //iteration counter
    }
    ie.close();
   
    double futi = cpu.stop();

    cout << "Solution "<< U << endl;

    bool sfullgraphs = ZeldDat.get_datum<bool>("showfullgraphs");

    if(sfullgraphs){
      GNUplot show;
      show("set key box right; set xlabel '$t$';set ylabel '$Y_k, k=1,\\\\ldots,7$'; plot 'Zeldovich.dat' using 1:2 with lines lw 2 title '$y_{\\\\ce{O}}$', 'Zeldovich.dat' using 1:3 with lines lw 2 title '$y_{\\\\ce{N2}}$', 'Zeldovich.dat' using 1:4 with lines lw 2 title 'y_{\\\\ce{NO}}', 'Zeldovich.dat' using 1:5 with lines lw 2 title 'y_{\\\\ce{N}}', 'Zeldovich.dat' using 1:6 with lines lw 2 title 'y_{\\\\ce{O2}}', 'Zeldovich.dat' using 1:7 with lines lw 2 title 'y_{\\\\ce{OH}}', 'Zeldovich.dat' using 1:8 with lines lw 2 title 'y_{\\\\ce{H}}'");
    
      GNUplot showtemp;
      showtemp("set ylabel '$T$';set xlabel '$t$';plot 'Zeldovich.dat' using 1:9 with lines lw 2 notitle");

    }

     cout<<          "=============================================="<<endl;
    cout << endl << "================ REDUCED MODEL ===============" <<endl;
    cout<<          "=============================================="<<endl;
  

    cpu.reset(); //reset timer


    Zeldovich<double> Zeld(8);
    cout << "x0 = "<< x0 << endl;
    cout << "full(·) = "<< Zeld(x0) << endl << endl;

#if ZELDOVICH_NUM_RPV == 3
    MyVec<double> rpv(4);
    rpv <<= x0[1], x0[2], x0[4], //species 
      x0[7];  //temperature 
    ReducedZeldovich<double> redfun(4); //4 th is temperature
    cout << "f(·) = "<<redfun(rpv) << endl;
#endif

#if ZELDOVICH_NUM_RPV == 2
    MyVec<double> rpv(3);
    rpv <<= x0[1], x0[4], //species 
      x0[7];  //temperature 
    TwoSpeciesReducedZeldovich<double> redfun(4);
#endif

    time = t0;  //reuse it
  
    ofstream of("redzeldo.dat");
    of << "# time     y_N2    y_NO    y_O2    T"<<endl;
    cpu.start();
    unsigned cred(1);
    //USE EXPLICIT EULER -- equidistant
    while(time < tend){
      of << " " <<time << "     ";
      for(size_t l = 0; l < rpv.size(); ++l)
	of << rpv[l] << "   ";
      of << endl;
     
      cout << cred << ".)     t = "<<time << "  dt = "<< hEstim << endl;
      rpv += hEstim*redfun(rpv);
      time += hEstim;
      ++cred; //counter
    }
    of.close();
    double reti = cpu.stop();

    cout << "Full time (after "<<c<<" iters): "; time_in_human_readable_form(futi);
    cout << "Red. time (after "<<cred<<" iters): ";   time_in_human_readable_form(reti);
    //slow_down_rate(reti,futi);
    
    cout << "Solution (FULL) = "<< U << endl;
    cout << "Solution (RED)  = "<<rpv <<endl;

    
    cout << "l_2 (species) at tend = "<<tend<< setprecision(12) <<": "<< 
#if ZELDOVICH_NUM_RPV == 3
sqrt(ntimes<2>((Abs(U[1] - rpv[0]))) + ntimes<2>((Abs(U[2] - rpv[1]))) + ntimes<2>((Abs(U[4] - rpv[2]))))<< endl;
#endif
#if ZELDOVICH_NUM_RPV == 2
    sqrt(ntimes<2>((Abs(U[1] - rpv[0])))  + ntimes<2>((Abs(U[4] - rpv[1]))))<< endl;
#endif

    is_faster(reti,futi); //check if first argument is smaller,equal,greather than second
    
    
    string outtex,outtex2,
      texfile = PD.get_datum<string>("outfilename"),
      texfile2 = PD.get_datum<string>("outfilename2");
    if(totex){  //change it in dat/solversettings.dat
      outtex = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile+".tex'; set format xy \"$%g$\";";
       outtex2 = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile2+".tex'; set format xy \"$%g$\";";
     }
    
    string xscale = "set xtics 0,.75e-06,2.e-06;";
    GNUplot bothspec;
#if ZELDOVICH_NUM_RPV == 3
    bothspec(outtex+xscale+"set key bottom right; set mxtics 4; set mytics 4;  set grid xtics ytics mxtics mytics; set xlabel \"$t$\"; plot 'Zeldovich.dat' using 1:3 with lines lw 2 lc rgbcolor \"#DC143C\" title \"$Y_{\\\\ce{N2}}$ (full)\", 'Zeldovich.dat' using 1:4 with lines lw 2 lc rgbcolor \"#32CD32\" title \"$Y_{\\\\ce{NO}}$ (full)\",'Zeldovich.dat' using 1:6 with lines lw 2 lc rgbcolor \"#5F9EA0\" title \"$Y_{\\\\ce{O2}}$ (full)\",'redzeldo.dat' using 1:2 with lines lw 2 title \"$Y_{\\\\ce{N2}}$ (red)\",'redzeldo.dat' using 1:3 with lines lw 2 title \"$Y_{\\\\ce{NO}}$ (red)\", 'redzeldo.dat' using 1:4 with lines lw 2 lc rgb 'blue' title \"$Y_{\\\\ce{O2}}$ (red)\"");
#endif
#if ZELDOVICH_NUM_RPV == 2
 bothspec(outtex+xscale+"set key center right; set mxtics 4; set mytics 4;  set grid xtics ytics mxtics mytics; set xlabel \"$t$\"; plot 'Zeldovich.dat' using 1:3 with lines lw 2 lc rgbcolor \"#DC143C\" title \"$Y_{\\\\ce{N2}}$ (full)\",'Zeldovich.dat' using 1:6 with lines lw 2 lc rgbcolor \"#5F9EA0\" title \"$Y_{\\\\ce{O2}}$ (full)\",'redzeldo.dat' using 1:2 with lines lw 2 title \"$Y_{\\\\ce{N2}}$ (red)\",'redzeldo.dat' using 1:3 with lines lw 2 title \"$Y_{\\\\ce{O2}}$ (red)\"");
#endif


    //plot temperature
    GNUplot bothtmp;
#if ZELDOVICH_NUM_RPV == 3
    // bothtmp(outtex2+"set key bottom right; unset ylabel; unset ytics; set xlabel \"$t$\";set y2tics mirror;set mxtics 4; set my2tics 4; set grid xtics y2tics mxtics my2tics; set xlabel \"$t$\"; plot 'Zeldovich.dat' using 1:9 with lines lw 2 title \"$T$ (full)\", 'redzeldo.dat' using 1:5 with lines lw 2 title \"$T$ (red)\" axes x1y2");
    bothtmp(outtex2+xscale+"set key center right; set xlabel \"$t$\";set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"$t$\"; plot 'Zeldovich.dat' using 1:9 with lines lw 2 lc rgb 'black' title \"$T$ (full)\", 'redzeldo.dat' using 1:5 with lines lw 2 title \"$T$ (red)\" ");
#endif
#if ZELDOVICH_NUM_RPV == 2
 bothtmp(outtex2+xscale+"set key center right; set xlabel \"$t$\";set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"$t$\"; plot 'Zeldovich.dat' using 1:9 with lines lw 2 lc rgb 'black' title \"$T$ (full)\", 'redzeldo.dat' using 1:4 with lines lw 2 title \"$T$ (red)\" ");
#endif


  }
  else if (executable == 9){
    cout << "********************************************************"<<endl;
    cout << "*************** 1D H2C6 REVISITED_ *********************"<<endl;
    cout << "********************************************************"<<endl;

    ParameterData PD;
    PD.read_from_file("dat/h2c6.dat");

    int nx = PD.get_datum<int>("Nx");
    double hx = (PD.get_datum<double>("b")-PD.get_datum<double>("a"))/(nx-1),
      t0 = PD.get_datum<double>("t0"),
      tend = PD.get_datum<double>("tend"),
      hEstim = PD.get_datum<double>("kn"),
      hmin = PD.get_datum<double>("hmin"),
      hmax = PD.get_datum<double>("hmax");


    MyVec<double> ustar(6);
    ustar[0] = 0.27;      // H2
    ustar[1] = 0.05;
    ustar[2] = 0.135;
    ustar[3] = 0.02;
    ustar[4] = 0.7;       // H2O
    ustar[5] = 0.01;
    //! start value
    StartingPoint<MyVec<double> > startpoint(6);
    
    MyVec<double> U0 = startpoint.get_point();

    MyVec<double> IV(6*nx);

    MyVec<double> u0(6), u1(6);
    
    u0 <<= 0.455, 0.779, 0.237, 0.363, 0.148, 0.015;

    u1 <<= 0.2, 0.95, 0.31, 0.03, 0.3, 0.05;

    double x = 0.;
    ofstream ivof("IV0.dat");
    MyVec<MyVec<double> > Tabularization(nx);

    //initial conditions -- constant at each grid point x_i
    for(int i = 0; i < nx ;++i){
      (Tabularization[i]).resize(6);
      x = i*hx;  //CAUTION: previously: x += i*hx
      if(0.002 <= x && x <= 0.008){
      	for(int k = 0; k < 6; ++k)
      	  IV[k*nx+i] = u0[k];
      }
      if(x < 0.002 || x > 0.008){
      	for(int k = 0; k < 6; ++k)
      	  IV[k*nx+i] = u1[k];
      }
      
      for(int k = 0; k < 6; ++k){
	Tabularization[i][k] = IV[k*nx+i];
      }
      //cout << "tabulated val = "<< Tabularization[i] <<endl;
      // // !older version
      // if( (x <= 0.002) || (x > 0.008) ){
      // 	IV[i] = ustar[0];
      // 	IV[nx+i] = ustar[1];
      // 	IV[2*nx+i] = ustar[2];
      // 	IV[3*nx+i] = ustar[3];
      // 	IV[4*nx+i] = ustar[4];
      // 	IV[5*nx+i] = ustar[5];
      // }
      // else{
      // 	IV[i] = U0[0];
      // 	IV[nx+i] = U0[1];
      // 	IV[2*nx+i] = U0[2];
      // 	IV[3*nx+i] = U0[3];
      // 	IV[4*nx+i] = U0[4];
      // 	IV[5*nx+i] = U0[5];
      // }
      ivof << x <<"  ";
      for(int k = 0; k < 6; ++k)
      	ivof << IV[k*nx+i] << "  ";
      ivof << endl;

    }
    ivof.close();
 
    
    //=============== test with valgrind ==============
    // int reddim = 2*nx;
    // MyVec<double> eval(reddim);
    // for(int i = 0; i < nx; ++i){
    //   eval[i] = IV[i];
    //   eval[1*nx+i] = IV[4*nx+i];
    // }
    // ReducedMOLH2C61D<double> fun(2*nx);
    // //fun(eval);

    // MyVec<double> zM(6);
    // fun.reconstruct_full_state(zM,nx-1,eval);
    //=================================================

#if FULL_MODEL == 1   
    GNUplot ivsee;
    ivsee("set title \"H2C6 initial conditions\"; plot 'IV0.dat' using 1:2 with lines lw 2 title \"[H2]\", 'IV0.dat' using 1:4 with lines lw 2 title \"[H2O]\"");
#endif

    string fname = PD.get_datum<string>("filenamefull");
    CPUTime<double> cpu;
    
#if FULL_MODEL == 1
    cout<<          "=============================================="<<endl;
    cout << endl << "================== FULL MODEL ================" <<endl;
    cout<<          "=============================================="<<endl;
    
   
    typedef ODE<double,MOLH2C61D,1,true,false> FullModelType;
    typedef FullModelType::IndexVecType IndexVecType;
    IndexVecType nothing(2);
    nothing <<= 0,4; //plot only species 0 and 4
    FullModelType ivp(t0,tend,IV);

    cpu.start();
    MyVec<double> SolFull = ivp.implicit_euler(10000000,7,1.e-08,1.e-08,fname,hEstim,1.11111,hmin,hmax,nothing,"dat/h2c6.dat",13);
    double tfull = cpu.stop();
#endif
  
#if SWITCH_ON_REDUCTION == 1
    cout<<          "=============================================="<<endl;
    cout << endl << "================ REDUCED MODEL ===============" <<endl;
    cout<<          "=============================================="<<endl;
    
    cpu.reset();
    
    fname = PD.get_datum<string>("filenameredu");
   
    x = 0;
    ofstream ivredof("rediv.dat");
    MyVec<double> RV(2*nx);
    for(int i = 0; i < nx; ++i){
      x = i*hx;
      RV[i] = IV[i];
      RV[1*nx+i] = IV[4*nx+i];
    
      ivredof << x << "  "<< RV[i] << "   " << RV[1*nx+i] << endl;
    }
    ivredof.close();
    GNUplot sri;
    sri("set title \"IV for the reduced model\";plot 'rediv.dat' using 1:2 with lines lw 2 lc rgbcolor  \"#DC143C\" title \"[H2]\", 'rediv.dat' using 1:3 with lines lw 2 lc rgbcolor \"#32CD32\" title \"[H2O]\"");
   

    cout << "RV = "<< RV << endl;

    typedef ODE<double,ReducedMOLH2C61D,1,true,false> REDModelType;
    typedef REDModelType::IndexVecType IndexVecType;
    IndexVecType noth;
    REDModelType redivp(t0,tend,RV);
    
    cpu.start();
    MyVec<double> SolRed = redivp.implicit_euler(10000000,7,1.e-08,1.e-08,fname,hEstim,1.11111,hmin,hmax,noth,"dat/h2c6.dat",13);
    double redfull = cpu.stop();
  
 
#if FULL_MODEL == 1
    cout << "time_FULL = "<< tfull << endl << "time_RED = "<< redfull <<endl;
    is_faster(redfull,tfull);
   
    MyVec<int> rpvix(2);
    rpvix <<= 0,4;
    cout << "rpvix = "<< rpvix << endl << "SolFull.size(): "<<SolFull.size() << "  SolRed.size() = "<< SolRed.size() << endl;

    cout << "l2_error at tend = "<< tend << ": ";
    //cout << l_2_error(SolFull,SolRed,&rpvix[0]) <<endl;
    double l2 = 0.;
    for(int i = 0; i < nx; ++i)
      l2 += ntimes<2>(Abs(SolFull[0*nx+i]-SolRed[0*nx + i])) + ntimes<2>(Abs(SolFull[4*nx+i]-SolRed[1*nx + i]));

    adonis_assert(l2 > 0.);
    l2 = sqrt(l2);
    cout << l2 << endl;
      

#endif
//   

#endif

      }
  
  //default argument 
  else{
    ADONIS_ERROR(MainProgramError, "Example '"<<executable<<"' hasn't been declared for \""<<argv[0]<<"\" yet. \n");

  }

    return 0;
}

/* ----- ----- eof ----- ----- */
