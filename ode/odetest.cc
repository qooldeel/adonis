#include "../want_2_use.h"  //put it at the very beginning of the main file!

#include "ode.hh"
#include "examples/testcases.hh"
#include "tol.hh"
#include "../misc/misctmps.hh"
#include "../common/error.hh"
#include "../common/numericaltraits.hh"
#include "../common/xmath.hh"

//++++++++++++++ begin ren pope ++++++++++++++++
#include "examples/diffusionexample.hh"
#include "examples/firstapproximation.hh"
#include "examples/withcorrection.hh"
#include "examples/secapproxCPA.hh"

#include "examples/extendedDavisSkodje.hh"  //source term only without diffusion
//++++++++++++++ ren pope end ++++++++++++++++++

#include "../common/tempusfugit.hh"

#include<cmath>
#include<map>

#include "sourceterms.hh"
#include "../misc/misctmps.hh"
#include "../dunexternal/dunextensions.hh"
#include "../containers/fdiagonal.hh"
#include "../common/smartassign.hh"
#include "../misc/useful.hh"
#include "rpparameter.hh"
#include "../graphics/printer.hh"
#include "../common/numerictypechecker.hh"
#include "../optimization/linesearch.hh"
#include "../containers/densesymmetricmatrix.hh"
#include "../io/readinparameters.hh"
#include "stiffnesschecker.hh"
#include "../noderivativemethods/finitedifferences.hh"
#include "examples/simpleparabolicexample.hh"

#include "../expressiontemplates/exprvec.hh"

// TEST SUITE //
#include "examples/Test_Set_4_IVP_Solvers/testset.hh"
////////////////


//===== thermodynamics ====

//new stuff
#include "../massactionkinetics/thermochemistry.hh"
#include "../massactionkinetics/stoichiometry.hh"
#include "../massactionkinetics/reactionrates.hh"
#include "../massactionkinetics/buildchemicalrhs.hh"
#include "../massactionkinetics/quantity.hh"
#include "../massactionkinetics/data/thermochemicaldata.hh"
#include "../massactionkinetics/indexinjector.hh"
//=========================

//==== automatically computed sources =============
#include "examples/automaticallygeneratedsourceterms/h2gri.hh"
#include "examples/automaticallygeneratedsourceterms/h2grimarc.hh"
#include "examples/automaticallygeneratedsourceterms/h2c6.hh"
#include "examples/automaticallygeneratedsourceterms/o3.hh"

#include "examples/automaticallygeneratedsourceterms/conaireH2.hh"
#include "examples/automaticallygeneratedsourceterms/chemkinH2default.hh"

#include "examples/automaticallygeneratedsourceterms/2Ddiscr/navierstokes.hh"
//=================================================


//======= transport stuff ==========================
#include "../moleculartransport/prefactors.hh"
#include "../moleculartransport/diffusion.hh"
#include "../moleculartransport/multicomponenttransport.hh"
#include "../moleculartransport/viscosity.hh"
#include "../moleculartransport/blockalgorithms.hh"
//==================================================

//========== fancy TMPs ==========================
#include "../templatemetaprograms/matrixunroller.hh"
//================================================

#include "../containers/staticarray.hh"

//======== SPECIAL TESTING  ===============================================
#include "examples/peculiarfun/sillyfun.hh"


#include "../optimization/augmentedlagrangian.hh"

#include "../optimization/mysbnlpexamples.hh"
#include "../optimization/steihaug.hh"
#include "../optimization/hessianapproximations.hh"
#include "../optimization/trustregion.hh"
#include "../optimization/derivauglag.hh"

#include "../regularization/guaranteepositivity.hh"


#include "../preconditioner/precond.hh"
#include "../statistics/probability.hh"

#include "additional/somestuff.hh"

#include "../sparse/umftypetraits.hh"
#include "../sparse/sparsematrix.hh"

#include "examples/geodynamics/transientheatequation2D.hh"
#include "examples/michaelismenten.hh"
#include "examples/lamgoussismech.hh"

#include "../marcvecmatrix/MAVEC.h"
#include "../modred/indexmanipulation.hh"
#include "explicit.hh"
#include "../modred/massbalanceinvariant.hh"

#include "../dunexternal/tangentspace.hh"

#include "../random/lecuyer32.hh"
#include "../random/lcgminstd.hh"

#include "../derivatives/sparsecppaddriver.hh"

#include "../fdm/nseh2conaire.hh"

#include "../io/readinmatrix.hh"

#include "../fdm/gensettings.hh"

//MOL examples
#include "examples/Schiesser_MOL/burgers.hh"
#include "examples/Schiesser_MOL/cdr_ch13.hh"

//

//=========================================================================



using namespace std;
using namespace Adonis;
using namespace ExprTmpl;
using namespace my_function_collection;

template<class V, class IT>
typename V::value_type phi_func(size_t k, size_t j,const V& mw, const IT& eta){
  return 1./sqrt(8.)*pow((1+mw[k]/mw[j]),-0.5)*pow((1 + sqrt(eta[k]/eta[j])*pow(mw[j]/mw[k],0.25)),2.);
}

inline int symm_offset(int i, int j, int n){
  return (j>=i ? i*n-(i+1)*i/2+j : symm_offset(j,i,n));
}



template<class T>
inline T exact_integral(const T& a, const T& b){
  return (-cos(b) + cos(a) + b*b - a*a + 2.*(ntimes<4>(b) - ntimes<4>(a)));
}





//! just for test cases
template<class T>
class ReducedH2MechIn6Species{
public:
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  typedef T value_type;
  typedef std::size_t SizeType;
  
  enum{H2=0,H=1,O2=2,O=3,H2O=4,OH=5};

  
  ReducedH2MechIn6Species(SizeType n = 2):v_(n){
    k1_ = 2.;     km1_ = 216.;
    k3_ = 1.;     km3_ = 1400.;
    k4_ = 1000.;  km4_ = 10800.;
    k6_ = 100.;   km6_ = 0.7714;
  }

  SizeType domain_dim() const {return 6;}
  SizeType dim() const {return 2;}

  template<class X>
  VType& operator()(const X& x){
    v_[0] = -k1_*x[H2] + km1_*Adonis::ntimes<2>(x[H]) - k4_*x[H2]*x[O] + km4_*x[H]*x[OH] - k6_*x[H2]*x[O] + km6_*x[H2O];
    v_[1] = -k3_*x[H2O] + km3_*x[H]*x[OH] + k6_*x[H2]*x[O] - km6_*x[H2O];
    return v_;
  }

private:
  VType v_;
  T k1_, km1_, k3_, km3_, k4_, km4_, k6_, km6_;
};


template<class T>
class TestFunction{
public:
  typedef T value_type;
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  size_t dim() {return 3;}
  size_t dom_dim() {return 4;}

  TestFunction():v_(3){}
  
  template<class X>
  VType& operator()(const X& x) {
    v_[0] = sin(x[0]*x[2]*x[2]) + 3.5*(x[3]) - x[1]*x[0];
    v_[1] = x[0]*x[1]*x[2]*x[3];
    v_[2] = cos(x[0] + x[1] - 2.45*x[3]) - x[2]*x[3];
    
    return v_;
  }

  private:
  VType v_;

};


template<class T>
class SomeFunctional{
public:
  typedef T value_type;
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  size_t dim() {return 1;}
  size_t dom_dim() {return 4;}

  template<class X>
  T operator()(const X& x){
    return sin(x[0]*x[1] - 1-5*x[2]) + x[3]*(x[3] - x[1]*x[0]*cos(x[0]));
  }
};


template<class T>
class Function2Test{
public:
  typedef T value_type;
  typedef Adonis::ExprTmpl::MyVec<T> VType;
  size_t dim() {return 2;}
  size_t domain_dim() {return 2;}

  Function2Test():v_(2){}
  
  template<class X>
  VType& operator()(const X& x) {
    v_[0] = sin(x[0]/x[1]) + x[0]*x[0];
    v_[1] = ntimes<2>(cos(x[0]*x[1])) - 2 + sqrt(x[1]*x[1]);  //x[0]*x[1] - 1.5*x[1]*x[1];
    
    return v_;
  }

  private:
  VType v_;

};




template<class FUN>
class ForwardFDJacobian{
public:
  typedef typename FUN::value_type value_type;
  typedef matrix<value_type> MatrixType;

  ForwardFDJacobian(FUN& f):fun_(f),jac_((int)f.dim(),(int)f.domain_dim()){}

  template<class X, class V>
  MatrixType& evaluate(const X& x, const V& fx , const value_type& varepsilon){
    V c(fun_.dim());
     for(size_t i = 0; i < fun_.domain_dim(); ++i){
       V p(fun_.domain_dim());
       p[i] = varepsilon;
       c = fun_(x + p);
       c -= fx;
       c /= varepsilon;
       for(size_t j = 0; j < fun_.dim(); ++j){
	 jac_[j][i] = c[j];  //watch out for right matrix access!!!
       }
     }
     return jac_;
  }
  
  const MatrixType& get_jacobian() const {return jac_;} 

private:
  FUN& fun_;
  MatrixType jac_;
};


template<class FUN>
class CentralFDJacobian{
public:
  typedef typename FUN::value_type value_type;
  typedef matrix<value_type> MatrixType;
  typedef Adonis::ExprTmpl::MyVec<value_type> VType;

  CentralFDJacobian(FUN& f):fun_(f),jac_((int)f.dim(),(int)f.domain_dim()){}

  template<class X>
  MatrixType& evaluate(const X& x, const value_type& varepsilon){
    VType c(fun_.dim());
     for(size_t i = 0; i < fun_.domain_dim(); ++i){
       VType p(fun_.domain_dim());
       p[i] = varepsilon;
       c = fun_(x + p);
       c -= fun_(x - p);
       c /= (2.*varepsilon);
       for(size_t j = 0; j < fun_.dim(); ++j){
	 jac_[j][i] = c[j];  //watch out for right matrix access!!!
       }
     }
     return jac_;
  }
  
  const MatrixType& get_jacobian() const {return jac_;} 

private:
  FUN& fun_;
  MatrixType jac_;
};



void change_some_boolean(bool& change){
  change = true;
}

template <class T1, class T2>
T2 frouzakis_profile(const T1& x, const T1& xignite, const T2& tmin, const T2& tmax){
  return tmin*(1.1*cos(PhysicalConstants<T1>::Pi*(x/0.5+1)) +2.1);
}

template<class T1, class T2>
T2 hyp_tang_T_distribution(const T1& x, const T1& xignite, const T2& tmin, const T2& tmax, const T1& a = 1.75, const T1& eps = 0.0000345){
  //return tmin + (tmax-tmin)*tanh((x-xignite));
  // return ( tmin + 0.5*(tmax-tmin)*(1+tanh(2*a*(x-xignite)/(tmax-tmin))) );
  //return ( tmin + 0.5*(tmax-tmin)*(1+tanh(a*(x-xignite))) );
  //return ( tmin + (tmax-tmin)*pow((tanh(a*(x-xignite))),3) );

  // return (1 + (tmax-tmin)*(1 + tanh((x-xignite)/a)) ); //Frouzakis

  //return ( (x <= xignite) ? tmin*(1.1*cos(PhysicalConstants<T1>::Pi*(x/0.5+1)) +2.1) : tmax ); 

  //MARC
  T2 res;
  if(x <= xignite)
    res = frouzakis_profile(x,xignite,tmin,tmax);
  else if ( (xignite < x) && (x <= xignite + eps))
    res = tmax;
  else if( x > xignite +eps)
    res = frouzakis_profile(-x,xignite,tmin,tmax); //-x
  return res;
}



///////////////////////////////////////////////////////////////////////////////
////////////////////////////// MAIN PROGRAM  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv){
 
  if(argc != 2){
    ADONIS_ERROR(MainProgramError, "\""<< argv[0] <<"\" takes exactely 1 argument, pal.");
  }

  int example = atoi(argv[1]);


  if(example == -1){
    
    double x(0.), hy = 1./22;
    for(int j = 0; j < 23; ++j){
      cout << " x = "<<x << "       veloprofile = " << InletVelocity<double,double,double,'c'>::velocity(0.00025, 0.0005, 0.75) << endl;   //first 2 args are dummies
      x += hy;
    }
    
  }
  
  else if(example == 1){   
  
    cout << "°°°°° NONLINEAR STIFF system due to Ralson & Rabinowitz, p.232:"<<endl;
  double ral[] = {0,0};
  MyVec<double> RR_0(ral,ral+2);

  ODE<double,RalstonRabinowitz> IVP(0.,100., RR_0);
  
  size_t n_sps = 1000,//1000, //uniform h = .1
  newton_max = 13;
  double tolera = 1.e-12;
  

  cout<< "Implicit Euler = "<< IVP.implicit_euler(n_sps,newton_max,tolera,0.1);
  // cout<< "Explicit Euler = "<< IVP.explicit_euler(n_sps);

  cout<<"°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"  <<endl;
 
  }
  else if (example == 2){
 
    cout<< "a stiff example due to Hairer & Wanner,Vol II, p.3: " <<endl; 
   double hai[] = {1,0,0};
  MyVec<double> Hai_0(hai,hai+3);


  const int PRINTFLAG = 0; //no method of lines printing
  const bool SOLVERFLAG = true; //false;  //works for true and false
  //                                               (sparse LA) (dense LA)
  const bool SCALE = true; 
  const bool PATTERNCREATIONVIAAD = true; //true = AD, false = FD
  ODE<double,HairerRobertson,PRINTFLAG,SOLVERFLAG,SCALE,PATTERNCREATIONVIAAD> Stiff(.0,40.,Hai_0);
  
  cout<<endl<<"=========== a STIFF problem due to Hairer, Vol II, p.3 ========"<<std::endl;
  
  string name = "RobertsonReaction.dat";
  size_t n_steps = 1000000; //205; // 209
  size_t iter_max = 5;
  double TOL = 1e-6,
    rtol = 2.5e-6,
    hEstim = 1.e-03;

  cout<< endl << " Implicit Euler = "<<Stiff.implicit_euler(n_steps, iter_max, TOL, rtol,name,hEstim) <<endl;
  
  cout << "Note: to compare with Fig. 1.3 in Vol II, p.4 in Hairer/Wanner you should plot only the second species via '...using 1:3'"<< endl;
  cout<<"===================================================================="<<endl;
  GNUplot plotme;
  plotme("set key box left; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \"time t\"; set ylabel \"Y(t =T)\"; plot \"" + name + "\" using 1:2 with lines lw 3 title \"y_1\", \""+name+"\" using 1:3 with lines lw 3 title \"y_2\",\""+name+"\" using 1:4 with lines lw 3 title \"y_3\"");
 
  /*//works
    cout<<endl<< "Parameter"<<endl;
  Parameter<double,7> P("Parameter.dat");
  cout << P <<endl;
  */
 
  /*  //works
    cout<<endl<<"Test ntimes again = ";
  double ar = 1.5;
  double brd[] = {1., -0.5};
  cout<<ntimes<2>(brd[0]/(1+ar*brd[1]))<<endl;
  */
  }
  else if(example == 3){

   
    cout << "The oscillatory brusselator: "<<endl;
    double ibrus[] = {1,1};//{1,3};
  MyVec<double> y0_BRUSS(ibrus,ibrus+2);

  ODE<double,Brusselator> NonStiff(0.,20.,y0_BRUSS);//but oscillatory!!
  cout << endl<< "Solution:"<<endl;
 
  size_t  t_iter = 1000;//time steps

  //! the following 3 lines are used by the implicit Euler
  // size_t t_iter_impl = 50000;
  //double ntol = 1.e-3,  //0.3965389207; //8669200; //
  //  rtol = 2.5e-2;
  
  cout<<endl<< "================== SOLUTION BY DIFFERENT METHODS ============="<<endl;
  cout << setprecision(14) << endl;
    //double step = Abs(NonStiff.tend() - NonStiff.tstart())/t_iter;
  //cout <<endl << "Implicit Euler = " <<NonStiff.implicit_euler(t_iter_impl, 7, ntol, rtol) << endl;
  // cout << "Euler = "<<NonStiff.explicit_euler(t_iter) << endl;

  cout << "Heun = "<<NonStiff.heun(t_iter) << endl;

  // cout <<"PC Euler = "<<  NonStiff.pc_euler(t_iter) <<endl;
  
  //cout<<endl<<"M_PI = "<<M_PI<<endl; //works
  cout << endl;

  which_addition_is_used();
}

  
  else if (example == 4){
    cout << "An extended Davis Skodje Model:"<<endl; 
    
    typedef Constant<double> ConstantType;

    const size_t nos = ConstantType::numberOfSpecies,
      spacepts = ConstantType::spacePoints,
      fulldim = nos*spacepts;
    
    adonis_assert((nos == 2));  //for this example


    //============= Set time horizon ================
    double t0 = 0.,
      tf = 1.;
    
    //===============================================


    //FULL FULL FULL FULL FULL FULL FULL FULL FULL  FULL FULL FULL FULL FULL 
    Parameter<double,7> Para("Parameter.dat");

    MyVec<double> Ini(fulldim);
    
    const double a = Para[0],
      hx = 1./(spacepts-1); //1e-2;

    for(size_t i = 0; i < spacepts; ++i){
      Ini[i] = i*hx;
      
      Ini[spacepts+i] = Ini[i]/(1 + a* Ini[i]);  //unrepresented lies on SIM
    }
    //cout<< "Initial value vector = "<<Ini <<endl;
    
    ODE<double,DiffusionDevidedDifferences> FiDi(t0,tf,Ini);
    
    //FULL FULL FULL FULL FULL FULL FULL FULL FULL  FULL FULL FULL FULL FULL
    
    
    
    //1st APPROXIMATION 1st APPROXIMATION 1st APPROXIMATION 1st APPROXIMATION 
    MyVec<double> IFirst(spacepts);
    
    for(size_t i = 0; i < spacepts; ++i)
      IFirst[i] = i*hx;
    
    ODE<double, FirstApproximation> FirstApprox(t0,tf,IFirst);

    //1st APPROXIMATION 1st APPROXIMATION 1st APPROXIMATION 1st APPROXIMATION  


    //WITH CORRECTION WITH CORRECTION WITH CORRECTION WITH CORRECTION 
    MyVec<double> ICorr(spacepts);
    
    for(size_t i = 0; i < spacepts; ++i)
      ICorr[i] = i*hx;
    
    ODE<double, WithCorrection>        //original 
      //ODE<double, RPSecApproxCPA>    //well, does virtually the same
      Corrected(t0,tf,ICorr);

    //WITH CORRECTION WITH CORRECTION WITH CORRECTION WITH CORRECTION 

    

    // ==================== implicit Euler settings =========================
    size_t timest = 100, //originally: 100
      mxi = 6;
    
    double ntoler = 1e-4; 
    //=======================================================================
    

 
    cout << "***********************************************************"<<endl;
    CPUTime<double> cpuFull, cpu1st, cpuWC;
    WallTime<double> wallFull, wall1st, wallWC;

   
    
    
    wallFull.start();
    //time_t tm1 = time(0);
    cpuFull.start();
    wallFull.start();
    MyVec<double> IEU = FiDi.implicit_euler(timest, mxi, ntoler,1.e-2);
    cout<<" Implicit Euler = "<< IEU<<endl;
    cpuFull.stop();
    wallFull.stop();
    //time_t tm2 = time(0);
    //cout<< "Wall time prime = "<<difftime(tm2,tm1)<<endl;

    cpu1st.start();
    wall1st.start();
    MyVec<double> FAPP = FirstApprox.implicit_euler(timest, mxi, ntoler,1.e-2);
    cout<<endl<<" Backward Euler (1st approx.) = "<< FAPP <<endl;
    cpu1st.stop();
    wall1st.stop();
   
   
    cpuWC.start();
    wallWC.start();
    MyVec<double> CORR = Corrected.implicit_euler(timest, mxi, ntoler,1.e-2);
    cout<<endl<<" Backward Euler (with correction) = "<< CORR <<endl;
    cpuWC.stop();
    wallWC.stop();
   

    cout << endl << "Computational time(s):"<<endl;
    
    cout << "Full:   cpu = "<<cpuFull.get_time()<< "  wall = "<< wallFull.get_time() <<endl;
    cout << "1st:   cpu = "<<cpu1st.get_time()<< "  wall = "<< wall1st.get_time() <<endl;
    cout << "corr:   cpu = "<<cpuWC.get_time()<< "  wall = "<< wallWC.get_time() <<endl;
    
    
    
    //PLOT THE FOLLOWING FILES!!!!!!!!!!!!!!!
    //TO COMPARE WITH THE RESULTS FROM "REN/POPE, COMBUSTION AND FLAME 147 (2006),PP.  253" and for case 2, e.g., plot column 1 against 2 in each file 
    //cout<<"Solution of y_1 and y_2 at endtime t_end = "<< FiDi.tend() <<" along the spacial domain:"<<endl;
    
    //OUTPUT
    std::ofstream myo("Solution.dat",std::ios_base::out);
    std::ofstream first("1stapprox.dat",std::ios_base::out);
    std::ofstream wc("WithCorrection.dat",std::ios_base::out);
    std::ofstream z2file("z2calc.dat", std::ios_base::out);

    myo << "# space x           z_1            z_2"<<endl;

    first << "# space x           z_1"<<endl;
    
    wc << "# space x           z_1"<<endl;
    
    z2file <<" # space x           z_2"<<endl;

    // ==== needed  for delta u =====================
    double fprime = 0,
      x = 0,
      diffuse1 = 0,
      delta_u = 0;
    
    const double b = Para[1];
    const double c = Para[2];
    const double d = Para[3];
    const double e = Para[4];
    const double D1 = Para[5];
    const double eps = Para[6];
    
    DiffusiveTransport<double> D2(D1,e);
    //===============================================
    
    for(size_t i = 0; i < spacepts; ++i){
      x = i*hx;
      
      myo << " "<< x << "    "<<IEU[i] << "    "<<IEU[i+spacepts] <<endl;
      
      first <<  " "<< x << "    "<<FAPP[i] << endl; 
      
      wc <<  " "<< x << "    "<< CORR[i] 
	//Corrected.original_function()[i]  //this does not work and gives crap
	 <<endl; 
      

     
	

      
      //=========== needed for determination of delta u ====================
      //note that in the case of left bdy (i = 0) and right bdy (i = spacepts-1)
      //we have zero conditions, thus delta u is zero
      if((i > 0) && ( i < spacepts-1)){ 
	fprime =  1./(ntimes<2>(1 + a*CORR[i]));
      
	delta_u = eps/(c*fprime + 1)*( d*CORR[i]*fprime -  CORR[i]/(ntimes<2>(1+b*CORR[i])) - fprime*diffuse1 + ( D2(x+0.5*hx)*(CORR[i+1]/(1.+a*CORR[i+1]) - CORR[i]/(1.+a*CORR[i])) - D2(x-0.5*hx)*(CORR[i]/(1.+a*CORR[i]) - CORR[i-1]/(1.+a*CORR[i-1])) )/(ntimes<2>(hx)) );
      }
      //======================================================================

      cout << "delta_u = "<< delta_u << endl;


      z2file << " "<< x << "   "<< //CORR[i]/(1+a*CORR[i])
	//Corrected.original_function().z2_index(i) //this does not work and gives crap, but in conjunction with the crap version of 'wc', the two crap results aggree
	
	//!z2(x,t) = z2^M(x,t) + delta z2
	( CORR[i]/(1+a*CORR[i]) + delta_u )
	     << endl; 
    }

    //int pcase = 2;
    //string toeps = "set terminal postscript eps color enhanced; set output 'cpa_with_z2_case_" + my_function_collection::num2str(pcase) + ".eps'";

    GNUplot pres;
    pres(//toeps + ";"+
	 "set key box left; set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"x\"; set ylabel \"z(x,t=1)\"; set title \"1st and 2nd approx. versus full model\";plot  \"./1stapprox.dat\" using 1:2 with lines  lc rgb \"violet\" lw 2 title \"1st approx.\", \"./WithCorrection.dat\" using 1:2 with points lc rgb \"blue\" lw 1 pt 7 title \"2nd approx.\", \"./Solution.dat\" using 1:2 with lines lc rgb \"red\" lw 2 title \"z_1 full\", \"z2calc.dat\" using 1:2 with points lc rgb \"orange\" lw 1 pt 7 title \"z_2 recovered from reduction\" , \"./Solution.dat\" using 1:3 with lines lc rgb \"green\" lw 2 title \"z_2 full\"");

    myo.close();
   
    first.close();
    wc.close();
    z2file.close();

    cout << "***********************************************************"<<endl;
  
    cout<<endl<<"Stiffness of source term:"<<endl;
    RPDavisSkodje<double> RPSource(2);
    StiffnessDetector<double,RPDavisSkodje> StiffRP(RPSource);
  
    MyVec<double> eval(2);
    eval[0] = 0.64;   eval[1] = 0.39023681084;
   
    cout << "Ren Pope Stiffness rate = "<<StiffRP.check(0.,1.,eval)<<endl;
    cout<< "Parameter: "; 
    RPSource.parameter();
    cout<< endl;

  }


  else if(example == 5){
    cout << " WITH CORRECTION WITH CORRECTION WITH CORRECTION WITH CORRECTION"<<endl;
    
    typedef Constant<double> ConstantType;

    const size_t spacepts = ConstantType::spacePoints;
    
    const double hx = 1./(spacepts-1); 


    MyVec<double> ICorr(spacepts);
    
    for(size_t i = 0; i < spacepts; ++i)
      ICorr[i] = i*hx;
    
    ODE<double, WithCorrection> Corrected(0.,1,ICorr);

    
    size_t timest = 100, //originally: 100
      mxi = 6;
    
    double ntoler = 1e-4; 
    
    CPUTime<double> cpuIE;
    cpuIE.start();
    MyVec<double> COR = Corrected.implicit_euler(timest,mxi,ntoler,1.e-2);
    cout<<endl<<" Backward Euler (with correction) = "<< COR <<endl;
    cpuIE.stop();
    
    std::ofstream wc("WithCorrection.dat",std::ios_base::out);
    
    
    wc << "# space x           z_1"<<endl;
  
    for(size_t i = 0; i < spacepts; ++i){
      wc <<  " "<< i*hx << "    "<<COR[i] << endl; 
    }

    wc.close();
    

    cout << " WITH CORRECTION WITH CORRECTION WITH CORRECTION WITH CORRECTION"<<endl;

  }

  
  else if(example == 6){
    //================= TODO: select implicit EULER ===================
    bool ByHand = false;   //true: hand-coded Euler, false: ode.hh Euler
    //=================================================================
    //initial value 
    MyVec<double> S(2);
    S[0] = 1; 
    S[1] = 0;  //time

    if(!ByHand){
      const int PRINTFLAG = 0; //no method of lines printing
      const bool SOLVERFLAG = true; //false;  //works for true and false
      //                                               (sparse LA) (dense LA)
      const bool SCALE = true; 
      const bool PATTERNCREATIONVIAAD = true; //true = AD, false = FD

      ODE<double, Non2Autonomous,PRINTFLAG,SOLVERFLAG,SCALE,PATTERNCREATIONVIAAD> Petz(0.,3.075,S);
    
      //cout<<"Order of Method "<<Constant<double>::accuracyOfMethod<<endl;
      CPUTime<double> time;
      time.start();
      string name;
      //cout << setprecision(16) << endl;
    
    
      name = "Impl_Euler_AscherPetzold.dat";
    
      size_t steps_impl = 200000000,
	maxi = 12;
      double rto = 1e-4,
	stpzctrl = 1.e-3,
	hestim = 0.000125, //0.01;   //see Ascher/Petzold 100·h_n <= 1 
	Cscal = 1.111111;
    
    
      cout<< "IMplicit method = "<< setprecision(9)<< Petz.implicit_euler(steps_impl,maxi,rto,stpzctrl,name,hestim,Cscal) 
	  <<endl;
      cout << "Implicit Euler -- calculation time: "
	   << time.stop() <<endl;
    
      GNUplot plotme;
      plotme("set key box left; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \"time t\"; set ylabel \"Y(t =T)\"; plot \"" + name + "\" using 1:2 with lines lw 3");

      which_addition_is_used();
  
    }
    else{
      int maxNewt = 5;
      double dt = 0.000125,
	tol = 1.e-04,
	tend(3.075),
	time(0.);
  
      MyVec<double> U(S),Uprev(S), b(2), Delta(2);

      Non2Autonomous<double> F(2);
      Dune::FieldMatrix<double,2,2> DG, Inverse;
      ofstream of("ByHand.dat");
      of << "# time       U[0]     U[1]"<<endl;
      int count(0);
      
      CPUTime<double> cpu;
      cpu.start();
      while(time < tend){
	cout << count << ".)   time: "<< time << "   dt = "<< dt << endl;
	of << setprecision(12) << time << "  "<< U[0] << "   "<< U[1] << endl; 
	//Newton iteration
	U = Uprev;
	for(int i = 0; i < maxNewt; ++i){
	  b = -(U - Uprev - dt*F(U));
	
	  DG[0][0] = 1 - dt*(-100); DG[0][1] = -dt*100*cos(U[1]);
	  DG[1][0] = 0;             DG[1][1] = 1.;

	  Dune::FMatrixHelp::invertMatrix(DG,Inverse);
	
	  matrix_vector_product(Delta,Inverse,b);
	  U += Delta;

	  if(Norm<'2',double>::norm(Delta) <= tol){
	    break;
	  }

	  if(i == maxNewt-1){
	    ADONIS_INFO(Information, "Newton iteration did not converge in "<< maxNewt << " iterations.");
	  }
	}//end Newton iteration
	Uprev = U;

	time += dt;
	count++;
      }
      cpu.stop();
      cout << "HAND-CODED EQUIDISTANT IMPLICIT EULER" <<endl;

      GNUplot show;
      show("set xlabel \"Time \"; set ylabel \"U[1]\";set title \"Ascher Petzold example direct IE, equidistant dt = "+num2str(dt)+" \";plot 'ByHand.dat' using 1:2 with lines lw 2 notitle");
    
    }

    // //!works
    // MyVec<double> eval(2);
    // eval <<= 0.5, 0.75;

    // Non2Autonomous<double> Fu(2);

    
    // cout << "F(x) = "<< Fu(eval) << "   ||F(x)||_2 = "<< Norm<'2',double>::norm(Fu(eval)) << endl;
    
    /* //try explicit Euler
       size_t steps = 100, //!to use implicit trapezoidal increase stepsize (e.g.  steps = 200), otherwise method does not converge  
      // works very good with, say 3950 or 5000 steps (result o.k.)
       time.reset();
    time.start();
    cout<< "Explicit Euler = "<< Petz.explicit_euler(steps) << "     (truely horrendous!) "<<endl;
    cout<<"Note: Plot y[0] against time and compare to Fig. 3.4, p. 49 of [Asch/Petz] "<<endl;
  
    cout << endl;

    cout << "Explicit Euler with "<<steps<<" equidist. steps"  << time.stop() <<endl;
    */
    
    /*
     //… WEIRD … WEIRD … WEIRD … WEIRD … WEIRD … WEIRD … WEIRD … WEIRD …
    //… WEIRD … WEIRD … WEIRD … WEIRD … WEIRD … WEIRD … WEIRD … WEIRD …
    cout << endl << "WEIRD EXPLICIT METHOD 4 STIFF ODEs:"<<endl<<endl;
    name = "Ex4Stiff_Petzold.dat";
    //these setting perform relatively well, i.e. y(t=T) ~ [0.079589, 3.07192] at the cost of 559 iterations
    //============ TODO ======
    double TOL = 1.e-05,   
      tol = 1.e-05,         
      hEstim = 0.0125,     
      St = 0.3,
      c = 0.98;  //close to one

    size_t maxit = 4;
    
    //bounds on stepsize
    double kmin = 1.25e-06,
      kmax = 1.25;
      //=======================


     time.reset();
     time.start();
     cout << "Method of Eriksson-Johnson-Logg = "<< Petz.explicit_4_stiff(TOL,tol,hEstim,St,c,maxit,name,kmin,kmax) << endl;
     cout << "time: " << time.stop() <<endl;
*/    

     
    /* //works
    cout << endl<< "Test Norm"<<endl;
    typedef ChooseNorm<double,ExprTmpl::MyVec,double,Constant<double>::whatNorm>::NormType NormType;
    MyVec<double> v(3);
    v[0] = 0.5; v[1] = -1.5; v[2] = 2.75;
    cout << "1-Norm = "<< OneNorm<double,ExprTmpl::MyVec,double>::norm(v) <<endl;
    cout<<"2-Norm = "<< NormType::norm(v) <<endl;
    cout<<"inf-Norm = "<< InfinityNorm<double,ExprTmpl::MyVec,double>::norm(v) <<endl;
    */
  }

  else if (example == 7){
    MyVec<double> S(2);
    S[0] = 1; 
    S[1] = 1;  //time

    const int PRINTFLAG = 0; //no method of lines printing
    const bool SOLVERFLAG = true; //false;  //works for true and false
    //                                               (sparse LA) (dense LA)
    const bool SCALE = true; 
    const bool PATTERNCREATIONVIAAD = false; //true = AD, false = FD

    double tend = 1.0;
    ODE<double,SimpleStiffSys,PRINTFLAG,SOLVERFLAG,SCALE,PATTERNCREATIONVIAAD> ivp(0.,tend,S);
    
    //cout<<"Order of Method "<<Constant<double>::accuracyOfMethod<<endl;
    CPUTime<double> time;
    time.start();
    string name;
    //cout << setprecision(16) << endl;
    
    
    name = "Simple nonlinear stiff test system";
    
    size_t steps_impl = 200000000,
      maxi = 12;
    double rto = 1e-4,
      stpzctrl = 1.e-3,
      hestim = 0.000125, //0.01;    
      Cscal = 1.111111;
    
    MyVec<double> numsol = ivp.implicit_euler(steps_impl,maxi,rto,stpzctrl,name,hestim,Cscal),
      exactsol(2);
    exactsol <<= exp(-2*tend),exp(-tend);

    cout<< "IMplicit method = "<< setprecision(9)<<  numsol
	<<endl << "exact solution =  "<< exactsol << endl;
    MyVec<double> err(2);
    err = numsol - exactsol;
    cout << "||Y(tend)-y(tend)|| = "<< Norm<'2',double>::norm(err) << endl;
    cout << "Implicit Euler -- calculation time: "
	 << time.stop() <<endl;
    

    GNUplot plotme;
    plotme("set key box right; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \"time t\"; set ylabel \"Y(t)\"; plot \"" + name + "\" using 1:2 with lines lw 3 title \"Y_1\", \"" + name + "\" using 1:3 with lines lw 3 title \"Y_2\"");

    which_addition_is_used();
  
    /*
    cout << "Simple nonlinear test mechanism" << endl;
    
    double t0 = 0.,
      tend = 300;
    double hEstim = 0.75, //0.015;
      Cscal = 1.111,
      hmin = 0.5,
      hmax = 1.;

    MyVec<double> x0(3);
    x0[0] = 0.35;
    x0[1] = 0.25;    //x[0] + x[1] + x[2] = 1.0
    x0[2] = 0.4;

    ODE<double,ToyMechanism> ode(t0,tend,x0);

    MyVec<double> Sol = ode.implicit_euler(200000,7,1.e-05,1.e-05,"ToyMechanism.dat",hEstim,Cscal,hmin,hmax);

    cout << "Solution = "<< Sol << endl;

    GNUplot show;
    show("plot 'ToyMechanism.dat' using 1:2 with lines lw 2 title 'y1', 'ToyMechanism.dat' using 1:3 with lines lw 2 title 'y2', 'ToyMechanism.dat' using 1:4 with lines lw 2 title 'y3'");
    
    GNUplot plot;
    plot("splot \"ToyMechanism.dat\" using 2:3:4");
*/
    
  }
  else if (example == 8){
    cout << "Stiff example taken from Burden Faires, § 5.11, p. 335:"<<std::endl;

    const int dim = 3; 
    double t_0 = 0.,
      t_end = 1.;

    MyVec<double> Y_0(dim);
    Y_0[0] = 4./3.;
    Y_0[1] = 2./3.;
    Y_0[2] = t_0;

    
    ODE<double, StiffBurdenFaires> BF(t_0, t_end, Y_0);

    CPUTime<double> watch;
    watch.start();
    cout << "Impl. Euler = "<<BF.implicit_euler(20, 7, 1e-5,1.e-2);
    watch.stop();
  
  }

  //=================== TEST FULL AND REDUCED H2 MODEL ========================
  else if (example == 9){
    MyVec<double> v(6);

    v[0] = -0.05; v[1] = -1.5; v[2] = 2.75; v[3] = -0.75; v[4] = 3.085; v[5] = -1.125;

    cout << "v = "<< v <<endl;
    cout << "\"Absify\" v:" <<endl;
    absify(v);
    cout << "v = "<< v <<endl;
    
  }

  
  else if (example == 10){
  
    cout << "THERMODYNAMIC STUFF" << endl;
    
    

    //NEW · NEW · NEW
    
    //O
    double thermo[] = { 2.56942078E+00,-8.59741137E-05, 4.19484589E-08,-1.00177799E-11, 1.22833691E-15, 2.92175791E+04, 4.78433864E+00, 3.16826710E+00,-3.27931884E-03, 6.64306396E-06,-6.12806624E-09, 2.11265971E-12, 2.91222592E+04, 2.05193346E+0,
  //O2
			3.28253784E+00,1.48308754E-03,-7.57966669E-07, 2.09470555E-10,-2.16717794E-14, -1.08845772E+03, 5.45323129E+00, 3.78245636E+00,-2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12,-1.06394356E+03, 3.65767573E+00,

//H
			2.50000001E+00,-2.30842973E-11, 1.61561948E-14,-4.73515235E-18, 4.98197357E-22,    
			2.54736599E+04,-4.46682914E-01, 2.50000000E+00, 7.05332819E-13,-1.99591964E-15,    
			2.30081632E-18,-9.27732332E-22, 2.54736599E+04,-4.46682853E-01,

//OH
3.09288767E+00, 5.48429716E-04, 1.26505228E-07,-8.79461556E-11, 1.17412376E-14,    
			 3.85865700E+03, 4.47669610E+00, 3.99201543E+00,-2.40131752E-03, 4.61793841E-06,    
			-3.88113333E-09, 1.36411470E-12, 3.61508056E+03,-1.03925458E-01,
 
//H2	
3.33727920E+00,-4.94024731E-05, 4.99456778E-07,-1.79566394E-10, 2.00255376E-14,   
			  -9.50158922E+02,-3.20502331E+00, 2.34433112E+00, 7.98052075E-03,-1.94781510E-05,    
			2.01572094E-08,-7.37611761E-12,-9.17935173E+02, 6.83010238E-01,

//HO2
    4.01721090E+00, 2.23982013E-03,-6.33658150E-07, 1.14246370E-10,-1.07908535E-14,    
		      1.11856713E+02, 3.78510215E+00, 4.30179801E+00,-4.74912051E-03, 2.11582891E-05,    
			-2.42763894E-08, 9.29225124E-12, 2.94808040E+02, 3.71666245E+00,

//H2O2
4.16500285E+00, 4.90831694E-03,-1.90139225E-06, 3.71185986E-10,-2.87908305E-14,    
			    -1.78617877E+04, 2.91615662E+00, 4.27611269E+00,-5.42822417E-04, 1.67335701E-05,    
			-2.15770813E-08, 8.62454363E-12,-1.77025821E+04, 3.43505074E+00,

 //H2O
 3.03399249E+00, 2.17691804E-03,-1.64072518E-07,-9.70419870E-11, 1.68200992E-14,    
			   -3.00042971E+04, 4.96677010E+00, 4.19864056E+00,-2.03643410E-03, 6.52040211E-06,    
			-5.48797062E-09, 1.77197817E-12,-3.02937267E+04,-8.49032208E-01,

 //N2
			0.02926640E+02, 0.14879768E-02,-0.05684760E-05, 0.10097038E-09,-0.06753351E-13,    
			 -0.09227977E+04, 0.05980528E+02, 0.03298677E+02, 0.14082404E-02,-0.03963222E-04,    
			 0.05641515E-07,-0.02444854E-10,-0.10208999E+04, 0.03950372E+02   };

   
    double tlh[] = {200,3500,1000,  //O
		    200,3500,1000,  //O2
		    200,3500,1000,  //H
		    200,3500,1000,  //OH
		    200,3500,1000,  //H2
		    200,3500,1000,  //HO2
		    200,3500,1000,  //H2O2
		    200,3500,1000,  //H2O
		    300,5000,1000   //N2
    };
    
    MyVec<double> temps(tlh, tlh+3*9);
    
  
    double lt = 856.135;
  
    NASA7CoefficientPolynomial<const double*,MyVec<double>::const_iterator> def;

    typedef NASA7CoefficientPolynomial<const double*,MyVec<double>::const_iterator> ThermoType;
    
    typedef ThermoData4Mechanism<double,9> DataType;


    

    size_t species = 3;
    ThermoType tc(DataType::nspec,thermo,temps.begin());
    cout << endl<<"Low temperatures:"<<endl;
    cout << "heat capacity C_p = "<<tc.heat_capacity(species,lt) <<endl;
    cout << "enthalpy H_T = " << tc.enthalpy(species,lt) << endl;
    cout << "entropy S_T = "<< tc.entropy(species,lt) << endl;
    
    double ht = 1765.78;

    cout << endl<<"High temperatures:"<<endl;
    cout << "heat capacity C_p = "<<tc.heat_capacity(species,ht) <<endl;
    cout << "enthalpy H_T = " << tc.enthalpy(species,ht) << endl;
    cout << "entropy S_T = "<< tc.entropy(species,ht) << endl;
    
    /* taken from '/massactionkinetics/data/thermochemicaldata.hh'
    //store only ni' 
    //O O2 H OH  H2  HO2  H2O2  H2O  N2
    size_t sm[] = {2,0,0,0,0,0,0,0,0,
		   0,1,0,0,0,0,0,0,0,
		   1,0,1,0,0,0,0,0,0,
		   0,0,0,1,0,0,0,0,0,
		   1,0,0,0,1,0,0,0,0,
		   0,0,1,1,0,0,0,0,0,
		   1,0,0,0,0,1,0,0,0,
		   0,1,0,1,0,0,0,0,0,
		   1,0,0,0,0,0,1,0,0,
		   0,0,0,1,0,1,0,0,0,
		   0,1,1,0,0,0,0,0,0,
		   0,0,0,0,0,1,0,0,0,
		   0,2,1,0,0,0,0,0,0,
		   0,1,0,0,0,1,0,0,0,
		   0,1,1,0,0,0,0,1,0,
		   0,0,0,0,0,1,0,1,0,
		   0,1,1,0,0,0,0,0,1,
		   0,0,0,0,0,1,0,0,1,
		   0,1,1,0,0,0,0,0,0,
		   1,0,0,1,0,0,0,0,0,
		   0,0,2,0,0,0,0,0,0,
		   0,0,0,0,1,0,0,0,0,
		   0,0,2,0,1,0,0,0,0,
		   0,0,0,0,2,0,0,0,0,
		   0,0,2,0,0,0,0,1,0,
		   0,0,0,0,1,0,0,1,0,
		   0,0,1,1,0,0,0,0,0,
		   0,0,0,0,0,0,0,1,0,
		   0,0,1,0,0,1,0,0,0,
		   1,0,0,0,0,0,0,1,0,
		   0,0,1,0,0,1,0,0,0,
		   0,1,0,0,1,0,0,0,0,
		   0,0,1,0,0,1,0,0,0,
		   0,0,0,2,0,0,0,0,0,
		   0,0,1,0,0,0,1,0,0,
		   0,0,0,0,1,1,0,0,0,
		   0,0,1,0,0,0,1,0,0,
		   0,0,0,1,0,0,0,1,0,
		   0,0,0,1,1,0,0,0,0,
		   0,0,1,0,0,0,0,1,0,
		   0,0,0,2,0,0,0,0,0,
		   0,0,0,0,0,0,1,0,0,
		   0,0,0,2,0,0,0,0,0,
		   1,0,0,0,0,0,0,1,0,
		   0,0,0,1,0,1,0,0,0,
		   0,1,0,0,0,0,0,1,0,
		   0,0,0,1,0,0,1,0,0,
		   0,0,0,0,0,1,0,1,0,
		   0,0,0,1,0,0,1,0,0,
		   0,0,0,0,0,1,0,1,0,
		   0,0,0,0,0,2,0,0,0,
		   0,1,0,0,0,0,1,0,0,
		   0,0,0,0,0,2,0,0,0,
		   0,1,0,0,0,0,1,0,0,
		   0,0,0,1,0,1,0,0,0,
		   0,1,0,0,0,0,0,1,0 
    };
*/
    
    //works
    //typedef StoichiometricMatrix<size_t*> StoichType;
    //StoichType SM(nos,nor,sm);
   
    typedef StoichiometricMatrix<const DataType::index_type*> StoichType;
    
    StoichType SM(DataType::nspec,DataType::nreac,DataType::stoichiometric_matrix());


    cout << SM.nu_forward(2,2) << "   " << SM.nu_reverse(2,2) << endl<<
      SM.nu_forward(2,3) << "   " << SM.nu_reverse(2,3) << endl;

    
   
   
    cout << endl<< "FORWARD REACTIONS:"<<endl;
    
    double AbetaEa[] = {1.200E+17, -1.000,     0.00,
			5.000E+17, -1.000,     0.00, 
			3.870E+04, 2.700,   26.1918, 
			2.000E+13,   .000,     0.00,
			9.630E+06,  2.000,   16.736,
			2.800E+18, -0.860,     0.00,        
			2.080E+19, -1.240,     0.00,        
			11.26E+18, -.760,      0.00,        
			2.600E+19, -1.240,     0.00,        
			2.650E+16, -0.6707,  71.2995,       
			1.000E+18, -1.000,     0.00,        
			9.000E+16, -0.600,     0.00,        
			6.000E+19, -1.250,     0.00,        
			2.200E+22, -2.000,     0.00,        
			3.970E+12,  0.000,   2.80746,        
			4.480E+13,  0.000,   4.46951,        
			0.840E+14,  0.000,   2.65684,        
			1.210E+07,  2.000,   21.7568,        
			1.000E+13,  0.000,   15.0624,        
			2.160E+08,  1.510,   14.3511,        
			7.400E+13, -0.370,     0.00,             
			3.570E+04,  2.400,   -8.8282,        
			1.450E+13,  0.000,   -2.029,        
			2.000E+12,  0.000,   1.78657,        
			
			1.700E+18,  0.000,   123.051,       
   
			1.300E+11,  0.000,   -6.8199,        
			
			4.200E+14,  0.000,   50.208,       
			
			0.500E+16, 0.000,   72.5087 
};

    /*  //taken from '/massactionkinetics/data/thermochemicaldata.hh'
    int troewtb[] = {0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			1,
			0,
			0,
			0,
			0,
			0,
			0,
			0, 
//now comes information about 3rd bodies:
			0,
			1,
			-1,
			-1,
			-1,
			2,
			-1,
			-1,
			-1,
			-1,
			3,
			-1,
			-1,
			4,
			-1,
			-1,
			-1,
			-1,
			-1,
			-1,
			5,
			-1,
			-1,
			-1,
			-1,
			-1,
			-1,
			-1
};
*/
   
    const bool Aalter = true,
      Ealter = true;
    const char storageForm = 'a';

    cout << endl<< "Test Transformation:"<<endl;
   
    //TROE Reactions -- stored LOW|TROE
   
    //must be transformed into right unit first
    double troecoeff[DataType::ntroereac*7] ={2.300E+18, -0.900,   -7.1128,      //LOW
				     //alpha   T***   T*       T**
				     0.7346,   94.0,  1756.0,  21.68};  //TROE


    typedef TransformArrheniusDataIntoRightUnits<double,const double*,const int*,StoichType,const double*,Aalter,Ealter,storageForm> ConverterType;
    
   

     
    
    //take from '/massactionkinetics/data/thermochemicaldata.hh'
    /*
    double collisionEfficiencies[] = {1,
				    1,
				    1,
				    1,
				    2.4,
				    1,
				    1,
				    15.4,
				    1,
				    1,
				    1,
				    1,
				    1,
				    2,
				    1,
				    1,
				    6,
				    1,
				    1,
				    0,
				    1,
				    1,
				    1,
				    1,
				    1,
				    0,
				    1,
				    1,
				    1,
				    1,
				    1,
				    0,
				    1,
				    1,
				    0,
				    1,
				    1,
				    1,
				    1,
				    1,
				    0.73,
				    1,
				    1,
				    3.65,
				    1,
				    1,
				    1,
				    1,
				    1,
				    2,
				    1,
				    1,
				    6,
				    1

    };
    */

     //================
    typedef ThermoData4Mechanism<double,9> DataType;
    //typedef DataType::int_type_pointer int_type_pointer;
    //================

    cout << endl<< "REVERSE reaction rates:"<<endl;

    //create troe index, e.g. at reaction i = 20 there is the 1st index 0
    //                                    i=  37              2nd       1... 
    TroeIndex TIndex;
    TIndex.create(DataType::nreac,DataType::troewtb());

    cout <<endl<< "TIndex = "<< TIndex << endl;
    cout << "Number of TROE reactions: "<< TIndex.number_of_TROE_reactions() << endl;
  cout << "TIndex[20] = "<< TIndex[20] << endl;
  //cout << "TIndex[6] = "<< TIndex[6] << endl;


  typedef DataType::value_type d_type;
  typedef DataType::int_type i_type;
  typedef ForwardReactionRates<ConverterType::ContainerType::const_iterator,const d_type*, const i_type*, const d_type*,storageForm> ForwardRatesType;
  ForwardRatesType FRR(DataType::nspec,DataType::nreac,ConverterType::instance(DataType::nreac,AbetaEa,DataType::troewtb(),SM,DataType::ntroereac,troecoeff).convert().begin(),DataType::temperature_bounds(), DataType::ntroereac, DataType::troewtb(),DataType::collision_efficiencies(),&TIndex);

  double temper =  1500.; //1235.067;
  size_t reacatinterest = 5;
  cout << "k_f,"<< reacatinterest << " = "<< FRR.k_forward(reacatinterest,temper) <<endl;
   
  MyVec<double> conc(DataType::nspec);
    random_container(conc);
    cout << "random concentrations = "<< conc <<endl;

    cout <<endl<< "----FORWARD RATES:"<<endl << "  T = "<<temper << " K"<<endl;
    for(size_t i = 0; i< DataType::nreac; ++i)
      cout << "k_f,"<<i<< " = "<< FRR.k_forward(i,temper,conc) <<endl;



    //CppAD::AD<double> 
    double tprt =  temper; //906.78;
    typedef ReverseReactionRates<ForwardRatesType,StoichType,ThermoType> ReverseRatesType;

    ReverseRatesType RRR(&FRR,&SM,&tc);

    //MyVec<CppAD::AD<double> > 
    MyVec<double> concentration(DataType::nspec);
    concentration[0] = 0.75;  concentration[1] = 12.5;
    concentration[2] = 2.45;  concentration[3] = 0.5;
    concentration[4] = 0.;  concentration[5] = 3.85;
    concentration[6] = 1.9;  concentration[7] = 0;
    concentration[8] = 5.22;

    cout <<endl<< "----REVERSE RATE:"<<endl << "  T = "<<tprt << " K"<<endl;
    for(int i = 0; i< DataType::nreac; ++i){
      cout << "k_{rev,"<<i<<"} = "<< RRR.k_reverse(i,tprt,concentration) <<endl;
    }

   
   
    cout <<endl << "Test reverse TROE reaction:"<<endl;

    size_t troe_reac = 20;
   

    

    cout << "k_{f,uni} = "<<FRR.k_forward(troe_reac,tprt,concentration) << endl;
 
  

    cout << "++++ nspec = "<< SM.number_of_species() << "   nreac = "<< FRR.number_of_forward_reactions() << endl;
    
    //========================================================
    cout<<endl<<endl<<"BUILD CHEMICAL RHS"<<endl;
    
    typedef double Type;            //CppAD::AD<double> 
      
    MyVec<Type> Uappx(DataType::nspec+1); //nos+1 is the temperature

    random_container(Uappx);
    Uappx[DataType::nspec] = 867.19;

    size_t reac = 20;

    cout << "Uappx = "<< Uappx << endl;
    

    //ThermoType nas;
    //nas = tc;
    //ForwardRatesType fwddef;
    //fwddef = FRR;

    //const ForwardRatesType* fwdptr = &FRR; //o.k.
    
    ReverseRatesType revdefault;
    
    revdefault = RRR;
  
      
    //works -- original implementation
    //BuildChemicalSourceTerm<ForwardRatesType,ReverseRatesType> BCS(&FRR,&RRR);
    CommonIndexInjector Cii(DataType::nspec,DataType::nreac); 
    BuildChemicalSourceTerm<ForwardRatesType,ReverseRatesType,CommonIndexInjector> BCS;
    
    BCS.initialize(&FRR,&RRR,&Cii);
    
    cout << "BCS: nspec = "<< BCS.number_of_species()<<std::endl;
    
    cout << "High precision output:"<< setprecision(16) << endl;
    cout << "eff3rd body = "<< BCS.effective_3rd_body_total_mixture_concentration(reac,Uappx) <<endl;


    size_t spec = 3;
    cout << "q_"<<reac<<" = "<< BCS.rate_of_progress(reac,Uappx[DataType::nspec],Uappx) <<endl;
    //                                                   species  temperature  u
    cout <<"omega_"<<spec<< " = "<< BCS.net_formation_rate(spec,Uappx[DataType::nspec],Uappx) << endl; 
    
    


    //gives some output
    cout << endl<<" TEST automatically built chemical source term:"<<endl;
    H2Gri<Type> Fun(DataType::nspec+1);
    cout << "dim = "<< Fun.dim()<<endl;
    // cout<< "rev-address = "<< Fun.rev_address()<<endl;
    cout << "nspec: "<<Fun.nspec() <<endl; //works after {static...}
    //start values in specific moles mol·kg
    double st[] = {0.0,7.0702437,0.0,0.0,14.140487,0.0,0.0,0.0,26.584116,//spec
		   1500 };   //temp
    
    MyVec<Type> Start(st,st+DataType::nspec+1);
    cout << "Input in specific moles:"<<endl;
    cout << "Start = "<< Start << endl; 
    

    cout<< "Fun(Start) = "<< Fun(Start) <<endl;
    
    //================= ODE ================================
    cout << "ODE integration:"<<endl;
    MyVec<double> Init(st,st+DataType::nspec+1);

    //END TIME
    double tf = 3.0e-4;

    ODE<double,H2Gri> ivp(0.,tf,Init); 

    MyVec<double> Sol;

    string ieFile; 

    CPUTime<double> cpu;
    cpu.start();
    
    string nme;

    //============= TODO: select integration method =====================
    bool useIE = true;   //true = use IE, false = try expl. method
    //===================================================================

    if(useIE){
      nme = "Implicit Euler";
      //implicit euler
      ieFile = "IE_H2Gri.dat";
      
      size_t msteps = 2000000,
	mNewt = 5;
      double ntol = 1.e-8,    //accept Newton as converged
	atol = 1.e-8,         //tol for stszctrl
	hEstim = 0.0000045,
	Cscal = 1.;
      
      //these are also worthwhile to play with
      double hmin = 1.125e-07,  //maybe method behaves 'quasi-equidistant'
	hmax = 4.15e-05;

      Sol = ivp.implicit_euler(msteps,mNewt,ntol,atol,ieFile,hEstim,Cscal,hmin,hmax);
    }
    else{
      
      //========================================================================
      //===== WEIRD · WEIRD · WEIRD · WEIRD · WEIRD · WEIRD · WEIRD ============
      //========================================================================
      //CAUTION : 
      //---------
      //         does not work very reliably  !!!!!!!!!! 
     
      nme = "Expl. method 4 stiff ODEs";
      ieFile = "ex4stiff_H2Gri.dat";
      double errTol = 1.e-8,
	fixedPtTol = 1.e-8,
	k0 = 0.0000012,
	St = 1.,
	cs = 0.9995;    //~ 1
      
      size_t mxit = 4;         //max number of function iterations

      double kmin = 1.12e-08,
	kmax = 1.e-06;
 
    
      Sol = ivp.explicit_4_stiff(errTol,fixedPtTol,k0,St,cs,mxit,ieFile,kmin,kmax);

      /*
      //test EXPLICIT EULER -- for the sake of completeness
      //works with relatively small stepsize, but the result is, kindly spoken, really bad
      nme = "explicit euler";
      ieFile = "ExEuler_H2Gri.dat";
      
      Sol = ivp.explicit_euler(1000000,ieFile);
      */
    }


    

    //print solution in formatted output -- note that we don't consider 
    // N2 here (the 9th species) for it is 
    GNUplot plotMePal;
    plotMePal("scale1(t) = 1.e-2*t; scale2(t) = 1.e-03*t; set yrange [0:15]; set title \""+ nme + " -- Order of method: " + num2str(ivp.order_of_method()) + "\" ; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \" time t \"; set ylabel \" U(t=T)\"; plot for[i=2:9] \"" + ieFile + "\" using 1:i with lines lw 1 notitle, \""+ ieFile + "\" using 1:(scale1($10)) with line lw 1 notitle, \"" + ieFile + "\" using 1:(scale2($11)) with lines lw 2 lc rgb \"dark-orange\" title \"T(t) (x 1.e-03)\"");
    

    double elti = cpu.stop();
   


    time_in_human_readable_form(elti);

    cout << endl<< "Solution = "<< Sol <<endl;
}



  

  else if (example == 11){
    
    cout << "TEST CLASSES"<<endl;

    typedef double Type;
    typedef ThermoData4Mechanism<Type,9> DataType;


    
    Type t3_high[]={ 3.09288767E+00, 5.48429716E-04, 1.26505228E-07,-8.79461556E-11, 1.17412376E-14,    
		  3.85865700E+03, 4.47669610E+00};

    Type t3_low[]={3.99201543E+00,-2.40131752E-03, 4.61793841E-06,    
		      -3.88113333E-09, 1.36411470E-12, 3.61508056E+03,-1.03925458E-01};

    typedef NASA7CoefficientPolynomial<const Type*, const Type*> NasaPolyType; 
    NasaPolyType nasa(DataType::nspec,DataType::thermo(), DataType::temperature_bounds());
    
    size_t spec = 3;
    Type T_low = 857.012;
    cout<<"____________________________________"<<endl;
    cout << "LOW temp: T_low = "<< T_low<<endl;
    cout<<"____________________________________"<<endl;
    cout << "heat_capa_low = "<< nasa.heat_capacity(spec,T_low)<<endl;
    cout<< "cross check: "<<heat_capacity(T_low,t3_low) <<endl;

    cout << "enthalpy_low = "<<nasa.enthalpy(spec,T_low)<<endl;
    cout<< "cross check: "<<enthalpy(T_low,t3_low) <<endl;

    cout << "entropy_low = "<<nasa.entropy(spec,T_low)<<endl;
    cout<< "cross check: "<<entropy(T_low,t3_low) << endl;

    Type T_high = 1478.456;
    cout<<"____________________________________"<<endl;
    cout << "HIGH temp: T_high = "<< T_high <<endl;
    cout<<"____________________________________"<<endl;
    cout << "heat_capa_low = "<< nasa.heat_capacity(spec,T_high)<<endl;
    cout<< "cross check: "<<heat_capacity(T_high,t3_high) <<endl;

    cout << "enthalpy_low = "<<nasa.enthalpy(spec,T_high)<<endl;
    cout<< "cross check: "<<enthalpy(T_high,t3_high) <<endl;

    cout << "entropy_low = "<<nasa.entropy(spec,T_high)<<endl;
    cout<< "cross check: "<<entropy(T_high,t3_high) << endl;
    
    
    cout << endl<< "2.) STOICHIOMETRY:"<<endl;
    typedef StoichiometricMatrix<DataType::index_type*> StoichType;
    StoichType SM(DataType::nspec,DataType::nreac,DataType::stoichiometric_matrix());

    cout << "nu' and nu'' (in this order); species are columns --> are even numbered rows, <-- odd ones"<<endl;
     for(size_t i = 0; i < 2*DataType::nreac; ++i){
       for(size_t k = 0; k < DataType::nspec; ++k){
	 cout << SM.nu_forward(k,i) << "  "<< SM.nu_reverse(k,i)<< "  ";
       }
       cout << endl;
     }
      

     cout << endl<< "3.) FORWARDRATES:"<<endl;
     TroeIndex TIndex;
     TIndex.create(DataType::nreac,DataType::troewtb());

     typedef DataType::value_type d_type;
     typedef DataType::int_type i_type;
     typedef ForwardReactionRates<const d_type*,const d_type*, const i_type*, const d_type*,'a'> FwdType;
     FwdType FRR(DataType::nspec,DataType::nreac,DataType::arrhenius_values(),DataType::temperature_bounds(), DataType::ntroereac, DataType::troewtb(),DataType::collision_efficiencies(),&TIndex);


     size_t reac = 7,
       r1 = 9;
     Type temperature = 1403.65;
     cout << "k_f,"<<reac<<" = "<< FRR.k_forward(reac,temperature) << endl;
     cout << "cross check: "<< arrhenius_formula(temperature,1.126000e+07,   -7.600000e-01,   0.000000e+00) <<endl;
     cout << "k_f,"<<r1<<" = "<< FRR.k_forward(r1,temperature) << endl;
     cout << "cross check: "<< arrhenius_formula(temperature,2.650000e+10,   -6.707000e-01,   7.129950e+04)<<endl;


     typedef ReverseReactionRates<FwdType,StoichType,NasaPolyType> RevType;
     
     
     
     RevType RRR(&FRR,&SM,&nasa);

    cout << "k_r,"<<endl << reac<< " = "<< RRR.k_reverse(spec,temperature) << endl;

    
    size_t tindex = 5;
    cout <<"alpha_{k," << tindex <<"} = " <<endl;
    for(size_t k = 0; k < DataType::nspec; ++k)
      cout<< DataType::collision_efficiencies()[tindex*DataType::nspec+k]<<endl;

    
    cout << "Test A_low,beta_low, Ea_low:"<<endl;
    size_t tix = 0;
    
    cout<< "A_low = "<< DataType::arrhenius_values()[FRR.index_A_low(tix)]<<"   beta_low = "<<  DataType::arrhenius_values()[FRR.index_beta_low(tix)] << "   Ea_low = "<< DataType::arrhenius_values()[FRR.index_Ea_low(tix)] << endl;
    cout << "alpha = "<< DataType::arrhenius_values()[FRR.index_troe_alpha(tix)]<<"   T*** = "<< DataType::arrhenius_values()[FRR.index_troe_t3star(tix)] << "   T* = "<< DataType::arrhenius_values()[FRR.index_troe_t1star(tix)] << "   T** = "<<  DataType::arrhenius_values()[FRR.index_troe_t2star(tix)]<<endl;



    typedef double ADType; 
    cout << "specific moles to concentrations:"<<endl;
     ADType st[] = {0.0,7.0702437,0.0,0.0,14.140487,0.0,0.0,0.0,26.584116,//spec
		   1500 };   //temp
    
    MyVec<ADType> Start(st,st+DataType::nspec+1);
    cout << "specific moles = "<< Start<<endl;
    typedef Quantity<'s'> QType; 
    MyVec<ADType> Concent;
    
    QType::resize(DataType::nspec+1,Concent);
    Concent[DataType::nspec] = Start[DataType::nspec];

    
    QType::concentration(DataType::nspec,Concent,Start,DataType::molar_masses(),DataType::density_of_mixture());
    cout << "Concent = "<< Concent<<endl;

    cout << "Test forward_rates of troe reaction:"<<endl;
    size_t troe_reac = 20;
   

    cout << "k_f,troe = "<< scientific << FRR.k_forward(troe_reac,Concent[DataType::nspec],Concent) <<endl;
    
    cout << "Forward Rates at T = "<< Concent[DataType::nspec]<<endl<<"____________________________________________"<<endl;
    for(size_t i = 0; i < DataType::nreac; ++i){
      cout << setprecision(4)<<FRR.k_forward(i,Concent[DataType::nspec],Concent)<< "  ";
    }
    cout<<endl;

    cout << endl<< "Reverse Rates at T = "<<  Concent[DataType::nspec]<<endl<<"____________________________________________"<<endl;
    for(size_t i = 0; i < DataType::nreac; ++i){
      cout << setprecision(4)<<RRR.k_reverse(i,Concent[DataType::nspec],Concent)<< "  ";
    }
    cout<<endl;

    
    //needed for a Matlab(r) application
    //sum alpha_ki*[X_k]
    size_t ntb = 6;
    MyVec<ADType> M(ntb);
    
    for(size_t n = 0; n < ntb; ++n){
      for(size_t k = 0; k < DataType::nspec; ++k){
        M[n] += DataType::collision_efficiencies()[n*DataType::nspec + k]*Concent[k];
      }
    }

    cout << endl<< "mass = "<< M <<endl;
    
    cout << endl<< "EVALUATE RHS:"<<endl;
    typedef H2Gri<Type> FunType;
    FunType Fun(DataType::nspec+1);
    Type start[] = {0.0,7.0702437,0.0,0.0,14.140487,0.0,0.0,0.0,26.584116,//spec
		   1500 };   //temp
    
    MyVec<Type> Init(start,start+DataType::nspec+1);
    cout << "Input in specific moles:"<<endl;
    cout << "Start = "<< Start << endl; 
    
    cout<< "Fun(Start) = "<< setprecision(15) << Fun(Init) <<endl;
    
    cout << endl << "Cross check with another value:"<<endl;
    Type oth[] = {6.5241105350736750e-01,
	     1.5037422689287090e+00,
	     3.2983414409901290e+00,
	     4.6981446488606887e-01,
	     2.2521264284480416e+00,
	     2.2121065214370252e-03,
	     3.1765653531269466e-03,
	     1.0000000000000016e+01,
	     2.6584116000000002e+01,
	     1.5000000000000000e+03};

    MyVec<Type> Other(oth,oth+DataType::nspec+1);

    cout << "Fun(Other) = "<< Fun(Other) <<endl;

    Type yet[] = {6.0760354976002506e-01, 
		  9.2729150282555606e-01,
		  1.4679590700772216e+00,
		  1.7371731637159893e+00,
		  2.5972674429795659e+00,
		  3.0739933389068818e-04,
		  1.2741748038786796e-05,
		  9.9999999987085442e+00,
		  2.6585499999999996e+01 ,
		  2.9054713921909702e+03};

    MyVec<Type> YetAnother(yet,yet+DataType::nspec+1);

    cout << "Fun(YetAnother) = "<< Fun(YetAnother) <<endl;


    cout<<endl<< "CONVERT SANDIA TRANSPORT DATA 2 SI UNITS"<<endl;
    Type transportdata[] = {
      //O
      0,    80.000,     2.750,     0.000,     0.000,     0.000,
      //O2
      1,   107.400,     3.458,     0.000,     1.600,     3.800,
      //H
      0,   145.000,     2.050,     0.000,     0.000,     0.000,
      //OH                
      1,    80.000,     2.750,     0.000,     0.000,     0.000,
      //H2                 
      1,    38.000,     2.920,     0.000,     0.790,   280.000,
      //HO2                
      2,   107.400,     3.458,     0.000,     0.000,     1.000,
      //H2O2               
      2,   107.400,     3.458,     0.000,     0.000,     3.800,
      //H2O                
      2,   572.400,     2.605,     1.844,     0.000,     4.000,
      //N2                 
      1,    97.530,     3.621,     0.000,     1.760,     4.000
    };

   
    typedef TransformTransportDataIntoRightUnits<Type,const Type*,'I'> TrafoTransType; 
    typedef TrafoTransType::ContainerType TTTContainer;

    TTTContainer v = TrafoTransType::instance(DataType::nspec,transportdata).convert();

    /* //works
    //H2  H  O2  O H2O  OH 
    Type stoih2c6[] = {
      1,0,0,0,0,0,
      0,2,0,0,0,0,
      0,0,1,0,0,0,
      0,0,0,2,0,0,
      0,0,0,0,1,0,
      0,1,0,0,0,1,
      1,0,0,1,0,0,
      0,1,0,0,0,1,
      0,1,1,0,0,0,
      0,0,0,1,0,1,
      1,0,0,1,0,0,
      0,0,0,0,1,0
     };
    
    StoichiometricMatrix<const Type*> Stoi(6,6,stoih2c6);
    cout << Stoi << endl;

    cout << endl<< "Test index injector:"<<"_________________________"<<endl;
    size_t k = 1,
      j = 2;
    
    CommonIndexInjector DCII; //default
    
    CommonIndexInjector CII(6,8);

    cout<<"common SPEC index at k = "<< k <<": "<< CII.s_ix(k) <<endl;
    cout<<"common REAC index at j = "<< j <<": "<< CII.r_ix(j) <<endl;

    cout<< "Test copy functionalities:"<<endl;
    cout << "copy constructor:"<<endl;
    CommonIndexInjector CC(CII);
    cout<< "CC -- species = " << CC.number_of_species() <<endl;
    cout<< "CC -- reactions =" << CC.number_of_forward_reactions() <<endl;
    
    cout << "copy assignment:"<<endl;
    DCII = CII;
    cout<< "DCII -- species = " << DCII.number_of_species() <<endl;
    cout<< "DCII -- reactions =" << DCII.number_of_forward_reactions() <<endl;


    cout<<endl<<"reduced indexer:"<<endl;
    MyVec<size_t> Rs(2);
    Rs[0] = 0; Rs[1] = 4;
    
    size_t Rr[] = {0,2,3,5};
 
    ReducedIndexInjector<MyVec<size_t>::const_iterator,const size_t*> RDef;

    ReducedIndexInjector<MyVec<size_t>::const_iterator,const size_t*> Red(2,4,Rs.begin(),Rr);

    cout<<"Reduced SPEC index at k = "<< k <<": "<< Red.s_ix(k)<<endl;
    cout<<"Reduced REAC index at j = "<< j <<": "<< Red.r_ix(j) <<endl;

    cout << "copy constructor:"<<endl;
    cout<<endl<<"Test copy functionality:"<<endl;
    ReducedIndexInjector<MyVec<size_t>::const_iterator,const size_t*>  CR(Red);
    cout<< "CR -- species = " << CR.number_of_species() <<endl;
    cout<< "CR -- reactions = " << CR.number_of_forward_reactions() <<endl;
    
    cout << "copy assignment:"<<endl;
    RDef = Red;
    cout<< "RDef -- species = " << RDef.number_of_species() <<endl;
    cout<< "RDef -- reactions =" << RDef.number_of_forward_reactions() <<endl;
    cout<<"Reduced SPEC index at k = "<< k <<": "<< RDef.s_ix(k)<<endl;
    cout<<"Reduced REAC index at j = "<< j <<": "<< RDef.r_ix(j) <<endl;
*/

   
    

    /* //works  
    cout<< "even (-->)/odd (<--) reactions"<<endl;
    for(size_t i = 0; i < 56; ++i){
      cout << i << ".)   "<< ((i%2 == 0) ? true : false) <<endl;
    }
*/
    /* //works
      cout<< "Test dot product" <<endl;
    const size_t dim = 3;

    MyVec<double> v(dim), w(dim);

    v[0] = 1.0038197461437246;
    v[1] = 0.5923111987435578;
    v[2] = -0.4152183943948576;

    w[0] = -0.5321635143211294;
    w[1] = 3.0021735453454532;
    w[2] = 2.7382466454787999;

    cout << "dot(v,w) = "<< setprecision(16) << dot(v,w) <<endl;

    int ex = -1;
    double b = 3.25;

    complex<double> zn(3.,4.);

    cout << endl<< "natural_power(R) = "<< natural_power(b,ex) <<endl;
    cout << endl<< "natural_power(C) = "<< natural_power(zn,ex) <<endl;
    */
   
  }

  else if ( example == 12){
  
    double val = -5.e-24;

    if(val < 0.)
      cout << "val = "<< val << " is NEGATIVE" << endl;
    else 
      cout << "val = "<< val << " is positive" << endl;


    MyVec<double> Yfrac(10);

    Yfrac <<= 1.9558687605162e-14,   0.35529993880904,   1.6848276600126e-19,   7.0664824696754e-19,   0.022400000121411,   1.691498162069e-12,   7.0566137620673e-17,   4.786745991972e-16,   0.62230004752831,   -1.2549002078591e-23;

    cout << "sum of entries = "<< UnrollLoop<0,10>::sum<1>(Yfrac.begin()) << endl;

    MyVec<double> Yfrac_corr(Yfrac);
    Yfrac[9] = 0.;
    
    cout << "sum of entries = "<< UnrollLoop<0,10>::sum<1>(Yfrac_corr.begin()) << endl;
    cout << setprecision(24) << "Y_AR = "<< 1- (Yfrac_corr[0] + Yfrac_corr[1] + Yfrac_corr[2] + Yfrac_corr[3] + Yfrac_corr[4] + Yfrac_corr[5] + Yfrac_corr[6] + Yfrac_corr[7] +  Yfrac_corr[8]) << endl;
    

  }  
  else if ( example == 13){
    cout<<" EXpression examples -- only explicit methods"<<endl;
    double xnu[] = {4., 1.25}; //initial
    double att = 1.5;

    MyVec<double> Xnu(xnu,xnu+2);
     
     ODE<double,KincaidCheney> ivp1(0.,att,Xnu);
     cout<<"NUMERICAL solution of 'KincaidCheney' = "<<
       ivp1.rk4(5000)          //CLASSICAL RK4
       //ivp1.rk3(5000)            //3-stage RK (a.k.a. RK3)
	 <<endl;

     
     //exact solution
     KincaidCheney<double> Ex1(2);
     double a = 2., b = -1.;
     cout<<" EXACT solution ="<<Ex1.exact_solution(att,a, b)<<endl;
   
     cout << endl;

     which_addition_is_used();
     
  }

  else if (example == 14){
   
    ParameterData PD;
    PD.read_from_file("datafiles/ODEsolverSettings.dat");
    double t0 = PD.get_datum<double>("t0");
    double tend = PD.get_datum<double>("tend");
    unsigned maxsteps = PD.get_datum<unsigned>("maxsteps");
    unsigned maxNewt = PD.get_datum<unsigned>("maxNewt");
    double ntol = PD.get_datum<double>("ntol");
    double atol = PD.get_datum<double>("atol");
    double hEstim = PD.get_datum<double>("hEstim");
    double Cscal = PD.get_datum<double>("Cscal");
    double hmin = PD.get_datum<double>("hmin");
    double hmax = PD.get_datum<double>("hmax");

    

    //initial value vector
    MyVec<double> ystart(4);
    double spec[3] = {0.,0.8,0.2}; // Braack  //{0.,0.75, 0.25};  //Warnatz  
    //
    

    ystart[0] = spec[0];            //O
    ystart[1] = spec[1];          //O2
    ystart[2] = spec[2];          //O3
    ystart[3] = PD.get_datum<double>("T0"); //1500.;         //T
    
    ////!EXPLICIT
    // unsigned Msteps = PD.get_datum<unsigned>("Msteps");
    // double hEst =  PD.get_datum<double>("hEst");
    // double TOL = PD.get_datum<double>("TOL");

    typedef ODE<double,O3Decomposition> DiffEqType;
    //typedef DiffEqType::IndexVecType IndexVecType;
  
    DiffEqType ivp(t0,tend,ystart);
   
    //IndexVecType sv(2);
    //sv[0] = 1; sv[1] = 3; 
    MyVec<double> Sol = ivp.implicit_euler(maxsteps,maxNewt,ntol,atol,"IE_O3.dat",hEstim,Cscal,hmin,hmax //,sv
);
    //ivp.rkf45(Msteps,hEst,TOL,"RKF45_O3.dat");

    bool makePicture = PD.get_datum<bool>("toeps");
    string mp;
    
    if(makePicture == true){
      mp = "set terminal postscript eps color enhanced; set output 'O3_ODE.eps'";
      mp +=  ";";
    }
   
    string liwi = my_function_collection::num2str(PD.get_datum<double>("lw"));

    GNUplot graphixs;
    graphixs(mp + "scaleme(t) = 1.e-3*t;set key box right; set ytics nomirror; set y2tics; set tics out; set autoscale y; set autoscale y2; set y2label \"Thermal energy in kK\"; set mxtics 2; set mytics 2; set grid xtics ytics mxtics mytics; set xlabel \"time (sec)\"; set ylabel \"Species mass fractions \"; set title \"O3 Decomposition\";plot \"IE_O3.dat\" using 1:2 with lines lw + " + liwi + " title 'O', \"IE_O3.dat\" using 1:3 with lines lw " + liwi + " title 'O2',\"IE_O3.dat\" using 1:4 with lines lw " + liwi +" title 'O3', \"IE_O3.dat\" using 1:(scaleme($5)) with lines lw "+ liwi +" axes x2y2 title 'T'");
  

    cout << "SOLUTION: "<< Sol << endl;

    // GNUplot smol;
    // smol("splot 'IE_O3.dat' using 2:3:4 title \"phase space\"");

     // cout<<"Test nice gnuplot stuff..."<<endl;
    // GNUplot yeah;
    // yeah.plot_from_file("fancy.gplot");
    
  }
  else if (example == 15){
 
    cout << "Conaire H2 mechanism:"<<endl;
    ParameterData PD;
    PD.read_from_file("datafiles/conH2.dat");
    double t0 = PD.get_datum<double>("t0");
    double tend = PD.get_datum<double>("tend");
    unsigned maxsteps = PD.get_datum<unsigned>("maxsteps");
    unsigned maxNewt = PD.get_datum<unsigned>("maxNewt");
    double ntol = PD.get_datum<double>("ntol");
    double atol = PD.get_datum<double>("atol");
    double hEstim = PD.get_datum<double>("hEstim");
    double Cscal = PD.get_datum<double>("Cscal");
    double hmin = PD.get_datum<double>("hmin");
    double hmax = PD.get_datum<double>("hmax");

    

    //initial value vector
    MyVec<double> ystart(11);
    ystart <<=  PD.get_datum<double>("O"), 
      PD.get_datum<double>("O2"),
      PD.get_datum<double>("H"),
      PD.get_datum<double>("OH"),
      PD.get_datum<double>("H2"),
      PD.get_datum<double>("HO2"),
      PD.get_datum<double>("H2O2"),
      PD.get_datum<double>("H2O"),
      PD.get_datum<double>("N2"),
      PD.get_datum<double>("AR"),
      PD.get_datum<double>("T0");  //temperature
    
    ////!EXPLICIT
    // unsigned Msteps = PD.get_datum<unsigned>("Msteps");
    // double hEst =  PD.get_datum<double>("hEst");
    // double TOL = PD.get_datum<double>("TOL");

    const bool SPARSE = true;
    const bool SCALEDLINSYS = true;
    typedef ODE<double,ConaireH2,0,SPARSE,SCALEDLINSYS> DiffEqType;
    //typedef DiffEqType::IndexVecType IndexVecType;
  
    DiffEqType ivp(t0,tend,ystart);
   
   
    // ConaireH2<double> fun(11); //works
    // cout << "fun(ystart) = "<< fun(ystart) << endl;

    
    MyVec<double> Sol = ivp.implicit_euler(maxsteps,maxNewt,ntol,atol,"H2_conaire.dat",hEstim,Cscal,hmin,hmax //,sv
);
  
    bool makePicture = PD.get_datum<bool>("toeps");
    string mp;
    
    if(makePicture == true){
      mp = "set terminal postscript eps color enhanced; set output 'H2_conaire.eps'";
      mp +=  ";";
    }
   
    string liwi = my_function_collection::num2str(PD.get_datum<double>("lw"));

    string fac = PD.get_datum<string>("scalesomespecies");

    string formatOfTics = " set format x \"%.1e\";";

    GNUplot graphixs;
    graphixs(mp + "scaleme(t) = "+fac+"*t; set autoscale xy; set mxtics 2; set mytics 2;"+formatOfTics+" set grid xtics ytics mxtics mytics; set xlabel \"time (sec)\"; set ylabel \"Species mass fractions \"; set title \"H2/O2 Conaire mech\";plot \"H2_conaire.dat\" using 1:3 with lines lw 2 title \"O2\",  \"H2_conaire.dat\" using 1:(scaleme($5)) with lines lw 2 title \"OH (x "+fac+")\",  \"H2_conaire.dat\" using 1:6 with lines lw 2 title \"H2\", \"H2_conaire.dat\" using 1:9 with lines lw 2 title \"H2O\"");
  

    //plot for [i=1:10] \"H2_conaire.dat\"  using 1:i with lines lw + " + liwi +"

    GNUplot TemperDisp;
    TemperDisp(formatOfTics+"set mxtics 2; set mytics 2; set autoscale xy; set grid xtics ytics mxtics mytics; set xlabel \"time (sec)\"; set ylabel \"Temperature in K \"; set title \"H2/O2 Conaire mech\";plot \"H2_conaire.dat\" using 1:12 with lines lw 2 title \"T\"");

    cout << "SOLUTION of CONAIRE H2 mech: "<< Sol << endl;



   // cout << "Test FDiagonal"<<endl;
    // FieldDiagonal<double,4> D1(4.2);
    
    // FieldDiagonal<double,5> D2;
    // D2 = 3;

    // MyVec<double> M(3);
    // M[0] = 1; M[1] = 3; M[2] = 5;

    // FieldDiagonal<double,3> D3(M.begin(), M.end());

    // cout<< D1<<endl;
    // cout << D2 <<endl;
    // cout << D3 << endl;
  
    // cout << "Test copy and assignment"<<endl;
    // FieldDiagonal<double,3> D4(D3);
    // cout<<D4<<endl;
    // FieldDiagonal<double,3> D5 = D3;
    // cout << D5<<endl;
  
    // int idx = 2;
    // cout << D3[idx] <<endl<< D3[1]<<endl;
    // D3 *= -1.5;
    // cout << D3<<endl;
    // D3 /= 2.;
    // D3 /= -0.0005;
    // cout<<D3<<endl;

    // typedef FieldDiagonal<double,3>::iterator iter;
    // typedef FieldDiagonal<double,3>::const_iterator c_iter;
 
    // for(iter it = D4.begin(); it != D4.end(); ++it){
    //   *it *= 2.;
    //   cout << *it << " ";
    // }
    // cout << endl;
    // for(c_iter it = D5.begin(); it != D5.end(); ++it){
    //   cout << *it << " ";
    // }
    // cout <<endl;
  

    // cout <<endl<<"det = "<<D5.determinant()<<endl;
    
    // D5.invert();
    // cout << "Inverse = "<<D5 <<endl;
  
    // cout<< "Is singular = "<<D5.is_singular()<<endl;

    // FieldDiagonal<double,3> Id;
    // Id.identity();
    // cout<< Id<<endl;

    // FieldDiagonal<double,3> Hi;
    // Hi[0] = 0.5; Hi[1] = -2.5; Hi[2] = 1.75;

    // Dune::FieldMatrix<double,3,3> W;
    // W[0][0] = 3;  W[0][1] = 5;  W[0][2] = 1;
    // W[1][0] = 2;  W[1][1] = 4;  W[1][2] = 5;
    // W[2][0] = 1;  W[2][1] = 2;  W[2][2] = 2;

    // cout << W <<endl;

    // cout << "Test operations:"<<endl;
    // cout<<"1) Addition"<<endl;
    // cout << Hi +W <<endl;
    // cout<< W + Hi <<endl;

    // cout <<"2) Subtraction"<<endl;
    // cout << Hi- W <<endl;
    // cout<< W - Hi <<endl;


    // cout << "Multiplication:"<<endl;
    // cout << Hi*W <<endl;
    // cout<< W*Hi<<endl;
  
    // cout<<"Diag vec multiplication and vice versa"<<endl;
    // MyVec<double> v1(3),v2(3);
    // v1[0] = 1; v1[1] = 2; v1[2] = 3;
    // v2[0] = 1; v2[1] = 2; v2[2] = 3;

    // cout<< Hi*v1<<endl;
    // cout<<v2*Hi<<endl;
  }
  
  else if (example == 16) {
    //H2 GRI 3.0 part for H2; exactly same settings as in example 16! (Only difference: here, no Ar is used, i.e. we have only 9 chem species)
    ParameterData PD;
    PD.read_from_file("datafiles/conH2.dat");
    double t0 = PD.get_datum<double>("t0");
    double tend = PD.get_datum<double>("tend");
    unsigned maxsteps = PD.get_datum<unsigned>("maxsteps");
    unsigned maxNewt = PD.get_datum<unsigned>("maxNewt");
    double ntol = PD.get_datum<double>("ntol");
    double atol = PD.get_datum<double>("atol");
    double hEstim = PD.get_datum<double>("hEstim");
    double Cscal = PD.get_datum<double>("Cscal");
    double hmin = PD.get_datum<double>("hmin");
    double hmax = PD.get_datum<double>("hmax");

    

    //initial value vector
    MyVec<double> ystart(10);
    ystart <<=  PD.get_datum<double>("O"), 
      PD.get_datum<double>("O2"),
      PD.get_datum<double>("H"),
      PD.get_datum<double>("OH"),
      PD.get_datum<double>("H2"),
      PD.get_datum<double>("HO2"),
      PD.get_datum<double>("H2O2"),
      PD.get_datum<double>("H2O"),
      PD.get_datum<double>("N2"),
      PD.get_datum<double>("T0");  //temperature
    
    ////!EXPLICIT
    // unsigned Msteps = PD.get_datum<unsigned>("Msteps");
    // double hEst =  PD.get_datum<double>("hEst");
    // double TOL = PD.get_datum<double>("TOL");

    const bool SPARSE = true;
    const bool SCALEDLINSYS = true;
    typedef ODE<double,H2Gri,0,SPARSE,SCALEDLINSYS> DiffEqType;
    //typedef DiffEqType::IndexVecType IndexVecType;
  
    DiffEqType ivp(t0,tend,ystart);
   
  
    
    MyVec<double> Sol = ivp.implicit_euler(maxsteps,maxNewt,ntol,atol,"H2_Gri.dat",hEstim,Cscal,hmin,hmax //,sv
);
  
    bool makePicture = PD.get_datum<bool>("toeps");
    string mp;
    
    if(makePicture == true){
      mp = "set terminal postscript eps color enhanced; set output 'H2_Gri.eps'";
      mp +=  ";";
    }
   
    string liwi = my_function_collection::num2str(PD.get_datum<double>("lw"));

    string fac = PD.get_datum<string>("scalesomespecies");

    GNUplot graphixs;
    graphixs(mp + "scaleme(t) = "+fac+"*t; set mxtics 2; set mytics 2; set grid xtics ytics mxtics mytics; set xlabel \"time (sec)\"; set ylabel \"Species mass fractions \"; set title \"H2 Gri 3.0 mech\";plot \"H2_Gri.dat\" using 1:3 with lines lw 2 title \"O2\",  \"H2_Gri.dat\" using 1:(scaleme($5)) with lines lw 2 title \"OH (x "+fac+")\",  \"H2_Gri.dat\" using 1:6 with lines lw 2 title \"H2\", \"H2_Gri.dat\" using 1:9 with lines lw 2 title \"H2O\"");
  

    //plot for [i=1:10] \"H2_conaire.dat\"  using 1:i with lines lw + " + liwi +"

    GNUplot TemperDisp;
    TemperDisp("set mxtics 2; set mytics 2; set grid xtics ytics mxtics mytics; set xlabel \"time (sec)\"; set ylabel \"Temperature in K \"; set title \"H2 Gri 3.0 mech\";plot \"H2_Gri.dat\" using 1:11 with lines lw 2 title \"T\"");

    cout << "SOLUTION: "<< Sol << endl;

   

    // map<string, string> Map;

    // Map["a"] = "hallo";
    // Map["b"] = "world";
    // Map["c"] = "i'm here...";

    // cout << Map.find("b")->second << endl;
    // //cout << Map.find("x")->second << endl; //not there, gives SEG.FAULT


    // cout << "Test something -- numerical critical:"<<endl;
    // double mH2O = 0.0180158;
    // cout << "molar mass [kg/mol] of H2O = "<< mH2O <<endl;
    // cout << "One molecule weighs [kg] = "<< mH2O/PhysicalConstants<double>::AvogadrosConstant <<endl;


    // //  cout <<endl<< "Prefactor of SELF diffusion [Curtiss, Hirschfelder] = "<< 3./8.*std::sqrt(M_PI*PhysicalConstants<double>::AvogadrosConstant*ntimes<3>(PhysicalConstants<double>::kB))/M_PI << endl; 


    // double Pi = 3.141592653589793238462643383279;

    // cout << "SELF diff--cross check = "<< 0.375/(Pi)*sqrt(1./2.*Pi*6.02214*ntimes<3>(1.38065))*1.e-23/1.e-15 << endl;
    
    // cout << "BIN. Diff = "<< 3./16.*sqrt(Pi*ntimes<3>(1.38065) * 6.02214)/(Pi) *1.e-8<< endl;
   
  }
  
  else if(example == 17){
    cout << "CHEMKIN default H2 mechanism"<<endl;
    ParameterData PD;
    PD.read_from_file("datafiles/conH2.dat");
    double t0 = PD.get_datum<double>("t0");
    double tend = PD.get_datum<double>("tend");
    unsigned maxsteps = PD.get_datum<unsigned>("maxsteps");
    unsigned maxNewt = PD.get_datum<unsigned>("maxNewt");
    double ntol = PD.get_datum<double>("ntol");
    double atol = PD.get_datum<double>("atol");
    double hEstim = PD.get_datum<double>("hEstim");
    double Cscal = PD.get_datum<double>("Cscal");
    double hmin = PD.get_datum<double>("hmin");
    double hmax = PD.get_datum<double>("hmax");

    

    //initial value vector
    MyVec<double> ystart(10);
    ystart <<=  PD.get_datum<double>("O"), 
      PD.get_datum<double>("O2"),
      PD.get_datum<double>("H"),
      PD.get_datum<double>("OH"),
      PD.get_datum<double>("H2"),
      PD.get_datum<double>("HO2"),
      PD.get_datum<double>("H2O2"),
      PD.get_datum<double>("H2O"),
      PD.get_datum<double>("N2"),
      PD.get_datum<double>("T0");  //temperature
    
    cout << "T0 = "<< PD.get_datum<double>("T0") << endl;

    // ChemkinH2Default<double> fun(10);
    // cout << "fun(ystart) = "<< fun(ystart) << endl;

    const bool SPARSE = true;
    const bool SCALEDLINSYS = true;
    typedef ODE<double,ChemkinH2Default,0,SPARSE,SCALEDLINSYS> DiffEqType;
    //typedef DiffEqType::IndexVecType IndexVecType;
  
    DiffEqType ivp(t0,tend,ystart);
   
  
    
    MyVec<double> Sol = ivp.implicit_euler(maxsteps,maxNewt,ntol,atol,"Chemkin_H2.dat",hEstim,Cscal,hmin,hmax //,sv
);
  
    bool makePicture = PD.get_datum<bool>("toeps");
    string mp;
    
    if(makePicture == true){
      mp = "set terminal postscript eps color enhanced; set output 'Chemkin_H2.eps'";
      mp +=  ";";
    }
   
    string liwi = my_function_collection::num2str(PD.get_datum<double>("lw"));

    string fac = PD.get_datum<string>("scalesomespecies");

    GNUplot graphixs;
    graphixs(mp + "scaleme(t) = "+fac+"*t; set mxtics 2; set mytics 2; set grid xtics ytics mxtics mytics; set xlabel \"time (sec)\"; set ylabel \"Species mass fractions \"; set title \"Chemkin H2 default mechanism \";plot \"Chemkin_H2.dat\" using 1:3 with lines lw 2 title \"O2\",  \"Chemkin_H2.dat\" using 1:(scaleme($5)) with lines lw 2 title \"OH (x "+fac+")\",  \"Chemkin_H2.dat\" using 1:6 with lines lw 2 title \"H2\", \"Chemkin_H2.dat\" using 1:9 with lines lw 2 title \"H2O\"");
  

    

    GNUplot TemperDisp;
    TemperDisp("set mxtics 2; set mytics 2; set grid xtics ytics mxtics mytics; set xlabel \"time (sec)\"; set ylabel \"Temperature in K \"; set title \"Chemkin H2 default mechanism\";plot \"Chemkin_H2.dat\" using 1:11 with lines lw 2 title \"T\"");

    cout << "SOLUTION: "<< Sol << endl;
    
   
    

    // cout << "AKZO-NOBEL PROBLEM"<<endl;
    // // initial and time horizon taken from "Eriksson et al., "Explicit Time-stepping for stiff ODEs""
    // double iv[6] = {0.437, 0.00123, 0., 0., 0., 0.367};
    // MyVec<double> u0(iv,iv+6);
    
    // double t0 = 0.,
    //   tend = 180.;
    
    // double time1, time2;

    // ODE<double,AkzoNobelProblem> ode(t0,tend,u0);
    
    
    // CPUTime<double> cpu;
    
    // string fname = "ImplEul_AkzoNoble.dat";
    // cpu.start();

    
    // MyVec<double> Simpl = ode.implicit_euler(10000, 7, 1.e-04, 1.e-03, fname, 0.00001);
    
    // time1 = cpu.stop();
    // cout<<"IMPLICIT EULER =  "<<  Simpl << endl;
   
    // cpu.reset();
    
    // //string latexstuff = "set terminal latex; set output \"AkzoNobel.tex\"; set format xy \"$%g$\"; ";

    // GNUplot plotter;
    // plotter("set title \"Impl. Euler\"; set yrange[0:0.425]; set mxtics 2; set mytics 2; set grid xtics ytics mxtics mytics; set xlabel \"t \"; set ylabel \" U(t) \"; plot for[i=2:" + num2str(Simpl.size()) + "] \"" +fname + "\" using 1:i with lines lw 2 notitle\"" );

    
   
    //================ tols for weird ====================
    
    // double TOL = 1.e-6,
    //   tol = 1.e-6,
    //   hEstim =  0.0001, //0.0000125,
    //   St = 0.35,
    //   c = 0.98;  //close to 1 

    // double kmin = 0.000005,
    //   kmax = 1.;
    
    // size_t maxit = 5;
    // //====================================================

    // string e4sname = "AkzoNobel_withMethodofErikssonetal.dat";

    // cpu.start();

    // MyVec<double> Result = ode.explicit_4_stiff(TOL,tol,hEstim,St,c,maxit,e4sname,kmin,kmax);

    // time2 = cpu.stop();

    // string ltx = "\n set terminal epslatex \n set output \"akzoNobel_expl_for_stiff.tex\" \n  replot";
    // //for using gnuplot and latex, see ~/books/Graphics/Gnuplot....pdf, p. 214
    // GNUplot grafix;
    // grafix("set title \"Expl. 4 stiff ODEs\"; set yrange[0:0.425]; set mxtics 2; set mytics 2; set grid xtics ytics mxtics mytics; set xlabel \"$t$\"; set ylabel \"$U(t)$\"; plot for[i=2:" + num2str(Result.size()) + "] \"" +e4sname + "\" using 1:i with lines lw 2 notitle\"" + ltx);


    // cout<<"My Implicit Euler     =  "<< Simpl << endl;
    // cout<<"ERIKKSON-JOHNSON-LOGG =  "<< Result << endl;
    
    // MyVec<double> Diff(6);
    // Diff = Simpl - Result;

    // cout << endl << "||IE - EX4Stiff|| = " << Norm<'2',double>::norm(Diff) << endl;
    // cout << "speed up: "<<time_ratio(time1,time2) << endl;
    
  }
  else if (example == 18){
   
    cout << "1D scalar heat equation:"<<endl;
    unsigned numberOfNodes = Constant<double>::discretizationPointsInSpace+2; 

    MyVec<double> u0(numberOfNodes);  //u(x,0) = 0, x \in [0,1]
    
    //u0 = 0.;
    
    ODE<double,HeatEquation1D> ivp(0.,1.,u0);

    CPUTime<double> cpu;
    cpu.start();

    //==============
    double TOL = 1.e-04,
      tol = 1.e-04,
      k1 = 0.000124,
      S_T = 0.45,
      c = 0.98;

    size_t maxit = 3;
    
    double kmin = 1.e-06,
      kmax = 12.35;
    //==============

    MyVec<double> Sol = ivp.explicit_4_stiff(TOL,tol,k1,S_T,c,maxit,"Heat_equation_1D.dat",kmin,kmax);
      //ivp.implicit_euler(30000, 7, 1.e-05, 1.e-04,"IE_Scalar1DHeat.dat",0.000124);
    cpu.stop();
      
  }

  else if (example == 19){
#if USE_CPPAD
    // // //WORKS as desired
    // typedef CppAD::AD<double> Type; 
    // SBEQProblem1<Type> prob;
    
    // MyVec<double> redu(2);
    // redu[0] = 0.345;
    // redu[1] = 0.885;
    

    // prob.set_additional_bounds(redu);

    // Type xf[] = {0.1, 0.2, 0.05, 0.3, 0.05, 0.3};
    // MyVec<Type> x(xf,xf+6),
    //   eq;

    // eq = prob.equality_constraints(x);

    // cout << "eq = "<< eq << endl;

    // ==================================================================
    
    cout << "Test globally convergent algorithm..."<<endl;

    //! fm \f$ \in \f$ box constraints
    //double fm[] = {0.1, 1.25, 0.68, 0.75, 0.98, 0.125};
    //except constant temperature
    double x_future[] ={//H2
      2.9999999999999999e-01,
      //H
      1.8946083116906379e-01,
      //O2
      1.4367564782438827e-01,
      //O
      1.0210953552028719e-01,
      //H2O
      5.9999999999999998e-01,
      //OH
      1.0539168830936280e-02};

    
    // double x_past[] = {// H2
    //   3.0004541410314495e-01,
    //   //H
    //   1.8944296646551151e-01,
    //   //O2
    //   1.4354979042038563e-01,
    //   //O
    //   1.0238618155887168e-01,
    //   //H2O
    //   5.9995196772784154e-01,
    //   //OH
    //   1.0562269872515486e-02
    // };

    const char mustrategy = 'n'; //'n' = no mu update
    typedef GloballyConvergentAugmentedLagrangianMethod<double,SBEQProblem1,mustrategy> AugLagMethType;
    const int eqdim = AugLagMethType::neq;
    double ld[] = {1.234, 0.55, -0.03, 0.75};
    MyVec<double> x(x_future,x_future+6),   //starting point
      lambda(ld,ld+eqdim),
      Sdiag(eqdim,1.),  //scales the constraints, default = 1.
      Mdiag(6,1.); //scales the gradient, default = 1,
    
    Sdiag[0] = 20.3; Sdiag[1] = 3.04; Sdiag[2] = 1.75; Sdiag[2] = 2.08;
    
    Mdiag[0] = 0.45; Mdiag[1] = 0.035; Mdiag[2] = 0.09; Mdiag[3] = 0.35;
    Mdiag[4] = 0.995; Mdiag[5] = 0.035;

    MyVec<double> reduced(2);
    reduced[0] = 0.39;   //0.345;
    reduced[1] = 0.61;   //0.885;
    
    
    CPUTime<double> cpu;
    
    AugLagMethType GCALM(x,lambda,Sdiag,Mdiag,"../optimization/Data4auglagrmethod.dat",reduced);

    cout << "#eq = "<< AugLagMethType::neq << endl;

    cpu.start();

    MyVec<double> xStar = GCALM.optimum();
    cout << "Alledged solution found after "<< GCALM.number_of_iterations() << " iterations with OUTRAGEOUS value:"<< endl;
    cout << "x* ~" << xStar  << std::endl;

    cout << endl<< "check if point fulfills equality constraints"<<endl << "Should be approx. ZERO...:"<<endl;
    cout << GCALM.equality_constraints(xStar) << endl;
    
    cpu.stop();
#endif
  }
  else if (example == 20){ 
#if USE_CPPAD   
 //================= TODO: choose your example here, pal ===================
    int prob = 3;  //1 = lin. obj., nonlin. constraint, 2 = convex obj., lin constraint, 3 = MIT, convex. obj., lin. constraints (1 slacked)
    //=========================================================================

    if(prob == 1){
      cout << " min{x_1 + x_2 | x_1² + x_2² -2}" << endl
	   << "--------------------------------"<<endl;
    
    
      MyVec<double> x(2), lambda(1),
	Sdiag(1),  //scales constraints
	Mdiag(2);  //scales gradient
      x[0] = 1.67;
      x[1] = 0.65;
      
      lambda[0] = -0.4;
      //everything 1, i.e. no scaling
      Sdiag[0] = 1;
      Mdiag[0] = 1;  Mdiag[1] = 1;  
    

    
    
      GloballyConvergentAugmentedLagrangianMethod<double,EasyNLPProblem> Algo(x,lambda,Sdiag,Mdiag,"../optimization/Data4auglagrmethod.dat");

      MyVec<double> xStar = Algo.optimum();
      cout << "===================================================================="<<endl;
      cout << "Alledged solution found after "<< Algo.number_of_iterations() << " iterations with OUTRAGEOUS value:"<< endl;
      cout << "x* ~"  << xStar  << std::endl;
      
      cout << endl<< "check if point fulfills equality constraints"<<endl << "Should be approx. ZERO...:"<<endl;
      cout <<Algo.equality_constraints(xStar) << endl;
   
    }
    
    if(prob == 2){
      cout << " min{6x1² + 4x1x2 + 3x2² | x_1 + x_2 -5}" << endl
	   << "--------------------------------"<<endl;
      
    
      MyVec<double> x(2), lambda(1),
	Sdiag(1),  //scales constraints
	Mdiag(2);  //scales gradient
      x[0] = -0.5;
      x[1] = 0.35;
      
      lambda[0] = 0.;
      //everything 1, i.e. no scaling
      Sdiag[0] = 1;
      Mdiag[0] = 1;  Mdiag[1] = 1;  
    

    
    
      GloballyConvergentAugmentedLagrangianMethod<double,SnymanNLPProblem> Algo(x,lambda,Sdiag,Mdiag,"../optimization/Data4auglagrmethod.dat");

      MyVec<double> xStar = Algo.optimum();
      cout << "===================================================================="<<endl;
      cout << "Alledged solution found after "<< Algo.number_of_iterations() << " iterations with OUTRAGEOUS value:"<< endl;
      cout << "x* ~"  << xStar  << std::endl;

      cout << endl<< "check if point fulfills equality constraints"<<endl << "Should be approx. ZERO...:"<<endl;
      cout <<Algo.equality_constraints(xStar) << endl;
   
    }

    
     if(prob == 3){
       cout << " min{2x1² + x1x2 + 2x2² -6x1 - 6x2| -2x1 + 2x2 + 1 = 0, "<< endl << "x1+ 2x2 -5 + slack = 0," << endl << "x1<= 1.75, x2 <= 2, x3 >= 0} "
	   << "--------------------------------"<<endl;
      
    
      MyVec<double> x(3), lambda(2),
	Sdiag(2,1.),  //scales constraints
	Mdiag(3,1.);  //scales gradient
      x[0] = 0.5;
      x[1] = 0.;
      x[2] = 0.015; //SLACK >= 0
      

      lambda[0] = 0.; //0.0023;
      lambda[1] = 0.; //-0.04;
    

    
    
      GloballyConvergentAugmentedLagrangianMethod<double,MITEasyNLPProblem> Algo(x,lambda,Sdiag,Mdiag,"../optimization/Data4auglagrmethod.dat");

      MyVec<double> xStar = Algo.optimum();
      cout << "===================================================================="<<endl;
      cout << "Alledged solution found after "<< Algo.number_of_iterations() << " iterations with OUTRAGEOUS value:"<< endl;
      cout << "x* ~"  << xStar  << std::endl;

      cout << endl<< "check if point fulfills equality constraints"<<endl << "Should be approx. ZERO...:"<<endl;
      cout <<Algo.equality_constraints(xStar) << endl;
   
    }

    /* //works 
    cout << "matrix_vector_multiplication:"<<endl;   
    double a[4][3] = {{0.5, -0.75, 1.2},
		      {-1.55, 0.65, -2.15},
		      {3.5, -1.35, 2.25},
		      {-1.45, -5.15, 6.65}};
    Dune::FieldMatrix<double,4,3> A;
    for(int i = 0; i < 4; ++i)
      for(int j = 0; j < 3; ++j)
	A[i][j] = a[i][j];
    


    std::vector<double> w(3);
    w[0] = -1; w[1] = 0.45; w[2] = 1.85;

    cout << "Amat = "<< endl << A <<endl;
    cout << "w = "; print_all(w);
    

    Dune::FieldVector<double,4> res;
    
    matrix_vector_product(res,A,w);
    
    cout << "res = "<<res <<endl;

    cout<<"Another example"<< endl;
    double b[2][6] = { {1, 0.5, -4.5, 3., -2., 1.5},
		       {0.5, 5, 2.5, -1, -2, 3.5}};
    Dune::FieldMatrix<double,2,6> B;
    for(int i = 0; i < 2; ++i)
      for(int j = 0; j < 6; ++j)
	B[i][j] = b[i][j];
    
    Dune::FieldVector<double,6> z;

    //MyVec<double> z(6);   //werks
    z[0] = 9; z[1] = 2; z[2] = -3; z[3] = 4; z[4] = -2; z[5] = 6;
   
     MyVec<double> res2(2);

     matrix_vector_product(res2,B,z);
    
    cout << "res2 = "<<res2 <<endl;

    cout << "Test vector-vector-operation:"<<endl;
    
    MyVec<double> u1(3), v1(3), w1(3);
      
      for(unsigned i = 0; i < 3; ++i){
        u1[i] = (i+1);
        v1[i] = 0.5*(i+1);
        w1[i] = std::pow((i+1)/2.,2.);
      }
      
   
      vector_vector_operation<SubtractBasicElements<double> >(u1,v1,w1);

      std::cout << "u1 = " << u1 << std::endl;
    
      vector_vector_addition(u1,v1,w1);

      cout<<"alternative; u1 = "<<u1 <<endl;
      */

      /*
      cout << "Test:"<<endl;
      std::vector<double> v(3);
     // v = 2.75; ERR: no match for ‘operator=’ in ‘v = 2.75e+0’
      */

    /*
    MyVec<double> test(6);
    for(int i= 0; i< 6; ++i)
      test[i] = (i+1);

    cout << "test·z = "  << // DotProduct<6,double>::result(test,z) <<endl; 
	 dot_product<6>(test,z) <<endl;

    //works
   typedef ValueTypeWrapper<IsADuneFieldContainer<Dune::FieldVector<double,6> >::Value, Dune::FieldVector<double,6> >::value_type value_type;

   value_type d = 1.765;
   
   cout << "d = "<< d << "  type: "<<typeid(value_type).name() <<endl;

    //works
    cout << "IsADuneFieldContainer ? : "<<IsADuneFieldContainer<Dune::FieldMatrix<double,5,13> >::Value <<endl;
    cout << "IsADuneFieldContainer ? : "<<IsADuneFieldContainer<std::vector<double> >::Value <<endl;

    */
#endif
  }

  
  else if (example == 21){
#if USE_CPPAD 
   cout << "Example form LANCELOT_SIMPLE"<<endl;

    MyVec<double> x(3), lambda(2),
      Sdiag(2),  //scales constraints
      Mdiag(3);  //scales gradient
    
    //Start value for 
    x[0] = -1.2;
    x[1] = 1.0;

    double slack = 0.01;  //1.56; //0.135;
    x[2] = slack;

    //starting values for lambda
    lambda[0] = 0.5; 
    lambda[1] = -0.75;
    //everything 1, i.e. no scaling
    Sdiag[0] = 1; Sdiag[1] = 1.;
    Mdiag[0] = 1;  Mdiag[1] = 1; Mdiag[2] = 1.;  
    
    
    GloballyConvergentAugmentedLagrangianMethod<double,LancSimpleNLPProblem> Algo(x,lambda,Sdiag,Mdiag,"../optimization/Data4auglagrmethod.dat");

    MyVec<double> xStar = Algo.optimum();
    cout << "===================================================================="<<endl;
    cout << "Alledged solution found after "<< Algo.number_of_iterations() << " iterations with OUTRAGEOUS value:"<< endl;
    cout << "x* ~"  << xStar  << std::endl;

    cout << endl<< "check if point fulfills equality constraints"<<endl << "Should be approx. ZERO...:"<<endl;
    cout <<Algo.equality_constraints(xStar) << endl;
   

   #endif
    
  }

  else if(example == 22){
#if USE_CPPAD
    cout << "Minimize rosenbrock_function(x):"<< std::endl
	 << "--------------------------------"<<std::endl;

     MyVec<double> x(2), lambda,
      Sdiag,  //scales constraints
      Mdiag(2);  //scales gradient
    x[0] = -1.2;
    x[1] = 1.0;

    //everything 1, i.e. no scaling
    Mdiag[0] = 1;  Mdiag[1] = 1;  
    
    
    GloballyConvergentAugmentedLagrangianMethod<double,RosenNLPProblem> Algo(x,lambda,Sdiag,Mdiag,"../optimization/Data4auglagrmethod.dat");

    MyVec<double> xStar = Algo.optimum();
    cout << "===================================================================="<<endl;
    cout << "Alledged solution found after "<< Algo.number_of_iterations() << " iterations with OUTRAGEOUS value:"<< endl;
    cout << "x* ~"  << xStar  << std::endl;

    cout << endl<< "check if point fulfills equality constraints"<<endl << "Should be approx. ZERO...:"<<endl;
    cout <<Algo.equality_constraints(xStar) << endl;
   #endif
  }
  
  else if(example == 23){
#if USE_CPPAD 
    //////////// TODO -- choose your example //////////////
    int problem = 2;  //2 = Hock-Schittkowsky no. 71 (Ipopt)
    //////////////////////////////////////////////////////

    if(problem == 1){
    cout << "min{exp(x1·x2·x3·x4·x5) - 0.5*(x1³ + x2³ +1)²| "<< std::endl
	 << "    x1²+ x2² + x3² + x4² + x5² -10 = 0 "<< std::endl
	 << "    x2·x3 - 5x4·x5 = 0"<< std::endl
	 << "    x1³ + x2³ +1 = 0}"<<std::endl;
      

     MyVec<double> x(5), lambda(3),
      Sdiag(3),  //scales constraints
      Mdiag(5);  //scales gradient
    
     //!given on page 562
    x[0] = -1.71; 
    x[1] = 1.59; 
    x[2] = 1.82; 
    x[3] = -0.763; 
    x[4] = -0.763;
    

    lambda[0] = 1.75; 
    lambda[1] = -2.35;
    lambda[2] = 4.5;

    //everything 1, i.e. no scaling
    Sdiag[0] = Sdiag[1] = Sdiag[2] = 1.;
    Mdiag[0] = Mdiag[1] = Mdiag[2] = Mdiag[3] = Mdiag[4] = 1;

    
    
    
    GloballyConvergentAugmentedLagrangianMethod<double,TestNLPProblem> Algo(x,lambda,Sdiag,Mdiag,"../optimization/Data4auglagrmethod.dat");

    MyVec<double> xStar = Algo.optimum();
    cout << "===================================================================="<<endl;
    cout << "Alledged solution found after "<< Algo.number_of_iterations() << " iterations with OUTRAGEOUS value:"<< endl;
    cout << "x* ~"  << xStar  << std::endl;

    cout << endl<< "check if point fulfills equality constraints"<<endl << "Should be approx. ZERO...:"<<endl;
    cout <<Algo.equality_constraints(xStar) << endl;
    }
    
    ////////////////////////////////////////////////////////
    if(problem == 2){

      cout << "min{x1x4(x1 + x2 + x3)+x3| "<< std::endl
	 << "    x1²+ x2² + x3² + x4² - 40 = 0 "<< std::endl
	 << "    -x1·x2·x3·x4 + 25 <= 0"<< std::endl
	 << "     1 <= x1, x2, x3, x4 <= 5}"<<std::endl;
      

      MyVec<double> x(5), lambda(2),
	Sdiag(2,1.),  //scales constraints
	Mdiag(5,1.);  //scales gradient
      
      //!given on page 562
      x[0] = 1.; 
      x[1] = 5.; 
      x[2] = 5.; 
      x[3] = 1.; 
      
      x[4] = 0.025;    //slack
    
      
      lambda[0] = -0.3; 
      lambda[1] = 0.435;

      //everything 1, i.e. no scaling
      Sdiag[0] = 0.1;
      Sdiag[1] = 0.01;
      
      Mdiag[0] = 1.e-01;
      // Mdiag[1] = 1.;
      // Mdiag[2] = 1.;
      // Mdiag[3] = 1.; 
      // Mdiag[4] = 1.;
     
      //! spd matrix -- fixed
      // double cp[] = {1.7392,0.7142,1.8934,1.4695,1.5499,
      // 		     0.5878,0.7385,0.5888,0.6450,
      // 		     2.1655,1.6312,1.6179,
      // 		     1.5914,1.5848,
      // 		     1.9289};

      MyVec<double> Cprec; //(cp,cp+gauss_sum(x.size()));
      
      
      
    
      GloballyConvergentAugmentedLagrangianMethod<double,HockSchittkowsky71> Algo(x,lambda,Sdiag,Mdiag,"../optimization/Data4auglagrmethod.dat");
    
      CPUTime<double> cpu;

      cpu.start();
      MyVec<double> xStar = Algo.optimum(Cprec);
      cout << "===================================================================="<<endl;
      cout << "Alledged solution found after "<< Algo.number_of_iterations() << " iterations with OUTRAGEOUS value:"<< endl;
      cout << "x* ~"  << xStar  << std::endl;
      
      cout << endl<< "check if point fulfills equality constraints"<<endl << "Should be approx. ZERO...:"<<endl;
      cout <<Algo.equality_constraints(xStar) << endl;
      cpu.stop();
    }
   #endif
  }



  else if (example == 24){
    cout << "Test BLOCK MATRIX ALGORITHM:"<<endl;
    
   
    cout << endl<< "SOLVE BLOCK DIAGONAL SYSTEM"<<endl;
    MyVec<double> D1(9),D2(9),D3(9);

    //diag blocks
    double d1[] = { 4,-1, 0,
		    -1, 4,-1,
		    0,-1,4};
    
   
    double d2[] = {1,2,3,
		   0,5,6,
		   4,9,2};

    
    double d3[] = {3,4,5,
		   1,9,1,
		   2,0,1};

    //subdiag blocks
    double l1[] = {1,0,-4,
		   1,3,3,
		   2,7,5};

    double l2[] = {5,1,2,
		   9,1,2,
		   1,4,3};

    //superdiag blocks
    double u1[] = {2,1,9,
		   1,4,-2,
		   1,3,5};

    double u2[] = {1, 3, 4,
		   2,-1, -1,
		   2, -5, 6};

    

    for(int i = 0; i < 9; ++i){
      D1[i] = d1[i];
      D2[i] = d2[i];
      D3[i] = d3[i];
    }
      
    
    //transpose 
    double Dia1[9],
      Dia2[9],
      Dia3[9],
      E1[9],
      E2[9],
      F1[9],
      F2[9];

    //transposed form
    
    MatrixUnroller<0,0,3,3>::transpose(&Dia1[0],&D1[0]);
    MatrixUnroller<0,0,3,3>::transpose(&Dia2[0],&D2[0]);
    //cout << "Transposed stuff: "; print_all(Dia2,Dia2+9);
    MatrixUnroller<0,0,3,3>::transpose(&Dia3[0],&D3[0]);

    MatrixUnroller<0,0,3,3>::transpose(&F1[0],&u1[0]);
    MatrixUnroller<0,0,3,3>::transpose(&F2[0],&u2[0]);
    
    MatrixUnroller<0,0,3,3>::transpose(&E1[0],&l1[0]);
    MatrixUnroller<0,0,3,3>::transpose(&E2[0],&l2[0]);
    


    //right hand side blocks (3 x 1 blocks)
    double b1[] = {1, 1, 2};
    double b2[] = {4, 5, 7};
    double b3[] = {2, 3, 5};

    // ============ This is passed to the block triang solver ==============
    // typedef MyVec<double>::iterator MyVIterType;
    //MyVec<MyVIterType> diag(3);
    //diag[0] = D1.begin(); diag[1] = D2.begin(); diag[2] = D3.begin();
    //double* diag[3];
    MyVec<double*> diag(3);
    //diag[0] = &D1[0]; diag[1] = &D2[0]; diag[2] = &D3[0];
    diag[0] = &Dia1[0]; diag[1] = &Dia2[0]; diag[2] = &Dia3[0]; //tranpose

    double* up[2];
    //up[0] = &u1[0];  up[1] = &u2[0];
    up[0] = &F1[0];  up[1] = &F2[0];//tranpose

    double* low[2];
    //low[0] = &l1[0]; low[1] = &l2[0];
    low[0] = &E1[0]; low[1] = &E2[0];                     //tranpose
    

    vector<double*> r(3);
    r[0] = &b1[0]; r[1] = &b2[0]; r[2] = &b3[0];
    
  

    cout << endl<< "Solution (Golub/van Loan) = " << endl;
    BlockTridiagonalAlgorithm<3,3,1,double,'t'>::solve(&low[0], diag.begin(), &up[0], r.begin());

    cout << "x-block 1 = "; print_all(r[0],r[0]+3); cout<<endl;
    cout << "x-block 2 = "; print_all(r[1],r[1]+3); cout<<endl;
    cout << "x-block 3 = "; print_all(r[2],r[2]+3); cout<<endl;
    


    
    cout << endl<< "Solve via square solve:"<<endl;
    double blkmtx[] = {
      4,    -1 ,    0 ,    2  ,   1 ,    9,     0 ,    0 ,    0,
      -1 ,    4 ,   -1 ,    1 ,    4 ,   -2 ,    0  ,   0  ,   0,
      0,    -1,     4 ,    1,     3,     5,     0 ,    0,     0,
      1 ,    0,    -4,     1 ,    2 ,    3 ,    1 ,    3,     4,
      1 ,    3,     3,     0 ,    5 ,    6  ,   2 ,   -1,    -1,
      2 ,    7 ,    5,     4 ,    9 ,    2 ,    2 ,   -5,     6,
      0 ,    0 ,    0,     5 ,    1 ,    2 ,    3  ,   4,     5,
      0 ,    0 ,    0,     9 ,    1 ,    2 ,    1 ,    9,     1,
      0 ,    0 ,    0 ,    1 ,    4 ,    3 ,    2,     0,     1};
    
    MyVec<double> Blkmtx(blkmtx,blkmtx+81);
    double bv[] = {1,     1,     2,     4,     5 ,    7 ,    2 ,    3 ,    5};
    MyVec<double> Rhs(bv,bv+9);
    
    good_square_solve(Blkmtx,9, Rhs, 1,true);

    cout << "Solution = "<< Rhs << endl;
  }
  
  else if (example == 25){
    cout << "REN-POPE-DAVIS-SKODJE"<<endl;
    cout << "---------------------"<<endl;
    cout << "a) source term" << endl;
    
    ParameterData PD;
    PD.read_from_file("inputfiles/renpopeParams.dat");


    double t0 = 0.,
      tend = PD.get_datum<double>("tend");

    // MyVec<double> z0(2);  //initial condition
    // z0 <<= 0.5, 1.75;
    
    RPDavisSkodje<double> Fun(2);

    //Implicit Euler the handwritten way
    double time(t0),
      step = 0.01;
    size_t maxNewt = 5;
   
    int order = 2;
    MyVec<double> Uprev(2),U(2),b(2),IterMatrix(ntimes<2>(order)),sim(2);
    const double Tol = 1.e-08;

    

    //some initial values 
    double z1_0 = PD.get_datum<double>("z1_0"),   //0.5
      z2_0 = PD.get_datum<double>("z2_0");
    
    double incr = PD.get_datum<double>("incr");

    const size_t ctMax = PD.get_datum<size_t>("niv");  //10
    
    // string fname = "Ren_Pope_Impl_Euler.dat";
    // ofstream of(fname.c_str());
    //of << "#  t   z[0]     z[1]     zM[0]     zM[1]" << endl;
    
    FancyMessages FM;

   
    int rows = int((tend-t0)/step),
      cols = 1+4*int(ctMax);
    
    cout << "rows = "<< rows << "  cols = "<< cols << endl;

    matrix<double> mtx(rows,cols);

    for(size_t p = 0; p < ctMax; ++p){  
      time = t0;  //reset
      
      //rule to generate different initial values for [z1,z2]
      z1_0 += p*incr;
      z2_0 -= p*incr;

      Uprev <<= z1_0, z2_0;
      U = Uprev;


      //terminal output
      FM.nice_output("Initial Value no. "+my_function_collection::num2str(p));
      //Numerical solution: BACKWARD EULER -- easy implementation
      size_t count(0);
      while(time < tend){
    	sim = Fun.exact_SIM(U);  //calculate exact SIM

    	mtx[count][0] = time;  //first is time 
    	mtx[count][1+4*p] = U[0];
    	mtx[count][2+4*p] = U[1];
    	mtx[count][3+4*p] = sim[0];
    	mtx[count][4+4*p] = sim[1];

    	// of << setprecision(14)  <<  time << "  " << U[0] << "   "<< U[1] << "  "
    	//    <<  sim[0] << "   " << sim[1] << endl;

    	cout << count++ << ".)   dt = "<< step << "   time = " << time << endl;
     
    	//Newton iteration
    	U = Uprev;   //start guess
    	for(size_t l = 1; l <= maxNewt; ++l){
    	  b = -(U-Uprev-step*Fun(U));
    	  IterMatrix = -step*Fun.exact_Jacobian(U);
    	  update_diagonal<AddBasicElements>(IterMatrix,order,1.);
    	  good_square_solve(IterMatrix,order,b,1); //b contains solution now
    	  //b = matrix_vector_product(inverse(IterMatrix,order),order,b);
    	  U += b;
    	  if(Norm<'2',double>::norm(b) <= Tol)
    	    break;

    	  if(l == maxNewt)
    	    ADONIS_INFO(Information,"Newton iteration did not converge within "<< maxNewt << " iterations for Tol = "<< Tol << ".");
    	} //end Newton iteraion
    	Uprev = U;    //store current approximation
    	time += step; //increment time 
      }//end time iteration
      

    } //end loop over all initial values
    

    string fina = "variousinitalvalues.dat";
    ofstream ofi(fina.c_str());
    //write matrix to file
    for(int i = 0; i < rows; ++i){
      for(int j = 0; j < cols; ++j){
    	ofi << setprecision(14) << mtx[i][j] << "   ";
      }
      ofi << endl;
    }
    ofi.close();
    
    //!GNUPLOT stuff
    //create appropriate gnuplot command since  plot "…" for [i=1:10] using (2+4*i):(3+4*i) does not to seem to work with my current version of gnuplot
    string plotcomm,
      style = " with lines lw 2 notitle",
      xrange = "set xrange[ "+PD.get_datum<string>("xfrom")+":"+PD.get_datum<string>("xto")+"];",
      yrange = "set yrange[ "+PD.get_datum<string>("yfrom")+":"+PD.get_datum<string>("yto")+"];",
      tolatex,
      texfile = PD.get_datum<string>("TeXfile"),
      end,
      rotate, // = "rotate by 90;";
      plotunit =  PD.get_datum<string>("plotunit"),
      plotsize = "size "+PD.get_datum<string>("width")+plotunit+","+PD.get_datum<string>("height")+plotunit,
      ylab, //= PD.get_datum<string>("ylab"),
      maintitle;//"set title '"+ PD.get_datum<string>("maintitle")+" ($\\varepsilon = \\ $" + PD.get_datum<string>("eps")+ ")'";
    

    if(PD.get_datum<bool>("TeXme")){
      tolatex = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile+".tex'; set format xy \"$%g$\";";
      //tolatex = "set terminal latex; set output \""+texfile+"\"; set format xy \"$%g$\";";
    }

    //!backslash is output as '\\'
    plotcomm = tolatex + xrange + yrange + "set xlabel '$"+PD.get_datum<string>("xlab")+"$'; set ylabel '$"+ylab+"$'"+rotate+maintitle+"; plot \""+fina+"\" using 2:3" + style + ",";
    
    string separator =",";
    for(size_t i = 1; i < ctMax; ++i){
      if(i == ctMax-1)
	separator="";
      plotcomm += (" \""+fina+"\" using "+my_function_collection::num2str(2+4*i)+":"+my_function_collection::num2str(3+4*i)+style+separator);
    }
    plotcomm += (", \""+fina+"\" using 4:5 with points pt 7 title \""+PD.get_datum<string>("addstr")+"\""+end);
    
    // cout << plotcomm << endl;

    GNUplot phsp;
    
    phsp(plotcomm);
    // cout << texfile << endl;
    // cout <<"MOVING *.tex and *.eps TO /home/mfein/INAUGURALDISSERTATION"<<endl;
    // int isys = system(("mv " +texfile + "* /home/mfein/INAUGURALDISSERTATION").c_str());
    // cout << "system status: "<< isys << endl;
   
    //of.close();
    // GNUplot show;
    // show("set xlabel \"time\"; set ylabel \"z(t)\"; plot '"+fname+"' using 1:2 with lines lw 2 title \"z_1\", '"+fname+"' using 1:3 with lines lw 2 title \"z_2\"");
    // GNUplot phasespace;
    // phasespace("set title \"phase space\";set xlabel \"z_1\";set ylabel \"z_2\"; plot '"+fname+"' using 2:3 with lines lw 2 notitle, '"+fname+"' using 4:5 with lines lw 3 title \"exact SIM\"");

  }
  else if (example == 26){
    cout << "Try QR decomposition" <<endl;
    double a[] = {2., 1., 0., 0., 2., 1.,
		  0., 0., 2., 1., 1., 1.};

    MyVec<double> A(a,a+2*6),
      AT = transpose(A,6),      //Fortran style (enters routine)
      ATrans = AT;

    int m = 2,
      n = 6,
      lwork = std::max(1,n),
      lda = std::max(1,m),
      tdim = std::min(m,n),
      info;

    double* tau = new double [tdim];                 
    double* work = new double [lwork];                
    

    F77Type<TypeTraits<double>::Value>::GEQRF(&m,&n,&AT[0],&lda,&tau[0],&work[0], &lwork,&info);
    
    cout << "A = QR: "<< AT << endl;
    
    cout << endl << "Orthonormalize:"<< endl;
    
    // int k = n;
    // F77Type<TypeTraits<double>::Value>::ORGQR(&m, &m, &k, &AT[0], &m, &tau[0], &work[0], &lwork, &info);

    delete[] tau;
    delete[] work;

    //DORGQR generates an M-by-N real matrix Q with orthonormal columns
    cout << "othonormal columns " << AT << endl;

    Dune::FieldMatrix<double,6,2> Fat;

    for(int i = 0; i < n; ++i)
      for(int j = 0; j < m; ++j)
     	Fat[i][j] = ATrans[i*2+j];


    cout << "Fat = "<< Fat << endl;
    cout << orthonormalize(Fat) << endl;

    
    


    // cout<<"A \"finer\" grid for the Ren Pope Davis Skodjie:"<<endl;
    // double iend = 1., istart = 0.;
    // unsigned numofpoints = 200;
    // double step = (iend-istart)/numofpoints;
    // cout<<" I = ["<<istart<<","<<iend<<"],  number of points = "<<numofpoints<<"   step = "<<step<<endl;

    // string name = "grid_"+num2str(numofpoints)+"_points.dat";
    // ofstream file(name.c_str(),ios_base::out);
    // file<<"z_1"<<endl;
    // for(double i = iend; i >= istart; i-= step){
    //   //cout<< i << "  ";
    //   file<<i<<endl;
    // }
    // file<<"eof"<<endl;
    
    // string jochendir = "~/JOCHEN/JOCHEN++/MoRe_examples/RenDS/START";
    // system(("mv "+name+" "+jochendir).c_str())
      ;
  }

  else if(example == 27){
    cout << "Test Choleky on symm. non-pd matrix:"<<endl;
    // double npd[] = {-1,0,2,   ////nearly s. p.d.
    // 		    0,5,3,
    // 		    2,3,-8};

    
    double npd[] = {-0.05, -3, 4, 
		    -3, 100, 6, 
		    4, 6, -0.04};


    MyVec<double> Nonspd(npd,npd+9),
      A(Nonspd);
    int n = 3,
      lda = std::max(1,n),
      info;
    char uplo = 'L';  //U = upper triangle is stored
   
   

    F77Type<'d'>::POTRF(&uplo,&n,&Nonspd[0],&lda,&info);
    if(info == 0)
      cout << "A = L^TL successfully computed :)"<<endl;
    else if (info < 0) 
      cout << info << "th entry is illegal"<<endl;
    else
      cout << "Matrix is NOT positive definite :(" <<endl;

    
    cout << endl<< "Make matrix p.d."<<endl;
    CholeskyWithAddedMultipleOfIdentity<MyVec<double> > CWAMI(A);

    CWAMI.compute(1.e-03,10.,'L',9,true);
    cout << "number of iterations to reach p.d.: "<<CWAMI.number_of_iterations() << endl;
    
    cout << "updated A = "<< CWAMI.updated_matrix() << endl;

  }
  
  else if(example == 28){
    cout<< "···························································"<<endl;
    cout<< "················· NAVIER-STOKES ···························" <<endl;
    cout<< "···························································"<<endl;

    string additionalFile = "examples/automaticallygeneratedsourceterms/2Ddiscr/2Dsettings.dat";
    
    ParameterData PD;
    PD.read_from_file(additionalFile);
    
    size_t Nx = PD.get_datum<size_t>("Nx"),
      Ny = PD.get_datum<size_t>("Ny"),
      numOfPoints = Nx*Ny;
    
    cout << "Nx: "<< Nx << "   Ny: "<< Ny << "    # total points: "<< numOfPoints << endl;


    //O3 mechanism initial values: rho, v1, v2, T, O, O2, O3: 7 variables
    //
    size_t totDim = numOfPoints*7;
    MyVec<double> y0(totDim);

    double v1 = PD.get_datum<double>("v1_in"), 
      v2 = PD.get_datum<double>("v2_in"),
      p0 = PD.get_datum<double>("pconst"),
      a = PD.get_datum<double>("a"),
      b = PD.get_datum<double>("b"),
      c = PD.get_datum<double>("c"),
      d = PD.get_datum<double>("d"),
      temp1 = PD.get_datum<double>("T_min"),
      temp2 = PD.get_datum<double>("T_ignition");
    
    double hx = (b-a)/(Nx-1),
      halfdiam = (d-c)/2.,
      hy = (d-c)/(Ny-1);

    typedef ThermoData4Mechanism<double,3> MechType;
    MyVec<double> Y0(MechType::nspec);
    Y0[0] = 0.;  Y0[1] = 0.8;  Y0[2] = 0.2;

    double Wbar = 0.;

    UnrollLoop<0,MechType::nspec>::mean_molecular_weight_Y(Wbar,Y0.begin(),MechType::molar_masses());

   
    //boundary temperature distribution at TOP and BOTTOM bdy 
    BoundaryTemperatureDistribution<double> BTD;
    BTD.initialize(temp1,temp2,
		   PD.get_datum<double>("hottest_x"));

    // int Nx = PD.get_datum<int>("Nx");
    // double a = PD.get_datum<double>("a"),
    //   b = PD.get_datum<double>("b"),
    //   length = b-a,
    //   hx = length/(Nx-1),
    //   x = a, y = 0.;

    // string fname = "ex28_bdyTemp.dat";
    // ofstream of(fname.c_str(),ios_base::out);
    // for(int i = 0; i < Nx; ++i){
    //   x = i*hx;
    //   of << x << "   " << BTD(x,y)<< endl;
      
    // }
    // of.close();
    

    
    //ofstream swframe("IV.dat",ios_base::out);
    ////! maybe hardly any reaction
    double temper, Tin = PD.get_datum<double>("T_in");
    double x = 0., y = 0., 
      radius;
    for(size_t i = 0; i < Nx; ++i){
      x = i*hx;

      radius = -halfdiam;
      for(size_t j = 0; j < Ny; ++j){
    	y = j*hy;

	if(j == 0 || j == Ny-1){
	  temper = BTD(x,y);
	  //cout << "UP/LOW = "<< temper << endl;
	}
	else
	  temper = Tin;

	y0[i+Nx*j] = p0*Wbar/(PhysicalConstants<double>::IdealGasConstant*temper);  //rho
	y0[i+Nx*j + numOfPoints] = zero(parabolic_inlet_velocity(radius,halfdiam,v1));
	radius += hy;
	y0[i+Nx*j + 2*numOfPoints] = v2;
	y0[i+Nx*j + 3*numOfPoints] = temper;
	y0[i+Nx*j + 4*numOfPoints] = Y0[0];       //O
	y0[i+Nx*j + 5*numOfPoints] = Y0[1];       //O2
	y0[i+Nx*j + 6*numOfPoints] = Y0[2];       //O3

	// swframe << x << "  "<< y << "  "<<y0[i+Nx*j] << "  " << y0[i+Nx*j + numOfPoints] << "  " << y0[i+Nx*j + 2*numOfPoints] << "  " << y0[i+Nx*j + 3*numOfPoints] << "  " << y0[i+Nx*j + 4*numOfPoints] << "  "<< y0[i+Nx*j + 5*numOfPoints] << "  "<< y0[i+Nx*j + 6*numOfPoints] << endl;
      }
      // swframe << endl;
    }
    
    //swframe.close(); abort();

    // double tmprture = PD.get_datum<double>("T_in");
    // for(size_t t = 0; t < numOfPoints; ++t){
    //   //rho  (for dry air at 300 K: ~ 1.161 kg/m³)
    //   y0[t] = p0*Wbar/(PhysicalConstants<double>::IdealGasConstant*tmprture);   
    //   y0[numOfPoints + t] = v1;        //v1    
    //   y0[2*numOfPoints +t] = v2;       //v2
    //   y0[3*numOfPoints +t] = tmprture; //T
    //   y0[4*numOfPoints +t] = Y0[0];       //O
    //   y0[5*numOfPoints +t] = Y0[1];      //O2
    //   y0[6*numOfPoints +t] = Y0[2];      //O3
    // }


    //cout << "y0 = "<< y0 << endl;

    const int TWOD = 2;
    
    const bool SPARSE = true; //true = use sparse lin. alg. (strongly recommended here!!)

    ParameterData Set;
    Set.read_from_file("datafiles/ODEsolverSettings.dat");
    unsigned maxsteps = Set.get_datum<unsigned>("maxsteps"),
      maxNewt = Set.get_datum<unsigned>("maxNewt");
    
    //time horizon
    double t0 =  Set.get_datum<double>("t0"),
      tend =  Set.get_datum<double>("tend"),
      
      ntol = Set.get_datum<double>("ntol"),
      atol = Set.get_datum<double>("atol"),
      hEstim = Set.get_datum<double>("hEstim"),
      Cscal = Set.get_datum<double>("Cscal"),
      hmin = Set.get_datum<double>("hmin"),
      hmax = Set.get_datum<double>("hmax");
   


    //!works
    // cout << "TEST NSE FUNCTOR:"<<endl;
    // typedef CppAD::AD<double> Type;
    // ReactiveNavierStokes2D<Type> RNS(y0.size());
    //cout << "operator(): "  << RNS(y0) << endl;

    typedef  ODE<double,ReactiveNavierStokes2D,TWOD,SPARSE> DiffEqType;
    typedef DiffEqType::IndexVecType IndexVecType;
    
    DiffEqType ode(t0,tend,y0);
    
    cout << "Dimension of ode: "<< ode.dimension() << endl;

   // int sb[] = {0,1,2,3,4,5,6};
    IndexVecType selectVariables; //(sb, sb+7); //empty species selection vector
    
    int printprec = Set.get_datum<int>("printprec");

    CPUTime<double> cpu;
    cpu.start();
    ode.implicit_euler(maxsteps,maxNewt,ntol,atol,"NSE_OZONE_",hEstim,Cscal,hmin,hmax,selectVariables,additionalFile,printprec);
    double elapsedTime = cpu.stop();

    time_in_human_readable_form(elapsedTime);
    
  
    
  }
  else if (example == 29){
    CompressedSparseStorage<double,'C',int> CSCdef;

    int n = 5;
    int nz = 12;

    int Ap[] = { 0, 2, 5, 9, 10, 12 };
    // int Ap[] = {0, 2, 5, 9, 10};
    int Ai[] =    { 0, 1, 0,  2, 4, 1,  2, 3, 4, 2, 1, 4 };
    double Ax[] = { 2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1 };

    MyVec<double> Av(Ax,Ax+12);
    MyVec<int> indexv(Ai,Ai+12), pv(Ap,Ap+n+1); 
    

    CompressedSparseStorage<double,'C',int> CSC1(n,n,nz,indexv,pv,Av);

    cout << "CSC1 = "<< CSC1 << endl;
    
    CompressedSparseStorage<double,'C',int> CSC2(n,n,nz,Ai,Ai+12,Ap,Ap+6,Ax,Ax+12);
    cout << "CSC2 = "<< CSC2 << endl;
    

    cout << "Copy constructor:"<<endl;
    CompressedSparseStorage<double,'C',int> Copy(CSC1);
    cout << "Copy" << Copy << endl;

    cout << "Assignment:"<<endl;
    

    CSCdef = CSC2;
    cout << "CSCddef = " << CSCdef<< endl;

    CSC1 = CSC2;
    cout << "CSC1 =" << CSC1 << endl;

    //cout << "compression scheme:"<< CompressedSparseStorage<double,'C',int>::Value << endl;
    cout << "compression scheme:"<< CSC1.value() <<endl;

    
    cout << "read in from dense row-wise matrix"<<endl;
    double afull[] = {2,   3,   0,   0,   0,
		      3,   0,   4,   0,   6,
		      0,  -1,  -3,   2,   0,
		      0,   0,   1,   0,   0,
		      0,   4,   2,   0,   1};

    MyVec<double> Afull(afull,afull+25);
  

    cout << "RECTANGULAR example:"<<endl;

    ////!works
    // double afu[] = {4,1,0,2,0,
    // 		    0,3,0,0,0,
    // 		    0,0,5,3,6,
    // 		    0,2,1,0,7};
    double afu[] = {4,1,0,
		    0,3,0,
		    0,0,5,
		    0,2,1};
   
    
    //MyVec<double> Afu(afu,afu+20);
    MyVec<double> Afu(afu,afu+12);
    // int r = 4,
    //   c = 5;
    
     int r = 4,
      c = 3;


     //example from SAAD, Ex. 3.7, p.84
     double ma[] = {1,0,0,2,0,
		    3,4,0,5,0,
		    6,0,7,8,9,
		    0,0,10,11,0,
		    0,0,0,0,12};

     MyVec<double> Ma(ma,ma+25);

   
     int ro = 3,
       co = 5;

     double cab[] = {2,0,5,6,0,
		     0,4,0,0,1,
		     3,0,8,0,0};
     
     MyVec<double> Cab(cab,cab+15);

     cout << endl;
     cout << "##########################################################"<<endl;
     cout << "########## TEST READING FROM FULL MTX-ARRAY ##############"<<endl;
     cout << "##########################################################"<<endl;

    CompressedSparseStorage<double,'c',int> CoSp;
    CoSp.fill_from_row_wise_dense(5,5,Afull);
    cout << "CoSp: "<< CoSp << endl;
    
    CompressedSparseStorage<double,'r',int> CoROWSp;
    CoROWSp.fill_from_row_wise_dense(5,5,Afull);
    cout << "CoROWSp: "<< CoROWSp << endl;
    

    cout << "RECTANGULAR example: rows GREATER cols"<<endl;
    CompressedSparseStorage<double,'c',int> CoSp1;
    CoSp1.fill_from_row_wise_dense(r,c,Afu);
    cout << "CoSp1: "<< CoSp1 << endl;
    
    CompressedSparseStorage<double,'r',int> CoROWSp1;
    CoROWSp1.fill_from_row_wise_dense(r,c,Afu);
    cout << "CoROWSp1: "<< CoROWSp1 << endl;
    
    cout << "SAAD, 2nd ed., Ex. 3.7, p. 84:"<<endl;
    CompressedSparseStorage<double,'c',int> CoSp2;
    CoSp2.fill_from_row_wise_dense(5,5,Ma);
    cout << "CoSp2: "<< CoSp2 << endl;
    
    CompressedSparseStorage<double,'R',int> CoROWSp2;
    CoROWSp2.fill_from_row_wise_dense(5,5,Ma);
    cout << "CoROWSp2: "<< CoROWSp2 << endl;
    
    cout << "RECTANGULAR example: rows SMALLER cols:"<<endl;
    CompressedSparseStorage<double,'c',int> CoSp3;
    CoSp3.fill_from_row_wise_dense(ro,co,Cab);
    cout << "CoSp3: "<< CoSp3 << endl;
    
    CompressedSparseStorage<double,'R',int> CoROWSp3;
    CoROWSp3.fill_from_row_wise_dense(ro,co,Cab);
    cout << "CoROWSp3: "<< CoROWSp3 << endl;
  }

  else if (example == 30){
    double a[] = {1.44,  -7.84,  -4.39,   4.53,
		  -9.96,  -0.28,  -3.24,   3.83,
		  -7.55,   3.24,   6.27,  -6.64,
		  8.34,   8.09,   5.28,   2.06,
		  7.08,   2.52,   0.74,  -2.47,
		  -5.45,  -5.70,  -1.19,   4.70};

    double b[] = {8.58,   9.35,
		  8.26,  -4.43,
		  8.48,  -0.70,
		  -5.28,  -0.26,
		  5.72,  -7.36,
		  8.93,  -2.52};

    MyVec<double> A(24), B(12); //B(b,b+12);
    
    //store both matrices in F77 style (i.e. column major)
    unsigned k = 0;
    for(int i= 0; i < 6; ++i)
       for(int j= 0; j < 4; ++j)
	 A[ColumnMajor::offset(i,j,6)] = a[k++];

    
     unsigned l = 0;
    for(int i= 0; i < 6; ++i)
       for(int j= 0; j < 2; ++j)
	 B[ColumnMajor::offset(i,j,6)] = b[l++];
    
    cout<<"B = "<<B <<endl;

    //cout<<"A = "<<A<<endl;
    solve(A,6,B,2,false);
    cout<<"B = "<<B<<endl<<endl;
    cout<<"B = "<< extract_solution(B,6,4,2) <<endl;
    // print_in_matrix_style(B,2,"",2);

    cout<<endl<<"Another example"<<endl;
    double a1[] = {-0.57,  -1.28,  -0.39,   0.25,
		 -1.93,   1.08,  -0.31,  -2.14,
		 2.30,   0.24,   0.40,  -0.35,
		 -1.93,   0.64,  -0.66,   0.08,
		 0.15,   0.30,   0.15,  -2.13,
		 -0.02,   1.03,  -1.43,   0.50};

    double b1[] = {-2.67,-0.55,3.34,-0.77,0.48,4.10   };
    
    MyVec<double> A1(24), B1(b1,b1+6);
    
    unsigned h = 0;
    for(int i= 0; i < 6; ++i)    //store in F77 style
       for(int j= 0; j < 4; ++j)
	 A1[ColumnMajor::offset(i,j,6)] = a1[h++];

    solve(A1,6,B1,1,false);

    cout<<"B1 = "<<B1<<endl;
    cout<<"B1 = "<< extract_solution(B1,6,4,1) <<endl;
    

    cout<<"A square example:"<<endl;
    double a2[] = {1.0000,      0,   -2.0000,
		   3.0000,    0.5000,   -1.0000,
		   0,    2.0000,    2.0000};

    double b2[] = {0.5,-0.5,6};

    MyVec<double> A2(a2,a2+9), B2(b2,b2+3);
    
    solve(A2,3,B2); //same as 'solve(A2,3,B2,1,true)'
    cout<<"B2 = " << B2<<endl;

  }
  else if (example == 31){
    
    cout << "test transpose of a matrix :"<<endl;

    double tot[] = {1,2,4,
		    8,3,6,
		    0,1,5};

    MyVec<double> Mtx(tot,tot+9), TransMtx(9);

    MatrixUnroller<0,0,3,3>::transpose(TransMtx.begin(),Mtx.begin());
    
    print_in_matrix_style(TransMtx,3);
  }


  else if (example == 32){
    cout << "Moore-Penrose pseudoinverse:"<<endl;
    //example taken from http://www.mathematics-online.org/inhalt/beispiel/beispiel565/
    Dune::FieldMatrix<double,4,3> A;
    double a[4][3] = {{2,-4,5},
		      {6,0,3},
		      {2,-4,5},
		      {6,0,3}};
    
    for(int i = 0; i < 4; ++i)
      for(int j = 0; j < 3; ++j)
	A[i][j] = a[i][j];
    
    cout << "A = "<<endl<< A <<endl;
    const int mx = 4;
    Dune::FieldMatrix<double,mx,mx> B;
    for(int i = 0; i < mx; ++i){
      for(int j = 0; j < mx; ++j){
	if (i == j)
	  B[i][j] = 1;
	else 
	  B[i][j] = 0;
      }
    }
    //cout << "B = "<<endl<< B << endl;
    //cout<< solve_rank_deficient_ls(A,B,-1,'c')<<endl;

    cout<< "Moore-Penrose = "<<endl<< pseudo_inverse(A) <<endl;
    
    cout<<endl<<"Next example under-determined (can differ from Matlab) "<< endl;
    
    double a1[2][6] = {{1, 0., 2., 5., 7., -1.},
		       {0., 0.5, -2.35, 1., 2., 0.75}};
    
    Dune::FieldMatrix<double,2,6> A1;
    for(int i = 0; i < 2; ++i)
      for(int j =0; j < 6; ++j)
	A1[i][j] = a1[i][j];
    cout << "A1 = "<<endl<< A1<<endl;
    
    
    cout<< "Moore-Penrose = "<<endl<< pseudo_inverse(A1) <<endl;
    

    cout << "Test special cases formula for pseudo-inverse:"<<endl;
    double v1[] = {1,2,3,4,5,6};
    double v2[] = {2,-1,1,0.5,1.4,-2};
    MyVec<double> V1(v1,v1+6), V2(v2,v2+6);
    normalize(V1);
    normalize(V2);

    cout<< "V1 = "<< V1 << endl<< "V2 = "<< V2<<endl;;
    

    double dp = 0.;
    for(int i = 0; i < 6; ++i)
      dp += V1[i]*V2[i];

    cout<< "dotproduct(v1,v2) = "<<dp <<endl;
    
    Dune::FieldMatrix<double,2,6> NT;
    for(int j = 0; j < 6; ++j){
      NT[0][j] = V1[j];
      NT[1][j] = V2[j];
    }

    cout << "NT = "<< endl<< NT <<endl;

    cout<<"pseudo_inverse_full_row_rank(NT) = "<<endl<<pseudo_inverse_full_row_rank(NT)<<endl<<endl;
    
    
    cout << "Test transpose: "<<endl;
    Dune::FieldMatrix<double,2,3> K;
    unsigned ct = 1;
    for(int i = 0; i < 2; ++i)
      for(int j = 0; j < 3; ++j)
	K[i][j] = ct++;

    cout << "K = "<<endl<< K<<endl;
    
    cout<< "transpose = "<<endl;
    for(int j =0; j < 3; ++j){
      for(int i = 0; i < 2; ++i){
	cout<< K[i][j] << " ";
      }
      cout<<endl;
    }

   
  
    //cout << endl<<transpose(K)*K<<endl<<endl;
    MyVec<double> td(3);
    td[0] = -0.5; td[1] = 1.25; td[2] = 0.75;

    cout << endl << "Test this hunky chunk..."<<endl;
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	cout<< N_i_times_N_T_j(i,j,K)<< " ";
      }
      cout<<endl;
    }

   
    cout << endl<< "Test template meta program DOT PRODUCT:"<<endl;
    cout<< DotProduct<3,double>::result(K[0],td) <<endl;

    cout<<"another dot product example"<<endl;
    Dune::FieldVector<double,6> d1;
    for(int i = 0; i < 6; ++i)
      d1[i] = (i+1)*0.5;

    double df[] = {-1, 0, 0.5, -2., 4, -0.4};
    MyVec<double> d2(df,df+6);
   
    cout << "dot_product = "<<dot_product<6>(d1,d2) <<endl;

    // cout<< "Size = "<<container_size(K[1])<<endl; //works
      
    //cout<< "Size = "<<std::distance(K[1].begin(),K[1].end())<<endl; //works

    

    double ef[2][6] = {{2., -0.5, 1., -1.5, 4., 3.25},
		       {1., 0., -0.75, 2.45, 0.065, -1.25}}; 

    Dune::FieldMatrix<double,2,6> E;
    for(int i = 0; i < 2; ++i)
      for(int j = 0; j < 6; ++j)
	E[i][j] = ef[i][j];

    
    Dune::FieldVector<double,2> g;
    g[0] = 1; g[1] = -3;
    
    cout<< "SOlve rank deficient ls (overloaded function)"<<endl;
    cout << solve_rank_deficient_ls(E,g) <<endl;

    cout<<endl;
    
    }
  else if (example == 33){
     cout << "Test Line search:"<<endl;
     LineSearchWithBackTracking<double,DennisSchnabelExa6_3_1> LS(1.,0.1,0.1);
     //============= TEST ===================
     MyVec<double>  xc(2), grad(2), pk(2); 
     xc[0] = 1.; xc[1] = 1.;
     grad[0] = 6; grad[1] = 2;
     pk = -1.*grad;
     
     
     cout<< "alpha_k = "<< LS.determine_step(xc,grad,pk) <<endl;

     double f2[] = {1,1};
     MyVec<double> x0(f2, f2+2);
     //======================================

     double tolerance = 1.e-6;
     cout<< "Minimum (Dennis Schnabel) = "<<LS.loop(1301,x0,tolerance)<<endl;


     cout << endl<< "Rosenbrock's function:"<<endl;
     
     double st[] = {-3.75, 1.89};  // some shitty start point  //{0,0}
     MyVec<double> Start(st,st+2);
     LineSearchWithBackTracking<double,Rosenbrock> RB(1.,0.1,0.1);
     cout << "Minimum of Rosenbrock's fct = "<< RB.loop(4000,Start,tolerance) <<endl;


     double zt[] = {0,0};
     MyVec<double> B(zt,zt+2);
     LineSearchWithBackTracking<double,FEx2_2> E2(1.,0.1,0.1);
     cout << "Minimum of non p.d. fct. = "<< E2.loop(3000,B,tolerance) << endl;
     

     /*//deprecated worked
       BackTrackingLS<double,DennisSchnabelExa6_3_1> BTLS(1.,0.1,0.1);
     
     DennisSchnabelExa6_3_1<double> fun(1);
     cout<< "alpha_k = "<<BTLS.determine_step(fun,xc,grad,pk) <<endl;
     */
     }
  
  
  else if (example == 34){

   cout<<endl<<"Test dense sym. matrix:"<<endl;
   
   double fis[] = {1,4,6,2,5,3};
   DenseSymmetricMatrix<double> DSM(3), Def,
     ItD(3,fis,fis+6);

   //works
   cout << "ItD = "<< ItD << endl;
   cout << "copy constr.:";
   DenseSymmetricMatrix<double> DSC(ItD);
   cout<< DSC <<endl;
   cout << "copy assign.:";
   Def = ItD;
   cout << Def <<endl;

   cout << "ItD.size() = "<< ItD.size() << endl;

   cout<<endl<<"Test mat-vec multiplication:"<<endl;
   double vf[] = {0.5, -1, 2.75};
   MyVec<double> v(vf, vf+3);

   cout<< ItD*v <<endl;
   //cout<<"...and vice versa:"<<endl<< v*ItD <<endl; //works
   cout << "another example:"<<endl;
   double symbd[] = {0.5, -7, 0.75, 1, 
		     6, 5.2, -4.15,
                     -1.5,3,
                      2};

   DenseSymmetricMatrix<double> Four(4,symbd, symbd + 10);
   MyVec<double> r(4);
   r[0] = 1; r[1] = 0.095; r[2] = -2.345; r[3] = 1.56;

   cout << r*Four <<endl;


   
   cout << endl<< "Do something..."<<endl;
   
   MyVec<double> rhs1(3);
   rhs1[0] = -0.3; rhs1[1] = 1.87; rhs1[2] = -1.54;
   cout<<"rhs1 = "<<rhs1<<endl;
   
   MyVec<double> r1(rhs1);

   ItD.solve(rhs1,1);

   cout << "Solution = "<<rhs1 <<endl;

   

   cout << endl<< "Another example cf. \n      http://www.nag.co.uk/lapack-ex/node19.html"<<endl;
   
   double a1[] = {-1.81, 1.15,-0.21, 2.07,
		  2.06, 1.87, 3.87,
		  0.63, 4.20,
		  -1.15};
   
   DenseSymmetricMatrix<double> A1(4,a1, a1+10);
   double b1[] = {0.96,6.07,8.38,9.50};
   MyVec<double> B1(b1,b1+4);
     
   A1.solve(B1,1);
   cout<< "Solution = "<<B1 <<endl;

   cout<<endl<< "A 3rd example"<<endl;
   double a2[] = {-5.86, 4.46, -8.52, 3.72, 9.83, 
		  3.99, 2.58, 8.57, 8.07,
		  -5.93, 4.42, 7.69,
		  -2.82,4.61, 
		  7.69};
   
   /*double b2[] = {1.32,  -6.33,  -8.77,
		  2.22,   1.69,  -8.33,
		  0.12,  -1.56,   9.54,
		  -6.41,  -9.49,   9.56,
		  6.33,  -3.67,   7.48};*/

   //column-wise storage !!!
   double b2[] = {1.32, 2.22, 0.12, -6.41, 6.33,
                  -6.33, 1.69, -1.56, -9.49, -3.67,
		  -8.77,-8.33, 9.54, 9.56, 7.48};
   
   
   DenseSymmetricMatrix<double> A2(5,a2, a2+15);
   MyVec<double> B2(b2,b2+15);
   
   A2.solve(B2,3,true);
   cout<< "Solution = " << B2 <<endl;
   print_in_matrix_style(B2,3);
   
   /* //works with 'full symmetric matrix'
   double a3[] = {-5.86,   3.99,  -5.93,  -2.82,   7.69,
		  3.99,   4.46,   2.58,   4.42,   4.61,
		  -5.93,   2.58,  -8.52,   8.57,   7.69,
		  -2.82,   4.42,   8.57,   3.72,   8.07,
		  7.69,   4.61,   7.69,   8.07,   9.83};

   MyVec<double> A3(a3,a3+25);
   solve_symmetric_ls(A3,5,B2,3,'L');
   MyVec<double> B2T = transpose(B2,5);
   //cout<<"B2T (Sol.) = "<<B2T<<endl;
   print_in_matrix_style(B2T,3);
   */ 
   /*MyVec<double> K1(5*5);
   A2.fill_me_into_a_full_matrix_given_as_rac<ColumnMajor>(K1);
   print_in_matrix_style(K1,5);

   solve_symmetric_ls(K1,5,B2,3,'U');*/
   
    
   
   bool testDuneVers = false;

   if(testDuneVers){
     cout << "Test solving a symmetric system involving DUNE field types" <<endl;
    double sys[5][5] = {{-5.86, 3.99, -5.93, -2.82, 7.69},
			{3.99, 4.46, 2.58, 4.42, 4.61},
			{-5.93, 2.58, -8.52, 8.57, 7.69},
			{-2.82, 4.42,  8.57, 3.72, 8.07},
			{7.69,  4.61,  7.69, 8.07, 9.83}};

    Dune::FieldMatrix<double,5,5> SymMa;
    for(int i = 0; i < 5; ++i)
      for(int j = 0; j < 5; ++j)
	SymMa[i][j] = sys[i][j];
	
    cout<< "SymMa ="<<endl<<SymMa<<endl;

    double b4s[5][3] = {{1.32, -6.33, -8.77},
			{2.22,  1.69, -8.33},
			{0.12,  -1.56, 9.54},
			{-6.41, -9.49, 9.56},
			{6.33,  -3.67, 7.48}};
    
    Dune::FieldMatrix<double,5,3> B4S, SolSym;
    for(int i = 0; i < 5; ++i)
      for(int j = 0; j < 3; ++j)
	B4S[i][j] = b4s[i][j];
    
    cout<<"B4S = "<<endl<<B4S<<endl;
    
    solve_symmetric_ls(SymMa,B4S);
    cout <<endl<< "Solution of symmetric system = "<<endl<<B4S<<endl;
   }

   
    //Since DSC is not p.d., no cholesky fac. exists
    /*cout<<endl << "DSC = " << DSC <<endl; 
   cout << "r1 = "<<r1 <<endl;
   cout<< "Is pos def ?:"<< DSC.is_positive_definite() <<endl;
   DSC.solve_pos_def(r1,1);
   cout<< "Is pos def (after cholesky)?:"<< DSC.is_positive_definite() <<endl;
    */
   //DSC.update_diagonal<AddBasicElements<double> >(0.021);
   //cout << "DSC (diag updated) = "<< DSC <<endl; //works
   
   


   cout<<endl<<"Test p.d. system, cf. \n     http://software.intel.com/sites/products/documentation/hpc/mkl/lapack/mkl_lapack_examples/dposv.htm"<<endl;
   
   double ddhgs[] = {3.14, 0.79, 4.53, 5.32, 1.98,
		     0.17, 0.83, -3.7, -1.37,
		     -0.9, -0.65, 1.6,
		     1.65, 0.28,
		     -0.72};

   /*double pr[] = {-7.29,  6.11,  0.59,
            9.25,  2.90,  8.88,
            5.99, -5.05,  7.57,
           -1.94, -3.80,  5.57,
           -8.30,  9.66, -1.67};*/

   //column-wise storage
   double pr[] = {7.29,  9.25,  5.99, -1.94, -8.30,
            6.11,  2.90, -5.05, -3.80,  9.66,
            0.59,  8.88,  7.57,  5.57, -1.67};

   MyVec<double> rP(pr, pr+15);
   DenseSymmetricMatrix<double> PD(5, ddhgs, ddhgs+15);
   
   /*MyVec<double> P(5*5);
   PD.fill_me_into_a_full_matrix_given_as_rac<RowMajor>(P);
   solve_positive_definite_ls(P,5,rP,3,'U',false);
   */
   PD.regularization(rP,3,true,0.1,true);
   cout<< "rP = " << rP <<endl;
   

   print_in_matrix_style(rP,3);


   const unsigned ydim = 2;
   DenseSymmetricMatrix<double> Symm(ydim);
   Symm[0] = 6; Symm[1] = 2; Symm[2] = 0;

   cout << "Symm = "<<Symm <<endl;

   MyVec<double> g(ydim);
   g[0] = 6; g[1] = 2;

   Symm.solve(g,1);

   cout << "g = "<<g<<endl;


   cout << "Test compution of symmetric ev's:"<<endl;
   cout << "DSC = "<<DSC <<endl;

    MyVec<double> ev(3);
    DSC.eigenvalues(ev);
    cout << "Eigen values = " << ev <<endl;
    
    cout << DSC.is_positive_definite() << endl;
    //cout << "eigenvalues = "<< DSC.eigenvalues(ev) <<endl; //works
    
    double pascal[] = {1, 2, 6, 20, 70,
		       1, 3, 10, 35,
		       1, 4, 15,
		       1, 5,
		       1};
    

    DenseSymmetricMatrix<double> Pascal(5,pascal,pascal+gauss_sum(5));
    MyVec<double> Ei(5);
    Pascal.eigenvalues(Ei);
    cout << "eig(Pascal) = "<<Ei <<endl;
    cout << "Pos. def. ?: "<<Pascal.is_positive_definite() <<endl;
    

    cout<< "Example taken from \n   http://software.intel.com/sites/products/documentation/hpc/mkl/lapack/mkl_lapack_examples/lapacke_dsyev_row.c.htm"<<endl;
    double itl[] = {1.96, 3.8, 4.17, 5.7, -7.1,
		    -6.49, -6.39, -1.51, 1.8,
		    -0.47, 1.5, 2.67,
		    -7.2, -6.34,
		    -0.65};

    DenseSymmetricMatrix<double> ITL(5,itl,itl+gauss_sum(5));
    
    ITL.eigenvalues(Ei);
    cout<< "Eigenvalues = "<<Ei <<endl;

   

  }
  else if (example == 35){
     //works 
    cout << "Test Projection onto Box Constraints:"<<endl;

    double lw[] = {0.5,  -1.5, -0.75};
    double u[] =  {0.65,  2.,   1.25};
    MyVec<double> low(lw,lw+3), up(u,u+3);
    
    ProjectionOntoBoxConstraints<MyVec<double> > PBXC(low,up);

    MyVec<double> Diag(3),
      proj(3),
      X(3),
      grad(3);
    Diag[0] = 2.0; Diag[1] = 0.5; Diag[2] = 1.5;
    grad[0] = 0.45; grad[1] = 0.15; grad[2] = 0.55;
    X[0] = -7;
    X[1] = 0.76;
    X[2] = 1.4;

    cout << "projection = "<< PBXC.projection(proj,X-Diag*grad) << endl;


    MyVec<double> Y(3);
    Y[0] = -2;
    Y[1] = -0.44;
    Y[2] = 1.56;

    //note, just for testing but not efficient since X-Y is formed X.size()-times ;)
    cout << "Check operator()(i,expr):"<<endl;
    for(int i = 0; i< 3; ++i)
      cout << PBXC(i,X-Y) << endl;

  }

  else if (example == 36){
    typedef MyVec<double> SymmMtxType;
    size_t dim = 2;
    		   
    // double tf[] = {3,1,
    // 		   1,5};
    
    // SymmMtxType C(tf,tf+ntimes<2>(dim));
    SymmMtxType B; 
   
    Rosenbrock<double> rosenbrock(1);
    double xe[] = {0.45, -0.75};
    
    MyVec<double> x(xe,xe+2),
      gradf = rosenbrock.gradient(x);
   
   

     cout << "gradf = "<< gradf<<endl<< "Hessian";
     B = rosenbrock.complete_hessian(x);
     //print_in_matrix_style(B,2);
     double tria[] = {3,1,5};
     SymmMtxType B0, Cspd(tria,tria+gauss_sum(dim));
     full_matrix_2_symmetric_matrix(B0,B);
     
     
     cout << endl<< "STEIHAUG:"<<endl;
     SteihaugPCG<SymmMtxType> stei(dim);

     MyVec<double> direction = stei.approximation(B0,gradf,0.0435, 1.e-04,Cspd); 
     
     cout << "p_approx = "<< direction << endl;

    


     //works
    double cf[] = {3,1,2,
		   1,5,4,
		   2,4,8};
    
    double ve[] = {0.5,-0.25, 0.75};
    MyVec<double> v(ve,ve+3), S1(cf,cf+ntimes<2>(3));
    
    //cout << matrix_vector_product(S1,v) << endl;

    cout << "||v||_C = " << scaled_norm(v,S1) << endl;
   

     // works
    cout << "SIGNUM"<<endl;

    cout << "Sgn(-1.5) = "<<  Sgn(-1.5)<<endl;
    cout << "Sgn(0.0) = "<<  Sgn(0.0)<<endl;
    cout << "Sgn(2.25) = "<<  Sgn(2.25)<<endl;

    cout <<endl<<"complex version:"<<endl;
    complex<double> z0, z1(-3,-4);
    cout << "Sgn(z0) = "<<  Sgn(z0)<<endl;
    cout << "Sgn(z1) = "<<  Sgn(z1)<<endl;

    cout <<endl<<"Test number conversion:"<<endl;
    cout << "real: "<<convert_number(3.75) <<endl;
    CppAD::AD<double> ad(1.5);
    cout << "AD: "<<convert_number(ad) <<endl;
     
    cout << "Heaviside:"<<endl;
    cout << "neg: "<<Heaviside(-1.75) <<endl;
    cout << "pos: "<<Heaviside(0.75) <<endl;
    

    /*
    //works well
    double drei[] = {1,2,3};
    MyVec<double> fucker(drei,drei+3),
      empty;
    
    cout<< "test resize_me;" <<endl;
    resize_me_when_empty(4,fucker,1.75);
    resize_me_when_empty(4,empty,1.75);
  
    cout << "fucker = "<<fucker << endl;
    cout << "empty = "<< empty << endl;
    */
  }

  else if (example == 37){
#if USE_CPPAD == 1   
    typedef CppAD::AD<double> ADType;

    double NaN = 0./0.;
    double a = 0.5,
      b = NaN;

    MyVec<double> eval(2);
    eval <<= a, b;
    MyVec<ADType> X(2), Y(2);
    CppAD::Independent(X);
    Y[0] = X[0]*X[0] + sin(X[1]);
    Y[1] = 3.5*X[1] - X[0]*X[1];

    CppAD::ADFun<double> adfun;
    adfun.Dependent(X,Y);
    cout << "start evaluation with a NaN value " << endl;
    MyVec<double> adJac = adfun.Jacobian(eval);
    print_in_matrix_style(adJac,2,"CppAD Jacobian:",12);
#endif

   
  }

  else if (example == 38){
    

   
    
    
    //works well
    double iv[] = {1,1};
    double t0 = 0,
      tend = 20;
    
    ODE<double, Brusselator> ode(t0,tend,iv,iv+2);
    cout << "Sol (adapt. RKF) = "<<ode.adaptive_rkf(1000,0.025,1e-6) <<endl;
    //cout << "Sol (Implicit Euler) = "<<ode.implicit_euler(1000,7,1e-4,2.5e-02)<<endl;
    
    

    /* //That's rather a PATHOLOCICAL example
    double init[] = {-12,6};
     ODE<double, Lambert> ode(0,10,init,init+2);
     cout << "Num. sol (Lambert) = "<<ode.adaptive_rkf(300,0.01,1e-7) <<endl;

     Lambert<double> Lamb(2);
     cout<<"Exact (Lambert) = "<<Lamb.exact_solution(10) <<endl;
    */

    /*
    double start[] = {3,1};
    ODE<double, KincaidCheney_Ex3> ode(1,-2,start,start+2);
    cout<< "adaptive_rkf (KincaidCheney_Ex3) = " << ode.adaptive_rkf(300, -0.1, 1e-4) << endl;;
    */
    
    
    
  }

  else if (example == 39){
   cout<< "CHECK STIFFNESS OF GIVEN PROBLEM:"<<endl;
   
   //automatically generated source term
   H2MechIn6Species<double> rhs(6);
    double t0 =0.,
      tend = 2.;

    const int stiff = 2;  //1 = computes stiffness due to [?], stiff (> 1.e-08)                          // 2 = due to Ascher/Petzold: stiff ( < -1)   

    StiffnessDetector<double,H2MechIn6Species,false,stiff> IPS(rhs);
    double pkt[] = {0.579301,  0.126241,  0.29082,  0.0505531,  0.34735,  0.0204567};
    MyVec<double> Punkt(pkt,pkt+6);
    cout << "H2C6 stiffness-rate = "<<IPS.check(t0, tend, Punkt) <<endl;
   
    which_derivative_is_used();
  }

  else if (example == 40){

    cout << "Test solve_rank_deficient_ls(A,B):" << endl;
    double a1[6][2] = {{-3.0000,    9.0000},
		       {-1.0000,    3.0000},
		       {0.5000,   -1.5000},
		       {0.2500,   -0.7500},
		       {-2.0000,    6.0000},
		       {-0.1000,    0.3000}};

    Dune::FieldMatrix<double,6,2> A1;
    for(int i = 0; i < 6; ++i)
      for(int j = 0; j < 2; ++j)
	A1[i][j] = a1[i][j];

    cout << "A1 = "<<endl << A1 <<endl;


    double b1[6][3] = {{-1.0000,    2.0000,    3.0000},
		       {0.7000,    0.6000,   -2.0000},
		       {0,    1.0000,   -0.9500},
		       {-2.0000,    4.0000,   -3.0000},
		       {-1.5000,   -2.3000,    0.0560},
		       {1.4500,    5.0000,   -8.0000}};

    Dune::FieldMatrix<double,6,3> B1;
    for(int i = 0; i < 6; ++i)
      for(int j = 0; j < 3; ++j)
	B1[i][j] = b1[i][j];

    cout << "B1 = "<<endl << B1 <<endl;

    Dune::FieldMatrix<double,2,3> Res1 = solve_rank_deficient_ls(A1,B1);

    cout << "Res1 = "<< Res1 << endl;


    cout << endl << "Example from http://software.intel.com/sites/products/documentation/hpc/mkl/lapack/mkl_lapack_examples/lapacke_dgelsd_row.c.htm"<<endl;


    double a2[4][5] ={{0.12, -8.19,   7.69,  -2.26,  -4.71},
		      {-6.91,   2.22,  -5.12,  -9.08,   9.96},
		      {-3.33,  -8.94,  -6.72,  -4.40,  -9.98},
		      {3.97,   3.33,  -2.74,  -7.92,  -3.20}};

    Dune::FieldMatrix<double,4,5> A2;
    for(int i = 0 ; i < 4; ++i)
      for(int j = 0; j < 5; ++j)
	A2[i][j] = a2[i][j];

    
    double b2[4][3] = {{7.30,   0.47,  -6.28},
		       {1.33,   6.58, -3.42},
		       {2.68,  -1.71,   3.46},
		       {-9.62,  -0.79,   0.41}};

    Dune::FieldMatrix<double,4,3> B2;

    for(int i = 0; i < 4; ++i)
      for(int j = 0; j < 3; ++j)
	B2[i][j] = b2[i][j];

    
    Dune::FieldMatrix<double,5,3> Res2;
    
    Res2 = solve_rank_deficient_ls(A2,B2,-1,'c');
    // //!Minimum norm solution as reported:
    // //!  -0.69  -0.24   0.06
    // //!  -0.80  -0.08   0.21
    // //!   0.38   0.12  -0.65
    // //!   0.29  -0.24   0.42
    // //!   0.29   0.35  -0.30

    cout << "Res2 = "<< Res2 << endl;

    
    cout << endl << "  OVERDETERMINED example: Usually an overdetermined system (more equations than unknowns) has no exact solution:" << endl;
    double a[6][4] = {{0.8147,    0.2785,    0.9572,    0.7922},
		      {0.9058,    0.5469,    0.4854,    0.9595},
		      {0.1270,    0.9575,    0.8003,    0.6557},
		      {0.9134,    0.9649,    0.1419,    0.0357},
		      {0.6324,    0.1576,    0.4218,    0.8491},
		      {0.0975,    0.9706,    0.9157,    0.9340}};

    Dune::FieldMatrix<double,6,4> A;
    for(int i = 0; i < 6; ++i)
      for(int j = 0; j < 4; ++j)
	A[i][j] = a[i][j];
    
    // cout << "A = "<< endl << A << endl;

    double b[6][5] = {{0.6787,    0.7060,    0.6948,    0.7655,    0.7094},
		      {0.7577,    0.0318,    0.3171,    0.7952,    0.7547},
		      {0.7431,    0.2769,    0.9502,    0.1869,    0.2760},
		      {0.3922,    0.0462,    0.0344,    0.4898,    0.6797},
		      {0.6555,    0.0971,    0.4387,    0.4456,    0.6551},
		      {0.1712,    0.8235,    0.3816,    0.6463,    0.1626}};

    Dune::FieldMatrix<double,6,5> B;
    for(int i = 0; i < 6; ++i)
      for(int j = 0; j < 5; ++j)
	B[i][j] = b[i][j];
    //cout << "B = "<< endl << B << endl;
    
    Dune::FieldMatrix<double,4,5> X;
    X = solve_rank_deficient_ls(A,B,-1,'c');
    cout << setprecision(5) << "X = "<< endl << X << endl;

    // cout << endl << "TEST: A·X = B?:"<< endl << A*X << endl;
    //cout << "Res = "<< Res << endl;
    

    cout << endl << "Conservation matrix C of H2C6 mech:"<<endl;
    double c[2][6] =  {{2,1,0,0,2,1},
		     {0,0,2,1,1,1}};

    Dune::FieldMatrix<double,2,6> C;
    for(int i = 0; i < 2; ++i)
      for(int j = 0; j < 6; ++j)
	C[i][j] = c[i][j];

    Dune::FieldVector<double,2> brhs;
    brhs[0] = 2.; brhs[1] = 1.;

    Dune::FieldVector<double,6> u = solve_rank_deficient_ls(C,brhs);
    cout << "u = " << u <<endl;

    cout << endl<< "Check: C*u = "<< C*u << endl;


    /* //works as expected
    cout << " TEST NORMAL SPACE" << endl;
    Dune::FieldMatrix<double,6,2> T;

    T[0][0] = 1;        T[0][1] = -3;
    T[1][0] = 0.5000 ;  T[1][1] = 2.4500;
    T[2][0] = 0.7500;   T[2][1] = -0.6450;
    T[3][0] = 3.2000;   T[3][1] = 4.8000;
    T[4][0] = 0.6000;   T[4][1] = 1.;
    T[5][0] = -2.1500;  T[5][1] = 1.7000;

    cout << "Tangent space = "<< endl << T << endl;

    // Dune::FieldMatrix<double,6,6> Q;
    // Q = orthonormalize(T);
    // cout << "Q matrix "<< endl << Q << endl;
  

    Dune::FieldMatrix<double,4,6> NT;
    // NT = extract_normalspace_transposed<2>(Q);
    // cout << "NT = "<< NT << endl;

    // //! cross check
    // Dune::FieldMatrix<double,6,4> N;
    // N = extract_normalspace<2>(Q);
    // cout << "N = "<< N << endl;
    //cout << "NT*T = "<<transpose(N)*T <<endl;
    
    //! all in one 
    NT = normalspace_transposed<2>(T);
    
    cout << "NT = "<< NT << endl;
    cout << "This has to be a zero matrix:"<<endl;
    cout << "NT*T = "<< NT*T << endl;
    */
  }

  else if (example == 41){
    
    // 5 Body Problem
    MyVec<double> Start(30);
    Start[0] = 3.42947415189; Start[1] = -.557160570446;
    Start[2] = 3.35386959711; Start[3] =  .505696783289;
    Start[4] = 1.35494901715; Start[5] =  .230578543901;
    
    Start[6] = 6.64145542550; Start[7] = -.415570776342;
    Start[8] = 5.97156957878; Start[9] =  .365682722812;
    Start[10] = 2.18231499728; Start[11] = .169143213293;
    
    Start[12] = 11.2630437207; Start[13] = -.325325669158;
    Start[14] = 14.6952576794; Start[15] =  .189706021964;
    Start[16] = 6.27960525067; Start[17] =  .0877265322780;
    
    Start[18] = -30.1552268759; Start[19] = -.0240476254170;
    Start[20] =  1.65699966404; Start[21] = -.287659532608;
    Start[22] = 1.43785752721;  Start[23] = -.117219543175;

    Start[24] = -21.1238353380; Start[25] = -.176860753121;
    Start[26] =  28.4465098142; Start[27] = -.216393453025;
    Start[28] =  15.3882659679; Start[29] = -.0148647893090;
    
    
    //===================================
    unsigned maxsteps = 7000;
    double hEstim = 0.075,
      tol = 1.e-07;
    //===================================
    

    ODE<double, FiveBodyProblem> ivp(0., 2.65, Start);
    
    CPUTime<double> cpu;
    cpu.start();
    cout<< "adaptive RKF = " << ivp.adaptive_rkf(maxsteps, hEstim, tol) << endl;
    
    //cout << "implicit Euler = "<<ivp.implicit_euler(maxsteps,7,tol,1e-02,"FiveBodyProblem.dat",hEstim) <<endl; //works
    cpu.stop();
    

    /*
    //works
    double ts = 1.789;
    FiveBodyProblem<double> FBP(30);
    cout<< FBP(ts, Start) <<endl;
    */
  }
  else if (example == 42){
    cout << "TRUST REGION "<<endl
	 << "------------" << endl;
    
    TRRadius<double,'2'> trRadius(0.0005,0.45,0.135);

    double ratio1 = 0.12432,
      ratio2 = 0.88621,
      delta = 0.073;
    double pf[] = {0.5, -0.75, 1.25};
    MyVec<double> p(pf,pf+3);
    
    cout << "delta_{i-1} = " << delta << endl;
    cout << "delta_{i} =   " << trRadius.update(delta,ratio1,p);
    cout << "     nu_ = "<<  trRadius.is_direction_changed() << endl;
    cout << "delta_{i+1} = " << trRadius.update(delta,ratio2,p);
    cout << "     nu_ = "<< trRadius.is_direction_changed() << endl;
    

  }
  else if (example == 43){
    cout << "VARIATIONAL DIFFERENTIAL EQUATION"<<endl;
    cout << "... DAHLQUIST'S TEST Equation y' = ay"<<endl;
    
    int maxit = 500;
    double t0 = 0,
      tend = 1.25;
    
    double kn = Abs(tend-t0)/maxit;

    double u0 = 1., 
      u(u0);
  
    double ascal = 1.745; 
    
   
    double v = 0,  
      w = 1,          //t = 0 then 1
      t = 0;
    //explicit Euler
    for(int i = 0; i < maxit; ++i){
      t = t0 + i*kn;
      cout << "time: "<< t << endl;
      //! derivatives w.r.t. parameter a
      u += kn*ascal*u;     //cout << "u = "<< u << endl;
      v += kn*(ascal*v+u);       //kn*(ascal*v + exp(ascal*t)); //with exact
    
      //! derivatives w.r.t. y
      w += kn*(ascal*w);
      
    }
    cout << "Solution of ivp = "<< u << endl << " analytical solution = "<< exp(ascal*tend) << endl;
    
   

    cout << " Psi (numerical.) = "<< v << endl << " d_a y(t) (anal.) = "<< tend*exp(ascal*tend) << endl;
    
    cout<< " deriv. w.r.t. y = "<< w << "   analyt. = "<< exp(ascal*tend) <<endl;

   
    // //! Works: complex scalar product
    // complex<double> z1(3,4), z2(-0.5,-2), z3(-0.75,0.25),
    //   w1(-4,-2.5),w2(0.15,0.25),w3(1.5,2.15);

    // //! 
    // MyVec<complex<double> > vz(3), wz(3);
    // vz[0] = z1; vz[1] = z2; vz[2] = z3;
    // wz[0] = w1; wz[1] = w2; wz[2] = w3;

    // cout << "vz^T·wz = "<< dot(vz,wz) << endl;

  }


  else if (example == 44){
    cout << "TR ALGORITHM test"<<endl;
    size_t dim = 2;
    Rosenbrock<double> rosenbrock(1);
    double xe[] = {0.45, -0.75};
    //double xe[] = {-1.2,1.0}; //TROUBLESOME STARTING POINT -- working

    MyVec<double> x0(xe,xe+dim);
     

    //============== create AD sequence ======================
    //FD may also be considered here 
#if USE_CPPAD 
   typedef CppAD::AD<double> ADType;
    MyVec<ADType> X(dim),Y(1);
    CppAD::Independent(X);
    
    Rosenbrock<ADType> rbad(1); 
    Y = rbad(X);
    CppAD::ADFun<double> dx;
    dx.Dependent(X,Y);
    dx.optimize();
#else
    FiniteDifference<Rosenbrock<double> > dx(rosenbrock);
#endif

    //========================================================
   
   
    double tria[] = {3,1,5};
    MyVec<double> B0, Cprec(tria,tria+gauss_sum(dim));;
    
    symm_identity(B0,dim);
    
   
    
    typedef SteihaugPCG<MyVec<double> > SolverType;
    TrustRegion<SolverType,'2','b'> TR(0.5, 1., 0.0075, 1.e-04); 
    
    MyVec<double> sol = TR.method(x0,B0,Cprec,rosenbrock,dx);

    cout << "Approximate solution found afer "<< TR.number_of_iterations() << " iterations:"<<endl;
    cout << "sol = "<< sol << endl; 

  }

  else if(example == 45){
   
    ExtendedRosenbrock<double> fun(1);
    size_t n = static_cast<size_t>(ExtendedRosenbrock<double>::dodim);
    cout << "n = "<< n << endl;

    MyVec<double> x0(n);

    // random_container(x0);
    for(size_t i = 0; i < n; ++i){
      x0[i] = -1;  //see NOCEDAL,WRIGHT, p. 191
      // if(i%2 == 0)  //DIXON, MILLS, p. 176
      // 	x0[i] = -1.2;
      // else
      // 	x0[i] = 1.0;
    }

    cout << "x0 = "<< x0 << endl;

    cout << endl<< setprecision(7)<<"fun(x0) = "<< fun(x0) << endl;
    

    //============== create AD sequence ======================
    //FD may also be considered here 
#if USE_CPPAD
    typedef CppAD::AD<double> ADType;
    MyVec<ADType> X(n),Y(1);
    CppAD::Independent(X);
    
    ExtendedRosenbrock<ADType> nocad(1); 
    Y = nocad(X);
    CppAD::ADFun<double> dx;
    dx.Dependent(X,Y);
    dx.optimize();
#else
    FiniteDifference<ExtendedRosenbrock<double> > dx(fun);
#endif  

  
    //========================================================
   
    MyVec<double> cpf(gauss_sum(n));
    random_container(cpf);  //randomly generated symmetric precond. matrix

    MyVec<double> B0, Cprec(cpf);;
    
    symm_identity(B0,n);
    
    typedef SteihaugPCG<MyVec<double> > SolverType;
    TrustRegion<SolverType,'2','b'> TR(0.5, 1., 0.0075, 1.e-02); 
    
    MyVec<double> sol = TR.method(x0,B0,Cprec,fun,dx);

    cout << "Approximate solution found afer "<< TR.number_of_iterations() << " iterations:"<<endl;
    cout << "sol = "<< sol << endl; 

    //int n1 = 9; cout << n1/2 << endl;
  }
  else if (example == 46){
    
    double ft[] = {1,3,7,6,
		   3,2,1,5,
		   7,1,8,4, 
		   6,5,4,9};

    double rf[] = {0.5, -0.25, 1.25, 0.75};
    MyVec<double> A(ft,ft+16), r(rf,rf+4);
    
    

    cout << "SSOR for full mtx = "<< Preconditioner<1,'f'>::prec(A,r) <<endl;
  
    cout <<endl << "SYMMETRIC Case:"<<endl;
    MyVec<double> Sym;
    full_matrix_2_symmetric_matrix(Sym,A);
    cout << "Sym = "<< Sym << endl;

    cout <<endl<< "SSOR for dense symm. = "<< Preconditioner<1,'y'>::prec(Sym,r) <<endl;


    
    double sing[] = {1,2,3,4,6,7};
    MyVec<double> SingSym(sing,sing+6), rhs(3);
    rhs[0] = 0.5; rhs[1] = 0.75; rhs[2] = 1.25;

    cout << endl<< "SSOR for SINGULAR,dense symm. = "<< Preconditioner<1,'y'>::prec(SingSym,rhs) << endl;
  }

  else if (example == 47){
#if USE_CPPAD
    MyVec<double> lambda(1), Sdiag(1);
    lambda[0] = -0.4;
    Sdiag[0] = 0.075;

    MyVec<double> x(2);
    x[0] = 1.67;
    x[1] = 0.65;

    double mu = 0.1;

    DerivativesOfAugmentedLagrangianFunction<double,EasyNLPProblem> DLaug;
    DLaug.set();
    

    cout << "With varying arguments"<<endl;
    MyVec<double> lf(1),ran(2);
    lf[0] = 0.25;
    
    for(int i = 0; i < 4; ++i){
      cout <<"INDEX "<< i<<".)"<<endl;
      cout <<"---------------"<<endl;
      
    
      random_container(ran);

      if(i%2 == 0){
    	lambda += lf;
    	x -= ran; 
      }
      else{
    	lambda -= lf;
    	x += ran;
      }

      mu *= 0.5; //decrease mu

      // cout <<"Perform it a second time:"<<endl;
      // cout << "DLaug= "<< DLaug.gradient(x,lambda,Sdiag,mu) << endl;
      cout << "DLaug = "<< DLaug.gradient(x,lambda,Sdiag,mu) << endl;

      //cout << endl<< "AD derivative:"<<endl;
      AugmentedLagrangianFunction<CppAD::AD<double>,EasyNLPProblem> NablaLaug;
      MyVec<CppAD::AD<double> > X(2),Y(1);
      CppAD::ADFun<double> ads;
      CppAD::Independent(X);
      Y = NablaLaug(X,lambda,Sdiag,mu);
      ads.Dependent(X,Y);
      ads.optimize();
      
      MyVec<double> gL = ads.Jacobian(x), gH = ads.Hessian(x,0);
    
      cout << "gradLaug = "<< gL << endl;

      cout << "D2Laug = "<< DLaug.hessian(x,lambda,Sdiag,mu) << endl;
      cout << "hessLaug = "<< gH << endl;

  } //end for
  #endif
  }
  
  else if(example == 48){
    cout << "check lls solver (A may be rank deficient):"<<endl;
    double a[] = {2,1,0,0,2,1,
		  0,0,2,1,1,1};
    
    MyVec<double> A(a,a+12);

    MyVec<double> b(2);
    b[0] = 2.; b[1] = 1.;

    cout << "X = " << solve_rank_deficient_ls(A,2,b,1) << endl;


    cout << endl<< "2nd example:" << endl;
    double a1[] = {0.12,  -8.19,   7.69,  -2.26,  -4.71,
		   -6.91,   2.22,  -5.12,  -9.08,   9.96,
		   -3.33,  -8.94,  -6.72,  -4.40,  -9.98,
		   3.97,   3.33,  -2.74,  -7.92,  -3.20};

    MyVec<double> A1(a1,a1+20);
    double b1[] = {7.30,   0.47,  -6.28,
		   1.33,   6.58,  -3.42,
		   2.68,  -1.71,   3.46,
		   -9.62,  -0.79,   0.41};
    MyVec<double> B1(b1,b1+12);
    //B1 = transpose(B1,3);

    MyVec<double> Res = solve_rank_deficient_ls(A1,4,B1,3);
    cout << "Res = " <<  Res << endl;

  
    cout << endl << "CHECK: A1*Res = B: "<<endl;
    MyVec<double> Ch(4*3);

    matrix_matrix_multiplication(Ch, A1, 4,Res,3);

    cout << Ch << endl;

    cout << endl<< "B^T d = -B^Tu +r:"<<endl;
    double bt[] = {1,0,0,0,0,0,
		   0,0,0,0,1,0};

    MyVec<double> BT(bt,bt+2*6);

    double xf[] = {0.035,0.01,0.015,0.04, 0.55, 0.35};
    MyVec<double> u(xf,xf+6), r(2), g(2);
    r[0] = 0.61; r[1] = 0.39;
    g = -matrix_vector_product(BT,u)+r;
    cout << "g = "<< g << endl;

    MyVec<double> Ser = solve_rank_deficient_ls(BT,2,g,1);

    cout << "DELTA = "<< Ser << endl;


    //! need to cut out appropriate part
    // cout << endl << "OVERDETERMINED example:"<< endl;
    // double od[] = {1,   2,  0, -1.3,
    // 		   0.5, 1.5, -1,  4,
    // 		   0.2, 0.1, 0.3, 1.7,
    // 		   1.3, 1, -4, -0.4,
    // 		   -0.45, -1.4, 3., 2.13,
    // 		   1., -2, 0.75, 3.1};

    // MyVec<double> OD(od,od+6*4);

    // double bd[] = {0.65, -1.25,
    // 		   9, -2,
    // 		   -3, 4,
    // 		   0.04, -1.24,
    // 		   3.2, 0.75,
    // 		   -0.67, -3.22};
    // MyVec<double> BD(bd,bd+6*2);

    // MyVec<double> ROver = solve_rank_deficient_ls(OD,6,BD,2);

    // cout << "ROver = "<< ROver << endl;

    // MyVec<double> C(12);
    // cout << matrix_matrix_multiplication(C,A1,4,solve_rank_deficient_ls(A1,4,B1,3),3) << endl;

    // double df = -4.35;
    // std::complex<double> c(3.66, df);

    // CppAD::AD<double> ad(3.75);

    // cout << "floor(df) = "<< Floor(df) << endl;
    // cout << "floor(c) = "<< Floor(c) << endl;
    // cout << "floor(ad) = "<< Floor(ad) << endl;
    
     // WORKS: changed MyVec operator= for arbitrary container 
    // double x_future[] ={2.9999999999999999e-01,1.8946083116906379e-01,1.4367564782438827e-01,1.0210953552028719e-01,5.9999999999999998e-01,1.0539168830936280e-02};
    // const int lamdim = 4;
    // double ld[] = {1.234, 0.55, 0.12, 0.75};
    // MyVec<double> x(x_future,x_future+6),   //starting point
    //   lambda(ld,ld+lamdim),
    //   Sdiag(lamdim,1.);  //scales constraints
    
    // MyVec<double> reduced(2);
    // reduced[0] = 0.345;
    // reduced[1] = 0.885;
    
    // double mu = 10.;

    // typedef DerivativesOfAugmentedLagrangianFunction<double,SBEQProblem1> DerivType;
    // DerivType DLaug;
    
    // DLaug.get_problem().set_additional_bounds(reduced);
    
    // DLaug.get_AD_problem().set_additional_bounds(reduced); 

    // cout << "get_red() = "<< DLaug.get_problem().get_red() << endl;
    // cout << "get_red() (2x) = "<< DLaug.get_problem().get_red() << endl;
    
    
    // cout << "get_AD_red() = "<< DLaug.get_AD_problem().get_red() << endl;

    // DLaug.set();

    // // typedef DerivType::ProblemType ProblemType;
    // // //typedef SBEQProblem1<double> ProblemType;
    // // ProblemType fox;
    // // fox.set_additional_bounds(reduced);
    // // cout << fox.equality_constraints(x) << endl;
   

    // // cout << "red = "<< DLaug.get_problem().get_red() <<endl;
    // //cout << "AD red = "<< DLaug.get_AD_problem().get_red() <<endl;

    // cout << "eq(4 entries) = "<<DLaug.get_AD_problem().equality_constraints(x) << endl;

    // MyVec<double> grad = DLaug.gradient(x,lambda,Sdiag,mu);
    // cout << "grad Laug = "<< grad << endl;
  }


  else if (example == 49){
    double wf[] = {2.25, 1.96, -1.69};
    double uf[] = {0.5, -1.3, -4, 2.5, 1.75};
    Dune::FieldVector<double,5> u;
    for(int i= 0; i < 5; ++i)
      u[i] = uf[i];
    MyVec<double> v;
    std::vector<double> w(wf,wf+3);
    
    cout << join(v,w,u) << endl;
    
    //! works
    // cout<< "Test me pal:"<<endl;
     // complex<double> z(2.,4.);

     // cout << "complex arccosine = "<< Acos(z) << endl;
     // cout << "Inverse hyperbolic cosine = "<< Acosh(z) << endl;

     // //works
     // cout << deg2rad(45.) << " rad"<< endl;
     // cout << rad2deg(0.) << "°"<< endl;
     // cout << rad2deg(UniversalConstants<double>::Pi/6.) << "°"<< endl;
     // cout << rad2deg(UniversalConstants<double>::Pi/4.) << "°"<< endl;
     // cout << rad2deg(UniversalConstants<double>::Pi/3.) << "°"<< endl;
     // cout << rad2deg(UniversalConstants<double>::Pi/2.) << "°"<< endl;
     // cout << rad2deg(UniversalConstants<double>::Pi) << "°"<< endl;
     // cout << rad2deg(3*UniversalConstants<double>::Pi/2.) << "°"<< endl;
     // cout << rad2deg(2*UniversalConstants<double>::Pi) << "°"<< endl;
     // //Wiki examples
     // cout << rad2deg(1.) << "°"<< endl;
     // cout << rad2deg(2.5) << " rad"<< endl;
     // cout << deg2rad(1.) << " rad"<< endl;
     // cout << deg2rad(23.) << " rad"<< endl;

     
     //! o.k. throw error since not defined for complex numbers
     //cout << rad2deg(complex<double>(0.75,1.5)) << endl;

     ////works
     // cout << " complex?: " << NumericDataTypeChecker<double>::IsComplex << endl;
    // cout << " integral type?: "<< NumericDataTypeChecker<double>::IsIntegralType << endl;
    // cout <<  " complex?: " << NumericDataTypeChecker<complex<double> >::IsComplex << endl;
    // cout << " integral type?: "<< NumericDataTypeChecker<complex<double> >::IsIntegralType << endl;
    // cout << " complex?: " <<  NumericDataTypeChecker<long>::IsComplex << endl;
    // cout << " integral type?: "<<NumericDataTypeChecker<long>::IsIntegralType << endl;
    // cout << " complex?: " <<  NumericDataTypeChecker<complex<short> >::IsComplex << endl;
    // cout << " integral type?: "<<NumericDataTypeChecker<complex<short> >::IsIntegralType << endl;
    // cout << " complex?: " <<  NumericDataTypeChecker<int>::IsComplex << endl;
    // cout << " integral type?: "<<NumericDataTypeChecker<int>::IsIntegralType << endl;
    //  cout << " complex?: " <<  NumericDataTypeChecker<complex<long double> >::IsComplex << endl;
    // cout << " integral type?: "<<NumericDataTypeChecker<complex<long double> >::IsIntegralType << endl;

    // complex<double> z1(3.,4.), z2(-0.5, 1.25), z3(-0.75,-2.15),
    //   z4(0,2), z5(-3.5,-0.5), z6(1.,2.45);

    // MyVec<complex<double> > ZV1(3), ZV2(3);
    // ZV1[0] = z1; ZV1[1] = z2; ZV1[2] = z3;
    // ZV2[0] = z4; ZV2[1] = z5; ZV2[2] = z6;

    // cout << angle(ZV1,ZV2) << endl;

    // double af[] = {1, 0.5, -2.5, 3.15};
    // double bf[] = {0.75, -1.5, 0.5, 2.35};
    // MyVec<double> u(af,af+4), 
    //   v(bf,bf+4);
    // cout << "PROJECTIONS:"<< endl;
    // cout << "scalar proj. = "<< scalar_projection(v,u) <<endl;
    // cout << "vector proj. = "<< vector_projection(v,u) <<endl;
    // cout << "orth. proj. mtx  = "<< orthogonal_projection_onto_span_of_vector(v) <<endl;
    // cout << "angle between u and v = "<< angle(u,v) << endl;

    ////!not defined for complex numbers ;)
    // complex<double> z(3,4);
    // cout << std::acos(z) << endl;

    // //works: concurs with Matlab(R)
    // double c[] = {2,1,0,0,2,1,
    // 		  0,0,2,1,1,1};
    // double f[] = {0.035,0.01,0.015,0.04, 0.55, 0.35};
    // MyVec<double> C(c,c+2*6), b(2), u(f,f+6);
    // b[0] = 2.; b[1] = 1.;
    
    
    
    // MyVec<double> g(2);
    // //for(int i = 0; i < 3;++i){
    //   g = -matrix_vector_product(C,u) + b;
    //   cout << "-g(u) = "<< g << endl;
    //   //}
  
  }
  else if (example == 50){
    cout << "Hodgkin-Huxley nerve conduction model for Loligo vulgaris:"<<endl;
   
    CPUTime<double> cpu;

    double t0 = 0.,
      tend = 130.0;   //20.

    MyVec<double> x0(4);
    x0[0] = 0.;
    x0[0] = 0.;
    x0[0] = 0.;
    x0[0] = 0.;
    
    
    double hEstim = 0.0625,  //0.025;
      tol = 0.001;
    
    string outfile = "hodgkinhuxley.dat";

    cpu.start();

    ODE<double, HodgkinHuxley> ode(t0,tend,x0);
    //cout << "Sol (adapt. RKF) = "<<ode.adaptive_rkf(1000,hEstim,tol,outfile) <<endl;
    cout << "Sol (RKF45) = "<<ode.rkf45(1000,hEstim,tol,outfile) <<endl;
    //cout << "Sol (Implicit Euler) = "<<ode.implicit_euler(1000,7,1e-4,2.5e-02,outfile)<<endl;

    cpu.stop();

    GNUplot myplot;

    myplot("set key box left; set mxtics 4; set mytics 4; set grid xtics ytics mxtics mytics; set xlabel \"time (msec)\"; set ylabel \"membrane potential (mV)\"; set title \"Hodgkin Huxley model\"; plot for[i=2:5] \""+ outfile+ "\" using 1:i with lines lw 2 notitle");
    
  }
  else if (example == 51){
    typedef ThermoData4Mechanism<double,9> DataType;

    cout << "ODE integration:"<<endl;
    //double st[] = {0,0.15,0,0,0.65,0,0,0,0.2,1500};
    //double st[] = {0.03,0.06,0.01,0.15,0.05,0.09,0.21,0.34,0.06,1500};

    //mass frac = spec.mole * molar mass 
    double y1 = 7.0702437   * 3.200000e-02,
      y2 =      14.140487   * 2.015800e-03,
      y3 =      26.584116   * 2.802000e-02;

    cout << "sum mass fractions (must be approx. 1) = "<< y1 + y2 + y3 << endl;

    double st[] = {0.,         // O
		   y1,         // O2    <--
		   0.,         // H
		   0.,         // OH
		   y2,         // H2    <--
		   0.,         // HO2
		   0.,         // H2O2
		   0.,         // H2O
		   y3,         // N2    <--
		   1500.};

    MyVec<double> Init(st,st+DataType::nspec+1);
    cout << "Init.size() = "<< Init.size() << endl;

    //END TIME
    double tf = 2.e-04; //3.0e-4;

    ODE<double,H2GriMech3> ivp(0.,tf,Init); 

    MyVec<double> Sol;

    string ieFile; 

    CPUTime<double> cpu;
    cpu.start();
    
    string nme;

    //============= TODO: select integration method =====================
    bool useIE = true;   //true = use IE, false = try expl. method
    //===================================================================

    if(useIE){
      nme = "Implicit Euler";
      //implicit euler
      ieFile = "IE_H2GriMech3.dat";
      
      ParameterData PDat;
      PDat.read_from_file("ODEsolverSettings.dat");

      size_t msteps = PDat.get_datum<size_t>("maxsteps"),
    	mNewt = PDat.get_datum<size_t>("maxNewt");
      double ntol = PDat.get_datum<double>("ntol"), //accept Newton as converged
    	atol = PDat.get_datum<double>("atol"),         //tol for stszctrl
    	hEstim =  PDat.get_datum<double>("hEstim"),
    	Cscal = PDat.get_datum<double>("Cscal");
      
      //these are also worthwhile to play with
      double hmin = PDat.get_datum<double>("hmin"),  //maybe method behaves 'quasi-equidistant'
    	hmax =  PDat.get_datum<double>("hmax");

      Sol = ivp.implicit_euler(msteps,mNewt,ntol,atol,ieFile,hEstim,Cscal,hmin,hmax);
    }
    else{
      
      //========================================================================
      //===== WEIRD · WEIRD · WEIRD · WEIRD · WEIRD · WEIRD · WEIRD ============
      //========================================================================
      //CAUTION : 
      //---------
      //         does not work very reliably  !!!!!!!!!! 
     
      nme = "Expl. method 4 stiff ODEs";
      ieFile = "ex4stiff_H2GriMech3.dat";
      double errTol = 1.e-8,
    	fixedPtTol = 1.e-8,
    	k0 = 0.0000012,
    	St = 1.,
    	cs = 0.9995;    //~ 1
      
      size_t mxit = 4;         //max number of function iterations

      double kmin = 1.12e-08,
    	kmax = 1.e-06;
 
    
      Sol = ivp.explicit_4_stiff(errTol,fixedPtTol,k0,St,cs,mxit,ieFile,kmin,kmax);

    }//end else

     //print solution in formatted output
    GNUplot plotMePal;
    // plotMePal("scale1(t) = 1.e-2*t; scale2(t) = 1.e-03*t; set yrange [0:15]; set title \""+ nme + " -- Order of method: " + num2str(ivp.order_of_method()) + "\" ; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \" time t \"; set ylabel \" U(t=T)\"; plot for[i=2:9] \"" + ieFile + "\" using 1:i with lines lw 1 notitle, \""+ ieFile + "\" using 1:(scale1($10)) with line lw 1 notitle, \"" + ieFile + "\" using 1:(scale2($11)) with lines lw 2 lc rgb \"dark-orange\" title \"T(t) (x 1.e-03)\"");
    
    string setlog; //("set log y2;");

    //species i = 10 (i.e. the 9th species is N2 and is constant throughout),
    //hence I won't plot it
    plotMePal("set title \""+ nme + " -- Order of method: " + num2str(ivp.order_of_method()) + "\" ;set y2tics 1000,500; set y2range [1450: 3050]; set xtics 0,4e-5; set ytics nomirror; set xtics nomirror; set y2label \"Thermal energy in K\" textcolor rgb \"orange\"; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \" time t \"; set ylabel \" U(t) in mass fractions\";" + setlog + "plot for[i=2:9] \"" + ieFile + "\" using 1:i with lines lw 2 notitle, \"" +  ieFile + "\" using 1:11 with lines lw 2 lc rgb \"orange\" axis x1y2 notitle");


    //plot temperature
    // GNUplot tempplot;
    // tempplot("set title \""+ nme + " -- Order of method: " + num2str(ivp.order_of_method()) + "\" ; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \" time t in sec \"; set ylabel \" T(t=tend) in K\"; plot \"" + ieFile + "\" using 1:11 with lines lw 2 title \"thermal energy\"");

    double elti = cpu.stop();
   
    time_in_human_readable_form(elti);

    cout << endl<< "Solution = "<< Sol <<endl;
  }
  
  else if (example == 52){
    /*   double c = 0.0025, 
      d = 0.005,
      v1_max = 0.35;

    double R = d -c, //no the diameter is the radius
      l = 0;

    int ny = 9;
    double hy = R/(ny-1);
    
    string fname = "half_parabolicinlet.dat";
     ofstream of(fname.c_str(),ios_base::out);
     cout << "j        distance(f. mid)   velocity in x-dir."<<endl;
     for(int j = 0; j < ny; ++j){
       double velo = zero(parabolic_inlet_velocity(l,R,v1_max));
       cout << fixed << j << ".)      l = "<< l <<"    "   << velo << " m/s" << endl;
       of << j << "  " << velo << endl;
       l += hy;
     }
     of.close();

    
     GNUplot show;
     show("plot '"+fname+"' using 2:1 with linespoints lw 2");
*/
        //! WORKS 
  double c = 0., d = 0.005;
    int Ny = 9;

    double hy = (d-c)/(Ny-1),
      v1max = 0.6;
   
    double R = d/2., r = -R;
    string fname = "parabolicinlet.dat";
    ofstream of(fname.c_str(),ios_base::out);
    for(int j = 0; j < Ny; ++j){
      cout << " r = "<< r << "   velocity: "<< zero(parabolic_inlet_velocity(r,R,v1max)) << " m/s" << endl;
      of << j << "  " << zero(parabolic_inlet_velocity(r,R,v1max)) << endl;
      r += hy;
    }
    of.close();


    bool wannaEps = false;
    
    string toEps;
    if(wannaEps)
      toEps = "set terminal postscript eps color enhanced \"Helvetica\" 20 dashed; set output \"inletvelocityprofile.eps\";";
    
    

    //set style arrow 1 head back filled linetype 1 linecolor rgb \"dark-violet\"  linewidth 2.000 size screen 0.025,30.000,45.000;
    GNUplot plot;
    plot(toEps + "set label \"R\" font \"Times-Italic,24\" at first 0.72, first 6.35;set xtics nomirror; set ytics nomirror; set xrange [0:0.8]; set ylabel \"inlet boundary \";\
set xlabel \"velocity in x-direction (m/s)\";\
set style arrow 3 head back filled linetype 1 linecolor rgb \"black\"  linewidth 2.000 size screen 0.030,15.000,45.000;\
set arrow 1 from 0,1 to  0.2625,1 lw 2 arrowstyle 3;\
          set arrow 2 from 0,2 to  0.45,2 lw 2 arrowstyle 3; \
          set arrow 3 from 0,3 to  0.5625,3 lw 2 arrowstyle 3; \
          set arrow 4 from 0,4 to  0.6,4 lw 2 arrowstyle 3; \
set arrow 5 from 0,5 to  0.5625,5 lw 2 arrowstyle 3; \
set arrow 6 from 0,6 to  0.45,6 lw 2 arrowstyle 3; \
set arrow 7 from 0,7 to  0.2625,7 lw 2 arrowstyle 3; \
set arrow 8 from 0.7,4 to 0.7,8;\
set arrow 9 to 0.7,4 from 0.7,8 lw 1;\
set arrow 10 from 0,4 to 0.8,4 nohead lw 1 lt 0;\
          plot \"" + fname + "\" using 2:1 with lines lw 2 notitle");

   
  }
  
  else if (example == 53){
    
    // string additionalFile = "examples/automaticallygeneratedsourceterms/2Ddiscr/2Dsettings.dat";
    string additionalFile = "../fdm/data/conaireh2.dat";
    ParameterData PD;
    PD.read_from_file(additionalFile);
    double a = PD.get_datum<double>("a"),
      b = PD.get_datum<double>("b");
    double xignition = PD.get_datum<double>("x_ignition");  //(b-a)/(PD.get_datum<int>("channel_length_frac"));
    cout << "xignition: "<< xignition << " m"<<endl;
    int Nx = PD.get_datum<int>("Nx"); //1200; //PD.get_datum<int>("Nx");
    double length = b-a,
      hx = length/(Nx-1),
      x = 0.; // y = 0.;

    //============== TODO: select favorite temperature profile ===========
    const char TP = FDMSettings::WallTemperatureType; //'M'; //'S'; //'F' //'M'
    bool useMorePointsInLoop = false;
    int mxL = 1001;
    //====================================================================

    int maxLoop = Nx;
    double step = hx;
    string linefeat = "linespoints pt 3";
    if(useMorePointsInLoop){
      maxLoop = mxL;
      step = length/(maxLoop-1);
      linefeat = "lines";
    }
    string pfl;
    if(TP == 'e' || TP == 'E')
      pfl = "exponential";
    else if (TP == 'F' || TP == 'f')
      pfl = "Frouzakis";
    else if (TP == 'M' || TP == 'm')
      pfl = "Marc";
    else if (TP == 'H' || TP == 'h')
      pfl = "hyperbolic";
    else if (TP == 'S' || TP == 's')
      pfl = "natural cubic spline-like in 4 points";
    else if (TP == 'A' || TP == 'a')
      pfl = "'A'-shaped distribution";
    else
      pfl = "temp. profile NOT specified";

    pfl += (", Nx = " + Num2str(Nx)+ " ");

    double Tmin = PD.get_datum<double>("T_min"),
      Tmax = PD.get_datum<double>("T_ignition");
    typedef BoundaryTemperatureDistribution<double,TP> BoundaryWallTemperatureType;
    BoundaryWallTemperatureType BTD;
    //BTD.initialize(Tmin,Tmax,xignition,(((BoundaryWallTemperatureType::Value == 's') || (BoundaryWallTemperatureType::Value == 'S')) ? PD.get_datum<double>("cubic_spline_from_ignition")*hx : PD.get_datum<double>("plateau_width")*hx),(((BoundaryWallTemperatureType::Value == 'a') || (BoundaryWallTemperatureType::Value == 'A')) ? hx : length));
    BTD.initialize(Tmin,Tmax,xignition,(((BoundaryWallTemperatureType::Value == 's') || (BoundaryWallTemperatureType::Value == 'S')) ? PD.get_datum<double>("cubic_spline_from_ignition")*hx : PD.get_datum<double>("plateau_width")*hx),(((BoundaryWallTemperatureType::Value == 'a') || (BoundaryWallTemperatureType::Value == 'A')) ? 2*hx : length));

    string fname = "bdytempdistribution.dat";
    ofstream of(fname.c_str(),ios_base::out);
    for(int i = 0; i <  maxLoop; ++i){
      x = i*step;
      of << setprecision(12) << x << "   " << BTD(x)<< endl;
    }
    of.close();
    
    bool wannaEps = false;  //true;
    
    string toEps;
    if(wannaEps)
      toEps = "set terminal postscript eps color enhanced \"Helvetica\" 20; set output \"bdytemp.eps\";";
    GNUplot plot;
    plot(toEps+ "set xrange["+ Num2str(a)+":"+ Num2str(b)+"]; set title \"bdy temperature profile:  " + pfl +"\"; set ylabel \"temperature (K)\"; set xlabel \"boundary in x-direction (m)\"; plot '"+ fname + "' using 1:2 with "+linefeat+" lw 2 notitle;");
    

    cout << "T(0,·) = " << BTD(0.) << endl;
    cout << "step = "<< step << endl;
    cout << "domain length = "<< length << endl;
  }
   else if (example == 54){
     const int nspec = 3;
     MyVec<double> v(nspec+1);
     v[0] = 0.1; v[1] = 0.2; v[2] = 0.745;
     v[3] = 1425.;

     v[2] = 1. - UnrollLoop<0,nspec-1>::sum<1>(v);
     cout << "v[2] = "<< v[2] <<endl;
  }
  
   else if (example == 55){
     MyVec<double> v(2);
     
     v <<= 0.5, 0.75;
     cout << "v = "<< v << endl;
     
     double pi = UniversalConstants<double>::Pi;
     MyVec<double> vtest(4);
     vtest <<= -pi*2.25, 0.5*pi, 1.5, 0.75+0.75*pi; 
      
     cout << "vtest = " << vtest << endl;
     cout << entrywise_function_evaluation<double,std::sin>(vtest)<<endl;

     ADONIS_ERROR(DerivedError,"Hey Marc");
   }

   else if (example == 56){
     ParameterData PD;
     PD.read_from_file("examples/geodynamics/data2readin.dat");
     size_t Nx = PD.get_datum<size_t>("Nx"),
       Ny = PD.get_datum<size_t>("Ny");
     //size_t  npt = Nx*Ny;
     
     double a = PD.get_datum<double>("a"),
       b = PD.get_datum<double>("b");
     double c = PD.get_datum<double>("c"),
       d = PD.get_datum<double>("d"),
       Tbot = PD.get_datum<double>("Tbot"),
       Wplume = PD.get_datum<double>("Wplume"),
       Tplume = PD.get_datum<double>("Tplume");
     
     double width = b-a;       //L
     double depth = d-c;            //H
     double hx = width/(Nx-1);
     double hy = depth/(Ny-1);
     

     cout << "PARAMETERS:" << endl<<  "L: "<< width << endl<< "H: "<< depth << endl << "Tbot: "<< Tbot << endl<< "Tsurf: "<< PD.get_datum<double>("Tsurf") << endl << "Tplume: "<< Tplume << endl << "kappa: "<< PD.get_datum<double>("kappa") <<endl << "Wplume: "<< Wplume << endl; 

     
     cout << "hx = "<< hx << "  hy = "<< hy <<endl;

     MyVec<double> zvec, xvec;
     for(double j = -depth; j <= 0; j += hy)
       zvec.push_back(j);  
   
     MyVec<double> one(Nx,1.), oneX(Ny,1.);
    
     for(double i = -width/2; i <= width/2; i += hx)
       xvec.push_back(i);

     MyVec<double> z2d = dyad(zvec,one), x2d = dyad(oneX,xvec), //both 51x101
       T0;
     
     //cout << "z2dmarc = "<< z2d << endl;
     //cout << "x2dmarc = "<< x2d << endl; 
     //! Done!: z2d and x2d agree with matlab


     for(size_t k = 0; k < z2d.size(); ++k)
       T0.push_back(Abs(z2d[k]/depth)*Tbot);

    
     //*************** ***********************************************
     //Imping plume beneath lithoshere
     MyVec<size_t> ind;
     for(size_t j = 0; j < Nx; ++j){
       if(Abs(x2d[RowMajor::offset(0,j,Nx)]) <= Wplume/2)
	 ind.push_back(j);
     }
     cout << "ind = "<<ind <<endl;  //agrees with matlab

     for(size_t k = 0; k < ind.size(); ++k){
       T0[RowMajor::offset(0,ind[k],Nx)] = Tplume;
     }
         
     double t0 = 0.,
       tend = PD.get_datum<double>("tend");
     
     unsigned nt = PD.get_datum<unsigned>("nt"),
       maxNewt = PD.get_datum<unsigned>("maxNewt");
     double ntol = PD.get_datum<double>("ntol"),
       atol = PD.get_datum<double>("atol"),
       hEstim = PD.get_datum<double>("hEstim"),
       Cscal = PD.get_datum<double>("Cscal"),
       hmin = PD.get_datum<double>("hmin"),
       hmax = PD.get_datum<double>("hmax");
     
     
     cout << "# of max. steps: "<< nt << endl << "hEstim: "<< hEstim << endl << "hmin: "<< hmin<< endl; 

     const int TWOD = 2;
     const bool SPARSE = true;
     
     typedef  ODE<double,GeoTransientHeatEquation,TWOD,SPARSE> DiffEqType;
     typedef DiffEqType::IndexVecType IndexVecType;
     
     IndexVecType selectVariables; //empty

     string additionalFile = "examples/geodynamics/data2readin.dat";
     int printprec = PD.get_datum<int>("printprec");

     DiffEqType ivp(t0,tend,T0);
     CPUTime<double> cpu;
     cpu.start();
     
     //ivp.explicit_4_stiff(1e-04,1e-04,hEstim,1.,0.9675,13,"EX4STIFFGEO",hmin,hmax,selectVariables,additionalFile,printprec);
     //ivp.explicit_euler(nt,"ExplicitGEO",selectVariables,additionalFile,printprec);
     ivp.implicit_euler(nt,maxNewt,ntol,atol,"GEOHEAT",hEstim,Cscal,hmin,hmax,selectVariables,additionalFile,printprec);
     double elapsedTime = cpu.stop();
     
     time_in_human_readable_form(elapsedTime); 
  
   }
  
   else if (example == 57){
     double a[] = {2, 0.5, -1.5, 0.75};
     double b[] = {-0.65, 0.45, 2.24, -3.5};

     std::vector<double> u(4), v(a,a+4), w(b,b+4);

    
     print_all(vector_operation<AddBasicElements>(u,v,w));

   }
   else if (example == 58){
     
   
    MichaelisMentenOMalley<double> Fun(2);

    cout << " Check exact against ad jacobian "<<endl; 
    MyVec<double> seed(2);
    seed[0] = 0.54; seed[1] = -1.45;

#if USE_CPPAD == 1   
    typedef CppAD::AD<double> ADType;

    MichaelisMentenOMalley<ADType> FunAD(2);
    MyVec<ADType> X(2), Y(2);
    CppAD::Independent(X);
    Y = FunAD(X);
    CppAD::ADFun<double> adfun;
    adfun.Dependent(X,Y);

    MyVec<double> adJac = adfun.Jacobian(seed);
    print_in_matrix_style(adJac,2,"CppAD Jacobian:",12);

    cout << endl;
    MyVec<double> exJac = Fun.exact_Jacobian(seed);
    print_in_matrix_style(exJac,2,"Exact Jacobian:",12);

    cout << endl<< "Reduced derivative:"<<endl;
    MyVec<ADType> ZX(2), ZY(1);
    CppAD::Independent(ZX);
    ZY = FunAD.reduced_dynamics(ZX);
    
    CppAD::ADFun<double> af;
    af.Dependent(ZX,ZY);
    cout << "dx_ad = "<< af.Jacobian(seed) << endl;
    cout << "dx_ex =  [ "<< Fun.reduced_dynamics_Jacobian(seed) << " ]" << endl;
    cout << "... end of derivative checking :)"<<endl<< endl;
#endif

    ParameterData PD;
    PD.read_from_file("inputfiles/michaelismenten.dat");
    //Implicit Euler the handwritten way
    double t0 = 0.,
      tend = PD.get_datum<double>("tend");

    double time(t0),
      step = PD.get_datum<double>("step");
    size_t maxNewt = PD.get_datum<size_t>("maxNewt");
   
    int order = 2;
    MyVec<double> Uprev(2),U(2),b(2),IterMatrix(ntimes<2>(order)),sim(2);
    double Tol = PD.get_datum<double>("breakTOL");

    //some initial values 
    double z1_0 = PD.get_datum<double>("y1_0"),   //fixed
      z2_0 = PD.get_datum<double>("y2_0");
    
  
    string fname = "MichaelisMenten.dat";
     ofstream of(fname.c_str());
     of << "#  t   y[0]     y[1]     zM[0]     zM[1]   y_red" << endl;
    
      Uprev <<= z1_0, z2_0;
      U = Uprev;

      double yprev = z1_0;
      double yn = yprev;
 
      //Numerical solution: BACKWARD EULER -- easy implementation
      size_t count(0);
      while(time < tend){
    	sim = Fun.exact_SIM(U);  //calculate exact SIM

    	of << setprecision(14)  <<  time << "  " << U[0] << "   "<< U[1] << "  "
    	   <<  sim[0] << "   " << sim[1] << "   "<< yn << endl;

    	cout << count++ << ".)   dt = "<< step << "   time = " << time << endl;
     
    	//Newton iteration
    	U = Uprev;   //start guess
    	for(size_t l = 1; l <= maxNewt; ++l){
    	  b = -(U-Uprev-step*Fun(U));
    	  IterMatrix = -step*Fun.exact_Jacobian(U);
    	  update_diagonal<AddBasicElements>(IterMatrix,order,1.);
    	  good_square_solve(IterMatrix,order,b,1); //b contains solution now
    	  U += b;
    	  if(Norm<'2',double>::norm(b) <= Tol)
    	    break;

    	  if(l == maxNewt)
    	    ADONIS_INFO(Information,"Newton iteration did not converge within "<< maxNewt << " iterations for Tol = "<< Tol << ".");
    	} //end Newton iteraion
    	Uprev = U;    //store current approximation
    

    	//scalar version
    	yn = yprev;
    	double deltscal;
	  
	// this is for the Newton iteration for the reduced scalar system 
    	for(size_t l = 1; l <= maxNewt; ++l){
    	  deltscal  = -(yn - yprev -step*Fun.reduced_rhs(yn))/(1-step*Fun.reduced_rhs_prime(yn));
    	  yn += deltscal;
    	  if(Abs(deltscal) <= Tol)
    	    break;
    	  if(l == maxNewt)
    	    ADONIS_INFO(Information,"SCALAR Newton iteration did not converge within "<< maxNewt << " iterations for Tol = "<< Tol << ".");
    	}
    	yprev = yn;

    	time += step; //increment time (independent of Newton iteration)

      }//end time iteration
      
      of.close();

      string tolatex1, tolatex2;
      
      bool inTex = PD.get_datum<bool>("inTeX"); 
    
      string plotsize = " size "+PD.get_datum<string>("width")+PD.get_datum<string>("unit")+","+PD.get_datum<string>("height")+PD.get_datum<string>("unit")+" ",
    	texfile1 = "michaelis_menten_malley_evol"+PD.get_datum<string>("stradd"),
    	texfile2 = "michaelis_menten_malley_phase_space"+PD.get_datum<string>("stradd"),
    	setrange1, 
    	setrange2;
      bool enforceRange = PD.get_datum<bool>("enforcerange");

      cout << plotsize << endl;

      if(enforceRange){
	setrange1  = "set xrange["+ PD.get_datum<string>("xrangeStart_1")+":"+  PD.get_datum<string>("xrangeEnd_1")+"] ;" + "set yrange["+ PD.get_datum<string>("yrangeStart_1")+":"+  PD.get_datum<string>("yrangeEnd_1")+"] ;";
	setrange2 = "set xrange["+ PD.get_datum<string>("xrangeStart_2")+":"+  PD.get_datum<string>("xrangeEnd_2")+"] ;" + "set yrange["+ PD.get_datum<string>("yrangeStart_2")+":"+  PD.get_datum<string>("yrangeEnd_2")+"] ;";
      }
      
      if(inTex){
    	tolatex1 = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile1+".tex'; set format xy \"$%g$\";";
    	tolatex2 = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile2+".tex'; set format xy \"$%g$\";";
      }

      //moves (or not) tics of y-axes to the right while deleting number on the left
      string y2tics; // = "set y2tics mirror; unset ytics;";
      GNUplot show;
      show(tolatex1+setrange1+"set key box left bottom; set xlabel \"time\"; set ylabel \"$y(t)$\"; plot '"+fname+"' using 1:2 with lines lw 2 lc rgb \"black\" title \"$y_1$\", '"+fname+"' using 1:3 with lines lw 2 lc rgb \"violet\" title \"$y_2$\", '"+fname+"' using 1:6 with lines lw 3 lc rgbcolor \"#DC143C\" title \"$y_1^{\\\\mathrm{red}}$\"");
      GNUplot phasespace;
      phasespace(tolatex2+setrange2+" set key box left bottom;set xlabel \"$y_1$\";set y2label \"$y_2$\";"+y2tics+ " plot '"+fname+"' using 2:3 axes x1y2 with lines lw 2 notitle, '"+fname+"' using 4:5 axes x1y2 with lines lw 3 lc rgb \"blue\" title \"$h(y_1)$\", '"+fname+"' using 6:3 axes x1y2 with lines lw 2 lc rgb \"black\" title \"red.\"");

   }

   else if(example == 59){

   }
   else if (example == 60){
     MichaelisMentenOMalley<double> F(2);
     
     ParameterData PD;
     PD.read_from_file("inputfiles/michaelismenten.dat");


     int MAX = PD.get_datum<int>("MAX"), 
       infty = PD.get_datum<int>("pmax"),
       m = PD.get_datum<int>("m"),
       n = PD.get_datum<int>("n");

     double h = F.epsilon()/100.,
       H = n*h,
       delta(0),
       simtol = PD.get_datum<double>("simtol"),
       entry0 = 1.; 

     double vbar(PD.get_datum<double>("vattract")),
       u0 = entry0,  //reduced model ...
       vk, vkp1;
     MyVec<double> W(2);

     MyVec<double> Full(2);
     Full <<= entry0, 0.;   //and full model have same init. val. for 1st comp.

     double t(0);

     double exact(entry0);

     ofstream of("michiment.dat");
     for(int l = 0; l < MAX; ++l){  //INTEGRATION
       cout << l << ".)   time: "<< t << "   dt: "<< h << endl; 
       of << t << "   "<< u0 << "   " << Full[0]<<  "   " << exact << "  " << Abs(Full[0]-u0) << "   "<< Abs(u0-exact)   << endl;
       Full += h*F(Full);
       exact += h*F.reduced_rhs(exact);
       
       //PROJECT ONTO SIM
       int p = 1;
       while(++p <= infty){
	 cout << p << ".)  ";
	 W[0] = u0;      //fixed
	 W[1] = vbar;  
       
	 delta = 0; //reset
       
	 //int mplu1 = m+1;
	 //for(int k = 0; k <= mplu1; ++k){
	 for(int k = 0; k <= m; ++k){
	   vk = W[1]; 
	   W += H*F(W);
	   vkp1 = W[1];
	   // delta += NaturalPow(-1,mplu1-p)*binomial_coefficient(mplu1,p)*(vkp1-vk); 
	   delta += NaturalPow(-1.,k)*(vkp1-vk);
	 }
	 double adelt = Abs(delta);
	 cout << "|delta| = "<< adelt << endl;
	 if(adelt <= simtol){
	   FancyMessages().nice_output("YEAH, CONSTRAINT ITERATION HAS CONVERGED FOR GIVEN TOL");
	   break;
	 }
	 else{
	   vbar += delta;  //update vbar
	 }   
	 if(p == infty)
	   ADONIS_INFO(Information, "Max number of outer iters reached ( p = "<<infty << ").");
        
       } //end p- loop
       //integration
       double redfun = F.redfun(W);
       if(is_zero(redfun)){
	 cout << "F.redfun(W) ~ 0: Fixed point achieved"<<endl;
	 break;
       }
       u0 += h*redfun;
       t += h;
     }
     of.close();


      bool inTex = PD.get_datum<bool>("inTeX"); 
    
      string plotsize = " size "+PD.get_datum<string>("width")+PD.get_datum<string>("unit")+","+PD.get_datum<string>("height")+PD.get_datum<string>("unit")+" ",
    	texfile1 = "mime_gear",         //NO SUFFIX
    	texfile2 = "mime_gear_error",
	tolatex1,
	tolatex2;
      cout << plotsize << endl;

      
      if(inTex){
    	tolatex1 = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile1+".tex'; set format xy \"$%g$\";";
    	tolatex2 = "set terminal epslatex "+plotsize+" color colortext; set output '"+texfile2+".tex'; set format xy \"$%g$\";";
      }


     GNUplot draw;
     draw(tolatex1+"set style line 1 lw 2 pt 6 pi 150 ps 1; set xlabel \"$t$\"; set ylabel \"$y_1(t)$\";plot 'michiment.dat' using 1:2 with linespoints ls 1 title \"red. (approx.)\", 'michiment.dat' using 1:3 with lines lw 2 lc rgb '#C42F56' title \"full\", 'michiment.dat' using 1:4 with lines lw 2 title \"red. (exact)\"");


// plot1
// datafile with lines not,
// datafile every 5 with points not,
// 1/0 with linespoints t "legend"
     string linestyle1 = "set style line 1 lc rgb '#C42F56' lw 2 pt 3 pi 150 ps 1;",
       linestyle2 = "set style line 2 lc rgb 'blue' lw 2 pt 4 pi 150 ps 1;"; 
     GNUplot err;
     err(tolatex2+linestyle1 + linestyle2+"set xlabel\"$t$\";plot 'michiment.dat' using 1:5 with linespoints ls 1 notitle, 'michiment.dat' using 1:6 with linespoints ls 2 notitle");

    
   }
   else if(example == 61){
     MyVec<double> sim(2);

     RenPopeDavisSkodjie<double> RP(2);

     MyVec<double> eval(2);
     eval <<= 0.5, 0.325;

     double h = 1.e-09; //pow(std::numeric_limits<double>::epsilon(),1./3); //1.e-09

     MyVec<double> hvec(2,h);

     cout << "hvec = "<< hvec << endl;
     MyVec<double> pert(2), zM(2), zMpert(2), fddiff(2);
     pert = eval+hvec;
     zM = RP.sim(eval); 
     zMpert = RP.sim(pert);
     
     //fddiff = (zMpert - zM)/h;
     fddiff[0] = (zMpert[0] - zM[0])/h;
     fddiff[1] = (zMpert[1] - zM[1])/h;
     cout << setprecision(15) << fixed << "FD diff = "<< fddiff << "analyt. = "<< RP.tangent(eval) << endl; 
     
     
     
   }
   else if (example == 62){
     LamGoussisMech<double> LGM(7);

     MyVec<double> init(7);
     init <<= 2e-8, 0.01e-8, 0, 0, 50e-8, 0, 0;
     cout << "y(0) = "<< init << endl;
     cout << "LGM(x) = "<< LGM(init) << endl;
   }
   else if (example == 63){

#if USE_CPPAD
     const int dim = 3;
     MyVec<double> Id(dim*dim);
     Id <<= 1,0,0,
       0,1,0,
       0,0,1;

     double a = 1.,
       b = 1.,
       h = 0.00025;

     
     typedef CppAD::AD<double> ADType;
     typedef MyVec<ADType> VecADType;
     typedef CppAD::ADFun<double> ADFunType;

     VecADType X(dim), Y(dim);
     ADFunType seq;

     CppAD::Independent(X);
     
     //F(U)
     Y[0] = X[0]*X[2];
     Y[1] = sin(X[1]);
     Y[2] = 2*X[0] - X[1]*X[1];

     seq.Dependent(X,Y);
     seq.optimize();

     MyVec<double> eval(3);
     eval <<= 1.5, -0.75, 2.35;

     cout << "F'(U) = "<<endl;
     MyVec<double> Fprime = seq.Jacobian(eval);
     print_in_matrix_style(Fprime,dim);
     cout << endl 
	  << "G'(U) = a·I -bh·F'(U) = "<<endl;
     
     VecADType V(dim), W(dim), Uprev(dim);
     ADFunType qes;

     Uprev <<= -0.8, -1.5, 2.;

     CppAD::Independent(V);
     
     //G(U) := aU - Uprev - b*h*F(U)
     W[0] = a*V[0] - Uprev[0] - b*h*V[0]*V[2];
     W[1] = a*V[1] - Uprev[1] - b*h*sin(V[1]);
     W[2] = a*V[2] - Uprev[2] - b*h*2*V[0] - V[1]*V[1];

     qes.Dependent(V,W);
     qes.optimize();
     MyVec<double> Gprime = qes.Jacobian(eval);
     print_in_matrix_style(Gprime,dim);
     
     cout << endl;
     cout << "Sparsity pattern of F'(U):"<<endl;
     print_matrix_pattern(Fprime,dim);
     cout << endl << "Sparsity pattern of G'(U):"<<endl;
     print_matrix_pattern(Gprime,dim);
     cout << endl;
     cout << "Are patterns equal?:  "<< yes_no(compare_matrix_pattern(Fprime,Gprime,dim)) << endl;

#endif


      
   }
   else if (example == 64){
    
     
   }
   else if (example == 65){
     
      Dune::FieldMatrix<double,4,3> A;
      A[0][0] = 1.;   A[0][1] = 2.; A[0][2] = 3.;
      A[1][0] = -3.;  A[1][1] = 2.; A[1][2] = 1.;
      A[2][0] = 2.;   A[2][1] = 0.; A[2][2] = -1.;
      A[3][0] = 3.;  A[3][1] = -1.; A[3][2] = 2.;

      std::cout << orthonormalize(A) << std::endl; 


      Dune::FieldMatrix<double,4,4> ToInv;
      double a[4][4] = {{1, 4, 3, 2},
		      {-3,6,6,8},
		      {-6,0,5,-1},
		      {-1,-7,9,2}};

      for(int i = 0; i < 4; ++i)
	for(int j = 0; j < 4; ++j)
	  ToInv[i][j] = a[i][j];


      inversion(ToInv);

      cout << "Inv = "<< endl << ToInv << endl;
      
      cout << endl << "Test Tangent space computation"<<endl;
      
      Dune::FieldMatrix<double,6,2> B;
      B = 0.;

      B[0][0] = 1.; B[4][1] = 1.;

      //cout << "B = "<<endl<<B <<endl;

      MyVec<double> r(2),zM(6);
      r <<= 0.5, 1.235;
      zM <<= 0.25,1.5,0.75,0.15,2.0,1.5;


     

      TangentSpaceRobust<double,2,6> TS(B);
      cout << "B = "<<endl << TS.get_B() << endl;
      Dune::FieldMatrix<double, 6,2> Tr = TS.compute(r,zM);
      cout << "T(r) = " << endl <<  Tr <<endl;

      Dune::FieldMatrix<double, 4,6> NTr = normalspace_transposed<2>(Tr);
      cout << "N^T(r) = "<< endl << NTr<<endl;

      cout << "Nt·T = "<< endl << NTr*Tr <<endl;

      
      cout << endl << "Max, Min of two complex numbers, i.e. z1<z2 :<==> |z1| < |z2|:"<<endl;
      complex<double> z1(3,-4), z2(-0.5,6);

      cout << "Min(z1,z2) = "<< Min(z1,z2) << endl;
      cout << "Max(z1,z2) = "<< Max(z1,z2) << endl;

      // //!works
      // cout << endl << "Diagonal matrix *field matrix:"<<endl;
      // MyVec<double> v3(3);
      // v3 <<= 0.65, -1.5, 0.75;
      // double fin[3][4] = {{1, 2, 3, 4},
      // 			  {5, 6, 7, 8},
      // 			  {9, 10, 11, 12}};
      // Dune::FieldMatrix<double, 3, 4> fx;
      // for(int i = 0; i < 3; ++i)
      // 	for(int j = 0; j < 4;  ++j)
      // 	  fx[i][j] = fin[i][j];

      // MyVec<double> v4(4);
      // v4 <<= -1.45, 0.5, -1.65, 2.5;
      // cout << diagonal_matrix_multiplication(v3,fx) <<endl
      // 	   << diagonal_matrix_multiplication(fx,v4) << endl;
      		  

     // MyVec<double> v(5*3);
     // v <<= 1, 0.5, -3,
     //   -2, 0.75,-0.65,
     //   4, -5, -3,
     //   -0.15, 0.35, 1.75,
     //   -2.5, 1.45, 7;

     // Dune::FieldMatrix<double,5,3> Mtx = fill_from_rac<5,3>(v);
     
     // cout << "Mtx "<< endl<< Mtx << endl;

     // MyVec<double> v1(3), v2(4);

     // v1 <<= 1,2,3;
     // v2 <<= -2,3,-1,1;
     // Dune::FieldMatrix<double,3,4> dyad = dyadic_product<3,4>(v1,v2);
     // cout << "dyad = "<< endl<< dyad << endl;
     

    // const int dim = 5;
    // MyVec<double> N(dim);
    // N <<= 0.45, 12.35, 2.65, 1.55, 4.5;
    // double sum = UnrollLoop<0,dim>::sum<1>(N);
    
    // matrix<double> M(dim,dim,1.);  //initialize with 1's
    
    // for(int i = 0; i < dim; ++i)
    //   M[i][i] = (N[i]-sum)/N[i];

    // M.print_in_matlab_style("H");
    // cout << "    det(M) = "<<det(M) <<endl;

  }
   else if (example == 66){
   
     //time(NULL) gives current time
     RNGlfsr113<> Lfsr113; //(time(NULL));default seeded with time(0)
     //Lfsr113.set_seed(time(NULL)); 
     //RNGlfsr258<> Lfsr258(time(NULL));
     LCGminstd<>  Minstd; //(time(NULL));
     Minstd.set_seed(time(NULL));
     cout<< "SEED of LFSR113 = "<< Lfsr113.get_seed() <<endl;
     //cout<< "SEED of LFSR258 = "<< Lfsr258.get_seed() <<endl;
     cout<< "SEED of min std = "<< Minstd.get_seed() << endl;
     for(int i = 1; i <= 200; ++i){ 
      
       cout << "RNG due to L'Ecuyer 32 bit  = "<< Lfsr113.draw_number() << endl;
       //cout << "RNG due to L'Ecuyer 64 bit  = "<< Lfsr258.draw_number() << endl;
     
       cout << "minimal standard gen. = "<< Minstd.draw_number() << endl;
   

     }
     cout << endl << "DRAW NUMBER FROM RANGE:" <<endl;
     for(int i = 0; i < 25; ++i){
       cout << "range num = "<< setprecision(15) << Lfsr113.draw_number_from_range(2.e-13,5.e-12) << endl;
       cout << "range num (Minstd) = "<< setprecision(15) << Minstd.draw_number_from_range(2.e-13,5.e-12) << endl;
     }
     // cout << endl << "over a loop of show me the random numbers:"<<endl;
     // for(int i = 0; i < 10; ++i){
     //   cout << "random num = "<< Lfsr113.draw_number() << endl;
     // }
     
     cout << "size of type of lfrs113: " << sizeof(RNGlfsr113<>::IntType) <<endl;
     // cout << "size of type of lfrs258: " << sizeof(RNGlfsr258<>::IntType) <<endl;
     //cout << "RNGlfsr113: "<< (typeid(RNGlfsr113<>::IntType).name()) <<endl;
     //cout << "RNGlfsr258: "<< (typeid(RNGlfsr258<>::IntType).name()) <<endl;
   
     // ui64 sixtyfour = 18446744073709551614LU;
     // ui64 b = (((3 << 1) ^ 3) >> 53);
     // sixtyfour = ((3 & sixtyfour) << 10) ^ b;
     // cout << "sixtyfour = "<< sixtyfour << endl;
   }
   else if (example == 67){
#if USE_CPPAD 
     typedef CppAD::AD<double> Type;
     //typedef double Type;
     //typedef complex<double> Type;

     Type a1 = 0.0/0.0,
       a2 = 4.5,
       b1 = 1./0.0,
       b2 = -2.5,
       b3 = -1./0.,
       small = 1.e-24;
     complex<double> z1(0.0/0.0,2.),
       z2(-3, 1./0.),
       z3(3,-4),
       z4(0./0.,1./0.);
     cout << "NaN?: "<< IsNan(a1) <<endl;
     cout << "NaN?: "<< IsNan(a2) <<endl;
     cout << "Inf?: "<< IsInf(b1) <<endl;
     cout << "Inf?: "<< IsInf(b2) <<endl;
     cout << "Inf?: "<< IsInf(b3) <<endl;
     
     cout << endl << "=========================="<<endl;
     cout << "Check well-definedness of some numbers:"<<endl;
     cout << "Is "<< a1 << " well defined? : "<< is_well_defined_value(a1) << endl;
     cout << "Is "<< a2 << " well defined? : "<< is_well_defined_value(a2) << endl;
     cout << "Is "<< b1 << " well defined? : " << is_well_defined_value(b1) << endl;
     cout << "Is "<< b2 << " well defined? : "<< is_well_defined_value(b2) << endl;
     cout << "Is "<< b3 << " well defined? : "<< is_well_defined_value(b3) << endl;
     cout << "Is "<< z1 << " well defined? : "<< is_well_defined_value(z1) << endl;
     cout << "Is "<< z2 << " well defined? : "<< is_well_defined_value(z2) << endl;
     cout << "Is "<< z3 << " well defined? : "<< is_well_defined_value(z3) << endl;
     cout << "Is "<< z4 << " well defined? : "<< is_well_defined_value(z4) << endl;
     cout << "Is "<<  small << " well defined? : "<< is_well_defined_value(small) << endl;
     

     
     cout << "Test this:"<<endl;
     CppAD::AD<double> NaN = 0./0.,
       Inf = 1./0.,
       goodvalue1 = -12.34,
       goodvalue2 = 388195;

     double nan = convert_number(NaN),
       infty = convert_number(Inf),
       dn = 0./0.;

     
     cout << "NaN = "<< NaN << "   Inf = "<< Inf << endl;
     cout << "nan = "<< nan << "    infty = " << infty << endl;
     cout << "dn = "<< dn << endl;
     double gv1 = convert_number(goodvalue1),
       gv2 = convert_number(goodvalue2);
     cout << "double(goodvalue1) = "<< gv1 << " double(goodvalue2) = "<< gv2 << endl;


     double arj;
     smart_assign(arj,goodvalue1);
     cout << "arj = "<< arj << endl;
     smart_assign(arj,NaN);
     cout << "arj (NaN CppAD) = "<< arj << endl;
     smart_assign(arj,0./0.);
     cout << "arj (nan) = "<< arj << endl;
     
     // //! works
     // cout << endl << "----------- PLAUSIBILITY CHECK ON VALUES ------------"<<endl;
     // Type a = sqrt(-3.5);

     // cout << "WELL DEFINED?: "<< is_well_defined_value(a) <<endl;

     // cout << "Is value ok?: " << endl;
     // value_plausibility_check(a); 

     // Type NaN = 0./0.;
     // MyVec<Type> v(5);
     // v <<= 1.0, -3.65, NaN, 0.6, 1.2;
     // cout << "Is vector ok?:"<<endl; 
     // value_plausibility_check(v); 
#endif
   }
   else if (example == 68){
     typedef MyVec<double> VType;
     VType w1(5*6), w2(2*6), g(2*6);

     for(int i = 0; i < 30; ++i)
       w1[i] = i+1;

     w2 <<= 2,1,0,0,2,1,
       0,0,2,1,1,1;

     MyVec<double> w3(4);
     cout << "w3 = "<<w3 << endl; 

     g <<= -3,9,4,2,1,4,
       -5,1,3,7,2,9;

     w3 <<= 4,5,
       6,7;

     cout << "w3 (assigned) = "<< w3 << endl;

     Dune::FieldMatrix<double,2,6> FM;
     FM = 0;  //important
     FM[0][0] = 1;
     FM[1][4] = 1;


     Dune::FieldMatrix<double,2,3> MF;
     MF = 0;
     MF[0][0] = -14;  MF[0][1] = 15;  MF[0][2] = 16;
     MF[1][0] = -17;  MF[1][1] = 18;  MF[1][2] = 19;

     ConcatenateMatrix<2,6,VType> CM,  CM2;
     CM.columnwise(5,6,w1,FM);
     
     

     cout << "COL concatenated matrix: "<< CM.data() << endl;
     
     CM2.columnwise(FM,2,6,w2);
     cout <<endl<< "COL concatenated matrix: "<< CM2.data() << endl;
     
     //CM.reset_flag();
     
     CM.rowwise(2,6,g,FM);
     cout << endl<< "CM.size() = "<< CM.size() << endl;
     cout <<endl<< "ROW concatenated matrix: "<< CM.data() << endl;

     ConcatenateMatrix<2,3,VType> MC;
     MC.rowwise(2,6,g,MF);
     cout <<endl<< "ROW concatenated matrix: "<< MC.data() << endl;

     //MC.reset_flag();
     MC.rowwise(MF,2,6,g);
     cout <<endl<< "ROW concatenated matrix: "<< MC.data() << endl;
     
     cout <<endl<< "MATRIX for 6 spec mech:"<<endl;
     ConcatenateMatrix<2,6,VType> CB;
     CB.columnwise(2,6,w2,FM);
     cout <<endl<< "COL concatenated matrix: "<< CB.data() << endl;
     

     MyVec<double> u(7);
     MyVec<double> v(3),w(4);

     v <<= 0.5, 2, -0.15;
     w <<= 1.75, -0.65, 3.2, -3.5;

     concatenate(u,v,w);
     cout << "u (after concat.) = "<< u << endl;

   }

   else if(example == 69){
    //example form http://www.nag.co.uk/numeric/fl/nagdoc_fl22/xhtml/F08/f08kcf.xml

    MyVec<double> A(5*6), b(5), res;

    A <<= -0.09, -1.56,	-1.48,	-1.09,	0.08,	-1.59,
      0.14,	0.20,	-0.43,	0.84,	0.55,  -0.72,
      -0.46,	0.29,	0.89,	0.77,	-1.13,	1.06,
      0.68,	1.09,	-0.71,	2.11,	0.14,	1.24,
      1.29,	0.51,	-0.96,	-1.27,	1.74,	0.34;

    cout << endl;
    print_in_matrix_style(A,6," A  ",12,true);

    b<<= 7.4, 4.3,-8.1, 1.8, 8.7;
    cout << "b = "<< b << endl;
    
    res = solve_rank_deficient_ls(A,5,b,1);

    cout << "res = "<< res << endl;
   

    cout <<endl<< "overdetermined system from [GLOLUB, LOAN, Ex. 5.3.2, p. 239]:"<<endl;
    MyVec<double> A1(3*2), b1(3);
    A1 <<= 1, 1,
      1.e-03, 0, 
      0, 1.e-03;
    b1 <<= 2,1.e-03,1.e-03;
    
    MyVec<double> ser = solve_rank_deficient_ls(A1,3,b1,1);
    cout <<" compare solution with that of Golub: x_LS = [1,1]^T: "<<endl;
    cout << "x = "<< ser << endl;
  }

  else if (example == 70){
   
  }

  else if (example == 71){
    std::cout << setprecision(13) << endl;
    typedef MyVec<double> VType;
    
    typedef ThermoData4Mechanism<double,6> DataType;
    typedef StoichiometricMatrix<typename DataType::index_type*> StoichType;

    typedef NASA7CoefficientPolynomial<typename DataType::value_type*,typename DataType::value_type*> ThermoType;
    
    typedef ForwardReactionRates<typename DataType::value_type*,typename DataType::value_type*, typename DataType::int_type*, typename DataType::value_type*, 'a'> FwdType;

  
    typedef IsoThermalReverseReactionRates<VType,StoichType> RevType;

    //=========================== DIFFERENCE ==================================
    typedef ReducedIndexInjector<const size_t*, const size_t*> InjectorType;
    //=========================================================================

    typedef BuildChemicalSourceTerm<FwdType,RevType,InjectorType> BuilderType;
    //typedef typename BuilderType::RevPointerType RevPointerType;
  

    static ThermoType nasa(DataType::nspec,DataType::thermo(),DataType::temperature_bounds());
	

    
    static StoichType stoichmatrix(DataType::nspec,DataType::nreac,DataType::stoichiometric_matrix());

    
    static TroeIndex TIndex;
    TIndex.create(DataType::nreac,DataType::troewtb());
	
    //! everything stated as  pointers	
    static FwdType forwardRates(DataType::nspec,DataType::nreac,DataType::arrhenius_values(),DataType::temperature_bounds(), DataType::ntroereac, DataType::troewtb(),DataType::collision_efficiencies(),&TIndex);

    //! everything stated as  pointers	
    //static RevType reverseRates(&forwardRates,&stoichmatrix,&nasa);
    //!NOTE: mechanism is isothermal!!!
    static double rr[] = {216., 337.5, 1400., 10800., 33750., 0.7714}; //!these are fixed reverse rates for this mechanism
    static RevType reverseRates(rr,rr+DataType::nspec,&stoichmatrix);
    
    //=========================== DIFFERENCE ==============================
    static InjectorType Inj(2,DataType::rednreac,DataType::rpv_index(), DataType::rpv_reaction_index());
    //=====================================================================
    BuilderType BCS_;
    BCS_.initialize(&forwardRates,&reverseRates,&Inj);
	

    //Full state estimate
    MyVec<double> z(6);
    random_container(z);

    const int RED = 2;
    MyVec<double> dotomega(RED);
    for(int k = 0; k < RED; ++k){
      
      dotomega[k] = BCS_.net_formation_rate(BCS_.species_index(k),double(1000.),z);
     
    }

    cout << "AUTOMATICALLY:"<<endl;
    cout << "dotomega = "<<dotomega << endl;

    cout << endl << "HANDCODED:"<<endl;

    double pk[12];
    double rV[12];

     pk[0] = 2.0;
      pk[1] = 1.0;
      pk[2] = 1.0;
      pk[3] = 1000.0;
      pk[4] = 1000.0;
      pk[5] = 100.0;
      pk[6] = 216.0;
      pk[7] = 337.5;
      pk[8] = 1400.0;
      pk[9] = 10800.0;
      pk[10] = 33750.0;
      pk[11] = 0.7714;

      rV[0]  =  pk[0]*z[0];
      rV[1]  =  pk[6]*z[1]*z[1];
      
      rV[4]  =  pk[2]*z[4];
      rV[5]  =  pk[8]*z[1]*z[5];
      rV[6]  =  pk[3]*z[0]*z[3];
      rV[7]  =  pk[9]*z[1]*z[5];
      
      rV[10] =  pk[5]*z[0]*z[3];
      rV[11] =  pk[11]*z[4];

      MyVec<double> dotHand(2);

      dotHand[0] = - rV[0]  + rV[1]              //H2
      - rV[6]  + rV[7]
      - rV[10] + rV[11];

      dotHand[1] = - rV[4]  + rV[5]               //H2O
      + rV[10] - rV[11];


      cout << "dotHand = "<< dotHand << endl;

    // //works
    // HodgkinHuxley<double> HoHux(4);
    // MyVec<double> v(4);
    // v <<= 0.75, 0.235, 0.95, 1.75;

    // MyVec<double> feval = HoHux(v);
    //  cout << "||feval||_2 = "<< Norm<'2',double>::norm(feval) << endl; 

    // cout << "||F(v)||_2 = "<< Norm<'2',double>::norm(HoHux(v)) << endl; 


    // bool george = false;
    // cout << "george (before) = "<< nice_bool_out(george) <<endl;
    // change_some_boolean(george);

    // cout << "george (after) = "<< nice_bool_out(george) <<endl;

  }

  else if (example == 72){
    cout << "preparing the Conaire H2..."<<endl;
    const int nreac = 19;
    const int nspec = 10;
    const int twicenreac = nreac*2;

    int S[twicenreac][nspec];

    for(int i = 0; i < twicenreac; ++i)
      for(int j = 0; j < nspec; ++j)
	S[i][j] = 0;

    enum{O,O2,H,OH,H2,HO2,H2O2,H2O,N2,AR};

    const string specNames[nspec] = {"O","O2","H","OH","H2","HO2","H2O2","H2O","N2","AR"};

    //1.)
    S[0][H] = 1;     S[0][O2] = 1;   //fwd    
    S[1][O] = 1;     S[1][OH] = 1;   //rev
    
    //2.)
    S[2][O] = 1;     S[2][H2] = 1;   //fwd
    S[3][H] = 1;     S[3][OH] = 1;   //rev

    //3.)
    S[4][OH] = 1;    S[4][H2] = 1;  //fwd
    S[5][H] = 1;     S[5][H2O] = 1;  //rev

    //4.)
    S[6][O] = 1;     S[6][H2O] = 1;  //fwd
    S[7][OH] = 2;                     //rev


    //5.)
    S[8][H2] = 1;                     //3RD BODY: 1
    S[9][H] = 2;
 
    //6.)
    S[10][O2] = 1;                    //3RD BODY: 2
    S[11][O] = 2;

    //7.)
    S[12][OH] = 1;
    S[13][O] = 1;    S[13][H] = 1;     //3RD BODY: 3

    //8.)
    S[14][H2O] = 1;
    S[15][H] = 1;    S[15][OH] = 1;     //3RD BODY: 4

    //9.)
    S[16][H] = 1;    S[16][O2] = 1;     //TROE: 1
    S[17][HO2] = 1; 
 
    //10.)
    S[18][HO2] = 1;   S[18][H] = 1;
    S[19][H2] = 1;    S[19][O2] = 1;
    
    //11.)
    S[20][HO2] = 1;   S[20][H] = 1;
    S[21][OH] = 2;
			
    //12.)
    S[22][HO2] = 1;   S[22][O] = 1;
    S[23][OH] = 1;    S[23][O2] = 1;

    //13.)
    S[24][HO2] = 1;   S[24][OH] = 1;
    S[25][H2O] = 1;   S[25][O2] = 1;

    //14.)
    S[26][H2O2] = 1;  S[26][O2] = 1;  //DUPLICATE; sum of the 2 rate expressions
    S[27][HO2] = 2;                   //canceled by 2 

    //15.)
    S[28][H2O2] = 1;                 //TROE: 2
    S[29][OH] = 2;  
    
    //16.)
    S[30][H2O2] = 1;  S[30][H] = 1;
    S[31][H2O] = 1;   S[31][OH] = 1;
    
    //17.)
    S[32][H2O2] = 1;  S[32][H] = 1;
    S[33][H2] = 1;    S[33][HO2] = 1;

    //18.)
    S[34][H2O2] = 1;  S[34][O] = 1; 
    S[35][OH]  = 1;   S[35][HO2] = 1;
    
    //19.) 
    S[36][H2O2] = 1;  S[36][OH] = 1; //DUPLICATE; sum of the 2 rate expressions
    S[37][H2O] = 1;   S[37][HO2] = 1; //canceled by 2
    

    //print stuff 
    cout << "STOICHIOMETRIC MATRIX for H2 mechanism of Conaire et al.:"<<endl;
    for(int k = 0; k < nspec; ++k)
      cout << specNames[k] << "  ";
    cout << endl;
    for(int i = 0; i < twicenreac; ++i){
      for(int j = 0; j < nspec; ++j){
	cout << S[i][j] << ",  ";
      }
      cout << endl;
    }
  
    //    A        b        Ea     
 static double fwdOrig[3*19] = {
        1.915E+14,   0.00,  1.644E+04,
        5.080E+04,   2.67,  6.292E+03,
        2.160E+08,   1.51,  3.430E+03, 
        2.970E+06,   2.02,  1.340E+04,
        4.577E+19,  -1.40,  1.044E+05,
        4.515E+17,  -0.64,  1.189E+05,
        9.880E+17,  -0.74,  1.021E+05,
        1.912E+23,  -1.83,  1.185E+05,
        1.475E+12, 0.60, 0.000E+00,      //TROE
        1.660E+13,   0.00,  8.230E+02,
        7.079E+13,   0.00,  2.950E+02,
        3.250E+13,   0.00,  0.000E+00,
        2.890E+13,   0.00, -4.970E+02,
        4.634E+16+1.434E+13,  -0.35+-0.35,  5.067E+04+3.706E+04, //DUP: sum up
        2.951E+14,   0.00,  4.843E+04,    //TROE
        2.410E+13,   0.00,  3.970E+03,
        6.025E+13,   0.00,  7.950E+03,
        9.550E+06,   2.00,  3.970E+03,
        1.000E+12+5.800E+14,   0.00+0.00,  0.000E+00+9.557E+03 //DUP: sum up
        
 };
 cout << endl << "COEFFICIENTS FOR REACTION RATES:"<<endl;
 cout << "     A               b               Ea"<<endl;
    //    A        b        Ea     
 static double revOrig[3*19] = {
     5.481E+11,   0.39, -2.930E+02,
     2.667E+04,   2.65,  4.880E+03,
     2.298E+09,   1.40,  1.832E+04,
     1.465E+05,   2.11, -2.904E+03,
     1.146E+20,  -1.68,  8.200E+02,
     6.165E+15,  -0.50,  0.000E+00,
     4.714E+18,  -1.00,  0.000E+00,
     4.500E+22,  -2.00,  0.000E+00,
     3.090E+12, 0.53, 4.887E+04,      //TROE -- commented
     3.164E+12,   0.35, 5.551E+04,
     2.027E+10,   0.72,  3.684E+04,
     3.252E+12,   0.33,  5.328E+0,  //CHANGED3.250E+13,   0.00,  0.000E+00,
     5.861E+13,   0.24,  6.908E+04,
     4.200E+14+1.300E+11,   0.00+0.00,  1.198E+04+-1.629E+03, //DUP: sum up
     3.656E+08,   1.14, -2.584E+03,     //TROE: commented
     1.269E+08,   1.31,  7.141E+04,
     1.041E+11,   0.70,  2.395E+04,
     8.660E+03,   2.68,  1.856E+04,
     1.838E+10+1.066E+13,   0.59+0.59,  3.089E+04+4.045E+04  //DUP: sum up
     
     
 };



//          fwd rev
  int molec[2*19] = { 
            2,2,  //1.)
            2,2,
            2,2,
            2,2,
            2,3,
            2,3,
            2,3,
            2,3,
            3,2,  //9.) TROE
            2,2,
            2,2,
            2,2,
            2,2,
            2,2, //14.)  DUP
            2,3, //15.) TROE
            2,2,
            2,2,
            2,2,
            2,2  //19.) DUP
    };

  int troe[2] = {8,14};
   
  // const double lowPressTroe[2*3] = {
  //   3.4820E+16, -4.1100E-01, -1.1150E+03,
  //   1.202E+17,  0.00, 45500.
  // }; 

   cout << endl << "FORWARD RATES" << scientific << endl;
   double val; 
   for(int i = 0; i < 19; ++i){
        
      
      if(molec[2*i] == 2){
	if(i == troe[0] || i == troe[1])
	  val = fwdOrig[3*i];
	else 
	  val = 1.0e-06*fwdOrig[3*i];
	cout << val;
      }
      if(molec[2*i] == 3){
        if(i == troe[0] || i == troe[1])
	  val = 1.0e-06*fwdOrig[3*i];
	else
	  val =1.0e-12*fwdOrig[3*i];
	cout << val;
      }  
      cout << ",   " << fwdOrig[3*i+1] << ",   "<< cal2joule(fwdOrig[3*i+2]) << ","<< endl;
        
    }
   cout << endl<< "REVERSE RATES" << endl; 
    for(int i = 0; i < 19; ++i){
        
      
      if(molec[2*i+1] == 2){
	if(i == troe[0] || i == troe[1])
	  val = revOrig[3*i];
	else 
	  val = 1.0e-06*revOrig[3*i];
	cout << val;
      }
      if(molec[2*i+1] == 3){
        if(i == troe[0] || i == troe[1])
	  val = 1.0e-06*revOrig[3*i];
	else
	  val =1.0e-12*revOrig[3*i];
	cout << val;
      }  
      cout << ",   " << revOrig[3*i+1] << ",   "<< cal2joule(revOrig[3*i+2]) << ","<< endl;
        
        
    }
   
    cout << endl<< "TROE index: "<< endl;
    
    typedef ThermoData4Mechanism<double,10> DataType;
    TroeIndex troeidx;
    troeidx.create(DataType::nreac,DataType::troewtb());
    cout << troeidx << endl;


    cout << endl << "DUPLICATES (all reacs bimolecular T-dependence):" <<endl;
    const int numOfDupReacs = 2;
    const double dups[numOfDupReacs*(3*4)] = {
      //h2o2+o2 = ho2+ho2
      4.634E+16,  -0.35,  5.067E+04,  //fwd
      4.200E+14,   0.00,  1.198E+04,  //rev
      1.434E+13,  -0.35,  3.706E+04,  //fwd1
      1.300E+11,   0.00, -1.629E+03,  //rev1
      //h2o2+oh = h2o+ho2
      1.000E+12,   0.00,  0.000E+00,  //fwd
      1.838E+10,   0.59,  3.089E+04, //rev
      5.800E+14,   0.00,  9.557E+03,  //fwd1
      1.066E+13,   0.59,  4.045E+04,  //rev1 
 };
    
    //all bimolecular
    string notify;
    for(int l = 0; l < 8; ++l){
      cout << (1.0e-06*dups[3*l]) << ",  " << dups[3*l+1] << ",  " << cal2joule(dups[3*l+2]) << ",    //" <<( ((l%2==0) ? ("fwd") : ("rev")) ) << endl;
    }
  }

  else if (example == 73){
    cout << "CHEMKIN default H2 mechanism:"<<endl;
    enum{O, O2, H, OH, H2, HO2, H2O2, H2O, N2}; //keep order 

    const int nspec = 9;
    const int nreac = 20;
    const int twicereac = 2*nreac;
    
    int S[twicereac][nspec];
    for(int i = 0; i < twicereac; ++i)
      for(int j = 0; j < nreac; ++j)
	S[i][j] = 0;

    //! split every = reaction into two forward ones and count the stoich. coeffs
    //1.)
    S[0][H] = 1;  S[0][O2] = 1;      //3rd body: 1
    S[1][HO2] = 1;

    //2.)
    S[2][H] = 2;  
    S[3][H2] = 1;
    
    //3.)
    S[4][H] = 2; S[4][H2] = 1;
    S[5][H2] = 2;

    //4.)
    S[6][H] = 2; S[6][H2O] = 1;
    S[7][H2] = 1;  S[7][H2O] = 1;

    // 5.)
    S[8][H] = 1;  S[8][OH] = 1;   //3rd body: 2
    S[9][H2O] = 1;

    //6.)
    S[10][H] = 1; S[10][O] = 1;   //3rd body: 3
    S[11][OH] = 1;
    
    //7.)
    S[12][O] = 2;
    S[13][O2] = 1;

    //8.)
    S[14][H2O2] = 1;
    S[15][OH] =2;
    
    //9.)
    S[16][H2] = 1;  S[16][O2] = 1;
    S[17][OH] = 2;
    
    //10.)
    S[18][OH] = 1;  S[18][H2] = 1;
    S[19][H2O] = 1;  S[19][H] = 1;
    
    //11.)
    S[20][O] = 1; S[20][OH] = 1;
    S[21][O2] = 1; S[21][H] = 1;

    //12.)
    S[22][O] = 1; S[22][H2] = 1;
    S[23][OH] = 1; S[23][H] = 1;

    //13.) 
    S[24][OH] = 1; S[24][HO2] = 1;
    S[25][H2O] = 1; S[25][O2] = 1;

    //14.)
    S[26][H] = 1;  S[26][HO2] = 1;
    S[27][OH] = 2;
    
    //15.)
    S[28][O] = 1; S[28][HO2] = 1;
    S[29][O2] = 1; S[29][OH] = 1;

    //16.)
    S[30][OH] = 2;  
    S[31][O] = 1;  S[31][H2O] = 1;
    
    //17.)
    S[32][H] = 1; S[32][HO2] = 1;
    S[33][H2] = 1; S[33][O2] = 1;
    
    //18.)
    S[34][HO2] = 2;
    S[35][H2O2] = 1; S[35][O2] = 1;

    //19.)
    S[36][H2O2] = 1;  S[36][H] = 1;
    S[37][HO2] = 1;  S[37][H2] = 1;

    //20.)
    S[38][H2O2] = 1; S[38][OH] = 1;
    S[39][H2O] = 1; S[39][HO2] = 1;
    
    cout << "STOICHIOMETRIC MATRIX for H2 mechanism CHEMKIN:"<<endl;
    for(int i = 0; i < twicereac; ++i){
      for(int j = 0; j < nspec; ++j){
	cout << S[i][j] << ",  ";
      }
      cout << endl;
    }

    double fwdOrig[3*nreac] = {
      3.61E17,  -0.72,       0.,  //3rd body:1

      1.0E18,   -1.0,        0.,
      9.2E16,   -0.6,        0.,
      6.0E19,   -1.25,       0.,
      1.6E22,   -2.0,        0., //3rd body:2

      6.2E16,   -0.6,        0.,   //3rd body:3

      1.89E13,   0.0,    -1788.,
      1.3E17,    0.0,    45500.,
      1.7E13,    0.0,    47780.,
      1.17E9,    1.3,     3626.,
      3.61E14,  -0.5,        0.,
      5.06E4,    2.67,    6290.,
      7.5E12,    0.0,      0.0,
      1.4E14,    0.0,     1073.,
      1.4E13,    0.0,     1073.,
      6.0E+8,    1.3,        0.,
      1.25E13,   0.0,        0.,
      2.0E12,    0.0,        0.,
      1.6E12,    0.0,     3800.,
      1.0E13,    0.0,     1800.
    };
    

    //! this time we only consider forward molecularities since the reverse
    //! rates are computed using the equilibrium constant
    //! NOTE: count third bodies as well!
    int molec[nreac] = { 
      3,
      3,
      3,
      3,
      3,
      3,
      3,
      2,
      2,
      2,
      2,
      2,
      2,
      2,
      2,
      2,
      2,
      2,
      2,
      2
    };

    bool statedincalories = true;

    if(statedincalories)
      cout << endl<<"CAL/mol, i.e. value x 4.184"<<endl<<endl;
    else
      cout << endl<<"kJ/mol,  i.e. value x 1e+03"<<endl<<endl; 

    double val;
    cout << scientific <<  "#  A      beta      Ea"<<endl;
    for(int i = 0; i < nreac; ++i){
      if(molec[i] == 2){
	val = 1.0e-06*fwdOrig[3*i];
	cout << val;
      }
      if(molec[i] == 3){
	val =1.0e-12*fwdOrig[3*i];
	cout << val;
      }                                        //1.e+03   //! in kJ/mole ????
      cout << ",   " << fwdOrig[3*i+1] << ",   "<< ((statedincalories==true) ?  cal2joule(fwdOrig[3*i+2]):val*=1.e+03) << ","<< endl;
    }
  }  

    else if (example == 74){
      cout << "H2/O2 FUEL-AIR EQUIVALENCE RATIO:"<<endl;
      cout << "n mole of H2 vs n mole of O2:"<<endl;
      //stoich: H2 + 1/2 O2  ==> H2O
      //same mole number cancel out
      double mH2devmO2 = (2*1.007900e-03)/(2*1.600000e-02);
      double mH2devmO2stoich = (1*2*1.007900e-03)/(0.5*2*1.600000e-02);
      double phi = mH2devmO2/mH2devmO2stoich;

      cout << "mH2/mH2 = "<< mH2devmO2 << "    (mH2/O2)_stoich = "<<  mH2devmO2stoich << endl;

      cout << "phi = "<< phi << endl;
      //Since for phi = 0.5 the number of moles of H2 and O2 is equal, so are their mass mass fractions, i.e. X_O2 = X_H2.
      //H2/air mixture: H2/O2/N2

      double Xfrac_same = 0.3; 

      double W_H2 = 2.015800e-03,
	W_O2 = 3.200000e-02,
	W_N2 = 2.802000e-02;
      
	double Wbar = Xfrac_same*(W_H2 + W_O2) + (1-2*Xfrac_same)*W_N2;
	
	double Y_H2 = (Xfrac_same*W_H2)/Wbar,
	  Y_O2 = (Xfrac_same*W_O2)/Wbar;

	cout << "H2:O2 (by mass) = 1 : "<< Y_O2/Y_H2 << endl;  

	cout << "Y_H2 = "<< Y_H2 << endl;
	cout << "Y_O2 = "<< Y_O2 << endl;
	cout << "Y_N2 = "<< 1 - (Y_H2 + Y_O2) << endl;


    }

    else if (example == 75){
      double x = 0,
	xignite = 0.00025,
	hx =  0.005/(59-1),
	tmin = 300.,
	tmax = 960;

      BoundaryTemperatureDistribution<double,'M'> BWT(tmin,tmax,xignite,3*hx);
      ofstream of("hyptang.dat",ios_base::out);
      for (int i = 0; i < 59; ++i){
	of << x << "   " <<  //hyp_tang_T_distribution(x,xignite,tmin,tmax,1.,2*hx) << endl;
	  BWT(x) << endl;
	x += hx;
      }
      of.close();
    
      GNUplot pres;
      pres("set xlabel \"x in m\"; set ylabel \"temperature in K\"; plot 'hyptang.dat' using 1:2 with linespoints lw 3");
    }
  
    else if (example == 76){
      double val = 3.5;
      complex<double> zval(-3,-4);

      cout << "test Sqrt:"<<endl;
      cout << Sqrt(val) << endl;
      cout << Sqrt(zval) << endl;
      cout << Sqrt(0./0.) << endl;
      cout << Sqrt(1./0.) << endl;
      int tos = 16;
       cout << Sqrt(tos) << endl;
      //cout << Sqrt(-4) << endl;  //catch error

#if USE_CPPAD
      CppAD::AD<double> valad = 3.5,
	negvalad = -4.;
      CppAD::AD<complex<double> > zvalad = zval;
      
      cout << "test Sqrt:"<<endl;
      cout << Sqrt(valad) << endl;
      cout << Sqrt(zvalad) << endl;
      CppAD::AD<double> nanad = 0./0.,
	infad = 1./0.;
      cout << Sqrt(nanad) << endl;
      cout << Sqrt(infad) << endl;
      //CppAD::AD<int> tosad = 16;
      //cout << Sqrt(tosad) << endl;
      //cout << Sqrt(negvalad) << endl; //catch error
      

      // complex<int> zi(3,4);
      // typedef complex<double>::value_type cplxD_Type;
#endif

    }
    else if (example == 77){
#if USE_CPPAD
      CppAD::AD<double> val1 = 3500,
	val2 = -8.5;
      double fac = 100;
      
      cout << " |3500| >> |-8.5|*100 ?: "<<  yes_no(much_larger_than_in_abs(val1,val2,fac)) <<endl;

      double v1 = -4,
	v2 = 3;

       cout << " |-4| >> |3|*100 ?: "<<  yes_no(much_larger_than_in_abs(v1,v2,fac)) <<endl;
      


       cout << endl << "Check for good values:"<<endl;
       typedef CppAD::AD<double> Type;
       MyVec<Type> vec(4);
       Type badval = numeric_limits<double>::infinity()/numeric_limits<double>::infinity();
       cout << "badval = "<< badval << endl;
       vec <<= 1.5, -6.5, badval, 2.45;

       cout << "is vec ok?: "<< yes_no(is_plausible_value(vec)) << endl;
    
       Type lavdab = 0./0.;
       cout << "is scalar ok?: "<< yes_no(is_plausible_value(lavdab)) << endl;

#endif
}
    else if (example == 78){
#if USE_CPPAD 
      //cp. http://www.coin-or.org/CppAD/Doc/sparse_jacobian.cpp.htm
      //    http://www.coin-or.org/CppAD/Doc/cppad_sparse_jacobian.cpp.htm

      //! both types work
      typedef std::set<std::size_t> SimpleType;
      //typedef bool SimpleType;
      
      typedef std::vector<SimpleType> SVecType;
      //typedef VectorSet<SimpleType> VSetType;

      SVecType p;
      
      typedef CppAD::AD<double> ADType;
      int m = 3, 
	n = 3;
      MyVec<ADType> X(n),Y(m);
      CppAD::Independent(X);
      Y[0] = X[0]*X[0]*X[1];
      Y[1] = X[2];
      Y[2] = X[1] + X[2]*X[2];
    
      CppAD::ADFun<double> f;
      f.Dependent(X,Y);
      f.optimize();

      //VSetType::sparsity(p,f);

      SparsityPattern<SimpleType> SP(m,n);
      SP.calc_sparsity(f);
      cout << "Sparsity pattern" << endl << SP << endl;
      cout << "size = "<< SP.size() << "   #(nnz) = " << SP.number_of_nonzeros() << endl;

      cout << "Sparse jacobian evaluation:"<<endl;
      MyVec<double> eval1(3);
      eval1 <<= -0.5, 0.75, -1.5;

      MyVec<double> Jac1 = f.SparseJacobian(eval1,SP.get_p());
      print_in_matrix_style(Jac1,n);
      cout << endl<< "exact Jacobian = "<<endl;
      MyVec<double> Exa1(m*n);
      Exa1 <<= 2*eval1[0]*eval1[1], eval1[0]*eval1[0],0,
	0,  0, 1,
	0, 1, 2*eval1[2];
      print_in_matrix_style(Exa1,n);
      cout << endl;


      // for(int i =0; i < n; ++i){   //Works
      // 	SP.get_p()[i].insert(i);
      // }
      SP.diagonal_never_zero();
      cout << "Sparsity pattern NEW (diag)" << endl << SP << endl;
      cout << "size = "<< SP.size() << "   #(nnz) = " << SP.number_of_nonzeros() << endl;
      
      //2nd example -- reverse
      m = 3;
      n = 4;
      MyVec<ADType> a_x(n), a_y(m);
      CppAD::Independent(a_x);
      a_y[0] = a_x[0] + a_x[1];
      a_y[1] = a_x[2] + a_x[3];
      a_y[2] = a_x[0] + a_x[1] + a_x[2] + a_x[3] * a_x[3] / 2.;
      CppAD::ADFun<double> adseq;
      adseq.Dependent(a_x,a_y);
      adseq.optimize();
      
      SparsityPattern<SimpleType> SP2(m,n);
      SP2.calc_sparsity(adseq);
      cout << "Sparsity pattern 2" << endl << SP2 << endl;
      cout << "size = "<< SP2.size() <<"   #(nnz) = " << SP2.number_of_nonzeros() << endl;
      //3rd example
      m = 4;
      n = 2;
      MyVec<ADType> D(n), R(m);
      CppAD::Independent(D);
      R[0] = D[0] + D[1];
      R[1] = D[0]*D[1];
      R[2] = D[1]*D[1];
      R[3] = D[1]*1.5;

      CppAD::ADFun<double> qws;
      qws.Dependent(D,R);
      qws.optimize();
      
      SparsityPattern<SimpleType> SP3(m,n);
      SP3.calc_sparsity(qws);
      cout << "Sparsity pattern 3" << endl << SP3 << endl;
      cout << "size = "<< SP3.size() <<"   #(nnz) = " << SP3.number_of_nonzeros() << endl;

      cout << endl << "simple 2D example's sparsity pattern"<<endl;
      typedef TwoDParabolic<ADType> Model;
      const int nx = Model::NX;
      const int ny = Model::NY;

      int dim = nx*ny;
      Model TwoD(dim);

      MyVec<ADType> V(dim), W(dim);
      CppAD::Independent(V);

      W = TwoD(V);
      CppAD::ADFun<double> sq;
      sq.Dependent(V,W);
      sq.optimize();

      
      typedef SparsityPattern<SimpleType> twdpat_type;
      typedef twdpat_type::SVecType StdVecBoolType;
      twdpat_type Sp2d(dim,dim);
      
      StdVecBoolType smtx = Sp2d.calc_sparsity(sq); //sparsity of jacobian
      Sp2d.diagonal_never_zero(); //no empty sets after invocation
      cout << "# nonempty elems = " << Sp2d.number_of_nonempty_elements() << endl;
      //cout << "Sparsity pattern 2D" << endl << Sp2d << endl;
      cout << "size = "<< Sp2d.size() <<"   #(nnz) = " << Sp2d.number_of_nonzeros() << " (out of "<< Sp2d.size()<< " x "<< Sp2d.size() << " = "<< ntimes<2>(Sp2d.size())<< ")" << endl;
      
      string fname = "ExampleFileSparseMtx",
	mtxdatafile = fname+".dat";
      //print_matrix_pattern_2_file(smtx,dim,mtxdatafile);
      Sp2d.print_sparsity_pattern_2_file(mtxdatafile);

      string pltmtxpath = "../fdm/sh/./plotmatrix.sh " + mtxdatafile;
      system(pltmtxpath.c_str());
      string epsfile = fname+".eps";
      system(("gimp "+epsfile).c_str());

      cout << "Jacobian = "<<endl;
      MyVec<double> eval2(dim);
      random_container(eval2);
      Jac1 = sq.SparseJacobian(eval2,Sp2d);
      
#endif
    }

    else if (example == 79){
#if USE_CPPAD
        //! both types work
      typedef std::set<std::size_t> SimpleType;
      //typedef bool SimpleType;
      
      typedef std::vector<SimpleType> SVecType;
      //typedef VectorSet<SimpleType> VSetType;

      SVecType p;
      
      typedef CppAD::AD<double> ADType;
      int m = 3, 
	n = 3;
      MyVec<ADType> X(n),Y(m);
      CppAD::Independent(X);
      Y[0] = X[0]*X[0]*X[1];
      Y[1] = X[2];
      Y[2] = X[1] + X[2]*X[2];
    
      CppAD::ADFun<double> f;
      f.Dependent(X,Y);
      f.optimize();

      //VSetType::sparsity(p,f);
      
      SparsityPattern<SimpleType> SP(m,n);
      SP.calc_sparsity(f);
      cout << "Sparsity pattern" << endl << SP << endl;
      cout << "size = "<< SP.size() << "   #(nnz) = " << SP.number_of_nonzeros() << endl;

      cout << endl << "SPARSE CPPAD DRIVER "<<endl;
      typedef SparsityPattern<SimpleType>::SVecType PattType;
      typedef MyVec<double> MyDouVecType;
      SparseCppADDriver<double,PattType,MyDouVecType> SDrive(f,SP.get_p());
      
      cout << "Sparse jacobian evaluation:"<<endl;
      MyVec<double> eval1(3);
      eval1 <<= -0.5, 0.75, -1.5;

      cout << "Jac = "<< SDrive.sparse_jacobian(eval1) << endl;
#endif
    }
  
else if (example == 80){
  cout << "binary diffusion prefactor = "<< PracticalPrefactor<double>::prefac_diffcoeff_g_per_mol << endl;

  double temperature = 1350;      //K
  double pressure = 101325;       //Pa
  cout << "T = "<<temperature << " K" <<endl;
  cout << "p = "<<pressure << " Pa = " << ChooseRightPrefactor<double,true>::in_bar(pressure) << " bar" << endl;

  cout << "Binary diffusion:"<<endl;
  PureSpeciesDiffusion<ThermoData4Mechanism<double,10>,false> PSD(ThermoData4Mechanism<double,10>::transport());
  PSD.compute_binary_diffusion_matrix(pressure,temperature);

  //                dim·(dim+1)/2
  for(int i = 0; i < 55; ++i)
    cout << PSD.get_binary_diffusion_coefficient(i) << "  ";
  cout << endl;

  
  cout << endl << "MIXTURE-AVERAGED DIFFUSION COEFFICIENTS" <<endl;
  MyVec<double> Xfrac(10);
  Xfrac[1] = 0.3; Xfrac[4] = 0.1; Xfrac[8] = 0.6; 
  MixtureAveragedDiffusion<ThermoData4Mechanism<double,10>,false> MAD(ThermoData4Mechanism<double,10>::transport());
  MAD.compute_mixture_averaged_diffusion_coefficients(pressure,temperature,Xfrac);

  for(int i= 0; i < 10; ++i)
    cout << MAD[i] << "  ";
  cout << endl;
  
  cout << endl << "COMPUTE D_{O2,H2O} BY HAND:" <<endl;
  double epsilon_kj = sqrt(107.4*572.4);   //K
  double Tstar = temperature/epsilon_kj;   //unitless
  double sigma_kj = 0.5*(3.45800+2.60500); //Angström
  double Omega_11_star = 1.0548*pow(Tstar,-0.15504) + pow((Tstar+0.55909),-2.1705); //unitless, //delta^*max = 0 ==> f^(1,1)* = 1

  const double fac = 3./8.*sqrt(1000)*1.e+15*sqrt(M_PI*ntimes<3>(PhysicalConstants<double>::Rgas)/ntimes<2>(PhysicalConstants<double>::AvogadrosConstant))/M_PI; 
  double Wk =  32.;      //g/mol
  double Wj =  18.0158;  //g/mol
  double pInBar = pressure*1.e-05; //bar
  double D_kj = fac*sqrt(ntimes<3>(temperature)*(Wk+Wj)/(2*Wk*Wj))/(pInBar*ntimes<2>(sigma_kj)*Omega_11_star); //Wk,Wj enter in g/mol, p in bar and sigma_kj in Angström for well-scaleness!!!! The overall bin. diffusion coeff is in m²/s (SI)

  cout << "D_kj (hand coded)  = "<< D_kj << " m²/s"<< endl;
  cout << "D_kj (tmp version) = "<< PSD.get_binary_diffusion_coefficient(SymmetricDenseMatrixAccess<1,7,10>::offset()) << " m²/s" << endl; //I = O2 = 1, J = H2O = 7


  // //here, nothing was changed in the source code
//   cout << "Viscosity:"<<endl;
//   cout << "mu = "<< PracticalPrefactor<double>::prefac_viscoeff_g_per_mol << endl;
//   PureSpeciesViscosity<ThermoData4Mechanism<double,10>,false> PSV(ThermoData4Mechanism<double,10>::transport());
//   PSV.compute_viscosities(temperature);
//   for(int i= 0; i < 10; ++i)
//     cout << PSV.get_viscosity(i) << "  ";
//   cout << endl;

//   cout << "N_A = "<< PhysicalConstants<double>::AvogadrosConstant <<endl;
//   double mprefac = 3./8.*sqrt(1000)*1.e+15*sqrt(M_PI*ntimes<3>(PhysicalConstants<double>::Rgas)/ntimes<2>(PhysicalConstants<double>::AvogadrosConstant))/M_PI; //3./8.*(sqrt(M_PI*ntimes<3>(PhysicalConstants<double>::Rgas)/ntimes<2>(PhysicalConstants<double>::AvogadrosConstant)))/M_PI*1.e+15*sqrt(1000);
//   cout << "Dbin-prefac by hand = " << mprefac << endl;
//   cout << "mu-prefac by hand   =" << (5./16.*sqrt(1./1000.)*sqrt(M_PI*PhysicalConstants<double>::Rgas/ntimes<2>(PhysicalConstants<double>::AvogadrosConstant))*1.e+20/M_PI) //(5./16.*sqrt(M_PI*PhysicalConstants<double>::Rgas/ntimes<2>(PhysicalConstants<double>::AvogadrosConstant))*sqrt(1./1000.)*1.e+20/M_PI) 
// <<endl;
}
 else if (example == 81){

   int nx = 7, ny = 5;
   
   MyVec<char> vec(nx*ny);
   

   for(int i = 0; i < nx; ++i){
     for(int j = 0; j < ny; ++j){
       int off = i +nx*j;
        if(j == 0){ 
      	  vec[off] = 'D';
      	}
      	else if(j == ny-1){ 
      	  vec[off] = 'U';
      	}
      	else if(i == 0){   
	  vec[off] = 'L';          
      }
      	else if(i == nx-1){  
	  vec[off] = 'R';
      }
      	else 
      	  vec[off] = 'x';
     }
   } 
     
   // print_in_matrix_style(vec.begin(),vec.end(),nx);

   cout << vec <<endl;
   
 }

 else if (example == 82){
  cout << "Conaire H2 mechanism:"<<endl;
    ParameterData PD;
    PD.read_from_file("../fdm/data/conaireh2.dat");
    
    int nx = PD.get_datum<int>("Nx");
    //int ny = PD.get_datum<int>("Ny");
    double a = PD.get_datum<double>("a"),
      b = PD.get_datum<double>("b");
     
     double length = b-a;   //length in x-direction of rectangular domain
     double hx = length/(nx-1);
     cout << "hx = "<<hx << endl;
    //initial mass fraction vector Y
    MyVec<double> Y(10);
    Y <<=  PD.get_datum<double>("O"), 
      PD.get_datum<double>("O2"),
      PD.get_datum<double>("H"),
      PD.get_datum<double>("OH"),
      PD.get_datum<double>("H2"),
      PD.get_datum<double>("HO2"),
      PD.get_datum<double>("H2O2"),
      PD.get_datum<double>("H2O"),
      PD.get_datum<double>("N2"),
      PD.get_datum<double>("AR");


    cout << "Y_O2 = "<< Y[1] << ",   Y_H2 = "<< Y[4] << ",   Y_N2 =" << Y[8] << endl;

    double p0 = 101325;   //in Pascal 
    
    typedef ThermoData4Mechanism<double,10> DataType;

    double wbar(0.);
    double temperatureIn(300.);

    for(int k = 0; k < 10; ++k){
      wbar += Y[k]/DataType::molar_masses()[k];
    }
    wbar = 1./wbar;

    cout << "Wbar = "<< wbar << endl;

    double rhoIn = (p0*wbar)/(PhysicalConstants<double>::Rgas*temperatureIn);

    double p1 = 120019; //(rhoIn*PhysicalConstants<double>::Rgas*temperatureIn)/wbar;

    cout << "p1 = "<< p1 << endl;
    cout << "p_ANDERSON = rho·Rgas·T = "<< (rhoIn*PhysicalConstants<double>::Rgas* temperatureIn) << endl; 
    double dx_p = (p1 - p0)/hx;
    cout << "dx_p = "<< dx_p << endl;

    
    

 }

 else if (example == 83){
   const int I = 17;
   const int J = 18;

   const int I1 = 23;
   const int J1 = 23;

   cout << "delta(I,J) = "<< (!Kronecker<I,J>::delta() )<< endl;
   cout << "delta(I1,J1) = "<< (!Kronecker<I1,J1>::delta() )<< endl;
 }
 
 else if (example == 84){
   cout << "Test the fucking matrix read in "<< endl;

   ReadInMatrix<double> RIM;
   RIM.read_from_file("matrix.dat");

   cout << "RIM.rows() = "<< RIM.rows() << endl;
   cout << "RIM.cols() = "<< RIM.cols() << endl;

   RIM.print_matrix("A");
   cout << endl << "I read in " << RIM.number_of_read_lines() << " lines in total."<<endl;
	
   cout << endl << "OMEGA u BDY:" <<endl;
   RIM.clear();
   RIM.read_from_file("discretized2Dphysdomain.dat");
   RIM.print_matrix("Omega");

   cout << endl << "try to access elements in 2D MOL fashion:" << endl;
   //========================= CHANGE THESE VALUES APPROPRIATLY =======================
   const int Nx = 6;
   const int Ny = 4;
   const int numberOfQuantities = 3;
   //==================================================================================
	
    
	
   bool withPressure = false;
   
   cout << endl << "ASSIGN:" << endl;
   
   int nofquant = (withPressure==true) ?  numberOfQuantities : (numberOfQuantities-1);
   vector<double> v(Nx*Ny*nofquant);

   RIM.assign_2_2D_MOL_slide(v,Nx,Ny,numberOfQuantities,withPressure);

   double hx = 0.00015;
   double hy = 0.0001;

   //usual output as in print MOL data
   ofstream ofile;
   ofile.open("RewrittenOmega.dat");
   cout << "write to file" << endl;
   double x(0), y(0);
   for(int i = 0; i < Nx; ++i){
     x = i*hx;
     y = 0;
     for(int j = 0; j < Ny; ++j){
       y = j*hy;
       ofile <<" " << x << "  " << y << "  ";
       for(int k = 0; k < nofquant; ++k){
	 ofile << v[i + Nx*j + Nx*Ny*k] << "  ";	
       }
       ofile << endl;
     }
    
   }
   cout << "Done!" << endl;

   


 }
 else if (example == 85){
   // //! use NaturalPow only for integers powers!!!
   // cout << "27^(1./3) = "<< NaturalPow(27.,-3) << endl; 
   // cout << "pow = "<< pow(27,1./3.) << endl;
 
   cout << endl << endl << "=================================================================" << endl;
   cout <<               "======== INITIAL CONCENTRATIONS of H2-Air combustion ============" << endl;
   cout << "=================================================================" << endl;
      
   //!equivalence ratio phi = 0.5 means 17.4 % H2 by volume, cf. [VARMA et al., "Studies of premixed laminar hydrogen-air flames using elementary and global kinetics models", Combustion and Flame 64, 1986, pp. 233--236]
   //!from this follows that there are 82.6 % air
   //!We assume that air is composed by volume of
   //! 78.08 % N2, 20.95 % O2 and 0.97 % Ar. The volumetric parts of N2 and O2 are rounded while Ar also contains the remaining atmospheric species by volumne.
   //!Based on the fact that we 82.6 % air in the mixture which is composed of the aforementioned volumetric parts, we have:
   //! H2: 17.4 %
   //! N2: 64.49408 %
   //! O2: 17.30470 %
   //! Ar:  0.80122 %
   //! That means the mixture on a molar basis is given by
   //! 0.174 H2 + 0.6449408 N2 + 0.1730470 O2 + 0.0080122 Ar.
   //! The mole fraction can now be calculated using the molar masses (either in gram or kilogram), viz.
   const double molarMass_H2 = 2.015800e-03,
     molarMass_N2 = 2.802000e-02,
     molarMass_O2 = 3.200000e-02,
     molarMass_Ar = 3.9948e-02;

   double O2_percent_by_vol = 0.18473710, //0.1730470,
     H2_percent_by_vol = 0.11820000,//0.11800, //0.174, //you can also write here 17.4%, etc
     N2_percent_by_vol = 0.68850944, //0.6449408,
     Ar_percent_by_vol = 0.00855346; //0.0080122;

   cout << setprecision(10)<< "Sum of % parts = "<< (H2_percent_by_vol+N2_percent_by_vol+O2_percent_by_vol+Ar_percent_by_vol) << endl;
   
   //!cf. http://homepages.cae.wisc.edu/~hessel/faqs/what_are_approximate_mole_and_mass_ratiosOfN2AndO2ForAir.htm
   double mass_H2 = H2_percent_by_vol*molarMass_H2,
     mass_N2 = N2_percent_by_vol*molarMass_N2,
     mass_O2 = O2_percent_by_vol*molarMass_O2,
     mass_Ar = Ar_percent_by_vol*molarMass_Ar;
 
   double totalMass = mass_H2 + mass_N2 + mass_O2 + mass_Ar;
 
   double Y_H2 = mass_H2/totalMass,
     Y_N2 = mass_N2/totalMass,
     Y_O2 = mass_O2/totalMass,
     Y_Ar = mass_Ar/totalMass;
            
   cout << "Mass fraction for phi = 0.5 of H2-air mixture:"<<endl;
   cout << "Y_O2 = "<< Y_O2 << endl;
   cout << fixed << "Y_H2 = "<< Y_H2 << endl;
   cout << "Y_N2 = "<< Y_N2 << endl;
   cout << "Y_Ar = "<< Y_Ar << endl;
            
   cout << endl << "Sum_k Y_k = " << Y_H2+Y_N2+Y_O2+Y_Ar << endl;


   //                    N2     O2     Ar 
   double airByVol[3] = {0.7808, 0.2095, 0.0097};
   double H2ByVol(0.); //this is in percent [0,1.]
   double AirVol, O2ByVol;
   enum{N2 = 0, O2, Ar};
   
   ofstream ofile("fuel_fuel_equivalence_ratio_Phi.dat", ios_base::out);
   ofile << "# O2          H2            N2          Ar" << endl;
   for( ; H2ByVol <= 1.; H2ByVol += 0.0001){
     AirVol = 1. - H2ByVol;
     O2ByVol = AirVol*airByVol[O2];
     ofile << fixed << setprecision(8) << O2ByVol << "   "<< H2ByVol << "   " <<   (AirVol*airByVol[N2]) << "   "<< (AirVol*airByVol[Ar]) << "       Phi = "<<  (0.5*(H2ByVol/O2ByVol)) <<endl;
   }

   cout<< "check sum = "<< (0.18471615 +  0.11830000 + 0.68843136 +  0.00855249) <<endl;
 }

 else if (example == 86){
  
   enum{BrO3=0,Br,HBrO2,P,CeIV};

   MyVec<double> S(5);
    S[BrO3] = 0.68; 
    S[Br] = 1.e-09; 
    S[HBrO2] = 1.e-10;
    S[P] = 2.45;
    S[CeIV] = 9.4;
   cout << endl << "STIFF OREGONATOR:" << endl;
    const int PRINTFLAG = 0; //no method of lines printing
      const bool SOLVERFLAG = false; //false;  //works for true and false
      //                                               (sparse LA) (dense LA)
      const bool SCALE = false; 
      const bool PATTERNCREATIONVIAAD = true; //true = AD, false = FD

      ODE<double, Oregonator,PRINTFLAG,SOLVERFLAG,SCALE,PATTERNCREATIONVIAAD> Oregon(0.,25.,S);
    
      //cout<<"Order of Method "<<Constant<double>::accuracyOfMethod<<endl;
      CPUTime<double> time;
      time.start();
      string name;
      //cout << setprecision(16) << endl;
    
    
      name = "Oregonator.dat";
    
      size_t steps_impl = 200000000,
	maxi = 5;
      double rto = 1.e-8,
	stpzctrl = 1.e-8,
	hestim = 1.e-7, //0.01;   //see Ascher/Petzold 100·h_n <= 1 
	Cscal = 1.111111;
    
    
      cout<< "IMplicit method = "<< setprecision(9)<< Oregon.implicit_euler(steps_impl,maxi,rto,stpzctrl,name,hestim,Cscal) 
	  <<endl;
      cout << "Implicit Euler -- calculation time: "
	   << time.stop() <<endl;
    
      GNUplot plotme;
      plotme("set key box left; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \"time t\"; set ylabel \"Y(t =T)\"; plot \"" + name + "\" using 1:2 with lines lw 3");

      // which_addition_is_used();
  

   //     CppAD::AD<double> d1 = 1.75, d2 = -1.75, d3 = 1.75000000000000034;
   // complex<double> z1(3.,-4), z2(3.0000003,-4.00000000998), z3(3,4);
 
   // cout << "Check for approximated equality:"<<endl;
   // cout << " double:"<<endl;
   // cout << is_approximately_equal(d1,d2) << endl;
   // cout << is_approximately_equal(d1,d3) << endl;
   
   // cout << " cplx:"<<endl;
   // cout << is_approximately_equal(z1,z2,1.e-6) << endl;
   // cout << is_approximately_equal(z1,z3) << endl;

 }
 else if (example == 87){
   FancyMessages FM;
   FM.nice_caution("Check this out, pal");

   cout << endl << "Sparsity pattern of Gprime:"<<endl;
   int dim = 40*40*14;  //grid dim * num of phys specs
   
   MyVec<double> Y(dim);
   ReadInMatrix<double> RIM;
   RIM.read_from_file("CONAIRE_H2_2D_FULL0.gnu");
   RIM.assign_2_2D_MOL_slide(Y,40,40,15,false); //plus pressure
   
   ConaireH2ReactiveNS2D<double> fun(dim);

   cout << "Compute G':"<<endl;
   typedef JacobianMatrixCalculator<true,double,ConaireH2ReactiveNS2D,ExprTmpl::MyVec,'C',int> JacobianType;  //!last argument must be int here (coz of UMFPACK)
   JacobianType Gprime(true,fun);   //true = calculate G'
   Gprime.initialize(dim,dim,Y);    //proper initialization with Y_, set up G'
 
   //! not needed for pattern generation of G'
   // Gprime.evaluate_jacobian(Y);
   // Gprime.compute_G_prime(5.e-12,1.,1.);
   
   cout << "Calculate sparsity pattern of G':"<<endl;
   typedef JacobianType::PatternType PatternType;
   PatternType Sp2d = Gprime.get_pattern();

   cout << "# nonempty elems = " << Sp2d.number_of_nonempty_elements() << endl;
   //cout << "Sparsity pattern 2D" << endl << Sp2d << endl;
   cout << "size = "<< Sp2d.size() <<"   #(nnz) = " << Sp2d.number_of_nonzeros() << " (out of "<< Sp2d.size()<< " x "<< Sp2d.size() << " = "<< ntimes<2>(Sp2d.size())<< ")" << endl;
      
   string fname = "SparsityPattern_Gprime_CONAIRE_H2",
     mtxdatafile = fname+".dat";
   //print_matrix_pattern_2_file(smtx,dim,mtxdatafile);
   Sp2d.print_sparsity_pattern_2_file(mtxdatafile);

   string pltmtxpath = "../fdm/sh/./plotmatrix.sh " + mtxdatafile;
   system(pltmtxpath.c_str());
   string epsfile = fname+".eps";
   system(("gimp "+epsfile).c_str());

   
 }

 else if(example == 88){
#ifndef NONCONSERVATIVE_FORM
   ADONIS_ERROR(ConfigurationError,"You must uncomment '#define NONCONSERVATIVE_FORM' in file ../fdm/gensettings.hh");
#endif
   
   ParameterData PM;
   PM.read_from_file("datafiles/burgers.dat");

   int Nx(PM.get_datum<int>("Nx"));
   double xl(PM.get_datum<double>("a")), xu(PM.get_datum<double>("b")), vis(PM.get_datum<double>("vis"));
   double delta_x = (xu-xl)/(Nx-1);
   double t0(PM.get_datum<double>("t0")), tf(PM.get_datum<double>("tf")),
     hEstim(PM.get_datum<double>("hEstim"));


   //set up IC, cf. [SCHIESSER and GRIFFITHS, "A Compendium of PDEs Models - Method of Lines Analysis with Matlab", chap. 5]
   MyVec<double> u0(Nx);
   double x(0);
   for(int i = 0; i < Nx; ++i){
     x = xl + i*delta_x;
     //analytical solution
     //assign
     u0[i] = analytical_burger_s_equation(x,t0,vis);
   }

    MyVec<int> index;

    //! print every 0.1 secs until tf = 1. s
    //Original Burgers equation -- use same data but different file names
    ODE<double,OriginalBurgersEquationSchiesser,1,true,false> OriginalMol1DBurger(t0,tf,u0);
    OriginalMol1DBurger.implicit_euler(20000000,3,1.e-05,1.e-05,"SCHIESSER_Burger_1D_",hEstim,1.11111, 1.e-03, 0.5,index, "datafiles/burgers.dat",12,0,0,-1,0.1);
   

   //My version of Burgers equation
   ODE<double,BurgersEquationSchiesser,1,true,false> Mol1DBurger(t0,tf,u0);
   Mol1DBurger.implicit_euler(20000000,3,1.e-05,1.e-05,"MF_Burger_1D_",hEstim,1.11111, 1.e-03, 0.5,index, "datafiles/burgers.dat",12,0,0,-1,0.1);

   GNUplot plotme;
   plotme("schiesserfile(n) = sprintf(\"SCHIESSER_Burger_1D_%d.gnu\", n); marcfile(n) = sprintf(\"MF_Burger_1D_%d.gnu\", n);  set title \"Burger's equation, t~0.,0.1, 0.2,...,1.\"; set xlabel \"x\"; set ylabel \"u(x,t)\"; plot for [i=0:10] schiesserfile(i) using 1:2 with linespoints lw 2 pt 7, for [i=0:10] marcfile(i) using 1:2 with linespoints lw 2 pt 4");

 }
 else if (example == 89){
   cout << "Van der Pol equations" << endl;
   double t0 = 0.,
     tf = 3000.;
   MyVec<double> U0(2);
   U0 <<= 2., 0.;
   MyVec<int> index;
   //Original equation
   ODE<double,VanderPolOscillator> VanPol(t0,tf,U0);
   VanPol.implicit_euler(20000000,3,1.e-03,1.e-06,"VanPol.dat",15.,1.11111,1.e-04,50.,index);
   GNUplot gp;
   gp("set autoscale xy; set yrange [-2.5:2.0]; set xlabel \"t\"; plot 'VanPol.dat' using 1:2 with linespoints lw 2 pt 6 title \"y_1\"");
   
   // typedef double Type1;
   // typedef CppAD::AD<double> Type2;
   
   // cout << yes_no(are_types_equal<Type1,Type2>(),true) << endl;
   // cout << yes_no(are_types_equal<complex<int>, complex<int> >(),false) << endl;
   // cout << are_types_equal<complex<float>, complex<int> >() << endl;
   // cout <<  yes_no(are_types_equal<complex<float>, CppAD::AD<Type1> >(),true) << endl; 
   // cout << are_types_equal<Type1, CppAD::AD<typename TypeAdapter<Type2>::BaseType> >() << endl;

   // cout << endl << "check whether a file exists" << endl;
   // cout << yes_no(does_file_exist("Oregonator.dat")) << endl;
   // cout << yes_no(does_file_exist("Some_rather_sillyFile.dat")) << endl;
   // cout << yes_no(does_file_exist("H2_Gri.dat")) << endl;
 }
 else if (example == 90){

   
   //! is this now defined???
   //! do not print primitive variables since we don't solve any eqs in
   //! conservative form
   //! does not yield any results
// #ifndef NONCONSERVATIVE_FORM
// #define NONCONSERVATIVE_FORM
//    std::cout << "Conservative form of equations assumed...NOT OK since you print nonsense otherwise. \nTry to blend via defining the macro now..." <<std::endl;
// #else
//    std::cout <<"Non-conservative form used. OK"<<std::endl;
// #endif

   //=========================== TODO ======================================
   bool myFunctor = true;  //true: use Marc CDR, false: use Schiesser CDR
   //=======================================================================

#ifndef NONCONSERVATIVE_FORM
   ADONIS_ERROR(ConfigurationError, "You must use '#define NONCONSERVATIVE_FORM'.\n    To be UNcommented in file ../fdm/gensettings.hh");
#endif
   
   ParameterData PM;
   PM.read_from_file("datafiles/cylindrical_CDR.dat");
   int nz = PM.get_datum<int>("Nx");
   int nr = PM.get_datum<int>("Ny");
   double ca0 = PM.get_datum<double>("ca0"),
     Tk0 = PM.get_datum<double>("Tk0"),
     hEstim = PM.get_datum<double>("hEstim"),
     hMin = PM.get_datum<double>("hMin"),
     hMax = PM.get_datum<double>("hMax");

   cout << "ca0 = "<< ca0 << ",  Tk0 = "<< Tk0 << endl;
     
   //initial conditions
   MyVec<double> U0(2*nr*nz);
   int ij, npt = nz*nr;
   for(int i = 0; i < nz; ++i){
     for(int j = 0; j < nr; ++j){
       ij = i+nz*j;
       U0[ij]=ca0;         //CONC
       U0[ij+npt]=Tk0;     //TEMP 
     }
   }

   cout << "U0 = "<< U0 << endl;
   cout << "U0.size() = "<<U0.size()<<endl; 
   
   
   //ODE integration
   double t0 = 0.,
     tf = 200.;
   MyVec<int> index;

   string outfilename = "CDR_CYL_SCHIESSER_2D_";
   if(myFunctor){
     cout << "Using MARC's Convection-diffusion-reaction eqs" << endl;
     outfilename = "MARC_CDR_CYL_2D_";
   }
   //Original equation
   ODE<double,CDRCylindrical2DSchiesser,2,true,false> CDRSchiesser(t0,tf,U0);
   //my functor
   ODE<double,MarcCDRCylindrical2D,2,true,false> MarcCDR(t0,tf,U0);

   if(myFunctor == false)
     CDRSchiesser.implicit_euler(20000000,3,1.e-04,1.e-04,outfilename,hEstim,1.11111,hMin,hMax,index, "datafiles/cylindrical_CDR.dat",12,0,0,-1,50.);
   else{
     MarcCDR.implicit_euler(20000000,3,1.e-04,1.e-04,outfilename,hEstim,1.11111,hMin,hMax,index, "datafiles/cylindrical_CDR.dat",12,0,0,-1,50.);
     cout << endl<<"============================="<<endl;
     cout << "  MARC's CDR functor used."<<endl;
     cout << "============================="<<endl;
   }
  
  
  
 }
 else if(example == 91){
   cout << "Print slides from 'example == 90':"<<endl;

   //! ============================= TODO ==============================
   bool printFirst = true;  //true: fixed z-axis; false: fixed r-axis
   //! =================================================================

   
   string yrange1 = "set yrange [305:450]; ",
     yrange2 = "set yrange [0:0.01];";
   ReadInMatrix<double> RIM;
   ReadInMatrix<double> MIR; //this is for cross-checking with my functor
   
   if(printFirst){
     cout << "fixed z-axis..."<<endl;
     RIM.read_from_file("CDR_CYL_SCHIESSER_2D_1.gnu");
     RIM.write_2_file("Out_x_0_t_50.dat","test this");
     //! CAUTION: x and y are just reversed (as in the Matlab code examples)

     MIR.read_from_file("MARC_CDR_CYL_2D_1.gnu");
     MIR.write_2_file("Marc_Out_x_0_t_50.dat","test this");
     
     //! r = 0, t = 50, 
     RIM.fixed_y_2D(0., "x_0_t_50.dat");
     MIR.fixed_y_2D(0., "marc_x_0_t_50.dat");
     GNUplot plot1;
     plot1("set xlabel \"z\"; set ylabel \"c(r=0,z,t=50)\"; plot \"x_0_t_50.dat\" using 1:3 with linespoints lw 2 pt 4 title \"Schiesser\", \"marc_x_0_t_50.dat\" using 1:3 with linespoints lw 2 pt 7 title \"Marc\"");
     GNUplot plot2;
     plot2(yrange1+"set xlabel \"z\"; set ylabel \"T(r=0,z,t=50)\"; plot \"x_0_t_50.dat\" using 1:4 with linespoints lw 2 pt 4 title \"Schiesser\", \"marc_x_0_t_50.dat\" using 1:4 with linespoints lw 2 pt 7 title \"Marc\"");


     //r = 0, t = 100
     RIM.clear(); //prepare for next read in
     RIM.read_from_file("CDR_CYL_SCHIESSER_2D_2.gnu");
     RIM.fixed_y_2D(0., "x_0_t_100.dat");

     MIR.clear(); //prepare for next read in
     MIR.read_from_file("MARC_CDR_CYL_2D_2.gnu");
     MIR.fixed_y_2D(0., "marc_x_0_t_100.dat");
     
     GNUplot plot3;
     plot3(yrange2+"set xlabel \"z\"; set ylabel \"c(r=0,z,t=100)\"; plot \"x_0_t_100.dat\" using 1:3 with linespoints pt 4 lw 2 title \"Schiesser\", \"marc_x_0_t_100.dat\" using 1:3 with linespoints pt 7 lw 2 title \"Marc\"");
     GNUplot plot4;
     plot4(yrange1+"set xlabel \"z\"; set ylabel \"T(r=0,z,t=100)\"; plot \"x_0_t_100.dat\" using 1:4 with linespoints pt 4 lw 2 title \"Schiesser\", \"marc_x_0_t_100.dat\" using 1:4 with linespoints pt 7 lw 2 title \"Marc\"");
   
     //! r = 0, t = 150
     RIM.clear(); //prepare for next read in
     RIM.read_from_file("CDR_CYL_SCHIESSER_2D_3.gnu");
     RIM.fixed_y_2D(0., "x_0_t_150.dat");

     MIR.clear(); //prepare for next read in
     MIR.read_from_file("MARC_CDR_CYL_2D_3.gnu");
     MIR.fixed_y_2D(0., "marc_x_0_t_150.dat");
     
     GNUplot plot5;
     plot5(yrange2+"set xlabel \"z\"; set ylabel \"c(r=0,z,t=150)\"; plot \"x_0_t_150.dat\" using 1:3 with linespoints pt 4 lw 2 title \"Schiesser\", \"marc_x_0_t_150.dat\" using 1:3 with linespoints pt 7 lw 2 title \"Marc\"");
     GNUplot plot6;
     plot6(yrange1+"set xlabel \"z\"; set ylabel \"T(r=0,z,t=150)\"; plot \"x_0_t_150.dat\" using 1:4 with linespoints pt 4 lw 2 title \"Schiesser\", \"marc_x_0_t_150.dat\" using 1:4 with linespoints pt 7 lw 2 title \"Marc\"");

     //! r = 0, t = 200
     RIM.clear(); //prepare for next read in
     RIM.read_from_file("CDR_CYL_SCHIESSER_2D_4.gnu");
     RIM.fixed_y_2D(0., "x_0_t_200.dat");

     MIR.clear(); //prepare for next read in
     MIR.read_from_file("MARC_CDR_CYL_2D_4.gnu");
     MIR.fixed_y_2D(0., "marc_x_0_t_200.dat");
     
     GNUplot plot7;
     plot7(yrange2+"set xlabel \"z\"; set ylabel \"c(r=0,z,t=200)\"; plot \"x_0_t_200.dat\" using 1:3 with linespoints pt 4 lw 2 title \"Schiesser\", \"marc_x_0_t_200.dat\" using 1:3 with linespoints pt 7 lw 2 title \"Marc\"");
     GNUplot plot8;
     plot8(yrange1+"set xlabel \"z\"; set ylabel \"T(r=0,z,t=200)\"; plot \"x_0_t_200.dat\" using 1:4 with linespoints pt 4 lw 2 title \"Schiesser\", \"marc_x_0_t_200.dat\" using 1:4 with linespoints pt 7 lw 2 title \"Marc\" ");
     
   }
   else{ //!print r-AXIS FIXED
     cout << "r-axis fixed..."<< endl;
     //! z = 5, t = 200
     RIM.clear(); //prepare for next read in
     RIM.read_from_file("CDR_CYL_SCHIESSER_2D_4.gnu");
     RIM.fixed_x_2D(5.26315789474, "y_5_t_200.dat");

     MIR.clear(); //prepare for next read in
     MIR.read_from_file("MARC_CDR_CYL_2D_4.gnu");
     MIR.fixed_x_2D(5.26315789474, "marc_y_5_t_200.dat");
     
     GNUplot plot9;
     plot9(yrange2+"set xlabel \"r\"; set ylabel \"c(r,z=5.263,t=200)\"; plot \"y_5_t_200.dat\" using 2:3 with linespoints pt 4 lw 2 lt rgb \"violet\" title \"Schiesser\",  \"marc_y_5_t_200.dat\" using 2:3 with linespoints pt 7 lw 2 lt rgb \"cyan\" title \"Marc\"");
     GNUplot plot10;
     plot10(yrange1+"set xlabel \"r\"; set ylabel \"T(r,z=5.263,t=200)\"; plot \"y_5_t_200.dat\" using 2:4 with linespoints pt 4 lw 2 lt rgb \"violet\" title \"Schiesser\", \"marc_y_5_t_200.dat\" using 2:4 with linespoints pt 7 lw 2 lt rgb \"cyan\" title \"Marc\"");

     
     //! z = 35, t = 200
     RIM.clear(); //prepare for next read in
     RIM.read_from_file("CDR_CYL_SCHIESSER_2D_4.gnu");
     RIM.fixed_x_2D(36.8421052632, "y_35_t_200.dat");

     MIR.clear(); //prepare for next read in
     MIR.read_from_file("MARC_CDR_CYL_2D_4.gnu");
     MIR.fixed_x_2D(36.8421052632, "marc_y_35_t_200.dat");
     
     GNUplot plot11;
     plot11(yrange2+"set xlabel \"r\"; set ylabel \"c(r,z=36.842.263,t=200)\"; plot \"y_35_t_200.dat\" using 2:3 with linespoints pt 4 lw 2 lt rgb \"violet\" title \"Schiesser\", \"marc_y_35_t_200.dat\" using 2:3 with linespoints pt 7 lw 2 lt rgb \"cyan\" title \"Marc\"");
     GNUplot plot12;
     plot12(yrange1+"set xlabel \"r\"; set ylabel \"T(r,z=36.842,t=200)\"; plot \"y_35_t_200.dat\" using 2:4 with linespoints pt 4 lw 2 lt rgb \"violet\" title \"Schiesser\", \"marc_y_35_t_200.dat\" using 2:4 with linespoints pt 7 lw 2 lt rgb \"cyan\" title \"Marc\"");

      //! z = 70, t = 200
     RIM.clear(); //prepare for next read in
     RIM.read_from_file("CDR_CYL_SCHIESSER_2D_4.gnu");
     RIM.fixed_x_2D(68.4210526316, "y_70_t_200.dat");

     MIR.clear(); //prepare for next read in
     MIR.read_from_file("MARC_CDR_CYL_2D_4.gnu");
     MIR.fixed_x_2D(68.4210526316, "marc_y_70_t_200.dat");
     
     GNUplot plot13;
     plot13(yrange2+"set xlabel \"r\"; set ylabel \"c(r,z=68.421.263,t=200)\"; plot \"y_70_t_200.dat\" using 2:3 with linespoints pt 4 lw 2 lt rgb \"violet\" title \"Schiesser\", \"marc_y_70_t_200.dat\" using 2:3 with linespoints pt 7 lw 2 lt rgb \"cyan\" title \"Marc\"");
     GNUplot plot14;
     plot14(yrange1+"set xlabel \"r\"; set ylabel \"T(r,z=68.421,t=200)\"; plot \"y_70_t_200.dat\" using 2:4 with linespoints pt 4 lw 2 lt rgb \"violet\" title \"Schiesser\", \"marc_y_70_t_200.dat\" using 2:4 with linespoints pt 7 lw 2 lt rgb \"cyan\" title \"Marc\"");

      //! z = 100, t = 200
     RIM.clear(); //prepare for next read in
     RIM.read_from_file("CDR_CYL_SCHIESSER_2D_4.gnu");
     RIM.fixed_x_2D(100., "y_100_t_200.dat");

     MIR.clear(); //prepare for next read in
     MIR.read_from_file("MARC_CDR_CYL_2D_4.gnu");
     MIR.fixed_x_2D(100., "marc_y_100_t_200.dat");
     
     GNUplot plot15;
     plot15(yrange2+"set xlabel \"r\"; set ylabel \"c(r,z=100,t=200)\"; plot \"y_100_t_200.dat\" using 2:3 with linespoints pt 4 lw 2 lt rgb \"violet\" title \"Schiesser\", \"marc_y_100_t_200.dat\" using 2:3 with linespoints pt 7 lw 2 lt rgb \"cyan\" title \"Marc\"");
     GNUplot plot16;
     plot16(yrange1+"set xlabel \"r\"; set ylabel \"T(r,z=100,t=200)\"; plot \"y_100_t_200.dat\" using 2:4 with linespoints pt 4 lw 2 lt rgb \"violet\" title \"Schiesser\", \"marc_y_100_t_200.dat\" using 2:4 with linespoints pt 7 lw 2 lt rgb \"cyan\" title \"Marc\"");
     
   }
 }
  
 else if (example == 92){
   typedef MarcCDRCylindrical2D<CppAD::AD<double> > Functor1;
   typedef HodgkinHuxley<double> Functor2;

   Functor1 fun1(2*20*7);
   Functor2 fun2(4);
   
   MyVec<double> u1(2*20*7),
     u2(4);

   cout << "Is MarcCDRCylindrical2D a MOL functor?: "<<IsMOLFunctor<Functor1>::Value<< endl;
   cout << "Is HodgkinHuxley a MOL functor?: "<< IsMOLFunctor<Functor2>::Value << endl;
   BoundaryFromFunctor<Functor1,IsMOLFunctor<Functor1>::Value>::set_boundary(u1,fun1);
   BoundaryFromFunctor<Functor2,IsMOLFunctor<Functor2>::Value>::set_boundary(u2,fun2);
   
   /*
   double delta_t = 0.8*(1./( 200./2.e-07 + 300./1.e-07 + 342.59* sqrt(1./Sqr(2.e-07) + 1./Sqr(1.e-07)) + 2*(4./3.*1.e-05*(1.e-05/0.71)/0.925)*(1./Sqr(2.e-07) + 1./Sqr(1.e-07))) );
   cout <<  "delta_t = "<< delta_t << endl;
   */
}
 else{
   ADONIS_ERROR(MainProgramError, "Example '"<<example<<"' hasn't been declared for \""<<argv[0]<<"\" yet. \n");
 }

  
 

  return 0;
}
