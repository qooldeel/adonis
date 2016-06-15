//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "../want_2_use.h" //MUST PUT HERE IF U WANNA USE AD      !!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "navierstokesviamethodoflines.hh"
#include "lightnse.hh"
#include "1dlightnse.hh"

#include "mol2D.hh"

#include "../moleculartransport/cfdadds.hh"

#include "fulldotomega.hh"
#include "reconstructspecies.hh"
#include "../ode/examples/automaticallygeneratedsourceterms/o3.hh" //test

#include "nseh2gri.hh"
#include "nseh2conaire.hh"  //primary reference fuel mechanism
#include "o3nse.hh" //ozone mechanism with full navier stokes

#include "../ode/examples/automaticallygeneratedsourceterms/grimech30exerpt.hh"

#include "../ode/examples/cstr.hh"
#include "../ode/solvers.hh"

//REDUCED FUNCTORS
#include "../ode/examples/automaticallygeneratedsourceterms/redo3.hh"
#include "1dREDlightnse.hh"
#include "REDlightnse.hh"
#include "REDnseh2conaire.hh"

#include "../io/readinmatrix.hh"

using namespace Adonis;
using namespace ExprTmpl;
using namespace std;
 


template<class T>
class SpecialFunctional{
public:
  typedef ExprTmpl::MyVec<T> VType;
  
  SpecialFunctional():cost_(0.),x_(4),q_(2),Q_(4),R_(2), tempx_(4), tempq_(2), xS_(4),qS_(2),x1_(4),q1_(2), wS_(6), A_(6), tempw_(6),w1_(6),grad_(6), hess_(6*6){
    Q_ <<= 0.2, 1.0, 0.5, 0.2;
    R_ <<= 0.5, 5.e-07;
    xS_ <<= 2.1402, 1.0903, 114.19, 112.91;
    qS_ <<= 14.19, -1113.5;
    concatenate(wS_,xS_,qS_);
    concatenate(A_,Q_,R_);

    for(int k = 0; k < 6; ++k)
      hess_[RowMajor::offset(k,k,6)] = A_[k];

  }
  
 
  //! compute L(x,q) := (x-x_S)^TQ(x-x_S) + (q-q_S)^TR(q-q_S)
   T operator()(const VType& x, const VType& q){
     tempx_ = x - xS_;
     x1_ = Q_*tempx_;

     tempq_ = q - qS_;
     q1_ = R_*tempq_;

     return (dot(tempx_,x1_) + dot(tempq_,q1_));
   }

 

  T compute_L(const VType& w){
    tempw_ = w - wS_;
    w1_ = A_*tempw_;
    return dot(tempw_,w1_);
  }

  const T& cost() const {return cost_;}

  //Riemann integral
  T one_segment(const VType& w, const VType& wprev, const T& h){
    return h*compute_L(w); //0.5*h*((*this).compute_L(w)+ (*this).compute_L(wprev));
  }

 

  //use trapezoidal rule to integrate cost functional
  T operator()(const VType& w, const VType& wprev, const T& h){
    cost_ += one_segment(w,wprev,h);
    return cost_;
  }

   VType& gradient(const T& h, const VType& w){
     grad_ = 2.*h*A_*(w-wS_);  //h*A*(w-wS) when trapezoidal rule is used
    return grad_;
  }


private:
  T cost_;
  VType x_, q_, Q_, R_, tempx_, tempq_, xS_, qS_, x1_, q1_, wS_, A_,tempw_,w1_;
  VType grad_, hess_;

};

template<class V>
class ProjectOntoBounds{
public:
  typedef V VType;
  typedef typename V::value_type value_type;

  ProjectOntoBounds(const V& l, const V& u):low_(l),up_(u), proj_(l.size()){adonis_assert(l.size() == u.size());}

  V& operator()(const V& w){
    proj_ = Max(low_,Min(w,up_));
    return proj_;
  }
  
  //! assume w is already projected
  V& grad_projection(const V& w, const V& grad){
    adonis_assert(w.size() == grad.size());
    for(size_t k = 0; k < w.size(); ++k){
      if((low_[k] < w[k]) && (w[k] < up_[k]))
	proj_[k] = grad[k];         
      if(is_zero(Abs(w[k]-low_[k]), 1e+02)) //factor to multiply machine prec.
	proj_[k] = min(grad[k],0.);
      if(is_zero(Abs(w[k]-up_[k]), 1e+02))
	proj_[k] = max(grad[k],0.);
    }
    return proj_;
  }

private:
  const V& low_, up_;
  V proj_;
};




int main(int argc, char** argv){

  if(argc != 2) //argv[0] = Program name, argv[1] = number of executable
    ADONIS_ERROR(MainProgramError, "\""<< argv[0] <<"\" takes exactely 1 argument, pal.");

  int example = atoi(argv[1]);

  //!START EXAMPLES HERE
  if(example == 1){
     cout << "USE METHOD OF LINES FOR 1D RECTANGULAR DOMAIN" << endl;
   
#if FULL_MODEL == 1
     ParameterData PDf;
     PDf.read_from_file(ChooseFDSetting::FILE2READIN);
  
     //================== TODO ============================
     //=== change model type here
     typedef OneDIMLightNSE<double> ModelType;
     //====================================================
 
     int nptf = PDf.get_datum<int>("Nx"),
       totf = ModelType::nprim*nptf;
     cout << "totf = "<<totf << endl;
  
     //……………………………………………… ODE …………………………………………………………………………
     //initial value 
     MyVec<double> u0f(totf);
  
 
     int nxf = PDf.get_datum<int>("Nx");

     double af = PDf.get_datum<double>("a"),
       bf = PDf.get_datum<double>("b"),
       hxf = (bf-af)/(nxf-1),
       T_inletf = PDf.get_datum<double>("T_min"),
       T_burnerf = PDf.get_datum<double>("T_ignition");
  
     cout << "T_burnerf = "<<T_burnerf <<endl;

     cout << "hxf = " << hxf << endl; 

  
     enum{Oxyf=1,Oxy2f=2,Oxy3f=3};

     for(int i = 0; i < nxf; ++i){
   
       u0f[i] = T_inletf;

       if(i == 0)
	 u0f[i] = T_burnerf; //at left boundary, burner temperature

       u0f[i+Oxyf*nptf] = 0.0;
       u0f[i+Oxy2f*nptf] = 0.8;
       u0f[i+Oxy3f*nptf] = 0.2;
   
     }

     //cout << u0 << endl;
     const bool SPARSEF = true;
     const bool SCALESYSF = true; //try preventing ill-conditioning of matrix

     MOL<1,double,OneDIMLightNSE,SPARSEF,SCALESYSF> molf;
     molf.solve(u0f,"2Dsettings.dat","O3_1D_FULL"); //write in files 
                                                    //"O3_1D_FULL*"
  #endif


#if SWITCH_ON_REDUCTION == 1
     cout << "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"<<endl;
     cout << "RR        REDUCED SYSTEM               RR"<<endl;
     cout << "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"<<endl;
     ParameterData PD;
     PD.read_from_file(ChooseFDSetting::FILE2READIN);
  
     //================== TODO ============================
     //=== change model type here
     typedef ReducedOneDIMLightNSE<double> RedModelType; //===RED
     //====================================================
 
     int npt = PD.get_datum<int>("Nx"),
       tot = RedModelType::nprim*npt;
     cout << "tot = "<<tot << endl;
  
     //……………………………………………… ODE …………………………………………………………………………
     //initial value 
     MyVec<double> u0(tot);
  
 
     int nx = PD.get_datum<int>("Nx");

     double a = PD.get_datum<double>("a"),
       b = PD.get_datum<double>("b"),
       hx = (b-a)/(nx-1),
       T_inlet = PD.get_datum<double>("T_min"),
       T_burner = PD.get_datum<double>("T_ignition");
  
     cout << "T_burner = "<<T_burner <<endl;

     cout << "hx = " << hx << endl; 

  
     enum{Oxy=1};  //===RED

     for(int i = 0; i < nx; ++i){
   
       u0[i] = T_inlet;

       if(i == 0)
	 u0[i] = T_burner; //at left boundary, burner temperature

       u0[i+Oxy*npt] = 0.0;  //!only O is considered
   
     }

     cout << "u0 = "<< u0 << endl;
     cout <<  "u0.size() = "<< u0.size() << endl;

     //cout << u0 << endl;
     const bool SPARSE = true;
     const bool SCALESYS = true; //true; //try preventing ill-conditioning of matrix

     MOL<1,double,ReducedOneDIMLightNSE,SPARSE,SCALESYS> mol;
     mol.solve(u0,"2Dsettings.dat","O3_1D_RED"); //ensure that data written in
                                                 //file "O3_1D_RED"
    
     //! test only
     //typedef ReducedOneDIMLightNSE<CppAD::AD<double> > Model;
     //Model redfun(2*npt); //1(temp)+1(#RPVs)
     // MyVec<CppAD::AD<double> > ev(u0.size());
  
     // //Left
     // u0[0] = 1000;
     // u0[0+nx*1]= 0;
     // u0[0+nx*2]= 0.8;
     // u0[0+nx*3]= 0.2;

     // //right
     // u0[nx-1] = u0[nx-2];
     // for(int k = 0; k < 3; ++k)
     //   u0[nx-1 + (1+k)*nx] = u0[nx-2 + (1+k)*nx];
     

     // // for(size_t i = 0; i < u0.size(); ++i)
     // //   ev[i] = u0[i];
     // // MyVec<CppAD::AD<double> > feval = redfun(ev);
     
     // JacS<double,ReducedOneDIMLightNSE,ExprTmpl::MyVec> nablachem;
     // nablachem.set_with_init(2*npt,2*npt,u0);
     // MyVec<double> Jac = nablachem.jacobian(u0);
     //cout << "Jac = "<< Jac << endl;
     
     
       

#endif

  }
  
  //·········································································
  //·················· 2DIM ·················································
  //·········································································
  else if(example == 2){
    cout << "USE METHOD OF LINES FOR 2D RECTANGULAR DOMAIN" << endl;
   
#if FULL_MODEL == 1
    ParameterData PD;
    PD.read_from_file(ChooseFDSetting::FILE2READIN);
  
    //================== TODO ============================
    //=== change model type here
    typedef LightNSE<double> ModelType;
    //====================================================
 
    int npt = PD.get_datum<int>("Nx")*PD.get_datum<int>("Ny"),
      tot = ModelType::nprim*npt;
    cout << "tot = "<<tot << endl;
  
    //……………………………………………… ODE …………………………………………………………………………
    //initial value 
    MyVec<double> u0(tot);
  
 
    int nx = PD.get_datum<int>("Nx"),
      ny = PD.get_datum<int>("Ny");
    double a = PD.get_datum<double>("a"),
      b = PD.get_datum<double>("b"),
      c = PD.get_datum<double>("c"),
      diam = PD.get_datum<double>("d") - c,
      halfdiam = diam/2.,
      hy = diam/(ny-1),
      hx = (b-a)/(nx-1),
      T_inlet = PD.get_datum<double>("T_min");
  
    cout << "T_inlet = "<<T_inlet <<endl;

    cout << "hx = " << hx << "   hy = "<< hy << "  halfdiam = "<< halfdiam<< endl;

    BoundaryTemperatureDistribution<double> BTD(PD.get_datum<double>("T_min"),PD.get_datum<double>("T_ignition"),PD.get_datum<double>("hottest_x"));

   
    enum{Oxy=1,Oxy2=2,Oxy3=3};

    for(int i = 0; i < nx; ++i){
      double x = a + i*hx; 
      for(int j = 0; j < ny; ++j){
	int ij = i + nx*j;

	u0[ij] = T_inlet;

        if(j == ny-1){
	  u0[ij] = BTD(x); //parts of upper bdy = ignition zone
	}
	

	u0[ij+Oxy*npt] = 0.0;
	u0[ij+Oxy2*npt] = 0.8;
	u0[ij+Oxy3*npt] = 0.2;
      }
    }

    
    const bool SPARSE = true;
    const bool SCALESYS = true; //try preventing ill-conditioning of matrix

    //MOL2D<double,ReactiveNS2D,SPARSE,SCALESYS> mol;
    MOL<2,double,LightNSE,SPARSE,SCALESYS> mol;
    mol.solve(u0);
#endif  


#if SWITCH_ON_REDUCTION == 1
     cout << "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"<<endl;
     cout << "RR     2D REDUCED SYSTEM               RR"<<endl;
     cout << "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"<<endl;
     ParameterData PDr;
     PDr.read_from_file(ChooseFDSetting::FILE2READIN);
     
     //=================================================
     typedef ReducedLightNSE<double> ReducedModelType;
     //=================================================

     int nptr = PDr.get_datum<int>("Nx")*PDr.get_datum<int>("Ny"),
      totr = ReducedModelType::nprim*nptr;
     cout << "tot = "<<totr << endl;
  
    //……………………………………………… ODE …………………………………………………………………………
    //initial value 
    MyVec<double> u0r(totr);
    
    int nxr = PDr.get_datum<int>("Nx"),
      nyr = PDr.get_datum<int>("Ny");
    double ar = PDr.get_datum<double>("a"),
      br = PDr.get_datum<double>("b"),
      cr = PDr.get_datum<double>("c"),
      diamr = PDr.get_datum<double>("d") - cr,
      halfdiamr = diamr/2.,
      hyr = diamr/(nyr-1),
      hxr = (br-ar)/(nxr-1),
      T_inletr = PDr.get_datum<double>("T_min");
  
    cout << "T_inlet = "<<T_inletr <<endl;

    cout << "hx = " << hxr << "   hy = "<< hyr << "  halfdiam = "<< halfdiamr<< endl;

    BoundaryTemperatureDistribution<double> BTDr(PDr.get_datum<double>("T_min"),PDr.get_datum<double>("T_ignition"),PDr.get_datum<double>("hottest_x"));

   
    enum{Oxyr=1};

    for(int i = 0; i < nxr; ++i){
      double x = ar + i*hxr; 
      for(int j = 0; j < nyr; ++j){
	int ij = i + nxr*j;

	u0r[ij] = T_inletr;

        if(j == nyr-1){
	  u0r[ij] = BTDr(x); //parts of upper bdy = ignition zone
	}
	

	u0r[ij+Oxyr*nptr] = 0.0;
      }
    }

    
    const bool SPARSEr = true;
    const bool SCALESYSr = true; //try preventing ill-conditioning of matrix

    MOL<2,double,ReducedLightNSE,SPARSEr,SCALESYSr> molred;
    molred.solve(u0r,"2Dsettings.dat","NSE_O3_RED_");


#endif


  }

  else if (example == 3){
    enum{ozone=3};
    typedef ReconstructSpecies<double,3,1,ozone> ChemType;
    ChemType RS;
    typedef ChemType::DataType DataType; 
    MyVec<double> zM0(3+1);

    zM0 <<= 0.,  //O 
      0.8,       //O2
      0.2,       //O3
      1000.;     //T
    
    MyVec<double> C(3*1), b(1);
    C <<= 1.,1.,1.;
    b <<= 1.;

    
    RS.initialize(zM0,5e-06,1.e-08,'E');


    MyVec<double> red(1+1);
    red <<= 0.005,   //1 RPV
      1000;          // TEMP

    RS.assign(red,zM0);
    RS.enforce_linearized_constraints_wrt_chemistry(C,b);
    MyVec<double> low(3+1),up(3+1);
    low <<= 0,0,0,DataType::temperature_bounds()[0];
    up <<= 1.,1.,1.,DataType::temperature_bounds()[1];
    cout << "low = "<< low << endl << "up = "<< up << endl;
   

    RS.project_onto_bounds(low,up);

    MyVec<double> Ynew = RS.get_z();
    cout << "Ynew = "<< Ynew << endl;
    
    //cout << "Ynew = "<< Ynew << endl;
    //now convert it into concentrations
    // MyVec<double> Ct(RS.size());
    // double wbar = 0;
    // for(int k = 0; k < 3; ++k)  //in full chemistry here, coz of full concentr.
    //   wbar += Ynew[k]/DataType::molar_masses()[k];
    // wbar = 1./wbar;

    // double rho = (101325* wbar)/(PhysicalConstants<double>::Rgas*Ynew[3]);
    // for(int k = 0; k < 3; ++k)
    //   Ct[k] = rho*Ynew[k]/DataType::molar_masses()[k];

    // Ct[3] = Ynew[3]; //temperature
    //cout<< "concentrations + temperature = "<< Ct << endl;
    //RS.time_step(Ct);
    

    RS.time_step(Ynew,101325);
    cout << "which time stepping method used: "<< RS.which_time_stepper() << endl;
    

    cout << "RPV index = "<< RS.rpv_index() << endl;
    cout << "unrep index = "<< RS.unrep_index() << endl;

  }
  else if (example == 4){
    
    MyVec<double> U0(4);
    U0 <<= 0.00333333, 0.798333, 0.198333, 1000; //0.,0.8,0.2,1000;

   
    DotOmegaWithTemperatureFull<double,3> Fun;
    cout << "Fun(U0) = " << Fun(U0) << endl;

    O3Decomposition<double> o3test(4);

    cout << "O3 test = "<< o3test(U0) << endl;

      // double t0 = 0.,
      // 	tend = 1.e-03,
      // 	kn = 5.e-06;
    // JacobianOfDotOmegaWithTemperatureFull<double,3> Jac;
    // Jac.initialize(U0);
    // //BACKWARD EULER -- easy implementation
    
    // double time(t0);
    // MyVec<double> b(4);
    // MyVec<double> Uprev(U0), U(U0), IterMatrix(4*4);
    // int count(0);
    
    // ofstream of("O3_SOURCE.dat");
    // while(time < tend){
    //   cout << count << ". )        t = "<<time << "   dt = "<< kn << endl; 
    //   of << time << "  ";
    //   for(int k = 0; k < 4; ++k)
    // 	of << U[k] << "  ";
    //   of << endl;

    //   U = Uprev;   //start guess
    //   int check(0);
    //   for(size_t l = 1; l <= 5; ++l){ //Newton iter
    // 	check = l;
    // 	b = -(U - Uprev - kn*Fun(U));
    // 	IterMatrix = -kn*Jac.jacobian(U);
    // 	update_diagonal<AddBasicElements>(IterMatrix,4,1.);
    // 	good_square_solve(IterMatrix,4,b,1);
    // 	U += b;
    // 	if(Norm<'2',double>::norm(b) <= 1.e-09)
    // 	  break;
    // 	if(check == 5)
    // 	  ADONIS_INFO(Information,"Newton iteration did not converge within 5 iterations");
    //   } //end NEWTON
    //   Uprev = U;

    //   time += kn;
    //   count++;
    // }
    // cout << "U(t="<<tend<<") = "<< U << endl;
  
    // GNUplot show;
    // show("plot 'O3_SOURCE.dat' using 1:2 with lines lw 2 title \"Y_{O}\",'O3_SOURCE.dat' using 1:3 with lines lw 2 title \"Y_{O2}\", 'O3_SOURCE.dat' using 1:4 with lines lw 2 title \"Y_{O3}\"");
    // GNUplot showtemp;
    // showtemp("plot 'O3_SOURCE.dat' using 1:5 with lines lw 2 title \"T\"");
  }


  else if (example == 5){
    cout << "PRIMARY REFERENCE FUEL H2 MECHANISM DUE TO CONAIRE, CURRAN, SIMMIE, PITZ and WESTBROOK:"<<endl;

    ParameterData PD;
    PD.read_from_file("data/conaireh2.dat");
    //================== TODO ============================
    //=== change model type here
    typedef ConaireH2ReactiveNS2D<double> ModelType;
    typedef ModelType::DataType DataType;
    typedef ModelType::InVeloType InletVelocityType;


    //=== Settings for MOL
    const bool SPARSE = true;
    //! caution: scaling in conjunction with Rosenbrock method can fail!!! 
    const bool SCALESYS = false; //false; // normally, UMFPACK internally scales as well //true; //true: try preventing ill-conditioning of matrix
    const bool PATTERNCREATIONVIAAD = true; //true: AD, false: FD
    const bool REPAIR = true;
    const char REPAIRTYPE = 'f';
    const bool PRINTPRESSURE = true;

    string cpuInfoOutFile = "H2mechanism_CPU_Info.cpu";

    int excessSpec = DataType::N2;
    //======================================================
 
    int npt = PD.get_datum<int>("Nx")*PD.get_datum<int>("Ny"),
      tot = ModelType::nprim*npt;
    cout << "tot = "<<tot << endl;

  
    //……………………………………………… ODE …………………………………………………………………………
    //initial value 
    MyVec<double> u0(tot);
    
 
    int nx = PD.get_datum<int>("Nx"),
      ny = PD.get_datum<int>("Ny");
    double a = PD.get_datum<double>("a"),
       b = PD.get_datum<double>("b"),
      length = b-a,
      hx = length/
#ifndef GHOST_POINTS_INCLUDED
      (nx-1)
#else
      (nx-3)
#endif
      ,
      c = PD.get_datum<double>("c"),
      d = PD.get_datum<double>("d"),
      diam = d - c,
      hy = (diam)/
#ifndef GHOST_POINTS_INCLUDED
      (ny-1)
#else
      (ny-3)
#endif 
      ,
      p0 = PD.get_datum<double>("pconst"),
      yhottest = PD.get_datum<double>("hottest_y"),
      Tignition = PD.get_datum<double>("T_ignition"),
      Tmin = PD.get_datum<double>("T_min"),
      Tin = Tmin,
      halfdiam = 0.5*diam,
      Lfrac = PD.get_datum<double>("channel_length_frac"),
      v1_in = PD.get_datum<double>("v1_in"),
      v2_in = PD.get_datum<double>("v2_in"),
      v1_wall = PD.get_datum<double>("v1_wall"),
      v2_wall = PD.get_datum<double>("v2_wall"),
      burnerLength = length/Lfrac; //PD.get_datum<double>("burner_length");
    cout << "hx = "<< hx << "  hy = "<< hy << "    yhottest = "<< yhottest << "  burner length = "<< burnerLength << endl;
    
    // BoundaryTemperatureDistribution<double> BTD(Tin,Tignition,yhottest);

    double burnerSpikeLength_y = PD.get_datum<double>("y_burner_spike_factor");
    cout << "burnerSpikeLength_y = " << burnerSpikeLength_y << endl;


    std::string testOut;

#ifdef GHOST_POINTS_INCLUDED
      ModelType Mt(u0.size());
#endif
      
      int fromLastSlideOn = PD.get_datum<int>("fromLastStepOn");

      if(fromLastSlideOn == 0){ //start from t = 0
	//either constant or parabolic
	MyVec<double> veloprofile(
#ifndef GHOST_POINTS_INCLUDED
				  ny
#else
				  ny-2
#endif
				  );
	cout << "veloprofile.size() = "<< veloprofile.size() << endl;
    
	double r = -halfdiam;  //R = d/2, r = -R

	for(int j = 0; j < (int)veloprofile.size(); ++j)
	  {
	
	    veloprofile[j] = InletVelocityType::velocity(r,halfdiam,v1_in); 
	    r += hy;
	    //cout << "veloprofile["<<j<<"] = "<< veloprofile[j] << endl;
	  }

 
	//O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR
	double air[10] = {PD.get_datum<double>("O"),    
			  PD.get_datum<double>("O2"),  
			  PD.get_datum<double>("H"),  
			  PD.get_datum<double>("OH"),   
			  PD.get_datum<double>("H2"),    
			  PD.get_datum<double>("HO2"),   
			  PD.get_datum<double>("H2O2"),     
			  PD.get_datum<double>("H2O"),     
			  PD.get_datum<double>("N2"),
			  PD.get_datum<double>("AR")
	};  

	//print_all(air,air+DataType::nspec,5,false,"initial vector");
	MyVec<double> Icomp(air,air+10);
	cout << "INITIAL COMPOSITION: "<< Icomp<< endl;

	double wbar(0.);
	for(int k = 0; k < DataType::nspec; ++k)
	  wbar += Icomp[k]/DataType::molar_masses()[k]; //initially filled with air
	wbar = 1./wbar;

#ifdef USE_OPENMP
	cout  <<endl;
	FancyMessages().nice_output("OpenMP will be used", 36); //33 = yello
	cout << endl;
#endif


	//-----------------------Reynolds number---------------------------------
	MixtureAveragedViscosity<DataType,false> Mav(DataType::transport());
	StaticArray<double,DataType::nspec> Xfrac;
	for(int k = 0; k < DataType::nspec; ++k){
	  Xfrac[k] = (Icomp[k]*wbar)/DataType::molar_masses()[k];
	}
    
	double mu = Mav.compute_mixture_averaged_viscosity(Tin,Xfrac);
	double rho0 = (p0*wbar)/(PhysicalConstants<double>::Rgas*Tin);
	cout << "mixture viscosity = "<< mu <<  "   char. mixture density = "<< rho0 << endl;
	double Re = Reynolds_number(rho0,v1_in,diam,mu);
	cout << endl << "REYNOLDS NUMBER = "<< Re << endl;
	//-----------------------------------------------------------------------
    
	//-------------------------Mach number-----------------------------------
	double Ma = Mach_number<DataType::ENCOD>(Tin,v1_in,Icomp);
	cout  << "MACH NUMBER = "<< Ma << endl; 
	//-----------------------------------------------------------------------

	cout << endl << "inlet velocity type   = "<< InletVelocityType::Value << endl;

	typedef BoundaryTemperatureDistribution<double,FDMSettings::WallTemperatureType> BoundaryWallTemperatureType;
	cout  << "wall temperature type = "<< ascii(BoundaryWallTemperatureType::Value) << endl;
  

	double x_ignition = PD.get_datum<double>("x_ignition");//(b-a)/PD.get_datum<double>("channel_length_frac");
	cout << "x_ignition = " << x_ignition << endl;


	BoundaryWallTemperatureType BWT(Tmin,Tignition,x_ignition,(((BoundaryWallTemperatureType::Value == 's') || (BoundaryWallTemperatureType::Value == 'S')) ? PD.get_datum<double>("cubic_spline_from_ignition")*hx : PD.get_datum<double>("plateau_width")*hx),(((BoundaryWallTemperatureType::Value == 'a') || (BoundaryWallTemperatureType::Value == 'A')) ? 2*hx : length));
	MyVec<double> walltempprofile(
#ifndef GHOST_POINTS_INCLUDED
				      nx
#else
				      nx-2
#endif
				      );
	double xdir(0.);


	for(int i = 0; i < (int)walltempprofile.size(); ++i)
	  {
	    xdir = i*hx;
	    walltempprofile[i] = BWT(xdir);
	    cout << "T_wall("<<xdir<<") = "<< walltempprofile[i] << endl;
	  }


	cout << "wbar = "<< wbar << endl;
	double troom = Tin;
	double rho = (p0*wbar)/(PhysicalConstants<double>::Rgas*troom);
	cout << "rho = "<< rho << endl;

	//initally domain (Omega with boundary) is filled with air at room temperature
	
	double temper = troom;
	double v1(v1_in), v2(v2_in);

#ifndef IGNITION_AT_HOT_WALLS
	double xdist(0.);   //ignition at inlet
#endif
	double y(0.);
	cout << "y-direction start: "<< y << endl;
	for(int i = 0; i < nx; ++i){
	  y = 0; //reset y-direction
     
#ifndef IGNITION_AT_HOT_WALLS
	  //xdist += hx;
	  xdist = a + i*hx; //a = 0, 
#endif
	  for(int j = 0; j < ny; ++j){
	    y = j*hy; 
	    int ij = i + nx*j;

#ifndef GHOST_POINTS_INCLUDED	    
	    v1 = veloprofile[j];
#else //!ghost points used
	    if((j >= 1) && (j <= ny-2)){
	      v1 = veloprofile[j-1];
	    }
#endif
	    v2 = v2_in;
	 
#ifdef IGNITION_AT_HOT_WALLS
#ifndef GHOST_POINTS_INCLUDED
	    if((j == 0) || (j == ny-1) )
#else
	      if((j==1) || (j == ny-2) )  
#endif
		{ //DOWN or UP boundary (wall)
		  v1 = v1_wall;  //no slip 
		  v2 = v2_wall;  //no slip
#ifndef GHOST_POINTS_INCLUDED	    
		  temper =  walltempprofile[i];
#else //!ghost points used
		  if((i >= 1) && (i <= nx-2)){
		    temper = walltempprofile[i-1];
		  }
#endif
		}
	      else{
		temper = troom;
	      }
	

	    //!=============================================================

#else //ignition at left side
	    if(
#ifndef GHOST_POINTS_INCLUDED
	       (i == 0)
#else
	       (i == 1)
#endif
	       && (is_contained(yhottest,y-hy,y+hy))){
	      temper = Tignition;
	    } 
	    else
	      temper = troom;
#endif 
	 
	    u0[ij] =  (p0*wbar)/(PhysicalConstants<double>::Rgas*temper);

	    //velocities
	    u0[ij+1*npt] = 
#ifdef NONCONSERVATIVE_FORM
	      v1;  //no slip at wall
#else
	    v1*u0[ij];
#endif

	    u0[ij+2*npt] = 
#ifdef NONCONSERVATIVE_FORM
	      v2;
#else
	    v2*u0[ij];
#endif

	    u0[ij+3*npt] = 
#ifdef NONCONSERVATIVE_FORM
	      temper;
#else
	    temper*u0[ij];
#endif
	 

	    for(int k = 0; k < DataType::nspec; ++k){
	      u0[ij+(4+k)*npt] = 
#ifdef NONCONSERVATIVE_FORM
		Icomp[k]; // mixture everywhere
#else
	      Icomp[k]*u0[ij];
#endif
	    }

	  } //end domain
	}

#ifdef GHOST_POINTS_INCLUDED
       Mt.set_ghost_points(u0); //Mt declared above
#endif

	cout << endl<< "REYNOLDS NUMBER = "<< Re << endl;
	cout  << "MACH NUMBER = "<< Ma << endl;

	testOut = "NormalOut";
      } //fromLastSlide == 0
      
      else if(fromLastSlideOn > 0){ //start from available record
      cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"<<endl;
      cout << "P         PRINT FROM LAST SLIDE " << fromLastSlideOn << " ON    P" << endl;
      cout << "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"<<endl;
      
      cout << "REARRANGE u0:"<<endl;
      ReadInMatrix<double> RIM;
      RIM.read_from_file(PD.get_datum<string>("lastSlide"));
      bool withPressure(true);
      int numberOfQuantities =  ModelType::nprim;
      if(PRINTPRESSURE){
	withPressure=false;  //exclude pressure
	numberOfQuantities =  ModelType::nprim + 1; //total quants with p
      }
	
      RIM.assign_2_2D_MOL_slide(u0,nx,ny,numberOfQuantities,withPressure);
      cout << "offs = "<< RIM.get_offset() << endl;

#ifdef GHOST_POINTS_INCLUDED
       Mt.set_ghost_points(u0); //Mt declared above

      //! DEBUGGING
      // for(size_t i = 0; i < u0.size(); ++i){
      // 	if(!is_well_defined_value(u0[i])){
      // 	  ADONIS_ERROR(ValueError, "Unphysical value at index = "<< i << ". Correct initialization applied??");
      // 	}
      
#endif

       testOut = "CheckOut";
   
         
      }
      else{
	ADONIS_ERROR(IndexError,"fromLastSlideOn = "<<fromLastSlideOn<<" must be POSITIVE!");
      }
      

      //! JUST FOR DEBUGGING
      // MyVec<int> ixtpe;
      // PrintSolution PS;
      // PrintMOLData<2,double,PRINTPRESSURE> PMD;
      // PMD.init("data/conaireh2.dat",testOut,PS,12);
      // PMD.write_2_file(5.e-12,700,1.3440234375e-08,u0,testOut,PS,ixtpe,DataType::molar_masses(),PD.get_datum<double>("elapsed_time_so_far"),PD.get_datum<double>("number_lasttimestep"),0);
     
      
       //MUST BE commented in late on
    cout << "Excess species: "<< DataType::species_names()[excessSpec] << endl; 
    //MOL2D<double,ReactiveNS2D,SPARSE,SCALESYS> mol;
    MOL<2,double,ConaireH2ReactiveNS2D,SPARSE,SCALESYS,PATTERNCREATIONVIAAD,REPAIR,REPAIRTYPE,PRINTPRESSURE> mol;
    mol.solve(u0,"data/conaireh2.dat","",cpuInfoOutFile,DataType::molar_masses(),fromLastSlideOn,excessSpec);
      



  } //end example

  //=========================================================================
  //=========================================================================
  //=========================================================================


  
  else if(example == 6){
    FancyMessages().nice_output("REDUCED H2 CONAIRE",35);

    ParameterData PD;
    PD.read_from_file("data/conaireh2.dat");
      //================== TODO ============================
    //=== change model type here
    typedef ConaireH2ReactiveNS2D<double> ModelType;
    typedef ModelType::DataType DataType;
    typedef ModelType::InVeloType InletVelocityType;
    //======================================================
 
    int npt = PD.get_datum<int>("Nx")*PD.get_datum<int>("Ny"),
      tot = ModelType::nprim*npt;
    cout << "tot = "<<tot << endl;

  
    //……………………………………………… ODE …………………………………………………………………………
    //initial value 
    MyVec<double> u0(tot);
  
 
    int nx = PD.get_datum<int>("Nx"),
      ny = PD.get_datum<int>("Ny");
    double a = PD.get_datum<double>("a"),
       b = PD.get_datum<double>("b"),
      hx = (b-a)/(nx-1),
      c = PD.get_datum<double>("c"),
      d = PD.get_datum<double>("d"),
      diam = d - c,
      hy = (diam)/(ny-1),
      p0 = PD.get_datum<double>("pconst"),
      yhottest = PD.get_datum<double>("hottest_y"),
      Tignition = PD.get_datum<double>("T_ignition"),
      Tmin = PD.get_datum<double>("T_min"),
      Tin = Tmin,
      halfdiam = 0.5*diam,
      Lfrac = PD.get_datum<double>("channel_length_frac"),
      v1_in = PD.get_datum<double>("v1_in"),
      v2_in = PD.get_datum<double>("v2_in"),
      burnerLength = (b-a)/Lfrac; //PD.get_datum<double>("burner_length");
    cout << "hx = "<< hx << "  hy = "<< hy << "    yhottest = "<< yhottest << "  burner length = "<< burnerLength << endl;
    
    // BoundaryTemperatureDistribution<double> BTD(Tin,Tignition,yhottest);

    
    //either constant or parabolic
    MyVec<double> veloprofile(ny);
    double r = -halfdiam;  //R = d/2, r = -R
      for(int j = 0; j < ny; ++j){
        veloprofile[j] = InletVelocityType::velocity(r,halfdiam,v1_in); 
	r += hy;
      }

    //O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR
    double air[] = {PD.get_datum<double>("O"),    
		    PD.get_datum<double>("O2"),  
		    PD.get_datum<double>("H"),  
		    PD.get_datum<double>("OH"),   
		    PD.get_datum<double>("H2"),    
		    PD.get_datum<double>("HO2"),   
		    PD.get_datum<double>("H2O2"),     
		    PD.get_datum<double>("H2O"),     
		    PD.get_datum<double>("N2"),
		    PD.get_datum<double>("AR")
    };  

    //print_all(air,air+DataType::nspec,5,false,"initial vector");
    MyVec<double> Icomp(air,air+10);
    cout << "INITIAL COMPOSITION: "<< Icomp<< endl;

    double wbar(0.);
    for(int k = 0; k < DataType::nspec; ++k)
      wbar += Icomp[k]/DataType::molar_masses()[k]; //initially filled with air
    wbar = 1./wbar;

#ifdef USE_OPENMP
    cout  <<endl;
    FancyMessages().nice_output("OpenMP will be used", 36); //33 = yello
    cout << endl;
#endif


    //-----------------------Reynolds number---------------------------------
    MixtureAveragedViscosity<DataType,false> Mav(DataType::transport());
    StaticArray<double,DataType::nspec> Xfrac;
    for(int k = 0; k < DataType::nspec; ++k){
      Xfrac[k] = (Icomp[k]*wbar)/DataType::molar_masses()[k];
    }
    
    double mu = Mav.compute_mixture_averaged_viscosity(Tin,Xfrac);
    double rho0 = (p0*wbar)/(PhysicalConstants<double>::Rgas*Tin);
    cout << "mixture viscosity = "<< mu <<  "   char. mixture density = "<< rho0 << endl;
    double Re = Reynolds_number(rho0,v1_in,diam,mu);
    cout << endl << "REYNOLDS NUMBER = "<< Re << endl;
    //-----------------------------------------------------------------------
    
    //-------------------------Mach number-----------------------------------
    double Ma = Mach_number<DataType::ENCOD>(Tin,v1_in,Icomp);
    cout  << "MACH NUMBER = "<< Ma << endl; 
    //-----------------------------------------------------------------------

    cout << endl << "inlet velocity type   = "<< InletVelocityType::Value << endl;

    typedef BoundaryTemperatureDistribution<double,FDMSettings::WallTemperatureType> BoundaryWallTemperatureType;
    cout  << "wall temperature type = "<< ascii(BoundaryWallTemperatureType::Value) << endl;
    
    double x_ignition = (b-a)/PD.get_datum<double>("channel_length_frac");
    cout << "x_ignition = " << x_ignition << endl;

    BoundaryWallTemperatureType BWT(Tmin,Tignition,x_ignition,PD.get_datum<double>("plateau_width")*hx);
    MyVec<double> walltempprofile(nx);
    double xdir(0.);
    for(int i = 0; i < nx; ++i){
      walltempprofile[i] = BWT(xdir);
      xdir += hx;
    }
    

    cout << "wbar = "<< wbar << endl;
    double troom = Tin;
    double rho = (p0*wbar)/(PhysicalConstants<double>::Rgas*troom);
    cout << "rho = "<< rho << endl;

    //initally domain is filled with air at room temperature
   
#ifndef IGNITION_AT_HOT_WALLS
    double xdist = 0.; 
#endif
      double temper = troom;
    double y(0.);
    cout << "y-direction start: "<< y << endl;
    for(int i = 0; i < nx; ++i){
      y = 0; //reset y-direction
   
#ifndef IGNITION_AT_HOT_WALLS
   //xdist += hx;
      xdist = a + i*hx; //a = 0, 
#endif
      for(int j = 0; j < ny; ++j){
	y = j*hy; 
	int ij = i + nx*j;

	 
	//u0[ij] = rho;
	
	 u0[ij+1*npt] = veloprofile[j];

	 u0[ij+2*npt] = v2_in;
	 
#ifdef IGNITION_AT_HOT_WALLS
	 //!==== TODO: (un)comment
	 if((j == 0) || (j == ny-1)){ //DOWN or UP boundary (wall)
	   u0[ij+1*npt] = 0.;  //no slip at wall
	   u0[ij+2*npt] = 0.;
	   temper = walltempprofile[i];
	 }
	 else
	 temper = troom;

	  //!hot region inside domain
	  // if(is_contained(xdist,x_ignition-hx,x_ignition+hx)){
	 //   if((y <= (burnerSpikeLength_y*hy)) || (y >= (diam-burnerSpikeLength_y*hy))){
	 //     temper = Tignition;
	 //   }
	 //}
	 
	 //!=============================================================
#else
	 //!==== TODO: burner placed inside domain (1st argument of if)
	 if((xdist <= burnerLength) &&
	    (is_contained(yhottest,y-hy,y+hy))){
	   temper = Tignition;
	 } 
	 else
	   temper = troom;
#endif 
	 
	 u0[ij] =  (p0*wbar)/(PhysicalConstants<double>::Rgas*temper);

	 

	 u0[ij+3*npt] = temper;
	 

	 for(int k = 0; k < DataType::rednspec; ++k) //===RED
	   u0[ij+(4+k)*npt] = Icomp[DataType::rpv_index()[k]]; // mixture everywhere

      } //end domain
    }




    const bool SPARSE = true;
    const bool SCALESYS = true; //try preventing ill-conditioning of matrix
    const bool PATTERNCREATIONVIAAD = true; //true: AD, false: FD

    
    MOL<2,double, ReducedConaireH2ReactiveNS2D,SPARSE,SCALESYS,PATTERNCREATIONVIAAD> mol;
    mol.solve(u0,"data/conaireh2.dat","REDUCED_CONAIRE_H2_2D_");
    
    cout << endl<< "REYNOLDS NUMBER = "<< Re << endl;
     cout  << "MACH NUMBER = "<< Ma << endl;
    /*
  ParameterData PD;
    PD.read_from_file("data/h2gri.dat");

    //O, O2, H, OH, H2, HO2, H2O2, H2O, N2
    double air[] = {0.,    
		    0.18, //0.21,  //nothing works for pure air without H2!! 
		    0.,
		    0.,
		    0.36, //0.29, //5.5e-07, //0.29,
		    0.,
		    0.,
		    0., //0.00999945, //0.01,
		    0.46  //0.78
    };

    double massbl = 0;
    for(int k = 0; k < 9; ++k)
      massbl += air[k];
    cout << "massbl = "<< massbl << endl;
    adonis_assert(is_zero(1.-massbl));

    MyVec<double> u0(10);
    for(int k = 0; k < 9; ++k)
      u0[k] = air[k];
    u0[9] = PD.get_datum<double>("T_ignition");

     MOL<0,double,GriMech30Exerpt,true,false> mol;
     MyVec<double> SolFull = mol.solve(u0,"data/h2gri.dat");
      
     GNUplot showspec;
     //showspec("set ylabel \"Y\"; set xlabel \"t\";plot for[i=2:10] 'Gri30excerpt' using 1:i with lines lw 2 notitle");
     showspec("set xlabel \"t\"; plot 'Gri30excerpt' using 1:3 with lines title \"O2\", 'Gri30excerpt' using 1:6 with lines title \"H2\",  'Gri30excerpt' using 1:10 with lines title \"N2\"");

     GNUplot show;
     show("plot 'Gri30excerpt' using 1:11 with lines lw 2 title \"$T$\"");
*/
  }


   else if(example == 7){

#if FULL_MODEL == 1
     ParameterData PDf;
     PDf.read_from_file("2Dsettings.dat");
     string fname = PDf.get_datum<string>("filename"); 

     MyVec<double> u0f(4);
     u0f <<= 0.0, 0.8, 0.2, PDf.get_datum<double>("T_ignition");
     

     MOL<0,double,O3Decomposition,true,false> molf;
     MyVec<double> SolFull = molf.solve(u0f,"2Dsettings.dat");
     
     GNUplot showY0f;
     showY0f("set xlabel \"t\"; plot '"+fname+"' using 1:2 with lines lw 2 title \"Y_0\"");
      GNUplot showTf;
     showTf("set xlabel \"t\"; plot '"+fname+"' using 1:5 with lines lw 2 title \"T\"");

#endif

#if SWITCH_ON_REDUCTION == 1
     ParameterData PD;
     PD.read_from_file("2Dsettings.dat");

     MyVec<double> u0(2);
     u0 <<= 0.0, PD.get_datum<double>("T_ignition"); //just entry from full initial vector
     
     string outname = "O3_red";
     MOL<0,double,ReducedO3Decomposition,true,false> mol;
     MyVec<double> SolRed = mol.solve(u0,"2Dsettings.dat",outname);
   
     GNUplot showY0;
     showY0("set title \"Y_0 (RED)\";set xlabel \"t\"; plot '"+outname+"' using 1:2 with lines lw 2 title \"Y_0\"");
      GNUplot showT;
     showT("set title \"T (RED)\";set xlabel \"t\"; plot '"+outname+"' using 1:3 with lines lw 2 title \"T\"");

#if FULL_MODEL == 1
     cout << "TIME: "<< endl<< "  full: "<< molf.elapsed_time() << endl << "  red:  "<< mol.elapsed_time() << endl;
     is_faster(mol.elapsed_time(),molf.elapsed_time());
     //calculate l2 error at tend
     cout << "|Y_O - Y_O_red| = "<< Abs(SolFull[0]-SolRed[0]) << endl;
     cout << "|T - T_red| = "<< Abs(SolFull[3]-SolRed[1]) << endl;
#endif

#endif
   }

   else if (example == 8){
     
    
     cout << "CONTINUOUS STIRRED TANK REACTOR" <<endl;
     cout << "does not work properly..."<<endl;
     //==========Settings =============================
     // ParameterData PD;
     // PD.read_from_file("data/putzi.dat");
     // double t0(PD.get_datum<double>("t0"));
     // double tend(PD.get_datum<double>("tend"));
     // double h = PD.get_datum<double>("hEstim"); 
     

     // //projection algorithm
     // double alp = PD.get_datum<double>("gamma"),
     //   beta = PD.get_datum<double>("beta"),
     //   tol = PD.get_datum<double>("tol");
     // int maxit = PD.get_datum<int>("maxit"), //2000;
     //   maxarmijo = PD.get_datum<int>("maxarmijo"); 

     // //===========end settings ========================

     // typedef CSTR<double> ModelType;
     // ModelType Fun(4);

     // Fun.set_controls(21.,0.);

     // ClassicalRK<ModelType> rk4(Fun);
     
     // MyVec<double> x(4), q(2);

     // cout << "FIXED SYSTEM QUANTITIES:"<<endl 
     //      << "------------------------"<<endl;
     // MyVec<double> Q(4), R(2), A(4+2);
     // Q <<= 0.2, 1.0, 0.5, 0.2;
     // R <<= 0.5, 5.e-07;
    
     // concatenate(A,Q,R);
     // cout << "A = "<< A << endl;  //turns out to be the Hessian! (always p.d.)

     // //reference values
     // MyVec<double> xS(4), qS(2), wS(4+2), wlow(6), wup(6),
     //   w(6);  //w = [x,q]
     // xS <<= 2.1402, 1.0903, 114.19, 112.91;
     // qS <<= 14.19, -1113.5;

     // concatenate(wS,xS,qS);
     // cout << "wS = "<< wS << endl; 
     // //first bounds on x then on q
     // wlow <<= 0.,0., 60., 60.,   3.0, -9000.;  //temperature not become too small
     // wup <<= 1e+20, 1e+20, 1e+20, 1e+20,  35., 0.;
     // cout << "wlow = "<< wlow << endl;
     // cout << "wup = "<< wup << endl;

     // cout << "OPTIMAL CONTROL PROBLEM:" << endl
     //      << "------------------------"<<endl;
     // //initial values
     // x <<= PD.get_datum<double>("cA"), PD.get_datum<double>("cB"),
     //   PD.get_datum<double>("theta"), PD.get_datum<double>("thetaK");
     // q <<= PD.get_datum<double>("q1"), PD.get_datum<double>("q2");
     
     // concatenate(w,x,q);  //starting value 

     //  //impl. Euler
     // JacS<double,CSTR,MyVec> Jac;
     // Jac.set_with_init(4,4,x);
     // MyVec<double> xprev(x),IterMatrix(4*4),b(4); 
       


     // ProjectOntoBounds<MyVec<double> > P_S(wlow,wup);

     // cout << "w = "<< w << endl;
     // SpecialFunctional<double> Cost;
   
     
     // ofstream of("cstr.dat");

     // double time(t0);
     // CPUTime<double> cpu;
     // cpu.start();
     // MyVec<double> wprev(w),grad(6), gradnext(6), g0(6), ray(6), s(6), y(6), p(6), eval(6), gradnew(6), pjg(6), tmp(6), tmp2(6); //needed for PABB stepsizes
     // MyVec<double> hess(6), pgc(6), xt(6), pl(6); 
     // hess = h*A;    //note Hessian is diagonal, const. and p.d.
    
     // MyVec<double> w1(6),w2(6);


     // typedef Norm<'2',double> NormType;
     // int count(0);
     // while(time < tend){
     //   wprev = w;   //save old value
     //   // xprev = x;

     //   of << setprecision(12) << time << "    "; 
     //   for(int k = 0; k < 6; ++k) of << w[k] << "   ";
     //   of << endl;
     //   cout << count << ".)      time = "<< time << "    h = "<< h << endl;
     //   cout << "STATE: "<< x << endl; 
     //   cout << "CONTROLS:  q1 = "<< w[4] << "   g2 = "<< w[5] << endl;
     //   //! SOLVE ODE: FIRST DISCRETIZE
     //   Fun.set_controls(w[4],w[5]);  //controls are assigned
      
     //   //===TODO
     //   for(int k = 0; k < 4; ++k)   //assign values of states
     //   	 x[k] = w[k];   
       
     //    //CSTR equations integrated via time stepper 
     //   // x = rk4.step(h,x);       //CLASSICAL RUNGE-KUTTA
     //   xprev = x;
       
     //   x = xprev;   //start guess
     //   size_t l;
     //   for(l = 1; l <= 7; ++l){
     //   	 b = -(x-xprev-h*Fun(x));
     //   	 IterMatrix = -h*Jac.jacobian(x);
     //   	 update_diagonal<AddBasicElements>(IterMatrix,4,1.);
     //   	 good_square_solve(IterMatrix,4,b,1); //b contains solution now
     //   	 x += b;
     //   	 if(Norm<'2',double>::norm(b) <= 1.e-06)
     //   	   break;
     //   } //end Newton iteraion
     //   if(l == 7)
     //   	 ADONIS_INFO(Information,"Newton iteration did not converge within given number of iterations");
     //   xprev = x;    //store current approximation
     //   //! end IMPL. EULER
       

     //   cout << "ODE: x = "<< x << endl;

     //   //!THEN OPTIMIZE
     //   //START: OPTIMIZE · OPTIMIZE · OPTIMIZE · OPTIMIZE · OPTIMIZE · OPTIMIZE
       
     //   for(int k = 0; k < 4; ++k)
     // 	 w[k] = x[k];            //assign back new values from integration

     //   cout << "w (before minim) = "<< w << endl;

     //   //! ++++ compare with KELLY's code
     //   //%%% put initial iterate in feasible set
     //   tmp = w- P_S(w);   // w and xc are the same
     //   if(NormType::norm(tmp) > 0.){
     // 	 cout << "initial iterate not feasible. Project onto bounds..."<<endl;
     // 	 w = P_S(w);
     //   }
     
       
     //   //! f(w) = 0.5*h*(L(w)+L(wprev)) + K
     //   grad = Cost.gradient(h,w);  //h*A*(w-wS);   //f'(w)
     //   tmp = w - grad;
     //   pgc = w - P_S(tmp);

     //   double fc = Cost(w,wprev,h);
     //   cout << "COST(w) = "<< fc << endl;

     //   //LOOP: PROJECTED GRAD ITERATION 
     //   int itc(1), numf(1);
     //   double lambda(1.);
     //   while(NormType::norm(pgc) > tol && itc <= maxit){
     // 	 itc++;
	
     // 	 //-------------- BARZILAI-BORWEIN -----------------------
     // 	 w1 = w;
     // 	 grad = Cost.gradient(h,w1);
     // 	 //stopping criterion
     // 	 tmp = w1-grad;
     // 	 pgc = P_S(tmp)-w1;
	 
	 
     // 	 tmp = w - lambda*grad;
     // 	 w = P_S(tmp);               //update w
     // 	 tmp2 = Cost.gradient(h,w);
     // 	 s = w-w1;
     // 	 y = tmp2 - grad;
     // 	 double denom = dot(s,y);
     // 	 if(is_zero(denom))
     // 	   ADONIS_ERROR(ZeroDivision,"denom = 0 in iteration i = "<<itc << " at t = "<< time << ".");
     // 	 double aBB1 = dot(s,s)/denom;
     // 	 lambda = aBB1;

     // 	 // //alternating Barzilai-Borwein
     // 	 // double denom1 = dot(y,y);
     // 	 // if(is_zero(denom1))
     // 	 //   ADONIS_ERROR(ZeroDivision,"dot(y,y) = 0 in iteration i = "<<itc << " at t = "<< time << ".");
     // 	 // double aBB2 = dot(s,y)/denom1;

     // 	 // if(!is_even(count))  //odd iteration index
     // 	 //   lambda = aBB1;
     // 	 // else                //even iteration index
     // 	 //   lambda = aBB2;
     // 	 //----------------------------------------------------------

     // 	 // lambda = 1;
     // 	 // tmp =  w - lambda*grad;
     // 	 // xt = P_S(tmp); 
     // 	 // double ft = Cost.cost() + Cost.one_segment(xt,wprev,h);  //(xt,wprev,h);  //compute f(w), note wprev = w(t_i)
     //     //                               //and w = w(t_i+1)
	
     // 	 // numf++;


     // 	 // int iarm(0);
     // 	 // pl = w - xt;
     // 	 // double fgoal = fc - dot(pl,pl)*(alp/lambda);

     // 	 // //%%% simple line search due to Armijo
     // 	 // cout << " ARMIJO: "<< endl;
     // 	 // cout << " -------" << endl;
     // 	 // while(ft > fgoal){
     // 	 //   cout << "ft = "<< ft << "   fgoal = "<< fgoal << endl;
     // 	 //   lambda *= beta;
     // 	 //   //cout << "lambda = "<< lambda << endl;
     // 	 //   iarm++;
     // 	 //   tmp = w-lambda*grad;
     // 	 //   xt = P_S(tmp);
     // 	 //   pl = w-xt;
     // 	 //   ft = Cost.cost() + Cost.one_segment(xt,wprev,h); //Cost(xt,wprev,h); 
     // 	 //   numf++;
     // 	 //   if(iarm > maxarmijo)
     // 	 //     ADONIS_ERROR(IterationError,"Armijo error in gradient projection");
     // 	 //   fgoal = fc -dot(pl,pl)*(alp/lambda);
     // 	 // }
     // 	 // w = xt;  //if everything has worked well, that's the output
     // 	 // fc = Cost.cost() + Cost.one_segment(w,wprev,h); //Cost(w,wprev,h); 
     // 	 // numf++; 
     // 	 // grad = Cost.gradient(h,w);      //compute f'(w)  	 
     // 	 // tmp  = w - grad;
     // 	 // pgc = w - P_S(tmp);
       
     // 	 // // cout << "||pgc|| = "<< NormType::norm(pgc) << endl;
     //   } //end loop
     //   if(itc == maxit)
     // 	 ADONIS_ERROR(IterationError, "Maximum number of optimization iterations reached; maxit = "<< maxit<<".");
     //   else
     // 	 cout << "#iters: "<< itc << "."<< endl;
     //   //END: OPTIMIZE · OPTIMIZE · OPTIMIZE · OPTIMIZE · OPTIMIZE · OPTIMIZE

     //     cout << "w (after minim) = "<< w << endl;

     //   time += h; //increment time
     //   count++;
     // }
     // of.close();

     // double elapsed = cpu.stop();
     // time_in_human_readable_form(elapsed);
     
     // GNUplot show1;
     // show1("plot 'cstr.dat' using 1:2 with lines lw 2 title \"$x_1$\", 'cstr.dat' using 1:3 with lines lw 2 title \"$x_2$\"");
     // GNUplot show2;
     // show2("plot 'cstr.dat' using 1:4 with lines lw 2 title \"$T$\"");
     // GNUplot show3;
     // show3("plot 'cstr.dat' using 1:6 with lines lw 2 title \"$q_1$\"");
     // GNUplot show4;
     // show4("plot 'cstr.dat' using 1:7 with lines lw 2 title \"$q_2$\"");



   }
   else if (example == 9){
     MyVec<double> v(5);
     v <<= 1,2,4,7,14;

     int chs = 2;
     double sm(0.);
     for(int l = 0; l < 5; ++l){
       if(l != chs)
	 sm += v[l];
     }
     cout << "sum = "<< sm << endl;
   }
   else if(example == 10){
       double c = 0., d = 0.005;
    int Ny = 9;

    double hy = (d-c)/(Ny-1),
      v1max = 0.35;
   
    double R = d/2., r = -R;
    string fname = "parabolicinlet.dat";
    ofstream of(fname.c_str(),ios_base::out);
    for(int j = 0; j < Ny; ++j){
      cout << " r = "<< r << "   velocity: "<< zero(parabolic_inlet_velocity(r,R,v1max)) << " m/s" << endl;
      of << j << "  " << zero(parabolic_inlet_velocity(r,R,v1max)) << endl;
      r += hy;
    }
    of.close();


    bool wannaEps = false; //true; //false;
    
    string toEps;
    if(wannaEps)
      toEps = "set terminal postscript eps color enhanced \"Helvetica\" 20 dashed; set output \"inletvelocityprofile.eps\";";
    
    

    string add = "unset ytics;";//"set style arrow 1 head back filled linetype 1 linecolor rgb \"dark-violet\"  linewidth 2.000 size screen 0.025,30.000,45.000;";
    GNUplot plot;
    // gnuplot coordinates (x1,y1) to (x2,y2)  
    plot(toEps + add+"set label \"R\" font \"Times-Italic,24\" at first 0.72, first 6.35;set xtics nomirror; set ytics nomirror; unset ytics; set xrange [0:0.8]; set ylabel \"inlet boundary \";\
set xlabel \"velocity in x-direction (m/s)\";\
set style arrow 3 head back filled linetype 1 linecolor rgb \"black\"  linewidth 2.000 size screen 0.030,15.000,45.000;\
set arrow 1 from 0,1 to  0.153125,1 lw 2 arrowstyle 3;\
set arrow 2 from 0,2 to  0.2625,2 lw 2 arrowstyle 3;   \
set arrow 3 from 0,3 to  0.328125,3 lw 2 arrowstyle 3; \
set arrow 4 from 0,4 to  0.35,4 lw 2 arrowstyle 3; \
set arrow 5 from 0,5 to  0.328125,5 lw 2 arrowstyle 3; \
set arrow 6 from 0,6 to  0.2625,6 lw 2 arrowstyle 3; \
set arrow 7 from 0,7 to  0.153125,7 lw 2 arrowstyle 3; \
set arrow 8 from 0.7,4 to 0.7,8;\
set arrow 9 to 0.7,4 from 0.7,8 lw 1;\
set arrow 10 from 0,4 to 0.8,4 nohead lw 1 lt 0;\
          plot \"" + fname + "\" using 2:1 with lines lw 2 lc rgb \"blue\" notitle");

   
   
   }
   else if (example == 11){
     cout << "Test O3 mechanism with full Navier-Stokes:"<<endl;
   
     typedef ThermoData4Mechanism<double,3> DataType;
     ParameterData PD;
    PD.read_from_file("2Dsettings.dat");
  
    //================== TODO ============================
    //=== change model type here
    typedef OzoneReactiveNS2D<double> ModelType;
    //====================================================
 
    int npt = PD.get_datum<int>("Nx")*PD.get_datum<int>("Ny"),
      tot = ModelType::nprim*npt;
    cout << "tot = "<<tot << endl;
  
    //……………………………………………… ODE …………………………………………………………………………
    //initial value 
    MyVec<double> u0(tot);
  
 
    int nx = PD.get_datum<int>("Nx"),
      ny = PD.get_datum<int>("Ny");
    double a = PD.get_datum<double>("a"),
      b = PD.get_datum<double>("b"),
      c = PD.get_datum<double>("c"),
      diam = PD.get_datum<double>("d") - c,
      halfdiam = diam/2.,
      hy = diam/(ny-1),
      hx = (b-a)/(nx-1),
      T_inlet = PD.get_datum<double>("T_min"),
      v1 = PD.get_datum<double>("v1_in");;
  
    cout << "T_inlet = "<<T_inlet <<endl;

    cout << "hx = " << hx << "   hy = "<< hy << "  halfdiam = "<< halfdiam<< endl;

    BoundaryTemperatureDistribution<double> BTD(PD.get_datum<double>("T_min"),PD.get_datum<double>("T_ignition"),PD.get_datum<double>("hottest_x"));

    double wbar = 1./(0.8/DataType::molar_masses()[1] + 0.2/DataType::molar_masses()[2]); 

    double temp;
    for(int i = 0; i < nx; ++i){
      double x = a + i*hx; 
      double radius(0.);
      for(int j = 0; j < ny; ++j){
	int ij = i + nx*j;

	
	u0[ij+1*npt] =  zero(parabolic_inlet_velocity(radius,diam,v1));
	radius += hy;

	u0[ij+2*npt] = 0.;

	temp = T_inlet;
	u0[ij+3*npt] = temp;

        if(j == ny-1){
	  temp = BTD(x);
	  u0[ij+3*npt] = temp; //parts of upper bdy = ignition zone
	}
	

	u0[ij+4*npt] = 0.0;
	u0[ij+5*npt] = 0.8;
	u0[ij+6*npt] = 0.2;

	u0[ij] = (101325.*wbar)/(PhysicalConstants<double>::Rgas*temp);
      }
    }

    
    const bool SPARSE = true;
    const bool SCALESYS = true; //try preventing ill-conditioning of matrix

    //MOL2D<double,ReactiveNS2D,SPARSE,SCALESYS> mol;
    MOL<2,double,OzoneReactiveNS2D,SPARSE,SCALESYS> mol;
    mol.solve(u0,"2Dsettings.dat","O3_2D_COMPLETE_NSE_FULL","elapsed_time_of_full_O3.dat");
     
   }

   else if (example == 12){
     cout << "Reynolds number = "<< Reynolds_number(0.25,0.75,0.001,5.65e-05) << endl;
   }
   else if (example == 13){
     #if USE_CPPAD
      ParameterData PD;
      PD.read_from_file("data/conaireh2.dat");

      int nx = PD.get_datum<int>("Nx"),
	ny = PD.get_datum<int>("Ny"),
	quantities = 14,
	npt = nx*ny,
	dim = quantities*npt           ;
      
      typedef CppAD::AD<double> ADType;
      typedef ConaireH2ReactiveNS2D<ADType> Model;
      Model TwoD(dim);

      
      MyVec<ADType> V(dim), W(dim);
      //!initialize V to take valid values for the mechanism
      //O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR
      double air[] = {PD.get_datum<double>("O"),    
		      PD.get_datum<double>("O2"),  
		      PD.get_datum<double>("H"),  
		      PD.get_datum<double>("OH"),   
		      PD.get_datum<double>("H2"),    
		      PD.get_datum<double>("HO2"),   
		      PD.get_datum<double>("H2O2"),     
		      PD.get_datum<double>("H2O"),     
		      PD.get_datum<double>("N2"),
		      PD.get_datum<double>("AR")
      };  
      for(int i= 0; i < nx; ++i){
	for(int j = 0; j < ny; ++j){
	  V[i + nx*j + 0*npt] = 2e-06;
	  V[i + nx*j + 1*npt] = 0.75;
	  V[i + nx*j + 2*npt] = 0;
	  V[i + nx*j + 3*npt] = 300.;
	  
	  for(int k = 0; k < 10; ++k)
	    V[i + nx*j + (4+k)*npt] = air[k];
	  
	}
      }
      CppAD::Independent(V);

      W = TwoD(V);
      CppAD::ADFun<double> sq;
      sq.Dependent(V,W);
      sq.optimize();

      typedef SparsityPattern<std::set<std::size_t> > twdpat_type;
      typedef typename std::set<std::size_t>::iterator SetIterType;
      typedef twdpat_type::SVecType StdVecBoolType;
      twdpat_type Sp2d(dim,dim);
      
      StdVecBoolType smtx = Sp2d.calc_sparsity(sq); //sparsity of jacobian
      //cout << "Sparsity pattern 2D" << endl << Sp2d << endl;
      Sp2d.diagonal_never_zero();

      cout << "size = "<< Sp2d.size() <<"   #(nnz) = " << Sp2d.number_of_nonzeros() << " (out of "<< Sp2d.size()<< " x "<< Sp2d.size() << " = "<< ntimes<2>(Sp2d.size())<< ")" << endl;
      
      string fname = "SparsPattConaireH2",
	mtxdatafile = fname+".dat";
      //print_matrix_pattern_2_file(smtx,dim,mtxdatafile);
      std::ofstream ofs(mtxdatafile.c_str(),std::ios_base::out);
      for(size_t i = 0; i < Sp2d.size(); ++i){
	SetIterType setIt = Sp2d[i].begin();
	for(size_t j = 0; j < Sp2d.size(); ++j){
	  if((j == (*setIt)) && (Sp2d[i].size() != 0)){
	    ofs << 1 << " ";
	    setIt++;
	  }
	  else{
	    ofs << 0 << " ";
	  }
	}
	ofs<<endl;
      }

      ofs.close();
      
      string pltmtxpath = "sh/./plotmatrix.sh " + mtxdatafile;
      system(pltmtxpath.c_str());
      string epsfile = fname+".eps";
      system(("gimp "+epsfile).c_str());
#endif
      
   }
   else if (example == 14){
     
      ParameterData PD;
      PD.read_from_file("data/conaireh2.dat");

      int nx = PD.get_datum<int>("Nx"),
	ny = PD.get_datum<int>("Ny"),
	quantities = 14,
	nchem = 10,
	npt = nx*ny,
	dim = quantities*npt;

      double Tlow = PD.get_datum<double>("Tlow"),
	Tup = PD.get_datum<double>("Tup");

      double hMin = PD.get_datum<double>("hmin");
      cout << "hMin = "<< hMin << endl;

      ReadInMatrix<double> RIM;
      
      MyVec<double> u(dim), unext(dim), ueul(dim);
      RIM.clear();
      //note that u already contains the primitive variables
      RIM.read_from_file("CONAIRE_H2_2D_FULL41.gnu");
      RIM.assign_2_2D_MOL_slide(u,nx,ny,quantities+1,false); //pressure
      RIM.print_matrix("Printed matrix");
      // RIM.read_from_file("CONAIRE_H2_2D_FULL42.gnu");
      // RIM.assign_2_2D_MOL_slide(unext,nx,ny,quantities+1,false); //pressure
      typedef ConaireH2ReactiveNS2D<double> ModelType;
      typedef ModelType::DataType DataType;
      typedef MyVec<int> IndexVecType;
      ModelType fun(dim);
      
      CPUTime<double> cpu;
      cpu.start(false);
    
    
      
      int excessSpec = DataType::N2;
     
      int off;
      cout << "show me odd Y_k data:"<<endl;
      int ctodd(0);
      for(int i = 0; i < nx; ++i){
	for(int j = 0; j < ny; ++j){
	  for(int k = 0; k < nchem; ++k){
	    off = i + nx*j + (4+k)*npt;
	    if((u[off] > 1.) ||
	       (u[off] 	< 0.)
	       ){
	      ctodd++;
	      cout << "("<<i<<", "<<j<<"):  " << endl << "-------- ";
	      //print out the violations
	      for(int l = 0; l < nchem; ++l){
		 
		cout << u[i + nx*j + (4+l)*npt] << ", ";
	      }
	      cout << endl;
	      break;
	     }//end if
	  }

	  if(
	     (u[i + nx*j + 3*npt] < Tlow) || (u[i + nx*j + 3*npt] > Tup)
	     ){
	    ctodd++;
	    cout << "("<<i<<", "<<j<<"):  " << endl << "-------- ";
	    cout << u[i + nx*j + 3*npt] << ", ";
	    
	      cout << endl;
	      break;
	  }
	 
	}
      }
      if(ctodd != 0){
	cout << "There were "<< ctodd << " obvious violations in primitive data found" << endl;
      }
      else
	cout << "No obvious violations found in Y_k and T."<<endl;
      

      MyVec<double> cp(u);  //old value
     //  cout << endl << " test nearby_functionality in 2D:" << endl;
     //  double avg(0);
     //  int idx = 40,
     // 	jdx = 35,
     // 	spec = 8;
      
     // nearby_value_2D(avg,u,idx,jdx,nx,ny,spec);
     
     // cout << "\"avg\" ("<<idx<<","<<jdx<<","<<DataType::species_names()[spec]<<") = " << avg << endl; 


      cout << endl << "REPAIR SLICE before proceeding"<<endl;
      typedef RepairPhysicalQuantity<double,2,true,'f'> RepairType;
      RepairType Repairer;
      Repairer.init("data/conaireh2.dat",u);
      cout << "------"<<endl;
      Repairer.repair(u,excessSpec); //repairs now u
      Repairer.info(41,0);
      
      int countdiff(0);
      for(int i = 0; i < nx; ++i){
	for(int j = 0; j < ny; ++j){
	  for(int quant = 0; quant < 14; ++quant){
	    off = i + nx*j + quant*npt;
	  
	    if(!is_zero(Abs(u[off]-cp[off]))){
	      if(quant == 3){ //only temperature
		cout << "("<<i<<","<<j<<"): u_i,j_NEW = "<< u[off] << "  u_i,j_OLD = "<<  cp[off] << endl; 

	      }
	      countdiff++;
	    }
	  }
	}

      }
      cout << "# of diffs: "<< countdiff << " (out of "<< nx*ny*14 << ")."<<endl;

      double tm = cpu.stop(false);

      IndexVecType select;
      PrintSolution PS;  
      PrintMOLData<2,double,true> PMD, DMP;
      PMD.init("data/conaireh2.dat","Something",PS,15);
      PMD.write_2_file(hMin,0,2.21875e-09,u,"Something",PS,select,DataType::molar_masses(),tm,0);

      PrintSolution SP;  
      DMP.init("data/conaireh2.dat","OrigSlide",SP,15);
      DMP.write_2_file(hMin,0,2.21875e-09,cp,"OrigSlide",SP,select,DataType::molar_masses(),PD.get_datum<double>("elapsed_time_so_far"),0);


      // cpu.reset();
      // cpu.start(false);

      // cout << "EXPLICIT EULER:"<<endl;
      // ueul = u + hMin*fun(u);
      // tm = cpu.stop(false);

      // PrintSolution PS1;  
      // PrintMOLData<2,double,true> PMD1;
      // PMD1.init("data/conaireh2.dat","Eul",PS1,14);

      // PMD1.write_2_file(0,2.8571875e-09,ueul,"Eul",PS1,select,DataType::molar_masses(),tm,0);
  
   }
   else if (example == 15){
     cout << endl << "Test this:" <<endl;
     ParameterData Param;
     Param.read_from_file("../fdm/data/conaireh2.dat");

     int Nx = Param.get_datum<int>("Nx");

     cout << "Nx (from file) = "<< Nx << endl;

     cout << endl << " Extract slide" <<endl;
     ReadInMatrix<double> RIM;
     RIM.read_from_file("CONAIRE_H2_2D_FULL55.gnu");
     cout << "rows = " << RIM.rows() << "   cols = "<< RIM.cols() << endl;
     
     double x = 7.6923076923077e-07;
     string slideoutfile1 = "Slide_"+Num2str(x)+".dat";
     RIM.fixed_x_2D(x,slideoutfile1);
     cout << "count = "<<RIM.count() << endl;
     cout << "Gnuplot:"<<endl;
     GNUplot plotme;
     plotme("set key box left; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \"y\"; set ylabel \"value(x_FIX,y,.)\"; plot \""+slideoutfile1+"\" using 2:17 with lines lw 3 title \"pressure\"");

     double y = 3.36e-07;
     string slideoutfile2 = "Slide_"+Num2str(y)+".dat";
     RIM.fixed_y_2D(y,slideoutfile2);
     GNUplot plotme2;
     plotme2("set key box left; set mxtics \"2\";set mytics \"2\"; set grid xtics ytics mxtics mytics; set xlabel \"x\"; set ylabel \"value(x,y_FIX,.)\"; plot \""+slideoutfile2+"\" using 1:17 with lines lw 3 title \"pressure\"");
   }

   else if(example == 16){
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
     string slide = "CONAIRE_H2_2D_FULL190.gnu";
     //=======================================================================
     cout << "Nx (from file) = "<< Nx << "   Ny (from file) = "<< Ny <<endl;
     //correct hx, hy for ghostpoints
     cout << "hx = "<< hx << "   hy = "<< hy << endl;

     cout << endl << " Extract slide" <<endl;
     ReadInMatrix<double> RIM;
     RIM.read_from_file(slide);
     cout << "rows = " << RIM.rows() << "   cols = "<< RIM.cols() << endl;

     RIM.write_2_file("Test_read_in_file.dat");
     
     RIM.show_2D_MOL_boundaries(Nx,Ny,hx,hy,quantity,numberOfQuantities);
   }
  
  
   else if(example == 17){
     cout << "if you can read this, you're silly and ugly" <<endl;
     int sum(0);
#ifdef USE_OPENMP
#pragma omp parallel for num_threads(NUMTHREADS)
#endif
     for(int i = 0; i < 16; ++i){
       sum+=i;
     }
     cout << "sum = "<< sum << endl;
   }

   //default argument 
  else{
    ADONIS_ERROR(MainProgramError, "Example '"<<example<<"' hasn't been declared for \""<<argv[0]<<"\" yet. \n");

  }
  return 0;
} //END MAIN
