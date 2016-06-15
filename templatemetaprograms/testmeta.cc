#include <iostream>
#include <string>

#include <cppad/cppad.hpp>

//!HAS TO BE PUT AT THE VERY BEGINNING !!!!!!!!!!!!!!!!!!!
#include "../want_2_use.h"  
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "../common/globalfunctions.hh"
#include "simpledotproduct.hh"

#include "unrollloop.hh"

#include "commonfunctions.hh"

#include "conditionalstructures.hh"

#include "../expressiontemplates/exprvec.hh"
#include "../expressiontemplates/xsettings.hh"


//JUST TEST
#include "matrixunroller.hh"
#include "threeloopsunrolled.hh"


//test
#include "../massactionkinetics/data/thermochemicaldata.hh"
#include "../massactionkinetics/physicalconstants.hh"
#include "../moleculartransport/diffusion.hh"
#include "../moleculartransport/viscosity.hh"
#include "../moleculartransport/multicomponenttransport.hh"
#include "../moleculartransport/conductivity.hh"

#include "../moleculartransport/thermaldiffusion.hh"

#include "../linalg/linearsystemsolvers.hh"

#include "compiletimearray.hh"

#include "../graphics/gnuplot.hh"
#include "../massactionkinetics/eos.hh"

#include "functionalities4containersandscalars.hh"

using namespace std;
using namespace Adonis;
using namespace ExprTmpl;

template<class T, int N, int M>
class TinyMatrix{
public:
  TinyMatrix(){}
  
  TinyMatrix(const T& d){
    for(int i = 0; i < N; ++i)
      for(int j = 0; j < M; ++j)
	m_[i][j] = d;
  }

  T* operator[](int i) {return m_[i];}
  const T* operator[](int i) const {return m_[i];}

  T& operator()(int i, int j) {
    return m_[i][j];
  }

  const T& operator()(int i, int j) const{
    return m_[i][j];
  }


  friend std::ostream& operator<<(std::ostream& os, const TinyMatrix& m){
    for(int i = 0; i < N; ++i){
      for(int j = 0; j < M; ++j){
	os << m[i][j] << " ";
      }
      os << std::endl;
    }
    return os;
  }

private:
  T m_[N][M];
};


//fake class 
template<class T, class ITER>
class TinyTransport{
public:
  
  TinyTransport(ITER it = ITER()):it_(it){}

  T delta(int i, int j) const{
    return 0.5*it_[6*i+3]*it_[6*j+3];   //generally this is zero, hence this form
  }

private:
  ITER it_;
};

//brute force kronecker
inline int kron_delta(int i, int j){
  return ((i==j)? 1 : 0);
}


template<class V>
typename V::value_type test_formula(int i, int j, int k, const V& M){
  return (M[i]/(M[j]*(pow(M[i]+M[k],2)))*(15./2*M[j]*M[j] + (25./4 - 3)*M[k]*M[k]) );
}

template<class V, class IT>
typename V::value_type phi(size_t k, size_t j,const V& mw, const IT& eta){
  return 1./sqrt(8.)*pow((1+mw[k]/mw[j]),-0.5)*pow((1 + sqrt(eta[k]/eta[j])*pow(mw[j]/mw[k],0.25)),2.);
}


template<class T>
class TestMe{
public:
  enum{kaitain = 3, giediprime = 5};
};


template<class T1, class T2, class T>
void assign_prop_zero(bool cool, T1& x1, const T2& x2, const T& molec){

  (cool == false) ? smart_assign(x1,x2) : smart_assign(x1,((x2 <= molec) ? 3.5 : 0));
}

// MAIN PROGRAM · MAIN PROGRAM · MAIN PROGRAM · MAIN PROGRAM ·  MAIN PROGRAM
int main(){
  
  const size_t fac = 7;
  cout<<"  The "<<fac<<"-th fibonacci number is "<< Fibonacci<7>::fib <<"."<<endl;

  cout<<endl<<endl;
  cout<<"Simple metaprogramme to execute the dot product in am efficent way:"<<endl;
  
 
   MetaVec<3, double> a, b;
   
   a[0] = 1; a[1] = 2; a[2] = 3;
   b[0] = 5; b[1] = 6; b[2] = 7;

   cout << a <<endl;
   cout << b <<endl;
   cout << " Dot product <a,b> = " << dot(a,b) <<endl; 
   
   cout<<endl<<"Test some loop unwinding" <<endl;
   
   UnrollLoop<0,10>::loop();
   
   

   cout<<endl<<endl;
   
   
   
   cout<<endl<<"Test IfElse class:"<<endl;
   const int x =1456;


   IfElse<x == 1456>::statement();

   cout<<endl<<"Test switch Metaprogram:"<<endl<<endl;
   Switch<0>::s();
   Switch<1>::s();
   Switch<2>::s();
   cout<<endl;
   Switch<16>::s(); //default statement
   
  

   cout<<endl<<"Test factorial the Marc-way:"<<endl;
   const int tof = 7;
   cout<< " "<<tof<<"! = "<< Factorial<tof>::result <<"."<<endl;


   cout << endl<< "Loop again:"<<endl;
  
   const int dim = 7;
   typedef ExprTmpl::MyVec<double> VType;
   double av[] = {0.5, -1.35, 2.45, -1.5, 1.75,0.15,3.85};
   VType v1(dim),
     v2(av,av+dim);
     
   //cout << "v1 = "<< v1 << endl;
   cout << "BEFORE: v1 = "<< v1 << endl;
   
   
   //static const int d[2] = {0,7}; //note that this cannot appear in constant expression!
  

   typedef UnrollLoop<0,7> LoopType;
   
   //!not very promising, though
   //LoopType::assign(v1[LoopType::index()],v2[LoopType::index()]);
   
  
   //LoopType::assign(v1.begin(),v2.begin()); //works
   LoopType::assign_rac2rac(v1,v2);

   cout << "AFTER:  v1 = "<< v1 <<endl;


   cout << endl<< "Assign only special values from a container"<<endl;
   double tofill[] = {1,2,3,4,5,6,7};
   ExprTmpl::MyVec<double> w1(tofill,tofill+dim),
     w2(v2);
   cout << "BEFORE: w1 = "<< w1 << endl;

   const int reddim = 3;
   size_t idx[] = {1,2,4};

   UnrollLoop<0,reddim>::assign_rac2rac(w1,w2,idx);

   cout << "AFTER: w1 = "<< w1 <<endl;


   cout<<"test this shit:"<<endl;
   double onmani[] = {0.4549999999992859, 0.7780585738931507, 0.2366143850825262, 0.3628298037265891, 0.1479999999999196, 0.01594142610843904};
    
   
   ExprTmpl::MyVec<double> IV(onmani,onmani+6);
    
  
   int index[] = {0,4};

   ExprTmpl::MyVec<double> RV(2);

   UnrollLoop<0,2>::assign_greater_2_smaller_rac(RV,IV,index);

   cout << "RV = "<< RV<<endl;

   std::vector<double> VR(2);
   VR[0] = 2.5; VR[1] = 3.25;

   UnrollLoop<0,2>::assign_smaller_2_greater_rac(IV,VR,index);
   
   cout<<" New IV = "<<IV <<endl;


   std::vector<double> vw(6);

   for(size_t i = 0; i < 6; ++i)
     vw[i] = (i+1)*0.5;

   ExprTmpl::MyVec<double> zag(6);

   UnrollLoop<0,6>::assign_rac2rac(zag,vw);

   cout << "zag = "<< zag << endl;
   

   cout << "Matrix style stuff: "<<endl;
   UnrollLoop<0,4>::matrix_style<3>(); 

   cout <<" test tuned floating point summation:"<<endl;
   const int Dim = 5;
   ExprTmpl::MyVec<double> Acc(Dim);

   Acc[0] =  1.001234987123452;
   Acc[1] = -0.671840823915887;
   Acc[2] =  4.992143473861539;
   Acc[3] = -1.501204324051362;
   Acc[4] =  0.071234945105499;
  
   double s = 0,
     c = 0;

   cout<<setprecision(16)<<endl;
   cout << "tuned fl.pt. sum = "<< UnrollLoop<0,Dim>::tuned_addition<ExprTmpl::XPCTSettings::floatingPointAdditionMethod>(Acc,s,c) <<endl;

   cout << "conventional sum = "<< conventional_sum(Acc) << endl;
   cout << "unr. covtl.  sum = "<< UnrollLoop<0,Dim>::conventional_addition(Acc) <<endl;


   cout << endl << "TEST Unrolling matrix stuff "<<endl;
   TinyMatrix<double,5,4> TM(0.);
   ExprTmpl::MyVec<double> spec(5), mw(5);

   spec[0] = 0.5; spec[1] = 2.; spec[2] = -0.75; spec[3] = 1.5; spec[4] = 1.25;
   mw[0] = -0.75; mw[1] = 0.35; mw[2] = 1.15; mw[3] = -1.5; mw[4] = -0.65;

   MatrixUnroller<0,0,5,4>::assign(TM,spec,mw);
   cout << setprecision(5) <<endl;
   cout << "TM = "<< endl << TM << endl;


   cout << endl << "3 Loops unroller"<<endl;
     

   cout << endl << endl<< "MATRIX-MATRIX-PRODUCT"<<endl;

   TinyMatrix<double,5,4> A;
   TinyMatrix<double,4,3> B;
   //result of c = a*b (not vice versa!!)
   TinyMatrix<double,5,3> C(0.);   //must initialized to be zero

   A[0][0] = 0.5; A[0][1] = 1; A[0][2] = -0.75; A[0][3] = -4.35;
   A[1][0] = 0; A[1][1] = -3.; A[1][2] = 0.15; A[1][3] = 0;
   A[2][0] = -6.05; A[2][1] = 0.05; A[2][2] = 3.95; A[2][3] = 9.5;
   A[3][0] = 4.265; A[3][1] = 1.4; A[3][2] = -2.965; A[3][3] = -1.25;
   A[4][0] = -3.25; A[4][1] = 5.65; A[4][2] = -4.2; A[4][3] = 3.5;

   cout << "A = "<< endl<<A<<endl;

   B[0][0] = -2.35; B[0][1] = 0; B[0][2] = 10.75;
   B[1][0] = 4.5; B[1][1] = 1.85; B[1][2] = -0.85;
   B[2][0] = 0.0035; B[2][1] = 0.1145; B[2][2] = 5.5;
   B[3][0] = 3.45; B[3][1] = -1.45; B[3][2] = -0.25;

   cout << endl<<"B = "<< endl<<B<<endl;

   

   //ThreeNestedLoopsUnroller<0,0,0,5,3,4>::matrix_matrix_prod(C,A,B);
   // cout << "C = " << endl << C << endl;
   cout << "C = " << endl << ThreeNestedLoopsUnroller<0,0,0,5,3,4>::matrix_matrix_prod(C,A,B) << endl;
   cout << endl<< "u might wanna cross-check it with Matlab(c)....or in conventional style"<<endl;
   
   TinyMatrix<double,5,3> CC(0.);   //must initialized to be zero
   for(int i=0; i< 5; i++){                //A.row
     for(int j=0; j< 3; j++){           //B.column
       for(int k=0; k< 4; k++){         //A.column
	 CC[i][j] +=A[i][k]*B[k][j];
	 
       }
     }
   }  
   cout << "CC = "<<endl << CC <<endl;
   

   cout << endl<< "Kronecker Delta"<<endl;
   const int ix = 3;
   const int jx = 3;

   cout << "delta_{"<<ix<<","<<jx<<"} = "<< Kronecker<ix,jx>::delta() <<endl;
   cout << "Should be negative (minus one): "<< (Kronecker<0,3>::delta() - Kronecker<2,2>::delta()) <<endl;

   const int sdim = 3;
   double sym[sdim*(sdim+1)/2] = {1,2,3,4,5,6};

   cout << "a_{0,2} = "<< sym[SymmetricDenseMatrixAccess<0,2,sdim>::offset()]<<endl; 
   cout << "a_{2,0} = "<< sym[SymmetricDenseMatrixAccess<2,0,sdim>::offset()]<<endl; 
   cout << "a_{1,1} = "<< sym[SymmetricDenseMatrixAccess<1,1,sdim>::offset()]<<endl; 



   //********************************************************************
   //*********** check TRANSPORT STUFF **********************************
   //********************************************************************
   cout << "Check handcoded stuff:"<<endl;
   //molar masses in g/mol, i.e. cgs units
   //H2,H,O2,O,H2O, OH
   double molma[] = {2.0158, 1.0079, 32., 16., 18.0158, 17.0079};
   ExprTmpl::MyVec<double> MolarMasses(molma,molma+6);

   ExprTmpl::MyVec<double> RedMolarMasses(2);
   RedMolarMasses[0] = MolarMasses[0];  RedMolarMasses[1] = MolarMasses[4]; 

   //mole fractions
   //double xfrac[] = {0.35, 0.165, 0.75, 0.95, 0.8, 0.5};
   double xfrac[] = {0.035,0.01,0.015,0.04, 0.55, 0.35};
   ExprTmpl::MyVec<double> Xfrac(xfrac,xfrac+6);

   

   TinyMatrix<double,6,6> BinDiff;

   double bd[6][6] = {
     {0.00158318,  0.00240802,  0.000891541,  0.00118851,  0.00105789,  0.00118456},  
     {0.00240802,  0.00365453,  0.00149239,  0.00204499,  0.00180367,  0.00204139},  
     {0.000891541,  0.00149239,  0.000239911,  0.000373651,  0.000315866,  0.000366196},  
     {0.00118851,  0.00204499,  0.000373651,  0.000563017,  0.00048635,  0.000554613},  
     {0.00105789,  0.00180367,  0.000315866,  0.00048635,  0.000390551,  0.000478657},  
     {0.00118456,  0.00204139,  0.000366196,  0.000554613,  0.000478657,  0.00054608}  
   };

   for(int i = 0 ; i < 6; ++i)
     for(int j = 0; j < 6; ++j)
       BinDiff[i][j] = bd[i][j];

   cout << "BinDiff = "<< endl<< BinDiff<<endl;

   TinyMatrix<double,6,6> L0000;
   
   double temperature = 1250.0,
     pressure = 1.;
     
   
   double innersum = 0;
   for(int i = 0; i < 6; ++i){
     for(int j = 0; j < 6; ++j){
       innersum = 0;              //reset!!
       for(int k = 0; k < 6; ++k){
	 innersum +=  Xfrac[k]/(MolarMasses[i]*BinDiff[i][k])*(MolarMasses[j]*Xfrac[j]*(1-kron_delta(i,k)) - MolarMasses[i]*Xfrac[j]*( kron_delta(i,j) - kron_delta(j,k) ) ); 
	 //cout << "D_{ik} = "<< BinDiff[i][k] <<endl;
       }
       L0000[i][j] = 16./(25.*1.e5)*(temperature)/pressure*innersum;
     }
   }
   
    cout << " L0000 = "<< endl<< L0000 <<endl;

    // cout <<endl<< "UNROLL ME:"<<endl;
    // double l0000[6*6];

    const bool FULL = false;
    const bool REDUCED = true;
    
    typedef ThermoData4Mechanism<double,6> DataType;
    PureSpeciesDiffusion<DataType,FULL> SpecDiff(DataType::transport());
    
    //! test copy assignment
    PureSpeciesDiffusion<DataType,FULL> PSD;
    PSD = SpecDiff; 
    
     //works   
    PSD.compute_binary_diffusion_matrix(pressure,temperature);
     // cout << "Smart computation of binary diff = "<<endl;
    // for(int i = 0; i < 6*7/2; ++i){
    //   cout << PSD.get_matrix_ptr()[i] << " ";
    // }
   
    cout << "COMPUTED Bin diff = "<< endl;
    print_all(PSD.coefficients(),PSD.coefficients()+gauss_sum(6));
    
    //reduced
    PureSpeciesDiffusion<DataType,REDUCED> RedPSD(DataType::transport());
    RedPSD.compute_binary_diffusion_matrix(pressure,temperature);
    cout << "REDUCED COMPUTED Bin diff = "<< endl;
    print_all(RedPSD.coefficients(),RedPSD.coefficients()+gauss_sum(2));
    


    CppAD::AD<double> yfrac[] = {0.01, 0.02, 0.01, 0.05, 0.6, 0.3};

    cout <<endl<<"Inbetween test:"<<endl;
    double mmw1 = 0,
      mmw2 = 0;

    ComputeTransportProperties<FULL>::mean_molecular_weight_Y<DataType>(mmw1,&yfrac[0],DataType::molar_masses());
    cout << "mmw_Y = "<< mmw1 << endl;

    ComputeTransportProperties<FULL>::mean_molecular_weight_X<DataType>(mmw2,Xfrac.begin(),DataType::molar_masses());
    cout << "mmw_X = "<< mmw2 << endl;

    
   
    mmw1 = mmw2 = 0; //reset!!

    double yfrac_red[] = {0.65, 0.35};
    double xfrac_red[] = {0.24, 0.76};

    ComputeTransportProperties<REDUCED>::mean_molecular_weight_Y<DataType>(mmw1,&yfrac_red[0],DataType::molar_masses());
    cout << "REDUCED mmw_Y = "<< mmw1 << endl;


    ComputeTransportProperties<REDUCED>::mean_molecular_weight_X<DataType>(mmw2,&xfrac_red[0],DataType::molar_masses());
    cout << "REDUCED mmw_X = "<< mmw2 << endl;


    //mixed average diffusion coeffs
    cout << "p = "<< pressure << "   T = "<< temperature << endl; 
    cout << "Xfrac = "<< Xfrac << endl;
    MixtureAveragedDiffusion<DataType,FULL> MAD(DataType::transport());
    MAD.compute_mixture_averaged_diffusion_coefficients(pressure,temperature,Xfrac);
    cout<< "Mix. avg. diff. coeffs = " << endl;
    print_all(MAD.coefficients(),MAD.coefficients()+6);
    
    //reduced stuff
    MixtureAveragedDiffusion<DataType,REDUCED> RedMAD(DataType::transport());
    MyVec<CppAD::AD<double> > RedXfrac(2);
    RedXfrac[0] = 0.44;//Xfrac[0];//
    RedXfrac[1] = 0.56;//Xfrac[4]; //
    
    RedMAD.compute_mixture_averaged_diffusion_coefficients(pressure,temperature,RedXfrac);
    cout<< "+++++ REDUCED Mix. avg. diff. coeffs = " << endl;
    print_all(RedMAD.coefficients(),RedMAD.coefficients()+2);
    
    //same as 
    for(int k = 0; k < 2; ++k)
      cout << " D^mix_"<<k<<" = "<< RedMAD[k] << endl; 

    

    cout << endl<< "*********************"<< "Test from the scratch "<<endl;
    double xfc[] = {0.227848, 0.0421941, 0.113924, 0.0168776, 0.590717, 0.00843882};
    MyVec<double> Xfc(xfc,xfc+6);
    double druck = 8314.4,
      hitze = 1000;

    MixtureAveragedDiffusion<DataType,false> ScratchMAD;
    ScratchMAD.initialize(DataType::transport());
    
    ScratchMAD.compute_mixture_averaged_diffusion_coefficients(druck,hitze,Xfc);

    cout << "mix. avg. diff coeffs = " <<endl;
    print_all(ScratchMAD.coefficients(),ScratchMAD.coefficients()+6);
    
    
    cout << endl << "REDUCTION * REDUCTION* REDUCTION  " << endl;
    cout << endl << "----------------------------------" << endl;
    MyVec<double> RdXfr(2);
    RdXfr[0] =  0.278351;
    RdXfr[1] =  0.721649;

    
    MixtureAveragedDiffusion<DataType,true> REDScratchMAD;
    REDScratchMAD.initialize(DataType::transport());
    REDScratchMAD.compute_mixture_averaged_diffusion_coefficients(druck,hitze,RdXfr);

    cout << "REDUCED mix. avg. diff coeffs = " <<endl;
    print_all(REDScratchMAD.coefficients(),REDScratchMAD.coefficients()+2);
    

    cout << endl << "VVVVVVVVV VISCOSITIES VVVVVVV"<< endl;
    cout << "temperature = "<< temperature << endl;
    PureSpeciesViscosity<DataType,FULL> PSV(DataType::transport());
    
    PSV.compute_viscosities(temperature);
    
    cout << "FULL vis.:"<<endl;
    for(int i = 0; i < 6; ++i)
      cout << PSV.get_viscosity(i) << "  ";
    cout << endl;
    
    cout << "REDUCED vis.:"<<endl;
    PureSpeciesViscosity<DataType,REDUCED> RedPSV(DataType::transport());
    RedPSV.compute_viscosities(temperature);
    
    //cout << "DIM = "<< PureSpeciesViscosity<DataType,REDUCED>::DIM << endl;
    for(int i = 0; i < 2; ++i)
      cout << RedPSV.get_viscosity(i) << "  ";
    cout << endl;
    
    cout << endl << "MMMMMMMM  MULTICOMPONENT DIFFUSION  MMMMMMM"<< endl;
    cout << "pressure = "<< pressure << endl;
    cout << "Xfrac = "<< Xfrac << endl;
    
    MultiComponentTransport<DataType,FULL> MCT(DataType::transport());

    MCT.properties(pressure,temperature,Xfrac);

    print_in_matrix_style(MCT.get_L0000_ptr(),MCT.get_L0000_ptr()+36,6);
    

    cout << endl << "...a REDUCED description:"<<endl;
    typedef double NumType;
    MyVec<NumType> ReducedXfrac(2);
    ReducedXfrac[0] = 0.21;
    ReducedXfrac[1] = 0.79;

    MultiComponentTransport<DataType,REDUCED> RedMCT(DataType::transport());
    RedMCT.properties(pressure,temperature,ReducedXfrac);

    print_in_matrix_style(RedMCT.get_L0000_ptr(),RedMCT.get_L0000_ptr()+4,2);



    cout << endl << "CCCCCCC  CONDUCTIVITIES  CCCCCCC"<<endl;
    double rho = 0.3;
    PureSpeciesConductivity<DataType,FULL> Cond(DataType::transport());
    CppAD::AD<double> temperatureAD = temperature;
    Cond.compute_conductivities(rho,pressure,temperatureAD);

    print_all(Cond.coefficients(),Cond.coefficients()+6);

    cout << "REDUCED:"<<endl;
    PureSpeciesConductivity<DataType,REDUCED> RedCond(DataType::transport());
    RedCond.compute_conductivities(rho,pressure,temperatureAD);
    print_all(RedCond.coefficients(),RedCond.coefficients()+2);
    
    cout << endl << "MIX. AVERAGED CONDUCTIVITY" << endl;
    MixtureAveragedConductivity<DataType,FULL> MAC;
    cout << "test copy constructor"<<endl;
    MixtureAveragedConductivity<DataType,FULL> MACY(DataType::transport());
    MAC = MACY;
    cout << "mix. avg. cond (full) = "<< MAC.compute_mixture_averaged_conductivity(rho,pressure,temperature,Xfrac) << endl;
    cout << endl << "...a REDUCED description:"<<endl;
    cout << "test copy cstr"<<endl;
    MixtureAveragedConductivity<DataType,REDUCED> RdefMAC(DataType::transport());
    MixtureAveragedConductivity<DataType,REDUCED> RedMAC(RdefMAC);
    cout << "mix. avg. cond (red.) = "<< RedMAC.compute_mixture_averaged_conductivity(rho,pressure,temperature,ReducedXfrac) << endl;
    

   
    double sys[] = {0.5, 1., 0.75, -2, 0.65};
    MyVec<double> Ys(sys,sys+5);
    cout << "sum_k=1^4 = "<< UnrollLoop<0,4>::sum<1>(Ys.begin()) << endl; //s = 0.25 


    cout << endl<< "EOS:"<<endl;
    CppAD::AD<double> p0 = 101315., wbar = 0.865;
    double tp = 1345.75;
    cout << "Eos = "<< EquationOfState<'i'>::density(p0,wbar,tp) << endl;
    
    
    //// works: IMPORTANT TO KNOW, gnuplot dates latex ;)
    UnrollLoop<0,TestMe<double>::giediprime>::print_cool_stuff();

    ////! just for plotting -- done
    cout << "compute mixture viscosity over time"<<endl;
    typedef ThermoData4Mechanism<double,3> O3Mech;
    MixtureAveragedViscosity<O3Mech,FULL> Mav(O3Mech::transport());
    double xo3[] = {0., 1., 0.}; //pure O2 
    MyVec<double> XO3(xo3,xo3+3);
    
    MixtureAveragedConductivity<O3Mech,FULL> MavgCond(O3Mech::transport());
    double p_atm = 101325.;  //Pa
    double Wbar = 0., rhoNow = 0.;
    for(int k = 0; k < O3Mech::nspec; ++k)
      Wbar += (XO3[k]*O3Mech::molar_masses()[k]);
   
    cout << "Wbar = "<< Wbar << endl;

    MixtureAveragedDiffusion<O3Mech,FULL> MixDiff(O3Mech::transport());

    size_t o2 = 1;  
    ofstream of;
    of.open("Mu.dat",ios_base::out);
    of<<"# T     "<<  "mu(T,X)    "<< "lambda(rho,p,T,X)       D^mix(p,T,X)" <<endl;
    // double mavgcond;
    for(double tpure = 300; tpure <= 3500; tpure += 10){
      rhoNow = p_atm*Wbar/(PhysicalConstants<double>::IdealGasConstant*tpure);
      //cout << "rho = "<< rhoNow << endl; 
      // mavgcond = MavgCond.compute_mixture_averaged_conductivity(rhoNow,p_atm,tpure,XO3);
      //adonis_assert(is_well_defined(mavgcond));
      //! not that for a single species ,this diffusion coeff is undefined
      MixDiff.compute_mixture_averaged_diffusion_coefficients(p_atm,tpure,XO3);
      //print_all(MixDiff.coefficients(),MixDiff.coefficients()+3);


      of << tpure << "  " << Mav.compute_mixture_averaged_viscosity(tpure,XO3) << "  " << MavgCond.compute_mixture_averaged_conductivity(rhoNow,p_atm,tpure,XO3)<< "   " << MixDiff.coefficients()[o2] << endl;
    }
    of.close();
    
    int out = 2;  //1 = eps, 2 = tex
    string togx, togxlambda, togxdiff, cl;
    
    if(out == 1)
      togx = "set terminal postscript eps color enhanced \"Helvetica\" 20; set output \"mixavgviscosity.eps\"; set xlabel \"temperature in K\"; set ylabel \"mixture viscosity {/Symbol m}(T)  in kg/(ms)\";";
    //!CAUTION: use single quotes ('') here to propagate the backslash.
    //!         with double quotes ("") use \\ for it (NOT working for me)  
    //!         FOR ME: use '' in conjunction with '' !!!
    if(out == 2){
      togx = "set mxtics 2; set xtics (500,1500,2500,3500);set size 0.625,0.6; set terminal latex; set output \"mixavgviscosity.tex\"; set format xy \"$%g$\"; set xlabel '$T \\ (\\si{K})$'; set ylabel '\\rotatebox[origin=c]{90}{$\\mu(T) \\ (\\si{kg.(m.s)^{-1}})$}'; unset key;";

       togxlambda = "set xtics (500,1500,2500,3500);set size 0.6,0.6; set terminal latex; set output \"mixavgconductivity.tex\"; set format xy \"$%g$\"; set xlabel '$T \\ (\\si{K})$'; set ylabel '\\rotatebox[origin=c]{90}{$\\lambda(T) \\ (\\si{W.(m.K)^{-1}})$}'; unset key;";

       togxdiff = "set xtics (500,1500,2500,3500);set size 0.6,0.6; set terminal latex; set output \"mixavgconductivity.tex\"; set format xy \"$%g$\"; set xlabel '$T \\ (\\si{K})$'; set ylabel '\\rotatebox[origin=c]{90}{$\\mathcal{D}^{\\mathrm{mix}}(T) \\ (\\si{m.(s)^{-1}})$}'; unset key;";

      cl = ";set output; set terminal pop;";
    }

    GNUplot plotmu;
    //set xtic 300,1000,3500;set xtics rotate by -45;
    plotmu(togx + "set xrange[300:3500]; plot \"Mu.dat\" using 1:2 with lines lw 3 notitle"+cl);

    GNUplot plotlambda;
    plotlambda(togxlambda + "set xrange[300:3500]; plot \"Mu.dat\" using 1:3 with lines lw 3 notitle\"" +cl);

    // GNUplot plotdiffusion;
    // plotdiffusion(togxdiff + "set xrange[300:3500]; plot \"Mu.dat\" using 1:4 with lines lw 3 notitle\"" +cl);

    GNUplot sq;
    sq("plot \"Mu.dat\" using 1:(sqrt($1)) with lines lw 3");
 

    MyVec<double> vxz(8), vyz(7);
    vxz[7] = 2200;
    
    vyz <<= 0.7, 0.14, 0.001, 0.12, 0.432, 0.89, 0.115;

    UnrollLoop<0,7>::assign(vxz.begin(),vyz.begin());
    cout << "vxz = "<< vxz << endl;
    

    MyVec<double> xtest(3);
    xtest <<= 1,2,3;
    double pert = 0.5;
    double xscalar = 3.25;
    
    cout << "\"perturbation\"   = " << perturbation(xtest,pert) <<endl;
    cout << "\"scalar perturb\" = " <<perturbation(xscalar,pert) << endl;

    cout << endl << "THERMAL DIFFUSION RATIO"<<endl;
    
    const bool WHXH = false;
    
    typedef ThermoData4Mechanism<double,10> ConH2MechType;
    //!CAUTION: always assign the right transport() array ;)
    ThermalDiffusionRatio<ConH2MechType,FULL,WHXH> Theta(ConH2MechType::transport());
    
    typedef CppAD::AD<double> SomeType;
    
    SomeType tpur = 1500;
    
    string str[10] = {"O", "O2", "H", "OH", "H2", "HO2", "H2O2", "H2O", "N2","AR"};
    SomeType xfuck[10] = {0., 0.232, 0., 0., 0., 0., 0., 0.001, 0.755, 0.012};
    MyVec<SomeType> Xfuck(xfuck,xfuck+10);
    Theta.compute_thermal_diffusion_ratios(tpur,Xfuck);

    cout << setprecision(7);

    for(int k = 0; k < 10; ++k)
      cout  << "Theta_"<<str[k] << " = " << Theta[k] << endl;
    

    cout << endl<< "---------------- REDUCED FORM -- only for the sake of completeness ---------"<<endl;
    MyVec<SomeType> RedXfuck(ConH2MechType::rednspec);
    RedXfuck <<= 0.3, 0, 0.7, 0.1;
    ThermalDiffusionRatio<ConH2MechType,REDUCED,WHXH> RedTheta(ConH2MechType::transport());
    RedTheta.compute_thermal_diffusion_ratios(tpur,RedXfuck);

    for(int k =0 ; k < ConH2MechType::rednspec; ++k)
      cout << "RedTheta_"<<k << " = " << RedTheta.ratios()[k] << endl; 


    cout << endl <<  "TestMe:" << endl;
    bool Cool[6];

   UnrollLoop<0,6>::assign_value(Cool,false);

   for(int i = 0; i < 6; ++i)
     cout << Cool[i] << ", "; 
    
    // //!works
    // cout << endl << "TEst:"<<endl;
    // double xc1 = 0, xc2 = 1.75, molec = 5., molec2 = 0.45;
    // assign_prop_zero(false, xc1,xc2,molec);
    // cout << "xc1 = "<< xc1 <<endl;
    // assign_prop_zero(true, xc1,xc2,molec);
    // cout << "xc1 = "<< xc1 <<endl;
    // cout << "same with above molec:"<<endl;
    //  assign_prop_zero(false, xc1,xc2,molec2);
    // cout << "xc1 = "<< xc1 <<endl;
    // assign_prop_zero(true, xc1,xc2,molec2);
    // cout << "xc1 = "<< xc1 <<endl;
      
//RedTheta[k] << endl; //works, of course
      


    // int idxy = 0;
    // cout << endl << "ekB("<<idxy << ") = "<< Theta.ekB(idxy) << endl;

    ////works
    // double xb1 = 5.;
    // double xb2 = 12.5;
    // cout << "x^(-1) = "<< ntimes<-1>(xb1) << endl;
    // cout << "x^(-1) = "<< ntimes<-1>(xb2) << endl;
    
    // double yd1[] = {0.5, 1.2, -4.6};
    // double yd2[] = {-1.5, 0.095, 0.75};
    
    // MyVec<double> Yd1(yd1,yd1+3);
    // MyVec<double> Yd2(yd2,yd2+3);

    // double su1 = 0;
    // for(int i = 0; i < 3; ++i)
    //   su1 += yd1[i]/yd2[i];

    // cout << "su1 = "<< su1 << endl;
    // cout << UnrollLoop<0,3>::sum_product_xtended<1,-1>(Yd1.begin(),Yd2.begin()) << endl;

    //********************************************************************
   //*************************** END *************************************
   //*********************************************************************
    

  
    
    // cout <<endl<< "Test compile time array:"<<endl;
    // typedef Element<0, Element<3, Element<6, Element<7> > > > Array;
    // cout << AccessElement<Array, 0>::result << endl;
    // cout << AccessElement<Array, 1>::result << endl;
    // cout << AccessElement<Array, 2>::result << endl;
    // cout << AccessElement<Array, 3>::result << endl;


    // cout<< endl<< "Test sum<P>:"<<endl;
    // double sqs[] = {0.5, 1.25, -0.75, 3., -2.5};
    // ExprTmpl::MyVec<double> sqsuvec(sqs,sqs+5);

    // cout << "sum_i=1^5 x_i² = " << setprecision(5) << UnrollLoop<0,5>::sum<2>(sqsuvec) << endl;   

   return 0;
}
