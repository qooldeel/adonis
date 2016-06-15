#include <iostream>
#include <cmath>
#include <complex>

#include <typeinfo>

#include "../want_2_use.h"  //include it at the very top

#include "exprvec.hh"
#include "../common/tempusfugit.hh"
#include "../misc/misctmps.hh"

#include "../dunestuff/fmatrix.hh"

#include "../common/globalfunctions.hh"
#include "../common/universalconstants.hh"
#include "../common/numericaltraits.hh"


//#include <cppad/cppad.hpp>

//test operations the FORTRAN,Java,C,Perl,Matlab(c),...style
template<class V>
inline V operator*(const V& v, const V& w){
  adonis_assert(v.size() == w.size());
  
  V res(v.size());                   //temporary object 
  
  for(size_t i = 0; i < v.dim; ++i)    //for loop
    res[i] = v[i]*w[i];

  return res;
}


//a function object 
template<class T>
class FunnyFunction{
public:
  typedef Adonis::ExprTmpl::MyVec<T> VType;

  FunnyFunction(int n = 0):dest_(n){}
  
  size_t size() const{return dest_.size();}

  template<class W>
  VType& operator()(const W& arg){

    //==== vector valued function ====
    dest_[0] = arg[0] + arg[1]*arg[0];
    dest_[1] = arg[0] + arg[1]*arg[2];
    dest_[2] = arg[0] + 2.5*arg[1] - 0.5*arg[2];

    //================================
    
    return dest_;
  }

private:
  VType dest_;

};


template<class T>
class FullFunc{
public:
   typedef Adonis::ExprTmpl::MyVec<T> VType;

  FullFunc(unsigned n = 0):f_(n){}

  unsigned size() const{return f_.size();}
  unsigned dim() const{return f_.size();}
  unsigned domain_dim() const{return 6;}

  template<class X>
  VType& operator()(const T& time, const X& x){
    //========================================================================
    f_[0] = x[0] + std::sin(x[1]*x[2] - x[3])*x[4] - x[5]*x[5]*x[3];
    f_[1] = x[1] + 2.*x[5] + std::cos(3.*x[0]);
    f_[2] = x[2]*x[5];
    f_[3] = std::cos(x[4]) - x[0];
    f_[4] = -x[5] + 3.5*x[1]*x[3];
    f_[5] = -x[4] + x[2]*4.75;
    //========================================================================
    return f_;
  }

private:
  VType f_;
};

//! a reduced function object \f[ f = \begin{bmatrix}  f_1(x_1, x_2, \ldots, x_6) \\ f_2(x_1, x_2, \ldots, x_6) \end{bmatrix}. \]
template<class T>
class RedFunc{
public:
   typedef Adonis::ExprTmpl::MyVec<T> VType;

  RedFunc(unsigned n = 0):r_(n){}

  unsigned size() const{return r_.size();}
  unsigned dim() const{return r_.size();}
  unsigned domain_dim() const{return 6;}

  template<class X>
  VType& operator()(const T& time, const X& x){
    //=========================================================================
    r_[0] = x[0] + std::sin(x[1]*x[2] - x[3])*x[4] - x[5]*x[5]*x[3];
    r_[1] = -x[5] + 3.5*x[1]*x[3];
    //=========================================================================
    return r_;
  }

private:
  VType r_;
};


/* 
 * \brief Just for checking the order of evaluation, i.e. from which side
 * an expression of the following form  is evaluated:
 *  \f[  x = a + b + c + d + e +f \]
 *
 * Apparently, this is done from <I> left to right </I>
 *
 *\code
    cout << endl<< "Check order of evaluation:"<< endl;
    OrderOfEvaluation<PrecType> ooe(Pa);
 
    PrecType rslt = ooe[0] + ooe[1] + ooe[2] + ooe[3] + ooe[4] + ooe[5];

    //likewise
    PrecType vor = PrecType();
    vor += (ooe[0] + ooe[1] + ooe[2] + ooe[3] + ooe[4] + ooe[5]);

    cout << endl<<"Example that should produce the same order:"<<endl;
    PrecType ltsr = PrecType();
    for(size_t i = 0; i < ooe.size(); ++i)
      ltsr += ooe[i];
 *\endcode
 *
 * The result looks like evalutate w.r.t. index 0
 *                       evalutate w.r.t. index 1
 *                       evalutate w.r.t. index 2
 *                                 .....
 * valgrind doesn't change anything ;)
*/
template<class T>
class OrderOfEvaluation{
public:
  typedef Adonis::ExprTmpl::MyVec<T> Type;

  OrderOfEvaluation(size_t dim = 0):v_(dim), dim_(dim){};

  OrderOfEvaluation(const Type& w):v_(w.size()), dim_(w.size()){
    for(size_t i = 0; i < w.size(); ++i)
      v_[i] = w[i];
  }

  size_t size() const {return dim_;}

  T& operator[](int i) {
    std::cout << "evalutate w.r.t. index "<<i<< std::endl;
    return v_[i];
  }

private:
    Type v_;
    size_t dim_;
};


//============= MAIN PROGRAM ==============================
using namespace std;
using namespace Adonis;
using namespace ExprTmpl;


int main(){

  double a[] = { 1, 2, 3, 4, 5};
  double b[] = { 1, 2, 3, 4, 5};
  double c[] = { 1, 2, 3, 4, 5};
  
  const size_t dim = 5;
  MyVec<double> X(dim), Y(dim), Z(dim), W(dim), Def;
  cout << "X.size() = "<< X.size() << endl;
  cout << "Some output..."<<endl;
  //cout << "X = "<< X <<endl;

  // //fill the vectors
  for(size_t i = 0; i < dim; ++i){
     cout << i << ") "<<endl;
     //cout << "X["<<i<<"] = "<< X[i] <<endl;
     X[i] = a[i];
     Y[i] = b[i];
     Z[i] = c[i];
  }
  
  cout<<"Test copy-constructor:"<<endl;
  MyVec<double> CC(X);
  for(size_t i = 0; i< dim; i++)
    cout<<CC[i]<<"  ";
 

  cout<<endl<<"Test normal assignment:"<<endl;
  MyVec<double> CA;
  CA = Y;
  for(size_t i = 0; i< dim; i++)
    cout<<CA[i]<<"  ";
 
  cout<<endl<<"Test iterator based constructor:"<<endl;
  std::vector<double> stlvec(c, c +dim);
  MyVec<double> FromIt(stlvec.begin(), stlvec.end());
  for(size_t i = 0; i< dim; i++)
    cout<<FromIt[i]<<"  ";
  cout<<endl;
  
  cout<<"ITERATORS of MyVec:"<<endl;
  cout<<"Non-const case:"<<endl;
  for(MyVec<double>::iterator it = Z.begin(); it != Z.end(); ++it)
    cout<<*it<<"  ";
  cout<<endl<<"Const case:"<<endl;
  for(MyVec<double>::const_iterator cit = Z.begin(); cit != Z.end(); ++cit)
    cout<<*cit<<"  ";
  cout<<endl;

  cout<<endl;
  cout<<"=================== TEST OPERATIONS: ======================"<<endl;
  
  W = X + Y + Z; //works ?
  cout<<W<<endl;
  
 
  W=(X+Y+Z)*3.; //works
  cout<<W<<endl;
  W = X + Y + Z; //works
  cout<<W<<endl;
  W = 1.*X + 0.5*Y; //works
  cout<<W<<endl;
  W = 1.*X + 0.5*Y + Z;  //works
  cout<<W<<endl;
  W = 1.*X + 0.5*Y + Z + 2.0*Y; //works

  
 
  W = X*0.75 - 4.5*Y - X*Y*Z + 1.35*X*X;     //works
  cout<<W<<endl;
  
  W += X + 0.5*Y;   //works
  cout<<W<<endl;
  W -= 1.*X;    //works
  cout<<W<<endl;
  W += X;      //works
   cout<<W<<endl;
  W -= Y+2.*Z;
   cout<<W<<endl;
  W *= 10.75*(X+Y);  //works
   cout<<W<<endl;
  
  W *= (X+Y)*1.;    //works
  cout<<W<<endl;
  
                    
  W = X + 0.5*Y;        //works
  cout<<W<<endl;
  W *= -3.*(X+Y)*(Y+Z);
  cout<<W<<endl;
  W *=1.;
  cout<<W<<endl;
  W = -3.75*(X+Y);
  cout<<W<<endl;
  W *= -1.*(X+Y);
  cout<<W<<endl;
  W += 0.5*X + Z; //works
  cout<<W<<endl;
  W += -1.75*(X +1.0*Y*Z); //works
  cout<<W<<endl;
  W = -3.*(X+Y) - Z*-1.5;
  cout<<W<<endl;

  
  cout<<endl<<endl<<"Test size and capacity:"<<endl;
  MyVec<double> V;  //empty vec
  cout<<" Empty size = "<<V.size()<<endl;
  const size_t dres = 7;
  V.reserve(dres);     //reserve space
  cout<<" Size after reserve = "<<V.size()<<endl;
  
  for(int i = 0; i < int(dres); ++i){
    V.push_back(i*0.2+3.075);
    cout<<"#of elements which can be pushed without causing memory reallocation = "<<V.capacity()-V.size()<<endl;
  }

  cout<<endl<<" V = "<<V<<endl;
  cout<< " size() = "<<V.size()<<endl;
  cout<<" capacity() = "<<V.capacity()<<endl;

  cout<<"O.k., now operator[] without range violation"<<endl;
  cout<< W[3]<<endl;
  
  
 //works when DNDEBUG is commented in your Makefile
  //cout<<"...and with obvious violation:"<<endl;
  //int fail = 18;
  //cout<<W[fail]<<endl;
  
  cout << "test operator ==" <<endl;
  MyVec<double> mc1(3), mc2(4), mc3(3),mc4(3);

  mc1[0] = 1; mc1[1] = 2; mc1[2] = -1; 
  mc2[0] = 1; mc2[1] = 2; mc2[2] = -1; mc2[3] = 0.5; 
  mc3[0] = 1; mc3[1] = 2; mc3[2] = -1; 
  mc4[0] = 1; mc4[1] = 0.5; mc4[2] = -1.5; 
 

  cout<<endl<<"are vec's equal?:"<<endl;
  cout << (mc1 == mc2) << " (NO -- different sizes)"<<endl;
  cout << (mc1 == mc3) << " (YES -- same entries)"<<endl;
  cout << (mc1 == mc4) <<  " (NO -- different entries)" <<endl;
  cout << (mc2 == mc2) << " (YES -- identity)" <<endl;
 
  //The bad stuff! Expressions will have different sizes!! 
  cout<<endl<<"---NOTE: Operations on differently sized MyVec's will fail, e.g. "<<endl<<endl;
  size_t b1 = 6, b2 = 3;
  double bad[] = {-0.5, 0.5, 1, 1.5, -1.75, 2.75};
  double bad2[] = {-1., 3.5, -0.6};

  MyVec<double> Bad(bad,bad + b1), Bad2(bad2, bad2 + b2);
  cout<<" Bad = "<<Bad<<endl;
  cout<<" Bad2 = "<<Bad2<<endl;
 
  cout << endl<< "===== TEST MATRIX-VECTOR / VECTOR-MATRIX EXPRESSION TEMPLATE OPERATIONS ====" <<endl;
  double smurf[] = {1,3,9,6};
  MyVec<double> v(smurf, smurf+4);
  
  Dune::FieldMatrix<double,3,4> M;
  M[0][0] = 2; M[0][1] = 5;  M[0][2] = -3; M[0][3] = 1;
  M[1][0] = 0; M[1][1] = 2;  M[1][2] = 3;  M[1][3] = 2;
  M[2][0] = 1; M[2][1] = -1; M[2][2] = -3; M[2][3] = 6;
  
  cout <<"M = "<<endl<< M << endl;

  cout <<"v = "<<endl<< v <<endl;

  MyVec<double> vtimesM(3);  //necessary to initialise!
  
  vtimesM = M*v;

  cout << "M*v = "<<endl;
  cout<<vtimesM <<endl;
  
  double garg[]={-4,2,-11};
  MyVec<double> w(garg, garg+3);
  cout<<"w="<<endl<< w <<endl;
  

  MyVec<double> Mtimesw(4);
  Mtimesw = w*M;
  cout << "w*M = "<<endl;
  cout << Mtimesw <<endl;

  cout<<endl<<"2nd example:"<<endl;
  double cfi = 1;
  Dune::FieldMatrix<double,3,3> Trois;
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Trois[i][j] = cfi++;
    }
  }

  cout<<Trois<<endl;
  MyVec<double> Bo(3), Ob(3);

  Bo = Trois*w;
  Ob = w*Trois;

  cout << "Trois*w = "<<endl<<Bo<<endl;
  cout << "w*Trois = "<<endl<< Ob<<endl; 
  cout << endl << "====== END MATRIX-VECTOR ========"<<endl;

  cout<<endl<<"An update the swift way:" <<endl;
  double yp[] = {0.5, -0.5, 1};
  double yn[] = {0.75, 0., -2.5};
  double fev[] = {-3.45, 2.55, 4.5};

  MyVec<double> Yp(3),Yn(3),Fev(3),G(3);
  for(unsigned i = 0; i<3; ++i){
    Yp[i] = yn[i];
    Yn[i] = yp[i];
    Fev[i] = fev[i];
  };
  double h = 0.05;

  G = -Yn + Yp + h*Fev;
  cout<<"G = "<<G<<endl;
    
  cout<< "G(1) = "<<G(1)<<endl;

  cout << endl << "Test with function object that posses an operator() returning a MyVec<T>& :"<<endl;
  const int sz = 3;
  FunnyFunction<double> FF(sz);
  
  double fi[] = {2,-1,3};
  MyVec<double> arg(fi, fi+sz);
  cout << FF(arg) <<endl;

  cout<<"FF.size() = "<<FF.size()<<endl;
  MyVec<double> funres(sz);
  funres = FF(arg)*Trois;
  cout << "FF(arg)*Trois = "<<endl<< funres <<endl;

 

  cout<< endl<< "FUNCTION EXPRESSIONS"<<endl;
  cout<< " -with differently sized input"<<endl;
  typedef MyVec<double> ExprVecType;
  ExprVecType  x1(6), y1(6), z1(6), 
    rres(2), 
    fres(6);
  for(unsigned i = 0; i < 6; ++i){
    x1[i] = (i+1);
    y1[i] = 0.5*(i+1);
    z1[i] = 0.1*(i+1);
  }


  //Adonis::ExprTmpl::Xpress<double,ExprVecType::const_iterator> x1pr(x1.begin()), y1pr(y1.begin()), z1pr(z1.begin());
  
  double time = 2.65;

  //work both
  RedFunc<double> rfun(2);
  FullFunc<double> ffun(6);

  rres = rfun(time,x1 + 2.*y1 + x1*z1*0.75
	      //x1pr + 2.*y1pr
	      );

  fres = ffun(time, x1 + 2.*y1
	      //x1pr + 2.*y1pr

	      );
  
  cout<<"rres = "<<rres<<endl;
  cout<<"fres = "<<fres<<endl;

  Dune::FieldMatrix<double,2,6> B_T(0.);
  B_T[0][0] = 1.; B_T[1][4] = 1.;

  cout << "B_T = "<< endl<< B_T <<endl;
  
  MyVec<double> Mprod(2);
  Mprod =  B_T*ffun(time, x1 + 2.*y1
 		    //x1pr + 2.*y1pr
		    );

  cout<< "B^T·ffun(time,expression) = "<< Mprod <<endl;
			   

  //cout<<endl<<"Further investigations in 'Xpress:'"<<endl;
  MyVec<double> K1(3), K2(3), K3(3), FK(3);

  for(unsigned i = 0; i < 3; ++i){
    K1[i] = i+1;
    K2[i] = 0.25*(i+1);
    K3[i] = -1.*(i+1);
  }

  cout<<"K1 = "<<K1 <<endl;
  cout<<"K2 = "<<K2 <<endl;
  cout<<"K3 = "<<K3 <<endl;
  // Adonis::ExprTmpl::Xpress<double,ExprVecType::const_iterator> k1(K1.begin()), k2(K2.begin()), k3(K3.begin());

  FK = K1 - 2.*K2 + K3*5.;
  cout << "FK = "<<FK <<endl;


  cout << "Test operator-:"<<endl;
  MyVec<double> ov1(3), ov2(3), ov3(3), Ov(3);

  for(unsigned i = 0; i < 3; ++i){
    ov1[i] = i+1;
    ov2[i] = 0.5*(i+1);
    ov3[i] = -0.25*(i+1);
  }

  Ov = -ov1 + 2.5*ov2*ov3;

  cout << "Ov = "<<Ov <<endl;

  cout<<"======================================================================"<<endl;
  
  const int prec = 17;

  cout<< setprecision(prec) <<endl; //PRECISION

  typedef double PrecType;
  cout << "sizeof(PrecType) = "<< sizeof(PrecType) <<endl;

  cout <<endl<< " KAHAN COMPENSATES SUMMATION TEST:" <<endl; 

  const unsigned pdim = 3; 
  MyVec<PrecType> XK(pdim),YK(pdim),ZK(pdim),WK(pdim),
    AK(pdim), BK(pdim), CK(pdim), DK(pdim);

  XK[0] =  3.50000000000089278;
  XK[1] = -1.00000000000000265;
  XK[2] =  0.75000000000000115;

  YK[0] =  2.00000006000158661;
  YK[1] =  0.00000000000000909;
  YK[2] = -0.25000000000000128;

  ZK[0] =  0.78932890129012190;
  ZK[1] = -4.00671007821947582;
  ZK[2] = -0.17631801229931225;
  
  AK[0] = -4.37819345347801721;
  AK[1] = -0.00321745123924099;
  AK[2] =  0.83129471500099913;

  BK[0] =  0.16780000125420912;
  BK[1] = -0.71367647361311109;
  BK[2] =  2.00782331846111754;

  CK[0] = -1.00471235870621079;
  CK[1] = -0.67460048287360248;
  CK[2] =  2.64238889009914552;


   PrecType fac1 = 0.3728011324891249,
    fac2 = -1.6910892499788139,
    fac3 = 0.07180003214110991;
  

  cout << "Test stuff MyVec = Expression:"<<endl;
  MyVec<PrecType> Te(6);
  Te[0] = fac1*XK[0]*YK[0]; 
  Te[1] = YK[0]; 
  Te[2] = ZK[0];
  Te[3] = AK[0]*ntimes<2>(fac2); 
  Te[4] = fac3*BK[0]; 
  Te[5] = CK[0]; 

  //unrolled loop for KAHAN
  PrecType sum = 0.,
    corr = 0.,
    y, t;

  
  sum = Te[0];
  // cout << "sum = " << sum << "   corr = "<< corr << endl; 

  y = Te[1] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;
  
  //cout << "sum = " << sum << "   corr = "<< corr << endl; 

  y = Te[2] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;


  y = Te[3] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;

  y = Te[4] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;

  y = Te[5] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;


  PrecType unrolledkah = sum;

  //unrolled loop for BABUSKA-KAHAN
  sum = Te[0];
  corr = 0.;

  t = sum + Te[1];
  if(Abs(sum) >= Abs(Te[1]))
    corr += ((sum-t)+Te[1]);
  else
    corr += ((Te[1]-t)+sum);

  sum = t;

  t = sum + Te[2];
  if(Abs(sum) >= Abs(Te[2]))
    corr += ((sum-t)+Te[2]);
  else
    corr += ((Te[2]-t)+sum);

  sum = t;

  t = sum + Te[3];
  if(Abs(sum) >= Abs(Te[3]))
    corr += ((sum-t)+Te[3]);
  else
    corr += ((Te[3]-t)+sum);

  sum = t;

  t = sum + Te[4];
  if(Abs(sum) >= Abs(Te[4]))
    corr += ((sum-t)+Te[4]);
  else
    corr += ((Te[4]-t)+sum);

  sum = t;

   t = sum + Te[5];
  if(Abs(sum) >= Abs(Te[5]))
    corr += ((sum-t)+Te[5]);
  else
    corr += ((Te[5]-t)+sum);

  sum = t;

  sum += corr; //update sum at the very end

  PrecType unrolledbabusk = sum;

  cout << "normal = "<< conventional_sum(Te) <<endl;
  cout << "direct = " << Te[0] + Te[1] + Te[2] + Te[3] + Te[4] + Te[5] <<endl;
  cout << "KAHAN  = "<< kahan_algorithm(Te) <<endl;
  cout << "dikah  = "<< unrolledkah << endl;
  cout << "BABUSK = "<< kahan_babuska_algorithm(Te) <<endl;
  cout << "dirBAB = "<< unrolledbabusk << endl;

  
  //WK[0] = -5.072183721909000278;
  //WK[1] =  0.123280012381290012;
  //WK[2] =  1.777903431151378993;
    

  WK = fac1*XK*YK + YK + ZK + AK*ntimes<2>(fac2) + fac3*BK + CK;

  cout << endl<< " WK[0] = "<< WK[0] <<endl;

  
  cout << endl << "Test 'MyVec += Expression':"<<endl;
  MyVec<PrecType> Ka(5);
  
  Ka[0] = AK[0];
  Ka[1] = BK[0];
  Ka[2] = CK[0];
  Ka[3] = XK[0];
  Ka[4] = YK[0];
  
  MyVec<PrecType> Sa(pdim);
  Sa[0] =  0.93812492374239479;
  Sa[1] = -0.53245239081439753;
  Sa[2] =  1.47123784836145554;
 
  MyVec<PrecType> Pa(6);
  
  Pa[0] = Sa[0];
  for(int i = 1; i < 6; ++i)
    Pa[i] = Ka[i-1];
  
 
  sum = Pa[0];
  corr = 0.;
  

  y = Pa[1] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;

  y = Pa[2] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;

  y = Pa[3] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;

  y = Pa[4] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;

  y = Pa[5] - corr;
  t = sum + y;
  corr = (t - sum) - y;
  sum = t;


  unrolledkah = sum;

  //unrolled loop for BABUSKA-KAHAN
  sum = Pa[0];
  corr = 0.;

  t = sum + Pa[1];
  if(Abs(sum) >= Abs(Pa[1]))
    corr += ((sum-t)+Pa[1]);
  else
    corr += ((Pa[1]-t)+sum);

  sum = t;

  t = sum + Pa[2];
  if(Abs(sum) >= Abs(Pa[2]))
    corr += ((sum-t)+Pa[2]);
  else
    corr += ((Pa[2]-t)+sum);

  sum = t;

  t = sum + Pa[3];
  if(Abs(sum) >= Abs(Pa[3]))
    corr += ((sum-t)+Pa[3]);
  else
    corr += ((Pa[3]-t)+sum);

  sum = t;

  t = sum + Pa[4];
  if(Abs(sum) >= Abs(Pa[4]))
    corr += ((sum-t)+Pa[4]);
  else
    corr += ((Pa[4]-t)+sum);

  sum = t;

   t = sum + Pa[5];
  if(Abs(sum) >= Abs(Pa[5]))
    corr += ((sum-t)+Pa[5]);
  else
    corr += ((Pa[5]-t)+sum);

  sum = t;

  sum += corr; //update sum at the very end

  unrolledbabusk = sum;

  cout << "normal = "<< conventional_sum(Pa) <<endl;
  cout << "direct = "<< Pa[0] + Pa[1] + Pa[2] + Pa[3] + Pa[4] + Pa[5]  << endl;
  cout << "KAHAN  = "<< kahan_algorithm(Pa) <<endl;
  cout << "dirKAH = "<< unrolledkah  <<endl;
  cout << "BABUSK = "<< kahan_babuska_algorithm(Pa) <<endl;
  cout << "dirBAB = "<< unrolledbabusk <<endl<<endl<<endl;

  //cout << Sa[0] + AK[0] + BK[0] + CK[0] + XK[0] + YK[0] << endl;
  Sa += (AK + BK + CK + XK + YK);    cout << " Sa[0] = "<< Sa[0] <<endl;
  

  cout << endl<< "Check order of evaluation:"<<endl;
  OrderOfEvaluation<PrecType> ooe(Pa);

  PrecType rslt = PrecType();
  rslt += (ooe[0] + ooe[1] + ooe[2] + ooe[3] + ooe[4] + ooe[5]);


  
   // o.k.
  cout << endl << "add 2 numbers: ";
  double ad = 1.5, 
    hst = 0.075;
  cout << IntelligentArithmeticOperation<XPCTSettings::floatingPointAdditionMethod>::add(ad,hst) << endl;
  IntelligentArithmeticOperation<XPCTSettings::floatingPointAdditionMethod>::pluseq(ad,hst);
  cout << "ad (after pluseq) = "<<ad <<endl;


  cout <<endl<< "Test with FunnyFunction"<<endl;
  FunnyFunction<PrecType> fun(pdim);

  MyVec<PrecType> nleq(pdim);
  PrecType step = 0.00125;

  nleq = -1.*AK + BK + step*fun(0.5*(AK + BK));

  cout << "nleq = "<< nleq <<endl;
  
  cout << "Brute force evaluation:"<<endl;
  MyVec<PrecType> farr(pdim), karr(pdim);
  farr = AK + BK;
  farr *= 0.5;
  karr = -1.*AK + BK + step*fun(farr);
  cout << "karr = "<< karr <<endl;
  
  
  cout<<endl << "Test -() "<<endl;
  double v1[] = {1,2,3};
  double v2[] = {-2,0.5,1.5};
  double v3[] = {1.75, -0.5, 6};
  MyVec<double> V1(v1,v1+3), V2(v2,v2+3), V3(v3,v3+3), VRES(3);
  
  VRES = -(0.75*V1 - V2*1.5 + V3);

  cout << "VRES = "<< VRES << endl;
  


  cout<<endl << "Test more than two arguments at the same time:"<<endl;
  cout << typeid(NumericalTypePromotion<long double, NumericalTypePromotion<double,std::complex<float> >::ReturnType>::ReturnType).name() <<endl;
  cout << typeid(NumericalTypePromotion<long double, NumericalTypePromotion<CppAD::AD<double>,std::complex<float> >::ReturnType>::ReturnType).name() <<endl;
  

  cout << endl<< "Test fancy assignment operator:"<<endl;
  MyVec<double> vtest(5);
  cout << "vtest.size() = "<<vtest.size()<<endl; 
  vtest <<= 0.5, -0.75, 
    0., -1.5, 2;
  
  cout << "vtest = " << vtest << endl;

  MyVec<double> vtest1(4);
  cout << "vtest1.size() = "<<vtest1.size()<<endl; 
  vtest1 <<= 4,6,
    7,8;

  cout << "vtest1 = "<<vtest1<<endl;
  
  cout<< "1 entry example:"<<endl;
  MyVec<double> vtest2(1);
  vtest2 <<= 9.65;
  cout << "vtest2 = "<< vtest2 << endl;

  cout << endl<< "Test more difficult example:"<<endl;
  double pi = UniversalConstants<double>::Pi;

  MyVec<double> vtest3(3);

  vtest3 <<= -pi/2., pi, sqrt(0.75*pi);

  cout << "vtest3 = "<<vtest3 <<endl;

  ////! works
  // cout << "Test base 2 logarithm:" <<endl;
  // double dval = 2.45;
  // cout << "std::log2(dval) = "<< log2(dval) << endl;
  // cout << "Log2(dval) = "<< Log2(dval) << endl;
  // complex<double> z(-3,-4), zzero(0,0);
  // cout << "Log2(z) = " << Log2(z) << endl;
  // cout << "Log2(zzero) = " << Log2(zzero) << endl;
  // CppAD::AD<double> adval = 0.5;
  // cout << "Log2(adval) = " << Log2(adval)<<endl;
  // cout << "log2(0.5) = "<< log2(0.5) << endl;
  //cout << typeid(NumericalTypePromotion<CppAD::AD<double>,std::complex<float> >::ReturnType).name() <<endl;


  /* //works as far as implemented............
  cout<<endl<<"Check NumericalTraits:"<<endl;
  cout << typeid(NumericalTypePromotion<float,float>::ReturnType).name() <<endl;
 
  cout << typeid(NumericalTypePromotion<double,float>::ReturnType).name() <<endl;
  
  cout << typeid(NumericalTypePromotion<float,long double>::ReturnType).name() <<endl;
  cout << typeid(NumericalTypePromotion<double,long double>::ReturnType).name() <<endl;
  cout << typeid(NumericalTypePromotion<long double,double>::ReturnType).name() <<endl;
    
  cout <<endl << "complex stuff:"<<endl;
  cout << typeid(NumericalTypePromotion<complex<float>,complex<float> >::ReturnType).name() <<endl;
		 
		 cout << typeid(NumericalTypePromotion<complex<double>,complex<float> >::ReturnType).name() <<endl;
  
		   cout << typeid(NumericalTypePromotion<complex<float>,complex<long double> >::ReturnType).name() <<endl;
		   cout << typeid(NumericalTypePromotion<complex<double>,complex<long double> >::ReturnType).name() <<endl;
		   cout << typeid(NumericalTypePromotion<complex<long double>,complex<double> >::ReturnType).name() <<endl;
    

				  cout<<endl<<"Mixed complex - normal types:"<<endl;
  cout << typeid(NumericalTypePromotion<complex<float>,double>::ReturnType).name() <<endl;
  
    cout << typeid(NumericalTypePromotion<complex<long double>,double>::ReturnType).name() <<endl;
     cout << typeid(NumericalTypePromotion<complex<float>,complex<float> >::ReturnType).name() <<endl;
    cout << typeid(NumericalTypePromotion<complex<float>,complex<double> >::ReturnType).name() <<endl;
  				  

				  
				  cout<<endl<<"CppAD stuff"<<endl;
    cout << typeid(NumericalTypePromotion<CppAD::AD<float>,double>::ReturnType).name() <<endl;
  
    cout << typeid(NumericalTypePromotion<CppAD::AD<long double>,double>::ReturnType).name() <<endl;
     cout << typeid(NumericalTypePromotion<CppAD::AD<float>,CppAD::AD<float> >::ReturnType).name() <<endl;
		    cout << typeid(NumericalTypePromotion<complex<CppAD::AD<float> >,CppAD::AD<double> >::ReturnType).name() <<endl;
		    cout << typeid(NumericalTypePromotion<complex<CppAD::AD<double> >,CppAD::AD<float> >::ReturnType).name() <<endl;
		    cout << typeid(NumericalTypePromotion<CppAD::AD<complex<float> > , CppAD::AD<complex<float> > >::ReturnType).name() <<endl;	   
    cout << typeid(NumericalTypePromotion<long double , CppAD::AD<complex<float> > >::ReturnType).name() <<endl;	 
		    cout << typeid(NumericalTypePromotion<std::complex<double> , CppAD::AD<complex<float> > >::ReturnType).name() <<endl;
	 
		    cout << typeid(NumericalTypePromotion<CppAD::AD<complex<float> >, double>::ReturnType).name() <<endl;	    
		    cout << typeid(NumericalTypePromotion<CppAD::AD<complex<float> >, complex<long double> >::ReturnType).name() <<endl;	    
				   cout << typeid(NumericalTypePromotion<complex<CppAD::AD<float> > , complex<long double> >::ReturnType).name() <<endl;	    
				   cout << typeid(NumericalTypePromotion<complex<CppAD::AD<float> > , CppAD::AD<long double> >::ReturnType).name() <<endl;
				   //  cout << typeid(NumericalTypePromotion<complex<CppAD::AD<float> > , CppAD::AD<complex<long double> > >::ReturnType).name() <<endl; //no idea what this means but it is silly to go further
*/

  which_addition_is_used();
  //okay
  /*
  cout<<endl<<"Test operator+=(MyVec) -- only addition of 2 nums:"<<endl;
  MyVec<PrecType> Vo(pdim);
  Vo[0] = 0.12849324551440012;
  Vo[1] = 2.00239964361001287;

  MyVec<PrecType> Noml(2);
  Noml[0] = Vo[0]; Noml[1] = AK[0];

  cout << "normal = "<< conventional_sum(Noml) <<endl;
  cout << "KAHAN  = "<< kahan_algorithm(Noml) <<endl;
  cout << "BABUSK = "<< kahan_babuska_algorithm(Noml) <<endl;
  

  Vo += AK;
  cout << "Vo = "<< Vo <<endl;
*/

  /*CPUTime<PrecType> cpu;
  cpu.start();
  WK = XK + 2.*YK + 0.5*XK*YK + 0.0008123371284*ZK + (YK*1.0004327 + 0.45*XK)*ZK + ZK+ YK;
  cpu.stop();

  cout << setprecision(18) << "WK = " << WK <<endl;
  */
/*
  cout<<endl<<"Test operator* for std::vector:"<<endl;
  vector<double> v1stl(dim), v2stl(dim), v3stl(dim);
  for(size_t i = 0; i < dim; ++i){
    v1stl[i] = a[i];
    v2stl[i] = b[i];
    v3stl[i] = c[3];
  }
  */
  
  
  



 
  /* cout<<endl<<"THIS WILL CAUSE SEVERE ERRORS DUE TO THE DIFFERENT SIZES OF THE VEC's!!"<<endl;
  MyVec<double> CoolVectorGonnaBad(dim);
  CoolVectorGonnaBad = X + 0.5*Bad2;   //undersized
  // CoolVectorGonnaBad= 1.*X + Bad;      //oversized
  // CoolVectorGonnaBad+= -2.05*X+(X*Bad); 
  cout<<CoolVectorGonnaBad <<endl;
  */
  
  
 
  
}
