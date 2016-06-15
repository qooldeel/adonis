#include "bandmatrix.hh"
#include "../expressiontemplates/exprvec.hh"

using namespace std;
using namespace Adonis;
using namespace ExprTmpl;
main(){
  
  const int dim = 5;
  MyVec<double> mtx(dim*dim);

  mtx <<= 1, 6, 10, 0,0,
    13,2,-1,11,0,
    0,14,3,8,12,
    0,0,-2,4,9,
    0,0,0,16,5;

  print_in_matrix_style(mtx,dim);

  BandMtx<double> Bdef;
  cout << "Bdef (default) = " << Bdef << " left:" << Bdef.left_band_width() << "   right:" << Bdef.right_band_width() << endl;

  BandMtx<double> B1(dim,1,2,mtx);

  cout << "B1 = "<< B1 <<endl;

  cout << "copy:"<<endl;
  BandMtx<double> B2(B1);
  cout << "B2 = "<< B2 << endl;
  
  cout << "assignment:"<<endl;
  Bdef = B1;
  cout << "Bdef = "<< Bdef << endl;
  
  MyVec<double> v(dim);

  v <<= 0.5,-2.15,4,-0.5,1.5;

  cout << "B1*v = "<< B1*v << endl;

  MyVec<double> rhs(5);

  rhs <<= -0.45,3.64,-1.3,1.2,-3;

  B1.solve(rhs);

  cout << "solution = "<< rhs << endl; 
}
