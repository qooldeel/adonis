#include <iostream>
#include "MAVEC.h"

using namespace std;

int main(){

  matrix<double> A(3,3);

  double a[3][3] = {{1,3,-2},
		    {3,5,6},
		    {2,4,3}};

  vec<double> b(3);
  b[0] = 5; b[1] = 7; b[2] = 8;

  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      A[i][j] = a[i][j];

  cout << "solution = "<< endl<<A.gausselimpartpivot(b) << endl;
 
  //cout << "Rank(A) = "<< A.rank() << endl;

}
