#include "../want_2_use.h"  //put it at the very beginning of the main file!

#include <iostream> 
#include <complex>

#include "../common/fancymessages.hh"
#include "../expressiontemplates/exprvec.hh"

#include "sparsematrix.hh"
#include "sparseutilities.hh"


#include "umfpack.h"
#include "umfdriver.hh"

#include "sparsitypattern.hh"

#include <cppad/cppad.hpp>

#include "../derivatives/derivativecalculator.hh"
#include "../ode/sourceterms.hh"
#include "sparsetmp.hh"
#include "../misc/xscaling.hh"

#include "squarelinearsystem.hh"



using namespace std;
using namespace Adonis;
using namespace ExprTmpl;

int main(){
  cout << "COMPLEX EXAMPLES"<<endl << "taken from \"http://help.scilab.org/docs/5.3.0/en_US/umfpack.html\"" << endl;
  int nz = 12;
  typedef complex<double> CType;
  CType a1(2,1),a2(3,-1), a3(3,2),a4(-1,1),a5(4,1),a6(-3,6),a7(1,-5),a8(2,-1),
    a9(6,-3);

  //rhs 
  CType b1(3,13), b2(58,32), b3(-19,13), b4(18,-12), b5(22,16);
  
  CType Acom[] = {a1,a2,a3,a4,4.,a5,a6,a7,a8,a8,a9,1.};

  CType bcom[] = {b1,b2,b3,b4,b5};
  MyVec<CType> bcplxrhs(bcom,bcom+5);
  //solution must be x [1+i; 2+2i; 3+3i; 4 + 4i; 5+5i]


  double are[12];
  MyVec<double> aim(nz);

  UmfpackUsesComplexNumbers<NumericDataTypeChecker<CType>::IsComplex>::split_complex_container_into_real_values(nz,Acom,are,aim.begin());

  cout << "Split into 2 arrays:"<<endl;
  cout << "real:"<<endl;
  print_all(are,are+nz);
  cout << "imag:"<<endl;
  print_all(aim);

  cout << endl<< "Test resize function:"<<endl;
  MyVec<double> tores(3);
  UmfpackUsesComplexNumbers<NumericDataTypeChecker<CType>::IsComplex>::resize(7,tores);
  cout << "tores = "<< tores <<endl;

  cout << endl<< "Test signature:" <<endl;
  typedef CompressedSparseStorage<double,'C',int> UmfStyleMatrixType;
  UmfStyleMatrixType CSC;

  cout << "compression style: "<< UmfStyleMatrixType::compressionType << endl;

  
  cout << endl<< "========================= UMFPACK DRIVER==========================="<<endl;
  int nrows = 5;
	//double x[5];

	//! Numeric void pointer contains the numerical factorization,
	//! that is used by Umfpack's solve routine
	void *Symbolic, *Numeric; 


	/* cumulative count of entries, as matrix is scanned columnwise */
	int Ap[] = { 0, 2, 5, 9, 10, 12 };

	/* row indices of entries, as matrix is scanned columnwise */
	int Ai[] =    { 0, 1, 0,  2, 4, 1,  2, 3, 4, 2, 1, 4 };

	/* matrix entries */
	double Ax[] = { 2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1 };

	/* the right hand side */
	double b[] = { 8, 45, -3, 3, 19 };

	MyVec<double> Av(Ax,Ax+12), sol(nrows);
	//	MyVec<double> bv(b,b+n);
	MyVec<int> indexv(Ai,Ai+12), pv(Ap,Ap+nrows+1); 

	//double* Info = 0;
	//double* az = 0;  //or just use 0 as argument in the fcts below 
	double Info[UMFPACK_INFO];
	double* Control = 0; 
	int sys = 0;
	
	MyVec<double> ImagAx(1);

#if USE_UMFPACK
	//!invoke umfpack routines to solve a real system A·x = b
	UmfpackDriver<double,int>::symbolic_analysis(nrows,nrows,&Ap[0],Ai,&Ax[0],&ImagAx[0],&Symbolic,Control,Info);
	UmfpackDriver<double,int>::lu_decomposition(Ap,Ai,Ax,0,Symbolic,&Numeric,Control,Info);

	UmfpackDriver<double,int>::solve(sys,Ap,Ai,Ax,0,&sol[0],0,b,0,Numeric,Control,Info);
	UmfpackDriver<double,int>::free_symbolic(&Symbolic);
	UmfpackDriver<double,int>::free_numeric(&Numeric);

	cout << endl<< "Solution = "<< sol << endl;

	// umfpack_di_symbolic(nrows,nrows,Ap,Ai,Ax,&Symbolic,0,Info);
	// umfpack_di_free_symbolic(&Symbolic);


	//! !invoke umfpack routines to solve a complex system A·x = b
	cout << endl << "  COMPLEX version"<<endl;
	int order = 5;
	MyVec<double> Re(nz), Im(nz), xx(order), xz(order), bx(order), bz(order);
	UmfpackUsesComplexNumbers<NumericDataTypeChecker<CType>::IsComplex>::split_complex_container_into_real_values(nz,Acom,Re.begin(),Im.begin());
	UmfpackUsesComplexNumbers<NumericDataTypeChecker<CType>::IsComplex>::split_complex_container_into_real_values(order,bcom,bx.begin(),bz.begin());	
	
	int ixcx[] = {0,1, 0,2,4, 1,2,3,4, 2, 1, 4};
	int pcx[] = {0,2,5,9,10,nz};
	MyVec<int> Icx(ixcx,ixcx+12), Pcx(pcx,pcx+order+1);
	UmfpackDriver<CType,int>::symbolic_analysis(order,order,&Pcx[0],&Icx[0],&Re[0],&Im[0],&Symbolic,Control,Info);
	UmfpackDriver<CType,int>::lu_decomposition(&Pcx[0],&Icx[0],&Re[0],&Im[0],Symbolic,&Numeric,Control,Info);
	UmfpackDriver<CType,int>::free_symbolic(&Symbolic);
	UmfpackDriver<CType,int>::solve(sys,&Pcx[0],&Icx[0],&Re[0],&Im[0],&xx[0],&xz[0],&bx[0],&bz[0],Numeric,Control,Info);
	UmfpackDriver<CType,int>::free_numeric(&Numeric);
	
	MyVec<CType> Csol(order);
	UmfpackUsesComplexNumbers<NumericDataTypeChecker<CType>::IsComplex>::join(order,&Csol[0],&xx[0],&xz[0]);

	//!MATLAB(R) solution: 
	//! 1.0000 + 1.0000i
	//! 2.0000 + 2.0000i
	//! 3.0000 + 3.0000i
	//! 4.0000 + 4.0000i
	//! 5.0000 + 5.0000i
	cout << "Cplx sol = "<< Csol << endl;
#endif


	////////////////////////////////////////////////////////////////
	FancyMessages FM;
	FM.display_in_bold("Hello Fucker this is in bold so that you don't miss reading it ;).");
	cout << "Now normal style again..."<<endl;
	

	const int rw = 5;
	const int cw = 3;
	double fmx[] = {1, 0, 3,
			0, 0, 2,
			0, 4, 0,
			1, 2, 4,
			6, 0, 2};
	
	MyVec<double> fullMtx(fmx,fmx+rw*cw);

	SparsityPattern<bool> SP(rw,cw);
	//SparsityPattern<set<size_t> > SP(rw,cw);
	SP.create(fullMtx);
	cout << "Status: #nonz = "<< SP.number_of_nonzeros() << "   is dense: "<< SP.is_dense() << "   is strictly dense: "<< SP.is_strictly_dense()<< endl << endl;
	cout << SP << endl;
	//print_in_matrix_style(SP.get_p(),cw);
	
	
	cout << endl<< "Check with example from \"http://www.coin-or.org/CppAD/Doc/sparse_jacobian.cpp.xml\"" <<endl;

	typedef CppAD::AD<double> ADType;
	const int domdim = 4, randim = 3;
	MyVec<ADType> X(domdim), Y(randim); //this is a randim x domdim Jacobian

	CppAD::Independent(X);
	
	Y[0] = X[0] + X[1];
	Y[1] = X[2] + X[3];
	Y[2] = X[0] + X[1] + X[2] + X[3]*X[3]/2.;

	//! a more difficult example
	// Y[0] = X[0]*X[1] + X[2];
	// Y[1] = X[2] + X[3];
	// Y[2] = X[0] + X[1]*(X[0]+2.75*X[1]) + X[2] + X[3]*X[3]/2.;


	CppAD::ADFun<double> adseq;
	adseq.Dependent(X,Y);
	MyVec<double> eval(domdim);
	for(int i = 0; i < domdim; ++i) 
	  eval[i] = double(i);
	cout << "eval = " << eval << endl << endl;
	MyVec<double> jacDense = adseq.SparseJacobian(eval);

	print_in_matrix_style(jacDense,domdim, "a dense Jacobian",1);
	
	
	//typedef bool SiType;
	typedef set<size_t> SiType; 
	typedef SparsityPattern<SiType> PatternType;
	

	//! note Jacobian defines a randim x domdim matrix
	PatternType Pattern(randim,domdim);
	Pattern.create(jacDense);
	cout << "Sparsity pattern of jacobian:" << endl << Pattern << endl;
	
	cout << "Compute Jacobian with known sparsity pattern:"<<endl;
	cout << "Patter.size(): "<< Pattern.size() << "  #(nonz): "<< Pattern.number_of_nonzeros() << endl;
	
	MyVec<double> jacSparse = ((Pattern.number_of_nonzeros() == 0) ? adseq.SparseJacobian(eval) : adseq.SparseJacobian(eval,Pattern.get_p()));

	print_in_matrix_style(jacSparse,domdim, "densly computed Jac",1);
	

	cout << setprecision(3) << endl<< "Test reciprocal of number:"<<endl;
	double d1 = -0.25;
	complex<double> z1(3.,-4.);


	cout << "double:  " << Reciprocal(d1) << endl;
	cout << "complex: " << Reciprocal(z1) << endl;


	int rows = 4,
	  cols = 5,
	  numzeros = 7;
	double val[] = {10,12,11,13,16,11,13};
	int colInd[] = {0,3,2,4,1,2,4};
	int rowPtr[] = {0,2,4,5,numzeros};

	CompressedSparseStorage<double,'R',int> CSRwiki(rows,cols,numzeros,colInd,colInd+numzeros,rowPtr,rowPtr+rows+1,val,val+numzeros);

	cout << "CSRwiki:" << CSRwiki << endl;

	CompressedSparseStorage<double,'R',int> DoesItWork;
	DoesItWork =  CompressedSparseStorage<double,'R',int>(rows,cols,numzeros,colInd,colInd+numzeros,rowPtr,rowPtr+rows+1,val,val+numzeros);
	cout << "DoesItWork = "<< DoesItWork << endl;
	

	cout << endl << "FILL from container"<<endl;
	
	MyVec<double> tf(30);
	tf <<= 0,1,2,0,3,4,
	  5,0,0,0,6,0,
	  0,0,7,8,0,0,
	  0,9,0,10,0,11,
	  12,13,0,0,14,0;
	
	CompressedSparseStorage<double,'C',int> CoSpSt;
	CoSpSt.fill_from_row_wise_dense(5,6,tf);
	cout << "CoSpSt = "<< CoSpSt <<endl; 
	

	cout << setprecision(5) << endl<< "Test sparse matrix operations:"<<endl;
	typedef MyVec<set<size_t> > VecSetType;
        int Dim = 4,
	  nnZ = 8;
	MyVec<double> diag(Dim);
	diag <<= 0.5, 0.75, 1.5, -1;
	VecSetType p_s(Dim);
	p_s[0].insert(0); p_s[0].insert(3);
	p_s[1].insert(2); p_s[1].insert(3);
	p_s[2].insert(0); p_s[2].insert(1); p_s[2].insert(3);
	p_s[3].insert(2);
	
	MyVec<double> Amx(nnZ), res(nnZ);
	Amx <<= 1.,3., 4.,5., 6.,7.,8., 9.;
	
	right_diagonal_sparse_matrix_multiplication(res,Amx,diag,p_s);
	cout << "D1*A = "<< res << endl;
	left_diagonal_sparse_matrix_multiplication(res,Amx,diag,p_s);
	cout << "A*D1 = "<< res << endl;
	
	cout << endl << "---------- Test 2----------" << endl;
	Dim = 4;
	nnZ = 8;
	MyVec<double> diag1(Dim), diag2(3);
	diag1 <<= 12,-20,30,-9;
	diag2 <<= 2.5, -10, -13;
	p_s.clear();
	p_s.resize(Dim);
	p_s[0].insert(0); p_s[0].insert(2);
	p_s[1].insert(0); p_s[1].insert(1); p_s[1].insert(2);
	p_s[2].insert(2);
	p_s[3].insert(0); p_s[3].insert(1);
	
	MyVec<double> Amx1(nnZ), res1(nnZ);
	Amx1 <<= 1,3, -2, 0.5,0.75, -4, 0.65,3.45;
	
	right_diagonal_sparse_matrix_multiplication(res1,Amx1,diag1,p_s);
	cout << "D1*A = "<< res1 << endl;
	left_diagonal_sparse_matrix_multiplication(res1,Amx1,diag2,p_s);
	cout << "A*D1 = "<< res1 << endl;
	
	cout << endl << "---------- Test 3----------" << endl;
	Dim = 3;
	nnZ = 6;
	p_s.clear();
	p_s.resize(Dim);
	p_s[0].insert(1); p_s[0].insert(2); p_s[0].insert(3);
	p_s[2].insert(0); p_s[2].insert(1); p_s[2].insert(2);
	MyVec<double> Amx2(nnZ), res2(nnZ);

	Amx2 <<= 2.3, 1.33, -8.9, 2.5, -3.8, -10;

	right_diagonal_sparse_matrix_multiplication(res2,Amx2,diag2,p_s);
	cout << "D2*A = "<< res2 << endl;
	left_diagonal_sparse_matrix_multiplication(res2,Amx2,diag1,p_s);
	cout << "A*D1 = "<< res2 << endl;
	

	
	cout << endl 
	     << "======================================================="<<endl              << "      Test new Jacobian generation object..."<<endl
	     << "======================================================="<<endl;

	//******* TODO change format (and if required format) *********
	const bool ISSPARSE = false;
	const char MTX_OD = 'r'; //no meaning for dense calculations
	bool computeGprime = true;
	//**************************************************************


	if(ISSPARSE==true){
	  cout << endl << endl;;
	  cout << "==================================="<<endl;
	  cout << "=       SPARSE Jacobian and G'     ="<<endl;
	  cout << "==================================="<<endl;
	  
	}
	else{
	  cout << endl << endl;;
	  cout << "==================================="<<endl;
	  cout << "=       DENSE Jacobian and G'     ="<<endl;
	  cout << "==================================="<<endl;
	}

	//COOL :) reset setprecision flags
	cout << resetiosflags( ios::fixed | ios::showpoint ) << endl;	

	typedef JacobianMatrixCalculator<ISSPARSE,double,MechWithSparseJac,MyVec,MTX_OD> SpMtxType; 
	MechWithSparseJac<double> fun1(3);
	SpMtxType JacMtx(computeGprime,fun1);
	//for ODE these two objects already exist
	MyVec<double> init;
	//MechWithSparseJac<double> mfun(3);
	//JacMtx.initialize(3,3,init,mfun);
	JacMtx.initialize(3,3,init);

	MyVec<double> eval1(3);
	eval1 <<= -0.5, 0.75, -1.5;
	MyVec<double> exctJac(5);
	exctJac <<= 2*eval1[0]*eval1[1], eval1[0]*eval1[0],
	  -2.5, 1, 2*eval1[2];
	cout << "J_exct <R> = " << exctJac;

	MyVec<double> Jres = JacMtx.evaluate_jacobian(eval1);
	
	double h = 0.0025, acoeff = 1, bcoeff = 1;
	JacMtx.compute_G_prime(h,acoeff,bcoeff);
	
	MyVec<double> exactGpr(6);  
	exactGpr <<= acoeff*1-h*bcoeff*2*eval1[0]*eval1[1], -h*bcoeff*eval1[0]*eval1[0], 
	  1, 2.5*h*bcoeff, 
	  -h*bcoeff, acoeff*1 - h*bcoeff*2*eval1[2];

	cout << "G'_exct  <R> = "<< exactGpr << endl; 
	MyVec<double> cGpr(6); //uses mapped order of <R>
	cGpr <<= exactGpr[0], exactGpr[1], exactGpr[2], exactGpr[4], exactGpr[3], exactGpr[5];

	cout << "G'_exct  <C> = " << cGpr << endl;
	cout << "G'_num   <"<<JacMtx.format_info() << "> =  ";
	JacMtx.reorder_values(); //reorder values if column compressed 
	print_all(JacMtx.value_pointer(),JacMtx.value_pointer() + JacMtx.dimension_of_values());

	

	cout << endl << "-------------- Test 2 -------------------"<<endl;
	typedef JacobianMatrixCalculator<ISSPARSE,double,MechWithSparseJac_1,MyVec,MTX_OD> JacType1;
	MechWithSparseJac_1<double> fun2(3);
	JacType1 JT2(computeGprime,fun2);
	
	MyVec<double> init2(4);
	//domaindim, rangedim
	JT2.initialize(4,3,init2);
	
	MyVec<double> eval2(4);
	eval2 <<= -1.65, 0.55, 0.24, 3.75;
	cout << "Jac_cppad_ex = "<< JT2.evaluate_jacobian(eval2) << endl;

	JT2.reorder_values();
	cout << "J <"<<JT2.format_info() <<"> = "; print_all(JT2.value_pointer(),JT2.value_pointer()+JT2.dimension_of_values());

	cout << endl << "-------------- Test 3 -------------------"<<endl;
	typedef JacobianMatrixCalculator<ISSPARSE,double,MechWithSparseJac_2,MyVec,MTX_OD> JacType3;
	MechWithSparseJac_2<double> fun3(3);
	JacType3 JT3(computeGprime,fun3);
	
	MyVec<double> init3(2);
	//domaindim, rangedim
	JT3.initialize(2,3,init3);
	
	MyVec<double> eval3(2);
	eval3 <<= -0.75, 1.8;
	cout << "Jac_cppad_ex_3 = "<< JT3.evaluate_jacobian(eval3) << endl;

	JT3.reorder_values();
	cout << "J <"<<JT3.format_info() <<"> = "; print_all(JT3.value_pointer(),JT3.value_pointer()+JT3.dimension_of_values());

	cout << endl << "-------------- Test 4 -------------------"<<endl;
	typedef JacobianMatrixCalculator<ISSPARSE,double,MechWithSparseJac_3,MyVec,MTX_OD> JacType4;
	MechWithSparseJac_3<double> fun4(5);
	JacType4 JT4(computeGprime,fun4);
	
	MyVec<double> init4(4);
	//domaindim, rangedim
	JT4.initialize(4,5,init4);
	
	MyVec<double> eval4(4);
	eval4 <<= 0.5, 2.15, 4.5, 1.5;

	cout << "Jac_cppad_ex_4 = "<< JT4.evaluate_jacobian(eval4) << endl;

	JT4.reorder_values();
	cout << "J <"<<JT4.format_info() <<"> = "; print_all(JT4.value_pointer(),JT4.value_pointer()+JT4.dimension_of_values());

	cout << endl << "-------------- Test 5 -------------------"<<endl;
	typedef JacobianMatrixCalculator<ISSPARSE,double,MechWithSparseJac_4,MyVec,MTX_OD> JacType5;
	MechWithSparseJac_4<double> fun5(4);
	JacType5 JT5(computeGprime,fun5);
	
	MyVec<double> init5(5);
	//domaindim, rangedim
	JT5.initialize(5,4,init5);
	
	MyVec<double> eval5(5);
	eval5 <<= 0.15, -3.0, 2.5, -0.75, 1.2;

	cout << "Jac_cppad_ex_5 = "<< JT5.evaluate_jacobian(eval5) << endl;

	JT5.reorder_values();
	cout << "J <"<<JT5.format_info() <<"> = "; print_all(JT5.value_pointer(),JT5.value_pointer()+JT5.dimension_of_values());


	cout << endl << "=================== G' TEST SECTION ================="<<endl;
	typedef JacobianMatrixCalculator<ISSPARSE,double,MechWithSparseJac_5,MyVec,MTX_OD> GPType1;
	MechWithSparseJac_5<double> fun6(5);
	GPType1 GP1(computeGprime,fun6);
	
	MyVec<double> init6(5);
	//domaindim, rangedim
	GP1.initialize(5,5,init6);
      
	GP1.evaluate_jacobian(eval5);         //1.) compute Jacobian
	GP1.compute_G_prime(h,acoeff,bcoeff); //2.) compute G'
	GP1.reorder_values();                 //3.) reorder values if necessary
	
        
	MyVec<double> exactGp1(8), J_exact_R(5);
	J_exact_R <<= eval5[2]*eval5[2], 2*eval5[1]*eval5[2],1, eval5[4], eval5[3];
	exactGp1 <<= acoeff*1, acoeff*1-h*bcoeff*eval5[2]*eval5[2], -h*bcoeff*2*eval5[1]*eval5[2], 1, -h*bcoeff, acoeff*1-h*bcoeff*eval5[4], -h*bcoeff*eval5[3],1;
	cout << "J<R>_exact  = "<< J_exact_R << endl;
	cout << "G'<R>_exact = " << exactGp1 << endl;
	MyVec<double> exactGp1_C(8);
	exactGp1_C <<= exactGp1[0], exactGp1[1], exactGp1[4], exactGp1[2], exactGp1[3], exactGp1[5], exactGp1[6], exactGp1[7];
	cout << "G'<C>_exact = " << exactGp1_C << endl;
	cout << "G'<"<<GP1.format_info()<<">_num   =    ";
	print_all(GP1.value_pointer(),GP1.value_pointer() + GP1.dimension_of_values());

	cout << endl;
	cout << "====================================================="<<endl;
	cout << "=           SCALE SQUARE LINEAR SYSTEM              ="<<endl;
	cout << "====================================================="<<endl;

	const bool SCALE =  true;
	const bool SPARSE = true;


	MyVec<double> valsOfMtx;

	if(SPARSE){
	  valsOfMtx.resize(11); //nnz
	  //! values in any case in CSR format because of pattern!!!
	  valsOfMtx <<= -3.75, 0.8, 23.4, 
				-0.05, -100.7, 
				0.00035, 
				-5000.8, 1.9, 
	    10212.4,-0.034, 42.41;
	}
	else{ //dense variant
	  valsOfMtx.resize(5*5);
	  valsOfMtx <<= 0,-3.75, 0, 0.8, 23.4,
	                -0.05,0,0,0,-100.7,
			 0,0,0.00035,0,0,
			 0,-5000.8,1.9,0,0,
			 0,10212.4,-0.034,0,42.14;
	}
	
	//just dummy or not
	SparsityPattern<set<size_t> > patt(5,5);
	patt[0].insert(1); patt[0].insert(3); patt[0].insert(4);
	patt[1].insert(0); patt[1].insert(4);
	patt[2].insert(2);
	patt[3].insert(1); patt[3].insert(2); 
	patt[4].insert(1);patt[4].insert(2); patt[4].insert(4);
	
	cout << "Sparsity Pattern = "<< endl<< patt << endl;
	

	MyVec<double> y_n(5), brisi(5);
	y_n <<= 22.35, 0.1, -1.75, -0.024, 301.07;
	brisi <<= 0.00012, +3004.1, 1.75, 0.0085, -72.345;
	double tol1 = 1.e-04,
	  tol2 = 1.e-03;
	cout << "y_n = "<< y_n << endl;
	ScaleSquareLinearSystem<SCALE,SPARSE,MyVec<double>,SparsityPattern<set<size_t> > > SLSsprs(patt,tol1,y_n.size());

	//for(int i = 0; i < 3; ++i){ // //ok iterator stuff ok
	SLSsprs.scale(valsOfMtx.begin(),brisi.begin(),y_n,tol1,tol2);

	cout << setprecision(7) << "values = "; 
	if(SPARSE)
	  print_all(valsOfMtx);
	else{
	  print_in_matrix_style(valsOfMtx,y_n.size());
	  // for(size_t i = 0; i < y_n.size(); ++i){
	  //   for(size_t j = 0; j < y_n.size(); ++j){
	  //     cout << valsOfMtx[i*y.size() + j] << " ";
	  //   }
	  //   cout << endl;
	  // }
	}

	
	  
	cout <<endl<< "b = "<< brisi << endl;
	// 	}

	cout << "==============================================="<<endl;
	cout << "=          TEST SQUARE LINEAR SOLVER          ="<<endl;
	cout << "==============================================="<<endl;

	//! NOTE: to test with SCILAB, use inv(A)*b, where b = [b1; b2;...;bn]
       
	int sysType = 0;
#if USE_UMFPACK
	sysType = UMFPACK_A;
#endif
	const bool UMF_PATTERN_ONCE = true;

	size_t idx1[] = {0,4, 1, 0,2, 3, 1,4};
	size_t pos1[] = {0,   2, 3,   5, 6,   8};
	
	bool make_singular = false;
	
	double dgent1 = 2.5;
	
	if(make_singular)  //provoke singular matrix by zero diag entry
	  dgent1 = 0.; 

	double val1[] = {0.5,-4.5, dgent1, 8.75,0.75, -3.85, -3.15,12.5};  

        typedef CompressedSparseStorage<double,'c',int> MatType;
	MatType S1(5,5,8, idx1,idx1+8, pos1, pos1+6,val1,val1+8);

	cout << "S1 = "<< S1 << endl;
	MyVec<double> br1(5);
	br1 <<= 2.35, -6.25, 0.25, 13.345, -0.45;


	SquareLinearSystem<MatType,MyVec<double>,CompressedSparseStorageIdentity<MatType::compressionType>::Value,UMF_PATTERN_ONCE,int> SQLS1(S1,br1,1,sysType);
	
	SQLS1.ls_solver_settings(0,2);

	SQLS1.solve();
	

	cout << "solution = " << br1 << endl;
	cout << "("<< SQLS1.system_to_solve() << ")"<<endl;

	cout << endl << "Complex square system"<<endl;
	complex<double> zv1(-3,4), zv2(0.5,-0.75), zv3(0.345,0), zv4(-2.5,-3),
	  zv5(-12,8), zv6(4,5), zv7(0.25,-0.25), zv8(23.5, -2), zv9(0.995,1.5),
	  bv1(0,13), bv2(-5,-6.5), bv3(0.45,-2.65), bv4(-0.55,30.5);

	int idx2[] = {0,2,3, 0,1,3, 2, 1,3};
	int pos2[] = {0,     3,     6, 7,   9};
	complex<double> val2[] = //{zv1,zv2,zv3,zv4,zv5,zv6,zv7,zv8,zv9}; 
{zv1,zv5,zv7,zv2,zv3,zv8,zv6,zv4,zv9};
	
	typedef CompressedSparseStorage<complex<double>,'c',int> CplxMatType;
	CplxMatType S2(4,4,9,idx2,idx2+9,pos2,pos2+5,val2,val2+9);
	MyVec<complex<double> > br2(4);
	br2 <<= bv1,bv2,bv3,bv4;
	cout << "br2 = "<< br2 << endl;
	cout << "S2 = "<< S2 << endl;
	
	SquareLinearSystem<CplxMatType,MyVec<complex<double> >,'C',UMF_PATTERN_ONCE,int> SQLS2(S2,br2,1,sysType);
	//SQLS2.system_to_solve(2);
	SQLS2.solve();
	
	cout << "cplx solution = "<< br2 << endl;
	cout << "("<< SQLS2.system_to_solve() << ")"<<endl;
								 

	cout << "test above systems with DENSE ls solver" <<endl;
	MyVec<double> S3(5*5), br3(5);
	S3 <<= 0.5, 0, 8.75, 0, 0,
	  0, 2.5, 0, 0, -3.15,
	  0,0,0.75,0,0,
	  0,0,0,-3.85,0,
	  -4.5,0,0,0,12.5;
	
	br3 <<= 2.35,  - 6.25,    0.25,    13.345,  - 0.45;

	typedef SquareLinearSystem<MyVec<double>,MyVec<double> > DenseMatType;

	DenseMatType DM1(S3,br3,1);
	DM1.solve();
	cout << "dense solution = "<< br3 << endl;
	cout << "("<< DM1.system_to_solve() << ")"<<endl;
	
	cout << endl<< "complex "<<endl;
	MyVec<complex<double> > S4(4*4), br4(4);
	S4 <<= zv1,zv2,0,0,
	  0,zv3,0,zv4,
	  zv5,0,zv6,0,
	  zv7,zv8,0,zv9;

	//S4 = transpose_array(S4,4);

	//cout << "transposed array = "<<endl;
	print_in_matrix_style(S4,4);
	
	br4 <<= bv1,bv2,bv3,bv4;
	
	SquareLinearSystem<MyVec<complex<double> >,MyVec<complex<double> >,'D'> DM2(S4,br4,1);
	DM2.solve();
	cout << "cplx dense solution = "<< br4 << endl;
	cout << "("<< DM2.system_to_solve() << ")"<<endl;

	cout << endl << "-----------------------------------------"<<endl;
	cout << endl << " SINGULAR sparse matrix" << endl;
	size_t idx5[] = {0,1, 0,1, 1,2};
	size_t pos5[] = {0,   2,   4,  6};

	double val5[] = {2,4, 1,2, -1,3};  

	MatType S5(3,3,6, idx5,idx5+6, pos5, pos5+4,val5,val5+6);

	cout << "S5 = "<< S5 << endl;
	MyVec<double> br5(3);
	br5 <<= -1,-7, 4;


	SquareLinearSystem<MatType,MyVec<double>,CompressedSparseStorageIdentity<MatType::compressionType>::Value,UMF_PATTERN_ONCE,int> SQLS5(S5,br5,1,sysType);
	

	SQLS5.solve();
	

	cout << "solution = "<< br5 << endl;
	cout << "("<< SQLS5.system_to_solve() << ")"<<endl;

	cout << endl << "FURTHER SINGULAR MATRIX" << endl;
	size_t idx6[] = {1,2,3,4,5,  0,1,3,4,  0,1,3,4, 5, 0,1,4,5, 0,1,4,5};
	size_t pos6[] = {0,          5,        9,       13,14,      18,     22};   
	double val6[] = {-1,6,-2,-3,4,   2,-16,4,-2,  4,-32,8,-4,  -2,  1,1,-6,2,  -3,-3, 18,-6};

	MatType S6(6,6,22, idx6,idx6+22, pos6, pos6+7,val6,val6+22);

	cout << "S6 = "<< S6 << endl;
	MyVec<double> br6(6);
	br6 <<= -0.5,2.75,-0.75, -3.45, 0.0035, 19.75;

	
	SquareLinearSystem<MatType,MyVec<double>,CompressedSparseStorageIdentity<MatType::compressionType>::Value,UMF_PATTERN_ONCE,int> SQLS6(S6,br6,1,sysType);
	
	SQLS6.solve();
	

	cout << "solution = "<< br6 << endl;
	cout << "("<< SQLS6.system_to_solve() << ")"<<endl;



	cout << endl << "GPRIME EXAMPLE i.e. G' = a*I - b*S'(x):"<<endl;
	cout <<         "---------------"<<endl;  
	double hstep = 0.25,
	  afac = 1.,
	  bfac = 1.;

	MyVec<double> y(5);
	y <<= 0.5, -1.75, 3., -2.65, -0.95;
	MyVec<double> Gtan(5*5);   //exact Jacobian S'(y)
	Gtan <<= 0, 1, 2*y[2], 0, 0,
	           0, 12*y[1]*y[1], -y[4],0,-y[2],
	           -y[3], 0,0,-y[0],0,
	  0,2,0,0,0,
	  0,2*y[1]*y[3], 0, y[1]*y[1], 0;

	Gtan *= (-bfac*hstep);        // A := b*h*S'(y)
	for(int i = 0; i < 5; ++i)
	  Gtan[i*5+i] += afac;        //update diagonal of A to yield G'
	cout << fixed << endl;
	print_in_matrix_style(Gtan,5,"G' (exact)");
	cout <<endl;
	FunWithSparseJac5<double> fun(5);

	typedef JacobianMatrixCalculator<true,double,FunWithSparseJac5,MyVec,'C',int> JabberwookT; 
	
	JabberwookT Griever(true,fun);
	 Griever.initialize(5,5,y);
	 Griever.evaluate_jacobian(y);
	 Griever.compute_G_prime(hstep,afac,bfac);
	 
	 Griever.reorder_values();
	cout << "G' (via AD) = "<<Griever.format_info() <<"> = "; print_all(Griever.value_pointer(),Griever.value_pointer()+Griever.dimension_of_values());
	  
	
	
	// cout << endl << "++++++++++++++++++++++++++ TEST ++++++++++++++++++++++++++" << endl;

	// ADONIS_WARNING(Warning,"This is a warning message....");

	/* //works
	cout << endl <<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	MyVec<double> MvE(5);
	MvE <<= 0.67, -3.4, -6.8, 0.56, -2.3;
	
	double* mveit = &MvE[0];
	double* mvestart = mveit;
	for(int i = 0; i < 5; ++i){
	  *mveit *= 2;
	  mveit++;
	}
	mveit = mvestart; //reset to beginning!!
	for(int i = 0; i < 5; ++i){
	  *mveit += 1.5;
	  mveit++;
	}
	cout << "MvE = "<< MvE << endl;
	*/
	  
	/*
	  //! working so far
	cout << endl 
	     << "======================================================="<<endl              << "      Test Sparse Matrix selection...                  "<<endl
	     << "======================================================="<<endl;
	typedef VecSetType::value_type SetType;
	typedef SetType::const_iterator SetIteratorType;
	 size_t m = 4,
          n = 5,
          nnz = 7;
	 
	 const char Char = 'c';

	 typedef SparseMatrixSelector<Char,double> MtxSelectorType;

	  VecSetType p(4);
	  p[0].insert(0); p[0].insert(3);
	  p[1].insert(2); p[1].insert(4);
	  p[2].insert(1);
	  p[3].insert(2); p[3].insert(4);
   
	  SetIteratorType sit;
	  cout << "p = " << endl;
	  for(size_t i = 0; i < m; ++i){
	    for(sit = p[i].begin(); sit != p[i].end(); ++sit){
	      cout << *sit << " ";
	    }
	    cout << endl;
	  }
   
	  
	  MyVec<size_t> mapperC(nnz),rowC(MtxSelectorType::dim_row_array(m,n,nnz)),colptrC(MtxSelectorType::dim_col_array(m,n,nnz));
	  MtxSelectorType::create_format_data(m,n,rowC.begin(),colptrC.begin(),mapperC.begin(),p);

	   cout << endl;
	   print_all(rowC,9, true,"row");
	   print_all(colptrC,9,true,"col_ptr");
	   print_all(mapperC, 9, true,"mapper vec");

	   MyVec<double> valR1(nnz), valC1(nnz);
	   valR1 <<= -2.5, 0.75, 1.5, -3.15, -9.75, -3.05, 12.5;
	   MtxSelectorType::fill_in_values(nnz,valC1.begin(),valR1.begin(),mapperC.begin());
	   print_all(valC1, 9, true,"VAL");


	   cout << endl << "TEST 7 -- 2 empty col" <<endl;
	   m = 4;
	   n = 4;
	   nnz = 6;
	   VecSetType p7(m);
	   p7[0].insert(0); p7[0].insert(3);
	   p7[1].insert(0); p7[1].insert(3);
	   p7[2].insert(3);
	   p7[3].insert(0);
	   
	   cout << "p7 = " << endl;
	   for(size_t i = 0; i < m; ++i){
	     for(sit = p7[i].begin(); sit != p7[i].end(); ++sit){
	       cout << *sit << " ";
	     }
	     cout << endl;
	   }
   
	   //clear entries
	   mapperC.clear();
	   rowC.clear();
	   colptrC.clear();
	   mapperC.resize(nnz);
	   rowC.resize(MtxSelectorType::dim_row_array(m,n,nnz));
	   colptrC.resize(MtxSelectorType::dim_col_array(m,n,nnz));
    
   
    
	   MtxSelectorType::create_format_data(m,n,rowC.begin(),colptrC.begin(),mapperC.begin(),p7);
	   cout << endl;
	   print_all(rowC,9, true,"row");
	   print_all(colptrC,9,true,"col_ptr");
	   print_all(mapperC, 9, true,"mapper vec");
	   
	   
	   cout << endl << "TEST 2 -- empty col" <<endl;
   m = 3;
   n = 3;
   nnz = 4;
   VecSetType p1(m);
   p1[0].insert(0);
   p1[1].insert(0); p1[1].insert(1);
   p1[2].insert(1);
   cout << "p1 = " << endl;
   for(size_t i = 0; i < m; ++i){
      for(sit = p1[i].begin(); sit != p1[i].end(); ++sit){
        cout << *sit << " ";
     }
     cout << endl;
    }
   
    //clear entries
    mapperC.clear();
    rowC.clear();
    colptrC.clear();
    mapperC.resize(nnz);
    rowC.resize(MtxSelectorType::dim_row_array(m,n,nnz));
    colptrC.resize(MtxSelectorType::dim_col_array(m,n,nnz));
    
   
    
    MtxSelectorType::create_format_data(m,n,rowC.begin(),colptrC.begin(),mapperC.begin(),p1);
   cout << endl;
   print_all(rowC,9, true,"row");
   print_all(colptrC,9,true,"col_ptr");
   print_all(mapperC, 9, true,"mapper vec");

   cout << endl << "TEST 3 -- empty row" <<endl;
   m = 4; 
   n = 4; 
   nnz = 6;
   
   VecSetType p2(m);
   p2[0].insert(1); p2[0].insert(3);
   p2[2].insert(0); p2[2].insert(1); p2[2].insert(2);
   p2[3].insert(1);
   cout << "p2 = " << endl;
   for(size_t i = 0; i < m; ++i){
      for(sit = p2[i].begin(); sit != p2[i].end(); ++sit){
        cout << *sit << " ";
     }
     cout << endl;
    }
   
    //clear entries
    mapperC.clear();
    rowC.clear();
    colptrC.clear();
    mapperC.resize(nnz);
    rowC.resize(MtxSelectorType::dim_row_array(m,n,nnz));
    colptrC.resize(MtxSelectorType::dim_col_array(m,n,nnz));
    MtxSelectorType::create_format_data(m,n,rowC.begin(),colptrC.begin(),mapperC.begin(),p2);
   cout << endl;
   print_all(rowC,9, true,"row");
   print_all(colptrC,9,true,"col_ptr");
   print_all(mapperC, 9, true,"mapper vec");
   
   MyVec<double> valR3(nnz), valC3(nnz);
   valR3 <<= 0.5, 1.95, -2.4, 1.9, 4.15, -2.8;
   MtxSelectorType::fill_in_values(nnz,valC3.begin(),valR3.begin(),mapperC.begin());
   print_all(valC3, 9, true,"VAL");
   
   cout << endl << "TEST 4 " <<endl;
   m = 4;
   n = 2;
   nnz = 5;
   
   VecSetType p3(m);
   p3[0].insert(0);
   p3[1].insert(1);
   p3[2].insert(0); p3[2].insert(1);
   p3[3].insert(1);
   
   cout << "p3 = " << endl;
   for(size_t i = 0; i < m; ++i){
      for(sit = p3[i].begin(); sit != p3[i].end(); ++sit){
        cout << *sit << " ";
     }
     cout << endl;
    }
   
   //clear entries
    mapperC.clear();
    rowC.clear();
    colptrC.clear();
    mapperC.resize(nnz);
    rowC.resize(MtxSelectorType::dim_row_array(m,n,nnz));
    colptrC.resize(MtxSelectorType::dim_col_array(m,n,nnz));
     MtxSelectorType::create_format_data(m,n,rowC.begin(),colptrC.begin(),mapperC.begin(),p3);
   cout << endl;
   print_all(rowC,9, true,"row");
   print_all(colptrC,9,true,"col_ptr");
   print_all(mapperC, 9, true,"mapper vec");
   
   
   
   cout << endl << "TEST 5 " <<endl;
   m = 6;
   n = 6;
   nnz = 11;
   
   VecSetType p4(m);
   p4[0].insert(1);  p4[0].insert(3);  p4[0].insert(4); 
   p4[1].insert(1);  p4[1].insert(5);
   p4[2].insert(3); 
   p4[3].insert(2);  p4[3].insert(4);
   p4[4].insert(3); 
   p4[5].insert(0);  p4[5].insert(2);
   
   cout << "p4 = " << endl;
   for(size_t i = 0; i < m; ++i){
      for(sit = p4[i].begin(); sit != p4[i].end(); ++sit){
        cout << *sit << " ";
     }
     cout << endl;
    }
   
   //clear entries
    mapperC.clear();
    rowC.clear();
    colptrC.clear();
    mapperC.resize(nnz);
    rowC.resize(MtxSelectorType::dim_row_array(m,n,nnz));
    colptrC.resize(MtxSelectorType::dim_col_array(m,n,nnz));
    // if((Char == 'R') || (Char == 'r')){
    //   rowC.resize(m+1);
    //   colptrC.resize(nnz);
    // }
    MtxSelectorType::create_format_data(m,n,rowC.begin(),colptrC.begin(),mapperC.begin(),p4);
   cout << endl;
   print_all(rowC,9, true,"row");
   print_all(colptrC,9,true,"col_ptr");
   print_all(mapperC, 9, true,"mapper vec");

    cout << endl << "TEST 6 " <<endl;
   m = 5;
   n = 4;
   nnz = 7;
   
   VecSetType p5(m);
  
   p5[0].insert(0);  p5[0].insert(2);  p5[0].insert(3); 
   
   p5[2].insert(0); p5[2].insert(3); 
   p5[3].insert(3);  
   p5[4].insert(2); 
 
 
   cout << "p5 = " << endl;
   for(size_t i = 0; i < m; ++i){
      for(sit = p5[i].begin(); sit != p5[i].end(); ++sit){
        cout << *sit << " ";
     }
     cout << endl;
    }
 
   //clear entries
    mapperC.clear();
    rowC.clear();
    colptrC.clear();
 
   
    mapperC.resize(nnz);
    rowC.resize(MtxSelectorType::dim_row_array(m,n,nnz));
    colptrC.resize(MtxSelectorType::dim_col_array(m,n,nnz));
    // if((Char == 'R') || (Char == 'r')){
    //   rowC.resize(m+1);
    //   colptrC.resize(nnz);
    // }
    MtxSelectorType::create_format_data(m,n,rowC.begin(),colptrC.begin(),mapperC.begin(),p5);
   cout << endl;
   print_all(rowC,9, true,"row");
   print_all(colptrC,9,true,"col_ptr");
   print_all(mapperC, 9, true,"mapper vec");

 cout << endl << "TEST 8 " <<endl;
   m = 3;
   n = 4;
   nnz = 3;
   
   VecSetType p8(m);
  
   p8[1].insert(1);  p8[1].insert(3);  
   p8[2].insert(2);  
  
 
   cout << "p8 = " << endl;
   for(size_t i = 0; i < m; ++i){
      for(sit = p8[i].begin(); sit != p8[i].end(); ++sit){
        cout << *sit << " ";
     }
     cout << endl;
    }
 
   //clear entries
    mapperC.clear();
    rowC.clear();
    colptrC.clear();
 
   
    mapperC.resize(nnz);
    rowC.resize(MtxSelectorType::dim_row_array(m,n,nnz));
    colptrC.resize(MtxSelectorType::dim_col_array(m,n,nnz));
    // if((Char == 'R') || (Char == 'r')){
    //   rowC.resize(m+1);
    //   colptrC.resize(nnz);
    // }
    MtxSelectorType::create_format_data(m,n,rowC.begin(),colptrC.begin(),mapperC.begin(),p8);
   cout << endl;
   print_all(rowC,9, true,"row");
   print_all(colptrC,9,true,"col_ptr");
   print_all(mapperC, 9, true,"mapper vec");

*/


   /*
   cout << endl<< "=============== TEST EXTENSIONS FOR COMPRESSED MATRIX ==========="<<endl;
   CompressedSparseStorage<double, 'C', size_t> SPM;
   m = 3;
   n = 4;
   nnz = 5;
// 'C'-case 
    size_t rowC0[] = {0,1,1,2,3};
    size_t colC0[] = {0,2,3,4,  5};
   vector<size_t> row(rowC0,rowC0+5),
    col(colC0,colC0+5);
   
   SPM.resize(m,n,nnz);
   SPM.fill_index(rowC0);

   print_all(SPM.index(),SPM.index()+nnz);
   cout << "----------------- TEST 2 -----------------"<<endl; 	
   m = 4;
   n = 2;
   nnz = 4;
   size_t rowC1[] = {0,1, 1, 0};
   //size_t colC1[] = {0, 2, 3,   4};
   
   SPM.resize(m,n,nnz);
   SPM.fill_index(rowC1);
   print_all(SPM.index(),SPM.index()+nnz);
   
   cout << "----------------- TEST 3 -----------------"<<endl; 
   m = 4;
   n = 4;
   nnz = 5;
   size_t rowC2[] = {1,2, 0,2, 3};
   //size_t colC2[] = {0, 2, 4,   nnz};
   
   SPM.resize(m,n,nnz);
   SPM.fill_index(rowC2);
   print_all(SPM.index(),SPM.index()+nnz);
*/



	//complex<double> z2; cout << "complex: " << Reciprocal(z2) << endl;
	////! works
	// cout << endl << "Check sparsity pattern with that obtained from Cppad:"<< endl << "-------------------------------------------------------------" <<endl;
	// int m = randim,
	//   n = domdim;

	// //using packed boolean sparsity patterns
	// std::vector<bool> s_b(m * m), p_b(m * n);
	// for(int i = 0; i < m; i++)
	// {	for(int ell = 0; ell < m; ell++)
	// 		s_b[i * m + ell] = false;
	// 	s_b[i * m + i] = true;
	// }
	// p_b   = adseq.RevSparseJac(m, s_b);

	// print_in_matrix_style(p_b,n);
	// cout <<endl;

	// std::vector< std::set<size_t> >  s_s(m), p_s(m);

	// for(int i = 0; i < m; i++)
	//   s_s[i].insert(i);
	// p_s   = adseq.RevSparseJac(m, s_s);

	// typedef std::set<size_t>::const_iterator setiter_type;
	// for(int i = 0; i < m; ++i){
	//   for(setiter_type cit = p_s[i].begin(); cit != p_s[i].end(); ++cit){
	//     cout << *cit << " ";
	//   }
	//   cout << endl;
	// }


	// double mv[] = {0, -2, 0, 1, 0,-5,
	// 	       0,0,1,4,-9,0,
	// 	       -3,0,2,-3,0,-0.5};

	// MyVec<double> DsMtx(mv,mv+3*6);
	// print_in_matrix_style(DsMtx,6, "Fullmatrix",2);
	// Pattern.resize(3,6);

	// Pattern.create(DsMtx);

	// cout << "#nonz: "<< Pattern.number_of_nonzeros() << " m: "<< Pattern.rows() << "  n: "<< Pattern.cols()<<  endl;
	// cout << "New pattern:" << endl<< Pattern<<endl;

	// PatternType Pdef;

	// Pdef.resize(3,6);
	// Pdef.create(DsMtx);
	// cout << "created from scratch:"<<endl << Pdef << endl;

	////! copy cstrc and copy assign. work
	// PatternType P2(Pattern);
	// cout << "#nonz: "<< P2.number_of_nonzeros() << " m: "<< P2.rows() << "  n: "<< P2.cols()<<  endl << P2 << endl;
	
	
	// PatternType P3;
	// P3 = Pattern;
	// cout << "#nonz: "<< P3.number_of_nonzeros() << " m: "<< P3.rows() << "  n: "<< P3.cols()<<  endl << P3 << endl;

  ////! works
  // int rows = 4,
  //   cols = 5,
  //   nz = 7;
  // double val[] = {10,12,11,13,16,11,13};
  // int colInd[] = {0,3,2,4,1,2,4};
  // int rowPtr[] = {0,2,4,5,nz};

  // CompressedSparseStorage<double,'R',int> CSRwiki(rows,cols,nz,colInd,colInd+nz,rowPtr,rowPtr+rows+1,val,val+nz);

  // cout << "CSRwiki:" << CSRwiki << endl;

  // MyVec<double> a1(nz);
  // MyVec<int> i1(nz), p1(cols+1);

  // cout << "Test converter function:"<<endl;
  // row2column(rows,cols,nz,colInd,rowPtr,val, i1.begin(),p1.begin(),a1.begin());
  // cout << endl << "Test 2nd example:"<<endl;
  

  // int r1 = 4,
  //   c1 = 3,
  //   nz1 = 6;
  // double v1[] = {4,  1,  3,  5,  2,  1};
  // int cI1[] = {0,  1,  1,  2,  1,  2};
  // int cp1[] = {0,  2,  3,  4,  nz1}; 

  // MyVec<double> a2(nz1);
  // MyVec<int> i2(nz1), p2(c1+1);

  // row2column(r1,c1,nz1,cI1,cp1,v1, i2.begin(),p2.begin(),a2.begin());

  // cout << endl << "Example from Davis:" <<endl;
 
  // int ro = 4,
  //   co = 4,
  //   nonz = 10;
  // double ava[] = {4.5, 3.2,   3.1, 2.9, 0.9,   1.7, 3.,   3.5, 0.4, 1.};
 
  // int cix[] = {0, 2,    0,1,3,  1,2,  0,1,3};
  // MyVec<int> Cx(cix,cix+nonz);
  // int pix[] = {0,2,5,7, nonz};

  // MyVec<double> a3(nonz);
  // MyVec<double> ix3(nonz), p3(co+1);

  // row2column(ro,co,nonz,&Cx[0],pix,ava, &ix3[0],p3.begin(),&a3[0]);
  
  // cout <<endl<<endl<< "COLUMN 2 ROW:"<<endl << "---------"<<endl;
  // double awi[] = {10,  16,  11,  11,  12,  13,  13};
  // int iwi[] = {0,  2,  1,  3,  0,  1,  3};
  // int pwi[] = {0,  1,  2,  4,  5,  nz};
  
  // MyVec<double> Aw(nz);
  // MyVec<int> Iw(nz), Pw(rows+1);

  // column2row(rows,cols,nz,iwi,pwi,awi,&Iw[0],Pw.begin(),Aw.begin());
  // cout << endl;

  // double mya[] = {4,  1,  3,  2,  5,  1};
  // int myi[] = {0,  0,  1,  3,  2,  3};
  // int myp[] = {0,  1,  4,  6};

  // MyVec<double> Mya(nz1);
  // MyVec<double> Myi(nz1), Myp(r1+1);

  // column2row(r1,c1,nz1,myi,myp,mya,&Myi[0],&Myp[0],&Mya[0]);

  // cout << endl;
  
  // double xv[] = {4.5,  3.1,  3.5,  2.9,  1.7,  0.4,  3.2,  3.,  0.9,  1.};
  // int ix[] = {0,  1,  3,  1,  2,  3,  0,  2,  1,  3}; 
  // int pi[] = {0,  3,  6,  8,  nonz};

  // MyVec<double> ave(nonz);
  // MyVec<int> ive(nonz), pve(ro+1);

  // column2row(ro,co,nonz,ix,pi,xv,ive.begin(),pve.begin(),ave.begin());




  return 0;
}
