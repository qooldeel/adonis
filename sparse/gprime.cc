#include <iostream>#include <vector>#include <set>#include <cassert>#include <string>/*** \brief Compute G' = aI - bhS'* This is only for square matrices. Actually, the sizes of the sparsity* patterns for both the Jacobian S' and G' are equal; the only difference lies in the fact* that G' posseses non-vanishing diagonal entries, thus the sets of each row differ by* at most one entry.* */template<class PATT, class VEC, class T>inline void create_G_prime(VEC& Gprime, const VEC& Jac, const PATT& pattNewt, const PATT& pattJac, const T& h, const T& a = 1, const T& b = 1){    assert(pattNewt.size() == pattJac.size());        typedef std::size_t SizeType;    typedef typename PATT::value_type SetType;    typedef typename SetType::const_iterator ConstSetIterType;        bool sameSize = false;    T bh = b*h;        SizeType k = 0; //positioin in G'    SizeType diffsz = 0;        ConstSetIterType pattJacIt,      setIt;            SizeType idx;    for(SizeType i = 0; i < pattNewt.size(); ++i){      pattJacIt = pattJac[i].begin();      (pattNewt[i].size() != pattJac[i].size()) ? sameSize = false : sameSize = true;                           for(setIt = pattNewt[i].begin(); setIt != pattNewt[i].end(); ++setIt){	if(i == *setIt){	  diffsz += pattNewt[i].size() - pattJac[i].size(); //either zero or one	  if(!sameSize){	    Gprime[k] = a*1;	  }          	  else{	    Gprime[k] = a*1 - bh*Jac[k-diffsz];	  }	}	else{	  Gprime[k] = -bh*Jac[k-diffsz];    	}          	k++;      }        }          }template<class SVEC>inline void print_pattern(const SVEC& p, const std::string& s = std::string()){     typedef typename SVEC::value_type SetType;     typedef typename SetType::const_iterator SetItType;         if(s.size() != 0)      std::cout << "Sparsity pattern \""<< s << "\":"<< std::endl;     for(std::size_t i = 0; i < p.size(); ++i){        std::cout << i << ".)  ";        for(SetItType it = p[i].begin(); it != p[i].end(); ++it){         std::cout << *it << " ";        }             std::cout << std::endl;    }    }template<class OBJ>inline void print_me(const OBJ& obj, const std::string& str = std::string()){    if(str.size() != 0)      std::cout << str << " = ";    std::cout << "[ "; for(typename OBJ::const_iterator it = obj.begin(); it != obj.end(); ++it)   std::cout << *it << " "; std::cout<<"]"<< std::endl;}template<class SVEC>inline void form_diag(SVEC& p){    for(std::size_t i = 0; i < p.size(); ++i){        p[i].insert(i);  //no double insertion possible            }    }using namespace std;int main(){   cout << "Test Me:" << endl;    int m = 4,   nnzJac = 7,   nnz = 8; //in G'      typedef vector<double> ValVecType;   typedef vector<set<size_t> > VecSetType;      VecSetType pattJac1(4), pattNewt;   pattJac1[0].insert(0); pattJac1[0].insert(1);   pattJac1[2].insert(1); pattJac1[2].insert(2); pattJac1[2].insert(3);   pattJac1[3].insert(1); pattJac1[3].insert(3);      cout << endl;   print_pattern(pattJac1,"pattJac1");   pattNewt = pattJac1;   form_diag(pattNewt);   print_pattern(pattNewt,"pattNewt");   double a[] = {4.5,-1.5, 0.5, 2.95, -0.75, 2.65, -0.95};   ValVecType jac(a,a+nnzJac),   Gprime(nnz);   double h = 1./4;   create_G_prime(Gprime, jac, pattNewt,pattJac1,h);      double tst[] = {1-h*4.5, -h*(-1.5), 1, -h*0.5, 1-h*2.95, -h*(-0.75), -h*2.65,1-h*(-0.95)};   ValVecType test1(tst,tst+nnz);   print_me(Gprime, "G'");   print_me(test1,  "a ");      cout << endl << "------------------ TEST 2 --------------------"<<endl;   m = 4;   nnzJac = 6;   nnz = 8; //in G'            VecSetType pattJac2(4), pattNewt2;   pattJac2[0].insert(1);   pattJac2[2].insert(1); pattJac2[2].insert(2); pattJac2[2].insert(3);   pattJac2[3].insert(1); pattJac2[3].insert(3);      cout << endl;   print_pattern(pattJac2,"pattJac2");   pattNewt2 = pattJac2;   form_diag(pattNewt2);   print_pattern(pattNewt2,"pattNewt2");   double a2[] = {-1.5, 0.5, 2.95, -0.75, 2.65, -0.95};   ValVecType jac2(a2,a2+nnzJac),   Gprime2(nnz);      create_G_prime(Gprime2, jac2, pattNewt2,pattJac2,h);      double tst2[] = {1,-h*(-1.5), 1, -h*0.5, 1-h*2.95, -h*(-0.75), -h*2.65,1-h*(-0.95)};    ValVecType test2(tst2,tst2+8);   print_me(Gprime2, "G'");   print_me(test2,  "a ");           cout << endl << "------------------ TEST 3 --------------------"<<endl;   m = 6;   nnzJac = 13;   nnz = 17; //in G'            VecSetType pattJac3(m), pattNewt3;   pattJac3[0].insert(2); pattJac3[0].insert(4);     pattJac3[1].insert(1); pattJac3[1].insert(3);   pattJac3[2].insert(0); pattJac3[2].insert(3); pattJac3[2].insert(4);pattJac3[2].insert(5);   pattJac3[3].insert(1); pattJac3[3].insert(2);   pattJac3[4].insert(4);   pattJac3[5].insert(0); pattJac3[5].insert(2);   cout << endl;   print_pattern(pattJac3,"pattJac3");   pattNewt3 = pattJac3;   form_diag(pattNewt3);   print_pattern(pattNewt3,"pattNewt3");   double a3[] = {0.45, 3.5, -5.6, 8.5, 0.75, 0.005, -0.75, 6.23, 9.01, -3.4, -0.7, 1.4, 0.3};   ValVecType jac3(a3,a3+nnzJac),   Gprime3(nnz);      create_G_prime(Gprime3, jac3, pattNewt3,pattJac3,h);       double tst3[] = {1, -h*a3[0], -h*a3[1],    1-h*a3[2], -h*a3[3],     -h*a3[4], 1, -h*a3[5], -h*a3[6], -h*a3[7],    -h*a3[8],-h*a3[9],1,    1-h*a3[10],    -h*a3[11], -h*a3[12], 1};     ValVecType test3(tst3,tst3+nnz);    print_me(Gprime3, "G'");    print_me(test3,  "a ");   cout << endl << "------------------ TEST 4 --------------------"<<endl;   m = 4;   nnzJac = 4;   nnz = 6; //in G'            VecSetType pattJac4(m), pattNewt4;   pattJac4[0].insert(0); pattJac4[0].insert(1);     pattJac4[3].insert(0); pattJac4[3].insert(3);   cout << endl;   print_pattern(pattJac4,"pattJac4");   pattNewt4 = pattJac4;   form_diag(pattNewt4);   print_pattern(pattNewt4,"pattNewt4");   double a4[] = {-2.5, 14.23, -0.44, 3.05};   ValVecType jac4(a4,a4+nnzJac),   Gprime4(nnz);      create_G_prime(Gprime4, jac4, pattNewt4,pattJac4,h);      double tst4[] = {1-h*a4[0],-h*a4[1],		    1,		    1,		    -h*a4[2], 1-h*a4[3]};     ValVecType test4(tst4,tst4+nnz);    print_me(Gprime4, "G'");    print_me(test4,  "a ");         cout << endl << "------------------ TEST 5 --------------------"<<endl;   m = 5;   nnzJac = 6;   nnz = 9; //in G'            VecSetType pattJac5(m), pattNewt5;   pattJac5[1].insert(1); pattJac5[1].insert(4);     pattJac5[3].insert(0); pattJac5[3].insert(1); pattJac5[3].insert(3); pattJac5[3].insert(4);   cout << endl;   print_pattern(pattJac5,"pattJac5");   pattNewt5 = pattJac5;   form_diag(pattNewt5);   print_pattern(pattNewt5,"pattNewt5");   double a5[] = {0.07, 1.32, -0.91, 0.31, 11.9, -5.32};   ValVecType jac5(a5,a5+nnzJac),   Gprime5(nnz);      create_G_prime(Gprime5, jac5, pattNewt5,pattJac5,h);      double tst5[] = {1,		    1-h*a5[0], -h*a5[1],		    1,		    -h*a5[2], -h*a5[3], 1-h*a5[4],-h*a5[5],		    1};     ValVecType test5(tst5,tst5+nnz);    print_me(Gprime5, "G'");    print_me(test5,  "a ");   return 0;}