#ifndef SCALING_AND_REGULARISATION_HH
#define SCALING_AND_REGULARISATION_HH

#include <cmath>
#include <typeinfo>
#include "../common/adonisassert.hh"
#include "../common/machinebase.hh"
#include "../common/globalfunctions.hh"
#include "normalspace.hh"
#include "dunextensions.hh"
#include "../dunestuff/fmatrix.hh" 
#include "../containers/fdiagonal.hh"


#include "../common/elementaryoperations.hh"

namespace Adonis{


 
  /**
   * \brief preconditioning by row equilibration. Multiply row i by \f$ r_i = \frac{1}{\max|a_{ij}|}.\f$ In practice \f$r_i\f$ is taken to be a number of the form  \f$ 2^m \f$ to avoid additional roundoff errors.
   *
   * CAUTION: Scaling has always been surrounded by a mystery and it  even might occur that a badly scaled problem becomes even worse after scaling!
   *
   * See e.g. [KINCAID/CHENEY Ch. 4] or [GOLUB/VANLOAN, § 3.5.2, pp. 125]
   */
  template<class T, int M, int N>
  inline void row_equilibration(Dune::FieldMatrix<T,M,N>& A){
    T r_i = T(); 
    T twolgthm = T();
    int m; 
    for(int i = 0; i < M; ++i){
      r_i = T();                     //reset
      for(int j = 0; j < N; ++j){      
	r_i = std::max(r_i,Abs(A[i][j]));  //determine max_{0<=j<N}(|a_ij|)
      }
      //adonis_assert(r_i != T());  //since log_2 0 isn't defined!
      if(r_i == T())
	twolgthm = T();
      else 
	twolgthm = log2(r_i);

      //std::cout << "max_"<<"|a|_ij = " << r_i << std::endl;

      //log_2 1/max(|a_ij|) = log_2 1 - log_2 max(|a_ij|) = - log_2 max(|a_ij|)
      //Hence m ~ -log_2 max(|a_ij|)
      m = round(-log2(r_i));
      r_i = pow_of_machine_base<double,2>(m);
	
      //std::cout << "r_i = 1./max|a_ij| ~ 2^m = "<<r_i << std::endl;
      
      for(int j = 0; j < N; ++j){      // multiply row i by r 
	A[i][j] *= r_i;
      }
    }
  }


  /**
   *\brief Don't alter input matrix. Rather store r_i's in a FieldDiagonal for later use.
   */
  template<class T, int M, int N>
  inline FieldDiagonal<T,M>& row_equilibration(FieldDiagonal<T,M>& Ri, const Dune::FieldMatrix<T,M,N>& A){
    T r_i = T(); 
    T twolgthm = T();
    int m; 
    for(int i = 0; i < M; ++i){
      r_i = T();                     //reset
      for(int j = 0; j < N; ++j){      
	r_i = std::max(r_i,Abs(A[i][j]));  //determine max_{0<=j<N}(|a_ij|)
      }
      //adonis_assert(r_i != T());  //since log_2 0 isn't defined!
      if(r_i == T())
	twolgthm = T();
      else 
	twolgthm = log2(r_i);

      m = round(-twolgthm);
      r_i = pow_of_machine_base<double,2>(m);
      
      Ri[i] = r_i;
    }
    return Ri;
  }
 

  /**
   *\brief Preconditioning by column-equilibration. This is similar to row equilibration except that the columns are regarded.
   */
  template<class T, int M, int N>
  inline void column_equilibration(Dune::FieldMatrix<T,M,N>& A){
    T r_j = T(); 
    T twolgthm = T();
    int m; 
    for(int j = 0; j < N; ++j){
      r_j = T();                     //reset
      for(int i = 0; i < M; ++i){      
	r_j = std::max(r_j,Abs(A[i][j]));  //determine max_{0<=j<N}(|a_ij|)
      }
      //adonis_assert(r_j != T());  //since log_2 0 isn't defined!
      
      if(r_j == T())
	twolgthm = T();
      else
	twolgthm = log2(r_j);

      //std::cout << "max_"<<"|a|_ij = " << r_j << std::endl;

      //log_2 1/max(|a_ij|) = log_2 1 - log_2 max(|a_ij|) = - log_2 max(|a_ij|)
      //Hence m ~ -log_2 max(|a_ij|)
      m = round(-twolgthm);
      r_j = pow_of_machine_base<double,2>(m);
	
      //std::cout << "r_j = 1./max|a_ij| ~ 2^m = "<<r_j << std::endl;
      
      for(int i = 0; i < M; ++i){      // multiply column i by r 
	A[i][j] *= r_j;
      }
    }
  }


  template<class T, int M, int N>
  inline  FieldDiagonal<T,N>& column_equilibration(FieldDiagonal<T,N>& Rj, const Dune::FieldMatrix<T,M,N>& A){
    T r_j = T(); 
    T twolgthm = T();
    int m; 
    for(int j = 0; j < N; ++j){
      r_j = T();                     //reset
      for(int i = 0; i < M; ++i){      
	r_j = std::max(r_j,Abs(A[i][j]));  //determine max_{0<=j<N}(|a_ij|)
      }
      // adonis_assert(r_j != T());  //since log_2 0 isn't defined!
      
      if(r_j == T())
	twolgthm = T();
      else
	twolgthm = log2(r_j);

      //std::cout << "max_"<<"|a|_ij = " << r_j << std::endl;
      
      //log_2 1/max(|a_ij|) = log_2 1 - log_2 max(|a_ij|) = - log_2 max(|a_ij|)
      //Hence m ~ -log_2 max(|a_ij|)
      m = round(-log2(r_j));
      r_j = pow_of_machine_base<double,2>(m);
      
      Rj[j] = r_j;
    }
    
    return Rj;
  }


  /**
     \brief Do row-column equilibration
   */
  template<class T, int M, int N>
  inline void row_column_equilibration(Dune::FieldMatrix<T,M,N>& A){
    row_equilibration(A);
    column_equilibration(A);
  }

  /**
     \brief Do column-row equilibration
   */
  template<class T, int M, int N>
  inline void column_row_equilibration(Dune::FieldMatrix<T,M,N>& A){
    column_equilibration(A);
    row_equilibration(A);
  }


  /**
   * \brief Get the diagonal matrices R and C which store the row and column scalings 
   */
  template<class T, int M, int N>
  inline std::pair<FieldDiagonal<T,M>, FieldDiagonal<T,N> > get_R_and_C(const Dune::FieldMatrix<T,M,N>& A){
    std::pair<FieldDiagonal<T,M>, FieldDiagonal<T,N> > p;
    FieldDiagonal<T,M> R;
    FieldDiagonal<T,N> C;
    
    row_equilibration(R,A);
    column_equilibration(C,A);

    p.first = R;
    p.second = C;
    
    return p;
  }


  /** 
   * \brief Variable scaling \f$ x \leftarrow w\cdot x + s, \f$ where \f$w,s \in \mathbb{R}^n \f$ are the (vector) weights and shifts respectively.  
   */
  template<class X>
  inline void variable_scaling(X& x, const X& weights, const X& shifts){
    x = weights*x +shifts;
  }
  
  template<class X>
  inline void variable_scaling(X& x, const X& weights){
    x = weights*x;
  }
  
  /**
   * \brief When a post scaling of the variable is desired, just recover \f$x \f$ by \f$ x \leftarrow w(x-s) \f$.
   */
  template<class X>
  inline void post_scaling(X& x, X& weights, const X& shifts){
    invert(weights);
    x = weights*(x-shifts);
  }

  template<class X>
  inline void post_scaling(X& x, X& weights){
    invert(weights);
    x = weights*x;
  }


  /**
   * \brief Change the diagonal of a matrix of order N
   */
  template<template<class S> class OP, class T, int N>
  inline void upgrade_diagonal(Dune::FieldMatrix<T,N,N>& A, const T& val){
    for(int i = 0; i < N; ++i){
      OP<T>::apply(A[i][i],val);
    }
  }



  /**
   * \brief 
   */
  template<class T, int K, int NRHS>
  class DuneRightHandSideSelector{
  public:
    typedef Dune::FieldMatrix<T,K,NRHS> RHSType;
  };
  
  //!partial specialisation
  template<class T, int K>
  class DuneRightHandSideSelector<T,K,1>{
  public:
    typedef Dune::FieldVector<T,K> RHSType;
  };


  /**
   * \brief scale linear system \f$ A\cdot x = b, \quad A \in R^{m \times n}, b \in R^m  f$ via \f[ (D_1^{-1}AD_2)y = D_1^{-1}b\f] using e.g. Gaussion elimination and then setting \f$ x = D_2y\f$. Here \f$ D_1, D_2\f$ are two appropriate nonsingular diagonal matrices, e.g. \f[ D_1 = \diag{\beta^{r_i}_{i = 1}^{m}, D_2 = \diag{\beta^{c_j}_{j = 1}^{n}, \f] and \f$ \beta\f$ denotes the machine precision, see [1]
   *
   * References: 
   *
   *  [1] [GOLUB,VAN LOAN, <I> Matrix Computations </I> 3rd, § 3.5.2]
   *  [2] [HIGHAM, <I> Accuracy and Stability of Numerical Algorithms </I>  2nd, § 9.8]
   */
  template<class T, int M, int N>
  inline void scale_linear_system( Dune::FieldMatrix<T,M,N>& Anew, Dune::FieldVector<T,M>& bnew, const Dune::FieldMatrix<T,M,N>& A, const Dune::FieldVector<T,M>& b, const FieldDiagonal<T,M>& D1, const FieldDiagonal<T,N>& D2){
    
    FieldDiagonal<T,M> D1_minus1(D1);
    D1_minus1.invert();
    
    Anew = D1_minus1*A*D2;
    bnew = D1_minus1*b;
  }

  
  template<class T, int N>
  inline Dune::FieldVector<T,N>& retrieve_solution_after_scaling(Dune::FieldVector<T,N>& x, const Dune::FieldVector<T,N>& y, const FieldDiagonal<T,N>& D2){
    return (x = D2*y);
  }


  /**
   *\brief Tikhonov (Tychonov, Tychonoff) Regularisation (sometimes referred to as ridge regression). Let \f[ AX = B \f] a system where \f$ A \in \mathbb{R}^{m \times n}\f$ is <B>rank-decicient</B>. Then the simples Tikhonov regularisation reads as  \f[ \min \{\| AX - B\|^2 + \lambda^2\| LX\|^2 \}. \f]
   * This is equivalent to the following system (here set the regularistion parameter to be the identity, i.e. \f$ L := I\f$) \f[A^TA+\lambda^2I)x_{\lambda} = A^TB. \f]
   *
   * See [HANSEN, "Rank-Deficient & Discrete Ill-Posed Problems", SIAM, 1987, Ch. 5]
   *
   * Let \f$ A \f$ be a matrix of order \f$ n, \f$ where \f$ A \f$ is rank-deficient. Then perform \f[ A^T\cdot A + \lambdaI = A^T\cdot B\f] to yield a system which becomes p.d.
   * The system then is always of the form \f$ KY = F, \f$ with \f$ K \in \mathbb{R}^{N \times N}, X,F \in \mathbb{R}^{N\times NRHS}.\f$ 
   *
   * NOTE 1: Determining the optimal regularisation parameter is another issue, but for the standard Tikhonov regularisation it can be taken to be the corner of the generic L-curve (\f$\lambda \in \mathcal{O}(1.e-2) \f$), 
   * cf. [HANSEN, "The L-curve and its use in the numerical treatment of inverse problems"]
   * NOTE 2: in case of the regularisation operator  \f$ L \f$ not being the identity \f$ I \f$, the most efficent way to compute the solution of the Tikhonov regularisation is to use Eld{\'e}n's algorithm, 
   * cf. [ELDEN, "Algorithms for the regularisation of ill-conditioned least squares problems"] 
   * When \f$\lambda \f$ is appropriately chosen, the solution to Tikhonov system should <B>coincide quite well</B> with the solution obtained from LAPACK's xgelsd_.
   */
  template<class T, int M, int N, int NRHS, 
	   template<class S> class OP = AddBasicElements,
	   class DetType = double>
  class TikhonovRegularisation{
  public:
    typedef Dune::FieldMatrix<T,M,N> LHSType;
    typedef Dune::FieldMatrix<T,N,M> TransposedLHSType;
    typedef Dune::FieldMatrix<T,N,N> TikhSystemType;
    typedef Dune::FieldMatrix<T,((M > N) ? M : N),NRHS> RHSType;
    typedef Dune::FieldMatrix<T,N,NRHS> SolutionType;
    

    TikhonovRegularisation(const LHSType& A, const RHSType& B, const T& lambda = 0.021568):A_(A), B_(B), lambda_(lambda), isRegularized_(false){}

    const LHSType& A() {return A_;}
    const RHSType& B() {return B_;}

    T& lambda() {return lambda_;}
    const T& lambda() const {return lambda_;}

    TikhSystemType& tikhonov_lhs_system() {
      adonis_assert(isRegularized_);
      return Sy_;
    }
    
    const TikhSystemType& tikhonov_lhs_system() const{
      adonis_assert(isRegularized_);
      return Sy_;
    }

    SolutionType& tikhonov_rhs_system() {
      adonis_assert(isRegularized_);
      return X_;
    }
    
    const SolutionType& tikhonov_rhs_system() const{
      adonis_assert(isRegularized_);
      return X_;
    }
    
    SolutionType& solution() {return tikhonov_rhs_system();}
    const SolutionType& solution() const {return tikhonov_rhs_system();}

    inline void regularize(){
      TransposedLHSType ATrans = transpose(A_);
      Sy_ = ATrans*A_;                         //A^T·A
      upgrade_diagonal<OP>(Sy_,lambda_);       //A^T·A +lambda·I
      X_ = ATrans*B_;                          //update rhs type 
      isRegularized_ = true;                   //now system is regularized
    }

    inline T determinant() const {
      adonis_assert(isRegularized_);
      return Sy_.determinant();  
    }

    //! system is symmetric an should be non-singular and can be solved using
    //! LAPACK's xSYSV routine
    inline void solve(bool isPd = true, char rc = 'c'){
      //adonis_assert(isRegularized_ && !is_zero(Sy_.determinant()));
      //if A p.d. ==> det(A) > 0
      adonis_assert(isRegularized_ && is_greater_than_zero(Sy_.determinant()));

      if(!isPd)
	solve_symmetric_ls(Sy_,X_,'U',rc);  //X_ contains solution x_lambda now
      else{
	std::cout << std::endl<< "Start solving p.d. system via "<<typeid(T).name()<<"posv_:" << std::endl;
	solve_positive_definite_ls(Sy_,X_,'U',rc);
      }
    }
    
  private:
    const LHSType& A_;
    const RHSType& B_;
    T lambda_;
    TikhSystemType Sy_;
    SolutionType X_;
    mutable bool isRegularized_;
  };
  

}//end of namespace 

#endif
