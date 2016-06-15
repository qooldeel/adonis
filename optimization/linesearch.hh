#ifndef LINE_SEARCH_HH
#define LINE_SEARCH_HH

#include "../common/adonisassert.hh"
#include "../common/globalfunctions.hh"
#include "../expressiontemplates/expressionvectortraits.hh"
#include "../expressiontemplates/scalar.hh"
#include "../containers/densesymmetricmatrix.hh"
#include "../common/elementaryoperations.hh"
#include "../common/typeadapter.hh"
#include "../misc/operations4randomaccesscontainers.hh"

namespace Adonis{

  /**
   * \ brief <B>Idea of Line Search:</B> given a descent direction \f$ p_k \f$, take a step in that direction that yields an "acceptable" \f$ x_{k+1}\f$, i.e.
   *
   * \f$ \nu = 0, 1, 2, \ldots\f$
   *
   *    (1) Calculate search direction \f$ p_k\f$
   *
   *    (2) \f$ x_{k+1} \leftarrow x_k + \alpha_k p_k\f$ for some \f$ \alpha_k > 0\f$ that makes \f$ x_{k+1} \f$ an acceptable next iterate (see backtracking)
   *
   * Usage: \f$ \min_x f(x), \ f: \mathbb{R}^n \longrightarrow \mathbb{R}\f$
   *      for systems, e.g. \f$ F(x) = 0, \ F: \mathbb{R}^n \longrightarrow \mathbb{R}^n \f$ use \f$ f(x) := \frac{1}{2}\|F(x) \|_2^2 = \frac{1}{2}F(x)^t\cdot F(x).\f$
   *
   * Convergence test (aka stopping criterion): In [FLETCHER, "Practical methods of optimization", 2nd ed, p 22/23] it is suggested to use \f$ \|\nabla f(x_k)\| \leq TOL\f$ or \f$ f(x_k) - f(x_{k+1}) \leq TOL, \f$ where \f$ TOL > 0\f$ is a  <I> user-defined</I> tolerance.
    NOTE: Backtracking unsing first Wolfe condition. This is used for <I>globalisation </I>
    *
    * See 1.) [NOCEDAL/WRIGHT, "Numerical Optimization", 2nd ed, p. 37]
    *     2.) [DENNIS/SCHNABEL, "Numerical Methods for Unconstrained Optimization and Nonlinear Equations", p. 126]
  */
  template<class T, template<class S> class FUN,
	   class NORM = Norm<'2',typename TypeAdapter<T>::BaseType> >
  class LineSearchWithBackTracking{
  public:
    typedef FUN<T> FunType;
    typedef typename ExprVecTraits<T>::ExprVecType ExprVecType;
    typedef DenseSymmetricMatrix<T,ExprTmpl::MyVec,SymmetricStorageByDiagonal<unsigned> > SymmType;
    
    typedef NORM NormType;

    LineSearchWithBackTracking(const T& alpha = 1., const T& rho = 0.1, const T& c = 1e-4):alpha_(alpha), rho_(rho), c_(c), fun_(1){
      adonis_assert(rho > 0 && rho < 1 && c > 0 && c < 1);
    }

    //! backtracking; often sophisticated in determining a nonfixed rho and difficult to code -- ARMIJO 
    T determine_step(const ExprVecType& x, const ExprVecType& gradf, const ExprVecType& p){
      T alpha = alpha_;
      
      while(fun_(x + alpha*p) > fun_(x) + c_*alpha*dot(gradf,p)){

	//std::cout << "alpha = "<<alpha << "    f(x+alpha*p) = "<< f(xxpress + alpha*pxpress) << std::endl;


	alpha *= rho_;   //!for some \f$ \rho in (0,1)\f$ 
	//rho is (ideally) chosen anew each time by the line search
      }
      return alpha;
    }

    

    /**
     * \brief At step \f$ \f$ calculate \f$ \nu \f$ for making \f$\nabla^2 f \f$ p.d.
     * This is needed when the Hessian is not p.d. 
     *
     * Direction \f$ p \f$ is altered
     */
    void determine_nu(const T& nu0, ExprVecType& x, SymmType& hessf, const ExprVecType& gradf, const T& lambda, ExprVecType& p){
      adonis_assert(nu0 >= T());

      //ExprVecType p(x.size());

      T nu = nu0;
      if(nu == 0) //otherwise we would loop without stopping
	nu = 0.025;  
      
      //(ii) !UNFURTUNATELY, \f$ \nabla^2 f\f$ must be positive definite!
      while(!hessf.is_positive_definite()){
	hessf.update_diagonal<AddBasicElements<T> >(nu); // G + nuI
	nu *= 4;
      }
      
      //(iii)
      p = -1.*gradf;
      hessf.solve_pos_def(p,1);

      //(iv)
      T r = (fun_(x) - fun_(x+lambda*p))/(fun_(x) - 0.5*dot(p*hessf,p) + dot(gradf,p));
      
      //(v)
      if(r < 0.25)
	nu *= 4;
      if(r > 0.75)
	nu *= 0.5;

      //(vi)
      //if(r > 0)
      //x += lambda*p;  //update
    }
    
    /**
     *\brief the line search with (simple) backtracking
     *
     * NOTE: To work with \f$ p_k : = -\nabla^2 f(x_k)^{-1}\cdot \nabla f(x_k)\f$, the Hessian \f$ \nabla^2 f(x_k) \f$ has to be positive definite, i.e. all its eigenvalues must be strictly positive!  
     */
    ExprVecType loop(unsigned maxiteration, const ExprVecType& x0, const T& tol, const T& nu = T()){
      adonis_assert(tol > 0 && tol < 1); //a meaningful tolerance

      unsigned dim = x0.size();
      // std::cout<< "dim = "<< dim << std::endl;

      adonis_assert(dim == fun_.domain_dim());
      
       ExprVecType x(x0),  //iterate
	 gradf(dim),       //gradient of fun_ 
	 p(dim);               //search direction
      SymmType hessf(dim); //hessian of f
      
      T lambda = 1.;       //step size for line search

      unsigned itercount = 0;
      for(unsigned k = 0; k < maxiteration; ++k){
	
	//(1) compute search direction pk
	gradf = fun_.gradient(x);
     
	determine_nu(nu,x,fun_.hessian(x),gradf,lambda,p);

	//p = -1.*gradf;
	//(fun_.hessian(x)).solve(p,1); 
	
	//(2) calculate "acceptable" step... 
	lambda = determine_step(x,gradf,p);
	std::cout << k <<".)      lambda_k = "<< lambda << std::endl; 
	
	//...and update along the line
	x += lambda*p;

	if(NormType::norm(gradf) <= tol){
	  itercount = k+1;
	  break;  //leave loop
	}
	if(k == maxiteration-1)
	  ADONIS_INFO(Information,"max. number of iterations reached without convergence");
      }
    
      std::cout << "Global minimum achieved within "<< itercount << " iterations." << std::endl;
      
      return x;
    }

  private:
    T alpha_,
      rho_,
      c_;
    FunType fun_;
  };


} //end namespace 

#endif
