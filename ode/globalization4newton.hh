#ifndef GLOBALIZATION_STRATEGY_4_NEWTON_S_METHOD
#define GLOBALIZATION_STRATEGY_4_NEWTON_S_METHOD

#include "../common/globalfunctions.hh"
#include "../common/typeadapter.hh"
#include "../misc/operations4randomaccesscontainers.hh"
#include "../common/fancymessages.hh"
#include "../misc/misctmps.hh"
#include "../random/randomnumbergenerators.hh"

namespace Adonis{

  /**
   * \brief Globalization is done via an Armijo line search, cf. [1] and [2].

   * References:
   * [1]  [KELLY, "Iterative methods for linear and nonlinear equations, Ch. 8"]
   * [2] The programs accompanying the book are freely available on the web:
   <A HREF="http://de.mathworks.com/matlabcentral/fileexchange/2198-iterative-methods-for-linear-and-nonlinear-equations/content/kelley/">Programs accompanying the book</A>
   *
   * \brief globalization for Newton with Armijo line search
   * forward declaration
   * See [1, Ch. 8, Algo. 8.2.1, step 2. (b), p. 139] or 
   *     [2, program nsola.m, l. 140 -- 182]
   */
  template<class V, class REPAIRER, bool GLOBALIZEME> class Globalization4Newton;

  template<class V, class REPAIRER>
  class Globalization4Newton<V,REPAIRER,true>{
  public:
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef Norm<'2',BaseType> NormType;
    typedef std::size_t SizeType;
    
    Globalization4Newton(REPAIRER& Repairer, SizeType dim, 
			 const BaseType& alpha = 1.e-04, 
			 const BaseType& sigma0 = 0.1, 
			 const BaseType& sigma1 = 0.5, 
			 unsigned maxarm = 50, 
			 const BaseType& lambdaMin = 0.25
			 ):Repairer_(Repairer), yOld_(dim), xt_(dim), ft_(dim), f0_(dim), alpha_(alpha), sigma0_(sigma0), sigma1_(sigma1), lambda_(1.), lamm_(1.), lamc_(1.), lambdaMin_(lambdaMin), nft_(0.), nf0_(0.), ff0_(0.), ffc_(0.), ffm_(0.), p3res_(1.), c2_(0.), iarm_(0), maxarm_(maxarm), isSet_(false){
      FancyMessages FM;
      FM.nice_output("GLOBALIZATION 4 NEWTON'S METHOD APPLIED...",36);
    }

    void get_f0(int timestep, const V& g){
      if(timestep == 0){
	f0_ = -g;
      } 
      isSet_ = true;
    }

    /**
     * \brief Invoke the function inside the Newton iteration
     * \tparam NLEQ, the nonlinear function for \f$ G(y_n) = 0 \f$
     * \param y_n [input, output] current iterate
     * \param step, solution of Newton system \f$G'(y_n)·step = -G(y_n)\f$
     */
    template<class NLEQ>
    void armijo_line_search(int timestep, V& y_n, const V& step, NLEQ& G, int excessSpecIndex, bool& isStationary){
      adonis_assert(isSet_==true);
      bool flag(false);
      unsigned countRefinement(0);

      //TEST TEST TEST
      // std::cout << "update that SHIT" << std::endl;
      // y_n += step;
      //TEST TEST TEST
      
      yOld_ = y_n;
      //reset 
      lambda_ = lamm_ = 1.;
      lamc_ = lambda_;
      iarm_ = 0;
      xt_ = y_n + lambda_*step;

      //repair action to be done (only for 1D and 2D MOL relevant
      Repairer_.repair(xt_,excessSpecIndex);
      Repairer_.info(timestep, -1, " Armijo LS globalization: \n");
      
      ft_ = G(xt_, isStationary);
      nft_ = NormType::norm(ft_); nf0_ = NormType::norm(f0_); 
      ff0_ = ntimes<2>(nf0_); ffc_ = ntimes<2>(nft_); ffm_ = ntimes<2>(nft_);
      
      while((nft_ >= (1 - alpha_*lambda_)*nf0_)){
	// apply the three point parabolic model
	if (iarm_ == 0)
            lambda_ *= sigma1_; 
        else
	  lambda_ = parab3p(lamc_, lamm_, ff0_, ffc_, ffm_); 
        

	//update xt_; keep the books on lambda
	xt_ = y_n + lambda_*step;
        lamm_ = lamc_;
        lamc_ = lambda_;

	//keep the books on the function norms
	//repair action to be done (only for 1D and 2D MOL relevant
	Repairer_.repair(xt_,excessSpecIndex);
	Repairer_.info(timestep, -1, (" Armijo LS globalization in lambda refinement step"+Num2str(countRefinement)));
	ft_ = G(xt_, isStationary);
	ffm_ = ffc_;
        ffc_ = ntimes<2>(nft_);
        iarm_++;
	if(iarm_ > maxarm_){
	  ADONIS_CAUTION(Caution, "Armijo failure, too many reductions ( maxarm = "<<maxarm_+1<<").");
	  flag = true;
	  break; // reset and leave loop 
	}
	countRefinement++;  //increment refinement counter
      } //end while
      
      if(flag == true){
	y_n = yOld_;
      }
      else{
	y_n = xt_;
	f0_ = ft_;
      }
      
      std::cout << "    Globalization lambda = " << lambda_ << "  (#refinements = "<< countRefinement << ")." << std::endl;
    } //end of line search



  private:
    REPAIRER& Repairer_;
    V yOld_, xt_, ft_, f0_;
    BaseType alpha_, sigma0_, sigma1_,
      lambda_, lamm_, lamc_, lambdaMin_, nft_, nf0_, ff0_, ffc_, ffm_, 
      p3res_, c2_;
    unsigned iarm_, maxarm_;
    bool isSet_;
    
    RandomNumberGenerator<BaseType,32,'L'> rng_;
    
    //! private function
    BaseType& parab3p(const BaseType& lambdac, const BaseType& lambdam, const BaseType& ff0, const BaseType& ffc, const BaseType& ffm){
      //compute coefficients of interpolation polynomial, 
      
      c2_ = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
      if(c2_ >= 0){
	p3res_ = sigma1_*lambdac;
      }
      else{ //execute rest of function
	p3res_ = -(lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0))*.5/c2_;
	if (p3res_ < sigma0_*lambdac)
	  p3res_=sigma0_*lambdac;
	if (p3res_ > sigma1_*lambdac) 
	  p3res_=sigma1_*lambdac;
      }
      return p3res_;
    }


  };

  // don't apply any Newton globalization via LS -- Do nothing
   template<class V, class REPAIRER>
  class Globalization4Newton<V,REPAIRER,false>{
  public:
    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef Norm<'2',BaseType> NormType;
    typedef std::size_t SizeType;
    
     Globalization4Newton(REPAIRER& Repairer, SizeType dim, const BaseType& alpha = 1.e-04, const BaseType& sigma0 = 0.1, const BaseType& sigma1 = 0.5, unsigned maxarm = 50){}

     void get_f0(int timestep, const V& g){}

     template<class NLEQ>
     void armijo_line_search(int timestep, V& y_n, const V& step, NLEQ& G, int excessSpecIndex, bool& isStationary){}

   };

} //end namespace

#endif
