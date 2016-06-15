#ifndef DAMPED_NEWTON_PROCEDURES_HH
#define DAMPED_NEWTON_PROCEDURES_HH

#include "../common/globalfunctions.hh"
#include "../io/readinparameters.hh"
#include "../fdm/gensettings.hh"
#include "constants.hh"

namespace Adonis{

  /**
   * \brief Calculate Newton step with lambda, starting from 1. Lower it if 
   *  mass fraction exceeds [0,1]
   * \tparam V Vector tpye
   * \tparam DIM dimension of MOL problem 
   * \tparam DAMPME false = don't damp, true = damping on
   */
  template<class V, int DIM, bool DAMPME>  //! default situation: do nothing
  class SimpleDampedNewtonMassFraction{
  public:
    typedef std::string StringType;
    typedef typename V::value_type value_type;

    SimpleDampedNewtonMassFraction(const value_type& lambda = value_type(1), const value_type& lambdaMin = value_type(1.e-10)){}
    
    value_type damping(bool isViolated) {return value_type(1);}
    void reset_lambda() {} //do nothing so far
    void info() const {} //do nothing 
    
    V& update_Newton_iterate(V& y_n, const V& g){
      return (y_n += g);
    }
  private:
  };


  

  //! partial specialization

  //! use damping 
  //!Only for 2D MOL (2) AND when damping is switch on (true)
  template<class V>                  
  class SimpleDampedNewtonMassFraction<V,2,true>{
  public:
    typedef std::string StringType;
    typedef typename V::value_type value_type;

    SimpleDampedNewtonMassFraction(const value_type& lambda = value_type(1), const value_type& lambdaMin = value_type(1.e-10)):lambda_(lambda),lambdaMin_(lambdaMin){}
    

    value_type& damping(bool isViolated){
      //some physical quants, mainly T and Y_k, were violated
      if(isViolated==true){
	lambda_ /= Constant<unsigned int>::dampFac;  //half lambda
      }
      
      lambda_ = Max(lambda_,lambdaMin_); //don't make lambda too small
      std::cout << "   lambda = "<< lambda_ << std::endl; //print info
      return lambda_;
    }

    void reset_lambda() {lambda_ = value_type(1);} //set back to 1

    void info() const{
      std::cout << "****DAMPED NEWTON for 2D MOL used for Y_k < 0 OR Y_k > 1. " << std::endl;
    }

    V& update_Newton_iterate(V& y_n, const V& g){
      return (y_n += lambda_*g);
    }

  private:
    value_type lambda_, lambdaMin_;
  };


 
  /* template<class V, int DIM, bool DAMPME> 
  class SimpleDampedNewtonMassFraction{
  public:
    typedef std::string StringType;
    typedef typename V::value_type value_type;

    SimpleDampedNewtonMassFraction(const StringType& molSettingFile, V& y, const value_type& lambda = value_type(1), const value_type& lambdaMin = value_type(1.e-10)){}
    
    value_type damping() {return value_type(1);}

  private:
  };


  

  //! partial specialization

  //! use damping 
  template<class V>                   //2D case 
  class SimpleDampedNewtonMassFraction<V,2,true>{
  public:
    typedef std::string StringType;
     typedef typename V::value_type value_type;

    SimpleDampedNewtonMassFraction(const StringType& molSettingFile, V& y, const value_type& lambda = value_type(1), const value_type& lambdaMin = value_type(1.e-10)):y_(y),lambda_(lambda),lambdaMin_(lambdaMin){

      ParameterData Param;
      Param.read_from_file(molSettingFile);
      Nx_ = Param.get_datum<int>("Nx");
      Ny_ = Param.get_datum<int>("Ny");
      nspec_ = (int)y_.size()/(Nx_*Ny_)-4;
    }
    

    value_type& damping(){
      bool YkBeyondBounds(false);
      for(int i = 0; i < Nx_; ++i){
	for(int j  = 0; j < Ny_; ++j){
	  for(int k = 0; k < nspec_; ++k){
	    if(//if Y_K(i,j) is not in [0.,1.]
#ifdef NONCONSERVATIVE_FORM
	       is_contained( y_[i + Nx_*j + (4+k)*Nx_*Ny_],0.,1.) == false 
#else
	       is_contained((y_[i + Nx_*j + (4+k)*Nx_*Ny_]/y[i + Nx_*j]),0.,1.) == false
#endif
	       ){
	      YkBeyondBounds = true;
	      break;
	    }
	    else
	      YkBeyondBounds = false;
		 
	  }
	}
      }
      if(YkBeyondBounds)
	lambda_ /= 2;  //half lambda
      lambda_ = Max(lambda_,lambdaMin_); //don't make lambda too small
      
      return lambda_;
    }

  private:
    V& y_;  //reference 
    value_type lambda_, lambdaMin_;
      int Nx_, Ny_, nspec_;
  };
 */

} //end namespace

#endif
