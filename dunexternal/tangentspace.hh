#include "dunextensions.hh"
#include "../misc/operations4randomaccesscontainers.hh"

namespace Adonis{

  /**
   * \brief Compute tangent space according to [1, Appendix 1, pp. 847]
   *
   * Reference:
   *
   * [1] [Ren, Goldin, Hiremath and Pope, "Reduced description of reactive flows with tabulation of chemistry", Combustion Theory and Modelling, 15, 2011, (827--848)]
   */
  template<class T, int REDDIM, int FULLDIM>
  class TangentSpaceRobust{
  public:
    typedef Dune::FieldMatrix<T,FULLDIM,REDDIM> BType;
    typedef Dune::FieldMatrix<T,REDDIM,REDDIM> RedMtxType;
    typedef Dune::FieldMatrix<T,FULLDIM,FULLDIM> FullMtxType;

    TangentSpaceRobust(const BType& B):B_(B){}

    const BType& get_B() const {return B_;}

    template<class V, class W>
    BType compute(const V& r, const W& zM){
      //! compute M, i.e. [1,eq. (16)]
      RedMtxType rrt = dyadic_product<REDDIM,REDDIM>(r,r),
	I,
	A;
      T d = dot(r,r);
      adonis_assert(!is_zero(d));

      rrt /= d;
     
      I = T();

      for(int i = 0; i < REDDIM; ++i)
	I[i][i] = 1.;

      A = I - rrt;
        
      //! we don't set up Z but rather calculate it more qickly using zM instead
      BType Mtx = dyadic_product<FULLDIM,REDDIM>(zM,r) + diagonal_matrix_multiplication(zM,B_)*A;

      //! W as QR decomposition of Mtx
      BType Wtx  = extract_Q1<REDDIM>(orthonormalize(Mtx));
      RedMtxType D = transpose(B_)*Wtx;
      inversion(D);  //!now D is overwritten by its inverse \f$(B^TW)^{-1}\f$
      return Wtx*D;  //!formula (20) from [1]  
    }

   
  private:
    const BType& B_;
  };



} //end namespace
