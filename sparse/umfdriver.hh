#ifndef UMFPACK_PRIMARY_ROUTINES_HH
#define UMFPACK_PRIMARY_ROUTINES_HH

#include <string>
#include <complex>

#if USE_UMFPACK
#include "umfpack.h"
#include "umftypetraits.hh"
// #else
// #error "U need USE_UMFPACK set to 1"
#endif

#include "../common/numerictypechecker.hh"


namespace Adonis{

#if USE_UMFPACK
  /**
     \brief Primary routines of UMFPACK, i.e. routines that can be solved to 
     solve a nonsymmetric square system \f$ A\cdot = b.\f$
     
     UMFPACK has 4 different such primary routines, viz.:
     1.) umfpack_di_fun     (double, int)
     2.) umfpack_dl_fun     (double, SuiteSparse_long)
     3.) umfpack_zi_fun     (complex double, int)
     4.) umfpack_zl_fun     (complex double, SuiteSparse_long),
     
     where fun is <I>any</I> routine from the UMFPACK package. 
  */
  template<class T, class INT> class UmfpackDriver;

  // //!specializations
  //!1.) 
  template<>
  class UmfpackDriver<double,int>{
  public:
    typedef double value_type; 
    typedef double BaseType;
    typedef int IntegerType;

    static const IntegerType signature = 1;

    static inline std::string routine_type() {return "di";} 

    //! for real version, az is a dummy argument and you can just pass 0 
    //! as function argument
    static inline IntegerType symbolic_analysis(IntegerType rows, IntegerType cols, const IntegerType* p, const IntegerType* ix, const double* a, const double* az, void** symbolic, const double* ctrl, double* info){ 
      return umfpack_di_symbolic(rows,cols,p,ix,a,symbolic,ctrl,info);
    }
 
    static inline void free_symbolic(void** symbolic){
      umfpack_di_free_symbolic(symbolic);
    }
  
    
    static inline IntegerType lu_decomposition(const IntegerType* p, const IntegerType* ix, const double* ax, const double* az, void* symbolic, void** numeric, const double* ctrl, double* info){  //az is a dummy here. see above
      return umfpack_di_numeric(p,ix,ax,symbolic,numeric,ctrl,info);
    }
  
    static inline void free_numeric(void** numeric){
      umfpack_di_free_numeric(numeric);
    }

    //! the placeholders for the imag. part az, xz and bz are meaningless here
    static inline IntegerType solve(IntegerType sys, const IntegerType* p, const IntegerType* ix, const double* ax, const double* az, double* x, const double* xz, const double* b, const double* bz, void* numeric, const double* ctrl, double* info){
      return umfpack_di_solve(sys,p,ix,ax,x,b,numeric,ctrl,info);
    }

    //! right dimensions for wsolve additional arguments
    static inline IntegerType w_dim(IntegerType n, IntegerType sys, const double* Control, const double* Info){
      IntegerType d = 0;
      if( ((sys == UMFPACK_A) || (sys == UMFPACK_At) || (sys == UMFPACK_Aat)) && (Control[UMFPACK_IRSTEP] > 0) //&& (Info[UMFPACK_STATUS] !=  UMFPACK_WARNING_singular_matrix) 
){ //iterative refinement is performed
	d = 5*n;
	//std::cout << "5*n"<<std::endl;
      }
      else {//no iterative refinement
	d = n;
	//std::cout << "n" << std::endl;
      }
      return d;
    }

    //! solution can be speeded up for iterative procedures
    static inline IntegerType wsolve(IntegerType sys, const IntegerType* p, const IntegerType* ix, const double* ax, const double* az, double* x, const double* xz, const double* b, const double* bz, void* numeric, const double* ctrl, double* info, IntegerType* Wi, double* W){
      return umfpack_di_wsolve(sys,p,ix,ax,x,b,numeric,ctrl,info,Wi,W);
    }

    //! load UMFPACK's control parameters. You can then modify individual
    //! settings by changing specific entries in the Control array. 
    static inline void defaults(double* Control){
      umfpack_di_defaults(Control);
    }

    //!get dimension of LU factor 
    static inline IntegerType get_lunz(IntegerType* lnz, IntegerType* unz, IntegerType* n_row, IntegerType* n_col, IntegerType* nz_udiag, void* Numeric){
      return umfpack_di_get_lunz(lnz,unz,n_row,n_col,nz_udiag,Numeric);
    }

    //! contains arrays with dimensions determined by get_lunz
    static inline IntegerType get_numeric(IntegerType* Lp, IntegerType* Lj, double* Lx, double* Lz, IntegerType* Up, IntegerType* Ui, double* Ux, double* Uz, IntegerType* P, IntegerType* Q, double* Dx, double* Dz, IntegerType* do_recip, double* Rs, void* Numeric){
      return umfpack_di_get_numeric(Lp,Lj,Lx,Up,Ui,Ux,P,Q,Dx,do_recip,Rs,Numeric);
    }
    
  };

  //2.)
  template<>
  class UmfpackDriver<double,SuiteSparse_long>{
  public:
    typedef double value_type; 
    typedef double BaseType;
    typedef SuiteSparse_long IntegerType;

    static const SuiteSparse_long signature = 2;

    static inline std::string routine_type() {return "dl";} 

    //! for real version, az is a dummy argument and you can just pass 0 
    //! as function argument
    static inline IntegerType symbolic_analysis(IntegerType rows, IntegerType cols, const IntegerType* p, const IntegerType* ix, const double* a, const double* az, void** symbolic, const double* ctrl, double* info){ 
      return umfpack_dl_symbolic(rows,cols,p,ix,a,symbolic,ctrl,info);
    }
 
    static inline void free_symbolic(void** symbolic){
      umfpack_dl_free_symbolic(symbolic);
    }
  
    
    static inline IntegerType lu_decomposition(const IntegerType* p, const IntegerType* ix, const double* ax, const double* az, void* symbolic, void** numeric, const double* ctrl, double* info){  //az is a dummy here. see above
      return umfpack_dl_numeric(p,ix,ax,symbolic,numeric,ctrl,info);
    }
  
    static inline void free_numeric(void** numeric){
      umfpack_dl_free_numeric(numeric);
    }

    //! the complex placeholders az, xz and bz are meaningless here
    static inline IntegerType solve(IntegerType sys, const IntegerType* p, const IntegerType* ix, const double* ax, const double* az, double* x, const double* xz, const double* b, const double* bz, void* numeric, const double* ctrl, double* info){
      return umfpack_dl_solve(sys,p,ix,ax,x,b,numeric,ctrl,info);
    }

    //! right dimensions for wsolve additional arguments
    static inline IntegerType w_dim(IntegerType n, IntegerType sys, const double* Control, const double* Info){
      IntegerType d = 0;
      if( ((sys == UMFPACK_A) || (sys == UMFPACK_At) || (sys == UMFPACK_Aat)) && (Control[UMFPACK_IRSTEP] > 0) 
	  //&& (Info[UMFPACK_STATUS] !=  UMFPACK_WARNING_singular_matrix) 
){ //iterative refinement is performed
	d = 5*n;
      }
      else //no iterative refinement
	d = n;
      return d;
    }

    //! solution can be speeded up for iterative procedures
    static inline IntegerType wsolve(IntegerType sys, const IntegerType* p, const IntegerType* ix, const double* ax, const double* az, double* x, const double* xz, const double* b, const double* bz, void* numeric, const double* ctrl, double* info, IntegerType* Wi, double* W){
      return umfpack_dl_wsolve(sys,p,ix,ax,x,b,numeric,ctrl,info,Wi,W);
    }
    
    static inline void defaults(double* Control){
      umfpack_dl_defaults(Control);
    }

    //!get dimension of LU factor 
    static inline IntegerType get_lunz(IntegerType* lnz, IntegerType* unz, IntegerType* n_row, IntegerType* n_col, IntegerType* nz_udiag, void* Numeric){
      return umfpack_dl_get_lunz(lnz,unz,n_row,n_col,nz_udiag,Numeric);
    }

    //! contains arrays with dimensions determined by get_lunz
    static inline IntegerType get_numeric(IntegerType* Lp, IntegerType* Lj, double* Lx, double* Lz, IntegerType* Up, IntegerType* Ui, double* Ux, double* Uz, IntegerType* P, IntegerType* Q, double* Dx, double* Dz, IntegerType* do_recip, double* Rs, void* Numeric){
      return umfpack_dl_get_numeric(Lp,Lj,Lx,Up,Ui,Ux,P,Q,Dx,do_recip,Rs,Numeric);
    }

  };
 
  
  //3.) COMPLEX VERSIONS
  template<>
  class UmfpackDriver<std::complex<double>,int>{
  public:
    typedef std::complex<double> value_type; 
    typedef double BaseType;
    typedef int IntegerType;

    static const IntegerType signature = 3;

    static inline std::string routine_type() {return "zi";} 

    //! NO dummy arguments any more
    static inline IntegerType symbolic_analysis(IntegerType rows, IntegerType cols, const IntegerType* p, const IntegerType* ix, const BaseType* a, const BaseType* az, void** symbolic, const BaseType* ctrl, BaseType* info){ 
      return umfpack_zi_symbolic(rows,cols,p,ix,a,az,symbolic,ctrl,info);
    }
 
    static inline void free_symbolic(void** symbolic){
      umfpack_zi_free_symbolic(symbolic);
    }
  
    
    static inline IntegerType lu_decomposition(const IntegerType* p, const IntegerType* ix, const BaseType* ax, const BaseType* az, void* symbolic, void** numeric, const BaseType* ctrl, BaseType* info){ 
      return umfpack_zi_numeric(p,ix,ax,az,symbolic,numeric,ctrl,info);
    }
  
    static inline void free_numeric(void** numeric){
      umfpack_zi_free_numeric(numeric);
    }

    //! the  placeholders for the imag. part az, xz and bz become IMPORTANT now
    static inline IntegerType solve(IntegerType sys, const IntegerType* p, const IntegerType* ix, const BaseType* ax, const BaseType* az, BaseType* x,  BaseType* xz, const BaseType* b, const BaseType* bz, void* numeric, const BaseType* ctrl, BaseType* info){
      return umfpack_zi_solve(sys,p,ix,ax,az,x,xz,b,bz,numeric,ctrl,info);
    }

     //! right dimensions for wsolve additional arguments
    static inline IntegerType w_dim(IntegerType n, IntegerType sys, const double* Control, const double* Info){
      IntegerType d = 0;
      if( ((sys == UMFPACK_A) || (sys == UMFPACK_At) || (sys == UMFPACK_Aat)) && (Control[UMFPACK_IRSTEP] > 0) //&& (Info[UMFPACK_STATUS] !=  UMFPACK_WARNING_singular_matrix) 
	  ){ //iterative refinement is performed
	d = 10*n;
      }
      else //no iterative refinement
	d = 4*n;
      return d;
    }

    //! solution can be speeded up for iterative procedures
    static inline IntegerType wsolve(IntegerType sys, const IntegerType* p, const IntegerType* ix, const double* ax, const double* az, double* x, double* xz, const double* b, const double* bz, void* numeric, const double* ctrl, double* info, IntegerType* Wi, double* W){
      return umfpack_zi_wsolve(sys,p,ix,ax,az,x,xz,b,bz,numeric,ctrl,info,Wi,W);
    }
    
    static inline void defaults(double* Control){
      umfpack_zi_defaults(Control);
    }

    //!get dimension of LU factor 
    static inline IntegerType get_lunz(IntegerType* lnz, IntegerType* unz, IntegerType* n_row, IntegerType* n_col, IntegerType* nz_udiag, void* Numeric){
      return umfpack_zi_get_lunz(lnz,unz,n_row,n_col,nz_udiag,Numeric);
    }

    //! contains arrays with dimensions determined by get_lunz
    static inline IntegerType get_numeric(IntegerType* Lp, IntegerType* Lj, double* Lx, double* Lz, IntegerType* Up, IntegerType* Ui, double* Ux, double* Uz, IntegerType* P, IntegerType* Q, double* Dx, double* Dz, IntegerType* do_recip, double* Rs, void* Numeric){
      return umfpack_zi_get_numeric(Lp,Lj,Lx,Lz,Up,Ui,Ux,Uz,P,Q,Dx,Dz,do_recip,Rs,Numeric);
    }

  };


  //!4.)
  template<>
  class UmfpackDriver<std::complex<double>,SuiteSparse_long>{
  public:
    typedef std::complex<double> value_type; 
    typedef double BaseType;
    typedef SuiteSparse_long IntegerType;

    static const IntegerType signature = 4;

    static inline std::string routine_type() {return "zl";} 

    //! here NO dummy arguments any more
    static inline IntegerType symbolic_analysis(IntegerType rows, IntegerType cols, const IntegerType* p, const IntegerType* ix, const BaseType* a, const BaseType* az, void** symbolic, const BaseType* ctrl, BaseType* info){ 
      return umfpack_zl_symbolic(rows,cols,p,ix,a,az,symbolic,ctrl,info);
    }
 
    static inline void free_symbolic(void** symbolic){
      umfpack_zl_free_symbolic(symbolic);
    }
  
    
    static inline IntegerType lu_decomposition(const IntegerType* p, const IntegerType* ix, const BaseType* ax, const BaseType* az, void* symbolic, void** numeric, const BaseType* ctrl, BaseType* info){ 
      return umfpack_zl_numeric(p,ix,ax,az,symbolic,numeric,ctrl,info);
    }
  
    static inline void free_numeric(void** numeric){
      umfpack_zl_free_numeric(numeric);
    }

    //! the placeholders for the imag. part az, xz and bz have now a meaning
    static inline IntegerType solve(IntegerType sys, const IntegerType* p, const IntegerType* ix, const BaseType* ax, const BaseType* az, BaseType* x, BaseType* xz, const BaseType* b, const BaseType* bz, void* numeric, const BaseType* ctrl, BaseType* info){
      return umfpack_zl_solve(sys,p,ix,ax,az,x,xz,b,bz,numeric,ctrl,info);
    }

     //! right dimensions for wsolve additional arguments
    static inline IntegerType w_dim(IntegerType n, IntegerType sys, const double* Control, const double* Info){
      IntegerType d = 0;
      if( ((sys == UMFPACK_A) || (sys == UMFPACK_At) || (sys == UMFPACK_Aat)) && (Control[UMFPACK_IRSTEP] > 0) //&& (Info[UMFPACK_STATUS] !=  UMFPACK_WARNING_singular_matrix) 
){ //iterative refinement is performed
	d = 10*n;
      }
      else //no iterative refinement
	d = 4*n;
      return d;
    }

    //! solution can be speeded up for iterative procedures
    static inline IntegerType wsolve(IntegerType sys, const IntegerType* p, const IntegerType* ix, const double* ax, const double* az, double* x, double* xz, const double* b, const double* bz, void* numeric, const double* ctrl, double* info, IntegerType* Wi, double* W){
      return umfpack_zl_wsolve(sys,p,ix,ax,az,x,xz,b,bz,numeric,ctrl,info,Wi,W);
    }

    static inline void defaults(double* Control){
      umfpack_zl_defaults(Control);
    }
    
    //!get dimension of LU factor 
    static inline IntegerType get_lunz(IntegerType* lnz, IntegerType* unz, IntegerType* n_row, IntegerType* n_col, IntegerType* nz_udiag, void* Numeric){
      return umfpack_zl_get_lunz(lnz,unz,n_row,n_col,nz_udiag,Numeric);
    }

    //! contains arrays with dimensions determined by get_lunz
    static inline IntegerType get_numeric(IntegerType* Lp, IntegerType* Lj, double* Lx, double* Lz, IntegerType* Up, IntegerType* Ui, double* Ux, double* Uz, IntegerType* P, IntegerType* Q, double* Dx, double* Dz, IntegerType* do_recip, double* Rs, void* Numeric){
      return umfpack_zl_get_numeric(Lp,Lj,Lx,Lz,Up,Ui,Ux,Uz,P,Q,Dx,Dz,do_recip,Rs,Numeric);
    }

  };

#endif //end USE_UMFPACK 

} //end namespace 

#endif
