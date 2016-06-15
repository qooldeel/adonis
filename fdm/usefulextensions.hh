#ifndef USE_FUL_EXTENSIONS_4_2D_REAC_FLOWS_HH
#define USE_FUL_EXTENSIONS_4_2D_REAC_FLOWS_HH

namespace Adonis{
  
  template<int OFS(int,int,int), int> class FiDi2D;

  template<int OFS(int,int,int)> //well uses a function pointer here
  class FiDi2D<OFS,1>{           //! first order FD derivative
  public:
    
    template<class V, class H>
    static inline typename V::value_type fd_downstream_x(int i, int j, const V& u, int s, const H& h){
      return ( (u[OFS(i+1,j,s)] -u[OFS(i,j,s)])/h ); 
    }


    template<class V, class H>
    static inline typename V::value_type fd_downstream_y(int i, int j, const V& u, int s, const H& h){
      return ( (u[OFS(i,j+1,s)] -u[OFS(i,j,s)])/h ); 
    }
  };

   //!  2nd order FD derivative
  template<int OFS(int,int,int)> //well uses a function pointer here
  class FiDi2D<OFS,2>{          
  public:
    
    template<class V, class H>
    static inline typename V::value_type fd_downstream_x(int i, int j, const V& u, int s, const H& h){
      return ( (u[OFS(i+1,j,s)] -u[OFS(i-1,j,s)])/(2.*h) ); 
    }


    template<class V, class H>
    static inline typename V::value_type fd_downstream_y(int i, int j, const V& u, int s, const H& h){
      return ( (u[OFS(i,j+1,s)] -u[OFS(i,j-1,s)])/(2.*h) ); 
    }
  };

}
#endif
