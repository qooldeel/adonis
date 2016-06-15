//! test double pointer
#include <iostream>
#include <string>


typedef struct{
 double rl;
 double im;
}tDCplx;


inline void init_cplx(tDCplx* z){
    z->rl = 0.0;
    z->im = 0.0;
}

inline void print_cplx(tDCplx* z, const std::string& str = std::string()){
    std::cout << ((str.size() == 0) ? "z" : str)<< " = "<< z->rl << " + " << z->im << "i" << std::endl;
    
}

//! double void pointer
inline void set_new_value(void** num, double a, double b){
    //conversion to particular object
    tDCplx* ptr = (tDCplx*)num; //pointer to pointer
    ptr->rl = a;
    ptr->im = b;
}

using namespace std;

int main()
{
   //local 
    tDCplx mycplx;
    mycplx.rl = 3;
    mycplx.im = 4;
    
    print_cplx(&mycplx);
    
    tDCplx* ptr   = &mycplx;
    //tDCplx** pptr = &ptr; 
    
	//test it: pointer cast to void
    void* numptr = (void*)ptr;
    void** numpptr = &numptr;
    
    set_new_value(numpptr,-2.5, 9.5);
    
    print_cplx(&mycplx,"z_new");
   
   return 0;
}

/*#include <iostream>
#include <string>


typedef struct{
 double rl;
 double im;
}tDCplx;


inline void init_cplx(tDCplx* z){
    z->rl = 0.0;
    z->im = 0.0;
}

inline void print_cplx(tDCplx* z, const std::string& str = std::string()){
    std::cout << ((str.size() == 0) ? "z" : str)<< " = "<< z->rl << " + " << z->im << "i" << std::endl;
    
}

inline void set_new_value(tDCplx** num, double a, double b){
    tDCplx* ptr = *num; //pointer to pointer
    ptr->rl = a;
    ptr->im = b;
}

using namespace std;

int main()
{
   //local 
    tDCplx mycplx;
    mycplx.rl = 3;
    mycplx.im = 4;
    
    print_cplx(&mycplx);
    
    tDCplx* ptr   = &mycplx;
    tDCplx** pptr = &ptr; 
    
    set_new_value(pptr,-2.5, 9.5);
    
    print_cplx(&mycplx,"z_new");
   
   return 0;
}
*/



//! UMFPACK stuff
//!square mtx only
//!umfpack_numeric.h
//!umf_internal.h: contains NumericType structure
some_function(...void **NumericHandle...){
NumericType *Numeric; //local obj // is this needed?
//equal to: status = UMFPACK__WARNING_singular_matrix
if((Numeric->nnzpiv < Symbolic->n_row) ||(SCALAR_IS_ZERO (Numeric->rcond)) || (SCALAR_IS_NAN (Numeric->rcond))){
	std::cout << "Singular matrix" << std::endl;
	//! try desingularization
	if(desingularize){
		//! singularity only occurs if there are zeros and/or NaN's on the diagonal of U
		for(Int i = 0; i < Symbolic->n_row; ++i){
			if(SCALAR_IS_ZERO(Numeric->D[i]) || SCALAR_IS_NAN(Numeric->D[i])){
				Numeric->D[i] = (Entry)smallPosValue;
				}
	    }	
	     Numeric->min_udiag = (double)smallPosValue;
		 Numeric->max_udiag = (double)smallPosValue;
	}	
}

//assign to function argument
*NumericHandle = (void*)Numeric; // is this needed?
}


