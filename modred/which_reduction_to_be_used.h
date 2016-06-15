
#define FULL_MODEL 1  // 1 = simulate full model, 0 = no

#define SWITCH_ON_REDUCTION 1  //1 = switched on, 0 = switched off

#define RED_METH_TO_APPLY 'M'  //'L' = Lebiedz, 'M' = Marc, 
                              //'F' = my very simple approach (not so good) 

//========== switch on reduction ==========================================

 //1 = yes (apply some reduction(s)), 0 = no
#define APPLY_ANY_REDUCTION 1 

//! 1 = 1st approx., 2 = 2nd approx.
#define REDUCTION_TO_BE_USED 1        

//plot both reductions at the same time: 1 = yes, 0 = no
#define PLOT_BOTH_REDUCTIONS_SIMULTANEOUSLY 0  
//===========================================================================

//'i' = implicit Euler, 
//'j' = method due to Eriksson, Johnson and Logg
//
#define SOLVER_TO_BE_USED 'i'  


#define IS_MECHANISM_ISOTHERMAL 1           //for H2-GRI


//ZELDOVICH MECHANISM
#define ZELDOVICH_NUM_RPV 2
