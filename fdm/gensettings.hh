#ifndef SOME_GENERAL_SETTINGS_FOR_USED_4_CHEMICALLY_REACTING_FLOWS_HH
#define SOME_GENERAL_SETTINGS_FOR_USED_4_CHEMICALLY_REACTING_FLOWS_HH

//no namespace here


//!===========================================================================
//!=============== TODO: comment or uncomment the following defines ==========
//!===========================================================================


//! DEBUGGING: show some output which can be used for debugging
//#define SHOW_IT_2_ME

//! use ghost points (grid assumed to include ghost points)
#define GHOST_POINTS_INCLUDED
#define REFLECTIVE_BDY_CONDITIONS
#define ONLY_CONSIDER_NEUMANN_BDY_CONDITIONS

//! do not solve for the excess species (this is set to 0 and will be recovered
//! via Y_excess = 1 - sum_{k != excess} Y_k when needed. This means that
//! the Jacobian's entries for the excess species (e.g. N2) are 0
#define SOLVE_FOR_K_MINUS_1_CONSERVATION_EQS

//! central or forward/backward FD depending on direction of flow field speed
//#define ORDER_2_CONVECTIVE_TERM

//! calculate d/d_x p(x,y,·) and d/d_y p(x,y,·) using corresponding 1st order
//! FDs, depending on the flow field direction
#define PRESSURE_GRADIENT_ORDER_1_FD

//!cf. e.g. Anderson, Chap. 10
//#define EXTRAPOLATION_BDY


//! uncommenting the following line switches correction velocity off
//#define NO_CORRECTION_DIFFUSION_VELOCITY

//! compute low mach number formulation of density change (uncomment)
//#define LOW_MACH_NUMBER_FORMULATION

//! by uncommenting the next line, ignition at the hot walls is considered.
//! Conversley, by uncommenting the next line, we consider ignition at the inlet

#define IGNITION_AT_HOT_WALLS

//#define REPAIR_ACTION //repair phys quantities INSIDE functor
//#define VALUE_PLAUSIBILITY_CHECK
//#define CHECK_4_TOO_LARGE_VALUES

//! commenting in this line uses the non-conservative form of the CFD equations,
//! i.e. if commented out you'll use the CONSERVATIVE form of the eqs.
//! this might be uncommented for special CDR equations (cf. Schiesser, etc.)
//#define NONCONSERVATIVE_FORM


//! used for non-differentiable density and/or velocity 
//#define HYDRO_PRESSURE

//#define PRESSURE_CONST   


//#define V1_CONSTANT //1st velocity component is const, i.e. d_t v1 = 0
//#define V2_CONSTANT //2nd velocity component is const, i.e. d_t v2 = 0
//#define RHO_ALGEBRAIC



//#define USE_OPENMP
#define NUMTHREADS 4

class FDMSettings{
public:


  //! false: compute $\Theta_k$ for all species, 
  //! true: only light-weight species i.e. species with m_k < 5.e-03 kg/mol, are considered, e.g. H2
  static const bool ThermalDiffusion4LightWeightSpeciesOnly = true; //false;
  
  static const char InletVelocityType = 'p';  // c/C = constant, p/P parabolic
 
  static const char WallTemperatureType = 'M'; //e/E = exponential, h/H hyperbolic tangent, f/F Frouzakis profile, m/M Marc profile, s/S natural cubic spline cut off profile, a/A profile roughly resembling an 'A'

  static const int orderOfBdyCond = 1;
  
  static const bool adjust2ndOrderNeumann = true;

  static const bool throwErrorWhenSingularMatrixDetected = false; //true = yes (error is thrown), false = no (only warning is thrown)
};

//!===========================================================================
//!===========================================================================
//!===========================================================================

#endif
