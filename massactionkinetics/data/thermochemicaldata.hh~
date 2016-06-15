#ifndef SOME_THERMO_CHEMICAL_DATA_HH
#define SOME_THERMO_CHEMICAL_DATA_HH

#include <iostream>
#include <string>
#include <vector> 

#include "mechutils.h"

#include "../../templatemetaprograms/compiletimearray.hh"

namespace Adonis{
 
  /**
   * \brief These are data for thermochemical reaction mechanisms. They <B>  
   * be stated in SI or SI-derived units!</B>, unless you accept badly scaled or even wrong results. These data are intrinsic to a particular combustion (<B>this should be made an exception for transport stuff!!</BB>).
   * 
   * mechanism and shouldn't be used frivolously for other combustion scenarios.
   * \tparam T the precision
   * \tparam N reference integer encoding for the mechanism
   * \code
       enum{h2c6 = 6, h2_gri = 9}
       const int k = 3;
       cout<< ThermoData4Mechanism<double,h2_gri> TD4M::thermo[9*k+2] << endl;
  *\endcode
  */
  template<class T, int N> class ThermoData4Mechanism;



    /**
   * \brief A simple ozone mechanism in 3 species due to [1]
   *
   *  References: 
   *    [1] [Maas and Pope, "Simulation of Thermal Ignition Processes in Two-Dimnsional Geometries", Zeitschrift fuer Physikalische Chemie, Neue Folge, Bd. 161, 1989, p. 72]
   */
  template<class T> 
  class ThermoData4Mechanism<T,3>{   //3 = O3
  public:
    typedef size_t index_type;
    typedef index_type* index_type_pointer;
    typedef int int_type;
    typedef int_type* int_type_pointer;
    typedef T value_type;
    typedef T* value_type_pointer; 
    //this is the order of which the species are to be considered
    enum{O,
	 O2,       
	 O3
    };

    enum{
      ENCOD=3
    };

    enum{
      nspec = 3,        //number of species involved in the mechanism 
      nreac = 3,       // total number of <--> reactions
      ntroereac = 0,    //number of TROE reactions
      ntb = 1,           //number of 3rd bodies
      nelem = 1         //number of chem. elements in mechanism
    };
    
     static inline std::string* species_names() {
      static std::string names[] = {
	"O","O2","O3"
      };
      return names;
     }

    static inline T* default_values(){
      static T dfv[] = {0.,0.8,0.2};
      return dfv;
    }

     static inline T density_of_mixture() {
       return 1.;  //default value
     }

    static inline T* thermo(){
      static T th[] ={
	//O
	0.02542059E+02, -0.02755061E-03, -0.03102803E-07, 0.04551067E-10, -0.04368051E-14,
	0.02923080E+06, 0.04920308E+02, 0.02946428E+02,-0.16381665E-02, 0.02421031E-04,    
	-0.16028431E-08, 0.03890696E-11, 0.02914764E+06, 0.02963995E+02,                   
	//O2   
	0.03697578E+02, 0.06135197E-02,-0.12588420E-06, 0.01775281E-09,-0.11364354E-14,    
	-0.12339301E+04, 0.03189165E+02, 0.03212936E+02, 0.11274864E-02,-0.05756150E-05,    
	0.13138773E-08,-0.08768554E-11,-0.10052490E+04, 0.06034737E+02,                   
	//O3    
	0.05429371E+02, 0.01820380E-01,-0.07705607E-05, 0.14992929E-09,-0.10755629E-13,    
	0.15235267E+05,-0.03266386E+02, 0.02462608E+02, 0.09582781E-01,-0.07087359E-04,    
	0.13633683E-08, 0.02969647E-11, 0.16061522E+05, 0.12141870E+02    
      };
      return th;
    }

    //set some temperature bounds
     static inline T* temperature_bounds(){
       static T tbds[] ={
	 273., 5000., 1000.,       //previously: T_low = 300 K    
	 273., 5000., 1000.,
	 273., 5000., 1000.
       };

       return tbds;
     }
    

    static inline index_type* stoichiometric_matrix(){
      static index_type sm[] = {
	2,0,0,
	0,1,0,
	0,0,1,
	1,1,0,
	1,0,1,
	0,2,0
      };
      return sm;
    }

    //! the values for the forward reactions
    //! NOTE: these are already in SI format
    static inline T* arrhenius_values(){
      static T arrh[] ={
	2.9e+05, -1., 0.,
	9.5e+08, 0., 9.5e+04,
	5.2e+06, 0., 1.74e+04
      };
      return arrh;
    }


    
     static inline int_type_pointer troewtb(){
      static int tw[] = {
	//troe info
	0,
	0,
	0,
	//3rd body info
	0,              //M(0)
	0,              //M(0) 
	-1    //third bidirectional reaction doesn't have 3rd body
      };
      return tw;
     }

     //! number of third bodies x number of species
    static inline T* collision_efficiencies(){
      static T ce[] ={  //no collision coefficients given
	//M(0)
	1.14,    //alpha_O
	0.40,    //alpha_O2
	0.92     //alpha_O3
      };
      return ce;
    }
  
    //molecular weights
    //! from periodic table -- here in SI units!
    inline static T* molar_masses(){
      static T mm[] ={
	//!"O","O2","O3"
	1.600000e-02, 3.200000e-02, 4.800000e-02
      };
      return mm;
    }


    inline static T* transport(){
      static T tpt[] ={
 //Geometry    epsilon/kappa_B   sigma      mu               alpha      Z_rot
//Units:--           K            A         Debye             A³         --    
//O   
	0,        80.0,        2.75,        0.,                0.,        0.,
//O2	
	1,        107.4,       3.458,       0.,               1.6e-30,    3.8,
//O3
	2,        180.0,       4.1,         0.,               0.,        2.0
      };
      return tpt;
    }

    
     //! # of bidirec reactions
    //! when \rev is given, compute reverse reaction via Arrhenius law
    // ! 1 = yes, 0 = no (then compute it like always)
     static inline bool* explicit_reverse_reaction_coefficients(){
      static bool revcoef[] = {
	0,
	0,
	0
      };
      return revcoef;
     }


    
    //! tells us if species k participates in any (at least one) reaction
    //! true if it does, false if not
    static inline bool* is_species_reactive(){
      static bool rtive[] = {
	true,   //O
	true,   //O2
	true    //O3
      };
      return rtive;
     }

    //====================== TODO: The following  might be changed  ===========
    //=====================        from time to time                 ===========
    enum{
      rednspec = 1,  //reduced dimension specifications for RPV O
      rednreac = 3   //
    };
    
    static inline index_type* rpv_index() {
      static index_type rpv[] = {0};   //O
      
      return rpv;
    }
  
    static inline index_type* rpv_reaction_index() {
      static index_type rri[] = {0,1,2};  //when O serves as RPV (all reactions)
      
      return rri;
    }
    

    // //!compile time version of <TT> rpv_index()</TT>. Looks quite bizarre, but
    // //!we get compile-time performance
    typedef Element<0> RpvIndexType;
    
    // //!compile time version of <TT> rpv_reaction_index() </TT>
    typedef Element<0,Element<1,Element<2> > > RpvReactionIndexType;

    
    inline static T* mass_balance_matrix(){
      static T mbm[] ={1.,1.,1.};   //! sum of mass fractions equals one
      return mbm;
    }

     inline static T* mass_sum(){
      static T ms[] ={1.};
      return ms;
     }

  }; //end of class specialization





  /**
   * \brief a ficticious H2 combustion mechanism in 6 species
   */
  template<class T> 
  class ThermoData4Mechanism<T,6>{ 
  public:
    typedef size_t index_type;
    typedef index_type* index_type_pointer;
    typedef int int_type;
    typedef int_type* int_type_pointer;
    typedef T value_type;
    typedef T* value_type_pointer; 
    //this is the order of which the species are to be considered
    enum{H2,
	 H,       
	 O2,      
	 O,       
	 H2O,     
	 OH};

    enum{
      ENCOD=6
    };

    enum{
      nspec = 6,        //number of species involved in the mechanism 
      nreac = 6,       // total number of <--> reactions
      ntroereac = 0,    //number of TROE reactions
      ntb = 0,           //number of 3rd bodies
      nelem = 2        //number of chem. elements in mechanism
    };

     static inline std::string* species_names() {
      static std::string names[] = {
	"H2","H","O2","O","H2O","OH"
      };
      return names;
     }

     static inline T density_of_mixture() {
       return 1.;
     }

    static inline T* thermo(){
      static T th[] ={
	//H2
	0.02991423E+02, 0.07000644E-02,-0.05633828E-06,-0.09231578E-10, 0.15827519E-14,    
	-0.08350340E+04,-0.13551101E+01, 0.03298124E+02, 0.08249441E-02,-0.08143015E-05,    
	-0.09475434E-09, 0.04134872E-11,-0.10125209E+04,-0.03294094E+02,
	//H
	0.02500000E+02, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,    
	0.02547162E+06,-0.04601176E+01, 0.02500000E+02, 0.00000000E+00, 0.00000000E+00,    
	0.00000000E+00, 0.00000000E+00, 0.02547162E+06,-0.04601176E+01,
	//O2
	0.03697578E+02, 0.06135197E-02,-0.12588420E-06, 0.01775281E-09,-0.11364354E-14,    
	-0.12339301E+04, 0.03189165E+02, 0.03212936E+02, 0.11274864E-02,-0.05756150E-05,    
	0.13138773E-08,-0.08768554E-11,-0.10052490E+04, 0.06034737E+02,
	//O
	0.02542059E+02,-0.02755061E-03,-0.03102803E-07, 0.04551067E-10,-0.04368051E-14,    
	0.02923080E+06, 0.04920308E+02, 0.02946428E+02,-0.16381665E-02, 0.02421031E-04,    
	-0.16028431E-08, 0.03890696E-11, 0.02914764E+06, 0.02963995E+02,
	//H2O
	0.02672145E+02, 0.03056293E-01,-0.08730260E-05, 0.12009964E-09,-0.06391618E-13,    
	-0.02989921E+06, 0.06862817E+02, 0.03386842E+02, 0.03474982E-01,-0.06354696E-04,    
	0.06968581E-07,-0.02506588E-10,-0.03020811E+06, 0.02590232E+02,
	//OH
	0.02882730E+02, 0.10139743E-02,-0.02276877E-05, 0.02174683E-09,-0.05126305E-14,    
	0.03886888E+05, 0.05595712E+02, 0.03637266E+02, 0.01850910E-02,-0.16761646E-05,    
	0.02387202E-07,-0.08431442E-11, 0.03606781E+05, 0.13588605E+01  
      };
      
      return th;
    }

    //set some temperature bounds
     static inline T* temperature_bounds(){
       static T tbds[] ={
	 300., 5000., 1000.,
	 300., 5000., 1000.,
	 300., 5000., 1000.,
	 300., 5000., 1000.,
	 300., 5000., 1000.,
	 300., 5000., 1000.
       };

       return tbds;
     }
    

    static inline index_type* stoichiometric_matrix(){
      static index_type sm[] = {
	1,0,0,0,0,0,
	0,2,0,0,0,0,
	0,0,1,0,0,0,
	0,0,0,2,0,0,
	0,0,0,0,1,0,
	0,1,0,0,0,1,
	1,0,0,1,0,0,
	0,1,0,0,0,1,
	0,1,1,0,0,0,
	0,0,0,1,0,1,
	1,0,0,1,0,0,
	0,0,0,0,1,0
      };
      return sm;
    }

    //! the values for the forward reactions
    static inline T* arrhenius_values(){
      static T arrh[] ={
	2., 0, 0,
	1., 0, 0,
	1., 0, 0,
	1000., 0, 0,
	1000., 0, 0,
	100., 0, 0,
      };
      return arrh;
    }


    
     static inline int_type_pointer troewtb(){
      static int tw[] = {
	//troe info
	0,
	0,
	0,
	0,
	0,
	0,
	//wtb info
	-1,
	-1,
	-1,
	-1,
	-1,
	-1,
      };
      return tw;
     }

    static inline T* collision_efficiencies(){
      static T ce[1];  //no collision coefficients given
      return ce;
    }
  
    //molecular weights
    //! from periodic table -- here in SI units!
    inline static T* molar_masses(){
      static T mm[] ={
	//!"H2","H","O2","O","H2O","OH"
	2.015800e-03, 1.007900e-03, 3.200000e-02, 1.600000e-02, 1.801580e-02, 1.700790e-02
	//2.0158, 1.0079, 32., 16., 18.0158, 17.0079
      };
      return mm;
    }


    inline static T* transport(){
      static T tpt[] ={
 //Geometry    epsilon/kappa_B   sigma      mu               alpha      Z_rot
//Units:--           K            A         Debye             A³         --    
//H2   
1.00000e+00,   3.80000e+01,   2.92000,   0.00000e+00,   7.90000e-31,   2.80000e+02,  
//H   
0.00000e+00,   1.45000e+02,   2.05000,   0.00000e+00,   0.00000e+00,   0.00000e+00, 
//O2
1.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   1.60000e-30,   3.80000e+00,
//O
0.00000e+00,   8.00000e+01,   2.75000,   0.00000e+00,   0.00000e+00,   0.00000e+00,   
//H2O
2.00000e+00,   5.72400e+02,   2.60500,   1.84400e-28,   0.00000e+00,   4.00000e+00,

//OH  
1.00000e+00,   8.00000e+01,   2.75000,   0.00000e+00,   0.00000e+00,   0.00000e+00
 
      };
      return tpt;
    }

 //! # of bidirec reactions
    //! when \rev is given, compute reverse reaction via Arrhenius law
    // ! 1 = yes, 0 = no (then compute it like always)
     static inline bool* explicit_reverse_reaction_coefficients(){
      static bool revcoef[] = {
	0,
	0,
	0,
	0,
	0,
	0
      };
      return revcoef;
     }


     //! tells us if species k participates in any (at least one) reaction
    //! true if it does, false if not
    static inline bool* is_species_reactive(){
      static bool rtive[] = {
	true,   
	true,   
	true,
	true,
	true,
	true
      };
      return rtive;
     }

    //====================== TODO: The following  might be changed  ===========
    //=====================        from time to time                 ===========
    enum{
      rednspec = 2,  //reduced dimension
      rednreac = 4
    };
    
    static inline index_type* rpv_index() {
      static index_type rpv[] = {0,4};
      
      return rpv;
    }
  
    //! README:
    //! take each bidirectional reaction (here we have 6 of such directions)
    //! numbering consecutively, beginning with 0 (here we have 0,...,5)
    //! then mark each direction in which H2 and/or H2O (the two rpvs) occur
    //! (H2 in 0, 3, 5  and H2O in 2 and 5, therefore, form the union of
    //! both sets, i.e. [0,2,3,5]
    static inline index_type* rpv_reaction_index() {
      static index_type rri[] = {0,2,3,5};
      
      return rri;
    }
    

    //!compile time version of <TT> rpv_index()</TT>. Looks quite bizarre, but
    //!we get compile-time performance
    typedef Element<0,Element<4> > RpvIndexType;
    
    //!compile time version of <TT> rpv_reaction_index() </TT>
    typedef Element<0,Element<2,Element<3,Element<5> > > > RpvReactionIndexType;

    
    inline static T* mass_balance_matrix(){
      static T mbm[] ={2., 1., 0., 0., 2., 1.,
		       0., 0., 2., 1., 1., 1.};
      return mbm;
    }

     inline static T* mass_sum(){
      static T ms[] ={2., 1.};
      return ms;
    }

  }; //end of class specialization







  /**
   * \brief Hydrogen combustion part of the GRI-Mech 3.0 mechanism
   *
   * References:
   * [1] <a href="http://www.me.berkeley.edu/gri_mech/"> GRI-Mech 3.0 </a>
   *
   * The full mechanism and its thermo data for the NASA polynoms can be found at 
   <a href="http://me.berkeley.edu/gri_mech/version30/files30/grimech30.dat"> Complete mechanism, including Arrhenius data </a>
   *
   * <a href="http://me.berkeley.edu/gri_mech/version30/files30/thermo30.dat"> thermodynamics in NASA polynomial format for CHEMKIN-II>
   *
   * The transport data are given at 
   *
   * <a href="http://me.berkeley.edu/gri_mech/version30/files30/transport.dat"> Transport specifications </a>
   * 
   * Note: each species has 6 numbers which represent the following quantities:
   * geometry, Lennard Jones potential well depth [K], Lennard Jones collision diameter [\f$10^{-10}\f$m], dipole moment [\f$ 10^{-28}\f$kg\f$\cdot\f$m\f$^{\frac{1}{2}}/s \f$]m, polarizability [10\f$^{-30}\f$m] and rotational collision number at 298K  
   *
   *
   * This mechanism can be found in [1] with corresponding data
   *
   * Ref.:
   * [1]  <a href="http://www.princeton.edu/~lam/documents/LSTEP/CHEMKIN/Hydrogen"> H2 Gri </a>
   */
  template<class T> 
  class ThermoData4Mechanism<T,9>{ 
  public:
    typedef size_t index_type;
    typedef index_type* index_type_pointer;
    typedef int int_type;
    typedef int_type* int_type_pointer;
    typedef T value_type;
    typedef T* value_type_pointer; 
    //this is the order of which the species are to be considered
    enum{O, O2, H, OH, H2, HO2, H2O2, H2O, N2};


    enum{
      ENCOD=9
    };

    enum{
      nspec = 9,
      nreac = 28,
      ntroereac = 1,
      ntb = 6,           //number of 3rd bodies
      nelem = 3         //number of chem. elements in mechanism
    };

    static inline std::string* species_names() {
      static std::string names[] = {
	"O","O2","H","OH", "H2","HO2","H2O2","H2O","N2"
      };
      return names;
    }

    
    static inline T density_of_mixture() {
      return 0.3;    //in kg/m³
    }

    static inline T* thermo(){
      static T th[] = {
	//O
	2.56942078E+00,-8.59741137E-05, 4.19484589E-08,-1.00177799E-11, 1.22833691E-15, 2.92175791E+04, 4.78433864E+00, 3.16826710E+00,-3.27931884E-03, 6.64306396E-06,-6.12806624E-09, 2.11265971E-12, 2.91222592E+04, 2.05193346E+0,
  //O2
			3.28253784E+00,1.48308754E-03,-7.57966669E-07, 2.09470555E-10,-2.16717794E-14, -1.08845772E+03, 5.45323129E+00, 3.78245636E+00,-2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12,-1.06394356E+03, 3.65767573E+00,

//H
			2.50000001E+00,-2.30842973E-11, 1.61561948E-14,-4.73515235E-18, 4.98197357E-22,    
			2.54736599E+04,-4.46682914E-01, 2.50000000E+00, 7.05332819E-13,-1.99591964E-15,    
			2.30081632E-18,-9.27732332E-22, 2.54736599E+04,-4.46682853E-01,

//OH
3.09288767E+00, 5.48429716E-04, 1.26505228E-07,-8.79461556E-11, 1.17412376E-14,    
			 3.85865700E+03, 4.47669610E+00, 3.99201543E+00,-2.40131752E-03, 4.61793841E-06,    
			-3.88113333E-09, 1.36411470E-12, 3.61508056E+03,-1.03925458E-01,
 
//H2	
3.33727920E+00,-4.94024731E-05, 4.99456778E-07,-1.79566394E-10, 2.00255376E-14,   
			  -9.50158922E+02,-3.20502331E+00, 2.34433112E+00, 7.98052075E-03,-1.94781510E-05,    
			2.01572094E-08,-7.37611761E-12,-9.17935173E+02, 6.83010238E-01,

//HO2
    4.01721090E+00, 2.23982013E-03,-6.33658150E-07, 1.14246370E-10,-1.07908535E-14,    
		      1.11856713E+02, 3.78510215E+00, 4.30179801E+00,-4.74912051E-03, 2.11582891E-05,    
			-2.42763894E-08, 9.29225124E-12, 2.94808040E+02, 3.71666245E+00,

//H2O2
4.16500285E+00, 4.90831694E-03,-1.90139225E-06, 3.71185986E-10,-2.87908305E-14,    
			    -1.78617877E+04, 2.91615662E+00, 4.27611269E+00,-5.42822417E-04, 1.67335701E-05,    
			-2.15770813E-08, 8.62454363E-12,-1.77025821E+04, 3.43505074E+00,

 //H2O
 3.03399249E+00, 2.17691804E-03,-1.64072518E-07,-9.70419870E-11, 1.68200992E-14,    
			   -3.00042971E+04, 4.96677010E+00, 4.19864056E+00,-2.03643410E-03, 6.52040211E-06,    
			-5.48797062E-09, 1.77197817E-12,-3.02937267E+04,-8.49032208E-01,

 //N2
			0.02926640E+02, 0.14879768E-02,-0.05684760E-05, 0.10097038E-09,-0.06753351E-13,    
			 -0.09227977E+04, 0.05980528E+02, 0.03298677E+02, 0.14082404E-02,-0.03963222E-04,    
			 0.05641515E-07,-0.02444854E-10,-0.10208999E+04, 0.03950372E+02   };


      return th;
    }

    

    static inline T* temperature_bounds(){
      static T tbds[] = {273,3500,1000,  //O   //previously 
			 273,3500,1000,  //O2
			 273,3500,1000,  //H
			 273,3500,1000,  //OH
			 273,3500,1000,  //H2
			 273,3500,1000,  //HO2
			 273,3500,1000,  //H2O2
			 273,3500,1000,  //H2O
			 273,5000,1000   //N2
      };
    
      return tbds;
    }


    static inline index_type* stoichiometric_matrix(){
      //store only ni' 
    //each column corresponds to O, O2, H, OH,  H2,  HO2,  H2O2,  H2O and  N2,
      //respectively. Each <--> reaction is split into a --> and a <-- reaction
      static index_type sm[] = {2,0,0,0,0,0,0,0,0,
		   0,1,0,0,0,0,0,0,0,
		   1,0,1,0,0,0,0,0,0,
		   0,0,0,1,0,0,0,0,0,
		   1,0,0,0,1,0,0,0,0,
		   0,0,1,1,0,0,0,0,0,
		   1,0,0,0,0,1,0,0,0,
		   0,1,0,1,0,0,0,0,0,
		   1,0,0,0,0,0,1,0,0,
		   0,0,0,1,0,1,0,0,0,
		   0,1,1,0,0,0,0,0,0,
		   0,0,0,0,0,1,0,0,0,
		   0,2,1,0,0,0,0,0,0,
		   0,1,0,0,0,1,0,0,0,
		   0,1,1,0,0,0,0,1,0,
		   0,0,0,0,0,1,0,1,0,
		   0,1,1,0,0,0,0,0,1,
		   0,0,0,0,0,1,0,0,1,
		   0,1,1,0,0,0,0,0,0,
		   1,0,0,1,0,0,0,0,0,
		   0,0,2,0,0,0,0,0,0,
		   0,0,0,0,1,0,0,0,0,
		   0,0,2,0,1,0,0,0,0,
		   0,0,0,0,2,0,0,0,0,
		   0,0,2,0,0,0,0,1,0,
		   0,0,0,0,1,0,0,1,0,
		   0,0,1,1,0,0,0,0,0,
		   0,0,0,0,0,0,0,1,0,
		   0,0,1,0,0,1,0,0,0,
		   1,0,0,0,0,0,0,1,0,
		   0,0,1,0,0,1,0,0,0,
		   0,1,0,0,1,0,0,0,0,
		   0,0,1,0,0,1,0,0,0,
		   0,0,0,2,0,0,0,0,0,
		   0,0,1,0,0,0,1,0,0,
		   0,0,0,0,1,1,0,0,0,
		   0,0,1,0,0,0,1,0,0,
		   0,0,0,1,0,0,0,1,0,
		   0,0,0,1,1,0,0,0,0,
		   0,0,1,0,0,0,0,1,0,
		   0,0,0,2,0,0,0,0,0,
		   0,0,0,0,0,0,1,0,0,
		   0,0,0,2,0,0,0,0,0,
		   1,0,0,0,0,0,0,1,0,
		   0,0,0,1,0,1,0,0,0,
		   0,1,0,0,0,0,0,1,0,
		   0,0,0,1,0,0,1,0,0,
		   0,0,0,0,0,1,0,1,0,
		   0,0,0,1,0,0,1,0,0,
		   0,0,0,0,0,1,0,1,0,
		   0,0,0,0,0,2,0,0,0,
		   0,1,0,0,0,0,1,0,0,
		   0,0,0,0,0,2,0,0,0,
		   0,1,0,0,0,0,1,0,0,
		   0,0,0,1,0,1,0,0,0,
		   0,1,0,0,0,0,0,1,0 
    };

      return sm;
    }

    
    //! NOTE: these are already in SI format, i.e. they were transformed 
    //! accordingly once before
    static inline T* arrhenius_values(){
      static T arrh[] ={
1.200000e+05,   -1.000000e+00,   0.000000e+00,   
5.000000e+05,   -1.000000e+00,   0.000000e+00,   
3.870000e-02,   2.700000e+00,   2.619180e+04,   
2.000000e+07,   0.000000e+00,   0.000000e+00,   
9.630000e+00,   2.000000e+00,   1.673600e+04,   
2.800000e+06,   -8.600000e-01,   0.000000e+00,   
2.080000e+07,   -1.240000e+00,   0.000000e+00,   
1.126000e+07,   -7.600000e-01,   0.000000e+00,   
2.600000e+07,   -1.240000e+00,   0.000000e+00,   
2.650000e+10,   -6.707000e-01,   7.129950e+04,   
1.000000e+06,   -1.000000e+00,   0.000000e+00,   
9.000000e+04,   -6.000000e-01,   0.000000e+00,   
6.000000e+07,   -1.250000e+00,   0.000000e+00,   
2.200000e+10,   -2.000000e+00,   0.000000e+00,   
3.970000e+06,   0.000000e+00,   2.807460e+03,   
4.480000e+07,   0.000000e+00,   4.469510e+03,   
8.400000e+07,   0.000000e+00,   2.656840e+03,   
1.210000e+01,   2.000000e+00,   2.175680e+04,   
1.000000e+07,   0.000000e+00,   1.506240e+04,   
2.160000e+02,   1.510000e+00,   1.435110e+04,   
	//TROE reaction at index 20 with molecularity 3
7.400000e+07,   -3.700000e-01,   0.000000e+00,   
3.570000e-02,   2.400000e+00,   -8.828200e+03,   
1.450000e+07,   0.000000e+00,   -2.029000e+03,   
2.000000e+06,   0.000000e+00,   1.786570e+03,   
1.700000e+12,   0.000000e+00,   1.230510e+05,   
1.300000e+05,   0.000000e+00,   -6.819900e+03,   
4.200000e+08,   0.000000e+00,   5.020800e+04,   
 5.000000e+09,   0.000000e+00,   7.250870e+04, 
     
	//!starting from 3*nreac, the troe coefficients are stored in 7-entries
	//!blocks -- each containing A_low, beta_low, Ea_low, alpha,T***,T* and 
	//! T** in this order
//A_low            beta_low          Ea_low
2.300000e+06,   -9.000000e-01,   -7.112800e+03,
//alpha   T***   T*       T**
0.7346,   94.0,  1756.0,  21.68
	
 };

      return arrh;
    }

    
    

    //! troe with third bodies information
    static inline int_type_pointer troewtb(){
      static int tw[] = {0,
			0,
			 0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			0,
			1,
			0,
			0,
			0,
			0,
			0,
			0,
			0, 
//now comes information about 3rd bodies:
			0,
			1,
			-1,
			-1,
			-1,
			2,
			-1,
			-1,
			-1,
			-1,
			3,
			-1,
			-1,
			4,
			-1,
			-1,
			-1,
			-1,
			-1,
			-1,
			5,
			-1,
			-1,
			-1,
			-1,
			-1,
			-1,
			-1
      };

      return tw;
    }

  
    //! 6 THIRD BODIES x 9 SPECIES = 54 entries
    static inline T* collision_efficiencies(){
      
      static T ce[] = {
	//M(0)
	1,
	1,
	1,
	1,
	2.4,                  //H2
	1,
	1,
	15.4,                 //H2O
	1,
	//M(1)
	1,
	1,
	1,
	1,
	2,
	1,
	1,
	6,
	1,
	//M(2)
	1,
	0,
	1,
	1,
	1,
	1,
	1,
	0,
	1,
	//M(3)
	1,
	1,
	1,
	1,
	0,
	1,
	1,
	0,
	1,
	//M(4)
	1,
	1,
	1,
	1,
	0.73,
	1,
	1,
	3.65,
	1,
	//M(5)
	1,
	1,
	1,
	1,
	2,
	1,
	1,
	6,
	1

      };
      return ce;
    }
    

    //! note that the molar masses (also know as molecular weights) 
    //! which can be computed from the periodic table are stated in g/mol. 
    //! Here we use SI units, i.e. kg/mol. Therefore, the
    //! molar masses \f$M_k, \ k=1, \ldots, K\f$ must be devided by 1000, after
    //! having been calculated from the periodic table, as done below.
    inline static T* molar_masses(){
      //O, O2, H, OH, H2, HO2, H2O2, H2O, N2
      static T mm[] ={1.600000e-02, 3.200000e-02, 1.007900e-03, 1.700790e-02, 2.015800e-03, 3.300790e-02, 3.401580e-02, 1.801580e-02, 2.802000e-02 };
      
      return mm;
    }


    inline static T* transport(){
      static T tpt[] ={
 //Geometry    epsilon/kappa_B      sigma      mu               alpha      Z_rot
 //O
0.00000e+00,   8.00000e+01,   2.75000,   0.00000e+00,   0.00000e+00,   0.00000e+00,   
//O2
1.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   1.60000e-30,   3.80000e+00,
//H   
0.00000e+00,   1.45000e+02,   2.05000,   0.00000e+00,   0.00000e+00,   0.00000e+00, 
//OH  
1.00000e+00,   8.00000e+01,   2.75000,   0.00000e+00,   0.00000e+00,   0.00000e+00,
//H2   
1.00000e+00,   3.80000e+01,   2.92000,   0.00000e+00,   7.90000e-31,   2.80000e+02,   
//HO2
2.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   0.00000e+00,   1.00000e+00,   
//H2O2
2.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   0.00000e+00,   3.80000e+00,   
//H2O
2.00000e+00,   5.72400e+02,   2.60500,   1.84400e-28,   0.00000e+00,   4.00000e+00,
//N2   
1.00000e+00,   9.75300e+01,   3.62100,   0.00000e+00,   1.76000e-30,   4.00000e+00 
      };

      return tpt;
    }

    //! # of bidirec reactions
    //! when \rev is given, compute reverse reaction via Arrhenius law
    // ! 1 = yes, 0 = no (then compute it like always)
     static inline bool* explicit_reverse_reaction_coefficients(){
      static bool revcoef[] = {
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0
      };
      return revcoef;
     }

    //! tells us if species k participates in any (at least one) reaction
    //! true if it does, false if not
    static inline bool* is_species_reactive(){
      static bool rtive[] = {
	true,   
	true,   
	true,
	true,
	true,
	true,
	true,
	true,
	true
      };
      return rtive;
     }

     //====================== TODO: The following  might be changed  ===========
    //=====================        from time to time                 ===========
    
    /*
    //! 3 rpvs involved
    enum{
      rednspec = 3,  //reduced dimension
      rednreac = 24
    };
    
    static inline index_type* rpv_index() {
      static index_type rpv[] = {1,4,7};
      
      return rpv;
    }
  
    static inline index_type* rpv_reaction_index() {
      static index_type rri[] = {0,2,3,5,6,7,8,9,10,11,12,13,14,15,17,18,19,21,22,23,24,25,26,27};
      
      return rri;
    }

    //! actually not needed since the construction can be obtained via
    //! species_names() and rpv_index()
    */

    
      //! 2 rpvs involved -- here you loose more information
    enum{
      rednspec = 2,  //reduced dimension
      rednreac = 16
    };
    
    static inline index_type* rpv_index() {
      static index_type rpv[] = {4,7};
      
      return rpv;
    }
  
    static inline index_type* rpv_reaction_index() {
      static index_type rri[] = {2,7,10,11,12,13,14,15,17,18,19,21,22,23,24,27};
      
      return rri;
    }

     
    //!compile-time versions
    typedef Element<4,Element<7> > RpvIndexType;
     
    typedef Element<2,Element<7,Element<10,Element<11, Element<12, Element<13, Element<14, Element<15, Element<17, Element<18, Element<19, Element<21, Element<22, Element<23, Element<24, Element<27> > > > > > > > > > > > > > > > RpvReactionIndexType;

  };  //end of class 



  
  /**
   *
   * \brief CHEMKIN H2 default mechanism.
   *
   * The mechanism along with thermochemical data can be found on the web [1]
   *
   * Ref.:
   * [1]  <a href="http://www.engr.colostate.edu/~marchese/combustion11/chemkin.html"> CHEMKIN H2 default mechanism </a>
   */
  template<class T> 
  class ThermoData4Mechanism<T,91>{  //mechanism 9.1
  public:
    typedef size_t index_type;
    typedef index_type* index_type_pointer;
    typedef int int_type;
    typedef int_type* int_type_pointer;
    typedef T value_type;
    typedef T* value_type_pointer; 
    //this is the order of which the species are to be considered
    enum{O, O2, H, OH, H2, HO2, H2O2, H2O, N2};

    enum{
      ENCOD=91
    };

    enum{
      nspec = 9,
      nreac = 20,
      ntroereac = 0,
      ntb = 6,           //number of 3rd bodies
      nelem = 3         //number of chem. elements in mechanism
    };

    static inline std::string* species_names() {
      static std::string names[] = {
	"O","O2","H","OH", "H2","HO2","H2O2","H2O","N2"
      };
      return names;
    }

    
    static inline T density_of_mixture() {
      return 0.3;    //in kg/m³
    }

    static inline T* thermo(){
      static T th[] = {
	//O
	2.54206000e+00,-2.75506200e-05,-3.10280300e-09, 4.55106700e-12,-4.36805200e-16,    
	2.92308000e+04, 4.92030800e+00, 2.94642900e+00,-1.63816600e-03, 2.42103200e-06,    
	-1.60284300e-09, 3.89069600e-13, 2.91476400e+04, 2.96399500e+00,
	// O2
	3.69757800e+00, 6.13519700e-04,-1.25884200e-07, 1.77528100e-11,-1.13643500e-15,    
	-1.23393000e+03, 3.18916600e+00, 3.21293600e+00, 1.12748600e-03,-5.75615000e-07,    
	1.31387700e-09,-8.76855400e-13,-1.00524900e+03, 6.03473800e+00,  
	//H
	2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,    
	2.54716300e+04,-4.60117600e-01, 2.50000000e+00, 0.00000000e+00, 0.00000000e+00,    
	0.00000000e+00, 0.00000000e+00, 2.54716300e+04,-4.60117600e-01,
	// OH 
	2.88273000e+00, 1.01397400e-03,-2.27687700e-07, 2.17468400e-11,-5.12630500e-16,    
	3.88688800e+03, 5.59571200e+00, 3.63726600e+00, 1.85091000e-04,-1.67616500e-06,    
	2.38720300e-09,-8.43144200e-13, 3.60678200e+03, 1.35886000e+00,
	//H2 
	2.99142300e+00, 7.00064400e-04,-5.63382900e-08,-9.23157800e-12, 1.58275200e-15,    
	-8.35034000e+02,-1.35511000e+00, 3.29812400e+00, 8.24944200e-04,-8.14301500e-07,    
	-9.47543400e-11, 4.13487200e-13,-1.01252100e+03,-3.29409400e+00, 
	// HO2 
	4.07219100e+00, 2.13129600e-03,-5.30814500e-07, 6.11226900e-11,-2.84116500e-15,    
	-1.57972700e+02, 3.47602900e+00, 2.97996300e+00, 4.99669700e-03,-3.79099700e-06,    
	2.35419200e-09,-8.08902400e-13, 1.76227400e+02, 9.22272400e+00, 
	// H2O2
	4.57316700e+00, 4.33613600e-03,-1.47468900e-06, 2.34890400e-10,-1.43165400e-14,    
	-1.80069600e+04, 5.01137000e-01, 3.38875400e+00, 6.56922600e-03,-1.48501300e-07,    
	-4.62580600e-09, 2.47151500e-12,-1.76631500e+04, 6.78536300e+00, 
	//H2O
	2.67214600e+00, 3.05629300e-03,-8.73026000e-07, 1.20099600e-10,-6.39161800e-15,    
	-2.98992100e+04, 6.86281700e+00, 3.38684200e+00, 3.47498200e-03,-6.35469600e-06,    
	6.96858100e-09,-2.50658800e-12,-3.02081100e+04, 2.59023300e+00,
	//N2
	2.92664000e+00, 1.48797700e-03,-5.68476100e-07, 1.00970400e-10,-6.75335100e-15,    
	-9.22797700e+02, 5.98052800e+00, 3.29867700e+00, 1.40824000e-03,-3.96322200e-06,    
	5.64151500e-09,-2.44485500e-12,-1.02090000e+03, 3.95037200e+00
 };


      return th;
    }

    

    static inline T* temperature_bounds(){
      static T tbds[] = {273,5000,1000,  //O   //previously 
			 273,5000,1000,  //O2
			 273,5000,1000,  //H
			 273,5000,1000,  //OH
			 273,5000,1000,  //H2
			 273,5000,1000,  //HO2
			 273,5000,1000,  //H2O2
			 273,5000,1000,  //H2O
			 273,5000,1000   //N2
      };
    
      return tbds;
    }


    static inline index_type* stoichiometric_matrix(){
      //store only ni' 
    //each column corresponds to O, O2, H, OH,  H2,  HO2,  H2O2,  H2O and  N2,
      //respectively. Each <--> reaction is split into a --> and a <-- reaction
      static index_type sm[] = {
	0,  1,  1,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  1,  0,  0,  0,  
	0,  0,  2,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  1,  0,  0,  0,  0,  
	0,  0,  2,  0,  1,  0,  0,  0,  0,  
	0,  0,  0,  0,  2,  0,  0,  0,  0,  
	0,  0,  2,  0,  0,  0,  0,  1,  0,  
	0,  0,  0,  0,  1,  0,  0,  1,  0,  
	0,  0,  1,  1,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  1,  0,  
	1,  0,  1,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  1,  0,  0,  0,  0,  0,  
	2,  0,  0,  0,  0,  0,  0,  0,  0,  
	0,  1,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  1,  0,  0,  
	0,  0,  0,  2,  0,  0,  0,  0,  0,  
	0,  1,  0,  0,  1,  0,  0,  0,  0,  
	0,  0,  0,  2,  0,  0,  0,  0,  0,  
	0,  0,  0,  1,  1,  0,  0,  0,  0,  
	0,  0,  1,  0,  0,  0,  0,  1,  0,  
	1,  0,  0,  1,  0,  0,  0,  0,  0,  
	0,  1,  1,  0,  0,  0,  0,  0,  0,  
	1,  0,  0,  0,  1,  0,  0,  0,  0,  
	0,  0,  1,  1,  0,  0,  0,  0,  0,  
	0,  0,  0,  1,  0,  1,  0,  0,  0,  
	0,  1,  0,  0,  0,  0,  0,  1,  0,  
	0,  0,  1,  0,  0,  1,  0,  0,  0,  
	0,  0,  0,  2,  0,  0,  0,  0,  0,  
	1,  0,  0,  0,  0,  1,  0,  0,  0,  
	0,  1,  0,  1,  0,  0,  0,  0,  0,  
	0,  0,  0,  2,  0,  0,  0,  0,  0,  
	1,  0,  0,  0,  0,  0,  0,  1,  0,  
	0,  0,  1,  0,  0,  1,  0,  0,  0,  
	0,  1,  0,  0,  1,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  2,  0,  0,  0,  
	0,  1,  0,  0,  0,  0,  1,  0,  0,  
	0,  0,  1,  0,  0,  0,  1,  0,  0,  
	0,  0,  0,  0,  1,  1,  0,  0,  0,  
	0,  0,  0,  1,  0,  0,  1,  0,  0,  
	0,  0,  0,  0,  0,  1,  0,  1,  0
    };

      return sm;
    }

    
    //! NOTE: these are already in SI format, i.e. they were transformed 
    //! accordingly once before
    //!      A              beta             Ea
    static inline T* arrhenius_values(){
      static T arrh[] ={
	3.610000e+05,   -7.200000e-01,   0.000000e+00,
	1.000000e+06,   -1.000000e+00,   0.000000e+00,
	9.200000e+04,   -6.000000e-01,   0.000000e+00,
	6.000000e+07,   -1.250000e+00,   0.000000e+00,
	1.600000e+10,   -2.000000e+00,   0.000000e+00,
	6.200000e+04,   -6.000000e-01,   0.000000e+00,
	1.890000e+01,   0.000000e+00,   -7.480992e+03,
	1.300000e+11,   0.000000e+00,   1.903720e+05,
	1.700000e+07,   0.000000e+00,   1.999115e+05,
	1.170000e+03,   1.300000e+00,   1.517118e+04,
	3.610000e+08,   -5.000000e-01,   0.000000e+00,
	5.060000e-02,   2.670000e+00,   2.631736e+04,
	7.500000e+06,   0.000000e+00,   0.000000e+00,
	1.400000e+08,   0.000000e+00,   4.489432e+03,
	1.400000e+07,   0.000000e+00,   4.489432e+03,
	6.000000e+02,   1.300000e+00,   0.000000e+00,
	1.250000e+07,   0.000000e+00,   0.000000e+00,
	2.000000e+06,   0.000000e+00,   0.000000e+00,
	1.600000e+06,   0.000000e+00,   1.589920e+04,
	1.000000e+07,   0.000000e+00,   7.531200e+03
	
 };

      return arrh;
    }

    
    

    //! troe with third bodies information
    static inline int_type_pointer troewtb(){
      static int tw[] = {  //no fall-off reaction at all
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
//now comes information about 3rd bodies:
	0,
	1,
	-1,
	-1,
	2,
	3,
	4,
	5,
	-1,
	-1,
	-1,
	-1
	-1,
	-1,
	-1,
	-1,
	-1,
	-1,
	-1,
	-1
      };

      return tw;
    }

  
    //! 6 THIRD BODIES x 9 SPECIES = 54 entries
    static inline T* collision_efficiencies(){
      
      static T ce[] = {
	//M(0)
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	2.86,  //H2
	1,  //HO2
	1,  //H2O2
	18.6,  //H2O
	1,  //N2

	//M(1)   //all species in mixture contribute equally
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	1,  //H2
	1,  //HO2
	1,  //H2O2
	1,  //H2O
	1,  //N2
	
	//M(2)
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	1,  //H2
	1,  //HO2
	1,  //H2O2
	5.,  //H2O
	1,  //N2
	
	//M(3)
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	1,  //H2
	1,  //HO2
	1,  //H2O2
	5.,  //H2O
	1,  //N2

	//M(4)      //all species in mixture contribute equally
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	1,  //H2
	1,  //HO2
	1,  //H2O2
	1,  //H2O
	1,  //N2
	
	//M(5)    //all species in mixture contribute equally
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	1,  //H2
	1,  //HO2
	1,  //H2O2
	1,  //H2O
	1  //N2

      };
      return ce;
    }
    

    //! note that the molar masses (also know as molecular weights) 
    //! which can be computed from the periodic table are stated in g/mol. 
    //! Here we use SI units, i.e. kg/mol. Therefore, the
    //! molar masses \f$M_k, \ k=1, \ldots, K\f$ must be devided by 1000, after
    //! having been calculated from the periodic table, as done below.
    inline static T* molar_masses(){
      //O, O2, H, OH, H2, HO2, H2O2, H2O, N2
      static T mm[] ={1.600000e-02, 3.200000e-02, 1.007900e-03, 1.700790e-02, 2.015800e-03, 3.300790e-02, 3.401580e-02, 1.801580e-02, 2.802000e-02 };
      
      return mm;
    }


    inline static T* transport(){
      static T tpt[] ={
 //Geometry    epsilon/kappa_B      sigma      mu               alpha      Z_rot
 //O
0.00000e+00,   8.00000e+01,   2.75000,   0.00000e+00,   0.00000e+00,   0.00000e+00,   
//O2
1.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   1.60000e-30,   3.80000e+00,
//H   
0.00000e+00,   1.45000e+02,   2.05000,   0.00000e+00,   0.00000e+00,   0.00000e+00, 
//OH  
1.00000e+00,   8.00000e+01,   2.75000,   0.00000e+00,   0.00000e+00,   0.00000e+00,
//H2   
1.00000e+00,   3.80000e+01,   2.92000,   0.00000e+00,   7.90000e-31,   2.80000e+02,   
//HO2
2.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   0.00000e+00,   1.00000e+00,   
//H2O2
2.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   0.00000e+00,   3.80000e+00,   
//H2O
2.00000e+00,   5.72400e+02,   2.60500,   1.84400e-28,   0.00000e+00,   4.00000e+00,
//N2   
1.00000e+00,   9.75300e+01,   3.62100,   0.00000e+00,   1.76000e-30,   4.00000e+00 
      };

      return tpt;
    }

    //! # of bidirec reactions
    //! when \rev is given, compute reverse reaction via Arrhenius law
    // ! 1 = yes, 0 = no (then compute it like always)
     static inline bool* explicit_reverse_reaction_coefficients(){
      static bool revcoef[] = {
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
	0
      };
      return revcoef;
     }

    //! tells us if species k participates in any (at least one) reaction
    //! true if it does, false if not
    static inline bool* is_species_reactive(){
      static bool rtive[] = {
	true,   
	true,   
	true,
	true,
	true,
	true,
	true,
	true,
	false   //N2 does not participate in reactions (only as third body)
      };
      return rtive;
     }


     //====================== TODO: The following  might be changed  ===========
    //=====================        from time to time                 ===========
    
    /*
    //! 3 rpvs involved
    enum{
      rednspec = 3,  //reduced dimension
      rednreac = 24
    };
    
    static inline index_type* rpv_index() {
      static index_type rpv[] = {1,4,7};
      
      return rpv;
    }
  
    static inline index_type* rpv_reaction_index() {
      static index_type rri[] = {0,2,3,5,6,7,8,9,10,11,12,13,14,15,17,18,19,21,22,23,24,25,26,27};
      
      return rri;
    }

    //! actually not needed since the construction can be obtained via
    //! species_names() and rpv_index()
    */

    
      //! 2 rpvs involved -- here you loose more information
    /*enum{
      rednspec = 2,  //reduced dimension
      rednreac = 16
    };
    
    static inline index_type* rpv_index() {
      static index_type rpv[] = {4,7};
      
      return rpv;
    }
  
    static inline index_type* rpv_reaction_index() {
      static index_type rri[] = {2,7,10,11,12,13,14,15,17,18,19,21,22,23,24,27};
      
      return rri;
    }

     
    //!compile-time versions
    typedef Element<4,Element<7> > RpvIndexType;
     
    typedef Element<2,Element<7,Element<10,Element<11, Element<12, Element<13, Element<14, Element<15, Element<17, Element<18, Element<19, Element<21, Element<22, Element<23, Element<24, Element<27> > > > > > > > > > > > > > > > RpvReactionIndexType;
    */

  };  //end of class 




 /**
   * \brief HAND-CODED H2/O2 mechanism due to [1], [2]
   *
   * 
   *
   * References:
   * [1] <a href="https://www-pls.llnl.gov/data/docs/science_and_technology/chemistry/combustion/h2_v1b_mech.txt"> H2 mech </a>
   *
   * [2] [M. O'. CONAIRE, H. J. CURRAN, J. M. SIMMIE,  W. J. PITZ and C.K. WESTBROOK, "A Comprehensive Modeling Study of Hydrogen Oxidation", Internationa Journal of Chemical Kinetics, 36 (11), 2004, pp. 603--622] 
   *
   * Note: 1.) Activation energies have been originally stated in calories (1 cal = 4.184 J)
           2.) There are 2 duplicate reactions. According to [2, p. 9, footnote h)], these reactions are expressed as sum of the corresponding rate expressions, thus yielding 19 (instead of 21) reactions in total.
	   3.) Except for the fall-off reactions, reverse rate coefficients are calculated by tabulated values for A_rev,i, b_rev,i, Ea_rev,i. 
	   4.) N2 and Ar don't participate in the reactions. Ar plays a role as
  third body in several bi- and trimolecular reactions.
   */
  template<class T> 
  class ThermoData4Mechanism<T,10>{ 
  public:
    typedef size_t index_type;
    typedef index_type* index_type_pointer;
    typedef int int_type;
    typedef int_type* int_type_pointer;
    typedef T value_type;
    typedef T* value_type_pointer; 
    //this is the order of which the species are to be considered
    enum{O, O2, H, OH, H2, HO2, H2O2, H2O, N2, AR};

    enum{
      ENCOD=10
    };
    
    enum{
      nspec = 10,
      nreac = 
#ifdef COUNT_DUPLICATE_REACTIONS
      21,
#else
      19,  //no duplicates
#endif     
      ntroereac = 2,
      ntb = 6,           //number of 3rd bodies
      nelem = 3         //number of chem. elements in mechanism
    };

    static inline std::string* species_names() {
      static std::string names[] = {
	"O","O2","H","OH", "H2","HO2","H2O2","H2O","N2","Ar"
      };
      return names;
    }

    
    static inline T density_of_mixture() {
      return 0.3;    //in kg/m³
    }

    static inline T* thermo(){
      static T th[] = {
	//O 
	0.02542060e+02,-0.02755062e-03,-0.03102803e-07, 0.04551067e-10,-0.04368052e-14,    
	0.02923080e+06, 0.04920308e+02, 0.02946429e+02,-0.01638166e-01, 0.02421032e-04,    
	-0.01602843e-07, 0.03890696e-11, 0.02914764e+06, 0.02963995e+02,              //O2
	0.03697578e+02, 0.06135197e-02,-0.01258842e-05, 0.01775281e-09,-0.01136435e-13,    
	-0.01233930e+05, 0.03189166e+02, 0.03212936e+02, 0.01127486e-01,-0.05756150e-05,    
	0.01313877e-07,-0.08768554e-11,-0.01005249e+05, 0.06034738e+02,
	//H
	0.02500000e+02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,    
	0.02547163e+06,-0.04601176e+01, 0.02500000e+02, 0.00000000e+00, 0.00000000e+00,    
	0.00000000e+00, 0.00000000e+00, 0.02547163e+06,-0.04601176e+01,
	//OH
	2.85376040e+00, 1.02994334e-03,-2.32666477e-07, 1.93750704e-11,-3.15759847e-16,    
	3.69949720e+03, 5.78756825e+00, 3.41896226e+00, 3.19255801e-04,-3.08292717e-07,    
	3.64407494e-10,-1.00195479e-13, 3.45264448e+03, 2.54433372e+00,
	//H2
	0.02991423e+02, 0.07000644e-02,-0.05633829e-06,-0.09231578e-10, 0.01582752e-13,    
	-0.08350340e+04,-0.01355110e+02, 0.03298124e+02, 0.08249442e-02,-0.08143015e-05,    
	-0.09475434e-09, 0.04134872e-11,-0.01012521e+05,-0.03294094e+02,
	//HO2
	4.01721090e+00, 2.23982013e-03,-6.33658150e-07, 1.14246370e-10,-1.07908535e-14,    
	1.11856713e+02, 3.78510215e+00, 4.30179801e+00,-4.74912051e-03, 2.11582891e-05,    
	-2.42763894e-08, 9.29225124e-12, 2.94808040e+02, 3.71666245e+00,
	//H2O2
	0.04573167e+02, 0.04336136e-01,-0.01474689e-04, 0.02348904e-08,-0.01431654e-12,    
	-0.01800696e+06, 0.05011370e+01, 0.03388754e+02, 0.06569226e-01,-0.01485013e-05,    
	-0.04625806e-07, 0.02471515e-10,-0.01766315e+06,0.06785363e+02,
	//H2O
	0.02672146e+02, 0.03056293e-01,-0.08730260e-05, 0.01200996e-08,-0.06391618e-13,    
	-0.02989921e+06, 0.06862817e+02, 0.03386842e+02, 0.03474982e-01,-0.06354696e-04,    
	0.06968581e-07,-0.02506588e-10,-0.03020811e+06, 0.02590233e+02,
	//N2
	0.02926640e+02, 0.01487977e-01,-0.05684761e-05, 0.01009704e-08,-0.06753351e-13,    
	-0.09227977e+04, 0.05980528e+02, 0.03298677e+02, 0.01408240e-01,-0.03963222e-04,    
	0.05641515e-07,-0.02444855e-10,-0.01020900e+05, 0.03950372e+02,
	//Ar
	0.02500000e+02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,    
	-0.07453750e+04, 0.04366001e+02, 0.02500000e+02, 0.00000000e+00, 0.00000000e+00,    
	0.00000000e+00, 0.00000000e+00,-0.07453750e+04, 0.04366001e+02 
   };


      return th;
    }

    

    static inline T* temperature_bounds(){
      static T tbds[] = {273,5000,1000,  //O   //previously 
			 273,5000,1000,  //O2
			 273,5000,1000,  //H
			 273,5000,1710,  //OH
			 273,5000,1000,  //H2
			 200,3500,1000,  //HO2
			 273,5000,1000,  //H2O2
			 273,5000,1000,  //H2O
			 273,5000,1000,  //N2
			 273,5000,1000   //Ar    
      };
    
      return tbds;
    }


    static inline index_type* stoichiometric_matrix(){
      //store only ni' 
      //each column corresponds to 
      //O,  O2, H,  OH, H2,HO2,H2O2,H2O,N2 and Ar
      //respectively. Each <--> reaction is split into a --> and a <-- reaction
      //3rd bodies (+m) are not count 
      static index_type sm[] = {
	0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  //h+o2 = o+oh
	1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  
	1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  //h+o2 = o+oh
	0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  //oh+h2 = h+h2o
	0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  
	1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  //o+h2o = oh+oh
	0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  //h2+m = h+h+m
	0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  
	0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  //o2+m = o+o+m
	2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  //oh+m = o+h+m
	1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  //h2o+m = h+oh+m
	0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  
	0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  //h+o2(+m) = ho2(+m)  TROE
	0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  
	0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  //ho2+h = h2+o2
	0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  
	0,  0,  1,  0,  0,  1,  0,  0,  0,  0,  //ho2+h = oh+oh
	0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  
	1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  //ho2+o = oh+o2
	0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  
	0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  //ho2+oh = h2o+o2 
	0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  
	0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  //h2o2+o2 = ho2+ho2
	0,  0,  0,  0,  0,  2,  0,  0,  0,  0, 
#ifdef COUNT_DUPLICATE_REACTIONS
	0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  //h2o2+o2 = ho2+ho2 DUPLICATE
	0,  0,  0,  0,  0,  2,  0,  0,  0,  0, 
#endif 
	0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  //h2o2(+m) = oh+oh(+m) TROE
	0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  
	0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  //h2o2+h = h2o+oh
	0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  
	0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  //h2o2+h = h2+ho2
	0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  
	1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  //h2o2+o = oh+ho2
	0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  
	0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  //h2o2+oh = h2o+ho2
	0,  0,  0,  0,  0,  1,  0,  1,  0,  0
#ifdef COUNT_DUPLICATE_REACTIONS
       ,0,  0,  0,  1,  0,  0,  1,  0,  0,  0, //h2o2+oh = h2o+ho2 DUPLICATE
	0,  0,  0,  0,  0,  1,  0,  1,  0,  0
#endif 
    };

      return sm;
    }

    //TODO: insert well-scaled parameters
    //! NOTE: these are already in SI format, i.e. they were transformed 
    //! accordingly once before
    static inline T* arrhenius_values(){
      static T arrh[] ={ //FORWARD
	1.915000e+08,   0.000000e+00,   6.878496e+04, //h+o2 = o+oh
	5.080000e-02,   2.670000e+00,   2.632573e+04, //o+h2 = h+oh 
	2.160000e+02,   1.510000e+00,   1.435112e+04, //oh+h2 = h+h2
	2.970000e+00,   2.020000e+00,   5.606560e+04, //o+h2o = oh+oh
	4.577000e+13,   -1.400000e+00,   4.368096e+05, // h2+m = h+h+m
	4.515000e+11,   -6.400000e-01,   4.974776e+05, //o2+m = o+o+m 
	9.880000e+11,   -7.400000e-01,   4.271864e+05, //oh+m = o+h+m
	1.912000e+17,   -1.830000e+00,   4.958040e+05, //h2o+m = h+oh+m
	1.475000e+06,   6.000000e-01,   0.000000e+00, //TROE h+o2(+m) = ho2(+m)
	1.660000e+07,   0.000000e+00,   3.443432e+03, //ho2+h = h2+o2
	7.079000e+07,   0.000000e+00,   1.234280e+03, //ho2+h = oh+oh 
	3.250000e+07,   0.000000e+00,   0.000000e+00, // ho2+o = oh+o2
	2.890000e+07,   0.000000e+00,   -2.079448e+03,//ho2+oh = h2o+o2 
#ifdef COUNT_DUPLICATE_REACTIONS  //DUP h2o2+o2 = ho2+ho2
        4.634000e+10,  -3.500000e-01,  2.120033e+05,    //fwd
	1.434000e+07,  -3.500000e-01,  1.550590e+05,    //fwd                   
#else
	4.635434e+10,   -7.000000e-01,   3.670623e+05, //DUPsingle h2o2+o2 = ho2+ho2
#endif
	2.951000e+14,   0.000000e+00,   2.026311e+05, //h2o2(+m) = oh+oh(+m)
	2.410000e+07,   0.000000e+00,   1.661048e+04, //h2o2+h = h2o+oh
	6.025000e+07,   0.000000e+00,   3.326280e+04, //h2o2+h = h2+ho2
	9.550000e+00,   2.000000e+00,   1.661048e+04, //h2o2+o = oh+ho2
#ifdef COUNT_DUPLICATE_REACTIONS  //DUP	h2o2+oh = h2o+ho2
	1.000000e+06,  0.000000e+00,  0.000000e+00,    //fwd
	5.800000e+08,  0.000000e+00,  3.998649e+04,    //fwd
#else
	5.810000e+08,   0.000000e+00,   3.998649e+04, //DUPsingle h2o2+oh = h2o+ho2
#endif
     
	//!starting from 3*nreac, the troe coefficients are stored in 7-entries
	//!blocks -- each containing A_low, beta_low, Ea_low, alpha,T***,T* and 
	//! T** in this order
//TROE1: A_low            beta_low          Ea_low
	3.4820E+04, -4.1100E-01, -4.66516E+03, //low forward
//alpha   T***   T*       T**
	0.5,  1.0000E-30,  1.0000E+30,  1.0000E+100,
//TROE2: A_low            beta_low          Ea_low
	//1.202E+12,  
	1.202E+11, 0.00,  190372,
//alpha   T***   T*       T**	
	0.5, 1.0e-30, 1.0e+30, 1.0e+100,



	//REVERSE RATE: EXPLICITELY STATED
	5.481000e+05,   3.900000e-01,   -1.225912e+03,  //h+o2 = o+oh
	2.667000e-02,   2.650000e+00,   2.041792e+04,   //o+h2 = h+oh
	2.298000e+03,   1.400000e+00,   7.665088e+04,   //oh+h2 = h+h2
	1.465000e-01,   2.110000e+00,   -1.215034e+04,  //o+h2o = oh+oh
	1.146000e+08,   -1.680000e+00,   3.430880e+03,  // h2+m = h+h+
	6.165000e+03,   -5.000000e-01,   0.000000e+00,  //o2+m = o+o+m 
	4.714000e+06,   -1.000000e+00,   0.000000e+00,  //oh+m = o+h+m
	4.500000e+10,   -2.000000e+00,   0.000000e+00, //h2o+m = h+oh+m
	3.090000e+12,   5.300000e-01,   2.044721e+05, //! Troe h+o2(+m) = ho2(+m)
	3.164000e+06,   3.500000e-01,   2.322538e+05, //ho2+h = h2+o2
	2.027000e+04,   7.200000e-01,   1.541386e+05, //ho2+h = oh+oh 
	3.252000e+06,   3.300000e-01,   2.229235e+01,// ho2+o = oh+o2 3.250000e+07,   0.000000e+00,   0.000000e+00,
	5.861000e+07,   2.400000e-01,   2.890307e+05, //ho2+oh = h2o+o2 
#ifdef COUNT_DUPLICATE_REACTIONS
        4.200000e+08,  0.000000e+00,  5.012432e+04,     //rev
	1.300000e+05,  0.000000e+00,  -6.815736e+03,    //rev
#else
	4.201300e+08,   0.000000e+00,   4.330858e+04,  //DUP h2o2+o2 = ho2+ho2
#endif
	3.656000e+02,   1.140000e+00,   -1.081146e+04, //! TROE  h2o2(+m) = oh+oh(+m)
	1.269000e+02,   1.310000e+00,   2.987794e+05, //h2o2+h = h2o+oh
	1.041000e+05,   7.000000e-01,   1.002068e+05, //h2o2+h = h2+ho2
	8.660000e-03,   2.680000e+00,   7.765504e+04, //h2o2+o = oh+ho2
#ifdef COUNT_DUPLICATE_REACTIONS
	1.838000e+04,  5.900000e-01,  1.292438e+05,    //rev
	1.066000e+07,  5.900000e-01,  1.692428e+05,    //rev
#else
	1.067838e+07,   1.180000e+00,   2.984866e+05,  //DUPsingle h2o2+oh = h2o+ho2
#endif

	//======= TODO: check if this is correct; I just take the low coeffs
	//              and, according to the molecularity adjust them.
	//              Here, however it can be skipped since the Troe reacs
	//              compute reaction rates via k_fwd,i/K_c,i.

	//TROE1: A_low            beta_low          Ea_low
	3.4820E+10, -4.1100E-01, -4.66516E+03,
//alpha   T***   T*       T**
	0.5,  1.0000E-30,  1.0000E+30,  1.0000E+100,
//TROE2: A_low            beta_low          Ea_low
	1.202E+05,  0.00,  190372,
//alpha   T***   T*       T**	
	0.5, 1.0e-30, 1.0e+30, 1.0e+100
 };

      return arrh;
    }

    
    

    //! troe with third bodies information
    static inline int_type_pointer troewtb(){
      static int tw[] = {0, //# of bidirec. reacs
			 0,   //TROE information part
			 0,
			0,
			0,
			0,
			0,
			0,
			1,
			0,
			0,
			0,
			0,
			0,
#ifdef COUNT_DUPLICATE_REACTIONS
			 0,          //no 3rd body in duplicate reaction 
#endif
			1,
			0,
			0,
			0,
			0, 
#ifdef COUNT_DUPLICATE_REACTIONS
			 0,          //no 3rd body in duplicate reaction 
#endif

			 //now comes information about 3rd bodies: # of 3rd bodies, -1 when no 3rd body is involved
			 -1,
			 -1,
			 -1,
			 -1,
			 0,
			 1,
			 2,
			 3,
			 4,
			 -1,
			 -1,
			 -1,
			 -1,
			 -1,
#ifdef COUNT_DUPLICATE_REACTIONS
			 -1,          //no 3rd body in duplicate reaction 
#endif
			 5,
			 -1,
			 -1,
			 -1,
			 -1
#ifdef COUNT_DUPLICATE_REACTIONS
			 ,-1          //no 3rd body in duplicate reaction 
#endif
			 
      };

      return tw;
    }

  
    //! 6 THIRD BODIES x 10 SPECIES = 60 entries; default: 1
    static inline T* collision_efficiencies(){
      
      static T ce[] = {
	//M(0)
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	2.5,  //H2
	1,  //HO2
	1,  //H2O2
	12.,  //H2O
	1,  //N2
	1,  //Ar
	//M(1)
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	2.5,  //H2
	1,  //HO2
	1,  //H2O2
        12.,  //H2O
	1,  //N2
	0.83,  //Ar
	//M(2)
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	2.5,  //H2
	1,  //HO2
	1,  //H2O2
	12.,  //H2O
	1,  //N2
	0.75,  //Ar
	//M(3)
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	0.73,  //H2
	1,  //HO2
	1,  //H2O2
	12.,  //H2O
	1,  //N2
	0.38,  //Ar
	//M(4)              //TROE
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	1.3,  //H2
	1,  //HO2
	1,  //H2O2
	14.,  //H2O
	1,  //N2
	0.67,  //Ar
	//M(5)             //TROE
	1,  //O
	1,  //O2
	1,  //H
	1,  //OH
	2.5,  //H2
	1,  //HO2
	1,  //H2O2
	12.,  //H2O
	1,  //N2
	0.64  //Ar

      };
      return ce;
    }
    

    //! note that the molar masses (also know as molecular weights) 
    //! which can be computed from the periodic table are stated in g/mol. 
    //! Here we use SI units, i.e. kg/mol. Therefore, the
    //! molar masses \f$M_k, \ k=1, \ldots, K\f$ must be devided by 1000, after
    //! having been calculated from the periodic table, as done below.
    inline static T* molar_masses(){
      //O, O2, H, OH, H2, HO2, H2O2, H2O, N2, Ar
      static T mm[] ={1.600000e-02, 3.200000e-02, 1.007900e-03, 1.700790e-02, 2.015800e-03, 3.300790e-02, 3.401580e-02, 1.801580e-02, 2.802000e-02, 3.9948e-02 };
      
      return mm;
    }


    inline static T* transport(){
      static T tpt[] ={
 //Geometry    epsilon/kappa_B      sigma      mu               alpha      Z_rot
 //O
0.00000e+00,   8.00000e+01,   2.75000,   0.00000e+00,   0.00000e+00,   0.00000e+00,   
//O2
1.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   1.60000e-30,   3.80000e+00,
//H   
0.00000e+00,   1.45000e+02,   2.05000,   0.00000e+00,   0.00000e+00,   0.00000e+00, 
//OH  
1.00000e+00,   8.00000e+01,   2.75000,   0.00000e+00,   0.00000e+00,   0.00000e+00,
//H2   
1.00000e+00,   3.80000e+01,   2.92000,   0.00000e+00,   7.90000e-31,   2.80000e+02,   
//HO2
2.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   0.00000e+00,   1.00000e+00,   
//H2O2
2.00000e+00,   1.07400e+02,   3.45800,   0.00000e+00,   0.00000e+00,   3.80000e+00,   
//H2O
2.00000e+00,   5.72400e+02,   2.60500,   1.84400e-28,   0.00000e+00,   4.00000e+00,
//N2   
1.00000e+00,   9.75300e+01,   3.62100,   0.00000e+00,   1.76000e-30,   4.00000e+00,
//Ar
0,             136.500,       3.330,     0.000,         0.000,         0.000
      };

      return tpt;
    }


    //! # of bidirec reactions
    //! when "rev /" is given, compute reverse reaction via Arrhenius law
    // ! 1 = yes, 0 = no (then compute it like always)
    //! NOTE: As soon as a '1' appears you must provide arrhenius coeffs, and 
    //! , possibly, coefficients for TROE-fall-off reactions in 'arrhenius_values()', too.
     static inline bool* explicit_reverse_reaction_coefficients(){
      static bool revcoef[] = {
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	1,
	0,  //9.) TROE -- "rev /" is commented out (!)
	1,
	1,
	1,
	1,
	1,
#ifdef COUNT_DUPLICATE_REACTIONS
	1,          //rev is given 
#endif
	0,  //15.) resp. 16.) for DUP reaction TROE -- "rev /" is commented out (!)
	1,
	1,
	1,
	1
#ifdef COUNT_DUPLICATE_REACTIONS
	,1          //rev is given 
#endif	
      };
      return revcoef;
     }


     //! tells us if species k participates in any (at least one) reaction
    //! true if it does, false if not
    static inline bool* is_species_reactive(){
      static bool rtive[] = {
	true,   
	true,   
	true,
	true,
	true,
	true,
	true,
	true,
	false,   //N2 does not participate in reactions (only as third body)
	false    //... neither does Argon (except as third body)
      };
      return rtive;
     }

    //===== TODO · TODO · TODO ... if it is every working....
    //===== TODO · TODO · TODO ... if it is every working....
    //===== TODO · TODO · TODO ... if it is every working....
     //====================== TODO: The following  might be changed  ===========
    //=====================        from time to time                 ===========
    enum{
      rednspec = 4,  //reduced dimension
      rednreac = 19   //all reactions involved
    };
    
    static inline index_type* rpv_index() {
      static index_type rpv[] = {1,3,4,7};
      
      return rpv;
    }
  
    static inline index_type* rpv_reaction_index() {
      static index_type rri[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
      
      return rri;
    }

     
    //!compile-time versions
    typedef Element<1,Element<3,Element<4,Element<7> > > > RpvIndexType;
     
    typedef Element<0,Element<1,Element<2,Element<3, Element<4, Element<5, Element<6, Element<7, Element<8, Element<9, Element<10, Element<11, Element<12, Element<13, Element<14, Element<15, Element<16, Element<17, Element<18> > > > > > > > > > > > > > > > > > >  RpvReactionIndexType;


    inline static T* mass_balance_matrix(){
      static T mbm[] ={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};   //! sum of mass fractions equals one
      return mbm;
    }

    inline static T* mass_sum(){
      static T ms[] ={1.};
      return ms;
    }

  };  //end of class 

 
}  //end namespace 

#endif
