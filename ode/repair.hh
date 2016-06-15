#ifndef REPAIR_PHYSICAL_QUANTITIES_HH
#define REPAIR_PHYSICAL_QUANTITIES_HH

#include "constants.hh"
#include "../common/fancymessages.hh"
#include "../common/date.hh"
#include "../common/globalfunctions.hh"
#include "../io/readinparameters.hh"
#include "../templatemetaprograms/unrollloop.hh"
#include "../containers/staticarray.hh"
#include "../expressiontemplates/exprvec.hh"

#include "../massactionkinetics/physicalconstants.hh"
#include "../fdm/primitivevars.hh"

namespace Adonis{

  /**
   * \brief Choose an appropriate number when a mass fraction is beyond 1
   */
  template<class V> 
  inline typename V::value_type& nearby_value_2D(typename V::value_type& avg, const V& u, int i, int j, int nx, int ny, int k){
    //acces right index only
    adonis_assert((i>= 0) && (i <= nx-1) && (j>= 0) && (j <= ny-1));
    adonis_assert(4+k >= 4); //only meant for chem. species
    
    avg = 0; //reset 
    int m = 0;

    //interior
    if(((i >= 1) && (i <= nx-2)) && ((j >= 1) && (j <= ny-2))){
      
      for(int s = -1; s <= +1; ++s){
	for(int t = -1; t <= +1; ++t){ 
	  if(is_properly_contained(FlowVariables<V>::primitive(i+s,j+t,u,4+k,nx,ny),0.,1.)){
	    avg += FlowVariables<V>::primitive(i+s,j+t,u,4+k,nx,ny);
	    m++;
	  }
	}
      }
    }
    
    //left
    if( (i==0) && ((j >=1) && (j <= ny-2))){
      //(0,j+1)
      if(is_properly_contained(FlowVariables<V>::primitive(0,j+1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(0,j+1,u,4+k,nx,ny);
	m++;
      }
      //(1,j+1)
      if(is_properly_contained(FlowVariables<V>::primitive(1,j+1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(1,j+1,u,4+k,nx,ny);
	m++;
      }
      //(1,j)
      if(is_properly_contained(FlowVariables<V>::primitive(1,j,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(1,j,u,4+k,nx,ny);
	m++;
      }
      //(1,j-1)
      if(is_properly_contained(FlowVariables<V>::primitive(1,j-1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(1,j-1,u,4+k,nx,ny);
	m++;
      }
      //(0,j-1)
      if(is_properly_contained(FlowVariables<V>::primitive(0,j-1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(0,j-1,u,4+k,nx,ny);
	m++;
      }
    }
    
    //right
    if((i == nx-1) && ((j >= 1) && (j <= ny-2))){
      //std::cout << "RIGHT"<<std::endl;
      //(nx-1,j+1)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-1,j+1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-1,j+1,u,4+k,nx,ny);
	m++;
      }
      //(nx-2,j+1)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-2,j+1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-2,j+1,u,4+k,nx,ny);
	m++;
      }
      //(nx-2,j)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-2,j,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-2,j,u,4+k,nx,ny);
	m++;
      }
      //(nx-2,j-1)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-2,j-1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-2,j-1,u,4+k,nx,ny);
	m++;
      }
      //(nx-1,j-1)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-1,j-1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-1,j-1,u,4+k,nx,ny);
	m++;
      }
    }
    //up
    if(((i >= 1) && (i <= nx-2)) && (j == ny-1)){
      //(i-1,ny-1)
      if(is_properly_contained(FlowVariables<V>::primitive(i-1,ny-1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i-1,ny-1,u,4+k,nx,ny);
	m++;
      }
      //(i-1,ny-2)
      if(is_properly_contained(FlowVariables<V>::primitive(i-1,ny-2,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i-1,ny-2,u,4+k,nx,ny);
	m++;
      }
      //(i,ny-2)
      if(is_properly_contained(FlowVariables<V>::primitive(i,ny-2,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i,ny-2,u,4+k,nx,ny);
	m++;
      }
      //(i+1,ny-2)
      if(is_properly_contained(FlowVariables<V>::primitive(i+1,ny-2,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i+1,ny-2,u,4+k,nx,ny);
	m++;
      }
      //(i+1,ny-1)
      if(is_properly_contained(FlowVariables<V>::primitive(i+1,ny-1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i+1,ny-1,u,4+k,nx,ny);
	m++;
      }
    }
    //down
    if(((i >= 1) && (i <= nx-2)) && (j == 0)){
      //(i-1,1)
      if(is_properly_contained(FlowVariables<V>::primitive(i-1,1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i-1,1,u,4+k,nx,ny);
	m++;
      }
      //(i-1,0)
      if(is_properly_contained(FlowVariables<V>::primitive(i-1,0,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i-1,0,u,4+k,nx,ny);
	m++;
      }
      //(i,1)
      if(is_properly_contained(FlowVariables<V>::primitive(i,1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i,1,u,4+k,nx,ny);
	m++;
      }
      //(i+1,1)
      if(is_properly_contained(FlowVariables<V>::primitive(i+1,1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i+1,1,u,4+k,nx,ny);
	m++;
      }
      //(i+1,0)
      if(is_properly_contained(FlowVariables<V>::primitive(i+1,0,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(i+1,0,u,4+k,nx,ny);
	m++;
      }
    }

    //CORNERS: special treatment
    //up-left corner
    if((i == 0) && (j == ny-1)){
      //(0,ny-2)
      if(is_properly_contained(FlowVariables<V>::primitive(0,ny-2,u,4+k,nx,ny),0.,1.)){
	 avg += FlowVariables<V>::primitive(0,ny-2,u,4+k,nx,ny);
	 m++;
      }
      //(1,ny-2)
      if(is_properly_contained(FlowVariables<V>::primitive(1,ny-2,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(1,ny-2,u,4+k,nx,ny);
	m++;
      }
      //(1,ny-1)
      if(is_properly_contained(FlowVariables<V>::primitive(1,ny-1,u,4+k,nx,ny),0.,1.)){
	 avg += FlowVariables<V>::primitive(1,ny-1,u,4+k,nx,ny);
	 m++;
      }
    }
    //up_right
    if((i == nx-1) && (j == ny-1)){
      //(nx-2,ny-1)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-2,ny-1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-2,ny-1,u,4+k,nx,ny);
	m++;
      }
      //(nx-2,ny-2)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-2,ny-2,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-2,ny-2,u,4+k,nx,ny);
	m++;
      }
      //(nx-1,ny-2)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-1,ny-2,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-1,ny-2,u,4+k,nx,ny);
	m++;
      }
    }
    //left-down
    if((i == 0) && (j == 0)){
      //(1,0)
      if(is_properly_contained(FlowVariables<V>::primitive(1,0,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(1,0,u,4+k,nx,ny);
	m++;
      }
      //(1,1)
      if(is_properly_contained(FlowVariables<V>::primitive(1,1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(1,1,u,4+k,nx,ny);
	m++;
      }
      //(0,1)
      if(is_properly_contained(FlowVariables<V>::primitive(0,1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(0,1,u,4+k,nx,ny);
	m++;
      }
    }
    //right-down
    if((i == nx-1) && (j == 0)){
      //(nx-1,1)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-1,1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-1,1,u,4+k,nx,ny);
	m++;
      }
      //(nx-2,1)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-2,1,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-2,1,u,4+k,nx,ny);
	m++;
      }
      //(nx-2,0)
      if(is_properly_contained(FlowVariables<V>::primitive(nx-2,0,u,4+k,nx,ny),0.,1.)){
	avg += FlowVariables<V>::primitive(nx-2,0,u,4+k,nx,ny);
	m++;
      }
    }
    
    //!calculate average value of reasonable values
    if( m != 0)
      avg /= m;
    if((m == 0) || (avg >= 1)){
      std::cout << "nearby_value_2D: no reasonable value found in the vicinity of current point. Take random number from [0,1) instead."<< std::endl;
      avg = (typename V::value_type)rand()/(typename V::value_type)(RAND_MAX-1); //take value from [0,1)
    }
    if(m > 8){
      ADONIS_ERROR(DerivedError,"Oooops, m = "<< m << ". Hasn't m been reset?");
    }
  
    return avg;
  }

 

  /**
   * \brief Imposing physical ranges to quantities of interest reduces
   * the risk of undesirable results during the computation, such as negative
   * density values, negative mass fractions, etc.
   * default: do nothing
   */
  template<class K, int DIM, bool REPAIR, char REPAIRTYPE> 
  class RepairPhysicalQuantity{
  public:
    typedef K value_type;
    typedef typename TypeAdapter<K>::BaseType BaseType;
    typedef std::vector<int> InfoVecType;
   
    
    RepairPhysicalQuantity():numOfSpec_(0){
      UnrollLoop<0,20>::assign_value(ViolatedVars_,false);
    }

    template<class V>
    void init(const std::string& file, const V& v, const std::string& fn = std::string()){
      UnrollLoop<0,20>::assign_value(ViolatedVars_,false);
      msg_.clear();
      dim_ = 0;
    } //do nothing

    template<class V>
    void repair(V& u, int excessSpecIndex){}  //do nothing here

    void info(int timeStep = 0, int NewtonIter = 0, const std::string& str = std::string()){

      if(ViolatedVars_[0] == true)
	msg_ += " Y_k < 0.0 ";
      if(ViolatedVars_[1] == true)
	msg_ += " Y_k > 1.0 ";
      if(ViolatedVars_[2] == true)
	msg_ += " Y_k = NaN ";
      
      if(msg_.size() != 0){
	//additional information provided
	if(str.size() != 0){
	  std::cout << "#### "<<str << std::endl;
	}
	msg_.insert(0,"VARIABLE BOUNDS VIOLATED in time step "+Num2str(timeStep)+ ", in Newton iter " +Num2str(NewtonIter)+ " : \n");
	FM_.nice_output(msg_,35);
      }
	
      //reset
      msg_.clear();
      UnrollLoop<0,20>::assign_value(ViolatedVars_,false);
    }

    bool are_phys_quantities_violated() const {return false;} //not violated
    
    const StaticArray<int,3>& space_dimensions() const {return dim_;}
    
  private:
    int numOfSpec_;
    bool ViolatedVars_[20]; //Y_k < 0, Y_k >0, Y_NaN, T<Tmin, T>Tmax, T_NaN
    std::string msg_;
    FancyMessages FM_;

    StaticArray<int,3> dim_;
  };

   //! partial specialization for rectangular FD 1d MOL
    //! do nothing for 1d MOL so far. 
    //! TODO: if you want to use it write your own specialization
  template<class K, bool REPAIR, char REPAIRTYPE> 
  class RepairPhysicalQuantity<K,1,REPAIR,REPAIRTYPE>{
  public:
    typedef K value_type;
      
    typedef typename TypeAdapter<K>::BaseType BaseType;
    typedef std::vector<int> InfoVecType;

    RepairPhysicalQuantity(){
      dim_ = 0;
    }
			     
    template<class V>
    void init(const std::string& file, const V& v, const std::string& fn = std::string()){} //do nothing
      
    template<class V>
    void repair(V& u, int excessSpecIndex){} //do nothing

    void info(int timeStep = 0, int NewtonIter = 0, const std::string& str = std::string()){}

    bool are_phys_quantities_violated() const {return false;} //not violated 

    const StaticArray<int,3>& space_dimensions() const {return dim_;}
  private:
    StaticArray<int,3> dim_;
  };

  // 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D
  // 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D
  // 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D · 2D
  //!partial specialization for rectangular FD 2d MOL
  template<class K>
  class RepairPhysicalQuantity<K,2,true,'f'>{
  public:
    typedef K value_type;
    typedef typename TypeAdapter<K>::BaseType BaseType;
    typedef std::vector<int> InfoVecType;

    RepairPhysicalQuantity():nx_(0),ny_(0),npt_(0),numOfSpec_(0),rhoSmall_(1.e-05), v1MaxAbs_(0.), v2MaxAbs_(0.), Tmin_(280.), Tmax_(2000.), SmallSpec_(0.0), sm_(0),violation_(false),maxk_(0), count0_(0), count1_(0), countLarge_(0), countOk_(0), iterCount_(0),sum_(0), partsum_(0),tmp_(0), mean_(0), approxEqual_(1.e-09),avg_(0){
       
      UnrollLoop<0,20>::assign_value(ViolatedVars_,false);
	
    }

    ~RepairPhysicalQuantity(){of_.close();}

    template<class V>
    void init(const std::string& file, const V& v, const std::string& outfile = "Recorded_Repair_History.dat"){
      //at the very beginning we assume a physically relevant composition
      //uGood_ = v; 
      
      FancyMessages FMSG;
      FMSG.nice_output("REPAIR PHYSICAL QUANTITIES FOR 2D RECTANGULAR GRID WITH MOL APPLIED INSIDE ODE SOLVER...",35);

      ParameterData PD;
      PD.read_from_file(file);

     
      if(PD.get_datum<unsigned>("fromLastStepOn") > 0){
	if(does_file_exist(outfile) == true){
	  of_.open(outfile.c_str(), std::ios_base::app); //append to existing file
	  of_<< std::endl << "#@@@@@ Append from here..."<<std::endl;
	  of_<< std::endl << " @@@@@  " << show_time() << std::endl;
	}
	else
	  of_.open(outfile.c_str(), std::ios_base::out);
      }
      else{
	of_.open(outfile.c_str(), std::ios_base::out);
      }

      if(outfile.size() != 0){
	std::cout << "REPAIRER 2D MOL: record variable violations in file '"<<outfile << "'."<< std::endl;
      }
      
      //! these parameters must not be renamed or left out in the corresponding
      //! input files!
      nx_ = PD.get_datum<int>("Nx");
      ny_ = PD.get_datum<int>("Ny");
      npt_ = nx_*ny_;
      rhoSmall_ = PD.get_datum<K>("rho_min");
      v1MaxAbs_ = PD.get_datum<K>("v1_max");
      v2MaxAbs_ = PD.get_datum<K>("v2_max");
      Tmin_ = PD.get_datum<K>("Tlow");
      Tmax_ = PD.get_datum<K>("Tup");
      //SmallSpec_ = PD.get_datum<K>("smallSpec");
      approxEqual_ = PD.get_datum<BaseType>("approx_equal");

      adonis_assert((v1MaxAbs_ > 0) && (v2MaxAbs_ > 0));
		
      msg2screen_ = addInfo_ = std::string(); //reset
      UnrollLoop<0,20>::assign_value(ViolatedVars_,false);

      dim_[0] = nx_; dim_[1] = ny_; dim_[2] = 0;
    }

	
    void info(int timeStep = 0, int NewtonIter = -1, const std::string& str = std::string()){
      bool vioDetected(false); //no violation detected by default
      (*this).check4Violation(0,"rho < rhoSmall = "+Num2str(rhoSmall_)+". ",vioDetected);
      (*this).check4Violation(1,"rho is NaN. ",vioDetected);
      (*this).check4Violation(2,"rho is Inf. ",vioDetected);
      (*this).check4Violation(3,"|v_1| > |v1_max|. ",vioDetected);
      (*this).check4Violation(4,"v1 is NaN. ",vioDetected);
      (*this).check4Violation(5,"|v_2| > |v2_max|. ",vioDetected);
      (*this).check4Violation(6,"v2 is NaN. ",vioDetected);
      (*this).check4Violation(7,"T < Tmin. ",vioDetected);
      (*this).check4Violation(8,"T > Tmax. ",vioDetected);
      (*this).check4Violation(9,"T is NaN. ",vioDetected);
      (*this).check4Violation(10,"Some chemical species < 0. ",vioDetected);
      (*this).check4Violation(11,(("Some chemical species > 1. \n")+addInfo_),vioDetected);
      (*this).check4Violation(12,"Some chemical species is NaN. ",vioDetected);
      (*this).check4Violation(13,"Failed to perserve mass. ",vioDetected);
      (*this).check4Violation(14,"After mass preservation: some species < 0.0. ",vioDetected);
      (*this).check4Violation(15,"Mass balance: sm ~ 0 in (1-sm). ",vioDetected);
     
      //insert something at the very beginning of the message string		
      if(vioDetected){
	//additional information provided
	if(str.size() != 0){
	  std::cout << "#### "<<str << std::endl;
	}

	if(NewtonIter >= 0){
	  tmpStr_ = ", in Newton iteration "+Num2str(NewtonIter)+" (corrections already performed): \n ";  
	}
	else{
	  tmpStr_ = " "; 
	}
	msg2screen_.insert(0,"VARIABLES OUT OF BOUNDS in time step "+Num2str(timeStep)+" "+tmpStr_);
	of_ << //("VARIABLES OUT OF BOUNDS in time step "+Num2str(timeStep)+tmpStr_) << 
	  msg2screen_ << std::endl;  //write violations also to a default file
	FM_.nice_output(msg2screen_,35);  
      }
      //reset 
      msg2screen_ = addInfo_ = std::string(); 
      UnrollLoop<0,20>::assign_value(ViolatedVars_,false);
      
    }	
	
    //! try to impose reasonable values for all computational quantities
    template<class V>
    void repair(V& u, int excessSpecIndex){

      violation_ = false; //reset; no violation at all 

      if(numOfSpec_ == 0){ //compute #species only once
	numOfSpec_ = ((int)u.size())/npt_ - 4; //1=rho,2=v1,3=v2,T=4
	std::cout << "Repairer: numOfSpec_ = "<< numOfSpec_ << std::endl;
	// resize Ivec_ only once
	Ivec_.resize(numOfSpec_);
      }
      adonis_assert((excessSpecIndex > 0) && (excessSpecIndex < numOfSpec_));
      
      msg2screen_ = addInfo_ = std::string(); //reset to empty string
      UnrollLoop<0,20>::assign_value(ViolatedVars_,false); //reset; no violation
      

      //! It turned out that a treatment of boundary values is also reasonable,
      //! in case these are subject to changes such as Y_k and rho
      for(int i = 0; i < nx_; ++i){     //formerly: 1...nx_-1 
	for(int j = 0; j < ny_; ++j){   //formerly: 1...ny_-1
	  
	  //RHO
	  if(FlowVariables<V>::primitive(i,j,u,0,nx_,ny_) < rhoSmall_){ //<0.0
	    FlowVariables<V>::assign(i,j,u,0,nx_,ny_,rhoSmall_);
	    ViolatedVars_[0] = true;
	  }

	  if(IsNan(FlowVariables<V>::primitive(i,j,u,0,nx_,ny_))){
	    FlowVariables<V>::assign(i,j,u,0,nx_,ny_,rhoSmall_);
	    ViolatedVars_[1] = true;
	  }
	  //...if rho is a nan, make it something more reasonable
	  if(IsInf(FlowVariables<V>::primitive(i,j,u,0,nx_,ny_))){ //Wbar value is just an arguable value here
	    FlowVariables<V>::assign(i,j,u,0,nx_,ny_, (101325*1.75)/(PhysicalConstants<typename TypeAdapter<typename V::value_type>::BaseType>::Rgas*Tmax_));
	    ViolatedVars_[2] = true;
	  }


	  //V
	    if(Abs(FlowVariables<V>::primitive(i,j,u,1,nx_,ny_)) > v1MaxAbs_){
	      //! take signum(current v1)*v1MaxAbs_
	      FlowVariables<V>::assign(i,j,u,1,nx_,ny_,Sgn(FlowVariables<V>::primitive(i,j,u,1,nx_,ny_))*v1MaxAbs_);
	    ViolatedVars_[3] = true;
	  }
	 
	  if(IsNan(FlowVariables<V>::primitive(i,j,u,1,nx_,ny_))){
	    FlowVariables<V>::assign(i,j,u,1,nx_,ny_,0.0);
	    ViolatedVars_[4] = true;
	  }

	  
	  if(Abs(FlowVariables<V>::primitive(i,j,u,2,nx_,ny_)) > v2MaxAbs_){
	    //! take signum(current v2)*v2MaxAbs_
	    FlowVariables<V>::assign(i,j,u,2,nx_,ny_,Sgn(FlowVariables<V>::primitive(i,j,u,2,nx_,ny_))*v2MaxAbs_);
	    ViolatedVars_[5] = true;
	  }
	  
	  if(IsNan(FlowVariables<V>::primitive(i,j,u,2,nx_,ny_))){
	    FlowVariables<V>::assign(i,j,u,2,nx_,ny_,0.0);
	    ViolatedVars_[6] = true;
	  }


	  //T
	  if(FlowVariables<V>::primitive(i,j,u,3,nx_,ny_) < Tmin_){
	    FlowVariables<V>::assign(i,j,u,3,nx_,ny_,Tmin_);
	    ViolatedVars_[7] = true;
	    violation_ |= true;
	  }
	  
	  if(FlowVariables<V>::primitive(i,j,u,3,nx_,ny_) > Tmax_){
	    FlowVariables<V>::assign(i,j,u,3,nx_,ny_,Tmax_);
	    ViolatedVars_[8] = true;
	    violation_ |= true;
	  }
	  
	  if(IsNan(FlowVariables<V>::primitive(i,j,u,3,nx_,ny_))){
	    FlowVariables<V>::assign(i,j,u,3,nx_,ny_,Tmin_);
	    ViolatedVars_[9] = true;
	    violation_ |= true;
	  }

	  //=============================================================
	  //smart balance to unity with species mass fractions
	  //SPECIES
	  // 0.0 <= Y_k <= 1.0, sum_k Y_k = 1
	  count0_ = count1_ = countLarge_ = countOk_ = iterCount_ =  0; //reset
	  sum_ = partsum_ = 0;  
	  for(int k = 0; k < numOfSpec_; ++k){
	    if(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) < 0.0){
	      FlowVariables<V>::assign(i,j,u,4+k,nx_,ny_,0.0); //clip to 0
	      count0_++;
	      Ivec_[k] = -1;
	      ViolatedVars_[10] = true;
	      violation_ |= true;
	    }
	    else if (FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) >= 1.0){
	      addInfo_ += "("+Num2str(i)+", "+Num2str(j)+"): Y_["+Num2str(k)+"] = "+ Num2str(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_))+"   §§  ";
	      
	      //! uncomment this if you want a waring also on screen
	      //ThrowErrorOrNot<Constant<K>::setErrorAndStop,MassBalanceError>::info(addInfo_);
	      FlowVariables<V>::assign(i,j,u,4+k,nx_,ny_,1.0); //clip to 1
	      count1_++;
	      Ivec_[k] = 2;
	      ViolatedVars_[11] = true;
	      violation_ |= true;
	    }
	    else if (IsNan(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_))){
	      FlowVariables<V>::assign(i,j,u,4+k,nx_,ny_,0.0); //clip to 0
	      count0_++;
	      Ivec_[k] = -1;
	      ViolatedVars_[12] = true;
	      violation_ |= true;
	    }
	    else if((*this).is_large_value(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_))){           //large value is just a definition
	      countLarge_++;
	      Ivec_[k] = 1; 
	    }
	    else if((FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) >= 0.) && (FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) < 1.)){  //values ok
	      Ivec_[k] = 0;
	      countOk_++;
	    }
	    else {
	      Ivec_[k] = 0; //ok by default
	    }
	  
	    sum_ +=  FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_);
	  }
	  //std::cout << "sum_ = "<< sum_ << std::endl;

	  preserve_mass(excessSpecIndex,i,j,u); 
	} //end for j
      }   //end for i

    }

    bool are_phys_quantities_violated() const {return violation_;}

    const StaticArray<int,3>& space_dimensions() const {return dim_;}
    
  private:
    int nx_, ny_, npt_, numOfSpec_;
    K rhoSmall_, v1MaxAbs_, v2MaxAbs_, Tmin_, Tmax_, SmallSpec_;
    K sm_;
    bool violation_;
    int maxk_, count0_, count1_, countLarge_, countOk_, iterCount_;
    K sum_, partsum_,tmp_, mean_, approxEqual_, avg_;
    InfoVecType Ivec_;
    ExprTmpl::MyVec<K> uGood_;

    StaticArray<int,3> dim_;

    //! private functions
    template<class V>
    bool does_Y_k_exceed_1(const V& u){
      bool res(false);
      for(int i = 0; i < nx_; ++i){
	for(int j = 0; j < ny_; ++j){
	  for(int k = 0; k < numOfSpec_; ++k){
	    if(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) >= 1.0){
	      res |= true;
	      ViolatedVars_[11] = true;
	      violation_ |= true;
	      //break; //leave loop
	    }
	    else if(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) < 0.0){
	      res |= false; //not greater than 1
	      ViolatedVars_[10] = true;
	      violation_ |= true;
	    }
	    else{
	      res |= false; //not greater than 1
	    }
	    
	  }
	}
      }
      return res;
    }
    
    //       0     1      2      3        4      5       6    7    8     9   10  11  12     13                 14     15  
  
    bool ViolatedVars_[20]; //elements in this order:{rho<0,rhoNan,rhoInf,v1MaxAbs,v1Nan,v2MaxAbs,v2Nan,Tmin,Tmax,Tnan,Y<0,Y>1,Y_NaN, preserve_mass fault,excess spec < 0,sm < 0}
    std::string msg2screen_, addInfo_, tmpStr_;
    FancyMessages FM_;
    std::ofstream of_;
    

    int OFF(int i, int j, int spec) const{
      return (i + nx_*j + spec*npt_);
    }

    void check4Violation(unsigned i, const std::string& msg, bool& vioDetected){
      if(ViolatedVars_[i] == true){
	msg2screen_ += msg;
	vioDetected = true;
      }
    }
	
    bool is_large_value(const value_type& val, const BaseType& big = 0.8){
      if((val >= big) && (val < 1.0))
	return true;
      else
	return false;
    }

    //! TODO · TODO · TODO · TODO 
    //! assumes that u contains corrected values 
    //! excess corresponds to index of excess species
    template<class V>
    bool preserve_mass(int excess, int i, int j, V& u){
     bool isOK = is_approximately_equal(sum_,1.,approxEqual_);
     iterCount_ = 0; //rest
     partsum_ = 0;
     int m(0);

     while(!isOK){
       isOK = is_approximately_equal(sum_,1.,approxEqual_);
       iterCount_++;

       m = 0; //reset
       partsum_ = sum_ = mean_ = 0;
       maxk_ = 0;
      
       //! mass fraction of first species, i.e. Y[0]
       tmp_ = FlowVariables<V>::primitive(i,j,u,4+0,nx_,ny_); //for max determination
       for(int k = 0; k < numOfSpec_; ++k){
	 if(Ivec_[k] == 2){
	   if(k != excess){
	     //take nearby value. Otherwise take random number between 0 and 1
	     FlowVariables<V>::assign(i,j,u,4+k,nx_,ny_, nearby_value_2D(avg_,u,i,j,nx_,ny_,k));
	     count1_--; //decrease by one
	   }
	 }
	 if((FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) >= 0) && (!is_approximately_equal(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_),1.,approxEqual_))){
	   Ivec_[k] = 0;
	   countOk_++;
	 }
	 if((*this).is_large_value(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_))){
	   Ivec_[k] = 1;
	 }
	 if (is_approximately_equal(FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_),1.,approxEqual_)){
	   Ivec_[k] = 2;
	  }
	
	
	 if(Ivec_[k] != 2){
	   if(k != excess){ //excess species will be recovered later on
	     partsum_ += FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_);
	     m++;
	   }
	 }

	 if(m != 0){
	   mean_ = partsum_/m;  //mean value
	 }
	 else{
	   mean_ = 0;
	 }

	 sum_ += FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_);  //just register sum as well

	 //! find max value that is not equal to 1
	 if((FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) > tmp_) && (Ivec_[k] != 2)){
	   tmp_ = FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_);
	   maxk_ = k;
	 }
	 
       }//end for

       //std::cout << iterCount_<<".)    maxk_ = "<< maxk_ << "." << std::endl;

       //
       if(partsum_ >= 1){
	 partsum_ = 0; //Reset again
	 for(int k = 0; k < numOfSpec_; ++k){
	   if((Ivec_[k] == 1) || (k == maxk_)){  //if large value encountered
	     if(iterCount_ < 10){
	       FlowVariables<V>::assign(i,j,u,4+k,nx_,ny_,0.9*FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_)); //reduce value by 10%
	     }
	     else if((iterCount_ >=10) && (iterCount_ < 25)){
	       FlowVariables<V>::assign(i,j,u,4+k,nx_,ny_,0.6*FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_)); // reduce value by 40% now
	     }
	     else{
	        FlowVariables<V>::assign(i,j,u,4+k,nx_,ny_,0.4*FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_)); // reduce value by 60% now
	     }
	   }
	   if(k != excess){
	     partsum_ += FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_);
	   }
	 }
       }  //don't use if-else here
       if(partsum_ < 1){
	 //std::cout << "partsum_ = "<< partsum_;
	 FlowVariables<V>::assign(i,j,u,4+excess,nx_,ny_,1. - partsum_);
	 sum_ = 0;
	 for(int k = 0; k < numOfSpec_; ++k)
	   sum_ += FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_);
	 
       }

      
       // if(iterCount_ >= 10){
       // 	std::cout << "#### INFO: Too many iterations required..."<<std::endl;
       // }

       if(iterCount_ >= 50){
	 std::cout << "#### WARNING: Too many iterations required..."<<std::endl;
	 isOK = false;
	 break;
       }
      
     } //end while

     if(isOK==false){
       std::cout << "****CAUTION: Smart balance to unity: Sum Y_k("<<i<<","<<j<<") = "<< sum_ << " does not add to 1 yet (partsum_ = "<< partsum_ << ")." << std::endl;
       std::cout << "Y_{"<<i<<","<<j<<"} = ";
       for(int k = 0; k < numOfSpec_; ++k){
	 std::cout << FlowVariables<V>::primitive(i,j,u,4+k,nx_,ny_) << ", ";
       }
       std::cout << std::endl;
     }
    
     return isOK;
    }

  };


  /**
   * \brief basically the same procedure as above ecept that it directly applies
   * to container Yfrac
   */
  template<class V>
  class CorrectFractions{
  public:

    typedef typename V::value_type value_type;
    typedef typename TypeAdapter<value_type>::BaseType BaseType;
    typedef std::vector<int> InfoVecType;
    
    CorrectFractions(int nspec = 0):numOfSpec_(nspec),Ivec_(nspec), sum_(0), partsum_(0),tmp_(0), mean_(0),maxk_(0) ,count0_(0),count1_(0),countLarge_(0), countOk_(0), iterCount_(0){}


    void initialize(int nspec){
      numOfSpec_ = nspec;
      Ivec_.resize(nspec);
    }

    void smart_balance_to_unity(V& Y, int excess, const BaseType& eps = 1.e-09){
      adonis_assert(excess < numOfSpec_);
      adonis_assert((int)Y.size() == numOfSpec_);
 
      count0_ = count1_ = countLarge_ = countOk_ = iterCount_ =  0; //reset
      sum_ = partsum_ = 0;
      for(int k = 0; k < numOfSpec_; ++k){
	if(Y[k] >= 1.0){
	  count1_++;
	  Y[k] = 1.;
	  Ivec_[k] = 2; 
	}
	else if(Y[k] < 0.){
	  count0_++;
	  Y[k] = 0.;
	  Ivec_[k] = -1; 
	}  //large value is just a definition
	else if((*this).is_large_value(Y[k])){
	  countLarge_++;
	  Ivec_[k] = 1;
	}
	else if((Y[k] >= 0.) && (Y[k] < 1.)){ //values ok
	  Ivec_[k] = 0;
	  countOk_++;
	}
	else{
	  Ivec_[k] = 0; //ok
	}

     
	sum_ += Y[k];
      }

      //std::cout << "sum_ = "<< sum_ << std::endl;

      //! now there all values are from [0,1]
      //! Check sum. If it happens that sum is already ~1 then stop
      bool isOK = is_approximately_equal(sum_,1.,eps);
      iterCount_ = 0; //rest
      partsum_ = 0;
      int m(0);
      while(!isOK){
	isOK = is_approximately_equal(sum_,1.,eps);
	iterCount_++;

	m = 0; //reset
	partsum_ = sum_ = mean_ = 0;
	maxk_ = 0;
      
	tmp_ = Y[0]; //for max determination
	for(int k = 0; k < numOfSpec_; ++k){
	  if(Ivec_[k] == 2){
	    if(k != excess){
	      //TODO: nearby value (whatever this is) or just take random number between 0 and 1. 
	      Y[k] = (value_type)rand()/(value_type)(RAND_MAX-1); //from [0,1)
	      count1_--; //decrease by one
	    }
	  }
	  if((Y[k] >= 0) && (!is_approximately_equal(Y[k],1.,eps))){ //(Y[k] < 1)){
	    Ivec_[k] = 0;
	    countOk_++;
	  }
	  if((*this).is_large_value(Y[k])){
	    Ivec_[k] = 1;
	  }
	  if (is_approximately_equal(Y[k],1.,eps)){
	    Ivec_[k] = 2;
	  }
	
	
	  if(Ivec_[k] != 2){
	    if(k != excess){ //excess species will be recovered later on
	      partsum_ += Y[k];
	      m++;
	    }
	  }

	  if(m != 0){
	    mean_ = partsum_/m;  //mean value
	  }
	  else{
	    mean_ = 0;
	  }

	  sum_ += Y[k];  //just register sum as well

	  //! find max value that is not equal to 1
	  if((Y[k] > tmp_) && (Ivec_[k] != 2)){
	    tmp_ = Y[k];
	    maxk_ = k;
	  }
	 
	}//end for

	//std::cout << iterCount_<<".)    maxk_ = "<< maxk_ << "." << std::endl;
	//(*this).print("Y");

	//
	if(partsum_ >= 1){
	  partsum_ = 0; //Reset again
	  for(int k = 0; k < numOfSpec_; ++k){
	    if((Ivec_[k] == 1) || (k == maxk_)){  //if large value encountered
	      if(iterCount_ < 10)
		Y[k] *= 0.9;   //reduce value by 10%
	      else  if((iterCount_ >=10) && (iterCount_ < 25))
		Y[k] *= 0.6; // reduce value by 40% now
	      else
		Y[k] *= 0.4;
	    }
	    if(k != excess){
	      partsum_ += Y[k];
	    }
	  }
	
	}  //don't use if-else here
	if(partsum_ < 1){
	  //std::cout << "partsum_ = "<< partsum_;
	  Y[excess] = 1. - partsum_;
	  //std::cout << ",   Y[excess] = "<< Y[excess] << "   sum = partsum_+Y[excess] = "<< std::setprecision(12) << (Y[excess] + partsum_) << std::endl;
	  //! just for output
	  sum_ = 0;
	  for(int k = 0; k < numOfSpec_; ++k)
	    sum_ += Y[k];
	}

      
	// if(iterCount_ >= 10){
	// 	std::cout << "#### INFO: Too many iterations required..."<<std::endl;
	// }

	if(iterCount_ >= 50){
	  std::cout << "#### WARNING: Too many iterations required..."<<std::endl;
	  isOK = false;
	  break;
	}
      
      } //end while

      // if(isOK){
      // 	std::cout << ":) SUCCESS: sum = "<< sum_ << std::endl;
      // }
      // else{
      // 	std::cout << ":( FAILURE: sum = "<< sum_ << std::endl;
      // }
      if(!isOK){
	std::cout << ":( FAILURE: sum = "<< sum_ << std::endl;
      }
    }
    
    


    const InfoVecType& get_value_info() const{ return Ivec_;}
    const int index_of_max_value() const{return maxk_;}
    const int num_of_iterations() const {return iterCount_;}
    
    
  private:
    int numOfSpec_;
    InfoVecType Ivec_;
    value_type sum_, partsum_, tmp_, mean_;
    int maxk_,count0_, count1_, countLarge_, countOk_, iterCount_;


    bool is_large_value(const value_type& val, const BaseType& big = 0.8){
    if((val >= big) && (val < 1.0))
      return true;
    else
      return false;
    }
  };

} //end namespace 

#endif
