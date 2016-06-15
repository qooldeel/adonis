#ifndef READ_IN_DATA_GIVEN_AS_A_DENSE_MATRIX_HH
#define READ_IN_DATA_GIVEN_AS_A_DENSE_MATRIX_HH

#include "../expressiontemplates/exprvec.hh"
#include "../common/adonisassert.hh"
#include "../common/error.hh"
#include "../common/typeadapter.hh"
#include "../marcvecmatrix/myfunctions.h"
#include "../fdm/gensettings.hh"
#include "../graphics/gnuplot.hh"

namespace Adonis{


  template<class T>
  class ReadInMatrix{
  public:
    typedef std::string StringType;
    typedef T value_type;
    typedef typename TypeAdapter<T>::BaseType BaseType;
    typedef ExprTmpl::MyVec<T> VType;
    typedef std::size_t SizeType;


    ReadInMatrix(const BaseType& eps = 1.e-08, int prec = 12):eps_(eps),prec_(prec),rows_(0),cols_(0),lc_(0),count_(0),isEmpty_(true),offs_(
#ifndef GHOST_POINTS_INCLUDED
																	    0
#else //ghost points included
																	    1
#endif
																	    ), idx_(-1)

    {}

  
    void read_from_file(const StringType& InfileName){
      if(isEmpty_ == false){
	ADONIS_ERROR(ClearError,"Before reading in again, you must invoke member function 'clear' in order to avoid possible spurious results");
      }
      inName_ = InfileName;
      std::ifstream infile(InfileName.c_str(),std::ios_base::in);
      rows_ = cols_ = lc_ = 0;  //reset before reading in
      int ncols(0);
      bool flag(false);   //if a line is meaningful then increase rows_
      lc_ = 0;          //line counter (every line)
      SizeType idx(0);
      if(!infile.is_open()){
	ADONIS_ERROR(FileError,"File \""<<InfileName<<"\" not found");
      }

      while(!infile.eof()){ //if end of file hasn't reached yet
	ncols = 0;         //set back
	line_.clear();
	ds_.clear();        //set back
	flag = false;       //reset

	// read line to end. String 'line' contains now everything except '\n'
	getline(infile,line_,'\n'); 
	//increase line number after each read in line. This corresponds to the 
	//exact line number of the opened file in your editor
	lc_++;                           
	//std::cout << lc_ << ".)   line_.size() = "<< line_.size() << std::endl;
	
	for(SizeType k = 0; k < line_.size(); ++k){ //iterate over read in line

	  //here we read a meaningful content of string 'line'
	  if(!my_function_collection::is_delimiter(line_[k]) && !my_function_collection::is_comment(line_[k])){
	    ds_.push_back(line_[k]);
	  }

	  //ds contains some content and the current char is a white char
	  if( (ds_.size() != 0) && ((my_function_collection::is_delimiter(line_[k])) || (my_function_collection::is_comment(line_[k])) || (k == line_.size()-1) ) ){
	    //data_.push_back(atof(ds_.c_str()));
	    data_.push_back(OutputTypeAdapter<T>::convert(ds_));
	    //std::cout << data_[idx++] << " ";
	    ncols++;  
	    ds_.clear();
	    flag = true;
	  }
					
	  //rest of line is comment so break here
	  if(my_function_collection::is_comment(line_[k])){
	    break;
	  }
	} //end for-loop
				
	if(flag==true){
	  rows_++;    //increase row number by one
	  if(rows_ == 1){
	    cols_ = ncols;  //only assign once		 
	  }
	  //check if number of columns stays equal in each meaningful row. 
	  //If not, no valid matrix format has been detected.
	  if(rows_ > 1){
	    if(cols_ != ncols){ 
	      ADONIS_ERROR(IODimensionError, "NO VALID MATRIX FORMAT: in row "<< rows_ << " (line no. "<<lc_<<" of file '" << InfileName <<"')," << std::endl << " I found "<< ncols << " columns compared to "<< cols_ << " in previous rows.");
	    }
	  }
	}

      }  //end of while not eof
      isEmpty_= false;
    }

    const int rows() const {return rows_;}
    const int cols() const {return cols_;}
    const bool is_empty() const {return isEmpty_;}

    //set back everything. This can be used prior to reading in a new file
    void clear() {
      rows_ = cols_ = lc_ = 0;
      line_.clear();
      ds_.clear();
      data_.clear();
      isEmpty_ = true;
    }

    const int number_of_read_lines() const {return lc_;}


    /**
     * \brief Print what has been read from file and stored 
     */
    void print_matrix(const StringType& s = StringType()) const{
      if(s.size() != 0)
	std::cout <<  s << " = "<< std::endl;
      for(int i = 0; i < rows_; ++i){
	for(int j = 0; j < cols_; ++j){
	  std::cout << data_[i*cols_+j] << "  ";
	}
	std::cout << std::endl;
      }
      //JUST for debugging
      // if(s.size() != 0)
      // 	std::cout <<  s << " = "<< std::endl;
      // for(int i = 0; i < rows_; ++i){
      // 	for(int j = 0; j < cols_; ++j){
      // 	  //std::cout << data_[i*cols_+j] << "  ";
      // 	  if(data_[i*cols_+j] > 3200 && data_[i*cols_+j] < 8000){
      // 	    std::cout << "<"<<i<<","<<j<<"> = "<< data_[i*cols_+j] << "; ";
      // 	    break;
      // 	  }
      // 	}
      // 	//std::cout << std::endl;
      // }
    }

    //! write data_ to file without any whitelines inbetween, i.e. write
    //! matrix to file
    void write_2_file(const StringType& OutfileName, const StringType& furtherComment = "##"){
      if(data_.size() != 0){
	std::ofstream ofile(OutfileName.c_str(),std::ios_base::out);
	ofile << "## Printed data_ created from file '"<< inName_ << "'"<<std::endl;
	ofile << "#"<<furtherComment << std::endl;
	for(int i = 0; i < rows_; ++i){
	  for(int j = 0; j < cols_; ++j){
	    ofile << std::setprecision(prec_) << data_[i*cols_+j] << "  ";
	  }
	  ofile << std::endl;
	}
      }
      else{
	ADONIS_INFO(Information,"data_ field is empty. Apparently, nothing has been read in or something has been cleared");
      }
    }

    T& operator()(int i, int j) {
      adonis_assert((i>=0) && (i<rows_) && (j>=0) && (j<cols_));
      return data_[i*cols_+j];
    }

    const T& operator()(int i, int j) const{
      adonis_assert((i>=0) && (i<rows_) && (j>=0) && (j<cols_));
      return data_[i*cols_+j];
    }

    const int get_offset() const {return offs_;}

    //! assign read in 2D MOL slide 
    //! \param v vector to be overwritten with read in result.Must be of size
    //!        Nx·Ny·numberOfQuantities if withPressure=true otherwise of size
    //!        Nx·Ny·(numberOfQuantities-1) 
    //! \param Nx number of grid points in x-direction
    //! \param Ny number of grid points in y-direction
    //! \param numberOfQuantities total # of physical quantities (including p)
    //! \param withPressure consider pressure (true) or not (false)
    template<class V>
    void assign_2_2D_MOL_slide(V& v, int Nx, int Ny, int numberOfQuantities, bool withPressure = true){ 
      if(withPressure){
	adonis_assert((int)v.size() == Nx*Ny*numberOfQuantities);
      }
      else{  //pressure not considered
	adonis_assert((int)v.size() == Nx*Ny*(numberOfQuantities-1));
	std::cout << "***** Pressure (last column of 2D MOL file) is omitted." << std::endl;
      }
      int npt = Nx*Ny;
#ifndef GHOST_POINTS_INCLUDED
      for(int j = 0; j < Ny; ++j)
#else
       for(int j = 1; j < Ny-1; ++j)	
#endif
	{
#ifndef GHOST_POINTS_INCLUDED
	for(int i = 0; i < Nx; ++i)
#else
	for(int i = 1; i < Nx-1; ++i)  
#endif
	  {
	  for(int k = 0; k < numberOfQuantities; ++k){
	    if(withPressure){
	      v[i + Nx*j + npt*k] = 
#ifdef NONCONSERVATIVE_FORM
		(*this).get_2D_MOL_entry((i-offs_),(j-offs_),Ny,k,numberOfQuantities);
#else //conservative form
	      (((k == 0) || (k == numberOfQuantities-1)) ? ((*this).get_2D_MOL_entry((i-offs_),(j-offs_),Ny,0,numberOfQuantities)) : ((*this).get_2D_MOL_entry((i-offs_),(j-offs_),Ny,k,numberOfQuantities)*(*this).get_2D_MOL_entry((i-offs_),(j-offs_),Ny,0,numberOfQuantities))); //multiplication by rho_i,j (except when rho or pressure is considered)
#endif

	    }
	    else{ //exclude pressure
	      if(k < numberOfQuantities-1){
		v[i + Nx*j + npt*k] = 
#ifdef NONCONSERVATIVE_FORM
		  (*this).get_2D_MOL_entry((i-offs_),(j-offs_),Ny,k,numberOfQuantities);
#else //conservative form
		((k == 0) ? ((*this).get_2D_MOL_entry((i-offs_),(j-offs_),Ny,k,numberOfQuantities)) : ((*this).get_2D_MOL_entry((i-offs_),(j-offs_),Ny,k,numberOfQuantities)*(*this).get_2D_MOL_entry((i-offs_),(j-offs_),Ny,0,numberOfQuantities))); //multiplication by rho_i,j (except when rho is considered)
#endif
	      }
	    }
	      
	  }
	}
      }
      std::cout << "Assignment DONE!" << std::endl;
    }

    //! extract slide in x-direction
    void fixed_x_2D(const BaseType& xval, const StringType& outfile){
      std::ofstream of(outfile.c_str(),std::ios_base::out);
      count_ = 0;
      for(int i = 0; i < rows_; ++i){
	if((Abs(xval-(*this)(i,0)) <= eps_) //|| (is_equal(xval,(BaseType)round(xval)))
	   ){
	  count_++;
	  for(int j = 0; j < cols_; ++j){
	    of << std::setprecision(prec_)<< (*this)(i,j) << " ";
	  }
	  of << std::endl;
	}
      }
      if(count_ == 0){
	ADONIS_WARNING(Warning,"count_ = 0. Are you sure xval = "<<xval<<" is the right value?...\n Are you sure eps_ = "<< eps_ << " has been chosen appropriately?");
      }
      of.close();
    }

    //! extract slide in y-direction
    void fixed_y_2D(const BaseType& yval, const StringType& outfile){
      std::ofstream of(outfile.c_str(),std::ios_base::out);
      count_ = 0;
      for(int i = 0; i < rows_; ++i){
	if((Abs(yval-(*this)(i,1)) <= eps_)// || (is_equal(yval,(BaseType)round(yval)))
	   ){
	  count_++;
	  for(int j = 0; j < cols_; ++j){
	    of << std::setprecision(prec_)<< (*this)(i,j) << " ";
	  }
	  of << std::endl;
	}
      }
      if(count_ == 0){
	ADONIS_WARNING(Warning,"count_ = 0. Are you sure yval = "<<yval<<" is the right value?...\n Are you sure eps_ = "<< eps_ << " has been chosen appropriately?");
      }
      of.close();
    }

    //! shouldn't be invoked too much due to too many local variables
    //! NOTE: the printed slide is always WITHOUT ghost points 
    void show_2D_MOL_boundaries(int Nx, int Ny, const BaseType& hx, const BaseType& hy, int k, int numberOfQuantities) const{
      adonis_assert(isEmpty_ != true); //file must be read in
      std::string leftBdyOut = "left_bdy_2_plot.txt",
	upBdyOut = "up_bdy_2_plot.txt", 
	downBdyOut = "down_bdy_2_plot.txt",
	rightBdyOut = "right_bdy_2_plot.txt";

      std::ofstream left(leftBdyOut.c_str(), std::ios_base::out);
      std::ofstream up(upBdyOut.c_str(), std::ios_base::out);
      std::ofstream down(downBdyOut.c_str(), std::ios_base::out);
      std::ofstream right(rightBdyOut.c_str(), std::ios_base::out);

      int nx =
#ifndef GHOST_POINTS_INCLUDED
	Nx
#else
	(Nx-2)
#endif
	;
     int ny =
#ifndef GHOST_POINTS_INCLUDED
	Ny
#else
	(Ny-2)
#endif
	;     

     //! read in file without ghost points. nx and ny calculated correctly, hence
     //! no special indexing required any more
     BaseType xaxis(0);
     //LEFT + RIGHT

      for(int j = 1; j < ny-1; ++j){
	xaxis = (j*hy);
	left << std::setprecision(prec_) << xaxis << "  "<<
	  //data_[0*((2+numberOfQuantities)*ny) + (2+numberOfQuantities)*j + (2+k)]
	  (*this). get_2D_MOL_entry_no_ghost_points(0,j,ny,k,numberOfQuantities)
	     << std::endl;

	right << std::setprecision(prec_) << xaxis << "  "<<
	  //data_[(ny-1)*((2+numberOfQuantities)*ny) + (2+numberOfQuantities)*j + (2+k)]
	  (*this). get_2D_MOL_entry_no_ghost_points(nx-1,j,ny,k,numberOfQuantities)
	      << std::endl;
      }
 
	  
      //UP + DOWN
      for(int i = 0; i < nx; ++i){
	xaxis = (i*hx);
	up << std::setprecision(prec_) << xaxis << "  "<<
	  //data_[i*((2+numberOfQuantities)*ny) + (2+numberOfQuantities)*(ny-1) + (2+k)]
	  (*this). get_2D_MOL_entry_no_ghost_points(i,ny-1,ny,k,numberOfQuantities)
	   << std::endl;

	down <<  std::setprecision(prec_) << xaxis << "  "<<
	  //data_[i*((2+numberOfQuantities)*ny) + (2+numberOfQuantities)*0 + (2+k)]
	  (*this). get_2D_MOL_entry_no_ghost_points(i,0,ny,k,numberOfQuantities)
	     << std::endl;

      }

      //! start plotting
      GNUplot p1, p2, p3, p4;

      std::string intro = "set key box right; set autoscale fix; set ylabel \"quantity_"+Num2str(k)+"\"; ";
      p1(intro+"set xlabel \"y\"; plot \""+leftBdyOut+"\" using 1:2 with lines lw 3 title \"left bdy\"");
      p2(intro+"set xlabel \"y\"; plot \""+rightBdyOut+"\" using 1:2 with lines lw 3 title \"right bdy\"");
       p3(intro+"set xlabel \"x\"; plot \""+upBdyOut+"\" using 1:2 with lines lw 3 title \"up bdy\"");
      p4(intro+"set xlabel \"x\"; plot \""+downBdyOut+"\" using 1:2 with lines lw 3 title \"down bdy\"");
      
    }

    //! NOTE again: file is always printed without ghost points
    void show_min_max_value_2D_MOL(int Nx, int Ny, int k, int numberOfQuantities){
      T valmx, valmn, valprev;
      valmx = valmn = (*this).get_2D_MOL_entry_no_ghost_points(0,0,Ny,k,numberOfQuantities);
      int nx =
#ifndef GHOST_POINTS_INCLUDED
	Nx
#else
	(Nx-2)
#endif
	;
     int ny =
#ifndef GHOST_POINTS_INCLUDED
	Ny
#else
	(Ny-2)
#endif
	;     

      int imx(0), jmx(0), imn(0), jmn(0);
      for(int i = 0; i < nx; ++i){
	for(int j = 0; j < ny; ++j){
	  //get max
	  valprev = valmx;
	  valmx = Max(valmx,(*this).get_2D_MOL_entry_no_ghost_points(i,j,ny,k,numberOfQuantities));
	  if(valmx > valprev){
	    imx = i; jmx = j;
	  }
	  //get min
	  valprev = valmn;
	  valmn = Min(valmn,(*this).get_2D_MOL_entry_no_ghost_points(i,j,ny,k,numberOfQuantities));
	  if(valmn < valprev){
	    imn = i; jmn = j;
	  }
	}
      }
   
      std::cout << std::endl<< "quant_MAX_"<< k << "("<<imx<<", "<<jmx<<") = "<< valmx << std::endl;
      std::cout << "quant_MIN_"<< k << "("<<imn<<", "<<jmn<<") = "<< valmn << std::endl<<std::endl;
    } 
    
    const int count() const {return count_;}
    int count() {return count_;}
    
  private:
    BaseType eps_;
    int prec_, rows_, cols_, lc_, count_;
    mutable bool isEmpty_;
    int offs_;
    mutable int idx_; //can be changed in const-member
    VType data_;
    StringType line_, ds_, inName_;
    

    //forbid to invoke this member directly in order to prevent user of extracting a wrong order
    //use member function <I> assign </I> instead.
    //first two columns contain x- and y-axis values
    const T& get_2D_MOL_entry(int i, int j, int Ny, int k, int numberOfQuantities) const{ //without x and y
      adonis_assert((j >= 0) && (j < Ny) );
      idx_ = i*((2+numberOfQuantities)*
#ifndef GHOST_POINTS_INCLUDED
	Ny
#else
		(Ny-2)   //don't forget the parantheses here, pal ;)
#endif
		) + (2+numberOfQuantities)*j + (2+k);
      adonis_assert(idx_ < (int)data_.size());
      return (data_[idx_]);//(data_[2 + i + Nx*j + npt*quantity]);
    }

    //! NOTE: a file is always printed WITHOUT ghostpoints. Hence, unlike the previous
    //!       fct, no special treatment regarding ghost points is necessary
     const T& get_2D_MOL_entry_no_ghost_points(int i, int j, int Ny, int k, int numberOfQuantities) const{ //without x and y
      adonis_assert((j >= 0) && (j < Ny) );
      idx_ = i*((2+numberOfQuantities)*Ny) + (2+numberOfQuantities)*j + (2+k);
      adonis_assert(idx_ < (int)data_.size());
      return (data_[idx_]);//(data_[2 + i + Nx*j + npt*quantity]);
    }


  };


} //end namespace 

#endif 
