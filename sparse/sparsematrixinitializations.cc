#include <iostream>
#include <vector>
#include <cassert>

template<class IT>
inline void print_me_yeah(std::size_t dim, IT it){
    for(std::size_t i = 0; i < dim; ++i){
        std::cout << *it << " ";
        ++it;
    }
    std::cout << std::endl;
}

//! Tentative class to test some compressed sparse features
template<char CHAR, class T>
class SparseMatrix{
  public:
    typedef std::size_t SizeType;
   
   //default construction
    SparseMatrix():m_(0),n_(0),nonz_(0),dimrowA_(0),dimcolA_(0),rowA_(0),colA_(0),valA_(0),isInitialized_(false){}
    
    
    void resize(SizeType m, SizeType n, SizeType nnz){
        m_ = m;
        n_ = n;
        nonz_ = nnz;
        dimrowA_ = ( ((CHAR=='c') || (CHAR=='C')) ? nonz_ : ( ((CHAR=='r') || (CHAR=='R')) ? (m_+1) : 0 ));
        dimcolA_ = ( ((CHAR=='c') || (CHAR=='C')) ? n_+1 : ( ((CHAR=='r') || (CHAR=='R')) ? (nonz_) : 0 ));
        
        if((rowA_ != 0) && (colA_ != 0) && (valA_ != 0)){
            delete[] rowA_;
            delete[] colA_;
            delete[] valA_;
        }
        rowA_ = new SizeType[dimrowA_];
        colA_ = new SizeType[dimcolA_];
        valA_ = new T[nonz_];
        
        //initialize arrays --- better it is to prevent spurious results ;)
        for(SizeType l = 0; l < dimrowA_; ++l)
           rowA_[l] = 0;
        for(SizeType l = 0; l < dimcolA_; ++l)
           colA_[l] = 0;
        for(SizeType l = 0; l < nonz_; ++l)   
           valA_[l] = T(0);
        //std::cout << "rowA_["<<l<<"] = "<< rowA_[l] << std::endl;
        
         //std::cout << "***"<<std::endl;
        isInitialized_ = true;
    }
   
   //destructor
   ~SparseMatrix(){
        delete[] rowA_;
        delete[] colA_;
        delete[] valA_;
       
   }
   
  void  print_row_array(){
      assert(isInitialized_);
       for(SizeType l = 0; l < dimrowA_;++l){
        std::cout << l << ".)  "<< rowA_[l] << " ";
       }
       std::cout << std::endl;
   }
   
   SizeType* get_row_address() {return rowA_;}
   
   const SizeType dim_of_row_array() const {return dimrowA_;}
   const SizeType dim_of_col_array() const {return dimcolA_;}
   
   //! CAUTION: no size check, it's up to you to provide properly sized
   //!          input data!
    template<class IT>
    void fill_row_array(IT rowDataIt){
        assert(isInitialized_);
        for(SizeType l = 0; l < dimrowA_; ++l){
            rowA_[l] = *(rowDataIt++); //rowDataIt[l];
        
        }
        
    }
    
    
  private:
    SizeType m_,n_,nonz_,
        dimrowA_, dimcolA_;  //proper dimension in dependence of compressed
                             // storage format
   
    SizeType* rowA_;    
    SizeType* colA_;
    T* valA_;
    mutable bool isInitialized_; 
};

using namespace std;

int main()
{
    size_t m = 3,
        n = 4,
        nnz = 5;
// 'C'-case 
    size_t rowC[] = {0,1,1,2,3};
    size_t colC[] = {0,2,3,4,  5};
   vector<size_t> row(rowC,rowC+5),
    col(colC,colC+5);
   
   cout << "Test me..." << endl; 
   SparseMatrix<'C',double> SPM;
   
   SPM.resize(m,n,nnz);
   SPM.fill_row_array(rowC);
   //SPM.print_row_array();
   cout << "Print stuff..."<<endl;
   print_me_yeah(SPM.dim_of_row_array(),SPM.get_row_address());
   
   cout << "----------------- TEST 2 -----------------"<<endl; 
   m = 4;
   n = 2;
   nnz = 4;
   size_t rowC1[] = {0,1, 1, 0};
   size_t colC1[] = {0, 2, 3,   4};
   
   SPM.resize(m,n,nnz);
   SPM.fill_row_array(rowC1);
   cout << "Print stuff..."<<endl;
   print_me_yeah(SPM.dim_of_row_array(),SPM.get_row_address());
   
    cout << "----------------- TEST 3 -----------------"<<endl; 
   m = 4;
   n = 4;
   nnz = 5;
   size_t rowC2[] = {1,2, 0,2, 3};
   size_t colC2[] = {0, 2, 4,   nnz};
   
   SPM.resize(m,n,nnz);
   SPM.fill_row_array(rowC2);
   cout << "Print stuff..."<<endl;
   print_me_yeah(SPM.dim_of_row_array(),SPM.get_row_address());
   return 0;
}
