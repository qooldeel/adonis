#ifndef INDEX_MANIPULATION_HH
#define INDEX_MANIPULATION_HH

#include<vector>
#include "../common/adonisassert.hh"

namespace Adonis{
  
  /**
   *\brief given an index of reaction progress variables (rvp), this class constructs the unrepresented index, i.e. the index of non-rvps.
   *
   * \tparam FULL dimension of the full PDE system
   * \tparam RED (< FULL) dimenson of the reduced system
   * \code
      IndexManipulator<F,R> IM;
      IM.create_unrepresented_index(ices);
      std::cout << IM << std::endl;
   * \endcode
   */
  template<int FULL, int RED, 
	   template<class D, class A = std::allocator<D> > 
	   class V = std::vector>
  class IndexManipulator{
  public:
    typedef V<int> VType;
    
    IndexManipulator(){ adonis_assert(FULL > RED);}
    
    VType& get_index(){
      adonis_assert(idx_.size() != 0);
      return idx_;
    }
    const VType& get_index() const {
      adonis_assert(idx_.size() != 0);
      return idx_;
    }

    friend  std::ostream& operator<<(std::ostream& os, const IndexManipulator& c){
      adonis_assert(c.idx_.size() != 0);

      os << " Index of unrepresented variables = "<<std::endl;
      const int nu = FULL-RED;
      os << " ";
      for(int i = 0; i < nu; ++i){
	os << c.idx_[i] <<" ";
      }
      os<<std::endl;
      return os;
    }


    /**
     * \brief This creates the index of unrepresented species and stores it in a STL-compliant random access container 'idx_' 
     * \param RandomAccessContainer can be any vector like object, e.g. a Dune::FieldVector as well as a std::vector
     * 
     * NOTE: no size-checking on input is done here
     */
    template<class RandomAccessContainer>
    inline void create_unrepresented_index(const RandomAccessContainer& red){
      bool isEqual = false;
      
      int securityCount = 1;
      const int nu = FULL-RED;
      idx_.reserve(nu);
      
      for(int i = 0; i < FULL; ++i){
	isEqual = false;   //reset
	for (int j = 0; j < RED; ++j){
	  if(i == static_cast<int>(red[j])){
	    isEqual = true;
	    break;
	  }
	}
	if(isEqual==false){
	  idx_.push_back(i);
	  adonis_assert(securityCount <= nu);
	  //std::cout<< "securityCount == "<<securityCount<<std::endl;
	  securityCount++;
	}
      }
    }
 

  private:
    VType idx_;

  };

}

#endif 
