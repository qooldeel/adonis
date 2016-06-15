#ifndef PERMUTATION_FIELD_MATRIX_HH
#define PERMUTATION_FIELD_MATRIX_HH

#include "../dunestuff/fvector.hh"
#include "../common/adonisassert.hh"


namespace Adonis{

  /**
  * \brief A permutation matrix class based on Dune::FieldVector.
  *
  * Permutation matrices are the unity matrix with columns being interchanged.
  *
  * \tparam N order of the matrix. 
  * A Permutation matrix is an identity matrix whose columns have been swapped.
  * 
  * Storage format: store column index representing the position where the sole  1s in P's kth row stands. See GOLUB/VAN LOAN, Matrix computations, § 4.4, p. 109
  */
  template<int N>
  class PermutationFieldMatrix{
  public:
    typedef Dune::FieldVector<int,N> VType;
    typedef PermutationFieldMatrix<N> ThisType;
    typedef int value_type; //simple: int is value_type 
    

    typedef int* iterator;
    typedef const int* const_iterator;
    
    //! default construction
    PermutationFieldMatrix(){}

    //!construct from given FieldVector
    PermutationFieldMatrix(const VType& v){
      for(int i = 0; i < N; ++i)
	idx_[i] = v[i];
    }

    //! construction from Iterators
    template<class Iter>
    PermutationFieldMatrix(Iter it1, Iter it2){
      int x = 0;
      for(Iter it = it1; it != it2; ++it)
	idx_[x++] = static_cast<int>(*it);
    }

    //! accessible via 'P.dimension' and 'P.order'
    enum {dimension = N, order = N};

    //!copy constructor
    PermutationFieldMatrix(const PermutationFieldMatrix& P){
       for(int i = 0; i < N; ++i)
	idx_[i] = P[i];
    }

    //!assignment operator
    PermutationFieldMatrix& operator=(const PermutationFieldMatrix& P){
      if(this != &P){
	for(int i = 0; i < N; ++i)
	  idx_[i] = P[i];
      }
      return *this;
    }

    const int& operator[](int i) const{
      adonis_assert((i >= 0) && (i < N));
      return idx_[i];
    }

    int& operator[](int i){
      adonis_assert((i >= 0) && (i < N));
      return idx_[i];
    }

    //!assign entries from given FieldVector
    PermutationFieldMatrix& operator=(const VType& v){
      for(int i = 0; i < N; ++i)
	idx_[i] = v[i];
      return *this;
    }

    
    //!assign entries from given STL-compliant random access container 
    template<class V>
    PermutationFieldMatrix& operator=(const V& v){
      adonis_assert((int)v.size() == N);
      for(int i = 0; i < N; ++i)
	idx_[i] = static_cast<int>(v[i]);
      return *this;
    }

    //! the unity matrix
    inline void unity() {
      for(int i = 0; i < N; ++i)
	idx_[i] = i;
    }

    
    //! tiny iterator interface
    iterator begin(){
      return &(idx_[0]);  //if you use VTpye as field: idx_.begin();
    }

    iterator end(){
      return  &(idx_[N]);           //idx_.end();
    }

    const_iterator begin() const{
      return &(idx_[0]);             //idx_.begin();
    }

    const_iterator end() const{
      return &(idx_[N]);             //idx_.end();
    }


    friend std::ostream& operator<<(std::ostream& os, const ThisType& P){
      os <<"A permutation matrix "<<std::endl; 
      os << "Human readable format -- it is NOT stored like that ;) "<<std::endl;
      for(int i = 0; i < N; ++i){
	for(int j = 0; j < N; ++j){
	  if(j == P[i])
	    os << 1 << "  ";
	  else 
	    os << 0 << "  ";
	}
	os << std::endl;
      }
	
      os << std::endl << "Storage format, i.e. column index (per row) that contains '1':"<<std::endl;
      for(int i = 0; i < N; ++i){
	os << P[i] << " ";
      }
      os << std::endl;
       
      return os;
    }

  private:
    int idx_[N];     //very simple: an in-built array
 
  };




  //operatations with Permutation matrices. 
  /**
   * \brief Right multiplication with permutation matrix, i.e. \f$ P*M \f$ 
   */
  template<class K, int M, int N>
  inline Dune::FieldMatrix<K,M,N> operator*(const PermutationFieldMatrix<M>& P, const Dune::FieldMatrix<K,M,N>& F){

    Dune::FieldMatrix<K,M,N> A; 

    for(int i = 0; i < M; ++i){
      for(int j = 0; j < N; ++j){
	A[i][j] = F[P[i]][j];
      }
    }
    return A;
  }
  


  /**
   * \brief Transpose of given permutation matrix. Needed for left multiplication.
   *
   * column index is stored in idx_ for row 0,1,.... . The indices are completed by the row index and then swapped. Finally, those index pairs are sorted by 1st entries (row indices) and the so rearranged column index corresponds to the transposed matrix. 
   *
   * Costs: Dominated by the sorting routine: \f$ \mathcal{O}(NlogN) \f$
   */
  template<int N> 
  inline PermutationFieldMatrix<N> transpose(const PermutationFieldMatrix<N>& P){
    // Dune::FieldVector<std::pair<int, int>, N> p;
    std::pair<int, int> p[N];

    for(int i = 0; i < N; ++i){                            //O(N)
      p[i].first = P[i];       //already in swapped form
      p[i].second = i; 
    }

    //sort using merge sort to sort array p
    sort_me(p,N);   //sort_me(p); //works                  //O(NlogN)

    PermutationFieldMatrix<N> Trans;
    
    for(int i = 0; i < N; ++i)                             //O(N)
      Trans[i] = p[i].second;

    return Trans;
  }


  /**
   * \brief Left multiplication with permutation matrix, i.e. \f$ M*P \f$
   *
   * Needed: Transposed Permutation matrix
   * Overall Costs: \f$ \mathcal{O}(MN) \approx \mathcal{O}(N^2), \ \textrm{for} \ M \approx N.\f$ 
   */
  template<class K, int M, int N>
  inline Dune::FieldMatrix<K,M,N> operator*(const Dune::FieldMatrix<K,M,N>& F, const PermutationFieldMatrix<N>& P ){

    Dune::FieldMatrix<K,M,N> A; 
   
    PermutationFieldMatrix<N> PTrans(transpose(P)); //transpose of P 
      
    for(int j = 0; j < N; ++j){
      for(int i = 0; i < M; ++i){
	A[i][j] = F[i][PTrans[j]];  //A[i][j] = F[i][P[j]]; //if P^T is known
      }
    }
    return A;
  }
  

  /**
   * \brief Multiplication between two permutations
   *
   * Costs: \f$ \mathcal{O}(N^2)\f$
   */
  template<int N>
  inline PermutationFieldMatrix<N> operator*(const PermutationFieldMatrix<N>& P1, const PermutationFieldMatrix<N>& P2){

    adonis_assert(N > 0);  //makes no sense when N is 0

    PermutationFieldMatrix<N> 
      Product;
    
    if(N == 1) 
      Product[0] = 0;  //i.e. the only entry can be found in the sole col 
    
    else{
      
      int x = -1;  //initialize it with some stupid value 
      for(int i = 0; i < N; ++i){
	for(int j = 0; j < N; ++j){          
	  if(P1[i] == j){                             //N² (at most N per idx)
	    x = j;
	    //std::cout << "j = "<< j << std::endl; 
	    break; //leave inner loop
	  }
	}
	Product[i] = P2[x]; //x;
      } 
    
    }                                             //_______________
    return Product;                                  
  }                                               //     O(N²)


}//end namespace 

#endif 
