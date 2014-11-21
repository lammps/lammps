// A class for wrapping matrix operations with dependency information to speed up execution

#ifndef DEPENDENCY_MANAGER_H
#define DEPENDENCY_MANAGER_H

#include "ATC_TypeDefs.h"
#include "MatrixLibrary.h"
#include "ATC_Error.h"
#include "MPI_Wrappers.h"

namespace ATC {

  class InterscaleManager;

  /** memory type */
  enum MemoryType
  {
    TEMPORARY = 0,
    PERSISTENT
  };

  /**
   *  @class  DependencyManager 
   *  @brief  Base class for defining objects that manage the dependencies of various objects
   */

  class DependencyManager {

  public:

    // used as a friend so it can perform a depth-first search to have safe deletions of managed dependencies
    friend class InterscaleManager;
    
    // constructor
    DependencyManager() : needReset_(true), isFixed_(false), memoryType_(TEMPORARY), dfsFound_(false) {};
    
    // destructor
    virtual ~DependencyManager() {};
    
    /** registration by other PerAtomQuantity objects */
    void register_dependence(DependencyManager * dependentQuantity)
    {dependentQuantities_.insert(dependentQuantity);};

    /** removes dependencies from the set */
    void remove_dependence(DependencyManager * dependentQuantity)
    {dependentQuantities_.erase(dependentQuantity);};

    /** check if a reset is required */
    bool need_reset() const {return needReset_ && !isFixed_;};
    
    /** propagate need to reset to to dependencies */
    void propagate_reset()
    {
      if (!isFixed_) {
        std::set<DependencyManager *>::iterator it;
        for (it = dependentQuantities_.begin(); it != dependentQuantities_.end(); it++)
          (*it)->force_reset();
      }
    };

    /** actions associated with indicating this quantity requires a reset */
    void set_reset()
    {
      needReset_ = true;
    }

    /** flip this object to needing a reset, and get dependencies */
    void force_reset()
    {
      set_reset();
      propagate_reset();
    }

    /** force quantity to be held fixed, enables dependent quantity to be used as persistent storage */
    void fix_quantity() {isFixed_ = true;};

    /** unfix the quantity */
    void unfix_quantity()
    {
      if (isFixed_) {
        isFixed_ = false;
        if (needReset_) propagate_reset();
      }
    };

    /** check on the memory type of the quantity */
    MemoryType memory_type() const {return memoryType_;};

    /** set the memory type of the quantity */
    void set_memory_type(MemoryType memoryType) {memoryType_ = memoryType;};


  protected:
    
    /** list of dependent atomic quantities */
    std::set<DependencyManager * > dependentQuantities_;
    
    /** flag for needing a recent */
    // mutable is applied because there can be internal updates because we update when needed rather than when pushed
    mutable bool needReset_;

    /** flag for if quantity is being held fixed */
    bool isFixed_;

    /** flag for if the quantity is temporary (per-run) */
    MemoryType memoryType_;

    /** flag for if the node has been found in depth-first search */
    bool dfsFound_;
    
  };
  
  /**
   *  @class  MatrixDependencyManager
   *  @brief  Class for defining objects that manage the dependencies of matrices
   */

  // Matrix class T, underlying type U
  template <template <typename> class T, typename U>
  class MatrixDependencyManager : public DependencyManager {

  public:

    MatrixDependencyManager() {};
    MatrixDependencyManager(int nRows, int nCols) : quantity_(nRows,nCols) {};
    virtual ~MatrixDependencyManager() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual T<U> & set_quantity() {propagate_reset(); return get_quantity();};
    
    /** access to a constant dense matrix of the quantity, indexed by AtC atom counts */
    virtual const T<U> & quantity() const {return get_quantity();};

    /** number of rows in quantity */
    virtual int nRows() const {return (this->quantity()).nRows();};

    /** number of columns in quantity */
    virtual int nCols() const {return (this->quantity()).nCols();};

    /** size of the matrix */
    virtual int size() const {return (this->quantity()).size();};

    /** reset the quantities size */
    void reset(INDEX nRows, INDEX nCols) {get_quantity().reset(nRows,nCols); propagate_reset();};

    /** resize the quantities size */
    void resize(INDEX nRows, INDEX nCols) {get_quantity().resize(nRows,nCols); propagate_reset();};

    /** sets the quantity to a given value */
    virtual void operator=(const T<U> & target) {get_quantity()=target; propagate_reset();};
    /** sets the quantity to a given constant value */
    virtual void operator=(U target) {get_quantity()=target; propagate_reset();};

    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const T<U> & addition) {get_quantity()+=addition; propagate_reset();};
    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator+=(U addition) {get_quantity()+=addition; propagate_reset();};

    /** adds the given data to the Lammps quantity */
    virtual void operator-=(const T<U> & subtraction) {get_quantity()-=subtraction; propagate_reset();};
    /** adds the scalar data to the Lammps quantity for AtC atoms */
    virtual void operator-=(U subtraction) {get_quantity()-=subtraction; propagate_reset();};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const T<U> & multiplier) {get_quantity()*=multiplier; propagate_reset();};
    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(U multiplier) {get_quantity()*=multiplier; propagate_reset();};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const T<U> & divisor) {get_quantity()/=divisor; propagate_reset();};
    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(U divisor) {get_quantity()/=divisor; propagate_reset();};

    // I have no idea why these won't compile (JAT, 04/07/11)
    /** adds the given data to the Lammps quantity */
    virtual void operator+=(const MatrixDependencyManager<T,U> & addition) {get_quantity()+=addition.quantity(); propagate_reset();};

    /** adds the given data to the Lammps quantity */
    virtual void operator-=(const MatrixDependencyManager<T,U> & subtraction) {get_quantity()-=subtraction.quantity(); propagate_reset();};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator*=(const MatrixDependencyManager<T,U> & multiplier) {get_quantity()*=multiplier.quantity(); propagate_reset();};

    /** multiples the Lammps quantity by the given data, input is indexed in AtC atom counts */
    virtual void operator/=(const MatrixDependencyManager<T,U> & divisor) {get_quantity()/=divisor.quantity(); propagate_reset();};

    /** execute the matrix print command */
    virtual void print(const std::string &name) const {get_quantity().print(name);};

  protected:

    // This getter can be overridden by derived classes if they need e.g. a
    // differently-constructed quantity, but would like to use the same operators.
    virtual T<U> &get_quantity() const { return quantity_; }

    /** matrix */
    // mutable is applied because there can be internal updates because we update when needed rather than when pushed
    mutable T<U> quantity_;

  };


  /**
   *  @class MatrixDependencyManager<ParSparseMatrix, T>
   *  @brief Class for defining objects that manage the dependencies of parallelized sparse matrices
   */
  template<typename T>
  class MatrixDependencyManager<ParSparseMatrix, T> :
    public MatrixDependencyManager<SparseMatrix, T>
  {
  public:

    MatrixDependencyManager(MPI_Comm comm) : 
      MatrixDependencyManager<SparseMatrix, T>(), quantity_(comm) {};
    
    MatrixDependencyManager(MPI_Comm comm, int nRows, int nCols) :
      MatrixDependencyManager<SparseMatrix, T>(), quantity_(comm, nRows, nCols) {};
    
    virtual ~MatrixDependencyManager() {};

  protected:

    // Let the superclass's operators work on our ParSparseMatrix.
    virtual ParSparseMatrix<T> &get_quantity() const { return quantity_; }

    mutable ParSparseMatrix<T> quantity_;

  };


  /**
   *  @class MatrixDependencyManager<ParDiagonalMatrix, T>
   *  @brief Class for defining objects that manage the dependencies of parallelized diagonal matrices
   */
  template<typename T>
  class MatrixDependencyManager<ParDiagonalMatrix, T> :
    public MatrixDependencyManager<DiagonalMatrix, T>
  {
  public:

    MatrixDependencyManager(MPI_Comm comm) : 
      MatrixDependencyManager<DiagonalMatrix, T>(), quantity_(comm) {};
    
    MatrixDependencyManager(MPI_Comm comm, int nRows, int nCols) :
      MatrixDependencyManager<DiagonalMatrix, T>(), quantity_(comm, nRows, nCols) {};
    
    virtual ~MatrixDependencyManager() {};

  protected:

    // Let the superclass's operators work on our ParDiagonalMatrix.
    virtual ParDiagonalMatrix<T> &get_quantity() const { return quantity_; }

    mutable ParDiagonalMatrix<T> quantity_;

  };

  /**
   *  @class  SetDependencyManager
   *  @brief  Class for defining objects that manage the dependencies of standard library sets
   */
  template <typename T>
  class SetDependencyManager : public DependencyManager {

  public:

    // constructor
    SetDependencyManager() :
      DependencyManager(), quantity_() {};
    
    // destructor
    virtual ~SetDependencyManager() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual std::set<T> & set_quantity() {propagate_reset(); return quantity_;};
    
    /** access to a constant dense matrix of the quantity, indexed by AtC atom counts */
    virtual const std::set<T> & quantity() const {return quantity_;};

    /** size of the set */
    virtual int size() const {return (this->quantity()).size();};

  protected:

    /** underlying set */
    // mutable is applied because there can be internal updates because we update when needed rather than when pushed
    mutable std::set<T> quantity_;

  };

  /**
   *  @class  VectorDependencyManager
   *  @brief  Class for defining objects that manage the dependencies of standard library vectors
   */
  template <typename T>
  class VectorDependencyManager : public DependencyManager {

  public:

    // constructor
    VectorDependencyManager() :
      DependencyManager(), quantity_() {};
    
    // destructor
    virtual ~VectorDependencyManager() {};

    /** returns a non-const version for manipulations and changes, resets dependent quantities */
    virtual std::vector<T> & set_quantity() {propagate_reset(); return quantity_;};
    
    /** access to a constant dense matrix of the quantity, indexed by AtC atom counts */
    virtual const std::vector<T> & quantity() const {return quantity_;};

    /** size of the set */
    virtual int size() const {return (this->quantity()).size();};

  protected:

    /** underlying set */
    // mutable is applied because there can be internal updates because we update when needed rather than when pushed
    mutable std::vector<T> quantity_;

  };

};

#endif
