// A class for defining atomic quantities for interscale operations

#ifndef PER_PAIR_QUANTITY_H
#define PER_PAIR_QUANTITY_H

// ATC_Method headers
#include "LammpsInterface.h"
#include "DependencyManager.h"
#include "PerAtomQuantity.h"
#include <map>

using std::map;
using std::pair;

namespace ATC {
  /** mapping of atomic pairs to pair index value */
  typedef map< pair< int,int >,int > PAIR_MAP;
  typedef map< pair< int,int >,int >::const_iterator PAIR_MAP_ITERATOR;
  typedef pair< pair< int,int >,int > ATOM_PAIR;

  /**
   *  @class  PairMap 
   *  @brief  Base class maps of pair indices to a single index
   */

  class PairMap : public DependencyManager {
  public:
    PairMap(LammpsInterface * lammpsInterface, int groupbit);
    virtual ~PairMap(void);
    virtual bool need_reset(void) const = 0; 
    virtual void reset(void) const = 0;
    void quantity() {}; 
    void set_quantity() { throw ATC_Error("inappropriate access to pair map");}
    // iterator interface
    int size(void) const { 
      if (need_reset()) reset(); 
      return nPairs_+nBonds_; 
    }
    int num_bonds(void) const { return nBonds_; }
    virtual ATOM_PAIR  start() const = 0; // const reset / non-const call propagte reset
    virtual ATOM_PAIR  next()  const = 0;
    //ATOM_PAIR& operator++ ()    {return next();} // prefix ++: no parameter, returns a reference
    ATOM_PAIR  operator++ (int) const {return next();} // postfix ++: dummy parameter, returns a value
    virtual bool       finished()  const = 0;
  protected:
    LammpsInterface * lammpsInterface_;
    int groupbit_;
    mutable int nPairs_;
    mutable int nBonds_;
  private:
    PairMap();// do not define
  };

  class PairMapNeighbor : public PairMap {
  public:
    PairMapNeighbor(LammpsInterface * lammpsInterface, int groupbit);
    virtual ~PairMapNeighbor(void) {};
    virtual bool need_reset(void) const;
    virtual void reset(void) const;
    virtual ATOM_PAIR  start(void) const {
      if (need_reset()) reset();
      iterator_ = pairMap_.begin(); return *iterator_;}
    virtual ATOM_PAIR  next(void)  const {
       iterator_++; return *iterator_;}
    virtual bool       finished()  const  { 
         return (iterator_==pairMap_.end());}
  protected:
    mutable PAIR_MAP pairMap_; 
  private:
    mutable PAIR_MAP_ITERATOR iterator_;
    PairMapNeighbor();// do not define
  };

  class PairMapBond : public PairMap {
  public:
    PairMapBond(LammpsInterface * lammpsInterface, int groupbit);
    virtual ~PairMapBond(void) {};
    virtual bool need_reset(void) const { return true;} 
    virtual void reset(void) const {nBonds_ = lammpsInterface_->bond_list_length(); };
    ATOM_PAIR  start() const {
      reset();
//    if (needs_reset()) propagate_reset()
      index_ = 0; return atom_pair(index_);} 
    ATOM_PAIR  next()  const {return atom_pair(++index_);}
    bool       finished()  const  { return index_==nBonds_; }
    ATOM_PAIR  atom_pair(int n) const {
      if ( !(n<nBonds_) ) {
        pair<int,int> pair_ij(-1,-1); // this is the "end" value
        ATOM_PAIR p(pair_ij,n);
        return p;
      }
      int * bond = (lammpsInterface_->bond_list())[n];
      pair<int,int> pair_ij(bond[0],bond[1]); 
      ATOM_PAIR p(pair_ij,n);
      return p;
    }
  private:
    mutable int index_;
    PairMapBond();// do not define
  };

  class PairMapBoth : public PairMapNeighbor {
  public:
    PairMapBoth(LammpsInterface * lammpsInterface, int groupbit);
    virtual ~PairMapBoth(void) {};
    virtual bool need_reset(void) const {
      nBonds_ = lammpsInterface_->bond_list_length(); 
      return PairMapNeighbor::need_reset();
    }
    virtual void reset(void) const { 
     nBonds_ = lammpsInterface_->bond_list_length(); 
     PairMapNeighbor::reset();
    }
    virtual ATOM_PAIR  start(void) const {
      if (need_reset()) reset();
      index_ = 0;
      iterator_ = pairMap_.begin(); 
      return atom_pair(index_); // start with bonds
    }
    virtual ATOM_PAIR  next(void)  const {
      ++index_;
      if (index_ < nBonds_) { return atom_pair(index_);} 
      else { if (index_>nBonds_) iterator_++; return *iterator_; }
    }
    ATOM_PAIR  atom_pair(int n) const {
      int * bond = (lammpsInterface_->bond_list())[n];
      pair<int,int> pair_ij(bond[0],bond[1]); 
      ATOM_PAIR p(pair_ij,n);
      return p;
    }
    virtual bool       finished()  const  { 
         return (iterator_==pairMap_.end());}
  private:
    mutable int index_;
    mutable PAIR_MAP_ITERATOR iterator_;
    PairMapBoth();// do not define
  };

  /**
   *  @class  DensePerPairQuantity 
   *  @brief  Base class for objects that manage pair/bond quantities 
   */
  class DensePerPairMatrix : public MatrixDependencyManager<DenseMatrix, double> {

  public:
    
    // constructor
    DensePerPairMatrix(LammpsInterface * lammpsInterface, 
       const PairMap & pairMap, 
       int nCols = 1);
    
    // destructor
    virtual ~DensePerPairMatrix(){};

    /** access to a constant dense matrix of the quantity */
    virtual const DenseMatrix<double> & quantity() const
      {this->reset(); return MatrixDependencyManager<DenseMatrix, double>::quantity();};

    /** access to a non-constant dens matrix of the quantity */
    virtual DenseMatrix<double> & set_quantity()
      {this->reset(); return MatrixDependencyManager<DenseMatrix, double>::set_quantity();}

    /** number of columns in quantity */
    INDEX nCols() const {return nCols_;};

    /** resets data, if necessary */
    virtual void reset() const = 0;

  protected:
    
    /** pointer to access Lammps data */
    LammpsInterface * lammpsInterface_;

    /** reference to pair map */
    const PairMap & pairMap_;

    /** number of columns of the per atom quantity -static */
    int nCols_;

  private:
    DensePerPairMatrix(); // do not define
  };

  /**
   *  @class  PairVirial 
   *  @brief  f_ab (x) x_ab where (ab) -> p
   */
  class PairVirial : public DensePerPairMatrix {

  public:
    // constructor
    PairVirial(LammpsInterface * lammpsInterface,
      const PairMap & pairMap, int nCols);

    // destructor
    virtual ~PairVirial(){};

    /** resets data, if necessary */
    virtual void reset() const = 0;
    
  private:
    PairVirial(void); // do not define
  };
  /**
   *  @class  PairVirial 
   *  @brief  f_ab (x) x_ab where (ab) -> p
   */
  class PairVirialEulerian : public PairVirial {

  public:
    // constructor
    PairVirialEulerian(LammpsInterface * lammpsInterface,
      const PairMap & pairMap);

    // destructor
    virtual ~PairVirialEulerian(){};

    /** resets data, if necessary */
    virtual void reset() const;
    
  private:
    PairVirialEulerian(void); // do not define
  };
  class PairVirialLagrangian : public PairVirial {

  public:
    // constructor
    PairVirialLagrangian(LammpsInterface * lammpsInterface,
      const PairMap & pairMap, 
      double ** xRef);
//    const PerAtomQuantity<double> * x0);

    // destructor
    virtual ~PairVirialLagrangian(){};

    /** resets data, if necessary */
    virtual void reset() const;
  protected:
    double ** xRef_; // note difficult to make a ** const
//  const PerAtomQuantity<double> * xRef_;
    
  private:
    PairVirialLagrangian(void); // do not define
  };

  /**`
   *  @class  PairPotentialHeatFlux 
   *  @brief  f_ab v_b where (ab) -> p
   */
  class PairPotentialHeatFlux : public DensePerPairMatrix {

  public:
    // constructor
    PairPotentialHeatFlux(LammpsInterface * lammpsInterface,
      const PairMap & pairMap);

    // destructor
    virtual ~PairPotentialHeatFlux(){};

    /** resets data, if necessary */
    virtual void reset() const =0;
    
  private:
    PairPotentialHeatFlux(void); // do not define
  };
  class PairPotentialHeatFluxEulerian : public PairPotentialHeatFlux {

  public:
    // constructor
    PairPotentialHeatFluxEulerian(LammpsInterface * lammpsInterface,
      const PairMap & pairMap);

    // destructor
    virtual ~PairPotentialHeatFluxEulerian(){};

    /** resets data, if necessary */
    virtual void reset() const;
    
  private:
    PairPotentialHeatFluxEulerian(void); // do not define
  };

  class PairPotentialHeatFluxLagrangian : public PairPotentialHeatFlux {

  public:
    // constructor
    PairPotentialHeatFluxLagrangian(LammpsInterface * lammpsInterface,
      const PairMap & pairMap, double ** xRef);

    // destructor
    virtual ~PairPotentialHeatFluxLagrangian(){};

    /** resets data, if necessary */
    virtual void reset() const;
  protected:
    double ** xRef_; // note difficult to make a ** const
    //const PerAtomQuantity<double> * x0_;
    
  private:
    PairPotentialHeatFluxLagrangian(void); // do not define
  };

  /**
   *  @class  SparsePerPairMatrix 
   *  @brief  Base class for objects that manage pair/bond quantities 
   */
  class SparsePerPairMatrix : public MatrixDependencyManager<SparseMatrix, double> {

  public:
    
    // constructor
    SparsePerPairMatrix(LammpsInterface * lammpsInterface, 
       const PairMap & pairMap); 
    
    // destructor
    virtual ~SparsePerPairMatrix(){};

    /** access to a constant dense matrix of the quantity */
    virtual const SparseMatrix<double> & quantity() const
      {reset(); return MatrixDependencyManager<SparseMatrix, double>::quantity();};

    /** access to a non-constant dens matrix of the quantity */
    virtual SparseMatrix<double> & set_quantity()
      {reset(); return MatrixDependencyManager<SparseMatrix, double>::set_quantity();}

    /** resets data, if necessary */
    virtual void reset() const = 0;

  protected:
    
    /** pointer to access Lammps data */
    LammpsInterface * lammpsInterface_;

    /** reference to pair map */
    const PairMap & pairMap_;

  private:
    SparsePerPairMatrix(); // do not define
  };

  /**
   *  @class  BondMatrix
   *  @brief  Hardy's B_Iab wher (ab) -> p
   */
  class BondMatrix : public SparsePerPairMatrix {

  public:
    // constructor
    BondMatrix(LammpsInterface * lammpsInterface,
      const PairMap & pairMap, double ** x_, const class FE_Mesh * feMesh);

    // destructor
    virtual ~BondMatrix(){};

    /** resets data, if necessary */
    virtual void reset() const = 0;

  protected:
    double ** x_;
    const class FE_Mesh * feMesh_;
    
  private:
    BondMatrix(void); // do not define
  };

  /**
   *  @class  BondMatrixKernel
   *  @brief  Hardy's B_Iab wher (ab) -> p
   */
  class BondMatrixKernel : public BondMatrix {

  public:
    // constructor
    BondMatrixKernel(LammpsInterface * lammpsInterface,
      const PairMap & pairMap, 
      double ** x,
      const class FE_Mesh * feMesh,
      const class KernelFunction * kernelFunction);

    // destructor
    virtual ~BondMatrixKernel(){};

    /** resets data, if necessary */
    virtual void reset() const;

  protected:
    const class KernelFunction * kernelFunction_;
    
  private:
    BondMatrixKernel(void); // do not define
  };

  /**
   *  @class  BondMatrixPartitionOfUnity
   *  @brief  Hardy's B_Iab wher (ab) -> p
   */
  class BondMatrixPartitionOfUnity : public BondMatrix {

  public:
    // constructor
    BondMatrixPartitionOfUnity(LammpsInterface * lammpsInterface,
      const PairMap & pairMap, 
      double ** x,
      const class FE_Mesh * feMesh,
      const DIAG_MAN * invVol);

    // destructor
    virtual ~BondMatrixPartitionOfUnity(){};

    /** resets data, if necessary */
    virtual void reset() const;
  protected:
    const DIAG_MAN * invVols_;
    static const int lineNgauss_ = 10;
    double lineXg_[lineNgauss_], lineWg_[lineNgauss_];
  private:
    BondMatrixPartitionOfUnity(void); // do not define
  };
}

#endif
