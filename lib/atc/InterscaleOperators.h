// An interscale operator class for sharing definitions of atomic quantities, e.g., temperature
// between different parts of the code

#ifndef INTERSCALE_MANAGER_H
#define INTERSCALE_MANAGER_H

#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"
#include "ATC_Error.h"
#include "LammpsInterface.h"
#include "PerAtomQuantity.h"
#include "PerPairQuantity.h"
#include "FundamentalAtomicQuantity.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>

namespace ATC {

  // forward declarations
  class ATC_Method;
  class SmallMoleculeSet;

  /**
   *  @class  InterscaleManager
   *  @brief  Handles definitions for atomistic quantities
   */

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class InterscaleManager
  //--------------------------------------------------------
  //--------------------------------------------------------

  class InterscaleManager {
  
  public:
  
    // constructor
    InterscaleManager(ATC_Method * atc);
    
    // destructor
    ~InterscaleManager();

    /** delete all allocated data */
    void clear();

    /** delete non-persistent data */
    void clear_temporary_data();

    /** set lammps data prefix */
    void set_lammps_data_prefix();
        
    /** parser/modifier */
    bool modify(int narg, char **arg);

    /** pre time integration */
    void initialize();

    // access to per atom quantity objects
    /** access to fundamental atomic quantities */
    FundamentalAtomQuantity * fundamental_atom_quantity(LammpsInterface::FundamentalAtomQuantity id,
                                                        AtomType atomType = INTERNAL);
    /** access to double per atom quantities */
    PerAtomQuantity<double> * per_atom_quantity(const std::string & tag);
    /** access to integer per atom quantities */
    PerAtomQuantity<int> * per_atom_int_quantity(const std::string & tag);
    /** access to double per atom diagonal matrices */
    PerAtomDiagonalMatrix<double> * per_atom_diagonal_matrix(const std::string & tag);
    /** access to double per atom sparse matrices */
    PerAtomSparseMatrix<double> * per_atom_sparse_matrix(const std::string & tag);
    /** access to pair maps */
    PairMap * pair_map(const std::string & tag);

    // addition of new atom quantities, note provider must allocate but the manager will clean-up
    /** addition of a double atomic quantity */
    void add_per_atom_quantity(PerAtomQuantity<double> * atomQuantity,
                               const std::string & tag);
    /** addition of an integer atomic quantity */
    void add_per_atom_int_quantity(PerAtomQuantity<int> * atomQuantity,
                                   const std::string & tag);
    /** addition of a double atomic diagonal matrix */
    void add_per_atom_diagonal_matrix(PerAtomDiagonalMatrix<double> * atomQuantity,
                                      const std::string & tag);
    /** addition of a double atomic sparse matrix */
    void add_per_atom_sparse_matrix(PerAtomSparseMatrix<double> * atomQuantity,
                                    const std::string & tag);
    /** addition of an pair map */
    void add_pair_map(PairMap * pairMap, const std::string & tag);

    /** access to dense matrices */
    DENS_MAN * dense_matrix(const std::string & tag);

    /** addition of dense matrices */
    void add_dense_matrix(DENS_MAN * denseMatrix,
                          const std::string & tag);

    /** access integer dense matrices */
    MatrixDependencyManager<DenseMatrix, int> * dense_matrix_int(const std::string & tag);

    /** addition of integer dense matrices */
    void add_dense_matrix_int(MatrixDependencyManager<DenseMatrix, int> * denseMatrix,
                               const std::string & tag);

    /** access boolean dense matrices */
    MatrixDependencyManager<DenseMatrix, bool> * dense_matrix_bool(const std::string & tag);

    /** addition of boolean dense matrices */
    void add_dense_matrix_bool(MatrixDependencyManager<DenseMatrix, bool> * denseMatrix,
                               const std::string & tag);

    /** access to sparse matrices */
    SPAR_MAN * sparse_matrix(const std::string & tag);

    /** addition of a sparse matrix */
    void add_sparse_matrix(SPAR_MAN * sparseMatrix,
                           const std::string & tag);

    /** access to diagonal matrices */
    DIAG_MAN * diagonal_matrix(const std::string & tag);

    /** addition of a diagonal matrix */
    void add_diagonal_matrix(DIAG_MAN * diagonalMatrix,
                             const std::string & tag);

    /** access to vectors of sparse matrices */
    VectorDependencyManager<SPAR_MAT * > * vector_sparse_matrix(const std::string & tag);

    /** addition of a vector of sparse matrices */
    void add_vector_sparse_matrix(VectorDependencyManager<SPAR_MAT * > * sparseMatrix,
                                  const std::string & tag);

    /** access to sets of ints */
    SetDependencyManager<int> * set_int(const std::string & tag);

    /** addition of a set of ints */
    void add_set_int(SetDependencyManager<int> * sparseMatrix,
                     const std::string & tag);

    /** access to molecule sets */
    SmallMoleculeSet * small_molecule_set(const std::string & tag);

    /** addition of a transfer operator */
    void add_small_molecule_set(SmallMoleculeSet * moleculeSet,
                                const std::string & tag);

    /** addition of exchange list object */
    void add_to_exchange_list(const std::string & tag);

    /** searches through all lists to see if a tag is registered */
    DependencyManager * find(const std::string & tag);

    /** schedules a quantity for deletion, if it exists */
    void remove(const std::string & tag);

    /** size communicated quantities initially */
    void size_comm_quantities();

    /** resets nlocal count of managed atomic quantities which do not perform parallel exchange */
    void reset_nlocal();

    /** resets specific lammps fundamental quantities data, as needed, to account for times when lammps can change quantities */
    void fundamental_force_reset(unsigned quantity);

    /** resets all lammps data, as needed, to account for times when lammps can change quantities */
    void lammps_force_reset();

    /** syncs lammps data to managed objects for parallel communication */
    void prepare_exchange();

    /** syncs managed objects to lammps data after parallel communication */
    void post_exchange();

    /** returns how much lammps memory is used in this function */
    int memory_usage() const;

    /** packs up data for parallel transfer, called from pack_exchange */
    int pack_exchange(int i, double *buffer);

    /** unpacks data after parallel transfer, called from unpack_exchange */
    int unpack_exchange(int i, double *buffer);

    /** packs up data for parallel transfer to ghost atoms on other processors */
    int pack_comm(int index, double *buf, 
                  int pbc_flag, int *pbc);

    /** unpacks data after parallel transfer to ghost atoms on other processors */
    int unpack_comm(int index, double *buf);

    /** changes size of temperary lammps storage data if transfer is being used */
    void grow_arrays(int nmax);

    /** rearrange memory of temporary lammps storage data, called from copy_array */
    void copy_arrays(int i, int j);

  protected:
  
    /** pointer to access ATC methods */
    ATC_Method * atc_;

    /** flag for if first initialization has happened */
    bool initialized_;

    /** containers for fundamental atom quantities, set on request */
    std::vector<std::vector<FundamentalAtomQuantity* > > fundamentalAtomQuantities_;

    /** container for per-atom quantities using dense matrices of doubles */
    std::map<std::string, PerAtomQuantity<double> * > perAtomQuantities_;

    /** container for integer atom quantities, set by AtC classes */
    std::map<std::string, PerAtomQuantity<int> * > perAtomIntQuantities_;

    /** container for per-atom quantities using diagonal matrices of doubles */
    std::map<std::string, PerAtomDiagonalMatrix<double> * > perAtomDiagonalMatrices_;

    /** container for per-atom quantities using sparse matrices of doubles */
    std::map<std::string, PerAtomSparseMatrix<double> * > perAtomSparseMatrices_;

    /** container for pair maps */
    std::map<std::string, PairMap * > pairMaps_;

    /** container for dense matrices */
    std::map<std::string, DENS_MAN * > denseMatrices_;

    /** container for dense matrices for integer quantities */
    std::map<std::string, MatrixDependencyManager<DenseMatrix, int> * > denseMatricesInt_;

    /** container for dense matrces for boolean quantities */
    std::map<std::string, MatrixDependencyManager<DenseMatrix, bool> * > denseMatricesBool_;

    /** container for sparse matrices */
    std::map<std::string, SPAR_MAN * > sparseMatrices_;

    /** container for diagonal matrices */
    std::map<std::string, DIAG_MAN * > diagonalMatrices_;

    /** container for vectors of vectors of sparse matrices */
    std::map<std::string, VectorDependencyManager<SPAR_MAT * > * > vectorSparMat_;

    /** container for sets of integer quantities */
    std::map<std::string, SetDependencyManager<int> * > setInt_;

    /** container for molecule sets */
    std::map<std::string, SmallMoleculeSet * > smallMoleculeSets_;

    /** container for atomic quantities which must be transfered when atoms cross processors */
    std::set<PerAtomQuantity<double> *> exchangeList_;

    /** container for atomic quantities which must be transfered to ghost atoms on other processors */
    std::vector<PerAtomQuantity<double> *> commList_;

    /** container for integer atomic quantities which must be transfered to ghost atoms on other processors */
    std::vector<PerAtomQuantity<int> *> commIntList_;

    /** container for atomic diagonal matrices which must be transfered to ghost atoms on other processors */
    std::vector<PerAtomDiagonalMatrix<double> *> commDmList_;

    /** container for atomic sparse matrices which must be transfered to ghost atoms on other processors */
    std::vector<PerAtomSparseMatrix<double> *> commSmList_;

    /** prefix for labeling associated lammps arrays */
    std::string prefix_;

    /** order of deletion list of managed quantities */
    std::vector<DependencyManager * > deletionList_;

    /** creates a reverse sorted depth-first search list for deleting managed quantities */
    void create_deletion_list();

    /** executes a depth-first search visit on a managed quantity */
    int dfs_visit(DependencyManager * quantity, const int index);

    /** helper function to access a data entry in a list */
    template <typename data>
    data * return_quantity(std::map<std::string,data * > & list, const std::string & tag)
    {
      typename std::map<std::string,data * >::iterator it = list.find(tag);
      if (it==list.end()) return NULL;
      return it->second;
    }

    /** helper function to add a data entry to a list */
    template <typename data>
    void add_quantity(std::map<std::string,data * > & list, data * quantity, const std::string & tag)
    {
      typename std::map<std::string,data * >::iterator it = list.find(tag);
      if (it!=list.end())
        throw ATC_Error("Tried to add another Quantity with tag "+tag+" in InterscaleManager::add_quantity");
      typename std::template pair<std::string,data * > myPair(tag,quantity);
      list.insert(myPair);
    }

    /** helper function to add a data entry to a list when it requires neighbor communication*/
    template <typename data>
    void add_comm_quantity(std::map<std::string,data * > & list, std::vector<data * > & commList, data * quantity, const std::string & tag)
    {
      add_quantity(list,quantity,tag);
      // allocate data for parallel communication
      quantity->grow_lammps_array(LammpsInterface::instance()->nmax(),prefix_+tag);
      if (quantity->atom_type() == PROC_GHOST) {
        commList.push_back(quantity);
      }
    }

     /** helper function to fina a data entry in a list */
    template <typename data>
    data * find_in_list(std::map<std::string,data * > & list, const std::string & tag)
    {
      typename std::map<std::string,data * >::iterator it = list.find(tag);
      if (it!=list.end()) return it->second;
      return NULL;
    }

    /** helper function to force the reset of all data in a list */
    template <typename data>
    void force_reset_loop(std::map<std::string,data * > & list)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it)
        (it->second)->force_reset();
    }

    /** helper function to set the memory type to temporary of a list */
    template <typename data>
    void set_memory_temporary(std::map<std::string,data * > & list)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it)
        (it->second)->set_memory_type(TEMPORARY);
    }

    /** helper function to perform intialization for dfs of a list */
    template <typename data>
    void dfs_prepare_loop(std::map<std::string,data * > & list)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it) {
        (it->second)->dfsFound_ = false;
      }
    }

    /** helper function to start the dfs visit for list */
    template <typename data>
    void dfs_visit_loop(std::map<std::string,data * > & list,
                        int & index)
    {
      typename std::map<std::string,data* >::iterator it = list.begin();
      while (it != list.end()) {
        if  (!((it->second)->dfsFound_)) index = dfs_visit(it->second,index);
        if ((it->second)->memory_type()==TEMPORARY) list.erase(it++);
        else ++it;
      }
    }

    // PAQ helper functions
    /** helper function to adjust local atom count for all data in a list before exchange, only valid with quantities that do that are aware of atom counts */
    template <typename data>
    void reset_nlocal_loop(std::map<std::string,data * > & list)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it) {
        (it->second)->reset_nlocal();
      }
    }

    /** helper function to indicate lammps data is stale for all data in a list before exchange, only valid with PAQs */
    template <typename data>
    void lammps_reset_loop(std::map<std::string,data * > & list)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it)
        (it->second)->lammps_force_reset();
    }

    /** helper function to size all data in a list, only valid with comm lists */
    template <typename data>
    void size_comm_loop(std::vector<data * > & list)
    {
      for (typename std::vector<data* >::iterator it = list.begin(); it != list.end(); ++it)
        (*it)->quantity();
    }

    /** helper function to pack all data in a list before exchange, only valid with quantities that do work before parallel communication */
    template <typename data>
    void prepare_exchange_loop(std::map<std::string,data * > & list)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it) {
        (it->second)->prepare_exchange();
      }
    }

    /** helper function to extract all data in a list after exchange, only valid with quantities that do work after parallel communication */
    template <typename data>
    void post_exchange_loop(std::map<std::string,data * > & list)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it) {
        (it->second)->post_exchange();
      }
    }

    /** helper function to determine memory usage of all data in a list, only valid with PAQs */
    template <typename data>
    void memory_usage_loop(const std::map<std::string,data * > & list, int & usage) const
    {
      for (typename std::map<std::string,data* >::const_iterator it = list.begin(); it != list.end(); ++it)
        usage += (it->second)->memory_usage();
    }

    /** helper function to pack arrays of all data before exchange in a list, only valid with PAQs */
    template <typename data>
    void pack_exchange_loop(std::map<std::string,data * > & list, int & index, int i, double *buffer)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it) {
        index += (it->second)->pack_exchange(i,&buffer[index]);
      }
    }

    /** helper function to unpack arrays of all data after exchange in a list, only valid with PAQs */
    template <typename data>
    void unpack_exchange_loop(std::map<std::string,data * > & list, int & index, int i, double *buffer)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it)
        index += (it->second)->unpack_exchange(i,&buffer[index]);
    }

    /** helper function to pack arrays of all data in a list, only valid with comm lists */
    template <typename data>
    void pack_comm_loop(std::vector<data * > & list, int & size, int index, double *buf, 
                        int pbc_flag, int *pbc)
    {
      for (typename std::vector<data* >::iterator it = list.begin(); it != list.end(); ++it)
        size += (*it)->pack_comm(index,&buf[size],pbc_flag,pbc);
    }

    /** helper function to unpack arrays of all data in a list, only valid with comm lists */
    template <typename data>
    void unpack_comm_loop(std::vector<data * > & list, int & size, int index, double *buf)
    {
      for (typename std::vector<data* >::iterator it = list.begin(); it != list.end(); ++it)
        size += (*it)->unpack_comm(index,&buf[size]);
    }

    /** helper function to grow arrays of all data in a list, only valid with PAQs */
    template <typename data>
    void grow_arrays_loop(std::map<std::string,data * > & list, int nmax)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it)
        (it->second)->grow_lammps_array(nmax,prefix_+it->first);
    }

    /** helper function to copy arrays of all data in a list, only valid with PAQs */
    template <typename data>
    void copy_arrays_loop(std::map<std::string,data * > & list, int i, int j)
    {
      for (typename std::map<std::string,data* >::iterator it = list.begin(); it != list.end(); ++it)
        (it->second)->copy_lammps_array(i,j);
    };

  private:

    
    InterscaleManager();
  
  };

}

#endif
