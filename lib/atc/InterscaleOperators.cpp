// ATC transfer headers
#include "InterscaleOperators.h"
#include "PerAtomQuantity.h"
#include "TransferOperator.h"
#include "MoleculeSet.h"
#include "ATC_Method.h"
//#include <typeinfo>

using std::set;
using std::map;
using std::string;

namespace ATC{

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class InterscaleManager
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  InterscaleManager::InterscaleManager(ATC_Method * atc) :
    atc_(atc),
    initialized_(false),
    prefix_(":fix_atc:")
  {
    fundamentalAtomQuantities_.resize(NUM_ATOM_TYPES);
    for (unsigned int i = 0; i < NUM_ATOM_TYPES; i++) {
      fundamentalAtomQuantities_[i].resize(LammpsInterface::NUM_FUNDAMENTAL_ATOM_QUANTITIES);
      for (unsigned int j = 0; j < LammpsInterface::NUM_FUNDAMENTAL_ATOM_QUANTITIES; j++)
        fundamentalAtomQuantities_[i][j] = NULL;
    }
  }

  //----------------------------------------------------
  //  Set_lammps_data_prefix
  //--------------------------------------------------------
  void InterscaleManager::set_lammps_data_prefix()
  {
    prefix_ = (atc_->lammps_interface())->fix_id() + prefix_;
  }
        
  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  InterscaleManager::~InterscaleManager()
  {
    clear();
  }

  //--------------------------------------------------------
  //  clear
  //--------------------------------------------------------
  void InterscaleManager::clear()
  {
    // set all memory types to temporary
    for (unsigned int i = 0; i < NUM_ATOM_TYPES; i++) {
      for (unsigned int j = 0; j < LammpsInterface::NUM_FUNDAMENTAL_ATOM_QUANTITIES; j++) {
        if (fundamentalAtomQuantities_[i][j]) {
          fundamentalAtomQuantities_[i][j]->set_memory_type(TEMPORARY);
        }
      }
    }
    set_memory_temporary(perAtomQuantities_);
    set_memory_temporary(perAtomIntQuantities_);
    set_memory_temporary(perAtomDiagonalMatrices_);
    set_memory_temporary(perAtomSparseMatrices_);
    set_memory_temporary(pairMaps_);
    set_memory_temporary(denseMatrices_);
    set_memory_temporary(denseMatricesInt_);
    set_memory_temporary(denseMatricesBool_);
    set_memory_temporary(sparseMatrices_);
    set_memory_temporary(diagonalMatrices_);
    set_memory_temporary(vectorSparMat_);
    set_memory_temporary(setInt_);
    set_memory_temporary(smallMoleculeSets_);

    // clean up maps and vectors
    clear_temporary_data();
  }

  //--------------------------------------------------------
  //  clear
  //--------------------------------------------------------
  void InterscaleManager::clear_temporary_data()
  {
    deletionList_.clear();
    int listSize = fundamentalAtomQuantities_.size()+perAtomQuantities_.size()+perAtomIntQuantities_.size()+perAtomDiagonalMatrices_.size()+perAtomSparseMatrices_.size();
    listSize += pairMaps_.size()+denseMatrices_.size()+denseMatricesInt_.size()+denseMatricesBool_.size()+sparseMatrices_.size()+diagonalMatrices_.size()+vectorSparMat_.size()+setInt_.size()+smallMoleculeSets_.size();
    deletionList_.reserve(listSize);
    create_deletion_list();
    for (unsigned int i = 0; i < deletionList_.size(); i++) {
      if (deletionList_[i]) {
        delete deletionList_[i];
      }
    }
  }

  //--------------------------------------------------------
  //  create_deletion_list
  //--------------------------------------------------------
  void InterscaleManager::create_deletion_list()
  {
    // set all quantities to unfound
    dfs_prepare_loop(perAtomQuantities_);
    dfs_prepare_loop(perAtomIntQuantities_);
    dfs_prepare_loop(perAtomDiagonalMatrices_);
    dfs_prepare_loop(perAtomSparseMatrices_);
    dfs_prepare_loop(perAtomQuantities_);
    dfs_prepare_loop(pairMaps_);
    dfs_prepare_loop(denseMatrices_);
    dfs_prepare_loop(denseMatricesInt_);
    dfs_prepare_loop(denseMatricesBool_);
    dfs_prepare_loop(sparseMatrices_);
    dfs_prepare_loop(diagonalMatrices_);
    dfs_prepare_loop(vectorSparMat_);
    dfs_prepare_loop(setInt_);
    dfs_prepare_loop(smallMoleculeSets_);

    // perform dfs, special case for fundamental atom quantities
    int index = 0;
    for (unsigned int i = 0; i < NUM_ATOM_TYPES; i++) {
      for (unsigned int j = 0; j < LammpsInterface::NUM_FUNDAMENTAL_ATOM_QUANTITIES; j++) {
        if (fundamentalAtomQuantities_[i][j]) {
          index = dfs_visit(fundamentalAtomQuantities_[i][j],index);
          if ((fundamentalAtomQuantities_[i][j])->memory_type()==TEMPORARY) {
            fundamentalAtomQuantities_[i][j] = NULL;
          }
        }
      }
    }
    // dfs for everything else
    dfs_visit_loop(perAtomQuantities_,index);
    dfs_visit_loop(perAtomIntQuantities_,index);
    dfs_visit_loop(perAtomDiagonalMatrices_,index);
    dfs_visit_loop(perAtomSparseMatrices_,index);
    dfs_visit_loop(pairMaps_,index);
    dfs_visit_loop(denseMatrices_,index);
    dfs_visit_loop(denseMatricesInt_,index);
    dfs_visit_loop(denseMatricesBool_,index);
    dfs_visit_loop(sparseMatrices_,index);
    dfs_visit_loop(diagonalMatrices_,index);
    dfs_visit_loop(vectorSparMat_,index);
    dfs_visit_loop(setInt_,index);
    dfs_visit_loop(smallMoleculeSets_,index);
  }

  //--------------------------------------------------------
  //  dfs_visit
  //--------------------------------------------------------
  int InterscaleManager::dfs_visit(DependencyManager * quantity, const int index)
  {
    int myIndex = index;
    set<DependencyManager * >::iterator it;
    bool isTemporary = (quantity->memory_type()==TEMPORARY);
    
    for (it = (quantity->dependentQuantities_).begin(); it != (quantity->dependentQuantities_).end(); it++) {
      // make sure that if quantity isn't persistent, none of it's dependencies are
      if ((*it)->memory_type()==PERSISTENT && isTemporary) {
        throw ATC_Error("InterscaleManager::dfs_visit - a persistent quantity has a temporary dependency");
      }

      if (!((*it)->dfsFound_)) {
          myIndex = dfs_visit(*it,myIndex);
        }
    }

    quantity->dfsFound_ = true;
    if (isTemporary)
      deletionList_.push_back(quantity);
    return ++myIndex;
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void InterscaleManager::initialize()
  {
    // force all existing objects to reset
    for (unsigned int i = 0; i < NUM_ATOM_TYPES; i++) {
      for (unsigned int j = 0; j < LammpsInterface::NUM_FUNDAMENTAL_ATOM_QUANTITIES; j++) {
        if (fundamentalAtomQuantities_[i][j]) {
          fundamentalAtomQuantities_[i][j]->force_reset();
        }
      }
    }
    force_reset_loop(perAtomQuantities_);
    force_reset_loop(perAtomIntQuantities_);
    force_reset_loop(perAtomDiagonalMatrices_);
    force_reset_loop(perAtomSparseMatrices_);
    force_reset_loop(perAtomQuantities_);
    force_reset_loop(pairMaps_);
    force_reset_loop(denseMatrices_);
    force_reset_loop(denseMatricesInt_);
    force_reset_loop(denseMatricesBool_);
    force_reset_loop(sparseMatrices_);
    force_reset_loop(diagonalMatrices_);
    force_reset_loop(vectorSparMat_);
    force_reset_loop(setInt_);
    force_reset_loop(smallMoleculeSets_);
  }

  // access methods for atomic quantities
  //--------------------------------------------------------
  //  fundamental_atom_quantity
  //--------------------------------------------------------
  FundamentalAtomQuantity * InterscaleManager::fundamental_atom_quantity(LammpsInterface::FundamentalAtomQuantity id,
                                                                             AtomType atomType)
  {
    if (!fundamentalAtomQuantities_[atomType][id]) { // create a new one if it doesn't exist
      if (id == LammpsInterface::ATOM_MASS) {
        double * mass = LammpsInterface::instance()->atom_mass();
        if (mass)
          fundamentalAtomQuantities_[atomType][id] = new AtomMass(atc_,atomType);
        else
          fundamentalAtomQuantities_[atomType][id] = new FundamentalAtomQuantity(atc_,id,atomType);
      }
      else
        fundamentalAtomQuantities_[atomType][id] = new FundamentalAtomQuantity(atc_,id,atomType);
      fundamentalAtomQuantities_[atomType][id]->set_memory_type(PERSISTENT);
    }

    return fundamentalAtomQuantities_[atomType][id];
  }

  //--------------------------------------------------------
  //  per_atom_quantity
  //--------------------------------------------------------
  PerAtomQuantity<double> * InterscaleManager::per_atom_quantity(const string & tag)
  {
    return return_quantity(perAtomQuantities_,tag);
  }

  //--------------------------------------------------------
  //  per_atom_int_quantity
  //--------------------------------------------------------
  PerAtomQuantity<int> * InterscaleManager::per_atom_int_quantity(const string & tag)
  {
    return return_quantity(perAtomIntQuantities_,tag);
  }

  //--------------------------------------------------------
  //  per_atom_diagonal_matrix
  //--------------------------------------------------------
  PerAtomDiagonalMatrix<double> * InterscaleManager::per_atom_diagonal_matrix(const string & tag)
  {
    return return_quantity(perAtomDiagonalMatrices_,tag);
  }

  //--------------------------------------------------------
  //  per_atom_sparse_matrix
  //--------------------------------------------------------
  PerAtomSparseMatrix<double> * InterscaleManager::per_atom_sparse_matrix(const string & tag)
  {
    return return_quantity(perAtomSparseMatrices_,tag);
  }

  //--------------------------------------------------------
  //  pair_map
  //--------------------------------------------------------
  PairMap * InterscaleManager::pair_map(const string & tag)
  {
    return return_quantity(pairMaps_,tag);
  }

  //--------------------------------------------------------
  //  add_per_atom_quantity
  //--------------------------------------------------------
  void InterscaleManager::add_per_atom_quantity(PerAtomQuantity<double> * atomQuantity,
                                                const string & tag)
  {
    add_comm_quantity(perAtomQuantities_,commList_,atomQuantity,tag);
  }

  //--------------------------------------------------------
  //  add_per_atom_int_quantity
  //--------------------------------------------------------
  void InterscaleManager::add_per_atom_int_quantity(PerAtomQuantity<int> * atomQuantity,
                                                   const string & tag)
  {
     add_comm_quantity(perAtomIntQuantities_,commIntList_,atomQuantity,tag);
  }

  //--------------------------------------------------------
  //  add_per_atom_diagonal_matrix
  //--------------------------------------------------------
  void InterscaleManager::add_per_atom_diagonal_matrix(PerAtomDiagonalMatrix<double> * atomQuantity,
                                                       const string & tag)
  {
     add_comm_quantity(perAtomDiagonalMatrices_,commDmList_,atomQuantity,tag);
  }

  //--------------------------------------------------------
  //  add_per_atom_sparse_matrix
  //--------------------------------------------------------
  void InterscaleManager::add_per_atom_sparse_matrix(PerAtomSparseMatrix<double> * atomQuantity,
                                                     const string & tag)
  {
    add_comm_quantity(perAtomSparseMatrices_,commSmList_,atomQuantity,tag);
  }

  //--------------------------------------------------------
  //  add_pair_map
  //--------------------------------------------------------
  void InterscaleManager::add_pair_map(PairMap * pairMap,
                                       const string & tag)
  {
    add_quantity(pairMaps_,pairMap,tag);
  }

  //--------------------------------------------------------
  //  dense_matrix
  //--------------------------------------------------------
  DENS_MAN * InterscaleManager::dense_matrix(const string & tag)
  {
    return return_quantity(denseMatrices_,tag);
  }

  //--------------------------------------------------------
  //  add_dense_matrix
  //--------------------------------------------------------
  void InterscaleManager::add_dense_matrix(DENS_MAN * denseMatrix,
                                           const string & tag)
  {
    add_quantity(denseMatrices_,denseMatrix,tag);
  }

  //--------------------------------------------------------
  //  dense_matrix_int
  //--------------------------------------------------------
  MatrixDependencyManager<DenseMatrix, int> * InterscaleManager::dense_matrix_int(const string & tag)
  {
    return return_quantity(denseMatricesInt_,tag);
  }

  //--------------------------------------------------------
  //  add_dense_matrix_int
  //--------------------------------------------------------
  void InterscaleManager::add_dense_matrix_int(MatrixDependencyManager<DenseMatrix, int> * denseMatrix,
                                               const string & tag)
  {
    add_quantity(denseMatricesInt_,denseMatrix,tag);
  }

  //--------------------------------------------------------
  //  dense_matrix_bool
  //--------------------------------------------------------
  MatrixDependencyManager<DenseMatrix, bool> * InterscaleManager::dense_matrix_bool(const string & tag)
  {
    return return_quantity(denseMatricesBool_,tag);
  }

  //--------------------------------------------------------
  //  add_dense_matrix_bool
  //--------------------------------------------------------
  void InterscaleManager::add_dense_matrix_bool(MatrixDependencyManager<DenseMatrix, bool> * denseMatrix,
                                                const string & tag)
  {
    add_quantity(denseMatricesBool_,denseMatrix,tag);
  }

  //--------------------------------------------------------
  //  sparse_matrix
  //--------------------------------------------------------
  SPAR_MAN * InterscaleManager::sparse_matrix(const string & tag)
  {
    return return_quantity(sparseMatrices_,tag);
  }

  //--------------------------------------------------------
  //  add_sparse_matrix
  //--------------------------------------------------------
  void InterscaleManager::add_sparse_matrix(SPAR_MAN * sparseMatrix,
                                            const string & tag)
  {
    add_quantity(sparseMatrices_,sparseMatrix,tag);
  }

  //--------------------------------------------------------
  //  diagonal_matrix
  //--------------------------------------------------------
  DIAG_MAN * InterscaleManager::diagonal_matrix(const string & tag)
  {
    return return_quantity(diagonalMatrices_,tag);
  }

  //--------------------------------------------------------
  //  add_sparse_matrix
  //--------------------------------------------------------
  void InterscaleManager::add_diagonal_matrix(DIAG_MAN * diagonalMatrix,
                                              const string & tag)
  {
    add_quantity(diagonalMatrices_,diagonalMatrix,tag);
  }

  //--------------------------------------------------------
  //  vector_spar_mat
  //--------------------------------------------------------
  VectorDependencyManager<SPAR_MAT * > * InterscaleManager::vector_sparse_matrix(const string & tag)
  {
    return return_quantity(vectorSparMat_,tag);
  }

  //--------------------------------------------------------
  //  add_vector_spar_mat
  //--------------------------------------------------------
  void InterscaleManager::add_vector_sparse_matrix(VectorDependencyManager<SPAR_MAT * > * vectorSparMat,
                                                   const string & tag)
  {
    add_quantity(vectorSparMat_,vectorSparMat,tag);
  }

  //--------------------------------------------------------
  //  set_int
  //--------------------------------------------------------
  SetDependencyManager<int> * InterscaleManager::set_int(const string & tag)
  {
   return return_quantity(setInt_,tag);
  }

  //--------------------------------------------------------
  //  add_set_int
  //--------------------------------------------------------
  void InterscaleManager::add_set_int(SetDependencyManager<int> * setInt,
                                      const string & tag)
  {
    add_quantity(setInt_,setInt,tag);
  }

  //--------------------------------------------------------
  //  molecule_set
  //--------------------------------------------------------
  SmallMoleculeSet * InterscaleManager::small_molecule_set(const string & tag)
  {
    return return_quantity(smallMoleculeSets_,tag);
  }

  //--------------------------------------------------------
  //  add_molecule_set
  //--------------------------------------------------------
  void InterscaleManager::add_small_molecule_set(SmallMoleculeSet * moleculeSet,
                                                 const string & tag)
  {
    add_quantity(smallMoleculeSets_,moleculeSet,tag);
    moleculeSet->initialize();
  }

  //--------------------------------------------------------
  //  find_tag
  //--------------------------------------------------------
  DependencyManager * InterscaleManager::find(const string & tag)
  {
    // REFACTOR add check for duplicate entries
    DependencyManager * quantity = NULL;
    
    quantity = find_in_list(perAtomQuantities_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(perAtomIntQuantities_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(perAtomDiagonalMatrices_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(perAtomSparseMatrices_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(pairMaps_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(denseMatrices_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(denseMatricesInt_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(denseMatricesBool_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(sparseMatrices_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(diagonalMatrices_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(vectorSparMat_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(setInt_,tag);
    if (quantity) return quantity;
    quantity = find_in_list(smallMoleculeSets_,tag);
    if (quantity) return quantity;

    return NULL;
  }

  //--------------------------------------------------------
  //  remove
  //--------------------------------------------------------
  void InterscaleManager::remove(const string & tag)
  {
    DependencyManager * toBeDeleted = find(tag);
    if (toBeDeleted) {
      toBeDeleted->set_memory_type(TEMPORARY);
    }
  }

  //--------------------------------------------------------
  //  reset_nlocal
  //--------------------------------------------------------
  void InterscaleManager::reset_nlocal()
  {
    reset_nlocal_loop(perAtomSparseMatrices_); // this goes in between diagonal matrices and molecule sets above
  }

  //--------------------------------------------------------
  //  computes_force_reset
  //--------------------------------------------------------
  void InterscaleManager::fundamental_force_reset(unsigned quantity)
  {
    for (unsigned int i = 0; i < NUM_ATOM_TYPES; i++) {
      if (fundamentalAtomQuantities_[i][quantity]) {
        fundamentalAtomQuantities_[i][quantity]->lammps_force_reset();
      }
    }
  }

  //--------------------------------------------------------
  //  computes_force_reset
  //--------------------------------------------------------
  void InterscaleManager::lammps_force_reset()
  {
    for (unsigned int j = 0; j < LammpsInterface::NUM_FUNDAMENTAL_ATOM_QUANTITIES; j++) {
      fundamental_force_reset(j);
    }
    lammps_reset_loop(perAtomQuantities_);
    lammps_reset_loop(perAtomIntQuantities_);
    lammps_reset_loop(perAtomDiagonalMatrices_);
    lammps_reset_loop(perAtomSparseMatrices_);
  }

  //--------------------------------------------------------
  //  size_comm_quantities
  //--------------------------------------------------------
  void InterscaleManager::size_comm_quantities()
  {
    size_comm_loop(commList_);
    size_comm_loop(commIntList_);
    size_comm_loop(commDmList_);
    size_comm_loop(commSmList_);
  }

  //--------------------------------------------------------
  //  prepare_exchange
  //--------------------------------------------------------
  void InterscaleManager::prepare_exchange()
  {
    prepare_exchange_loop(perAtomIntQuantities_);
    prepare_exchange_loop(perAtomQuantities_);
    prepare_exchange_loop(perAtomDiagonalMatrices_);
    prepare_exchange_loop(perAtomSparseMatrices_);
  }

  //--------------------------------------------------------
  //  post_exchange
  //--------------------------------------------------------
  void InterscaleManager::post_exchange()
  {
    post_exchange_loop(perAtomIntQuantities_);
    post_exchange_loop(perAtomQuantities_);
    post_exchange_loop(perAtomDiagonalMatrices_);
    post_exchange_loop(perAtomSparseMatrices_);
    post_exchange_loop(smallMoleculeSets_);
    post_exchange_loop(pairMaps_);
  }

  //--------------------------------------------------------
  //  memory_usage
  //--------------------------------------------------------
  int InterscaleManager::memory_usage() const
  {
    int usage = 0;

    memory_usage_loop(perAtomQuantities_,usage);
    memory_usage_loop(perAtomIntQuantities_,usage);
    memory_usage_loop(perAtomDiagonalMatrices_,usage);
    memory_usage_loop(perAtomSparseMatrices_,usage);

    return usage;
  }

  //--------------------------------------------------------
  //  pack_exchange
  //--------------------------------------------------------
  int InterscaleManager::pack_exchange(int i, double *buffer)
  {
    int index = 0;

    pack_exchange_loop(perAtomQuantities_,index,i,&buffer[index]);
    pack_exchange_loop(perAtomIntQuantities_,index,i,&buffer[index]);
    pack_exchange_loop(perAtomDiagonalMatrices_,index,i,&buffer[index]);
    pack_exchange_loop(perAtomSparseMatrices_,index,i,&buffer[index]);

    return index;
  }

  //--------------------------------------------------------
  //  unpack_exchange
  //--------------------------------------------------------
  int InterscaleManager::unpack_exchange(int i, double *buffer)
  {
    int index = 0;

    unpack_exchange_loop(perAtomQuantities_,index,i,&buffer[index]);
    unpack_exchange_loop(perAtomIntQuantities_,index,i,&buffer[index]);
    unpack_exchange_loop(perAtomDiagonalMatrices_,index,i,&buffer[index]);
    unpack_exchange_loop(perAtomSparseMatrices_,index,i,&buffer[index]);

    return index;
  }

  //--------------------------------------------------------
  //  pack_comm
  //--------------------------------------------------------
  int InterscaleManager::pack_comm(int index, double *buf, 
                                   int pbc_flag, int *pbc)
  {
    int size = 0;

    //pack_comm_loop(commList_,size,index,buf,pbc_flag,pbc);
    pack_comm_loop(commIntList_,size,index,buf,pbc_flag,pbc);
    //pack_comm_loop(commDmList_,size,index,buf,pbc_flag,pbc);
    //pack_comm_loop(commSmList_,size,index,buf,pbc_flag,pbc);

    return size;
  }

  //--------------------------------------------------------
  //  unpack_comm
  //--------------------------------------------------------
  int InterscaleManager::unpack_comm(int index, double *buf)
  {
    int size = 0;

    //unpack_comm_loop(commList_,size,index,buf);
    unpack_comm_loop(commIntList_,size,index,buf);
    //unpack_comm_loop(commDmList_,size,index,buf);
    //unpack_comm_loop(commSmList_,size,index,buf);

    return size;
  }

  //--------------------------------------------------------
  //  grow_array
  //--------------------------------------------------------
  void InterscaleManager::grow_arrays(int nmax)
  {
    grow_arrays_loop(perAtomQuantities_,nmax);
    grow_arrays_loop(perAtomIntQuantities_,nmax);
    grow_arrays_loop(perAtomDiagonalMatrices_,nmax);
    grow_arrays_loop(perAtomSparseMatrices_,nmax);
  }

  //--------------------------------------------------------
  //  copy_array
  //--------------------------------------------------------
  void InterscaleManager::copy_arrays(int i, int j)
  {
    copy_arrays_loop(perAtomQuantities_,i,j);
    copy_arrays_loop(perAtomIntQuantities_,i,j);
    copy_arrays_loop(perAtomDiagonalMatrices_,i,j);
    copy_arrays_loop(perAtomSparseMatrices_,i,j);
  }

};
