#ifndef PER_ATOM_QUANTITY_INL_H
#define PER_ATOM_QUANTITY_INL_H

#include <string>

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class PerAtomQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  template <typename T>
  PerAtomQuantity<T>::PerAtomQuantity(ATC_Method * atc,
                                      int nCols,
                                      AtomType atomType) :
    MatrixDependencyManager<DenseMatrix, T>(),
    atc_(atc,atomType),
    lammpsInterface_(LammpsInterface::instance()),
    atomType_(atomType),
    nCols_(nCols),
    quantityToLammps_(atc_.atc_to_lammps_map()),
    lammpsScalar_(NULL),
    lammpsVector_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  template <typename T>
  PerAtomQuantity<T>::~PerAtomQuantity()
  {
    if (lammpsScalar_) lammpsInterface_->destroy_array(lammpsScalar_);
    if (lammpsVector_) lammpsInterface_->destroy_array(lammpsVector_);
  }

  //--------------------------------------------------------
  //  set_lammps_to_quantity
  //--------------------------------------------------------
  template <typename T>
  void PerAtomQuantity<T>::set_lammps_to_quantity() const
  {
    const DenseMatrix<T> & myQuantity(this->quantity_); // necessary to access quantity_ this way because of templating
    if (myQuantity.nRows()>0) {
      // full matrix copy
      
      if (atomType_ == ALL || atomType_ == PROC_GHOST) {
        if (nCols_==1) { // scalar
          T * lammpsQuantity = this->lammps_scalar();
          
          for (int i = 0; i < atc_.nlocal_total(); i++)
            lammpsQuantity[i] = myQuantity(i,0);
        }
        else{ // vector
          T ** lammpsQuantity = this->lammps_vector();
          
          for (int i = 0; i < atc_.nlocal_total(); i++)
            for (int j = 0; j < nCols_; j++)
              lammpsQuantity[i][j] = myQuantity(i,j);
        }
      }
      // map quantities
      else {
        int atomIndex;
        const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
        if (nCols_==1) { // scalar
          T * lammpsQuantity = this->lammps_scalar();
          for (int i = 0; i < myQuantity.nRows(); i++) {
            atomIndex = quantityToLammps(i);
            lammpsQuantity[atomIndex] = myQuantity(i,0);
          }
        }
        else{ // vector
          T ** lammpsQuantity = this->lammps_vector();
          for (int i = 0; i < myQuantity.nRows(); i++) {
            atomIndex = quantityToLammps(i);
            for (int j = 0; j < nCols_; j++) {
              lammpsQuantity[atomIndex][j] = myQuantity(i,j);
            }
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //  set_quantity_to_lammps
  //--------------------------------------------------------
  template <typename T>
  void PerAtomQuantity<T>::set_quantity_to_lammps() const
  {
    DenseMatrix<T> & myQuantity(this->quantity_);
    if (myQuantity.nRows()>0) {
      // full matrix copy
      // in the case where processor ghosts are in the quantity, don't set them back
      if (atomType_ == ALL || atomType_ == PROC_GHOST) {
        if (nCols_==1) { // scalar
          const T * lammpsQuantity = this->lammps_scalar();
          for (int i = 0; i < myQuantity.nRows(); i++) {
            myQuantity(i,0) = lammpsQuantity[i];
          }
        }
        else {
          const T * const * lammpsQuantity = this->lammps_vector();
          for (int i = 0; i < myQuantity.nRows(); i++) {
            for (int j = 0; j < nCols_; j++) {
              myQuantity(i,j) = lammpsQuantity[i][j];
            }
          }
        }
      }
      // map quantities
      else {
        const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
        int atomIndex;
        
        if (nCols_==1) { // scalar
          const T * lammpsQuantity = this->lammps_scalar();
          for (int i = 0; i < myQuantity.nRows(); i++) {
            atomIndex = quantityToLammps(i);
            myQuantity(i,0) = lammpsQuantity[atomIndex];
          }
        }
        else {
          const T * const * lammpsQuantity = this->lammps_vector();
          for (int i = 0; i < myQuantity.nRows(); i++) {
            atomIndex = quantityToLammps(i);
            for (int j = 0; j < nCols_; j++) {
              myQuantity(i,j) = lammpsQuantity[atomIndex][j];
            }
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for exchange with another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomQuantity<T>::pack_exchange(int i, double *buffer)
  {
    if (nCols_ == 1)
      buffer[0] = static_cast<double>(this->lammps_scalar()[i]);
    else
      for (int j = 0; j < nCols_; j++) {
        T ** lammpsVector = this->lammps_vector();
        buffer[j] = static_cast<double>(lammpsVector[i][j]);
      }
    return nCols_;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomQuantity<T>::unpack_exchange(int i, double *buffer)
  {
    if (nCols_ == 1)
      this->lammps_scalar()[i] = static_cast<T>(buffer[0]);
    else
      for (int j = 0; j < nCols_; j++) {
        T ** lammpsVector = this->lammps_vector();
        lammpsVector[i][j] = static_cast<T>(buffer[j]);
      }
    return nCols_;
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for passing to ghosts on another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomQuantity<T>::pack_comm(int index, double *buf, 
                                    int pbc_flag, int *pbc)
  {
    if (this->need_reset()) this->reset();
    DenseMatrix<T> & myQuantity(this->quantity_);
    for (int k = 0; k < nCols_; k++) {
      buf[k] = static_cast<double>(myQuantity(index,k));
    }
    return nCols_;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays for passing to ghosts on another proc
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomQuantity<T>::unpack_comm(int index, double *buf)
  {
    DenseMatrix<T> & myQuantity(this->quantity_);
    for (int k = 0; k < nCols_; k++) {
      myQuantity(index,k) = static_cast<T>(buf[k]);
    }
    this->propagate_reset();
    return nCols_;
  }

  //-----------------------------------------------------------------
  // allocate local atom-based arrays 
  //-----------------------------------------------------------------
  template <typename T>
  void PerAtomQuantity<T>::grow_lammps_array(int nmax, const std::string & tag)
  {
    
    if (nCols_ == 1)
      this->lammpsScalar_ = lammpsInterface_->grow_array(this->lammpsScalar_,nmax,tag.c_str());
    else
      this->lammpsVector_ = lammpsInterface_->grow_array(this->lammpsVector_,nmax,nCols_,tag.c_str());
  }

  //-----------------------------------------------------------------
  // copy values within local atom-based arrays 
  //-----------------------------------------------------------------
  template <typename T>
  void PerAtomQuantity<T>::copy_lammps_array(int i, int j)
  {
    if (nCols_ == 1) {
      T * lammpsScalar = this->lammps_scalar();
      lammpsScalar[j] = lammpsScalar[i];
    }
    else {
      T ** lammpsVector = this->lammps_vector();
      for (int k = 0; k < nCols_; k++)
        lammpsVector[j][k] = lammpsVector[i][k];
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LammpsAtomQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  reset
  //    reset is only true if set_quantity was used,
  //    so this syncs the quantity back with lammps
  //--------------------------------------------------------
  template <typename T>
  void LammpsAtomQuantity<T>::reset() const
  {
    if (this->need_reset()) {
      PerAtomQuantity<T>::reset();
      this->set_quantity_to_lammps();
    }
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for passing to ghosts on another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int LammpsAtomQuantity<T>::pack_comm(int index, double *buf, 
                                       int pbc_flag, int *pbc)
  {
    if (this->need_reset()) this->reset();
    int bufIdx = 0;

    if (this->nCols_ == 1) {
      T * lammpsQuantity = this->lammps_scalar();
      buf[bufIdx++] = double(lammpsQuantity[index]);
    }
    else {
      T ** lammpsQuantity = this->lammps_vector();
      for (int k = 0; k < this->nCols_; k++)
        buf[bufIdx++] = double(lammpsQuantity[index][k]);
    }
    return bufIdx;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays for passing to ghosts on another proc
  //-----------------------------------------------------------------
  template <typename T>
  int LammpsAtomQuantity<T>::unpack_comm(int index, double *buf)
  {
    DenseMatrix<T> & myQuantity(this->quantity_);
    int bufIdx = 0;

    if (this->nCols_ == 1) {
      T * lammpsQuantity = this->lammps_scalar();
      myQuantity(index,0) = T(buf[bufIdx++]);
      lammpsQuantity[index] = myQuantity(index,0);
    }
    else {
      T ** lammpsQuantity = this->lammps_vector();
      for (int k = 0; k < this->nCols_; k++) {
        myQuantity(index,k) = T(buf[bufIdx++]);
        lammpsQuantity[index][k] = myQuantity(index,k);
      }
    }
    
    this->propagate_reset();
    return bufIdx;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ProtectedMappedAtomQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  set_lammps_to_quantity
  //--------------------------------------------------------
  template <typename T>
  void ProtectedMappedAtomQuantity<T>::set_lammps_to_quantity() const
  {
    this->reset();
    const DenseMatrix<T> & myQuantity(this->quantity_); // necessary to access quantity_ this way because of templating
    int nCols = myQuantity.nCols();
    const INT_ARRAY & atomMap(atomMap_->quantity());
    if (myQuantity.nRows()>0) {
      // full matrix copy
      
      if (PerAtomQuantity<T>::atomType_ == ALL || PerAtomQuantity<T>::atomType_ == PROC_GHOST) {
        if (nCols==1) { // scalar
          T * lammpsQuantity = ProtectedAtomQuantity<T>::lammps_scalar();
          
          for (int i = 0; i < PerAtomQuantity<T>::atc_.nlocal_total(); i++) {
            int idx = atomMap(i,0);
            if (idx > -1) {
              lammpsQuantity[i] = myQuantity(idx,0);
            }
          }
        }
        else{ // vector
          T ** lammpsQuantity = ProtectedAtomQuantity<T>::lammps_vector();
          
          for (int i = 0; i < PerAtomQuantity<T>::atc_.nlocal_total(); i++) {
            int idx = atomMap(i,0);
            if (idx > -1) {
              for (int j = 0; j < nCols; j++) {
                lammpsQuantity[i][j] = myQuantity(idx,j);
              }
            }
          }
        }
      }
      // map quantities
      else {
        int atomIndex;
        const Array<int> & quantityToLammps = (PerAtomQuantity<T>::atc_).atc_to_lammps_map();
        if (nCols==1) { // scalar
          T * lammpsQuantity = ProtectedAtomQuantity<T>::lammps_scalar();
          for (int i = 0; i < atomMap.nRows(); i++) {
            int idx = atomMap(i,0);
            if (idx > -1) {
              atomIndex = quantityToLammps(i);
              lammpsQuantity[atomIndex] = myQuantity(idx,0);
            }
          }
        }
        else{ // vector
          T ** lammpsQuantity = ProtectedAtomQuantity<T>::lammps_vector();
          for (int i = 0; i < atomMap.nRows(); i++) {
            int idx = atomMap(i,0);
            if (idx > -1) {
              atomIndex = quantityToLammps(i);
              for (int j = 0; j < nCols; j++) {
                lammpsQuantity[atomIndex][j] = myQuantity(idx,j);
              }
            }
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //  set_quantity_to_lammps
  //--------------------------------------------------------
  template <typename T>
  void ProtectedMappedAtomQuantity<T>::set_quantity_to_lammps() const
  {
    DenseMatrix<T> & myQuantity(this->quantity_);
    int nCols = myQuantity.nCols();
    const INT_ARRAY & atomMap(atomMap_->quantity());
    if (myQuantity.nRows()>0) {
      // full matrix copy
      // in the case where processor ghosts are in the quantity, don't set them back
      if (PerAtomQuantity<T>::atomType_ == ALL || PerAtomQuantity<T>::atomType_ == PROC_GHOST) {
        if (nCols==1) { // scalar
          const T * lammpsQuantity = ProtectedAtomQuantity<T>::lammps_scalar();
          for (int i = 0; i < PerAtomQuantity<T>::atc_.nlocal_total(); i++) {
            int idx = atomMap(i,0);
            if (idx > -1) {
              myQuantity(idx,0) = lammpsQuantity[i];
            }
          }
        }
        else {
          const T * const * lammpsQuantity = ProtectedAtomQuantity<T>::lammps_vector();
          for (int i = 0; i < PerAtomQuantity<T>::atc_.nlocal_total(); i++) {
            int idx = atomMap(i,0);
            if (idx > -1) {
              for (int j = 0; j < nCols; j++)
                myQuantity(idx,j) = lammpsQuantity[i][j];
            }
          }
        }
      }
      // map quantities
      else {
        int atomIndex;
        const Array<int> & quantityToLammps = (PerAtomQuantity<T>::atc_).atc_to_lammps_map();
        if (nCols==1) { // scalar
          const T * lammpsQuantity = ProtectedAtomQuantity<T>::lammps_scalar();
          for (int i = 0; i < atomMap.nRows(); i++) {
            int idx = atomMap(i,0);
            if (idx > -1) {
              atomIndex = quantityToLammps(i);
              myQuantity(idx,0) = lammpsQuantity[atomIndex];
            }
          }
        }
        else {
          const T * const * lammpsQuantity = ProtectedAtomQuantity<T>::lammps_vector();
          for (int i = 0; i < atomMap.nRows(); i++) {
            int idx = atomMap(i,0);
            if (idx > -1) {
              atomIndex = quantityToLammps(i);
              for (int j = 0; j < nCols; j++) {
                myQuantity(idx,j) = lammpsQuantity[atomIndex][j];
              }
            }
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class PerAtomDiagonalMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  template <typename T>
  PerAtomDiagonalMatrix<T>::PerAtomDiagonalMatrix(ATC_Method * atc,
                                                  AtomType atomType) :
    MatrixDependencyManager<DiagonalMatrix, T>(),
    atc_(atc,atomType),
    lammpsInterface_(LammpsInterface::instance()),
    atomType_(atomType),
    quantityToLammps_(atc_.atc_to_lammps_map()),
    lammpsScalar_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  template <typename T>
  PerAtomDiagonalMatrix<T>::~PerAtomDiagonalMatrix()
  {
    if (lammpsScalar_) lammpsInterface_->destroy_array(lammpsScalar_);
  }

  //--------------------------------------------------------
  //  set_lammps_to_quantity
  //--------------------------------------------------------
  template <typename T>
  void PerAtomDiagonalMatrix<T>::set_lammps_to_quantity() const
  {
    const DiagonalMatrix<T> & myQuantity(this->quantity_); // necessary to access quantity_ this way because of templating
    if (myQuantity.size()>0) {
      // full matrix copy
      
      if (atomType_ == ALL || atomType_ == PROC_GHOST) {
        T * lammpsQuantity = this->lammps_scalar();
        
        for (int i = 0; i < atc_.nlocal_total(); i++) {
          lammpsQuantity[i] = myQuantity(i,i);
        }
      }
      // map quantities
      else {
        int atomIndex;
        const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
        T * lammpsQuantity = this->lammps_scalar();
        for (int i = 0; i < myQuantity.nRows(); i++) {
          atomIndex = quantityToLammps(i);
          lammpsQuantity[atomIndex] = myQuantity(i,i);
        }
      }
    }
  }

  //--------------------------------------------------------
  //  set_quantity_to_lammps
  //--------------------------------------------------------
  template <typename T>
  void PerAtomDiagonalMatrix<T>::set_quantity_to_lammps() const
  {
    DiagonalMatrix<T> & myQuantity(this->quantity_);
    if (myQuantity.size()>0) {
      // full matrix copy
      // in the case where processor ghosts are in the quantity, don't set them back
      if (atomType_ == ALL || atomType_ == PROC_GHOST) {
        const T * lammpsQuantity = this->lammps_scalar();
        for (int i = 0; i < myQuantity.nRows(); i++) {
          myQuantity(i,i) = lammpsQuantity[i];
        }
      }
      // map quantities
      else {
        int atomIndex;
        const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
        const T * lammpsQuantity = this->lammps_scalar();
        for (int i = 0; i < myQuantity.nRows(); i++) {
          atomIndex = quantityToLammps(i);
          myQuantity(i,i) = lammpsQuantity[atomIndex];
        }
      }
    }
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for exchange with another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomDiagonalMatrix<T>::pack_exchange(int i, double *buffer)
  {
    buffer[0] = static_cast<double>(lammps_scalar()[i]);
    return 1;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomDiagonalMatrix<T>::unpack_exchange(int i, double *buffer)
  {
    lammps_scalar()[i] = static_cast<T>(buffer[0]);
    return 1;
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for passing to ghosts on another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomDiagonalMatrix<T>::pack_comm(int index, double *buf, 
                                          int pbc_flag, int *pbc)
  {
    if (this->need_reset()) this->reset();
    DiagonalMatrix<T> & myQuantity(this->quantity_);
    buf[0] = static_cast<double>(myQuantity(index,index));
    return 1;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays for passing to ghosts on another proc
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomDiagonalMatrix<T>::unpack_comm(int index, double *buf)
  {
    DiagonalMatrix<T> & myQuantity(this->quantity_);
    myQuantity(index,index) = static_cast<T>(buf[0]);
    this->propagate_reset();
    return 1;
  }

  //-----------------------------------------------------------------
  // allocate local atom-based arrays 
  //-----------------------------------------------------------------
  template <typename T>
  void PerAtomDiagonalMatrix<T>::grow_lammps_array(int nmax, const std::string & tag)
  {
    this->lammpsScalar_ = lammpsInterface_->grow_array(this->lammpsScalar_,nmax,tag.c_str());
  }

  //-----------------------------------------------------------------
  // copy values within local atom-based arrays 
  //-----------------------------------------------------------------
  template <typename T>
  void PerAtomDiagonalMatrix<T>::copy_lammps_array(int i, int j)
  {
    T * lammpsScalar = this->lammps_scalar();
    lammpsScalar[j] = lammpsScalar[i];
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class PerAtomSparseMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  template <typename T>
  PerAtomSparseMatrix<T>::PerAtomSparseMatrix(ATC_Method * atc,
                                              int nCols,
                                              int maxEntriesPerRow,
                                              AtomType atomType) :
    MatrixDependencyManager<SparseMatrix, T>(),
    atc_(atc,atomType),
    lammpsInterface_(LammpsInterface::instance()),
    atomType_(atomType),
    nCols_(nCols),
    maxEntriesPerRow_(maxEntriesPerRow),
    quantityToLammps_(atc_.atc_to_lammps_map()),
    lammpsVector_(NULL),
    lammpsColIndices_(NULL)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  template <typename T>
  PerAtomSparseMatrix<T>::~PerAtomSparseMatrix()
  {
    if (lammpsVector_) lammpsInterface_->destroy_array(lammpsVector_);
    if (lammpsColIndices_) lammpsInterface_->destroy_array(lammpsColIndices_);
  }

  //--------------------------------------------------------
  //  set_lammps_to_quantity
  //--------------------------------------------------------
  template <typename T>
  void PerAtomSparseMatrix<T>::set_lammps_to_quantity() const
  {
    const SparseMatrix<T> & myQuantity(this->quantity_); // necessary to access quantity_ this way because of templating
    if (myQuantity.nRows()>0) {
      // full matrix copy
      if (atomType_ == ALL || atomType_ == PROC_GHOST) {
        
        T ** lammpsQuantity = this->lammps_vector();
        int ** lammpsColIndices = this->lammps_column_indices();
        
        for (int i = 0; i < atc_.nlocal_total(); i++) {
          myQuantity.row(i,_values_,_colIndices_);
          for (int j = 0; j < _values_.size(); j++) {
            lammpsQuantity[i][j] = _values_(j);
            lammpsColIndices[i][j] = _colIndices_(j);
          }
          for (int j = _values_.size(); j < maxEntriesPerRow_; j++) {
            lammpsColIndices[i][j] = -1;
          }
        }
      }
      // map quantities
      else {
        int atomIndex;
        const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
        T ** lammpsQuantity = this->lammps_vector();
        int ** lammpsColIndices = this->lammps_column_indices();
        for (int i = 0; i < myQuantity.nRows(); i++) {
          myQuantity.row(i,_values_,_colIndices_);
          atomIndex = quantityToLammps(i);
          for (int j = 0; j < _values_.size(); j++) {
            lammpsQuantity[atomIndex][j] = _values_(j);
            lammpsColIndices[atomIndex][j] = _colIndices_(j);
          }
          for (int j = _values_.size(); j < maxEntriesPerRow_; j++) {
            lammpsQuantity[atomIndex][j] = 0;
            lammpsColIndices[atomIndex][j] = -1;
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //  set_quantity_to_lammps
  //--------------------------------------------------------
  template <typename T>
  void PerAtomSparseMatrix<T>::set_quantity_to_lammps() const
  {
    SparseMatrix<T> & myQuantity(this->quantity_);
    if (myQuantity.nRows()>0) {
      // full matrix copy
      // in the case where processor ghosts are in the quantity, don't set them back
      if (atomType_ == ALL || atomType_ == PROC_GHOST) {
        const T * const * lammpsQuantity = this->lammps_vector();
        const int * const * lammpsColIndices = this->lammps_column_indices();
        for (int i = 0; i < myQuantity.nRows(); i++) {
          for (int j = 0; j < maxEntriesPerRow_; j++) {
            if (lammpsColIndices[i][j] < 0) {
              break;
            }
            myQuantity.set(i,lammpsColIndices[i][j],lammpsQuantity[i][j]);
          }
        }
      }
      // map quantities
      else {
        int atomIndex;
        const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
        const T * const * lammpsQuantity = this->lammps_vector();
        const int * const * lammpsColIndices = this->lammps_column_indices();
        for (int i = 0; i < myQuantity.nRows(); i++) {
          atomIndex = quantityToLammps(i);
          for (int j = 0; j < maxEntriesPerRow_; j++) {
            if (lammpsColIndices[atomIndex][j] < 0) {
              break;
            }
            myQuantity.set(i,lammpsColIndices[atomIndex][j],lammpsQuantity[atomIndex][j]);
          }
        }
      }

      myQuantity.compress();
    }
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for exchange with another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomSparseMatrix<T>::pack_exchange(int i, double *buffer)
  {
    int idx = 0;
    T ** lammpsVector = this->lammps_vector();
    for (int j = 0; j < maxEntriesPerRow_; j++) {
      buffer[idx++] = static_cast<double>(lammpsVector[i][j]);
    }
    int ** lammpsColIndices = this->lammps_column_indices();
    for (int j = 0; j < maxEntriesPerRow_; j++) {
      buffer[idx++] = static_cast<double>(lammpsColIndices[i][j]);
    }
    return 2*maxEntriesPerRow_;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomSparseMatrix<T>::unpack_exchange(int i, double *buffer)
  {
    int idx = 0;
    T ** lammpsVector = this->lammps_vector();
    for (int j = 0; j < maxEntriesPerRow_; j++) {
      lammpsVector[i][j] = static_cast<T>(buffer[idx++]);
    }
    int ** lammpsColIndices = this->lammps_column_indices();
    for (int j = 0; j < maxEntriesPerRow_; j++) {
      lammpsColIndices[i][j] = static_cast<int>(buffer[idx++]);
    }
    return 2*maxEntriesPerRow_;
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for passing to ghosts on another proc 
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomSparseMatrix<T>::pack_comm(int index, double *buf, 
                                        int pbc_flag, int *pbc)
  {
    return 0;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays for passing to ghosts on another proc
  //-----------------------------------------------------------------
  template <typename T>
  int PerAtomSparseMatrix<T>::unpack_comm(int index, double *buf)
  {
    return 0;
  }

  //-----------------------------------------------------------------
  // allocate local atom-based arrays 
  //-----------------------------------------------------------------
  template <typename T>
  void PerAtomSparseMatrix<T>::grow_lammps_array(int nmax, const std::string & tag)
  {
    
    this->lammpsVector_ = lammpsInterface_->grow_array(this->lammpsVector_,nmax,maxEntriesPerRow_,tag.c_str());
    std::string myString(tag+std::string("Columns"));
    this->lammpsColIndices_ = lammpsInterface_->grow_array(this->lammpsColIndices_,nmax,maxEntriesPerRow_,myString.c_str());
  }

  //-----------------------------------------------------------------
  // copy values within local atom-based arrays 
  //-----------------------------------------------------------------
  template <typename T>
  void PerAtomSparseMatrix<T>::copy_lammps_array(int i, int j)
  {
    T ** lammpsVector = this->lammps_vector();
    int ** lammpsColIndices = this->lammps_column_indices();
    for (int k = 0; k < maxEntriesPerRow_; k++) {
      lammpsVector[j][k] = lammpsVector[i][k];
      lammpsColIndices[j][k] = lammpsColIndices[i][k];
    }
  }

}
#endif
