// ATC_Transfer headers
#include "ATC_Transfer.h"
#include "FE_Engine.h"
#include "Array.h"
#include "Array2D.h"
#include "ATC_Error.h"
#include "CG.h"
#include "XT_Function.h"
#include "PrescribedDataManager.h"
#include "LammpsInterface.h"

#include <fstream>

namespace ATC {

  //=================================================================
  // quadrature weight functions
  //=================================================================
  
  //-----------------------------------------------------------------
  void ATC_Transfer::set_atomic_weights(void)
  {
    switch (atomWeightType_) {
    case USER:
      reset_atomicWeightsGroup();
      break;
    case LATTICE:
      reset_atomicWeightsLattice();
      break;
    case ELEMENT:
      reset_atomicWeightsElement();
      break;
    case REGION:
      if (!initialized_) {
        // Uses region bounding box as total volume
        // note will *only* work if atoms exactly fill the region
        int nregion = lammpsInterface_->nregion();
        Array<double> volume(nregion);
    
        for (int i = 0; i < nregion; ++i) {
          volume(i) = 
            (lammpsInterface_->region_xhi(i)-lammpsInterface_->region_xlo(i))
            *(lammpsInterface_->region_yhi(i)-lammpsInterface_->region_ylo(i))
            *(lammpsInterface_->region_zhi(i)-lammpsInterface_->region_zlo(i));
        }
      
        DenseVector<int> localCount(nregion);
        DenseVector<int> globalCount(nregion);
      
        // loop over atoms
        localCount = 0;
        for (int i = 0; i < nLocal_; ++i) {
          for (int j = 0; j < nregion; ++j) {
            if (lammpsInterface_->region_match(j,
                                               atomicCoords_(0,i),
                                               atomicCoords_(1,i),
                                               atomicCoords_(2,i))) {
              localCount(j)++;
            }
          }
        }
        // communication to get total atom counts per group
        lammpsInterface_->int_allsum(localCount.get_ptr(),
                                     globalCount.get_ptr(),nregion);
        
        for (int i = 0; i < nregion; ++i) {
          if (globalCount(i) > 0)
            Valpha_[i] = volume(i)/globalCount(i);
          else
            Valpha_[i] = 0;
        }
      }
    
      reset_atomicWeightsRegion();
      break;
    case GROUP:
      if (!initialized_) {
        // Uses group bounding box as total volume
        // note will *only* work if atoms exactly fill the box
        int ngroup = lammpsInterface_->ngroup();
        int *mask = lammpsInterface_->atom_mask();
        double * bounds;
        map<int, double> volume;
      
        bounds = new double[6];
        for (int i = 0; i < ngroup; ++i) {
          lammpsInterface_->group_bounds(i, bounds);
          volume[lammpsInterface_->group_bit(i)] = 
            (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]);
        }
        delete [] bounds;
      
        DenseVector<int> localCount(ngroup);
        DenseVector<int> globalCount(ngroup);
      
        // loop over atoms
        localCount = 0;
        for (int i = 0; i < nLocal_; ++i) {
          for (int j = 0; j < ngroup; ++j) {
            if (mask[internalToAtom_(i)] & lammpsInterface_->group_bit(j))
              localCount(j)++;
          }
        }
      
        // communication to get total atom counts per group
        lammpsInterface_->int_allsum(localCount.get_ptr(),
                                     globalCount.get_ptr(),ngroup);
     
        for (int i = 0; i < ngroup; ++i) {
          int iGroupBit = lammpsInterface_->group_bit(i);
          if (globalCount(i) > 0)
            Valpha_[iGroupBit] = 
              volume[iGroupBit]/globalCount(i);
          else
            Valpha_[iGroupBit] = 0;
        }
      }
      reset_atomicWeightsGroup();
      break;
    case MULTISCALE:
      reset_atomicWeightsMultiscale(shpFcn_,atomicWeights_);
      break;
    }

    if( nLocalMask_ > 0 ) {
      atomicWeightsMask_.reset(nLocalMask_,nLocalMask_);
      int nmask = 0;
      for (int i = 0; i < nLocal_; ++i) {
        if (atomMask_(i)) {
          atomicWeightsMask_(nmask,nmask) = atomicWeights_(i,i);
          nmask++;
        }
      }
    }
  }

  //-----------------------------------------------------------------
  // check to see if shape function contains internal or ghost atoms
  //-----------------------------------------------------------------
  int ATC_Transfer::check_shape_function_type(int nodeIdx)
  {
   ShapeFunctionType thisShpType;
    // NOTE change this check to not use shpWeight
    if (shpWeight_(nodeIdx,nodeIdx) > 0.) { // NOTE need tolerance check
      DENS_VEC tmpWeights(nNodes_);
      DENS_VEC weights(nNodes_);
      if (nLocalGhost_>0)
        tmpWeights = shpFcnGhost_.col_sum();
      lammpsInterface_->allsum(tmpWeights.get_ptr(),weights.get_ptr(),nNodes_);
      
      if (weights(nodeIdx) > 0.) // NOTE need tolerance check
        thisShpType = BOUNDARY;
      else
        thisShpType = MD_ONLY;
    }
    else
      thisShpType = FE_ONLY;
    
    return thisShpType;
  }

  //-----------------------------------------------------------------
  // checks to see if an element is completely inside the MD domain
  //-----------------------------------------------------------------
  bool ATC_Transfer::check_internal(int eltIdx)
  {
    // NOTE we use a check only on the bounding boxes, not the actual shapes
    int nNode;
    DENS_MAT coords;
    bool notInMD = false;

    feEngine_->element_coordinates(eltIdx,coords);
    nNode = coords.nCols();

    if (regionID_ < 0) { // loop over groups to find bounds
      double xhi,xlo,yhi,ylo,zhi,zlo;
      double * tmpBounds = new double[6];
      set<int>::iterator iset; 
      iset = igroups_.begin(); // NOTE inefficient, use some bounds as initial values
      lammpsInterface_->group_bounds(*iset,tmpBounds);
      xlo = tmpBounds[0];
      xhi = tmpBounds[1];
      ylo = tmpBounds[2];
      yhi = tmpBounds[3];
      zlo = tmpBounds[4];
      zhi = tmpBounds[5];
      for (iset = igroups_.begin(); iset != igroups_.end(); iset++) {
        lammpsInterface_->group_bounds(*iset,tmpBounds);
        if (xlo>tmpBounds[0]) xlo = tmpBounds[0];
        if (xhi<tmpBounds[1]) xhi = tmpBounds[1];
        if (ylo>tmpBounds[2]) ylo = tmpBounds[2];
        if (yhi<tmpBounds[3]) yhi = tmpBounds[3];
        if (zlo>tmpBounds[4]) zlo = tmpBounds[4];
        if (zhi<tmpBounds[5]) zhi = tmpBounds[5];
      }
      bool checkx = true;
      bool checky = true;
      bool checkz = true;
      if (lammpsInterface_->xperiodic() == 1) checkx = false;
      if (lammpsInterface_->yperiodic() == 1) checky = false;
      if (lammpsInterface_->zperiodic() == 1) checkz = false;
      for (int i = 0; i < nNode; ++i) {
        if (((coords(0,i)<xlo)||(coords(0,i)>xhi))&&checkx) notInMD = true;
        if (((coords(1,i)<ylo)||(coords(1,i)>yhi))&&checky) notInMD = true;
        if (((coords(2,i)<zlo)||(coords(2,i)>zhi))&&checkz) notInMD = true;
      }
      delete [] tmpBounds;
    }
    else { // use region denoted by regionID_
      for (int i = 0; i < nNode; ++i) {
        if (!lammpsInterface_->region_match(regionID_,coords(0,i),coords(1,i),coords(2,i))) {
          notInMD = true;
        }
      }
    }
    
    return notInMD;
  }
  
  //-----------------------------------------------------------------
  // checks to see if an element intersects the ghost atoms
  //-----------------------------------------------------------------
  bool ATC_Transfer::intersect_ghost(int eltIdx)
  {
    // NOTE we use a check only on the bounding boxes, not the actual shapes
    int nNode;
    DENS_MAT coords;
    bool intersectGhost = false;

    feEngine_->element_coordinates(eltIdx,coords);
    nNode = coords.nCols();

    double xhi,xlo,yhi,ylo,zhi,zlo;
    double * tmpBounds;
    tmpBounds = new double[6];
    set<int>::iterator iset; 
    iset = igroupsGhost_.begin(); // NOTE inefficient, use some bounds as initial values
    lammpsInterface_->group_bounds(*iset,tmpBounds);
    xlo = tmpBounds[0];
    xhi = tmpBounds[1];
    ylo = tmpBounds[2];
    yhi = tmpBounds[3];
    zlo = tmpBounds[4];
    zhi = tmpBounds[5];
    for (iset = igroupsGhost_.begin(); iset != igroupsGhost_.end(); iset++) {
      lammpsInterface_->group_bounds(*iset,tmpBounds);
      if (xlo<tmpBounds[0]) xlo = tmpBounds[0];
      if (xhi>tmpBounds[1]) xhi = tmpBounds[1];
      if (ylo<tmpBounds[2]) ylo = tmpBounds[2];
      if (yhi>tmpBounds[3]) yhi = tmpBounds[3];
      if (zlo<tmpBounds[4]) zlo = tmpBounds[4];
      if (zhi>tmpBounds[5]) zhi = tmpBounds[5];
    }
    
    bool xlower = true;
    bool xhigher = true;
    bool ylower = true;
    bool yhigher = true;
    bool zlower = true;
    bool zhigher = true;
    for (int i = 0; i < nNode; ++i) {
      if (coords(0,i)>=xlo)
        xlower = false;
      if (coords(0,i)<=xhi)
        xhigher = false;
      if (coords(0,i)>=ylo)
        ylower = false;
      if (coords(0,i)<=yhi)
        yhigher = false;
      if (coords(0,i)>=zlo)
        zlower = false;
      if (coords(0,i)<=zhi)
        zhigher = false;
    }
    if (!(xlower||xhigher||ylower||yhigher||zlower||zhigher))
      intersectGhost = true;
    delete [] tmpBounds;

    return intersectGhost;
  }

  //=================================================================
  // memory management and processor information exchange
  //=================================================================

  //-----------------------------------------------------------------
  // memory usage of local atom-based arrays 
  //-----------------------------------------------------------------
  int ATC_Transfer::memory_usage()
  {
    int bytes = lammpsInterface_->nmax() * 3 * sizeof(double);
    return bytes;
  }

  //-----------------------------------------------------------------
  // allocate local atom-based arrays 
  //-----------------------------------------------------------------
  void ATC_Transfer::grow_arrays(int nmax)
  {
    xref_ =
      lammpsInterface_->grow_2d_double_array(xref_,nmax,3,"fix_atc:xref");
  }

  //-----------------------------------------------------------------
  // copy values within local atom-based arrays 
  //-----------------------------------------------------------------
  void ATC_Transfer::copy_arrays(int i, int j)
  {
    xref_[j][0] = xref_[i][0];
    xref_[j][1] = xref_[i][1];
    xref_[j][2] = xref_[i][2];
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays for exchange with another proc 
  //-----------------------------------------------------------------
  int ATC_Transfer::pack_exchange(int i, double *buf)
  {
    buf[0] = xref_[i][0]; 
    buf[1] = xref_[i][1]; 
    buf[2] = xref_[i][2]; 
    return 3;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  int ATC_Transfer::unpack_exchange(int nlocal, double *buf)
  {
    xref_[nlocal][0] = buf[0];
    xref_[nlocal][1] = buf[1];
    xref_[nlocal][2] = buf[2];
    return 3;
  }

  //-----------------------------------------------------------------
  // pack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  int ATC_Transfer::pack_comm(int n, int *list, double *buf, 
                              int pbc_flag, int *pbc)
  {
    int i,j,m;
    double dx,dy,dz;
  
    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = xref_[j][0];
        buf[m++] = xref_[j][1];
        buf[m++] = xref_[j][2];
      }
    } 
    else {
      if (lammpsInterface_->domain_triclinic() == 0) {
        dx = pbc[0]*Xprd_;
        dy = pbc[1]*Yprd_;
        dz = pbc[2]*Zprd_;
      } 
      else {
        dx = pbc[0]*Xprd_ + pbc[5]*XY_ + pbc[4]*XZ_;
        dy = pbc[1]*Yprd_ + pbc[3]*YZ_;
        dz = pbc[2]*Zprd_;
      }
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = xref_[j][0] + dx;
        buf[m++] = xref_[j][1] + dy;
        buf[m++] = xref_[j][2] + dz;
      }
    }
    return 3;
  }

  //-----------------------------------------------------------------
  // unpack values in local atom-based arrays from exchange with another proc 
  //-----------------------------------------------------------------
  void ATC_Transfer::unpack_comm(int n, int first, double *buf)
  {
    int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      xref_[i][0] = buf[m++];
      xref_[i][1] = buf[m++];
      xref_[i][2] = buf[m++];
    }
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  bool ATC_Transfer::get_region_bounds(const char * regionName,
                                       double & xmin, double & xmax,
                                       double & ymin, double & ymax,
                                       double & zmin, double & zmax,
                                       double & xscale,
                                       double & yscale,
                                       double & zscale)
  {
    int iRegion = lammpsInterface_->get_region_id(regionName);

    xscale = lammpsInterface_->region_xscale(iRegion);
    yscale = lammpsInterface_->region_yscale(iRegion);
    zscale = lammpsInterface_->region_zscale(iRegion);
    xmin = lammpsInterface_->region_xlo(iRegion); 
    xmax = lammpsInterface_->region_xhi(iRegion);
    ymin = lammpsInterface_->region_ylo(iRegion); 
    ymax = lammpsInterface_->region_yhi(iRegion);
    zmin = lammpsInterface_->region_zlo(iRegion); 
    zmax = lammpsInterface_->region_zhi(iRegion);

    if (strcmp(lammpsInterface_->region_style(iRegion),"block")==0) {
      return true;
    }
    else {
      return false;
    }
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Transfer::reset_nlocal()
  {
    nLocalTotal_ = lammpsInterface_->nlocal();
    int * mask = lammpsInterface_->atom_mask();
    nLocal_ = 0;
    nLocalGhost_ = 0;
    for (int i = 0; i < nLocalTotal_; ++i) {
      if (mask[i] & groupbit_) nLocal_++;
      if (mask[i] & groupbitGhost_) nLocalGhost_++;
    }
    int local_data[2] = {nLocal_, nLocalGhost_};
    int data[2] = {0, 0};
    lammpsInterface_->int_allsum(local_data,data,2);
    nInternal_ = data[0];
    nGhost_    = data[1];

    // set up internal & ghost maps & coordinates
    if (nLocal_>0) {
      // set map
      internalToAtom_.reset(nLocal_);
      int j = 0;
      // construct internalToAtom map 
      //  : internal index -> local lammps atom index
      for (int i = 0; i < nLocalTotal_; ++i) {
        if (mask[i] & groupbit_) internalToAtom_(j++) = i;
      }
      // construct reverse map
      for (int i = 0; i < nLocal_; ++i) {
        atomToInternal_[internalToAtom_(i)] = i;
      }
      atomicCoords_.reset(nsd_,nLocal_);
    }
    if (nLocalGhost_>0) {
      // set map
      ghostToAtom_.reset(nLocalGhost_);
      int j = 0;
      for (int i = 0; i < nLocalTotal_; ++i) {
        if (mask[i] & groupbitGhost_) ghostToAtom_(j++) = i;
      }
      ghostAtomCoords_.reset(nsd_,nLocalGhost_);
    }

    reset_coordinates();
  }

  void ATC_Transfer::reset_coordinates()
  {
    // reset atomic coordinates
    double ** xx = xref_;
    if (atomToElementMapType_ == EULERIAN) xx = lammpsInterface_->xatom();
    
    // fill coordinates
    if (nLocal_ > 0) {
      for (int i = 0; i < nLocal_; i++) {
        for (int j = 0; j < nsd_; j++)
          atomicCoords_(j,i) = xx[internalToAtom_(i)][j];
      }
    }
    
    // fill ghost coordinates
    if (nLocalGhost_ > 0) {
      for (int i = 0; i < nLocalGhost_; i++) {
        for (int j = 0; j < nsd_; j++)
          ghostAtomCoords_(j,i) = xx[ghostToAtom_(i)][j];
      }
    }
  }

  // shape functions are based on reference positions 
  // and only change when atoms cross processor boundaries
  // (x_ref migrates with atoms but does not change value)
  void ATC_Transfer::reset_shape_functions()
  {
    // FE_Engine allocates all required memory
    // assume initial atomic position is the reference position for now
    if (!initialized_) {
      shpFcnDerivs_ = vector<SPAR_MAT >(nsd_);
      shpFcnDerivsGhost_ = vector<SPAR_MAT >(nsd_);
      shpFcnDerivsMask_ = vector<SPAR_MAT >(nsd_);
    }

    // use FE_Engine to compute shape functions at atomic points
    if (nLocal_>0) {
      feEngine_->evaluate_shape_functions(atomicCoords_,
                                          shpFcn_,
                                          shpFcnDerivs_,
                                          atomToElementMap_);
    }
    if (nLocalGhost_>0) {
      feEngine_->evaluate_shape_functions(ghostAtomCoords_,
                                          shpFcnGhost_,
                                          shpFcnDerivsGhost_,
                                          ghostAtomToElementMap_);
    }
    
    // shpFcnMasked are the shape functions over atoms for which we need to do approximate
    // atomic quadrature to remove the FE contributions, i.e. all elements with 
    // internal atoms if atomQuadForInternal is true and all elements containing 
    // but not fully filled by atoms (e.g. has internal and ghost types) if 
    // atomQuadForInternal is false (this has *nothing* to do with thermostats!)
    // determine full/partial/empty elements
    if (nLocal_>0) {
      atomMask_.reset(nLocal_);
      if (atomQuadForInternal_) {
        atomMask_ = true;
        nLocalMask_ = nLocal_;
      }
      else {
        Array<bool> shapeFunctionMask(nNodes_);
        shapeFunctionMask = false;
        atomMask_ = false;
        nLocalMask_ = 0;
        for (int i = 0; i < nLocal_; ++i) {
          DENS_VEC coords(3);
          coords(0) = atomicCoords_(0,i);
          coords(1) = atomicCoords_(1,i);
          coords(2) = atomicCoords_(2,i);
          int eltID = feEngine_->get_feMesh()->map_to_element(coords);
          atomMask_(i) = elementMask_(eltID);
          if (atomMask_(i)) nLocalMask_++;
        }  
      }
      
      if( nLocalMask_ > 0 ) {
        atomicCoordsMask_.reset(nsd_,nLocalMask_);
        int nmask = 0;
        for (int i = 0; i < nLocal_; ++i) {
          if (atomMask_(i)) {
            for (int j = 0; j < nsd_; ++j) {
              atomicCoordsMask_(j,nmask) = atomicCoords_(j,i);
            }
            nmask++;
          }
        }
        
        // NOTE inefficient but easy to implement
        Array<int> tempArray;
        feEngine_->evaluate_shape_functions(atomicCoordsMask_,
                                            shpFcnMasked_,
                                            shpFcnDerivsMask_,
                                            tempArray);
        
      }
    }

    shpWeight_.reset(nNodes_,nNodes_);
    DENS_VEC tmpWeights(nNodes_);
    DENS_VEC weights(nNodes_);
    // Computed shape function weights
    // If weight is function of shape function index only, must be (sum_alpha
    // N_I^alpha)^-1
    if (nLocal_>0) tmpWeights = shpFcn_.col_sum();
    lammpsInterface_->allsum(tmpWeights.get_ptr(),weights.get_ptr(),nNodes_);

    for (int i = 0; i < nNodes_; i++) {
      if (weights(i) > 0.) {
        shpWeight_(i,i) = 1./weights(i);
      }
    }

    return;
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Transfer::reset_NhatOverlap()
  {
    // determine which atoms are being thermostatted, if we use a global lambda
    // it is all of them or if we used a local atom it is only those in the support
    // of shape functions on the boundary, i.e. those containing atoms in their support
    // but not having their support completely in the MD region (e.g. they contain both
    // internal and ghost atoms).  (this has *nothing* to do with FE integration!)
    if (nLocal_ > 0) {
      // form mask
      Array<bool> atomMask(nLocal_);
      if (!useLocalizedLambda_) {
        atomMask = true;
        nLocalLambda_ = nLocal_;
      }
      else {
        atomMask = false;
        nLocalLambda_ = 0;
        DENS_VEC coords(nsd_);
        for (int i = 0; i < nLocal_; ++i) {
          for (int j = 0; j < nsd_; ++j) {
            coords(j) = atomicCoords_(j,i);
          }
          int eltID = feEngine_->get_feMesh()->map_to_element(coords);
          Array<int> nodes;
          feEngine_->get_feMesh()->element_connectivity_unique(eltID,nodes);
          for (int j = 0; j < nodes.size(); ++j) {
            if (maskMat_(nodes(j),nodes(j)) > 0.) {
              atomMask(i) = true;
              nLocalLambda_++;
              break;
            }
          }
        }  
      }
      
      // compute shape functions and normalize
      NhatOverlap_.reset(nLocalLambda_,nNodeOverlap_); 
      if( nLocalLambda_ > 0 ) {
        DENS_MAT atomicCoords(nsd_,nLocalLambda_);
        int count = 0;
        for (int i = 0; i < nLocal_; ++i) {
          if (atomMask(i)) {
            for (int j = 0; j < nsd_; ++j) {
              atomicCoords(j,count) = atomicCoords_(j,i);
            }
            count++;
          }
        }
        // from shpFcnLambda get nz and map sparsity pattern
        SPAR_MAT shpFcnLambda;
        Array<int> tempArray;
        feEngine_->evaluate_shape_functions(atomicCoords,shpFcnLambda,tempArray);

        // NhatOverlap__ = shpFcnLambda*maskMat_
        DIAG_MAT scaling(maskMat_*shpWeight_);
        CLON_VEC scaleVector(scaling);
        NhatOverlapWeights_.reset(NhatOverlap_.nCols(),NhatOverlap_.nCols());
        for (int i = 0; i < shpFcnLambda.size(); ++i) {
          TRIPLET<double> myTriplet = shpFcnLambda.get_triplet(i);
          int myCol = nodeToOverlapMap_(myTriplet.j);
          if (myCol > -1) {
            NhatOverlap_.set(myTriplet.i,myCol,myTriplet.v);
            NhatOverlapWeights_(myCol) = scaleVector(myTriplet.j);
          }
        }
        NhatOverlap_.compress();

        // set up matrix that grabs only atoms to be used in the thermostat
        Trestrict_.reset(nLocalLambda_,nLocal_);
        int j = 0;
        for (int i = 0; i < nLocal_; ++i) {
          if (atomMask(i))  Trestrict_.set(j++,i,1.0);
        }
        Trestrict_.compress();
      }
    }
    else {
      nLocalLambda_ = 0;
    }

    // compute overlap ghost shape function
    // NOTE if all ghosts do not lie in the support of boundary nodes
    // this process will not be correct
    if (nLocalGhost_>0) {
      shpFcnGhostOverlap_ = shpFcnGhost_*maskMat_;
      
      DENS_VEC activeNodes(nNodes_);
      for (int i = 0; i < nNodes_; ++i) {
        activeNodes(i) = maskMat_(i,i);
      }
      DENS_VEC weights = shpFcnGhost_*activeNodes;
      DIAG_MAT weightMatrix(nLocalGhost_,nLocalGhost_);
      for (int i = 0; i < nLocalGhost_; ++i) {
        weightMatrix(i,i) = 1./weights(i);
      }
      shpFcnGhostOverlap_ = weightMatrix*shpFcnGhostOverlap_;
    }
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Transfer::reset_atomicWeightsGroup()
  {
    if (nLocal_>0) {
      atomicWeights_.reset(nLocal_,nLocal_);
      int *mask = lammpsInterface_->atom_mask();
      map<int, double>::const_iterator igroup;
      for (igroup = Valpha_.begin(); igroup != Valpha_.end(); igroup++) {
        int gid = igroup->first;
        double weight = igroup->second;
        for (int i = 0; i < nLocal_; ++i) {
          if (mask[internalToAtom_(i)] & gid) {
            atomicWeights_(i,i) = weight;
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Transfer::reset_atomicWeightsRegion()
  {
    if (nLocal_>0) {
      atomicWeights_.reset(nLocal_,nLocal_);
      map<int, double>::const_iterator iregion;
      for (iregion = Valpha_.begin(); iregion != Valpha_.end(); iregion++) {
        int rid = iregion->first;
        double weight = iregion->second;
        for (int i = 0; i < nLocal_; ++i) {
          if (lammpsInterface_->region_match(rid,
                                             atomicCoords_(0,i),
                                             atomicCoords_(1,i),
                                             atomicCoords_(2,i)))
            atomicWeights_(i,i) = weight;
        }
      }
    }
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Transfer::reset_atomicWeightsLattice()
  {
    // setup quadrature weights uniformly
    if (nLocal_>0) {
      atomicWeights_.reset(nLocal_,nLocal_);
      double volume_per_atom = lammpsInterface_->volume_per_atom();
      atomicWeights_ = volume_per_atom;
    }
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Transfer::reset_atomicWeightsElement()
  { 
    const FE_Mesh * feMesh = feEngine_->get_feMesh();
    int nElts = feMesh->get_nElements();
    DenseVector<int> eltAtomCtLocal(nElts);
    DenseVector<int> eltAtomCt(nElts);
    DenseVector<int> atomEltID;
    DENS_VEC eltAtomVol(nElts);

    for (int i = 0; i < nElts; ++i)
      eltAtomCtLocal(i) = 0;

    if (nLocal_>0) {
      atomicWeights_.reset(nLocal_,nLocal_);
      atomEltID.reset(nLocal_);

      // determine number of atoms in each element, partial sum
      const FE_Mesh * feMesh = feEngine_->get_feMesh();
      int nElts = feMesh->get_nElements();

      for (int i = 0; i < nLocal_; ++i) {
        DENS_VEC coords(3);
        coords(0) = xref_[internalToAtom_(i)][0];
        coords(1) = xref_[internalToAtom_(i)][1];
        coords(2) = xref_[internalToAtom_(i)][2];
        int thisElt = feMesh->map_to_element(coords);
        eltAtomCtLocal(thisElt) += 1;
        atomEltID(i) = thisElt;
      }
    }
  
    // mpi to determine total atoms
    lammpsInterface_->int_allsum(eltAtomCtLocal.get_ptr(),eltAtomCt.get_ptr(),nElts);

    // divide element volume by total atoms to get atomic volume
    if (nLocal_>0) {
      for (int i = 0; i < nElts; ++i) {
      
        double minx, maxx, miny, maxy, minz, maxz;
        DENS_MAT nodalCoords;
        feMesh->element_coordinates(i,nodalCoords);
        minx = nodalCoords(0,0); maxx = nodalCoords(0,0);
        miny = nodalCoords(1,0); maxy = nodalCoords(1,0);
        minz = nodalCoords(2,0); maxz = nodalCoords(2,0);
        for (int j = 1; j < nodalCoords.nCols(); ++j) {
          if (nodalCoords(0,j)<minx) minx = nodalCoords(0,j);
          if (nodalCoords(0,j)>maxx) maxx = nodalCoords(0,j);
          if (nodalCoords(1,j)<miny) miny = nodalCoords(1,j);
          if (nodalCoords(1,j)>maxy) maxy = nodalCoords(1,j);
          if (nodalCoords(2,j)<minz) minz = nodalCoords(2,j);
          if (nodalCoords(2,j)>maxz) maxz = nodalCoords(2,j);
        }
        double eltVol = (maxx-minx)*(maxy-miny)*(maxz-minz);
        if (eltVol<0) eltVol *= -1.;
        eltAtomVol(i) = eltVol/eltAtomCt(i);
      }

      for (int i = 0; i < nLocal_; ++i)
        atomicWeights_(i,i) = eltAtomVol(atomEltID(i));
    }
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Transfer::reset_atomicWeightsMultiscale(const SPAR_MAT & shapeFunctionMatrix,
                                                   DIAG_MAT & atomicVolumeMatrix)
  {
    // solve equation \sum_a N_Ia \sum_J N_Ja dV_J = \int_Omega N_I dV
    // form left-hand side
    SPAR_MAT lhs;
    compute_consistent_md_mass_matrix(shapeFunctionMatrix,lhs);
    // form right-hand side
    DIAG_MAT rhsMatrix;
    feEngine_->compute_lumped_mass_matrix(rhsMatrix);
    DENS_VEC rhs;
    rhs.copy(rhsMatrix.get_ptr(),rhsMatrix.size(),1);

    // change entries for all entries if no atoms in shape function support
    double totalVolume = rhs.sum();
    double averageVolume = totalVolume/nInternal_;
    DENS_VEC scale(nNodes_);
    for (int i = 0; i < nNodes_; i++) {
      if ((abs(lhs(i,i)) > 0.))
        scale(i) = 1.;
      else
        scale(i) = 0.;
    }
    lhs.row_scale(scale);
    for (int i = 0; i < nNodes_; i++) {
      if (scale(i) < 0.5) {
        lhs.set(i,i,1.);
        rhs(i) = averageVolume;
      }
    }
    lhs.compress();

    // solve equation
    DIAG_MAT preConditioner = lhs.get_diag();
    int myMaxIt = lhs.nRows(); // NOTE could also use a fixed parameter
    double myTol = 1.e-10; // NOTE could also use a fixed parameter

    DENS_VEC nodalAtomicVolume(nNodes_);
    int convergence = CG(lhs, nodalAtomicVolume, rhs, preConditioner, myMaxIt, myTol);

    if (convergence>0) // error if didn't converge
      throw ATC_Error(0,"CG solver did not converge in ATC_Transfer::reset_atomicWeightsMultiscale");

    // compute each atom's volume
    DENS_VEC atomicVolume(nLocal_);
    prolong(nodalAtomicVolume,atomicVolume);
    atomicVolumeMatrix.reset(atomicVolume);
  }

  //-----------------------------------------------------------------
  //
  //-----------------------------------------------------------------
  void ATC_Transfer::compute_consistent_md_mass_matrix(const SPAR_MAT & shapeFunctionMatrix,
                                                       SPAR_MAT & mdMassMatrix)
  {
    // NOTE could use sparse matrix and do communication broadcasting to exchange triplets
    int nRows = shapeFunctionMatrix.nRows();
    int nCols = shapeFunctionMatrix.nCols();
    DENS_MAT massMatrixLocal(nCols,nCols);
    DENS_MAT denseMdMassMatrix(nCols,nCols);
    if (nLocal_>0)
      massMatrixLocal = shapeFunctionMatrix.transMat(shapeFunctionMatrix);
    
    lammpsInterface_->allsum(massMatrixLocal.get_ptr(),
                             denseMdMassMatrix.get_ptr(),
                             denseMdMassMatrix.size());
    mdMassMatrix.reset(denseMdMassMatrix,1.e-10);
  }

  //=================================================================
  // Interscale operators
  //=================================================================
  //--------------------------------------------------------
  void ATC_Transfer::restrict(const MATRIX & atomData,
                              MATRIX & nodeData)
  {
    // computes nodeData = N*atomData where N are the normalized shape functions
    DENS_MAT workNodeArray(nodeData.nRows(),nodeData.nCols());
    if (nLocal_>0) {
      workNodeArray = shpWeight_*shpFcn_.transMat(atomData);
    }
    int count = nodeData.nRows()*nodeData.nCols();
    lammpsInterface_->allsum(workNodeArray.get_ptr(),nodeData.get_ptr(),count);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::restrict_unscaled(const MATRIX & atomData,
                                       MATRIX & nodeData)
  {
    // computes nodeData = N*DeltaVAtom*atomData where N are the shape functions
    DENS_MAT workNodeArray(nodeData.nRows(),nodeData.nCols());

    if (nLocal_>0) {
      workNodeArray = shpFcn_.transMat(atomicWeights_*atomData);
    }
    int count = nodeData.nRows()*nodeData.nCols();
    lammpsInterface_->allsum(workNodeArray.get_ptr(),nodeData.get_ptr(),count);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::restrict_volumetric_quantity(const MATRIX & atomData,
                                                  MATRIX & nodeData,
                                                  const SPAR_MAT & shpFcn)
  {
    // computes nodeData = N*DeltaVAtom*atomData where N are the shape functions
    DENS_MAT workNodeArray(nodeData.nRows(),nodeData.nCols());

    if (nLocal_>0) {
      workNodeArray = shpFcn.transMat(atomData);
    }
    int count = nodeData.nRows()*nodeData.nCols();
    lammpsInterface_->allsum(workNodeArray.get_ptr(),nodeData.get_ptr(),count);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::restrict_volumetric_quantity(const MATRIX & atomData,
                                                  MATRIX & nodeData)
  {
    restrict_volumetric_quantity(atomData,nodeData,shpFcn_);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::project(const MATRIX & atomData,
                             MATRIX & nodeData,
                             FieldName thisField)
  {
    // computes M*nodeData = N*DeltaVAtom*atomData where N are the shape functions
    restrict_unscaled(atomData,nodeData);

    apply_inverse_mass_matrix(nodeData,thisField);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::project_volumetric_quantity(const MATRIX & atomData,
                                                 MATRIX & nodeData,
                                                 FieldName thisField)
  {
    project_volumetric_quantity(atomData,nodeData,shpFcn_,thisField);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::project_volumetric_quantity(const MATRIX & atomData,
                                                 MATRIX & nodeData,
                                                 const SPAR_MAT & shpFcn,
                                                 FieldName thisField)
  {
    // computes nodeData = N*atomData where N are the shape functions
    restrict_volumetric_quantity(atomData,nodeData,shpFcn);

    apply_inverse_mass_matrix(nodeData,thisField);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::project_md(const MATRIX & atomData,
                                MATRIX & nodeData,
                                FieldName thisField)
  {
    // computes M*nodeData = N*DeltaVAtom*atomData where N are the shape functions
    restrict_unscaled(atomData,nodeData);

    apply_inverse_md_mass_matrix(nodeData,thisField);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::project_md_volumetric_quantity(const MATRIX & atomData,
                                                    MATRIX & nodeData,
                                                    FieldName thisField)
  {
    project_md_volumetric_quantity(atomData,nodeData,shpFcn_,thisField);
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::project_md_volumetric_quantity(const MATRIX & atomData,
                                                    MATRIX & nodeData,
                                                    const SPAR_MAT & shpFcn,
                                                    FieldName thisField)
  {
    // computes nodeData = N*atomData where N are the shape functions
    restrict_volumetric_quantity(atomData,nodeData,shpFcn);
    
    apply_inverse_md_mass_matrix(nodeData,thisField);
    
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::prolong(const MATRIX & nodeData,
                             MATRIX & atomData)
  {
    // computes the finite element interpolation at atoms atomData = N*nodeData
    if (nLocal_>0) {
      atomData = shpFcn_*nodeData;
    }
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::prolong_scaled(const MATRIX & nodeData,
                                    MATRIX & atomData)
  {
    // computes the finite element interpolation at atoms atomData = N*nodeData
    if (nLocal_>0) {
      atomData = shpFcn_*(shpWeight_*nodeData);
    }
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::prolong_ghost(const MATRIX & nodeData,
                                   MATRIX & atomData)
  {
    // computes the finite element interpolation at atoms atomData = N*nodeData
    if (nLocalGhost_>0) {
      atomData = shpFcnGhost_*nodeData;
    }
    return;
  }

  //=================================================================
  //  Interscale physics constructors
  //=================================================================
  void ATC_Transfer::compute_atomic_mass(MATRIX & atomicMasses)
  {
    atomicMasses.reset(nLocal_, 1);
    if (nLocal_>0) {
      int * type = LammpsInterface::instance()->atom_type();
      double * mass = LammpsInterface::instance()->atom_mass();
      double * rmass = LammpsInterface::instance()->atom_rmass();
      int atomIndex;
      double atomMass;
      for (int i = 0; i < nLocal_; i++) {
        atomIndex = internalToAtom_(i);
        if (mass) {
          atomMass = mass[type[atomIndex]];
        }
        else {
          atomMass = rmass[atomIndex];
        }
        atomicMasses[i] = atomMass;
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_charge(MATRIX & atomicCharges)
  {
    atomicCharges.reset(nLocal_, 1);
    if (nLocal_>0) {
      double * charge = LammpsInterface::instance()->atom_charge();
      //double coef = LammpsInterface::instance()->elemental_charge();
      int atomIndex;
      for (int i = 0; i < nLocal_; i++) {
        atomIndex = internalToAtom_(i);
        atomicCharges[i] = charge[atomIndex];
        //atomicCharges[i] = coef*charge[atomIndex];
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_temperature(MATRIX & T, 
                                                const double * const* v,
                                                double ** v2)
  {
    // v2 is for fractional step
    if (v2==NULL) v2 = const_cast<double**>(v);
    // T_alpha = 2/(3kB) * 1/2*m_alpha*(v_alpha \dot v2_alpha)
    T.reset(nLocal_, 1);
    if (nLocal_>0) {
      double ma;
      int * type = lammpsInterface_->atom_type();
      double * mass = lammpsInterface_->atom_mass();
      double * rmass = lammpsInterface_->atom_rmass();
      double kBoltzmann = lammpsInterface_->kBoltzmann();
      int atomIdx;
      double coef = 1./(nsd_*kBoltzmann);
      for (int i = 0; i < nLocal_; i++) {
        atomIdx = internalToAtom_(i);
        if (mass) {
          ma = mass[type[atomIdx]];
        }
        else {
          ma = rmass[atomIdx];
        }
        for (int j = 0; j < nsd_; j++) {
          T[i] += coef*ma*(v[atomIdx][j]*v2[atomIdx][j]);
        }
      }
    }
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_kinetic_energy(MATRIX & kineticEnergy, 
                                                   const double * const* v,
                                                   double ** v2)
  {
    // v2 is for fractional step
    if (v2==NULL) v2 = const_cast<double**>(v);
    // KE_alpha = 1/2*m_alpha*(v_alpha \dot v2_alpha)
    if (nLocal_>0) {
      kineticEnergy.reset(nLocal_, 1);
      double ma;
      int * type = lammpsInterface_->atom_type();
      double * mass = lammpsInterface_->atom_mass();
      double * rmass = lammpsInterface_->atom_rmass();
      int atomIdx;
      for (int i = 0; i < nLocal_; i++) {
        atomIdx = internalToAtom_(i);
        if (mass) {
          ma = mass[type[atomIdx]];
        }
        else {
          ma = rmass[atomIdx];
        }
        for (int j = 0; j < nsd_; j++) {
          kineticEnergy[i] += ma*(v[atomIdx][j]*v2[atomIdx][j]);
        }
      }
      kineticEnergy *= 0.5;
    }
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_temperature_roc(MATRIX & atomicPower,
                                                    const double * const* v,
                                                    const double * const* f)
  {
    // atom_power = 2/(3kB) * v_alpha * f
    atomicPower.reset(nLocal_,1);
    if (nLocal_>0) {
      double coef = 2./(nsd_*lammpsInterface_->boltz());
      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; j++) {
          atomicPower[i] += coef*(v[atomIdx][j]*f[atomIdx][j]);
        }
      }
    }
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_power(MATRIX & atomicPower,
                                          const double * const* v,
                                          const double * const* f)
  {
    // atom_power = v_alpha * f
    atomicPower.reset(nLocal_,1);
    if (nLocal_>0) {
      double coef = lammpsInterface_->ftm2v();
      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; j++) {
          atomicPower[i] += coef*(v[atomIdx][j]*f[atomIdx][j]);
        }
      }
    }
    return;
  }
  
  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_temperature_roc(MATRIX & atomicPower,
                                                    const double * const* v,
                                                    const MATRIX & f)
  {
    // atom_power = 2/(3kB) * v_alpha * f
    atomicPower.reset(nLocal_,1);
    if (nLocal_>0) {
      double coef = 2./(nsd_*(lammpsInterface_->kBoltzmann()));
      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; j++) {
          atomicPower[i] += coef*(v[atomIdx][j]*f(i,j));
        }
      }
    }
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_power(MATRIX & atomicPower,
                                          const double * const* v,
                                          const MATRIX & f)
  {
    // atom_power = v_alpha * f
    atomicPower.reset(nLocal_,1);
    if (nLocal_>0) {
      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; j++) {
          atomicPower[i] += v[atomIdx][j]*f(i,j);
        }
      }
    }
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_force_strength(MATRIX & atomicForceStrength,
                                                   const double * const* f)
  {
    // atom_force_strength = 2/(3kB) * f * f / ma
    atomicForceStrength.reset(nLocal_,1);
    if (nLocal_>0) {
      double mvv2e = lammpsInterface_->mvv2e();
      double coef = 2./(nsd_*(lammpsInterface_->kBoltzmann())) / mvv2e / mvv2e;
      int * type = lammpsInterface_->atom_type();
      double * mass = lammpsInterface_->atom_mass();
      double * rmass = lammpsInterface_->atom_rmass();
      double ma;

      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        if (mass) {
          ma = mass[type[atomIdx]];
        }
        else {
          ma = rmass[atomIdx];
        }
        for (int j = 0; j < nsd_; j++) {
          atomicForceStrength[i] += coef*(f[atomIdx][j]*f[atomIdx][j])/ma;
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_force_strength(MATRIX & atomicForceStrength,
                                                   const MATRIX & f)
  {
    // atom_force_strength = 1/(3kB) * f * f / ma
    atomicForceStrength.reset(nLocal_,1);
    if (nLocal_>0) {
      double coef = 1./(nsd_*(lammpsInterface_->kBoltzmann()));
      int * type = lammpsInterface_->atom_type();
      double * mass = lammpsInterface_->atom_mass();
      double * rmass = lammpsInterface_->atom_rmass();
      double ma;

      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        if (mass) {
          ma = mass[type[atomIdx]];
        }
        else {
          ma = rmass[atomIdx];
        }
        for (int j = 0; j < nsd_; j++) {
          atomicForceStrength[i] += coef*(f(i,j)*f(i,j))/ma;
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_force_dot(MATRIX & atomicForceDot,
                                              const double * const* f1,
                                              const MATRIX & f2)
  {
    // atom_force_strength = 2/(3kB) * f * f / ma
    atomicForceDot.reset(nLocal_,1);
    if (nLocal_>0) {
      double mvv2e = lammpsInterface_->mvv2e();
      double coef = 1./(nsd_*(lammpsInterface_->kBoltzmann())) / mvv2e;
      int * type = lammpsInterface_->atom_type();
      double * mass = lammpsInterface_->atom_mass();
      double * rmass = lammpsInterface_->atom_rmass();
      double ma;

      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        if (mass) {
          ma = mass[type[atomIdx]];
        }
        else {
          ma = rmass[atomIdx];
        }
        for (int j = 0; j < nsd_; j++) {
          atomicForceDot[i] += coef*(f1[atomIdx][j]*f2(i,j))/ma;
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_lambda_power(MATRIX & lambdaPower,
                                                 const MATRIX & lambda,
                                                 const double * const* v)
  {
    // atom_lambda_power = -(1/2) m_a * v_a^2 * lambda
    lambdaPower.reset(nLocal_,1);
    if (nLocal_>0) {
      DENS_VEC lambdaAtom(nLocal_);
      prolong_scaled(lambda,lambdaAtom);
    
      double ma;
      int * type = lammpsInterface_->atom_type();
      double * mass = lammpsInterface_->atom_mass();
      double * rmass = lammpsInterface_->atom_rmass();
      int atomIdx;
      double coef = 0.5;
      for (int i = 0; i < nLocal_; i++) {
        atomIdx = internalToAtom_(i);
        if (mass) {
          ma = mass[type[atomIdx]];
        }
        else {
          ma = rmass[atomIdx];
        }
        for (int j = 0; j < nsd_; j++) {
          lambdaPower[i] -= coef*ma*v[atomIdx][j]*v[atomIdx][j]*lambdaAtom(i);
        }
      }
    }
    return;
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_temperature_roc(MATRIX & atomicPower,
                                                    const MATRIX & force,
                                                    const double * const* v)
  {
    // atom_lambda_power = v_{\alpha} \dot force
    // NOTE account for missing mass 
    atomicPower.reset(nLocal_,1);
    if (nLocal_>0) {
      double ma;
      int * type = lammpsInterface_->atom_type();
      double * mass = lammpsInterface_->atom_mass();
      double * rmass = lammpsInterface_->atom_rmass();
      int atomIdx;
      double coef = 2./(nsd_*(lammpsInterface_->kBoltzmann()));
      for (int i = 0; i < nLocal_; i++) {
        atomIdx = internalToAtom_(i);
        if (mass) {
          ma = mass[type[atomIdx]];
        }
        else {
          ma = rmass[atomIdx];
        }
        for (int j = 0; j < nsd_; j++) {
          atomicPower[i] += coef*ma*v[atomIdx][j]*force(i,j);
        }
      }
    }
    return;
   }

  //--------------------------------------------------------
  // atomic kinematic quantities
  void ATC_Transfer::compute_atomic_displacement(DENS_MAT & atomicDisplacement,
                                                 const double * const* x)
  {
    atomicDisplacement.reset(nLocal_, nsd_);
    if (nLocal_>0) {
      for (int i = 0; i < nLocal_; ++i) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; ++j) {
          atomicDisplacement(i,j) = x[atomIdx][j] - xref_[atomIdx][j];
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_centerOfMass_displacement(DENS_MAT & atomicComD, 
                                                              const double * const* x)
  {
    atomicComD.reset(nLocal_, nsd_);
    if (nLocal_>0) {
      double *mass = lammpsInterface_->atom_mass();
      if (mass) {
        int *type = lammpsInterface_->atom_type();
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int j = 0; j < nsd_; j++)
            atomicComD(i,j) = mass[type[atomIdx]]*(x[atomIdx][j]- xref_[atomIdx][j]);
        }
      }
      else {
        double *rmass = lammpsInterface_->atom_rmass();
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int j = 0; j < nsd_; j++)
            atomicComD(i,j) = rmass[atomIdx]*(x[atomIdx][j]- xref_[atomIdx][j]);
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_position(DENS_MAT & atomicDisplacement,
                                             const double * const* x)
  {
    atomicDisplacement.reset(nLocal_, nsd_);
    if (nLocal_>0) {
      for (int i = 0; i < nLocal_; ++i) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; ++j)
          atomicDisplacement(i,j) = x[atomIdx][j];
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_centerOfMass(DENS_MAT & atomicCom, 
                                                 const double * const* x)
  {
    atomicCom.reset(nLocal_, nsd_);
    if (nLocal_>0) {
      double *mass = lammpsInterface_->atom_mass();
      if (mass) {
        int *type = lammpsInterface_->atom_type();
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int j = 0; j < nsd_; j++)
            atomicCom(i,j) = mass[type[atomIdx]]*x[atomIdx][j];
        }
      }
      else {
        double *rmass = lammpsInterface_->atom_rmass();
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int j = 0; j < nsd_; j++)
            atomicCom(i,j) = rmass[atomIdx]*x[atomIdx][j];
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_velocity(DENS_MAT & atomicVelocity,
                                             const double * const* v)
  {
    atomicVelocity.reset(nLocal_, nsd_);
    if (nLocal_>0) {
      for (int i = 0; i < nLocal_; ++i) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; ++j)
          atomicVelocity(i,j) = v[atomIdx][j];
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_momentum(DENS_MAT & atomicMomentum, 
                                             const double * const* v)
  {
    atomicMomentum.reset(nLocal_, nsd_);
    if (nLocal_>0) {
      double *mass = lammpsInterface_->atom_mass();
      if (mass) {
        int *type = lammpsInterface_->atom_type();
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int j = 0; j < nsd_; j++)
            atomicMomentum(i,j) = mass[type[atomIdx]]*v[atomIdx][j];
        }
      }
      else {
        double *rmass = lammpsInterface_->atom_rmass();
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int j = 0; j < nsd_; j++)
            atomicMomentum(i,j) = rmass[atomIdx]*v[atomIdx][j];
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_acceleration(DENS_MAT & atomicAcceleration, 
                                                 const double * const* f)
  {
    atomicAcceleration.reset(nLocal_, nsd_);
    if (nLocal_>0) {
      double coef = lammpsInterface_->ftm2v();
      double *mass = lammpsInterface_->atom_mass();
      if (mass) {
        int *type = lammpsInterface_->atom_type();
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int j = 0; j < nsd_; j++)
            atomicAcceleration(i,j) = coef*f[atomIdx][j]/mass[type[atomIdx]];
        }
      }
      else {
        double *rmass = lammpsInterface_->atom_rmass();
        for (int i = 0; i < nLocal_; i++) {
          int atomIdx = internalToAtom_(i);
          for (int j = 0; j < nsd_; j++)
            atomicAcceleration(i,j) = coef*f[atomIdx][j]/rmass[atomIdx];
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::compute_atomic_force(DENS_MAT & atomicForce, 
                                          const double * const* f)
  {
    atomicForce.reset(nLocal_, nsd_);
    if (nLocal_>0) {
      double coef = lammpsInterface_->ftm2v();
      for (int i = 0; i < nLocal_; i++) {
        int atomIdx = internalToAtom_(i);
        for (int j = 0; j < nsd_; j++)
          atomicForce(i,j) = coef*f[atomIdx][j];
      }
    }
  }

  //=================================================================
  // FE mapping operations
  //=================================================================
  // Mappings between overlap nodes and unique nodes
  void ATC_Transfer::map_unique_to_overlap(const MATRIX & uniqueData,
                                           MATRIX & overlapData)
  {
    for (int i = 0; i < nNodes_; i++) {
      if (nodeToOverlapMap_(i)!=-1) {
        for (int j = 0; j < uniqueData.nCols(); j++) {
          overlapData(nodeToOverlapMap_(i),j) = uniqueData(i,j);
        }
      }
    }
  }

  void ATC_Transfer::map_overlap_to_unique(const MATRIX & overlapData,
                                           MATRIX & uniqueData)
  {
    uniqueData.reset(nNodes_,overlapData.nCols());
    for (int i = 0; i < nNodeOverlap_; i++) {
      for (int j = 0; j < overlapData.nCols(); j++) {
        uniqueData(overlapToNodeMap_(i),j) = overlapData(i,j);
      }
    }
  }

  //--------------------------------------------------------
  // 
  //--------------------------------------------------------
  void ATC_Transfer::construct_prescribed_data_manager (void) {
    prescribedDataMgr_ = new PrescribedDataManager(feEngine_,fieldSizes_);
  }

  //========================================================
  // FE functions
  //========================================================

  //--------------------------------------------------------
  // handle nodesets
  //--------------------------------------------------------
  void ATC_Transfer::set_fixed_nodes()
  {
    // set fields
    prescribedDataMgr_->set_fixed_fields(simTime_,
                                         fields_,dot_fields_,ddot_fields_,dddot_fields_);


    // set related data
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          if (prescribedDataMgr_->is_fixed(inode,thisField,thisIndex)) {
            dot_fieldsMD_    [thisField](inode,thisIndex) = 0.;
            ddot_fieldsMD_   [thisField](inode,thisIndex) = 0.;
            dot_dot_fieldsMD_[thisField](inode,thisIndex) = 0.;
            rhs_[thisField](inode,thisIndex) = 0.;
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::set_initial_conditions()
  {
    // set fields 
    prescribedDataMgr_->set_initial_conditions(simTime_,
                                               fields_,dot_fields_,ddot_fields_,dddot_fields_);

    // set (all) related data
    map<FieldName,int>::const_iterator field;
    for (field = fieldSizes_.begin(); field!=fieldSizes_.end(); field++) {
      FieldName thisField = field->first;
      int thisSize = field->second;
      for (int inode = 0; inode < nNodes_ ; ++inode) {
        for (int thisIndex = 0; thisIndex < thisSize ; ++thisIndex) {
          dot_fieldsMD_    [thisField](inode,thisIndex) = 0.;
          ddot_fieldsMD_   [thisField](inode,thisIndex) = 0.;
          dot_dot_fieldsMD_[thisField](inode,thisIndex) = 0.;
          rhs_[thisField](inode,thisIndex) = 0.;
        }
      }
    }
  }

  //--------------------------------------------------------
  void ATC_Transfer::set_sources()
  {
    // set fields
    prescribedDataMgr_->set_sources(simTime_,sources_);

  }

  //--------------------------------------------------------
  void ATC_Transfer::output()
  {
    if (lammpsInterface_->comm_rank() == 0) {
      map< pair <string, FieldName>, double >::iterator iter;
      for (iter = nsetCompute_.begin(); iter != nsetCompute_.end();iter++){
        pair <string, FieldName> id = iter->first;
        string nsetName = id.first;
        FieldName field = id.second;
        double sum = 0.0;
        const set<int> nset = feEngine_->get_feMesh()->get_nodeset(nsetName);
        set< int >::const_iterator itr;
        for (itr = nset.begin(); itr != nset.end();itr++){
          int node = *itr;
          sum += fields_[field](node,0);
        }
        iter->second = sum;
        string name = nsetName + "_" + field_to_string(field);
        feEngine_->add_global(name, sum);
      }
    }
  }


  //--------------------------------------------------------
  // use part, return to OutputManager maps of data & geometry
  // base class atomic output or just data creation
  void ATC_Transfer::atomic_output()
  {
    int me = lammpsInterface_->comm_rank();
    int nTotal = (int) lammpsInterface_->natoms();
    double ** x = lammpsInterface_->xatom();
    double ** v = lammpsInterface_->vatom();
    double ** f = lammpsInterface_->fatom();
    int * type = lammpsInterface_->atom_type();
    double * mass = lammpsInterface_->atom_mass();
    double * rmass = lammpsInterface_->atom_rmass();
    // map for output data
    OUTPUT_LIST atom_data;

    // geometry
    // NOTE this does not properly handle the output of subset of the entire system
    DENS_MAT atomicCoordinatesGlobal(nsd_,nTotal);
    DENS_MAT atomicCoordinatesLocal(nsd_,nTotal);
    DENS_MAT atomicTypeLocal(nTotal,1);
    DENS_MAT atomicType(nTotal,1);
    int * tag = lammpsInterface_->atom_tag();
    if (nLocal_>0) {
      for (int i = 0; i < nLocal_; ++i) {
        //cout << "procID : " << me << ", iLocal_ : " << i << flush;
        int atom_index   = internalToAtom_(i);
        //cout << ", internalToAtom_ : " << atom_index << flush;
        int global_index = tag[atom_index]-1;
        //cout << ", global_index : " << global_index << flush;
        for (int j = 0; j < nsd_; j++) {
          atomicCoordinatesLocal(j,global_index) = x[atom_index][j];
          //cout << ", x : " << x[atom_index][j] << flush;
        }
        atomicTypeLocal(global_index,0) = type[atom_index];
        //cout << "\n";
      }
    }
    if (nLocalGhost_>0) {
      for (int i = 0; i < nLocalGhost_; ++i) {
        int atom_index = ghostToAtom_(i);
        int global_index = tag[atom_index]-1;
        atomicTypeLocal(global_index,0) = 1.0;
        for (int j = 0; j < nsd_; j++) {
          atomicCoordinatesLocal(j,global_index) = x[atom_index][j];
        }
      }
    }
    lammpsInterface_->allsum(atomicCoordinatesLocal.get_ptr(),
                             atomicCoordinatesGlobal.get_ptr(), nTotal*nsd_);
    if (me==0) mdOutputManager_.write_geometry(atomicCoordinatesGlobal);

    lammpsInterface_->allsum(atomicTypeLocal.get_ptr(),
                             atomicType.get_ptr(), nTotal);

    DENS_MAT atomicPositionGlobal(nTotal,nsd_);
    if (atomicOutputMask_.count("position") > 0) {
      DENS_MAT atomicPositionLocal(nTotal,nsd_);
      if (nLocal_>0) {
        for (int i = 0; i < nLocal_; ++i) {
          int atom_index   = internalToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicPositionLocal(global_index,j) = x[atom_index][j];
          }
        }
      }
      if (nLocalGhost_>0) {
        for (int i = 0; i < nLocalGhost_; ++i) {
          int atom_index = ghostToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicPositionLocal(global_index,j) = x[atom_index][j];
          }
        }
      }
      lammpsInterface_->allsum(atomicPositionLocal.get_ptr(),
                               atomicPositionGlobal.get_ptr(), nTotal*nsd_);
      if (me==0) atom_data["atomicPosition"] = & atomicPositionGlobal;
    }

    // point data, output: displacement, velocity, temperature
    if (me==0) atom_data["atomicType"] = & atomicType;

    // displacement
    // NOTE make this a stand alone function
#define CHEAP_UNWRAP
#ifdef CHEAP_UNWRAP
    double box_bounds[3][3];
    lammpsInterface_->get_box_bounds(box_bounds[0][0],box_bounds[1][0],
                                     box_bounds[0][1],box_bounds[1][1],
                                     box_bounds[0][2],box_bounds[1][2]);
    double box_length[3];
    for (int k = 0; k < 3; k++) {
      box_length[k] = box_bounds[1][k] - box_bounds[0][k];
    }
#endif
    bool latticePeriodicity[3];
    latticePeriodicity[0] = (bool) lammpsInterface_->xperiodic();   
    latticePeriodicity[1] = (bool) lammpsInterface_->yperiodic();
    latticePeriodicity[2] = (bool) lammpsInterface_->zperiodic();

    DENS_MAT atomicDisplacementGlobal(nTotal,nsd_);
    
    if (atomicOutputMask_.count("displacement") > 0) {
      DENS_MAT atomicDisplacementLocal(nTotal,nsd_);
      if (nLocal_>0) {
        for (int i = 0; i < nLocal_; ++i) {
          int atom_index   = internalToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicDisplacementLocal(global_index,j) 
              = x[atom_index][j] -xref_[atom_index][j];
            if (latticePeriodicity[j]) {
#ifdef CHEAP_UNWRAP
              // the cheap version of unwrapping atomic positions
              double u = atomicDisplacementLocal(global_index,j);
              if (u >=  0.5*box_length[j]) { u -= box_length[j]; }
              if (u <= -0.5*box_length[j]) { u += box_length[j]; }
              atomicDisplacementLocal(global_index,j) = u;
#else
              double x_a[3];
              lammpsInterface_->unwrap_coordinates(atom_index,x_a);
              atomicDisplacementLocal(global_index,j) 
                = x_a[j] -xref_[atom_index][j];
#endif
            }
          }
        }
      }
      if (nLocalGhost_>0) {
        for (int i = 0; i < nLocalGhost_; ++i) {
          int atom_index = ghostToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicDisplacementLocal(global_index,j) 
              = x[atom_index][j] -xref_[atom_index][j];
            if (latticePeriodicity[j]) {
#ifdef CHEAP_UNWRAP
              // the cheap version of unwrapping atomic positions
              double u = atomicDisplacementLocal(global_index,j);
              if (u >=  0.5*box_length[j]) { u -= box_length[j]; }
              if (u <= -0.5*box_length[j]) { u += box_length[j]; }
              atomicDisplacementLocal(global_index,j) = u;
#else
              double x_a[3];
              lammpsInterface_->unwrap_coordinates(atom_index,x_a);
              atomicDisplacementLocal(global_index,j) 
                = x_a[j] -xref_[atom_index][j];
#endif
            }
          }
        }
      }
      lammpsInterface_->allsum(atomicDisplacementLocal.get_ptr(),
                               atomicDisplacementGlobal.get_ptr(), nTotal*nsd_);
      if (me==0) atom_data["atomicDisplacement"] = & atomicDisplacementGlobal;
    }

    // velocity
    DENS_MAT atomicVelocityGlobal(nTotal,nsd_);
    if (atomicOutputMask_.count("velocity") > 0) {
      DENS_MAT atomicVelocityLocal(nTotal,nsd_);
      if (nLocal_>0) {
        for (int i = 0; i < nLocal_; ++i) {
          int atom_index   = internalToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicVelocityLocal(global_index,j) = v[atom_index][j];
          }
        }
      }
      if (nLocalGhost_>0) {
        for (int i = 0; i < nLocalGhost_; ++i) {
          int atom_index = ghostToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicVelocityLocal(global_index,j) = v[atom_index][j];
          }
        }
      }
      lammpsInterface_->allsum(atomicVelocityLocal.get_ptr(),
                               atomicVelocityGlobal.get_ptr(), nTotal*nsd_);
      if (me==0) atom_data["atomicVelocity"] = & atomicVelocityGlobal;
    }

    // temperature
    DENS_MAT atomicTemperatureGlobal(nTotal,1);
    if (atomicOutputMask_.count("temperature") > 0) {
      DENS_MAT atomicTemperatureLocal(nTotal,1);
      atomicTemperatureLocal = 0.0;
      double ma;
      double coef = 1./(nsd_*(lammpsInterface_->kBoltzmann()));
      if (nLocal_>0) {
        for (int i = 0; i < nLocal_; ++i) {
          int atom_index   = internalToAtom_(i);
          int global_index = tag[atom_index]-1;
          if (mass) ma = mass[type[atom_index]];
          else      ma = rmass[atom_index];
          for (int j = 0; j < nsd_; j++) {
            double vj = v[atom_index][j];
            atomicTemperatureLocal(global_index,0) += coef*ma*vj*vj;
          }
        }
      }
      if (nLocalGhost_>0) {
        for (int i = 0; i < nLocalGhost_; ++i) {
          int atom_index = ghostToAtom_(i);
          int global_index = tag[atom_index]-1;
          if (mass) ma = mass[type[atom_index]];
          else      ma = rmass[atom_index];
          for (int j = 0; j < nsd_; j++) {
            double vj = v[atom_index][j];
            atomicTemperatureLocal(global_index,0) += coef*ma*vj*vj;
          }
        }
      }
      lammpsInterface_->allsum(atomicTemperatureLocal.get_ptr(),
                               atomicTemperatureGlobal.get_ptr(), nTotal);
      if (me==0) atom_data["atomicTemperature"] = & atomicTemperatureGlobal;
    }

    // force
    DENS_MAT atomicForceGlobal(nTotal,nsd_);
    if (atomicOutputMask_.count("force") > 0) {
      DENS_MAT atomicForceLocal(nTotal,nsd_);
      if (nLocal_>0) {
        for (int i = 0; i < nLocal_; ++i) {
          int atom_index   = internalToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicForceLocal(global_index,j) = f[atom_index][j];
          }
        }
      }
      if (nLocalGhost_>0) {
        for (int i = 0; i < nLocalGhost_; ++i) {
          int atom_index = ghostToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicForceLocal(global_index,j) = f[atom_index][j];
          }
        }
      }
      lammpsInterface_->allsum(atomicForceLocal.get_ptr(),
                               atomicForceGlobal.get_ptr(), nTotal*nsd_);
      if (me==0) atom_data["atomicForce"] = & atomicVelocityGlobal;
    }


    // atomic stress 
    double ** atomStress = NULL;
    try{ 
      atomStress = lammpsInterface_->compute_vector_data("stress/atom");
    }
    catch(ATC::ATC_Error& atcError) {
      atomStress = NULL;
    }
    // if compute stress/atom isn't create, this code is skipped
    DENS_MAT atomicVirialGlobal(nTotal,6);
    if (atomicOutputMask_.count("virial") > 0 && atomStress) {
      DENS_MAT atomicVirialLocal(nTotal,6);
      if (nLocal_>0) {
        for (int i = 0; i < nLocal_; ++i) {
          int atom_index   = internalToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int k = 0; k < 6; k++) {
            atomicVirialLocal(global_index,k) = atomStress[atom_index][k];
          }
        }
      }
      if (nLocalGhost_>0) {
        for (int i = 0; i < nLocalGhost_; ++i) {
          int atom_index = ghostToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int k = 0; k < 6; k++) {
            atomicVirialLocal(global_index,k) = atomStress[atom_index][k]; 
          }
        }
      }
      lammpsInterface_->allsum(atomicVirialLocal.get_ptr(),
                               atomicVirialGlobal.get_ptr(), nTotal*6);
      //atomicVirialGlobal *= nktv2p;
      if (me==0) atom_data["atomicVirial"] = & atomicVirialGlobal;
    }

    // potential energy
    // if compute pe/atom isn't created, this code is skipped
    double * atomPE = NULL;
    try{ 
      atomPE = lammpsInterface_->compute_scalar_data("pe/atom");
    }
    catch(ATC::ATC_Error& atcError) {
      atomPE = NULL;
    }
    DENS_MAT atomicEnergyGlobal(nTotal,1);
    if (atomicOutputMask_.count("potential_energy") > 0 && atomPE) {
      DENS_MAT atomicEnergyLocal(nTotal,1);
      if (nLocal_>0) {
        for (int i = 0; i < nLocal_; ++i) {
          int atom_index   = internalToAtom_(i);
          int global_index = tag[atom_index]-1;
          atomicEnergyLocal(global_index,0) = atomPE[atom_index];
        }
      }
      if (nLocalGhost_>0) {
        for (int i = 0; i < nLocalGhost_; ++i) {
          int atom_index = ghostToAtom_(i);
          int global_index = tag[atom_index]-1;
          atomicEnergyLocal(global_index,0) = atomPE[atom_index];
        }
      }
      lammpsInterface_->allsum(atomicEnergyLocal.get_ptr(),
                               atomicEnergyGlobal.get_ptr(), nTotal);
      if (me==0) atom_data["atomicPotentialEnergy"] = & atomicEnergyGlobal;
    }

    // centrosymmetry 
    // NOTE if hardy doesn't create compute centro/atom this is skipped
    // we should move this to base class with a request from derived
    //double * atomCentro = lammpsInterface_->atomCentro_compute();
    double * atomCentro = NULL;
    try{ 
      atomCentro = lammpsInterface_->compute_scalar_data("centro/atom");
    }
    catch(ATC::ATC_Error& atcError) {
      atomCentro = NULL;
    }
    // if compute centro/atom isn't create, this code is skipped
    DENS_MAT atomicCentroGlobal(nTotal,1);
    if (atomicOutputMask_.count("centrosymmetry") > 0 && atomCentro) {
      DENS_MAT atomicCentroLocal(nTotal,1);
      if (nLocal_>0) {
        for (int i = 0; i < nLocal_; ++i) {
          int atom_index   = internalToAtom_(i);
          int global_index = tag[atom_index]-1;
          for (int j = 0; j < nsd_; j++) {
            atomicCentroLocal(global_index,0) = atomCentro[atom_index];
          }
        }
      }
      if (nLocalGhost_>0) {
        for (int i = 0; i < nLocalGhost_; ++i) {
          int atom_index = ghostToAtom_(i);
          int global_index = tag[atom_index]-1;
          atomicCentroLocal(global_index,0) = atomCentro[atom_index];
        }
      }
      lammpsInterface_->allsum(atomicCentroLocal.get_ptr(),
                               atomicCentroGlobal.get_ptr(), nTotal);
      if (me==0) atom_data["atomicCentrosymmetry"] = & atomicCentroGlobal;
    }

    if (me==0) mdOutputManager_.write_data(simTime_, & atom_data);
 
  }

  //=================================================================
  // FE interfaces
  //=================================================================
  void ATC_Transfer::compute_boundary_flux(const Array2D<bool> & rhsMask,
                                           const FIELDS & fields, 
                                           FIELDS & rhs) 
  {
    if (bndyIntType_ == FE_QUADRATURE) {
      //cout << "FE QUAD " << lammpsInterface_->comm_rank() << endl;
      feEngine_->compute_boundary_flux(rhsMask,
                                       fields,
                                       physicsModel_,
                                       elementToMaterialMap_,
                                       (* bndyFaceSet_),
                                       rhs);
    }
    else if (bndyIntType_ == FE_INTERPOLATION) {
      //cout << "FE INTERP " << lammpsInterface_->comm_rank() << endl;
      feEngine_->compute_boundary_flux(rhsMask,
                                       fields,
                                       physicsModel_,
                                       elementToMaterialMap_,
                                       atomMaterialGroups_,
                                       atomicWeights_,
                                       shpFcn_,
                                       shpFcnDerivs_,
                                       fluxMask_,
                                       rhs);
    }
    else if (bndyIntType_ == NO_QUADRATURE) {
      //cout << "RESET BOUNDARY FLUX " << lammpsInterface_->comm_rank() << endl;
      FIELDS::const_iterator field;
      for (field = fields.begin(); field != fields.end(); field++) {
	FieldName thisFieldName = field->first;
	if (rhsMask(thisFieldName,FLUX)) {
	  int nrows = (field->second).nRows();
	  int ncols = (field->second).nCols();
	  rhs[thisFieldName].reset(nrows,ncols);
	}
      }
    }
  }

  //-----------------------------------------------------------------
  void ATC_Transfer::compute_flux(const Array2D<bool> & rhsMask,
                                  const FIELDS & fields, 
                                  GRAD_FIELDS & flux,
                                  const PhysicsModel * physicsModel) 
  {
    if (! physicsModel) { physicsModel = physicsModel_; }
    feEngine_->compute_flux(rhsMask,
                            fields,
                            physicsModel,
                            elementToMaterialMap_,
                            flux);
  }

  //-----------------------------------------------------------------
  void ATC_Transfer::evaluate_rhs_integral(const Array2D<bool> & rhsMask,
                                           const FIELDS & fields, FIELDS & rhs,
                                           const IntegrationDomainType domain,
                                           const PhysicsModel * physicsModel)
  {
    if (!physicsModel) physicsModel = physicsModel_;

    // compute B and N flux NOTE linear or not
    if      (domain == FE_DOMAIN ) {
      feEngine_->compute_rhs_vector(rhsMask, 
                                    fields, 
                                    physicsModel,
                                    elementToMaterialMap_,
                                    rhs, 
                                    &elementMask_);
      if (nLocalMask_ > 0) {
        feEngine_->compute_rhs_vector(rhsMask,
                                      fields,
                                      physicsModel,
                                      atomMaterialGroups_,
                                      atomicWeightsMask_,
                                      shpFcnMasked_,
                                      shpFcnDerivsMask_,
                                      rhsAtomDomain_);
        for (FIELDS::const_iterator field = fields.begin(); 
             field != fields.end(); field++) {
          FieldName thisFieldName = field->first;
          rhs[thisFieldName] -= rhsAtomDomain_[thisFieldName];
        }
      }
    }
    else if (domain == ATOM_DOMAIN) {
      // NOTE should we let this actually do a full atom quadrature integration?
      if (nLocalMask_ > 0) {
        feEngine_->compute_rhs_vector(rhsMask,
                                    fields,
                                    physicsModel,
                                    atomMaterialGroups_,
                                    atomicWeightsMask_,
                                    shpFcnMasked_,
                                    shpFcnDerivsMask_,
                                    rhs);
      }
    }
    else { // domain == FULL_DOMAIN
      feEngine_->compute_rhs_vector(rhsMask,
                                    fields,
                                    physicsModel,
                                    elementToMaterialMap_,
                                    rhs);
    }
  }

  //-----------------------------------------------------------------
  void ATC_Transfer::compute_rhs_vector(const Array2D<bool> & rhsMask,
                                        const FIELDS & fields, FIELDS & rhs,
                                        const IntegrationDomainType domain,
                                        const PhysicsModel * physicsModel)
  {
    if (!physicsModel) physicsModel = physicsModel_;
    // NOTE what about FE_DOMAIN & Extrinsic???
    int index;
    // compute FE contributions
    evaluate_rhs_integral(rhsMask,fields,rhs,domain,physicsModel);
    
    // NOTE change massMask to rhsMask
    Array<FieldName> massMask;
    for (FIELDS::const_iterator field = fields.begin();
         field != fields.end(); field++)
      {
        FieldName thisFieldName = field->first;
      if (rhsMask(thisFieldName,PRESCRIBED_SOURCE)) {
        if (is_intrinsic(thisFieldName)) {
          rhs[thisFieldName] += fluxMaskComplement_*sources_[thisFieldName];
        } // NOTE perhaps replace fluxMask in future
        else {
          rhs[thisFieldName] +=                     sources_[thisFieldName];
        }
      }
      
      // add in sources from extrinsic models
      if (rhsMask(thisFieldName,EXTRINSIC_SOURCE)) {
        if (is_intrinsic(thisFieldName)) {
          rhs[thisFieldName] += fluxMaskComplement_*extrinsicSources_[thisFieldName];
        }
        else {
          rhs[thisFieldName] +=                     extrinsicSources_[thisFieldName];
        }
      }

      
      }
    
  }

  //-----------------------------------------------------------------
  void ATC_Transfer::compute_atomic_sources(const Array2D<bool> & fieldMask,
                                            const FIELDS & fields,
                                            FIELDS & atomicSources)
  {

    // NOTE change massMask to rhsMask
    Array<FieldName> massMask;
    for (FIELDS::const_iterator field = fields.begin(); 
         field != fields.end(); field++) {
      FieldName thisFieldName = field->first;
      if (is_intrinsic(thisFieldName)) {
        atomicSources[thisFieldName] = 0.;
        if (fieldMask(thisFieldName,FLUX)) {
          atomicSources[thisFieldName] = boundaryFlux_[thisFieldName];
        }
        if (fieldMask(thisFieldName,PRESCRIBED_SOURCE)) {
          atomicSources[thisFieldName] -= fluxMask_*sources_[thisFieldName];
        } // NOTE perhaps replace fluxMask in future
        
        
        // add in sources from extrinsic models
        if (fieldMask(thisFieldName,EXTRINSIC_SOURCE))
          atomicSources[thisFieldName] -= fluxMask_*extrinsicSources_[thisFieldName];
        
      }
    }
  }



  void ATC_Transfer::update_filter(MATRIX & filteredQuantity,
                                   const MATRIX & unfilteredQuantity,
                                   MATRIX & unfilteredQuantityOld,
                                   const double dt)
  {
    double tau = timeFilterManager_.get_filter_scale();
    filteredQuantity = 1./(1./dt+1./(2*tau))*( 1./(2*tau)*
                                               (unfilteredQuantity+unfilteredQuantityOld) +
                                               (1./dt-1./(2*tau))*filteredQuantity);
    unfilteredQuantityOld = unfilteredQuantity;
    return;
  }
  
  //--------------------------------------------------------
  void ATC_Transfer::update_filter_implicit(MATRIX & filteredQuantity,
                                            const MATRIX & unfilteredQuantity,
                                            const double dt)
  {
    double tau = timeFilterManager_.get_filter_scale();
    filteredQuantity = (1./(1.+dt/tau))*((dt/tau)*unfilteredQuantity + filteredQuantity);
    return;
  }
};
