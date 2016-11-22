// ATC transfer headers
#include "PerAtomQuantityLibrary.h"
#include "ATC_Transfer.h"
#include "FE_Engine.h"
#include "LammpsInterface.h"
#include <typeinfo>
#include <sstream>
#include <iostream>

using std::map;
using std::ifstream;
using std::stringstream;
using std::set;
using std::string;
using std::vector;

using ATC_Utility::to_string;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomToElementMap
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomToElementMap::AtomToElementMap(ATC_Method * atc,
                                     PerAtomQuantity<double> * atomPositions,
                                     AtomType atomType) :
    ProtectedAtomQuantity<int>(atc,1,atomType),
    atomPositions_(atomPositions)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomPositions_)
      atomPositions_ = interscaleManager.per_atom_quantity("AtomicCoarseGrainingPositions");
    if (!atomPositions_)
      throw ATC_Error("AtomToElementMap::AtomToElementMap - atom position quantity is undefined");

    atomPositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomToElementMap::~AtomToElementMap()
  {
    atomPositions_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomToElementMap::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<int>::reset();
      const DENS_MAT & position(atomPositions_->quantity());
      const FE_Mesh * feMesh = atc_.fe_engine()->fe_mesh();
      int nsd = atc_.nsd();
      DENS_VEC coords(nsd);
      for (int i = 0; i < quantity_.nRows(); i++) {
        for (int j = 0; j < nsd; j++) {
          coords(j) = position(i,j);
        }
        quantity_(i,0) = feMesh->map_to_element(coords);
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomInElementSet
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomInElementSet::AtomInElementSet(ATC_Method * atc,
                     PerAtomQuantity<int> * map,
                     ESET eset, int type):
  atc_(atc,INTERNAL), 
  map_(map),eset_(eset),type_(type),
  quantityToLammps_(atc_.atc_to_lammps_map())
  {
    map_->register_dependence(this);
    needReset_ = true;
    
  }
  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomInElementSet::~AtomInElementSet()
  {
    map_->remove_dependence(this);
  }
  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomInElementSet::reset() 
  {
    if (map_->need_reset() || needReset_) {
      list_.clear();
      INT_ARRAY map = map_->quantity();
      int * type = ATC::LammpsInterface::instance()->atom_type(); 
      ESET::const_iterator itr;
      const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
      for (int i = 0; i < map.size(); i++) {
        for (itr=eset_.begin(); itr != eset_.end(); itr++) {
          if (map(i,0) == *itr) {
            int a = quantityToLammps(i);
            if (type[a] == type_) {
              list_.push_back(ID_PAIR(i,a)); 
              break;
            }
          }
        }
      }
      needReset_ = false;
    }
  }
  //--------------------------------------------------------
  //  qauntity
  //--------------------------------------------------------
  const ID_LIST &  AtomInElementSet::quantity()
  {
    reset();
    return list_;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomVolumeUser
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomVolumeUser::AtomVolumeUser(ATC_Method * atc,
                                 map<int,double> & atomGroupVolume,
                                 AtomType atomType) :
    ProtectedAtomDiagonalMatrix<double>(atc,atomType),
    atomGroupVolume_(atomGroupVolume),
    lammpsInterface_(atc->lammps_interface()),
    atcToLammps_(atc->internal_to_atom_map())
  {
    // do nothing
  }
  
  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomVolumeUser::reset() const
  {
    if (need_reset()) {
      PerAtomDiagonalMatrix<double>::reset();
      const int *mask = lammpsInterface_->atom_mask();
      quantity_ = 0.;

      map<int, double>::const_iterator igroup;
      for (igroup = atomGroupVolume_.begin(); igroup != atomGroupVolume_.end(); igroup++) {
        int gid = igroup->first;
        double weight = igroup->second;
        for (int i = 0; i < quantity_.nRows(); ++i) {
          if (mask[atcToLammps_(i)] & gid) {
            quantity_(i,i) = weight;
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomVolumeGroup
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomVolumeGroup::AtomVolumeGroup(ATC_Method * atc,
                                   map<int,double> & atomGroupVolume,
                                   AtomType atomType) :
    AtomVolumeUser(atc,atomGroupVolume,atomType),
    atomGroupVolume_(atomGroupVolume),
    lammpsInterface_(atc->lammps_interface()),
    atcToLammps_(atc->internal_to_atom_map())
  {
    // Uses group bounding box as total volume
    // ASSUME will *only* work if atoms exactly fill the box
    int ngroup = lammpsInterface_->ngroup();
    const int *mask = lammpsInterface_->atom_mask();
    double * bounds;
      
    bounds = new double[6];
    for (int i = 0; i < ngroup; ++i) {
      lammpsInterface_->group_bounds(i, bounds);
      atomGroupVolume_[lammpsInterface_->group_bit(i)] = 
        (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]);
    }
    delete [] bounds;
      
    INT_VECTOR localCount(ngroup);
    INT_VECTOR globalCount(ngroup);
      
    // loop over atoms
    localCount = 0;
    for (int i = 0; i < atcToLammps_.size(); ++i) {
      for (int j = 0; j < ngroup; ++j) {
        if (mask[atcToLammps_(i)] & lammpsInterface_->group_bit(j))
          localCount(j)++;
      }
    }
      
    // communication to get total atom counts per group
    lammpsInterface_->int_allsum(localCount.ptr(),
                                 globalCount.ptr(),ngroup);
    
    for (int i = 0; i < ngroup; ++i) {
      int iGroupBit = lammpsInterface_->group_bit(i); 
      if (globalCount(i) > 0)
        atomGroupVolume_[iGroupBit] /= globalCount(i);
      else
        atomGroupVolume_[iGroupBit] = 0;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomVolumeLattice
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomVolumeLattice::AtomVolumeLattice(ATC_Method * atc,
                                       AtomType atomType) :
    ProtectedAtomDiagonalMatrix<double>(atc,atomType),
    lammpsInterface_(atc->lammps_interface())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomVolumeLattice::reset() const
  {
    if (need_reset()) {
      PerAtomDiagonalMatrix<double>::reset();
      quantity_ = lammpsInterface_->volume_per_atom();
    }
  }
  
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomVolumeElement
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomVolumeElement::AtomVolumeElement(ATC_Method * atc,
                                       PerAtomQuantity<int> * atomElement,
                                       AtomType atomType) :
    ProtectedAtomDiagonalMatrix<double>(atc,atomType),
    atomElement_(atomElement),
    lammpsInterface_(atc->lammps_interface()),
    feMesh_((atc->fe_engine())->fe_mesh())
  {
    if (!atomElement_) {
      InterscaleManager & interscaleManager(atc->interscale_manager());
      atomElement_ = interscaleManager.per_atom_int_quantity("AtomElement");
    }
    atomElement_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtomVolumeElement::~AtomVolumeElement()
  {
    atomElement_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomVolumeElement::reset() const
  {
    if (need_reset()) {
      PerAtomDiagonalMatrix<double>::reset();
      int nElts = feMesh_->num_elements();
      int nLocal = quantity_.nRows();
      _elementAtomCountLocal_.reset(nElts);
      _elementAtomCount_.resize(nElts);
      const INT_ARRAY & atomElement(atomElement_->quantity());
       _elementAtomVolume_.resize(nElts);
      
      // determine number of atoms in each element, partial sum
      for (int i = 0; i < nLocal; ++i) {
        _elementAtomCountLocal_(atomElement(i,0)) += 1;
      }
      
      // mpi to determine total atoms
      lammpsInterface_->int_allsum(_elementAtomCountLocal_.ptr(),_elementAtomCount_.ptr(),nElts);
      
      // divide element volume by total atoms to get atomic volume
      if (nLocal>0) {
        for (int i = 0; i < nElts; ++i) {
          
          double minx, maxx, miny, maxy, minz, maxz;
          feMesh_->element_coordinates(i,_nodalCoordinates_);
          minx = _nodalCoordinates_(0,0); maxx = _nodalCoordinates_(0,0);
          miny = _nodalCoordinates_(1,0); maxy = _nodalCoordinates_(1,0);
          minz = _nodalCoordinates_(2,0); maxz = _nodalCoordinates_(2,0);
          for (int j = 1; j < _nodalCoordinates_.nCols(); ++j) {
            if (_nodalCoordinates_(0,j)<minx) minx = _nodalCoordinates_(0,j);
            if (_nodalCoordinates_(0,j)>maxx) maxx = _nodalCoordinates_(0,j);
            if (_nodalCoordinates_(1,j)<miny) miny = _nodalCoordinates_(1,j);
            if (_nodalCoordinates_(1,j)>maxy) maxy = _nodalCoordinates_(1,j);
            if (_nodalCoordinates_(2,j)<minz) minz = _nodalCoordinates_(2,j);
            if (_nodalCoordinates_(2,j)>maxz) maxz = _nodalCoordinates_(2,j);
          }
          double eltVol = (maxx-minx)*(maxy-miny)*(maxz-minz);
          if (eltVol<0) eltVol *= -1.;
          _elementAtomVolume_(i) = eltVol/_elementAtomCount_(i);
        }
        
        for (int i = 0; i < nLocal; ++i)
          quantity_(i,i) = _elementAtomVolume_(atomElement(i,0));
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomVolumeRegion
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomVolumeRegion::AtomVolumeRegion(ATC_Method * atc,
                                     DENS_MAN * atomCoarseGrainingPositions,
                                     AtomType atomType) :
    ProtectedAtomDiagonalMatrix<double>(atc,atomType),
    lammpsInterface_(atc->lammps_interface())
  {
    if (!atomCoarseGrainingPositions) {
      InterscaleManager & interscaleManager(atc->interscale_manager());
      atomCoarseGrainingPositions_ = interscaleManager.per_atom_quantity("AtomicCoarseGrainingPositions");
    }
    atomCoarseGrainingPositions_->register_dependence(this);

    // compute volumes and atom counts in each region
    int nregion = lammpsInterface_->nregion();
    regionalAtomVolume_.resize(nregion);
    
    for (int i = 0; i < nregion; ++i) {
      regionalAtomVolume_(i) = 
        (lammpsInterface_->region_xhi(i)-lammpsInterface_->region_xlo(i))
        *(lammpsInterface_->region_yhi(i)-lammpsInterface_->region_ylo(i))
        *(lammpsInterface_->region_zhi(i)-lammpsInterface_->region_zlo(i));
    }
      
    INT_VECTOR localCount(nregion);
    INT_VECTOR globalCount(nregion);  
    // loop over atoms
    localCount = 0;
    const DENS_MAT atomicCoordinates(atomCoarseGrainingPositions->quantity());
    for (int i = 0; i < quantity_.nRows(); ++i) {
      for (int j = 0; j < nregion; ++j) {
        if (lammpsInterface_->region_match(j,
                                           atomicCoordinates(i,0),
                                           atomicCoordinates(i,1),
                                           atomicCoordinates(i,2))) {
          localCount(j)++;
        }
      }
    }
    // communication to get total atom counts per group
    lammpsInterface_->int_allsum(localCount.ptr(),
                                 globalCount.ptr(),nregion);
        
    for (int i = 0; i < nregion; ++i) {
      if (globalCount(i) > 0)
        regionalAtomVolume_(i) /= globalCount(i);
      else
        regionalAtomVolume_(i) = 0;
    }
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  AtomVolumeRegion::~AtomVolumeRegion()
  {
    atomCoarseGrainingPositions_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomVolumeRegion::reset() const
  {
    if (need_reset()) {
      PerAtomDiagonalMatrix<double>::reset();
      const DENS_MAT & atomicCoordinates(atomCoarseGrainingPositions_->quantity());  
      for (int i = 0; i < quantity_.nRows(); i++) {
        for (int iregion = 0; iregion < lammpsInterface_->nregion(); iregion++) {
          if (lammpsInterface_->region_match(iregion,
                                             atomicCoordinates(i,0),
                                             atomicCoordinates(i,0),
                                             atomicCoordinates(i,2)))
            quantity_(i,i) = regionalAtomVolume_(iregion);
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomVolumeFile
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomVolumeFile::AtomVolumeFile(ATC_Method * atc,
                                 const string & atomVolumeFile,
                                 AtomType atomType) :
    ProtectedAtomDiagonalMatrix<double>(atc,atomType),
    atomVolumeFile_(atomVolumeFile),
    lammpsInterface_(atc->lammps_interface())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomVolumeFile::reset() const
  {
    if (need_reset()) {
      PerAtomDiagonalMatrix<double>::reset();
      int nlocal = lammpsInterface_->nlocal();
      ifstream in;
      const int lineSize = 256;
      char line[lineSize];
      const char* fname = &atomVolumeFile_[0]; 

      // create tag to local id map for this processor
      map <int,int> tag2id;
      map <int,int>::const_iterator itr;
      const int * atom_tag = lammpsInterface_->atom_tag();
      for (int i = 0; i < nlocal; ++i) {
        tag2id[atom_tag[i]] = i;
      }

      // get number of atoms
      int natoms = 0;
      if (ATC::LammpsInterface::instance()->rank_zero()) {
        in.open(fname);
        string msg;
        string name = atomVolumeFile_;
        msg = "no "+name+" atomic weights file found";
        if (! in.good()) throw ATC_Error(msg);
        in.getline(line,lineSize); // header
        in.getline(line,lineSize); // blank line
        in.getline(line,lineSize); // number of atoms 
        stringstream inss (line,stringstream::in | stringstream::out);
        inss >> natoms; // number of atoms
        stringstream ss; 
        ss << " found " << natoms << " atoms in atomic weights file";
        lammpsInterface_->print_msg(ss.str());
        if (natoms != lammpsInterface_->natoms()) { 
          throw ATC_Error("Incorrect number of atomic weights read-in!");
        }
        in.getline(line,lineSize); // blank line
      }

      // read and assign atomic weights
      int nread = 0, tag = -1, count = 0;
      double atomic_weight;
      while (nread < natoms) {
        if (ATC::LammpsInterface::instance()->rank_zero()) {
          in.getline(line,lineSize);
          stringstream ss (line,stringstream::in | stringstream::out);
          ss >> tag >> atomic_weight; 
          nread++;
        }
        lammpsInterface_->int_broadcast(&nread);
        lammpsInterface_->int_broadcast(&tag);
        lammpsInterface_->broadcast(&atomic_weight);
        itr = tag2id.find(tag);
        if (itr != tag2id.end()) {
          int iatom = itr->second;
          quantity_(iatom,0) = atomic_weight;
          count++;
        }
      }
      if (lammpsInterface_->rank_zero()) {
        in.close();
        stringstream ss; 
        ss << " read  " << nread << " atomic weights";
        lammpsInterface_->print_msg(ss.str());
      }
      if (count != nlocal) 
        throw ATC_Error("reset "+to_string(count)+" atoms vs "+to_string(nlocal));
    }
  }

  // need to add capability to take in group bit (JAT, 04/02/11)
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomNumber
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomNumber::AtomNumber(ATC_Method * atc, AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,1,atomType),
    atc_(atc)  
  {
  }

  void AtomNumber::reset() const
  {
    int nlocal = atc_->nlocal();
    quantity_.reset(nlocal,1);
    quantity_ = 1;
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomTypeVector
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomTypeVector::AtomTypeVector(ATC_Method * atc, 
    vector<int> typeList, AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,
    typeList.size(),
    atomType),
    atc_(atc), ntypes_(ATC::LammpsInterface::instance()->ntypes()),
    typeList_(typeList)  
  {
    if (typeList_.size() == 0) throw ATC_Error("type list is empty");
    index_.resize(ntypes_,-1); 
    for (unsigned int i = 0; i < typeList_.size(); i++) {
      index_[typeList_[i]] = i;
    }
  }
  AtomTypeVector::AtomTypeVector(ATC_Method * atc, 
    vector<int> typeList, vector<int> groupList, AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,
    typeList.size()+groupList.size(),
    atomType),
    atc_(atc), ntypes_(ATC::LammpsInterface::instance()->ntypes()),
    typeList_(typeList),
    groupList_(groupList)  
  {
    if ((typeList_.size() == 0) && (groupList_.size() == 0)) throw ATC_Error("type/group lists are empty");
    // reverse map
    index_.resize(ntypes_,-1); 
    for (unsigned int i = 0; i < typeList_.size(); i++) {
      index_[typeList_[i]-1] = i;
    }
  }

  void AtomTypeVector::reset() const
  {
    if (need_reset()) {
      //PerAtomQuantity<double>::reset();
      int nlocal = atc_->nlocal();
      quantity_.reset(nlocal,typeList_.size()+groupList_.size());
      const Array<int> & quantityToLammps = (PerAtomQuantity<double>::atc_).atc_to_lammps_map();
      if (typeList_.size()) { 
        int * type = ATC::LammpsInterface::instance()->atom_type(); 
        for (int i = 0; i < nlocal; i++) {
          int a = quantityToLammps(i); 
          int index = index_[type[a]-1];
          if (index > -1) quantity_(i,index) = 1;
        }
      }
      int index = typeList_.size();
      if (groupList_.size()) { 
        const int * mask = ATC::LammpsInterface::instance()->atom_mask();
        for (unsigned int j = 0; j < groupList_.size(); j++) {
          int group = groupList_[j];
          for (int i = 0; i < nlocal; i++) {
            int a = quantityToLammps(i);
            if (mask[a] & group) quantity_(i,index) = 1;
          }
          index++;
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class XrefWrapper
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  XrefWrapper::XrefWrapper(ATC_Method * atc, AtomType atomType) :
    ProtectedClonedAtomQuantity<double>(atc,atc->nsd(),atomType),
    atc_(atc)
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  lammps_vector
  //--------------------------------------------------------
  double ** XrefWrapper::lammps_vector() const
  {
    return atc_->xref();
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicMassWeightedDisplacement
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomicMassWeightedDisplacement::AtomicMassWeightedDisplacement(ATC_Method * atc,
    PerAtomQuantity<double> * atomPositions,
    PerAtomQuantity<double> * atomMasses,
    PerAtomQuantity<double> * atomReferencePositions,
    AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,atc->nsd(),atomType),
    atomPositions_(atomPositions),
    atomMasses_(atomMasses),
    atomReferencePositions_(atomReferencePositions)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomPositions_)
      atomPositions_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_POSITION,
                                                                       atomType);
    if (!atomMasses_)
      atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                    atomType);
    if (!atomReferencePositions_)
      atomReferencePositions_ = interscaleManager.per_atom_quantity("AtomicCoarseGrainingPositions");

    atomPositions_->register_dependence(this);
    atomMasses_->register_dependence(this);
    atomReferencePositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomicMassWeightedDisplacement::~AtomicMassWeightedDisplacement()
  {
    atomPositions_->remove_dependence(this);
    atomMasses_->remove_dependence(this);
    atomReferencePositions_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomicMassWeightedDisplacement::reset() const
  {
    if (need_reset()) { 
      PerAtomQuantity<double>::reset();
      const DENS_MAT & position(atomPositions_->quantity());
      const DENS_MAT & mass(atomMasses_->quantity());
      const DENS_MAT & refPosition(atomReferencePositions_->quantity());

      // q = m * (x - xref)
      quantity_ = position;
      quantity_ -= refPosition;
      quantity_ *= mass;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicMomentum
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomicMomentum::AtomicMomentum(ATC_Method * atc,
                                 PerAtomQuantity<double> * atomVelocities,
                                 PerAtomQuantity<double> * atomMasses,
                                 AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,atc->nsd(),atomType),
    atomVelocities_(atomVelocities),
    atomMasses_(atomMasses)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomVelocities_)
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                    atomType);
    if (!atomMasses_)
      atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                    atomType);

    atomVelocities_->register_dependence(this);
    atomMasses_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomicMomentum::~AtomicMomentum()
  {
    atomVelocities_->remove_dependence(this);
    atomMasses_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomicMomentum::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & mass(atomMasses_->quantity());
      const DENS_MAT & velocity(atomVelocities_->quantity());

      // q = m * v
      quantity_ = velocity;
      quantity_ *= mass;    
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TwiceKineticEnergy
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  TwiceKineticEnergy::TwiceKineticEnergy(ATC_Method * atc,
                                         PerAtomQuantity<double> * atomVelocities,
                                         PerAtomQuantity<double> * atomMasses,
                                         AtomType atomType) :
    AtomicEnergyForTemperature(atc,atomType),
    atomVelocities_(atomVelocities),
    atomMasses_(atomMasses)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomVelocities_)
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                        atomType);
    if (!atomMasses_)
      atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                    atomType);

    atomVelocities_->register_dependence(this);
    atomMasses_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  TwiceKineticEnergy::~TwiceKineticEnergy()
  {
    atomVelocities_->remove_dependence(this);
    atomMasses_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void TwiceKineticEnergy::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & mass(atomMasses_->quantity());
      const DENS_MAT & velocity(atomVelocities_->quantity());
      
      // q = m * (v dot v)
      for (int i = 0; i < quantity_.nRows(); i++) {
          quantity_(i,0) = velocity(i,0)*velocity(i,0);
        for (int j = 1; j < velocity.nCols(); j++) {
          quantity_(i,0) += velocity(i,j)*velocity(i,j);
        }
        quantity_(i,0) *= mass(i,0);
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FluctuatingVelocity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  FluctuatingVelocity::FluctuatingVelocity(ATC_Method * atc,
    PerAtomQuantity<double> * atomVelocities,
    PerAtomQuantity<double> * atomMeanVelocities,
    AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,3,atomType),
    atomVelocities_(atomVelocities),
    atomMeanVelocities_(atomMeanVelocities)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomVelocities_) {
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(
        LammpsInterface::ATOM_VELOCITY, atomType);
    }
    if (!atomMeanVelocities_) {
      atomMeanVelocities_ = interscaleManager.per_atom_quantity(field_to_prolongation_name(VELOCITY));
    }

    atomVelocities_->register_dependence(this);
    atomMeanVelocities_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  FluctuatingVelocity::~FluctuatingVelocity()
  {
    atomVelocities_->remove_dependence(this);
    atomMeanVelocities_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void FluctuatingVelocity::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & velocity(atomVelocities_->quantity());
      const DENS_MAT & meanVelocity(atomMeanVelocities_->quantity());
      
      // q = m * (v dot v)
      for (int i = 0; i < quantity_.nRows(); i++) {
        for (int j = 0; j < velocity.nCols(); j++) {
          quantity_(i,j) = velocity(i,j) - meanVelocity(i,j);
        }
      }
    }
  }
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class ChargeVelocity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  ChargeVelocity::ChargeVelocity(ATC_Method * atc,
    PerAtomQuantity<double> * atomVelocities,
    FundamentalAtomQuantity * atomCharge,
    AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,3,atomType),
    fluctuatingVelocities_(atomVelocities),
    atomCharge_(atomCharge)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!fluctuatingVelocities_) {
      fluctuatingVelocities_ = interscaleManager.per_atom_quantity("AtomicFluctuatingVelocity");
    }
    if (!atomCharge_) {
      atomCharge_ = interscaleManager.fundamental_atom_quantity(
        LammpsInterface::ATOM_CHARGE, atomType);
    }
    fluctuatingVelocities_->register_dependence(this);
    atomCharge_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  ChargeVelocity::~ChargeVelocity()
  {
    fluctuatingVelocities_->remove_dependence(this);
    atomCharge_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void ChargeVelocity::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & velocity(fluctuatingVelocities_->quantity());
      const DENS_MAT & charge(atomCharge_->quantity());
      for (int i = 0; i < quantity_.nRows(); i++) {
        for (int j = 0; j < velocity.nCols(); j++) {
          quantity_(i,j) = charge(i,0)*velocity(i,j);
        }
      }
    }
  }
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SpeciesVelocity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  SpeciesVelocity::SpeciesVelocity(ATC_Method * atc,
    PerAtomQuantity<double> * fluctuatingVelocities,
    PerAtomQuantity<double> * atomTypeVector,
    AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,3*(atomTypeVector->nCols()),atomType),
    fluctuatingVelocities_(fluctuatingVelocities),
    atomTypeVector_(atomTypeVector)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!fluctuatingVelocities_) {
      fluctuatingVelocities_ = interscaleManager.per_atom_quantity("AtomicFluctuatingVelocity");
    }
    if (!atomTypeVector_) {
      atomTypeVector_ = interscaleManager.per_atom_quantity("AtomicTypeVector");
    }
    fluctuatingVelocities_->register_dependence(this);
    atomTypeVector_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  SpeciesVelocity::~SpeciesVelocity()
  {
    fluctuatingVelocities_->remove_dependence(this);
    atomTypeVector_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void SpeciesVelocity::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & velocity(fluctuatingVelocities_->quantity());
      const DENS_MAT & tv(atomTypeVector_->quantity());
      
      for (int i = 0; i < quantity_.nRows(); i++) {
        int m = 0;
        for (int j = 0; j < velocity.nCols(); j++) {
          for (int k = 0; j < tv.nCols(); j++) {
            quantity_(i,m++) = tv(i,k)*velocity(i,j);
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TwiceFluctuatingKineticEnergy
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  TwiceFluctuatingKineticEnergy::TwiceFluctuatingKineticEnergy(ATC_Method * atc,
    PerAtomQuantity<double> * atomVelocities,
    PerAtomQuantity<double> * atomMasses,
    PerAtomQuantity<double> * atomMeanVelocities,
    AtomType atomType) :
    AtomicEnergyForTemperature(atc,atomType),
    atomVelocities_(atomVelocities),
    atomMasses_(atomMasses),
    atomMeanVelocities_(atomMeanVelocities)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomVelocities_) {
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                    atomType);
    }
    if (!atomMasses_) {
      atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                atomType);
    }
    if (!atomMeanVelocities_) {
      atomMeanVelocities_ = interscaleManager.per_atom_quantity(field_to_prolongation_name(VELOCITY));
    }

    atomVelocities_->register_dependence(this);
    atomMasses_->register_dependence(this);
    atomMeanVelocities_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  TwiceFluctuatingKineticEnergy::~TwiceFluctuatingKineticEnergy()
  {
    atomVelocities_->remove_dependence(this);
    atomMasses_->remove_dependence(this);
    atomMeanVelocities_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void TwiceFluctuatingKineticEnergy::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & mass(atomMasses_->quantity());
      const DENS_MAT & velocity(atomVelocities_->quantity());
      const DENS_MAT & meanVelocity(atomMeanVelocities_->quantity());
      
      // q = m * (v dot v)
      double vRel;
      for (int i = 0; i < quantity_.nRows(); i++) {
        vRel = velocity(i,0) - meanVelocity(i,0);
        quantity_(i,0) = vRel*vRel;
        for (int j = 1; j < velocity.nCols(); j++) {
          vRel = velocity(i,j) - meanVelocity(i,j);
          quantity_(i,0) += vRel*vRel;;
        }
        quantity_(i,0) *= mass(i,0);
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class KineticTensor
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  KineticTensor::KineticTensor(ATC_Method * atc,
     PerAtomQuantity<double> * atomVelocities,
     PerAtomQuantity<double> * atomMasses,
     AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,6,atomType),
    atomVelocities_(atomVelocities),
    atomMasses_(atomMasses)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomVelocities_) {
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                    atomType);
    }
    if (!atomMasses_) {
      atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                atomType);
    }

    atomVelocities_->register_dependence(this);
    atomMasses_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  KineticTensor::~KineticTensor()
  {
    atomVelocities_->remove_dependence(this);
    atomMasses_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void KineticTensor::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & mass(atomMasses_->quantity());
      const DENS_MAT & velocity(atomVelocities_->quantity());
      
      // K = m * (v \otimes v)
      for (int i = 0; i < quantity_.nRows(); i++) {
        double m = mass(i,0);
        double v[3] = {velocity(i,0),velocity(i,1),velocity(i,2)};
        quantity_(i,0) -= m*v[0]*v[0];
        quantity_(i,1) -= m*v[1]*v[1];
        quantity_(i,2) -= m*v[2]*v[2];
        quantity_(i,3) -= m*v[0]*v[1];
        quantity_(i,4) -= m*v[0]*v[2];
        quantity_(i,5) -= m*v[1]*v[2];
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FluctuatingKineticTensor
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  FluctuatingKineticTensor::FluctuatingKineticTensor(ATC_Method * atc,
     PerAtomQuantity<double> * atomVelocities,
     PerAtomQuantity<double> * atomMasses,
     PerAtomQuantity<double> * atomMeanVelocities,
     AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,6,atomType),
    atomVelocities_(atomVelocities),
    atomMasses_(atomMasses)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomVelocities_) {
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                    atomType);
    }
    if (!atomMasses_) {
      atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                atomType);
    }
    if (!atomMeanVelocities_) {
      atomMeanVelocities_ = interscaleManager.per_atom_quantity(field_to_prolongation_name(VELOCITY));
    }

    atomVelocities_->register_dependence(this);
    atomMasses_->register_dependence(this);
    atomMeanVelocities_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  FluctuatingKineticTensor::~FluctuatingKineticTensor()
  {
    atomVelocities_->remove_dependence(this);
    atomMasses_->remove_dependence(this);
    atomMeanVelocities_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void FluctuatingKineticTensor::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & mass(atomMasses_->quantity());
      const DENS_MAT & velocity(atomVelocities_->quantity());
      const DENS_MAT & meanVelocity(atomMeanVelocities_->quantity());
      
      // K = m * (v \otimes v)
      for (int i = 0; i < quantity_.nRows(); i++) {
        double m = mass(i,0);
        double v[3] = {velocity(i,0)-meanVelocity(i,0),
                       velocity(i,1)-meanVelocity(i,0),
                       velocity(i,2)-meanVelocity(i,0)};
        quantity_(i,0) -= m*v[0]*v[0];
        quantity_(i,1) -= m*v[1]*v[1];
        quantity_(i,2) -= m*v[2]*v[2];
        quantity_(i,3) -= m*v[0]*v[1];
        quantity_(i,4) -= m*v[0]*v[2];
        quantity_(i,5) -= m*v[1]*v[2];
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class MixedKePeEnergy
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  MixedKePeEnergy::MixedKePeEnergy(ATC_Method * atc,
                                   double keMultiplier,
                                   double peMultiplier,
                                   PerAtomQuantity<double> * twiceKineticEnergy,
                                   PerAtomQuantity<double> * potentialEnergy,
                                   AtomType atomType) :
    AtomicEnergyForTemperature(atc,atomType),
    keMultiplier_(keMultiplier),
    peMultiplier_(peMultiplier),
    twiceKineticEnergy_(twiceKineticEnergy),
    potentialEnergy_(potentialEnergy)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!twiceKineticEnergy_) {
       twiceKineticEnergy_ = interscaleManager.per_atom_quantity("AtomicTwiceKineticEnergy");
    }
    if (!potentialEnergy_) {
      potentialEnergy_ = interscaleManager.per_atom_quantity("AtomicFluctuatingPotentialEnergy");
    }

    twiceKineticEnergy_->register_dependence(this);
    potentialEnergy_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  MixedKePeEnergy::~MixedKePeEnergy()
  {
    twiceKineticEnergy_->remove_dependence(this);
    potentialEnergy_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void MixedKePeEnergy::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & twoKe(twiceKineticEnergy_->quantity());
      const DENS_MAT & pe(potentialEnergy_->quantity());

      // q = peScale * pe + keScale * ke
      quantity_ = pe;
      quantity_ *= peMultiplier_;
      quantity_ += (keMultiplier_/2.)*twoKe;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class TotalEnergy
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  TotalEnergy::TotalEnergy(ATC_Method * atc,
                           PerAtomQuantity<double> * twiceKineticEnergy,
                           PerAtomQuantity<double> * potentialEnergy,
                           AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,atomType),
    twiceKineticEnergy_(twiceKineticEnergy),
    potentialEnergy_(potentialEnergy)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());

    if (!twiceKineticEnergy_) {
       twiceKineticEnergy_ = interscaleManager.per_atom_quantity("AtomicTwiceKineticEnergy");
    }
    if (!potentialEnergy_) {
      potentialEnergy_ = interscaleManager.per_atom_quantity("AtomicPotentialEnergy");
    }
    twiceKineticEnergy_->register_dependence(this);
    potentialEnergy_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  TotalEnergy::~TotalEnergy()
  {
    twiceKineticEnergy_->remove_dependence(this);
    potentialEnergy_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void TotalEnergy::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & twoKe(twiceKineticEnergy_->quantity());
      const DENS_MAT & pe(potentialEnergy_->quantity());
      quantity_ = pe;
      quantity_ += (0.5)*twoKe;
    }
  }
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class FluctuatingPotentialEnergy
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  FluctuatingPotentialEnergy::FluctuatingPotentialEnergy(ATC_Method * atc,
                                                         PerAtomQuantity<double> * potentialEnergy,
                                                         PerAtomQuantity<double> * referencePotential,
                                                         AtomType atomType) :
    AtomicEnergyForTemperature(atc,atomType),
    potentialEnergy_(potentialEnergy),
    referencePotential_(referencePotential)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!potentialEnergy_) {
      potentialEnergy_ = interscaleManager.per_atom_quantity("AtomicPotentialEnergy");
    }
    if (!referencePotential_) {
       referencePotential_ = interscaleManager.per_atom_quantity("AtomicReferencePotential");
    }
    
    potentialEnergy_->register_dependence(this);
    referencePotential_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  FluctuatingPotentialEnergy::~FluctuatingPotentialEnergy()
  {
    potentialEnergy_->remove_dependence(this);
    referencePotential_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void FluctuatingPotentialEnergy::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & pe(potentialEnergy_->quantity());
      const DENS_MAT & refPe(referencePotential_->quantity());

      quantity_ = pe;
      quantity_ -= refPe;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class DotTwiceKineticEnergy
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  DotTwiceKineticEnergy::DotTwiceKineticEnergy(ATC_Method * atc,
                                               PerAtomQuantity<double> * atomForces,
                                               PerAtomQuantity<double> * atomVelocities,
                                               AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,1,atomType),
    atomForces_(atomForces),
    atomVelocities_(atomVelocities)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomForces_)
      atomForces_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_FORCE,
                                                                    atomType);
    if (!atomVelocities_)
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                        atomType);

    atomForces_->register_dependence(this);
    atomVelocities_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  DotTwiceKineticEnergy::~DotTwiceKineticEnergy()
  {
    atomForces_->remove_dependence(this);
    atomVelocities_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void DotTwiceKineticEnergy::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & velocity(atomVelocities_->quantity());
      const DENS_MAT & force(atomForces_->quantity());

      for (int i = 0; i < quantity_.nRows(); i++) {
        quantity_(i,0) = velocity(i,0)*force(i,0);
        for (int j = 1; j < velocity.nCols(); j++) {
          quantity_(i,0) += velocity(i,j)*force(i,j);
        }
        quantity_(i,0) *= 2.;
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class VelocitySquared
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  VelocitySquared::VelocitySquared(ATC_Method * atc,
                                   PerAtomQuantity<double> * atomVelocities,
                                   AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,1,atomType),
    atomVelocities_(atomVelocities)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomVelocities_)
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                    atomType);
    atomVelocities_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  VelocitySquared::~VelocitySquared()
  {
    atomVelocities_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void VelocitySquared::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & velocity(atomVelocities_->quantity());

      for (int i = 0; i < quantity_.nRows(); i++) {
        quantity_(i,0) = velocity(i,0)*velocity(i,0);
        for (int j = 1; j < velocity.nCols(); j++) {
          quantity_(i,0) += velocity(i,j)*velocity(i,j);
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaSquared
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  LambdaSquared::LambdaSquared(ATC_Method * atc,
                               PerAtomQuantity<double> * atomMasses,
                               PerAtomQuantity<double> * atomVelocitiesSquared,
                               PerAtomQuantity<double> * atomLambdas,
                               AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,1,atomType),
    atomMasses_(atomMasses),
    atomVelocitiesSquared_(atomVelocitiesSquared),
    atomLambdas_(atomLambdas)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomMasses_) {
      atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                    atomType);
    }
    atomMasses_->register_dependence(this);

    if (!atomVelocitiesSquared) {
      atomVelocitiesSquared_ = interscaleManager.per_atom_quantity("LambdaSquared");
    }
    atomVelocitiesSquared_->register_dependence(this);

    if (!atomLambdas_) {
      atomLambdas_ = interscaleManager.per_atom_quantity("AtomLambdaEnergy");
    }
    atomLambdas_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  LambdaSquared::~LambdaSquared()
  {
    atomMasses_->remove_dependence(this);
    atomVelocitiesSquared_->remove_dependence(this);
    atomLambdas_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void LambdaSquared::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & mass(atomMasses_->quantity());
      const DENS_MAT & velocitySquared(atomVelocitiesSquared_->quantity());
      const DENS_MAT & lambda(atomLambdas_->quantity());

      for (int i = 0; i < quantity_.nRows(); i++) {
        quantity_(i,0) = lambda(i,0)*lambda(i,0)*velocitySquared(i,0)/mass(i,0);
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomToType
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomToType::AtomToType(ATC_Method * atc,
                         int type,
                         AtomType atomType) :
    LargeToSmallAtomMap(atc,atomType),
    type_(type)
  {
    // DO NOTHING
    
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomToType::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<int>::reset();
      quantity_ = -1.;
      size_ = 0;
      const int * type = lammpsInterface_->atom_type();

      const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
      for (int i = 0; i < quantity_.nRows(); ++i) {
        int atomIdx = quantityToLammps(i);
        if (type[atomIdx] == type_) {
          quantity_(i,0) = size_;
          ++size_;
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomToGroup
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomToGroup::AtomToGroup(ATC_Method * atc,
                           int group,
                           AtomType atomType) :
    LargeToSmallAtomMap(atc,atomType),
    group_(group)
  {
    // DO NOTHING
    
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomToGroup::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<int>::reset();
      quantity_ = -1.;
      size_ = 0;
      const int * mask = lammpsInterface_->atom_mask();

      const Array<int> & quantityToLammps = atc_.atc_to_lammps_map();
      for (int i = 0; i < quantity_.nRows(); ++i) {
        int atomIdx = quantityToLammps(i);
        if (mask[atomIdx] & group_) {
          quantity_(i,0) = size_;
          ++size_;
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomToNodeset
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomToNodeset::AtomToNodeset(ATC_Method * atc,
                               SetDependencyManager<int> * subsetNodes,
                               PerAtomQuantity<int> * atomElement,
                               AtomType atomType) :
    LargeToSmallAtomMap(atc,atomType),
    subsetNodes_(subsetNodes),
    atomElement_(atomElement),
    feMesh_((atc->fe_engine())->fe_mesh())
  {
    if (!atomElement_) {
      atomElement_ = (atc->interscale_manager()).per_atom_int_quantity("AtomElement");
    }
    if (atomElement_) {
      atomElement_->register_dependence(this);
    }
    else {
      throw ATC_Error("PerAtomQuantityLibrary::AtomToRegulated - No AtomElement provided");
    }

    subsetNodes_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  quantity
  //--------------------------------------------------------
  void AtomToNodeset::reset() const
  {
    
    //so it has been commented out.
    /*
    if (need_reset()) {
      PerAtomQuantity<int>::reset();
      const INT_ARRAY & atomElement(atomElement_->quantity());
      const set<int> & subsetNodes(subsetNodes_->quantity());
      int nLocal = atomElement.nRows();
      quantity_.resize(nLocal,1);
      quantity_ = -1;
      size_ = 0;

      for (int i = 0; i < nLocal; ++i) {
        feMesh_->element_connectivity_unique(atomElement(i,0),_nodes_);
        for (int j = 0; j < _nodes_.size(); ++j) {
          if (subsetNodes.find(_nodes_(j)) != subsetNodes.end()) {
            quantity_(i,0) = size_;
            size_++;
            break;
          }
        }
      }
    }
    */
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomToElementset
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  AtomToElementset::AtomToElementset(ATC_Method * atc,
                                     MatrixDependencyManager<DenseMatrix, bool> * elementMask,
                                     PerAtomQuantity<int> * atomElement,
                                     AtomType atomType) :
    LargeToSmallAtomMap(atc,atomType),
    elementMask_(elementMask),
    atomElement_(atomElement),
    feMesh_((atc->fe_engine())->fe_mesh())
  {
    if (!atomElement_) {
      atomElement_ = (atc->interscale_manager()).per_atom_int_quantity("AtomElement");
    }
    if (atomElement_) {
      atomElement_->register_dependence(this);
    }
    else {
      throw ATC_Error("PerAtomQuantityLibrary::AtomToRegulated - No AtomElement provided");
    }

    elementMask_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomToElementset::~AtomToElementset()
  {
    atomElement_->remove_dependence(this);
    elementMask_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  quantity
  //--------------------------------------------------------
  void AtomToElementset::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<int>::reset();
      const INT_ARRAY & atomElement(atomElement_->quantity());
      const DenseMatrix<bool> & elementMask(elementMask_->quantity());
      int nLocal = atomElement.nRows();
      quantity_.resize(nLocal,1);
      quantity_ = -1;
      size_ = 0;

      for (int i = 0; i < nLocal; ++i) {
        if (elementMask(atomElement(i,0),0)) {
          quantity_(i,0) = size_;
          size_++;
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class MappedAtomQuantity
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  MappedAtomQuantity::MappedAtomQuantity(ATC_Method * atc,
                                         PerAtomQuantity<double> * source,
                                         LargeToSmallAtomMap * map,
                                         AtomType atomType) :
    ProtectedMappedAtomQuantity<double>(atc,map,source->nCols(),atomType),
    source_(source),
    map_(map)
  {
    source_->register_dependence(this);
    map_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void MappedAtomQuantity::reset() const
  {
    if (needReset_) {
      ProtectedMappedAtomQuantity<double>::reset();
      const DENS_MAT & source(source_->quantity());
      const INT_ARRAY & map(map_->quantity());
      quantity_.resize(map_->size(),source.nCols());
      
      for (int i = 0; i < source.nRows(); i++) {
        int idx = map(i,0);
        if (idx > -1) {
          for (int j = 0; j < source.nCols(); j++) {
            quantity_(idx,j) = source(i,j);
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class VelocitySquaredMapped
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  VelocitySquaredMapped::VelocitySquaredMapped(ATC_Method * atc,
                                               MatrixDependencyManager<DenseMatrix, int> * atomMap,
                                               PerAtomQuantity<double> * atomVelocities,
                                               AtomType atomType) :
    ProtectedMappedAtomQuantity<double>(atc,atomMap,1,atomType),
    atomVelocities_(atomVelocities)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomVelocities_)
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                        atomType);
    atomVelocities_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  VelocitySquaredMapped::~VelocitySquaredMapped()
  {
    atomVelocities_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void VelocitySquaredMapped::reset() const
  {
    if (need_reset()) {
      ProtectedMappedAtomQuantity<double>::reset();
      const DENS_MAT & velocity(atomVelocities_->quantity());
      const INT_ARRAY & atomMap(atomMap_->quantity());

      for (int i = 0; i < atomMap.nRows(); i++) {
        int idx = atomMap(i,0);
        if (idx > -1) {
          quantity_(idx,0) = velocity(i,0)*velocity(i,0);
          for (int j = 1; j < velocity.nCols(); j++) {
            quantity_(idx,0) += velocity(i,j)*velocity(i,j);
          }
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaSquaredMapped
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  LambdaSquaredMapped::LambdaSquaredMapped(ATC_Method * atc,
                                           MatrixDependencyManager<DenseMatrix, int> * atomMap,
                                           PerAtomQuantity<double> * atomMasses,
                                           PerAtomQuantity<double> * atomVelocitiesSquared,
                                           PerAtomQuantity<double> * atomLambdas,
                                           AtomType atomType) :
    ProtectedMappedAtomQuantity<double>(atc,atomMap,1,atomType),
    atomMasses_(atomMasses),
    atomVelocitiesSquared_(atomVelocitiesSquared),
    atomLambdas_(atomLambdas)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomMasses_) {
      atomMasses_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                                    atomType);
    }
    atomMasses_->register_dependence(this);

    if (!atomVelocitiesSquared) {
      atomVelocitiesSquared_ = interscaleManager.per_atom_quantity("LambdaSquared");
    }
    atomVelocitiesSquared_->register_dependence(this);

    if (!atomLambdas_) {
      atomLambdas_ = interscaleManager.per_atom_quantity("AtomLambdaEnergy");
    }
    atomLambdas_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  LambdaSquaredMapped::~LambdaSquaredMapped()
  {
    atomMasses_->remove_dependence(this);
    atomVelocitiesSquared_->remove_dependence(this);
    atomLambdas_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void LambdaSquaredMapped::reset() const
  {
    if (need_reset()) {
      ProtectedMappedAtomQuantity<double>::reset();
      const DENS_MAT & mass(atomMasses_->quantity());
      const DENS_MAT & velocitySquared(atomVelocitiesSquared_->quantity());
      const DENS_MAT & lambda(atomLambdas_->quantity());
      const INT_ARRAY & atomMap(atomMap_->quantity());

      for (int i = 0; i < atomMap.nRows(); i++) {
        int idx = atomMap(i,0);
        if (idx > -1) {
          quantity_(idx,0) = lambda(i,0)*lambda(i,0)*velocitySquared(idx,0)/mass(i,0);
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class HeatCapacity
  //    computes the classical atomic heat capacity
  //--------------------------------------------------------
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  HeatCapacity::HeatCapacity(ATC_Method * atc, AtomType atomType) :
    ConstantQuantity<double>(atc,1,1,atomType)
  {
    constant_ = (atc->nsd())*(lammpsInterface_->kBoltzmann());
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicVelocityRescaleFactor
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomicVelocityRescaleFactor::AtomicVelocityRescaleFactor(ATC_Method * atc,
                                                           PerAtomQuantity<double> * atomLambdas,
                                                           AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,1,atomType),
    atomLambdas_(atomLambdas)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomLambdas) {
      atomLambdas_ = interscaleManager.per_atom_quantity("AtomLambdaEnergy");
    }
    atomLambdas_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomicVelocityRescaleFactor::~AtomicVelocityRescaleFactor()
  {
    atomLambdas_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomicVelocityRescaleFactor::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & lambda(atomLambdas_->quantity());

      for (int i = 0; i < quantity_.nRows(); i++) {
        quantity_(i,0) = sqrt(lambda(i,0));
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicFluctuatingVelocityRescaled
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomicFluctuatingVelocityRescaled::AtomicFluctuatingVelocityRescaled(ATC_Method * atc,
                                                                       PerAtomQuantity<double> * atomRescaleFactor,
                                                                       PerAtomQuantity<double> * atomFluctuatingVelocity,
                                                                       AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,atc->nsd(),atomType),
    atomRescaleFactor_(atomRescaleFactor),
    atomFluctuatingVelocity_(atomFluctuatingVelocity)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomRescaleFactor_) {
      atomRescaleFactor_ = interscaleManager.per_atom_quantity("AtomVelocityRescaling");
    }
    if (!atomFluctuatingVelocity_) {
      atomFluctuatingVelocity_ = interscaleManager.per_atom_quantity("AtomicFluctuatingVelocity");
    }

    atomRescaleFactor_->register_dependence(this);
    atomFluctuatingVelocity_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomicFluctuatingVelocityRescaled::~AtomicFluctuatingVelocityRescaled()
  {
    atomRescaleFactor_->remove_dependence(this);
    atomFluctuatingVelocity_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomicFluctuatingVelocityRescaled::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & factor(atomRescaleFactor_->quantity());
      const DENS_MAT & v(atomFluctuatingVelocity_->quantity());

      for (int i = 0; i < quantity_.nRows(); i++) {
        for (int j = 0; j < quantity_.nCols(); j++) {
          quantity_(i,j) = factor(i,0)*v(i,j);
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicCombinedRescaleThermostatError
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomicCombinedRescaleThermostatError::AtomicCombinedRescaleThermostatError(ATC_Method * atc,
                                                                             PerAtomQuantity<double> * atomFluctuatingMomentumRescaled,
                                                                             PerAtomQuantity<double> * atomMeanVelocity,
                                                                             PerAtomQuantity<double> * atomStreamingVelocity,
                                                                             PerAtomQuantity<double> * atomMass,
                                                                             AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,1,atomType),
    atomFluctuatingMomentumRescaled_(atomFluctuatingMomentumRescaled),
    atomMeanVelocity_(atomMeanVelocity),
    atomStreamingVelocity_(atomStreamingVelocity),
    atomMass_(atomMass)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomFluctuatingMomentumRescaled_) {
      atomFluctuatingMomentumRescaled_ = interscaleManager.per_atom_quantity("AtomFluctuatingMomentumRescaled");
    }
    if (!atomMeanVelocity_) {
      atomMeanVelocity_ = interscaleManager.per_atom_quantity(field_to_prolongation_name(VELOCITY));
    }
    if (!atomStreamingVelocity_) {
      atomStreamingVelocity_ = interscaleManager.per_atom_quantity("AtomLambdaMomentum");
    }
    if (!atomMass_) {
      atomMass_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                              atomType);
    }

    atomFluctuatingMomentumRescaled_->register_dependence(this);
    atomMeanVelocity_->register_dependence(this);
    atomStreamingVelocity_->register_dependence(this);
    atomMass_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomicCombinedRescaleThermostatError::~AtomicCombinedRescaleThermostatError()
  {
    atomFluctuatingMomentumRescaled_->remove_dependence(this);
    atomMeanVelocity_->remove_dependence(this);
    atomStreamingVelocity_->remove_dependence(this);
    atomMass_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomicCombinedRescaleThermostatError::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & m(atomMass_->quantity());
      const DENS_MAT & p(atomFluctuatingMomentumRescaled_->quantity());
      const DENS_MAT & v(atomMeanVelocity_->quantity());
      const DENS_MAT & s(atomStreamingVelocity_->quantity());

      double diff;
      for (int i = 0; i < quantity_.nRows(); i++) {
        diff = s(i,0)-v(i,0);
        quantity_(i,0) = 2.*p(i,0)*diff + m(i,0)*diff*diff;
        for (int j = 1; j < p.nCols(); j++) {
          diff = s(i,j)-v(i,j);
          quantity_(i,0) += 2.*p(i,j)*diff + m(i,0)*diff*diff;
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicThermostatForce
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomicThermostatForce::AtomicThermostatForce(ATC_Method * atc,
                                               PerAtomQuantity<double> * atomLambdas,
                                               PerAtomQuantity<double> * atomVelocities,
                                               AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,atc->nsd(),atomType),
    atomLambdas_(atomLambdas),
    atomVelocities_(atomVelocities)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomLambdas) {
      atomLambdas_ = interscaleManager.per_atom_quantity("AtomLambdaEnergy");
    }
    if (!atomVelocities_) {
      atomVelocities_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_VELOCITY,
                                                                    atomType);
    }
    
    atomLambdas_->register_dependence(this);
    atomVelocities_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomicThermostatForce::~AtomicThermostatForce()
  {
    atomLambdas_->remove_dependence(this);
    atomVelocities_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomicThermostatForce::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & v(atomVelocities_->quantity());
      const DENS_MAT & lambda(atomLambdas_->quantity());
      
      // force = -lambda*v
      quantity_ = v;
      quantity_ *= lambda;
      quantity_ *= -1.;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicKinetostatForceDisplacement
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomicKinetostatForceDisplacement::AtomicKinetostatForceDisplacement(ATC_Method * atc,
                                               PerAtomQuantity<double> * atomLambda,
                                               PerAtomQuantity<double> * atomMass,
                                               AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,atc->nsd(),atomType),
    atomLambda_(atomLambda),
    atomMass_(atomMass)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomLambda) {
      atomLambda_ = interscaleManager.per_atom_quantity("AtomLambdaMomentum");
    }
    if (!atomMass_) {
      atomMass_ = interscaleManager.fundamental_atom_quantity(LammpsInterface::ATOM_MASS,
                                                              atomType);
    }
    
    atomLambda_->register_dependence(this);
    atomMass_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomicKinetostatForceDisplacement::~AtomicKinetostatForceDisplacement()
  {
    atomLambda_->remove_dependence(this);
    atomMass_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomicKinetostatForceDisplacement::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & m(atomMass_->quantity());
      const DENS_MAT & lambda(atomLambda_->quantity());
      double timeFactor = time_step_factor(0.5*atc_.dt());
      
      //force = -lambda*m*(timestep factor)
      quantity_ = lambda;
      quantity_ *= m;
      quantity_ *= -1.*timeFactor;
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class AtomicKinetostatForceStress
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  AtomicKinetostatForceStress::AtomicKinetostatForceStress(ATC_Method * atc,
                                                           PerAtomQuantity<double> * atomLambda,
                                                           AtomType atomType) :
    ProtectedAtomQuantity<double>(atc,atc->nsd(),atomType),
    atomLambda_(atomLambda)
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomLambda_) {
      atomLambda_ = interscaleManager.per_atom_quantity("AtomLambdaMomentum");
    }
    
    atomLambda_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  AtomicKinetostatForceStress::~AtomicKinetostatForceStress()
  {
    atomLambda_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void AtomicKinetostatForceStress::reset() const
  {
    if (need_reset()) {
      PerAtomQuantity<double>::reset();
      const DENS_MAT & lambda(atomLambda_->quantity());
      
      // force = -lambda
      quantity_ = lambda;
      quantity_ *= -1.;
    }
  }
  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class PerAtomKernelFunction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  PerAtomKernelFunction::PerAtomKernelFunction(ATC_Method * atc,
                                               PerAtomQuantity<double> * atomPositions,
                                               AtomType atomType) :
    ProtectedAtomSparseMatrix<double>(atc,(atc->fe_engine())->num_nodes(),atc->accumulant_bandwidth(),atomType),
    atomPositions_(atomPositions),
    feEngine_(atc->fe_engine())
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomPositions_) {
      atomPositions_ = interscaleManager.per_atom_quantity("AtomicCoarseGrainingPositions");
    }
    atomPositions_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  PerAtomKernelFunction::~PerAtomKernelFunction()
  {
    atomPositions_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void PerAtomKernelFunction::reset() const
  {
    if (need_reset()) {
      PerAtomSparseMatrix<double>::reset();
      const DENS_MAT & positions(atomPositions_->quantity());
      if (positions.nRows() > 0) {
        feEngine_->evaluate_kernel_functions(positions,quantity_);
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class PerAtomShapeFunction
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  constructor
  //--------------------------------------------------------
  PerAtomShapeFunction::PerAtomShapeFunction(ATC_Method * atc,
                                             PerAtomQuantity<double> * atomPositions,
                                             PerAtomQuantity<int> * atomElements,
                                             AtomType atomType) :
    ProtectedAtomSparseMatrix<double>(atc,atc->num_nodes(),(atc->fe_engine())->num_nodes_per_element(),atomType),
    atomPositions_(atomPositions),
    atomElements_(atomElements),
    feEngine_(atc->fe_engine())
  {
    InterscaleManager & interscaleManager(atc->interscale_manager());
    if (!atomPositions_) {
      atomPositions_ = interscaleManager.per_atom_quantity("AtomicCoarseGrainingPositions");
    }
    if (!atomElements_) {
      atomElements_ = interscaleManager.per_atom_int_quantity("AtomElement");
    }

    atomPositions_->register_dependence(this);
    atomElements_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  destructor
  //--------------------------------------------------------
  PerAtomShapeFunction::~PerAtomShapeFunction()
  {
    atomPositions_->remove_dependence(this);
    atomElements_->remove_dependence(this);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void PerAtomShapeFunction::reset() const
  {
    if (need_reset()) {
      PerAtomSparseMatrix<double>::reset();
      const DENS_MAT & positions(atomPositions_->quantity());
      const INT_ARRAY & elements(atomElements_->quantity());
      if (positions.nRows() > 0) {
        feEngine_->evaluate_shape_functions(positions,
                                            elements,
                                            quantity_);
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LambdaCouplingMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  LambdaCouplingMatrix::LambdaCouplingMatrix(ATC_Method * atc,
                                             MatrixDependencyManager<DenseMatrix, int> * nodeToOverlapMap,
                                             SPAR_MAN * shapeFunction) :
    ProtectedMappedAtomSparseMatrix<double>(atc,nodeToOverlapMap->size(),
                                            (atc->fe_engine())->num_nodes_per_element(),
                                            INTERNAL),
    shapeFunction_(shapeFunction),
    nodeToOverlapMap_(nodeToOverlapMap)
  {
    if (!shapeFunction_) {
      shapeFunction_ = (atc->interscale_manager()).per_atom_sparse_matrix("Interpolant");;
    }
    if (!nodeToOverlapMap_) {
      nodeToOverlapMap_ = (atc->interscale_manager()).dense_matrix_int("NodeToOverlapMap");
    }
    
    shapeFunction_->register_dependence(this);
    nodeToOverlapMap_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void LambdaCouplingMatrix::reset() const
  {
    if (need_reset()) {
      PerAtomSparseMatrix<double>::reset();
      int nNodeOverlap = nodeToOverlapMap_->size();
      const SPAR_MAT & shapeFunction(shapeFunction_->quantity());
      quantity_.reset(shapeFunction.nRows(),nNodeOverlap); // number of atoms X number of nodes overlapping MD region
      const INT_ARRAY nodeToOverlapMap(nodeToOverlapMap_->quantity());
    
      for (int i = 0; i < shapeFunction.size(); ++i) {
        TRIPLET<double> myTriplet = shapeFunction.triplet(i);
        int myCol = nodeToOverlapMap(myTriplet.j,0);
        if (myCol > -1) {
          quantity_.set(myTriplet.i,myCol,myTriplet.v);
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class LocalLambdaCouplingMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  LocalLambdaCouplingMatrix::LocalLambdaCouplingMatrix(ATC_Method * atc,
                                                       MatrixDependencyManager<DenseMatrix, int> * lambdaAtomMap,
                                                       MatrixDependencyManager<DenseMatrix, int> * nodeToOverlapMap,
                                                       SPAR_MAN * shapeFunction) :
    LambdaCouplingMatrix(atc,nodeToOverlapMap,shapeFunction),
    lambdaAtomMap_(lambdaAtomMap)
  {
    if (!lambdaAtomMap_) {
      lambdaAtomMap_ = (atc->interscale_manager()).dense_matrix_int("LambdaAtomMap");
    }
    
    lambdaAtomMap_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void LocalLambdaCouplingMatrix::reset() const
  {
    if (need_reset()) {
      PerAtomSparseMatrix<double>::reset();
      int nNodeOverlap = nodeToOverlapMap_->size();
      int nLocalLambda = lambdaAtomMap_->size();

      quantity_.reset(nLocalLambda,nNodeOverlap); // number of regulated atoms X number of nodes containing them
      const SPAR_MAT & shapeFunction(shapeFunction_->quantity());
      const INT_ARRAY nodeToOverlapMap(nodeToOverlapMap_->quantity());
      const INT_ARRAY lambdaAtomMap(lambdaAtomMap_->quantity());
      
      for (int i = 0; i < shapeFunction.size(); ++i) {
        TRIPLET<double> myTriplet = shapeFunction.triplet(i);
        int myRow = lambdaAtomMap(myTriplet.i,0);
        int myCol = nodeToOverlapMap(myTriplet.j,0);
        if ((myRow > -1) && (myCol > -1)) {
          quantity_.set(myRow,myCol,myTriplet.v);
        }
      }
    }
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class GhostCouplingMatrix
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  GhostCouplingMatrix::GhostCouplingMatrix(ATC_Method * atc,
                                           SPAR_MAN * shapeFunction,
                                           SetDependencyManager<int> * subsetNodes,
                                           MatrixDependencyManager<DenseMatrix, int> * nodeToOverlapMap) :
    LambdaCouplingMatrix(atc,nodeToOverlapMap,shapeFunction),
    subsetNodes_(subsetNodes)
  {
    subsetNodes_->register_dependence(this);
  }

  //--------------------------------------------------------
  //  reset_quantity
  //--------------------------------------------------------
  void GhostCouplingMatrix::reset() const
  {
    if (need_reset()) {
      PerAtomSparseMatrix<double>::reset();
      const SPAR_MAT & shapeFunction(shapeFunction_->quantity());
      const set<int> & subsetNodes(subsetNodes_->quantity());
      int nNodes = nodeToOverlapMap_->nRows();
      int nLocalGhost = shapeFunction.nRows();
      quantity_.reset(nLocalGhost,nNodes);
      const INT_ARRAY nodeToOverlapMap(nodeToOverlapMap_->quantity());
      for (int i = 0; i < shapeFunction.size(); ++i) {
        TRIPLET<double> myTriplet = shapeFunction.triplet(i);
        int myCol = myTriplet.j;
        if (nodeToOverlapMap(myCol,0) > -1) {
          quantity_.set(myTriplet.i,myCol,myTriplet.v);
        }
      }
      quantity_.compress();
      
      //int nNodes = nodeToOverlapMap_->nRows();
      _activeNodes_.reset(nNodes);
      for (int i = 0; i < nNodes; ++i) {
        if (subsetNodes.find(i) != subsetNodes.end()) {
          _activeNodes_(i) = 1.;
        }
      }

      //const SPAR_MAT & shapeFunction(shapeFunction_->quantity());
      //int nLocalGhost = shapeFunction.nRows();
      _weights_ = shapeFunction*_activeNodes_;
      _weightMatrix_.resize(nLocalGhost,nLocalGhost);
      for (int i = 0; i < nLocalGhost; ++i) {
        _weightMatrix_(i,i) = 1./_weights_(i);
      }
      quantity_ = _weightMatrix_*quantity_;
    }
  }
};
