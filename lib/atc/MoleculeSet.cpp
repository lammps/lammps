// ATC_Method headers
#include "MoleculeSet.h"
#include "ATC_Method.h"
#include "LammpsInterface.h"
#include "ATC_Error.h"
#include <queue>
#include <utility>
#include <sstream>

using std::multimap;
using std::map;
using std::pair;
using std::set;
using std::stringstream;

namespace ATC {

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class MoleculeSet
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  MoleculeSet::MoleculeSet(ATC_Method * atc, int groupBit) :
    atc_(atc),
    groupBit_(groupBit),
    lammps_(atc->lammps_interface())
  {
    // do nothing
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  MoleculeSet::~MoleculeSet()
  {
    clear();
  }

  //--------------------------------------------------------
  //  clear
  //--------------------------------------------------------
  void MoleculeSet::clear()
  {
    moleculeToAtoms_.clear();
    localMoleculeToAtoms_.clear();
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void MoleculeSet::initialize(map<int,double> * globalAtomsPerMolecule)
  {
    // determine the total number of molecules in this group
    // essentially ripped Compute::molecules_in_group from lammps

    // find lo/hi molecule ID for any atom in group
    int i;
    int *molecule = lammps_->atom_to_molecule();
    const int *mask = lammps_->atom_mask();
    int nlocal = lammps_->nlocal();

    int lo = lammps_->natoms();
    int hi = -1;
    int flag = 0;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupBit_) {
        if (molecule[i] == 0) flag = 1;
        lo = MIN(lo,molecule[i]);
        hi = MAX(hi,molecule[i]);
      }
    }

    int flagall;
    lammps_->int_allsum(&flag,&flagall);
    if (flagall) throw ATC_Error("Atom with molecule ID = 0 included in atc molecule group");

    int globalLo, globalHi;
    lammps_->int_allmin(&lo,&globalLo);
    lammps_->int_allmax(&hi,&globalHi);
    if (globalLo == lammps_->natoms()) throw ATC_Error("MoleculeSet:initialize - no molecules correspond to the group");

    // molmap = vector of length nlen
    // set to 1 for IDs that appear in group across all procs, else 0

    int nlen = globalHi-globalLo+1;
    int * localCount = new int[nlen];
    for (i = 0; i < nlen; i++) localCount[i] = 0;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupBit_)
        localCount[molecule[i]-globalLo]++;

    int * globalCount = new int[nlen];
    lammps_->int_allsum(localCount,globalCount,nlen);

    // nmolecules = # of non-zero IDs in molmap
    // molmap[i] = index of molecule, skipping molecules not in group with -1

    nMoleculesTotal_ = 0;
    for (i = 0; i < nlen; i++)
      if (globalCount[i]) nMoleculesTotal_++;

    if (globalAtomsPerMolecule) {
      for (i = 0; i < nlen; i++)
        if (globalCount[i]) globalAtomsPerMolecule->insert(pair<int,double>(i+globalLo,double(globalCount[i])));
    }

    // deallocate storage
    delete [] localCount;
    delete [] globalCount;
  }

  //--------------------------------------------------------
  //  atoms_by_global_molecule
  //--------------------------------------------------------
  set<int> MoleculeSet::atoms_by_global_molecule(int id) const
  {
    if (need_reset()) reset();
    typedef multimap<int, set<int> >::const_iterator MMIT;
    typedef set<int>::const_iterator SIT;
    pair<MMIT,MMIT> mol = moleculeToAtoms_.equal_range(id);

    set<int> realAtoms;
    const set<int> * myAtoms;
    for (MMIT molIt = mol.first; molIt != mol.second; molIt++) {
      myAtoms = &(molIt->second);
        for (SIT atomIt = myAtoms->begin(); atomIt != myAtoms->end(); atomIt++)
          realAtoms.insert(*atomIt);
    }
    return realAtoms;
  }

  //--------------------------------------------------------
  //  atoms_by_local_molecule
  //--------------------------------------------------------
  const set<int> & MoleculeSet::atoms_by_local_molecule(int id) const
  {
    if (need_reset()) reset();
    return localMoleculeToAtoms_[id]->second;
  }

  //--------------------------------------------------------
  //  set_local_molecules_to_atoms
  //--------------------------------------------------------
  void MoleculeSet::set_local_molecules_to_atoms() const
  {
    localMoleculeToAtoms_.clear();
    localMoleculeToAtoms_.reserve(moleculeToAtoms_.size());
    multimap<int, set<int> >::const_iterator molecule;
    for (molecule = moleculeToAtoms_.begin(); molecule != moleculeToAtoms_.end(); molecule++)
      localMoleculeToAtoms_.push_back(molecule);
  }

  //--------------------------------------------------------
  //--------------------------------------------------------
  //  Class SmallMoleculeSet
  //--------------------------------------------------------
  //--------------------------------------------------------

  //--------------------------------------------------------
  //  Constructor
  //--------------------------------------------------------
  SmallMoleculeSet::SmallMoleculeSet(ATC_Method * atc, int groupBit,
                                     PerAtomQuantity<int> * bondList, PerAtomQuantity<int> * /* numBond */) :
    MoleculeSet(atc,groupBit),
    bondList_(bondList)
  {
  }

  //--------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------
  SmallMoleculeSet::~SmallMoleculeSet()
  {
  }

  //--------------------------------------------------------
  //  clear
  //--------------------------------------------------------
  void SmallMoleculeSet::clear()
  {
    MoleculeSet::clear();
  }

  //--------------------------------------------------------
  //  initialize
  //--------------------------------------------------------
  void SmallMoleculeSet::initialize(std::map<int, double> * globalAtomsPerMolecule)
  {
    // make sure newton_bond is off, otherwise use large molecule set
    if (lammps_->newton_bond())
      throw ATC_Error("Cannot use newton_bond with small molecules");

    MoleculeSet::initialize(globalAtomsPerMolecule);
  }

  //--------------------------------------------------------
  //  reset
  //--------------------------------------------------------
  void SmallMoleculeSet::reset() const
  {
    stringstream message;
    using std::queue;
    moleculeToAtoms_.clear();
    lammps_->forward_comm_fix();
    int * numBond = lammps_->num_bond();
    int ** bondAtom = lammps_->bond_atom();

    // add in real atoms for molecules
    int *molecule = lammps_->atom_to_molecule();
    const int *mask = lammps_->atom_mask();
    int nlocal = lammps_->nlocal();
    _atomFound_.resize(atc_->nproc_ghost());
    _atomFound_ = false;
    int nmol = 0;

    for (int i = 0; i < nlocal; i++) {
      queue<int> myQueue;
      if ((mask[i] & groupBit_) && !_atomFound_(i)) {
        set<int> myAtoms;
        myAtoms.insert(i);
        _atomFound_(i) = true;
        for (int j = 0; j < numBond[i]; j++) {
          int localIdx = lammps_->local_to_global_map(bondAtom[i][j]);
          if (!_atomFound_(localIdx)) {
            myQueue.push(localIdx);
            _atomFound_(localIdx) = true;
          }
        }
        while (!myQueue.empty()) {
          int myIdx = myQueue.front();
          myQueue.pop();
          myAtoms.insert(myIdx);
          for (int j = 0; j < numBond[myIdx]; j++) {
            int localIdx = lammps_->local_to_global_map(bondAtom[myIdx][j]);
            if (!_atomFound_(localIdx)) {
              myQueue.push(localIdx);
              _atomFound_(localIdx) = true;
            }
          }
        }
        nmol++;
        moleculeToAtoms_.insert(pair<int,set<int> >(molecule[i],myAtoms));
      }
    }
    // set local molecule order
    MoleculeSet::set_local_molecules_to_atoms();
    needReset_ = false;
  }

  //--------------------------------------------------------
  //  atoms_by_global_molecule
  //--------------------------------------------------------
  set<int> SmallMoleculeSet::atoms_by_global_molecule(int id) const
  {
    // take all atoms and prune out ghosts
    set<int> realAtoms = MoleculeSet::atoms_by_global_molecule(id);
    remove_proc_ghosts(realAtoms);
    return realAtoms;
  }

  //--------------------------------------------------------
  //  local_fraction
  //--------------------------------------------------------
  double SmallMoleculeSet::local_fraction(int id) const
  {
    if (need_reset()) reset();
    set<int> realAtoms = (localMoleculeToAtoms_[id])->second;
    int totalAtoms = realAtoms.size();
    remove_proc_ghosts(realAtoms);
    return double(realAtoms.size())/double(totalAtoms);
  }

  //--------------------------------------------------------
  //  remove_proc_ghosts
  //--------------------------------------------------------
  void SmallMoleculeSet::remove_proc_ghosts(set<int> & atomSet) const
  {
    int nlocalIdx = lammps_->nlocal() - 1;
    set<int>::const_iterator atomIt;
    for (atomIt = atomSet.begin(); atomIt != atomSet.end(); atomIt++) {
      if (*atomIt > nlocalIdx)
        atomSet.erase(atomIt);
    }
  }


};
