// A class for managing the data associated with sets of molecules

#ifndef PER_MOLECULE_SET_H
#define PER_MOLECULE_SET_H

// ATC_Method headers
#include "LammpsInterface.h"
#include "DependencyManager.h"
#include <map>
#include <set>
#include <vector>

namespace ATC {

  // forward declarations
  class ATC_Method;
  template <typename T> class PerAtomQuantity;

  /**
   *  @class  MoleculeSet
   *  @brief  A class for handling all the data associated with sets of molecules
   */

  class MoleculeSet : public DependencyManager {

  public:

    MoleculeSet(ATC_Method * atc, int groupBit);

    virtual ~MoleculeSet();

    /** reset all data */
    virtual void clear();

    /** initialize global data */
    virtual void initialize(std::map<int, double> * globalAtomsPerMolecule = NULL);

    /** reset the number of atoms/molecules on this processor */
    void reset_nlocal() {this->set_reset();};

    /** recompute data when atoms cross processors */
    void post_exchange() {this->set_reset();};

    /** access the number of total molecules */
    int global_molecule_count() const {return nMoleculesTotal_;};

    /** access the number of local molecules */
    int local_molecule_count() const {if (need_reset()) reset(); return moleculeToAtoms_.size();};

    /** access molecule atoms by lammps id */
    std::set<int> atoms_by_global_molecule(int id) const;

    /** access molecules by local indexing */
    const std::set<int> & atoms_by_local_molecule(int id) const;

    /** access fraction of a locally indexed molecule on this processor */
    virtual double local_fraction(int id) const = 0;

    /** use global index to get local index */
    //int global_to_local(int id) const;

    /** use local index to get global index */
    int local_to_global(int id) const {return (*localMoleculeToAtoms_[id]).first;};

  protected:

    /** pointer for access to atc data */
    ATC_Method * atc_;

    /** lammps group bit corresponding to desired molecules */
    int groupBit_;

    /** pointer to lammps interface */
    const LammpsInterface * lammps_;

    /** total number of molecules in this group */
    // see Compute::molecules_in_group
    int nMoleculesTotal_;

    /** multimap from lammps molecule id to ids of consituent atoms, all atoms are real */
    // multiple map to account for periodic images
    mutable std::multimap<int, std::set<int> > moleculeToAtoms_;

    /** vector in processor-local molecule order to constituent atom sets, atoms include ghosts */
    mutable std::vector< std::map<int, std::set<int> >::const_iterator > localMoleculeToAtoms_;

    /** resets the quantity based on the latest data */
    virtual void reset() const = 0;

    /** creates the ordered list of local molecules */
    void set_local_molecules_to_atoms() const;

  private:

    // do not define this
    MoleculeSet();

  };


  /**
   *  @class  SmallMoleculeSet
   *  @brief  A class for handling data for small molecules, i.e., molecules with maximum distance between atoms less than the lammps cutoff radius.  Atom ids are in [0,nlocalTotal-1].
   */

  class SmallMoleculeSet : public MoleculeSet {

  public:

    SmallMoleculeSet(ATC_Method * atc, int groupBit,
                     PerAtomQuantity<int> * bondList = NULL,
                     PerAtomQuantity<int> * numBond = NULL);

    virtual ~SmallMoleculeSet();
    
    /** reset all data */
    virtual void clear();

    /** initialize global data */
    virtual void initialize(std::map<int, double> * globalAtomsPerMolecule = NULL);

    /** access molecule atoms by lammps id */
    std::set<int> atoms_by_global_molecule(int id) const;

    /** access fraction of a locally indexed molecule on this processor */
    virtual double local_fraction(int id) const;

  protected:

    /** store the number of atoms in a molecule on this processor */
    //std::map<int, int> localAtomsPerMolecule_;

    /** resets the quantity based on the latest data */
    virtual void reset() const;

    /** data structure containing bond list information, forces parallel communication of bond lists */
    PerAtomQuantity<int> * bondList_;

    /** data structure containing bond list information, forces parallel communication of bond lists */
    PerAtomQuantity<int> * numBond_;

    /** removes processor ghosts from a set of atom ids */
    void remove_proc_ghosts(std::set<int> & atomSet) const;

    // workspace variable for determining if we've hit an internal atom already
    mutable Array<bool> _atomFound_;

  private:

    // do not define this
    SmallMoleculeSet();

  };


};

#endif
