/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_kokkos.h"

#include "atom_masks.h"
#include "atom_vec.h"
#include "atom_vec_kokkos.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "error.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "update.h"
#include "kokkos_base.h"
#include "modify.h"
#include "fix.h"
#include "fix_property_atom_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomKokkos::AtomKokkos(LAMMPS *lmp) : Atom(lmp),
mapBinner(1, 0.0, 1.0), // no default constructor, these values are not used
mapSorter(d_tag_sorted, 0, 1, mapBinner, true)
{
  avecKK = nullptr;

  k_error_flag = DAT::tdual_int_scalar("atom:error_flag");

  d_tag_min_max = t_tagint_2(Kokkos::NoInit("atom:tag_min_max"));
  h_tag_min_max = t_host_tagint_2(Kokkos::NoInit("atom:tag_min_max"));

  d_tag_min = Kokkos::subview(d_tag_min_max,0);
  d_tag_max = Kokkos::subview(d_tag_min_max,1);

  h_tag_min = Kokkos::subview(h_tag_min_max,0);
  h_tag_max = Kokkos::subview(h_tag_min_max,1);

  nprop_atom = 0;
  fix_prop_atom = nullptr;
}

/* ---------------------------------------------------------------------- */

AtomKokkos::~AtomKokkos()
{
  memoryKK->destroy_kokkos(k_tag, tag);
  memoryKK->destroy_kokkos(k_mask, mask);
  memoryKK->destroy_kokkos(k_type, type);
  memoryKK->destroy_kokkos(k_image, image);
  memoryKK->destroy_kokkos(k_molecule, molecule);

  memoryKK->destroy_kokkos(k_x, x);
  memoryKK->destroy_kokkos(k_v, v);
  memoryKK->destroy_kokkos(k_f, f);

  memoryKK->destroy_kokkos(k_mass, mass);
  memoryKK->destroy_kokkos(k_q, q);

  memoryKK->destroy_kokkos(k_radius, radius);
  memoryKK->destroy_kokkos(k_rmass, rmass);
  memoryKK->destroy_kokkos(k_omega, omega);
  memoryKK->destroy_kokkos(k_angmom, angmom);
  memoryKK->destroy_kokkos(k_torque, torque);

  memoryKK->destroy_kokkos(k_nspecial, nspecial);
  memoryKK->destroy_kokkos(k_special, special);
  memoryKK->destroy_kokkos(k_num_bond, num_bond);
  memoryKK->destroy_kokkos(k_bond_type, bond_type);
  memoryKK->destroy_kokkos(k_bond_atom, bond_atom);
  memoryKK->destroy_kokkos(k_num_angle, num_angle);
  memoryKK->destroy_kokkos(k_angle_type, angle_type);
  memoryKK->destroy_kokkos(k_angle_atom1, angle_atom1);
  memoryKK->destroy_kokkos(k_angle_atom2, angle_atom2);
  memoryKK->destroy_kokkos(k_angle_atom3, angle_atom3);
  memoryKK->destroy_kokkos(k_num_dihedral, num_dihedral);
  memoryKK->destroy_kokkos(k_dihedral_type, dihedral_type);
  memoryKK->destroy_kokkos(k_dihedral_atom1, dihedral_atom1);
  memoryKK->destroy_kokkos(k_dihedral_atom2, dihedral_atom2);
  memoryKK->destroy_kokkos(k_dihedral_atom3, dihedral_atom3);
  memoryKK->destroy_kokkos(k_dihedral_atom4, dihedral_atom4);
  memoryKK->destroy_kokkos(k_num_improper, num_improper);
  memoryKK->destroy_kokkos(k_improper_type, improper_type);
  memoryKK->destroy_kokkos(k_improper_atom1, improper_atom1);
  memoryKK->destroy_kokkos(k_improper_atom2, improper_atom2);
  memoryKK->destroy_kokkos(k_improper_atom3, improper_atom3);
  memoryKK->destroy_kokkos(k_improper_atom4, improper_atom4);

  AtomKokkos::map_delete();

  // SPIN package

  memoryKK->destroy_kokkos(k_sp, sp);
  memoryKK->destroy_kokkos(k_fm, fm);
  memoryKK->destroy_kokkos(k_fm_long, fm_long);

  // DPD-REACT package
  memoryKK->destroy_kokkos(k_uCond, uCond);
  memoryKK->destroy_kokkos(k_uMech, uMech);
  memoryKK->destroy_kokkos(k_uChem, uChem);
  memoryKK->destroy_kokkos(k_uCG, uCG);
  memoryKK->destroy_kokkos(k_uCGnew, uCGnew);
  memoryKK->destroy_kokkos(k_rho, rho);
  memoryKK->destroy_kokkos(k_dpdTheta, dpdTheta);
  memoryKK->destroy_kokkos(k_duChem, duChem);

  memoryKK->destroy_kokkos(k_dvector, dvector);
  dvector = nullptr;
  delete [] fix_prop_atom;
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::init()
{
  Atom::init();

  sort_classic = lmp->kokkos->sort_classic;
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::update_property_atom()
{
  nprop_atom = 0;
  std::vector<Fix *> prop_atom_fixes;
  for (auto &ifix : modify->get_fix_by_style("^property/atom")) {
    if (!ifix->kokkosable)
      error->all(FLERR, "KOKKOS package requires a Kokkos-enabled version of fix property/atom");

    ++nprop_atom;
    prop_atom_fixes.push_back(ifix);
  }

  delete[] fix_prop_atom;
  fix_prop_atom = new FixPropertyAtomKokkos *[nprop_atom];

  int n = 0;
  for (auto &ifix : prop_atom_fixes)
    fix_prop_atom[n++] = dynamic_cast<FixPropertyAtomKokkos *>(ifix);
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sync(const ExecutionSpace space, unsigned int mask)
{
  if (space == Device && lmp->kokkos->auto_sync) {
    avecKK->modified(Host, mask);
    for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->modified(Host, mask);
  }

  avecKK->sync(space, mask);
  for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->sync(space, mask);
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::modified(const ExecutionSpace space, unsigned int mask)
{
  avecKK->modified(space, mask);
  for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->modified(space, mask);

  if (space == Device && lmp->kokkos->auto_sync) {
    avecKK->sync(Host, mask);
    for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->sync(Host, mask);
  }
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sync_overlapping_device(const ExecutionSpace space, unsigned int mask)
{
  avecKK->sync_overlapping_device(space, mask);
  for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->sync_overlapping_device(space, mask);
}
/* ---------------------------------------------------------------------- */

void AtomKokkos::allocate_type_arrays()
{
  if (avec->mass_type == AtomVec::PER_TYPE) {
    k_mass = DAT::tdual_float_1d("Mass", ntypes + 1);
    mass = k_mass.h_view.data();
    mass_setflag = new int[ntypes + 1];
    for (int itype = 1; itype <= ntypes; itype++) mass_setflag[itype] = 0;
    k_mass.modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sort()
{
  // check if all fixes with atom-based arrays support sort on device

  if (!sort_classic) {
    int flag = 1;
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) {
      auto fix_iextra = modify->fix[atom->extra_grow[iextra]];
      if (!fix_iextra->sort_device) {
        flag = 0;
        break;
      }
    }
    if (!flag) {
      if (comm->me == 0) {
        error->warning(FLERR,"Fix with atom-based arrays not compatible with Kokkos sorting on device, "
                           "switching to classic host sorting");
      }
      sort_classic = true;
    }
  }

  if (sort_classic) {
    sync(Host, ALL_MASK);
    Atom::sort();
    modified(Host, ALL_MASK);
  } else sort_device();
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sort_device()
{
  // set next timestep for sorting to take place

  nextsort = (update->ntimestep / sortfreq) * sortfreq + sortfreq;

  // re-setup sort bins if needed

  if (domain->box_change) setup_sort_bins();
  if (nbins == 1) return;

  // for triclinic, atoms must be in box coords (not lamda) to match bbox

  if (domain->triclinic) domain->lamda2x(nlocal);

  auto d_x = k_x.d_view;
  sync(Device, X_MASK);

  // sort

  int max_bins[3];
  max_bins[0] = nbinx;
  max_bins[1] = nbiny;
  max_bins[2] = nbinz;

  using KeyViewType = DAT::t_x_array;
  using BinOp = BinOp3DLAMMPS<KeyViewType>;
  BinOp binner(max_bins, bboxlo, bboxhi);
  Kokkos::BinSort<KeyViewType, BinOp> Sorter(d_x, 0, nlocal, binner, false);
  Sorter.create_permute_vector(LMPDeviceType());

  avecKK->sort_kokkos(Sorter);

  if (atom->nextra_grow) {
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) {
      auto fix_iextra = modify->fix[atom->extra_grow[iextra]];
      KokkosBase *kkbase = dynamic_cast<KokkosBase*>(fix_iextra);

      kkbase->sort_kokkos(Sorter);
    }
  }

 //  convert back to lamda coords

 if (domain->triclinic) domain->x2lamda(nlocal);
}

/* ----------------------------------------------------------------------
   reallocate memory to the pointer selected by the mask
------------------------------------------------------------------------- */

void AtomKokkos::grow(unsigned int mask)
{
  if (mask & SPECIAL_MASK) {
    memoryKK->destroy_kokkos(k_special, special);
    sync(Device, mask);
    modified(Device, mask);
    memoryKK->grow_kokkos(k_special, special, nmax, maxspecial, "atom:special");
    avec->grow_pointers();
    sync(Host, mask);
  }
}

/* ----------------------------------------------------------------------
   add a custom variable with name of type flag = 0/1 for int/double
   assumes name does not already exist
   return index in ivector or dvector of its location
------------------------------------------------------------------------- */

int AtomKokkos::add_custom(const char *name, int flag, int cols, int ghost)
{
  int index = -1;

  if (flag == 0 && cols == 0) {
    index = nivector;
    nivector++;
    ivname = (char **) memory->srealloc(ivname, nivector * sizeof(char *), "atom:ivname");
    ivname[index] = utils::strdup(name);
    ivghost = (int *) memory->srealloc(ivghost,nivector * sizeof(int),"atom:ivghost");
    ivghost[index] = ghost;
    ivector = (int **) memory->srealloc(ivector, nivector * sizeof(int *), "atom:ivector");
    memory->create(ivector[index], nmax, "atom:ivector");

  } else if (flag == 1 && cols == 0) {
    index = ndvector;
    ndvector++;
    dvname = (char **) memory->srealloc(dvname, ndvector * sizeof(char *), "atom:dvname");
    dvname[index] = utils::strdup(name);
    dvghost = (int *) memory->srealloc(dvghost, ndvector * sizeof(int), "atom:dvghost");
    dvghost[index] = ghost;
    dvector = (double **) memory->srealloc(dvector, ndvector * sizeof(double *), "atom:dvector");
    this->sync(Device, DVECTOR_MASK);
    memoryKK->grow_kokkos(k_dvector, dvector, ndvector, nmax, "atom:dvector");
    this->modified(Device, DVECTOR_MASK);

  } else if (flag == 0 && cols) {
    index = niarray;
    niarray++;
    ianame = (char **) memory->srealloc(ianame, niarray * sizeof(char *), "atom:ianame");
    ianame[index] = utils::strdup(name);
    iaghost = (int *) memory->srealloc(iaghost, niarray * sizeof(int), "atom:iaghost");
    iaghost[index] = ghost;
    iarray = (int ***) memory->srealloc(iarray, niarray * sizeof(int **), "atom:iarray");
    memory->create(iarray[index], nmax, cols, "atom:iarray");

    icols = (int *) memory->srealloc(icols, niarray * sizeof(int), "atom:icols");
    icols[index] = cols;

  } else if (flag == 1 && cols) {
    index = ndarray;
    ndarray++;
    daname = (char **) memory->srealloc(daname, ndarray * sizeof(char *), "atom:daname");
    daname[index] = utils::strdup(name);
    daghost = (int *) memory->srealloc(daghost, ndarray * sizeof(int), "atom:daghost");
    daghost[index] = ghost;
    darray = (double ***) memory->srealloc(darray, ndarray * sizeof(double **), "atom:darray");
    memory->create(darray[index], nmax, cols, "atom:darray");

    dcols = (int *) memory->srealloc(dcols, ndarray * sizeof(int), "atom:dcols");
    dcols[index] = cols;
  }

  if (index < 0)
    error->all(FLERR,"Invalid call to AtomKokkos::add_custom()");

  return index;
}

/* ----------------------------------------------------------------------
   remove a custom variable of type flag = 0/1 for int/double at index
   free memory for vector/array and name and set ptrs to a null pointer
   these lists never shrink
------------------------------------------------------------------------- */

void AtomKokkos::remove_custom(int index, int flag, int cols)
{
  if (flag == 0 && cols == 0) {
    memory->destroy(ivector[index]);
    ivector[index] = nullptr;
    delete[] ivname[index];
    ivname[index] = nullptr;

  } else if (flag == 1 && cols == 0) {
    dvector[index] = nullptr;
    delete[] dvname[index];
    dvname[index] = nullptr;

  } else if (flag == 0 && cols) {
    memory->destroy(iarray[index]);
    iarray[index] = nullptr;
    delete[] ianame[index];
    ianame[index] = nullptr;

  } else if (flag == 1 && cols) {
    memory->destroy(darray[index]);
    darray[index] = nullptr;
    delete[] daname[index];
    daname[index] = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::deallocate_topology()
{
  memoryKK->destroy_kokkos(k_bond_type, bond_type);
  memoryKK->destroy_kokkos(k_bond_atom, bond_atom);

  memoryKK->destroy_kokkos(k_angle_type, angle_type);
  memoryKK->destroy_kokkos(k_angle_atom1, angle_atom1);
  memoryKK->destroy_kokkos(k_angle_atom2, angle_atom2);
  memoryKK->destroy_kokkos(k_angle_atom3, angle_atom3);

  memoryKK->destroy_kokkos(k_dihedral_type, dihedral_type);
  memoryKK->destroy_kokkos(k_dihedral_atom1, dihedral_atom1);
  memoryKK->destroy_kokkos(k_dihedral_atom2, dihedral_atom2);
  memoryKK->destroy_kokkos(k_dihedral_atom3, dihedral_atom3);
  memoryKK->destroy_kokkos(k_dihedral_atom4, dihedral_atom4);

  memoryKK->destroy_kokkos(k_improper_type, improper_type);
  memoryKK->destroy_kokkos(k_improper_atom1, improper_atom1);
  memoryKK->destroy_kokkos(k_improper_atom2, improper_atom2);
  memoryKK->destroy_kokkos(k_improper_atom3, improper_atom3);
  memoryKK->destroy_kokkos(k_improper_atom4, improper_atom4);
}

/* ---------------------------------------------------------------------- */

AtomVec *AtomKokkos::new_avec(const std::string &style, int trysuffix, int &sflag)
{
  // check if avec already exists, if so this is a hybrid substyle

  int hybrid_substyle_flag = (avec != nullptr);

  AtomVec *avec = Atom::new_avec(style, trysuffix, sflag);
  if (!avec->kokkosable) error->all(FLERR, "KOKKOS package requires a Kokkos-enabled atom_style");

  if (!hybrid_substyle_flag)
    avecKK = dynamic_cast<AtomVecKokkos*>(avec);

  return avec;
}
