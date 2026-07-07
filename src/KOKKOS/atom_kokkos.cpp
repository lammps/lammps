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
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "update.h"
#include "kokkos_base.h"
#include "modify.h"
#include "fix.h"
#include "fix_property_atom_kokkos.h"

#include <map>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomKokkos::AtomKokkos(LAMMPS *lmp) : Atom(lmp)
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
  hybrid_flag = 0;
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
  memoryKK->destroy_kokkos(k_mu, mu);

  memoryKK->destroy_kokkos(k_radius, radius);
  memoryKK->destroy_kokkos(k_rmass, rmass);
  memoryKK->destroy_kokkos(k_omega, omega);
  memoryKK->destroy_kokkos(k_angmom, angmom);
  memoryKK->destroy_kokkos(k_torque, torque);
  memoryKK->destroy_kokkos(k_ellipsoid, ellipsoid);

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

  // ivector/dvector are single contiguous Kokkos views, with the legacy
  // ivector[i]/dvector[i] pointers aliasing into them.  Free the view data and
  // null those aliases, but leave the row-pointer arrays themselves in place:
  // the base Atom destructor then safely memory->destroy()s the (null) aliases
  // and sfree()s the row-pointer arrays via the standard nullptr-safe path, so
  // no Kokkos-specific handling is needed in ~Atom.

  k_dvector = DAT::ttransform_kkfloat_2d();
  for (int i = 0; i < ndvector; i++) dvector[i] = nullptr;

  k_ivector = DAT::tdual_int_2d_lr();
  for (int i = 0; i < nivector; i++) ivector[i] = nullptr;

  // views-of-views: destroy each inner DualView (this frees its row-pointer
  // array and nulls the legacy iarray[i]/darray[i] alias, so the base Atom
  // destructor's nullptr-safe memory->destroy() is a no-op for them), then
  // release the outer DualView, which frees the inner views' data.

  for (int i = 0; i < niarray; i++)
    memoryKK->destroy_kokkos(k_iarray.view_host()[i].k_view, iarray[i]);
  k_iarray = tdual_struct_tdual_int_2d_1d();

  for (int i = 0; i < ndarray; i++)
    memoryKK->destroy_kokkos(k_darray.view_host()[i].k_view, darray[i]);
  k_darray = tdual_struct_tdual_double_2d_1d();

  delete [] fix_prop_atom;
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::init()
{
  Atom::init();

  sort_legacy = lmp->kokkos->sort_legacy;
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

/* ----------------------------------------------------------------------
   return a pointer to a per-atom property by name (used by the library
   interface).  Override the base-class version so that data which lives on
   the device is first synced back to the host.  Without this, library or
   Python calls to extract_atom() that are not aligned with an output step
   (e.g. issued from the LAMMPS GUI or a "python" command during a run) may
   hand out stale host data when running with KOKKOS on a GPU.  See issue #3945.
------------------------------------------------------------------------- */

void *AtomKokkos::extract(const char *name)
{
  // map the public extract name to the KOKKOS data mask of the dual view that
  // holds it.  Names whose data is not device-resident are simply absent here
  // (their host copy is always current) and fall through to Atom::extract().

  static const std::map<std::string, uint64_t> extract_mask = {
      {"id", TAG_MASK}, {"type", TYPE_MASK}, {"mask", MASK_MASK}, {"image", IMAGE_MASK},
      {"x", X_MASK}, {"v", V_MASK}, {"f", F_MASK}, {"q", Q_MASK}, {"mu", MU_MASK},
      {"omega", OMEGA_MASK}, {"angmom", ANGMOM_MASK}, {"torque", TORQUE_MASK},
      {"radius", RADIUS_MASK}, {"rmass", RMASS_MASK}, {"ellipsoid", ELLIPSOID_MASK},
      {"molecule", MOLECULE_MASK}, {"nspecial", SPECIAL_MASK}, {"special", SPECIAL_MASK},
      {"num_bond", BOND_MASK}, {"bond_type", BOND_MASK}, {"bond_atom", BOND_MASK},
      {"num_angle", ANGLE_MASK}, {"angle_type", ANGLE_MASK},
      {"angle_atom1", ANGLE_MASK}, {"angle_atom2", ANGLE_MASK}, {"angle_atom3", ANGLE_MASK},
      {"num_dihedral", DIHEDRAL_MASK}, {"dihedral_type", DIHEDRAL_MASK},
      {"dihedral_atom1", DIHEDRAL_MASK}, {"dihedral_atom2", DIHEDRAL_MASK},
      {"dihedral_atom3", DIHEDRAL_MASK}, {"dihedral_atom4", DIHEDRAL_MASK},
      {"num_improper", IMPROPER_MASK}, {"improper_type", IMPROPER_MASK},
      {"improper_atom1", IMPROPER_MASK}, {"improper_atom2", IMPROPER_MASK},
      {"improper_atom3", IMPROPER_MASK}, {"improper_atom4", IMPROPER_MASK},
      {"sp", SP_MASK}, {"dpdTheta", DPDTHETA_MASK}};

  const auto it = extract_mask.find(name);
  if (it != extract_mask.end()) {
    sync(Host, it->second);
  } else if (utils::strmatch(name, "^[id]2?_")) {
    // custom per-atom data (fix property/atom). each prefix maps to its own
    // data mask:
    //   i_  -> ivector (IVECTOR_MASK)    d_  -> dvector (DVECTOR_MASK)
    //   i2_ -> iarray  (IARRAY_MASK)     d2_ -> darray  (DARRAY_MASK)
    // all four custom data types are device-resident, so each sync pulls the
    // requested property back to the host before extract() returns it.
    const bool dbl = (name[0] == 'd');
    const bool arr = (name[1] == '2');
    uint64_t cmask;
    if (!dbl && !arr) cmask = IVECTOR_MASK;
    else if (dbl && !arr) cmask = DVECTOR_MASK;
    else if (!dbl && arr) cmask = IARRAY_MASK;
    else cmask = DARRAY_MASK;
    sync(Host, cmask);
  }

  return Atom::extract(name);
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sync(const ExecutionSpace space, uint64_t mask)
{
  if ((space == Device || space == HostKK) && lmp->kokkos->auto_sync) {

    // sync HostKK -> Host if needed

    avecKK->sync(Host, mask);
    for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->sync(Host, mask);

    avecKK->modified(Host, mask);
    for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->modified(Host, mask);
  }

  avecKK->sync(space, mask);
  for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->sync(space, mask);
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::modified(const ExecutionSpace space, uint64_t mask)
{
  avecKK->modified(space, mask);
  for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->modified(space, mask);

  if ((space == Device || space == HostKK) && lmp->kokkos->auto_sync) {
    avecKK->sync(Host, mask);
    for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->sync(Host, mask);
  }
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sync_pinned(const ExecutionSpace space, uint64_t mask, int async_flag)
{
  avecKK->sync_pinned(space, mask, async_flag);
  for (int n = 0; n < nprop_atom; n++) fix_prop_atom[n]->sync_pinned(space, mask, async_flag);
}
/* ---------------------------------------------------------------------- */

void AtomKokkos::allocate_type_arrays()
{
  if (avec->mass_type == AtomVec::PER_TYPE) {
    memoryKK->create_kokkos(k_mass,mass,ntypes + 1,"atom::mass");
    mass_setflag = new int[ntypes + 1];
    for (int itype = 1; itype <= ntypes; itype++) mass_setflag[itype] = 0;
    k_mass.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sort()
{
  // check if all fixes with atom-based arrays support sort on device

  if (!sort_legacy) {
    int flag = 1;
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++) {
      auto fix_iextra = modify->fix[atom->extra_grow[iextra]];
      if (!fix_iextra->sort_device) {
        flag = 0;
        if (comm->me == 0)
          error->warning(FLERR,"Fix {} not (yet) compatible with Kokkos sorting on device", fix_iextra->style);
        break;
      }
    }
    if (!flag) {
      if (comm->me == 0) {
        error->warning(FLERR,"Fix with atom-based arrays not (yet) compatible with Kokkos sorting on device, "
                           "switching to legacy host sorting");
      }
      sort_legacy = true;
    }

    int bonus_flag = (ellipsoid_flag || line_flag || tri_flag || body_flag);

    if (bonus_flag) {
      if (comm->me == 0) {
        error->warning(FLERR,"Atom bonus data not (yet) compatible with Kokkos sorting on device, "
                           "switching to legacy host sorting");
      }
      sort_legacy = true;
    }

    if (hybrid_flag) {
      if (comm->me == 0) {
        error->warning(FLERR,"Atom_style hybrid not (yet) compatible with Kokkos sorting on device, "
                           "switching to legacy host sorting");
      }
      sort_legacy = true;
    }
  }

  if (sort_legacy) {
    sync(Host, ALL_MASK);
    int prev_auto_sync = lmp->kokkos->auto_sync;
    lmp->kokkos->auto_sync = 1;
    Atom::sort();
    lmp->kokkos->auto_sync = prev_auto_sync;
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

  auto d_x = k_x.view_device();
  sync(Device, X_MASK);

  // sort

  int max_bins[3];
  max_bins[0] = nbinx;
  max_bins[1] = nbiny;
  max_bins[2] = nbinz;

  using KeyViewType = DAT::t_kkfloat_1d_3_lr;
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
    this->sync(Device, IVECTOR_MASK);
    memoryKK->grow_kokkos(k_ivector, ivector, nivector, nmax, "atom:ivector");
    this->modified(Device, IVECTOR_MASK);

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
    iarray[index] = nullptr;

    icols = (int *) memory->srealloc(icols, niarray * sizeof(int), "atom:icols");
    icols[index] = cols;

    // grow the outer view-of-views by one inner DualView for this property.
    // SequentialHostInit is required: the struct elements wrap a DualView
    // (non-trivial ctor + atomic refcount), so the default parallel host init
    // would race constructing/copying the inner DualViews.
    k_iarray.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit), niarray);
    memoryKK->create_kokkos(k_iarray.view_host()[index].k_view, iarray[index],
                            nmax, cols, "atom:iarray");
    k_iarray.modify_host();
    k_iarray.sync_device();

  } else if (flag == 1 && cols) {
    index = ndarray;
    ndarray++;
    daname = (char **) memory->srealloc(daname, ndarray * sizeof(char *), "atom:daname");
    daname[index] = utils::strdup(name);
    daghost = (int *) memory->srealloc(daghost, ndarray * sizeof(int), "atom:daghost");
    daghost[index] = ghost;
    darray = (double ***) memory->srealloc(darray, ndarray * sizeof(double **), "atom:darray");
    darray[index] = nullptr;

    dcols = (int *) memory->srealloc(dcols, ndarray * sizeof(int), "atom:dcols");
    dcols[index] = cols;

    // see iarray branch above for why SequentialHostInit is required here
    k_darray.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit), ndarray);
    memoryKK->create_kokkos(k_darray.view_host()[index].k_view, darray[index],
                            nmax, cols, "atom:darray");
    k_darray.modify_host();
    k_darray.sync_device();
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
  // the per-atom data is Kokkos-managed (k_ivector/k_dvector are contiguous, and
  // k_iarray/k_darray are views-of-views), so do NOT memory->destroy() it here --
  // that would free pointers that alias into a Kokkos view.  For ivector/dvector
  // the data lives in a shared contiguous view and cannot be freed per index, so
  // just drop the legacy alias; for iarray/darray we free both the row-pointer
  // array and the inner DualView's host+device storage.

  if (flag == 0 && cols == 0) {
    ivector[index] = nullptr;
    delete[] ivname[index];
    ivname[index] = nullptr;

  } else if (flag == 1 && cols == 0) {
    dvector[index] = nullptr;
    delete[] dvname[index];
    dvname[index] = nullptr;

  } else if (flag == 0 && cols) {
    // destroy_kokkos receives the inner view by value, so it only frees the
    // row-pointer array; reset the stored inner DualView to release its
    // host+device data, then push the emptied slot to the device-side view
    memoryKK->destroy_kokkos(k_iarray.view_host()[index].k_view, iarray[index]);
    k_iarray.view_host()[index].k_view = DAT::tdual_int_2d_lr();
    k_iarray.modify_host();
    k_iarray.sync_device();
    iarray[index] = nullptr;
    delete[] ianame[index];
    ianame[index] = nullptr;

  } else if (flag == 1 && cols) {
    memoryKK->destroy_kokkos(k_darray.view_host()[index].k_view, darray[index]);
    k_darray.view_host()[index].k_view = DAT::tdual_double_2d_lr();
    k_darray.modify_host();
    k_darray.sync_device();
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

  hybrid_flag = (avec != nullptr);

  AtomVec *avec = Atom::new_avec(style, trysuffix, sflag);
  if (!avec->kokkosable) error->all(FLERR, "KOKKOS package requires a Kokkos-enabled atom_style");

  if (!hybrid_flag)
    avecKK = dynamic_cast<AtomVecKokkos*>(avec);

  return avec;
}
