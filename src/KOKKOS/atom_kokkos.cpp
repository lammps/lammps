/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomKokkos::AtomKokkos(LAMMPS *lmp) : Atom(lmp)
{
  k_error_flag = DAT::tdual_int_scalar("atom:error_flag");
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
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sync(const ExecutionSpace space, unsigned int mask)
{
  if (space == Device && lmp->kokkos->auto_sync) ((AtomVecKokkos *) avec)->modified(Host, mask);

  ((AtomVecKokkos *) avec)->sync(space, mask);
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::modified(const ExecutionSpace space, unsigned int mask)
{
  ((AtomVecKokkos *) avec)->modified(space, mask);

  if (space == Device && lmp->kokkos->auto_sync) ((AtomVecKokkos *) avec)->sync(Host, mask);
}

void AtomKokkos::sync_overlapping_device(const ExecutionSpace space, unsigned int mask)
{
  ((AtomVecKokkos *) avec)->sync_overlapping_device(space, mask);
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
  int i, m, n, ix, iy, iz, ibin, empty;

  // set next timestep for sorting to take place

  nextsort = (update->ntimestep / sortfreq) * sortfreq + sortfreq;

  // re-setup sort bins if needed

  if (domain->box_change) setup_sort_bins();
  if (nbins == 1) return;

  // reallocate per-atom vectors if needed

  if (atom->nmax > maxnext) {
    memory->destroy(next);
    memory->destroy(permute);
    maxnext = atom->nmax;
    memory->create(next, maxnext, "atom:next");
    memory->create(permute, maxnext, "atom:permute");
  }

  // insure there is one extra atom location at end of arrays for swaps

  if (nlocal == nmax) avec->grow(0);

  sync(Host, ALL_MASK);
  modified(Host, ALL_MASK);

  // bin atoms in reverse order so linked list will be in forward order

  for (i = 0; i < nbins; i++) binhead[i] = -1;

  HAT::t_x_array_const h_x = k_x.view<LMPHostType>();
  for (i = nlocal - 1; i >= 0; i--) {
    ix = static_cast<int>((h_x(i, 0) - bboxlo[0]) * bininvx);
    iy = static_cast<int>((h_x(i, 1) - bboxlo[1]) * bininvy);
    iz = static_cast<int>((h_x(i, 2) - bboxlo[2]) * bininvz);
    ix = MAX(ix, 0);
    iy = MAX(iy, 0);
    iz = MAX(iz, 0);
    ix = MIN(ix, nbinx - 1);
    iy = MIN(iy, nbiny - 1);
    iz = MIN(iz, nbinz - 1);
    ibin = iz * nbiny * nbinx + iy * nbinx + ix;
    next[i] = binhead[ibin];
    binhead[ibin] = i;
  }

  // permute = desired permutation of atoms
  // permute[I] = J means Ith new atom will be Jth old atom

  n = 0;
  for (m = 0; m < nbins; m++) {
    i = binhead[m];
    while (i >= 0) {
      permute[n++] = i;
      i = next[i];
    }
  }

  // current = current permutation, just reuse next vector
  // current[I] = J means Ith current atom is Jth old atom

  int *current = next;
  for (i = 0; i < nlocal; i++) current[i] = i;

  // reorder local atom list, when done, current = permute
  // perform "in place" using copy() to extra atom location at end of list
  // inner while loop processes one cycle of the permutation
  // copy before inner-loop moves an atom to end of atom list
  // copy after inner-loop moves atom at end of list back into list
  // empty = location in atom list that is currently empty

  for (i = 0; i < nlocal; i++) {
    if (current[i] == permute[i]) continue;
    avec->copy(i, nlocal, 0);
    empty = i;
    while (permute[empty] != i) {
      avec->copy(permute[empty], empty, 0);
      empty = current[empty] = permute[empty];
    }
    avec->copy(nlocal, empty, 0);
    current[empty] = permute[empty];
  }

  // sanity check that current = permute

  //int flag = 0;
  //for (i = 0; i < nlocal; i++)
  //  if (current[i] != permute[i]) flag = 1;
  //int flagall;
  //MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  //if (flagall) errorX->all(FLERR,"Atom sort did not operate correctly");
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

int AtomKokkos::add_custom(const char *name, int flag, int cols)
{
  int index;

  if (flag == 0 && cols == 0) {
    index = nivector;
    nivector++;
    ivname = (char **) memory->srealloc(ivname,nivector*sizeof(char *),
					"atom:ivname");
    int n = strlen(name) + 1;
    ivname[index] = new char[n];
    strcpy(ivname[index],name);
    ivector = (int **) memory->srealloc(ivector,nivector*sizeof(int *),
                                        "atom:ivector");
    memory->create(ivector[index],nmax,"atom:ivector");

  } else if (flag == 1 && cols == 0) {
    index = ndvector;
    ndvector++;
    dvname = (char **) memory->srealloc(dvname,ndvector*sizeof(char *),
					"atom:dvname");
    int n = strlen(name) + 1;
    dvname[index] = new char[n];
    strcpy(dvname[index],name);
    dvector = (double **) memory->srealloc(dvector,ndvector*sizeof(double *),
                                           "atom:dvector");
    this->sync(Device,DVECTOR_MASK);
    memoryKK->grow_kokkos(k_dvector,dvector,ndvector,nmax,
                        "atom:dvector");
    this->modified(Device,DVECTOR_MASK);

  } else if (flag == 0 && cols) {
    index = niarray;
    niarray++;
    ianame = (char **) memory->srealloc(ianame,niarray*sizeof(char *),
					"atom:ianame");
    int n = strlen(name) + 1;
    ianame[index] = new char[n];
    strcpy(ianame[index],name);
    iarray = (int ***) memory->srealloc(iarray,niarray*sizeof(int **),
					"atom:iarray");
    memory->create(iarray[index],nmax,cols,"atom:iarray");

    icols = (int *) memory->srealloc(icols,niarray*sizeof(int),"atom:icols");
    icols[index] = cols;

  } else if (flag == 1 && cols) {
    index = ndarray;
    ndarray++;
    daname = (char **) memory->srealloc(daname,ndarray*sizeof(char *),
					"atom:daname");
    int n = strlen(name) + 1;
    daname[index] = new char[n];
    strcpy(daname[index],name);
    darray = (double ***) memory->srealloc(darray,ndarray*sizeof(double **),
					   "atom:darray");
    memory->create(darray[index],nmax,cols,"atom:darray");

    dcols = (int *) memory->srealloc(dcols,ndarray*sizeof(int),"atom:dcols");
    dcols[index] = cols;
  }

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
    ivector[index] = NULL;
    delete [] ivname[index];
    ivname[index] = NULL;

  } else if (flag == 1 && cols == 0) {
    dvector[index] = NULL;
    delete [] dvname[index];
    dvname[index] = NULL;

  } else if (flag == 0 && cols) {
    memory->destroy(iarray[index]);
    iarray[index] = NULL;
    delete [] ianame[index];
    ianame[index] = NULL;

  } else if (flag == 1 && cols) {
    memory->destroy(darray[index]);
    darray[index] = NULL;
    delete [] daname[index];
    daname[index] = NULL;
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

/* ----------------------------------------------------------------------
   perform sync and modify for each of 2 masks
   called by individual styles to override default sync/modify calls
     done at higher levels (Verlet,Modify,etc)
------------------------------------------------------------------------- */

void AtomKokkos::sync_modify(ExecutionSpace execution_space, unsigned int datamask_read,
                             unsigned int datamask_modify)
{
  sync(execution_space, datamask_read);
  modified(execution_space, datamask_modify);
}

AtomVec *AtomKokkos::new_avec(const std::string &style, int trysuffix, int &sflag)
{
  AtomVec *avec = Atom::new_avec(style, trysuffix, sflag);
  if (!avec->kokkosable) error->all(FLERR, "KOKKOS package requires a kokkos enabled atom_style");
  return avec;
}
