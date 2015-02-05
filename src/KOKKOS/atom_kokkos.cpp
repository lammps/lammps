/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "atom_kokkos.h"
#include "atom_vec.h"
#include "atom_vec_kokkos.h"
#include "comm_kokkos.h"
#include "update.h"
#include "domain.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomKokkos::AtomKokkos(LAMMPS *lmp) : Atom(lmp)
{
  // set CommKokkos pointer to Atom class, since CommKokkos allocated first

  ((CommKokkos *) comm)->atomKK = this;
}

/* ---------------------------------------------------------------------- */

AtomKokkos::~AtomKokkos()
{
  memory->destroy_kokkos(k_tag, tag);
  memory->destroy_kokkos(k_mask, mask);
  memory->destroy_kokkos(k_type, type);
  memory->destroy_kokkos(k_image, image);
  memory->destroy_kokkos(k_molecule, molecule);

  memory->destroy_kokkos(k_x, x);
  memory->destroy_kokkos(k_v, v);
  memory->destroy_kokkos(k_f, f);

  memory->destroy_kokkos(k_mass, mass);
  memory->destroy_kokkos(k_q, q);

  memory->destroy_kokkos(k_radius, radius);
  memory->destroy_kokkos(k_rmass, rmass);
  memory->destroy_kokkos(k_omega, omega);
  memory->destroy_kokkos(k_torque, torque);

  memory->destroy_kokkos(k_nspecial, nspecial);
  memory->destroy_kokkos(k_special, special);
  memory->destroy_kokkos(k_num_bond, num_bond);
  memory->destroy_kokkos(k_bond_type, bond_type);
  memory->destroy_kokkos(k_bond_atom, bond_atom);
  memory->destroy_kokkos(k_num_angle, num_angle);
  memory->destroy_kokkos(k_angle_type, angle_type);
  memory->destroy_kokkos(k_angle_atom1, angle_atom1);
  memory->destroy_kokkos(k_angle_atom2, angle_atom2);
  memory->destroy_kokkos(k_angle_atom3, angle_atom3);
  memory->destroy_kokkos(k_num_dihedral, num_dihedral);
  memory->destroy_kokkos(k_dihedral_type, dihedral_type);
  memory->destroy_kokkos(k_dihedral_atom1, dihedral_atom1);
  memory->destroy_kokkos(k_dihedral_atom2, dihedral_atom2);
  memory->destroy_kokkos(k_dihedral_atom3, dihedral_atom3);
  memory->destroy_kokkos(k_dihedral_atom4, dihedral_atom4);
  memory->destroy_kokkos(k_num_improper, num_improper);
  memory->destroy_kokkos(k_improper_type, improper_type);
  memory->destroy_kokkos(k_improper_atom1, improper_atom1);
  memory->destroy_kokkos(k_improper_atom2, improper_atom2);
  memory->destroy_kokkos(k_improper_atom3, improper_atom3);
  memory->destroy_kokkos(k_improper_atom4, improper_atom4);
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sync(const ExecutionSpace space, unsigned int mask)
{
  ((AtomVecKokkos *) avec)->sync(space,mask);
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::modified(const ExecutionSpace space, unsigned int mask)
{
  ((AtomVecKokkos *) avec)->modified(space,mask);
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::allocate_type_arrays()
{
  if (avec->mass_type) {
    k_mass = DAT::tdual_float_1d("Mass",ntypes+1);
    mass = k_mass.h_view.ptr_on_device();
    mass_setflag = new int[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) mass_setflag[itype] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::sort()
{
  int i,m,n,ix,iy,iz,ibin,empty;

  // set next timestep for sorting to take place

  nextsort = (update->ntimestep/sortfreq)*sortfreq + sortfreq;

  // re-setup sort bins if needed

  if (domain->box_change) setup_sort_bins();
  if (nbins == 1) return;

  // reallocate per-atom vectors if needed

  if (nlocal > maxnext) {
    memory->destroy(next);
    memory->destroy(permute);
    maxnext = atom->nmax;
    memory->create(next,maxnext,"atom:next");
    memory->create(permute,maxnext,"atom:permute");
  }

  // insure there is one extra atom location at end of arrays for swaps

  if (nlocal == nmax) avec->grow(0);

  sync(Host,ALL_MASK);
  modified(Host,ALL_MASK);

  // bin atoms in reverse order so linked list will be in forward order

  for (i = 0; i < nbins; i++) binhead[i] = -1;

  HAT::t_x_array_const h_x = k_x.view<LMPHostType>();
  for (i = nlocal-1; i >= 0; i--) {
    ix = static_cast<int> ((h_x(i,0)-bboxlo[0])*bininvx);
    iy = static_cast<int> ((h_x(i,1)-bboxlo[1])*bininvy);
    iz = static_cast<int> ((h_x(i,2)-bboxlo[2])*bininvz);
    ix = MAX(ix,0);
    iy = MAX(iy,0);
    iz = MAX(iz,0);
    ix = MIN(ix,nbinx-1);
    iy = MIN(iy,nbiny-1);
    iz = MIN(iz,nbinz-1);
    ibin = iz*nbiny*nbinx + iy*nbinx + ix;
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
    avec->copy(i,nlocal,0);
    empty = i;
    while (permute[empty] != i) {
      avec->copy(permute[empty],empty,0);
      empty = current[empty] = permute[empty];
    }
    avec->copy(nlocal,empty,0);
    current[empty] = permute[empty];
  }

  // sanity check that current = permute

  //int flag = 0;
  //for (i = 0; i < nlocal; i++)
  //  if (current[i] != permute[i]) flag = 1;
  //int flagall;
  //MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  //if (flagall) error->all(FLERR,"Atom sort did not operate correctly");
}

/* ----------------------------------------------------------------------
   reallocate memory to the pointer selected by the mask
------------------------------------------------------------------------- */

void AtomKokkos::grow(unsigned int mask){

  if (mask & SPECIAL_MASK){
    memory->destroy_kokkos(k_special, special);
    sync(Device, mask);
    modified(Device, mask);
    memory->grow_kokkos(k_special,special,nmax,maxspecial,"atom:special");
    avec->grow_reset();
    sync(Host, mask);
  }

}

/* ---------------------------------------------------------------------- */

void AtomKokkos::deallocate_topology()
{
  memory->destroy_kokkos(k_bond_type, bond_type);
  memory->destroy_kokkos(k_bond_atom, bond_atom);

  memory->destroy_kokkos(k_angle_type, angle_type);
  memory->destroy_kokkos(k_angle_atom1, angle_atom1);
  memory->destroy_kokkos(k_angle_atom2, angle_atom2);
  memory->destroy_kokkos(k_angle_atom3, angle_atom3);

  memory->destroy_kokkos(k_dihedral_type, dihedral_type);
  memory->destroy_kokkos(k_dihedral_atom1, dihedral_atom1);
  memory->destroy_kokkos(k_dihedral_atom2, dihedral_atom2);
  memory->destroy_kokkos(k_dihedral_atom3, dihedral_atom3);
  memory->destroy_kokkos(k_dihedral_atom4, dihedral_atom4);

  memory->destroy_kokkos(k_improper_type, improper_type);
  memory->destroy_kokkos(k_improper_atom1, improper_atom1);
  memory->destroy_kokkos(k_improper_atom2, improper_atom2);
  memory->destroy_kokkos(k_improper_atom3, improper_atom3);
  memory->destroy_kokkos(k_improper_atom4, improper_atom4);
}

/* ----------------------------------------------------------------------
   perform sync and modify for each of 2 masks
   called by individual styles to override default sync/modify calls
     done at higher levels (Verlet,Modify,etc)
------------------------------------------------------------------------- */

void AtomKokkos::sync_modify(ExecutionSpace execution_space,
                             unsigned int datamask_read, 
                             unsigned int datamask_modify)
{
  sync(execution_space,datamask_read);
  modified(execution_space,datamask_modify);
}
