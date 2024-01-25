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

#include "reset_atoms_id.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory.h"
#include "modify.h"
#include "special.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#if defined(LMP_QSORT)
// allocate space for static class variable
// prototype for non-class function
ResetAtomsID::AtomRvous *ResetAtomsID::sortrvous;
static int compare_coords(const void *, const void *);
#else
// prototype for non-class function
static int compare_coords(const int, const int, void *);
#endif

static constexpr int PERBIN = 10;
static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

ResetAtomsID::ResetAtomsID(LAMMPS *lmp) : Command(lmp)
{
  binlo = binhi = -1;
}

/* ---------------------------------------------------------------------- */

void ResetAtomsID::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "Reset_atoms id command before simulation box is defined");
  if (atom->tag_enable == 0) error->all(FLERR, "Cannot use reset_atoms id unless atoms have IDs");

  for (const auto &ifix : modify->get_fix_list())
    if (ifix->stores_ids)
      error->all(FLERR, "Cannot use reset_atoms id with a fix {} storing atom IDs", ifix->style);

  if (comm->me == 0) utils::logmesg(lmp, "Resetting atom IDs ...\n");

  // process args

  int sortflag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "sort") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "reset_atoms id", error);
      sortflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown reset_atoms id keyword: {}", arg[iarg]);
  }

  // create an atom map if one doesn't exist already

  int mapflag = 0;
  if (atom->map_style == Atom::MAP_NONE) {
    mapflag = 1;
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // initialize system since comm->borders() will be invoked

  lmp->init();

  // setup domain, communication
  // acquire ghosts - is that really necessary?
  // exchange will clear map, borders will reset
  // this is the map needed to lookup current global IDs for bond topology

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  comm->exchange();
  comm->borders();
  if (domain->triclinic) domain->lamda2x(atom->nlocal + atom->nghost);

  // oldIDs = copy of current owned IDs

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  tagint *oldIDs;
  memory->create(oldIDs, nlocal, "reset_atom_ids:oldIDs");

  for (int i = 0; i < nlocal; i++) {
    oldIDs[i] = tag[i];
    tag[i] = 0;
  }

  // assign new contiguous IDs to owned atoms
  // if sortflag = no: use fast tag_extend()
  // if sortflag = yes: use slower full spatial sort plus rendezvous comm

  if (sortflag == 0)
    atom->tag_extend();
  else
    sort();

  // newIDs = copy of new IDs
  // restore old IDs, consistent with existing atom map
  // forward_comm_array acquires new IDs for ghost atoms

  double **newIDs;
  memory->create(newIDs, nall, 1, "reset_atom_ids:newIDs");

  for (int i = 0; i < nlocal; i++) {
    newIDs[i][0] = ubuf(tag[i]).d;
    tag[i] = oldIDs[i];
  }

  comm->forward_comm_array(1, newIDs);

  // loop over bonds, angles, etc and reset IDs in stored topology arrays
  // only necessary for molecular = Atom::MOLECULAR, not molecular = Atom::TEMPLATE
  // badcount = atom IDs that could not be found

  int badcount = 0;

  if (atom->molecular == Atom::MOLECULAR) {
    int j, m;
    tagint oldID;

    if (atom->avec->bonds_allow) {
      int *num_bond = atom->num_bond;
      tagint **bond_atom = atom->bond_atom;
      for (int i = 0; i < nlocal; i++) {
        for (j = 0; j < num_bond[i]; j++) {
          oldID = bond_atom[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            bond_atom[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;
        }
      }
    }

    if (atom->avec->angles_allow) {
      int *num_angle = atom->num_angle;
      tagint **angle_atom1 = atom->angle_atom1;
      tagint **angle_atom2 = atom->angle_atom2;
      tagint **angle_atom3 = atom->angle_atom3;
      for (int i = 0; i < nlocal; i++) {
        for (j = 0; j < num_angle[i]; j++) {
          oldID = angle_atom1[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            angle_atom1[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;

          oldID = angle_atom2[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            angle_atom2[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;

          oldID = angle_atom3[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            angle_atom3[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;
        }
      }
    }

    if (atom->avec->dihedrals_allow) {
      int *num_dihedral = atom->num_dihedral;
      tagint **dihedral_atom1 = atom->dihedral_atom1;
      tagint **dihedral_atom2 = atom->dihedral_atom2;
      tagint **dihedral_atom3 = atom->dihedral_atom3;
      tagint **dihedral_atom4 = atom->dihedral_atom4;
      for (int i = 0; i < nlocal; i++) {
        for (j = 0; j < num_dihedral[i]; j++) {
          oldID = dihedral_atom1[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            dihedral_atom1[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;

          oldID = dihedral_atom2[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            dihedral_atom2[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;

          oldID = dihedral_atom3[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            dihedral_atom3[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;

          oldID = dihedral_atom4[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            dihedral_atom4[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;
        }
      }
    }

    if (atom->avec->impropers_allow) {
      int *num_improper = atom->num_improper;
      tagint **improper_atom1 = atom->improper_atom1;
      tagint **improper_atom2 = atom->improper_atom2;
      tagint **improper_atom3 = atom->improper_atom3;
      tagint **improper_atom4 = atom->improper_atom4;
      for (int i = 0; i < nlocal; i++) {
        for (j = 0; j < num_improper[i]; j++) {
          oldID = improper_atom1[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            improper_atom1[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;

          oldID = improper_atom2[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            improper_atom2[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;

          oldID = improper_atom3[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            improper_atom3[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;

          oldID = improper_atom4[i][j];
          m = atom->map(oldID);
          if (m >= 0)
            improper_atom4[i][j] = (tagint) ubuf(newIDs[m][0]).i;
          else
            badcount++;
        }
      }
    }
  }

  // error check

  int all;
  MPI_Allreduce(&badcount, &all, 1, MPI_INT, MPI_SUM, world);
  if (all)
    error->all(FLERR,
               "reset_atoms id is missing {} bond topology atom IDs - use comm_modify cutoff", all);

  // reset IDs and atom map for owned atoms

  atom->map_clear();
  atom->nghost = 0;
  for (int i = 0; i < nlocal; i++) tag[i] = (tagint) ubuf(newIDs[i][0]).i;
  atom->map_init();
  atom->map_set();

  // need to update exclusions with new atom IDs

  if (atom->molecular == Atom::MOLECULAR) {
    Special special(lmp);
    special.build();
  }

  // delete temporary atom map

  if (mapflag) {
    atom->map_delete();
    atom->map_style = Atom::MAP_NONE;
  }

  // clean up

  memory->destroy(oldIDs);
  memory->destroy(newIDs);
}

/* ----------------------------------------------------------------------
   spatial sort of atoms followed by rendezvous comm to reset atom IDs
------------------------------------------------------------------------- */

void ResetAtomsID::sort()
{
  double mylo[3], myhi[3], bboxlo[3], bboxhi[3];

  int me = comm->me;
  int nprocs = comm->nprocs;
  int dim = domain->dimension;

  // bboxlo,bboxhi = bounding box on all atoms in system
  // expanded by 0.01 percent
  // bbox should work for orthogonal or triclinic system

  double **x = atom->x;
  int nlocal = atom->nlocal;

  mylo[0] = mylo[1] = mylo[2] = BIG;
  myhi[0] = myhi[1] = myhi[2] = -BIG;

  for (int i = 0; i < nlocal; i++) {
    mylo[0] = MIN(mylo[0], x[i][0]);
    mylo[1] = MIN(mylo[1], x[i][1]);
    mylo[2] = MIN(mylo[2], x[i][2]);
    myhi[0] = MAX(myhi[0], x[i][0]);
    myhi[1] = MAX(myhi[1], x[i][1]);
    myhi[2] = MAX(myhi[2], x[i][2]);
  }

  if (dim == 2) mylo[2] = myhi[2] = 0.0;

  // must ensure that bounding box volume is > 0.0

  for (int i = 0; i < 3; ++i) {
    if (mylo[i] == myhi[i]) {
      mylo[i] -= 0.5;
      myhi[i] += 0.5;
    }
  }

  MPI_Allreduce(mylo, bboxlo, 3, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(myhi, bboxhi, 3, MPI_DOUBLE, MPI_MAX, world);

  bboxlo[0] -= 0.0001 * (bboxhi[0] - bboxlo[0]);
  bboxlo[1] -= 0.0001 * (bboxhi[1] - bboxlo[1]);
  bboxlo[2] -= 0.0001 * (bboxhi[2] - bboxlo[2]);
  bboxhi[0] += 0.0001 * (bboxhi[0] - bboxlo[0]);
  bboxhi[1] += 0.0001 * (bboxhi[1] - bboxlo[1]);
  bboxhi[2] += 0.0001 * (bboxhi[2] - bboxlo[2]);

  // nbin_estimate = estimate of total number of bins, each with PERBIN atoms
  // binsize = edge length of a cubic bin
  // nbin xyz = bin count in each dimension

  bigint nbin_estimate = atom->natoms / PERBIN + 1;

  double vol;
  if (dim == 2)
    vol = (bboxhi[0] - bboxlo[0]) * (bboxhi[1] - bboxlo[1]);
  else
    vol = (bboxhi[0] - bboxlo[0]) * (bboxhi[1] - bboxlo[1]) * (bboxhi[2] - bboxlo[2]);
  double binsize = pow(vol / nbin_estimate, 1.0 / dim);

  int nbinx = static_cast<int>((bboxhi[0] - bboxlo[0]) / binsize) + 1;
  int nbiny = static_cast<int>((bboxhi[1] - bboxlo[1]) / binsize) + 1;
  int nbinz = static_cast<int>((bboxhi[2] - bboxlo[2]) / binsize) + 1;

  double invx = 1.0 / (bboxhi[0] - bboxlo[0]);
  double invy = 1.0 / (bboxhi[1] - bboxlo[1]);
  double invz;
  if (dim == 2)
    invz = 0.0;
  else
    invz = 1.0 / (bboxhi[2] - bboxlo[2]);

  // nbins = actual total number of bins
  // bins are split evenly across procs
  // lo-numbered procs may have own one bin less than rest of procs
  // nlo = # of bins/proc for lowest procs
  // nhi = nlo+1 = # of bins/proc for highest procs
  // nplo = # of procs in low group
  // nbinlo = # of bins owned by low group procs
  // binlo to binhi-1 = bin indices this proc will own in Rvous decomp
  // bins are numbered from 0 to Nbins-1

  bigint nbins = (bigint) nbinx * nbiny * nbinz;
  bigint nlo = nbins / nprocs;
  bigint nhi = nlo + 1;
  bigint nplo = nprocs - (nbins % nprocs);
  bigint nbinlo = nplo * nlo;

  if (me < nplo) {
    binlo = me * nlo;
    binhi = (me + 1) * nlo;
  } else {
    binlo = nbinlo + (me - nplo) * nhi;
    binhi = nbinlo + (me + 1 - nplo) * nhi;
  }

  // fill atombuf with info on my atoms
  // ibin = which bin the atom is in
  // proclist = proc that owns ibin

  int *proclist;
  memory->create(proclist, nlocal, "special:proclist");
  auto atombuf =
      (AtomRvous *) memory->smalloc((bigint) nlocal * sizeof(AtomRvous), "resetIDs:idbuf");

  int ibinx, ibiny, ibinz, iproc;
  bigint ibin;

  for (int i = 0; i < nlocal; i++) {
    ibinx = static_cast<int>((x[i][0] - bboxlo[0]) * invx * nbinx);
    ibiny = static_cast<int>((x[i][1] - bboxlo[1]) * invy * nbiny);
    ibinz = static_cast<int>((x[i][2] - bboxlo[2]) * invz * nbinz);
    ibin = (bigint) ibinz * nbiny * nbinx + (bigint) ibiny * nbinx + ibinx;

    if (ibin < nbinlo)
      iproc = ibin / nlo;
    else
      iproc = nplo + (ibin - nbinlo) / nhi;
    proclist[i] = iproc;

    atombuf[i].ibin = ibin;
    atombuf[i].proc = me;
    atombuf[i].ilocal = i;
    atombuf[i].x[0] = x[i][0];
    atombuf[i].x[1] = x[i][1];
    atombuf[i].x[2] = x[i][2];
  }

  // perform rendezvous operation, send atombuf to other procs

  char *buf;
  int nreturn = comm->rendezvous(1, nlocal, (char *) atombuf, sizeof(AtomRvous), 0, proclist,
                                 sort_bins, 0, buf, sizeof(IDRvous), (void *) this);
  auto outbuf = (IDRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(atombuf);

  // set new ID for all owned atoms

  int ilocal;

  for (int i = 0; i < nreturn; i++) {
    ilocal = outbuf[i].ilocal;
    atom->tag[ilocal] = outbuf[i].newID;
  }

  memory->sfree(outbuf);
}

/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N AtomRvous datums
   outbuf = list of N IDRvous datums, sent back to sending proc
------------------------------------------------------------------------- */

int ResetAtomsID::sort_bins(int n, char *inbuf, int &flag, int *&proclist, char *&outbuf, void *ptr)
{
  int i, ibin, index;

  auto rptr = (ResetAtomsID *) ptr;
  Memory *memory = rptr->memory;
  Error *error = rptr->error;
  MPI_Comm world = rptr->world;
  bigint binlo = rptr->binlo;
  bigint binhi = rptr->binhi;

  if ((binlo < 0) || (binhi < 0))
    error->one(FLERR, "Called sort_bins without previous setup of bins");

  // nbins = my subset of bins from binlo to binhi-1
  // initialize linked lists of my Rvous atoms in each bin

  int *next, *head, *last, *count;

  int nbins = binhi - binlo;
  memory->create(next, n, "resetIDs:next");
  memory->create(head, nbins, "resetIDs:head");
  memory->create(last, nbins, "resetIDs:last");
  memory->create(count, nbins, "resetIDs:count");

  for (i = 0; i < n; i++) next[i] = -1;

  for (ibin = 0; ibin < nbins; ibin++) {
    head[ibin] = last[ibin] = -1;
    count[ibin] = 0;
  }

  auto in = (AtomRvous *) inbuf;

  for (i = 0; i < n; i++) {
    if (in[i].ibin < binlo || in[i].ibin >= binhi) {
      //printf("BAD me %d i %d %d ibin %d binlohi %d %d\n",
      //             rptr->comm->me,i,n,in[i].ibin,binlo,binhi);
      error->one(FLERR, "Bad spatial bin assignment in reset_atoms id sort");
    }
    ibin = in[i].ibin - binlo;
    if (head[ibin] < 0) head[ibin] = i;
    if (last[ibin] >= 0) next[last[ibin]] = i;
    last[ibin] = i;
    count[ibin]++;
  }

  int maxcount = 0;
  for (ibin = 0; ibin < nbins; ibin++) maxcount = MAX(maxcount, count[ibin]);

  int *order;
  memory->create(order, maxcount, "resetIDs:order");

  // sort atoms in each bin by their spatial position
  // recreate linked list for each bin based on sorted order

  for (ibin = 0; ibin < nbins; ibin++) {
    index = head[ibin];
    for (i = 0; i < count[ibin]; i++) {
      order[i] = index;
      index = next[index];
    }

#if defined(LMP_QSORT)
    sortrvous = in;
    qsort(order, count[ibin], sizeof(int), compare_coords);
#else
    utils::merge_sort(order, count[ibin], (void *) in, compare_coords);
#endif

    head[ibin] = last[ibin] = -1;
    for (i = 0; i < count[ibin]; i++) {
      if (head[ibin] < 0) head[ibin] = order[i];
      if (last[ibin] >= 0) next[last[ibin]] = order[i];
      last[ibin] = order[i];
    }
  }

  // MPI_Scan() to find out where my atoms begin in globally sorted list

  tagint ntag = n;
  tagint nprev;
  MPI_Scan(&ntag, &nprev, 1, MPI_LMP_TAGINT, MPI_SUM, world);
  nprev -= n;

  // loop over bins and sorted atoms in each bin, reset ids to be consecutive
  // pass outbuf of IDRvous datums back to comm->rendezvous

  int nout = n;
  memory->create(proclist, nout, "resetIDs:proclist");
  auto out = (IDRvous *) memory->smalloc(nout * sizeof(IDRvous), "resetIDs:out");

  tagint one = nprev + 1;
  for (ibin = 0; ibin < nbins; ibin++) {
    index = head[ibin];
    for (i = 0; i < count[ibin]; i++) {
      proclist[index] = in[index].proc;
      out[index].newID = one++;
      out[index].ilocal = in[index].ilocal;
      index = next[index];
    }
  }

  outbuf = (char *) out;

  // clean up
  // Comm::rendezvous will delete proclist and out (outbuf)

  memory->destroy(next);
  memory->destroy(head);
  memory->destroy(last);
  memory->destroy(count);
  memory->destroy(order);

  // flag = 2: new outbuf

  flag = 2;
  return nout;
}

#if defined(LMP_QSORT)

/* ----------------------------------------------------------------------
   comparison function invoked by qsort()
   accesses static class member sortrvous, set before call to qsort()
------------------------------------------------------------------------- */

int compare_coords(const void *iptr, const void *jptr)
{
  int i = *((int *) iptr);
  int j = *((int *) jptr);
  ResetAtomsID::AtomRvous *rvous = ResetAtomsID::sortrvous;
  double *xi = rvous[i].x;
  double *xj = rvous[j].x;
  if (xi[0] < xj[0]) return -1;
  if (xi[0] > xj[0]) return 1;
  if (xi[1] < xj[1]) return -1;
  if (xi[1] > xj[1]) return 1;
  if (xi[2] < xj[2]) return -1;
  if (xi[2] > xj[2]) return 1;
  return 0;
}

#else

/* ----------------------------------------------------------------------
   comparison function invoked by merge_sort()
   void pointer contains sortrvous
------------------------------------------------------------------------- */

int compare_coords(const int i, const int j, void *ptr)
{
  auto rvous = (ResetAtomsID::AtomRvous *) ptr;
  double *xi = rvous[i].x;
  double *xj = rvous[j].x;
  if (xi[0] < xj[0]) return -1;
  if (xi[0] > xj[0]) return 1;
  if (xi[1] < xj[1]) return -1;
  if (xi[1] > xj[1]) return 1;
  if (xi[2] < xj[2]) return -1;
  if (xi[2] > xj[2]) return 1;
  return 0;
}

#endif
