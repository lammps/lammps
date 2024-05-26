// clang-format off
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

#include "special.h"

#include "accelerator_kokkos.h"  // IWYU pragma: export
#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "comm.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"

using namespace LAMMPS_NS;

static constexpr int RVOUS = 1;   // 0 for irregular, 1 for all2all

/* ---------------------------------------------------------------------- */

Special::Special(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  onetwo = onethree = onefour = onefive = nullptr;
}

/* ---------------------------------------------------------------------- */

Special::~Special()
{
  memory->destroy(onetwo);
  memory->destroy(onethree);
  memory->destroy(onefour);
  memory->destroy(onefive);
}

/* ----------------------------------------------------------------------
   create 1-2, 1-3, 1-4 lists of topology neighbors, 1-5 list is optional
   store in onetwo, onethree, onefour, onefive for each atom
   store first 3 counters in nspecial[i], and 4th in nspecial15[i]
------------------------------------------------------------------------- */

void Special::build()
{
  MPI_Barrier(world);
  double time1 = platform::walltime();

  if (me == 0) {
    const double * const special_lj   = force->special_lj;
    const double * const special_coul = force->special_coul;
    utils::logmesg(lmp, "Finding 1-2 1-3 1-4 neighbors ...\n"
                   "  special bond factors lj:    {:<8} {:<8} {:<8}\n"
                   "  special bond factors coul:  {:<8} {:<8} {:<8}\n",
                   special_lj[1],special_lj[2],special_lj[3],
                   special_coul[1],special_coul[2],special_coul[3]);
  }

  // set onefive_flag if special_bonds command set it

  onefive_flag = force->special_onefive;

  // initialize nspecial counters to 0

  int **nspecial = atom->nspecial;
  int *nspecial15 = atom->nspecial15;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    nspecial[i][0] = 0;
    nspecial[i][1] = 0;
    nspecial[i][2] = 0;
  }

  if (onefive_flag) {
    for (int i = 0; i < nlocal; i++) nspecial15[i] = 0;
  }

  // setup atomIDs and procowner vectors in rendezvous decomposition

  atom_owners();

  // tally nspecial[i][0] = # of 1-2 neighbors of atom i
  // create onetwo[i] = list of 1-2 neighbors for atom i

  if (force->newton_bond) onetwo_build_newton();
  else onetwo_build_newton_off();

  // print max # of 1-2 neighbors

  if (me == 0)
    utils::logmesg(lmp,"{:>6} = max # of 1-2 neighbors\n",maxall);

  // done if special_bond weights for 1-3, 1-4 are set to 1.0
  // onefive_flag must also be off, else 1-4 is needed to create 1-5

  if (!onefive_flag &&
      force->special_lj[2] == 1.0 && force->special_coul[2] == 1.0 &&
      force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) {
    dedup();
    combine();
    fix_alteration();
    memory->destroy(procowner);
    memory->destroy(atomIDs);
    timer_output(time1);
    return;
  }

  // tally nspecial[i][1] = # of 1-3 neighbors of atom i
  // create onethree[i] = list of 1-3 neighbors for atom i

  onethree_build();

  // print max # of 1-3 neighbors

  if (me == 0)
    utils::logmesg(lmp,"{:>6} = max # of 1-3 neighbors\n",maxall);

  // done if special_bond weights for 1-4 are set to 1.0
  // onefive_flag must also be off, else 1-4 is needed to create 1-5

  if (!onefive_flag &&
      force->special_lj[3] == 1.0 && force->special_coul[3] == 1.0) {
    dedup();
    if (force->special_angle) angle_trim();
    combine();
    fix_alteration();
    memory->destroy(procowner);
    memory->destroy(atomIDs);
    timer_output(time1);
    return;
  }

  // tally nspecial[i][2] = # of 1-4 neighbors of atom i
  // create onefour[i] = list of 1-4 neighbors for atom i

  onefour_build();

  // print max # of 1-4 neighbors

  if (me == 0)
    utils::logmesg(lmp,"{:>6} = max # of 1-4 neighbors\n",maxall);

  // optionally store 1-5 neighbors
  // tally nspecial15[i] = # of 1-5 neighbors of atom i
  // create onefive[i] = list of 1-5 neighbors for atom i

  if (onefive_flag) {
    onefive_build();
    if (me == 0)
      utils::logmesg(lmp,"{:>6} = max # of 1-5 neighbors\n",maxall);
  }

  // finish processing the onetwo, onethree, onefour, onefive lists

  dedup();
  if (force->special_angle) angle_trim();
  if (force->special_dihedral) dihedral_trim();
  combine();
  fix_alteration();
  memory->destroy(procowner);
  memory->destroy(atomIDs);
  timer_output(time1);
}

/* ----------------------------------------------------------------------
   setup atomIDs and procowner
------------------------------------------------------------------------- */

void Special::atom_owners()
{
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int *proclist;
  memory->create(proclist,nlocal,"special:proclist");
  auto idbuf = (IDRvous *) memory->smalloc((bigint) nlocal*sizeof(IDRvous),"special:idbuf");

  // setup input buf for rendezvous comm
  // one datum for each owned atom: datum = owning proc, atomID
  // each proc assigned every 1/Pth atom

  for (int i = 0; i < nlocal; i++) {
    proclist[i] = tag[i] % nprocs;
    idbuf[i].me = me;
    idbuf[i].atomID = tag[i];
  }

  // perform rendezvous operation

  char *buf;
  comm->rendezvous(RVOUS,nlocal,(char *) idbuf,sizeof(IDRvous),0,proclist,
                   rendezvous_ids,0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->sfree(idbuf);
}

/* ----------------------------------------------------------------------
   onetwo build when newton_bond flag on
   uses rendezvous comm
------------------------------------------------------------------------- */

void Special::onetwo_build_newton()
{
  int i,j,m;

  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  // nsend = # of my datums to send

  int nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < num_bond[i]; j++) {
      m = atom->map(bond_atom[i][j]);
      if (m < 0 || m >= nlocal) nsend++;
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"special:proclist");
  auto inbuf = (PairRvous *) memory->smalloc((bigint) nsend*sizeof(PairRvous),"special:inbuf");

  // setup input buf to rendezvous comm
  // one datum for each unowned bond partner: bond partner ID, atomID
  // owning proc for each datum = bond partner ID % nprocs

  nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < num_bond[i]; j++) {
      m = atom->map(bond_atom[i][j]);
      if (m >= 0 && m < nlocal) continue;
      proclist[nsend] = bond_atom[i][j] % nprocs;
      inbuf[nsend].atomID = bond_atom[i][j];
      inbuf[nsend].partnerID = tag[i];
      nsend++;
    }
  }

  // perform rendezvous operation

  char *buf;
  int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(PairRvous), 0,proclist,
                                 rendezvous_pairs,0,buf,sizeof(PairRvous), (void *) this);
  auto outbuf = (PairRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set nspecial[0] and onetwo for all owned atoms
  // based on owned info plus rendezvous output info
  // output datums = pairs of atoms that are 1-2 neighbors

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < num_bond[i]; j++) {
      nspecial[i][0]++;
      m = atom->map(bond_atom[i][j]);
      if (m >= 0 && m < nlocal) nspecial[m][0]++;
    }
  }

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    nspecial[i][0]++;
  }

  int max = 0;
  for (i = 0; i < nlocal; i++)
    max = MAX(max,nspecial[i][0]);

  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
  memory->create(onetwo,nlocal,maxall,"special:onetwo");

  for (i = 0; i < nlocal; i++) nspecial[i][0] = 0;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < num_bond[i]; j++) {
      onetwo[i][nspecial[i][0]++] = bond_atom[i][j];
      m = atom->map(bond_atom[i][j]);
      if (m >= 0 && m < nlocal) onetwo[m][nspecial[m][0]++] = tag[i];
    }
  }

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    onetwo[i][nspecial[i][0]++] = outbuf[m].partnerID;
  }

  memory->sfree(outbuf);
}

/* ----------------------------------------------------------------------
   onetwo build with newton_bond flag off
   no need for rendezvous comm
------------------------------------------------------------------------- */

void Special::onetwo_build_newton_off()
{
  int i,j;

  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  int max = 0;
  for (i = 0; i < nlocal; i++)
    max = MAX(max,num_bond[i]);

  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
  memory->create(onetwo,nlocal,maxall,"special:onetwo");

  // nsend = # of my datums to send
  // include nlocal datums with owner of each atom

  for (i = 0; i < nlocal; i++) {
    nspecial[i][0] = num_bond[i];
    for (j = 0; j < num_bond[i]; j++)
      onetwo[i][j] = bond_atom[i][j];
  }
}

/* ----------------------------------------------------------------------
   onethree build
   uses rendezvous comm
------------------------------------------------------------------------- */

void Special::onethree_build()
{
  int i,j,k,m,proc;

  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  // nsend = # of my datums to send

  int nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][0]; j++) {
      m = atom->map(onetwo[i][j]);
      if (m < 0 || m >= nlocal) nsend += nspecial[i][0]-1;
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"special:proclist");
  auto inbuf = (PairRvous *) memory->smalloc((bigint) nsend*sizeof(PairRvous),"special:inbuf");

  // setup input buf to rendezvous comm
  // datums = pairs of onetwo partners where either is unknown
  //          these pairs are onethree neighbors
  //          datum = onetwo ID, onetwo ID
  // owning proc for each datum = onetwo ID % nprocs

  nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][0]; j++) {
      m = atom->map(onetwo[i][j]);
      if (m >= 0 && m < nlocal) continue;
      proc = onetwo[i][j] % nprocs;
      for (k = 0; k < nspecial[i][0]; k++) {
        if (j == k) continue;
        proclist[nsend] = proc;
        inbuf[nsend].atomID = onetwo[i][j];
        inbuf[nsend].partnerID = onetwo[i][k];
        nsend++;
      }
    }
  }

  // perform rendezvous operation

  char *buf;
  int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(PairRvous), 0,proclist,
                                 rendezvous_pairs,0,buf,sizeof(PairRvous), (void *) this);
  auto outbuf = (PairRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set nspecial[1] and onethree for all owned atoms
  // based on owned info plus rendezvous output info
  // output datums = pairs of atoms that are 1-3 neighbors

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][0]; j++) {
      m = atom->map(onetwo[i][j]);
      if (m >= 0 && m < nlocal) nspecial[m][1] += nspecial[i][0]-1;
    }
  }

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    nspecial[i][1]++;
  }

  int max = 0;
  for (i = 0; i < nlocal; i++)
    max = MAX(max,nspecial[i][1]);

  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
  memory->create(onethree,nlocal,maxall,"special:onethree");

  for (i = 0; i < nlocal; i++) nspecial[i][1] = 0;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][0]; j++) {
      m = atom->map(onetwo[i][j]);
      if (m < 0 || m >= nlocal) continue;
      for (k = 0; k < nspecial[i][0]; k++) {
        if (j == k) continue;
        onethree[m][nspecial[m][1]++] = onetwo[i][k];
      }
    }
  }

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    onethree[i][nspecial[i][1]++] = outbuf[m].partnerID;
  }

  memory->sfree(outbuf);
}

/* ----------------------------------------------------------------------
   onefour build
   uses rendezvous comm
------------------------------------------------------------------------- */

void Special::onefour_build()
{
  int i,j,k,m,proc;

  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;

  // nsend = # of my datums to send

  int nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][1]; j++) {
      m = atom->map(onethree[i][j]);
      if (m < 0 || m >= nlocal) nsend += nspecial[i][0];
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"special:proclist");
  auto inbuf = (PairRvous *) memory->smalloc((bigint) nsend*sizeof(PairRvous),"special:inbuf");

  // setup input buf to rendezvous comm
  // datums = pairs of onethree and onetwo partners where onethree is unknown
  //          these pairs are onefour neighbors
  //          datum = onethree ID, onetwo ID
  // owning proc for each datum = onethree ID % nprocs

  nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][1]; j++) {
      m = atom->map(onethree[i][j]);
      if (m >= 0 && m < nlocal) continue;
      proc = onethree[i][j] % nprocs;
      for (k = 0; k < nspecial[i][0]; k++) {
        proclist[nsend] = proc;
        inbuf[nsend].atomID = onethree[i][j];
        inbuf[nsend].partnerID = onetwo[i][k];
        nsend++;
      }
    }
  }

  // perform rendezvous operation

  char *buf;
  int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(PairRvous), 0,proclist,
                                 rendezvous_pairs,0,buf,sizeof(PairRvous), (void *) this);
  auto outbuf = (PairRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set nspecial[2] and onefour for all owned atoms
  // based on owned info plus rendezvous output info
  // output datums = pairs of atoms that are 1-4 neighbors

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][1]; j++) {
      m = atom->map(onethree[i][j]);
      if (m >= 0 && m < nlocal) nspecial[m][2] += nspecial[i][0];
    }
  }

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    nspecial[i][2]++;
  }

  int max = 0;
  for (i = 0; i < nlocal; i++)
    max = MAX(max,nspecial[i][2]);

  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
  memory->create(onefour,nlocal,maxall,"special:onefour");

  for (i = 0; i < nlocal; i++) nspecial[i][2] = 0;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][1]; j++) {
      m = atom->map(onethree[i][j]);
      if (m < 0 || m >= nlocal) continue;
      for (k = 0; k < nspecial[i][0]; k++) {
        onefour[m][nspecial[m][2]++] = onetwo[i][k];
      }
    }
  }

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    onefour[i][nspecial[i][2]++] = outbuf[m].partnerID;
  }

  memory->sfree(outbuf);
}

/* ----------------------------------------------------------------------
   optional onefive build
   uses rendezvous comm
------------------------------------------------------------------------- */

void Special::onefive_build()
{
  int i,j,k,m,proc;

  int **nspecial = atom->nspecial;
  int *nspecial15 = atom->nspecial15;
  int nlocal = atom->nlocal;

  // nsend = # of my datums to send

  int nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][2]; j++) {
      m = atom->map(onefour[i][j]);
      if (m < 0 || m >= nlocal) nsend += nspecial[i][0];
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"special:proclist");
  PairRvous *inbuf = (PairRvous *)
    memory->smalloc((bigint) nsend*sizeof(PairRvous),"special:inbuf");

  // setup input buf to rendezvous comm
  // datums = pairs of onefour and onetwo partners where onefour is unknown
  //          these pairs are onefive neighbors
  //          datum = onefour ID, onetwo ID
  // owning proc for each datum = onefour ID % nprocs

  nsend = 0;
  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][2]; j++) {
      m = atom->map(onefour[i][j]);
      if (m >= 0 && m < nlocal) continue;
      proc = onefour[i][j] % nprocs;
      for (k = 0; k < nspecial[i][0]; k++) {
        proclist[nsend] = proc;
        inbuf[nsend].atomID = onefour[i][j];
        inbuf[nsend].partnerID = onetwo[i][k];
        nsend++;
      }
    }
  }

  // perform rendezvous operation

  char *buf;
  int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(PairRvous),
                                 0,proclist,
                                 rendezvous_pairs,0,buf,sizeof(PairRvous),
                                 (void *) this);
  PairRvous *outbuf = (PairRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set nspecial15 and onefive for all owned atoms
  // based on owned info plus rendezvous output info
  // output datums = pairs of atoms that are 1-5 neighbors

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][2]; j++) {
      m = atom->map(onefour[i][j]);
      if (m >= 0 && m < nlocal) nspecial15[m] += nspecial[i][0];
    }
  }

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    nspecial15[i]++;
  }

  int max = 0;
  for (i = 0; i < nlocal; i++)
    max = MAX(max,nspecial15[i]);

  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
  memory->create(onefive,nlocal,maxall,"special:onefive");

  for (i = 0; i < nlocal; i++) nspecial15[i] = 0;

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < nspecial[i][2]; j++) {
      m = atom->map(onefour[i][j]);
      if (m < 0 || m >= nlocal) continue;
      for (k = 0; k < nspecial[i][0]; k++) {
        onefive[m][nspecial15[m]++] = onetwo[i][k];
      }
    }
  }

  for (m = 0; m < nreturn; m++) {
    i = atom->map(outbuf[m].atomID);
    onefive[i][nspecial15[i]++] = outbuf[m].partnerID;
  }

  memory->sfree(outbuf);
}

/* ----------------------------------------------------------------------
   remove duplicates within each of onetwo, onethree, onefour individually
   also dedup onefive if enabled
------------------------------------------------------------------------- */

void Special::dedup()
{
  int i,j;
  tagint m;

  // clear map so it can be used as scratch space

  atom->map_clear();

  // use map to cull duplicates
  // exclude original atom explicitly
  // adjust onetwo, onethree, onefour, onefive values to remove duplicates
  // must unset map for each atom

  int **nspecial = atom->nspecial;
  int *nspecial15 = atom->nspecial15;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int unique;

  // dedup onetwo

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);
    for (j = 0; j < nspecial[i][0]; j++) {
      m = onetwo[i][j];
      if (atom->map(m) < 0) {
        onetwo[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][0] = unique;
    atom->map_one(tag[i],-1);
    for (j = 0; j < unique; j++) atom->map_one(onetwo[i][j],-1);
  }

  // dedup onethree

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);
    for (j = 0; j < nspecial[i][1]; j++) {
      m = onethree[i][j];
      if (atom->map(m) < 0) {
        onethree[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][1] = unique;
    atom->map_one(tag[i],-1);
    for (j = 0; j < unique; j++) atom->map_one(onethree[i][j],-1);
  }

  // dedup onefour

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);
    for (j = 0; j < nspecial[i][2]; j++) {
      m = onefour[i][j];
      if (atom->map(m) < 0) {
        onefour[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][2] = unique;
    atom->map_one(tag[i],-1);
    for (j = 0; j < unique; j++) atom->map_one(onefour[i][j],-1);
  }

  // dedup onefive

  if (onefive_flag) {
    for (i = 0; i < nlocal; i++) {
      unique = 0;
      atom->map_one(tag[i],0);
      for (j = 0; j < nspecial15[i]; j++) {
        m = onefive[i][j];
        if (atom->map(m) < 0) {
          onefive[i][unique++] = m;
          atom->map_one(m,0);
        }
      }
      nspecial15[i] = unique;
      atom->map_one(tag[i],-1);
      for (j = 0; j < unique; j++) atom->map_one(onefive[i][j],-1);
    }
  }

  // re-create map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();
}

/* ----------------------------------------------------------------------
   concatenate onetwo, onethree, onefour into master atom->special list
   remove duplicates between 3 lists, leave dup in first list it appears in
   convert nspecial[0], nspecial[1], nspecial[2] into cumulative counters
   if 1-5 is enabled, reset nspecial15/special15 to remove dups with 1-2,1-3,1-4
------------------------------------------------------------------------- */

void Special::combine()
{
  int i,j;
  tagint m;

  int me;
  MPI_Comm_rank(world,&me);

  int **nspecial = atom->nspecial;
  int *nspecial15 = atom->nspecial15;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // ----------------------------------------------------
  // compute culled maxspecial = max # of special neighs of any atom
  // ----------------------------------------------------

  // clear map so it can be used as scratch space

  atom->map_clear();

  // unique = # of unique nspecial neighbors of one atom
  // unique15 = ditto for 1-5 interactions
  // cull duplicates using map to check for them
  // exclude original atom explicitly
  // must unset map for each atom

  int unique,unique15;
  int maxspecial = 0;
  int maxspecial15 = 0;

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);

    for (j = 0; j < nspecial[i][0]; j++) {
      m = onetwo[i][j];
      if (atom->map(m) < 0) {
        unique++;
        atom->map_one(m,0);
      }
    }
    for (j = 0; j < nspecial[i][1]; j++) {
      m = onethree[i][j];
      if (atom->map(m) < 0) {
        unique++;
        atom->map_one(m,0);
      }
    }
    for (j = 0; j < nspecial[i][2]; j++) {
      m = onefour[i][j];
      if (atom->map(m) < 0) {
        unique++;
        atom->map_one(m,0);
      }
    }

    maxspecial = MAX(maxspecial,unique);

    if (onefive_flag) {
      unique15 = 0;
      for (j = 0; j < nspecial15[i]; j++) {
        m = onefive[i][j];
        if (atom->map(m) < 0) {
          unique15++;
          atom->map_one(m,0);
        }
      }
      maxspecial15 = MAX(maxspecial15,unique15);
    }

    atom->map_one(tag[i],-1);
    for (j = 0; j < nspecial[i][0]; j++) atom->map_one(onetwo[i][j],-1);
    for (j = 0; j < nspecial[i][1]; j++) atom->map_one(onethree[i][j],-1);
    for (j = 0; j < nspecial[i][2]; j++) atom->map_one(onefour[i][j],-1);
    if (onefive_flag)
      for (j = 0; j < nspecial15[i]; j++) atom->map_one(onefive[i][j],-1);
  }

  // if atom->maxspecial has been updated before,
  //   make certain it is not reset to a smaller value
  // since atom->maxspecial is initialized to 1,
  //   this ensures that it stays larger than zero

  maxspecial = MAX(atom->maxspecial,maxspecial);

  // compute global maxspecial
  // add in extra factor from special_bonds command
  // allocate correct special array with same nmax, new maxspecial
  // previously allocated one must be destroyed
  // must make AtomVec class update its ptr to special

  MPI_Allreduce(&maxspecial,&atom->maxspecial,1,MPI_INT,MPI_MAX,world);
  atom->maxspecial += force->special_extra;

  // add force->special_extra only once

  force->special_extra = 0;

  if (me == 0)
    utils::logmesg(lmp,"{:>6} = max # of special neighbors\n",atom->maxspecial);

  if (lmp->kokkos) {
    auto  atomKK = dynamic_cast<AtomKokkos*>(atom);
    atomKK->modified(Host,SPECIAL_MASK);
    atomKK->sync(Device,SPECIAL_MASK);
    auto  memoryKK = dynamic_cast<MemoryKokkos*>(memory);
    memoryKK->grow_kokkos(atomKK->k_special,atom->special,
                        atom->nmax,atom->maxspecial,"atom:special");
    atomKK->modified(Device,SPECIAL_MASK);
    atomKK->sync(Host,SPECIAL_MASK);
    atom->avec->grow_pointers();
  } else {
    memory->destroy(atom->special);
    memory->create(atom->special,atom->nmax,atom->maxspecial,"atom:special");
  }

  // if 1-5 is enabled, similarly compute global maxspecial15 and reallocate

  if (onefive_flag) {
    maxspecial15 = MAX(atom->maxspecial15,maxspecial15);
    MPI_Allreduce(&maxspecial15,&atom->maxspecial15,1,MPI_INT,MPI_MAX,world);
    memory->destroy(atom->special15);
    memory->create(atom->special15,atom->nmax,atom->maxspecial15,"atom:special15");
  }

  // ----------------------------------------------------
  // fill special array with 1-2, 1-3, 1-4 neighs for each atom
  // optionally fill special15 array with 1-5 neighs
  // ----------------------------------------------------

  tagint **special = atom->special;
  tagint **special15 = atom->special15;

  // again use map to cull duplicates
  // exclude original atom explicitly
  // adjust nspecial[i] values to reflect removed duplicates
  // nspecial[i][1] and nspecial[i][2] now become cumulative counters

  for (i = 0; i < nlocal; i++) {
    unique = 0;
    atom->map_one(tag[i],0);

    for (j = 0; j < nspecial[i][0]; j++) {
      m = onetwo[i][j];
      if (atom->map(m) < 0) {
        special[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][0] = unique;

    for (j = 0; j < nspecial[i][1]; j++) {
      m = onethree[i][j];
      if (atom->map(m) < 0) {
        special[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][1] = unique;

    for (j = 0; j < nspecial[i][2]; j++) {
      m = onefour[i][j];
      if (atom->map(m) < 0) {
        special[i][unique++] = m;
        atom->map_one(m,0);
      }
    }
    nspecial[i][2] = unique;

    if (onefive_flag) {
      unique15 = 0;
      for (j = 0; j < nspecial15[i]; j++) {
        m = onefive[i][j];
        if (atom->map(m) < 0) {
          special15[i][unique15++] = m;
          atom->map_one(m,0);
        }
      }
      nspecial15[i] = unique15;
    }

    atom->map_one(tag[i],-1);
    for (j = 0; j < nspecial[i][2]; j++) atom->map_one(special[i][j],-1);
    if (onefive_flag)
      for (j = 0; j < nspecial15[i]; j++) atom->map_one(special15[i][j],-1);
  }

  // re-create map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();
}

/* ----------------------------------------------------------------------
   trim list of 1-3 neighbors by checking defined angles
   delete a 1-3 neigh if they are not end atoms of a defined angle
     and if they are not 1,3 or 2,4 atoms of a defined dihedral
   uses rendezvous comm
------------------------------------------------------------------------- */

void Special::angle_trim()
{
  int i,j,k,m;

  int *num_angle = atom->num_angle;
  int *num_dihedral = atom->num_dihedral;
  tagint **angle_atom1 = atom->angle_atom1;
  tagint **angle_atom2 = atom->angle_atom2;
  tagint **angle_atom3 = atom->angle_atom3;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom3 = atom->dihedral_atom3;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // stats on old 1-3 neighbor counts

  double onethreecount = 0.0;
  for (i = 0; i < nlocal; i++) onethreecount += nspecial[i][1];
  double allcount;
  MPI_Allreduce(&onethreecount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0)
    utils::logmesg(lmp,"  {} = # of 1-3 neighbors before angle trim\n",
                   allcount);

  // if angles or dihedrals are defined
  // rendezvous angle 1-3 and dihedral 1-3,2-4 pairs

  if ((num_angle && atom->nangles) || (num_dihedral && atom->ndihedrals)) {

    // nsend = # of my datums to send
    // latter is only for angles or dihedrlas where I own atom2 (newton bond off)

    int nsend = 0;
    for (i = 0; i < nlocal; i++) {
      if (num_angle) {
        for (j = 0; j < num_angle[i]; j++) {
          if (tag[i] != angle_atom2[i][j]) continue;
          m = atom->map(angle_atom1[i][j]);
          if (m < 0 || m >= nlocal) nsend++;
          m = atom->map(angle_atom3[i][j]);
          if (m < 0 || m >= nlocal) nsend++;
        }
      }

      if (num_dihedral) {
        for (j = 0; j < num_dihedral[i]; j++) {
          if (tag[i] != dihedral_atom2[i][j]) continue;
          m = atom->map(dihedral_atom1[i][j]);
          if (m < 0 || m >= nlocal) nsend++;
          m = atom->map(dihedral_atom3[i][j]);
          if (m < 0 || m >= nlocal) nsend++;
          m = atom->map(dihedral_atom4[i][j]);
          if (m < 0 || m >= nlocal) nsend++;
        }
      }
    }

    int *proclist;
    memory->create(proclist,nsend,"special:proclist");
    auto inbuf = (PairRvous *) memory->smalloc((bigint) nsend*sizeof(PairRvous),"special:inbuf");

    // setup input buf to rendezvous comm
    // datums = pairs of onetwo partners where either is unknown
    //          these pairs are onethree neighbors
    //          datum = onetwo ID, onetwo ID
    // owning proc for each datum = onetwo ID % nprocs

    nsend = 0;
    for (i = 0; i < nlocal; i++) {
      if (num_angle) {
        for (j = 0; j < num_angle[i]; j++) {
          if (tag[i] != angle_atom2[i][j]) continue;

          m = atom->map(angle_atom1[i][j]);
          if (m < 0 || m >= nlocal) {
            proclist[nsend] = angle_atom1[i][j] % nprocs;
            inbuf[nsend].atomID = angle_atom1[i][j];
            inbuf[nsend].partnerID = angle_atom3[i][j];
            nsend++;
          }

          m = atom->map(angle_atom3[i][j]);
          if (m < 0 || m >= nlocal) {
            proclist[nsend] = angle_atom3[i][j] % nprocs;
            inbuf[nsend].atomID = angle_atom3[i][j];
            inbuf[nsend].partnerID = angle_atom1[i][j];
            nsend++;
          }
        }
      }

      if (num_dihedral) {
        for (j = 0; j < num_dihedral[i]; j++) {
          if (tag[i] != dihedral_atom2[i][j]) continue;

          m = atom->map(dihedral_atom1[i][j]);
          if (m < 0 || m >= nlocal) {
            proclist[nsend] = dihedral_atom1[i][j] % nprocs;
            inbuf[nsend].atomID = dihedral_atom1[i][j];
            inbuf[nsend].partnerID = dihedral_atom3[i][j];
            nsend++;
          }

          m = atom->map(dihedral_atom3[i][j]);
          if (m < 0 || m >= nlocal) {
            proclist[nsend] = dihedral_atom3[i][j] % nprocs;
            inbuf[nsend].atomID = dihedral_atom3[i][j];
            inbuf[nsend].partnerID = dihedral_atom1[i][j];
            nsend++;
          }

          m = atom->map(dihedral_atom4[i][j]);
          if (m < 0 || m >= nlocal) {
            proclist[nsend] = dihedral_atom4[i][j] % nprocs;
            inbuf[nsend].atomID = dihedral_atom4[i][j];
            inbuf[nsend].partnerID = dihedral_atom2[i][j];
            nsend++;
          }
        }
      }
    }

    // perform rendezvous operation

    char *buf;
    int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(PairRvous), 0,proclist,
                                   rendezvous_pairs,0,buf,sizeof(PairRvous), (void *) this);
    auto outbuf = (PairRvous *) buf;

    memory->destroy(proclist);
    memory->sfree(inbuf);

    // flag all onethree atoms to keep

    int max = 0;
    for (i = 0; i < nlocal; i++)
      max = MAX(max,nspecial[i][1]);
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

    int **flag;
    memory->create(flag,nlocal,maxall,"special:flag");

    for (i = 0; i < nlocal; i++)
      for (j = 0; j < nspecial[i][1]; j++)
        flag[i][j] = 0;

    // reset nspecial[1] and onethree for all owned atoms based on output info
    // based on owned info plus rendezvous output info
    // output datums = pairs of atoms that are 1-3 neighbors

    for (i = 0; i < nlocal; i++) {
      if (num_angle) {
        for (j = 0; j < num_angle[i]; j++) {
          if (tag[i] != angle_atom2[i][j]) continue;

          m = atom->map(angle_atom1[i][j]);
          if (m >= 0 && m < nlocal) {
            for (k = 0; k < nspecial[m][1]; k++)
              if (onethree[m][k] == angle_atom3[i][j]) {
                flag[m][k] = 1;
                break;
              }
          }

          m = atom->map(angle_atom3[i][j]);
          if (m >= 0 && m < nlocal) {
            for (k = 0; k < nspecial[m][1]; k++)
              if (onethree[m][k] == angle_atom1[i][j]) {
                flag[m][k] = 1;
                break;
              }
          }
        }
      }

      if (num_dihedral) {
        for (j = 0; j < num_dihedral[i]; j++) {
          if (tag[i] != dihedral_atom2[i][j]) continue;

          m = atom->map(dihedral_atom1[i][j]);
          if (m >= 0 && m < nlocal) {
            for (k = 0; k < nspecial[m][1]; k++)
              if (onethree[m][k] == dihedral_atom3[i][j]) {
                flag[m][k] = 1;
                break;
              }
          }

          m = atom->map(dihedral_atom3[i][j]);
          if (m >= 0 && m < nlocal) {
            for (k = 0; k < nspecial[m][1]; k++)
              if (onethree[m][k] == dihedral_atom1[i][j]) {
                flag[m][k] = 1;
                break;
              }
          }

          m = atom->map(dihedral_atom4[i][j]);
          if (m >= 0 && m < nlocal) {
            for (k = 0; k < nspecial[m][1]; k++)
              if (onethree[m][k] == dihedral_atom2[i][j]) {
                flag[m][k] = 1;
                break;
              }
          }
        }
      }
    }

    for (m = 0; m < nreturn; m++) {
      i = atom->map(outbuf[m].atomID);
      for (k = 0; k < nspecial[i][1]; k++)
        if (onethree[i][k] == outbuf[m].partnerID) {
          flag[i][k] = 1;
          break;
        }
    }

    memory->destroy(outbuf);

    // use flag values to compress onefour list for each atom

    for (i = 0; i < nlocal; i++) {
      j = 0;
      while (j < nspecial[i][1]) {
        if (flag[i][j] == 0) {
          onethree[i][j] = onethree[i][nspecial[i][1]-1];
          flag[i][j] = flag[i][nspecial[i][1]-1];
          nspecial[i][1]--;
        } else j++;
      }
    }

    memory->destroy(flag);

    // if no angles or dihedrals are defined, delete all 1-3 neighs

  } else {
    for (i = 0; i < nlocal; i++) nspecial[i][1] = 0;
  }

  // stats on new 1-3 neighbor counts

  onethreecount = 0.0;
  for (i = 0; i < nlocal; i++) onethreecount += nspecial[i][1];
  MPI_Allreduce(&onethreecount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0)
    utils::logmesg(lmp,"  {} = # of 1-3 neighbors after angle trim\n",
                   allcount);
}

/* ----------------------------------------------------------------------
   trim list of 1-4 neighbors by checking all defined dihedrals
   delete a 1-4 neigh if they are not end atoms of a defined dihedral
   uses rendezvous comm
------------------------------------------------------------------------- */

void Special::dihedral_trim()
{
  int i,j,k,m;

  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom1 = atom->dihedral_atom1;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  tagint **dihedral_atom4 = atom->dihedral_atom4;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // stats on old 1-4 neighbor counts

  double onefourcount = 0.0;
  for (i = 0; i < nlocal; i++) onefourcount += nspecial[i][2];
  double allcount;
  MPI_Allreduce(&onefourcount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0)
    utils::logmesg(lmp,"  {} = # of 1-4 neighbors before dihedral trim\n",
                   allcount);

  // if dihedrals are defined, rendezvous the dihedral 1-4 pairs

  if (num_dihedral && atom->ndihedrals) {

    // nsend = # of my datums to send

    int nsend = 0;
    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < num_dihedral[i]; j++) {
        if (tag[i] != dihedral_atom2[i][j]) continue;
        m = atom->map(dihedral_atom1[i][j]);
        if (m < 0 || m >= nlocal) nsend++;
        m = atom->map(dihedral_atom4[i][j]);
        if (m < 0 || m >= nlocal) nsend++;
      }
    }

    int *proclist;
    memory->create(proclist,nsend,"special:proclist");
    auto inbuf = (PairRvous *) memory->smalloc((bigint) nsend*sizeof(PairRvous),"special:inbuf");

    // setup input buf to rendezvous comm
    // datums = pairs of onefour atom IDs in a dihedral defined for my atoms
    //          only dihedrals where I own atom2 (in case newton_bond off)
    //          datum = atom1 ID and atom4 ID
    //          send the datum twice, to owner of atom1 ID and atom4 ID
    // owning procs for each datum = atom1 or atom4 ID % nprocs

    nsend = 0;
    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < num_dihedral[i]; j++) {
        if (tag[i] != dihedral_atom2[i][j]) continue;

        m = atom->map(dihedral_atom1[i][j]);
        if (m < 0 || m >= nlocal) {
          proclist[nsend] = dihedral_atom1[i][j] % nprocs;
          inbuf[nsend].atomID = dihedral_atom1[i][j];
          inbuf[nsend].partnerID = dihedral_atom4[i][j];
          nsend++;
        }

        m = atom->map(dihedral_atom4[i][j]);
        if (m < 0 || m >= nlocal) {
          proclist[nsend] = dihedral_atom4[i][j] % nprocs;
          inbuf[nsend].atomID = dihedral_atom4[i][j];
          inbuf[nsend].partnerID = dihedral_atom1[i][j];
          nsend++;
        }
      }
    }

    // perform rendezvous operation

    char *buf;
    int nreturn = comm->rendezvous(RVOUS,nsend,(char *) inbuf,sizeof(PairRvous), 0,proclist,
                                   rendezvous_pairs,0,buf,sizeof(PairRvous), (void *) this);
    auto outbuf = (PairRvous *) buf;

    memory->destroy(proclist);
    memory->sfree(inbuf);

    // flag all of my onefour IDs to keep

    int max = 0;
    for (i = 0; i < nlocal; i++)
      max = MAX(max,nspecial[i][2]);
    MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

    int **flag;
    memory->create(flag,nlocal,maxall,"special:flag");

    for (i = 0; i < nlocal; i++)
      for (j = 0; j < nspecial[i][2]; j++)
        flag[i][j] = 0;

    for (i = 0; i < nlocal; i++) {
      for (j = 0; j < num_dihedral[i]; j++) {
        if (tag[i] != dihedral_atom2[i][j]) continue;

        m = atom->map(dihedral_atom1[i][j]);
        if (m >= 0 && m < nlocal) {
          for (k = 0; k < nspecial[m][2]; k++)
            if (onefour[m][k] == dihedral_atom4[i][j]) {
              flag[m][k] = 1;
              break;
            }
        }

        m = atom->map(dihedral_atom4[i][j]);
        if (m >= 0 && m < nlocal) {
          for (k = 0; k < nspecial[m][2]; k++)
            if (onefour[m][k] == dihedral_atom1[i][j]) {
              flag[m][k] = 1;
              break;
            }
        }
      }
    }

    for (m = 0; m < nreturn; m++) {
      i = atom->map(outbuf[m].atomID);
      for (k = 0; k < nspecial[i][2]; k++)
        if (onefour[i][k] == outbuf[m].partnerID) {
          flag[i][k] = 1;
          break;
        }
    }

    memory->destroy(outbuf);

    // use flag values to compress onefour list for each atom

    for (i = 0; i < nlocal; i++) {
      j = 0;
      while (j < nspecial[i][2]) {
        if (flag[i][j] == 0) {
          onefour[i][j] = onefour[i][nspecial[i][2]-1];
          flag[i][j] = flag[i][nspecial[i][2]-1];
          nspecial[i][2]--;
        } else j++;
      }
    }

    memory->destroy(flag);

  // if no dihedrals are defined, delete all 1-4 neighs

  } else {
    for (i = 0; i < nlocal; i++) nspecial[i][2] = 0;
  }

  // stats on new 1-4 neighbor counts

  onefourcount = 0.0;
  for (i = 0; i < nlocal; i++) onefourcount += nspecial[i][2];
  MPI_Allreduce(&onefourcount,&allcount,1,MPI_DOUBLE,MPI_SUM,world);

  if (me == 0)
    utils::logmesg(lmp,"  {} = # of 1-4 neighbors after dihedral trim\n",
                   allcount);
}

/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N IDRvous datums
   no outbuf
------------------------------------------------------------------------- */

int Special::rendezvous_ids(int n, char *inbuf, int &flag, int *& /*proclist*/, char *& /*outbuf*/, void *ptr)
{
  auto sptr = (Special *) ptr;
  Memory *memory = sptr->memory;

  int *procowner;
  tagint *atomIDs;

  memory->create(procowner,n,"special:procowner");
  memory->create(atomIDs,n,"special:atomIDs");

  auto in = (IDRvous *) inbuf;

  for (int i = 0; i < n; i++) {
    procowner[i] = in[i].me;
    atomIDs[i] = in[i].atomID;
  }

  // store rendezvous data in Special class

  sptr->nrvous = n;
  sptr->procowner = procowner;
  sptr->atomIDs = atomIDs;

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}


/* ----------------------------------------------------------------------
   process data for atoms assigned to me in rendezvous decomposition
   inbuf = list of N PairRvous datums
   outbuf = same list of N PairRvous datums, routed to different procs
------------------------------------------------------------------------- */

int Special::rendezvous_pairs(int n, char *inbuf, int &flag, int *&proclist,
                              char *&outbuf, void *ptr)
{
  auto sptr = (Special *) ptr;
  Atom *atom = sptr->atom;
  Memory *memory = sptr->memory;

  // clear atom map so it can be used here as a hash table
  // faster than an STL map for large atom counts

  atom->map_clear();

  // hash atom IDs stored in rendezvous decomposition

  int nrvous = sptr->nrvous;
  tagint *atomIDs = sptr->atomIDs;

  for (int i = 0; i < nrvous; i++)
    atom->map_one(atomIDs[i],i);

  // proclist = owner of atomID in caller decomposition

  auto in = (PairRvous *) inbuf;
  int *procowner = sptr->procowner;
  memory->create(proclist,n,"special:proclist");

  int m;
  for (int i = 0; i < n; i++) {
    m = atom->map(in[i].atomID);
    proclist[i] = procowner[m];
  }

  outbuf = inbuf;

  // re-create atom map

  atom->map_init(0);
  atom->nghost = 0;
  atom->map_set();

  // flag = 1: outbuf = inbuf

  flag = 1;
  return n;
}

/* ----------------------------------------------------------------------
   allow fixes to alter special list
   currently, only fix drude does this
     so that both the Drude core and electron are same level of neighbor
------------------------------------------------------------------------- */

void Special::fix_alteration()
{
  for (const auto &ifix : modify->get_fix_list())
    if (ifix->special_alter_flag) ifix->rebuild_special();
}

/* ----------------------------------------------------------------------
   print timing output
------------------------------------------------------------------------- */

void Special::timer_output(double time1)
{
  if (comm->me == 0)
    utils::logmesg(lmp,"  special bonds CPU = {:.3f} seconds\n",
                   platform::walltime()-time1);
}
