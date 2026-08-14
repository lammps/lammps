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

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (Sandia, tnshan@sandia.gov)
------------------------------------------------------------------------- */

#include "fix_reaxff_bonds.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "graphics.h"
#include "memory.h"
#include "neigh_list.h"
#include "update.h"

#include "pair_reaxff.h"
#include "reaxff_api.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace ReaxFF;

/* ---------------------------------------------------------------------- */

FixReaxFFBonds::FixReaxFFBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), neighid(nullptr), abo(nullptr), fp(nullptr), lists(nullptr),
  reaxff(nullptr), list(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg != 5)
    error->all(FLERR, Error::NOPOINTER, "Fix reaxff/bonds expected 5 arguments but got {}", narg);

  nmax = atom->nmax;
  compressed = 0;
  multifile = 0;
  padflag = 0;
  first_flag = true;
  numobjs = 0;
  dynamic_group_allow = 1;     // applies only to FixReaxFFBonds::image()

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix reaxff/bonds nevery value {}", nevery);
  global_freq = nevery;

  filename = arg[4];
  if (filename.rfind('*') != std::string::npos) multifile = 1;
  if (platform::has_compress_extension(filename)) compressed = 1;

  if (atom->tag_consecutive() == 0)
    error->all(FLERR, Error::NOLASTLINE, "Fix reaxff/bonds requires consecutive atom-IDs");

  allocate();
}

/* ---------------------------------------------------------------------- */

FixReaxFFBonds::~FixReaxFFBonds()
{
  destroy();

  if (fp) fclose(fp);

  // clean up dump image data
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixReaxFFBonds::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::setup(int /*vflag*/)
{
  if (atom->natoms > MAXSMALLINT)
    error->all(FLERR, Error::NOLASTLINE, "Too many atoms for fix {}", style);

  // only print output during setup() at the very beginning
  // to avoid duplicate outputs when using multiple run statements
  if (first_flag) end_of_step();
  first_flag = false;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::init()
{
  reaxff = dynamic_cast<PairReaxFF *>(force->pair_match("^reax..",0));
  if (reaxff == nullptr)
    error->all(FLERR, Error::NOLASTLINE, "Cannot use fix reaxff/bonds without "
               "pair_style reaxff, reaxff/kk, or reaxff/omp");
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::end_of_step()
{
  Output_ReaxFF_Bonds();
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::Output_ReaxFF_Bonds()

{
  int i, j;
  int nbuf, nbuf_local;
  int nlocal_max, numbonds, numbonds_max;
  double *buf;

  int nlocal = atom->nlocal;
  int nlocal_tot = static_cast<int>(atom->natoms);

  if (atom->nmax > nmax) {
    destroy();
    nmax = atom->nmax;
    allocate();
  }

  for (i = 0; i < nmax; i++) {
    numneigh[i] = 0;
    for (j = 0; j < MAXREAXBOND; j++) {
      neighid[i][j] = 0;
      abo[i][j] = 0.0;
    }
  }
  numbonds = FindBond();

  // no file output with NULL file name.
  if (filename == "NULL") return;

  // allocate a temporary buffer for the snapshot info
  MPI_Allreduce(&numbonds,&numbonds_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nlocal,&nlocal_max,1,MPI_INT,MPI_MAX,world);

  nbuf = 1+(numbonds_max*2+10)*nlocal_max;
  memory->create(buf,nbuf,"reaxff/bonds:buf");
  for (i = 0; i < nbuf; i ++) buf[i] = 0.0;

  // Pass information to buffer
  PassBuffer(buf, nbuf_local);

  // Receive information from buffer for output
  RecvBuffer(buf, nbuf, nbuf_local, nlocal_tot, numbonds_max);

  memory->destroy(buf);
}

/* ---------------------------------------------------------------------- */

int FixReaxFFBonds::FindBond()
{
  int *ilist, i, ii, inum;
  int j, pj, nj;
  tagint jtag;
  double bo_tmp,bo_cut;

  inum = reaxff->list->inum;
  ilist = reaxff->list->ilist;
  bond_data *bo_ij;
  bo_cut = reaxff->api->control->bg_cut;

  tagint *tag = atom->tag;
  int numbonds = 0;
  numobjs = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    nj = 0;

    for (pj = Start_Index(i, reaxff->api->lists); pj < End_Index(i, reaxff->api->lists); ++pj) {
      bo_ij = &(reaxff->api->lists->select.bond_list[pj]);
      j = bo_ij->nbr;
      jtag = tag[j];
      bo_tmp = bo_ij->bo_data.BO;

      if (bo_tmp > bo_cut) {
        neighid[i][nj] = jtag;
        abo[i][nj] = bo_tmp;
        nj ++;
      }
    }
    numneigh[i] = nj;
    numobjs += nj;
    if (nj > numbonds) numbonds = nj;
  }
  return numbonds;
}
/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::PassBuffer(double *buf, int &nbuf_local)
{
  int i, j, k, numbonds;
  int nlocal = atom->nlocal;

  j = 2;
  buf[0] = nlocal;
  for (i = 0; i < nlocal; i++) {
    buf[j-1] = atom->tag[i];
    buf[j+0] = atom->type[i];
    buf[j+1] = reaxff->api->workspace->total_bond_order[i];
    buf[j+2] = reaxff->api->workspace->nlp[i];
    buf[j+3] = atom->q[i];
    buf[j+4] = numneigh[i];
    numbonds = std::lround(buf[j+4]);

    for (k = 5; k < 5+numbonds; k++) {
      buf[j+k] = neighid[i][k-5];
    }
    j += (5+numbonds);

    if (atom->molecule == nullptr) buf[j] = 0.0;
    else buf[j] = atom->molecule[i];
    j ++;

    for (k = 0; k < numbonds; k++) {
      buf[j+k] = abo[i][k];
    }
    j += (1+numbonds);
  }
  nbuf_local = j - 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::RecvBuffer(double *buf, int nbuf, int nbuf_local,
                               int natoms, int maxnum)
{
  int i, j, k, itype;
  int inode, nlocal_tmp, numbonds;
  tagint itag;
  int nlocal = atom->nlocal;
  bigint ntimestep = update->ntimestep;
  double sbotmp, nlptmp, avqtmp;

  double cutof3 = reaxff->api->control->bg_cut;
  MPI_Request irequest, irequest2;

  if (comm->me == 0) {
    std::string myfile = filename;
    if (multifile) myfile = utils::star_subst(myfile, update->ntimestep, padflag);
    if (multifile && fp) {
      fclose(fp);
      fp = nullptr;
    }
    if (!fp) {
      if (compressed) {
        fp = platform::compressed_write(myfile);
      } else {
        fp = fopen(myfile.c_str(), "w");
      }
      if (!fp)
        error->one(FLERR, Error::NOLASTLINE,
                   "Cannot open fix reaxff/bonds file {}: {}", myfile, utils::getsyserror());
    }

    utils::print(fp,"# Timestep {}\n#\n",ntimestep);
    utils::print(fp,"# Number of particles {}\n#\n",natoms);
    utils::print(fp,"# Max number of bonds per atom {} with coarse bond order cutoff {:5.3f}\n",
               maxnum,cutof3);
    utils::print(fp,"# Particle connection table and bond orders\n"
               "# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q\n");
  }

  j = 2;
  if (comm->me == 0) {
    for (inode = 0; inode < comm->nprocs; inode ++) {
      if (inode == 0) {
        nlocal_tmp = nlocal;
      } else {
        MPI_Irecv(&buf[0],nbuf,MPI_DOUBLE,inode,0,world,&irequest);
        MPI_Wait(&irequest,MPI_STATUS_IGNORE);
        nlocal_tmp = std::lround(buf[0]);
      }
      j = 2;
      for (i = 0; i < nlocal_tmp; i ++) {
        itag = static_cast<tagint> (buf[j-1]);
        itype = std::lround(buf[j+0]);
        sbotmp = buf[j+1];
        nlptmp = buf[j+2];
        avqtmp = buf[j+3];
        numbonds = std::lround(buf[j+4]);

        auto mesg = fmt::format(" {} {} {}",itag,itype,numbonds);
        for (k = 5; k < 5+numbonds; k++)
          mesg += " " + std::to_string(static_cast<tagint> (buf[j+k]));
        j += (5+numbonds);

        mesg += " " + std::to_string(static_cast<tagint> (buf[j]));
        j ++;

        for (k = 0; k < numbonds; k++) mesg += fmt::format("{:14.3f}",buf[j+k]);
        j += (1+numbonds);

        mesg += fmt::format("{:14.3f}{:14.3f}{:14.3f}\n",sbotmp,nlptmp,avqtmp);
        utils::print(fp, mesg);
      }
    }
  } else {
    MPI_Isend(&buf[0],nbuf_local,MPI_DOUBLE,0,0,world,&irequest2);
    MPI_Wait(&irequest2,MPI_STATUS_IGNORE);
  }
  if (fp) {
    fputs("# \n",fp);
    fflush(fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::destroy()
{
  memory->destroy(abo);
  memory->destroy(neighid);
  memory->destroy(numneigh);
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBonds::allocate()
{
  memory->create(abo,nmax,MAXREAXBOND,"reaxff/bonds:abo");
  memory->create(neighid,nmax,MAXREAXBOND,"reaxff/bonds:neighid");
  memory->create(numneigh,nmax,"reaxff/bonds:numneigh");
}

/* ---------------------------------------------------------------------- */

int FixReaxFFBonds::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0], "pad") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "fix_modify pad", error);
    padflag = utils::inumeric(FLERR, arg[1], false, lmp);
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

double FixReaxFFBonds::memory_usage()
{
  double bytes;

  bytes = 3.0*nmax*sizeof(double);
  bytes += (double)nmax*sizeof(int);
  bytes += (double)1.0*nmax*MAXREAXBOND*sizeof(double);
  bytes += (double)1.0*nmax*MAXREAXBOND*sizeof(int);

  return bytes;
}

/* ---------------------------------------------------------------------- */

int FixReaxFFBonds::image(int *&objs, double **&parms)
{
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, Error::NOLASTLINE,
               "Using fix reaxff/bonds with dump image requires an atom map");

  if (!numobjs) return 0;
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
  memory->create(imgobjs, numobjs, "reaxff/bonds:imgobjs");
  memory->create(imgparms, numobjs, 8, "reaxff/bonds:imgparms");

  const int *type = atom->type;
  const int *mask = atom->mask;
  const double * const * const x = atom->x;

  int inum = reaxff->list->inum;
  int *ilist = reaxff->list->ilist;

  int n = 0;
  for (int ii = 0; ii < inum; ++ii) {
    int i = ilist[ii];
    if (mask[i] & groupbit) {
      for (int jj = 0; jj < numneigh[i]; ++jj) {
        int j = atom->map(neighid[i][jj]);
        j = domain->closest_image(i,j);
        if (j < 0) continue;
        if (mask[j] & groupbit) {
          imgobjs[n] = Graphics::BOND;
          imgparms[n][0] = type[i];
          imgparms[n][1] = type[j];
          imgparms[n][2] = x[i][0];
          imgparms[n][3] = x[i][1];
          imgparms[n][4] = x[i][2];
          imgparms[n][5] = x[j][0];
          imgparms[n][6] = x[j][1];
          imgparms[n][7] = x[j][2];
          ++n;
        }
      }
    }
  }

  objs = imgobjs;
  parms = imgparms;
  return n;
}
