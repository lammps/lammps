/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Christoph Scherer (MPIP Mainz)
   scherer@mpip-mainz.mpg.de
------------------------------------------------------------------------- */

#include "pair_threebody_table.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"
#include "table_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

static constexpr int DELTA = 4;

/* ---------------------------------------------------------------------- */

PairThreebodyTable::PairThreebodyTable(LAMMPS *lmp) :
    Pair(lmp), params(nullptr), neighshort(nullptr)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  maxshort = 10;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairThreebodyTable::~PairThreebodyTable()
{
  if (copymode) return;

  for (int m = 0; m < nparams; m++) free_param(&params[m]);    // free_param will call free_table
  memory->sfree(params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
  }
}

/* ---------------------------------------------------------------------- */

void PairThreebodyTable::compute(int eflag, int vflag)
{
  int i, j, k, ii, jj, kk, inum, jnum, jnumm1;
  int itype, jtype, ktype, ijparam, ijkparam;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl;
  double rsq, rsq1, rsq2;
  double delr1[3], delr2[3], fi[3], fj[3], fk[3];
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp, fytmp, fztmp;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (rsq >= params[ijparam].cutsq) {
        continue;
      } else {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort / 2;
          memory->grow(neighshort, maxshort, "pair:neighshort");
        }
      }

      jtag = tag[j];
      if (itag > jtag) {
        if ((itag + jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag + jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      //two-body interactions are not computed
    }

    jnumm1 = numshort - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];

      double fjxtmp, fjytmp, fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj + 1; kk < numshort; kk++) {
        k = neighshort[kk];
        ktype = map[type[k]];
        ijkparam = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];

        threebody(&params[ijkparam], rsq1, rsq2, delr1, delr2, fi, fj, fk, eflag, evdwl);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i, j, k, evdwl, 0.0, fj, fk, delr1, delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairThreebodyTable::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(neighshort, maxshort, "pair:neighshort");
  map = new int[np1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairThreebodyTable::settings(int narg, char ** /*arg*/)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairThreebodyTable::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  map_element2type(narg - 3, arg + 3);

  // read potential file and initialize potential parameters

  if (params) {
    for (int m = 0; m < nparams; m++) free_param(&params[m]);    // free_param will call free_table
    memory->sfree(params);
    params = nullptr;
  }
  read_file(arg[2]);
  setup_params();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairThreebodyTable::init_style()
{
  if (atom->tag_enable == 0) error->all(FLERR, "Pair style threebody/table requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style threebody/table requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairThreebodyTable::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairThreebodyTable::read_file(char *file)
{
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "threebody", unit_convert_flag);
    char *line;

    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement, kelement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;
        for (kelement = 0; kelement < nelements; kelement++)
          if (kname == elements[kelement]) break;
        if (kelement == nelements) continue;

        // load up parameter settings and error check their values

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params, maxparam * sizeof(Param), "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA * sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        // if jelement = kelement, symmetric is true, if not then it is false
        params[nparams].symmetric = false;
        if (params[nparams].jelement == params[nparams].kelement) params[nparams].symmetric = true;
        // read cut off
        params[nparams].cut = values.next_double();

        // read parameters of angle table
        std::string name = values.next_string();
        params[nparams].tablenamelength = name.length() + 1;
        params[nparams].tablename = utils::strdup(name);

        name = values.next_string();
        params[nparams].keywordlength = name.length() + 1;
        params[nparams].keyword = utils::strdup(name);

        name = values.next_string();
        if (name != "linear") error->all(FLERR, "Unknown table style {} in threebody table", name);
        params[nparams].tablength = values.next_int();

      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      if (params[nparams].cut < 0.0 || params[nparams].tablength < 0.0)
        error->one(FLERR, "Illegal threebody/table parameters");

      nparams++;
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0)
    params = (Param *) memory->srealloc(params, maxparam * sizeof(Param), "pair:params");
  MPI_Bcast(params, maxparam * sizeof(Param), MPI_BYTE, 0, world);

  // for each set of parameters, broadcast table name and keyword and read threebody table
  for (int m = 0; m < nparams; ++m) {
    if (comm->me != 0) {
      params[m].tablename = new char[params[m].tablenamelength];
      params[m].keyword = new char[params[m].keywordlength];
    }

    MPI_Bcast(params[m].tablename, params[m].tablenamelength, MPI_CHAR, 0, world);
    MPI_Bcast(params[m].keyword, params[m].keywordlength, MPI_CHAR, 0, world);

    // initialize threebodytable
    memory->create(params[m].mltable, 1, "param:threebodytable");
    null_table(params[m].mltable);

    //call read_table to read corresponding tabulated threebody file (only called by process 0)
    if (comm->me == 0) {
      read_table(params[m].mltable, params[m].tablename, params[m].keyword, params[m].symmetric);
    }

    // broadcast read in threebodytable to all processes
    bcast_table(params[m].mltable, params[m].symmetric);

    // error check on table parameters
    if (params[m].mltable->ninput <= 1) error->one(FLERR, "Invalid threebody table length");
  }
}

/* ---------------------------------------------------------------------- */

void PairThreebodyTable::setup_params()
{
  int i, j, k, m, n;
  double rtmp;

  // set elem3param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem3param);
  memory->create(elem3param, nelements, nelements, nelements, "pair:elem3param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement && k == params[m].kelement) {
            if (n >= 0)
              error->all(FLERR, "Potential file has a duplicate entry for: {} {} {}", elements[i],
                         elements[j], elements[k]);
            n = m;
          }
        }
        if (n < 0)
          error->all(FLERR, "Potential file is missing an entry for: {} {} {}", elements[i],
                     elements[j], elements[k]);
        elem3param[i][j][k] = n;
      }

  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  // (cut = a)

  for (m = 0; m < nparams; m++) {
    rtmp = params[m].cut;
    params[m].cutsq = rtmp * rtmp;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void PairThreebodyTable::read_table(Table *tb, char *file, char *keyword, bool symmetric)
{
  TableFileReader reader(lmp, file, "threebodytable");

  char *line = reader.find_section_start(keyword);

  if (!line) error->one(FLERR, "Did not find keyword in table file");

  // read args on 2nd line of section
  // allocate table arrays for file values

  line = reader.next_line();
  param_extract(tb, line);

  // if it is a symmetric threebody interaction, less table entries are required
  if (symmetric) {
    memory->create(tb->r12file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:r12file");
    memory->create(tb->r13file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:r13file");
    memory->create(tb->thetafile, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:thetafile");
    memory->create(tb->f11file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f11file");
    memory->create(tb->f12file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f12file");
    memory->create(tb->f21file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f21file");
    memory->create(tb->f22file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f22file");
    memory->create(tb->f31file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f31file");
    memory->create(tb->f32file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f32file");
    memory->create(tb->efile, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:efile");
  }
  // else, more (full) table entries are required
  else {
    memory->create(tb->r12file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:r12file");
    memory->create(tb->r13file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:r13file");
    memory->create(tb->thetafile, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:thetafile");
    memory->create(tb->f11file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f11file");
    memory->create(tb->f12file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f12file");
    memory->create(tb->f21file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f21file");
    memory->create(tb->f22file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f22file");
    memory->create(tb->f31file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f31file");
    memory->create(tb->f32file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f32file");
    memory->create(tb->efile, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:efile");
  }

  // read threebody table values from file

  int cerror = 0;
  reader.skip_line();
  // if it is a symmetric threebody interaction, less table entries are required
  if (symmetric) {
    for (int i = 0; i < tb->ninput * tb->ninput * (tb->ninput + 1); i++) {
      line = reader.next_line(11);
      try {
        ValueTokenizer values(line);
        values.next_int();
        tb->r12file[i] = values.next_double();
        tb->r13file[i] = values.next_double();
        tb->thetafile[i] = values.next_double();
        tb->f11file[i] = values.next_double();
        tb->f12file[i] = values.next_double();
        tb->f21file[i] = values.next_double();
        tb->f22file[i] = values.next_double();
        tb->f31file[i] = values.next_double();
        tb->f32file[i] = values.next_double();
        tb->efile[i] = values.next_double();
      } catch (TokenizerException &) {
        ++cerror;
      }
    }
  } else {
    for (int i = 0; i < 2 * tb->ninput * tb->ninput * tb->ninput; i++) {
      line = reader.next_line(11);
      try {
        ValueTokenizer values(line);
        values.next_int();
        tb->r12file[i] = values.next_double();
        tb->r13file[i] = values.next_double();
        tb->thetafile[i] = values.next_double();
        tb->f11file[i] = values.next_double();
        tb->f12file[i] = values.next_double();
        tb->f21file[i] = values.next_double();
        tb->f22file[i] = values.next_double();
        tb->f31file[i] = values.next_double();
        tb->f32file[i] = values.next_double();
        tb->efile[i] = values.next_double();
      } catch (TokenizerException &) {
        ++cerror;
      }
    }
  }

  // warn if data was read incompletely, e.g. columns were missing

  if (cerror)
    error->warning(FLERR, "{} of {} lines in table incomplete or could not be parsed", cerror,
                   tb->ninput);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ theta0
   N is required, other params are optional

   only called by read_table, only called by proc 0
------------------------------------------------------------------------- */

void PairThreebodyTable::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->rmin = 0.0;
  tb->rmax = 0.0;

  try {
    ValueTokenizer values(line);

    while (values.has_next()) {
      std::string word = values.next_string();

      if (word == "N") {
        tb->ninput = values.next_int();
      } else if (word == "rmin") {
        tb->rmin = values.next_double();
      } else if (word == "rmax") {
        tb->rmax = values.next_double();
      } else {
        error->one(FLERR, "Invalid keyword {} in angle table parameters", word);
      }
    }
  } catch (TokenizerException &e) {
    error->one(FLERR, e.what());
  }

  if (tb->ninput == 0) error->one(FLERR, "threebodytable parameters did not set N");
  if (tb->rmin == 0.0) error->one(FLERR, "threebodytable parameters did not set rmin");
  if (tb->rmax == 0.0) error->one(FLERR, "threebodytable parameters did not set rmax");
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,afile,efile,ffile,fpflag,fplo,fphi,theta0
------------------------------------------------------------------------- */

void PairThreebodyTable::bcast_table(Table *tb, bool symmetric)
{
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);

  int me;
  MPI_Comm_rank(world, &me);
  if (me > 0) {
    // if it is a symmetric threebody interaction, less table entries are required
    if (symmetric) {
      memory->create(tb->r12file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:r12file");
      memory->create(tb->r13file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:r13file");
      memory->create(tb->thetafile, tb->ninput * tb->ninput * (tb->ninput + 1),
                     "mltable:thetafile");
      memory->create(tb->f11file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f11file");
      memory->create(tb->f12file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f12file");
      memory->create(tb->f21file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f21file");
      memory->create(tb->f22file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f22file");
      memory->create(tb->f31file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f31file");
      memory->create(tb->f32file, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:f32file");
      memory->create(tb->efile, tb->ninput * tb->ninput * (tb->ninput + 1), "mltable:efile");
    }
    // else, more (full) table entries are required
    else {
      memory->create(tb->r12file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:r12file");
      memory->create(tb->r13file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:r13file");
      memory->create(tb->thetafile, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:thetafile");
      memory->create(tb->f11file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f11file");
      memory->create(tb->f12file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f12file");
      memory->create(tb->f21file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f21file");
      memory->create(tb->f22file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f22file");
      memory->create(tb->f31file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f31file");
      memory->create(tb->f32file, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:f32file");
      memory->create(tb->efile, 2 * tb->ninput * tb->ninput * tb->ninput, "mltable:efile");
    }
  }

  // if it is a symmetric threebody interaction, less table entries are required
  if (symmetric) {
    MPI_Bcast(tb->r12file, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->r13file, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->thetafile, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f11file, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f12file, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f21file, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f22file, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f31file, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f32file, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->efile, tb->ninput * tb->ninput * (tb->ninput + 1), MPI_DOUBLE, 0, world);
  }
  // else, more (full) table entries are required
  else {
    MPI_Bcast(tb->r12file, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->r13file, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->thetafile, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f11file, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f12file, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f21file, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f22file, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f31file, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->f32file, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
    MPI_Bcast(tb->efile, 2 * tb->ninput * tb->ninput * tb->ninput, MPI_DOUBLE, 0, world);
  }

  MPI_Bcast(&tb->rmin, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&tb->rmax, 1, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairThreebodyTable::free_param(Param *pm)
{
  // call free_table to destroy associated threebodytables
  free_table(pm->mltable);
  // then destroy associated threebodytable
  delete[] pm->tablename;
  delete[] pm->keyword;
  memory->sfree(pm->mltable);
}

/* ---------------------------------------------------------------------- */

void PairThreebodyTable::free_table(Table *tb)
{
  memory->destroy(tb->r12file);
  memory->destroy(tb->r13file);
  memory->destroy(tb->thetafile);
  memory->destroy(tb->f11file);
  memory->destroy(tb->f12file);
  memory->destroy(tb->f21file);
  memory->destroy(tb->f22file);
  memory->destroy(tb->f31file);
  memory->destroy(tb->f32file);
  memory->destroy(tb->efile);
}

/* ---------------------------------------------------------------------- */

void PairThreebodyTable::null_table(Table *tb)
{
  tb->r12file = tb->r13file = tb->thetafile = nullptr;
  tb->f11file = tb->f12file = nullptr;
  tb->f21file = tb->f22file = nullptr;
  tb->f31file = tb->f32file = nullptr;
  tb->efile = nullptr;
}

/* ----------------------------------------------------------------------
   calculate potential u and force f at angle x
------------------------------------------------------------------------- */

void PairThreebodyTable::uf_lookup(Param *pm, double r12, double r13, double theta, double &f11,
                                   double &f12, double &f21, double &f22, double &f31, double &f32,
                                   double &u)
{
  int i, itable, nr12, nr13, ntheta;
  double dr, dtheta;
  dr = (pm->mltable->rmax - pm->mltable->rmin) / (pm->mltable->ninput - 1);
  dtheta = (180.0 - 0.0) / (pm->mltable->ninput * 2);

  //lookup scheme

  // if it is a symmetric threebody interaction, less table entries are required
  if (pm->symmetric) {
    nr12 = (r12 - pm->mltable->rmin + 0.5 * dr - 0.00000001) / dr;
    if (r12 == (pm->mltable->rmin - 0.5 * dr)) nr12 = 0;
    nr13 = (r13 - pm->mltable->rmin + 0.5 * dr - 0.00000001) / dr;
    if (r13 == (pm->mltable->rmin - 0.5 * dr)) nr13 = 0;
    nr13 -= nr12;
    ntheta = (theta - 0.00000001) / dtheta;
    if (theta >= 180.0) ntheta = (pm->mltable->ninput * 2) - 1;
    itable = 0;
    for (i = 0; i < nr12; i++) itable += (pm->mltable->ninput - i);
    itable += nr13;
    itable *= (pm->mltable->ninput * 2);
    itable += ntheta;
  } else {
    // else, more (full) table entries are required
    nr12 = (r12 - pm->mltable->rmin + 0.5 * dr - 0.00000001) / dr;
    if (r12 == (pm->mltable->rmin - 0.5 * dr)) nr12 = 0;
    nr13 = (r13 - pm->mltable->rmin + 0.5 * dr - 0.00000001) / dr;
    if (r13 == (pm->mltable->rmin - 0.5 * dr)) nr13 = 0;
    ntheta = (theta - 0.00000001) / dtheta;
    if (theta >= 180.0) ntheta = (pm->mltable->ninput * 2) - 1;
    itable = nr12 * (pm->mltable->ninput);
    itable += nr13;
    itable *= (pm->mltable->ninput * 2);
    itable += ntheta;
  }

  f11 = pm->mltable->f11file[itable];
  f12 = pm->mltable->f12file[itable];
  f21 = pm->mltable->f21file[itable];
  f22 = pm->mltable->f22file[itable];
  f31 = pm->mltable->f31file[itable];
  f32 = pm->mltable->f32file[itable];
  u = pm->mltable->efile[itable];
}

/* ---------------------------------------------------------------------- */

void PairThreebodyTable::threebody(Param *paramijk, double rsq1, double rsq2, double *delr1,
                                   double *delr2, double *fi, double *fj, double *fk, int eflag,
                                   double &eng)
{
  double r12, r13, theta, rinv, cs;

  double f11, f12, f21, f22, f31, f32, u, temp;
  bool swapped;
  double dr;
  dr = (paramijk->mltable->rmax - paramijk->mltable->rmin) / (paramijk->mltable->ninput - 1);

  //if swap indices or not
  swapped = false;

  r12 = sqrt(rsq1);
  r13 = sqrt(rsq2);
  rinv = 1.0 / (r12 * r13);
  cs = (delr1[0] * delr2[0] + delr1[1] * delr2[1] + delr1[2] * delr2[2]) * rinv;
  //compute angle between r12 and r13 in degrees
  theta = acos(cs) * 180.0 / MY_PI;

  //if r12 > r13 swap them, as in lookup table always r13 > r12 do to symmetry reasons
  if (r12 > r13) {
    temp = r12;
    r12 = r13;
    r13 = temp;
    swapped = true;
  }

  //look up forces and energy in table belonging to parameter set paramijk

  //only do lookup and add three-body interactions if r12 and r13 are both between rmin and rmax

  if ((r12 >= (paramijk->mltable->rmin - 0.5 * dr)) &&
      (r13 <= (paramijk->mltable->rmax + 0.5 * dr)) &&
      (r13 >= (paramijk->mltable->rmin - 0.5 * dr)) &&
      (r13 <= (paramijk->mltable->rmax + 0.5 * dr))) {
    uf_lookup(paramijk, r12, r13, theta, f11, f12, f21, f22, f31, f32, u);
  } else {
    f11 = f12 = f21 = f22 = f31 = f32 = u = 0.0;
  }

  // if the indices have been swapped, swap them back
  if (swapped) {
    temp = r12;
    r12 = r13;
    r13 = temp;
    temp = f11;
    f11 = f12;
    f12 = temp;
    temp = f21;
    f21 = f31;
    f31 = temp;
    temp = f22;
    f22 = -f32;
    f32 = -temp;
  }

  fi[0] = delr1[0] * f11 + delr2[0] * f12;
  fi[1] = delr1[1] * f11 + delr2[1] * f12;
  fi[2] = delr1[2] * f11 + delr2[2] * f12;

  fj[0] = delr1[0] * f21 + (delr2[0] - delr1[0]) * f22;
  fj[1] = delr1[1] * f21 + (delr2[1] - delr1[1]) * f22;
  fj[2] = delr1[2] * f21 + (delr2[2] - delr1[2]) * f22;

  fk[0] = delr2[0] * f31 + (delr2[0] - delr1[0]) * f32;
  fk[1] = delr2[1] * f31 + (delr2[1] - delr1[1]) * f32;
  fk[2] = delr2[2] * f31 + (delr2[2] - delr1[2]) * f32;

  if (eflag) eng = u;
}
