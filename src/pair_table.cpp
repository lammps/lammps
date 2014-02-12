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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "pair_table.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ,BMP};

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairTable::PairTable(LAMMPS *lmp) : Pair(lmp)
{
  ntables = 0;
  tables = NULL;
}

/* ---------------------------------------------------------------------- */

PairTable::~PairTable()
{
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(tabindex);
  }
}

/* ---------------------------------------------------------------------- */

void PairTable::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,factor_lj,fraction,value,a,b;
  int *ilist,*jlist,*numneigh,**firstneigh;
  Table *tb;

  union_int_float_t rsq_lookup;
  int tlm1 = tablength - 1;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        tb = &tables[tabindex[itype][jtype]];
        if (rsq < tb->innersq)
          error->one(FLERR,"Pair distance < table inner cutoff");

        if (tabstyle == LOOKUP) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          fpair = factor_lj * tb->f[itable];
        } else if (tabstyle == LINEAR) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          fraction = (rsq - tb->rsq[itable]) * tb->invdelta;
          value = tb->f[itable] + fraction*tb->df[itable];
          fpair = factor_lj * value;
        } else if (tabstyle == SPLINE) {
          itable = static_cast<int> ((rsq - tb->innersq) * tb->invdelta);
          if (itable >= tlm1)
            error->one(FLERR,"Pair distance > table outer cutoff");
          b = (rsq - tb->rsq[itable]) * tb->invdelta;
          a = 1.0 - b;
          value = a * tb->f[itable] + b * tb->f[itable+1] +
            ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
            tb->deltasq6;
          fpair = factor_lj * value;
        } else {
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & tb->nmask;
          itable >>= tb->nshiftbits;
          fraction = (rsq_lookup.f - tb->rsq[itable]) * tb->drsq[itable];
          value = tb->f[itable] + fraction*tb->df[itable];
          fpair = factor_lj * value;
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (tabstyle == LOOKUP)
            evdwl = tb->e[itable];
          else if (tabstyle == LINEAR || tabstyle == BITMAP)
            evdwl = tb->e[itable] + fraction*tb->de[itable];
          else
            evdwl = a * tb->e[itable] + b * tb->e[itable+1] +
              ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) *
              tb->deltasq6;
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairTable::allocate()
{
  allocated = 1;
  const int nt = atom->ntypes + 1;

  memory->create(setflag,nt,nt,"pair:setflag");
  memory->create(cutsq,nt,nt,"pair:cutsq");
  memory->create(tabindex,nt,nt,"pair:tabindex");

  memset(&setflag[0][0],0,nt*nt*sizeof(int));
  memset(&cutsq[0][0],0,nt*nt*sizeof(double));
  memset(&tabindex[0][0],0,nt*nt*sizeof(int));
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTable::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // new settings

  if (strcmp(arg[0],"lookup") == 0) tabstyle = LOOKUP;
  else if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
  else if (strcmp(arg[0],"bitmap") == 0) tabstyle = BITMAP;
  else error->all(FLERR,"Unknown table style in pair_style command");

  tablength = force->inumeric(FLERR,arg[1]);
  if (tablength < 2) error->all(FLERR,"Illegal number of pair table entries");

  // optional keywords
  // assert the tabulation is compatible with a specific long-range solver

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ewald") == 0) ewaldflag = 1;
    else if (strcmp(arg[iarg],"pppm") == 0) pppmflag = 1;
    else if (strcmp(arg[iarg],"msm") == 0) msmflag = 1;
    else if (strcmp(arg[iarg],"dispersion") == 0) dispersionflag = 1;
    else if (strcmp(arg[iarg],"tip4p") == 0) tip4pflag = 1;
    else error->all(FLERR,"Illegal pair_style command");
    iarg++;
  }

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(tabindex);
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTable::coeff(int narg, char **arg)
{
  if (narg != 4 && narg != 5) error->all(FLERR,"Illegal pair_coeff command");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int me;
  MPI_Comm_rank(world,&me);
  tables = (Table *)
    memory->srealloc(tables,(ntables+1)*sizeof(Table),"pair:tables");
  Table *tb = &tables[ntables];
  null_table(tb);
  if (me == 0) read_table(tb,arg[2],arg[3]);
  bcast_table(tb);

  // set table cutoff

  if (narg == 5) tb->cut = force->numeric(FLERR,arg[4]);
  else if (tb->rflag) tb->cut = tb->rhi;
  else tb->cut = tb->rfile[tb->ninput-1];

  // error check on table parameters
  // insure cutoff is within table
  // for BITMAP tables, file values can be in non-ascending order

  if (tb->ninput <= 1) error->one(FLERR,"Invalid pair table length");
  double rlo,rhi;
  if (tb->rflag == 0) {
    rlo = tb->rfile[0];
    rhi = tb->rfile[tb->ninput-1];
  } else {
    rlo = tb->rlo;
    rhi = tb->rhi;
  }
  if (tb->cut <= rlo || tb->cut > rhi)
    error->all(FLERR,"Invalid pair table cutoff");
  if (rlo <= 0.0) error->all(FLERR,"Invalid pair table cutoff");

  // match = 1 if don't need to spline read-in tables
  // this is only the case if r values needed by final tables
  //   exactly match r values read from file
  // for tabstyle SPLINE, always need to build spline tables

  tb->match = 0;
  if (tabstyle == LINEAR && tb->ninput == tablength &&
      tb->rflag == RSQ && tb->rhi == tb->cut) tb->match = 1;
  if (tabstyle == BITMAP && tb->ninput == 1 << tablength &&
      tb->rflag == BMP && tb->rhi == tb->cut) tb->match = 1;
  if (tb->rflag == BMP && tb->match == 0)
    error->all(FLERR,"Bitmapped table in file does not match requested table");

  // spline read-in values and compute r,e,f vectors within table

  if (tb->match == 0) spline_table(tb);
  compute_table(tb);

  // store ptr to table in tabindex

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      tabindex[i][j] = ntables;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Illegal pair_coeff command");
  ntables++;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTable::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  tabindex[j][i] = tabindex[i][j];

  return tables[tabindex[i][j]].cut;
}

/* ----------------------------------------------------------------------
   read a table section from a tabulated potential file
   only called by proc 0
   this function sets these values in Table:
     ninput,rfile,efile,ffile,rflag,rlo,rhi,fpflag,fplo,fphi,ntablebits
------------------------------------------------------------------------- */

void PairTable::read_table(Table *tb, char *file, char *keyword)
{
  char line[MAXLINE];

  // open file

  FILE *fp = open_potential(file);
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }

  // loop until section found with matching keyword

  while (1) {
    if (fgets(line,MAXLINE,fp) == NULL)
      error->one(FLERR,"Did not find keyword in table file");
    if (strspn(line," \t\n\r") == strlen(line)) continue;  // blank line
    if (line[0] == '#') continue;                          // comment
    char *word = strtok(line," \t\n\r");
    if (strcmp(word,keyword) == 0) break;           // matching keyword
    fgets(line,MAXLINE,fp);                         // no match, skip section
    param_extract(tb,line);
    fgets(line,MAXLINE,fp);
    for (int i = 0; i < tb->ninput; i++) fgets(line,MAXLINE,fp);
  }

  // read args on 2nd line of section
  // allocate table arrays for file values

  fgets(line,MAXLINE,fp);
  param_extract(tb,line);
  memory->create(tb->rfile,tb->ninput,"pair:rfile");
  memory->create(tb->efile,tb->ninput,"pair:efile");
  memory->create(tb->ffile,tb->ninput,"pair:ffile");

  // setup bitmap parameters for table to read in

  tb->ntablebits = 0;
  int masklo,maskhi,nmask,nshiftbits;
  if (tb->rflag == BMP) {
    while (1 << tb->ntablebits < tb->ninput) tb->ntablebits++;
    if (1 << tb->ntablebits != tb->ninput)
      error->one(FLERR,"Bitmapped table is incorrect length in table file");
    init_bitmap(tb->rlo,tb->rhi,tb->ntablebits,masklo,maskhi,nmask,nshiftbits);
  }

  // read r,e,f table values from file
  // if rflag set, compute r
  // if rflag not set, use r from file

  int itmp;
  double rtmp;
  union_int_float_t rsq_lookup;

  fgets(line,MAXLINE,fp);
  for (int i = 0; i < tb->ninput; i++) {
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg %lg %lg",&itmp,&rtmp,&tb->efile[i],&tb->ffile[i]);

    if (tb->rflag == RLINEAR)
      rtmp = tb->rlo + (tb->rhi - tb->rlo)*i/(tb->ninput-1);
    else if (tb->rflag == RSQ) {
      rtmp = tb->rlo*tb->rlo +
        (tb->rhi*tb->rhi - tb->rlo*tb->rlo)*i/(tb->ninput-1);
      rtmp = sqrt(rtmp);
    } else if (tb->rflag == BMP) {
      rsq_lookup.i = i << nshiftbits;
      rsq_lookup.i |= masklo;
      if (rsq_lookup.f < tb->rlo*tb->rlo) {
        rsq_lookup.i = i << nshiftbits;
        rsq_lookup.i |= maskhi;
      }
      rtmp = sqrtf(rsq_lookup.f);
    }

    tb->rfile[i] = rtmp;
  }

  // close file

  fclose(fp);
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,rfile,efile,ffile,rflag,rlo,rhi,fpflag,fplo,fphi
------------------------------------------------------------------------- */

void PairTable::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me > 0) {
    memory->create(tb->rfile,tb->ninput,"pair:rfile");
    memory->create(tb->efile,tb->ninput,"pair:efile");
    memory->create(tb->ffile,tb->ninput,"pair:ffile");
  }

  MPI_Bcast(tb->rfile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->efile,tb->ninput,MPI_DOUBLE,0,world);
  MPI_Bcast(tb->ffile,tb->ninput,MPI_DOUBLE,0,world);

  MPI_Bcast(&tb->rflag,1,MPI_INT,0,world);
  if (tb->rflag) {
    MPI_Bcast(&tb->rlo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->rhi,1,MPI_DOUBLE,0,world);
  }
  MPI_Bcast(&tb->fpflag,1,MPI_INT,0,world);
  if (tb->fpflag) {
    MPI_Bcast(&tb->fplo,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&tb->fphi,1,MPI_DOUBLE,0,world);
  }
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in Table: e2file,f2file
------------------------------------------------------------------------- */

void PairTable::spline_table(Table *tb)
{
  memory->create(tb->e2file,tb->ninput,"pair:e2file");
  memory->create(tb->f2file,tb->ninput,"pair:f2file");

  double ep0 = - tb->ffile[0];
  double epn = - tb->ffile[tb->ninput-1];
  spline(tb->rfile,tb->efile,tb->ninput,ep0,epn,tb->e2file);

  if (tb->fpflag == 0) {
    tb->fplo = (tb->ffile[1] - tb->ffile[0]) / (tb->rfile[1] - tb->rfile[0]);
    tb->fphi = (tb->ffile[tb->ninput-1] - tb->ffile[tb->ninput-2]) /
      (tb->rfile[tb->ninput-1] - tb->rfile[tb->ninput-2]);
  }

  double fp0 = tb->fplo;
  double fpn = tb->fphi;
  spline(tb->rfile,tb->ffile,tb->ninput,fp0,fpn,tb->f2file);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value R/RSQ/BITMAP lo hi FP fplo fphi
   N is required, other params are optional
------------------------------------------------------------------------- */

void PairTable::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->rflag = NONE;
  tb->fpflag = 0;

  char *word = strtok(line," \t\n\r\f");
  while (word) {
    if (strcmp(word,"N") == 0) {
      word = strtok(NULL," \t\n\r\f");
      tb->ninput = atoi(word);
    } else if (strcmp(word,"R") == 0 || strcmp(word,"RSQ") == 0 ||
               strcmp(word,"BITMAP") == 0) {
      if (strcmp(word,"R") == 0) tb->rflag = RLINEAR;
      else if (strcmp(word,"RSQ") == 0) tb->rflag = RSQ;
      else if (strcmp(word,"BITMAP") == 0) tb->rflag = BMP;
      word = strtok(NULL," \t\n\r\f");
      tb->rlo = atof(word);
      word = strtok(NULL," \t\n\r\f");
      tb->rhi = atof(word);
    } else if (strcmp(word,"FP") == 0) {
      tb->fpflag = 1;
      word = strtok(NULL," \t\n\r\f");
      tb->fplo = atof(word);
      word = strtok(NULL," \t\n\r\f");
      tb->fphi = atof(word);
    } else {
      printf("WORD: %s\n",word);
      error->one(FLERR,"Invalid keyword in pair table parameters");
    }
    word = strtok(NULL," \t\n\r\f");
  }

  if (tb->ninput == 0) error->one(FLERR,"Pair table parameters did not set N");
}

/* ----------------------------------------------------------------------
   compute r,e,f vectors from splined values
------------------------------------------------------------------------- */

void PairTable::compute_table(Table *tb)
{
  int tlm1 = tablength-1;

  // inner = inner table bound
  // cut = outer table bound
  // delta = table spacing in rsq for N-1 bins

  double inner;
  if (tb->rflag) inner = tb->rlo;
  else inner = tb->rfile[0];
  tb->innersq = inner*inner;
  tb->delta = (tb->cut*tb->cut - tb->innersq) / tlm1;
  tb->invdelta = 1.0/tb->delta;

  // direct lookup tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // e,f = value at midpt of bin
  // e,f are N-1 in length since store 1 value at bin midpt
  // f is converted to f/r when stored in f[i]
  // e,f are never a match to read-in values, always computed via spline interp

  if (tabstyle == LOOKUP) {
    memory->create(tb->e,tlm1,"pair:e");
    memory->create(tb->f,tlm1,"pair:f");

    double r,rsq;
    for (int i = 0; i < tlm1; i++) {
      rsq = tb->innersq + (i+0.5)*tb->delta;
      r = sqrt(rsq);
      tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
      tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
    }
  }

  // linear tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // rsq,e,f = value at lower edge of bin
  // de,df values = delta from lower edge to upper edge of bin
  // rsq,e,f are N in length so de,df arrays can compute difference
  // f is converted to f/r when stored in f[i]
  // e,f can match read-in values, else compute via spline interp

  if (tabstyle == LINEAR) {
    memory->create(tb->rsq,tablength,"pair:rsq");
    memory->create(tb->e,tablength,"pair:e");
    memory->create(tb->f,tablength,"pair:f");
    memory->create(tb->de,tlm1,"pair:de");
    memory->create(tb->df,tlm1,"pair:df");

    double r,rsq;
    for (int i = 0; i < tablength; i++) {
      rsq = tb->innersq + i*tb->delta;
      r = sqrt(rsq);
      tb->rsq[i] = rsq;
      if (tb->match) {
        tb->e[i] = tb->efile[i];
        tb->f[i] = tb->ffile[i]/r;
      } else {
        tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
      }
    }

    for (int i = 0; i < tlm1; i++) {
      tb->de[i] = tb->e[i+1] - tb->e[i];
      tb->df[i] = tb->f[i+1] - tb->f[i];
    }
  }

  // cubic spline tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // rsq,e,f = value at lower edge of bin
  // e2,f2 = spline coefficient for each bin
  // rsq,e,f,e2,f2 are N in length so have N-1 spline bins
  // f is converted to f/r after e is splined
  // e,f can match read-in values, else compute via spline interp

  if (tabstyle == SPLINE) {
    memory->create(tb->rsq,tablength,"pair:rsq");
    memory->create(tb->e,tablength,"pair:e");
    memory->create(tb->f,tablength,"pair:f");
    memory->create(tb->e2,tablength,"pair:e2");
    memory->create(tb->f2,tablength,"pair:f2");

    tb->deltasq6 = tb->delta*tb->delta / 6.0;

    double r,rsq;
    for (int i = 0; i < tablength; i++) {
      rsq = tb->innersq + i*tb->delta;
      r = sqrt(rsq);
      tb->rsq[i] = rsq;
      if (tb->match) {
        tb->e[i] = tb->efile[i];
        tb->f[i] = tb->ffile[i]/r;
      } else {
        tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r);
      }
    }

    // ep0,epn = dh/dg at inner and at cut
    // h(r) = e(r) and g(r) = r^2
    // dh/dg = (de/dr) / 2r = -f/2r

    double ep0 = - tb->f[0] / (2.0 * sqrt(tb->innersq));
    double epn = - tb->f[tlm1] / (2.0 * tb->cut);
    spline(tb->rsq,tb->e,tablength,ep0,epn,tb->e2);

    // fp0,fpn = dh/dg at inner and at cut
    // h(r) = f(r)/r and g(r) = r^2
    // dh/dg = (1/r df/dr - f/r^2) / 2r
    // dh/dg in secant approx = (f(r2)/r2 - f(r1)/r1) / (g(r2) - g(r1))

    double fp0,fpn;
    double secant_factor = 0.1;
    if (tb->fpflag) fp0 = (tb->fplo/sqrt(tb->innersq) - tb->f[0]/tb->innersq) /
      (2.0 * sqrt(tb->innersq));
    else {
      double rsq1 = tb->innersq;
      double rsq2 = rsq1 + secant_factor*tb->delta;
      fp0 = (splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,sqrt(rsq2)) /
             sqrt(rsq2) - tb->f[0] / sqrt(rsq1)) / (secant_factor*tb->delta);
    }

    if (tb->fpflag && tb->cut == tb->rfile[tb->ninput-1]) fpn =
      (tb->fphi/tb->cut - tb->f[tlm1]/(tb->cut*tb->cut)) / (2.0 * tb->cut);
    else {
      double rsq2 = tb->cut * tb->cut;
      double rsq1 = rsq2 - secant_factor*tb->delta;
      fpn = (tb->f[tlm1] / sqrt(rsq2) -
             splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,sqrt(rsq1)) /
             sqrt(rsq1)) / (secant_factor*tb->delta);
    }

    for (int i = 0; i < tablength; i++) tb->f[i] /= sqrt(tb->rsq[i]);
    spline(tb->rsq,tb->f,tablength,fp0,fpn,tb->f2);
  }

  // bitmapped linear tables
  // 2^N bins from inner to cut, spaced in bitmapped manner
  // f is converted to f/r when stored in f[i]
  // e,f can match read-in values, else compute via spline interp

  if (tabstyle == BITMAP) {
    double r;
    union_int_float_t rsq_lookup;
    int masklo,maskhi;

    // linear lookup tables of length ntable = 2^n
    // stored value = value at lower edge of bin

    init_bitmap(inner,tb->cut,tablength,masklo,maskhi,tb->nmask,tb->nshiftbits);
    int ntable = 1 << tablength;
    int ntablem1 = ntable - 1;

    memory->create(tb->rsq,ntable,"pair:rsq");
    memory->create(tb->e,ntable,"pair:e");
    memory->create(tb->f,ntable,"pair:f");
    memory->create(tb->de,ntable,"pair:de");
    memory->create(tb->df,ntable,"pair:df");
    memory->create(tb->drsq,ntable,"pair:drsq");

    union_int_float_t minrsq_lookup;
    minrsq_lookup.i = 0 << tb->nshiftbits;
    minrsq_lookup.i |= maskhi;

    for (int i = 0; i < ntable; i++) {
      rsq_lookup.i = i << tb->nshiftbits;
      rsq_lookup.i |= masklo;
      if (rsq_lookup.f < tb->innersq) {
        rsq_lookup.i = i << tb->nshiftbits;
        rsq_lookup.i |= maskhi;
      }
      r = sqrtf(rsq_lookup.f);
      tb->rsq[i] = rsq_lookup.f;
      if (tb->match) {
        tb->e[i] = tb->efile[i];
        tb->f[i] = tb->ffile[i]/r;
      } else {
        tb->e[i] = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        tb->f[i] = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
      }
      minrsq_lookup.f = MIN(minrsq_lookup.f,rsq_lookup.f);
    }

    tb->innersq = minrsq_lookup.f;

    for (int i = 0; i < ntablem1; i++) {
      tb->de[i] = tb->e[i+1] - tb->e[i];
      tb->df[i] = tb->f[i+1] - tb->f[i];
      tb->drsq[i] = 1.0/(tb->rsq[i+1] - tb->rsq[i]);
    }

    // get the delta values for the last table entries
    // tables are connected periodically between 0 and ntablem1

    tb->de[ntablem1] = tb->e[0] - tb->e[ntablem1];
    tb->df[ntablem1] = tb->f[0] - tb->f[ntablem1];
    tb->drsq[ntablem1] = 1.0/(tb->rsq[0] - tb->rsq[ntablem1]);

    // get the correct delta values at itablemax
    // smallest r is in bin itablemin
    // largest r is in bin itablemax, which is itablemin-1,
    //   or ntablem1 if itablemin=0

    // deltas at itablemax only needed if corresponding rsq < cut*cut
    // if so, compute deltas between rsq and cut*cut
    //   if tb->match, data at cut*cut is unavailable, so we'll take
    //   deltas at itablemax-1 as a good approximation

    double e_tmp,f_tmp;
    int itablemin = minrsq_lookup.i & tb->nmask;
    itablemin >>= tb->nshiftbits;
    int itablemax = itablemin - 1;
    if (itablemin == 0) itablemax = ntablem1;
    int itablemaxm1 = itablemax - 1;
    if (itablemax == 0) itablemaxm1 = ntablem1;
    rsq_lookup.i = itablemax << tb->nshiftbits;
    rsq_lookup.i |= maskhi;
    if (rsq_lookup.f < tb->cut*tb->cut) {
      if (tb->match) {
        tb->de[itablemax] = tb->de[itablemaxm1];
        tb->df[itablemax] = tb->df[itablemaxm1];
        tb->drsq[itablemax] = tb->drsq[itablemaxm1];
      } else {
            rsq_lookup.f = tb->cut*tb->cut;
        r = sqrtf(rsq_lookup.f);
        e_tmp = splint(tb->rfile,tb->efile,tb->e2file,tb->ninput,r);
        f_tmp = splint(tb->rfile,tb->ffile,tb->f2file,tb->ninput,r)/r;
        tb->de[itablemax] = e_tmp - tb->e[itablemax];
        tb->df[itablemax] = f_tmp - tb->f[itablemax];
        tb->drsq[itablemax] = 1.0/(rsq_lookup.f - tb->rsq[itablemax]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set all ptrs in a table to NULL, so can be freed safely
------------------------------------------------------------------------- */

void PairTable::null_table(Table *tb)
{
  tb->rfile = tb->efile = tb->ffile = NULL;
  tb->e2file = tb->f2file = NULL;
  tb->rsq = tb->drsq = tb->e = tb->de = NULL;
  tb->f = tb->df = tb->e2 = tb->f2 = NULL;
}

/* ----------------------------------------------------------------------
   free all arrays in a table
------------------------------------------------------------------------- */

void PairTable::free_table(Table *tb)
{
  memory->destroy(tb->rfile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->e2file);
  memory->destroy(tb->f2file);

  memory->destroy(tb->rsq);
  memory->destroy(tb->drsq);
  memory->destroy(tb->e);
  memory->destroy(tb->de);
  memory->destroy(tb->f);
  memory->destroy(tb->df);
  memory->destroy(tb->e2);
  memory->destroy(tb->f2);
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

void PairTable::spline(double *x, double *y, int n,
                       double yp1, double ypn, double *y2)
{
  int i,k;
  double p,qn,sig,un;
  double *u = new double[n];

  if (yp1 > 0.99e30) y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
  }
  for (i = 1; i < n-1; i++) {
    sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0) / p;
    u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
    u[i] = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
  }
  if (ypn > 0.99e30) qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
  }
  y2[n-1] = (un-qn*u[n-2]) / (qn*y2[n-2] + 1.0);
  for (k = n-2; k >= 0; k--) y2[k] = y2[k]*y2[k+1] + u[k];

  delete [] u;
}

/* ---------------------------------------------------------------------- */

double PairTable::splint(double *xa, double *ya, double *y2a, int n, double x)
{
  int klo,khi,k;
  double h,b,a,y;

  klo = 0;
  khi = n-1;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi]-xa[klo];
  a = (xa[khi]-x) / h;
  b = (x-xa[klo]) / h;
  y = a*ya[klo] + b*ya[khi] +
    ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;
  return y;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTable::write_restart(FILE *fp)
{
  write_restart_settings(fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTable::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTable::write_restart_settings(FILE *fp)
{
  fwrite(&tabstyle,sizeof(int),1,fp);
  fwrite(&tablength,sizeof(int),1,fp);
  fwrite(&ewaldflag,sizeof(int),1,fp);
  fwrite(&pppmflag,sizeof(int),1,fp);
  fwrite(&msmflag,sizeof(int),1,fp);
  fwrite(&dispersionflag,sizeof(int),1,fp);
  fwrite(&tip4pflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTable::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&tabstyle,sizeof(int),1,fp);
    fread(&tablength,sizeof(int),1,fp);
    fread(&ewaldflag,sizeof(int),1,fp);
    fread(&pppmflag,sizeof(int),1,fp);
    fread(&msmflag,sizeof(int),1,fp);
    fread(&dispersionflag,sizeof(int),1,fp);
    fread(&tip4pflag,sizeof(int),1,fp);
  }
  MPI_Bcast(&tabstyle,1,MPI_INT,0,world);
  MPI_Bcast(&tablength,1,MPI_INT,0,world);
  MPI_Bcast(&ewaldflag,1,MPI_INT,0,world);
  MPI_Bcast(&pppmflag,1,MPI_INT,0,world);
  MPI_Bcast(&msmflag,1,MPI_INT,0,world);
  MPI_Bcast(&dispersionflag,1,MPI_INT,0,world);
  MPI_Bcast(&tip4pflag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairTable::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  int itable;
  double fraction,value,a,b,phi;
  int tlm1 = tablength - 1;

  Table *tb = &tables[tabindex[itype][jtype]];
  if (rsq < tb->innersq) error->one(FLERR,"Pair distance < table inner cutoff");

  if (tabstyle == LOOKUP) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    fforce = factor_lj * tb->f[itable];
  } else if (tabstyle == LINEAR) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    fraction = (rsq - tb->rsq[itable]) * tb->invdelta;
    value = tb->f[itable] + fraction*tb->df[itable];
    fforce = factor_lj * value;
  } else if (tabstyle == SPLINE) {
    itable = static_cast<int> ((rsq-tb->innersq) * tb->invdelta);
    if (itable >= tlm1) error->one(FLERR,"Pair distance > table outer cutoff");
    b = (rsq - tb->rsq[itable]) * tb->invdelta;
    a = 1.0 - b;
    value = a * tb->f[itable] + b * tb->f[itable+1] +
      ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
      tb->deltasq6;
    fforce = factor_lj * value;
  } else {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    itable = rsq_lookup.i & tb->nmask;
    itable >>= tb->nshiftbits;
    fraction = (rsq_lookup.f - tb->rsq[itable]) * tb->drsq[itable];
    value = tb->f[itable] + fraction*tb->df[itable];
    fforce = factor_lj * value;
  }

  if (tabstyle == LOOKUP)
    phi = tb->e[itable];
  else if (tabstyle == LINEAR || tabstyle == BITMAP)
    phi = tb->e[itable] + fraction*tb->de[itable];
  else
    phi = a * tb->e[itable] + b * tb->e[itable+1] +
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) * tb->deltasq6;
  return factor_lj*phi;
}

/* ----------------------------------------------------------------------
   return the Coulomb cutoff for tabled potentials
   called by KSpace solvers which require that all pairwise cutoffs be the same
   loop over all tables not just those indexed by tabindex[i][j] since
     no way to know which tables are active since pair::init() not yet called
------------------------------------------------------------------------- */

void *PairTable::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") != 0) return NULL;
  if (ntables == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut_coul = tables[0].cut;
  for (int m = 1; m < ntables; m++)
    if (tables[m].cut != cut_coul)
      error->all(FLERR,
                 "Pair table cutoffs must all be equal to use with KSpace");
  dim = 0;
  return &tables[0].cut;
}
