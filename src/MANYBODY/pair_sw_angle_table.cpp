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
   Contributing author: Christoph Scherer (MPIP Mainz)
   scherer@mpip-mainz.mpg.de
------------------------------------------------------------------------- */

#include "pair_sw_angle_table.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "table_file_reader.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

using MathConst::DEG2RAD;
using MathConst::MY_PI;
using MathConst::RAD2DEG;

static constexpr int DELTA = 4;

enum { LINEAR, SPLINE };

static constexpr double TINY = 1.0e-10;

/* ---------------------------------------------------------------------- */

PairSWAngleTable::PairSWAngleTable(LAMMPS *lmp) : PairSW(lmp), table_params(nullptr)
{
  unit_convert_flag = utils::NOCONVERT;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairSWAngleTable::~PairSWAngleTable()
{
  if (copymode) return;

  for (int m = 0; m < nparams; m++) free_param(&table_params[m]); // free_param will call free_table
  memory->destroy(params);
  memory->destroy(table_params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
  }
}

/* ---------------------------------------------------------------------- */

void PairSWAngleTable::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;

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
      rsq = delx*delx + dely*dely + delz*delz;

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (rsq >= params[ijparam].cutsq) {
        continue;
      } else {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort/2;
          memory->grow(neighshort,maxshort,"pair:neighshort");
        }
      }

      jtag = tag[j];
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      twobody(&params[ijparam],rsq,fpair,eflag,evdwl);

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    jnumm1 = numshort - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      double fjxtmp,fjytmp,fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort[kk];
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        threebody_table(&params[ijparam],&params[ikparam],&table_params[ijkparam],
                        rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

        fxtmp -= fj[0] + fk[0];
        fytmp -= fj[1] + fk[1];
        fztmp -= fj[2] + fk[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
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

void PairSWAngleTable::read_file(char *file)
{
  if (params) {
    for (int m = 0; m < nparams; m++) free_param(&table_params[m]); // free_param will call free_table
    memory->destroy(params);
    memory->destroy(table_params);
    memory->destroy(elem3param);
  }

  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "sw", unit_convert_flag);
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
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                              "pair:params");
          table_params = (ParamTable *) memory->srealloc(table_params,maxparam*sizeof(ParamTable),
                                              "pair:table_params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA*sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        params[nparams].epsilon  = values.next_double();
        params[nparams].sigma    = values.next_double();
        params[nparams].littlea  = values.next_double();
        params[nparams].lambda   = values.next_double();
        params[nparams].gamma    = values.next_double();
        params[nparams].costheta = values.next_double();
        params[nparams].biga     = values.next_double();
        params[nparams].bigb     = values.next_double();
        params[nparams].powerp   = values.next_double();
        params[nparams].powerq   = values.next_double();
        params[nparams].tol      = values.next_double();

        // read parameters of angle table
        std::string name = values.next_string();
        table_params[nparams].tablenamelength = name.length()+1;
        table_params[nparams].tablename = utils::strdup(name);

        name = values.next_string();
        table_params[nparams].keywordlength = name.length()+1;
        table_params[nparams].keyword = utils::strdup(name);

        name = values.next_string();
        if (name == "linear") table_params[nparams].tabstyle = LINEAR;
        else if (name == "spline") table_params[nparams].tabstyle = SPLINE;
        else error->all(FLERR,"Unknown table style {} of angle table file", name);
        table_params[nparams].tablength = values.next_int();

      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      if (params[nparams].epsilon < 0.0 || params[nparams].sigma < 0.0 ||
          params[nparams].littlea < 0.0 || params[nparams].lambda < 0.0 ||
          params[nparams].gamma < 0.0 || params[nparams].biga < 0.0 ||
          params[nparams].bigb < 0.0 || params[nparams].powerp < 0.0 ||
          params[nparams].powerq < 0.0 || params[nparams].tol < 0.0 ||
          table_params[nparams].tabstyle < 0.0 || table_params[nparams].tablength < 0.0)
        error->one(FLERR,"Illegal Stillinger-Weber parameter");

      nparams++;
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");
    table_params = (ParamTable *) memory->srealloc(table_params,maxparam*sizeof(ParamTable), "pair:table_params");
  }

  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
  MPI_Bcast(table_params, maxparam*sizeof(ParamTable), MPI_BYTE, 0, world);

  // for each set of parameters, broadcast table name and keyword and read angle table
  for (int m = 0; m < nparams; ++m){
    if (comm->me != 0) {
      table_params[m].tablename = new char[table_params[m].tablenamelength];
      table_params[m].keyword = new char[table_params[m].keywordlength];
    }
    MPI_Bcast(table_params[m].tablename, table_params[m].tablenamelength, MPI_CHAR, 0, world);
    MPI_Bcast(table_params[m].keyword, table_params[m].keywordlength, MPI_CHAR, 0, world);

    // initialize angtable
    memory->create(table_params[m].angtable,1,"table_params:angtable");
    null_table(table_params[m].angtable);

    // call read_table to read corresponding tabulated angle file (only called by process 0)
    if (comm->me == 0) read_table(table_params[m].angtable,table_params[m].tablename,table_params[m].keyword);

    // broadcast read in angtable to all processes
    bcast_table(table_params[m].angtable);

    // the following table manipulations are done in all processes
    // error check on table parameters
    if (table_params[m].angtable->ninput <= 1) error->one(FLERR,"Invalid angle table length");

    // error check on parameter range of angle table
    double alo,ahi;
    alo = table_params[m].angtable->afile[0];
    ahi = table_params[m].angtable->afile[table_params[m].angtable->ninput-1];
    if (fabs(alo-0.0) > TINY || fabs(ahi-180.0) > TINY)
      error->all(FLERR,"Angle table must range from 0 to 180 degrees");

    // convert theta from degrees to radians
    for (int i = 0; i < table_params[m].angtable->ninput; ++i){
      table_params[m].angtable->afile[i] *= MY_PI/180.0;
      table_params[m].angtable->ffile[i] *= 180.0/MY_PI;
    }

    // spline read-in table and compute a,e,f vectors within table
    spline_table(table_params[m].angtable);
    // compute_table needs parameter params[m].length for this specific interaction as
    // read in value length from .sw file can be different from value in angle table file
    compute_table(table_params[m].angtable,table_params[m].tablength);
  }
}

/* ---------------------------------------------------------------------- */

void PairSWAngleTable::threebody_table(Param *paramij, Param *paramik, ParamTable *table_paramijk,
                                    double rsq1, double rsq2, double *delr1, double *delr2,
                                    double *fj, double *fk, int eflag, double &eng)
{
  double r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  double r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  double rinv12,cs,facexp;
  double ftheta,facradtable,frad1table,frad2table,var;
  double acosprime,gradj1,gradj2,gradk1,gradk2,fprimetheta;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - paramij->cut);
  gsrainv1 = paramij->sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  rainv2 = 1.0/(r2 - paramik->cut);
  gsrainv2 = paramik->sigma_gamma * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2;
  expgsrainv2 = exp(gsrainv2);

  facexp = expgsrainv1*expgsrainv2;

  rinv12 = 1.0/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;

  var = acos(cs);

  // look up energy (f(theta), ftheta) and force (df(theta)/dtheta, fprimetheta) at
  // angle theta (var) in angle table belonging to parameter set paramijk
  uf_lookup(table_paramijk, var, ftheta, fprimetheta);

  acosprime = 1.0 / (sqrt(1 - cs*cs ) );

  facradtable = facexp*ftheta;
  frad1table = facradtable*gsrainvsq1;
  frad2table = facradtable*gsrainvsq2;
  gradj1 = acosprime * cs * rinvsq1  * facexp * fprimetheta;
  gradj2 = acosprime * rinv12 * facexp * fprimetheta;
  gradk1 = acosprime * cs * rinvsq2  * facexp * fprimetheta;
  gradk2 = acosprime * rinv12 * facexp * fprimetheta;

  fj[0] = delr1[0]*(frad1table+gradj1)-delr2[0]*gradj2;
  fj[1] = delr1[1]*(frad1table+gradj1)-delr2[1]*gradj2;
  fj[2] = delr1[2]*(frad1table+gradj1)-delr2[2]*gradj2;

  fk[0] = delr2[0]*(frad2table+gradk1)-delr1[0]*gradk2;
  fk[1] = delr2[1]*(frad2table+gradk1)-delr1[1]*gradk2;
  fk[2] = delr2[2]*(frad2table+gradk1)-delr1[2]*gradk2;

  if (eflag) eng = facradtable;
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void PairSWAngleTable::read_table(Table *tb, char *file, char *keyword)
{
  TableFileReader reader(lmp, file, "angletable");

  char *line = reader.find_section_start(keyword);

  if (!line) { error->one(FLERR, "Did not find keyword in table file"); }

  // read args on 2nd line of section
  // allocate table arrays for file values

  line = reader.next_line();
  param_extract(tb, line);
  memory->create(tb->afile, tb->ninput, "angle:afile");
  memory->create(tb->efile, tb->ninput, "angle:efile");
  memory->create(tb->ffile, tb->ninput, "angle:ffile");

  // read a,e,f table values from file

  int cerror = 0;
  reader.skip_line();
  for (int i = 0; i < tb->ninput; i++) {
    line = reader.next_line(4);
    try {
      ValueTokenizer values(line);
      values.next_int();
      tb->afile[i] = values.next_double();
      tb->efile[i] = values.next_double();
      tb->ffile[i] = values.next_double();
    } catch (TokenizerException &) {
      ++cerror;
    }
  }

  // warn if data was read incompletely, e.g. columns were missing

  if (cerror)
    error->warning(FLERR, "{} of {} lines in table incomplete or could not be parsed", cerror,
                   tb->ninput);
}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file,f2file
------------------------------------------------------------------------- */

void PairSWAngleTable::spline_table(Table *tb)
{
  memory->create(tb->e2file, tb->ninput, "angle:e2file");
  memory->create(tb->f2file, tb->ninput, "angle:f2file");

  double ep0 = -tb->ffile[0];
  double epn = -tb->ffile[tb->ninput - 1];
  spline(tb->afile, tb->efile, tb->ninput, ep0, epn, tb->e2file);

  if (tb->fpflag == 0) {
    tb->fplo = (tb->ffile[1] - tb->ffile[0]) / (tb->afile[1] - tb->afile[0]);
    tb->fphi = (tb->ffile[tb->ninput - 1] - tb->ffile[tb->ninput - 2]) /
        (tb->afile[tb->ninput - 1] - tb->afile[tb->ninput - 2]);
  }

  double fp0 = tb->fplo;
  double fpn = tb->fphi;
  spline(tb->afile, tb->ffile, tb->ninput, fp0, fpn, tb->f2file);
}

/* ----------------------------------------------------------------------
   compute a,e,f vectors from splined values
------------------------------------------------------------------------- */

void PairSWAngleTable::compute_table(Table *tb,int length)
{
  // delta = table spacing in angle for N-1 bins

  int tlm1 = length - 1;
  tb->delta = MY_PI / tlm1;
  tb->invdelta = 1.0 / tb->delta;
  tb->deltasq6 = tb->delta * tb->delta / 6.0;

  // N-1 evenly spaced bins in angle from 0 to PI
  // ang,e,f = value at lower edge of bin
  // de,df values = delta values of e,f
  // ang,e,f are N in length so de,df arrays can compute difference

  memory->create(tb->ang, length, "angle:ang");
  memory->create(tb->e, length, "angle:e");
  memory->create(tb->de, length, "angle:de");
  memory->create(tb->f, length, "angle:f");
  memory->create(tb->df, length, "angle:df");
  memory->create(tb->e2, length, "angle:e2");
  memory->create(tb->f2, length, "angle:f2");

  double a;
  for (int i = 0; i < length; i++) {
    a = i * tb->delta;
    tb->ang[i] = a;
    tb->e[i] = splint(tb->afile, tb->efile, tb->e2file, tb->ninput, a);
    tb->f[i] = splint(tb->afile, tb->ffile, tb->f2file, tb->ninput, a);
  }

  for (int i = 0; i < tlm1; i++) {
    tb->de[i] = tb->e[i + 1] - tb->e[i];
    tb->df[i] = tb->f[i + 1] - tb->f[i];
  }
  // get final elements from linear extrapolation
  tb->de[tlm1] = 2.0 * tb->de[tlm1 - 1] - tb->de[tlm1 - 2];
  tb->df[tlm1] = 2.0 * tb->df[tlm1 - 1] - tb->df[tlm1 - 2];

  double ep0 = -tb->f[0];
  double epn = -tb->f[tlm1];
  spline(tb->ang, tb->e, length, ep0, epn, tb->e2);
  spline(tb->ang, tb->f, length, tb->fplo, tb->fphi, tb->f2);
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,afile,efile,ffile,fpflag,fplo,fphi,theta0
------------------------------------------------------------------------- */

void PairSWAngleTable::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);

  int me;
  MPI_Comm_rank(world, &me);
  if (me > 0) {
    memory->create(tb->afile, tb->ninput, "angle:afile");
    memory->create(tb->efile, tb->ninput, "angle:efile");
    memory->create(tb->ffile, tb->ninput, "angle:ffile");
  }

  MPI_Bcast(tb->afile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->efile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->ffile, tb->ninput, MPI_DOUBLE, 0, world);

  MPI_Bcast(&tb->fpflag, 1, MPI_INT, 0, world);
  if (tb->fpflag) {
    MPI_Bcast(&tb->fplo, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&tb->fphi, 1, MPI_DOUBLE, 0, world);
  }
  MPI_Bcast(&tb->theta0, 1, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairSWAngleTable::null_table(Table *tb)
{
  tb->afile = tb->efile = tb->ffile = nullptr;
  tb->e2file = tb->f2file = nullptr;
  tb->ang = tb->e = tb->de = nullptr;
  tb->f = tb->df = tb->e2 = tb->f2 = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairSWAngleTable::free_table(Table *tb)
{
  memory->destroy(tb->afile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->e2file);
  memory->destroy(tb->f2file);

  memory->destroy(tb->ang);
  memory->destroy(tb->e);
  memory->destroy(tb->de);
  memory->destroy(tb->f);
  memory->destroy(tb->df);
  memory->destroy(tb->e2);
  memory->destroy(tb->f2);
}

/* ---------------------------------------------------------------------- */

void PairSWAngleTable::free_param(ParamTable *pm)
{
  // call free_table to destroy associated angle table
  free_table(pm->angtable);
  // then destroy associated angle table
  delete[] pm->keyword;
  delete[] pm->tablename;
  memory->sfree(pm->angtable);
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ theta0
   N is required, other params are optional

   only called by read_table, only called by proc 0
------------------------------------------------------------------------- */

void PairSWAngleTable::param_extract(Table *tb, char *line)
{
  tb->ninput = 0;
  tb->fpflag = 0;
  tb->theta0 = MY_PI;

  try {
    ValueTokenizer values(line);

    while (values.has_next()) {
      std::string word = values.next_string();

      if (word == "N") {
        tb->ninput = values.next_int();
      } else if (word == "FP") {
        tb->fpflag = 1;
        tb->fplo = values.next_double();
        tb->fphi = values.next_double();
        tb->fplo *= RAD2DEG * RAD2DEG;
        tb->fphi *= RAD2DEG * RAD2DEG;
      } else if (word == "EQ") {
        tb->theta0 = DEG2RAD * values.next_double();
      } else {
        error->one(FLERR, "Invalid keyword in angle table parameters");
      }
    }
  } catch (TokenizerException &e) {
    error->one(FLERR, e.what());
  }

  if (tb->ninput == 0) error->one(FLERR, "Angle table parameters did not set N");
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

void PairSWAngleTable::spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
  int i, k;
  double p, qn, sig, un;
  double *u = new double[n];

  if (yp1 > 0.99e300)
    y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }
  for (i = 1; i < n - 1; i++) {
    sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }
  if (ypn > 0.99e300)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }
  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
  for (k = n - 2; k >= 0; k--) y2[k] = y2[k] * y2[k + 1] + u[k];

  delete[] u;
}
/* ---------------------------------------------------------------------- */

double PairSWAngleTable::splint(double *xa, double *ya, double *y2a, int n, double x)
{
  int klo, khi, k;
  double h, b, a, y;

  klo = 0;
  khi = n - 1;
  while (khi - klo > 1) {
    k = (khi + klo) >> 1;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }
  h = xa[khi] - xa[klo];
  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  y = a * ya[klo] + b * ya[khi] +
      ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
  return y;
}

/* ----------------------------------------------------------------------
   calculate potential u and force f at angle x
------------------------------------------------------------------------- */

void PairSWAngleTable::uf_lookup(ParamTable *pm, double x, double &u, double &f)
{
  if (!std::isfinite(x)) { error->one(FLERR, "Illegal angle in angle style table"); }

  double fraction,a,b;

  // invdelta is based on tablength-1
  int itable = static_cast<int>(x * pm->angtable->invdelta);
  if (itable < 0) itable = 0;
  if (itable >=  pm->tablength) itable =  pm->tablength - 1;

  if (pm->tabstyle == LINEAR) {
    fraction = (x - pm->angtable->ang[itable]) * pm->angtable->invdelta;
    u = pm->angtable->e[itable] + fraction*pm->angtable->de[itable];
    f = pm->angtable->f[itable] + fraction*pm->angtable->df[itable];
  } else if (pm->tabstyle == SPLINE) {
    fraction = (x - pm->angtable->ang[itable]) * pm->angtable->invdelta;

    b = (x - pm->angtable->ang[itable]) * pm->angtable->invdelta;
    a = 1.0 - b;
    u = a * pm->angtable->e[itable] + b * pm->angtable->e[itable+1] +
      ((a * a * a - a) * pm->angtable->e2[itable] + (b * b * b - b) * pm->angtable->e2[itable+1]) *
      pm->angtable->deltasq6;
    f = a * pm->angtable->f[itable] + b * pm->angtable->f[itable+1] +
      ((a * a * a - a) * pm->angtable->f2[itable] + (b * b * b - b) * pm->angtable->f2[itable+1]) *
      pm->angtable->deltasq6;
  }
}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSWAngleTable::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style sw/angle/table command");
}
