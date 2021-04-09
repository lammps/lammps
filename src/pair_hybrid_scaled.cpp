/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_hybrid_scaled.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "memory.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairHybridScaled::PairHybridScaled(LAMMPS *lmp)
  : PairHybrid(lmp), fsum(nullptr), scaleval(nullptr), scaleidx(nullptr)
{
  nmaxfsum = -1;
}

/* ---------------------------------------------------------------------- */

PairHybridScaled::~PairHybridScaled()
{
  memory->destroy(fsum);
  memory->destroy(scaleval);
}

/* ----------------------------------------------------------------------
  call each sub-style's compute() or compute_outer() function
  accumulate sub-style global/peratom energy/virial in hybrid
  for global vflag = VIRIAL_PAIR:
    each sub-style computes own virial[6]
    sum sub-style virial[6] to hybrid's virial[6]
  for global vflag = VIRIAL_FDOTR:
    call sub-style with adjusted vflag to prevent it calling
      virial_fdotr_compute()
    hybrid calls virial_fdotr_compute() on final accumulated f
------------------------------------------------------------------------- */

void PairHybridScaled::compute(int eflag, int vflag)
{
  int i,j,m,n;

  // update scale values from variables where needed

  const int nvars = scalevars.size();
  if (nvars > 0) {
    double *vals = new double[nvars];
    for (i = 0; i < nvars; ++i) {
      j = input->variable->find(scalevars[i].c_str());
      vals[i] = input->variable->compute_equal(j);
    }
    for (i = 0; i < nstyles; ++i) {
      if (scaleidx[i] >= 0)
        scaleval[i] = vals[scaleidx[i]];
    }
    delete[] vals;
  }

  // check if no_virial_fdotr_compute is set and global component of
  //   incoming vflag = VIRIAL_FDOTR
  // if so, reset vflag as if global component were VIRIAL_PAIR
  // necessary since one or more sub-styles cannot compute virial as F dot r

  if (no_virial_fdotr_compute && (vflag & VIRIAL_FDOTR))
    vflag = VIRIAL_PAIR | (vflag & ~VIRIAL_FDOTR);

  ev_init(eflag,vflag);

  // grow fsum array if needed, and copy existing forces (usually 0.0) to it.

  if (atom->nmax > nmaxfsum) {
    memory->destroy(fsum);
    nmaxfsum = atom->nmax;
    memory->create(fsum,nmaxfsum,3,"pair:fsum");
  }
  const int nall = atom->nlocal + atom->nghost;
  auto f = atom->f;
  for (i = 0; i < nall; ++i) {
    fsum[i][0] = f[i][0];
    fsum[i][1] = f[i][1];
    fsum[i][2] = f[i][2];
  }

  // check if global component of incoming vflag = VIRIAL_FDOTR
  // if so, reset vflag passed to substyle so VIRIAL_FDOTR is turned off
  // necessary so substyle will not invoke virial_fdotr_compute()

  int vflag_substyle;
  if (vflag & VIRIAL_FDOTR) vflag_substyle = vflag & ~VIRIAL_FDOTR;
  else vflag_substyle = vflag;

  double *saved_special = save_special();

  // check if we are running with r-RESPA using the hybrid keyword

  Respa *respa = nullptr;
  respaflag = 0;
  if (utils::strmatch(update->integrate_style,"^respa")) {
    respa = (Respa *) update->integrate;
    if (respa->nhybrid_styles > 0) respaflag = 1;
  }

  for (m = 0; m < nstyles; m++) {

    // clear forces

    memset(&f[0][0],0,nall*3*sizeof(double));

    set_special(m);

    if (!respaflag || (respaflag && respa->hybrid_compute[m])) {

      // invoke compute() unless compute flag is turned off or
      // outerflag is set and sub-style has a compute_outer() method

      if (styles[m]->compute_flag == 0) continue;
      if (outerflag && styles[m]->respa_enable)
        styles[m]->compute_outer(eflag,vflag_substyle);
      else styles[m]->compute(eflag,vflag_substyle);
    }

    // add scaled forces to global sum
    const double scale = scaleval[m];
    for (i = 0; i < nall; ++i) {
      fsum[i][0] += scale*f[i][0];
      fsum[i][1] += scale*f[i][1];
      fsum[i][2] += scale*f[i][2];
    }

    restore_special(saved_special);

    // jump to next sub-style if r-RESPA does not want global accumulated data

    if (respaflag && !respa->tally_global) continue;

    if (eflag_global) {
      eng_vdwl += scale * styles[m]->eng_vdwl;
      eng_coul += scale * styles[m]->eng_coul;
    }
    if (vflag_global) {
      for (n = 0; n < 6; n++) virial[n] += scale * styles[m]->virial[n];
    }
    if (eflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double *eatom_substyle = styles[m]->eatom;
      for (i = 0; i < n; i++) eatom[i] += scale * eatom_substyle[i];
    }
    if (vflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double **vatom_substyle = styles[m]->vatom;
      for (i = 0; i < n; i++)
        for (j = 0; j < 6; j++)
          vatom[i][j] += scale * vatom_substyle[i][j];
    }

    // substyles may be CENTROID_SAME or CENTROID_AVAIL

    if (cvflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      if (styles[m]->centroidstressflag == CENTROID_AVAIL) {
        double **cvatom_substyle = styles[m]->cvatom;
        for (i = 0; i < n; i++)
          for (j = 0; j < 9; j++)
            cvatom[i][j] += scale * cvatom_substyle[i][j];
      } else {
        double **vatom_substyle = styles[m]->vatom;
        for (i = 0; i < n; i++) {
          for (j = 0; j < 6; j++) {
            cvatom[i][j] += scale * vatom_substyle[i][j];
          }
          for (j = 6; j < 9; j++) {
            cvatom[i][j] += scale * vatom_substyle[i][j-3];
          }
        }
      }
    }
  }

  // copy accumulated forces to original force array

  for (i = 0; i < nall; ++i) {
    f[i][0] = fsum[i][0];
    f[i][1] = fsum[i][1];
    f[i][2] = fsum[i][2];
  }
  delete [] saved_special;

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   create one pair style for each arg in list
------------------------------------------------------------------------- */

void PairHybridScaled::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");
  if (lmp->kokkos && !utils::strmatch(force->pair_style,"^hybrid.*/kk$"))
    error->all(FLERR,fmt::format("Must use pair_style {}/kk with Kokkos",
                                 force->pair_style));

  // delete old lists, since cannot just change settings

  if (nstyles > 0) {
    for (int m = 0; m < nstyles; m++) {
      delete styles[m];
      delete [] keywords[m];
      if (special_lj[m])   delete [] special_lj[m];
      if (special_coul[m]) delete [] special_coul[m];
    }
    delete[] styles;
    delete[] keywords;
    delete[] multiple;
    delete[] special_lj;
    delete[] special_coul;
    delete[] compute_tally;
    delete[] scaleval;
    delete[] scaleidx;
    scalevars.clear();
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    memory->destroy(nmap);
    memory->destroy(map);
  }
  allocated = 0;

  // allocate list of sub-styles as big as possibly needed if no extra args

  styles = new Pair*[narg];
  keywords = new char*[narg];
  multiple = new int[narg];

  special_lj = new double*[narg];
  special_coul = new double*[narg];
  compute_tally = new int[narg];

  scaleval = new double[narg];
  scaleidx = new int[narg];
  scalevars.reserve(narg);

  // allocate each sub-style
  // allocate uses suffix, but don't store suffix version in keywords,
  //   else syntax in coeff() will not match
  // call settings() with set of args that are not pair style names
  // use force->pair_map to determine which args these are

  int iarg,jarg,dummy;

  iarg = 0;
  nstyles = 0;
  while (iarg < narg-1) {

    // first process scale factor or variable
    // idx < 0 indicates constant value otherwise index in variable name list

    double val = 0.0;
    int idx = -1;
    if (utils::strmatch(arg[iarg],"^v_")) {
      for (std::size_t i=0; i < scalevars.size(); ++i) {
        if (scalevars[i] == arg[iarg]+2) {
          idx = i;
          break;
        }
      }
      if (idx < 0) {
        idx = scalevars.size();
        scalevars.push_back(arg[iarg]+2);
      }
    } else {
      val = utils::numeric(FLERR,arg[iarg],false,lmp);
    }
    scaleval[nstyles] = val;
    scaleidx[nstyles] = idx;
    ++iarg;

    if (utils::strmatch(arg[iarg],"^hybrid"))
      error->all(FLERR,"Pair style hybrid cannot have hybrid as an argument");
    if (strcmp(arg[iarg],"none") == 0)
      error->all(FLERR,"Pair style hybrid cannot have none as an argument");

    styles[nstyles] = force->new_pair(arg[iarg],1,dummy);
    force->store_style(keywords[nstyles],arg[iarg],0);
    special_lj[nstyles] = special_coul[nstyles] = nullptr;
    compute_tally[nstyles] = 1;

    // determine list of arguments for pair style settings
    // by looking for the next known pair style name.

    jarg = iarg + 1;
    while ((jarg < narg)
           && !force->pair_map->count(arg[jarg])
           && !lmp->match_style("pair",arg[jarg])) jarg++;

    // decrement to account for scale factor except when last argument

    if (jarg < narg) --jarg;

    styles[nstyles]->settings(jarg-iarg-1,arg+iarg+1);
    iarg = jarg;
    nstyles++;
  }

  // multiple[i] = 1 to M if sub-style used multiple times, else 0

  for (int i = 0; i < nstyles; i++) {
    int count = 0;
    for (int j = 0; j < nstyles; j++) {
      if (strcmp(keywords[j],keywords[i]) == 0) count++;
      if (j == i) multiple[i] = count;
    }
    if (count == 1) multiple[i] = 0;
  }

  // set pair flags from sub-style flags

  flags();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHybridScaled::coeff(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // 3rd arg = pair sub-style name
  // 4th arg = pair sub-style index if name used multiple times
  // allow for "none" as valid sub-style name

  int multflag;
  int m;

  for (m = 0; m < nstyles; m++) {
    multflag = 0;
    if (strcmp(arg[2],keywords[m]) == 0) {
      if (multiple[m]) {
        multflag = 1;
        if (narg < 4) error->all(FLERR,"Incorrect args for pair coefficients");
        if (!isdigit(arg[3][0]))
          error->all(FLERR,"Incorrect args for pair coefficients");
        int index = utils::inumeric(FLERR,arg[3],false,lmp);
        if (index == multiple[m]) break;
        else continue;
      } else break;
    }
  }

  int none = 0;
  if (m == nstyles) {
    if (strcmp(arg[2],"none") == 0) none = 1;
    else error->all(FLERR,"Pair coeff for hybrid has invalid style");
  }

  // move 1st/2nd args to 2nd/3rd args
  // if multflag: move 1st/2nd args to 3rd/4th args
  // just copy ptrs, since arg[] points into original input line

  arg[2+multflag] = arg[1];
  arg[1+multflag] = arg[0];

  // invoke sub-style coeff() starting with 1st remaining arg

  if (!none) styles[m]->coeff(narg-1-multflag,&arg[1+multflag]);

  // set setflag and which type pairs map to which sub-style
  // if sub-style is none: set hybrid subflag, wipe out map
  // else: set hybrid setflag & map only if substyle setflag is set
  //       if sub-style is new for type pair, add as multiple mapping
  //       if sub-style exists for type pair, don't add, just update coeffs

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (none) {
        setflag[i][j] = 1;
        nmap[i][j] = 0;
        count++;
      } else if (styles[m]->setflag[i][j]) {
        int k;
        for (k = 0; k < nmap[i][j]; k++)
          if (map[i][j][k] == m) break;
        if (k == nmap[i][j]) map[i][j][nmap[i][j]++] = m;
        setflag[i][j] = 1;
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   we need to handle Pair::svector special for hybrid/scaled
------------------------------------------------------------------------- */

void PairHybridScaled::init_svector()
{
  // single_extra = list all sub-style single_extra
  // allocate svector

  single_extra = 0;
  for (int m = 0; m < nstyles; m++)
    single_extra += styles[m]->single_extra;

  if (single_extra) {
    delete [] svector;
    svector = new double[single_extra];
  }
}

/* ----------------------------------------------------------------------
   we need to handle Pair::svector special for hybrid/scaled
------------------------------------------------------------------------- */

void PairHybridScaled::copy_svector(int itype, int jtype)
{
  int n=0;
  Pair *this_style;

  // fill svector array.
  // copy data from active styles and use 0.0 for inactive ones
  for (int m = 0; m < nstyles; m++) {
    for (int k = 0; k < nmap[itype][jtype]; ++k) {
      if (m == map[itype][jtype][k]) {
        this_style = styles[m];
      } else {
        this_style = nullptr;
      }
    }
    for (int l = 0; l < styles[m]->single_extra; ++l) {
      if (this_style) {
        svector[n++] = this_style->svector[l];
      } else {
        svector[n++] = 0.0;
      }
    }
  }
}
