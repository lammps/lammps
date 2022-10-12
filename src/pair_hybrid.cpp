// clang-format off

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

#include "pair_hybrid.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "respa.h"
#include "suffix.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairHybrid::PairHybrid(LAMMPS *lmp) :
    Pair(lmp), styles(nullptr), cutmax_style(nullptr), keywords(nullptr), multiple(nullptr),
    nmap(nullptr), map(nullptr), special_lj(nullptr), special_coul(nullptr), compute_tally(nullptr)
{
  nstyles = 0;

  outerflag = 0;
  respaflag = 0;
}

/* ---------------------------------------------------------------------- */

PairHybrid::~PairHybrid()
{
  if (nstyles > 0) {
    for (int m = 0; m < nstyles; m++) {
      delete styles[m];
      delete[] keywords[m];
      delete[] special_lj[m];
      delete[] special_coul[m];
    }
  }
  delete[] styles;
  delete[] cutmax_style;
  delete[] keywords;
  delete[] multiple;

  delete[] special_lj;
  delete[] special_coul;
  delete[] compute_tally;

  delete[] svector;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    memory->destroy(nmap);
    memory->destroy(map);
  }
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

void PairHybrid::compute(int eflag, int vflag)
{
  int i,j,m,n;

  // check if no_virial_fdotr_compute is set and global component of
  //   incoming vflag = VIRIAL_FDOTR
  // if so, reset vflag as if global component were VIRIAL_PAIR
  // necessary since one or more sub-styles cannot compute virial as F dot r

  if (no_virial_fdotr_compute && (vflag & VIRIAL_FDOTR))
    vflag = VIRIAL_PAIR | (vflag & ~VIRIAL_FDOTR);

  ev_init(eflag,vflag);

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
    respa = dynamic_cast<Respa *>(update->integrate);
    if (respa->nhybrid_styles > 0) respaflag = 1;
  }

  for (m = 0; m < nstyles; m++) {

    set_special(m);

    if (!respaflag || (respaflag && respa->hybrid_compute[m])) {

      // invoke compute() unless compute flag is turned off or
      // outerflag is set and sub-style has a compute_outer() method

      if (styles[m]->compute_flag == 0) continue;
      if (outerflag && styles[m]->respa_enable)
        styles[m]->compute_outer(eflag,vflag_substyle);
      else styles[m]->compute(eflag,vflag_substyle);
    }

    restore_special(saved_special);

    // jump to next sub-style if r-RESPA does not want global accumulated data

    if (respaflag && !respa->tally_global) continue;

    if (eflag_global) {
      eng_vdwl += styles[m]->eng_vdwl;
      eng_coul += styles[m]->eng_coul;
    }
    if (vflag_global) {
      for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
    }
    if (eflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double *eatom_substyle = styles[m]->eatom;
      for (i = 0; i < n; i++) eatom[i] += eatom_substyle[i];
    }
    if (vflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double **vatom_substyle = styles[m]->vatom;
      for (i = 0; i < n; i++)
        for (j = 0; j < 6; j++)
          vatom[i][j] += vatom_substyle[i][j];
    }

    // substyles may be CENTROID_SAME or CENTROID_AVAIL

    if (cvflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      if (styles[m]->centroidstressflag == CENTROID_AVAIL) {
        double **cvatom_substyle = styles[m]->cvatom;
        for (i = 0; i < n; i++)
          for (j = 0; j < 9; j++)
            cvatom[i][j] += cvatom_substyle[i][j];
      } else {
        double **vatom_substyle = styles[m]->vatom;
        for (i = 0; i < n; i++) {
          for (j = 0; j < 6; j++) {
            cvatom[i][j] += vatom_substyle[i][j];
          }
          for (j = 6; j < 9; j++) {
            cvatom[i][j] += vatom_substyle[i][j-3];
          }
        }
      }
    }

  }

  delete[] saved_special;

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairHybrid::finish()
{
  for (int m = 0; m < nstyles; m++) styles[m]->finish();
}

/* ---------------------------------------------------------------------- */

void PairHybrid::add_tally_callback(Compute *ptr)
{
  for (int m = 0; m < nstyles; m++)
    if (compute_tally[m]) styles[m]->add_tally_callback(ptr);
}

/* ---------------------------------------------------------------------- */

void PairHybrid::del_tally_callback(Compute *ptr)
{
  for (int m = 0; m < nstyles; m++)
    if (compute_tally[m]) styles[m]->del_tally_callback(ptr);
}

/* ---------------------------------------------------------------------- */

void PairHybrid::compute_inner()
{
  for (int m = 0; m < nstyles; m++)
    if (styles[m]->respa_enable) styles[m]->compute_inner();
}

/* ---------------------------------------------------------------------- */

void PairHybrid::compute_middle()
{
  for (int m = 0; m < nstyles; m++)
    if (styles[m]->respa_enable) styles[m]->compute_middle();
}

/* ---------------------------------------------------------------------- */

void PairHybrid::compute_outer(int eflag, int vflag)
{
  outerflag = 1;
  compute(eflag,vflag);
  outerflag = 0;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairHybrid::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");

  memory->create(nmap,n+1,n+1,"pair:nmap");
  memory->create(map,n+1,n+1,nstyles,"pair:map");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      nmap[i][j] = 0;
}

/* ----------------------------------------------------------------------
   create one pair style for each arg in list
------------------------------------------------------------------------- */

void PairHybrid::settings(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");
  if (lmp->kokkos && !utils::strmatch(force->pair_style,"^hybrid.*/kk$"))
    error->all(FLERR,"Must use pair_style {}/kk with Kokkos",
                                 force->pair_style);

  // delete old lists, since cannot just change settings

  if (nstyles > 0) {
    for (int m = 0; m < nstyles; m++) {
      delete styles[m];
      delete[] keywords[m];
      delete[] special_lj[m];
      delete[] special_coul[m];
    }
    delete[] styles;
    delete[] keywords;
    delete[] multiple;
    delete[] special_lj;
    delete[] special_coul;
    delete[] compute_tally;
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

  styles = new Pair *[narg];
  cutmax_style = new double[narg];
  memset(cutmax_style, 0.0, narg*sizeof(double));
  keywords = new char *[narg];
  multiple = new int[narg];

  special_lj = new double*[narg];
  special_coul = new double*[narg];
  compute_tally = new int[narg];

  // allocate each sub-style
  // allocate uses suffix, but don't store suffix version in keywords,
  //   else syntax in coeff() will not match
  // call settings() with set of args that are not pair style names
  // use force->pair_map to determine which args these are

  int iarg,jarg,dummy;

  iarg = 0;
  nstyles = 0;
  while (iarg < narg) {
    if (utils::strmatch(arg[iarg],"^hybrid"))
      error->all(FLERR,"Pair style hybrid cannot have hybrid as an argument");
    if (strcmp(arg[iarg],"none") == 0)
      error->all(FLERR,"Pair style hybrid cannot have none as an argument");

    styles[nstyles] = force->new_pair(arg[iarg],1,dummy);
    keywords[nstyles] = force->store_style(arg[iarg],0);
    special_lj[nstyles] = special_coul[nstyles] = nullptr;
    compute_tally[nstyles] = 1;

    // determine list of arguments for pair style settings
    // by looking for the next known pair style name.

    jarg = iarg + 1;
    while ((jarg < narg)
           && !force->pair_map->count(arg[jarg])
           && !lmp->match_style("pair", arg[jarg])) jarg++;

    styles[nstyles]->settings(jarg-iarg-1,&arg[iarg+1]);
    iarg = jarg;
    nstyles++;
  }

  // multiple[i] = 1 to M if sub-style used multiple times, else 0

  int num_tip4p = 0, num_coul = 0; // count sub-styles with tip4p and coulomb

  for (int i = 0; i < nstyles; i++) {
    int count = 0;
    for (int j = 0; j < nstyles; j++) {
      if (strcmp(keywords[j],keywords[i]) == 0) count++;
      if (j == i) multiple[i] = count;
    }
    if (count == 1) multiple[i] = 0;

    if (utils::strmatch(keywords[i],"/tip4p/")) ++num_tip4p;
    if (utils::strmatch(keywords[i],"/coul/")
        || utils::strmatch(keywords[i],"^comb")
        || utils::strmatch(keywords[i],"^reax/c")) ++num_coul;
  }

  if ((num_tip4p > 1) && (comm->me == 0))
    error->warning(FLERR,"Using multiple tip4p sub-styles can result in "
                   "inconsistent calculation of coulomb interactions");

  if ((num_tip4p > 0) && (num_coul > 0) && (comm->me == 0))
    error->warning(FLERR,"Using a tip4p sub-style with other sub-styles "
                   "that include coulomb interactions can result in "
                   "inconsistent calculation of the coulomb interactions");

  // set pair flags from sub-style flags

  flags();
}

/* ----------------------------------------------------------------------
   set top-level pair flags from sub-style flags
------------------------------------------------------------------------- */

void PairHybrid::flags()
{
  int m;

  // set comm_forward, comm_reverse, comm_reverse_off to max of any sub-style

  for (m = 0; m < nstyles; m++) {
    if (styles[m]) comm_forward = MAX(comm_forward,styles[m]->comm_forward);
    if (styles[m]) comm_reverse = MAX(comm_reverse,styles[m]->comm_reverse);
    if (styles[m]) comm_reverse_off = MAX(comm_reverse_off,
                                          styles[m]->comm_reverse_off);
  }

  // single_enable = 1 if all sub-styles are set
  // respa_enable = 1 if all sub-styles are set
  // manybody_flag = 1 if any sub-style is set
  // born_matrix_enable = 1 if all sub-styles are set
  // no_virial_fdotr_compute = 1 if any sub-style is set
  // ghostneigh = 1 if any sub-style is set
  // ewaldflag, pppmflag, msmflag, dipoleflag, dispersionflag, tip4pflag = 1
  //   if any sub-style is set
  // compute_flag = 1 if any sub-style is set

  single_enable = 0;
  compute_flag = 0;
  respa_enable = 0;
  restartinfo = 0;
  born_matrix_enable = 0;

  for (m = 0; m < nstyles; m++) {
    if (styles[m]->single_enable) ++single_enable;
    if (styles[m]->respa_enable) ++respa_enable;
    if (styles[m]->restartinfo) ++restartinfo;
    if (styles[m]->born_matrix_enable) ++born_matrix_enable;
    if (styles[m]->manybody_flag) manybody_flag = 1;
    if (styles[m]->no_virial_fdotr_compute) no_virial_fdotr_compute = 1;
    if (styles[m]->ghostneigh) ghostneigh = 1;
    if (styles[m]->ewaldflag) ewaldflag = 1;
    if (styles[m]->pppmflag) pppmflag = 1;
    if (styles[m]->msmflag) msmflag = 1;
    if (styles[m]->dipoleflag) dipoleflag = 1;
    if (styles[m]->spinflag) spinflag = 1;
    if (styles[m]->dispersionflag) dispersionflag = 1;
    if (styles[m]->tip4pflag) tip4pflag = 1;
    if (styles[m]->compute_flag) compute_flag = 1;
    if (styles[m]->finitecutflag) finitecutflag = 1;
  }
  single_enable = (single_enable == nstyles) ? 1 : 0;
  respa_enable = (respa_enable == nstyles) ? 1 : 0;
  restartinfo = (restartinfo == nstyles) ? 1 : 0;
  born_matrix_enable = (born_matrix_enable == nstyles) ? 1 : 0;
  init_svector();

  // set centroidstressflag for pair hybrid
  // set to CENTROID_NOTAVAIL if any substyle is NOTAVAIL
  // else set to CENTROID_AVAIL if any substyle is AVAIL
  // else set to CENTROID_SAME if all substyles are SAME

  centroidstressflag = CENTROID_SAME;

  for (m = 0; m < nstyles; m++) {
    if (styles[m]->centroidstressflag == CENTROID_NOTAVAIL)
      centroidstressflag = CENTROID_NOTAVAIL;
    if (centroidstressflag == CENTROID_SAME &&
        styles[m]->centroidstressflag == CENTROID_AVAIL)
      centroidstressflag = CENTROID_AVAIL;
  }
}

/* ----------------------------------------------------------------------
   initialize Pair::svector array
------------------------------------------------------------------------- */

void PairHybrid::init_svector()
{
  // single_extra = list all sub-style single_extra
  // allocate svector

  single_extra = 0;
  for (int m = 0; m < nstyles; m++)
    single_extra = MAX(single_extra,styles[m]->single_extra);

  if (single_extra) {
    delete[] svector;
    svector = new double[single_extra];
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHybrid::coeff(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // 3rd arg = pair sub-style name
  // 4th arg = pair sub-style index if name used multiple times
  // allow for "none" as valid sub-style name

  int multflag = 0;
  int m;

  for (m = 0; m < nstyles; m++) {
    multflag = 0;
    if (strcmp(arg[2],keywords[m]) == 0) {
      if (multiple[m]) {
        multflag = 1;
        if (narg < 4) error->all(FLERR,"Incorrect args for pair coefficients");
        if (multiple[m] == utils::inumeric(FLERR,arg[3],false,lmp)) break;
        else continue;
      } else break;
    }
  }

  int none = 0;
  if (m == nstyles) {
    if (strcmp(arg[2],"none") == 0) none = 1;
    else error->all(FLERR,"Pair coeff for hybrid has invalid style: {}",arg[2]);
  }

  // move 1st/2nd args to 2nd/3rd args
  // if multflag: move 1st/2nd args to 3rd/4th args
  // just copy ptrs, since arg[] points into original input line

  arg[2+multflag] = arg[1];
  arg[1+multflag] = arg[0];

  // invoke sub-style coeff() starting with 1st remaining arg

  if (!none) styles[m]->coeff(narg-1-multflag,arg+1+multflag);

  // if sub-style only allows one pair coeff call (with * * and type mapping)
  // then unset setflag/map assigned to that style before setting it below
  // in case pair coeff for this sub-style is being called for 2nd time

  if (!none && styles[m]->one_coeff) {
    if ((strcmp(arg[0],"*") != 0) || (strcmp(arg[1],"*") != 0))
      error->all(FLERR,"Incorrect args for pair coefficients");
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (nmap[i][j] && map[i][j][0] == m) {
          setflag[i][j] = 0;
          nmap[i][j] = 0;
        }
  }

  // set setflag and which type pairs map to which sub-style
  // if sub-style is none: set hybrid setflag, wipe out map
  // else: set hybrid setflag & map only if substyle setflag is set
  //       previous mappings are wiped out

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (none) {
        setflag[i][j] = 1;
        nmap[i][j] = 0;
        count++;
      } else if (styles[m]->setflag[i][j]) {
        setflag[i][j] = 1;
        nmap[i][j] = 1;
        map[i][j][0] = m;
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairHybrid::init_style()
{
  int i,m,itype,jtype,used,istyle,skip;

  // error if a sub-style is not used

  int ntypes = atom->ntypes;

  for (istyle = 0; istyle < nstyles; istyle++) {
    used = 0;
    for (itype = 1; itype <= ntypes; itype++)
      for (jtype = itype; jtype <= ntypes; jtype++)
        for (m = 0; m < nmap[itype][jtype]; m++)
          if (map[itype][jtype][m] == istyle) used = 1;
    if (used == 0) error->all(FLERR,"Pair hybrid sub-style is not used");
  }

  // The GPU library uses global data for each pair style, so the
  // same style must not be used multiple times

  for (istyle = 0; istyle < nstyles; istyle++) {
    bool is_gpu = styles[istyle]->suffix_flag & Suffix::GPU;
    if (multiple[istyle] && is_gpu)
      error->all(FLERR,"GPU package styles must not be used multiple times");
  }

  // check if special_lj/special_coul overrides are compatible

  for (istyle = 0; istyle < nstyles; istyle++) {
    if (special_lj[istyle]) {
      for (i = 1; i < 4; ++i) {
        if (((force->special_lj[i] == 0.0) || (force->special_lj[i] == 1.0))
            && (force->special_lj[i] != special_lj[istyle][i]))
          error->all(FLERR,"Pair_modify special lj 1-{} setting for pair hybrid substyle {} "
                     "incompatible with global special_bonds setting", i+1, keywords[istyle]);
      }
    }

    if (special_coul[istyle]) {
      for (i = 1; i < 4; ++i) {
        if (((force->special_coul[i] == 0.0) || (force->special_coul[i] == 1.0))
            && (force->special_coul[i] != special_coul[istyle][i]))
          error->all(FLERR,"Pair_modify special coul 1-{} setting for pair hybrid substyle {} "
                     "incompatible with global special_bonds setting", i+1, keywords[istyle]);
      }
    }
  }

  // check beyond contact (set during pair coeff) before init style
  for (istyle = 0; istyle < nstyles; istyle++)
    if (styles[istyle]->beyond_contact) beyond_contact = 1;

  // each sub-style makes its neighbor list request(s)

  for (istyle = 0; istyle < nstyles; istyle++) styles[istyle]->init_style();

  // create skip lists inside each pair neigh request
  // any kind of list can have its skip flag set in this loop

  for (auto &request : neighbor->get_pair_requests()) {

    // istyle = associated sub-style for the request

    for (istyle = 0; istyle < nstyles; istyle++)
      if (styles[istyle] == request->get_requestor()) break;

    // allocate iskip and ijskip
    // initialize so as to skip all pair types
    // set ijskip = 0 if type pair matches any entry in sub-style map
    // set ijskip = 0 if mixing will assign type pair to this sub-style
    //   will occur if type pair is currently unassigned
    //   and both I,I and J,J are assigned to single sub-style
    //   and sub-style for both I,I and J,J match istyle
    // set iskip = 1 only if all ijskip for itype are 1

    int *iskip = new int[ntypes+1];
    int **ijskip;
    memory->create(ijskip,ntypes+1,ntypes+1,"pair_hybrid:ijskip");

    for (itype = 1; itype <= ntypes; itype++)
      for (jtype = 1; jtype <= ntypes; jtype++)
        ijskip[itype][jtype] = 1;

    for (itype = 1; itype <= ntypes; itype++)
      for (jtype = itype; jtype <= ntypes; jtype++) {
        for (m = 0; m < nmap[itype][jtype]; m++)
          if (map[itype][jtype][m] == istyle)
            ijskip[itype][jtype] = ijskip[jtype][itype] = 0;
        if (nmap[itype][jtype] == 0 &&
            nmap[itype][itype] == 1 && map[itype][itype][0] == istyle &&
            nmap[jtype][jtype] == 1 && map[jtype][jtype][0] == istyle)
          ijskip[itype][jtype] = ijskip[jtype][itype] = 0;
      }

    for (itype = 1; itype <= ntypes; itype++) {
      iskip[itype] = 1;
      for (jtype = 1; jtype <= ntypes; jtype++)
        if (ijskip[itype][jtype] == 0) iskip[itype] = 0;
    }

    // if any skipping occurs
    // set request->skip and copy iskip and ijskip into request
    // else delete iskip and ijskip
    // no skipping if pair style assigned to all type pairs

    skip = 0;
    for (itype = 1; itype <= ntypes; itype++)
      for (jtype = 1; jtype <= ntypes; jtype++)
        if (ijskip[itype][jtype] == 1) skip = 1;

    if (skip) {
      request->set_skip(iskip, ijskip);
    } else {
      delete[] iskip;
      memory->destroy(ijskip);
    }
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairHybrid::init_one(int i, int j)
{
  // if I,J is not set explicitly:
  // perform mixing only if I,I sub-style = J,J sub-style
  // also require I,I and J,J are both assigned to single sub-style

  if (setflag[i][j] == 0) {
    if (nmap[i][i] != 1 || nmap[j][j] != 1 || map[i][i][0] != map[j][j][0])
      error->one(FLERR,"All pair coeffs are not set");
    nmap[i][j] = 1;
    map[i][j][0] = map[i][i][0];
  }

  // call init/mixing for all sub-styles of I,J
  // set cutsq in sub-style just as Pair::init() does via call to init_one()
  // set cutghost for I,J and J,I just as sub-style does
  // sum tail corrections for I,J
  // return max cutoff of all sub-styles assigned to I,J
  // if no sub-styles assigned to I,J (pair_coeff none), cutmax = 0.0 returned

  double cutmax = 0.0;
  cutghost[i][j] = cutghost[j][i] = 0.0;
  if (tail_flag) etail_ij = ptail_ij = 0.0;

  nmap[j][i] = nmap[i][j];

  for (int k = 0; k < nmap[i][j]; k++) {
    map[j][i][k] = map[i][j][k];
    double cut = styles[map[i][j][k]]->init_one(i,j);
    if (styles[map[i][j][k]]->did_mix) did_mix = true;
    styles[map[i][j][k]]->cutsq[i][j] = styles[map[i][j][k]]->cutsq[j][i] = cut*cut;
    if (styles[map[i][j][k]]->ghostneigh)
      cutghost[i][j] = cutghost[j][i] = MAX(cutghost[i][j],styles[map[i][j][k]]->cutghost[i][j]);
    if (tail_flag) {
      etail_ij += styles[map[i][j][k]]->etail_ij;
      ptail_ij += styles[map[i][j][k]]->ptail_ij;
    }
    cutmax = MAX(cutmax,cut);

    int istyle;
    for (istyle = 0; istyle < nstyles; istyle++)
      if (styles[istyle] == styles[map[i][j][k]]) break;

    if (styles[istyle]->trim_flag) {

      if (cut > cutmax_style[istyle]) {
        cutmax_style[istyle] = cut;

        for (auto &request : neighbor->get_pair_requests()) {
          if (styles[istyle] == request->get_requestor() && styles[istyle]->trim_flag) {
            request->set_cutoff(cutmax_style[istyle]);
            break;
          }
        }
      }
    }
  }

  return cutmax;
}

/* ----------------------------------------------------------------------
   invoke setup for each sub-style
------------------------------------------------------------------------- */

void PairHybrid::setup()
{
  for (int m = 0; m < nstyles; m++) styles[m]->setup();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairHybrid::write_restart(FILE *fp)
{
  fwrite(&nstyles,sizeof(int),1,fp);

  // each sub-style writes its settings, but no coeff info

  fwrite(compute_tally,sizeof(int),nstyles,fp);

  int n;
  for (int m = 0; m < nstyles; m++) {
    n = strlen(keywords[m]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(keywords[m],sizeof(char),n,fp);
    styles[m]->write_restart_settings(fp);
    // write out per style special settings, if present
    n = (special_lj[m] == nullptr) ? 0 : 1;
    fwrite(&n,sizeof(int),1,fp);
    if (n) fwrite(special_lj[m],sizeof(double),4,fp);
    n = (special_coul[m] == nullptr) ? 0 : 1;
    fwrite(&n,sizeof(int),1,fp);
    if (n) fwrite(special_coul[m],sizeof(double),4,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairHybrid::read_restart(FILE *fp)
{
  int me = comm->me;
  if (me == 0) utils::sfread(FLERR,&nstyles,sizeof(int),1,fp,nullptr,error);
  MPI_Bcast(&nstyles,1,MPI_INT,0,world);

  // allocate list of sub-styles

  delete[] styles;
  delete[] keywords;
  delete[] multiple;
  delete[] special_lj;
  delete[] special_coul;
  delete[] compute_tally;

  styles = new Pair*[nstyles];
  cutmax_style = new double[nstyles];
  memset(cutmax_style, 0.0, nstyles*sizeof(double));
  keywords = new char*[nstyles];
  multiple = new int[nstyles];

  special_lj = new double*[nstyles];
  special_coul = new double*[nstyles];
  compute_tally = new int[nstyles];

  // each sub-style is created via new_pair()
  // each reads its settings, but no coeff info

  if (me == 0) utils::sfread(FLERR,compute_tally,sizeof(int),nstyles,fp,nullptr,error);
  MPI_Bcast(compute_tally,nstyles,MPI_INT,0,world);

  int n,dummy;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) utils::sfread(FLERR,&n,sizeof(int),1,fp,nullptr,error);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    keywords[m] = new char[n];
    if (me == 0) utils::sfread(FLERR,keywords[m],sizeof(char),n,fp,nullptr,error);
    MPI_Bcast(keywords[m],n,MPI_CHAR,0,world);
    styles[m] = force->new_pair(keywords[m],1,dummy);
    styles[m]->read_restart_settings(fp);
    // read back per style special settings, if present
    special_lj[m] = special_coul[m] = nullptr;
    if (me == 0) utils::sfread(FLERR,&n,sizeof(int),1,fp,nullptr,error);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n > 0) {
      special_lj[m] = new double[4];
      if (me == 0) utils::sfread(FLERR,special_lj[m],sizeof(double),4,fp,nullptr,error);
      MPI_Bcast(special_lj[m],4,MPI_DOUBLE,0,world);
    }
    if (me == 0) utils::sfread(FLERR,&n,sizeof(int),1,fp,nullptr,error);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n > 0) {
      special_coul[m] = new double[4];
      if (me == 0) utils::sfread(FLERR,special_coul[m],sizeof(double),4,fp,nullptr,error);
      MPI_Bcast(special_coul[m],4,MPI_DOUBLE,0,world);
    }
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
   call sub-style to compute single interaction
   error if sub-style does not support single() call
   since overlay could have multiple sub-styles, sum results explicitly
------------------------------------------------------------------------- */

double PairHybrid::single(int i, int j, int itype, int jtype,
                          double rsq, double factor_coul, double factor_lj,
                          double &fforce)
{
  if (nmap[itype][jtype] == 0)
    error->one(FLERR,"Invoked pair single on pair style none");

  double fone;
  fforce = 0.0;
  double esum = 0.0;

  for (int m = 0; m < nmap[itype][jtype]; m++) {
    if (rsq < styles[map[itype][jtype][m]]->cutsq[itype][jtype]) {
      if (styles[map[itype][jtype][m]]->single_enable == 0)
        error->one(FLERR,"Pair hybrid sub-style does not support single call");

      if ((special_lj[map[itype][jtype][m]] != nullptr) ||
          (special_coul[map[itype][jtype][m]] != nullptr))
        error->one(FLERR,"Pair hybrid single calls do not support"
                   " per sub-style special bond values");

      esum += styles[map[itype][jtype][m]]->
        single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fone);
      fforce += fone;
    }
  }

  if (single_extra) copy_svector(itype,jtype);
  return esum;
}

/* ----------------------------------------------------------------------
   call sub-style to compute born matrix interaction
   error if sub-style does not support born_matrix call
   since overlay could have multiple sub-styles, sum results explicitly
------------------------------------------------------------------------- */

void PairHybrid::born_matrix(int i, int j, int itype, int jtype, double rsq,
                             double factor_coul, double factor_lj,
                             double &dupair, double &du2pair)
{
  if (nmap[itype][jtype] == 0)
    error->one(FLERR,"Invoked pair born_matrix on pair style none");

  double du, du2;
  dupair = du2pair = 0.0;

  for (int m = 0; m < nmap[itype][jtype]; m++) {
    if (rsq < styles[map[itype][jtype][m]]->cutsq[itype][jtype]) {
      if (styles[map[itype][jtype][m]]->born_matrix_enable == 0)
        error->one(FLERR,"Pair hybrid sub-style does not support born_matrix call");

      if ((special_lj[map[itype][jtype][m]] != nullptr) ||
          (special_coul[map[itype][jtype][m]] != nullptr))
        error->one(FLERR,"Pair hybrid born_matrix calls do not support"
                   " per sub-style special bond values");

      du = du2 = 0.0;
      styles[map[itype][jtype][m]]->born_matrix(i,j,itype,jtype,rsq,factor_coul,factor_lj,du,du2);
      dupair += du;
      du2pair += du2;
    }
  }
}

/* ----------------------------------------------------------------------
   copy Pair::svector data
------------------------------------------------------------------------- */

void PairHybrid::copy_svector(int itype, int jtype)
{
  memset(svector,0,single_extra*sizeof(double));

  // there is only one style in pair style hybrid for a pair of atom types
  Pair *this_style = styles[map[itype][jtype][0]];

  for (int l = 0; this_style->single_extra; ++l) {
    svector[l] = this_style->svector[l];
  }
}

/* ----------------------------------------------------------------------
   modify parameters of the pair style and its sub-styles
------------------------------------------------------------------------- */

void PairHybrid::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal pair_modify command");

  // if 1st keyword is pair, apply other keywords to one sub-style

  if (strcmp(arg[0],"pair") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal pair_modify command");
    int m;
    for (m = 0; m < nstyles; m++)
      if (strcmp(arg[1],keywords[m]) == 0) break;
    if (m == nstyles) error->all(FLERR,"Unknown pair_modify hybrid sub-style: {}",arg[1]);
    int iarg = 2;

    if (multiple[m]) {
      if (narg < 3) error->all(FLERR,"Illegal pair_modify command");
      int multiflag = utils::inumeric(FLERR,arg[2],false,lmp);
      for (m = 0; m < nstyles; m++)
        if (strcmp(arg[1],keywords[m]) == 0 && multiflag == multiple[m]) break;
      if (m == nstyles)
        error->all(FLERR,"Unknown pair_modify hybrid sub-style: {}",arg[1]);
      iarg = 3;
    }

    // keywords "special" and "compute/tally" have to be listed directly
    // after "pair" but can be given multiple times

again:

    if (iarg < narg && strcmp(arg[iarg],"special") == 0) {
      if (narg < iarg+5) error->all(FLERR,"Illegal pair_modify special command");
      modify_special(m,narg-iarg,&arg[iarg+1]);
      iarg += 5;
      goto again;
    }

    // if 2nd keyword (after pair) is compute/tally:
    // set flag to register TALLY computes accordingly

    if (iarg < narg && strcmp(arg[iarg],"compute/tally") == 0) {
      if (narg < iarg+2) error->all(FLERR,"Illegal pair_modify compute/tally command");
      compute_tally[m] = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
      goto again;
    }

    // apply the remaining keywords to the base pair style itself and
    // the sub-style except for "pair" and "special" or "compute/tally"
    // and their arguments. the former is important for some keywords
    // like "tail" or "compute"

    if (narg-iarg > 0) {
      Pair::modify_params(narg-iarg,&arg[iarg]);
      styles[m]->modify_params(narg-iarg,&arg[iarg]);
    }

  // apply all keywords to pair hybrid itself and every sub-style

  } else {
    Pair::modify_params(narg,arg);
    for (int m = 0; m < nstyles; m++) styles[m]->modify_params(narg,arg);
  }

  // reset global compute_flag since there may have been changes
  // to any of the substyles
  compute_flag = 0;
  for (int m = 0; m < nstyles; m++)
    if (styles[m]->compute_flag) compute_flag = 1;
}

/* ----------------------------------------------------------------------
   store a local per pair style override for special_lj and special_coul
------------------------------------------------------------------------- */

void PairHybrid::modify_special(int m, int /*narg*/, char **arg)
{
  double special[4];
  int i;

  special[0] = 1.0;
  special[1] = utils::numeric(FLERR,arg[1],false,lmp);
  special[2] = utils::numeric(FLERR,arg[2],false,lmp);
  special[3] = utils::numeric(FLERR,arg[3],false,lmp);

  if (styles[m]->suffix_flag & (Suffix::INTEL|Suffix::GPU))
    error->all(FLERR,"Pair_modify special not compatible with suffix version of hybrid substyle");

  if (strcmp(arg[0],"lj/coul") == 0) {
    if (!special_lj[m]) special_lj[m] = new double[4];
    if (!special_coul[m]) special_coul[m] = new double[4];
    for (i = 0; i < 4; ++i)
      special_lj[m][i] = special_coul[m][i] = special[i];

  } else if (strcmp(arg[0],"lj") == 0) {
    if (!special_lj[m]) special_lj[m] = new double[4];
    for (i = 0; i < 4; ++i)
      special_lj[m][i] = special[i];

  } else if (strcmp(arg[0],"coul") == 0) {
    if (!special_coul[m]) special_coul[m] = new double[4];
    for (i = 0; i < 4; ++i)
      special_coul[m][i] = special[i];

  } else error->all(FLERR,"Illegal pair_modify special command");
}

/* ----------------------------------------------------------------------
   override global special bonds settings with per substyle values
------------------------------------------------------------------------- */

void PairHybrid::set_special(int m)
{
  int i;
  if (special_lj[m])
    for (i = 0; i < 4; ++i) force->special_lj[i] = special_lj[m][i];
  if (special_coul[m])
    for (i = 0; i < 4; ++i) force->special_coul[i] = special_coul[m][i];
}

/* ----------------------------------------------------------------------
   store global special settings
------------------------------------------------------------------------- */

double * PairHybrid::save_special()
{
  auto saved = new double[8];

  for (int i = 0; i < 4; ++i) {
    saved[i] = force->special_lj[i];
    saved[i+4] = force->special_coul[i];
  }
  return saved;
}

/* ----------------------------------------------------------------------
   restore global special settings from saved data
------------------------------------------------------------------------- */

void PairHybrid::restore_special(double *saved)
{
  for (int i = 0; i < 4; ++i) {
    force->special_lj[i] = saved[i];
    force->special_coul[i] = saved[i+4];
  }
}

/* ----------------------------------------------------------------------
   extract a ptr to a particular quantity stored by pair
   pass request thru to sub-styles
   return first non-nullptr result except for cut_coul request
   for cut_coul, insure all non-nullptr results are equal since required by Kspace
------------------------------------------------------------------------- */

void *PairHybrid::extract(const char *str, int &dim)
{
  void *cutptr = nullptr;
  void *ptr;
  double cutvalue = 0.0;
  int couldim = -1;

  for (int m = 0; m < nstyles; m++) {
    ptr = styles[m]->extract(str,dim);
    if (ptr && strcmp(str,"cut_coul") == 0) {
      if (couldim != -1 && dim != couldim)
        error->all(FLERR, "Coulomb styles of pair hybrid sub-styles do not match");
      auto p_newvalue = (double *) ptr;
      double newvalue = *p_newvalue;
      if (cutptr && (newvalue != cutvalue))
        error->all(FLERR, "Coulomb cutoffs of pair hybrid sub-styles do not match");
      if (dim == 0) {
        cutptr = ptr;
        cutvalue = newvalue;
      }
      couldim = dim;
    } else if (ptr) return ptr;
  }

  if (strcmp(str,"cut_coul") == 0) return cutptr;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

void PairHybrid::reset_dt()
{
  for (int m = 0; m < nstyles; m++) styles[m]->reset_dt();
}

/* ----------------------------------------------------------------------
   check if itype,jtype maps to sub-style
------------------------------------------------------------------------- */

int PairHybrid::check_ijtype(int itype, int jtype, char *substyle)
{
  for (int m = 0; m < nmap[itype][jtype]; m++)
    if (strcmp(keywords[map[itype][jtype][m]],substyle) == 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if substyles calculate self-interaction range of particle
------------------------------------------------------------------------- */

double PairHybrid::atom2cut(int i)
{
  double temp, cut;

  cut = 0.0;
  for (int m = 0; m < nstyles; m++) {
    if (styles[m]->finitecutflag) {
      temp = styles[m]->atom2cut(i);
      if (temp > cut) cut = temp;
    }
  }
  return cut;
}

/* ----------------------------------------------------------------------
   check if substyles calculate maximum interaction range for two finite particles
------------------------------------------------------------------------- */

double PairHybrid::radii2cut(double r1, double r2)
{
  double temp, cut;

 cut = 0.0;
  for (int m = 0; m < nstyles; m++) {
    if (styles[m]->finitecutflag) {
      temp = styles[m]->radii2cut(r1,r2);
      if (temp > cut) cut = temp;
    }
  }
  return cut;
}

/* ----------------------------------------------------------------------
   memory usage of each sub-style
------------------------------------------------------------------------- */

double PairHybrid::memory_usage()
{
  double bytes = (double)maxeatom * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  bytes += (double)maxcvatom*9 * sizeof(double);
  for (int m = 0; m < nstyles; m++) bytes += styles[m]->memory_usage();
  return bytes;
}
