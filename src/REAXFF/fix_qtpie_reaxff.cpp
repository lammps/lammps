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
   Contributing authors:
      Efstratios M Kritikos, California Institute of Technology
      (Implemented original version in LAMMMPS Aug 2019)
      Navraj S Lalli, Imperial College London
      (Reimplemented QTPIE as a new fix in LAMMPS Aug 2024 and extended functionality)
------------------------------------------------------------------------- */

#include "fix_qtpie_reaxff.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_efield.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "region.h"
#include "respa.h"
#include "text_file_reader.h"
#include "update.h"

#include "pair_reaxff.h"
#include "reaxff_api.h"

#include <cmath>
#include <cstring>
#include <exception>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double CONV_TO_EV = 14.4;
static constexpr double QSUMSMALL = 0.00001;
static constexpr double ANGSTROM_TO_BOHRRADIUS = 1.8897261259;

/* ---------------------------------------------------------------------- */

FixQtpieReaxFF::FixQtpieReaxFF(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), matvecs(0), pertype_option(nullptr), gauss_file(nullptr)
{
  // this fix returns a global scalar (the number of iterations)
  scalar_flag = 1;
  extscalar = 0;

  // this fix returns a per-atom vector (the effective electronegativity)
  peratom_flag = 1;
  size_peratom_cols = 0;

  imax = 200;
  maxwarn = 1;

  if ((narg < 9) || (narg > 12)) error->all(FLERR,"Illegal fix {} command", style);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix {} command", style);

  swa = utils::numeric(FLERR,arg[4],false,lmp);
  swb = utils::numeric(FLERR,arg[5],false,lmp);
  tolerance = utils::numeric(FLERR,arg[6],false,lmp);
  pertype_option = utils::strdup(arg[7]);
  gauss_file = utils::strdup(arg[8]);

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nowarn") == 0) maxwarn = 0;
    else if (strcmp(arg[iarg],"maxiter") == 0) {
      if (iarg+1 > narg-1)
        error->all(FLERR,"Illegal fix {} command", style);
      imax = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg++;
    } else error->all(FLERR,"Illegal fix {} command", style);
    iarg++;
  }
  shld = nullptr;

  nn = nt = n_cap = 0;
  nmax = 0;
  m_fill = m_cap = 0;
  pack_flag = 0;
  s = nullptr;
  t = nullptr;
  nprev = 4;

  Hdia_inv = nullptr;
  b_s = nullptr;
  chi_eff = nullptr;
  b_t = nullptr;
  b_prc = nullptr;
  b_prm = nullptr;

  // CG

  p = nullptr;
  q = nullptr;
  r = nullptr;
  d = nullptr;

  // H matrix

  H.firstnbr = nullptr;
  H.numnbrs = nullptr;
  H.jlist = nullptr;
  H.val = nullptr;

  comm_forward = comm_reverse = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  reaxff = dynamic_cast<PairReaxFF *>(force->pair_match("^reaxff",0));

  s_hist = t_hist = nullptr;
  atom->add_callback(Atom::GROW);
}

/* ---------------------------------------------------------------------- */

FixQtpieReaxFF::~FixQtpieReaxFF()
{
  if (copymode) return;

  delete[] pertype_option;
  delete[] gauss_file;

  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);

  memory->destroy(s_hist);
  memory->destroy(t_hist);

  FixQtpieReaxFF::deallocate_storage();
  FixQtpieReaxFF::deallocate_matrix();

  memory->destroy(shld);

  memory->destroy(gauss_exp);
  if (!reaxflag) {
    memory->destroy(chi);
    memory->destroy(eta);
    memory->destroy(gamma);
  }
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::post_constructor()
{
  grow_arrays(atom->nmax);
  for (int i = 0; i < atom->nmax; i++)
    for (int j = 0; j < nprev; ++j)
      s_hist[i][j] = t_hist[i][j] = 0;

  pertype_parameters(pertype_option);
}

/* ---------------------------------------------------------------------- */

int FixQtpieReaxFF::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::pertype_parameters(char *arg)
{
  const int nlocal = atom->nlocal;
  const int *mask = atom->mask;
  const int *type = atom->type;
  const int ntypes = atom->ntypes;

  // read gaussian orbital exponents
  memory->create(gauss_exp,ntypes+1,"qtpie/reaxff:gauss_exp");
  if (comm->me == 0) {
    gauss_exp[0] = 0.0;
    try {
      FILE *fp = utils::open_potential(gauss_file, lmp, nullptr);
      if (!fp) throw TokenizerException("Fix qtpie/reaxff: could not open gauss file", gauss_file);
      TextFileReader reader(fp,"qtpie/reaxff gaussian exponents");
      reader.ignore_comments = true;
      for (int i = 1; i <= ntypes; i++) {
        const char *line = reader.next_line();
        if (!line)
          throw TokenizerException("Fix qtpie/reaxff: Incorrect number of atom types in gauss file","");
        ValueTokenizer values(line);

        if (values.count() != 2)
          throw TokenizerException("Fix qtpie/reaxff: Incorrect number of values per line "
                                   "in gauss file",std::to_string(values.count()));

        int itype = values.next_int();
        if ((itype < 1) || (itype > ntypes))
          throw TokenizerException("Fix qtpie/reaxff: Invalid atom type in gauss file",
                                   std::to_string(itype));

        double exp = values.next_double();
        if (exp < 0)
          throw TokenizerException("Fix qtpie/reaxff: Invalid orbital exponent in gauss file",
                                   std::to_string(exp));
        gauss_exp[itype] = exp;
      }
      fclose(fp);
    } catch (std::exception &e) {
      error->one(FLERR,e.what());
    }
  }

  MPI_Bcast(gauss_exp,ntypes+1,MPI_DOUBLE,0,world);

  // define a cutoff distance (in atomic units) beyond which overlap integrals are neglected
  // in calc_chi_eff()
  const double exp_min = find_min_exp(gauss_exp,ntypes+1);
  const int olap_cut = 10; // overlap integrals are neglected if less than pow(10,-olap_cut)
  dist_cutoff = sqrt(2*olap_cut/exp_min*log(10.0));

  // read chi, eta and gamma

  if (utils::strmatch(arg,"^reaxff")) {
    reaxflag = 1;
    Pair *pair = force->pair_match("^reaxff",0);
    if (!pair) error->all(FLERR,"No reaxff pair style for fix qtpie/reaxff");

    int tmp, tmp_all;
    chi = (double *) pair->extract("chi",tmp);
    eta = (double *) pair->extract("eta",tmp);
    gamma = (double *) pair->extract("gamma",tmp);
    if ((chi == nullptr) || (eta == nullptr) || (gamma == nullptr))
      error->all(FLERR, "Fix qtpie/reaxff could not extract qtpie parameters from pair reaxff");
    tmp = tmp_all = 0;
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        if ((chi[type[i]] == 0.0) && (eta[type[i]] == 0.0) && (gamma[type[i]] == 0.0))
          tmp = type[i];
      }
    }
    MPI_Allreduce(&tmp, &tmp_all, 1, MPI_INT, MPI_MAX, world);
    if (tmp_all)
      error->all(FLERR, "No qtpie parameters for atom type {} provided by pair reaxff", tmp_all);
    return;
  } else if (utils::strmatch(arg,"^reax/c")) {
    error->all(FLERR, "Fix qtpie/reaxff keyword 'reax/c' is obsolete; please use 'reaxff'");
  } else if (platform::file_is_readable(arg)) {
    ; // arg is readable file. will read below
  } else {
    error->all(FLERR, "Unknown fix qtpie/reaxff keyword {}", arg);
  }

  reaxflag = 0;

  memory->create(chi,ntypes+1,"qtpie/reaxff:chi");
  memory->create(eta,ntypes+1,"qtpie/reaxff:eta");
  memory->create(gamma,ntypes+1,"qtpie/reaxff:gamma");

  if (comm->me == 0) {
    chi[0] = eta[0] = gamma[0] = 0.0;
    try {
      TextFileReader reader(arg,"qtpie/reaxff parameter");
      reader.ignore_comments = false;
      for (int i = 1; i <= ntypes; i++) {
        const char *line = reader.next_line();
        if (!line)
          throw TokenizerException("Fix qtpie/reaxff: Invalid param file format","");
        ValueTokenizer values(line);

        if (values.count() != 4)
          throw TokenizerException("Fix qtpie/reaxff: Incorrect format of param file","");

        int itype = values.next_int();
        if ((itype < 1) || (itype > ntypes))
          throw TokenizerException("Fix qtpie/reaxff: Invalid atom type in param file",
                                   std::to_string(itype));

        chi[itype] = values.next_double();
        eta[itype] = values.next_double();
        gamma[itype] = values.next_double();
      }
    } catch (std::exception &e) {
      error->one(FLERR,e.what());
    }
  }

  MPI_Bcast(chi,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(eta,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(gamma,ntypes+1,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::allocate_storage()
{
  nmax = atom->nmax;

  memory->create(s,nmax,"qtpie:s");
  memory->create(t,nmax,"qtpie:t");

  memory->create(Hdia_inv,nmax,"qtpie:Hdia_inv");
  memory->create(b_s,nmax,"qtpie:b_s");
  memory->create(chi_eff,nmax,"qtpie:chi_eff");
  vector_atom = chi_eff;
  memory->create(b_t,nmax,"qtpie:b_t");
  memory->create(b_prc,nmax,"qtpie:b_prc");
  memory->create(b_prm,nmax,"qtpie:b_prm");

  int size = nmax;

  memory->create(p,size,"qtpie:p");
  memory->create(q,size,"qtpie:q");
  memory->create(r,size,"qtpie:r");
  memory->create(d,size,"qtpie:d");
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::deallocate_storage()
{
  memory->destroy(s);
  memory->destroy(t);

  memory->destroy(Hdia_inv);
  memory->destroy(b_s);
  memory->destroy(b_t);
  memory->destroy(b_prc);
  memory->destroy(b_prm);
  memory->destroy(chi_eff);

  memory->destroy(p);
  memory->destroy(q);
  memory->destroy(r);
  memory->destroy(d);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::reallocate_storage()
{
  deallocate_storage();
  allocate_storage();
  init_storage();
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::allocate_matrix()
{
  int i,ii;
  bigint m;

  int mincap;
  double safezone;

  if (reaxflag && reaxff) {
    mincap = reaxff->api->system->mincap;
    safezone = reaxff->api->system->safezone;
  } else {
    mincap = REAX_MIN_CAP;
    safezone = REAX_SAFE_ZONE;
  }

  n_cap = MAX((int)(atom->nlocal * safezone), mincap);

  // determine the total space for the H matrix

  m = 0;
  for (ii = 0; ii < nn; ii++) {
    i = ilist[ii];
    m += numneigh[i];
  }
  bigint m_cap_big = (bigint)MAX(m * safezone, mincap * REAX_MIN_NBRS);
  if (m_cap_big > MAXSMALLINT)
    error->one(FLERR,"Too many neighbors in fix {}",style);
  m_cap = m_cap_big;

  H.n = n_cap;
  H.m = m_cap;
  memory->create(H.firstnbr,n_cap,"qtpie:H.firstnbr");
  memory->create(H.numnbrs,n_cap,"qtpie:H.numnbrs");
  memory->create(H.jlist,m_cap,"qtpie:H.jlist");
  memory->create(H.val,m_cap,"qtpie:H.val");
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::deallocate_matrix()
{
  memory->destroy(H.firstnbr);
  memory->destroy(H.numnbrs);
  memory->destroy(H.jlist);
  memory->destroy(H.val);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::reallocate_matrix()
{
  deallocate_matrix();
  allocate_matrix();
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix {} requires atom attribute q", style);

  if (group->count(igroup) == 0)
    error->all(FLERR,"Fix {} group has no atoms", style);

  // compute net charge and print warning if too large
  double qsum_local = 0.0, qsum = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit)
      qsum_local += atom->q[i];
  }
  MPI_Allreduce(&qsum_local,&qsum,1,MPI_DOUBLE,MPI_SUM,world);

  if ((comm->me == 0) && (fabs(qsum) > QSUMSMALL))
    error->warning(FLERR,"Fix {} group is not charge neutral, net charge = {:.8}", style, qsum);

  // get pointer to fix efield if present. there may be at most one instance of fix efield in use.
  efield = nullptr;
  auto fixes = modify->get_fix_by_style("^efield");
  if (fixes.size() == 1) efield = dynamic_cast<FixEfield *>(fixes.front());
  else if (fixes.size() > 1)
    error->all(FLERR, "There may be only one fix efield instance used with fix {}", style);

  // ensure that fix efield is properly initialized before accessing its data and check some settings
  if (efield) {
    efield->init();
    if (strcmp(update->unit_style,"real") != 0)
      error->all(FLERR,"Must use unit_style real with fix {} and external fields", style);

    if (efield->groupbit != 1){ // if efield is not applied to all atoms
      error->all(FLERR,"Must use group id all for fix efield when using fix {}", style);
    }

    if (efield->region){ // if efield is not applied to all atoms
      error->all(FLERR,"Keyword region not supported for fix efield when using fix {}", style);
    }

    if (efield->varflag == FixEfield::ATOM && efield->pstyle != FixEfield::ATOM)
      error->all(FLERR,"Atom-style external electric field requires atom-style "
                       "potential variable when used with fix {}", style);
  }

  // we need a half neighbor list w/ Newton off
  // built whenever re-neighboring occurs

  neighbor->add_request(this, NeighConst::REQ_NEWTON_OFF);

  init_shielding();
  init_taper();

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

double FixQtpieReaxFF::compute_scalar()
{
  return matvecs/2.0;
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::init_shielding()
{
  int i,j;
  int ntypes;

  ntypes = atom->ntypes;
  if (shld == nullptr)
    memory->create(shld,ntypes+1,ntypes+1,"qtpie:shielding");

  for (i = 1; i <= ntypes; ++i)
    for (j = 1; j <= ntypes; ++j)
      shld[i][j] = pow(gamma[i] * gamma[j], -1.5);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::init_taper()
{
  double d7, swa2, swa3, swb2, swb3;

  if (fabs(swa) > 0.01 && comm->me == 0)
    error->warning(FLERR,"Fix qtpie/reaxff has non-zero lower Taper radius cutoff");
  if (swb < 0)
    error->all(FLERR, "Fix qtpie/reaxff has negative upper Taper radius cutoff");
  else if (swb < 5 && comm->me == 0)
    error->warning(FLERR,"Fix qtpie/reaxff has very low Taper radius cutoff");

  d7 = pow(swb - swa, 7);
  swa2 = SQR(swa);
  swa3 = CUBE(swa);
  swb2 = SQR(swb);
  swb3 = CUBE(swb);

  Tap[7] =  20.0 / d7;
  Tap[6] = -70.0 * (swa + swb) / d7;
  Tap[5] =  84.0 * (swa2 + 3.0*swa*swb + swb2) / d7;
  Tap[4] = -35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3) / d7;
  Tap[3] = 140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3) / d7;
  Tap[2] =-210.0 * (swa3*swb2 + swa2*swb3) / d7;
  Tap[1] = 140.0 * swa3 * swb3 / d7;
  Tap[0] = (-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 -
            7.0*swa*swb3*swb3 + swb3*swb3*swb) / d7;
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::setup_pre_force(int vflag)
{
  if (reaxff) {
    nn = reaxff->list->inum;
    nt = reaxff->list->inum + reaxff->list->gnum;
    ilist = reaxff->list->ilist;
    numneigh = reaxff->list->numneigh;
    firstneigh = reaxff->list->firstneigh;
  } else {
    nn = list->inum;
    nt = list->inum + list->gnum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }

  deallocate_storage();
  allocate_storage();

  init_storage();

  deallocate_matrix();
  allocate_matrix();

  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::min_setup_pre_force(int vflag)
{
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::init_storage()
{
  calc_chi_eff();

  for (int ii = 0; ii < nn; ii++) {
    int i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      Hdia_inv[i] = 1. / eta[atom->type[i]];
      b_s[i] = -chi_eff[i];
      b_t[i] = -1.0;
      b_prc[i] = 0;
      b_prm[i] = 0;
      s[i] = t[i] = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  int n = atom->nlocal;

  if (reaxff) {
    nn = reaxff->list->inum;
    nt = reaxff->list->inum + reaxff->list->gnum;
    ilist = reaxff->list->ilist;
    numneigh = reaxff->list->numneigh;
    firstneigh = reaxff->list->firstneigh;
  } else {
    nn = list->inum;
    nt = list->inum + list->gnum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
  }

  // grow arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) reallocate_storage();
  if (n > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  calc_chi_eff();

  init_matvec();

  matvecs_s = CG(b_s, s);       // CG on s - parallel
  matvecs_t = CG(b_t, t);       // CG on t - parallel
  matvecs = matvecs_s + matvecs_t;

  calculate_Q();
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::pre_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::init_matvec()
{
  /* fill-in H matrix */
  compute_H();

  int ii, i;

  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {

      /* init pre-conditioner for H and init solution vectors */
      Hdia_inv[i] = 1. / eta[atom->type[i]];
      b_s[i]      = -chi_eff[i];
      b_t[i]      = -1.0;

      /* quadratic extrapolation for s & t from previous solutions */
      t[i] = t_hist[i][2] + 3 * (t_hist[i][0] - t_hist[i][1]);

      /* cubic extrapolation for s & t from previous solutions */
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
    }
  }

  pack_flag = 2;
  comm->forward_comm(this); //Dist_vector(s);
  pack_flag = 3;
  comm->forward_comm(this); //Dist_vector(t);
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::compute_H()
{
  int jnum;
  int i, j, ii, jj, flag;
  double dx, dy, dz, r_sqr;
  constexpr double EPSILON = 0.0001;

  int *type = atom->type;
  tagint *tag = atom->tag;
  double **x = atom->x;
  int *mask = atom->mask;

  // fill in the H matrix
  m_fill = 0;
  r_sqr = 0;
  for (ii = 0; ii < nn; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      jlist = firstneigh[i];
      jnum = numneigh[i];
      H.firstnbr[i] = m_fill;

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
        r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

        flag = 0;
        if (r_sqr <= SQR(swb)) {
          if (j < atom->nlocal) flag = 1;
          else if (tag[i] < tag[j]) flag = 1;
          else if (tag[i] == tag[j]) {
            if (dz > EPSILON) flag = 1;
            else if (fabs(dz) < EPSILON) {
              if (dy > EPSILON) flag = 1;
              else if (fabs(dy) < EPSILON && dx > EPSILON)
                flag = 1;
            }
          }
        }

        if (flag) {
          H.jlist[m_fill] = j;
          H.val[m_fill] = calculate_H(sqrt(r_sqr), shld[type[i]][type[j]]);
          m_fill++;
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }

  if (m_fill >= H.m)
    error->all(FLERR,"Fix qtpie/reaxff H matrix size has been exceeded: m_fill={} H.m={}\n",
               m_fill, H.m);
}

/* ---------------------------------------------------------------------- */

double FixQtpieReaxFF::calculate_H(double r, double gamma)
{
  double Taper, denom;

  Taper = Tap[7] * r + Tap[6];
  Taper = Taper * r + Tap[5];
  Taper = Taper * r + Tap[4];
  Taper = Taper * r + Tap[3];
  Taper = Taper * r + Tap[2];
  Taper = Taper * r + Tap[1];
  Taper = Taper * r + Tap[0];

  denom = r * r * r + gamma;
  denom = pow(denom,1.0/3.0);

  return Taper * CONV_TO_EV / denom;
}

/* ---------------------------------------------------------------------- */

int FixQtpieReaxFF::CG(double *b, double *x)
{
  int  i, j;
  double tmp, alpha, beta, b_norm;
  double sig_old, sig_new;

  int jj;

  pack_flag = 1;
  sparse_matvec(&H, x, q);
  comm->reverse_comm(this); //Coll_Vector(q);

  vector_sum(r , 1.,  b, -1., q, nn);

  for (jj = 0; jj < nn; ++jj) {
    j = ilist[jj];
    if (atom->mask[j] & groupbit)
      d[j] = r[j] * Hdia_inv[j]; //pre-condition
  }

  b_norm = parallel_norm(b, nn);
  sig_new = parallel_dot(r, d, nn);

  for (i = 1; i < imax && sqrt(sig_new) / b_norm > tolerance; ++i) {
    comm->forward_comm(this); //Dist_vector(d);
    sparse_matvec(&H, d, q);
    comm->reverse_comm(this); //Coll_vector(q);

    tmp = parallel_dot(d, q, nn);
    alpha = sig_new / tmp;

    vector_add(x, alpha, d, nn);
    vector_add(r, -alpha, q, nn);

    // pre-conditioning
    for (jj = 0; jj < nn; ++jj) {
      j = ilist[jj];
      if (atom->mask[j] & groupbit)
        p[j] = r[j] * Hdia_inv[j];
    }

    sig_old = sig_new;
    sig_new = parallel_dot(r, p, nn);

    beta = sig_new / sig_old;
    vector_sum(d, 1., p, beta, d, nn);
  }

  if ((i >= imax) && maxwarn && (comm->me == 0))
    error->warning(FLERR, "Fix qtpie/reaxff CG convergence failed after {} iterations at step {}",
                   i,update->ntimestep);
  return i;
}


/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::sparse_matvec(sparse_matrix *A, double *x, double *b)
{
  int i, j, itr_j;
  int ii;

  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      b[i] = eta[atom->type[i]] * x[i];
  }

  int nall = atom->nlocal + atom->nghost;
  for (i = atom->nlocal; i < nall; ++i)
      b[i] = 0;

  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      for (itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
        j = A->jlist[itr_j];
        b[i] += A->val[itr_j] * x[j];
        b[j] += A->val[itr_j] * x[i];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::calculate_Q()
{
  int i, k;
  double u, s_sum, t_sum;
  double *q = atom->q;

  int ii;

  s_sum = parallel_vector_acc(s, nn);
  t_sum = parallel_vector_acc(t, nn);
  u = s_sum / t_sum;

  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      q[i] = s[i] - u * t[i];

      /* backup s & t */
      for (k = nprev-1; k > 0; --k) {
        s_hist[i][k] = s_hist[i][k-1];
        t_hist[i][k] = t_hist[i][k-1];
      }
      s_hist[i][0] = s[i];
      t_hist[i][0] = t[i];
    }
  }

  pack_flag = 4;
  comm->forward_comm(this); //Dist_vector(atom->q);
}

/* ---------------------------------------------------------------------- */

int FixQtpieReaxFF::pack_forward_comm(int n, int *list, double *buf,
                                  int /*pbc_flag*/, int * /*pbc*/)
{
  int m;

  if (pack_flag == 1)
    for (m = 0; m < n; m++) buf[m] = d[list[m]];
  else if (pack_flag == 2)
    for (m = 0; m < n; m++) buf[m] = s[list[m]];
  else if (pack_flag == 3)
    for (m = 0; m < n; m++) buf[m] = t[list[m]];
  else if (pack_flag == 4)
    for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1)
    for (m = 0, i = first; m < n; m++, i++) d[i] = buf[m];
  else if (pack_flag == 2)
    for (m = 0, i = first; m < n; m++, i++) s[i] = buf[m];
  else if (pack_flag == 3)
    for (m = 0, i = first; m < n; m++, i++) t[i] = buf[m];
  else if (pack_flag == 4)
    for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
}

/* ---------------------------------------------------------------------- */

int FixQtpieReaxFF::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  for (m = 0, i = first; m < n; m++, i++) buf[m] = q[i];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::unpack_reverse_comm(int n, int *list, double *buf)
{
  for (int m = 0; m < n; m++) q[list[m]] += buf[m];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQtpieReaxFF::memory_usage()
{
  double bytes;

  bytes = (double)atom->nmax*nprev*2 * sizeof(double); // s_hist & t_hist
  bytes += (double)atom->nmax*11 * sizeof(double); // storage
  bytes += (double)n_cap*2 * sizeof(int); // matrix...
  bytes += (double)m_cap * sizeof(int);
  bytes += (double)m_cap * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

void FixQtpieReaxFF::grow_arrays(int nmax)
{
  memory->grow(s_hist,nmax,nprev,"qtpie:s_hist");
  memory->grow(t_hist,nmax,nprev,"qtpie:t_hist");
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

void FixQtpieReaxFF::copy_arrays(int i, int j, int /*delflag*/)
{
  for (int m = 0; m < nprev; m++) {
    s_hist[j][m] = s_hist[i][m];
    t_hist[j][m] = t_hist[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixQtpieReaxFF::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nprev; m++) buf[m] = s_hist[i][m];
  for (int m = 0; m < nprev; m++) buf[nprev+m] = t_hist[i][m];
  return nprev*2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixQtpieReaxFF::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nprev; m++) s_hist[nlocal][m] = buf[m];
  for (int m = 0; m < nprev; m++) t_hist[nlocal][m] = buf[nprev+m];
  return nprev*2;
}

/* ---------------------------------------------------------------------- */

double FixQtpieReaxFF::parallel_norm(double *v, int n)
{
  int  i;
  double my_sum, norm_sqr;

  int ii;

  my_sum = 0.0;
  norm_sqr = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_sum += SQR(v[i]);
  }

  MPI_Allreduce(&my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world);

  return sqrt(norm_sqr);
}

/* ---------------------------------------------------------------------- */

double FixQtpieReaxFF::parallel_dot(double *v1, double *v2, int n)
{
  int  i;
  double my_dot, res;

  int ii;

  my_dot = 0.0;
  res = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_dot += v1[i] * v2[i];
  }

  MPI_Allreduce(&my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

double FixQtpieReaxFF::parallel_vector_acc(double *v, int n)
{
  int  i;
  double my_acc, res;

  int ii;

  my_acc = 0.0;
  res = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_acc += v[i];
  }

  MPI_Allreduce(&my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::vector_sum(double* dest, double c, double* v,
                                double d, double* y, int k)
{
  int kk;

  for (--k; k>=0; --k) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] = c * v[kk] + d * y[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::vector_add(double* dest, double c, double* v, int k)
{
  int kk;

  for (--k; k>=0; --k) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] += c * v[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQtpieReaxFF::calc_chi_eff()
{
  memset(&chi_eff[0],0,atom->nmax*sizeof(double));

  const auto x = (const double * const *)atom->x;
  const int *type = atom->type;

  double dist,overlap,sum_n,sum_d,expa,expb,chia,chib,phia,phib,p,m;
  int i,j;

  // check ghost atoms are stored up to the distance cutoff for overlap integrals
  const double comm_cutoff = MAX(neighbor->cutneighmax,comm->cutghostuser);
  if(comm_cutoff < dist_cutoff/ANGSTROM_TO_BOHRRADIUS) {
    error->all(FLERR,"comm cutoff = {} Angstrom is smaller than distance cutoff = {} Angstrom "
               "for overlap integrals in {}. Increase comm cutoff with comm_modify",
               comm_cutoff, dist_cutoff/ANGSTROM_TO_BOHRRADIUS, style);
  }

  // efield energy is in real units of kcal/mol, factor needed for conversion to eV
  const double qe2f = force->qe2f;
  const double factor = 1.0/qe2f;

  if (efield) {
    if (efield->varflag != FixEfield::CONSTANT)
      efield->update_efield_variables();
  }

  // compute chi_eff for each local atom
  for (i = 0; i < nn; i++) {
    expa = gauss_exp[type[i]];
    chia = chi[type[i]];
    if (efield) {
      if (efield->varflag != FixEfield::ATOM) {
        phia = -factor*(x[i][0]*efield->ex  + x[i][1]*efield->ey + x[i][2]*efield->ez);
      } else { // atom-style potential from FixEfield
        phia = efield->efield[i][3];
      }
    }

    sum_n = 0.0;
    sum_d = 0.0;

    for (j = 0; j < nt; j++) {
      dist = distance(x[i],x[j])*ANGSTROM_TO_BOHRRADIUS; // in atomic units

      if (dist < dist_cutoff) {
        expb = gauss_exp[type[j]];
        chib = chi[type[j]];

        // overlap integral of two normalised 1s Gaussian type orbitals
        p = expa + expb;
        m = expa * expb / p;
        overlap = pow((4.0*m/p),0.75) * exp(-m*dist*dist);

        if (efield) {
          if (efield->varflag != FixEfield::ATOM) {
            phib = -factor*(x[j][0]*efield->ex  + x[j][1]*efield->ey + x[j][2]*efield->ez);
          } else { // atom-style potential from FixEfield
            phib = efield->efield[j][3];
          }
          sum_n += (chia - chib + phia - phib) * overlap;
        } else {
          sum_n += (chia - chib) * overlap;
        }
        sum_d += overlap;
      }
    }

    chi_eff[i] = sum_n / sum_d;
  }
}

/* ---------------------------------------------------------------------- */

double FixQtpieReaxFF::find_min_exp(const double *array, const int array_length)
{
  // index of first gaussian orbital exponent is 1
  double exp_min = array[1];
  for (int i = 2; i < array_length; i++)
  {
    if (array[i] < exp_min)
      exp_min = array[i];
  }
  return exp_min;
}

/* ---------------------------------------------------------------------- */

double FixQtpieReaxFF::distance(const double *posa, const double *posb)
{
  double dx, dy, dz;
  dx = posb[0] - posa[0];
  dy = posb[1] - posa[1];
  dz = posb[2] - posa[2];
  return sqrt(dx*dx + dy*dy + dz*dz);
}
