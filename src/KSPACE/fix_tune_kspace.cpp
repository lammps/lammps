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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "fix_tune_kspace.h"

#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "suffix.h"
#include "timer.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <limits>
#include <utility>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static constexpr double GOLD = 1.618034;

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixTuneKspace::FixTuneKspace(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  global_freq = 1;
  firststep = 0;
  niter = 0;
  niter_adjust_rcut = 0;
  keep_bracketing = true;
  first_brent_pass = true;
  converged = false;
  need_fd2_brent = false;

  // parse arguments

  if (narg < 4) utils::missing_cmd_args(FLERR, "fix tune/kspace", error);
  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix tune/kspace interval {}", nevery);

  msmflag = espflag = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "msm") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix tune/kspace msm", error);
      msmflag = utils::logical(FLERR, arg[iarg+1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "esp") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix tune/kspace esp", error);
      espflag = utils::logical(FLERR, arg[iarg+1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix tune/kspace keyword: {}", arg[iarg]);
    }
  }

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

int FixTuneKspace::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTuneKspace::init()
{
  if (!force->kspace)
    error->all(FLERR,"Cannot use fix tune/kspace without a kspace style");
  if (!force->pair)
    error->all(FLERR,"Cannot use fix tune/kspace without a pair style");
  if (strncmp(force->pair_style,"hybrid",6) == 0)
    error->all(FLERR,"Cannot use fix tune/kspace with a hybrid pair style");
  if (force->kspace->dispersionflag)
    error->all(FLERR,"Cannot use fix tune/kspace with long-range dispersion");
  if (force->kspace->tip4pflag)
    error->all(FLERR,"Cannot use fix tune/kspace with TIP4P water");
  if (force->kspace->dipoleflag)
    error->all(FLERR,"Cannot use fix tune/kspace with dipole long-range solver");

  // the MSM and ESP styles are not (yet) supported in combination with the OPENMP package

  if (msmflag || espflag) {
    bool ompflag = (force->pair->suffix_flag & Suffix::OMP) != 0;
    if (lmp->suffix_enable) {
      if (lmp->suffix && (strcmp(lmp->suffix, "omp") == 0)) ompflag = true;
      if (lmp->suffix2 && (strcmp(lmp->suffix2, "omp") == 0)) ompflag = true;
    }
    if (ompflag && msmflag)
      error->all(FLERR, "Fix tune/kspace keyword 'msm yes' is not supported with the OPENMP package");
    if (ompflag && espflag)
      error->all(FLERR, "Fix tune/kspace keyword 'esp yes' is not supported with the OPENMP package");
  }

  store_old_kspace_settings();
  double old_acc = force->kspace->accuracy/force->kspace->two_charge_force;
  acc_str = std::to_string(old_acc);

  int itmp;
  auto *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  pair_cut_coul = *p_cutoff;
}

/* ----------------------------------------------------------------------
   perform dynamic kspace parameter optimization
------------------------------------------------------------------------- */

void FixTuneKspace::pre_exchange()
{
  if (!nevery) return;
  if (!force->kspace) return;
  if (!force->pair) return;
  if (next_reneighbor != update->ntimestep) return;
  next_reneighbor = update->ntimestep + nevery;

  double time = get_timing_info();
  niter++;

  // determine the kspace styles to test.  Ewald and PPPM use the same pair style,
  // MSM and ESP are only tested if requested and if a matching pair style exists.

  if (niter == 1) {
    store_old_kspace_settings();
    test_styles = {"ewald", "pppm"};
    Info info(lmp);
    if (msmflag && info.has_style("pair", base_pair_style + "/msm")) test_styles.emplace_back("msm");
    if (espflag && info.has_style("pair", base_pair_style + "/esp")) test_styles.emplace_back("esp");
    test_times.assign(test_styles.size(), 0.0);
  }

  // the timing is for the style that was tested since the previous invocation

  const int ntest = test_styles.size();
  if ((niter > 1) && (niter <= ntest + 1)) test_times[niter - 2] = time;

  if (niter <= ntest) {
    // test the next kspace style
    switch_kspace_style(test_styles[niter - 1]);
  } else if (niter == ntest + 1) {
    // switch to the fastest kspace style
    int best = 0;
    for (int i = 0; i < ntest; ++i) {
      if (comm->me == 0) utils::logmesg(lmp, "{}_time = {}\n", test_styles[i], test_times[i]);
      if (test_times[i] < test_times[best]) best = i;
    }
    switch_kspace_style(test_styles[best]);
  } else {
    adjust_rcut(time);
  }

  last_spcpu = timer->elapsed(Timer::TOTAL);
}

/* ----------------------------------------------------------------------
   figure out CPU time per timestep since last time checked
------------------------------------------------------------------------- */

double FixTuneKspace::get_timing_info()
{
  double dvalue;
  double new_cpu;
  int new_step = update->ntimestep;

  if (firststep == 0) {
    new_cpu = 0.0;
    dvalue = 0.0;
    firststep = 1;
  } else {
    new_cpu = timer->elapsed(Timer::TOTAL);
    double cpu_diff = new_cpu - last_spcpu;
    int step_diff = new_step - last_step;
    if (step_diff > 0.0) dvalue = cpu_diff/step_diff;
    else dvalue = 0.0;
  }

  last_step = new_step;
  last_spcpu = new_cpu;

  return dvalue;
}

/* ----------------------------------------------------------------------
   store old kspace settings: style, accuracy, order, etc
------------------------------------------------------------------------- */

void FixTuneKspace::store_old_kspace_settings()
{
  // store the style names without accelerator suffix, since the
  // styles are re-created with the suffix applied, if needed

  kspace_style = utils::strip_style_suffix(force->kspace_style, lmp);
  pair_style = utils::strip_style_suffix(force->pair_style, lmp);

  std::size_t found;
  if (std::string::npos != (found = pair_style.rfind("/long")))
    base_pair_style = pair_style.substr(0,found);
  else if (std::string::npos != (found = pair_style.rfind("/msm")))
    base_pair_style = pair_style.substr(0,found);
  else if (std::string::npos != (found = pair_style.rfind("/esp")))
    base_pair_style = pair_style.substr(0,found);
  else base_pair_style = pair_style;

  old_differentiation_flag = force->kspace->differentiation_flag;
  old_slabflag = force->kspace->slabflag;
  old_slab_volfactor = force->kspace->slab_volfactor;
}

/* ----------------------------------------------------------------------
   switch to the given kspace style and the matching pair style
------------------------------------------------------------------------- */

void FixTuneKspace::switch_kspace_style(const std::string &new_kspace_style)
{
  store_old_kspace_settings();
  if (new_kspace_style == "msm")
    pair_style = base_pair_style + "/msm";
  else if (new_kspace_style == "esp")
    pair_style = base_pair_style + "/esp";
  else
    pair_style = base_pair_style + "/long";
  kspace_style = new_kspace_style;
  update_pair_style(pair_style, pair_cut_coul);
  update_kspace_style(kspace_style, acc_str);
}

/* ----------------------------------------------------------------------
   update the pair style if necessary, preserving the settings
------------------------------------------------------------------------- */

void FixTuneKspace::update_pair_style(const std::string &new_pair_style,
                                      double pair_cut_coul)
{
  int itmp;
  auto *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  *p_cutoff = pair_cut_coul;

  // check to see if we need to change pair styles
  // the current style name may include an accelerator suffix
  if (utils::strip_style_suffix(new_pair_style, lmp) ==
      utils::strip_style_suffix(force->pair_style, lmp))
    return;

  // create a temporary file to store current pair settings
  FILE *p_pair_settings_file;
  p_pair_settings_file = tmpfile();
  force->pair->write_restart(p_pair_settings_file);
  rewind(p_pair_settings_file);
  if (comm->me == 0)
    utils::logmesg(lmp,"Creating new pair style: {}\n",new_pair_style);

  // delete old pair style and create new one
  force->create_pair(new_pair_style,1);

  // restore current pair settings from temporary file
  force->pair->read_restart(p_pair_settings_file);

  auto *pcutoff = (double *) force->pair->extract("cut_coul",itmp);
  double current_cutoff = *pcutoff;
  if (comm->me == 0)
    utils::logmesg(lmp,"Coulomb cutoff for real space: {}\n",current_cutoff);

  // close temporary file
  fclose(p_pair_settings_file);
}

/* ----------------------------------------------------------------------
   update the kspace style if necessary
------------------------------------------------------------------------- */

void FixTuneKspace::update_kspace_style(const std::string &new_kspace_style,
                                        const std::string &new_acc_str)
{
  // delete old kspace style and create new one

  auto *tmp_acc_str = (char *)new_acc_str.c_str();
  force->create_kspace(new_kspace_style,1);
  force->kspace->settings(1,&tmp_acc_str);
  force->kspace->differentiation_flag = old_differentiation_flag;
  force->kspace->slabflag = old_slabflag;
  force->kspace->slab_volfactor = old_slab_volfactor;

  // initialize new kspace style, pair style, molecular styles

  force->init();

  // set up grid

  force->kspace->reset_grid();

  // the force styles have been replaced, so fix omp must update
  // its pointer to the style that does the force reduction

  for (auto &ifix : modify->get_fix_by_style("^OMP$")) ifix->init();

  // re-init neighbor list
  // probably only needed when redefining the pair style
  // should happen after pair->init() to get pair style
  //   neighbor list request registered

  neighbor->init();

  // Re-init computes to update pointers to virials, etc.

  for (const auto &icompute : modify->get_compute_list()) icompute->init();
}

/* ----------------------------------------------------------------------
   find the optimal real space coulomb cutoff
------------------------------------------------------------------------- */

void FixTuneKspace::adjust_rcut(double time)
{
  if (utils::strmatch(force->kspace_style,"^msm")) return;
  if (converged) return;

  constexpr double TINY = 1.0e-20;

  // get the current cutoff
  int itmp;
  auto *p_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  double current_cutoff = *p_cutoff;
  if (comm->me == 0)
    utils::logmesg(lmp,"Old Coulomb cutoff for real space: {}\n",current_cutoff);

  // use Brent's method from Numerical Recipes to find optimal real space cutoff

  // first time through, get ax_brent and fa_brent, and adjust cutoff
  if (keep_bracketing) {
    if (niter_adjust_rcut == 0) {
      pair_cut_coul /= 2;
    } else if (niter_adjust_rcut == 1) {
      ax_brent = current_cutoff;
      fa_brent = time;
      pair_cut_coul *= 2;

    // second time through, get bx_brent and fb_brent, and adjust cutoff
    } else if (niter_adjust_rcut == 2) {
      bx_brent = current_cutoff;
      fb_brent = time;
      if (fb_brent > fa_brent) {
        std::swap(ax_brent,bx_brent);
        std::swap(fb_brent,fa_brent);
        pair_cut_coul /= 4;
      } else {
        pair_cut_coul *= 2;
      }

    // third time through, get cx_brent and fc_brent, and adjust cutoff if needed
    } else if (niter_adjust_rcut == 3) {
      cx_brent = current_cutoff;
      fc_brent = time;
      if (fc_brent > fb_brent) keep_bracketing = false;
      else {
        double r = (bx_brent - ax_brent)*(fb_brent - fc_brent);
        double q = (bx_brent - cx_brent)*(fb_brent - fa_brent);
        dx_brent = bx_brent - ((bx_brent - cx_brent)*q - (bx_brent - ax_brent)*r)/
         (2.0*SIGN(MAX(fabs(q - r),TINY),q - r));
        pair_cut_coul = dx_brent;
      }

    // after third time through, bracket the minimum, and adjust cutoff
    } else if (niter_adjust_rcut > 3) {
      dx_brent = current_cutoff;
      if (need_fd2_brent) fd2_brent = time;
      else fd_brent = time;
      mnbrak();
      pair_cut_coul = dx_brent;
    }
  }

  if (!keep_bracketing) {
    dx_brent = current_cutoff;
    fd_brent = time;
    if (first_brent_pass) brent0();
    else brent2();
    brent1();
    pair_cut_coul = dx_brent;
  }

  niter_adjust_rcut++;

  if (pair_cut_coul <= 0.0) pair_cut_coul = fabs(MIN(ax_brent,MIN(bx_brent,(MIN(cx_brent,dx_brent))))/2.0) + TINY;

  if (pair_cut_coul != pair_cut_coul)
    error->all(FLERR,"Bad real space Coulomb cutoff in fix tune/kspace");

  // change the cutoff to pair_cut_coul
  *p_cutoff = pair_cut_coul;

  // report the new cutoff
  auto *new_cutoff = (double *) force->pair->extract("cut_coul",itmp);
  current_cutoff = *new_cutoff;
  if (comm->me == 0)
    utils::logmesg(lmp,"Adjusted Coulomb cutoff for real space: {}\n", current_cutoff);

  store_old_kspace_settings();
  update_pair_style(pair_style,pair_cut_coul);
  update_kspace_style(kspace_style,acc_str);
}

/* ----------------------------------------------------------------------
   bracket a minimum using parabolic extrapolation
------------------------------------------------------------------------- */

void FixTuneKspace::mnbrak()
{
  const double GLIMIT = 100.0, TINY = 1.0e-20;
  double r,q;
  r = (bx_brent - ax_brent)*(fb_brent - fc_brent);
  q = (bx_brent - cx_brent)*(fb_brent - fa_brent);
  dx_brent = bx_brent - ((bx_brent - cx_brent)*q - (bx_brent - ax_brent)*r)/
   (2.0*SIGN(MAX(fabs(q - r),TINY),q - r));
  dxlim = bx_brent + GLIMIT*(cx_brent - bx_brent);

  if ((bx_brent - dx_brent)*(dx_brent - cx_brent) > 0.0) {
    if (fd_brent < fc_brent) {
      ax_brent = bx_brent;
      bx_brent = dx_brent;
      fa_brent = fb_brent;
      fb_brent = fd_brent;
      keep_bracketing = false;
      return;
    } else if (fd_brent > fb_brent) {
      cx_brent = dx_brent;
      fc_brent = fd_brent;
      keep_bracketing = false;
      return;
    }
    dx_brent = cx_brent + GOLD*(cx_brent - bx_brent);
    if (need_fd2_brent) {
      fd_brent = fd2_brent;
      need_fd2_brent = false;
    } else {
      need_fd2_brent = true;
      return;
    }
  } else if ((cx_brent - dx_brent)*(dx_brent - dxlim) > 0.0) {
    if (fd_brent < fc_brent) {
      if (need_fd2_brent) {
        need_fd2_brent = false;
      } else {
        need_fd2_brent = true;
        dx_brent += GOLD*(dx_brent - cx_brent);
        return;
      }
      shft3(bx_brent,cx_brent,dx_brent,dx_brent + GOLD*(dx_brent - cx_brent));
      shft3(fb_brent,fc_brent,fd_brent,fd2_brent);
    }
  } else if ((dx_brent - dxlim)*(dxlim - cx_brent) >= 0.0) {
    dx_brent = dxlim;
    if (need_fd2_brent) {
      fd_brent = fd2_brent;
      need_fd2_brent = false;
    } else {
      need_fd2_brent = true;
      return;
    }
  } else {
    dx_brent = cx_brent + GOLD*(cx_brent - bx_brent);
    if (need_fd2_brent) {
      fd_brent = fd2_brent;
      need_fd2_brent = false;
    } else {
      need_fd2_brent = true;
      return;
    }
  }
  shft3(ax_brent,bx_brent,cx_brent,dx_brent);
  shft3(fa_brent,fb_brent,fc_brent,fd_brent);
}

/* ----------------------------------------------------------------------
   Brent's method from Numerical Recipes
------------------------------------------------------------------------- */

void FixTuneKspace::brent0()
{
  a_brent=(ax_brent < cx_brent ? ax_brent : cx_brent);
  b_brent=(ax_brent > cx_brent ? ax_brent : cx_brent);
  x_brent=w_brent=v_brent=bx_brent;
  fw_brent=fv_brent=fx_brent=fb_brent;
}

/* ----------------------------------------------------------------------
   Brent's method from Numerical Recipes
------------------------------------------------------------------------- */

void FixTuneKspace::brent1()
{
  constexpr double CGOLD=0.3819660;
  const double ZEPS=std::numeric_limits<double>::epsilon()*1.0e-3;
  double d=0.0,etemp;
  double p,q,r,tol1,tol2,xm;
  double e=0.0;
  double tol=0.001;

  xm=0.5*(a_brent+b_brent);
  tol2=2.0*(tol1=tol*fabs(x_brent)+ZEPS);
  if (fabs(x_brent-xm) <= (tol2-0.5*(b_brent-a_brent))) {
    converged = true;
    dx_brent = x_brent;
    return;
  }
  if (fabs(e) > tol1) {
    r=(x_brent-w_brent)*(fx_brent-fv_brent);
    q=(x_brent-v_brent)*(fx_brent-fw_brent);
    p=(x_brent-v_brent)*q-(x_brent-w_brent)*r;
    q=2.0*(q-r);
    if (q > 0.0) p = -p;
    q=fabs(q);
    etemp=e;
    e=d;
    if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a_brent-x_brent) || p >= q*(b_brent-x_brent))
      d=CGOLD*(e=(x_brent >= xm ? a_brent-x_brent : b_brent-x_brent));
    else {
      d=p/q;
      dx_brent=x_brent+d;
      if (dx_brent-a_brent < tol2 || b_brent-dx_brent < tol2)
        d=SIGN(tol1,xm-x_brent);
    }
  } else {
    d=CGOLD*(e=(x_brent >= xm ? a_brent-x_brent : b_brent-x_brent));
  }
  dx_brent=(fabs(d) >= tol1 ? x_brent+d : x_brent+SIGN(tol1,d));

  first_brent_pass = false;
}

/* ----------------------------------------------------------------------
   Brent's method from Numerical Recipes
------------------------------------------------------------------------- */

void FixTuneKspace::brent2()
{
  if (fd_brent <= fx_brent) {
    if (dx_brent >= x_brent) a_brent=x_brent; else b_brent=x_brent;
    shft3(v_brent,w_brent,x_brent,dx_brent);
    shft3(fv_brent,fw_brent,fx_brent,fd_brent);
  } else {
    if (dx_brent < x_brent) a_brent=dx_brent; else b_brent=dx_brent;
    if (fd_brent <= fw_brent || w_brent == x_brent) {
      v_brent=w_brent;
      w_brent=dx_brent;
      fv_brent=fw_brent;
      fw_brent=fd_brent;
    } else if (fd_brent <= fv_brent || v_brent == x_brent || v_brent == w_brent) {
      v_brent=dx_brent;
      fv_brent=fd_brent;
    }
  }
}

