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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "fix_nve_spin.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_langevin_spin.h"
#include "fix_precession_spin.h"
#include "fix_setforce_spin.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair_hybrid.h"
#include "pair_spin.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_nve_spin[] =
  "fix nve/spin command: https://doi.org/10.1016/j.jcp.2018.06.042\n\n"
  "@article{tranchida2018massively,\n"
  "title={Massively Parallel Symplectic Algorithm for Coupled Magnetic Spin "
  "   Dynamics and Molecular Dynamics},\n"
  "author={Tranchida, J and Plimpton, S J and Thibaudeau, P and Thompson, A P},\n"
  "journal={Journal of Computational Physics},\n"
  "volume={372},\n"
  "pages={406--425},\n"
  "year={2018},\n"
  "publisher={Elsevier}\n"
  "doi={10.1016/j.jcp.2018.06.042}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixNVESpin::FixNVESpin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  pair(nullptr), spin_pairs(nullptr), locklangevinspin(nullptr),
  locksetforcespin(nullptr), lockprecessionspin(nullptr),
  rsec(nullptr), stack_head(nullptr), stack_foot(nullptr),
  backward_stacks(nullptr), forward_stacks(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_nve_spin);

  if (narg < 4) error->all(FLERR,"Illegal fix/nve/spin command");

  time_integrate = 1;
  lattice_flag = 1;
  nlocal_max = 0;
  npairs = 0;
  npairspin = 0;

  // test nprec
  nprecspin = nlangspin = nsetspin = 0;

  // checking if map array or hash is defined

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR,"Fix nve/spin requires an atom map, see atom_modify");

  // define sector_flag

  sector_flag = (comm->nprocs > 1) ? 1 : 0;

  // no self-consistent iterations of the single spin updates by default

  sciter = 0;

  // defining lattice_flag

  // changing the lattice option, from (yes,no) -> (moving,frozen)
  // for now, (yes,no) still works (to avoid user's confusions).

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"lattice") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix/nve/spin command");
      const std::string latarg = arg[iarg+1];
      if ((latarg == "no") || (latarg == "off") || (latarg == "false") || (latarg == "frozen"))
        lattice_flag = 0;
      else if ((latarg == "yes") || (latarg == "on") || (latarg == "true") || (latarg == "moving"))
        lattice_flag = 1;
      else error->all(FLERR,"Illegal fix/nve/spin command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"selfconsistent") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix nve/spin selfconsistent", error);
      sciter = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (sciter < 0)
        error->all(FLERR, iarg+1, "Fix nve/spin selfconsistent must be >= 0, not {}", sciter);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix/nve/spin command");
  }

  // check if the atom/spin style is defined

  if (!atom->sp_flag)
    error->all(FLERR,"Fix nve/spin requires atom/spin style");

  // initialize the magnetic interaction flags

  pair_spin_flag = 0;
  long_spin_flag = 0;
  precession_spin_flag = 0;
  maglangevin_flag = 0;
  tdamp_flag = temp_flag = 0;
  setforce_spin_flag = 0;
}

/* ---------------------------------------------------------------------- */

FixNVESpin::~FixNVESpin()
{
  memory->destroy(rsec);
  memory->destroy(stack_head);
  memory->destroy(stack_foot);
  memory->destroy(forward_stacks);
  memory->destroy(backward_stacks);
  delete [] spin_pairs;
  delete [] locklangevinspin;
  delete [] lockprecessionspin;
}

/* ---------------------------------------------------------------------- */

int FixNVESpin::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_NEIGHBOR;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESpin::init()
{

  // set timesteps

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dts = 0.25 * update->dt;
  npairs = npairspin = 0;

  // set ptrs on Pair/Spin styles

  // loop 1: obtain # of Pairs, and # of Pair/Spin styles

  npairspin = 0;
  PairHybrid *hybrid = dynamic_cast<PairHybrid *>(force->pair_match("^hybrid",0));
  if (force->pair_match("^spin",0,0)) {        // only one Pair/Spin style
    pair = force->pair_match("^spin",0,0);
    if (hybrid == nullptr) npairs = 1;
    else npairs = hybrid->nstyles;
    npairspin = 1;
  } else if (force->pair_match("^spin",0,1)) { // more than one Pair/Spin style
    pair = force->pair_match("^spin",0,1);
    if (hybrid == nullptr) npairs = 1;
    else npairs = hybrid->nstyles;
    for (int i = 0; i<npairs; i++) {
      if (force->pair_match("^spin",0,i+1)) {
        npairspin ++;
      }
    }
  }

  // init length of vector of ptrs to Pair/Spin styles

  if (npairspin > 0) {
    delete[] spin_pairs;
    spin_pairs = new PairSpin*[npairspin];
  }

  // loop 2: fill vector with ptrs to Pair/Spin styles

  int count1 = 0;
  if (npairspin == 1) {
    count1 = 1;
    spin_pairs[0] = dynamic_cast<PairSpin *>(force->pair_match("^spin",0,0));
  } else if (npairspin > 1) {
    for (int i = 0; i<npairs; i++) {
      if (force->pair_match("^spin",0,i+1)) {
        spin_pairs[count1] = dynamic_cast<PairSpin *>(force->pair_match("^spin",0,i+1));
        count1++;
      }
    }
  }

  if (count1 != npairspin)
    error->all(FLERR,"Incorrect number of spin pair styles");

  // set pair/spin and long/spin flags

  if (npairspin >= 1) pair_spin_flag = 1;

  for (int i = 0; i<npairs; i++) {
    if (force->pair_match("spin/long",0,i)) {
      long_spin_flag = 1;
    }
  }

  // set ptrs for fix precession/spin styles
  // reset flags and cached pointers first, since the fixes
  // may have been deleted since a previous run

  precession_spin_flag = maglangevin_flag = setforce_spin_flag = 0;
  locksetforcespin = nullptr;

  auto precfixes = modify->get_fix_by_style("^precession/spin");
  nprecspin = (int) precfixes.size();
  if (nprecspin > 0) {
    precession_spin_flag = 1;
    delete[] lockprecessionspin;
    lockprecessionspin = new FixPrecessionSpin*[nprecspin];
    for (int i = 0; i < nprecspin; i++)
      lockprecessionspin[i] = dynamic_cast<FixPrecessionSpin *>(precfixes[i]);
  }

  // set ptrs for fix langevin/spin styles

  auto langfixes = modify->get_fix_by_style("^langevin/spin");
  nlangspin = (int) langfixes.size();
  if (nlangspin > 0) {
    maglangevin_flag = 1;
    delete[] locklangevinspin;
    locklangevinspin = new FixLangevinSpin*[nlangspin];
    for (int i = 0; i < nlangspin; i++) {
      locklangevinspin[i] = dynamic_cast<FixLangevinSpin *>(langfixes[i]);

      // the random field would be drawn again at each iteration

      if ((sciter > 0) && locklangevinspin[i]->temp_flag)
        error->all(FLERR, Error::NOLASTLINE, "Fix nve/spin selfconsistent cannot be used with "
                   "fix langevin/spin {} at a non-zero temperature", locklangevinspin[i]->id);
    }
  }

  // ptrs FixSetForceSpin classes

  auto setforcefixes = modify->get_fix_by_style("^setforce/spin");
  if (!setforcefixes.empty()) {
    setforce_spin_flag = 1;
    locksetforcespin = dynamic_cast<FixSetForceSpin *>(setforcefixes.back());
  }

  // setting the sector variables/lists

  nsectors = 0;
  memory->destroy(rsec);
  memory->create(rsec,3,"nve/spin:rsec");

  // perform the sectoring operation

  if (sector_flag) sectoring();

  // init. size of stacking lists (sectoring)

  nlocal_max = atom->nlocal;
  memory->grow(stack_head,nsectors,"nve/spin:stack_head");
  memory->grow(stack_foot,nsectors,"nve/spin:stack_foot");
  memory->grow(backward_stacks,nlocal_max,"nve/spin:backward_stacks");
  memory->grow(forward_stacks,nlocal_max,"nve/spin:forward_stacks");
}

/* ---------------------------------------------------------------------- */

void FixNVESpin::initial_integrate(int /*vflag*/)
{
  double dtfm;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  int *type = atom->type;
  int *mask = atom->mask;

  // update half v for all atoms

  if (lattice_flag) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (rmass) dtfm = dtf / rmass[i];
        else dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
    }
  }

  // update half s for all atoms

  if (sector_flag) {                            // sectoring seq. update
    for (int j = 0; j < nsectors; j++) {        // advance quarter s for nlocal
      comm->forward_comm();
      int i = stack_foot[j];
      while (i >= 0) {
        const int next = forward_stacks[i];
        if (mask[i] & groupbit) {
          UpdateSingleSpin(i);
        }
        i = next;
      }
    }
    for (int j = nsectors-1; j >= 0; j--) {     // advance quarter s for nlocal
      comm->forward_comm();
      int i = stack_head[j];
      while (i >= 0) {
        const int next = backward_stacks[i];
        if (mask[i] & groupbit) {
          UpdateSingleSpin(i);
        }
        i = next;
      }
    }
  } else {                                       // serial seq. update
    comm->forward_comm();                        // comm. positions of ghost atoms
    for (int i = 0; i < nlocal; i++) {           // advance quarter s for nlocal
      if (mask[i] & groupbit) {
        UpdateSingleSpin(i);
      }
    }
    for (int i = nlocal-1; i >= 0; i--) {        // advance quarter s for nlocal
      if (mask[i] & groupbit) {
        UpdateSingleSpin(i);
      }
    }
  }

  // update x for all particles

  if (lattice_flag) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  }

  // update half s for all particles

  if (sector_flag) {                            // sectoring seq. update
    for (int j = 0; j < nsectors; j++) {        // advance quarter s for nlocal
      comm->forward_comm();
      int i = stack_foot[j];
      while (i >= 0) {
        const int next = forward_stacks[i];
        if (mask[i] & groupbit) {
          UpdateSingleSpin(i);
        }
        i = next;
      }
    }
    for (int j = nsectors-1; j >= 0; j--) {     // advance quarter s for nlocal
      comm->forward_comm();
      int i = stack_head[j];
      while (i >= 0) {
        const int next = backward_stacks[i];
        if (mask[i] & groupbit) {
          UpdateSingleSpin(i);
        }
        i = next;
      }
    }
  } else {                                      // serial seq. update
    comm->forward_comm();                       // comm. positions of ghost atoms
    for (int i = 0; i < nlocal; i++) {          // advance quarter s for nlocal-1
      if (mask[i] & groupbit) {
        UpdateSingleSpin(i);
      }
    }
    for (int i = nlocal-1; i >= 0; i--) {       // advance quarter s for nlocal-1
      if (mask[i] & groupbit) {
        UpdateSingleSpin(i);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   setup pre_neighbor()
---------------------------------------------------------------------- */

void FixNVESpin::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ----------------------------------------------------------------------
   store in two linked lists the advance order of the spins (sectoring)
---------------------------------------------------------------------- */

void FixNVESpin::pre_neighbor()
{
  double **x = atom->x;
  int nlocal = atom->nlocal;

  if (nlocal_max < nlocal) {                    // grow linked lists if necessary
    nlocal_max = nlocal;
    memory->grow(backward_stacks,nlocal_max,"nve/spin:backward_stacks");
    memory->grow(forward_stacks,nlocal_max,"nve/spin:forward_stacks");
  }

  for (int j = 0; j < nsectors; j++) {
    stack_head[j] = -1;
    stack_foot[j] = -1;
  }

  int nseci;
  for (int j = 0; j < nsectors; j++) {          // stacking backward order
    for (int i = 0; i < nlocal; i++) {
      nseci = coords2sector(x[i]);
      if (j != nseci) continue;
      backward_stacks[i] = stack_head[j];
      stack_head[j] = i;
    }
  }
  for (int j = nsectors-1; j >= 0; j--) {       // stacking forward order
    for (int i = nlocal-1; i >= 0; i--) {
      nseci = coords2sector(x[i]);
      if (j != nseci) continue;
      forward_stacks[i] = stack_foot[j];
      stack_foot[j] = i;
    }
  }

}

/* ----------------------------------------------------------------------
   compute the magnetic torque for the spin ii
---------------------------------------------------------------------- */

void FixNVESpin::ComputeInteractionsSpin(int i)
{
  double spi[3], fmi[3];

  double **sp = atom->sp;
  double **fm = atom->fm;

  // force computation for spin i

  spi[0] = sp[i][0];
  spi[1] = sp[i][1];
  spi[2] = sp[i][2];

  fmi[0] = fmi[1] = fmi[2] = 0.0;

  // update magnetic pair interactions

  if (pair_spin_flag) {
    for (int k = 0; k < npairspin; k++) {
      spin_pairs[k]->compute_single_pair(i,fmi);
    }
  }

  // update magnetic precession interactions

  if (precession_spin_flag) {
    for (int k = 0; k < nprecspin; k++) {
      lockprecessionspin[k]->compute_single_precession(i,spi,fmi);
    }
  }

  // update langevin damping and random force

  if (maglangevin_flag) {               // mag. langevin
    for (int k = 0; k < nlangspin; k++) {
      locklangevinspin[k]->compute_single_langevin(i,spi,fmi);
    }
  }

  // update setforce of magnetic interactions

  if (setforce_spin_flag) {
    locksetforcespin->single_setforce_spin(i,fmi);
  }

  // replace the magnetic force fm[i] by its new value fmi

  fm[i][0] = fmi[0];
  fm[i][1] = fmi[1];
  fm[i][2] = fmi[2];
}

/* ----------------------------------------------------------------------
   divide each domain into 8 sectors
   for triclinic boxes, the sectors are defined in lamda coordinates
---------------------------------------------------------------------- */

void FixNVESpin::sectoring()
{
  // extract larger cutoff from PairSpin styles

  double rv, cutoff;
  rv = cutoff = 0.0;
  int dim = 0;
  for (int i = 0; i < npairspin ; i++) {
    cutoff = *((double *) spin_pairs[i]->extract("cut",dim));
    rv = MAX(rv,cutoff);
  }

  if (rv == 0.0)
   error->all(FLERR, Error::NOLASTLINE,
              "No suitable cutoff found for sectoring operation: rv = {}", rv);

  // extent of the subdomain and cutoff in each direction
  // for triclinic boxes, convert the cutoff to the lamda distance
  //   between parallel faces of the box, same as for the ghost cutoff in Comm

  double rs[3], rcut[3];
  if (domain->triclinic) {
    const double *h_inv = domain->h_inv;
    for (int d = 0; d < 3; d++) rs[d] = domain->subhi_lamda[d] - domain->sublo_lamda[d];
    rcut[0] = rv * sqrt(h_inv[0] * h_inv[0] + h_inv[5] * h_inv[5] + h_inv[4] * h_inv[4]);
    rcut[1] = rv * sqrt(h_inv[1] * h_inv[1] + h_inv[3] * h_inv[3]);
    rcut[2] = rv * h_inv[2];
  } else {
    for (int d = 0; d < 3; d++) {
      rs[d] = domain->subhi[d] - domain->sublo[d];
      rcut[d] = rv;
    }
  }

  // split a direction in two sectors only if each is at least as wide as the cutoff

  int sec[3];
  for (int d = 0; d < 3; d++) sec[d] = (rs[d] >= 2.0 * rcut[d]) ? 2 : 1;

  nsectors = sec[0]*sec[1]*sec[2];

  if (sector_flag && (nsectors != 8))
    error->all(FLERR, Error::NOLASTLINE,
               "Illegal sectoring operation resulting in {} sectors instead of 8", nsectors);

  for (int d = 0; d < 3; d++) rsec[d] = rs[d] / sec[d];
}

/* ----------------------------------------------------------------------
   define sector for an atom at a position x[i]
---------------------------------------------------------------------- */

int FixNVESpin::coords2sector(double *x)
{
  double lamda[3];
  double *xs = x;
  double *sublo = domain->sublo;

  if (domain->triclinic) {
    domain->x2lamda(x, lamda);
    xs = lamda;
    sublo = domain->sublo_lamda;
  }

  int seci[3];
  seci[0] = xs[0] > (sublo[0] + rsec[0]);
  seci[1] = xs[1] > (sublo[1] + rsec[1]);
  seci[2] = xs[2] > (sublo[2] + rsec[2]);

  return (seci[0] + 2*seci[1] + 4*seci[2]);
}

/* ----------------------------------------------------------------------
   rotate the spin s about the field fm during a time dts:
   g = s + dts/2 fm x (s + g) (Cayley transform), so that |g| = |s|
   single-spin propagator of Omelyan, Mryglod, and Folk, Phys Rev Lett 86,
   898 (2001), in the form of Eq. (15) of Tranchida et al. (2018)
---------------------------------------------------------------------- */

static void RotateSpin(const double *s, const double *fm, double dts, double *g)
{
  double fm2,energy,dts2;
  double cp[3];

  fm2 = (fm[0]*fm[0])+(fm[1]*fm[1])+(fm[2]*fm[2]);
  energy = (s[0]*fm[0])+(s[1]*fm[1])+(s[2]*fm[2]);
  dts2 = dts*dts;

  cp[0] = fm[1]*s[2] - fm[2]*s[1];
  cp[1] = fm[2]*s[0] - fm[0]*s[2];
  cp[2] = fm[0]*s[1] - fm[1]*s[0];

  g[0] = s[0] + cp[0]*dts;
  g[1] = s[1] + cp[1]*dts;
  g[2] = s[2] + cp[2]*dts;

  g[0] += (fm[0]*energy - 0.5*s[0]*fm2)*0.5*dts2;
  g[1] += (fm[1]*energy - 0.5*s[1]*fm2)*0.5*dts2;
  g[2] += (fm[2]*energy - 0.5*s[2]*fm2)*0.5*dts2;

  g[0] /= (1.0 + 0.25*fm2*dts2);
  g[1] /= (1.0 + 0.25*fm2*dts2);
  g[2] /= (1.0 + 0.25*fm2*dts2);
}

/* ----------------------------------------------------------------------
   compute the field on spin i and advance it of a timestep dts
   with sciter > 0, the field is evaluated at the midpoint (s + s')/2 of
   the update, found by sciter fixed-point iterations; this makes the
   update time-symmetric also when the field on spin i depends on s_i
   (anisotropy, biquadratic exchange, ...)
---------------------------------------------------------------------- */

void FixNVESpin::UpdateSingleSpin(int i)
{
  double **sp = atom->sp;
  double **fm = atom->fm;
  double s0[3],g[3];

  ComputeInteractionsSpin(i);

  if (sciter > 0) {
    s0[0] = sp[i][0];
    s0[1] = sp[i][1];
    s0[2] = sp[i][2];
    for (int k = 0; k < sciter; k++) {
      RotateSpin(s0,fm[i],dts,g);
      sp[i][0] = 0.5*(s0[0] + g[0]);
      sp[i][1] = 0.5*(s0[1] + g[1]);
      sp[i][2] = 0.5*(s0[2] + g[2]);
      ComputeInteractionsSpin(i);
    }
    sp[i][0] = s0[0];
    sp[i][1] = s0[1];
    sp[i][2] = s0[2];
  }

  AdvanceSingleSpin(i);
}

/* ----------------------------------------------------------------------
   advance the spin i of a timestep dts
---------------------------------------------------------------------- */

void FixNVESpin::AdvanceSingleSpin(int i)
{
  int j=0;
  int *sametag = atom->sametag;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double g[3];

  RotateSpin(sp[i],fm[i],dts,g);

  sp[i][0] = g[0];
  sp[i][1] = g[1];
  sp[i][2] = g[2];

  // renormalization (check if necessary)

  // msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
  // scale = 1.0/sqrt(msq);
  // sp[i][0] *= scale;
  // sp[i][1] *= scale;
  // sp[i][2] *= scale;

  // comm. sp[i] to atoms with same tag (for serial algo)

  if (!sector_flag) {
    if (sametag[i] >= 0) {
      j = sametag[i];
      while (j >= 0) {
        sp[j][0] = sp[i][0];
        sp[j][1] = sp[i][1];
        sp[j][2] = sp[i][2];
        j = sametag[j];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixNVESpin::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  int *type = atom->type;
  int *mask = atom->mask;

  // update half v for all particles

  if (lattice_flag) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (rmass) dtfm = dtf / rmass[i];
        else dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

double FixNVESpin::memory_usage()
{
  double bytes = (double) nlocal_max * 2 * sizeof(int);     // backward_stacks + forward_stacks
  bytes += (double) nsectors * 2 * sizeof(int);             // stack_head + stack_foot
  return bytes;
}
