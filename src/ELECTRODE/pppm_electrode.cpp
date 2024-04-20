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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meißner (TUHH)
------------------------------------------------------------------------- */

#include "pppm_electrode.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "boundary_correction.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "force.h"
#include "grid3d.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neighbor.h"
#include "pair.h"
#include "remap_wrap.h"
#include "slab_dipole.h"
#include "update.h"
#include "wire_dipole.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

static constexpr int MAXORDER = 7;
static constexpr int OFFSET = 16384;
static constexpr double EPS_HOC = 1.0e-7;

enum { REVERSE_RHO };
enum { FORWARD_IK, FORWARD_AD, FORWARD_IK_PERATOM, FORWARD_AD_PERATOM };

static constexpr FFT_SCALAR ZEROF = 0.0;

static const char cite_pppm_electrode[] =
    "kspace_style pppm/electrode command:\n\n"
    "@article{Ahrens2021,\n"
    "author = {Ahrens-Iwers, Ludwig J.V. and Mei{\\ss}ner, Robert H.},\n"
    "doi = {10.1063/5.0063381},\n"
    "title = {{Constant potential simulations on a mesh}},\n"
    "journal = {Journal of Chemical Physics},\n"
    "year = {2021}\n"
    "volume = {155},\n"
    "pages = {104104},\n"
    "}\n";

/* ---------------------------------------------------------------------- */

PPPMElectrode::PPPMElectrode(LAMMPS *lmp) :
    PPPM(lmp), electrolyte_density_brick(nullptr), electrolyte_density_fft(nullptr),
    boundcorr(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_pppm_electrode);

  group_group_enable = 0;
  electrolyte_density_brick = nullptr;
  electrolyte_density_fft = nullptr;
  compute_vector_called = false;
  last_source_grpbit = 1 << 0;    // initialize to "all"
  last_invert_source = false;     // not sure what to initialize here
}
/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

PPPMElectrode::~PPPMElectrode()
{
  if (copymode) return;

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();
  if (group_allocate_flag) deallocate_groups();
  memory->destroy(part2grid);
  memory->destroy(acons);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMElectrode::init()
{
  if (me == 0) utils::logmesg(lmp, "PPPM/electrode initialization ...\n");

  // error check
  if (slabflag == 3)
    error->all(FLERR, "Cannot (yet) use PPPM/electrode with 'kspace_modify slab ew2d'");

  triclinic_check();
  triclinic = domain->triclinic;
  if (triclinic) error->all(FLERR, "Cannot (yet) use PPPM/electrode with triclinic box ");
  if (domain->dimension == 2) error->all(FLERR, "Cannot use PPPM/electrode with 2d simulation");

  if (!atom->q_flag) error->all(FLERR, "KSpace style requires atom attribute q");

  if (slabflag == 0 && wireflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR, "Cannot use non-periodic boundaries with PPPM/electrode");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || domain->boundary[2][0] != 1 ||
        domain->boundary[2][1] != 1)
      error->all(FLERR, "Incorrect boundaries with slab PPPM/electrode");
  } else if (wireflag) {
    if (domain->zperiodic != 1 || domain->boundary[0][0] != 1 || domain->boundary[0][1] != 1 ||
        domain->boundary[1][0] != 1 || domain->boundary[1][1] != 1)
      error->all(FLERR, "Incorrect boundaries with wire PPPM/electrode");
  }

  if (order < 2 || order > MAXORDER)
    error->all(FLERR, "PPPM/electrode order cannot be < 2 or > {}", MAXORDER);

  // compute two charge force

  two_charge();

  // extract short-range Coulombic cutoff from pair style

  pair_check();

  int itmp = 0;
  double *p_cutoff = (double *) force->pair->extract("cut_coul", itmp);
  if (p_cutoff == nullptr) error->all(FLERR, "KSpace style is incompatible with Pair style");
  cutoff = *p_cutoff;

  // if kspace is TIP4P, extract TIP4P params from pair style
  // bond/angle are not yet init(), so ensure equilibrium request is valid

  qdist = 0.0;

  if (tip4pflag) {
    if (me == 0) utils::logmesg(lmp, "  extracting TIP4P info from pair style\n");

    double *p_qdist = (double *) force->pair->extract("qdist", itmp);
    int *p_typeO = (int *) force->pair->extract("typeO", itmp);
    int *p_typeH = (int *) force->pair->extract("typeH", itmp);
    int *p_typeA = (int *) force->pair->extract("typeA", itmp);
    int *p_typeB = (int *) force->pair->extract("typeB", itmp);
    if (!p_qdist || !p_typeO || !p_typeH || !p_typeA || !p_typeB)
      error->all(FLERR, "Pair style is incompatible with TIP4P KSpace style");
    qdist = *p_qdist;
    typeO = *p_typeO;
    typeH = *p_typeH;
    int typeA = *p_typeA;
    int typeB = *p_typeB;

    if (force->angle == nullptr || force->bond == nullptr || force->angle->setflag == nullptr ||
        force->bond->setflag == nullptr)
      error->all(FLERR, "Bond and angle potentials must be defined for TIP4P");
    if (typeA < 1 || typeA > atom->nangletypes || force->angle->setflag[typeA] == 0)
      error->all(FLERR, "Bad TIP4P angle type for PPPM/TIP4P");
    if (typeB < 1 || typeB > atom->nbondtypes || force->bond->setflag[typeB] == 0)
      error->all(FLERR, "Bad TIP4P bond type for PPPM/TIP4P");
    double theta = force->angle->equilibrium_angle(typeA);
    double blen = force->bond->equilibrium_distance(typeB);
    alpha = qdist / (cos(0.5 * theta) * blen);
  }

  // compute qsum & qsqsum and warn if not charge-neutral

  scale = 1.0;
  qqrd2e = force->qqrd2e;
  qsum_qsq();
  natoms_original = atom->natoms;

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0)
    accuracy = accuracy_absolute;
  else
    accuracy = accuracy_relative * two_charge_force;

  // free all arrays previously allocated

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();
  if (group_allocate_flag) deallocate_groups();

  // setup FFT grid resolution and g_ewald
  // normally one iteration thru while loop is all that is required
  // if grid stencil does not extend beyond neighbor proc
  //   or overlap is allowed, then done
  // else reduce order and try again

  gc = nullptr;
  int iteration = 0;

  while (order >= minorder) {
    if (iteration && me == 0)
      error->warning(FLERR,
                     "Reducing PPPM/electrode order b/c stencil extends "
                     "beyond nearest neighbor processor");

    if (stagger_flag && !differentiation_flag) compute_gf_denom();
    set_grid_global();
    set_grid_local();
    if (overlap_allowed) break;

    gc = new Grid3d(lmp, world, nx_pppm, ny_pppm, nz_pppm, nxlo_in, nxhi_in, nylo_in, nyhi_in,
                    nzlo_in, nzhi_in, nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);

    int tmp1, tmp2;
    gc->setup_comm(tmp1, tmp2);
    if (gc->ghost_adjacent()) break;
    delete gc;

    order--;
    iteration++;
  }

  if (order < minorder) error->all(FLERR, "PPPM/electrode order < minimum allowed order");
  if (!overlap_allowed && !gc->ghost_adjacent())
    error->all(FLERR, "PPPM/electrode grid stencil extends beyond nearest neighbor processor");
  if (gc) delete gc;

  // adjust g_ewald

  if (!gewaldflag) adjust_gewald();

  // calculate the final accuracy

  double estimated_accuracy = final_accuracy();

  // print stats

  int ngrid_max, nfft_both_max;
  MPI_Allreduce(&ngrid, &ngrid_max, 1, MPI_INT, MPI_MAX, world);
  MPI_Allreduce(&nfft_both, &nfft_both_max, 1, MPI_INT, MPI_MAX, world);

  if (me == 0) {
    std::string mesg = fmt::format("  G vector (1/distance) = {:.8g}\n", g_ewald);
    mesg += fmt::format("  grid = {} {} {}\n", nx_pppm, ny_pppm, nz_pppm);
    mesg += fmt::format("  stencil order = {}\n", order);
    mesg += fmt::format("  estimated absolute RMS force accuracy = {:.8g}\n", estimated_accuracy);
    mesg += fmt::format("  estimated relative force accuracy = {:.8g}\n",
                        estimated_accuracy / two_charge_force);
    mesg += "  using " LMP_FFT_PREC " precision " LMP_FFT_LIB "\n";
    mesg += fmt::format("  3d grid and FFT values/proc = {} {}\n", ngrid_max, nfft_both_max);
    utils::logmesg(lmp, mesg);
  }

  // allocate K-space dependent memory
  // don't invoke allocate peratom() or group(), will be allocated when needed

  allocate();

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  if (differentiation_flag == 1) compute_sf_precoeff();
  compute_rho_coeff();
  compute_step = -1;
}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void PPPMElectrode::setup()
{
  // perform some checks to avoid illegal boundaries with read_data

  if (slabflag == 0 && wireflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR, "Cannot use non-periodic boundaries with PPPM/electrode");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || domain->boundary[2][0] != 1 ||
        domain->boundary[2][1] != 1)
      error->all(FLERR, "Incorrect boundaries with slab PPPM/electrode");
  } else if (wireflag) {
    if (domain->zperiodic != 1 || domain->boundary[0][0] != 1 || domain->boundary[0][1] != 1 ||
        domain->boundary[1][0] != 1 || domain->boundary[1][1] != 1)
      error->all(FLERR, "Incorrect boundaries with wire PPPM/electrode");
  }

  int i, j, k, n;
  double *prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  prd = domain->prd;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;
  volume = xprd_wire * yprd_wire * zprd_slab;

  delxinv = nx_pppm / xprd_wire;
  delyinv = ny_pppm / yprd_wire;
  delzinv = nz_pppm / zprd_slab;

  delvolinv = delxinv * delyinv * delzinv;

  double unitkx = (MY_2PI / xprd_wire);
  double unitky = (MY_2PI / yprd_wire);
  double unitkz = (MY_2PI / zprd_slab);

  // fkx,fky,fkz for my FFT grid pts

  for (i = nxlo_fft; i <= nxhi_fft; i++) {
    int per = i - nx_pppm * (2 * i / nx_pppm);
    fkx[i] = unitkx * per;
  }

  for (i = nylo_fft; i <= nyhi_fft; i++) {
    int per = i - ny_pppm * (2 * i / ny_pppm);
    fky[i] = unitky * per;
  }

  for (i = nzlo_fft; i <= nzhi_fft; i++) {
    int per = i - nz_pppm * (2 * i / nz_pppm);
    fkz[i] = unitkz * per;
  }

  // virial coefficients

  double sqk, vterm;

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++) {
    for (j = nylo_fft; j <= nyhi_fft; j++) {
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        sqk = fkx[i] * fkx[i] + fky[j] * fky[j] + fkz[k] * fkz[k];
        if (sqk == 0.0) {
          vg[n][0] = 0.0;
          vg[n][1] = 0.0;
          vg[n][2] = 0.0;
          vg[n][3] = 0.0;
          vg[n][4] = 0.0;
          vg[n][5] = 0.0;
        } else {
          vterm = -2.0 * (1.0 / sqk + 0.25 / (g_ewald * g_ewald));
          vg[n][0] = 1.0 + vterm * fkx[i] * fkx[i];
          vg[n][1] = 1.0 + vterm * fky[j] * fky[j];
          vg[n][2] = 1.0 + vterm * fkz[k] * fkz[k];
          vg[n][3] = vterm * fkx[i] * fky[j];
          vg[n][4] = vterm * fkx[i] * fkz[k];
          vg[n][5] = vterm * fky[j] * fkz[k];
        }
        n++;
      }
    }
  }

  if (differentiation_flag == 1)
    compute_gf_ad();
  else
    compute_gf_ik();
}

/* ----------------------------------------------------------------------
   reset local grid arrays and communication stencils
   called by fix balance b/c it changed sizes of processor sub-domains
------------------------------------------------------------------------- */

void PPPMElectrode::reset_grid()
{
  // free all arrays previously allocated

  deallocate();
  if (peratom_allocate_flag) deallocate_peratom();
  if (group_allocate_flag) deallocate_groups();

  // reset portion of global grid that each proc owns

  set_grid_local();

  // reallocate K-space dependent memory
  // check if grid communication is now overlapping if not allowed
  // don't invoke allocate peratom() or group(), will be allocated when needed

  allocate();

  if (!overlap_allowed && !gc->ghost_adjacent())
    error->all(FLERR,
               "PPPM/electrode grid stencil extends "
               "beyond nearest neighbor processor");

  // pre-compute Green's function denomiator expansion
  // pre-compute 1d charge distribution coefficients

  compute_gf_denom();
  if (differentiation_flag == 1) compute_sf_precoeff();
  compute_rho_coeff();

  // pre-compute volume-dependent coeffs for portion of grid I now own

  setup();
}

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMElectrode::compute(int eflag, int vflag)
{
  int i, j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag, vflag);

  if (evflag_atom && !peratom_allocate_flag) allocate_peratom();

  // if atom count has changed, update qsum and qsqsum

  qsum_qsq();
  natoms_original = atom->natoms;

  // return if there are no charges

  start_compute();

  /*
  if (compute_vector_called && last_invert_source) {
    // electrolyte_density_brick is filled, so we can grab only electrode atoms.
    // Does not work for direct cg algorithm because electrode charges change after compute_vector.
    // Therefore, only when last_invert_source true.
    // TODO: this is dangerous now that compute_vector's interface has been
    // changed since a compute could call an arbitrary source, needs tightening
    make_rho_in_brick(last_source_grpbit, density_brick, !last_invert_source);
    gc->reverse_comm(Grid3d::KSPACE, this, REVERSE_RHO, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                     MPI_FFT_SCALAR);
    for (int nz = nzlo_out; nz <= nzhi_out; nz++)
      for (int ny = nylo_out; ny <= nyhi_out; ny++)
        for (int nx = nxlo_out; nx <= nxhi_out; nx++) {
          density_brick[nz][ny][nx] += electrolyte_density_brick[nz][ny][nx];
        }
  } else {
  */
  particle_map();
  make_rho();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  gc->reverse_comm(Grid3d::KSPACE, this, REVERSE_RHO, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                   MPI_FFT_SCALAR);
  //}

  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  poisson();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1)
    gc->forward_comm(Grid3d::KSPACE, this, FORWARD_AD, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                     MPI_FFT_SCALAR);
  else
    gc->forward_comm(Grid3d::KSPACE, this, FORWARD_IK, 3, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                     MPI_FFT_SCALAR);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      gc->forward_comm(Grid3d::KSPACE, this, FORWARD_AD_PERATOM, 6, sizeof(FFT_SCALAR), gc_buf1,
                       gc_buf2, MPI_FFT_SCALAR);
    else if (differentiation_flag == 0)
      gc->forward_comm(Grid3d::KSPACE, this, FORWARD_IK_PERATOM, 7, sizeof(FFT_SCALAR), gc_buf1,
                       gc_buf2, MPI_FFT_SCALAR);
  }

  // calculate the force on my particles

  fieldforce();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in volume-dependent term

  const double qscale = qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy, &energy_all, 1, MPI_DOUBLE, MPI_SUM, world);
    energy = energy_all;

    energy *= 0.5 * volume;
    energy -= g_ewald * qsqsum / MY_PIS + MY_PI2 * qsum * qsum / (g_ewald * g_ewald * volume);
    energy *= qscale;
  }

  // sum global virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial, virial_all, 6, MPI_DOUBLE, MPI_SUM, world);
    for (i = 0; i < 6; i++) virial[i] = 0.5 * qscale * volume * virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction
  // ntotal accounts for TIP4P tallying eatom/vatom for ghost atoms

  if (evflag_atom) {
    double *q = atom->q;
    int nlocal = atom->nlocal;
    int ntotal = nlocal;
    if (tip4pflag) ntotal += atom->nghost;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        eatom[i] -=
            g_ewald * q[i] * q[i] / MY_PIS + MY_PI2 * q[i] * qsum / (g_ewald * g_ewald * volume);
        eatom[i] *= qscale;
      }
      for (i = nlocal; i < ntotal; i++) eatom[i] *= 0.5 * qscale;
    }

    if (vflag_atom) {
      for (i = 0; i < ntotal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5 * qscale;
    }
  }

  boundcorr->compute_corr(qsum, eflag_atom, eflag_global, energy, eatom);
  compute_vector_called = false;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */
void PPPMElectrode::start_compute()
{
  if (compute_step < update->ntimestep) {
    if (compute_step == -1) setup();
    boxlo = domain->boxlo;
    // extend size of per-atom arrays if necessary
    if (atom->nmax > nmax) {
      memory->destroy(part2grid);
      nmax = atom->nmax;
      memory->create(part2grid, nmax, 3, "pppm/electrode:part2grid");
    }
    // find grid points for all my particles
    // map my particle charge onto my local 3d density grid
    particle_map();
    compute_step = update->ntimestep;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */
void PPPMElectrode::compute_vector(double *vec, int sensor_grpbit, int source_grpbit,
                                   bool invert_source)
{
  start_compute();

  // temporarily store and switch pointers so we can use brick2fft() for
  // electrolyte density (without writing an additional function)
  FFT_SCALAR ***density_brick_real = density_brick;
  FFT_SCALAR *density_fft_real = density_fft;
  particle_map();
  make_rho_in_brick(source_grpbit, electrolyte_density_brick, invert_source);
  density_brick = electrolyte_density_brick;
  density_fft = electrolyte_density_fft;
  gc->reverse_comm(Grid3d::KSPACE, this, REVERSE_RHO, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                   MPI_FFT_SCALAR);
  brick2fft();
  // switch back pointers
  density_brick = density_brick_real;
  density_fft = density_fft_real;

  // transform electrolyte charge density (r -> k) (complex conjugate)
  for (int i = 0, n = 0; i < nfft; i++) {
    work1[n++] = electrolyte_density_fft[i];
    work1[n++] = ZEROF;
  }
  fft1->compute(work1, work1, -1);

  // k->r FFT of Green's * electrolyte density = u_brick
  for (int i = 0, n = 0; i < nfft; i++) {
    work2[n] = work1[n] * greensfn[i];
    n++;
    work2[n] = work1[n] * greensfn[i];
    n++;
  }
  fft2->compute(work2, work2, 1);
  for (int k = nzlo_in, n = 0; k <= nzhi_in; k++)
    for (int j = nylo_in; j <= nyhi_in; j++)
      for (int i = nxlo_in; i <= nxhi_in; i++) {
        u_brick[k][j][i] = work2[n];
        n += 2;
      }
  gc->forward_comm(Grid3d::KSPACE, this, FORWARD_AD, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                   MPI_FFT_SCALAR);
  project_psi(vec, sensor_grpbit);
  compute_vector_called = true;
}

void PPPMElectrode::project_psi(double *vec, int sensor_grpbit)
{
  // project u_brick with weight matrix
  double **x = atom->x;
  int *mask = atom->mask;
  const bigint ngridtotal = (bigint) nx_pppm * ny_pppm * nz_pppm;
  const double scaleinv = 1.0 / ngridtotal;

  for (int i = 0; i < atom->nlocal; i++) {
    if (!(mask[i] & sensor_grpbit)) continue;
    double v = 0.;
    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // (dx,dy,dz) = distance to "lower left" grid pt
    // (mx,my,mz) = global coords of moving stencil pt
    int nix = part2grid[i][0];
    int niy = part2grid[i][1];
    int niz = part2grid[i][2];
    FFT_SCALAR dix = nix + shiftone - (x[i][0] - boxlo[0]) * delxinv;
    FFT_SCALAR diy = niy + shiftone - (x[i][1] - boxlo[1]) * delyinv;
    FFT_SCALAR diz = niz + shiftone - (x[i][2] - boxlo[2]) * delzinv;
    compute_rho1d(dix, diy, diz);
    for (int ni = nlower; ni <= nupper; ni++) {
      double iz0 = rho1d[2][ni];
      int miz = ni + niz;
      for (int mi = nlower; mi <= nupper; mi++) {
        double iy0 = iz0 * rho1d[1][mi];
        int miy = mi + niy;
        for (int li = nlower; li <= nupper; li++) {
          int mix = li + nix;
          double ix0 = iy0 * rho1d[0][li];
          v += ix0 * u_brick[miz][miy][mix];
        }
      }
    }
    vec[i] += v * scaleinv;
  }
}
/* ----------------------------------------------------------------------
-------------------------------------------------------------------------
*/

void PPPMElectrode::compute_matrix(bigint *imat, double **matrix, bool timer_flag)
{
  compute(1, 0);    // make sure density bricks etc. are set up

  // fft green's function k -> r (double)
  double *greens_real;
  memory->create(greens_real, nz_pppm * ny_pppm * nx_pppm, "pppm/electrode:greens_real");
  memset(greens_real, 0,
         (std::size_t) nz_pppm * (std::size_t) ny_pppm * (std::size_t) nx_pppm * sizeof(double));
  for (int i = 0, n = 0; i < nfft; i++) {
    work2[n++] = greensfn[i];
    work2[n++] = ZEROF;
  }
  fft2->compute(work2, work2, -1);
  for (int k = nzlo_in, n = 0; k <= nzhi_in; k++)
    for (int j = nylo_in; j <= nyhi_in; j++)
      for (int i = nxlo_in; i <= nxhi_in; i++) {
        greens_real[ny_pppm * nx_pppm * k + nx_pppm * j + i] = work2[n];
        n += 2;
      }
  MPI_Allreduce(MPI_IN_PLACE, greens_real, nz_pppm * ny_pppm * nx_pppm, MPI_DOUBLE, MPI_SUM, world);
  int const nlocal = atom->nlocal;
  int nmat = std::count_if(&imat[0], &imat[nlocal], [](int x) {
    return x >= 0;
  });
  MPI_Allreduce(MPI_IN_PLACE, &nmat, 1, MPI_INT, MPI_SUM, world);

  // gather x_ele
  double **x_ele;
  memory->create(x_ele, nmat, 3, "pppm/electrode:x_ele");
  memset(&(x_ele[0][0]), 0, nmat * 3 * sizeof(double));
  double **x = atom->x;
  for (int i = 0; i < nlocal; i++) {
    int ipos = imat[i];
    if (ipos < 0) continue;
    for (int dim = 0; dim < 3; dim++) x_ele[ipos][dim] = x[i][dim];
  }
  MPI_Allreduce(MPI_IN_PLACE, &(x_ele[0][0]), nmat * 3, MPI_DOUBLE, MPI_SUM, world);

  if (conp_one_step)
    one_step_multiplication(imat, greens_real, x_ele, matrix, nmat, timer_flag);
  else
    two_step_multiplication(imat, greens_real, x_ele, matrix, nmat, timer_flag);
  memory->destroy(greens_real);
  memory->destroy(x_ele);
}

/* ----------------------------------------------------------------------*/

void PPPMElectrode::one_step_multiplication(bigint *imat, double *greens_real, double **x_ele,
                                            double **matrix, int const nmat, bool timer_flag)
{
  // map green's function in real space from mesh to particle positions
  // with matrix multiplication 'W^T G W' in one steps. Uses less memory than
  // two_step_multiplication
  //
  int const nlocal = atom->nlocal;
  double **x = atom->x;
  MPI_Barrier(world);
  double step1_time = MPI_Wtime();

  // precalculate rho_1d for local electrode
  std::vector<int> j_list;
  for (int j = 0; j < nlocal; j++) {
    int jpos = imat[j];
    if (jpos < 0) continue;
    j_list.push_back(j);
  }
  int const nj_local = j_list.size();

  FFT_SCALAR ***rho1d_j;
  memory->create(rho1d_j, nj_local, 3, order, "pppm/electrode:rho1d_j");

  for (int jlist_pos = 0; jlist_pos < nj_local; jlist_pos++) {
    int j = j_list[jlist_pos];
    int njx = part2grid[j][0];
    int njy = part2grid[j][1];
    int njz = part2grid[j][2];
    FFT_SCALAR const djx = njx + shiftone - (x[j][0] - boxlo[0]) * delxinv;
    FFT_SCALAR const djy = njy + shiftone - (x[j][1] - boxlo[1]) * delyinv;
    FFT_SCALAR const djz = njz + shiftone - (x[j][2] - boxlo[2]) * delzinv;
    compute_rho1d(djx, djy, djz);
    for (int dim = 0; dim < 3; dim++) {
      for (int oi = 0; oi < order; oi++) { rho1d_j[jlist_pos][dim][oi] = rho1d[dim][oi + nlower]; }
    }
  }

  // nested loops over weights of electrode atoms i and j
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  int const order2 = order * order;
  int const order6 = order2 * order2 * order2;
  double *amesh;
  memory->create(amesh, order6, "pppm/electrode:amesh");
  for (int ipos = 0; ipos < nmat; ipos++) {
    double *_noalias xi_ele = x_ele[ipos];
    // new calculation for nx, ny, nz because part2grid available for nlocal,
    // only
    int nix = static_cast<int>((xi_ele[0] - boxlo[0]) * delxinv + shift) - OFFSET;
    int niy = static_cast<int>((xi_ele[1] - boxlo[1]) * delyinv + shift) - OFFSET;
    int niz = static_cast<int>((xi_ele[2] - boxlo[2]) * delzinv + shift) - OFFSET;
    FFT_SCALAR const dix = nix + shiftone - (xi_ele[0] - boxlo[0]) * delxinv;
    FFT_SCALAR const diy = niy + shiftone - (xi_ele[1] - boxlo[1]) * delyinv;
    FFT_SCALAR const diz = niz + shiftone - (xi_ele[2] - boxlo[2]) * delzinv;
    compute_rho1d(dix, diy, diz);
    int njx = -1;
    int njy = -1;
    int njz = -1;    // force initial build_amesh
    for (int jlist_pos = 0; jlist_pos < nj_local; jlist_pos++) {
      int j = j_list[jlist_pos];
      int ind_amesh = 0;
      int jpos = imat[j];
      if ((ipos < jpos) == !((ipos - jpos) % 2)) continue;
      double aij = 0.;
      if (njx != part2grid[j][0] || njy != part2grid[j][1] || njz != part2grid[j][2]) {
        njx = part2grid[j][0];
        njy = part2grid[j][1];
        njz = part2grid[j][2];
        build_amesh(njx - nix, njy - niy, njz - niz, amesh, greens_real);
      }
      for (int ni = nlower; ni <= nupper; ni++) {    // i's rho1d[dim] indexed from nlower to nupper
        FFT_SCALAR const iz0 = rho1d[2][ni];
        for (int nj = 0; nj < order; nj++) {    // j's rho1d_j[][dim] indexed from 0 to order-1
          FFT_SCALAR const jz0 = rho1d_j[jlist_pos][2][nj];
          for (int mi = nlower; mi <= nupper; mi++) {
            FFT_SCALAR const iy0 = iz0 * rho1d[1][mi];
            for (int mj = 0; mj < order; mj++) {
              FFT_SCALAR const jy0 = jz0 * rho1d_j[jlist_pos][1][mj];
              for (int li = nlower; li <= nupper; li++) {
                FFT_SCALAR const ix0 = iy0 * rho1d[0][li];
                double aij_xscan = 0.;
                for (int lj = 0; lj < order; lj++) {
                  aij_xscan += (double) amesh[ind_amesh] * rho1d_j[jlist_pos][0][lj];
                  ind_amesh++;
                }
                aij += (double) ix0 * jy0 * aij_xscan;
              }
            }
          }
        }
      }
      matrix[ipos][jpos] += aij / volume;
      if (ipos != jpos) matrix[jpos][ipos] += aij / volume;
    }
  }
  memory->destroy(amesh);
  memory->destroy(rho1d_j);
  MPI_Barrier(world);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(lmp, "Single step time: {:.4g} s\n", MPI_Wtime() - step1_time);
}

/* ----------------------------------------------------------------------*/

void PPPMElectrode::build_amesh(const int dx,    // = njx - nix
                                const int dy,    // = njy - niy
                                const int dz,    // = njz - niz
                                double *amesh, double *const greens_real)
{
  auto fmod = [](int x, int n) {    // fast unsigned mod
    int r = abs(x);
    while (r >= n) r -= n;
    return r;
  };
  int ind_amesh = 0;

  for (int iz = 0; iz < order; iz++)
    for (int jz = 0; jz < order; jz++) {
      int const mz = fmod(dz + jz - iz, nz_pppm) * nx_pppm * ny_pppm;
      for (int iy = 0; iy < order; iy++)
        for (int jy = 0; jy < order; jy++) {
          int const my = fmod(dy + jy - iy, ny_pppm) * nx_pppm;
          for (int ix = 0; ix < order; ix++)
            for (int jx = 0; jx < order; jx++) {
              int const mx = fmod(dx + jx - ix, nx_pppm);
              amesh[ind_amesh] = greens_real[mz + my + mx];
              ind_amesh++;
            }
        }
    }
}

/* ----------------------------------------------------------------------*/

void PPPMElectrode::two_step_multiplication(bigint *imat, double *greens_real, double **x_ele,
                                            double **matrix, int const nmat, bool timer_flag)
{
  // map green's function in real space from mesh to particle positions
  // with matrix multiplication 'W^T G W' in two steps. gw is result of
  // first multiplication.
  int const nlocal = atom->nlocal;
  MPI_Barrier(world);
  double step1_time = MPI_Wtime();
  int nx_ele = nxhi_out - nxlo_out + 1;    // nx_pppm + order + 1;
  int ny_ele = nyhi_out - nylo_out + 1;    // ny_pppm + order + 1;
  int nz_ele = nzhi_out - nzlo_out + 1;    // nz_pppm + order + 1;
  int nxyz = nx_ele * ny_ele * nz_ele;

  double **gw;
  memory->create(gw, nmat, nxyz, "pppm/electrode:gw");
  memset(&(gw[0][0]), 0, (std::size_t) nmat * (std::size_t) nxyz * sizeof(double));

  auto fmod = [](int x, int n) {    // fast unsigned mod
    int r = abs(x);
    while (r >= n) r -= n;
    return r;
  };

  // loops over weights of electrode atoms and weights of complete grid
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  for (int ipos = 0; ipos < nmat; ipos++) {
    double *_noalias xi_ele = x_ele[ipos];
    // new calculation for nx, ny, nz because part2grid available for
    // nlocal, only
    int nix = static_cast<int>((xi_ele[0] - boxlo[0]) * delxinv + shift) - OFFSET;
    int niy = static_cast<int>((xi_ele[1] - boxlo[1]) * delyinv + shift) - OFFSET;
    int niz = static_cast<int>((xi_ele[2] - boxlo[2]) * delzinv + shift) - OFFSET;
    FFT_SCALAR dx = nix + shiftone - (xi_ele[0] - boxlo[0]) * delxinv;
    FFT_SCALAR dy = niy + shiftone - (xi_ele[1] - boxlo[1]) * delyinv;
    FFT_SCALAR dz = niz + shiftone - (xi_ele[2] - boxlo[2]) * delzinv;
    compute_rho1d(dx, dy, dz);
    // DO NOT TOUCH THESE LOOPS!
    // Attempted optimizations (such as calculating parts of indices
    // in the outer loops) have been shown harmful upon benchmarking!
    for (int mjz = nzlo_out; mjz <= nzhi_out; mjz++) {
      for (int ni = nlower; ni <= nupper; ni++) {
        double const iz0 = rho1d[2][ni];
        int const mz = fmod(mjz - ni - niz, nz_pppm);
        for (int mjy = nylo_out; mjy <= nyhi_out; mjy++) {
          for (int mi = nlower; mi <= nupper; mi++) {
            double const iy0 = iz0 * rho1d[1][mi];
            int const my = fmod(mjy - mi - niy, ny_pppm);
            for (int mjx = nxlo_out; mjx <= nxhi_out; mjx++) {
              for (int li = nlower; li <= nupper; li++) {
                double const ix0 = iy0 * rho1d[0][li];
                int const mx = fmod(mjx - li - nix, nx_pppm);
                gw[ipos][nx_ele * ny_ele * (mjz - nzlo_out) + nx_ele * (mjy - nylo_out) +
                         (mjx - nxlo_out)] +=
                    ix0 * greens_real[mz * nx_pppm * ny_pppm + my * nx_pppm + mx];
              }
            }
          }
        }
      }
    }
  }
  MPI_Barrier(world);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(lmp, "step 1 time: {:.4g} s\n", MPI_Wtime() - step1_time);

  // nested loop over electrode atoms i and j and stencil of i
  // in theory could reuse make_rho1d_j here -- but this step is already
  // super-fast
  double step2_time = MPI_Wtime();
  double **x = atom->x;
  for (int i = 0; i < nlocal; i++) {
    int ipos = imat[i];
    if (ipos < 0) continue;
    int nix = part2grid[i][0];
    int niy = part2grid[i][1];
    int niz = part2grid[i][2];
    FFT_SCALAR dix = nix + shiftone - (x[i][0] - boxlo[0]) * delxinv;
    FFT_SCALAR diy = niy + shiftone - (x[i][1] - boxlo[1]) * delyinv;
    FFT_SCALAR diz = niz + shiftone - (x[i][2] - boxlo[2]) * delzinv;
    compute_rho1d(dix, diy, diz);
    for (int jpos = 0; jpos < nmat; jpos++) {
      double aij = 0.;
      for (int ni = nlower; ni <= nupper; ni++) {
        double iz0 = rho1d[2][ni];
        int miz = ni + niz;
        for (int mi = nlower; mi <= nupper; mi++) {
          double iy0 = iz0 * rho1d[1][mi];
          int miy = mi + niy;
          for (int li = nlower; li <= nupper; li++) {
            int mix = li + nix;
            double ix0 = iy0 * rho1d[0][li];
            int miz0 = miz - nzlo_out;
            int miy0 = miy - nylo_out;
            int mix0 = mix - nxlo_out;
            aij += ix0 * gw[jpos][nx_ele * ny_ele * miz0 + nx_ele * miy0 + mix0];
          }
        }
      }
      matrix[ipos][jpos] += aij / volume;
    }
  }
  MPI_Barrier(world);
  memory->destroy(gw);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(lmp, "step 2 time: {:.4g} s\n", MPI_Wtime() - step2_time);
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMElectrode::allocate()
{
  if (slabflag == 1) {
    // EW3Dc dipole correction
    boundcorr = new SlabDipole(lmp);
  } else if (wireflag == 1) {
    // EW3Dc wire correction
    boundcorr = new WireDipole(lmp);
  } else {
    // dummy BoundaryCorrection for ffield
    boundcorr = new BoundaryCorrection(lmp);
  }

  // ----------------------------------------------------------------------
  // code from PPPM::allocate(), altered to use different Grid3d constructor
  // ----------------------------------------------------------------------

  // create ghost grid object for rho and electric field communication
  // returns local owned and ghost grid bounds
  // setup communication patterns and buffers

  gc = new Grid3d(lmp, world, nx_pppm, ny_pppm, nz_pppm, nxlo_in, nxhi_in, nylo_in, nyhi_in,
                  nzlo_in, nzhi_in, nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out);

  gc->setup_comm(ngc_buf1, ngc_buf2);

  if (differentiation_flag)
    npergrid = 1;
  else
    npergrid = 3;

  memory->create(gc_buf1, npergrid * ngc_buf1, "pppm:gc_buf1");
  memory->create(gc_buf2, npergrid * ngc_buf2, "pppm:gc_buf2");

  // tally local grid sizes
  // ngrid = count of owned+ghost grid cells on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  //              same as count of owned grid cells
  // nfft = FFT points in x-pencil FFT decomposition on this proc
  // nfft_both = greater of nfft and nfft_brick

  ngrid = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) * (nzhi_out - nzlo_out + 1);

  nfft_brick = (nxhi_in - nxlo_in + 1) * (nyhi_in - nylo_in + 1) * (nzhi_in - nzlo_in + 1);

  nfft = (nxhi_fft - nxlo_fft + 1) * (nyhi_fft - nylo_fft + 1) * (nzhi_fft - nzlo_fft + 1);

  nfft_both = MAX(nfft, nfft_brick);

  // allocate distributed grid data

  memory->create3d_offset(density_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:density_brick");

  memory->create(density_fft, nfft_both, "pppm:density_fft");
  memory->create(greensfn, nfft_both, "pppm:greensfn");
  memory->create(work1, 2 * nfft_both, "pppm:work1");
  memory->create(work2, 2 * nfft_both, "pppm:work2");
  memory->create(vg, nfft_both, 6, "pppm:vg");

  if (triclinic == 0) {
    memory->create1d_offset(fkx, nxlo_fft, nxhi_fft, "pppm:fkx");
    memory->create1d_offset(fky, nylo_fft, nyhi_fft, "pppm:fky");
    memory->create1d_offset(fkz, nzlo_fft, nzhi_fft, "pppm:fkz");
  } else {
    memory->create(fkx, nfft_both, "pppm:fkx");
    memory->create(fky, nfft_both, "pppm:fky");
    memory->create(fkz, nfft_both, "pppm:fkz");
  }

  if (differentiation_flag == 1) {
    memory->create3d_offset(u_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                            "pppm:u_brick");

    memory->create(sf_precoeff1, nfft_both, "pppm:sf_precoeff1");
    memory->create(sf_precoeff2, nfft_both, "pppm:sf_precoeff2");
    memory->create(sf_precoeff3, nfft_both, "pppm:sf_precoeff3");
    memory->create(sf_precoeff4, nfft_both, "pppm:sf_precoeff4");
    memory->create(sf_precoeff5, nfft_both, "pppm:sf_precoeff5");
    memory->create(sf_precoeff6, nfft_both, "pppm:sf_precoeff6");

  } else {
    memory->create3d_offset(vdx_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                            "pppm:vdx_brick");
    memory->create3d_offset(vdy_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                            "pppm:vdy_brick");
    memory->create3d_offset(vdz_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                            "pppm:vdz_brick");
  }

  // summation coeffs

  order_allocated = order;
  if (!stagger_flag) memory->create(gf_b, order, "pppm:gf_b");
  memory->create2d_offset(rho1d, 3, -order / 2, order / 2, "pppm:rho1d");
  memory->create2d_offset(drho1d, 3, -order / 2, order / 2, "pppm:drho1d");
  memory->create2d_offset(rho_coeff, order, (1 - order) / 2, order / 2, "pppm:rho_coeff");
  memory->create2d_offset(drho_coeff, order, (1 - order) / 2, order / 2, "pppm:drho_coeff");

  // create 2 FFTs and a Remap
  // 1st FFT keeps data in FFT decomposition
  // 2nd FFT returns data in 3d brick decomposition
  // remap takes data from 3d brick to FFT decomposition

  int tmp;

  fft1 = new FFT3d(lmp, world, nx_pppm, ny_pppm, nz_pppm, nxlo_fft, nxhi_fft, nylo_fft, nyhi_fft,
                   nzlo_fft, nzhi_fft, nxlo_fft, nxhi_fft, nylo_fft, nyhi_fft, nzlo_fft, nzhi_fft,
                   0, 0, &tmp, collective_flag);

  fft2 = new FFT3d(lmp, world, nx_pppm, ny_pppm, nz_pppm, nxlo_fft, nxhi_fft, nylo_fft, nyhi_fft,
                   nzlo_fft, nzhi_fft, nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in, 0, 0,
                   &tmp, collective_flag);

  remap = new Remap(lmp, world, nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in, nxlo_fft,
                    nxhi_fft, nylo_fft, nyhi_fft, nzlo_fft, nzhi_fft, 1, 0, 0, FFT_PRECISION,
                    collective_flag);

  // ELECTRODE specific allocations

  memory->create3d_offset(electrolyte_density_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out,
                          nxlo_out, nxhi_out, "pppm/electrode:electrolyte_density_brick");
  memory->create(electrolyte_density_fft, nfft_both, "pppm/electrode:electrolyte_density_fft");

  if (differentiation_flag != 1)
    memory->create3d_offset(u_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                            "pppm:u_brick");
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
-------------------------------------------------------------------------
*/

void PPPMElectrode::deallocate()
{
  delete gc;
  gc = nullptr;
  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);

  if (boundcorr != nullptr) delete boundcorr;
  memory->destroy3d_offset(electrolyte_density_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy(electrolyte_density_fft);

  memory->destroy3d_offset(density_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(u_brick, nzlo_out, nylo_out, nxlo_out);
  if (differentiation_flag == 1) {
    memory->destroy(sf_precoeff1);
    memory->destroy(sf_precoeff2);
    memory->destroy(sf_precoeff3);
    memory->destroy(sf_precoeff4);
    memory->destroy(sf_precoeff5);
    memory->destroy(sf_precoeff6);
  } else {
    memory->destroy3d_offset(vdx_brick, nzlo_out, nylo_out, nxlo_out);
    memory->destroy3d_offset(vdy_brick, nzlo_out, nylo_out, nxlo_out);
    memory->destroy3d_offset(vdz_brick, nzlo_out, nylo_out, nxlo_out);
  }

  memory->destroy(density_fft);
  memory->destroy(greensfn);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(vg);

  memory->destroy1d_offset(fkx, nxlo_fft);
  memory->destroy1d_offset(fky, nylo_fft);
  memory->destroy1d_offset(fkz, nzlo_fft);

  memory->destroy(gf_b);
  if (stagger_flag) gf_b = nullptr;
  memory->destroy2d_offset(rho1d, -order_allocated / 2);
  memory->destroy2d_offset(drho1d, -order_allocated / 2);
  memory->destroy2d_offset(rho_coeff, (1 - order_allocated) / 2);
  memory->destroy2d_offset(drho_coeff, (1 - order_allocated) / 2);

  delete fft1;
  delete fft2;
  delete remap;
  fft1 = nullptr;
  fft2 = nullptr;
  remap = nullptr;
}

void PPPMElectrode::allocate_peratom()
{
  peratom_allocate_flag = 1;

  memory->create3d_offset(v0_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v0_brick");
  memory->create3d_offset(v1_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v1_brick");
  memory->create3d_offset(v2_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v2_brick");
  memory->create3d_offset(v3_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v3_brick");
  memory->create3d_offset(v4_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v4_brick");
  memory->create3d_offset(v5_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v5_brick");

  // use same GC ghost grid object for peratom grid communication
  // but need to reallocate a larger gc_buf1 and gc_buf2

  if (differentiation_flag)
    npergrid = 6;
  else
    npergrid = 7;

  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);
  memory->create(gc_buf1, npergrid * ngc_buf1, "pppm:gc_buf1");
  memory->create(gc_buf2, npergrid * ngc_buf2, "pppm:gc_buf2");
}

/* ----------------------------------------------------------------------
   set global size of PPPM grid = nx,ny,nz_pppm
   used for charge accumulation, FFTs, and electric field interpolation
-------------------------------------------------------------------------
*/

void PPPMElectrode::set_grid_global()
{
  // use xprd,yprd,zprd (even if triclinic, and then scale later)
  // adjust z dimension for 2d slab PPPM
  // 3d PPPM just uses zprd since slab_volfactor = 1.0

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;

  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  double h;
  bigint natoms = atom->natoms;

  if (!gewaldflag) {
    if (accuracy <= 0.0) error->all(FLERR, "KSpace accuracy must be > 0");
    if (q2 == 0.0) error->all(FLERR, "Must use kspace_modify gewald for uncharged system");
    g_ewald = accuracy * sqrt(natoms * cutoff * xprd * yprd * zprd) / (2.0 * q2);
    if (g_ewald >= 1.0)
      g_ewald = (1.35 - 0.15 * log(accuracy)) / cutoff;
    else
      g_ewald = sqrt(-log(g_ewald)) / cutoff;
  }

  // set optimal nx_pppm,ny_pppm,nz_pppm based on order and accuracy
  // nz_pppm uses extended zprd_slab instead of zprd
  // reduce it until accuracy target is met

  if (!gridflag) {
    if (differentiation_flag == 1 || stagger_flag) {
      h = h_x = h_y = h_z = 4.0 / g_ewald;
      int count = 0;
      while (true) {
        // set grid dimensions

        nx_pppm = static_cast<int>(xprd_wire / h_x);
        ny_pppm = static_cast<int>(yprd_wire / h_y);
        nz_pppm = static_cast<int>(zprd_slab / h_z);

        if (nx_pppm <= 1) nx_pppm = 2;
        if (ny_pppm <= 1) ny_pppm = 2;
        if (nz_pppm <= 1) nz_pppm = 2;

        // estimate KSpace force error

        double df_kspace = compute_df_kspace();

        // break loop if the accuracy has been reached or
        // too many loops have been performed

        count++;
        if (df_kspace <= accuracy) break;

        if (count > 500) error->all(FLERR, "Could not compute grid size");
        h *= 0.95;
        h_x = h_y = h_z = h;
      }

    } else {
      double err;
      h_x = h_y = h_z = 1.0 / g_ewald;

      nx_pppm = static_cast<int>(xprd_wire / h_x) + 1;
      ny_pppm = static_cast<int>(yprd_wire / h_y) + 1;
      nz_pppm = static_cast<int>(zprd_slab / h_z) + 1;

      err = estimate_ik_error(h_x, xprd_wire, natoms);
      while (err > accuracy) {
        err = estimate_ik_error(h_x, xprd_wire, natoms);
        nx_pppm++;
        h_x = xprd_wire / nx_pppm;
      }

      err = estimate_ik_error(h_y, yprd_wire, natoms);
      while (err > accuracy) {
        err = estimate_ik_error(h_y, yprd_wire, natoms);
        ny_pppm++;
        h_y = yprd_wire / ny_pppm;
      }

      err = estimate_ik_error(h_z, zprd_slab, natoms);
      while (err > accuracy) {
        err = estimate_ik_error(h_z, zprd_slab, natoms);
        nz_pppm++;
        h_z = zprd_slab / nz_pppm;
      }
    }
  }

  // boost grid size until it is factorable

  while (!factorable(nx_pppm)) nx_pppm++;
  while (!factorable(ny_pppm)) ny_pppm++;
  while (!factorable(nz_pppm)) nz_pppm++;

  h_x = xprd_wire / nx_pppm;
  h_y = yprd_wire / ny_pppm;
  h_z = zprd_slab / nz_pppm;

  if (nx_pppm >= OFFSET || ny_pppm >= OFFSET || nz_pppm >= OFFSET)
    error->all(FLERR, "PPPM/electrode grid is too large");
}

/* ----------------------------------------------------------------------
   compute estimated kspace force error
-------------------------------------------------------------------------
*/

double PPPMElectrode::compute_df_kspace()
{
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;
  bigint natoms = atom->natoms;
  double df_kspace = 0.0;
  if (differentiation_flag == 1 || stagger_flag) {
    double qopt = compute_qopt();
    df_kspace = sqrt(qopt / natoms) * q2 / (xprd_wire * yprd_wire * zprd_slab);
  } else {
    double lprx = estimate_ik_error(h_x, xprd_wire, natoms);
    double lpry = estimate_ik_error(h_y, yprd_wire, natoms);
    double lprz = estimate_ik_error(h_z, zprd_slab, natoms);
    df_kspace = sqrt(lprx * lprx + lpry * lpry + lprz * lprz) / sqrt(3.0);
  }
  return df_kspace;
}

/* ----------------------------------------------------------------------
   compute qopt
-------------------------------------------------------------------------
*/

double PPPMElectrode::compute_qopt()
{
  int k, l, m, nx, ny, nz;
  double argx, argy, argz, wx, wy, wz, sx, sy, sz, qx, qy, qz;
  double u1, u2, sqk;
  double sum1, sum2, sum3, sum4, dot2;

  double *prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double xprd_wire = xprd * wire_volfactor;
  const double yprd_wire = yprd * wire_volfactor;
  const double zprd_slab = zprd * slab_volfactor;
  volume = xprd_wire * yprd_wire * zprd_slab;

  const double unitkx = (MY_2PI / xprd_wire);
  const double unitky = (MY_2PI / yprd_wire);
  const double unitkz = (MY_2PI / zprd_slab);

  const int twoorder = 2 * order;

  // loop over entire FFT grid
  // each proc calculates contributions from every Pth grid point

  bigint ngridtotal = (bigint) nx_pppm * ny_pppm * nz_pppm;
  bigint nxy_pppm = (bigint) nx_pppm * ny_pppm;

  double qopt = 0.0;

  for (bigint i = me; i < ngridtotal; i += nprocs) {
    k = i % nx_pppm;
    l = (i / nx_pppm) % ny_pppm;
    m = i / nxy_pppm;

    const int kper = k - nx_pppm * (2 * k / nx_pppm);
    const int lper = l - ny_pppm * (2 * l / ny_pppm);
    const int mper = m - nz_pppm * (2 * m / nz_pppm);

    sqk = square(unitkx * kper) + square(unitky * lper) + square(unitkz * mper);
    if (sqk == 0.0) continue;

    sum1 = sum2 = sum3 = sum4 = 0.0;

    for (nx = -2; nx <= 2; nx++) {
      qx = unitkx * (kper + nx_pppm * nx);
      sx = exp(-0.25 * square(qx / g_ewald));
      argx = 0.5 * qx * xprd_wire / nx_pppm;
      wx = powsinxx(argx, twoorder);
      qx *= qx;

      for (ny = -2; ny <= 2; ny++) {
        qy = unitky * (lper + ny_pppm * ny);
        sy = exp(-0.25 * square(qy / g_ewald));
        argy = 0.5 * qy * yprd_wire / ny_pppm;
        wy = powsinxx(argy, twoorder);
        qy *= qy;

        for (nz = -2; nz <= 2; nz++) {
          qz = unitkz * (mper + nz_pppm * nz);
          sz = exp(-0.25 * square(qz / g_ewald));
          argz = 0.5 * qz * zprd_slab / nz_pppm;
          wz = powsinxx(argz, twoorder);
          qz *= qz;

          dot2 = qx + qy + qz;
          u1 = sx * sy * sz;
          u2 = wx * wy * wz;

          sum1 += u1 * u1 / dot2 * MY_4PI * MY_4PI;
          sum2 += u1 * u2 * MY_4PI;
          sum3 += u2;
          sum4 += dot2 * u2;
        }
      }
    }

    sum2 *= sum2;
    qopt += sum1 - sum2 / (sum3 * sum4);
  }

  // sum qopt over all procs

  double qopt_all;
  MPI_Allreduce(&qopt, &qopt_all, 1, MPI_DOUBLE, MPI_SUM, world);
  return qopt_all;
}

/* ----------------------------------------------------------------------
   set local subset of PPPM/FFT grid that I own
   n xyz lo/hi in = 3d brick that I own (inclusive)
   n xyz lo/hi out = 3d brick + ghost cells in 6 directions (inclusive)
   n xyz lo/hi fft = FFT columns that I own (all of x dim, 2d decomp in yz)
------------------------------------------------------------------------- */

void PPPMElectrode::set_grid_local()
{
  // global indices of PPPM grid range from 0 to N-1
  // nlo_in,nhi_in = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that I own without ghost cells
  // for slab PPPM, assign z grid as if it were not extended
  // both non-tiled and tiled proc layouts use 0-1 fractional sumdomain
  // info

  if (comm->layout != Comm::LAYOUT_TILED) {
    nxlo_in = static_cast<int>(comm->xsplit[comm->myloc[0]] * nx_pppm / wire_volfactor);
    nxhi_in = static_cast<int>(comm->xsplit[comm->myloc[0] + 1] * nx_pppm / wire_volfactor) - 1;

    nylo_in = static_cast<int>(comm->ysplit[comm->myloc[1]] * ny_pppm / wire_volfactor);
    nyhi_in = static_cast<int>(comm->ysplit[comm->myloc[1] + 1] * ny_pppm / wire_volfactor) - 1;

    nzlo_in = static_cast<int>(comm->zsplit[comm->myloc[2]] * nz_pppm / slab_volfactor);
    nzhi_in = static_cast<int>(comm->zsplit[comm->myloc[2] + 1] * nz_pppm / slab_volfactor) - 1;

  } else {
    nxlo_in = static_cast<int>(comm->mysplit[0][0] * nx_pppm / wire_volfactor);
    nxhi_in = static_cast<int>(comm->mysplit[0][1] * nx_pppm / wire_volfactor) - 1;

    nylo_in = static_cast<int>(comm->mysplit[1][0] * ny_pppm / wire_volfactor);
    nyhi_in = static_cast<int>(comm->mysplit[1][1] * ny_pppm / wire_volfactor) - 1;

    nzlo_in = static_cast<int>(comm->mysplit[2][0] * nz_pppm / slab_volfactor);
    nzhi_in = static_cast<int>(comm->mysplit[2][1] * nz_pppm / slab_volfactor) - 1;
  }

  // nlower,nupper = stencil size for mapping particles to PPPM grid

  nlower = -(order - 1) / 2;
  nupper = order / 2;

  // shift values for particle <-> grid mapping
  // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

  if (order % 2)
    shift = OFFSET + 0.5;
  else
    shift = OFFSET;

  if (order % 2)
    shiftone = 0.0;
  else
    shiftone = 0.5;

  // nlo_out,nhi_out = lower/upper limits of the 3d sub-brick of
  //   global PPPM grid that my particles can contribute charge to
  // effectively nlo_in,nhi_in + ghost cells
  // nlo,nhi = global coords of grid pt to "lower left" of
  // smallest/largest
  //           position a particle in my box can be at
  // dist[3] = particle position bound = subbox + skin/2.0 + qdist
  //   qdist = offset due to TIP4P fictitious charge
  //   convert to triclinic if necessary
  // nlo_out,nhi_out = nlo,nhi + stencil size for particle mapping
  // for slab PPPM, assign z grid as if it were not extended

  double *prd, *sublo, *subhi;

  prd = domain->prd;
  boxlo = domain->boxlo;
  sublo = domain->sublo;
  subhi = domain->subhi;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;

  double dist[3] = {0.0, 0.0, 0.0};
  double cuthalf = 0.5 * neighbor->skin + qdist;
  dist[0] = dist[1] = dist[2] = cuthalf;

  int nlo, nhi;
  nlo = nhi = 0;

  nlo = static_cast<int>((sublo[0] - dist[0] - boxlo[0]) * nx_pppm / xprd_wire + shift) - OFFSET;
  nhi = static_cast<int>((subhi[0] + dist[0] - boxlo[0]) * nx_pppm / xprd_wire + shift) - OFFSET;
  nxlo_out = nlo + nlower;
  nxhi_out = nhi + nupper;

  nlo = static_cast<int>((sublo[1] - dist[1] - boxlo[1]) * ny_pppm / yprd_wire + shift) - OFFSET;
  nhi = static_cast<int>((subhi[1] + dist[1] - boxlo[1]) * ny_pppm / yprd_wire + shift) - OFFSET;
  nylo_out = nlo + nlower;
  nyhi_out = nhi + nupper;

  nlo = static_cast<int>((sublo[2] - dist[2] - boxlo[2]) * nz_pppm / zprd_slab + shift) - OFFSET;
  nhi = static_cast<int>((subhi[2] + dist[2] - boxlo[2]) * nz_pppm / zprd_slab + shift) - OFFSET;
  nzlo_out = nlo + nlower;
  nzhi_out = nhi + nupper;

  if (stagger_flag) {
    nxhi_out++;
    nyhi_out++;
    nzhi_out++;
  }

  // for slab PPPM, change the grid boundary for processors at +z end
  //   to include the empty volume between periodically repeating slabs
  // for slab PPPM, want charge data communicated from -z proc to +z proc,
  //   but not vice versa, also want field data communicated from +z proc
  //   to -z proc, but not vice versa
  // this is accomplished by nzhi_in = nzhi_out on +z end (no ghost cells)
  // also ensure no other procs use ghost cells beyond +z limit
  // differnet logic for non-tiled vs tiled decomposition

  if (slabflag == 1) {
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[2] == comm->procgrid[2] - 1) nzhi_in = nzhi_out = nz_pppm - 1;
    } else {
      if (comm->mysplit[2][1] == 1.0) nzhi_in = nzhi_out = nz_pppm - 1;
    }
    nzhi_out = MIN(nzhi_out, nz_pppm - 1);
  }

  if (wireflag == 1) {
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[0] == comm->procgrid[0] - 1) nxhi_in = nxhi_out = nx_pppm - 1;
      if (comm->myloc[1] == comm->procgrid[1] - 1) nyhi_in = nyhi_out = ny_pppm - 1;
    } else {
      if (comm->mysplit[0][1] == 1.0) nxhi_in = nxhi_out = nx_pppm - 1;
      if (comm->mysplit[1][1] == 1.0) nyhi_in = nyhi_out = ny_pppm - 1;
    }
    nxhi_out = MIN(nxhi_out, nx_pppm - 1);
    nyhi_out = MIN(nyhi_out, ny_pppm - 1);
  }

  // x-pencil decomposition of FFT mesh
  // global indices range from 0 to N-1
  // each proc owns entire x-dimension, clumps of columns in y,z
  // dimensions npey_fft,npez_fft = # of procs in y,z dims if nprocs is
  // small enough, proc can own 1 or more entire xy planes,
  //   else proc owns 2d sub-blocks of yz plane
  // me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
  // nlo_fft,nhi_fft = lower/upper limit of the section
  //   of the global FFT mesh that I own in x-pencil decomposition

  int npey_fft, npez_fft;
  if (nz_pppm >= nprocs) {
    npey_fft = 1;
    npez_fft = nprocs;
  } else
    procs2grid2d(nprocs, ny_pppm, nz_pppm, &npey_fft, &npez_fft);

  int me_y = me % npey_fft;
  int me_z = me / npey_fft;

  nxlo_fft = 0;
  nxhi_fft = nx_pppm - 1;
  nylo_fft = me_y * ny_pppm / npey_fft;
  nyhi_fft = (me_y + 1) * ny_pppm / npey_fft - 1;
  nzlo_fft = me_z * nz_pppm / npez_fft;
  nzhi_fft = (me_z + 1) * nz_pppm / npez_fft - 1;

  // ngrid = count of PPPM grid pts owned by this proc, including ghosts

  ngrid = (nxhi_out - nxlo_out + 1) * (nyhi_out - nylo_out + 1) * (nzhi_out - nzlo_out + 1);

  // count of FFT grids pts owned by this proc, without ghosts
  // nfft = FFT points in x-pencil FFT decomposition on this proc
  // nfft_brick = FFT points in 3d brick-decomposition on this proc
  // nfft_both = greater of 2 values

  nfft = (nxhi_fft - nxlo_fft + 1) * (nyhi_fft - nylo_fft + 1) * (nzhi_fft - nzlo_fft + 1);
  int nfft_brick = (nxhi_in - nxlo_in + 1) * (nyhi_in - nylo_in + 1) * (nzhi_in - nzlo_in + 1);
  nfft_both = MAX(nfft, nfft_brick);
}

/* ----------------------------------------------------------------------
   pre-compute modified (Hockney-Eastwood) Coulomb Green's function
-------------------------------------------------------------------------
*/

void PPPMElectrode::compute_gf_ik()
{
  const double *const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double xprd_wire = xprd * wire_volfactor;
  const double yprd_wire = yprd * wire_volfactor;
  const double zprd_slab = zprd * slab_volfactor;
  const double unitkx = (MY_2PI / xprd_wire);
  const double unitky = (MY_2PI / yprd_wire);
  const double unitkz = (MY_2PI / zprd_slab);

  double snx, sny, snz;
  double argx, argy, argz, wx, wy, wz, sx, sy, sz, qx, qy, qz;
  double sum1, dot1, dot2;
  double numerator, denominator;
  double sqk;

  int k, l, m, n, nx, ny, nz, kper, lper, mper;

  const int nbx =
      static_cast<int>((g_ewald * xprd_wire / (MY_PI * nx_pppm)) * pow(-log(EPS_HOC), 0.25));
  const int nby =
      static_cast<int>((g_ewald * yprd_wire / (MY_PI * ny_pppm)) * pow(-log(EPS_HOC), 0.25));
  const int nbz =
      static_cast<int>((g_ewald * zprd_slab / (MY_PI * nz_pppm)) * pow(-log(EPS_HOC), 0.25));
  const int twoorder = 2 * order;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm * (2 * m / nz_pppm);
    snz = square(sin(0.5 * unitkz * mper * zprd_slab / nz_pppm));

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm * (2 * l / ny_pppm);
      sny = square(sin(0.5 * unitky * lper * yprd_wire / ny_pppm));

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm * (2 * k / nx_pppm);
        snx = square(sin(0.5 * unitkx * kper * xprd_wire / nx_pppm));

        sqk = square(unitkx * kper) + square(unitky * lper) + square(unitkz * mper);

        if (sqk != 0.0) {
          numerator = 12.5663706 / sqk;    // 4 pi / k^2
          denominator = gf_denom(snx, sny, snz);
          sum1 = 0.0;

          for (nx = -nbx; nx <= nbx; nx++) {
            qx = unitkx * (kper + nx_pppm * nx);
            sx = exp(-0.25 * square(qx / g_ewald));
            argx = 0.5 * qx * xprd_wire / nx_pppm;
            wx = powsinxx(argx, twoorder);

            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky * (lper + ny_pppm * ny);
              sy = exp(-0.25 * square(qy / g_ewald));
              argy = 0.5 * qy * yprd_wire / ny_pppm;
              wy = powsinxx(argy, twoorder);

              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz * (mper + nz_pppm * nz);
                sz = exp(-0.25 * square(qz / g_ewald));
                argz = 0.5 * qz * zprd_slab / nz_pppm;
                wz = powsinxx(argz, twoorder);

                dot1 = unitkx * kper * qx + unitky * lper * qy + unitkz * mper * qz;
                dot2 = qx * qx + qy * qy + qz * qz;
                sum1 += (dot1 / dot2) * sx * sy * sz * wx * wy * wz;
              }
            }
          }
          greensfn[n++] = numerator * sum1 / denominator;
        } else
          greensfn[n++] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute optimized Green's function for energy calculation
-------------------------------------------------------------------------
*/

void PPPMElectrode::compute_gf_ad()
{
  const double *const prd = domain->prd;

  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];
  const double xprd_wire = xprd * wire_volfactor;
  const double yprd_wire = yprd * wire_volfactor;
  const double zprd_slab = zprd * slab_volfactor;
  const double unitkx = (MY_2PI / xprd_wire);
  const double unitky = (MY_2PI / yprd_wire);
  const double unitkz = (MY_2PI / zprd_slab);

  double snx, sny, snz, sqk;
  double argx, argy, argz, wx, wy, wz, sx, sy, sz, qx, qy, qz;
  double numerator, denominator;
  int k, l, m, n, kper, lper, mper;

  const int twoorder = 2 * order;

  for (int i = 0; i < 6; i++) sf_coeff[i] = 0.0;

  n = 0;
  for (m = nzlo_fft; m <= nzhi_fft; m++) {
    mper = m - nz_pppm * (2 * m / nz_pppm);
    qz = unitkz * mper;
    snz = square(sin(0.5 * qz * zprd_slab / nz_pppm));
    sz = exp(-0.25 * square(qz / g_ewald));
    argz = 0.5 * qz * zprd_slab / nz_pppm;
    wz = powsinxx(argz, twoorder);

    for (l = nylo_fft; l <= nyhi_fft; l++) {
      lper = l - ny_pppm * (2 * l / ny_pppm);
      qy = unitky * lper;
      sny = square(sin(0.5 * qy * yprd_wire / ny_pppm));
      sy = exp(-0.25 * square(qy / g_ewald));
      argy = 0.5 * qy * yprd_wire / ny_pppm;
      wy = powsinxx(argy, twoorder);

      for (k = nxlo_fft; k <= nxhi_fft; k++) {
        kper = k - nx_pppm * (2 * k / nx_pppm);
        qx = unitkx * kper;
        snx = square(sin(0.5 * qx * xprd_wire / nx_pppm));
        sx = exp(-0.25 * square(qx / g_ewald));
        argx = 0.5 * qx * xprd_wire / nx_pppm;
        wx = powsinxx(argx, twoorder);

        sqk = qx * qx + qy * qy + qz * qz;

        if (sqk != 0.0) {
          numerator = MY_4PI / sqk;
          denominator = gf_denom(snx, sny, snz);
          greensfn[n] = numerator * sx * sy * sz * wx * wy * wz / denominator;
          sf_coeff[0] += sf_precoeff1[n] * greensfn[n];
          sf_coeff[1] += sf_precoeff2[n] * greensfn[n];
          sf_coeff[2] += sf_precoeff3[n] * greensfn[n];
          sf_coeff[3] += sf_precoeff4[n] * greensfn[n];
          sf_coeff[4] += sf_precoeff5[n] * greensfn[n];
          sf_coeff[5] += sf_precoeff6[n] * greensfn[n];
          n++;
        } else {
          greensfn[n] = 0.0;
          sf_coeff[0] += sf_precoeff1[n] * greensfn[n];
          sf_coeff[1] += sf_precoeff2[n] * greensfn[n];
          sf_coeff[2] += sf_precoeff3[n] * greensfn[n];
          sf_coeff[3] += sf_precoeff4[n] * greensfn[n];
          sf_coeff[4] += sf_precoeff5[n] * greensfn[n];
          sf_coeff[5] += sf_precoeff6[n] * greensfn[n];
          n++;
        }
      }
    }
  }

  // compute the coefficients for the self-force correction

  double prex, prey, prez;
  prex = prey = prez = MY_PI / volume;
  prex *= nx_pppm / xprd_wire;
  prey *= ny_pppm / yprd_wire;
  prez *= nz_pppm / zprd_slab;
  sf_coeff[0] *= prex;
  sf_coeff[1] *= prex * 2;
  sf_coeff[2] *= prey;
  sf_coeff[3] *= prey * 2;
  sf_coeff[4] *= prez;
  sf_coeff[5] *= prez * 2;

  // communicate values with other procs

  double tmp[6];
  MPI_Allreduce(sf_coeff, tmp, 6, MPI_DOUBLE, MPI_SUM, world);
  for (n = 0; n < 6; n++) sf_coeff[n] = tmp[n];
}

/* ----------------------------------------------------------------------
   create discretized "density" of a set of particles (c.f. make_rho())
   in a specified scratch space.
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including
ghosts) in global grid
-------------------------------------------------------------------------
*/

void PPPMElectrode::make_rho_in_brick(int source_grpbit, FFT_SCALAR ***scratch_brick,
                                      bool invert_source)
{
  int l, m, n, nx, ny, nz, mx, my, mz;
  FFT_SCALAR dx, dy, dz, x0, y0, z0;

  last_source_grpbit = source_grpbit;
  last_invert_source = invert_source;

  // clear 3d density array
  memset(&(scratch_brick[nzlo_out][nylo_out][nxlo_out]), 0, ngrid * sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    bool const i_in_source = !!(mask[i] & source_grpbit) != invert_source;
    if (!i_in_source) continue;
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx + shiftone - (x[i][0] - boxlo[0]) * delxinv;
    dy = ny + shiftone - (x[i][1] - boxlo[1]) * delyinv;
    dz = nz + shiftone - (x[i][2] - boxlo[2]) * delzinv;

    compute_rho1d(dx, dy, dz);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n + nz;
      y0 = z0 * rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m + ny;
        x0 = y0 * rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l + nx;
          scratch_brick[mz][my][mx] += x0 * rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   group-group interactions
 -------------------------------------------------------------------------
*/
void PPPMElectrode::compute_group_group(int /*groupbit_A*/, int /*groupbit_B*/, int /*AA_flag*/)
{
  error->all(FLERR, "group group interaction not implemented in pppm/electrode yet");
}

void PPPMElectrode::compute_matrix_corr(bigint *imat, double **matrix)
{
  boundcorr->matrix_corr(imat, matrix);
}

/* ----------------------------------------------------------------------
   compute b-vector EW3DC correction of constant potential approach
 -------------------------------------------------------------------------
*/

void PPPMElectrode::compute_vector_corr(double *vec, int sensor_grpbit, int source_grpbit,
                                        bool invert_source)
{
  boundcorr->vector_corr(vec, sensor_grpbit, source_grpbit, invert_source);
}
