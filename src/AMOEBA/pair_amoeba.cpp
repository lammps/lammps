// clang-format off
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

#include "pair_amoeba.h"

#include "amoeba_convolution.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "timer.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

using MathSpecial::powint;

enum{INDUCE,RSD,SETUP_AMOEBA,SETUP_HIPPO,KMPOLE,AMGROUP,PVAL};  // forward comm
enum{FIELD,ZRSD,TORQUE,UFLD};                                   // reverse comm
enum{ARITHMETIC,GEOMETRIC,CUBIC_MEAN,R_MIN,SIGMA,DIAMETER,HARMONIC,HHG,W_H};
enum{HAL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};
enum{MPOLE_GRID,POLAR_GRID,POLAR_GRIDC,DISP_GRID,INDUCE_GRID,INDUCE_GRIDC};
enum{MUTUAL,OPT,TCG,DIRECT};
enum{GEAR,ASPC,LSQR};

#define DELTASTACK 16
#define DEBUG_AMOEBA 0

/* ---------------------------------------------------------------------- */

PairAmoeba::PairAmoeba(LAMMPS *lmp) : Pair(lmp)
{
  amoeba = true;
  mystyle = "amoeba";

  // pair style settings

  one_coeff = 1;
  single_enable = 0;
  no_virial_fdotr_compute = 1;

  nextra = 6;
  pvector = new double[nextra];

  // force field settings

  nmax = 0;
  xaxis2local = yaxis2local = zaxis2local = nullptr;
  rpole = nullptr;
  tq = nullptr;

  red2local = nullptr;
  xred = nullptr;

  uind = uinp = udirp = nullptr;
  uopt = uoptp = nullptr;
  fopt = foptp = nullptr;
  field = fieldp = nullptr;
  ufld = dufld = nullptr;
  rsd = rsdp = nullptr;
  zrsd = zrsdp = nullptr;

  cmp = fmp = nullptr;
  cphi = fphi = nullptr;

  _moduli_array = nullptr;
  _moduli_bsarray = nullptr;
  _nfft_max = 0;

  poli = nullptr;
  conj = conjp = nullptr;
  vec = vecp = nullptr;
  udir = usum = usump = nullptr;

  fuind = fuinp = nullptr;
  fdip_phi1 = fdip_phi2 = fdip_sum_phi = nullptr;
  dipfield1 = dipfield2 = nullptr;

  fphid = fphip = nullptr;
  fphidp = cphidp = nullptr;

  bsordermax = 0;
  thetai1 = thetai2 = thetai3 = nullptr;
  bsmod1 = bsmod2 = bsmod3 = nullptr;
  bsbuild = nullptr;
  igrid = nullptr;
  m_kspace = p_kspace = pc_kspace = d_kspace = nullptr;
  i_kspace = ic_kspace = nullptr;

  numneigh_dipole = nullptr;
  firstneigh_dipole = nullptr;
  firstneigh_dipdip = nullptr;
  ipage_dipole = nullptr;
  dpage_dipdip = nullptr;

  numneigh_precond = nullptr;
  firstneigh_precond = nullptr;
  ipage_precond = nullptr;

  firstneigh_pcpc = nullptr;
  dpage_pcpc = nullptr;

  qfac = nullptr;
  gridfft1 = nullptr;

  initialize_type_class();
  initialize_vdwl();
  initialize_smallsize();

  forcefield = nullptr;

  id_pole = id_udalt = id_upalt = nullptr;

  memset(special_hal, 0 , sizeof(special_hal));
  memset(special_repel, 0 , sizeof(special_repel));
  memset(special_disp, 0 , sizeof(special_disp));
  memset(special_mpole, 0 , sizeof(special_mpole));
  memset(special_polar_pscale, 0 , sizeof(special_polar_pscale));
  memset(special_polar_piscale, 0 , sizeof(special_polar_piscale));
  memset(special_polar_wscale, 0 , sizeof(special_polar_wscale));

  nualt = 0;
  first_flag = 1;
  first_flag_compute = 1;

  // use Tinker value = 332.063713 (one extra digit)
  // LAMMPS value = 332.06371

  electric = 332.063713;
  //electric = force->qqr2e;

  // factors for FFT grid size

  nfactors = 3;
  factors = new int[nfactors];
  factors[0] = 2;
  factors[1] = 3;
  factors[2] = 5;
}

/* ---------------------------------------------------------------------- */

PairAmoeba::~PairAmoeba()
{
  delete[] pvector;

  // check nfix in case all fixes have already been deleted

  if (modify->nfix) {
    if (id_pole) modify->delete_fix(id_pole);
    if (id_udalt) modify->delete_fix(id_udalt);
    if (id_upalt) modify->delete_fix(id_upalt);
  }

  delete[] id_pole;
  delete[] id_udalt;
  delete[] id_upalt;

  memory->destroy(xaxis2local);
  memory->destroy(yaxis2local);
  memory->destroy(zaxis2local);
  memory->destroy(rpole);
  memory->destroy(tq);

  memory->destroy(red2local);
  memory->destroy(xred);

  memory->destroy(uind);
  memory->destroy(uinp);
  memory->destroy(udirp);
  memory->destroy(uopt);
  memory->destroy(uoptp);
  memory->destroy(fopt);
  memory->destroy(foptp);

  memory->destroy(field);
  memory->destroy(fieldp);
  memory->destroy(ufld);
  memory->destroy(dufld);
  memory->destroy(rsd);
  memory->destroy(rsdp);
  memory->destroy(zrsd);
  memory->destroy(zrsdp);

  memory->destroy(cmp);
  memory->destroy(fmp);
  memory->destroy(cphi);
  memory->destroy(fphi);

  memory->destroy(poli);
  memory->destroy(conj);
  memory->destroy(conjp);
  memory->destroy(vec);
  memory->destroy(vecp);
  memory->destroy(udir);
  memory->destroy(usum);
  memory->destroy(usump);

  memory->destroy(fuind);
  memory->destroy(fuinp);
  memory->destroy(fdip_phi1);
  memory->destroy(fdip_phi2);
  memory->destroy(fdip_sum_phi);
  memory->destroy(dipfield1);
  memory->destroy(dipfield2);

  memory->destroy(fphid);
  memory->destroy(fphip);
  memory->destroy(fphidp);
  memory->destroy(cphidp);

  memory->destroy(_moduli_array);
  memory->destroy(_moduli_bsarray);

  memory->destroy(thetai1);
  memory->destroy(thetai2);
  memory->destroy(thetai3);
  memory->destroy(igrid);

  memory->destroy(bsmod1);
  memory->destroy(bsmod2);
  memory->destroy(bsmod3);
  memory->destroy(bsbuild);

  memory->destroy(qfac);
  memory->destroy(gridfft1);

  delete m_kspace;
  delete p_kspace;
  delete pc_kspace;
  delete d_kspace;
  delete i_kspace;
  delete ic_kspace;

  memory->destroy(numneigh_dipole);
  memory->sfree(firstneigh_dipole);
  memory->sfree(firstneigh_dipdip);
  delete ipage_dipole;
  delete dpage_dipdip;

  memory->destroy(numneigh_precond);
  memory->sfree(firstneigh_precond);
  delete ipage_precond;

  memory->sfree(firstneigh_pcpc);
  delete dpage_pcpc;

  deallocate_type_class();
  if (amoeba) deallocate_vdwl();
  deallocate_smallsize();

  delete[] forcefield;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

  delete[] factors;
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  if (eflag_atom || vflag_atom)
    error->all(FLERR,"Cannot (yet) compute per-atom energy/virial with pair_style {}", mystyle);

  // zero energy/virial components

  ehal = erepulse = edisp = epolar = empole = eqxfer = 0.0;

  for (int i = 0; i < 6; i++) {
    virhal[i] = 0.0;
    virrepulse[i] = 0.0;
    virdisp[i] = 0.0;
    virpolar[i] = 0.0;
    virmpole[i] = 0.0;
    virqxfer[i] = 0.0;
  }

  // grow local vectors and arrays if necessary

  if (atom->nmax > nmax) grow_local();

  // set amtype/amgroup ptrs for rest of compute() to use it

  amtype = atom->ivector[index_amtype];
  amgroup = atom->ivector[index_amgroup];

  // -------------------------------------------------------------------
  // one-time initializations
  // can't do in init_style() b/c these operations require communication
  // -------------------------------------------------------------------

  // assignment of atoms to polarization groups

  if (first_flag_compute) assign_groups();

  // assigmment of multipole neighbors to each owned atom
  // sets xaxis,yaxis,zaxis
  // for HIPPO, also set pval for each atom, then ghost comm of pval

  if (first_flag_compute) {
    cfstyle = KMPOLE;
    comm->forward_comm(this);
    kmpole();

    if (!amoeba) {
      double *pval = atom->dvector[index_pval];
      double **pole = fixpole->astore;
      int nlocal = atom->nlocal;
      int itype,iclass;
      for (int i = 0; i < nlocal; i++) {
        itype = amtype[i];
        iclass = amtype2class[itype];
        pval[i] = pole[i][0] - pcore[iclass];
      }
      cfstyle = PVAL;
      comm->forward_comm(this);
    }
  }

  first_flag_compute = 0;

  // -------------------------------------------------------------------
  // end of one-time initializations
  // -------------------------------------------------------------------

  // initialize timers on first compute() call after setup

  if (update->ntimestep <= update->beginstep+1) {
    time_init = time_hal = time_repulse = time_disp = time_mpole = 0.0;
    time_induce = time_polar = time_qxfer = 0.0;

    time_mpole_rspace = time_mpole_kspace = 0.0;
    time_direct_rspace = time_direct_kspace = 0.0;
    time_mutual_rspace = time_mutual_kspace = 0.0;
    time_polar_rspace = time_polar_kspace = 0.0;

    time_grid_uind = time_fphi_uind = 0.0;
    if (ic_kspace) {
      ic_kspace->time_fft = 0.0;
    }
  }

  double time0,time1,time2,time3,time4,time5,time6,time7,time8;

  if (timer->has_sync()) MPI_Barrier(world);
  time0 = platform::walltime();

  // if reneighboring step:
  // augment neighbor list to include 1-5 neighbor flags
  // re-create red2local and xyz axis2local
  // re-create induce neighbor list

  if (neighbor->ago == 0) {
    add_onefive_neighbors();
    if (amoeba) find_hydrogen_neighbors();
    find_multipole_neighbors();
    if (poltyp == MUTUAL && pcgprec) precond_neigh();
  }

  // reset KSpace recip matrix if box size/shape change dynamically

  if (domain->box_change) lattice();

  // compute reduced H coords for owned atoms
  // needs to be computed before forward_comm with cfstyle = SETUP

  if (amoeba) {
    int j,iclass;
    double rdn;

    double **x = atom->x;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      j = red2local[i];
      iclass = amtype2class[amtype[i]];
      rdn = kred[iclass];
      xred[i][0] = rdn*(x[i][0]-x[j][0]) + x[j][0];
      xred[i][1] = rdn*(x[i][1]-x[j][1]) + x[j][1];
      xred[i][2] = rdn*(x[i][2]-x[j][2]) + x[j][2];
    }
  }

  // compute rpole for owned atoms

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    chkpole(i);
    rotmat(i);
    rotsite(i);
  }

  // communicate quantities needed by ghost atoms: xred, rpole
  // for xred, need to account for PBC

  if (amoeba) cfstyle = SETUP_AMOEBA;
  else cfstyle = SETUP_HIPPO;
  comm->forward_comm(this);

  if (amoeba) pbc_xred();
  time1 = platform::walltime();

  // ----------------------------------------
  // compute components of force field
  // ----------------------------------------

  // buffered 14-7 Vdwl, pairwise

  if (amoeba && hal_flag) hal();
  time2 = platform::walltime();

  // Pauli repulsion, pairwise

  if (!amoeba && repulse_flag) repulsion();
  time3 = platform::walltime();

  // Ewald dispersion, pairwise and long range

  if (!amoeba && (disp_rspace_flag || disp_kspace_flag)) dispersion();
  time4 = platform::walltime();

  // multipole, pairwise and long range

  if (mpole_rspace_flag || mpole_kspace_flag) multipole();
  time5 = platform::walltime();

  // induced dipoles, interative CG relaxation
  // communicate induce() output values needed by ghost atoms

  if (polar_rspace_flag || polar_kspace_flag) {
    induce();
    cfstyle = INDUCE;
    comm->forward_comm(this);
  }
  time6 = platform::walltime();

  // dipoles, pairwise and long range

  if (polar_rspace_flag || polar_kspace_flag) polar();
  time7 = platform::walltime();

  // charge transfer, pairwise

  if (!amoeba && qxfer_flag) charge_transfer();
  time8 = platform::walltime();

  // store energy components for output by compute pair command

  pvector[0] = ehal;
  pvector[1] = erepulse;
  pvector[2] = edisp;
  pvector[3] = empole;
  pvector[4] = epolar;
  pvector[5] = eqxfer;

  // energy & virial summations

  eng_vdwl = ehal + edisp;
  eng_coul = erepulse + empole + epolar + eqxfer;

  for (int i = 0; i < 6; i++)
    virial[i] = virhal[i] + virrepulse[i] + virdisp[i] +
      virpolar[i] + virmpole[i] + virqxfer[i];

  // accumulate timing information

  time_init    += time1 - time0;
  time_hal     += time2 - time1;
  time_repulse += time3 - time2;
  time_disp    += time4 - time3;
  time_mpole   += time5 - time4;
  time_induce  += time6 - time5;
  time_polar   += time7 - time6;
  time_qxfer   += time8 - time7;
}

/* ----------------------------------------------------------------------
   print out AMOEBA/HIPPO timing info at end of run
------------------------------------------------------------------------- */

void PairAmoeba::finish()
{
  double ave;
  MPI_Allreduce(&time_init,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_init = ave/comm->nprocs;

  MPI_Allreduce(&time_hal,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_hal = ave/comm->nprocs;

  MPI_Allreduce(&time_repulse,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_repulse = ave/comm->nprocs;

  MPI_Allreduce(&time_disp,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_disp = ave/comm->nprocs;

  MPI_Allreduce(&time_mpole,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_mpole = ave/comm->nprocs;

  MPI_Allreduce(&time_induce,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_induce = ave/comm->nprocs;

  MPI_Allreduce(&time_polar,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_polar = ave/comm->nprocs;

  MPI_Allreduce(&time_qxfer,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_qxfer = ave/comm->nprocs;

  #if DEBUG_AMOEBA
  // real-space/kspace breakdown
  MPI_Allreduce(&time_mpole_rspace,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_mpole_rspace = ave/comm->nprocs;

  MPI_Allreduce(&time_mpole_kspace,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_mpole_kspace = ave/comm->nprocs;

  MPI_Allreduce(&time_direct_rspace,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_direct_rspace = ave/comm->nprocs;

  MPI_Allreduce(&time_direct_kspace,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_direct_kspace = ave/comm->nprocs;

  MPI_Allreduce(&time_mutual_rspace,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_mutual_rspace = ave/comm->nprocs;

  MPI_Allreduce(&time_mutual_kspace,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_mutual_kspace = ave/comm->nprocs;

  MPI_Allreduce(&time_polar_rspace,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_polar_rspace = ave/comm->nprocs;

  MPI_Allreduce(&time_polar_kspace,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_polar_kspace = ave/comm->nprocs;

  MPI_Allreduce(&time_grid_uind,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_grid_uind = ave/comm->nprocs;

  MPI_Allreduce(&time_fphi_uind,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_fphi_uind = ave/comm->nprocs;

  double time_mutual_fft = 0;
  if (ic_kspace) time_mutual_fft = ic_kspace->time_fft;
  MPI_Allreduce(&time_mutual_fft,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_mutual_fft = ave/comm->nprocs;
  #endif // DEBUG_AMOEBA

  double time_total = (time_init + time_hal + time_repulse + time_disp +
                       time_mpole + time_induce + time_polar + time_qxfer) / 100.0;

  if (comm->me == 0) {
    utils::logmesg(lmp,"\n{} timing breakdown:\n", utils::uppercase(mystyle));
    utils::logmesg(lmp,"  Init    time: {:<12.6g} {:6.2f}%\n", time_init, time_init/time_total);
    if (amoeba) {
      utils::logmesg(lmp,"  Hal     time: {:<12.6g} {:6.2f}%\n", time_hal, time_hal/time_total);
    } else { // hippo
      utils::logmesg(lmp,"  Repulse time: {:<12.6g} {:6.2f}%\n", time_repulse, time_repulse/time_total);
      utils::logmesg(lmp,"  Disp    time: {:<12.6g} {:6.2f}%\n", time_disp, time_disp/time_total);
    }
    utils::logmesg(lmp,"  Mpole   time: {:<12.6g} {:6.2f}%\n", time_mpole, time_mpole/time_total);
    utils::logmesg(lmp,"  Induce  time: {:<12.6g} {:6.2f}%\n", time_induce, time_induce/time_total);
    utils::logmesg(lmp,"  Polar   time: {:<12.6g} {:6.2f}%\n", time_polar, time_polar/time_total);
    if (!amoeba)
      utils::logmesg(lmp,"  Qxfer   time: {:.6g} {:.6g}\n", time_qxfer, time_qxfer/time_total);
    utils::logmesg(lmp,"  Total   time: {:.6g}\n",time_total * 100.0);

    #if DEBUG_AMOEBA
    double rspace_time = time_mpole_rspace + time_direct_rspace + time_mutual_rspace + time_polar_rspace;
    double kspace_time = time_mpole_kspace + time_direct_kspace + time_mutual_kspace + time_polar_kspace;

    utils::logmesg(lmp,"    Real-space timing breakdown: {:.3g}%\n", rspace_time/time_total);
    utils::logmesg(lmp,"      Mpole  time: {:.6g} {:.3g}%\n", time_mpole_rspace, time_mpole_rspace/time_total);
    utils::logmesg(lmp,"      Direct time: {:.6g} {:.3g}%\n", time_direct_rspace, time_direct_rspace/time_total);
    utils::logmesg(lmp,"      Mutual time: {:.6g} {:.3g}%\n", time_mutual_rspace, time_mutual_rspace/time_total);
    utils::logmesg(lmp,"      Polar  time: {:.6g} {:.3g}%\n", time_polar_rspace, time_polar_rspace/time_total);
    utils::logmesg(lmp,"    K-space timing breakdown   : {:.3g}%\n", kspace_time/time_total);
    utils::logmesg(lmp,"      Mpole  time: {:.6g} {:.3g}%\n", time_mpole_kspace, time_mpole_kspace/time_total);
    utils::logmesg(lmp,"      Direct time: {:.6g} {:.3g}%\n", time_direct_kspace, time_direct_kspace/time_total);
    utils::logmesg(lmp,"      Mutual time: {:.6g} {:.3g}%\n", time_mutual_kspace, time_mutual_kspace/time_total);
    utils::logmesg(lmp,"       - Grid    : {:.6g} {:.3g}%\n", time_grid_uind, time_grid_uind/time_total);
    utils::logmesg(lmp,"       - FFT     : {:.6g} {:.3g}%\n", time_mutual_fft, time_mutual_fft/time_total);
    utils::logmesg(lmp,"       - Interp  : {:.6g} {:.3g}%\n", time_fphi_uind, time_fphi_uind/time_total);
    utils::logmesg(lmp,"      Polar  time: {:.6g} {:.3g}%\n", time_polar_kspace, time_polar_kspace/time_total);
    #endif
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairAmoeba::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
   NOTE: these undocumented args are only for debugging
------------------------------------------------------------------------- */

void PairAmoeba::settings(int narg, char **arg)
{
  // turn on all FF components by default
  // first 4 lines are non-bonded terms
  // last 2 lines are bonded terms

  hal_flag = repulse_flag = qxfer_flag = 1;
  disp_rspace_flag = disp_kspace_flag = 1;
  polar_rspace_flag = polar_kspace_flag = 1;
  mpole_rspace_flag = mpole_kspace_flag = 1;
  bond_flag = angle_flag = dihedral_flag = improper_flag = 1;
  urey_flag = pitorsion_flag = bitorsion_flag = 1;

  int newvalue = -1;

  // include only specified FF components

  if (narg && (strcmp(arg[0],"include") == 0)) {
    newvalue = 1;
    hal_flag = repulse_flag = qxfer_flag = 0;
    disp_rspace_flag = disp_kspace_flag = 0;
    polar_rspace_flag = polar_kspace_flag = 0;
    mpole_rspace_flag = mpole_kspace_flag = 0;
    bond_flag = angle_flag = dihedral_flag = improper_flag = 0;
    urey_flag = pitorsion_flag = bitorsion_flag = 0;

  // exclude only specified FF components

  } else if (narg && (strcmp(arg[0],"exclude") == 0)) {
    newvalue = 0;

  } else if (narg) error->all(FLERR,"Illegal pair_style command");

  if (narg == 0) return;

  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // toggle components to include or exclude

  for (int iarg = 1; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"hal") == 0) hal_flag = newvalue;
    else if (strcmp(arg[iarg],"repulse") == 0) repulse_flag = newvalue;
    else if (strcmp(arg[iarg],"qxfer") == 0) qxfer_flag = newvalue;
    else if (strcmp(arg[iarg],"disp") == 0)
      disp_rspace_flag = disp_kspace_flag = newvalue;
    else if (strcmp(arg[iarg],"disp/rspace") == 0) disp_rspace_flag = newvalue;
    else if (strcmp(arg[iarg],"disp/kspace") == 0) disp_kspace_flag = newvalue;
    else if (strcmp(arg[iarg],"polar") == 0)
      polar_rspace_flag = polar_kspace_flag = newvalue;
    else if (strcmp(arg[iarg],"polar/rspace") == 0) polar_rspace_flag = newvalue;
    else if (strcmp(arg[iarg],"polar/kspace") == 0) polar_kspace_flag = newvalue;
    else if (strcmp(arg[iarg],"mpole") == 0)
      mpole_rspace_flag = mpole_kspace_flag = newvalue;
    else if (strcmp(arg[iarg],"mpole/rspace") == 0) mpole_rspace_flag = newvalue;
    else if (strcmp(arg[iarg],"mpole/kspace") == 0) mpole_kspace_flag = newvalue;
    else if (strcmp(arg[iarg],"bond") == 0) bond_flag = newvalue;
    else if (strcmp(arg[iarg],"angle") == 0) angle_flag = newvalue;
    else if (strcmp(arg[iarg],"dihedral") == 0) dihedral_flag = newvalue;
    else if (strcmp(arg[iarg],"improper") == 0) improper_flag = newvalue;
    else if (strcmp(arg[iarg],"urey") == 0) urey_flag = newvalue;
    else if (strcmp(arg[iarg],"pitorsion") == 0) pitorsion_flag = newvalue;
    else if (strcmp(arg[iarg],"bitorsion") == 0) bitorsion_flag = newvalue;
    else error->all(FLERR,"Illegal pair_style command");
  }

  // cannot disable bond and dihedral terms b/c those classes not in AMOEBA pkg

  if ((bond_flag == 0 || dihedral_flag == 0) && comm->me == 0)
    error->warning(FLERR,"Cannot disable AMOEBA bonds or dihedrals - "
                   "use bond_style or dihedral_style none instead");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAmoeba::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  if ((narg < 3) || (narg > 4)) error->all(FLERR,"Incorrect args for pair coefficients");

  // set setflag since coeff() is only called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 1;

  // read force field PRM file and optional KEY file

  set_defaults();
  read_prmfile(arg[2]);
  if (narg == 3) read_keyfile(nullptr);
  else read_keyfile(arg[3]);

  // compute Vdwl mixing rules, only for AMOEBA

  if (amoeba) {
    allocate_vdwl();
    mix();
  }

  // allocate arrays that depend on optorder or maxualt values from keyfile

  allocate_smallsize();

  // set copt and comp values, now that allocated to 0:optorder

  for (i = 0; i <= optorder; i++)
    copt[i] = copm[i] = 0.0;

  if (optorder == 1) {
    copt[0] = 0.530;
    copt[1] = 0.604;
  } else if (optorder == 2) {
    copt[0] = 0.042;
    copt[1] = 0.635;
    copt[2] = 0.414;
  } else if (optorder == 3) {
    copt[0] = -0.132;
    copt[1] = 0.218;
    copt[2] = 0.637;
    copt[3] = 0.293;
  } else if (optorder == 4) {
    copt[0] = -0.071;
    copt[1] = -0.096;
    copt[2] = 0.358;
    copt[3] = 0.587;
    copt[4] = 0.216;
  } else if (optorder == 5) {
    copt[0] = -0.005;
    copt[1] = -0.129;
    copt[2] = -0.026;
    copt[3] = 0.465;
    copt[4] = 0.528;
    copt[5] = 0.161;
  } else if (optorder == 6) {
    copt[0] = 0.014;
    copt[1] = -0.041;
    copt[2] = -0.172;
    copt[3] = 0.073;
    copt[4] = 0.535;
    copt[5] = 0.467;
    copt[6] = 0.122;
  }

  for (i = 0; i <= optorder; i++)
    for (j = optorder; j >= i; j--)
      copm[i] += copt[j];
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAmoeba::init_style()
{
  // error checks
  if (strcmp(update->unit_style,"real") != 0)
    error->all(FLERR, "Pair style {} requires real units", mystyle);
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style {} requires newton pair on", mystyle);
  if (domain->dimension == 2)
    error->all(FLERR, "Pair style {} requires a 3d system", mystyle);
  if (domain->triclinic)
    error->all(FLERR, "Pair style {} does not (yet) support triclinic systems", mystyle);

  int nperiodic = domain->xperiodic + domain->yperiodic + domain->zperiodic;
  if ((nperiodic != 0) && (nperiodic != 3))
    error->all(FLERR,"Pair style {} requires a fully periodic or fully non-periodic system",
               mystyle);

  if (!atom->q_flag)
    error->all(FLERR,"Pair style {} requires atom attribute q", mystyle);
  //if (!force->special_onefive)
  //  error->all(FLERR,"Pair style amoeba/hippo requires special_bonds one/five be set");

  // b/c polar uses mutipole virial terms

  if (apewald == aeewald && polar_kspace_flag && !mpole_kspace_flag)
    error->all(FLERR, "Pair {} with apewald = aeewald requires mpole and polar together", mystyle);

  // check if all custom atom arrays were set via fix property/atom

  int flag,cols;

  index_amtype = atom->find_custom("amtype",flag,cols);
  if (index_amtype < 0 || flag || cols)
    error->all(FLERR,"Pair {} amtype is not defined", mystyle);
  index_amgroup = atom->find_custom("amgroup",flag,cols);
  if (index_amgroup < 0 || flag || cols)
    error->all(FLERR,"Pair {} amgroup is not defined", mystyle);

  index_redID = atom->find_custom("redID",flag,cols);
  if (index_redID < 0 || !flag || cols)
    error->all(FLERR,"Pair {} redID is not defined", mystyle);
  index_xyzaxis = atom->find_custom("xyzaxis",flag,cols);
  if (index_xyzaxis < 0 || !flag || cols == 0)
    error->all(FLERR,"Pair {} xyzaxis is not defined", mystyle);

  index_polaxe = atom->find_custom("polaxe",flag,cols);
  if (index_polaxe < 0 || flag || cols)
    error->all(FLERR,"Pair {} polaxe is not defined", mystyle);
  index_pval = atom->find_custom("pval",flag,cols);
  if (index_pval < 0 || !flag || cols)
    error->all(FLERR,"Pair {} pval is not defined", mystyle);

  // -------------------------------------------------------------------
  // one-time initializations
  // can't do earlier b/c need all atoms to exist
  // -------------------------------------------------------------------

  // creation of per-atom storage
  // create a new fix STORE style for each atom's pole vector
  // id = "AMOEBA_pole", fix group = all

  // TODO: shouldn't there be an instance_me added to the identifier
  // in case there would be multiple pair style instances in a hybrid pair style?

  Fix *myfix;
  if (first_flag) {
    id_pole = utils::strdup("AMOEBA_pole");
    myfix = modify->add_fix(fmt::format("{} {} STORE/ATOM 13 0 0 1",id_pole,group->names[0]));
    fixpole = dynamic_cast<FixStoreAtom *>(myfix);
  }

  // creation of per-atom storage
  // create 2 new fix STORE styles for each atom's induced dipole history info
  // id = "AMOEBA_udalt", fix group = all
  // id = "AMOEBA_upalt", fix group = all
  // only if using preconditioner

  if (first_flag && use_pred) {
    id_udalt = utils::strdup("AMOEBA_udalt");
    myfix = modify->add_fix(fmt::format("{} {} STORE/ATOM {} 3 0 1",
                                        id_udalt, group->names[0], maxualt));
    fixudalt = dynamic_cast<FixStoreAtom *>(myfix);

    id_upalt = utils::strdup("AMOEBA_upalt");
    myfix = modify->add_fix(fmt::format("{} {} STORE/ATOM {} 3 0 1",
                                        id_upalt, group->names[0], maxualt));
    fixupalt = dynamic_cast<FixStoreAtom *>(myfix);
  }

  // create pages for storing pairwise data:
  // dipole/dipole interactions and preconditioner values

  if (first_flag) {
    ipage_dipole = new MyPage<int>();
    dpage_dipdip = new MyPage<double>();
    ipage_dipole->init(neighbor->oneatom,neighbor->pgsize);
    dpage_dipdip->init(6*neighbor->oneatom,6*neighbor->pgsize);

    if (poltyp == MUTUAL && pcgprec) {
      ipage_precond = new MyPage<int>();
      dpage_pcpc = new MyPage<double>();
      ipage_precond->init(neighbor->oneatom,neighbor->pgsize);
      dpage_pcpc->init(6*neighbor->oneatom,6*neighbor->pgsize);
    }
  }

  // initialize KSpace Ewald settings and FFTs and parallel grid objects
  // Coulombic grid is used with two orders: bseorder and bsporder
  //   so need two Grid3d instantiations for ghost comm

  if (first_flag) {
    kewald();
    if (use_ewald) {
      m_kspace =
        new AmoebaConvolution(lmp,this,nefft1,nefft2,nefft3,bseorder,MPOLE_GRID);
      p_kspace =
        new AmoebaConvolution(lmp,this,nefft1,nefft2,nefft3,bsporder,POLAR_GRID);
      pc_kspace =
        new AmoebaConvolution(lmp,this,nefft1,nefft2,nefft3,bsporder,POLAR_GRIDC);
      i_kspace =
        new AmoebaConvolution(lmp,this,nefft1,nefft2,nefft3,bsporder,INDUCE_GRID);
      ic_kspace =
        new AmoebaConvolution(lmp,this,nefft1,nefft2,nefft3,bsporder,INDUCE_GRIDC);

      // qfac is shared by induce and polar
      // gridfft1 is copy of FFT grid used within polar

      int nmine = p_kspace->nfft_owned;
      memory->create(qfac,nmine,"ameoba/induce:qfac");
      memory->create(gridfft1,2*nmine,"amoeba/polar:gridfft1");
    }
    if (use_dewald) {
      d_kspace =
        new AmoebaConvolution(lmp,this,ndfft1,ndfft2,ndfft3,bsdorder,DISP_GRID);
    }
  }

  // set csixpr = sum of csix[i]*csix[j] for a double loop over all atoms
  // compute this efficiently as M^2 instead of N^2, where M = # of classes
  // csix_num[iclass] = # of atoms in class Iclass

  if (first_flag) {
    amtype = atom->ivector[index_amtype];
    int nlocal = atom->nlocal;

    int *csix_num_one = new int[n_amclass+1];
    for (int i = 0; i <= n_amclass; i++) csix_num_one[i] = 0;

    int itype,iclass;

    for (int i = 0; i < nlocal; i++) {
      itype = amtype[i];
      iclass = amtype2class[itype];
      csix_num_one[iclass]++;
    }

    int *csix_num = new int[n_amclass+1];
    MPI_Allreduce(csix_num_one,csix_num,n_amclass+1,MPI_INT,MPI_SUM,world);

    csixpr = 0.0;
    for (int i = 1; i <= n_amclass; i++) {
      for (int j = i+1; j <= n_amclass; j++) {
        csixpr += csix[i]*csix[j] * csix_num[i]*csix_num[j];
      }
    }
    csixpr *= 2.0;

    for (int i = 1; i <= n_amclass; i++)
      csixpr += csix[i]*csix[i] * csix_num[i]*csix_num[i];

    delete[] csix_num_one;
    delete[] csix_num;
  }

  // initialize peratom pval to zero
  // so that initial ghost comm will be valid
  // pval is not set until first call to compute(), and only for HIPPO

  if (first_flag) {
    double *pval = atom->dvector[index_pval];
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) pval[i] = 0.0;
  }

  // output FF settings to screen and logfile

  if (first_flag && (comm->me == 0)) print_settings();

  // all done with one-time initializations

  first_flag = 0;

  // -------------------------------------------------------------------
  // end of one-time initializations
  // -------------------------------------------------------------------

  // check for fixes which store persistent per-atom properties

  if (id_pole) {
    myfix = modify->get_fix_by_id(id_pole);
    if (!myfix)
      error->all(FLERR,"Could not find internal pair amoeba fix STORE/ATOM id {}", id_pole);
    fixpole = dynamic_cast<FixStoreAtom *>(myfix);

  }

  if (id_udalt) {
    myfix = modify->get_fix_by_id(id_udalt);
    if (!myfix)
      error->all(FLERR,"Could not find internal pair amoeba fix STORE/ATOM id {}", id_udalt);
    fixudalt = dynamic_cast<FixStoreAtom *>(myfix);

    myfix = modify->get_fix_by_id(id_upalt);
    if (!myfix)
      error->all(FLERR,"Could not find internal pair amoeba fix STORE/ATOM id {}", id_upalt);
    fixupalt = dynamic_cast<FixStoreAtom *>(myfix);
  }

  // assign hydrogen neighbors (redID) to each owned atom
  // only set if kred[i] is non-zero and I is bonded to a single atom
  // conceptually: non-zero if I is hydrogen bonded to another atom

  if (amoeba) {
    amtype = atom->ivector[index_amtype];
    double *redID = atom->dvector[index_redID];

    int **nspecial = atom->nspecial;
    tagint **special = atom->special;
    int nlocal = atom->nlocal;

    int itype,iclass;

    for (int i = 0; i < nlocal; i++) {
      itype = amtype[i];
      iclass = amtype2class[itype];
      if (kred[iclass] == 0.0) {
        redID[i] = 0.0;
      }
      else if (nspecial[i][0] != 1) {
        redID[i] = 0.0;
      }
      else {
        redID[i] = ubuf(special[i][0]).d;
      }
    }
  }

  // set KSpace recip matrix based on box size and shape

  lattice();

  // can now set comm size needed by this Pair
  // cfstyle KMPOLE is max # of 1-2 bond partners, smaller than comm_forward

  if (amoeba) comm_forward = 16; // xred, rpole
  else comm_forward = 13;        // just rpole
  int fsize = 6;
  if (poltyp == OPT) fsize += 6*optorder;
  //if (poltyp == TCG) fsize += 12*tcgnab;
  comm_forward = MAX(comm_forward,fsize);

  comm_reverse = 9;

  // request standard neighbor list

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   print settings to screen and logfile
------------------------------------------------------------------------- */

void PairAmoeba::print_settings()
{
  std::string mesg = utils::uppercase(mystyle) + " force field settings\n";

  if (amoeba) {
    choose(HAL);
    mesg += fmt::format("  hal: cut {} taper {} vscale {} {} {} {}\n", sqrt(off2),sqrt(cut2),
                        special_hal[1],special_hal[2],special_hal[3],special_hal[4]);
  } else {
    choose(REPULSE);
    mesg += fmt::format("  repulsion: cut {} taper {} rscale {} {} {} {}\n", sqrt(off2),sqrt(cut2),
                        special_repel[1],special_repel[2],special_repel[3],special_repel[4]);

    choose(QFER);
    mesg += fmt::format("  qxfer: cut {} taper {} mscale {} {} {} {}\n", sqrt(off2),sqrt(cut2),
                        special_mpole[1],special_mpole[2],special_mpole[3],special_mpole[4]);

    if (use_dewald) {
      choose(DISP_LONG);
      mesg += fmt::format("  dispersion: cut {} aewald {} bsorder {} FFT {} {} {} "
                          "dspscale {} {} {} {}\n", sqrt(off2),aewald,bsdorder,ndfft1,ndfft2,ndfft3,
                          special_disp[1],special_disp[2],special_disp[3],special_disp[4]);
    } else {
      choose(DISP);
      mesg += fmt::format("  dispersion: cut {} aewald {} dspscale {} {} {} {}\n",
                          sqrt(off2),aewald,special_disp[1],
                          special_disp[2],special_disp[3],special_disp[4]);
    }
  }

  if (use_ewald) {
    choose(MPOLE_LONG);
    mesg += fmt::format("  multipole: cut {} aewald {} bsorder {} FFT {} {} {} "
                        "mscale {} {} {} {}\n",
                        sqrt(off2),aewald,bseorder,nefft1,nefft2,nefft3,
                        special_mpole[1],special_mpole[2],special_mpole[3],special_mpole[4]);
  } else {
    choose(MPOLE);
    mesg += fmt::format("  multipole: cut {} aewald {} mscale {} {} {} {}\n", sqrt(off2),aewald,
                        special_mpole[1],special_mpole[2],special_mpole[3],special_mpole[4]);
  }

  if (use_ewald) {
    choose(POLAR_LONG);
    mesg += fmt::format("  polar: cut {} aewald {} bsorder {} FFT {} {} {}\n",
                        sqrt(off2),aewald,bsporder,nefft1,nefft2,nefft3);
    mesg += fmt::format("         pscale {} {} {} {} piscale {} {} {} {} "
                        "wscale {} {} {} {} d/u scale {} {}\n",
                        special_polar_pscale[1],special_polar_pscale[2],
                        special_polar_pscale[3],special_polar_pscale[4],
                        special_polar_piscale[1],special_polar_piscale[2],
                        special_polar_piscale[3],special_polar_piscale[4],
                        special_polar_wscale[1],special_polar_wscale[2],
                        special_polar_wscale[3],special_polar_wscale[4],
                        polar_dscale,polar_uscale);
  } else {
    choose(POLAR);
    mesg += fmt::format("  polar: cut {} aewald {}\n",sqrt(off2),aewald);
    mesg += fmt::format("         pscale {} {} {} {} piscale {} {} {} {} "
                        "wscale {} {} {} {} d/u scale {} {}\n",
                        special_polar_pscale[1],special_polar_pscale[2],
                        special_polar_pscale[3],special_polar_pscale[4],
                        special_polar_piscale[1],special_polar_piscale[2],
                        special_polar_piscale[3],special_polar_piscale[4],
                        special_polar_wscale[1],special_polar_wscale[2],
                        special_polar_wscale[3],special_polar_wscale[4],
                        polar_dscale,polar_uscale);
  }

  choose(USOLV);
  mesg += fmt::format("  precondition: cut {}\n",sqrt(off2));
  utils::logmesg(lmp, mesg);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAmoeba::init_one(int /*i*/, int /*j*/)
{
  double cutoff = 0.0;

  if (amoeba) {
    choose(HAL);
    cutoff = MAX(cutoff,sqrt(off2));
  } else {
    choose(REPULSE);
    cutoff = MAX(cutoff,sqrt(off2));
    if (use_dewald) choose(DISP_LONG);
    else choose(DISP);
    cutoff = MAX(cutoff,sqrt(off2));
  }

  if (use_ewald) choose(MPOLE_LONG);
  else choose(MPOLE);
  cutoff = MAX(cutoff,sqrt(off2));

  if (use_ewald) choose(POLAR_LONG);
  else choose(POLAR);
  cutoff = MAX(cutoff,sqrt(off2));

  if (!amoeba) {
    choose(QFER);
    cutoff = MAX(cutoff,sqrt(off2));
  }

  return cutoff;
}

/* ---------------------------------------------------------------------- */

int PairAmoeba::pack_forward_comm(int n, int *list, double *buf,
                                  int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m;

  m = 0;

  if (cfstyle == INDUCE) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = uind[j][0];
      buf[m++] = uind[j][1];
      buf[m++] = uind[j][2];
      buf[m++] = uinp[j][0];
      buf[m++] = uinp[j][1];
      buf[m++] = uinp[j][2];
    }

    if (poltyp == OPT) {
      for (i = 0; i < n; i++) {
        j = list[i];
        for (k = 0; k < optorder; k++) {
          buf[m++] = uopt[j][k][0];
          buf[m++] = uopt[j][k][1];
          buf[m++] = uopt[j][k][2];
          buf[m++] = uoptp[j][k][0];
          buf[m++] = uoptp[j][k][1];
          buf[m++] = uoptp[j][k][2];
        }
      }
    }

    /*
    if (poltyp == TCG) {
      for (i = 0; i < n; i++) {
        j = list[i];
        for (k = 0; k < tcgnab; k++) {
          buf[m++] = uad[k][j][0];
          buf[m++] = uad[k][j][1];
          buf[m++] = uad[k][j][2];
          buf[m++] = uap[k][j][0];
          buf[m++] = uap[k][j][1];
          buf[m++] = uap[k][j][2];
          buf[m++] = ubd[k][j][0];
          buf[m++] = ubd[k][j][1];
          buf[m++] = ubd[k][j][2];
          buf[m++] = ubp[k][j][0];
          buf[m++] = ubp[k][j][1];
          buf[m++] = ubp[k][j][2];
        }
      }
    }
    */

  } else if (cfstyle == RSD) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = rsd[j][0];
      buf[m++] = rsd[j][1];
      buf[m++] = rsd[j][2];
      buf[m++] = rsdp[j][0];
      buf[m++] = rsdp[j][1];
      buf[m++] = rsdp[j][2];
    }

  } else if (cfstyle == SETUP_AMOEBA) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = xred[j][0];
      buf[m++] = xred[j][1];
      buf[m++] = xred[j][2];
      for (k = 0; k < 13; k++)
        buf[m++] = rpole[j][k];
    }

  } else if (cfstyle == SETUP_HIPPO) {
    for (i = 0; i < n; i++) {
      j = list[i];
      for (k = 0; k < 13; k++)
        buf[m++] = rpole[j][k];
    }

  } else if (cfstyle == KMPOLE) {
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(nspecial[j][0]).d;
      for (k = 0; k < nspecial[j][0]; k++)
        buf[m++] = ubuf(special[j][k]).d;
    }

  } else if (cfstyle == AMGROUP) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(amgroup[j]).d;
    }

  } else if (cfstyle == PVAL) {
    double *pval = atom->dvector[index_pval];
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = pval[j];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::unpack_forward_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;

  if (cfstyle == INDUCE) {
    for (i = first; i < last; i++) {
      uind[i][0] = buf[m++];
      uind[i][1] = buf[m++];
      uind[i][2] = buf[m++];
      uinp[i][0] = buf[m++];
      uinp[i][1] = buf[m++];
      uinp[i][2] = buf[m++];
    }

    if (poltyp == OPT) {
      for (i = first; i < last; i++) {
        for (k = 0; k < optorder; k++) {
          uopt[i][k][0] = buf[m++];
          uopt[i][k][1] = buf[m++];
          uopt[i][k][2] = buf[m++];
          uoptp[i][k][0] = buf[m++];
          uoptp[i][k][1] = buf[m++];
          uoptp[i][k][2] = buf[m++];
        }
      }
    }

    /*
    if (poltyp == TCG) {
      for (i = first; i < last; i++) {
        for (k = 0; k < tcgnab; k++) {
          uad[k][i][0] = buf[m++];
          uad[k][i][1] = buf[m++];
          uad[k][i][2] = buf[m++];
          uap[k][i][0] = buf[m++];
          uap[k][i][1] = buf[m++];
          uap[k][i][2] = buf[m++];
          ubd[k][i][0] = buf[m++];
          ubd[k][i][1] = buf[m++];
          ubd[k][i][2] = buf[m++];
          ubp[k][i][0] = buf[m++];
          ubp[k][i][1] = buf[m++];
          ubp[k][i][2] = buf[m++];
        }
      }
    }
    */

  } else if (cfstyle == RSD) {
    for (i = first; i < last; i++) {
      rsd[i][0] = buf[m++];
      rsd[i][1] = buf[m++];
      rsd[i][2] = buf[m++];
      rsdp[i][0] = buf[m++];
      rsdp[i][1] = buf[m++];
      rsdp[i][2] = buf[m++];
    }

  } else if (cfstyle == SETUP_AMOEBA) {
    for (i = first; i < last; i++) {
      xred[i][0] = buf[m++];
      xred[i][1] = buf[m++];
      xred[i][2] = buf[m++];
      for (k = 0; k < 13; k++)
        rpole[i][k] = buf[m++];
    }

  } else if (cfstyle == SETUP_HIPPO) {
    for (i = first; i < last; i++) {
      for (k = 0; k < 13; k++)
        rpole[i][k] = buf[m++];
    }


  } else if (cfstyle == KMPOLE) {
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;
    for (i = first; i < last; i++) {
      nspecial[i][0] = (int) ubuf(buf[m++]).i;
      for (k = 0; k < nspecial[i][0]; k++)
        special[i][k] = (tagint) ubuf(buf[m++]).i;
    }

  } else if (cfstyle == AMGROUP) {
    for (i = first; i < last; i++) {
      amgroup[i] = (int) ubuf(buf[m++]).i;
    }

  } else if (cfstyle == PVAL) {
    double *pval = atom->dvector[index_pval];
    for (i = first; i < last; i++) {
      pval[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairAmoeba::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (crstyle == FIELD) {
    for (i = first; i < last; i++) {
      buf[m++] = field[i][0];
      buf[m++] = field[i][1];
      buf[m++] = field[i][2];
      buf[m++] = fieldp[i][0];
      buf[m++] = fieldp[i][1];
      buf[m++] = fieldp[i][2];
    }
  } else if (crstyle == ZRSD) {
    for (i = first; i < last; i++) {
      buf[m++] = zrsd[i][0];
      buf[m++] = zrsd[i][1];
      buf[m++] = zrsd[i][2];
      buf[m++] = zrsdp[i][0];
      buf[m++] = zrsdp[i][1];
      buf[m++] = zrsdp[i][2];
    }
  } else if (crstyle == TORQUE) {
    for (i = first; i < last; i++) {
      buf[m++] = tq[i][0];
      buf[m++] = tq[i][1];
      buf[m++] = tq[i][2];
    }
  } else if (crstyle == UFLD) {
    for (i = first; i < last; i++) {
      buf[m++] = ufld[i][0];
      buf[m++] = ufld[i][1];
      buf[m++] = ufld[i][2];
      buf[m++] = dufld[i][0];
      buf[m++] = dufld[i][1];
      buf[m++] = dufld[i][2];
      buf[m++] = dufld[i][3];
      buf[m++] = dufld[i][4];
      buf[m++] = dufld[i][5];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  if (crstyle == FIELD) {
    for (i = 0; i < n; i++) {
      j = list[i];
      field[j][0] += buf[m++];
      field[j][1] += buf[m++];
      field[j][2] += buf[m++];
      fieldp[j][0] += buf[m++];
      fieldp[j][1] += buf[m++];
      fieldp[j][2] += buf[m++];
    }
  } else if (crstyle == ZRSD) {
    for (i = 0; i < n; i++) {
      j = list[i];
      zrsd[j][0] += buf[m++];
      zrsd[j][1] += buf[m++];
      zrsd[j][2] += buf[m++];
      zrsdp[j][0] += buf[m++];
      zrsdp[j][1] += buf[m++];
      zrsdp[j][2] += buf[m++];
    }
  } else if (crstyle == TORQUE) {
    for (i = 0; i < n; i++) {
      j = list[i];
      tq[j][0] += buf[m++];
      tq[j][1] += buf[m++];
      tq[j][2] += buf[m++];
    }
  } else if (crstyle == UFLD) {
    for (i = 0; i < n; i++) {
      j = list[i];
      ufld[j][0] += buf[m++];
      ufld[j][1] += buf[m++];
      ufld[j][2] += buf[m++];
      dufld[j][0] += buf[m++];
      dufld[j][1] += buf[m++];
      dufld[j][2] += buf[m++];
      dufld[j][3] += buf[m++];
      dufld[j][4] += buf[m++];
      dufld[j][5] += buf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   subset of FFT grids assigned to each proc may have changed
   notify each instance of AmoebaConvolution class
   called by load balancer when proc subdomains are adjusted
------------------------------------------------------------------------- */

void PairAmoeba::reset_grid()
{
  if (use_ewald) {
    m_kspace->reset_grid();
    p_kspace->reset_grid();
    pc_kspace->reset_grid();
    i_kspace->reset_grid();
    ic_kspace->reset_grid();
  }
  if (use_dewald) d_kspace->reset_grid();

  // qfac is shared by induce and polar
  // gridfft1 is copy of FFT grid used within polar

  memory->destroy(qfac);
  memory->destroy(gridfft1);
  int nmine = p_kspace->nfft_owned;
  memory->create(qfac,nmine,"ameoba/induce:qfac");
  memory->create(gridfft1,2*nmine,"amoeba/polar:gridfft1");
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void PairAmoeba::pack_forward_grid(int which, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  if (which == MPOLE_GRID) {
    FFT_SCALAR *src = m_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (which == POLAR_GRID) {
    FFT_SCALAR *src = p_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (which == POLAR_GRIDC) {
    FFT_SCALAR *src = pc_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src[2*list[i]];
      buf[n++] = src[2*list[i]+1];
    }
  } else if (which == DISP_GRID) {
    FFT_SCALAR *src = d_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (which == INDUCE_GRID) {
    FFT_SCALAR *src = i_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (which == INDUCE_GRIDC) {
    FFT_SCALAR *src = ic_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src[2*list[i]];
      buf[n++] = src[2*list[i]+1];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void PairAmoeba::unpack_forward_grid(int which, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  if (which == MPOLE_GRID) {
    FFT_SCALAR *dest = m_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (which == POLAR_GRID) {
    FFT_SCALAR *dest = p_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (which == POLAR_GRIDC) {
    FFT_SCALAR *dest = pc_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      dest[2*list[i]] = buf[n++];
      dest[2*list[i]+1] = buf[n++];
    }
  } else if (which == DISP_GRID) {
    FFT_SCALAR *dest = d_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (which == INDUCE_GRID) {
    FFT_SCALAR *dest = i_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (which == INDUCE_GRIDC) {
    FFT_SCALAR *dest = ic_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      dest[2*list[i]] = buf[n++];
      dest[2*list[i]+1] = buf[n++];
    }
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void PairAmoeba::pack_reverse_grid(int which, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  if (which == MPOLE_GRID) {
    FFT_SCALAR *src = m_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (which == POLAR_GRID) {
    FFT_SCALAR *src = p_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (which == POLAR_GRIDC) {
    FFT_SCALAR *src = pc_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src[2*list[i]];
      buf[n++] = src[2*list[i]+1];
    }
  } else if (which == DISP_GRID) {
    FFT_SCALAR *src = d_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (which == INDUCE_GRID) {
    FFT_SCALAR *src = i_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (which == INDUCE_GRIDC) {
    FFT_SCALAR *src = ic_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src[2*list[i]];
      buf[n++] = src[2*list[i]+1];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void PairAmoeba::unpack_reverse_grid(int which, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  if (which == MPOLE_GRID) {
    FFT_SCALAR *dest = m_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (which == POLAR_GRID) {
    FFT_SCALAR *dest = p_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (which == POLAR_GRIDC) {
    FFT_SCALAR *dest = pc_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      dest[2*list[i]] += buf[n++];
      dest[2*list[i]+1] += buf[n++];
    }
  } else if (which == DISP_GRID) {
    FFT_SCALAR *dest = d_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (which == INDUCE_GRID) {
    FFT_SCALAR *dest = i_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (which == INDUCE_GRIDC) {
    FFT_SCALAR *dest = ic_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      dest[2*list[i]] += buf[n++];
      dest[2*list[i]+1] += buf[n++];
    }
  }
}

// ----------------------------------------------------------------------
// AMOEBA/HIPPO specific methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   assign atoms to polarization groups with unique IDs
------------------------------------------------------------------------- */

void PairAmoeba::assign_groups()
{
  int i,j,m,jtype,mtype;
  int nbond,ngroup,ibond,igroup,ghostmark,anyghostmark;
  tagint jglobal;

  int nstack = 0;
  int maxstack = 0;
  int *stack = nullptr;

  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // initially, groupID = atomID
  // communicate new groupIDs to ghost atoms

  for (i = 0; i < nlocal; i++) amgroup[i] = tag[i];
  cfstyle = AMGROUP;
  comm->forward_comm(this);

  // loop until no ghost atom groupIDs are reset

  while (true) {

    // loop over all atoms and their group neighborhoods

    ghostmark = 0;
    for (i = 0; i < nlocal; i++) {

      // push atom I on stack

      nstack = 0;
      if (nstack == maxstack) {
        maxstack += DELTASTACK;
        memory->grow(stack,maxstack,"amoeba:stack");
      }
      stack[nstack++] = i;

      // loop over I's group neighborhood until stack is empty

      while (nstack > 0) {

        // pop atom M off stack

        m = stack[nstack-1];
        nstack--;
        mtype = amtype[m];

        // loop over bond partners of atom M

        nbond = nspecial[m][0];
        for (ibond = 0; ibond < nbond; ibond++) {
          jglobal = special[m][ibond];
          j = atom->map(jglobal);
          if (j < 0)
            error->one(FLERR,"AMOEBA group assignment bond neighbor not found");
          jtype = amtype[j];

          // if amtype of bondpartner J is not in polgroup of atom M, continue

          ngroup = npolgroup[mtype];
          for (igroup = 0; igroup < ngroup; igroup++)
            if (jtype == polgroup[mtype][igroup]) break;
          if (igroup == ngroup) continue;

          // if groupID of atoms J and M are the same, continue
          // else set atom with larger groupID to smaller groupID
          // if changed atom is ghost, set ghostmark, else push atom on stack

          if (amgroup[m] == amgroup[j]) continue;

          if (amgroup[m] > amgroup[j]) {
            amgroup[m] = amgroup[j];
            if (nstack == maxstack) {
              maxstack += DELTASTACK;
              memory->grow(stack,maxstack,"amoeba:stack");
            }
            stack[nstack++] = m;
          } else {
            amgroup[j] = amgroup[m];
            if (j >= nlocal) ghostmark = 1;
            else {
              if (nstack == maxstack) {
                maxstack += DELTASTACK;
                memory->grow(stack,maxstack,"amoeba:stack");
              }
              stack[nstack++] = j;
            }
          }
        }
      }
    }

    // communicate new groupIDs to ghost atoms

    cfstyle = AMGROUP;
    comm->forward_comm(this);

    // done if no proc reset groupID of a ghost atom

    MPI_Allreduce(&ghostmark,&anyghostmark,1,MPI_INT,MPI_MAX,world);
    if (!anyghostmark) break;
  }

  memory->destroy(stack);

  // print group count

  int count = 0;
  for (i = 0; i < nlocal; i++)
    if (tag[i] == amgroup[i]) count++;
  bigint bcount = count;
  bigint allbcount;
  MPI_Allreduce(&bcount,&allbcount,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (comm->me == 0)
    utils::logmesg(lmp, "  {} group count: {}\n",utils::uppercase(mystyle), allbcount);
}

/* ----------------------------------------------------------------------
   adjust xred for ghost atoms due to PBC
------------------------------------------------------------------------- */

void PairAmoeba::pbc_xred()
{
  double prd,prd_half,delta;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (domain->xperiodic) {
    prd = domain->xprd;
    prd_half = domain->xprd_half;
    for (int i = nlocal; i < nall; i++) {
      delta = xred[i][0] - x[i][0];
      while (fabs(delta) > prd_half) {
        if (delta < 0.0) xred[i][0] += prd;
        else xred[i][0] -= prd;
        delta = xred[i][0] - x[i][0];
      }
    }
  }

  if (domain->yperiodic) {
    prd = domain->yprd;
    prd_half = domain->yprd_half;
    for (int i = nlocal; i < nall; i++) {
      delta = xred[i][1] - x[i][1];
      while (fabs(delta) > prd_half) {
        if (delta < 0.0) xred[i][1] += prd;
        else xred[i][1] -= prd;
        delta = xred[i][1] - x[i][1];
      }
    }
  }

  if (domain->zperiodic) {
    prd = domain->zprd;
    prd_half = domain->zprd_half;
    for (int i = nlocal; i < nall; i++) {
      delta = xred[i][2] - x[i][2];
      while (fabs(delta) > prd_half) {
        if (delta < 0.0) xred[i][2] += prd;
        else xred[i][2] -= prd;
        delta = xred[i][2] - x[i][2];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   build reduced-size preconditioner neigh list from master neigh list
------------------------------------------------------------------------- */

void PairAmoeba::precond_neigh()
{
  int i,j,ii,jj,n,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *neighptr;

  // NOTE: no skin added to cutoff for this shorter neighbor list
  // also note that Tinker (and thus LAMMPS) does not apply the
  //   distance cutoff in the CG iterations in induce.cpp,
  //   rather all interactions in the precond neigh list are
  //   used every step until the neighbor list is rebuilt,
  //   this means the cutoff distance is not exactly enforced,
  //   on later steps atoms outside may contribute, atoms inside may not

  choose(USOLV);

  // atoms and neighbor list

  double **x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all induce neighs of owned atoms within shorter cutoff
  // scan longer-cutoff neighbor list of I

  ipage_precond->reset();

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    n = 0;
    neighptr = ipage_precond->vget();

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK15;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < off2) neighptr[n++] = jlist[jj];
    }

    firstneigh_precond[i] = neighptr;
    numneigh_precond[i] = n;
    ipage_precond->vgot(n);
  }
}


/* ----------------------------------------------------------------------
   allocate Vdwl arrays
   note that n_amclass = # of classes in Tinker PRM file
   actual number of classes for atoms in simulation may be smaller
   this is determined by the AMOEBA types listed in Tinker xyz file
     and their mapping to AMOEBA classes
------------------------------------------------------------------------- */

void PairAmoeba::initialize_vdwl()
{
  radmin = radmin4 = epsilon = epsilon4 = nullptr;
}

void PairAmoeba::allocate_vdwl()
{
  memory->create(radmin,n_amclass+1,n_amclass+1,"amoeba:radmin");
  memory->create(radmin4,n_amclass+1,n_amclass+1,"amoeba:radmin4");
  memory->create(epsilon,n_amclass+1,n_amclass+1,"amoeba:epsilon");
  memory->create(epsilon4,n_amclass+1,n_amclass+1,"amoeba:epsilon4");
}

void PairAmoeba::deallocate_vdwl()
{
  memory->destroy(radmin);
  memory->destroy(radmin4);
  memory->destroy(epsilon);
  memory->destroy(epsilon4);
}

/* ----------------------------------------------------------------------
   allocate small-size arrays
------------------------------------------------------------------------- */

void PairAmoeba::initialize_smallsize()
{
  copt = copm = nullptr;
  a_ualt = ap_ualt = nullptr;
  b_ualt = bp_ualt = nullptr;
  c_ualt = cp_ualt = nullptr;
  bpred = bpredp = bpreds = bpredps = nullptr;
  gear = aspc = nullptr;
}

void PairAmoeba::allocate_smallsize()
{
  // note use of optorder+1

  copt = new double[optorder+1];
  copm = new double[optorder+1];

  a_ualt = new double[maxualt*(maxualt+1)/2];
  ap_ualt = new double[maxualt*(maxualt+1)/2];
  b_ualt = new double[maxualt];
  bp_ualt = new double[maxualt];
  memory->create(c_ualt,maxualt,maxualt,"amoeba:c_ualt");
  memory->create(cp_ualt,maxualt,maxualt,"amoeba:cp_ualt");
  bpred = new double[maxualt];
  bpredp = new double[maxualt];
  bpreds = new double[maxualt];
  bpredps = new double[maxualt];
  if (use_pred) {
    if (polpred == GEAR) gear = new double[maxualt];
    if (polpred == ASPC) aspc = new double[maxualt];
  }
}

void PairAmoeba::deallocate_smallsize()
{
  delete[] copt;
  delete[] copm;
  delete[] a_ualt;
  delete[] ap_ualt;
  delete[] b_ualt;
  delete[] bp_ualt;
  memory->destroy(c_ualt);
  memory->destroy(cp_ualt);
  delete[] bpred;
  delete[] bpredp;
  delete[] bpreds;
  delete[] bpredps;
  delete[] gear;
  delete[] aspc;
}

/* ----------------------------------------------------------------------
   set cutoffs, taper constants, PME params for a FF component
------------------------------------------------------------------------- */

void PairAmoeba::choose(int which)
{
  double off = 0.0;
  double cut = 0.0;

  // short-range only terms

  if (which == HAL) {
    off = vdwcut;
    cut = vdwtaper;
  } else if (which == REPULSE) {
    off = repcut;
    cut = reptaper;
  } else if (which == QFER) {
    off = ctrncut;
    cut = ctrntaper;
  } else if (which == DISP) {
    off = dispcut;
    cut = disptaper;
    aewald = 0.0;
  } else if (which == MPOLE) {
    off = mpolecut;
    cut = mpoletaper;
    aewald = 0.0;
  } else if (which == POLAR) {
    off = ewaldcut;
    cut = 0.99*off;   // not used
    aewald = 0.0;
  } else if (which == USOLV) {
    off = usolvcut;
    cut = 0.99*off;   // not used

  // short-range + long-range terms

  } else if (which == DISP_LONG) {
    off = dispcut;
    cut = 0.99*cut;   // not used
    aewald = adewald;
  } else if (which == MPOLE_LONG) {
    off = mpolecut;
    cut = 0.99*cut;   // not used
    aewald = aeewald;
  } else if (which == POLAR_LONG) {
    off = ewaldcut;
    cut = 0.99*off;   // not used
    aewald = apewald;
  }

  off2 = off*off;
  cut2 = cut*cut;

  // taper coeffs

  double denom = powint(off-cut,5);
  c0 = off*off2 * (off2 - 5.0*off*cut + 10.0*cut2) / denom;
  c1 = -30.0 * off2*cut2 / denom;
  c2 = 30.0 * (off2*cut+off*cut2) / denom;
  c3 = -10.0 * (off2 + 4.0*off*cut + cut2) / denom;
  c4 = 15.0 * (off+cut) / denom;
  c5 = -6.0 / denom;
}

/* ----------------------------------------------------------------------
   compute mixing rules for all pairwise params on a per-class basis
   override default mixing with VDWLPR entries in force field file
   no vdwl14 terms are used by AMOEBA or HIPPO force fields
------------------------------------------------------------------------- */

void PairAmoeba::mix()
{
  int i,j,m;
  double ei,ej,sei,sej,eij;
  double ri,rj,sri,srj,rij;

  double TWOSIX = pow(2.0,1.0/6.0);

  for (i = 1; i <= n_amclass; i++) {
    for (j = i; j <= n_amclass; j++) {

      ei = vdwl_eps[i];
      ej = vdwl_eps[j];
      ri = vdwl_sigma[i];
      rj = vdwl_sigma[j];

      if (radius_type == SIGMA) {
        ri *= TWOSIX;
        rj *= TWOSIX;
      }
      if (radius_size == DIAMETER) {
        ri *= 0.5;
        rj *= 0.5;
      }

      sri = sqrt(ri);
      ei = fabs(ei);
      sei = sqrt(ei);
      srj = sqrt(rj);
      ej = fabs(ej);
      sej = sqrt(ej);

      if (ri == 0.0 && rj == 0.0) {
        rij = 0.0;
      } else if (radius_rule == ARITHMETIC) {
        rij = ri + rj;
      } else if (radius_rule == GEOMETRIC) {
        rij = 2.0 * sri * srj;
      } else if (radius_rule == CUBIC_MEAN) {
        rij = 2.0 * (ri*ri*ri + rj*rj*rj) / (ri*ri + rj*rj);
      } else {
        rij = ri + rj;
      }

      if (ei == 0.0 && ej == 0.0) {
        eij = 0.0;
      } else if (epsilon_rule == ARITHMETIC) {
        eij = 0.5 * (ei + ej);
      } else if (epsilon_rule == GEOMETRIC) {
        eij = sei * sej;
      } else if (epsilon_rule == HARMONIC) {
        eij = 2.0 * (ei*ej) / (ei+ej);
      } else if (epsilon_rule == HHG) {
        eij = 4.0 * (ei*ej) / ((sei+sej)*(sei+sej));
      } else if (epsilon_rule == W_H) {
        eij = 2.0 * (sei*sej) * powint(ri*rj,3) / (powint(ri,6) + powint(rj,6));
      } else {
        eij = sei * sej;
      }

      radmin[j][i] = radmin[i][j] = rij;
      radmin4[j][i] = radmin4[i][j] = rij;
      epsilon[j][i] = epsilon[i][j] = eij;
      epsilon4[j][i] = epsilon4[i][j] = eij;
    }
  }

  // override with VDWPR pairwise entries from force field file

  for (m = 0; m < nvdwl_pair; m++) {
    i = vdwl_class_pair[m][0];
    j = vdwl_class_pair[m][1];
    rij = vdwl_sigma_pair[m];
    eij = vdwl_eps_pair[m];

    if (radius_type == SIGMA) rij *= TWOSIX;

    radmin[j][i] = radmin[i][j] = rij;
    radmin4[j][i] = radmin4[i][j] = rij;
    epsilon[j][i] = epsilon[i][j] = eij;
    epsilon4[j][i] = epsilon4[i][j] = eij;
  }
}

/* ---------------------------------------------------------------------- */

void *PairAmoeba::extract(const char *str, int &dim)
{
  dim = 0;

  if (strcmp(str,"amtype") == 0) return (void *) amtype;
  if (strcmp(str,"atomic_num") == 0) return (void *) atomic_num;

  if (strcmp(str,"bond_flag") == 0) return (void *) &bond_flag;
  if (strcmp(str,"angle_flag") == 0) return (void *) &angle_flag;
  if (strcmp(str,"dihedral_flag") == 0) return (void *) &dihedral_flag;
  if (strcmp(str,"improper_flag") == 0) return (void *) &improper_flag;
  if (strcmp(str,"urey_flag") == 0) return (void *) &urey_flag;
  if (strcmp(str,"pitorsion_flag") == 0) return (void *) &pitorsion_flag;
  if (strcmp(str,"bitorsion_flag") == 0) return (void *) &bitorsion_flag;

  if (strcmp(str,"opbend_cubic") == 0) return (void *) &opbend_cubic;
  if (strcmp(str,"opbend_quartic") == 0) return (void *) &opbend_quartic;
  if (strcmp(str,"opbend_pentic") == 0) return (void *) &opbend_pentic;
  if (strcmp(str,"opbend_sextic") == 0) return (void *) &opbend_sextic;

  return nullptr;
}

/* ----------------------------------------------------------------------
   peratom requests from FixPair
   return ptr to requested data
   also return ncol = # of quantites per atom
     0 = per-atom vector
     1 or more = # of columns in per-atom array
   return NULL if str is not recognized
---------------------------------------------------------------------- */

void *PairAmoeba::extract_peratom(const char *str, int &ncol)
{
  if (strcmp(str,"uind") == 0) {
    ncol = 3;
    return (void *) uind;
  } else if (strcmp(str,"uinp") == 0) {
    ncol = 3;
    return (void *) uinp;
  }

  return nullptr;
}

/* ----------------------------------------------------------------------
  grow local vectors and arrays if necessary
  keep them all atom->nmax in length even if ghost storage not needed
------------------------------------------------------------------------- */

void PairAmoeba::grow_local()
{
 // free vectors and arrays

  memory->destroy(xaxis2local);
  memory->destroy(yaxis2local);
  memory->destroy(zaxis2local);
  memory->destroy(rpole);
  memory->destroy(tq);

  if (amoeba) {
    memory->destroy(red2local);
    memory->destroy(xred);
  }

  memory->destroy(uind);
  memory->destroy(uinp);
  memory->destroy(udirp);
  if (poltyp == OPT) {
    memory->destroy(uopt);
    memory->destroy(uoptp);
    memory->destroy(fopt);
    memory->destroy(foptp);
  }

  memory->destroy(field);
  memory->destroy(fieldp);
  memory->destroy(ufld);
  memory->destroy(dufld);
  memory->destroy(rsd);
  memory->destroy(rsdp);
  memory->destroy(zrsd);
  memory->destroy(zrsdp);

  // multipole

  memory->destroy(cmp);
  memory->destroy(fmp);
  memory->destroy(cphi);
  memory->destroy(fphi);

  // induce

  memory->destroy(poli);
  memory->destroy(conj);
  memory->destroy(conjp);
  memory->destroy(vec);
  memory->destroy(vecp);
  memory->destroy(udir);
  memory->destroy(usum);
  memory->destroy(usump);

  memory->destroy(fuind);
  memory->destroy(fuinp);
  memory->destroy(fdip_phi1);
  memory->destroy(fdip_phi2);
  memory->destroy(fdip_sum_phi);
  memory->destroy(dipfield1);
  memory->destroy(dipfield2);

  // polar

  memory->destroy(fphid);
  memory->destroy(fphip);
  memory->destroy(fphidp);
  memory->destroy(cphidp);

  if (use_ewald || use_dewald) {
    memory->destroy(thetai1);
    memory->destroy(thetai2);
    memory->destroy(thetai3);
    memory->destroy(igrid);
  }

  // dipole and PCG neighbor lists

  memory->destroy(numneigh_dipole);
  memory->sfree(firstneigh_dipole);
  memory->sfree(firstneigh_dipdip);

  if (poltyp == MUTUAL && pcgprec) {
    memory->destroy(numneigh_precond);
    memory->sfree(firstneigh_precond);
    memory->sfree(firstneigh_pcpc);
  }

  // reset nmax

  nmax = atom->nmax;

  // re-allocate vectors and arrays

  memory->create(xaxis2local,nmax,"amoeba:xaxis2local");
  memory->create(yaxis2local,nmax,"amoeba:yaxis2local");
  memory->create(zaxis2local,nmax,"amoeba:zaxis2local");
  memory->create(rpole,nmax,13,"amoeba:rpole");
  memory->create(tq,nmax,3,"amoeba:tq");

  if (amoeba) {
    memory->create(red2local,nmax,"amoeba:red2local");
    memory->create(xred,nmax,3,"amoeba:xred");
  }

  // note use of optorder+1 for uopt and uoptp

  memory->create(uind,nmax,3,"amoeba:uind");
  memory->create(uinp,nmax,3,"amoeba:uinp");
  memory->create(udirp,nmax,3,"amoeba:uinp");
  if (poltyp == OPT) {
    memory->create(uopt,nmax,optorder+1,3,"amoeba:uopt");
    memory->create(uoptp,nmax,optorder+1,3,"amoeba:uopt");
    memory->create(fopt,nmax,optorder,10,"amoeba:fopt");
    memory->create(foptp,nmax,optorder,10,"amoeba:foptp");
  }

  memory->create(field,nmax,3,"amoeba:field");
  memory->create(fieldp,nmax,3,"amoeba:fieldp");
  memory->create(ufld,nmax,3,"amoeba:ufld");
  memory->create(dufld,nmax,6,"amoeba:dufld");
  memory->create(rsd,nmax,3,"amoeba:rsd");
  memory->create(rsdp,nmax,3,"amoeba:rsdp");
  memory->create(zrsd,nmax,3,"amoeba:zrsd");
  memory->create(zrsdp,nmax,3,"amoeba:zrsdp");

  // multipole

  memory->create(cmp,nmax,10,"ameoba/mpole:cmp");
  memory->create(fmp,nmax,10,"ameoba/mpole:fmp");
  memory->create(cphi,nmax,10,"ameoba/mpole:cphi");
  memory->create(fphi,nmax,20,"ameoba/mpole:fphi");

  // induce

  memory->create(poli,nmax,"ameoba/induce:poli");
  memory->create(conj,nmax,3,"ameoba/induce:conj");
  memory->create(conjp,nmax,3,"ameoba/induce:conjp");
  memory->create(vec,nmax,3,"ameoba/induce:vec");
  memory->create(vecp,nmax,3,"ameoba/induce:vecp");
  memory->create(udir,nmax,3,"ameoba/induce:udir");
  memory->create(usum,nmax,3,"ameoba/induce:usum");
  memory->create(usump,nmax,3,"ameoba/induce:usump");

  memory->create(fuind,nmax,3,"ameoba/induce:fuind");
  memory->create(fuinp,nmax,3,"ameoba/induce:fuinp");
  memory->create(fdip_phi1,nmax,10,"ameoba/induce:fdip_phi1");
  memory->create(fdip_phi2,nmax,10,"ameoba/induce:fdip_phi2");
  memory->create(fdip_sum_phi,nmax,20,"ameoba/induce:fdip_dum_phi");
  memory->create(dipfield1,nmax,3,"ameoba/induce:dipfield1");
  memory->create(dipfield2,nmax,3,"ameoba/induce:dipfield2");

  // polar

  memory->create(fphid,nmax,10,"polar:fphid");
  memory->create(fphip,nmax,10,"polar:fphip");
  memory->create(fphidp,nmax,20,"polar:fphidp");
  memory->create(cphidp,nmax,10,"polar:cphidp");

  if (use_ewald || use_dewald) {
    memory->create(thetai1,nmax,bsordermax,4,"amoeba:thetai1");
    memory->create(thetai2,nmax,bsordermax,4,"amoeba:thetai2");
    memory->create(thetai3,nmax,bsordermax,4,"amoeba:thetai3");
    memory->create(igrid,nmax,3,"amoeba:igrid");
  }

  memory->create(numneigh_dipole,nmax,"amoeba:numneigh_dipole");
  firstneigh_dipole = (int **)
    memory->smalloc(nmax*sizeof(int *),"induce:firstneigh_dipole");
  firstneigh_dipdip = (double **)
    memory->smalloc(nmax*sizeof(double *),"induce:firstneigh_dipdip");

  if (poltyp == MUTUAL && pcgprec) {
    memory->create(numneigh_precond,nmax,"amoeba:numneigh_precond");
    firstneigh_precond = (int **)
      memory->smalloc(nmax*sizeof(int *),"induce:firstneigh_precond");
    firstneigh_pcpc = (double **)
      memory->smalloc(nmax*sizeof(double *),"induce:firstneigh_pcpc");
  }

  memory->create(_moduli_array,bsordermax,"amoeba:_moduli_array");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays and FFTs
------------------------------------------------------------------------- */

double PairAmoeba::memory_usage()
{
  double bytes = 0.0;

  bytes += (double) 3 * nmax * sizeof(int);       // xyz axis2local
  bytes += (double) 13 * nmax * sizeof(double);   // rpole
  bytes += (double) 3 * nmax * sizeof(double);    // tq

  if (amoeba) {
    bytes += (double) nmax * sizeof(int);          // red22local
    bytes += (double) 3 * nmax * sizeof(double);   // xred
  }

  bytes += (double) 9 * nmax * sizeof(double);         // uind/uinp/udirp
  if (poltyp == OPT) {
    bytes += (double) 6 * (optorder+1) * nmax * sizeof(double);   // uopt/uoptp
    bytes += (double) 20 * optorder * nmax * sizeof(double);      // fopt/foptp
  }

  bytes += (double) 15 * nmax * sizeof(double);   // field/fieldp/ufld/dufld
  bytes += (double) 12 * nmax * sizeof(double);   // rsd/rsdp/zrsd/zrsdp

  bytes += (double) 50 * nmax * sizeof(double);   // cmp/fmp/cphi/fphi

  bytes += (double) nmax * sizeof(double);        // poli
  bytes += (double) 12 * nmax * sizeof(double);   // conj/conjp/vec/vecp
  bytes += (double) 9 * nmax * sizeof(double);    // udir/usum/usump

  bytes += (double) 6 * nmax * sizeof(double);    // fuind/fuinp
  bytes += (double) 20 * nmax * sizeof(double);   // fdip_phi1/fdip_phi1
  bytes += (double) 20 * nmax * sizeof(double);   // fdip_sum_phi
  bytes += (double) 6 * nmax * sizeof(double);    // dipfield1/dipfield2

  bytes += (double) 50 * nmax * sizeof(double);   // fphid/fphip/fphidp/cphidp

  if (use_ewald || use_dewald) {
    bytes += (double) 12 * bsordermax * nmax *sizeof(double);   // theta123
    bytes += (double) 3 * nmax *sizeof(int);                    // igrid
  }

  bytes += (double) nmax * sizeof(int);        // numneigh_dipole
  bytes += (double) nmax * sizeof(int *);      // firstneigh_dipole
  bytes += (double) nmax * sizeof(double *);   // firstneigh_dipdip
  for (int i = 0; i < comm->nthreads; i++) {   // 2 neighbor lists
    bytes += ipage_dipole[i].size();
    bytes += dpage_dipdip[i].size();
  }

  if (poltyp == MUTUAL && pcgprec) {
    bytes += (double) nmax * sizeof(int);        // numneigh_rpecond
    bytes += (double) nmax * sizeof(int *);      // firstneigh_precond
    bytes += (double) nmax * sizeof(double *);   // firstneigh_pcpc
    for (int i = 0; i < comm->nthreads; i++) {   // 2 neighbor lists
      bytes += ipage_precond[i].size();
      bytes += dpage_pcpc[i].size();
    }
  }

  return bytes;
}
