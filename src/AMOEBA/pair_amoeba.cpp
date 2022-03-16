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

#include "pair_amoeba.h"
#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include "amoeba_convolution.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "fix_store.h"
#include "fft3d_wrap.h"
#include "gridcomm.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

enum{INDUCE,RSD,SETUP_AMOEBA,SETUP_HIPPO,KMPOLE,AMGROUP,PVAL};  // forward comm
enum{FIELD,ZRSD,TORQUE,UFLD};                                   // reverse comm
enum{ARITHMETIC,GEOMETRIC,CUBIC_MEAN,R_MIN,SIGMA,DIAMETER,HARMONIC,HHG,W_H};
enum{HAL,REPULSE,QFER,DISP,MPOLE,POLAR,USOLV,DISP_LONG,MPOLE_LONG,POLAR_LONG};
enum{MPOLE_GRID,POLAR_GRID,POLAR_GRIDC,DISP_GRID,INDUCE_GRID,INDUCE_GRIDC};
enum{MUTUAL,OPT,TCG,DIRECT};
enum{GEAR,ASPC,LSQR};

#define DELTASTACK 16

#define UIND_DEBUG 0       // also in amoeba_induce.cpp

/* ---------------------------------------------------------------------- */

PairAmoeba::PairAmoeba(LAMMPS *lmp) : Pair(lmp)
{
  // error checks

  if (strcmp(update->unit_style,"real") != 0) 
    error->all(FLERR,"Pair style amoeba/hippo require real units");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style amoeba/hippo require newton pair on");
  if (domain->dimension == 2)
    error->all(FLERR,"Pair style amoeba/hippo requires 3d system");

  me = comm->me;
  nprocs = comm->nprocs;
  
  // force field settings

  one_coeff = 1;
  single_enable = 0;

  amoeba = 1;
  hippo = 0;

  optorder = 0;
  maxualt = 0;
  tcgnab = 0;
  
  nmax = 0;
  xaxis2local = yaxis2local = zaxis2local = NULL;
  rpole = NULL;
  tq = NULL;

  ired2local = NULL;
  xred = NULL;

  uind = uinp = udirp = NULL;
  uopt = uoptp = NULL;
  fopt = foptp = NULL;
  field = fieldp = NULL;
  ufld = dufld = NULL;
  rsd = rsdp = NULL;
  zrsd = zrsdp = NULL;

  bsordermax = 0;
  thetai1 = thetai2 = thetai3 = NULL;
  bsmod1 = bsmod2 = bsmod3 = NULL;
  bsbuild = NULL;
  igrid = NULL;

  m_kspace = p_kspace = pc_kspace = d_kspace = NULL;
  i_kspace = ic_kspace = NULL;
  
  numneigh_dipole = NULL;
  firstneigh_dipole = NULL;
  firstneigh_dipdip = NULL;
  ipage_dipole = NULL;
  dpage_dipdip = NULL;

  numneigh_precond = NULL;
  firstneigh_precond = NULL;
  ipage_precond = NULL;

  firstneigh_pcpc = NULL;
  dpage_pcpc = NULL;

  initialize_type_class();
  initialize_vdwl();
  initialize_smallsize();

  forcefield = NULL;

  id_pole = id_udalt = id_upalt = NULL;

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
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) {
    if (id_pole) modify->delete_fix(id_pole);
    if (id_udalt) modify->delete_fix(id_udalt);
    if (id_upalt) modify->delete_fix(id_upalt);
  }

  delete [] id_pole;
  delete [] id_udalt;
  delete [] id_upalt;

  memory->destroy(xaxis2local);
  memory->destroy(yaxis2local);
  memory->destroy(zaxis2local);
  memory->destroy(rpole);
  memory->destroy(tq);

  memory->destroy(ired2local);
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

  memory->destroy(thetai1);
  memory->destroy(thetai2);
  memory->destroy(thetai3);
  memory->destroy(bsmod1);
  memory->destroy(bsmod2);
  memory->destroy(bsmod3);
  memory->destroy(bsbuild);
  memory->destroy(igrid);

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

  delete [] forcefield;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

  delete [] factors;

  // DEBUG

  if (me == 0 && UIND_DEBUG) fclose(fp_uind);
}

/* ---------------------------------------------------------------------- */

void PairAmoeba::compute(int eflag, int vflag)
{
  // DEBUG timer init

  if (update->ntimestep <= 1) {
    time_init = time_hal = time_repulse = time_disp = time_mpole = 0.0;
    time_induce = time_polar = time_qxfer = 0.0;
  }

  double evdwl;

  evdwl = 0.0;
  ev_init(eflag,vflag);

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

  if (atom->nlocal + atom->nghost > nmax) grow_local();

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
    comm->forward_comm_pair(this);
    kmpole();
    
    if (hippo) {
      pval = atom->dvector[index_pval];
      double **pole = fixpole->astore;
      int nlocal = atom->nlocal;
      int itype,iclass;
      for (int i = 0; i < nlocal; i++) {
	itype = amtype[i];
	iclass = amtype2class[itype];
	pval[i] = pole[i][0] - pcore[iclass];
      }
      cfstyle = PVAL;
      comm->forward_comm_pair(this);
    }
  }

  first_flag_compute = 0;

  // -------------------------------------------------------------------
  // end of one-time initializations
  // -------------------------------------------------------------------

  double time0,time1,time2,time3,time4,time5,time6,time7,time8;

  MPI_Barrier(world);
  time0 = MPI_Wtime();

  // if reneighboring step:
  // augment neighbor list to include 1-5 neighbor flags
  // re-create xyz axis2local and ired2local
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
      j = ired2local[i];
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
  comm->forward_comm_pair(this);

  if (amoeba) pbc_xred();

  time1 = MPI_Wtime();

  // ----------------------------------------
  // compute components of force field
  // ----------------------------------------

  // buffered 14-7 Vdwl, pairwise

  if (amoeba && hal_flag) hal();
  time2 = MPI_Wtime();

  // Pauli repulsion, pairwise

  if (hippo && repulse_flag) repulsion();
  time3 = MPI_Wtime();

  // Ewald dispersion, pairwise and long range

  if (hippo && (disp_rspace_flag || disp_kspace_flag)) dispersion();
  time4 = MPI_Wtime();

  // multipole, pairwise and long range

  if (mpole_rspace_flag || mpole_kspace_flag) multipole();
  time5 = MPI_Wtime();

  // induced dipoles, interative CG relaxation
  // communicate induce() output values needed by ghost atoms

  if (polar_rspace_flag || polar_kspace_flag) {
    induce();
    cfstyle = INDUCE;
    comm->forward_comm_pair(this);
  }
  time6 = MPI_Wtime();

  // dipoles, pairwise and long range

  if (polar_rspace_flag || polar_kspace_flag) polar();
  time7 = MPI_Wtime();

  // charge transfer, pairwise

  if (hippo && qxfer_flag) charge_transfer();

  time8 = MPI_Wtime();

  // energy, force, virial summations

  eng_vdwl = ehal + erepulse + edisp + empole + epolar + eqxfer;

  double **f = atom->f;
  nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (int i = 0; i < 6; i++)
    virial[i] = virhal[i] + virrepulse[i] + virdisp[i] + 
      virpolar[i] + virmpole[i] + virqxfer[i];

  // virial computation
  // NOTE: how does this work for AMOEBA ?
  // do all terms get summed this way, or only pairwise
  // it is 2x what summed virial[6] above is
  
  // if (vflag_fdotr) virial_fdotr_compute();

  // timing information

  time_init    += time1 - time0;
  time_hal     += time2 - time1;
  time_repulse += time3 - time2;
  time_disp    += time4 - time3;
  time_mpole   += time5 - time4;
  time_induce  += time6 - time5;
  time_polar   += time7 - time6;
  time_qxfer   += time8 - time7;

  // timing output

  if (update->ntimestep < update->laststep) return;

  double ave;
  MPI_Allreduce(&time_init,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_init = ave/nprocs;

  MPI_Allreduce(&time_hal,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_hal = ave/nprocs;

  MPI_Allreduce(&time_repulse,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_repulse = ave/nprocs;

  MPI_Allreduce(&time_disp,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_disp = ave/nprocs;

  MPI_Allreduce(&time_mpole,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_mpole = ave/nprocs;

  MPI_Allreduce(&time_induce,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_induce = ave/nprocs;

  MPI_Allreduce(&time_polar,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_polar = ave/nprocs;

  MPI_Allreduce(&time_qxfer,&ave,1,MPI_DOUBLE,MPI_SUM,world);
  time_qxfer = ave/nprocs;

  double time_amtot = time_init + time_hal + time_repulse + time_disp +
    time_mpole + time_induce + time_polar + time_qxfer;

  if (me == 0) {
    utils::logmesg(lmp,"\nAMEOBA/HIPPO timing info:\n");
    utils::logmesg(lmp,"  Init   time: {:.6g} {:.6g}\n",
                   time_init,time_init/time_amtot);
    if (amoeba)
      utils::logmesg(lmp,"  Hal    time: {:.6g} {:.6g}\n",
                     time_hal,time_hal/time_amtot*100);
    if (hippo)
      utils::logmesg(lmp,"  Repls  time: {:.6g} {:.6g}\n",
                     time_repulse,time_repulse/time_amtot*100);
    if (hippo)
      utils::logmesg(lmp,"  Disp   time: {:.6g} {:.6g}\n",
                     time_disp,time_disp/time_amtot*100);
    utils::logmesg(lmp,"  Mpole  time: {:.6g} {:.6g}\n",
                   time_mpole,time_mpole/time_amtot*100);
    utils::logmesg(lmp,"  Induce time: {:.6g} {:.6g}\n",
                   time_induce,time_induce/time_amtot*100);
    utils::logmesg(lmp,"  Polar  time: {:.6g} {:.6g}\n",
                   time_polar,time_polar/time_amtot*100);
    if (hippo)
      utils::logmesg(lmp,"  Qxfer  time: {:.6g} {:.6g}\n",
                     time_qxfer,time_qxfer/time_amtot*100);
    utils::logmesg(lmp,"  Total  time: {:.6g}\n\n",time_amtot);
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
------------------------------------------------------------------------- */

void PairAmoeba::settings(int narg, char **arg)
{
  // turn on all FF components by default
  
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
    else if (strcmp(arg[iarg],"disp/kspace") == 0) disp_rspace_flag = newvalue;
    else if (strcmp(arg[iarg],"polar") == 0) 
      polar_rspace_flag = polar_kspace_flag = newvalue;
    else if (strcmp(arg[iarg],"polar/rspace") == 0) polar_rspace_flag = newvalue;
    else if (strcmp(arg[iarg],"polar/kspace") == 0) polar_rspace_flag = newvalue;
    else if (strcmp(arg[iarg],"mpole") == 0) 
      mpole_rspace_flag = mpole_kspace_flag = newvalue;
    else if (strcmp(arg[iarg],"mpole/rspace") == 0) mpole_rspace_flag = newvalue;
    else if (strcmp(arg[iarg],"mpole/kspace") == 0) mpole_rspace_flag = newvalue;
    else if (strcmp(arg[iarg],"bond") == 0) bond_flag = newvalue;
    else if (strcmp(arg[iarg],"angle") == 0) angle_flag = newvalue;
    else if (strcmp(arg[iarg],"dihedral") == 0) dihedral_flag = newvalue;
    else if (strcmp(arg[iarg],"improper") == 0) improper_flag = newvalue;
    else if (strcmp(arg[iarg],"urey") == 0) urey_flag = newvalue;
    else if (strcmp(arg[iarg],"pitorsion") == 0) pitorsion_flag = newvalue;
    else if (strcmp(arg[iarg],"bitorsion") == 0) bitorsion_flag = newvalue;
    else error->all(FLERR,"Illegal pair_style command");
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAmoeba::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  if (narg < 3 && narg > 4)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // set setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 1;

  // read force field PRM file and optional KEY file

  set_defaults();
  read_prmfile(arg[2]);
  if (narg == 3) read_keyfile(NULL);
  else read_keyfile(arg[3]);

  // required for now, until 1-5 settings are implemented in LAMMPS

  //if (special_hal[4] != 1.0 || special_repel[4] != 1.0 ||
  //    special_disp[4] != 1.0 || special_mpole[4] != 1.0 ||
  //    special_polar_pscale[4] != 1.0 || special_polar_piscale[4] != 1.0 ||
  //    special_polar_wscale[4] != 1.0)
  //  error->all(FLERR,"AMOEBA 1-5 weights must be 1.0 for now");
  
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

  if (!atom->q_flag)
    error->all(FLERR,"Pair style amoeba/hippo requires atom attribute q");
  //if (!force->special_onefive)
  //  error->all(FLERR,"Pair style amoeba/hippo requires special_bonds one/five be set");

  // b/c polar uses mutipole virial terms

  if (apewald == aeewald && polar_kspace_flag && !mpole_kspace_flag)
    error->all(FLERR,
	       "Pair amoeba with apewald = aeewald requires mpole and polar together");
    
  // check if all custom atom arrays were set via fix property/atom

  int flag,cols;

  index_amtype = atom->find_custom("amtype",flag,cols);
  if (index_amtype < 0 || flag || cols) 
    error->all(FLERR,"Pair amoeba amtype is not defined");
  index_amgroup = atom->find_custom("amgroup",flag,cols);
  if (index_amgroup < 0 || flag || cols) 
    error->all(FLERR,"Pair amoeba amgroup is not defined");
  
  index_ired = atom->find_custom("ired",flag,cols);
  if (index_ired < 0 || flag || cols) 
    error->all(FLERR,"Pair amoeba ired is not defined");
  index_xaxis = atom->find_custom("xaxis",flag,cols);
  if (index_xaxis < 0 || flag || cols) 
    error->all(FLERR,"Pair amoeba xaxis is not defined");
  index_yaxis = atom->find_custom("yaxis",flag,cols);
  if (index_yaxis < 0 || flag || cols) 
    error->all(FLERR,"Pair amoeba yaxis is not defined");
  index_zaxis = atom->find_custom("zaxis",flag,cols);
  if (index_zaxis < 0 || flag || cols) 
    error->all(FLERR,"Pair amoeba zaxis is not defined");
  
  index_polaxe = atom->find_custom("polaxe",flag,cols);
  if (index_polaxe < 0 || flag || cols) 
    error->all(FLERR,"Pair amoeba polaxe is not defined");
  index_pval = atom->find_custom("pval",flag,cols);
  if (index_pval < 0 || !flag || cols) 
    error->all(FLERR,"Pair amoeba pval is not defined");

  // -------------------------------------------------------------------
  // one-time initializations
  // can't do earlier b/c need all atoms to exist
  // -------------------------------------------------------------------
  
  // creation of per-atom storage
  // create a new fix STORE style for each atom's pole vector
  // id = "AMOEBA_pole", fix group = all

  if (first_flag) {
    int n = strlen("AMOEBA_pole") + 1;
    id_pole = new char[n];
    strcpy(id_pole,"AMOEBA_pole");

    char **newarg = new char*[6];
    newarg[0] = id_pole;
    newarg[1] = group->names[0];
    newarg[2] = (char *) "STORE";
    newarg[3] = (char *) "peratom";
    newarg[4] = (char *) "1";
    newarg[5] = (char *) "13";
    modify->add_fix(6,newarg);
    fixpole = (FixStore *) modify->fix[modify->nfix-1];
    delete [] newarg;
  }

  // creation of per-atom storage
  // create 2 new fix STORE styles for each atom's induced dipole history info
  // id = "AMOEBA_udalt", fix group = all
  // id = "AMOEBA_upalt", fix group = all
  // only if using preconditioner

  if (first_flag && use_pred) {
    char ualt[8];
    sprintf(ualt,"%d",maxualt);

    int n = strlen("AMOEBA_udalt") + 1;
    id_udalt = new char[n];
    strcpy(id_udalt,"AMOEBA_udalt");

    char **newarg = new char*[7];
    newarg[0] = id_udalt;
    newarg[1] = group->names[0];
    newarg[2] = (char *) "STORE";
    newarg[3] = (char *) "peratom";
    newarg[4] = (char *) "1";
    newarg[5] = ualt;
    newarg[6] = (char *) "3";
    modify->add_fix(7,newarg);
    fixudalt = (FixStore *) modify->fix[modify->nfix-1];
    delete [] newarg;

    n = strlen("AMOEBA_upalt") + 1;
    id_upalt = new char[n];
    strcpy(id_udalt,"AMOEBA_upalt");

    newarg = new char*[7];
    newarg[0] = id_upalt;
    newarg[1] = group->names[0];
    newarg[2] = (char *) "STORE";
    newarg[3] = (char *) "peratom";
    newarg[4] = (char *) "1";
    newarg[5] = ualt;
    newarg[6] = (char *) "3";
    modify->add_fix(7,newarg);
    fixupalt = (FixStore *) modify->fix[modify->nfix-1];
    delete [] newarg;
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
  //   so need two GridComm instantiations for ghost comm
  
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

    delete [] csix_num_one;
    delete [] csix_num;
  }

  // initialize peratom pval to zero
  // so that initial ghost comm will be valid
  // pval is not set until first call to compute(), and only for HIPPO

  if (first_flag) {
    pval = atom->dvector[index_pval];
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) pval[i] = 0.0;
  }
  
  // output FF settings to screen and logfile
  
  if (first_flag) {
    if (comm->me == 0) {
      if (screen) print_settings(screen);
      if (logfile) print_settings(logfile);
    }
  }

  // all done with one-time initializations
  
  first_flag = 0;

  // -------------------------------------------------------------------
  // end of one-time initializations
  // -------------------------------------------------------------------

  // check for fixes which store persistent per-atom properties

  if (id_pole) {
    int ifix = modify->find_fix(id_pole);
    if (ifix < 0) error->all(FLERR,"Could not find pair amoeba fix ID");
    fixpole = (FixStore *) modify->fix[ifix];
  }

  if (id_udalt) {
    int ifix = modify->find_fix(id_udalt);
    if (ifix < 0) error->all(FLERR,"Could not find pair amoeba fix ID");
    fixudalt = (FixStore *) modify->fix[ifix];
    ifix = modify->find_fix(id_upalt);
    if (ifix < 0) error->all(FLERR,"Could not find pair amoeba fix ID");
    fixupalt = (FixStore *) modify->fix[ifix];
  }

  // assign hydrogen neighbors (ired) to each owned atom
  // only set if kred[i] is non-zero and I is bonded to a single atom
  // conceptually: non-zero if I is a hydrogen bonded to another atom
  // NOTE: ired needs to be a tagint vector?  but atom ivector is not!

  if (amoeba) {
    amtype = atom->ivector[index_amtype];
    ired = atom->ivector[index_ired];

    int **nspecial = atom->nspecial;
    int **special = atom->special;
    int nlocal = atom->nlocal;

    int itype,iclass;

    for (int i = 0; i < nlocal; i++) {
      itype = amtype[i];
      iclass = amtype2class[itype];
      if (kred[iclass] == 0.0) {
        ired[i] = 0;
      }
      else if (nspecial[i][0] != 1) {
        ired[i] = 0;
      }
      else {
        ired[i] = special[i][0];
      }
    }
  }

  // set KSpace recip matrix based on box size and shape

  lattice();

  // can now set comm size needed by this Pair
  // cfstyle KMPOLE is max # of 1-2 bond partners, smaller than comm_forward

  if (amoeba) comm_forward = 16;         // xred, rpole
  else if (hippo) comm_forward = 13;     // just rpole
  int fsize = 6;
  if (poltyp == OPT) fsize += 6*optorder;
  //if (poltyp == TCG) fsize += 12*tcgnab;
  comm_forward = MAX(comm_forward,fsize);

  comm_reverse = 9;

  // request neighbor lists

  int irequest = neighbor->request(this,instance_me);
  // for DEBUGGING with GPU
  //neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;

  // open debug output files
  // names are hard-coded

  if (me == 0) {
    char fname[32];
    sprintf(fname,"tmp.uind.kspace.%d",nprocs);
    if (UIND_DEBUG) fp_uind = fopen(fname,"w");
  }
}

/* ----------------------------------------------------------------------
   print settings to screen and logfile
------------------------------------------------------------------------- */

void PairAmoeba::print_settings(FILE *fp)
{
  fprintf(fp,"AMOEBA/HIPPO force field settings\n");

  if (amoeba) {
    choose(HAL);
    fprintf(fp,"  hal: cut %g taper %g vscale %g %g %g %g\n",
            sqrt(off2),sqrt(cut2),
            special_hal[1],special_hal[2],special_hal[3],special_hal[4]);
  }

  if (hippo) {
    choose(REPULSE);
    fprintf(fp,"  repulsion: cut %g taper %g rscale %g %g %g %g\n",
            sqrt(off2),sqrt(cut2),
            special_repel[1],special_repel[2],special_repel[3],special_repel[4]);
  }
  
  if (hippo) {
    choose(QFER);
    fprintf(fp,"  qxfer: cut %g taper %g mscale %g %g %g %g\n",
            sqrt(off2),sqrt(cut2),
            special_mpole[1],special_mpole[2],special_mpole[3],special_mpole[4]);
  }
  
  if (hippo) {
    if (use_dewald) {
      choose(DISP_LONG);
      fprintf(fp,"  dispersion: cut %g aewald %g bsorder %d "
              "FFT %d %d %d dspscale %g %g %g %g\n",
              sqrt(off2),aewald,bsdorder,ndfft1,ndfft2,ndfft3,
              special_disp[1],special_disp[2],special_disp[3],special_disp[4]);
    } else {
      choose(DISP);
      fprintf(fp,"  dispersion: cut %g aewald %g dspscale %g %g %g %g\n",
              sqrt(off2),aewald,
              special_disp[1],special_disp[2],special_disp[3],special_disp[4]);
    }
  }
  
  if (use_ewald) {
    choose(MPOLE_LONG);
    fprintf(fp,"  multipole: cut %g aewald %g bsorder %d "
	    "FFT %d %d %d mscale %g %g %g %g\n",
	    sqrt(off2),aewald,bseorder,nefft1,nefft2,nefft3,
	    special_mpole[1],special_mpole[2],special_mpole[3],special_mpole[4]);
  } else {
    choose(MPOLE);
    fprintf(fp,"  multipole: cut %g aewald %g mscale %g %g %g %g\n",
	    sqrt(off2),aewald,
	    special_mpole[1],special_mpole[2],special_mpole[3],special_mpole[4]);
  }
  
  if (use_ewald) {
    choose(POLAR_LONG);
    fprintf(fp,"  polar: cut %g aewald %g bsorder %d FFT %d %d %d\n",
	    sqrt(off2),aewald,bsporder,nefft1,nefft2,nefft3); 
    fprintf(fp,"         pscale %g %g %g %g piscale %g %g %g %g "
	    "wscale %g %g %g %g d/u scale %g %g\n",
	    special_polar_pscale[1],special_polar_pscale[2],
	    special_polar_pscale[3],special_polar_pscale[4],
	    special_polar_piscale[1],special_polar_piscale[2],
	    special_polar_piscale[3],special_polar_piscale[4],
	    special_polar_wscale[1],special_polar_wscale[2],
	    special_polar_wscale[3],special_polar_wscale[4],
	    polar_dscale,polar_uscale);
  } else {
    choose(POLAR);
    fprintf(fp,"  polar: cut %g aewald %g\n",sqrt(off2),aewald);
    fprintf(fp,"         pscale %g %g %g %g piscale %g %g %g %g "
	    "wscale %g %g %g %g d/u scale %g %g\n",
	    special_polar_pscale[1],special_polar_pscale[2],
	    special_polar_pscale[3],special_polar_pscale[4],
	    special_polar_piscale[1],special_polar_piscale[2],
	    special_polar_piscale[3],special_polar_piscale[4],
	    special_polar_wscale[1],special_polar_wscale[2],
	    special_polar_wscale[3],special_polar_wscale[4],
	    polar_dscale,polar_uscale);
  }
  
  choose(USOLV);
  fprintf(fp,"  precondition: cut %g\n",sqrt(off2));
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAmoeba::init_one(int i, int j)
{
  double cutoff = 0.0;

  if (amoeba) {
    choose(HAL);
    cutoff = MAX(cutoff,sqrt(off2));
  }

  if (hippo) {
    choose(REPULSE);
    cutoff = MAX(cutoff,sqrt(off2));
  }

  if (hippo) {
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

  if (hippo) {
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
    int **special = atom->special;
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
    int **special = atom->special;
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
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void PairAmoeba::pack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  if (flag == MPOLE_GRID) {
    FFT_SCALAR *src = m_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == POLAR_GRID) {
    FFT_SCALAR *src = p_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == POLAR_GRIDC) {
    FFT_SCALAR *src = pc_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src[2*list[i]];
      buf[n++] = src[2*list[i]+1];
    }
  } else if (flag == DISP_GRID) {
    FFT_SCALAR *src = d_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == INDUCE_GRID) {
    FFT_SCALAR *src = i_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == INDUCE_GRIDC) {
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

void PairAmoeba::unpack_forward_grid(int flag, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  if (flag == MPOLE_GRID) {
    FFT_SCALAR *dest = m_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (flag == POLAR_GRID) {
    FFT_SCALAR *dest = p_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (flag == POLAR_GRIDC) {
    FFT_SCALAR *dest = pc_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      dest[2*list[i]] = buf[n++];
      dest[2*list[i]+1] = buf[n++];
    }
  } else if (flag == DISP_GRID) {
    FFT_SCALAR *dest = d_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (flag == INDUCE_GRID) {
    FFT_SCALAR *dest = i_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (flag == INDUCE_GRIDC) {
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

void PairAmoeba::pack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  if (flag == MPOLE_GRID) {
    FFT_SCALAR *src = m_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == POLAR_GRID) {
    FFT_SCALAR *src = p_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == POLAR_GRIDC) {
    FFT_SCALAR *src = pc_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src[2*list[i]];
      buf[n++] = src[2*list[i]+1];
    }
  } else if (flag == DISP_GRID) {
    FFT_SCALAR *src = d_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == INDUCE_GRID) {
    FFT_SCALAR *src = i_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == INDUCE_GRIDC) {
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

void PairAmoeba::unpack_reverse_grid(int flag, void *vbuf, int nlist, int *list)
{
  FFT_SCALAR *buf = (FFT_SCALAR *) vbuf;

  if (flag == MPOLE_GRID) {
    FFT_SCALAR *dest = m_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (flag == POLAR_GRID) {
    FFT_SCALAR *dest = p_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (flag == POLAR_GRIDC) {
    FFT_SCALAR *dest = pc_kspace->grid_brick_start;
    int n = 0;
    for (int i = 0; i < nlist; i++) {
      dest[2*list[i]] += buf[n++];
      dest[2*list[i]+1] += buf[n++];
    }
  } else if (flag == DISP_GRID) {
    FFT_SCALAR *dest = d_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (flag == INDUCE_GRID) {
    FFT_SCALAR *dest = i_kspace->grid_brick_start;
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[i];
  } else if (flag == INDUCE_GRIDC) {
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
  int *stack = NULL;

  int **special = atom->special;
  int **nspecial = atom->nspecial;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // initially, groupID = atomID
  // communicate new groupIDs to ghost atoms
  
  for (i = 0; i < nlocal; i++) amgroup[i] = tag[i];
  cfstyle = AMGROUP;
  comm->forward_comm_pair(this);

  // loop until no ghost atom groupIDs are reset

  while (1) {

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
    comm->forward_comm_pair(this);

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

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"  AMOEBA/HIPPO group count: " BIGINT_FORMAT "\n",allbcount);
    if (logfile)
      fprintf(logfile,"  AMOEBA/HIPPO group count: " BIGINT_FORMAT "\n",allbcount);
  }
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

  // set cutoffs and taper coeffs
  // NOTE: should this cutoff include skin = 2.0 ?
  // Josh is checking

  choose(USOLV);

  // atoms and neighbor list

  double **x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all induce neighs of owned atoms
  // scan full neighbor list of I

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
------------------------------------------------------------------------- */

void PairAmoeba::initialize_vdwl()
{
  radmin = radmin4 = epsilon = epsilon4 = NULL;
}

// NOTE: n_amclass may be much larger than actual atom classes ??
//       due to format of Tinker PRM file

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
  copt = copm = NULL;
  a_ualt = ap_ualt = NULL;
  b_ualt = bp_ualt = NULL;
  c_ualt = cp_ualt = NULL;
  bpred = bpredp = bpreds = bpredps = NULL;
  gear = aspc = NULL;
}

void PairAmoeba::allocate_smallsize()
{
  // NOTE: are optorder and maxualt always initialized ?
  //       maybe there should be if tests here

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
  delete [] copt;
  delete [] copm;
  delete [] a_ualt;
  delete [] ap_ualt;
  delete [] b_ualt;
  delete [] bp_ualt;
  memory->destroy(c_ualt);
  memory->destroy(cp_ualt);
  delete [] bpred;
  delete [] bpredp;
  delete [] bpreds;
  delete [] bpredps;
  delete [] gear;
  delete [] aspc;
}

/* ----------------------------------------------------------------------
   set cutoffs, taper constants, PME params for a FF component
------------------------------------------------------------------------- */

void PairAmoeba::choose(int which)
{
  double off,cut;

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

  double denom = pow(off-cut,5.0);
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
   NOTE: not yet processing explicit VDWL14 entries in force field file
------------------------------------------------------------------------- */

void PairAmoeba::mix()
{
  int i,j,m;
  double ei,ej,sei,sej,eij;
  double ri,rj,sri,srj,rij;

  double TWOSIX = pow(2.0,1.0/6.0);

  for (i = 1; i <= n_amclass; i++) {
    
    // printf("MIX i %d nclass %d eps %g sigma %g\n",
    //        i,n_amclass,vdwl_eps[i],vdwl_sigma[i]);
    
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
        eij = 2.0 * (sei*sej) * pow(ri*rj,3.0) / (pow(ri,6.0) + pow(rj,6.0));
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
  grow local vectors and arrays if necessary
  keep them atom->nmax in length
  NOTE: some of these do not need to grow unless nlocal > atom->nmax
        these are ones that never store ghost values
        could realloc them separately
        e.g. thetai,igrid,fopt
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
    memory->destroy(ired2local);
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

  memory->destroy(thetai1);
  memory->destroy(thetai2);
  memory->destroy(thetai3);
  memory->destroy(igrid);

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
  // note use of optorder+1 for uopt and uoptp

  memory->create(xaxis2local,nmax,"amoeba:xaxis2local");
  memory->create(yaxis2local,nmax,"amoeba:yaxis2local");
  memory->create(zaxis2local,nmax,"amoeba:zaxis2local");
  memory->create(rpole,nmax,13,"amoeba:rpole");
  memory->create(tq,nmax,3,"amoeba:tq");

  if (amoeba) {
    memory->create(ired2local,nmax,"amoeba:ired2local");
    memory->create(xred,nmax,3,"amoeba:xred");
  }

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
}

// ----------------------------------------------------------------------
// debug output methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   dump ID + 6 values from 2 (N,3) per-atom arrays
   only proc 0 can write
   file is already open
------------------------------------------------------------------------- */

void PairAmoeba::dump6(FILE *fp, const char *columns, double scale,
		       double **a, double **b)
{
  int i,j,m;
  MPI_Status status;
  MPI_Request request;
  
  // setup

  int size_one = 7;
  int nlocal = atom->nlocal;

  char boundstr[9];          // encoding of boundary flags
  domain->boundary_string(boundstr);

  int maxlocal;
  MPI_Allreduce(&nlocal,&maxlocal,1,MPI_INT,MPI_MAX,world);

  double *buf;
  memory->create(buf,maxlocal*size_one,"amoeba:buf");

  // pack my data

  tagint *tag = atom->tag;

  m = 0;
  for (i = 0; i < nlocal; i++) {
    buf[m++] = tag[i];
    buf[m++] = scale*a[i][0];
    buf[m++] = scale*a[i][1];
    buf[m++] = scale*a[i][2];
    buf[m++] = scale*b[i][0];
    buf[m++] = scale*b[i][1];
    buf[m++] = scale*b[i][2];
  }

  // write file
  
  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,BIGINT_FORMAT "\n",atom->natoms);
    fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[0],domain->boxhi[0]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[1],domain->boxhi[1]);
    fprintf(fp,"%-1.16e %-1.16e\n",domain->boxlo[2],domain->boxhi[2]);
    fprintf(fp,"ITEM: ATOMS %s\n",columns);
  }

  int nlines;
  double tmp;
  
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buf,maxlocal*size_one,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	nlines /= size_one;
      } else nlines = nlocal;
      
      m = 0;
      for (i = 0; i < nlines; i++) {
	for (j = 0; j < size_one; j++) {
	  if (j == 0) fprintf(fp,"%d",static_cast<tagint> (buf[m]));
	  else fprintf(fp," %g",buf[m]);
	  m++;
	}
	fprintf(fp,"\n");
      }
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(buf,nlocal*size_one,MPI_DOUBLE,0,0,world);
  }

  // clean up

  memory->destroy(buf);
}
