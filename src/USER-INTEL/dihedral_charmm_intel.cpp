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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "dihedral_charmm_intel.h"
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "error.h"

#include "suffix.h"
using namespace LAMMPS_NS;

#define PTOLERANCE (flt_t)1.05
#define MTOLERANCE (flt_t)-1.05
#define SMALL      (flt_t)0.001
typedef struct { int a,b,c,d,t;  } int5_t;

/* ---------------------------------------------------------------------- */

DihedralCharmmIntel::DihedralCharmmIntel(class LAMMPS *lmp)
  : DihedralCharmm(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

void DihedralCharmmIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    DihedralCharmm::compute(eflag, vflag);
    return;
  }
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers(),
                          force_const_single);
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    compute<double,double>(eflag, vflag, fix->get_double_buffers(),
                           force_const_double);
  else
    compute<float,float>(eflag, vflag, fix->get_single_buffers(),
                         force_const_single);
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void DihedralCharmmIntel::compute(int eflag, int vflag,
				  IntelBuffers<flt_t,acc_t> *buffers,
				  const ForceConst<flt_t> &fc)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = 0;

  // insure pair->ev_tally() will use 1-4 virial contribution

  if (weightflag && vflag_global == 2)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  if (evflag) {
    if (eflag) {
      if (force->newton_bond)
	eval<1,1,1>(vflag, buffers, fc);
      else
	eval<1,1,0>(vflag, buffers, fc);
    } else {
      if (force->newton_bond)
	eval<1,0,1>(vflag, buffers, fc);
      else
	eval<1,0,0>(vflag, buffers, fc);
    }
  } else {
    if (force->newton_bond)
      eval<0,0,1>(vflag, buffers, fc);
    else
      eval<0,0,0>(vflag, buffers, fc);
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_BOND, class flt_t, class acc_t>
void DihedralCharmmIntel::eval(const int vflag, 
			       IntelBuffers<flt_t,acc_t> *buffers,
			       const ForceConst<flt_t> &fc)

{
  const int inum = neighbor->ndihedrallist;
  if (inum == 0) return;

  ATOM_T * _noalias const x = buffers->get_x(0);
  flt_t * _noalias const q = buffers->get_q(0);
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  int f_stride;
  if (NEWTON_BOND) f_stride = buffers->get_stride(nall);
  else f_stride = buffers->get_stride(nlocal);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(0, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;

  acc_t oedihedral, ov0, ov1, ov2, ov3, ov4, ov5;
  acc_t oevdwl, oecoul, opv0, opv1, opv2, opv3, opv4, opv5;
  if (EVFLAG) {
    if (EFLAG)
      oevdwl = oecoul = oedihedral = (acc_t)0.0;
    if (vflag) {
      ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0.0;
      opv0 = opv1 = opv2 = opv3 = opv4 = opv5 = (acc_t)0.0;
    }
  }

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(f_start,f_stride,fc)		  \
    reduction(+:oevdwl,oecoul,oedihedral,ov0,ov1,ov2,ov3,ov4,ov5, \
	      opv0,opv1,opv2,opv3,opv4,opv5)
  #endif
  {
    int nfrom, nto, tid;
    IP_PRE_omp_range_id(nfrom, nto, tid, inum, nthreads);

    FORCE_T * _noalias const f = f_start + (tid * f_stride);
    if (fix->need_zero(tid))
      memset(f, 0, f_stride * sizeof(FORCE_T));

    const int5_t * _noalias const dihedrallist = 
      (int5_t *) neighbor->dihedrallist[0];
    const flt_t qqrd2e = force->qqrd2e;

    acc_t sedihedral, sv0, sv1, sv2, sv3, sv4, sv5;
    acc_t sevdwl, secoul, spv0, spv1, spv2, spv3, spv4, spv5;
    if (EVFLAG) {
      if (EFLAG)
	sevdwl = secoul = sedihedral = (acc_t)0.0;
      if (vflag) {
	sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0.0;
	spv0 = spv1 = spv2 = spv3 = spv4 = spv5 = (acc_t)0.0;
      }
    }

    for (int n = nfrom; n < nto; n++) {
      const int i1 = dihedrallist[n].a;
      const int i2 = dihedrallist[n].b;
      const int i3 = dihedrallist[n].c;
      const int i4 = dihedrallist[n].d;
      const int type = dihedrallist[n].t;

      // 1st bond

      const flt_t vb1x = x[i1].x - x[i2].x;
      const flt_t vb1y = x[i1].y - x[i2].y;
      const flt_t vb1z = x[i1].z - x[i2].z;
      const int itype = x[i1].w;

      // 2nd bond

      const flt_t vb2xm = x[i2].x - x[i3].x;
      const flt_t vb2ym = x[i2].y - x[i3].y;
      const flt_t vb2zm = x[i2].z - x[i3].z;

      // 3rd bond
      
      const flt_t vb3x = x[i4].x - x[i3].x;
      const flt_t vb3y = x[i4].y - x[i3].y;
      const flt_t vb3z = x[i4].z - x[i3].z;
      const int jtype = x[i4].w;

      // 1-4

      const flt_t delx = x[i1].x - x[i4].x;
      const flt_t dely = x[i1].y - x[i4].y;
      const flt_t delz = x[i1].z - x[i4].z;


      // c,s calculation

      const flt_t ax = vb1y*vb2zm - vb1z*vb2ym;
      const flt_t ay = vb1z*vb2xm - vb1x*vb2zm;
      const flt_t az = vb1x*vb2ym - vb1y*vb2xm;
      const flt_t bx = vb3y*vb2zm - vb3z*vb2ym;
      const flt_t by = vb3z*vb2xm - vb3x*vb2zm;
      const flt_t bz = vb3x*vb2ym - vb3y*vb2xm;

      const flt_t rasq = ax*ax + ay*ay + az*az;
      const flt_t rbsq = bx*bx + by*by + bz*bz;
      const flt_t rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
      const flt_t rg = sqrt(rgsq);

      flt_t rginv, ra2inv, rb2inv;
      rginv = ra2inv = rb2inv = (flt_t)0.0;
      if (rg > 0) rginv = (flt_t)1.0/rg;
      if (rasq > 0) ra2inv = (flt_t)1.0/rasq;
      if (rbsq > 0) rb2inv = (flt_t)1.0/rbsq;
      const flt_t rabinv = sqrt(ra2inv*rb2inv);

      flt_t c = (ax*bx + ay*by + az*bz)*rabinv;
      const flt_t s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

      // error check
      if (c > PTOLERANCE || c < MTOLERANCE) {
	int me = comm->me;

	if (screen) {
	  char str[128];
	  sprintf(str,"Dihedral problem: %d/%d " BIGINT_FORMAT " "
		  TAGINT_FORMAT " " TAGINT_FORMAT " "
		  TAGINT_FORMAT " " TAGINT_FORMAT,
		  me,tid,update->ntimestep,
		  atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
	  error->warning(FLERR,str,0);
	  fprintf(screen,"  1st atom: %d %g %g %g\n",
		  me,x[i1].x,x[i1].y,x[i1].z);
	  fprintf(screen,"  2nd atom: %d %g %g %g\n",
		  me,x[i2].x,x[i2].y,x[i2].z);
	  fprintf(screen,"  3rd atom: %d %g %g %g\n",
		  me,x[i3].x,x[i3].y,x[i3].z);
	  fprintf(screen,"  4th atom: %d %g %g %g\n",
		  me,x[i4].x,x[i4].y,x[i4].z);
	}
      }

      if (c > (flt_t)1.0) c = (flt_t)1.0;
      if (c < (flt_t)-1.0) c = (flt_t)-1.0;

      const flt_t tcos_shift = fc.bp[type].cos_shift;
      const flt_t tsin_shift = fc.bp[type].sin_shift;
      const flt_t tk = fc.bp[type].k;
      const int m = fc.bp[type].multiplicity;

      flt_t p = (flt_t)1.0;
      flt_t ddf1, df1;
      ddf1 = df1 = (flt_t)0.0;

      for (int i = 0; i < m; i++) {
	ddf1 = p*c - df1*s;
	df1 = p*s + df1*c;
	p = ddf1;
      }

      p = p*tcos_shift + df1*tsin_shift;
      df1 = df1*tcos_shift - ddf1*tsin_shift;
      df1 *= -m;
      p += (flt_t)1.0;
      
      if (m == 0) {
	p = (flt_t)1.0 + tcos_shift;
	df1 = (flt_t)0.0;
      }

      const flt_t fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
      const flt_t hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
      const flt_t fga = fg*ra2inv*rginv;
      const flt_t hgb = hg*rb2inv*rginv;
      const flt_t gaa = -ra2inv*rg;
      const flt_t gbb = rb2inv*rg;

      const flt_t dtfx = gaa*ax;
      const flt_t dtfy = gaa*ay;
      const flt_t dtfz = gaa*az;
      const flt_t dtgx = fga*ax - hgb*bx;
      const flt_t dtgy = fga*ay - hgb*by;
      const flt_t dtgz = fga*az - hgb*bz;
      const flt_t dthx = gbb*bx;
      const flt_t dthy = gbb*by;
      const flt_t dthz = gbb*bz;

      const flt_t df = -tk * df1;

      const flt_t sx2 = df*dtgx;
      const flt_t sy2 = df*dtgy;
      const flt_t sz2 = df*dtgz;

      flt_t f1x = df*dtfx;
      flt_t f1y = df*dtfy;
      flt_t f1z = df*dtfz;

      const flt_t f2x = sx2 - f1x;
      const flt_t f2y = sy2 - f1y;
      const flt_t f2z = sz2 - f1z;

      flt_t f4x = df*dthx;
      flt_t f4y = df*dthy;
      flt_t f4z = df*dthz;

      const flt_t f3x = -sx2 - f4x;
      const flt_t f3y = -sy2 - f4y;
      const flt_t f3z = -sz2 - f4z;

      if (EVFLAG) {
	flt_t deng;
	if (EFLAG) deng = tk * p;
	IP_PRE_ev_tally_dihed(EFLAG, eatom, vflag, deng, i1, i2, i3, i4, f1x, 
			      f1y, f1z, f3x, f3y, f3z, f4x, f4y, f4z, vb1x, 
			      vb1y, vb1z, -vb2xm, -vb2ym, -vb2zm, vb3x, vb3y, 
			      vb3z, sedihedral, f, NEWTON_BOND, nlocal,
			      sv0, sv1, sv2, sv3, sv4, sv5);
      }


      {
        if (NEWTON_BOND || i2 < nlocal) {
	  f[i2].x += f2x;
	  f[i2].y += f2y;
	  f[i2].z += f2z;
        }

        if (NEWTON_BOND || i3 < nlocal) {
	  f[i3].x += f3x;
	  f[i3].y += f3y;
	  f[i3].z += f3z;
        }
      }

      // 1-4 LJ and Coulomb interactions
      // tally energy/virial in pair, using newton_bond as newton flag

      const flt_t tweight = fc.weight[type];
      const flt_t rsq = delx*delx + dely*dely + delz*delz;
      const flt_t r2inv = (flt_t)1.0/rsq;
      const flt_t r6inv = r2inv*r2inv*r2inv;

      flt_t forcecoul;
      if (implicit) forcecoul = qqrd2e * q[i1]*q[i4]*r2inv;
      else forcecoul = qqrd2e * q[i1]*q[i4]*sqrt(r2inv);
      const flt_t forcelj = r6inv * (fc.ljp[itype][jtype].lj1*r6inv - 
				     fc.ljp[itype][jtype].lj2);
      const flt_t fpair = tweight * (forcelj+forcecoul)*r2inv;

      if (NEWTON_BOND || i1 < nlocal) {
	f1x += delx*fpair;
	f1y += dely*fpair;
	f1z += delz*fpair;
      }
      if (NEWTON_BOND || i4 < nlocal) {
	f4x -= delx*fpair;
	f4y -= dely*fpair;
	f4z -= delz*fpair;
      }

      if (EVFLAG) {
	flt_t ev_pre = (flt_t)0;
	if (NEWTON_BOND || i1 < nlocal)
	  ev_pre += (flt_t)0.5;
	if (NEWTON_BOND || i4 < nlocal)
	  ev_pre += (flt_t)0.5;

	if (EFLAG) {
	  flt_t ecoul, evdwl;
	  ecoul = tweight * forcecoul;
	  evdwl = tweight * r6inv * (fc.ljp[itype][jtype].lj3*r6inv - 
				     fc.ljp[itype][jtype].lj4);
	  secoul += ev_pre * ecoul;
	  sevdwl += ev_pre * evdwl;
	  if (eatom) {
	    evdwl *= (flt_t)0.5;
	    evdwl += (flt_t)0.5 * ecoul;
	    if (NEWTON_BOND || i1 < nlocal)
	      f[i1].w += evdwl;
	    if (NEWTON_BOND || i4 < nlocal)
	      f[i4].w += evdwl;
	  }
	}
	if (vflag) {                                                    
	  spv0 += ev_pre * delx * delx * fpair;                               
	  spv1 += ev_pre * dely * dely * fpair;                               
	  spv2 += ev_pre * delz * delz * fpair;                               
	  spv3 += ev_pre * delx * dely * fpair;                               
	  spv4 += ev_pre * delx * delz * fpair;                               
	  spv5 += ev_pre * dely * delz * fpair;                               
	}                                                                    
      }

      // apply force to each of 4 atoms
      {
        if (NEWTON_BOND || i1 < nlocal) {
	  f[i1].x += f1x;
	  f[i1].y += f1y;
	  f[i1].z += f1z;
        }

        if (NEWTON_BOND || i4 < nlocal) {
	  f[i4].x += f4x;
	  f[i4].y += f4y;
	  f[i4].z += f4z;
        }
      }
    } // for n
    if (EVFLAG) {
      if (EFLAG) {
	oedihedral += sedihedral;
	oecoul += secoul;
	oevdwl += sevdwl;
      }
      if (vflag) {
	ov0 += sv0; ov1 += sv1; ov2 += sv2; ov3 += sv3; ov4 += sv4; ov5 += sv5;
	opv0 += spv0; opv1 += spv1; opv2 += spv2; 
	opv3 += spv3; opv4 += spv4; opv5 += spv5;
      }
    }
  } // omp parallel

  if (EVFLAG) {
    if (EFLAG) {
      energy += oedihedral;
      force->pair->eng_vdwl += oevdwl;
      force->pair->eng_coul += oecoul;
    }
    if (vflag) {
      virial[0] += ov0; virial[1] += ov1; virial[2] += ov2;
      virial[3] += ov3; virial[4] += ov4; virial[5] += ov5;
      force->pair->virial[0] += opv0;
      force->pair->virial[1] += opv1;
      force->pair->virial[2] += opv2;
      force->pair->virial[3] += opv3;
      force->pair->virial[4] += opv4;
      force->pair->virial[5] += opv5;
    }
  }

  fix->set_reduce_flag();
}

/* ---------------------------------------------------------------------- */

void DihedralCharmmIntel::init_style()
{
  DihedralCharmm::init_style();

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  #ifdef _LMP_INTEL_OFFLOAD
  _use_base = 0;
  if (fix->offload_balance() != 0.0) {
    _use_base = 1;
    return;
  }
  #endif

  fix->bond_init_check();

  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    pack_force_const(force_const_single, fix->get_mixed_buffers());
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    pack_force_const(force_const_double, fix->get_double_buffers());
  else
    pack_force_const(force_const_single, fix->get_single_buffers());
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void DihedralCharmmIntel::pack_force_const(ForceConst<flt_t> &fc,
	                                   IntelBuffers<flt_t,acc_t> *buffers)
{

  const int tp1 = atom->ntypes + 1;
  const int bp1 = atom->ndihedraltypes + 1;
  fc.set_ntypes(tp1,bp1,memory);
  buffers->set_ntypes(tp1);

  for (int i = 0; i < tp1; i++) {
    for (int j = 0; j < tp1; j++) {
      fc.ljp[i][j].lj1 = lj14_1[i][j];
      fc.ljp[i][j].lj2 = lj14_2[i][j];
      fc.ljp[i][j].lj3 = lj14_3[i][j];
      fc.ljp[i][j].lj4 = lj14_4[i][j];
    }
  }

  for (int i = 0; i < bp1; i++) {
    fc.bp[i].multiplicity = multiplicity[i];
    fc.bp[i].cos_shift = cos_shift[i];
    fc.bp[i].sin_shift = sin_shift[i];
    fc.bp[i].k = k[i];
    fc.weight[i] = weight[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void DihedralCharmmIntel::ForceConst<flt_t>::set_ntypes(const int npairtypes,
            	                                        const int nbondtypes,
	                                                Memory *memory) {
  if (npairtypes != _npairtypes) {
    if (_npairtypes > 0)
      _memory->destroy(ljp);
    if (npairtypes > 0)
      memory->create(ljp,npairtypes,npairtypes,"fc.ljp");
  }

  if (nbondtypes != _nbondtypes) {
    if (_nbondtypes > 0) {
      _memory->destroy(bp);
      _memory->destroy(weight);
    }
    
    if (nbondtypes > 0) {
      _memory->create(bp,nbondtypes,"dihedralcharmmintel.bp");
      _memory->create(weight,nbondtypes,"dihedralcharmmintel.weight");
    }
  }
  _npairtypes = npairtypes;
  _nbondtypes = nbondtypes;
  _memory = memory;
}
