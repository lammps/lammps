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

#include "compute_pressure_cylinder.h"

#include <cmath>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "citeme.h"
#include "domain.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

static const char cite_compute_pressure_cylinder[] =
  "compute pressure/cylinder:\n\n"
  "@Article{Addington,\n"
  " author = {C. K. Addington, Y. Long, K. E. Gubbins},\n"
  " title = {The pressure in interfaces having cylindrical geometry},\n"
  " journal = {J.~Chem.~Phys.},\n"
  " year =    2018,\n"
  " volume =  149,\n"
  " pages =   {084109}\n"
  "}\n\n";

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Calculate the configurational components of the pressure tensor in
  cylindrical geometry, according to the formulation of Addington et al. (2018)
  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

ComputePressureCyl::ComputePressureCyl(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  Pr_temp(nullptr), Pr_all(nullptr), Pz_temp(nullptr), Pz_all(nullptr), Pphi_temp(nullptr),
  Pphi_all(nullptr), R(nullptr), Rinv(nullptr), R2(nullptr), PrAinv(nullptr), PzAinv(nullptr),
  R2kin(nullptr), density_temp(nullptr), invVbin(nullptr), density_all(nullptr),
  tangent(nullptr), ephi_x(nullptr), ephi_y(nullptr), binz(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_compute_pressure_cylinder);
  if (narg != 7) error->all(FLERR,"Illegal compute pressure/cylinder command");

  zlo=utils::numeric(FLERR,arg[3],false,lmp);
  zhi=utils::numeric(FLERR,arg[4],false,lmp);
  Rmax=utils::numeric(FLERR,arg[5],false,lmp);
  bin_width=utils::numeric(FLERR,arg[6],false,lmp);

  if ((bin_width <= 0.0) || (bin_width > Rmax))
    error->all(FLERR,"Illegal compute pressure/cylinder command");
  if ((zhi < zlo) || ((zhi-zlo) < bin_width))
    error->all(FLERR,"Illegal compute pressure/cylinder command");
  if ((zhi > domain->boxhi[2]) || (zlo < domain->boxlo[2]))
    error->all(FLERR,"Illegal compute pressure/cylinder command");

  nbins=(int)(Rmax/bin_width);
  nzbins=(int)((zhi-zlo)/bin_width);

  // NOTE: at 2^22 = 4.2M bins, we will be close to exhausting allocatable
  // memory on a 32-bit environment. so we use this as an upper limit.

  if ((nbins < 1) || (nzbins < 1) || (nbins > 2<<22) || (nzbins > 2<<22))
    error->all(FLERR,"Illegal compute pressure/cylinder command");

  array_flag=1;
  vector_flag=0;
  extarray=0;
  size_array_cols = 5;  // r, number density, Pr, Pphi, Pz
  size_array_rows = nbins;

  Pr_temp = new double[nbins];
  Pr_all = new double[nbins];
  Pz_temp = new double[nbins];
  Pz_all = new double[nbins];
  Pphi_temp = new double[nbins];
  Pphi_all = new double[nbins];
  R  = new double[nbins];
  R2 = new double[nbins];
  PrAinv = new double[nbins];
  PzAinv = new double[nbins];
  Rinv = new double[nbins];
  binz = new double[nzbins];

  R2kin = new double[nbins];
  density_temp = new double[nbins];
  invVbin = new double[nbins];
  density_all = new double[nbins];

  memory->create(array,nbins,5,"PN:array");

  nphi=360;
  tangent = new double[nphi];
  ephi_x = new double[nphi];
  ephi_y = new double[nphi];

  nktv2p = force->nktv2p;

}

/* ---------------------------------------------------------------------- */

ComputePressureCyl::~ComputePressureCyl()
{
  // count all of these for memory usage
  memory->destroy(array);
  delete [] R;
  delete [] Rinv;
  delete [] R2;
  delete [] R2kin;
  delete [] invVbin;
  delete [] density_temp;
  delete [] density_all;
  delete [] tangent;
  delete [] ephi_x;
  delete [] ephi_y;
  delete [] Pr_temp;
  delete [] Pr_all;
  delete [] Pz_temp;
  delete [] Pz_all;
  delete [] Pphi_temp;
  delete [] Pphi_all;
  delete [] PrAinv;
  delete [] PzAinv;
  delete [] binz;
}

/* ---------------------------------------------------------------------- */

void ComputePressureCyl::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"No pair style is defined for compute pressure/cylinder");
  if (force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute pressure/cylinder");

  double phi;

  for (int iphi = 0; iphi < nphi; iphi++) {
    phi=((double)iphi)*MY_PI/180.0;
    tangent[iphi]=tan(phi);
    ephi_x[iphi]=-sin(phi);
    ephi_y[iphi]=cos(phi);
  }

  for (int iq = 0; iq < nbins; iq++) {
    R[iq]=((double)iq+0.5)*bin_width;
    Rinv[iq]=1.0/R[iq];
    R2[iq]=R[iq]*R[iq];
    R2kin[iq]=(((double)iq)+1.0)*bin_width;
    R2kin[iq]*=R2kin[iq];
    PrAinv[iq]=1.0/(2.0*MY_PI*(zhi-zlo)*R[iq]);
  }
  PphiAinv=1.0/((zhi-zlo)*bin_width*2.0*(double)nphi);

  invVbin[0]=1.0/((zhi-zlo)*MY_PI*R2kin[0]);
  PzAinv[0]=1.0/(MY_PI*R2kin[0]*((double)nzbins));

  for (int jq = 1; jq < nbins; jq++) {
    invVbin[jq]=1.0/((zhi-zlo)*MY_PI*(R2kin[jq]-R2kin[jq-1]));
    PzAinv[jq]=1.0/(MY_PI*(R2kin[jq]-R2kin[jq-1])*((double)nzbins));
  }

  // need an occasional half neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;

  for (int zzz = 0; zzz < nzbins; zzz++) binz[zzz]=(((double)zzz)+0.5)*bin_width+zlo;

}

/* ---------------------------------------------------------------------- */

void ComputePressureCyl::init_list(int /* id */, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if flag is set, compute requested info about pair
------------------------------------------------------------------------- */

void ComputePressureCyl::compute_array()
{
  invoked_array = update->ntimestep;

  int ibin;

  // clear pressures
  for (ibin = 0; ibin < nbins; ibin++) {
    density_temp[ibin]=0.0;
    density_all[ibin]=0.0;
    Pr_temp[ibin]=0.0;
    Pr_all[ibin]=0.0;
    Pphi_temp[ibin]=0.0;
    Pphi_all[ibin]=0.0;
    Pz_temp[ibin]=0.0;
    Pz_all[ibin]=0.0;
  }

  // what processor am I?
  int me;
  MPI_Comm_rank(world,&me);

  int i,j,ii,jj,inum,jnum,itype,jtype;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // calculate number density (by radius)
  double temp_R2;
  for (i = 0; i < nlocal; i++) if ((x[i][2] < zhi) && (x[i][2] > zlo)) {
    temp_R2=x[i][0]*x[i][0]+x[i][1]*x[i][1];
    if (temp_R2 > R2kin[nbins-1]) continue; // outside of Rmax

    for (j = 0; j < nbins; j++) if (temp_R2 < R2kin[j]) break;

    density_temp[j]+=invVbin[j];
  }
  MPI_Allreduce(density_temp,density_all,nbins,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nbins; i++) array[i][1]=density_all[i]; // NEW

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for newton = 0 and J = ghost atom,
  //   need to insure I,J pair is only output by one proc
  //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()
  // for flag = 0, just count pair interactions within force cutoff
  // for flag = 1, calculate requested output fields

  Pair *pair = force->pair;
  double **cutsq = force->pair->cutsq;

  double r1=0.0;
  double r2=0.0;
  double risq,rjsq;
  double A,B,C,D;
  double alpha1,alpha2;
  double xi,yi,zi,dx,dy,dz;
  double xR,yR,zR,fn;
  double alpha,xL,yL,zL,L2,ftphi,ftz;
  double sqrtD;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    r1=x[i][0]*x[i][0]+x[i][1]*x[i][1];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      // itag = jtag is possible for long cutoffs that include images of self
      // do calculation only on appropriate processor
      if (newton_pair == 0 && j >= nlocal) {
        jtag = tag[j];
        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      r2=x[j][0]*x[j][0]+x[j][1]*x[j][1];

      // ri is smaller of r1 and r2
      if (r2 < r1) {
        risq=r2;
        rjsq=r1;
        xi=x[j][0];
        yi=x[j][1];
        zi=x[j][2];
        dx=x[i][0]-x[j][0];
        dy=x[i][1]-x[j][1];
        dz=x[i][2]-x[j][2];
      } else {
        risq=r1;
        rjsq=r2;
        xi=x[i][0];
        yi=x[i][1];
        zi=x[i][2];
        dx=x[j][0]-x[i][0];
        dy=x[j][1]-x[i][1];
        dz=x[j][2]-x[i][2];
      }

      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;

      pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

      A=dx*dx+dy*dy;
      B=2.0*(xi*dx+yi*dy);

      // normal pressure contribution P_rhorho
      for (ibin = 0; ibin < nbins; ibin++) {
        // completely inside of R
        if (rjsq < R2[ibin]) continue;

        C=risq-R2[ibin];
        D=B*B-4.0*A*C;

        // completely outside of R
        if (D < 0.0) continue;

        sqrtD=sqrt(D);
        alpha1=0.5*(-B+sqrtD)/A;
        alpha2=0.5*(-B-sqrtD)/A;

        if ((alpha1 > 0.0) && (alpha1 < 1.0)) {
          zR=zi+alpha1*dz;
          if ((zR < zhi) && (zR > zlo))
          {
            xR=xi+alpha1*dx;
            yR=yi+alpha1*dy;
            fn=fpair*fabs(xR*dx+yR*dy);

            Pr_temp[ibin]+=fn;
          }
        }
        if ((alpha2 > 0.0) && (alpha2 < 1.0)) {
          zR=zi+alpha2*dz;
          if ((zR < zhi) && (zR > zlo)) {
            xR=xi+alpha2*dx;
            yR=yi+alpha2*dy;
            fn=fpair*fabs(xR*dx+yR*dy);

            Pr_temp[ibin]+=fn;
          }
        }
      }

      // azimuthal pressure contribution (P_phiphi)
      for (int iphi = 0; iphi < nphi; iphi++) {
        alpha=(yi-xi*tangent[iphi])/(dx*tangent[iphi]-dy);

        // no intersection with phi surface
        if ((alpha >= 1.0) || (alpha <= 0.0)) continue;

        // no contribution (outside of averaging region)
        zL=zi+alpha*dz;
        if ((zL > zhi) || (zL < zlo)) continue;

        xL=xi+alpha*dx;
        yL=yi+alpha*dy;

        L2=xL*xL+yL*yL;

        // no intersection (outside of Rmax)
        if (L2 > R2kin[nbins-1]) continue;

        ftphi=fabs(dx*ephi_x[iphi]+dy*ephi_y[iphi])*fpair;

        // add to appropriate bin
        for (ibin = 0; ibin < nbins; ibin++) if (L2 < R2kin[ibin]) {
          Pphi_temp[ibin]+=ftphi;
          break;
        }
      }

      // z pressure contribution (P_zz)
      for (int zbin = 0; zbin < nzbins; zbin++) {
        // check if interaction contributes
        if ((x[i][2] > binz[zbin]) && (x[j][2] > binz[zbin])) continue;
        if ((x[i][2] < binz[zbin]) && (x[j][2] < binz[zbin])) continue;

        alpha=(binz[zbin]-zi)/dz;

        xL=xi+alpha*dx;
        yL=yi+alpha*dy;

        L2=xL*xL+yL*yL;

        if (L2 > R2kin[nbins-1]) continue;

        ftz=fabs(dz)*fpair;

        // add to appropriate bin
        for (ibin = 0; ibin < nbins; ibin++) if (L2 < R2kin[ibin]) {
          Pz_temp[ibin]+=ftz;
          break;
        }
      }
    }
  }

  // calculate pressure (force over area)
  for (ibin = 0; ibin < nbins; ibin++) {
    Pr_temp[ibin]*=PrAinv[ibin]*Rinv[ibin];
    Pphi_temp[ibin]*=PphiAinv;
    Pz_temp[ibin]*=PzAinv[ibin];
  }

  // communicate these values across processors
  MPI_Allreduce(Pr_temp,Pr_all,nbins,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(Pphi_temp,Pphi_all,nbins,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(Pz_temp,Pz_all,nbins,MPI_DOUBLE,MPI_SUM,world);

  // populate array
  for (ibin = 0; ibin < nbins; ibin++) {
    array[ibin][0]=R[ibin];
    array[ibin][2]=Pr_all[ibin]*nktv2p;
    array[ibin][3]=Pphi_all[ibin]*nktv2p;
    array[ibin][4]=Pz_all[ibin]*nktv2p;
  }

}

/* ----------------------------------------------------------------------
memory usage of data
------------------------------------------------------------------------- */

double ComputePressureCyl::memory_usage()
{
  double bytes =
  (3.0*(double)nphi + 16.0*(double)nbins+5.0*(double)nbins) * sizeof(double);
  return bytes;
}
