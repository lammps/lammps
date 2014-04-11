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
   Contributing authors: Christopher Weinberger (SNL), Stephen Foiles (SNL),
                         Chandra Veer Singh (Cornell)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_adp.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairADP::PairADP(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  nmax = 0;
  rho = NULL;
  fp = NULL;
  mu = NULL;
  lambda = NULL;

  setfl = NULL;

  frho = NULL;
  rhor = NULL;
  z2r = NULL;
  u2r = NULL;
  w2r = NULL;

  frho_spline = NULL;
  rhor_spline = NULL;
  z2r_spline = NULL;
  u2r_spline = NULL;
  w2r_spline = NULL;

  // set comm size needed by this Pair

  comm_forward = 10;
  comm_reverse = 10;

  single_enable = 0;
  one_coeff = 1;
  manybody_flag = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairADP::~PairADP()
{
  memory->destroy(rho);
  memory->destroy(fp);
  memory->destroy(mu);
  memory->destroy(lambda);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] type2frho;
    memory->destroy(type2rhor);
    memory->destroy(type2z2r);
    memory->destroy(type2u2r);
    memory->destroy(type2w2r);
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    memory->destroy(setfl->u2r);
    memory->destroy(setfl->w2r);
    delete setfl;
  }

  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(z2r);
  memory->destroy(u2r);
  memory->destroy(w2r);

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
  memory->destroy(u2r_spline);
  memory->destroy(w2r_spline);
}

/* ---------------------------------------------------------------------- */

void PairADP::compute(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double u2,u2p,w2,w2p,nu;
  double *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double delmux,delmuy,delmuz,trdelmu,tradellam;
  double adpx,adpy,adpz,fx,fy,fz;
  double sumlamxx,sumlamyy,sumlamzz,sumlamyz,sumlamxz,sumlamxy;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // grow local arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    memory->destroy(mu);
    memory->destroy(lambda);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
    memory->create(mu,nmax,3,"pair:mu");
    memory->create(lambda,nmax,6,"pair:lambda");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) {
      rho[i] = 0.0;
      mu[i][0] = 0.0; mu[i][1] = 0.0; mu[i][2] = 0.0;
      lambda[i][0] = 0.0; lambda[i][1] = 0.0; lambda[i][2] = 0.0;
      lambda[i][3] = 0.0; lambda[i][4] = 0.0; lambda[i][5] = 0.0;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      rho[i] = 0.0;
      mu[i][0] = 0.0; mu[i][1] = 0.0; mu[i][2] = 0.0;
      lambda[i][0] = 0.0; lambda[i][1] = 0.0; lambda[i][2] = 0.0;
      lambda[i][3] = 0.0; lambda[i][4] = 0.0; lambda[i][5] = 0.0;
    }
  }

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = u2r_spline[type2u2r[jtype][itype]][m];
        u2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        mu[i][0] += u2*delx;
        mu[i][1] += u2*dely;
        mu[i][2] += u2*delz;
        coeff = w2r_spline[type2w2r[jtype][itype]][m];
        w2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        lambda[i][0] += w2*delx*delx;
        lambda[i][1] += w2*dely*dely;
        lambda[i][2] += w2*delz*delz;
        lambda[i][3] += w2*dely*delz;
        lambda[i][4] += w2*delx*delz;
        lambda[i][5] += w2*delx*dely;

        if (newton_pair || j < nlocal) {
          // verify sign difference for mu and lambda
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          coeff = u2r_spline[type2u2r[itype][jtype]][m];
          u2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          mu[j][0] -= u2*delx;
          mu[j][1] -= u2*dely;
          mu[j][2] -= u2*delz;
          coeff = w2r_spline[type2w2r[itype][jtype]][m];
          w2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
          lambda[j][0] += w2*delx*delx;
          lambda[j][1] += w2*dely*dely;
          lambda[j][2] += w2*delz*delz;
          lambda[j][3] += w2*dely*delz;
          lambda[j][4] += w2*delx*delz;
          lambda[j][5] += w2*delx*dely;
        }
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (eflag) {
      phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      phi += 0.5*(mu[i][0]*mu[i][0]+mu[i][1]*mu[i][1]+mu[i][2]*mu[i][2]);
      phi += 0.5*(lambda[i][0]*lambda[i][0]+lambda[i][1]*
                  lambda[i][1]+lambda[i][2]*lambda[i][2]);
      phi += 1.0*(lambda[i][3]*lambda[i][3]+lambda[i][4]*
                  lambda[i][4]+lambda[i][5]*lambda[i][5]);
      phi -= 1.0/6.0*(lambda[i][0]+lambda[i][1]+lambda[i][2])*
        (lambda[i][0]+lambda[i][1]+lambda[i][2]);
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // communicate derivative of embedding function

  comm->forward_comm_pair(this);

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        r = sqrt(rsq);
        p = r*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // u2 = u
        // u2p = u'
        // w2 = w
        // w2p = w'
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

        coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = u2r_spline[type2u2r[itype][jtype]][m];
        u2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        u2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        coeff = w2r_spline[type2w2r[itype][jtype]][m];
        w2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        w2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r;
        phi = z2*recip;
        phip = z2p*recip - phi*recip;
        psip = fp[i]*rhojp + fp[j]*rhoip + phip;
        fpair = -psip*recip;

        delmux = mu[i][0]-mu[j][0];
        delmuy = mu[i][1]-mu[j][1];
        delmuz = mu[i][2]-mu[j][2];
        trdelmu = delmux*delx+delmuy*dely+delmuz*delz;
        sumlamxx = lambda[i][0]+lambda[j][0];
        sumlamyy = lambda[i][1]+lambda[j][1];
        sumlamzz = lambda[i][2]+lambda[j][2];
        sumlamyz = lambda[i][3]+lambda[j][3];
        sumlamxz = lambda[i][4]+lambda[j][4];
        sumlamxy = lambda[i][5]+lambda[j][5];
        tradellam = sumlamxx*delx*delx+sumlamyy*dely*dely+
          sumlamzz*delz*delz+2.0*sumlamxy*delx*dely+
          2.0*sumlamxz*delx*delz+2.0*sumlamyz*dely*delz;
        nu = sumlamxx+sumlamyy+sumlamzz;

        adpx = delmux*u2 + trdelmu*u2p*delx*recip +
          2.0*w2*(sumlamxx*delx+sumlamxy*dely+sumlamxz*delz) +
          w2p*delx*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*delx;
        adpy = delmuy*u2 + trdelmu*u2p*dely*recip +
          2.0*w2*(sumlamxy*delx+sumlamyy*dely+sumlamyz*delz) +
          w2p*dely*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*dely;
        adpz = delmuz*u2 + trdelmu*u2p*delz*recip +
          2.0*w2*(sumlamxz*delx+sumlamyz*dely+sumlamzz*delz) +
          w2p*delz*recip*tradellam - 1.0/3.0*nu*(w2p*r+2.0*w2)*delz;
        adpx*=-1.0; adpy*=-1.0; adpz*=-1.0;

        fx = delx*fpair+adpx;
        fy = dely*fpair+adpy;
        fz = delz*fpair+adpz;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
        }

        if (eflag) evdwl = phi;
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,
                                 fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairADP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n+1];
  memory->create(type2rhor,n+1,n+1,"pair:type2rhor");
  memory->create(type2z2r,n+1,n+1,"pair:type2z2r");
  memory->create(type2u2r,n+1,n+1,"pair:type2u2r");
  memory->create(type2w2r,n+1,n+1,"pair:type2w2r");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairADP::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read concatenated *.plt file
------------------------------------------------------------------------- */

void PairADP::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read ADP parameter file

  if (setfl) {
    for (i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    memory->destroy(setfl->u2r);
    memory->destroy(setfl->w2r);
    delete setfl;
  }
  setfl = new Setfl();
  read_file(arg[2]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < setfl->nelements; j++)
      if (strcmp(arg[i],setfl->elements[j]) == 0) break;
    if (j < setfl->nelements) map[i-2] = j;
    else error->all(FLERR,"No matching element in ADP potential file");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(i,setfl->mass[map[i]]);
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairADP::init_style()
{
  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairADP::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (setfl) cutmax = setfl->cut;
  cutforcesq = cutmax*cutmax;

  return cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */

void PairADP::read_file(char *filename)
{
  Setfl *file = setfl;

  // open potential file

  int me = comm->me;
  FILE *fp;
  char line[MAXLINE];

  if (me == 0) {
    fp = force->open_potential(filename);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open ADP potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read and broadcast header
  // extract element names from nelements line

  int n;
  if (me == 0) {
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  sscanf(line,"%d",&file->nelements);
  int nwords = atom->count_words(line);
  if (nwords != file->nelements + 1)
    error->all(FLERR,"Incorrect element names in ADP potential file");

  char **words = new char*[file->nelements+1];
  nwords = 0;
  strtok(line," \t\n\r\f");
  while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

  file->elements = new char*[file->nelements];
  for (int i = 0; i < file->nelements; i++) {
    n = strlen(words[i]) + 1;
    file->elements[i] = new char[n];
    strcpy(file->elements[i],words[i]);
  }
  delete [] words;

  if (me == 0) {
    fgets(line,MAXLINE,fp);
    sscanf(line,"%d %lg %d %lg %lg",
           &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
  }

  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  file->mass = new double[file->nelements];
  memory->create(file->frho,file->nelements,file->nrho+1,"pair:frho");
  memory->create(file->rhor,file->nelements,file->nr+1,"pair:rhor");
  memory->create(file->z2r,file->nelements,file->nelements,file->nr+1,
                 "pair:z2r");
  memory->create(file->u2r,file->nelements,file->nelements,file->nr+1,
                 "pair:u2r");
  memory->create(file->w2r,file->nelements,file->nelements,file->nr+1,
                 "pair:w2r");

  int i,j,tmp;
  for (i = 0; i < file->nelements; i++) {
    if (me == 0) {
      fgets(line,MAXLINE,fp);
      sscanf(line,"%d %lg",&tmp,&file->mass[i]);
    }
    MPI_Bcast(&file->mass[i],1,MPI_DOUBLE,0,world);

    if (me == 0) grab(fp,file->nrho,&file->frho[i][1]);
    MPI_Bcast(&file->frho[i][1],file->nrho,MPI_DOUBLE,0,world);
    if (me == 0) grab(fp,file->nr,&file->rhor[i][1]);
    MPI_Bcast(&file->rhor[i][1],file->nr,MPI_DOUBLE,0,world);
  }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fp,file->nr,&file->z2r[i][j][1]);
      MPI_Bcast(&file->z2r[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fp,file->nr,&file->u2r[i][j][1]);
      MPI_Bcast(&file->u2r[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fp,file->nr,&file->w2r[i][j][1]);
      MPI_Bcast(&file->w2r[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }

  // close the potential file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void PairADP::file2array()
{
  int i,j,m,n;
  int ntypes = atom->ntypes;

  // set function params directly from setfl file

  nrho = setfl->nrho;
  nr = setfl->nr;
  drho = setfl->drho;
  dr = setfl->dr;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of setfl elements + 1 for zero array

  nfrho = setfl->nelements + 1;
  memory->destroy(frho);
  memory->create(frho,nfrho,nrho+1,"pair:frho");

  // copy each element's frho to global frho

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nrho; m++) frho[i][m] = setfl->frho[i][m];

  // add extra frho of zeroes for non-ADP types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-ADP atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-ADP atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of setfl elements

  nrhor = setfl->nelements;
  memory->destroy(rhor);
  memory->create(rhor,nrhor,nr+1,"pair:rhor");

  // copy each element's rhor to global rhor

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nr; m++) rhor[i][m] = setfl->rhor[i][m];

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for setfl files, I,J mapping only depends on I
  // OK if map = -1 (non-APD atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = setfl->nelements * (setfl->nelements+1) / 2;
  memory->destroy(z2r);
  memory->create(z2r,nz2r,nr+1,"pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = setfl->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // OK if map = -1 (non-ADP atom in pair hybrid) b/c type2z2r not used

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) continue;
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }

  // ------------------------------------------------------------------
  // setup u2r arrays
  // ------------------------------------------------------------------

  // allocate u2r arrays
  // nu2r = N*(N+1)/2 where N = # of setfl elements

  nu2r = setfl->nelements * (setfl->nelements+1) / 2;
  memory->destroy(u2r);
  memory->create(u2r,nu2r,nr+1,"pair:u2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) u2r[n][m] = setfl->u2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // OK if map = -1 (non-ADP atom in pair hybrid) b/c type2z2r not used

  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) continue;
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2u2r[i][j] = n;
    }
  }

  // ------------------------------------------------------------------
  // setup w2r arrays
  // ------------------------------------------------------------------

  // allocate w2r arrays
  // nw2r = N*(N+1)/2 where N = # of setfl elements

  nw2r = setfl->nelements * (setfl->nelements+1) / 2;
  memory->destroy(w2r);
  memory->create(w2r,nw2r,nr+1,"pair:w2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) w2r[n][m] = setfl->w2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // OK if map = -1 (non-ADP atom in pair hybrid) b/c type2z2r not used

  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) continue;
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2w2r[i][j] = n;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairADP::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
  memory->destroy(u2r_spline);
  memory->destroy(w2r_spline);

  memory->create(frho_spline,nfrho,nrho+1,7,"pair:frho");
  memory->create(rhor_spline,nrhor,nr+1,7,"pair:rhor");
  memory->create(z2r_spline,nz2r,nr+1,7,"pair:z2r");
  memory->create(u2r_spline,nz2r,nr+1,7,"pair:u2r");
  memory->create(w2r_spline,nz2r,nr+1,7,"pair:w2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho,drho,frho[i],frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr,dr,rhor[i],rhor_spline[i]);

  for (int i = 0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],z2r_spline[i]);

  for (int i = 0; i < nu2r; i++)
    interpolate(nr,dr,u2r[i],u2r_spline[i]);

  for (int i = 0; i < nw2r; i++)
    interpolate(nr,dr,w2r[i],w2r_spline[i]);
}

/* ---------------------------------------------------------------------- */

void PairADP::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                    8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0*(spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void PairADP::grab(FILE *fp, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fp);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL," \t\n\r\f"))) list[i++] = atof(ptr);
  }
}

/* ---------------------------------------------------------------------- */

int PairADP::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = fp[j];
    buf[m++] = mu[j][0];
    buf[m++] = mu[j][1];
    buf[m++] = mu[j][2];
    buf[m++] = lambda[j][0];
    buf[m++] = lambda[j][1];
    buf[m++] = lambda[j][2];
    buf[m++] = lambda[j][3];
    buf[m++] = lambda[j][4];
    buf[m++] = lambda[j][5];
  }
  return 10;
}

/* ---------------------------------------------------------------------- */

void PairADP::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    fp[i] = buf[m++];
    mu[i][0] = buf[m++];
    mu[i][1] = buf[m++];
    mu[i][2] = buf[m++];
    lambda[i][0] = buf[m++];
    lambda[i][1] = buf[m++];
    lambda[i][2] = buf[m++];
    lambda[i][3] = buf[m++];
    lambda[i][4] = buf[m++];
    lambda[i][5] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairADP::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
  buf[m++] = rho[i];
  buf[m++] = mu[i][0];
  buf[m++] = mu[i][1];
  buf[m++] = mu[i][2];
  buf[m++] = lambda[i][0];
  buf[m++] = lambda[i][1];
  buf[m++] = lambda[i][2];
  buf[m++] = lambda[i][3];
  buf[m++] = lambda[i][4];
  buf[m++] = lambda[i][5];
  }
  return 10;
}

/* ---------------------------------------------------------------------- */

void PairADP::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
    mu[j][0] += buf[m++];
    mu[j][1] += buf[m++];
    mu[j][2] += buf[m++];
    lambda[j][0] += buf[m++];
    lambda[j][1] += buf[m++];
    lambda[j][2] += buf[m++];
    lambda[j][3] += buf[m++];
    lambda[j][4] += buf[m++];
    lambda[j][5] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairADP::memory_usage()
{
  double bytes = Pair::memory_usage();
  bytes += 21 * nmax * sizeof(double);
  return bytes;
}
