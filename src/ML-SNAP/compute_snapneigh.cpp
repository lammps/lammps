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

#include "compute_snapneigh.h"

#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

enum{SCALAR,VECTOR,ARRAY};

ComputeSnapneigh::ComputeSnapneigh(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cutsq(nullptr), list(nullptr), radelem(nullptr), wjelem(nullptr),
  sinnerelem(nullptr), dinnerelem(nullptr)
{

  array_flag = 1;
  //vector_flag = 1;
  extarray = 0;

  double rfac0, rmin0;
  int twojmax, switchflag, bzeroflag, bnormflag, wselfallflag;

  int ntypes = atom->ntypes;
  int nargmin = 6+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal compute snap command");

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;
  bikflag = 0;
  dbirjflag = 0;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  switchinnerflag = 0;
  nelements = 1;

  // process required arguments

  memory->create(radelem,ntypes+1,"snapneigh:radelem"); // offset by 1 to match up with types
  memory->create(wjelem,ntypes+1,"snapneigh:wjelem");
  rcutfac = atof(arg[3]);
  rfac0 = atof(arg[4]);
  twojmax = atoi(arg[5]);
  for (int i = 0; i < ntypes; i++)
    radelem[i+1] = atof(arg[6+i]);
  for (int i = 0; i < ntypes; i++)
    wjelem[i+1] = atof(arg[6+ntypes+i]);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,ntypes+1,ntypes+1,"snapneigh:cutsq");
  for (int i = 1; i <= ntypes; i++) {
    cut = 2.0*radelem[i]*rcutfac;
    //printf("cut: %f\n", cut);
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut*cut;
    for (int j = i+1; j <= ntypes; j++) {
      cut = (radelem[i]+radelem[j])*rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rmin0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      rmin0 = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bzeroflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      bzeroflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      switchflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quadraticflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      quadraticflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"chem") == 0) {
      if (iarg+2+ntypes > narg)
        error->all(FLERR,"Illegal compute snap command");
      chemflag = 1;
      memory->create(map,ntypes+1,"compute_snapneigh:map");
      nelements = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR,arg[iarg+2+i],false,lmp);
        if (jelem < 0 || jelem >= nelements)
          error->all(FLERR,"Illegal compute snap command");
        map[i+1] = jelem;
      }
      iarg += 2+ntypes;
    } else if (strcmp(arg[iarg],"bnormflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      bnormflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"wselfallflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      wselfallflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bikflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      bikflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dbirjflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      dbirjflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchinnerflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute snap command");
      switchinnerflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"sinner") == 0) {
      iarg++;
      if (iarg+ntypes > narg)
        error->all(FLERR,"Illegal compute snap command");
      memory->create(sinnerelem,ntypes+1,"snapneigh:sinnerelem");
      for (int i = 0; i < ntypes; i++)
        sinnerelem[i+1] = utils::numeric(FLERR,arg[iarg+i],false,lmp);
      sinnerflag = 1;
      iarg += ntypes;
    } else if (strcmp(arg[iarg],"dinner") == 0) {
      iarg++;
      if (iarg+ntypes > narg)
        error->all(FLERR,"Illegal compute snap command");
      memory->create(dinnerelem,ntypes+1,"snapneigh:dinnerelem");
      for (int i = 0; i < ntypes; i++)
        dinnerelem[i+1] = utils::numeric(FLERR,arg[iarg+i],false,lmp);
      dinnerflag = 1;
      iarg += ntypes;
    } else error->all(FLERR,"Illegal compute snap command");
  }

  if (switchinnerflag && !(sinnerflag && dinnerflag))
    error->all(FLERR,"Illegal compute snap command: switchinnerflag = 1, missing sinner/dinner keyword");

  if (!switchinnerflag && (sinnerflag || dinnerflag))
    error->all(FLERR,"Illegal compute snap command: switchinnerflag = 0, unexpected sinner/dinner keyword");

  /*
  snaptr = new SNA(lmp, rfac0, twojmax,
                   rmin0, switchflag, bzeroflag,
                   chemflag, bnormflag, wselfallflag,
                   nelements, switchinnerflag);
  */

  //ncoeff = snaptr->ncoeff;
  nperdim = ncoeff;
  if (quadraticflag) nperdim += (ncoeff*(ncoeff+1))/2;
  ndims_force = 3;
  ndims_virial = 6;
  yoffset = nperdim;
  zoffset = 2*nperdim;
  natoms = atom->natoms;
  bik_rows = 1;
  if (bikflag) bik_rows = natoms;
  //size_array_rows = bik_rows+ndims_force*natoms+ndims_virial;
  dbirj_rows = ndims_force*natoms;
  size_array_rows = bik_rows+dbirj_rows+ndims_virial;
  size_array_cols = nperdim*atom->ntypes+1;
  lastcol = size_array_cols-1;

  ndims_peratom = ndims_force;
  size_peratom = ndims_peratom*nperdim*atom->ntypes;

  nmax = 0;

}

/* ---------------------------------------------------------------------- */

ComputeSnapneigh::~ComputeSnapneigh()
{

  memory->destroy(neighs);
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  //delete snaptr;

  if (chemflag) memory->destroy(map);

  if (switchinnerflag) {
    memory->destroy(sinnerelem);
    memory->destroy(dinnerelem);
  }

  if (dbirjflag){
    //printf("dbirj_rows: %d\n", dbirj_rows);
    //memory->destroy(dbirj);
    memory->destroy(nneighs);
    memory->destroy(neighsum);
    memory->destroy(icounter);
    //memory->destroy(dbiri);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSnapneigh::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute snap requires a pair style be defined");

  if (cutmax > force->pair->cutforce){
    //printf("----- cutmax cutforce: %f %f\n", cutmax, force->pair->cutforce);
    error->all(FLERR,"Compute snap cutoff is longer than pairwise cutoff");
  }

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

}


/* ---------------------------------------------------------------------- */

void ComputeSnapneigh::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSnapneigh::compute_array()
{

  if (dbirjflag){
    //printf("----- dbirjflag true.\n");
    get_dbirj_length();
    //printf("----- got dbirj_length\n");
  }

  //printf("----- cutmax cutforce: %f %f\n", cutmax, force->pair->cutforce);
  //else{
  int ntotal = atom->nlocal + atom->nghost;

  invoked_array = update->ntimestep;

  // invoke full neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  const int inum = list->inum;
  const int* const ilist = list->ilist;
  const int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;
  int * const type = atom->type;

  // compute sna derivatives for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double** const x = atom->x;
  const int* const mask = atom->mask;
  //printf("----- inum: %d\n", inum);
  //printf("----- NEIGHMASK: %d\n", NEIGHMASK);
  int ninside;
  int numneigh_sum = 0;
  int dbirj_row_indx;
  int dbiri_indx=0;
  for (int ii = 0; ii < inum; ii++) {
    int irow = 0;
    if (bikflag) irow = atom->tag[ilist[ii] & NEIGHMASK]-1;
    //printf("----- i, itag: %d %d\n", ilist[ii] & NEIGHMASK, atom->tag[ilist[ii]]);
    const int i = ilist[ii];
    //printf("----- ii, i: %d %d\n", ii, i);
    //printf("----- mask[i] groupbit: %d %d\n", mask[i], groupbit);
    if (mask[i] & groupbit) {

      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int itype = type[i];
      int ielem = 0;
      if (chemflag)
        ielem = map[itype];
      const double radi = radelem[itype];
      const int* const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      const int typeoffset_local = ndims_peratom*nperdim*(itype-1);
      const int typeoffset_global = nperdim*(itype-1);

      // insure rij, inside, and typej  are of size jnum

      //snaptr->grow_rij(jnum);

      // rij[][3] = displacements between atom I and those neighbors
      // inside = indices of neighbors of I within cutoff
      // typej = types of neighbors of I within cutoff
      // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

      /*
      This loop assigns quantities in snaptr.
      snaptr is a SNA class instance, see sna.h

      */
      //int ninside = 0;
      ninside=0;
      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;

        const double delx = x[j][0] - xtmp;
        const double dely = x[j][1] - ytmp;
        const double delz = x[j][2] - ztmp;
        const double rsq = delx*delx + dely*dely + delz*delz;
        int jtype = type[j];
        int jelem = 0;
        if (chemflag)
          jelem = map[jtype];
        if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
          if (dbirjflag){
            //dbirj_row_indx = 3*neighsum[atom->tag[j]-1] + 3*icounter[atom->tag[j]-1] ; // THIS IS WRONG, SEE NEXT VAR.
            //dbirj_row_indx = 3*neighsum[atom->tag[i]-1] + 3*(atom->tag[i]-1) + 3*icounter[atom->tag[j]-1]; // 3*tagi is to leave space for dBi/dRi
            dbirj_row_indx = 3*neighsum[atom->tag[j]-1] + 3*icounter[atom->tag[j]-1] + 3*(atom->tag[j]-1); // THIS IS WRONG, SEE NEXT VAR.
            //printf("--- %d %d %d %d\n", dbirj_row_indx, 3*neighsum[atom->tag[i]-1], 3*(atom->tag[i]-1), jj);
            //printf("jtag, icounter, dbirj_row_indx: %d, %d, %d %d %d\n", atom->tag[j], icounter[atom->tag[j]-1], dbirj_row_indx+0, dbirj_row_indx+1, dbirj_row_indx+2);
            icounter[atom->tag[j]-1] += 1;

            neighs[dbirj_row_indx+0][0] = atom->tag[i];
            neighs[dbirj_row_indx+1][0] = atom->tag[i];
            neighs[dbirj_row_indx+2][0] = atom->tag[i];

            neighs[dbirj_row_indx+0][1] = atom->tag[j];
            neighs[dbirj_row_indx+1][1] = atom->tag[j];
            neighs[dbirj_row_indx+2][1] = atom->tag[j];

            neighs[dbirj_row_indx+0][2] = 0;
            neighs[dbirj_row_indx+1][2] = 1;
            neighs[dbirj_row_indx+2][2] = 2;

            dbiri_indx = dbiri_indx+3;
          }
        }
      }

      //printf("--- dbiri_indx: %d\n", dbiri_indx);
      // Put dBi/dRi in

      neighs[dbiri_indx+0][0] = atom->tag[i];
      neighs[dbiri_indx+1][0] = atom->tag[i];
      neighs[dbiri_indx+2][0] = atom->tag[i];

      neighs[dbiri_indx+0][1] = atom->tag[i];
      neighs[dbiri_indx+1][1] = atom->tag[i];
      neighs[dbiri_indx+2][1] = atom->tag[i];

      neighs[dbiri_indx+0][2] = 0;
      neighs[dbiri_indx+1][2] = 1;
      neighs[dbiri_indx+2][2] = 2;

      dbiri_indx = dbiri_indx+3;
    }

    numneigh_sum += ninside;
  } // for (int ii = 0; ii < inum; ii++)


  // Check icounter.
  /*
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      printf("icounter[i]: %d\n", icounter[i]);
    }
  }
  */


  // sum up over all processes
  // I'll need to do something like this...
  /*
  MPI_Allreduce(&snap[0][0],&snapall[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);
  // assign energy to last column
  for (int i = 0; i < bik_rows; i++) snapall[i][lastcol] = 0;
  int irow = 0;
  double reference_energy = c_pe->compute_scalar();
  snapall[irow][lastcol] = reference_energy;
  */

  /*
  for (int i=0; i<dbirj_rows; i++){
    printf("----- %d: %f %f\n", i, array[i][0], array[i][1]);
  }

  //printf("vector[0]: %d\n", vector[0]);
  printf("----- End of compute snapneigh.\n");
  */
}

/* ----------------------------------------------------------------------
   compute dbirj length
------------------------------------------------------------------------- */

void ComputeSnapneigh::get_dbirj_length()
{
  // invoke full neighbor list (will copy or build if necessary)
  neighbor->build_one(list);
  dbirj_rows = 0;
  const int inum = list->inum;
  const int* const ilist = list->ilist;
  const int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;
  int * const type = atom->type;
  const int* const mask = atom->mask;
  double** const x = atom->x;
  //printf("----- inum: %d\n", inum);
  memory->create(neighsum, atom->nlocal, "snapneigh:neighsum");
  memory->create(nneighs, atom->nlocal, "snapneigh:nneighs");
  memory->create(icounter, atom->nlocal, "snapneigh:icounter");
  //memory->create(dbiri, 3*inum,ncoeff, "snapneigh:dbiri");
  /*
  for (int ii=0; ii<inum; ii++){
    for (int icoeff=0; icoeff<ncoeff; icoeff++){
      dbiri[ii][icoeff]=0.0;
    }
  }
  */
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      icounter[i]=0;
      neighsum[i] = 0;
      nneighs[i] = 0;
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const int itype = type[i];
      const int* const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      //printf("-----   jnum: %d\n", jnum);
      int jnum_cutoff = 0;
      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;

        const double delx = x[j][0] - xtmp;
        const double dely = x[j][1] - ytmp;
        const double delz = x[j][2] - ztmp;
        const double rsq = delx*delx + dely*dely + delz*delz;
        int jtype = type[j];

        if (rsq < cutsq[itype][jtype]&&rsq>1e-20) {
          dbirj_rows += 1; //jnum + 1;
          jnum_cutoff += 1;
          nneighs[i]+=1;
        }
      }
      //printf("----- jnum_cutoff: %d\n", jnum_cutoff);
    }
  }

  dbirj_rows *= ndims_force;

  // Loop over all atoms again to calculate neighsum.
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      //printf("nneighs[i]: %d\n", nneighs[i]);
      //neighsum[i] = 0;
      //printf("i nneighs[i]: %d %d\n", i, nneighs[i]);
      if (i==0){
        neighsum[i]=0;
      }
      else{
        for (int jj=0; jj < ii; jj++){
          const int j = ilist[jj];
          if (mask[j] & groupbit) {
            //printf("  j nneighs[j-1]: %d %d\n", j, nneighs[j]);
            neighsum[i] += nneighs[j];
          }
        }
        //neighsum[i] += 1; // Add 1 for the self term dBi/dRi
      }
    }
    //printf("%d\n", neighsum[i]);
  }

  size_array_rows = dbirj_rows+(3*atom->nlocal);
  size_array_cols = 3;
  //memory->create(dbirj, dbirj_rows, ncoeff, "snapneigh:dbirj");
  memory->create(neighs, size_array_rows, size_array_cols, "snapneigh:neighs");

  //vector = neighs;
  array = neighs;
  // Set size array rows which now depends on dbirj_rows.
  //size_array_rows = bik_rows+dbirj_rows+ndims_virial;
  //printf("----- dbirj_rows: %d\n", dbirj_rows);
  //printf("----- end of dbirj length.\n");

}

/* ----------------------------------------------------------------------
   compute array length
------------------------------------------------------------------------- */

double ComputeSnapneigh::compute_scalar()
{
  if (dbirjflag) get_dbirj_length();
  return size_array_rows;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSnapneigh::memory_usage()
{

  double bytes = (double)size_array_rows*size_array_cols *
    sizeof(double);                                     // array

  return bytes;
}
