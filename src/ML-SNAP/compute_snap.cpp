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

#include "compute_snap.h"

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

ComputeSnap::ComputeSnap(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), cutsq(nullptr), list(nullptr), snap(nullptr),
  snapall(nullptr), snap_peratom(nullptr), radelem(nullptr), wjelem(nullptr),
  sinnerelem(nullptr), dinnerelem(nullptr), snaptr(nullptr)
{

  array_flag = 1;
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

  memory->create(radelem,ntypes+1,"snap:radelem"); // offset by 1 to match up with types
  memory->create(wjelem,ntypes+1,"snap:wjelem");
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
  memory->create(cutsq,ntypes+1,ntypes+1,"snap:cutsq");
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
      memory->create(map,ntypes+1,"compute_snap:map");
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
      memory->create(sinnerelem,ntypes+1,"snap:sinnerelem");
      for (int i = 0; i < ntypes; i++)
        sinnerelem[i+1] = utils::numeric(FLERR,arg[iarg+i],false,lmp);
      sinnerflag = 1;
      iarg += ntypes;
    } else if (strcmp(arg[iarg],"dinner") == 0) {
      iarg++;
      if (iarg+ntypes > narg)
        error->all(FLERR,"Illegal compute snap command");
      memory->create(dinnerelem,ntypes+1,"snap:dinnerelem");
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

  snaptr = new SNA(lmp, rfac0, twojmax,
                   rmin0, switchflag, bzeroflag,
                   chemflag, bnormflag, wselfallflag,
                   nelements, switchinnerflag);

  ncoeff = snaptr->ncoeff;
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
  if (dbirjflag) size_array_cols = nperdim;
  else size_array_cols = nperdim*atom->ntypes+1;
  lastcol = size_array_cols-1;

  ndims_peratom = ndims_force;
  size_peratom = ndims_peratom*nperdim*atom->ntypes;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeSnap::~ComputeSnap()
{

  memory->destroy(snap);
  memory->destroy(snapall);
  memory->destroy(snap_peratom);
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;

  if (chemflag) memory->destroy(map);

  if (switchinnerflag) {
    memory->destroy(sinnerelem);
    memory->destroy(dinnerelem);
  }

  if (dbirjflag){
    //printf("dbirj_rows: %d\n", dbirj_rows);
    //printf("----- destroy dbirj\n");
    memory->destroy(dbirj);
    //printf("----- 1-1-1-1-1-1\n");
    memory->destroy(nneighs);
    //printf("----- 2-1-1-1-1-1\n");
    memory->destroy(neighsum);
    //printf("----- 3-1-1-1-1-1\n");
    memory->destroy(icounter);
    //printf("----- 4-1-1-1-1-1\n");
    memory->destroy(dbiri);
    //printf("----- 5-1-1-1-1-1\n");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute snap requires a pair style be defined");

  if (cutmax > force->pair->cutforce){
    //printf("----- cutmax cutforce: %f %f\n", cutmax, force->pair->cutforce);
    error->all(FLERR,"Compute snap cutoff is longer than pairwise cutoff");
  }

  // need an occasional full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("snap").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute snap");
  snaptr->init();

  // allocate memory for global array

  //printf("----- dbirjflag: %d\n", dbirjflag);
  memory->create(snap,size_array_rows,size_array_cols,
                 "snap:snap");
  memory->create(snapall,size_array_rows,size_array_cols,
                 "snap:snapall");
  array = snapall;

  // find compute for reference energy

  std::string id_pe = std::string("thermo_pe");
  int ipe = modify->find_compute(id_pe);
  if (ipe == -1)
    error->all(FLERR,"compute thermo_pe does not exist.");
  c_pe = modify->compute[ipe];

  // add compute for reference virial tensor
  std::string id_virial = std::string("snap_press");
  std::string pcmd = id_virial + " all pressure NULL virial";
  modify->add_compute(pcmd);

  int ivirial = modify->find_compute(id_virial);
  if (ivirial == -1)
    error->all(FLERR,"compute snap_press does not exist.");
  c_virial = modify->compute[ivirial];

}


/* ---------------------------------------------------------------------- */

void ComputeSnap::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSnap::compute_array()
{

  //printf(" -----2 dbirjflag: %d\n", dbirjflag);
  if (dbirjflag){
    //printf("----- dbirjflag true.\n");
    get_dbirj_length();
    //printf("----- got dbirj_length\n");
  }
  //else{
  //  printf("----- dbirjflag false.\n");
  //}

  //printf("----- cutmax cutforce: %f %f\n", cutmax, force->pair->cutforce);
  //else{
  int ntotal = atom->nlocal + atom->nghost;

  invoked_array = update->ntimestep;

  // grow snap_peratom array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(snap_peratom);
    nmax = atom->nmax;
    memory->create(snap_peratom,nmax,size_peratom,
                   "snap:snap_peratom");
  }
  // clear global array
  //printf("size_array_rows: %d\n", size_array_rows);
  for (int irow = 0; irow < size_array_rows; irow++){
    for (int icoeff = 0; icoeff < size_array_cols; icoeff++){
      //printf("%d %d\n", irow, icoeff);
      snap[irow][icoeff] = 0.0;
    }
  }

  // clear local peratom array
  for (int i = 0; i < ntotal; i++)
    for (int icoeff = 0; icoeff < size_peratom; icoeff++) {
      snap_peratom[i][icoeff] = 0.0;
    }

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

      snaptr->grow_rij(jnum);

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
          //printf("cutsq: %f\n", cutsq[itype][jtype]);
          snaptr->rij[ninside][0] = delx;
          snaptr->rij[ninside][1] = dely;
          snaptr->rij[ninside][2] = delz;
          snaptr->inside[ninside] = j;
          snaptr->wj[ninside] = wjelem[jtype];
          snaptr->rcutij[ninside] = (radi+radelem[jtype])*rcutfac;
          if (switchinnerflag) {
            snaptr->sinnerij[ninside] = 0.5*(sinnerelem[itype]+sinnerelem[jtype]);
            snaptr->dinnerij[ninside] = 0.5*(dinnerelem[itype]+dinnerelem[jtype]);
          }
          if (chemflag) snaptr->element[ninside] = jelem;
          ninside++;
        }
      }

      /*
      Now that we have assigned neighbor quantities with previous loop, we are ready to compute things.
      Here we compute the wigner functions (U), Z is some other quantity, and bi is bispectrum.
      */
      snaptr->compute_ui(ninside, ielem);
      snaptr->compute_zi();
      snaptr->compute_bi(ielem);

      /*
      Looks like this loop computes derivatives.
      How does snaptr know what atom I we're dealing with?
      I think it only needs neighbor info, and then it goes from there.
      */
      //printf("----- Derivative loop - looping over neighbors j.\n");
      //printf("-----   ninside: %d\n", ninside); // numneighs of I within cutoff
      for (int jj = 0; jj < ninside; jj++) {
        //printf("-----   jj: %d\n", jj);
        const int j = snaptr->inside[jj];
        //printf("-----   jj, j, jtag: %d %d %d\n", jj, j, atom->tag[j]);
        //int dbirj_row_indx = 3*neighsum[i] + 3*jj ; // THIS IS WRONG, SEE NEXT LINE.
        //int dbirj_row_indx = 3*neighsum[j] + 3*i_indx; // need to get i_indx.
        // How to get i_indx?
        /*
        i_indx must start at zero and end at (nneighs[j]-1).
        We can start a counter for each atom j.
        Maybe this icounter can serve as an index for i as a neighbor of j.
        icounter starts at zero and ends at (nneighs[j]-1).
        */
        //icounter[atom->tag[j]-1] += 1;
        if (dbirjflag){
          dbirj_row_indx = 3*neighsum[atom->tag[j]-1] + 3*icounter[atom->tag[j]-1] ; // THIS IS WRONG, SEE NEXT VAR.
          //printf("jtag, icounter, dbirj_row_indx: %d, %d, %d %d %d\n", atom->tag[j], icounter[atom->tag[j]-1], dbirj_row_indx+0, dbirj_row_indx+1, dbirj_row_indx+2);
          icounter[atom->tag[j]-1] += 1;
        }
        //int dbirj_row_indx = 3*neighsum[atom->tag[j]-1] + 3*icounter[atom->tag[j]-1] ; // THIS IS WRONG, SEE NEXT VAR.
        //printf("jtag, icounter, dbirj_row_indx: %d, %d, %d %d %d\n", atom->tag[j], icounter[atom->tag[j]-1], dbirj_row_indx+0, dbirj_row_indx+1, dbirj_row_indx+2);
        //icounter[atom->tag[j]-1] += 1;
        /*
        j is an atom index starting from 0.
        Use atom->tag[j] to get the atom in the box (index starts at 1).
        Need to make sure that the order of these ij pairs is the same when multiplying by dE/dD later.
        */
        //printf("-----   jj, j, jtag: %d %d %d\n", jj, j, atom->tag[j]);
        snaptr->compute_duidrj(jj);
        snaptr->compute_dbidrj();

        // Accumulate dBi/dRi, -dBi/dRj

        /*
        snap_peratom[i] has type double * because each atom index has indices for descriptors.
        */
        double *snadi = snap_peratom[i]+typeoffset_local;
        double *snadj = snap_peratom[j]+typeoffset_local;
        //printf("-----   ncoeff: %d\n", ncoeff);
        //printf("snadi: %f %f %f %f %f\n", snadi[0], snadi[1], snadi[2], snadi[3], snadi[4]);
        //printf("----- typeoffset_local: %d\n", typeoffset_local);
        //printf("snadi: ");
        //for (int s=0; s<(ncoeff*3); s++){
        //  printf("%f ", snadi[s]);
        //}
        /*
        printf("snadj: ");
        for (int s=0; s<(ncoeff*3); s++){
          printf("%f ", snadj[s]);
        }
        */
        //printf("\n");
        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          //printf("-----     dblist[icoeff]: %f %f %f\n", snaptr->dblist[icoeff][0], snaptr->dblist[icoeff][1], snaptr->dblist[icoeff][2]);
          /*
          I think these are the descriptor derivatives.
          Desriptor derivatives wrt atom i.
          What exactly is being summed here?
          This is a loop over descriptors or coeff k.

          */
          snadi[icoeff] += snaptr->dblist[icoeff][0];
          snadi[icoeff+yoffset] += snaptr->dblist[icoeff][1];
          snadi[icoeff+zoffset] += snaptr->dblist[icoeff][2];
          /*
          Descriptor derivatives wrt atom j
          */
          snadj[icoeff] -= snaptr->dblist[icoeff][0];
          snadj[icoeff+yoffset] -= snaptr->dblist[icoeff][1];
          snadj[icoeff+zoffset] -= snaptr->dblist[icoeff][2];


          if (dbirjflag){
            dbirj[dbirj_row_indx+0][icoeff] = snaptr->dblist[icoeff][0];
            dbirj[dbirj_row_indx+1][icoeff] = snaptr->dblist[icoeff][1];
            dbirj[dbirj_row_indx+2][icoeff] = snaptr->dblist[icoeff][2];
            // Accumulate dBi/dRi = sum (-dBi/dRj) for neighbors j of if i.
            dbiri[3*(atom->tag[i]-1)+0][icoeff] -= snaptr->dblist[icoeff][0];
            dbiri[3*(atom->tag[i]-1)+1][icoeff] -= snaptr->dblist[icoeff][1];
            dbiri[3*(atom->tag[i]-1)+2][icoeff] -= snaptr->dblist[icoeff][2];
          }


        }

        if (quadraticflag) {
          const int quadraticoffset = ncoeff;
          snadi += quadraticoffset;
          snadj += quadraticoffset;
          int ncount = 0;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bi = snaptr->blist[icoeff];
            double bix = snaptr->dblist[icoeff][0];
            double biy = snaptr->dblist[icoeff][1];
            double biz = snaptr->dblist[icoeff][2];

            // diagonal elements of quadratic matrix

            double dbxtmp = bi*bix;
            double dbytmp = bi*biy;
            double dbztmp = bi*biz;

            snadi[ncount] +=         dbxtmp;
            snadi[ncount+yoffset] += dbytmp;
            snadi[ncount+zoffset] += dbztmp;
            snadj[ncount] -=         dbxtmp;
            snadj[ncount+yoffset] -= dbytmp;
            snadj[ncount+zoffset] -= dbztmp;

            ncount++;

            // upper-triangular elements of quadratic matrix

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
              double dbxtmp = bi*snaptr->dblist[jcoeff][0]
                + bix*snaptr->blist[jcoeff];
              double dbytmp = bi*snaptr->dblist[jcoeff][1]
                + biy*snaptr->blist[jcoeff];
              double dbztmp = bi*snaptr->dblist[jcoeff][2]
                + biz*snaptr->blist[jcoeff];

              snadi[ncount] +=         dbxtmp;
              snadi[ncount+yoffset] += dbytmp;
              snadi[ncount+zoffset] += dbztmp;
              snadj[ncount] -=         dbxtmp;
              snadj[ncount+yoffset] -= dbytmp;
              snadj[ncount+zoffset] -= dbztmp;

              ncount++;
            }
          }

        }
      } // for (int jj = 0; jj < ninside; jj++)
      //printf("---- irow after jj loop: %d\n", irow);

      // Accumulate Bi
      //printf("----- Accumulate Bi.\n");

      // linear contributions
      int k;
      if (dbirjflag) k = 0;
      else k = typeoffset_global;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++){
        //printf("----- %d %f\n", icoeff, snaptr->blist[icoeff]);
        snap[irow][k++] += snaptr->blist[icoeff];
      }

      // quadratic contributions

      if (quadraticflag) {
        for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
          double bveci = snaptr->blist[icoeff];
          snap[irow][k++] += 0.5*bveci*bveci;
          for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            double bvecj = snaptr->blist[jcoeff];
            snap[irow][k++] += bveci*bvecj;
          }
        }
      }
    }

    numneigh_sum += ninside;
  } // for (int ii = 0; ii < inum; ii++)

  //printf("----- bik_rows: %d\n", bik_rows);
  //printf("----- bikflag: %d\n", bikflag);

  // Check icounter.
  /*
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      printf("icounter[i]: %d\n", icounter[i]);
    }
  }
  */

  // Sum all the derivatives we calculated to check usual compute snap output.
  /*
  if (dbirjflag){
    fh_d = fopen("DEBUG", "w");
    int row_indx=0;
    double sum;
    for (int ii=0; ii<inum; ii++){
      printf("ii: %d\n", ii);
      for (int a=0; a<3; a++){
        for (int k=0; k<ncoeff; k++){
          //double sumx = 0.0;
          //double sumy = 0.0;
          //double sumz = 0.0;
          sum=0.0;
          for (int jj=0; jj<nneighs[ii]; jj++){
            row_indx=3*neighsum[ii]+3*jj;
            //sumx+=dbirj[row_indx+a][k];
            //sumy+=dbirj[row_indx+1][k];
            //sumz+=dbirj[row_indx+2][k];
            sum+=dbirj[row_indx+a][k];
          }
          // Add dBi/dRi
          //sumx+=dbiri[3*ii+0][k];
          //sumy+=dbiri[3*ii+1][k];
          //sumz+=dbiri[3*ii+2][k];

          sum+=dbiri[3*ii+a][k];

          //printf("%d %f %f %f\n", ii, sumx, sumy, sumz);
          fprintf(fh_d, "%f ", sum);
        }
        fprintf(fh_d, "\n");

      }
    }
    fclose(fh_d);
  }
  */

  if (dbirjflag){

    //printf(" Accumulate to global array.\n");
    //printf("ntypes: %d\n", atom->ntypes);
    //for (int itype = 0; itype < atom->ntypes; itype++) {
    int dbiri_indx;
    int irow;
    for (int itype=0; itype<1; itype++){
      const int typeoffset_local = ndims_peratom*nperdim*itype;
      const int typeoffset_global = nperdim*itype;
      //printf("----- nperdim: %d\n", nperdim);
      for (int icoeff = 0; icoeff < nperdim; icoeff++) {
        //printf("-----   icoeff: %d\n", icoeff);
        dbiri_indx=0;
        for (int i = 0; i < atom->nlocal; i++) {
            //printf("i: %d\n", i);

            //int dbiri_indx = 0;
            //int irow;
            for (int jj=0; jj<nneighs[i]; jj++){
              //printf("  jj: %d\n", jj);
              int dbirj_row_indx = 3*neighsum[atom->tag[i]-1] + 3*jj;
              int snap_row_indx = 3*neighsum[atom->tag[i]-1] + 3*(atom->tag[i]-1) + 3*jj;
              //printf("snap_row_indx: %d\n", snap_row_indx);
              //int irow = dbirj_row_indx+bik_rows;
              irow = snap_row_indx + bik_rows;
              //printf("  row_indx, irow: %d %d\n", dbirj_row_indx, irow);
              snap[irow++][icoeff+typeoffset_global] += dbirj[dbirj_row_indx+0][icoeff];
              //printf("    irow: %d\n", irow);
              snap[irow++][icoeff+typeoffset_global] += dbirj[dbirj_row_indx+1][icoeff];
              //printf("    irow: %d\n", irow);
              snap[irow][icoeff+typeoffset_global]   += dbirj[dbirj_row_indx+2][icoeff];
              dbiri_indx = dbiri_indx+3;
            }
            // Put dBi/dRi at end of each dBj/dRi chunk.
            //int dbiri_row_indx;
            irow = dbiri_indx + bik_rows;
            //printf("dbiri_indx: %d\n", dbiri_indx);
            //printf("dbiri: %f %f %f\n", dbiri[3*i+0][icoeff], dbiri[3*i+1][icoeff], dbiri[3*i+2][icoeff]);
            snap[irow++][icoeff+typeoffset_global] += dbiri[3*i+0][icoeff];
            //printf("    irow: %d\n", irow);
            snap[irow++][icoeff+typeoffset_global] += dbiri[3*i+1][icoeff];
            //printf("    irow: %d\n", irow);
            snap[irow][icoeff+typeoffset_global]   += dbiri[3*i+2][icoeff];

            dbiri_indx = dbiri_indx+3;
        }
      }
    }
    //printf(" Accumulated to global array.\n");

  }

  else{
    //printf("----- Accumulate bispecturm force contributions to global array.\n");
    // accumulate bispectrum force contributions to global array
    //printf("----- ntotal, nmax, natoms: %d %d %d\n", ntotal, nmax, atom->natoms);
    for (int itype = 0; itype < atom->ntypes; itype++) {
      const int typeoffset_local = ndims_peratom*nperdim*itype;
      const int typeoffset_global = nperdim*itype;
      //printf("----- nperdim: %d\n", nperdim);
      /*nperdim = ncoeff set previsouly*/
      for (int icoeff = 0; icoeff < nperdim; icoeff++) {
        //printf("-----   icoeff: %d\n", icoeff);
        for (int i = 0; i < ntotal; i++) {
          double *snadi = snap_peratom[i]+typeoffset_local;
          int iglobal = atom->tag[i];
          if (icoeff==4){
            if ( (snadi[icoeff] != 0.0) || (snadi[icoeff+yoffset] != 0.0) || (snadi[icoeff+zoffset] != 0.0) ){
              //printf("%d %d %f %f %f\n", i, iglobal, snadi[icoeff], snadi[icoeff+yoffset], snadi[icoeff+zoffset]);
            }
          }
          int irow = 3*(iglobal-1)+bik_rows;
          //printf("-----     snadi[icoeff]: %f\n", snadi[icoeff]);
          snap[irow++][icoeff+typeoffset_global] += snadi[icoeff];
          snap[irow++][icoeff+typeoffset_global] += snadi[icoeff+yoffset];
          snap[irow][icoeff+typeoffset_global] += snadi[icoeff+zoffset];
        }
      }
    }
  }

  //printf("----- Accumulate forces to global array.\n");
  /*
  These are the last columns.
  */
 // accumulate forces to global array
 if (dbirjflag){

 }
 else{
    for (int i = 0; i < atom->nlocal; i++) {
      int iglobal = atom->tag[i];
      int irow = 3*(iglobal-1)+bik_rows;
      //printf("---- irow: %d\n", irow);
      snap[irow++][lastcol] = atom->f[i][0];
      //printf("---- irow: %d\n", irow);
      snap[irow++][lastcol] = atom->f[i][1];
      //printf("---- irow: %d\n", irow);
      snap[irow][lastcol] = atom->f[i][2];
    }
  }

  // accumulate bispectrum virial contributions to global array

  //dbdotr_compute();

  // sum up over all processes
  MPI_Allreduce(&snap[0][0],&snapall[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);
  // assign energy to last column
  if (dbirjflag){

  }
  else{
    for (int i = 0; i < bik_rows; i++) snapall[i][lastcol] = 0;
    int irow = 0;
    double reference_energy = c_pe->compute_scalar();
    snapall[irow][lastcol] = reference_energy;
  }

  // assign virial stress to last column
  // switch to Voigt notation

  c_virial->compute_vector();
  /*
  irow += 3*natoms+bik_rows;
  snapall[irow++][lastcol] = c_virial->vector[0];
  snapall[irow++][lastcol] = c_virial->vector[1];
  snapall[irow++][lastcol] = c_virial->vector[2];
  snapall[irow++][lastcol] = c_virial->vector[5];
  snapall[irow++][lastcol] = c_virial->vector[4];
  snapall[irow][lastcol] = c_virial->vector[3];
  */

  //}// else
}

/* ----------------------------------------------------------------------
   compute global virial contributions via summing r_i.dB^j/dr_i over
   own & ghost atoms
------------------------------------------------------------------------- */

void ComputeSnap::dbdotr_compute()
{
  double **x = atom->x;

  int irow0;
  if (dbirjflag){
    irow0 = bik_rows+dbirj_rows+(3*natoms);
  }
  else{
    irow0 = bik_rows+ndims_force*natoms;
  }
  //int irow0 = bik_rows+ndims_force*natoms;

  // sum over bispectrum contributions to forces
  // on all particles including ghosts

  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++)
    for (int itype = 0; itype < atom->ntypes; itype++) {
      const int typeoffset_local = ndims_peratom*nperdim*itype;
      const int typeoffset_global = nperdim*itype;
      double *snadi = snap_peratom[i]+typeoffset_local;
      for (int icoeff = 0; icoeff < nperdim; icoeff++) {
        double dbdx = snadi[icoeff];
        double dbdy = snadi[icoeff+yoffset];
        double dbdz = snadi[icoeff+zoffset];
        int irow = irow0;
        snap[irow++][icoeff+typeoffset_global] += dbdx*x[i][0];
        snap[irow++][icoeff+typeoffset_global] += dbdy*x[i][1];
        snap[irow++][icoeff+typeoffset_global] += dbdz*x[i][2];
        snap[irow++][icoeff+typeoffset_global] += dbdz*x[i][1];
        snap[irow++][icoeff+typeoffset_global] += dbdz*x[i][0];
        snap[irow][icoeff+typeoffset_global] += dbdy*x[i][0];
      }
    }
}

/* ----------------------------------------------------------------------
   compute dbirj length
------------------------------------------------------------------------- */

void ComputeSnap::get_dbirj_length()
{

  memory->destroy(snap);
  memory->destroy(snapall);
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
  memory->create(neighsum, inum, "snap:neighsum");
  memory->create(nneighs, inum, "snap:nneighs");
  memory->create(icounter, inum, "snap:icounter");
  memory->create(dbiri, 3*atom->nlocal,ncoeff, "snap:dbiri");
  for (int ii=0; ii<3*atom->nlocal; ii++){
    for (int icoeff=0; icoeff<ncoeff; icoeff++){
      dbiri[ii][icoeff]=0.0;
    }
  }
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

  //printf("----- dbirj_rows: %d\n", dbirj_rows);

  memory->create(dbirj, dbirj_rows, ncoeff, "snap:dbirj");
  for (int i=0; i<dbirj_rows; i++){
    for (int j=0; j<ncoeff; j++){
      dbirj[i][j]=0.0;
    }
  }
  // Set size array rows which now depends on dbirj_rows.
  //size_array_rows = bik_rows+dbirj_rows+ndims_virial;
  size_array_rows = bik_rows+dbirj_rows+ndims_virial+3*atom->nlocal; // Add 3*N for dBi/dRi
  //printf("----- dbirj_rows: %d\n", dbirj_rows);
  //printf("----- end of dbirj length.\n");

  memory->create(snap,size_array_rows,size_array_cols, "snap:snap");
  memory->create(snapall,size_array_rows,size_array_cols, "snap:snapall");
  array = snapall;

}

/* ----------------------------------------------------------------------
   compute array length
------------------------------------------------------------------------- */

double ComputeSnap::compute_scalar()
{
  if (dbirjflag) get_dbirj_length();
  return size_array_rows;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeSnap::memory_usage()
{

  double bytes = (double)size_array_rows*size_array_cols *
    sizeof(double);                                     // snap
  bytes += (double)size_array_rows*size_array_cols *
    sizeof(double);                                     // snapall
  bytes += (double)nmax*size_peratom * sizeof(double);  // snap_peratom
  bytes += snaptr->memory_usage();                      // SNA object
  int n = atom->ntypes+1;
  bytes += (double)n*sizeof(int);        // map

  return bytes;
}
