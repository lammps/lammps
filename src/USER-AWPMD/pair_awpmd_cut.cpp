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
   Contributing author: Ilya Valuev (JIHT, Moscow, Russia)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_awpmd_cut.h"
#include "atom.h"
#include "update.h"
#include "min.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

#include "TCP/wpmd_split.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairAWPMDCut::PairAWPMDCut(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  nmax = 0;
  min_var = NULL;
  min_varforce = NULL;
  nextra = 4;
  pvector = new double[nextra];

  ermscale=1.;
  width_pbc=0.;
  wpmd= new AWPMD_split();

  half_box_length=0;
}

/* ---------------------------------------------------------------------- */

PairAWPMDCut::~PairAWPMDCut()
{
  delete [] pvector;
  memory->destroy(min_var);
  memory->destroy(min_varforce);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }

  delete wpmd;
}


struct cmp_x{
  double **xx;
  double tol;
  cmp_x(double **xx_=NULL, double tol_=1e-12):xx(xx_),tol(tol_){}
  bool operator()(const pair<int,int> &left, const pair<int,int> &right) const {
    if(left.first==right.first){
      double d=xx[left.second][0]-xx[right.second][0];
      if(d<-tol)
        return true;
      else if(d>tol)
        return false;
      d=xx[left.second][1]-xx[right.second][1];
      if(d<-tol)
        return true;
      else if(d>tol)
        return false;
      d=xx[left.second][2]-xx[right.second][2];
      if(d<-tol)
        return true;
      else
        return false;
    }
    else
      return left.first<right.first;
  }
};

/* ---------------------------------------------------------------------- */

void PairAWPMDCut::compute(int eflag, int vflag)
{

  // pvector = [KE, Pauli, ecoul, radial_restraint]
  for (int i=0; i<4; i++) pvector[i] = 0.0;

  if (eflag || vflag)
    ev_setup(eflag,vflag);
  else
    evflag = vflag_fdotr = 0; //??

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *spin = atom->spin;
  int *type = atom->type;
  int *etag = atom->etag;
  double **v = atom->v;

  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntot=nlocal+nghost;

  int newton_pair = force->newton_pair;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;




  // width pbc
  if(width_pbc<0)
    wpmd->Lextra=2*half_box_length;
  else
    wpmd->Lextra=width_pbc;

  wpmd->newton_pair=newton_pair;



# if 1
  // mapping of the LAMMPS numbers to the AWPMC numbers
  vector<int> gmap(ntot,-1);

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    // local particles are all there
    gmap[i]=0;
    Vector_3 ri=Vector_3(x[i][0],x[i][1],x[i][2]);
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      if(j>=nlocal){ // this is a ghost
        Vector_3 rj=Vector_3(x[j][0],x[j][1],x[j][2]);
        int jtype = type[j];
        double rsq=(ri-rj).norm2();
        if (rsq < cutsq[itype][jtype])
          gmap[j]=0; //bingo, this ghost is really needed

      }
    }
  }

# else  // old mapping
  // mapping of the LAMMPS numbers to the AWPMC numbers
  vector<int> gmap(ntot,-1);
  // map for filtering the clones out: [tag,image] -> id
  typedef  map< pair<int,int>, int, cmp_x >  map_t;
  cmp_x cmp(x);
  map_t idmap(cmp);
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    // local particles are all there
    idmap[make_pair(atom->tag[i],i)]=i;
    bool i_local= i<nlocal ? true : false;
    if(i_local)
      gmap[i]=0;
    else if(gmap[i]==0) // this is a ghost which already has been tested
      continue;
    Vector_3 ri=Vector_3(x[i][0],x[i][1],x[i][2]);
    int itype = type[i];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      pair<map_t::iterator,bool> res=idmap.insert(make_pair(make_pair(atom->tag[j],j),j));
      bool have_it=!res.second;
      if(have_it){ // the clone of this particle is already listed
        if(res.first->second!=j) // check that was not the very same particle
          gmap[j]=-1; // filter out
        continue;
      }

      bool j_local= j<nlocal ? true : false;
      if((i_local && !j_local) || (j_local && !i_local)){ // some of them is a ghost
        Vector_3 rj=Vector_3(x[j][0],x[j][1],x[j][2]);
        int jtype = type[j];
        double rsq=(ri-rj).norm2();
        if (rsq < cutsq[itype][jtype]){
          if(!i_local){
            gmap[i]=0; //bingo, this ghost is really needed
            break; // don't need to continue j loop
          }
          else
            gmap[j]=0; //bingo, this ghost is really needed
        }
      }
    }
  }
# endif
  // prepare the solver object
  wpmd->reset();

  map<int,vector<int> > etmap;
  // add particles to the AWPMD solver object
  for (int i = 0; i < ntot; i++) {
    //int i = ilist[ii];
    if(gmap[i]<0) // this particle was filtered out
      continue;
    if(spin[i]==0)  // this is an ion
      gmap[i]=wpmd->add_ion(q[i], Vector_3(x[i][0],x[i][1],x[i][2]),i<nlocal ? atom->tag[i] : -atom->tag[i]);
    else if(spin[i]==1 || spin[i]==-1){ // electron, sort them according to the tag
      etmap[etag[i]].push_back(i);
    }
    else
      error->all(FLERR,fmt("Invalid spin value (%d) for particle %d !",spin[i],i));
  }
  // ion force vector
  Vector_3 *fi=NULL;
  if(wpmd->ni)
    fi= new Vector_3[wpmd->ni];

  // adding electrons
  for(map<int,vector<int> >::iterator it=etmap.begin(); it!= etmap.end(); ++it){
    vector<int> &el=it->second;
    if(!el.size()) // should not happen
      continue;
    int s=spin[el[0]] >0 ? 0 : 1;
    wpmd->add_electron(s); // starts adding the spits
    for(size_t k=0;k<el.size();k++){
      int i=el[k];
      if(spin[el[0]]!=spin[i])
        error->all(FLERR,fmt("WP splits for one electron should have the same spin (at particles %d, %d)!",el[0],i));
      double m= atom->mass ? atom->mass[type[i]] : force->e_mass;
      Vector_3 xx=Vector_3(x[i][0],x[i][1],x[i][2]);
      Vector_3 rv=m*Vector_3(v[i][0],v[i][1],v[i][2]);
      double pv=ermscale*m*atom->ervel[i];
      Vector_2 cc=Vector_2(atom->cs[2*i],atom->cs[2*i+1]);
      gmap[i]=wpmd->add_split(xx,rv,atom->eradius[i],pv,cc,1.,atom->q[i],i<nlocal ? atom->tag[i] : -atom->tag[i]);
      // resetting for the case constraints were applied
      v[i][0]=rv[0]/m;
      v[i][1]=rv[1]/m;
      v[i][2]=rv[2]/m;
      atom->ervel[i]=pv/(m*ermscale);
    }
  }
  wpmd->set_pbc(NULL); // not required for LAMMPS
  wpmd->interaction(0x1|0x4|0x10,fi);

   // get forces from the AWPMD solver object
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if(gmap[i]<0) // this particle was filtered out
      continue;
    if (spin[i]==0) {  // this is an ion, copying forces
      int ion=gmap[i];
      f[i][0]=fi[ion][0];
      f[i][0]=fi[ion][1];
      f[i][0]=fi[ion][2];
    } else { // electron
      int iel=gmap[i];
      int s=spin[i] >0 ? 0 : 1;
      wpmd->get_wp_force(s,iel,(Vector_3 *)f[i],(Vector_3 *)(atom->vforce+3*i),atom->erforce+i,atom->ervelforce+i,(Vector_2 *)(atom->csforce+2*i));
    }
  }

  if(fi)
    delete [] fi;

  // update LAMMPS energy
  if (eflag_either) {
    if (eflag_global){
      eng_coul+= wpmd->get_energy();
      // pvector = [KE, Pauli, ecoul, radial_restraint]
      pvector[0] = wpmd->Ee[0]+wpmd->Ee[1];
      pvector[2] = wpmd->Eii+wpmd->Eei[0]+wpmd->Eei[1]+wpmd->Eee;
      pvector[1] = pvector[0] + pvector[2] - wpmd->Edk - wpmd->Edc - wpmd->Eii;  // All except diagonal terms
      pvector[3] = wpmd->Ew;
    }

    if (eflag_atom) {
      // transfer per-atom energies here
      for (int i = 0; i < ntot; i++) {
        if (gmap[i]<0) // this particle was filtered out
          continue;
        if (spin[i]==0) {
          eatom[i]=wpmd->Eiep[gmap[i]]+wpmd->Eiip[gmap[i]];
        } else {
          int s=spin[i] >0 ? 0 : 1;
          eatom[i]=wpmd->Eep[s][gmap[i]]+wpmd->Eeip[s][gmap[i]]+wpmd->Eeep[s][gmap[i]]+wpmd->Ewp[s][gmap[i]];
        }
      }
    }
  }
  if (vflag_fdotr) {
    virial_fdotr_compute();
    if (flexible_pressure_flag)
       virial_eradius_compute();
  }
}

/* ----------------------------------------------------------------------
   electron width-specific contribution to global virial
------------------------------------------------------------------------- */

void PairAWPMDCut::virial_eradius_compute()
{
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  double e_virial;
  int *spin = atom->spin;

  // sum over force on all particles including ghosts

  if (neighbor->includegroup == 0) {
    int nall = atom->nlocal + atom->nghost;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }

  // neighbor includegroup flag is set
  // sum over force on initial nfirst particles and ghosts

  } else {
    int nall = atom->nfirst;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }

    nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }
  }
}



/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairAWPMDCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ---------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
// the format is: pair_style awpmd/cut [<global_cutoff|-1> [command1] [command2] ...]
// commands:
// [hartree|dproduct|uhf]  -- quantum approximation level (default is hartree)
// [free|pbc <length|-1>|fix <w0|-1>|relax|harm <w0>] -- width restriction (default is free)
// [ermscale <number>]  -- scaling factor between electron mass and effective width mass (used for equations of motion only) (default is 1)
// [flex_press]  -- set flexible pressure flag
// -1 for length means default setting (L/2 for cutoff and L for width PBC)

void PairAWPMDCut::settings(int narg, char **arg){
  if (narg < 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  ermscale=1.;
  width_pbc=0.;

  for(int i=1;i<narg;i++){
    // reading commands
    if(!strcmp(arg[i],"hartree"))
      wpmd->approx=AWPMD::HARTREE;
    else if(!strcmp(arg[i],"dproduct"))
      wpmd->approx=AWPMD::DPRODUCT;
    else if(!strcmp(arg[i],"uhf"))
      wpmd->approx=AWPMD::UHF;
    else if(!strcmp(arg[i],"free"))
      wpmd->constraint=AWPMD::NONE;
    else if(!strcmp(arg[i],"fix")){
      wpmd->constraint=AWPMD::FIX;
      i++;
      if(i>=narg)
        error->all(FLERR,"Setting 'fix' should be followed by a number in awpmd/cut");
      wpmd->w0=force->numeric(FLERR,arg[i]);
    }
    else if(!strcmp(arg[i],"harm")){
      wpmd->constraint=AWPMD::HARM;
      i++;
      if(i>=narg)
        error->all(FLERR,"Setting 'harm' should be followed by a number in awpmd/cut");
      wpmd->w0=force->numeric(FLERR,arg[i]);
      wpmd->set_harm_constr(wpmd->w0);
    }
    else if(!strcmp(arg[i],"pbc")){
      i++;
      if(i>=narg)
        error->all(FLERR,"Setting 'pbc' should be followed by a number in awpmd/cut");
      width_pbc=force->numeric(FLERR,arg[i]);
    }
    else if(!strcmp(arg[i],"relax"))
      wpmd->constraint=AWPMD::RELAX;
    else if(!strcmp(arg[i],"ermscale")){
      i++;
      if(i>=narg)
        error->all(FLERR,"Setting 'ermscale' should be followed by a number in awpmd/cut");
      ermscale=force->numeric(FLERR,arg[i]);
    }
    else if(!strcmp(arg[i],"flex_press"))
      flexible_pressure_flag = 1;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
// pair settings are as usual
void PairAWPMDCut::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR,"Incorrect args for pair coefficients");

  /*if(domain->xperiodic == 1 || domain->yperiodic == 1 ||
    domain->zperiodic == 1) {*/
  double delx = domain->boxhi[0]-domain->boxlo[0];
  double dely = domain->boxhi[1]-domain->boxlo[1];
  double delz = domain->boxhi[2]-domain->boxlo[2];
  half_box_length = 0.5 * MIN(delx, MIN(dely, delz));
  //}
  if(cut_global<0)
    cut_global=half_box_length;

  if (!allocated) {
    allocate();
  } else {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double cut_one = cut_global;
  if (narg == 3) cut_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAWPMDCut::init_style()
{
  // error and warning checks

  if (!atom->q_flag || !atom->spin_flag ||
      !atom->eradius_flag || !atom->erforce_flag )  // TO DO: adjust this to match approximation used
    error->all(FLERR,"Pair awpmd/cut requires atom attributes "
               "q, spin, eradius, erforce");

  /*
  if(vflag_atom){ // can't compute virial per atom
    //warning->
    error->all(FLERR,"Pair style awpmd can't compute per atom virials");
  }*/

  // add hook to minimizer for eradius and erforce

  if (update->whichflag == 2)
    int ignore = update->minimize->request(this,1,0.01);

  // make sure to use the appropriate timestep when using real units

  /*if (update->whichflag == 1) {
    if (force->qqr2e == 332.06371 && update->dt == 1.0)
      error->all(FLERR,"You must lower the default real units timestep for pEFF ");
  }*/

  // need a half neigh list and optionally a granular history neigh list

  //int irequest = neighbor->request(this,instance_me);

  //if (atom->tag_enable == 0)
  //  error->all(FLERR,"Pair style reax requires atom IDs");

  //if (force->newton_pair == 0)
    //error->all(FLERR,"Pair style awpmd requires newton pair on");

  //if (strcmp(update->unit_style,"real") != 0 && comm->me == 0)
    //error->warning(FLERR,"Not using real units with pair reax");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->newton = 2;

  if(force->e_mass==0. || force->hhmrr2e==0. || force->mvh2r==0.)
    error->all(FLERR,"Pair style awpmd requires e_mass and conversions hhmrr2e, mvh2r to be properly set for unit system");

  wpmd->me=force->e_mass;
  wpmd->h2_me=force->hhmrr2e/force->e_mass;
  wpmd->one_h=force->mvh2r;
  wpmd->coul_pref=force->qqrd2e;

  wpmd->calc_ii=1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAWPMDCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairAWPMDCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) fwrite(&cut[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairAWPMDCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) fread(&cut[i][j],sizeof(double),1,fp);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairAWPMDCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairAWPMDCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   returns pointers to the log() of electron radius and corresponding force
   minimizer operates on log(radius) so radius never goes negative
   these arrays are stored locally by pair style
------------------------------------------------------------------------- */

void PairAWPMDCut::min_xf_pointers(int ignore, double **xextra, double **fextra)
{
  // grow arrays if necessary
  // need to be atom->nmax in length
  int nvar=atom->nmax*(3+1+1+2);  // w(1), vel(3),  pw(1), cs(2)

  if (nvar > nmax) {
    memory->destroy(min_var);
    memory->destroy(min_varforce);
    nmax = nvar;
    memory->create(min_var,nmax,"pair:min_var");
    memory->create(min_varforce,nmax,"pair:min_varforce");
  }

  *xextra = min_var;
  *fextra = min_varforce;
}

/* ----------------------------------------------------------------------
   minimizer requests the log() of electron radius and corresponding force
   calculate and store in min_eradius and min_erforce
------------------------------------------------------------------------- */

void PairAWPMDCut::min_xf_get(int ignore)
{
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  double **v=atom->v;
  double *vforce=atom->vforce;
  double *ervel=atom->ervel;
  double *ervelforce=atom->ervelforce;
  double *cs=atom->cs;
  double *csforce=atom->csforce;

  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (spin[i]) {
      min_var[7*i] = log(eradius[i]);
      min_varforce[7*i] = eradius[i]*erforce[i];
      for(int j=0;j<3;j++){
        min_var[7*i+1+3*j] = v[i][j];
        min_varforce[7*i+1+3*j] = vforce[3*i+j];
      }
      min_var[7*i+4] = ervel[i];
      min_varforce[7*i+4] = ervelforce[i];
      min_var[7*i+5] = cs[2*i];
      min_varforce[7*i+5] = csforce[2*i];
      min_var[7*i+6] = cs[2*i+1];
      min_varforce[7*i+6] = csforce[2*i+1];

    } else {
      for(int j=0;j<7;j++)
        min_var[7*i+j] = min_varforce[7*i+j] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   propagate the minimizer values to the atom values
------------------------------------------------------------------------- */

void PairAWPMDCut::min_x_set(int ignore)
{
  double *eradius = atom->eradius;
  double **v=atom->v;
  double *ervel=atom->ervel;
  double *cs=atom->cs;

  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (spin[i]){
      eradius[i]=exp(min_var[7*i]);
      for(int j=0;j<3;j++)
        v[i][j]=min_var[7*i+1+3*j];
      ervel[i]=min_var[7*i+4];
      cs[2*i]=min_var[7*i+5];
      cs[2*i+1]=min_var[7*i+6];
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairAWPMDCut::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}
