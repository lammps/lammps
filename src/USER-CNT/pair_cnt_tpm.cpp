/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#include "pair_cnt_tpm.h"
#include "cntlist.h"
#include "export_cnt.h"

#include <mpi.h>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

#include <cstring>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

using namespace LAMMPS_NS;

// the cutoff distance between walls of tubes
static const double TPBRcutoff  = 3.0*3.4;

/* ---------------------------------------------------------------------- */

PairCNTTPM::PairCNTTPM(LAMMPS *lmp) : Pair(lmp) {
  writedata=1;
  BendingMode = 0;  // Harmonic bending model
  TPMType = 0;      // Inter-tube segment-segment interaction
  tab_path = NULL;
  tab_path_length = 0;
  
  eatom_s = NULL;
  eatom_b = NULL;
  eatom_t = NULL;
}

/* ---------------------------------------------------------------------- */

PairCNTTPM::~PairCNTTPM()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    
    memory->destroy(eatom_s);
    memory->destroy(eatom_b);
    memory->destroy(eatom_t);
  }
  if (tab_path != NULL) memory->destroy(tab_path);
}

/* ---------------------------------------------------------------------- */

void PairCNTTPM::compute(int eflag, int vflag){
  int nlocal = atom->nlocal;  //number of local atoms at the node
  //total number of atoms in the node and ghost shell
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;  //check "newton command"
  if(!newton_pair) error->all(FLERR,"set newton_pair");

  double **x = atom->x;
  double **f = atom->f;
  double *r = atom->radius;
  double *m = atom->rmass;
  double *l = atom->length;
  int *buckling = atom->buckling;
  tagint *g_id = atom->tag;
  tagint **bonds = atom->bond_cnt;
  //generate bonds and chain nblist
  CNTList cntlist(atom, list, cut_global*cut_global);

  //reorder data to make it contiguous within tubes
  //and compatable with Fortran functions
  std::vector<double> x_sort(3*nall), f_sort(3*nall), s_sort(9*nall);
  std::vector<double> u_ts_sort(nall), u_tb_sort(nall), u_tt_sort(nall);
  std::vector<int> b_sort(nall);
  for (int i = 0; i < nall; i++){
    int idx = cntlist.get_idx(i);
    for (int j = 0; j < 3; j++) x_sort[3*i+j] = x[idx][j];
    b_sort[i] = buckling[idx];
  }

  //bending potential
  int n_triplets = cntlist.get_triplets().size();
  for (int i = 0; i < n_triplets; i++) {
    const array2003<int,3>& t = cntlist.get_triplets()[i];
    //idx of nodes of a triplet in sorted representation
    int idx_s0 = cntlist.get_idxb(t[0]);
    int idx_s1 = cntlist.get_idxb(t[1]);
    int idx_s2 = cntlist.get_idxb(t[2]);

    double* X1 = &(x_sort[3*idx_s0]);
    double* X2 = &(x_sort[3*idx_s1]);
    double* X3 = &(x_sort[3*idx_s2]);
    double& U1b = u_tb_sort[idx_s0];
    double& U2b = u_tb_sort[idx_s1];
    double& U3b = u_tb_sort[idx_s2];
    double* F1 = &(f_sort[3*idx_s0]);
    double* F2 = &(f_sort[3*idx_s1]);
    double* F3 = &(f_sort[3*idx_s2]);
    double* S1 = &(s_sort[9*idx_s0]);
    double* S2 = &(s_sort[9*idx_s1]);
    double* S3 = &(s_sort[9*idx_s2]);
    double& R123 = r[t[1]];
    double& L123 = l[t[1]];
    int& BBF2 = b_sort[idx_s1];

    TubeBendingForceField(U1b, U2b, U3b, F1, F2, F3, S1, S2, S3, X1, X2, X3,
     R123, L123, BBF2);
  }

  //segment-segment and segment-tube interactions
  int n_segments = cntlist.get_segments().size();
  double Lmax = 0.0, Rmax = 0.0;
  double RT = get_R();
  for (int i = 0; i < n_segments; i++) {
    const array2003<int,2>& s = cntlist.get_segments()[i];
    //idx of a segment end 1 in sorted representation
    int idx_s0 = cntlist.get_idxb(s[0]);
    //idx of a segment end 2 in sorted representation
    int idx_s1 = cntlist.get_idxb(s[1]);
    double* X1 = &(x_sort[3*idx_s0]);
    double* X2 = &(x_sort[3*idx_s1]);
    double length = std::sqrt(std::pow(X1[0]-X2[0],2) +
     std::pow(X1[1]-X2[1],2) + std::pow(X1[2]-X2[2],2));
    if (length > Lmax) Lmax = length;
    double& U1t = u_tt_sort[idx_s0];
    double& U2t = u_tt_sort[idx_s1];
    double& U1s = u_ts_sort[idx_s0];
    double& U2s = u_ts_sort[idx_s1];
    double* F1 = &(f_sort[3*idx_s0]);
    double* F2 = &(f_sort[3*idx_s1]);
    double* S1 = &(s_sort[9*idx_s0]);
    double* S2 = &(s_sort[9*idx_s1]);
    double R12 = r[s[0]]; if (R12 > Rmax) Rmax = R12;
    if (std::abs(R12 - RT) > 1e-3) 
        error->all(FLERR,"Inconsistent input and potential table");
    //assume that the length of the segment is defined by the node with
    //smallest global id
    double L12 = (g_id[s[0]] > g_id[s[1]]) ? l[s[1]] : l[s[0]];
    TubeStretchingForceField(U1s, U2s, F1, F2, S1, S2, X1, X2, R12, L12);

    for (int nc = 0; nc < cntlist.get_nbs()[i].size(); nc++){
      //id of the beginning and end of the chain in the sorted representation
      const array2003<int,2>& chain = cntlist.get_nbs()[i][nc]; 
      int N = chain[1] - chain[0] + 1;  //number of elements in the chain
      int end1 = cntlist.get_idx(chain[0]);  //chain ends (real representation)
      int end2 = cntlist.get_idx(chain[1]);
      double* X = &(x_sort[3*chain[0]]);
      double* Ut = &(u_tt_sort[chain[0]]);
      double* Us = &(u_ts_sort[chain[0]]);
      double* F = &(f_sort[3*chain[0]]);
      double* S = &(s_sort[9*chain[0]]);
      double R = r[end1];
      int* BBF = &(b_sort[chain[0]]);
      int E1 = cntlist.is_end(end1);
      int E2 = cntlist.is_end(end2);

      int Ee = 0;
      double* Xe = X; double* Fe = F; double* Se = S;
      if (!E1 && cntlist.get_triplet(end1)[0] != CNTList::domain_end && 
       cntlist.get_triplet(cntlist.get_triplet(end1)[0])[0] == 
       CNTList::cnt_end){
        Ee = 1;
        int idx = cntlist.get_idxb(cntlist.get_triplet(end1)[0]);
        Xe = &(x_sort[3*idx]);
        Fe = &(f_sort[3*idx]);
        Se = &(s_sort[9*idx]);
      }
      else if (!E2 && cntlist.get_triplet(end2)[2] != CNTList::domain_end && 
       cntlist.get_triplet(cntlist.get_triplet(end2)[2])[2] == 
       CNTList::cnt_end){
        Ee = 2;
        int idx = cntlist.get_idxb(cntlist.get_triplet(end2)[2]);
        Xe = &(x_sort[3*idx]);
        Fe = &(f_sort[3*idx]);
        Se = &(s_sort[9*idx]);
      }
            
      SegmentTubeForceField(U1t, U2t, Ut, F1, F2, F, Fe, S1, S2, S, Se, X1,
       X2, R12, N, X, Xe, BBF, R, E1, E2, Ee, TPMType);
    }
  }

  if(neighbor->old_nrequest > 0){ //check if cutoff is chosen correctly
    double Rcut_min = std::max(2.0*Lmax, std::sqrt(0.5*Lmax*Lmax +
     std::pow((2.0*Rmax + TPBRcutoff),2)));
    if (cut_global < Rcut_min){
      std::cout << "L_max: = " << Lmax << ", R_max = " << Rmax << ", Rc = " 
       << cut_global << ", Rcut_min = " << Rcut_min << std::endl;
      error->all(FLERR,"The selected cutoff is too small for the current system");
    }
  }

  // set per atom values and accumulators
  // reallocate per-atom arrays if necessary
  if (atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    memory->destroy(eatom);
    memory->create(eatom,comm->nthreads*maxeatom,"pair:eatom");
    memory->destroy(eatom_s);
    memory->create(eatom_s,comm->nthreads*maxeatom,"pair:eatom_s");
    memory->destroy(eatom_b);
    memory->create(eatom_b,comm->nthreads*maxeatom,"pair:eatom_b");
    memory->destroy(eatom_t);
    memory->create(eatom_t,comm->nthreads*maxeatom,"pair:eatom_t");
  }

  if (atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom,comm->nthreads*maxvatom,6,"pair:vatom");
  }

  // zero accumulators
  eng_vdwl = 0.0; energy_s = 0.0;
  energy_b = 0.0; energy_t = 0.0;
  for (int i = 0; i < 6; i++) virial[i] = 0.0;
  for (int i = 0; i < nall; i++){
    eatom[i] = 0.0; eatom_s[i] = 0.0;
    eatom_b[i] = 0.0; eatom_t[i] = 0.0;
  }
  for (int i = 0; i < nall; i++) 
  for (int j = 0; j < 6; j++) vatom[i][j] = 0.0;

  //convert from sorted representation
  for (int i = 0; i < nall; i++){
    int idx = cntlist.get_idx(i);
    for (int j = 0; j < 3; j++) f[idx][j] = f_sort[3*i+j];
    eatom_s[idx] = u_ts_sort[i];
    eatom_b[idx] = u_tb_sort[i];
    eatom_t[idx] = u_tt_sort[i];
    eatom[idx] = u_ts_sort[i] + u_tb_sort[i] + u_tt_sort[i];
    energy_s += u_ts_sort[i];
    energy_b += u_tb_sort[i];
    energy_t += u_tt_sort[i];
    vatom[idx][0] = s_sort[9*i+0]; //xx
    vatom[idx][1] = s_sort[9*i+4]; //yy
    vatom[idx][2] = s_sort[9*i+8]; //zz
    vatom[idx][3] = s_sort[9*i+1]; //xy
    vatom[idx][4] = s_sort[9*i+2]; //xz
    vatom[idx][5] = s_sort[9*i+5]; //yz
    for (int j = 0; j < 6; j++) virial[j] += vatom[idx][j];
    buckling[idx] = b_sort[i];
  }
  eng_vdwl = energy_s + energy_b + energy_t;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCNTTPM::allocate(){
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCNTTPM::settings(int narg, char **arg){
  if ((narg == 0) || (narg > 4)) 
    error->all(FLERR,"Illegal pair_style command");
  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        cut[i][j] = cut_global;
  }
  std::string TPMAFile;
  if (narg > 1) {
    std::string path = arg[1];
    if(path.back() != '/') path += '/';
    tab_path_length = path.length();
    memory->create(tab_path,tab_path_length,"pair:path");
    std::memcpy(tab_path, path.c_str(), tab_path_length);
    std::string TPMSSTPFile = path + "TPMSSTP.xrs";
    TPMAFile = path + "TPMA.xrs";
    SetTablePath(TPMSSTPFile.c_str(), TPMSSTPFile.length(), TPMAFile.c_str(),
    TPMAFile.length());
  }
  else TPMAFile = "TPMA.xrs";
  if (narg > 2) {
    BendingMode = force->numeric(FLERR,arg[2]);
    if ((BendingMode < 0) || (BendingMode > 1))
      error->all(FLERR,"Incorrect BendingMode");
  }
  if (narg > 3){
    TPMType = force->numeric(FLERR,arg[3]);
    if ((TPMType < 0) || (TPMType > 1))
      error->all(FLERR,"Incorrect TPMType");
  }
  
  TPBInit();
  int M, N;
  std::ifstream in(TPMAFile);
  if (!in.is_open()) error->all(FLERR,"Incorrect table path");
  in >> M >> N;
  in.close();
  TPMInit(M, N);
  InitCNTPotModule(1, 3, 0, BendingMode, get_R());
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCNTTPM::coeff(int narg, char **arg){
  if ((narg < 2) || (narg > 3))
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

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
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCNTTPM::init_one(int i, int j){
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCNTTPM::write_restart(FILE *fp){
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCNTTPM::read_restart(FILE *fp){
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCNTTPM::write_restart_settings(FILE *fp){
  fwrite(&BendingMode,sizeof(int),1,fp);
  fwrite(&TPMType,sizeof(int),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&tab_path_length,sizeof(int),1,fp);
  fwrite(tab_path,tab_path_length,1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCNTTPM::read_restart_settings(FILE *fp){
  int me = comm->me;
  if (me == 0) {
    fread(&BendingMode,sizeof(int),1,fp);
    fread(&TPMType,sizeof(int),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&tab_path_length,sizeof(int),1,fp);
  }
  MPI_Bcast(&BendingMode,1,MPI_INT,0,world);
  MPI_Bcast(&TPMType,1,MPI_INT,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&tab_path_length,1,MPI_INT,0,world);
  
  memory->create(tab_path,tab_path_length,"pair:path");
  if (me == 0) fread(tab_path,tab_path_length,1,fp);
  MPI_Bcast(tab_path,tab_path_length,MPI_CHAR,0,world);
  
  if (tab_path != NULL) {
    std::string TPMSSTPFile = std::string(tab_path) + "TPMSSTP.xrs";
    std::string TPMAFile = std::string(tab_path) + "TPMA.xrs";
    SetTablePath(TPMSSTPFile.c_str(), TPMSSTPFile.length(), TPMAFile.c_str(),
     TPMAFile.length());
  }
  
  std::string TPMAFile = std::string((tab_path == NULL) ? "" : tab_path) 
   + "TPMA.xrs";
  TPBInit();
  int M, N;
  std::ifstream in(TPMAFile);
  if (!in.is_open()) error->all(FLERR,"Incorrect table path");
  in >> M >> N;
  in.close();
  TPMInit(M, N);
  InitCNTPotModule(1, 3, 0, BendingMode, get_R());
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairCNTTPM::write_data(FILE *fp){
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\n",i);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairCNTTPM::write_data_all(FILE *fp){
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g\n",i,j,cut[i][j]);
}


/* ---------------------------------------------------------------------- */

void PairCNTTPM::init_style(){
  //make sure that a full list is created (including ghost nodes)
  int r = neighbor->request(this,instance_me);
  neighbor->requests[r]->half = false;
  neighbor->requests[r]->full = true;
  neighbor->requests[r]->ghost = true;
}
