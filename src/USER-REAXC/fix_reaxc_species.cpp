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
   Contributing authors: Ray Shan (Sandia, tnshan@sandia.gov)
   			 Oleg Sergeev (VNIIA, sergeev@vniia.ru)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "stdlib.h"
#include "math.h"
#include "atom.h"
#include "string.h"
#include "fix_ave_atom.h"
#include "fix_reaxc_species.h"
#include "domain.h"
#include "update.h"
#include "pair_reax_c.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "reaxc_list.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxCSpecies::FixReaxCSpecies(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix reax/c/species command");

  force_reneighbor = 0;

  vector_flag = 1;
  size_vector = 2;

  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  global_freq = nfreq = atoi(arg[5]);

  comm_forward = 1;
  
  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix reax/c/species command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix reax/c/species command");
  
  // Neighbor lists must stay unchanged during averaging of bonds, 
  // but may be updated when no averaging is performed.
  
  int rene_flag = 0;
  if (nfreq % neighbor->every != 0 || neighbor->every < nevery * nrepeat) {
    int newneighborevery = nevery * nrepeat;
    while (nfreq % newneighborevery != 0 && newneighborevery <= nfreq / 2)
      newneighborevery++;

    if (nfreq % newneighborevery != 0)
      newneighborevery = nfreq;

    neighbor->every = newneighborevery;
    rene_flag = 1;
  }

  if (neighbor->delay != 0 || neighbor->dist_check != 0) {
    neighbor->delay = 0;
    neighbor->dist_check = 0;
    rene_flag = 1;
  }

  if (me == 0 && rene_flag) {
    char str[128];
    sprintf(str,"Resetting reneighboring criteria for fix reax/c/species");
    error->warning(FLERR,str);
  }

  tmparg = NULL;
  memory->create(tmparg,4,4,"reax/c/species:tmparg");
  strcpy(tmparg[0],arg[3]);
  strcpy(tmparg[1],arg[4]);
  strcpy(tmparg[2],arg[5]);

  if (me == 0) {
    fp = fopen(arg[6],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix reax/c/species file %s",arg[6]);
      error->one(FLERR,str);
    }
  }

  x0 = NULL;
  PBCconnected = NULL;
  clusterID = NULL;

  int ntmp = 1;
  memory->create(x0,ntmp,"reax/c/species:x0");
  memory->create(PBCconnected,ntmp,"reax/c/species:PBCconnected");
  memory->create(clusterID,ntmp,"reax/c/species:clusterID");
  vector_atom = clusterID;

  BOCut = NULL;
  Name = NULL;
  MolName = NULL;
  MolType = NULL;
  NMol = NULL;
  nd = NULL;
  molmap = NULL;

  nmax = 0;
  setupflag = 0;

  // set default bond order cutoff
  int n, i, j, itype, jtype;
  double bo_cut;
  bg_cut = 0.30;
  n = ntypes+1;
  memory->create(BOCut,n,n,"reax/c/species:BOCut");
  for (i = 1; i < n; i ++)
    for (j = 1; j < n; j ++)
      BOCut[i][j] = bg_cut;

  // optional args
  eletype = NULL;
  ele = filepos = NULL;
  eleflag = posflag = padflag = 0;

  singlepos_opened = multipos_opened = 0;
  multipos = 0;
  posfreq = 0;

  int iarg = 7;
  while (iarg < narg) {

    // set BO cutoff
    if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix reax/c/species command");
      itype = atoi(arg[iarg+1]);
      jtype = atoi(arg[iarg+2]);
      bo_cut = atof(arg[iarg+3]);
      if (itype > ntypes || jtype > ntypes) 
      	error->all(FLERR,"Illegal fix reax/c/species command");
      if (itype <= 0 || jtype <= 0) 
      	error->all(FLERR,"Illegal fix reax/c/species command");
      if (bo_cut > 1.0 || bo_cut < 0.0)
      	error->all(FLERR,"Illegal fix reax/c/species command");

      BOCut[itype][jtype] = bo_cut; 
      BOCut[jtype][itype] = bo_cut;
      iarg += 4;

    // modify element type names
    } else if (strcmp(arg[iarg],"element") == 0) {
      if (iarg+ntypes+1 > narg) error->all(FLERR,"Illegal fix reax/c/species command");

      eletype = (char**) malloc(ntypes*sizeof(char*));
      for (int i = 0; i < ntypes; i ++) {
        eletype[i] = (char*) malloc(2*sizeof(char));
      	strcpy(eletype[i],arg[iarg+1+i]);
      }
      eleflag = 1;
      iarg += ntypes + 1;

    // position of molecules
    } else if (strcmp(arg[iarg],"position") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix reax/c/species command");
      posflag = 1;
      posfreq = atoi(arg[iarg+1]);
      if (posfreq < nfreq || (posfreq%nfreq != 0))
	error->all(FLERR,"Illegal fix reax/c/species command");

      filepos = new char[255];
      strcpy(filepos,arg[iarg+2]);
      if (strchr(filepos,'*')) {
        multipos = 1;
      } else {
        if (me == 0) {
          pos = fopen(filepos, "w");
          if (pos == NULL) error->one(FLERR,"Cannot open fix reax/c/species position file");
        }
      	singlepos_opened = 1;
      	multipos = 0;
      }
      iarg += 3;
    } else error->all(FLERR,"Illegal fix reax/c/species command");
  }

  if (!eleflag) {
    memory->create(ele,ntypes+1,"reax/c/species:ele");
    ele[0]='C';
    if (ntypes > 1)
      ele[1]='H';
    if (ntypes > 2)
      ele[2]='O';
    if (ntypes > 3)
      ele[3]='N';
  }

  vector_nmole = 0;
  vector_nspec = 0;

}

/* ---------------------------------------------------------------------- */

FixReaxCSpecies::~FixReaxCSpecies()
{
  memory->destroy(ele);
  memory->destroy(BOCut);
  memory->destroy(clusterID);
  memory->destroy(PBCconnected);
  memory->destroy(x0);

  memory->destroy(nd);
  memory->destroy(Name);
  memory->destroy(NMol);
  memory->destroy(MolType);
  memory->destroy(MolName);
  memory->destroy(tmparg);

  if (filepos)
    delete [] filepos;

  if (me == 0) fclose(fp);
  if (me == 0 && posflag && multipos_opened) fclose(pos);

  modify->delete_compute("SPECATOM");
  modify->delete_fix("SPECBOND");
}

/* ---------------------------------------------------------------------- */

int FixReaxCSpecies::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::setup(int vflag)
{
  ntotal = static_cast<int> (atom->natoms);
  memory->create(Name,ntypes,"reax/c/species:Name");

  post_integrate();
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix reax/c/species unless atoms have IDs");

  reaxc = (PairReaxC *) force->pair_match("reax/c",1);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/species without "
		  "pair_style reax/c");

  reaxc->fixspecies_flag = 1;
  nvalid = update->ntimestep+nfreq;

  // check if this fix has been called twice
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"reax/c/species") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one fix reax/c/species");

  if (!setupflag) {
    // create a compute to store properties
    create_compute();

    // create a fix to point to fix_ave_atom for averaging stored properties
    create_fix();

    setupflag = 1;
  }

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::create_compute()
{
  int narg;
  char **args;

  narg = 34;
  args = new char*[narg];
  args[0]  = (char *) "SPECATOM";
  args[1]  = (char *) "all";
  args[2]  = (char *) "spec/atom";
  args[3]  = (char *) "q";
  args[4]  = (char *) "x";
  args[5]  = (char *) "y";
  args[6]  = (char *) "z";
  args[7]  = (char *) "vx";
  args[8]  = (char *) "vy";
  args[9]  = (char *) "vz";
  args[10] = (char *) "abo01";
  args[11] = (char *) "abo02";
  args[12] = (char *) "abo03";
  args[13] = (char *) "abo04";
  args[14] = (char *) "abo05";
  args[15] = (char *) "abo06";
  args[16] = (char *) "abo07";
  args[17] = (char *) "abo08";
  args[18] = (char *) "abo09";
  args[19] = (char *) "abo10";
  args[20] = (char *) "abo11";
  args[21] = (char *) "abo12";
  args[22] = (char *) "abo13";
  args[23] = (char *) "abo14";
  args[24] = (char *) "abo15";
  args[25] = (char *) "abo16";
  args[26] = (char *) "abo17";
  args[27] = (char *) "abo18";
  args[28] = (char *) "abo19";
  args[29] = (char *) "abo20";
  args[30] = (char *) "abo21";
  args[31] = (char *) "abo22";
  args[32] = (char *) "abo23";
  args[33] = (char *) "abo24";
  modify->add_compute(narg,args);
  delete [] args;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::create_fix()
{
  int narg;
  char **args;

  narg = 37;
  args = new char*[narg];
  args[0]  = (char *) "SPECBOND";
  args[1]  = (char *) "all";
  args[2]  = (char *) "ave/atom";
  args[3]  = tmparg[0];
  args[4]  = tmparg[1];
  args[5]  = tmparg[2];
  args[6]  = (char *) "c_SPECATOM[1]";	 // q, array_atoms[i][0]
  args[7]  = (char *) "c_SPECATOM[2]";	 // x, 1
  args[8]  = (char *) "c_SPECATOM[3]";	 // y, 2
  args[9]  = (char *) "c_SPECATOM[4]";	 // z, 3
  args[10] = (char *) "c_SPECATOM[5]";	 // vx, 4
  args[11] = (char *) "c_SPECATOM[6]";	 // vy, 5
  args[12] = (char *) "c_SPECATOM[7]";	 // vz, 6
  args[13] = (char *) "c_SPECATOM[8]";	 // abo01, 7
  args[14] = (char *) "c_SPECATOM[9]";	
  args[15] = (char *) "c_SPECATOM[10]";	
  args[16] = (char *) "c_SPECATOM[11]"; 
  args[17] = (char *) "c_SPECATOM[12]";
  args[18] = (char *) "c_SPECATOM[13]";
  args[19] = (char *) "c_SPECATOM[14]"; 
  args[20] = (char *) "c_SPECATOM[15]";
  args[21] = (char *) "c_SPECATOM[16]";
  args[22] = (char *) "c_SPECATOM[17]";
  args[23] = (char *) "c_SPECATOM[18]";
  args[24] = (char *) "c_SPECATOM[19]"; // abo12, 18
  args[25] = (char *) "c_SPECATOM[20]";	
  args[26] = (char *) "c_SPECATOM[21]"; 
  args[27] = (char *) "c_SPECATOM[22]";
  args[28] = (char *) "c_SPECATOM[23]";
  args[29] = (char *) "c_SPECATOM[24]"; 
  args[30] = (char *) "c_SPECATOM[25]";
  args[31] = (char *) "c_SPECATOM[26]";
  args[32] = (char *) "c_SPECATOM[27]";
  args[33] = (char *) "c_SPECATOM[28]";
  args[34] = (char *) "c_SPECATOM[29]"; 
  args[35] = (char *) "c_SPECATOM[30]";
  args[36] = (char *) "c_SPECATOM[31]";
  modify->add_fix(narg,args);
  f_SPECBOND = (FixAveAtom *) modify->fix[modify->nfix-1];
  delete [] args;

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::post_integrate()
{
  Output_ReaxC_Bonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::Output_ReaxC_Bonds(bigint ntimestep, FILE *fp)

{
  int Nmole, Nspec;

  // point to fix_ave_atom
  f_SPECBOND->end_of_step();

  if (ntimestep != nvalid) return;

  nlocal = atom->nlocal;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(x0);
    memory->destroy(PBCconnected);
    memory->destroy(clusterID);
    memory->create(x0,nmax,"reax/c/species:x0");
    memory->create(PBCconnected,nmax,"reax/c/species:PBCconnected");
    memory->create(clusterID,nmax,"reax/c/species:clusterID");
    vector_atom = clusterID;
  }

  for (int i = 0; i < nmax; i++) {
    PBCconnected[i] = 0;
    x0[i].x = x0[i].y = x0[i].z = 0.0;
  }

  Nmole = Nspec = 0;

  FindMolecule();

  SortMolecule (Nmole);

  FindSpecies(Nmole, Nspec);

  vector_nmole = Nmole;
  vector_nspec = Nspec;

  if (me == 0 && ntimestep >= 0)
    WriteFormulas (Nmole, Nspec);

  if (posflag && ((ntimestep)%posfreq==0)) {
    WritePos(Nmole, Nspec);
    if (me == 0) fflush(pos);
  }

  nvalid += nfreq;
}

/* ---------------------------------------------------------------------- */

AtomCoord chAnchor(AtomCoord in1, AtomCoord in2)
{
  if (in1.x < in2.x) 
    return in1;
  return in2;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::FindMolecule ()
{
  int i,j,ii,jj,inum,jnum,n,itype,jtype,itag,jtag,loop,looptot;
  int change,done,anychange;
  int *mask = atom->mask;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double bo_tmp,bo_cut;
  double **spec_atom = f_SPECBOND->array_atom;

  inum = reaxc->list->inum;
  ilist = reaxc->list->ilist;
  numneigh = reaxc->list->numneigh;
  firstneigh = reaxc->list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      clusterID[i] = atom->tag[i];
      x0[i].x = spec_atom[i][1];
      x0[i].y = spec_atom[i][2];
      x0[i].z = spec_atom[i][3];
    }
    else clusterID[i] = 0.0;
  }

  loop = 0;
  while (1) {
    comm->forward_comm_fix(this);
    loop ++;

    change = 0;
    while (1) {
      done = 1;

      for (ii = 0; ii < inum; ii++) {
      	i = ilist[ii];
      	if (!(mask[i] & groupbit)) continue;

      	itype = atom->type[i];

        for (jj = 0; jj < MAXSPECBOND; jj++) {
      	  j = reaxc->tmpid[i][jj];

      	  if (j < i) continue;
      	  if (!(mask[j] & groupbit)) continue;

      	  if (clusterID[i] == clusterID[j] && PBCconnected[i] == PBCconnected[j] 
	    && x0[i].x == x0[j].x && x0[i].y == x0[j].y && x0[i].z == x0[j].z) continue;

          jtype = atom->type[j];
      	  bo_cut = BOCut[itype][jtype];
      	  bo_tmp = spec_atom[i][jj+7];

      	  if (bo_tmp > bo_cut) {
            clusterID[i] = clusterID[j] = MIN(clusterID[i], clusterID[j]);
            PBCconnected[i] = PBCconnected[j] = MAX(PBCconnected[i], PBCconnected[j]);
            x0[i] = x0[j] = chAnchor(x0[i], x0[j]);
            if ((fabs(spec_atom[i][1] - spec_atom[j][1]) > reaxc->control->bond_cut)
             || (fabs(spec_atom[i][2] - spec_atom[j][2]) > reaxc->control->bond_cut) 
             || (fabs(spec_atom[i][3] - spec_atom[j][3]) > reaxc->control->bond_cut))
              PBCconnected[i] = PBCconnected[j] = 1;
      	    done = 0;
      	  }
      	}
      }
      if (!done) change = 1;
      if (done) break;
    }
    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    if (!anychange) break;

    MPI_Allreduce(&loop,&looptot,1,MPI_INT,MPI_SUM,world);
    if (looptot >= 400*nprocs) break;

  }
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::SortMolecule(int &Nmole)
{
  memory->destroy(molmap);
  molmap = NULL;

  int m, n, idlo, idhi;
  int *mask =atom->mask;
  int lo = ntotal;
  int hi = -ntotal;
  int flag = 0;
  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    if (clusterID[n] == 0.0) flag = 1;
    lo = MIN(lo,nint(clusterID[n]));
    hi = MAX(hi,nint(clusterID[n]));
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && me == 0) 
    error->warning(FLERR,"Atom with cluster ID = 0 included in "
		    "fix reax/c/species group");
  MPI_Allreduce(&lo,&idlo,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&hi,&idhi,1,MPI_INT,MPI_MAX,world);
  if (idlo == ntotal)
    if (me == 0)
      error->warning(FLERR,"Atom with cluster ID = maxmol "
		    "included in fix reax/c/species group");

  int nlen = idhi - idlo + 1;
  memory->create(molmap,nlen,"reax/c/species:molmap");
  for (n = 0; n < nlen; n++) molmap[n] = 0;

  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    molmap[nint(clusterID[n])-idlo] = 1;
  }

  int *molmapall;
  memory->create(molmapall,nlen,"reax/c/species:molmapall");
  MPI_Allreduce(molmap,molmapall,nlen,MPI_INT,MPI_MAX,world);

  Nmole = 0;
  for (n = 0; n < nlen; n++) {
    if (molmapall[n]) molmap[n] = Nmole++;
    else molmap[n] = -1;
  }
  memory->destroy(molmapall);

  flag = 0;
  for (n = 0; n < nlocal; n++) {
    if (mask[n] & groupbit) continue;
    if (nint(clusterID[n]) < idlo || nint(clusterID[n]) > idhi) continue;
    if (molmap[nint(clusterID[n])-idlo] >= 0) flag = 1;
  }

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"One or more cluster has atoms not in group");

  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    clusterID[n] = molmap[nint(clusterID[n])-idlo] + 1;
  }

  memory->destroy(molmap);
  molmap = NULL;

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::FindSpecies(int Nmole, int &Nspec)
{
  int inum, *ilist;
  int i, j, k, l, m, n, itype, cid;
  int flag_identity, flag_mol, flag_spec;
  int flag_tmp;
  int *mask =atom->mask;
  int *Nameall, *NMolall;

  memory->destroy(MolName);
  MolName = NULL;
  memory->create(MolName,Nmole*(ntypes+1),"reax/c/species:MolName");

  memory->destroy(NMol);
  NMol = NULL;
  memory->create(NMol,Nmole,"reax/c/species:NMol");
  for (m = 0; m < Nmole; m ++)
    NMol[m] = 1;

  memory->create(Nameall,ntypes,"reax/c/species:Nameall");
  memory->create(NMolall,Nmole,"reax/c/species:NMolall");

  for (m = 1, Nspec = 0; m <= Nmole; m ++) {
    for (n = 0; n < ntypes; n ++) Name[n] = 0;
    for (n = 0, flag_mol = 0; n < nlocal; n ++) {
      if (!(mask[n] & groupbit)) continue;
      cid = nint(clusterID[n]);
      if (cid == m) {
        itype = atom->type[n]-1;
        Name[itype] ++;
        flag_mol = 1; 
      }
    }
    MPI_Allreduce(&flag_mol,&flag_tmp,1,MPI_INT,MPI_MAX,world);
    flag_mol = flag_tmp;

    MPI_Allreduce(Name,Nameall,ntypes,MPI_INT,MPI_SUM,world);
    for (n = 0; n < ntypes; n++) Name[n] = Nameall[n];

    if (flag_mol == 1) {
      flag_identity = 1;
      for (k = 0; k < Nspec; k ++) {
        flag_spec=0;
        for (l = 0; l < ntypes; l ++) 
          if (MolName[ntypes*k+l] != Name[l]) flag_spec = 1;
        if (flag_spec == 0) NMol[k] ++;  
        flag_identity *= flag_spec; 
      }
      if (Nspec == 0 || flag_identity == 1) {
        for (l = 0; l < ntypes; l ++)
          MolName[ntypes*Nspec+l] = Name[l];
        Nspec ++; 
      }
    }
  }
  memory->destroy(NMolall);
  memory->destroy(Nameall);

  memory->destroy(nd);
  nd = NULL;
  memory->create(nd,Nspec,"reax/c/species:nd");

  memory->destroy(MolType);
  MolType = NULL;
  memory->create(MolType,Nspec*(ntypes+2),"reax/c/species:MolType");
}

/* ---------------------------------------------------------------------- */

int FixReaxCSpecies::CheckExistence(int id, int ntypes)
{
  int i, j, molid, flag;

  for (i = 0; i < Nmoltype; i ++) {
    flag = 0;
    for (j = 0; j < ntypes; j ++) {
      molid = MolType[ntypes * i + j];
      if (molid != MolName[ntypes * id + j]) flag = 1;
    }
    if (flag == 0) return i; 
  }
  for (i = 0; i < ntypes; i ++)
    MolType[ntypes * Nmoltype + i] = MolName[ntypes *id + i];

  Nmoltype ++;
  return Nmoltype - 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::WriteFormulas(int Nmole, int Nspec)
{
  int i, j, itemp;
  bigint ntimestep = update->ntimestep;

  fprintf(fp,"# Timestep     No_Moles     No_Specs     ");

  Nmoltype = 0;
 
  for (i = 0; i < Nspec; i ++)
    nd[i] = CheckExistence(i, ntypes);

  for (i = 0; i < Nmoltype; i ++) {
    for (j = 0;j < ntypes; j ++) {
      itemp = MolType[ntypes * i + j];
      if (itemp != 0) {
        if (eletype) fprintf(fp,"%s",eletype[j]);
        else fprintf(fp,"%c",ele[j]);
        if (itemp != 1) fprintf(fp,"%d",itemp);
      }
    }
    fprintf(fp,"\t");
  }
  fprintf(fp,"\n");

  fprintf(fp,BIGINT_FORMAT,ntimestep);
  fprintf(fp,"%11d%11d\t",Nmole,Nspec);

  for (i = 0; i < Nmoltype; i ++)
    fprintf(fp," %d\t",NMol[i]);
  fprintf(fp,"\n");

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::OpenPos()
{
  char *filecurrent;
  bigint ntimestep = update->ntimestep;

  filecurrent = (char*) malloc((strlen(filepos)+16)*sizeof(char));
  char *ptr = strchr(filepos,'*');
  *ptr = '\0';
  if (padflag == 0) 
    sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
	filepos,ntimestep,ptr+1);
  else {
    char bif[8],pad[16];
    strcpy(bif,BIGINT_FORMAT);
    sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
    sprintf(filecurrent,pad,filepos,ntimestep,ptr+1);
  }
  *ptr = '*';

  if (me == 0) {
    pos = fopen(filecurrent, "w");
    if (pos == NULL) error->one(FLERR,"Cannot open fix reax/c/species position file");
  } else pos = NULL;
  multipos_opened = 1;

  free(filecurrent);
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::WritePos(int Nmole, int Nspec)
{
  int i, itype, cid;
  int count, count_tmp, m, n, k;
  int *Nameall;
  int *mask =atom->mask;
  double avq, avq_tmp, avx[3], avx_tmp, box[3], halfbox[3];
  double **spec_atom = f_SPECBOND->array_atom;

  if (multipos) OpenPos();
  
  box[0] = domain->boxhi[0] - domain->boxlo[0];
  box[1] = domain->boxhi[1] - domain->boxlo[1];
  box[2] = domain->boxhi[2] - domain->boxlo[2];

  for (int j = 0; j < 3; j++)
    halfbox[j] = box[j] / 2;

  if (me == 0) {
    fprintf(pos,"Timestep "BIGINT_FORMAT " NMole %d  NSpec %d  xlo %f  "
		"xhi %f  ylo %f  yhi %f  zlo %f  zhi %f\n",
	        update->ntimestep,Nmole, Nspec,
		domain->boxlo[0],domain->boxhi[0],
		domain->boxlo[1],domain->boxhi[1],
		domain->boxlo[2],domain->boxhi[2]);

    fprintf(pos,"ID\tAtom_Count\tType\tAve_q\t\tCoM_x\t\tCoM_y\t\tCoM_z\n");
  }

  Nameall = NULL;
  memory->create(Nameall,ntypes,"reax/c/species:Nameall");

  for (m = 1; m <= Nmole; m ++) {

    count = 0;
    avq = 0.0;
    for (n = 0; n < 3; n++) 
      avx[n] = 0.0;
    for (n = 0; n < ntypes; n ++) 
      Name[n] = 0;

    for (i = 0; i < nlocal; i ++) {
      if (!(mask[i] & groupbit)) continue;
      cid = nint(clusterID[i]);
      if (cid == m) {
        itype = atom->type[i]-1;
        Name[itype] ++;
      	count ++;
      	avq += spec_atom[i][0];
        if (PBCconnected[i]) {
          if ((x0[i].x - spec_atom[i][1]) > halfbox[0])
            spec_atom[i][1] += box[0];
          if ((spec_atom[i][1] - x0[i].x) > halfbox[0])
             spec_atom[i][1] -= box[0];
          if ((x0[i].y - spec_atom[i][2]) > halfbox[1])
              spec_atom[i][2] += box[1];
          if ((spec_atom[i][2] - x0[i].y) > halfbox[1])
             spec_atom[i][2] -= box[1];
          if ((x0[i].z - spec_atom[i][3]) > halfbox[2])
            spec_atom[i][3] += box[2];
          if ((spec_atom[i][3] - x0[i].z) > halfbox[2])
           spec_atom[i][3] -= box[2];
        }
        for (n = 0; n < 3; n++)
          avx[n] += spec_atom[i][n+1];
      }
    }

    avq_tmp = 0.0;
    MPI_Allreduce(&avq,&avq_tmp,1,MPI_DOUBLE,MPI_SUM,world);
    avq = avq_tmp;
    
    for (n = 0; n < 3; n++) {
      avx_tmp = 0.0;
      MPI_Reduce(&avx[n],&avx_tmp,1,MPI_DOUBLE,MPI_SUM,0,world);
      avx[n] = avx_tmp;
    }

    MPI_Reduce(&count,&count_tmp,1,MPI_INT,MPI_SUM,0,world);
    count = count_tmp;

    MPI_Reduce(Name,Nameall,ntypes,MPI_INT,MPI_SUM,0,world);
    for (n = 0; n < ntypes; n++) Name[n] = Nameall[n];

    if (me == 0) { 
      fprintf(pos,"%d\t%d\t",m,count);
      for (n = 0; n < ntypes; n++) { 
        if (Name[n] != 0) {
          if (eletype) fprintf(pos,"%s",eletype[n]);
          else fprintf(pos,"%c",ele[n]);
          if (Name[n] != 1) fprintf(pos,"%d",Name[n]);
        }
      }
      if (count > 0) {
        avq /= count;
        for (k = 0; k < 3; k++) { 
      	  avx[k] /= count;
          if (avx[k] >= domain->boxhi[k])
            avx[k] -= box[k];
          if (avx[k] < domain->boxlo[k])
            avx[k] += box[k];

      	  avx[k] -= domain->boxlo[k];
      	  avx[k] /= box[k];
        }
        fprintf(pos,"\t%.8f \t%.8f \t%.8f \t%.8f",
		    avq,avx[0],avx[1],avx[2]);
      }
      fprintf(pos,"\n");
    }
  }
  if (me == 0 && !multipos) fprintf(pos,"#\n");
  memory->destroy(Nameall);
}

/* ---------------------------------------------------------------------- */

double FixReaxCSpecies::compute_vector(int n)
{
  if (n == 0) 
    return vector_nmole;
  if (n == 1) 
    return vector_nspec;
  return 0.0;

}

/* ---------------------------------------------------------------------- */

int FixReaxCSpecies::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}

/* ---------------------------------------------------------------------- */

int FixReaxCSpecies::pack_comm(int n, int *list, double *buf, 
				  int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m] = clusterID[j];
    buf[m+1] = (double)PBCconnected[j];
    buf[m+2] = x0[j].x;
    buf[m+3] = x0[j].y;
    buf[m+4] = x0[j].z;
    m += 5;
  }
  return 5;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    clusterID[i] = buf[m];
    PBCconnected[i] = (int)buf[m+1];
    x0[i].x = buf[m+2];
    x0[i].y = buf[m+3];
    x0[i].z = buf[m+4];
    m += 5;
  }
}

/* ---------------------------------------------------------------------- */

double FixReaxCSpecies::memory_usage()
{
  double bytes;

  bytes = 5*nmax*sizeof(double);  // clusterID + PBCconnected + x0

  return bytes;
}

/* ---------------------------------------------------------------------- */
