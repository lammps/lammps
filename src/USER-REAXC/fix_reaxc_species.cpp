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
   Contributing author: Ray Shan (Sandia, tnshan@sandia.gov)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "stdlib.h"
#include "math.h"
#include "atom.h"
#include "string.h"
#include "fix_reaxc_species.h"
#include "update.h"
#include "domain.h"
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
#include "reaxc_types.h"
#include "reaxc_defs.h"

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

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  global_freq = nfreq = force->inumeric(FLERR,arg[5]);

  comm_forward = 1;
  
  if (nevery == 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix reax/c/species command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix reax/c/species command");

  if (neighbor->every != nfreq || neighbor->delay != 0 || neighbor->dist_check != 0){
    if (me == 0) {
      char str[128];
      sprintf(str,"Resetting reneighboring criteria for fix reax/c/species");
      error->warning(FLERR,str);
    }
    neighbor->every = nfreq;
    neighbor->delay = 0;
    neighbor->dist_check = 0;
  }

  if (me == 0) {
    fp = fopen(arg[6],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix reax/c/species file %s",arg[6]);
      error->one(FLERR,str);
    }
  }

  tmpq = NULL;
  tmpx = NULL;
  abo = NULL;
  BOCut = NULL;
  clusterID = NULL;

  Name = NULL;
  MolName = NULL;
  MolType = NULL;
  NMol = NULL;
  nd = NULL;
  molmap = NULL;

  singlepos_opened = multipos_opened = 0;
  multipos = 0;
  posfreq = 0;

  // initialize bond order cutoff
  int n, i, j;
  n = ntypes+1;
  memory->create(BOCut,n,n,"reaxc/c/species:BOCut");
  for (i = 1; i < n; i ++)
    for (j = 1; j < n; j ++)
      BOCut[i][j] = 0.0;

  // optional args
  eletype = NULL;
  ele = posspec = filepos = NULL;
  eleflag = posflag = padflag = 0;

  int iarg = 7;
  int itype, jtype;
  double bo_cut;

  while (iarg < narg) {
    // set BO cutoff
    if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix reax/c/species command");
      itype = force->inumeric(FLERR,arg[iarg+1]);
      jtype = force->inumeric(FLERR,arg[iarg+2]);
      bo_cut = force->numeric(FLERR,arg[iarg+3]);
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
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix species command");
      posflag = 1;
      posfreq = force->inumeric(FLERR,arg[iarg+1]);
      filepos = new char[n];
      strcpy(filepos,arg[iarg+2]);
      if (strchr(filepos,'*')) {
	multipos = 1;
      } else {
	if (me == 0) {
	  pos = fopen(filepos, "w");
	  if (pos == NULL) error->one(FLERR,"Cannot open fix species position file");
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
    ele[1]='H';
    ele[2]='O';
    ele[3]='N';
  }

  nmax = 0;
  vector_nmole = vector_nspec = 0;

  irepeat = 0;
  nvalid = nextvalid();

  memory->create(Name,ntypes,"reax/c/species:Name");

}

/* ---------------------------------------------------------------------- */

FixReaxCSpecies::~FixReaxCSpecies()
{
  memory->destroy(tmpq);
  memory->destroy(tmpx);
  memory->destroy(abo);
  memory->destroy(Name);
  memory->destroy(BOCut);
  memory->destroy(clusterID);

  if (me == 0) fclose(fp);
  if (me == 0 && posflag && multipos_opened) fclose(pos);
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

  post_integrate();
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix reax/c/specis unless atoms have IDs");

  reaxc = (PairReaxC *) force->pair_match("reax/c",1);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/specis without "
		  "pair_style reax/c");

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
  }

  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"reax/c/species") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one fix reax/c/species");

  // set default bond order cutoff
  int n, i, j;
  bg_cut = reaxc->control->bg_cut;
  n = ntypes+1;
  for (i = 1; i < n; i ++)
    for (j = 1; j < n; j ++)
      if (BOCut[i][j] == 0.0) BOCut[i][j] = bg_cut;

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::post_integrate()
{
  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  Output_ReaxC_Bonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::Output_ReaxC_Bonds(bigint ntimestep, FILE *fp)
{
  int i, j, ii,jj;
  int Nmole, Nspec;

  MPI_Barrier(world);
  nlocal = atom->nlocal;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(tmpq);
    memory->destroy(tmpx);
    memory->destroy(abo);
    memory->destroy(clusterID);
    memory->create(tmpq,nmax,"reax/c/species:tmpq");
    memory->create(tmpx,nmax,3,"reax/c/species:tmpx");
    memory->create(abo,nmax,nmax,"reax/c/species:abo");
    memory->create(clusterID,nmax,"reax/c/species:clusterID");
  }

  repeat = nrepeat;
  Nmole = Nspec = 0;
  vector_nmole = vector_nspec = 0;

  if (irepeat == 0)
    for (i = 0; i < nmax; i++) {
      tmpq[i] = 0.0;
      for (j = 0; j < nmax; j++)
        abo[i][j] = 0.0;
      for (j = 0; j < 3; j++)
        tmpx[i][j] = 0.0;
    }

  GatherBondOrder(lists);

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    return;
  }
  irepeat = 0;
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;

  FindMolecule();

  SortMolecule( Nmole);

  FindSpecies(Nmole, Nspec);

  vector_nmole = Nmole;
  vector_nspec = Nspec;

  if (me == 0) WriteFormulae( Nmole, Nspec);

  if (posflag && (ntimestep%posfreq==0)) WritePos(Nmole, Nspec);

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::GatherBondOrder(struct _reax_list *lists)
{
  int i, ii, j, jj, nj, rj, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double bo_tmp;
  bond_data *bo_ij;

  inum = reaxc->list->inum;
  ilist = reaxc->list->ilist;
  numneigh = reaxc->list->numneigh;
  firstneigh = reaxc->list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    tmpq[i] += atom->q[i]/repeat;
    for (jj = 0; jj < 3; jj++)
      tmpx[i][jj] += atom->x[i][jj]/repeat;
    
    jnum = numneigh[i];
    jlist = firstneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      for( nj = Start_Index(i, reaxc->lists); nj < End_Index(i, reaxc->lists); ++nj ) {
        bo_ij = &( reaxc->lists->select.bond_list[nj] );
        rj = bo_ij->nbr;

	if (atom->tag[j] == atom->tag[rj]) {
          bo_tmp = bo_ij->bo_data.BO;
	  abo[i][j] += bo_tmp/repeat;
	}
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::FindMolecule()
{
  int i,j,ii,jj,inum,jnum,n,itype,jtype;
  int change, done, anychange, loop, looptot;
  int *mask = atom->mask;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double bo_tmp, bo_cut;

  inum = reaxc->list->inum;
  ilist = reaxc->list->ilist;
  numneigh = reaxc->list->numneigh;
  firstneigh = reaxc->list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) clusterID[i] = atom->tag[i];
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
	jlist = firstneigh[i];
	jnum = numneigh[i];

	for (jj = 0; jj < jnum; jj++) {
	  j = jlist[jj];

	  j &= NEIGHMASK;
	  if (!(mask[j] & groupbit)) continue;
	  if (clusterID[i] == clusterID[j]) continue;

	  jtype = atom->type[j];
	  bo_cut = BOCut[itype][jtype];
	  bo_tmp = abo[i][j];

	  if (bo_tmp > bo_cut) {
	    clusterID[i] = clusterID[j] = MIN(clusterID[i],clusterID[j]);
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
    if (looptot >= 200*nprocs) break;
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
    error->warning(FLERR,"One or more reax/c cluster has atoms not in group");

  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    clusterID[n] = molmap[nint(clusterID[n])-idlo]+1;
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

    MPI_Reduce(Name,Nameall,ntypes,MPI_INT,MPI_SUM,0,world);
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

void FixReaxCSpecies::WriteFormulae(int Nmole, int Nspec)
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
    if (pos == NULL) error->one(FLERR,"Cannot open fix species position file");
  } else pos = NULL;
  multipos_opened = 1;

  delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */
void FixReaxCSpecies::WritePos(int Nmole, int Nspec)
{
  int i,itype,cid;
  int count, count_tmp, m, n, k;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double avq, avq_tmp, avx[3], avx_tmp, box[3];
  int *mask =atom->mask;
  int *Nameall;

  if (multipos) OpenPos();

  box[0] = domain->boxhi[0] - domain->boxlo[0];
  box[1] = domain->boxhi[1] - domain->boxlo[1];
  box[2] = domain->boxhi[2] - domain->boxlo[2];

  if (me == 0) {
    fprintf(pos,"Timestep" BIGINT_FORMAT "NMole %d  NSpec %d  xlo %f  "
		"xhi %f  ylo %f  yhi %f  zlo %f  zhi %f\n",
	        update->ntimestep,Nmole, Nspec,
		domain->boxlo[0],domain->boxhi[0],
		domain->boxlo[1],domain->boxhi[1],
		domain->boxlo[2],domain->boxhi[2]);

    fprintf(pos,"ID\tAtom_Count\tType\tAve_q\t\tCoM_x\t\tCoM_y\t\tCoM_z\n");
  }

  Nameall = NULL;
  memory->create(Nameall,ntypes,"species:Nameall");

  for (m = 1; m <= Nmole; m ++) {

    count = 0;
    avq = 0.0;
    for (n = 0; n < 3; n++) avx[n] = 0.0;
    for (n = 0; n < ntypes; n ++) Name[n] = 0;

    for (i = 0; i < nlocal; i ++) {
      if (!(mask[i] & groupbit)) continue;
      cid = nint(clusterID[i]);
      if (cid == m) {
        itype = atom->type[i]-1;
        Name[itype] ++;
	count ++;
	avq += tmpq[i];
	for (n = 0; n < 3; n++) avx[n] += tmpx[i][n];
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
      fprintf(pos,"%d\t%d\t\t",m,count);
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
  if (n == 0) {
    return vector_nmole;
  } else if (n == 1) {
    return vector_nspec;
  }
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

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixReaxCSpecies::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid-0;
}

/* ---------------------------------------------------------------------- */

int FixReaxCSpecies::pack_comm(int n, int *list, double *buf, 
				  int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = clusterID[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) clusterID[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

double FixReaxCSpecies::memory_usage()
{
  double bytes;

  bytes = 2.0*nmax*sizeof(double);
  bytes += nmax*nmax*sizeof(double);

  return bytes;
}
