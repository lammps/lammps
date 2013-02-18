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
#include "fix_ave_atom.h"
#include "fix_species.h"
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
#include "reaxc_defs.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSpecies::FixSpecies(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix species command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  global_freq = nfreq = atoi(arg[5]);

  comm_forward = 1;
  
  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix species command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix species command");

  tmparg = NULL;
  memory->create(tmparg,4,4,"species:tmparg");
  strcpy(tmparg[0],arg[3]);
  strcpy(tmparg[1],arg[4]);
  strcpy(tmparg[2],arg[5]);

  if (me == 0) {
    fp = fopen(arg[6],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix species file %s",arg[6]);
      error->one(FLERR,str);
    }
  }

  BOCut = NULL;
  clusterID = NULL;

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
  memory->create(BOCut,n,n,"/species:BOCut");
  for (i = 1; i < n; i ++)
    for (j = 1; j < n; j ++)
      BOCut[i][j] = bg_cut;

  // optional args
  eletype = NULL;
  ele = posspec = filepos = NULL;
  eleflag = posflag = padflag = 0;

  int iarg = 7;
  while (iarg < narg) {

    // set BO cutoff
    if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix species command");
      itype = atoi(arg[iarg+1]);
      jtype = atoi(arg[iarg+2]);
      bo_cut = atof(arg[iarg+3]);
      if (itype > ntypes || jtype > ntypes) 
	error->all(FLERR,"Illegal fix species command");
      if (itype <= 0 || jtype <= 0) 
	error->all(FLERR,"Illegal fix species command");
      if (bo_cut > 1.0 || bo_cut < 0.0)
	error->all(FLERR,"Illegal fix species command");

      BOCut[itype][jtype] = bo_cut; 
      BOCut[jtype][itype] = bo_cut;
      iarg += 4;

    // modify element type names
    } else if (strcmp(arg[iarg],"element") == 0) {
      if (iarg+ntypes+1 > narg) error->all(FLERR,"Illegal fix species command");

      int nchar = 2;
      eletype = (char**) malloc(ntypes*sizeof(char*));
      for (int i = 0; i < ntypes; i ++) {
	if (strlen(arg[iarg+1+i]) > nchar)
	  error->all(FLERR,"Illegal fix species command");
        eletype[i] = (char*) malloc(nchar*sizeof(char));
	strcpy(eletype[i],arg[iarg+1+i]);
      }
      eleflag = 1;
      iarg += ntypes + 1;

    } else error->all(FLERR,"Illegal fix species command");
  }

  if (!eleflag) {
    memory->create(ele,ntypes+1,"species:ele");
    ele[0]='C';
    ele[1]='H';
    ele[2]='O';
    ele[3]='N';
  }

}

/* ---------------------------------------------------------------------- */

FixSpecies::~FixSpecies()
{
  memory->destroy(ele);
  memory->destroy(BOCut);
  memory->destroy(clusterID);

  memory->destroy(nd);
  memory->destroy(Name);
  memory->destroy(NMol);
  memory->destroy(MolType);
  memory->destroy(MolName);

  if (me == 0) fclose(fp);

  modify->delete_compute("SPECATOM");
  modify->delete_fix("SPECBOND");
}

/* ---------------------------------------------------------------------- */

int FixSpecies::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpecies::setup(int vflag)
{
  ntotal = static_cast<int> (atom->natoms);
  memory->create(Name,ntypes,"species:Name");
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixSpecies::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix specis unless atoms have IDs");

  reaxc = (PairReaxC *) force->pair_match("reax/c",1);

  if (reaxc == NULL) error->all(FLERR,"Cannot use fix species without "
		  "pair_style reax/c");

  reaxc->fixspecies_flag = 1;

  nvalid = update->ntimestep+nfreq;

  // request neighbor list
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  // check if this fix has been called twice
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"species") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one fix species");

  if (!setupflag) {
    // create a compute to store properties
    create_compute();

    // create a fix to point to fix_ave_atom for averaging stored properties
    create_fix();

    setupflag = 1;
  }

}

/* ---------------------------------------------------------------------- */

void FixSpecies::create_compute()
{
  int narg;
  char **args;

  narg = 22;
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
  modify->add_compute(narg,args);
  delete [] args;
}

/* ---------------------------------------------------------------------- */

void FixSpecies::create_fix()
{
  int narg;
  char **args;

  narg = 25;
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
  modify->add_fix(narg,args);
  f_SPECBOND = (FixAveAtom *) modify->fix[modify->nfix-1];
  delete [] args;

}

/* ---------------------------------------------------------------------- */

void FixSpecies::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixSpecies::end_of_step()
{
  Output_ReaxC_Bonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixSpecies::Output_ReaxC_Bonds(bigint ntimestep, FILE *fp)

{
  int i, j, k, itype, jtype, itag, jtag;
  int b, nbuf, nbuf_local, inode;
  int nlocal_max, numbonds, numbonds_max, count;
  int Nmole, Nspec;
  double *bbuf;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int ii, jj, inum, jnum;

  // point to fix_ave_atom
  f_SPECBOND->end_of_step();

  if (ntimestep != nvalid) return;

  nlocal = atom->nlocal;
  nghost = atom->nghost;
  nall = nlocal + nghost;
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(clusterID);
    memory->create(clusterID,nmax,"species:clusterID");
  }

  Nmole = Nspec = count = 0;

  FindMolecule();

  SortMolecule (Nmole);

  FindSpecies(Nmole, Nspec);

  if (me == 0)
    WriteFormulas (Nmole, Nspec);

  nvalid += nfreq;
}

/* ---------------------------------------------------------------------- */

void FixSpecies::FindMolecule ()
{
  int i,j,ii,jj,inum,jnum,itype,jtype,jtag,k,ktag,itag,loop,ntot;
  int change,done,anychange;
  int *mask = atom->mask;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double bo_tmp,bo_cut;
  double **spec_atom = f_SPECBOND->array_atom;

  neighbor->build_one(list->index);
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) clusterID[i] = atom->tag[i];
    else clusterID[i] = 0;
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
	  if (clusterID[i] == clusterID[j]) continue;

          jtype = atom->type[j];
	  bo_cut = BOCut[itype][jtype];
	  bo_tmp = spec_atom[i][jj+7];

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

    MPI_Allreduce(&loop,&ntot,1,MPI_INT,MPI_SUM,world);
    if (ntot >= 20*nprocs) break;
  }
}

/* ---------------------------------------------------------------------- */

void FixSpecies::SortMolecule(int &Nmole)
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
    if (clusterID[n] == 0) flag = 1;
    lo = MIN(lo,nint(clusterID[n]));
    hi = MAX(hi,nint(clusterID[n]));
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && me == 0) 
    error->warning(FLERR,"Atom with cluster ID = 0 included in "
		    "fix species group");
  MPI_Allreduce(&lo,&idlo,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&hi,&idhi,1,MPI_INT,MPI_MAX,world);
  if (idlo == ntotal)
    if (me == 0)
      error->warning(FLERR,"Atom with cluster ID = maxmol "
		    "included in fix species group");

  int nlen = idhi - idlo + 1;
  memory->create(molmap,nlen,"species:molmap");
  for (n = 0; n < nlen; n++) molmap[n] = 0;

  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    molmap[nint(clusterID[n])-idlo] = 1;
  }

  int *molmapall;
  memory->create(molmapall,nlen,"species:molmapall");
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
    clusterID[n] = molmap[nint(clusterID[n])-idlo]+1;
  }
  memory->destroy(molmap);
  molmap = NULL;

}

/* ---------------------------------------------------------------------- */

void FixSpecies::FindSpecies(int Nmole, int &Nspec)
{
  int inum, *ilist;
  int i, j, k, l, m, n, itype, cid;
  int flag_identity, flag_mol, flag_spec;
  int flag_tmp;
  int *mask =atom->mask;
  int *Nameall, *NMolall;

  memory->destroy(MolName);
  MolName = NULL;
  memory->create(MolName,Nmole*(ntypes+1),"species:MolName");

  memory->destroy(NMol);
  NMol = NULL;
  memory->create(NMol,Nmole,"species:NMol");
  for (m = 0; m < Nmole; m ++)
    NMol[m] = 1;

  memory->create(Nameall,ntypes,"species:Nameall");
  memory->create(NMolall,Nmole,"species:NMolall");

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
  memory->create(nd,Nspec,"species:nd");

  memory->destroy(MolType);
  MolType = NULL;
  memory->create(MolType,Nspec*(ntypes+2),"species:MolType");
}

/* ---------------------------------------------------------------------- */

int FixSpecies::CheckExistence(int id, int ntypes)
{
  int i, j, molid, flag, num;

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

void FixSpecies::WriteFormulas(int Nmole, int Nspec)
{
  int i, j, k, l, jj, itemp;
  int inode;
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

  fprintf(fp,"%11d",ntimestep);
  fprintf(fp,"%11d%11d\t",Nmole,Nspec);

  for (i = 0; i < Nmoltype; i ++)
    fprintf(fp," %d\t",NMol[i]);
  fprintf(fp,"\n");

}

/* ---------------------------------------------------------------------- */

int FixSpecies::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}

/* ---------------------------------------------------------------------- */

int FixSpecies::pack_comm(int n, int *list, double *buf, 
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

void FixSpecies::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) clusterID[i] = nint(buf[m++]);
}

/* ---------------------------------------------------------------------- */

double FixSpecies::memory_usage()
{
  double bytes;

  bytes += nmax*sizeof(double);
  bytes += nmax*nall*sizeof(double);

  return bytes;
}

/* ---------------------------------------------------------------------- */
