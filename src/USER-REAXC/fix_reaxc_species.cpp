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
#include "string.h"
#include "fix_reaxc_species.h"
#include "atom.h"
#include "pair_reax_c.h"
#include "update.h"
#include "force.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "reaxc_list.h"

using namespace LAMMPS_NS;

#define Carbon 0
#define Hydrogen 1
#define Nitrogen 3
#define Oxygen 2

#define write_BL

/* ---------------------------------------------------------------------- */

FixReaxCSpecies::FixReaxCSpecies(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix reax/c/species command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  nmax = (atom->nmax+10) * nprocs;
  ntypes = atom->ntypes;

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  global_freq = nfreq = atoi(arg[5]);
  
  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix reax/c/species command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix reax/c/species command");

  if (me == 0) {
    fp = fopen(arg[6],"w");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open fix reax/c/species file %s",arg[6]);
      error->one(FLERR,str);
    }
  }
  multifile = padflag = bondflag = 0;

  // optional args
  mass = bond = read = NULL;
  eletype = NULL;

  int iarg = 7;
  while (iarg < narg) {

    // mass spectra 
    if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) 
	error->all(FLERR,"Illegal fix reax/c/species mass command");

      int n = strlen(arg[iarg+1]) + 1;
      filename = new char[n];
      strcpy(filename,arg[iarg+1]);
      if (strchr(filename,'*')) multifile = 1;

      OpenFile();
      iarg += 2;

    // reax/c/bond info
    } else if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+2 > narg) 
	error->all(FLERR,"Illegal fix reax/c/species bond command");
      
      if (me == 0) {
        bond = fopen(arg[iarg+1],"w");
	if (bond == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix reax/c/species bond file %s",
		  arg[iarg+1]);
          error->one(FLERR,str);
	}
      }
      bondflag = 1;
      iarg += 2;

    // read Cutoff.dic
    } else if (strcmp(arg[iarg],"read") == 0) {
      if (iarg+2 > narg)
	error->all(FLERR,"Illegal fix reax/c/species read command");
      
      if (me == 0) {
        read = fopen(arg[iarg+1],"r");
	if (read == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix reax/c/species read file %s",
		  arg[iarg+1]);
          error->one(FLERR,str);
	}
	fprintf(stderr,
		"	Fix Reax/C/Species: reading BO cutoffs from %s\n",
		arg[iarg+1]);
	
	ReadDict(read);
	fclose(read);
      }
      iarg += 2;

    // modify element type names
    } else if (strcmp(arg[iarg],"element") == 0) {
      if (iarg+ntypes+1 > narg) 
	error->all(FLERR,"Illegal fix reax/c/species element command");

      eletype = new char*[ntypes];
      for (int i = 0; i < ntypes; i ++) {
	eletype[i] = new char[2];
	strcpy(eletype[i],arg[iarg+1+i]);
      }
      
      iarg += ntypes + 1;

    } else error->all(FLERR,"Illegal fix reax/c/species command");
  }

  // allocate memory for averaging BO
  allocate();

  // nvalid = next step on which end_of_step does something
  nvalid = nextvalid();

#if defined(write_BL)
  blfp=fopen("BLout","w");
#endif
 
}

/* ---------------------------------------------------------------------- */

FixReaxCSpecies::~FixReaxCSpecies()
{
  memory->destroy(BO);
  memory->destroy(ele);
  memory->destroy(Mol);
  memory->destroy(tag);
  memory->destroy(NMol);
  memory->destroy(Name);
  memory->destroy(type);
  memory->destroy(MolID);
  memory->destroy(MolType);
  memory->destroy(MolName);
  memory->destroy(nd);
  memory->destroy(numcount);

  if(mass) {
    memory->destroy(MolMass);
    memory->destroy(masstype);
  }

  memory->destroy(sbo);
  memory->destroy(nlp);
  memory->destroy(avq);

  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixReaxCSpecies::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   Only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixReaxCSpecies::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::init()
{
  // ensure ReaxFF/C is defined
  reaxc = (PairReaxC *) force->pair_match("reax/c",1);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/species without "
		  "pair_style reax/c");

  // Notify pair_reax_c to calculation bonding information
  reaxc->fixspecies_flag = 1;

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
  }

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::end_of_step()
{

  // skip if not step that requires calculating bonds
  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  OutputReaxCBonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::ReadDict(FILE *read)
{
  int i, j, iatom, jatom, c, d;
  double cut;
  char buffer[BUFLEN];
  
  BOCut = NULL;
  BOCut = (double**) malloc((ntypes+2)*sizeof(double*));
  for (int m = 0; m <= ntypes; m ++)
    BOCut[m] = (double*) malloc((ntypes+2)*sizeof(double));

  for (i = 0; i < ntypes; i ++)
    for (j = 0; j < ntypes;j ++)
      BOCut[i][j] = 0.30;
  d = 0;

  while(fscanf(read,"%s",buffer)!=EOF) {
    if (strcmp(buffer,"#") == 0) {
      while ((c = getc(read)) != EOF)
      if (c == '\n') break;
      continue;
    }
    iatom = atoi(buffer);
    fscanf(read,"%s",buffer);
    jatom = atoi(buffer);
    fscanf(read,"%lf",&cut);
    BOCut[iatom][jatom]=cut;
    BOCut[jatom][iatom]=cut;
    d ++;
  }
  /*
  for(i=0;i<ntypes;i++) 
    for(j=0;j<ntypes;j++)
      fprintf(stderr," %d - %d: %f  (%d)\n",i, j, BOCut[i][j], d);
  */

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::OpenFile()
{
  char *filecurrent;
  if (multifile == 0) filecurrent = filename;
  else {
    filecurrent = new char[strlen(filename) + 16];
    char *ptr = strchr(filename,'*');
    *ptr = '\0';
    if (padflag == 0) 
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
		filename,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filename,update->ntimestep,ptr+1);
    }
    *ptr = '*';
  }
  if (me == 0) {
    mass = fopen(filecurrent, "w");
    if (mass == NULL) 
      error->one(FLERR,"Cannot open fix reax/c/species mass file");
  }
  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::OutputReaxCBonds(bigint ntimestep, FILE *fp)

{
  int i, j, k, itype, jtype, iatom, itag, jtag;
  int numbonds,nsbmax, Nmole, Nspec, count;
  int nlocal_max, nsbmax_max, Nmole_sum, Nspec_sum, count_sum;
  int b, nbuf, nbuf_local, inode;
  double *buf, BO_tmp;
 
  int nlocal = atom->nlocal;
  int nlocal_tot = static_cast<int> (atom->natoms);
  // int nmax = (nlocal+2)*nprocs;
  // nmax = (atom->nmax+10) * nprocs;
  double repeat = nrepeat;

  Nmole = Nspec = count = Nmole_sum = Nspec_sum = 0;

  // zero out average BO for next Nfreq
  if (irepeat == 0)
    for (i = 0; i < nmax; i++) {
      for (j = 0; j < nmax; j++) BO[i][j] = 0.0;
      sbo[i] = nlp[i] = avq[i] = 0.0;
    }

  // get maxval from all nodes
  nsbmax = reaxc->system->my_bonds;	// max bond for each atom
  MPI_Allreduce(&nlocal,&nlocal_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nsbmax,&nsbmax_max,1,MPI_INT,MPI_MAX,world);
 
  // allocate a temporary buffer for the snapshot info
  nbuf = 1+(2*nsbmax_max+10)*nlocal_max*nprocs;
  memory->create(buf,nbuf,"reax/c/species:buf");
  if (irepeat == 0)
    for (i = 0; i < nbuf; i ++) buf[i] = 0.0;

  // Pass information to buffer
  PassBuffer( system, buf, nbuf_local);

  // Receive information from buffer
  RecvBuffer( system, buf, nbuf, nbuf_local);

  // done if irepeat < nrepeat, else reset irepeat and nvalid
  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    return;
  }
  irepeat = 0;
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;

  // set arrays to initial values
  for (i = 0; i < nmax; i++) {
    MolID[i] = -1;
    NMol[i] = 1;
    for (j = 0; j < nmax; j++) {
      BO[i][j] /= repeat;
    }
    sbo[i] /= repeat;
    nlp[i] /= repeat;
    avq[i] /= repeat;
  }

  if (multifile) OpenFile();

#if defined(write_BL)
  for (int n = 0; n <= ntypes; n ++)
    for (int m = 0; m <= ntypes;m ++) {
      BLcount[n][m] = 0;
      BLsum[n][m] = 0.0;
    }
#endif

  if (me == 0) {

    // find bonds from averaged BO values
    FindBond( system, nlocal_tot);

#if defined(write_BL)
  fprintf(blfp," # Timestep %d \n",ntimestep);
  for (int n = 1; n <= ntypes; n ++)
    for (int m = 1; m <= ntypes;m ++) {
      fprintf(blfp," %d - %d: %f \n",n,m,BLsum[n][m]/BLcount[n][m]);
    }
  fprintf(blfp," # \n");
#endif

    // find molecules from found bonds
    FindMolecule( system, nlocal_tot, count);

    // find species from found molecules
    FindSpecies( system, nlocal_tot, ntypes, count, Nmole, Nspec);

    // find chemical formulae of found species
    FindFormulas( ntypes, Nmole, Nspec);

    // optional: plot mass spectra
    if (mass) PlotSpectra( system, nlocal_tot, ntypes, Nmole, Nspec);

    // optional: output bonding info in same format as reax/c/bonds
    if (bondflag) WriteBond( system, nlocal_tot, nsbmax_max, buf, nbuf);

  }
  memory->destroy(buf);
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::PassBuffer( reax_system *system, double *buf, 
		int &nbuf_local)
{
  int i, j, k, numbonds;
  int nlocal = atom->nlocal;

  j = 2;
  buf[0] = nlocal;
  for (i = 0; i < nlocal; i++) {
    buf[j-1] = atom->tag[i];
    buf[j+0] = atom->type[i];
    buf[j+1] = reaxc->system->my_atoms[i].numbonds; 
    buf[j+2] = reaxc->workspace->total_bond_order[i];
    buf[j+3] = reaxc->workspace->nlp[i];
    buf[j+4] = atom->q[i];
    numbonds = nint(buf[j+1]);

    for (k = 5; k < 5+numbonds; k++) {
      buf[j+k] = reaxc->system->my_atoms[i].nbr_id[k-4];
    }
    j += (5+numbonds);
    for (k = 0; k < numbonds; k++) {
      buf[j+k] = reaxc->system->my_atoms[i].nbr_bo[k+1];
    }
    j += (1+numbonds);
  }
  nbuf_local = j - 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::RecvBuffer( reax_system *system, double *buf, 
		int nbuf, int nbuf_local)
{
  int i, j, k, itype, itag, jtag;
  int inode, nlocal_tmp, count, numbonds;
  MPI_Request irequest, irequest2;
  MPI_Status istatus;
  int nlocal = atom->nlocal;

  j = 2;
  if (me == 0) {
    for (inode = 0; inode < nprocs; inode ++) {
      if (inode == 0) {
	nlocal_tmp = nlocal;
      } else {
	MPI_Irecv(&buf[0],nbuf,MPI_DOUBLE,inode,0,world,&irequest);
        MPI_Wait(&irequest,&istatus);
	nlocal_tmp = nint(buf[0]);
      }
      j = 2;
      for (i = 0; i < nlocal_tmp; i ++) {
	itag = nint(buf[j-1]);
	itype = type[itag] = nint(buf[j+0]);
	numbonds = nint(buf[j+1]);
	sbo[itag] += buf[j+2];
	nlp[itag] += buf[j+3];
	avq[itag] += buf[j+4];

	for (k = 5; k < 5+numbonds; k++) {
          reaxc->system->my_atoms[i].nbr_id[k-4] = nint(buf[j+k]);
	}
	j += (5+numbonds);
	for (k = 0; k < numbonds; k++) {
          jtag = reaxc->system->my_atoms[i].nbr_id[k+1];
	  BO[itag][jtag] += buf[j+k];
	}
	j += (1+numbonds);
      }
    }
  } else {
    // MPI_Rsend(&buf[0],nbuf_local,MPI_DOUBLE,0,0,world);
    MPI_Isend(&buf[0],nbuf_local,MPI_DOUBLE,0,0,world,&irequest2);
    MPI_Wait(&irequest2,&istatus);
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::FindBond( reax_system *system, int nlocal)
{
  int i, j, itag, jtag, itype, jtype, num;
  double BO_tmp, bg_cut;

#if defined(write_BL)
  double **x = atom->x;
  double r, rsq, delx, dely, delz;
#endif

  bg_cut = 0.3;  // default value if not reading from Cutoff.dic

  for (i = 1; i <= nlocal; i ++) {
    itag = atom->tag[i];
    itype = type[i];
    num = 0;
    for (j = 1; j <= nlocal; j ++) {
      jtag = atom->tag[j];
      jtype = type[j];
      BO_tmp = BO[i][j];
      if (read) bg_cut = BOCut[itype][jtype];
      if (BO_tmp >= bg_cut ) {
	num ++;
	reaxc->system->my_atoms[i].nbr_id[num] = j;
		// fprintf(stderr,"%d(%d) - %d(%d): %f cut(%f)\n",i,itype,j,jtype,BO_tmp,bg_cut);
#if defined(write_BL)
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r = sqrt(rsq);
      	 // fprintf(stderr,"%d-%d: %f \n",itype,jtype,r);
      BLsum[itype][jtype] += r;
      BLcount[itype][jtype] ++;
#endif
      }
    }
    reaxc->system->my_atoms[i].numbonds = num;
    	// fprintf(stderr," i: %d, num: %d \n",i,num);
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::FindMolecule(reax_system *system, int nlocal, int &count)
{
  int i, j, k, itag, iatom, jatom, molid, num;

  for (i = 1; i <= nlocal; i ++) {
    iatom = i;
    if(MolID[iatom] == -1) {
      MolID[iatom] = count;
      count++;
    }
    num = reaxc->system->my_atoms[iatom].numbonds;
    
    for (j = 1; j <= num; j ++) {
      jatom = reaxc->system->my_atoms[iatom].nbr_id[j];
      molid = MolID[jatom];
      if(molid != MolID[iatom]) {
        if(MolID[jatom] == -1) MolID[jatom] = MolID[iatom];
        else {
          for(k = 1; k <= nlocal; k ++)
            if(MolID[k] == molid) MolID[k] = MolID[iatom];
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::FindSpecies( reax_system *system, int nlocal,
		int ntypes, int count, int &Nmole, int &Nspec)
{
  int i, j, k, l, m, n, jtype;
  int flag_identity, flag_mol, flag_spec;

  for (i = 0, Nspec = 0, Nmole = 0; i < count; i ++) {
    for (j = 0; j < ntypes; j ++) Name[j] = 0;
    for (j = 1, flag_mol = 0; j <= nlocal; j ++) {
      if (MolID[j] == i) {
        jtype = type[j]-1;
        Name[jtype] ++;
        flag_mol = 1; 
      }
    }
    if (flag_mol == 1) {
      flag_identity = 1;
      for (k = 0; k < Nspec; k ++) {
        flag_spec=0;
        for (l = 0; l < ntypes; l ++) {
          if (MolName[ntypes*k+l] != Name[l]) flag_spec = 1;
	}
        if (flag_spec == 0) NMol[k] ++;  
        flag_identity *= flag_spec; 
      }
      if (Nspec == 0 || flag_identity == 1) {
        for (l = 0; l < ntypes; l ++) {
          MolName[ntypes*Nspec+l] = Name[l];
	}
        Nspec ++; 
      }
      Nmole ++;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixReaxCSpecies::CheckExistence(int id, int ntypes)
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

void FixReaxCSpecies::FindFormulas(int ntypes, int Nmole, int Nspec)
{
  int i, j, id, itemp, step, molnum, num;
  bigint ntimestep = update->ntimestep;
  
  fprintf(fp,"# Timestep     No_Moles     No_Specs     ");

  Nmoltype = 0;
  for (i = 0; i < Nspec; i ++) {
    nd[i] = CheckExistence(i, ntypes);
  }

  for (i = 0; i < Nmoltype; i ++) {
    for (j = 0;j < ntypes; j ++) {
      itemp = MolType[ntypes * i + j];
      if( itemp != 0) {
        if(eletype) fprintf(fp,"%s",eletype[j]);
        else fprintf(fp,"%c",ele[j]);
        if (itemp != 1) {
          fprintf(fp,"%d",itemp);
        }
      }
    }
    fprintf(fp,"\t");
    numcount[i] = 0;
  }
  fprintf(fp,"\n");

  fprintf(fp,"%11d", ntimestep);
  fprintf(fp,"%11d%11d\t",Nmole,Nspec);

  for (i = 0; i < Nmoltype; i ++) numcount[i]=0;
  for (i = 0; i < Nspec; i ++) {
    numcount[nd[i]]=NMol[i];
  }
  for (i = 0; i < Nmoltype; i ++) {
    fprintf(fp," %d\t",numcount[i]);
  }

  fprintf(fp,"\n");
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::PlotSpectra(reax_system *system, int natoms, 
		int ntypes, int Nmole, int Nspec)
{
  int i, j, id, itemp, step, molnum, num;

  for (i = 0; i < Nmole; i ++) MolMass[i] = 0.0;

  for (i = 0; i < Nspec; i ++) {
    nd[i] = CheckExistence(i, ntypes);
  }

  for (i = 0; i < Nmoltype; i ++) {
    for (j = 0;j < ntypes; j ++) {
      itemp = MolType[ntypes * i + j];
      if( itemp != 0) {
        if (j == Carbon) {
	  masstype[j] = 12.01070;
        }
        else if (j == Hydrogen) {
	  masstype[j] = 1.007940;
        }
        else if (j == Oxygen) {
	  masstype[j] = 15.99940;
        }
        else if (j == Nitrogen) {
	  masstype[j] = 14.00670;
        }
	MolMass[i] += masstype[j];
        if (itemp != 1) {
	  MolMass[i] += itemp * masstype[j];
        }
      }
    }
    numcount[i] = 0;
  }

  fprintf(mass,"# Timestep: %d\t No_Molecules: %d\t No_Species: %d\n",
		  update->ntimestep,Nmole,Nspec);

  for (i = 0; i < Nmoltype; i ++) numcount[i]=0;
  for (i = 0; i < Nspec; i ++) {
    numcount[nd[i]]=NMol[i];
  }

  for (i = 0; i < Nmoltype; i ++) {
    fprintf(mass," %f   %d\n",MolMass[i],numcount[i]);
  }

}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::WriteBond( reax_system *system, int natoms, 
		int maxnumbonds, double *buf, int nbuf)
{
  int i, j, k, n;
  int itag, jtag, itype, numbonds, nbuf_local, molid;
  
  int nlocal = atom->nlocal;
  int ntimestep = update->ntimestep;
  double cutof3 = reaxc->control->bg_cut;

  fprintf(bond,"# Timestep " BIGINT_FORMAT " \n",ntimestep);
  fprintf(bond,"# \n");
  fprintf(bond,"# Number of particles %d \n",natoms);
  fprintf(bond,"# \n");
  fprintf(bond,"# Max number of bonds per atom %d with "
	    "coarse bond order cutoff %5.3f \n",
	    maxnumbonds,cutof3);
  fprintf(bond,"# Particle connection table and bond orders \n");
  fprintf(bond,"# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q \n");

  for (i = 1; i <= natoms; i ++) {
    
    // print atom tag, atom type, no.bonds
    itag = atom->tag[i-1];
    itype = type[i];
    numbonds = reaxc->system->my_atoms[i].numbonds;
    fprintf(bond," %d %d %d",i,itype,numbonds);
    
    // print connection table
    for (k = 2; k < 2+numbonds; k ++) {
      jtag = reaxc->system->my_atoms[i].nbr_id[k-1];
      fprintf(bond," %d",jtag);
    }

    // molecule id
    if (atom->molecule == NULL )
      molid = 0;
    else
      molid = atom->molecule[i];
    fprintf(bond," %d",molid);
    
    // print bond orders
    for (k = 0; k < numbonds; k ++) {
      jtag = reaxc->system->my_atoms[i].nbr_id[k+1];
      fprintf(bond,"%14.3f",BO[i][jtag]);
    }
    
    // print sum of bond orders, no. of lone pairs, charge
    fprintf(bond,"%14.3f%14.3f%14.3f\n",sbo[i],nlp[i],avq[i]);
  }
  fprintf(bond,"# \n");
 
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
  return nvalid;
}

/* ---------------------------------------------------------------------- */

void FixReaxCSpecies::allocate()
{

  BO = NULL;
  ele = NULL;
  Mol = NULL;
  tag = NULL;
  NMol = NULL;
  Name = NULL;
  type = NULL;
  MolID = NULL;
  MolName = NULL;
  MolType = NULL;
  nd = NULL;
  numcount = NULL;
  MolMass = NULL;
  masstype = NULL;

  ele = new char[4];
  ele[0]='C';
  ele[1]='H';
  ele[2]='O';
  ele[3]='N';

  Mol = (int*) malloc(nmax*sizeof(int));
  tag = (int*) malloc(nmax*sizeof(int));
  NMol = (int*) malloc(nmax*sizeof(int));
  Name = (int*) malloc(ntypes*sizeof(int));
  type = (int*) malloc(nmax*sizeof(int));
  MolID = (int*) malloc(nmax*sizeof(int));
  MolName = (int*) malloc(nmax*ntypes*sizeof(int));
  MolType = (int*) malloc((ntypes+20)*MaxMolTypes*sizeof(int));
  nd = (int*) malloc(nmax*sizeof(int));
  numcount = (int*) malloc(nmax*sizeof(int));

  BO = (double**) malloc(nmax*sizeof(double*));
  for (int m = 0; m < nmax; m ++)
    BO[m] = (double*) malloc(nmax*sizeof(double));

  irepeat = Nmoltype = 0;
  for (int n = 0; n < nmax; n++) {
    NMol[n] = 1;
    for (int m = 0; m < nmax; m++)
      BO[n][m] = 0.0;
  }

  if (mass) {
    MolMass = (double*) malloc(nmax*sizeof(double));
    for (int n = 0; n < nmax; n++) MolMass[n] = 0.0;

    masstype = (double*) malloc(ntypes*sizeof(double));
  }

  sbo = (double*) malloc(nmax*sizeof(double));
  nlp = (double*) malloc(nmax*sizeof(double));
  avq = (double*) malloc(nmax*sizeof(double));

#if defined(write_BL)
  BLcount = NULL;
  BLsum = NULL;
  BLcount = (int**) malloc((ntypes+2)*sizeof(int*));
  BLsum = (double**) malloc((ntypes+2)*sizeof(double*));
  for (int m = 0; m <= ntypes; m ++) {
    BLcount[m] = (int*) malloc((ntypes+2)*sizeof(int));
    BLsum[m] = (double*) malloc((ntypes+2)*sizeof(double));
  }
  for (int n = 1; n <= ntypes; n ++)
    for (int m = 1; m <= ntypes;m ++) {
      BLcount[n][m] = 0;
      BLsum[n][m] = 0.0;
    }
#endif
 
}
