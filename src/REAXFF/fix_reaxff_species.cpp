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

/* ----------------------------------------------------------------------
   Contributing authors: Ray Shan (Sandia, tnshan@sandia.gov)
                         Oleg Sergeev (VNIIA, sergeev@vniia.ru)
------------------------------------------------------------------------- */

#include "fix_reaxff_species.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_ave_atom.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include "pair_reaxff.h"
#include "reaxff_defs.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxFFSpecies::FixReaxFFSpecies(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix reaxff/species command");

  force_reneighbor = 0;

  vector_flag = 1;
  size_vector = 2;
  extvector = 0;

  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;

  nvalid = -1;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  ntypes = atom->ntypes;

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  global_freq = nfreq = atoi(arg[5]);

  comm_forward = 4;

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix reaxff/species command");
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Illegal fix reaxff/species command");

  // Neighbor lists must stay unchanged during averaging of bonds,
  // but may be updated when no averaging is performed.

  int rene_flag = 0;
  if (nevery * nrepeat != 1 && (nfreq % neighbor->every != 0 || neighbor->every < nevery * nrepeat)) {
    int newneighborevery = nevery * nrepeat;
    while (nfreq % newneighborevery != 0 && newneighborevery <= nfreq / 2)
      newneighborevery++;

    if (nfreq % newneighborevery != 0)
      newneighborevery = nfreq;

    neighbor->every = newneighborevery;
    rene_flag = 1;
  }

  if (nevery * nrepeat != 1 && (neighbor->delay != 0 || neighbor->dist_check != 0)) {
    neighbor->delay = 0;
    neighbor->dist_check = 0;
    rene_flag = 1;
  }

  if (me == 0 && rene_flag) {
    error->warning(FLERR,"Resetting reneighboring criteria for fix reaxff/species");
  }

  tmparg = nullptr;
  memory->create(tmparg,4,4,"reaxff/species:tmparg");
  strcpy(tmparg[0],arg[3]);
  strcpy(tmparg[1],arg[4]);
  strcpy(tmparg[2],arg[5]);

  if (me == 0) {
    char *suffix = strrchr(arg[6],'.');
    if (suffix && strcmp(suffix,".gz") == 0) {
#ifdef LAMMPS_GZIP
      auto gzip = fmt::format("gzip -6 > {}",arg[6]);
#ifdef _WIN32
      fp = _popen(gzip.c_str(),"wb");
#else
      fp = popen(gzip.c_str(),"w");
#endif
#else
      error->one(FLERR,"Cannot open gzipped file");
#endif
    } else fp = fopen(arg[6],"w");

    if (!fp)
      error->one(FLERR,fmt::format("Cannot open fix reaxff/species file {}: "
                                   "{}",arg[6],utils::getsyserror()));
  }

  x0 = nullptr;
  clusterID = nullptr;

  int ntmp = 1;
  memory->create(x0,ntmp,"reaxff/species:x0");
  memory->create(clusterID,ntmp,"reaxff/species:clusterID");
  vector_atom = clusterID;

  BOCut = nullptr;
  Name = nullptr;
  MolName = nullptr;
  MolType = nullptr;
  NMol = nullptr;
  nd = nullptr;
  molmap = nullptr;

  nmax = 0;
  setupflag = 0;

  // set default bond order cutoff
  int n, i, j, itype, jtype;
  double bo_cut;
  bg_cut = 0.30;
  n = ntypes+1;
  memory->create(BOCut,n,n,"reaxff/species:BOCut");
  for (i = 1; i < n; i ++)
    for (j = 1; j < n; j ++)
      BOCut[i][j] = bg_cut;

  // optional args
  eletype = nullptr;
  ele = filepos = nullptr;
  eleflag = posflag = padflag = 0;

  singlepos_opened = multipos_opened = 0;
  multipos = 0;
  posfreq = 0;

  int iarg = 7;
  while (iarg < narg) {

    // set BO cutoff
    if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix reaxff/species command");
      itype = atoi(arg[iarg+1]);
      jtype = atoi(arg[iarg+2]);
      bo_cut = atof(arg[iarg+3]);
      if (itype > ntypes || jtype > ntypes)
        error->all(FLERR,"Illegal fix reaxff/species command");
      if (itype <= 0 || jtype <= 0)
        error->all(FLERR,"Illegal fix reaxff/species command");
      if (bo_cut > 1.0 || bo_cut < 0.0)
        error->all(FLERR,"Illegal fix reaxff/species command");

      BOCut[itype][jtype] = bo_cut;
      BOCut[jtype][itype] = bo_cut;
      iarg += 4;

      // modify element type names
    } else if (strcmp(arg[iarg],"element") == 0) {
      if (iarg+ntypes+1 > narg) error->all(FLERR,"Illegal fix reaxff/species command");

      eletype = (char**) malloc(ntypes*sizeof(char*));
      int len;
      for (int i = 0; i < ntypes; i ++) {
        len = strlen(arg[iarg+1+i])+1;
        eletype[i] = (char*) malloc(len*sizeof(char));
        strcpy(eletype[i],arg[iarg+1+i]);
      }
      eleflag = 1;
      iarg += ntypes + 1;

      // position of molecules
    } else if (strcmp(arg[iarg],"position") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix reaxff/species command");
      posflag = 1;
      posfreq = atoi(arg[iarg+1]);
      if (posfreq < nfreq || (posfreq%nfreq != 0))
        error->all(FLERR,"Illegal fix reaxff/species command");

      filepos = new char[255];
      strcpy(filepos,arg[iarg+2]);
      if (strchr(filepos,'*')) {
        multipos = 1;
      } else {
        if (me == 0) {
          pos = fopen(filepos, "w");
          if (pos == nullptr) error->one(FLERR,"Cannot open fix reaxff/species position file");
        }
        singlepos_opened = 1;
        multipos = 0;
      }
      iarg += 3;
    } else error->all(FLERR,"Illegal fix reaxff/species command");
  }

  if (!eleflag) {
    memory->create(ele,ntypes+1,"reaxff/species:ele");
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

FixReaxFFSpecies::~FixReaxFFSpecies()
{
  memory->destroy(ele);
  memory->destroy(BOCut);
  memory->destroy(clusterID);
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

int FixReaxFFSpecies::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::setup(int /*vflag*/)
{
  ntotal = static_cast<int> (atom->natoms);
  if (Name == nullptr)
    memory->create(Name,ntypes,"reaxff/species:Name");

  post_integrate();
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix reaxff/species unless atoms have IDs");

  reaxff = (PairReaxFF *) force->pair_match("^reax..",0);
  if (reaxff == nullptr) error->all(FLERR,"Cannot use fix reaxff/species without "
                                "pair_style reaxff, reaxff/kk, or reaxff/omp");

  reaxff->fixspecies_flag = 1;

  // reset next output timestep if not yet set or timestep has been reset
  if (nvalid != update->ntimestep)
    nvalid = update->ntimestep+nfreq;

  // check if this fix has been called twice
  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"reaxff/species") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one fix reaxff/species");

  if (!setupflag) {
    // create a compute to store properties
    modify->add_compute("SPECATOM all SPEC/ATOM q x y z vx vy vz abo01 abo02 abo03 abo04 "
                        "abo05 abo06 abo07 abo08 abo09 abo10 abo11 abo12 abo13 abo14 "
                        "abo15 abo16 abo17 abo18 abo19 abo20 abo21 abo22 abo23 abo24");

    // create a fix to point to fix_ave_atom for averaging stored properties
    auto fixcmd = fmt::format("SPECBOND all ave/atom {} {} {}",tmparg[0],tmparg[1],tmparg[2]);
    for (int i = 1; i < 32; ++i) fixcmd += " c_SPECATOM[" + std::to_string(i) + "]";
    f_SPECBOND = (FixAveAtom *) modify->add_fix(fixcmd);
    setupflag = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::post_integrate()
{
  Output_ReaxFF_Bonds(update->ntimestep,fp);
  if (me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::Output_ReaxFF_Bonds(bigint ntimestep, FILE * /*fp*/)

{
  int Nmole, Nspec;

  // point to fix_ave_atom
  f_SPECBOND->end_of_step();

  if (ntimestep != nvalid) return;

  nlocal = atom->nlocal;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(x0);
    memory->destroy(clusterID);
    memory->create(x0,nmax,"reaxff/species:x0");
    memory->create(clusterID,nmax,"reaxff/species:clusterID");
    vector_atom = clusterID;
  }

  for (int i = 0; i < nmax; i++) {
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

AtomCoord FixReaxFFSpecies::chAnchor(AtomCoord in1, AtomCoord in2)
{
  if (in1.x < in2.x)
    return in1;
  else if (in1.x == in2.x) {
    if (in1.y < in2.y)
      return in1;
    else if (in1.y == in2.y) {
      if (in1.z < in2.z)
        return in1;
    }
  }
  return in2;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::FindMolecule ()
{
  int i,j,ii,jj,inum,itype,jtype,loop,looptot;
  int change,done,anychange;
  int *mask = atom->mask;
  int *ilist;
  double bo_tmp,bo_cut;
  double **spec_atom = f_SPECBOND->array_atom;

  inum = reaxff->list->inum;
  ilist = reaxff->list->ilist;

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
          j = reaxff->tmpid[i][jj];

          if ((j == 0) || (j < i)) continue;
          if (!(mask[j] & groupbit)) continue;

          if (clusterID[i] == clusterID[j]
              && x0[i].x == x0[j].x
              && x0[i].y == x0[j].y
              && x0[i].z == x0[j].z) continue;

          jtype = atom->type[j];
          bo_cut = BOCut[itype][jtype];
          bo_tmp = spec_atom[i][jj+7];

          if (bo_tmp > bo_cut) {
            clusterID[i] = clusterID[j] = MIN(clusterID[i], clusterID[j]);
            x0[i] = x0[j] = chAnchor(x0[i], x0[j]);
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

void FixReaxFFSpecies::SortMolecule(int &Nmole)
{
  memory->destroy(molmap);
  molmap = nullptr;

  int n, idlo, idhi;
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
                   "fix reaxff/species group");
  MPI_Allreduce(&lo,&idlo,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&hi,&idhi,1,MPI_INT,MPI_MAX,world);
  if (idlo == ntotal)
    if (me == 0)
      error->warning(FLERR,"Atom with cluster ID = maxmol "
                     "included in fix reaxff/species group");

  int nlen = idhi - idlo + 1;
  memory->create(molmap,nlen,"reaxff/species:molmap");
  for (n = 0; n < nlen; n++) molmap[n] = 0;

  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    molmap[nint(clusterID[n])-idlo] = 1;
  }

  int *molmapall;
  memory->create(molmapall,nlen,"reaxff/species:molmapall");
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
  molmap = nullptr;

}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::FindSpecies(int Nmole, int &Nspec)
{
  int k, l, m, n, itype, cid;
  int flag_identity, flag_mol, flag_spec;
  int flag_tmp;
  int *mask =atom->mask;
  int *Nameall, *NMolall;

  memory->destroy(MolName);
  MolName = nullptr;
  memory->create(MolName,Nmole*(ntypes+1),"reaxff/species:MolName");

  memory->destroy(NMol);
  NMol = nullptr;
  memory->create(NMol,Nmole,"reaxff/species:NMol");
  for (m = 0; m < Nmole; m ++)
    NMol[m] = 1;

  memory->create(Nameall,ntypes,"reaxff/species:Nameall");
  memory->create(NMolall,Nmole,"reaxff/species:NMolall");

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
  nd = nullptr;
  memory->create(nd,Nspec,"reaxff/species:nd");

  memory->destroy(MolType);
  MolType = nullptr;
  memory->create(MolType,Nspec*(ntypes+2),"reaxff/species:MolType");
}

/* ---------------------------------------------------------------------- */

int FixReaxFFSpecies::CheckExistence(int id, int ntypes)
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

void FixReaxFFSpecies::WriteFormulas(int Nmole, int Nspec)
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

void FixReaxFFSpecies::OpenPos()
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
    if (pos == nullptr) error->one(FLERR,"Cannot open fix reaxff/species position file");
  } else pos = nullptr;
  multipos_opened = 1;

  free(filecurrent);
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::WritePos(int Nmole, int Nspec)
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
    fprintf(pos,"Timestep " BIGINT_FORMAT " NMole %d  NSpec %d  xlo %f  "
            "xhi %f  ylo %f  yhi %f  zlo %f  zhi %f\n",
            update->ntimestep,Nmole, Nspec,
            domain->boxlo[0],domain->boxhi[0],
            domain->boxlo[1],domain->boxhi[1],
            domain->boxlo[2],domain->boxhi[2]);

    fprintf(pos,"ID\tAtom_Count\tType\tAve_q\t\tCoM_x\t\tCoM_y\t\tCoM_z\n");
  }

  Nameall = nullptr;
  memory->create(Nameall,ntypes,"reaxff/species:Nameall");

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

double FixReaxFFSpecies::compute_vector(int n)
{
  if (n == 0)
    return vector_nmole;
  if (n == 1)
    return vector_nspec;
  return 0.0;

}

/* ---------------------------------------------------------------------- */

int FixReaxFFSpecies::nint(const double &r)
{
  int i = 0;
  if (r>0.0) i = static_cast<int>(r+0.5);
  else if (r<0.0) i = static_cast<int>(r-0.5);
  return i;
}

/* ---------------------------------------------------------------------- */

int FixReaxFFSpecies::pack_forward_comm(int n, int *list, double *buf,
                                       int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m] = clusterID[j];
    buf[m+1] = x0[j].x;
    buf[m+2] = x0[j].y;
    buf[m+3] = x0[j].z;
    m += 4;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    clusterID[i] = buf[m];
    x0[i].x = buf[m+1];
    x0[i].y = buf[m+2];
    x0[i].z = buf[m+3];
    m += 4;
  }
}

/* ---------------------------------------------------------------------- */

double FixReaxFFSpecies::memory_usage()
{
  double bytes;

  bytes = 4*nmax*sizeof(double);  // clusterID + x0

  return bytes;
}

/* ---------------------------------------------------------------------- */
