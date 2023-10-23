/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ray Shan (Sandia, tnshan@sandia.gov)
                         Oleg Sergeev (VNIIA, sergeev@vniia.ru)
                         Jacob Gissinger (NASA, jacob.r.gissinger@gmail.com), 'delete' keyword
------------------------------------------------------------------------- */

#include "fix_reaxff_species.h"

#include "atom.h"
#include "atom_vec.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_ave_atom.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "variable.h"

#include "pair_reaxff.h"
#include "reaxff_defs.h"

#include <algorithm>
#include <cstring>
#include <exception>
#include <random>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_reaxff_species_delete[] =
    "fix reaxff/species, 'delete' keyword: https://doi.org/10.1016/j.carbon.2022.11.002\n\n"
    "@Article{Gissinger23,\n"
    " author = {J. R. Gissinger, S. R. Zavada, J. G. Smith, J. Kemppainen, I. Gallegos, G. M. "
    "Odegard, E. J. Siochi, K. E. Wise},\n"
    " title = {Predicting char yield of high-temperature resins},\n"
    " journal = {Carbon},\n"
    " year =    2023,\n"
    " volume =  202,\n"
    " pages =   {336-347}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixReaxFFSpecies::FixReaxFFSpecies(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), Name(nullptr), MolName(nullptr), NMol(nullptr), nd(nullptr),
    MolType(nullptr), molmap(nullptr), mark(nullptr), Mol2Spec(nullptr), clusterID(nullptr),
    x0(nullptr), BOCut(nullptr), fp(nullptr), pos(nullptr), fdel(nullptr), delete_Tcount(nullptr),
    ele(nullptr), eletype(nullptr), filepos(nullptr), filedel(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix reaxff/species", error);

  force_reneighbor = 1;
  next_reneighbor = -1;

  vector_flag = 1;
  size_vector = 2;
  extvector = 0;

  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;

  compressed = 0;
  nvalid = -1;

  ntypes = atom->ntypes;

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  global_freq = nfreq = utils::inumeric(FLERR, arg[5], false, lmp);
  if (nrepeat == 1) dynamic_group_allow = 1;

  comm_forward = 4;

  if (nevery <= 0) error->all(FLERR, "Invalid fix reaxff/species nevery value {}", nevery);
  if (nrepeat <= 0) error->all(FLERR, "Invalid fix reaxff/species nrepeat value {}", nrepeat);
  if (nfreq <= 0) error->all(FLERR, "Invalid fix reaxff/species nfreq value {}", nfreq);
  if ((nfreq % nevery) || (nrepeat * nevery > nfreq))
    error->all(FLERR, "Incompatible fix reaxff/species nevery/nrepeat/nfreq settings");

  // Neighbor lists must stay unchanged during averaging of bonds,
  // but may be updated when no averaging is performed.

  int rene_flag = 0;
  if (nevery * nrepeat != 1 &&
      (nfreq % neighbor->every != 0 || neighbor->every < nevery * nrepeat)) {
    int newneighborevery = nevery * nrepeat;
    while (nfreq % newneighborevery != 0 && newneighborevery <= nfreq / 2) newneighborevery++;

    if (nfreq % newneighborevery != 0) newneighborevery = nfreq;

    neighbor->every = newneighborevery;
    rene_flag = 1;
  }

  if (nevery * nrepeat != 1 && (neighbor->delay != 0 || neighbor->dist_check != 0)) {
    neighbor->delay = 0;
    neighbor->dist_check = 0;
    rene_flag = 1;
  }

  if (comm->me == 0 && rene_flag)
    error->warning(FLERR,
                   "Resetting reneighboring criteria to 'delay {} every {} check no' "
                   "due to fix reaxff/species averaging of bond data",
                   neighbor->delay, neighbor->every);

  if (comm->me == 0) {
    if (platform::has_compress_extension(arg[6])) {
      fp = platform::compressed_write(arg[6]);
      compressed = 1;
      if (!fp) error->one(FLERR, "Cannot open compressed file");
    } else
      fp = fopen(arg[6], "w");

    if (!fp)
      error->one(FLERR, "Cannot open fix reaxff/species file {}: {}", arg[6], utils::getsyserror());
  }

  x0 = nullptr;
  clusterID = nullptr;

  int ntmp = 1;
  memory->create(x0, ntmp, "reaxff/species:x0");
  memory->create(clusterID, ntmp, "reaxff/species:clusterID");
  vector_atom = clusterID;

  nmax = 0;
  setupflag = 0;

  // set default bond order cutoff
  double bo_cut = 0.30;
  int np1 = ntypes + 1;
  memory->create(BOCut, np1, np1, "reaxff/species:BOCut");
  for (int i = 1; i < np1; i++)
    for (int j = 1; j < np1; j++) BOCut[i][j] = bo_cut;

  // optional args
  eletype = nullptr;
  ele = filepos = filedel = nullptr;
  eleflag = posflag = padflag = 0;
  delflag = specieslistflag = masslimitflag = 0;
  delete_Nlimit = delete_Nsteps = 0;

  singlepos_opened = multipos_opened = del_opened = 0;
  multipos = 0;
  posfreq = 0;

  int iarg = 7;
  while (iarg < narg) {
    // set BO cutoff
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix reaxff/species cutoff", error);
      int ilo, ihi, jlo, jhi;
      utils::bounds(FLERR, arg[iarg + 1], 1, atom->ntypes, ilo, ihi, error);
      utils::bounds(FLERR, arg[iarg + 2], 1, atom->ntypes, jlo, jhi, error);
      bo_cut = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if ((bo_cut > 1.0) || (bo_cut < 0.0))
        error->all(FLERR, "Fix reaxff/species invalid cutoff value: {}", bo_cut);

      for (int i = ilo; i <= ihi; ++i) {
        for (int j = MAX(jlo, i); j <= jhi; ++j) {
          BOCut[i][j] = bo_cut;
          BOCut[j][i] = bo_cut;
        }
      }
      iarg += 4;

      // modify element type names
    } else if (strcmp(arg[iarg], "element") == 0) {
      if (iarg + ntypes + 1 > narg)
        utils::missing_cmd_args(FLERR, "fix reaxff/species element", error);

      eletype = (char **) malloc(ntypes * sizeof(char *));
      int len;
      for (int i = 0; i < ntypes; i++) {
        len = strlen(arg[iarg + 1 + i]) + 1;
        eletype[i] = (char *) malloc(len * sizeof(char));
        strcpy(eletype[i], arg[iarg + 1 + i]);
      }
      eleflag = 1;
      iarg += ntypes + 1;

      // delete species
    } else if (strcmp(arg[iarg], "delete") == 0) {
      delflag = 1;
      delete[] filedel;
      filedel = utils::strdup(arg[iarg + 1]);
      if (comm->me == 0) {
        fdel = fopen(filedel, "w");
        if (!fdel)
          error->one(FLERR, "Cannot open fix reaxff/species delete file {}: {}", filedel,
                     utils::getsyserror());
      }
      del_opened = 1;

      if (strcmp(arg[iarg + 2], "masslimit") == 0) {
        if (iarg + 5 > narg) utils::missing_cmd_args(FLERR, "fix reaxff/species masslimit", error);
        masslimitflag = 1;
        massmin = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
        massmax = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
        iarg += 5;

      } else if (strcmp(arg[iarg + 2], "specieslist") == 0) {
        specieslistflag = 1;
        ndelspec = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
        if (iarg + ndelspec + 4 > narg)
          utils::missing_cmd_args(FLERR, "fix reaxff/species delete specieslist", error);

        del_species.resize(ndelspec);
        for (int i = 0; i < ndelspec; i++) del_species[i] = arg[iarg + 4 + i];

        if (comm->me == 0) {
          fprintf(fdel, "Timestep");
          for (int i = 0; i < ndelspec; i++) fprintf(fdel, "\t%s", del_species[i].c_str());
          fprintf(fdel, "\n");
          fflush(fdel);
        }
        iarg += ndelspec + 4;

      } else
        error->all(FLERR, "Unknown fix reaxff/species delete option: {}", arg[iarg]);
      // rate limit when deleting molecules
    } else if (strcmp(arg[iarg], "delete_rate_limit") == 0) {
      if (iarg + 3 > narg)
        utils::missing_cmd_args(FLERR, "fix reaxff/species delete_rate_limit", error);
      delete_Nlimit_varid = -1;
      if (strncmp(arg[iarg + 1], "v_", 2) == 0) {
        delete_Nlimit_varname = &arg[iarg + 1][2];
        delete_Nlimit_varid = input->variable->find(delete_Nlimit_varname.c_str());
        if (delete_Nlimit_varid < 0)
          error->all(FLERR, "Fix reaxff/species: Variable name {} does not exist",
                     delete_Nlimit_varname);
        if (!input->variable->equalstyle(delete_Nlimit_varid))
          error->all(FLERR, "Fix reaxff/species: Variable {} is not equal-style",
                     delete_Nlimit_varname);
      } else
        delete_Nlimit = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      delete_Nsteps = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      iarg += 3;
      // position of molecules
    } else if (strcmp(arg[iarg], "position") == 0) {
      if (iarg + 3 > narg) utils::missing_cmd_args(FLERR, "fix reaxff/species position", error);
      posflag = 1;
      posfreq = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (posfreq < nfreq || (posfreq % nfreq != 0))
        error->all(FLERR, "Incompatible fix reaxff/species postion frequency {}", posfreq);

      filepos = new char[255];
      strcpy(filepos, arg[iarg + 2]);
      if (strchr(filepos, '*')) {
        multipos = 1;
      } else {
        if (comm->me == 0) {
          pos = fopen(filepos, "w");
          if (pos == nullptr)
            error->one(FLERR, "Cannot open fix reaxff/species position file: {}",
                       utils::getsyserror());
        }
        singlepos_opened = 1;
        multipos = 0;
      }
      iarg += 3;
    } else
      error->all(FLERR, "Unknown fix reaxff/species keyword: {}", arg[iarg]);
  }

  if (!eleflag) {
    memory->create(ele, ntypes + 1, "reaxff/species:ele");
    ele[0] = 'C';
    if (ntypes > 1) ele[1] = 'H';
    if (ntypes > 2) ele[2] = 'O';
    if (ntypes > 3) ele[3] = 'N';
  }

  if (delflag && specieslistflag && masslimitflag)
    error->all(FLERR, "Incompatible combination fix reaxff/species command options");

  if (delete_Nsteps > 0) {
    if (lmp->citeme) lmp->citeme->add(cite_reaxff_species_delete);
    memory->create(delete_Tcount, delete_Nsteps, "reaxff/species:delete_Tcount");

    for (int i = 0; i < delete_Nsteps; i++) delete_Tcount[i] = -1;
    delete_Tcount[0] = 0;
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
  memory->destroy(Mol2Spec);
  memory->destroy(MolType);
  memory->destroy(MolName);
  memory->destroy(delete_Tcount);

  delete[] filepos;
  delete[] filedel;

  if (comm->me == 0) {
    if (compressed)
      platform::pclose(fp);
    else
      fclose(fp);
    if (posflag && multipos_opened) fclose(pos);
    if (fdel) fclose(fdel);
  }

  try {
    modify->delete_compute(fmt::format("SPECATOM_{}", id));
    modify->delete_fix(fmt::format("SPECBOND_{}", id));
  } catch (std::exception &) {
  }
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
  ntotal = static_cast<int>(atom->natoms);
  if (Name == nullptr) memory->create(Name, ntypes, "reaxff/species:Name");

  post_integrate();
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR, "Cannot use fix reaxff/species unless atoms have IDs");

  reaxff = dynamic_cast<PairReaxFF *>(force->pair_match("^reax..", 0));
  if (reaxff == nullptr)
    error->all(FLERR, "Cannot use fix reaxff/species without a reaxff pair_style");

  reaxff->fixspecies_flag = 1;

  // reset next output timestep if not yet set or timestep has been reset
  if (nvalid != update->ntimestep) nvalid = update->ntimestep + nfreq;

  if (!setupflag) {
    // create a compute to store properties
    modify->add_compute(fmt::format("SPECATOM_{} all SPEC/ATOM q x y z vx vy vz abo01 abo02 "
                                    "abo03 abo04 abo05 abo06 abo07 abo08 abo09 abo10 abo11 "
                                    "abo12 abo13 abo14 abo15 abo16 abo17 abo18 abo19 abo20 "
                                    "abo21 abo22 abo23 abo24",
                                    id));

    // create a fix to point to fix_ave_atom for averaging stored properties
    auto fixcmd = fmt::format("SPECBOND_{} all ave/atom {} {} {}", id, nevery, nrepeat, nfreq);
    for (int i = 1; i < 32; ++i) fixcmd += fmt::format(" c_SPECATOM_{}[{}]", id, i);
    f_SPECBOND = dynamic_cast<FixAveAtom *>(modify->add_fix(fixcmd));
    setupflag = 1;
  }

  // check for valid variable name for delete Nlimit keyword
  if (delete_Nsteps > 0) {
    delete_Nlimit_varid = input->variable->find(delete_Nlimit_varname.c_str());
    if (delete_Nlimit_varid < 0)
      error->all(FLERR, "Fix reaxff/species: Variable name {} does not exist",
                 delete_Nlimit_varname);
    if (!input->variable->equalstyle(delete_Nlimit_varid))
      error->all(FLERR, "Fix reaxff/species: Variable {} is not equal-style",
                 delete_Nlimit_varname);
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
  Output_ReaxFF_Bonds(update->ntimestep, fp);
  if (comm->me == 0) fflush(fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::Output_ReaxFF_Bonds(bigint ntimestep, FILE * /*fp*/)

{
  int Nmole, Nspec;

  // point to fix_ave_atom
  f_SPECBOND->end_of_step();

  if (ntimestep != nvalid) {
    // push back delete_Tcount on every step
    if (delete_Nsteps > 0)
      for (int i = delete_Nsteps - 1; i > 0; i--) delete_Tcount[i] = delete_Tcount[i - 1];
    return;
  }

  nlocal = atom->nlocal;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(x0);
    memory->destroy(clusterID);
    memory->create(x0, nmax, "reaxff/species:x0");
    memory->create(clusterID, nmax, "reaxff/species:clusterID");
    vector_atom = clusterID;
  }

  for (int i = 0; i < nmax; i++) { x0[i].x = x0[i].y = x0[i].z = 0.0; }

  Nmole = Nspec = 0;

  FindMolecule();
  SortMolecule(Nmole);
  FindSpecies(Nmole, Nspec);

  vector_nmole = Nmole;
  vector_nspec = Nspec;

  if (comm->me == 0 && ntimestep >= 0) WriteFormulas(Nmole, Nspec);

  if (posflag && ((ntimestep) % posfreq == 0)) {
    WritePos(Nmole, Nspec);
    if (comm->me == 0) fflush(pos);
  }

  if (delflag) DeleteSpecies(Nmole, Nspec);

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
      if (in1.z < in2.z) return in1;
    }
  }
  return in2;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::FindMolecule()
{
  int i, j, ii, jj, inum, itype, jtype, loop, looptot;
  int change, done, anychange;
  int *mask = atom->mask;
  int *ilist;
  double bo_tmp, bo_cut;
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
    } else
      clusterID[i] = 0.0;
  }

  loop = 0;
  while (true) {
    comm->forward_comm(this);
    loop++;

    change = 0;
    while (true) {
      done = 1;

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;

        itype = atom->type[i];

        for (jj = 0; jj < MAXSPECBOND; jj++) {
          j = reaxff->tmpid[i][jj];

          if ((j == 0) || (j < i)) continue;
          if (!(mask[j] & groupbit)) continue;

          if (clusterID[i] == clusterID[j] && x0[i].x == x0[j].x && x0[i].y == x0[j].y &&
              x0[i].z == x0[j].z)
            continue;

          jtype = atom->type[j];
          bo_cut = BOCut[itype][jtype];
          bo_tmp = spec_atom[i][jj + 7];

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
    MPI_Allreduce(&change, &anychange, 1, MPI_INT, MPI_MAX, world);
    if (!anychange) break;

    MPI_Allreduce(&loop, &looptot, 1, MPI_INT, MPI_SUM, world);
    if (looptot >= 400 * comm->nprocs) break;
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::SortMolecule(int &Nmole)
{
  memory->destroy(molmap);
  molmap = nullptr;

  int n, idlo, idhi;
  int *mask = atom->mask;
  int lo = ntotal;
  int hi = -ntotal;
  int flag = 0;
  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    if (clusterID[n] == 0.0) flag = 1;
    lo = MIN(lo, nint(clusterID[n]));
    hi = MAX(hi, nint(clusterID[n]));
  }
  int flagall;
  MPI_Allreduce(&lo, &idlo, 1, MPI_INT, MPI_MIN, world);
  MPI_Allreduce(&hi, &idhi, 1, MPI_INT, MPI_MAX, world);
  int nlen = idhi - idlo + 1;
  if (nlen <= 0) {    // no atoms in group
    Nmole = 0;
    return;
  }
  if (idlo == ntotal)
    if (comm->me == 0)
      error->warning(FLERR, "Atom with cluster ID = maxmol included in fix reaxff/species group {}",
                     group->names[igroup]);

  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall && comm->me == 0)
    error->warning(FLERR, "Atom with cluster ID = 0 included in fix reaxff/species group {}",
                   group->names[igroup]);

  memory->create(molmap, nlen, "reaxff/species:molmap");
  for (n = 0; n < nlen; n++) molmap[n] = 0;

  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    molmap[nint(clusterID[n]) - idlo] = 1;
  }

  int *molmapall;
  memory->create(molmapall, nlen, "reaxff/species:molmapall");
  MPI_Allreduce(molmap, molmapall, nlen, MPI_INT, MPI_MAX, world);

  Nmole = 0;
  for (n = 0; n < nlen; n++) {
    if (molmapall[n])
      molmap[n] = Nmole++;
    else
      molmap[n] = -1;
  }
  memory->destroy(molmapall);

  flag = 0;
  for (n = 0; n < nlocal; n++) {
    if (mask[n] & groupbit) continue;
    if (nint(clusterID[n]) < idlo || nint(clusterID[n]) > idhi) continue;
    if (molmap[nint(clusterID[n]) - idlo] >= 0) flag = 1;
  }

  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
  if (flagall && comm->me == 0) error->warning(FLERR, "One or more cluster has atoms not in group");

  for (n = 0; n < nlocal; n++) {
    if (!(mask[n] & groupbit)) continue;
    clusterID[n] = molmap[nint(clusterID[n]) - idlo] + 1;
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
  int *mask = atom->mask;
  int *Nameall, *NMolall;

  memory->destroy(MolName);
  MolName = nullptr;
  memory->create(MolName, Nmole * (ntypes + 1), "reaxff/species:MolName");

  memory->destroy(NMol);
  NMol = nullptr;
  memory->create(NMol, Nmole, "reaxff/species:NMol");
  for (m = 0; m < Nmole; m++) NMol[m] = 1;

  memory->create(Nameall, ntypes, "reaxff/species:Nameall");
  memory->create(NMolall, Nmole, "reaxff/species:NMolall");

  memory->destroy(Mol2Spec);
  Mol2Spec = nullptr;
  memory->create(Mol2Spec, Nmole, "reaxff/species:Mol2Spec");
  for (m = 0; m < Nmole; m++) Mol2Spec[m] = -1;

  for (m = 1, Nspec = 0; m <= Nmole; m++) {
    for (n = 0; n < ntypes; n++) Name[n] = 0;
    for (n = 0, flag_mol = 0; n < nlocal; n++) {
      if (!(mask[n] & groupbit)) continue;
      cid = nint(clusterID[n]);
      if (cid == m) {
        itype = atom->type[n] - 1;
        Name[itype]++;
        flag_mol = 1;
      }
    }
    MPI_Allreduce(&flag_mol, &flag_tmp, 1, MPI_INT, MPI_MAX, world);
    flag_mol = flag_tmp;

    MPI_Allreduce(Name, Nameall, ntypes, MPI_INT, MPI_SUM, world);
    for (n = 0; n < ntypes; n++) Name[n] = Nameall[n];

    if (flag_mol == 1) {
      flag_identity = 1;
      for (k = 0; k < Nspec; k++) {
        flag_spec = 0;
        for (l = 0; l < ntypes; l++)
          if (MolName[ntypes * k + l] != Name[l]) flag_spec = 1;
        if (flag_spec == 0) {
          NMol[k]++;
          Mol2Spec[m - 1] = k;
        }
        flag_identity *= flag_spec;
      }
      if (Nspec == 0 || flag_identity == 1) {
        for (l = 0; l < ntypes; l++) MolName[ntypes * Nspec + l] = Name[l];
        Mol2Spec[m - 1] = Nspec;
        Nspec++;
      }
    }
  }
  memory->destroy(NMolall);
  memory->destroy(Nameall);

  memory->destroy(nd);
  nd = nullptr;
  memory->create(nd, Nspec, "reaxff/species:nd");

  memory->destroy(MolType);
  MolType = nullptr;
  memory->create(MolType, Nspec * (ntypes + 2), "reaxff/species:MolType");
}

/* ---------------------------------------------------------------------- */

int FixReaxFFSpecies::CheckExistence(int id, int ntypes)
{
  int i, j, molid, flag;

  for (i = 0; i < Nmoltype; i++) {
    flag = 0;
    for (j = 0; j < ntypes; j++) {
      molid = MolType[ntypes * i + j];
      if (molid != MolName[ntypes * id + j]) flag = 1;
    }
    if (flag == 0) return i;
  }
  for (i = 0; i < ntypes; i++) MolType[ntypes * Nmoltype + i] = MolName[ntypes * id + i];

  Nmoltype++;
  return Nmoltype - 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::WriteFormulas(int Nmole, int Nspec)
{
  int i, j, itemp;
  bigint ntimestep = update->ntimestep;

  fprintf(fp, "#  Timestep    No_Moles    No_Specs");

  Nmoltype = 0;

  for (i = 0; i < Nspec; i++) nd[i] = CheckExistence(i, ntypes);

  for (i = 0; i < Nmoltype; i++) {
    std::string molname;
    for (j = 0; j < ntypes; j++) {
      itemp = MolType[ntypes * i + j];
      if (itemp != 0) {
        if (eletype)
          molname += eletype[j];
        else
          molname += ele[j];
        if (itemp != 1) molname += std::to_string(itemp);
      }
    }
    fmt::print(fp, " {:>11}", molname);
  }
  fputs("\n", fp);

  fmt::print(fp, "{:>11} {:>11} {:>11}", ntimestep, Nmole, Nspec);
  for (i = 0; i < Nmoltype; i++) fmt::print(fp, " {:>11}", NMol[i]);
  fputs("\n", fp);
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::OpenPos()
{
  if (comm->me == 0) {
    auto filecurrent = utils::star_subst(filepos, update->ntimestep, padflag);
    pos = fopen(filecurrent.c_str(), "w");
    if (pos == nullptr)
      error->one(FLERR, "Cannot open fix reaxff/species position file {}: {}", filecurrent,
                 utils::getsyserror());
  } else
    pos = nullptr;
  multipos_opened = 1;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::WritePos(int Nmole, int Nspec)
{
  int i, itype, cid;
  int count, count_tmp, m, n, k;
  int *Nameall;
  int *mask = atom->mask;
  double avq, avq_tmp, avx[3], avx_tmp, box[3], halfbox[3];
  double **spec_atom = f_SPECBOND->array_atom;

  if (multipos) OpenPos();

  box[0] = domain->boxhi[0] - domain->boxlo[0];
  box[1] = domain->boxhi[1] - domain->boxlo[1];
  box[2] = domain->boxhi[2] - domain->boxlo[2];

  for (int j = 0; j < 3; j++) halfbox[j] = box[j] / 2;

  if (comm->me == 0) {
    fmt::print(pos,
               "Timestep {} NMole {}  NSpec {}  xlo {:f}  "
               "xhi {:f}  ylo {:f}  yhi {:f}  zlo {:f}  zhi {:f}\n",
               update->ntimestep, Nmole, Nspec, domain->boxlo[0], domain->boxhi[0],
               domain->boxlo[1], domain->boxhi[1], domain->boxlo[2], domain->boxhi[2]);

    fprintf(pos, "ID\tAtom_Count\tType\tAve_q\t\tCoM_x\t\tCoM_y\t\tCoM_z\n");
  }

  Nameall = nullptr;
  memory->create(Nameall, ntypes, "reaxff/species:Nameall");

  for (m = 1; m <= Nmole; m++) {

    count = 0;
    avq = 0.0;
    for (n = 0; n < 3; n++) avx[n] = 0.0;
    for (n = 0; n < ntypes; n++) Name[n] = 0;

    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      cid = nint(clusterID[i]);
      if (cid == m) {
        itype = atom->type[i] - 1;
        Name[itype]++;
        count++;
        avq += spec_atom[i][0];
        if ((x0[i].x - spec_atom[i][1]) > halfbox[0]) spec_atom[i][1] += box[0];
        if ((spec_atom[i][1] - x0[i].x) > halfbox[0]) spec_atom[i][1] -= box[0];
        if ((x0[i].y - spec_atom[i][2]) > halfbox[1]) spec_atom[i][2] += box[1];
        if ((spec_atom[i][2] - x0[i].y) > halfbox[1]) spec_atom[i][2] -= box[1];
        if ((x0[i].z - spec_atom[i][3]) > halfbox[2]) spec_atom[i][3] += box[2];
        if ((spec_atom[i][3] - x0[i].z) > halfbox[2]) spec_atom[i][3] -= box[2];
        for (n = 0; n < 3; n++) avx[n] += spec_atom[i][n + 1];
      }
    }

    avq_tmp = 0.0;
    MPI_Allreduce(&avq, &avq_tmp, 1, MPI_DOUBLE, MPI_SUM, world);
    avq = avq_tmp;

    for (n = 0; n < 3; n++) {
      avx_tmp = 0.0;
      MPI_Reduce(&avx[n], &avx_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, world);
      avx[n] = avx_tmp;
    }

    MPI_Reduce(&count, &count_tmp, 1, MPI_INT, MPI_SUM, 0, world);
    count = count_tmp;

    MPI_Reduce(Name, Nameall, ntypes, MPI_INT, MPI_SUM, 0, world);
    for (n = 0; n < ntypes; n++) Name[n] = Nameall[n];

    if (comm->me == 0) {
      fprintf(pos, "%d\t%d\t", m, count);
      for (n = 0; n < ntypes; n++) {
        if (Name[n] != 0) {
          if (eletype)
            fprintf(pos, "%s", eletype[n]);
          else
            fprintf(pos, "%c", ele[n]);
          if (Name[n] != 1) fprintf(pos, "%d", Name[n]);
        }
      }
      if (count > 0) {
        avq /= count;
        for (k = 0; k < 3; k++) {
          avx[k] /= count;
          if (avx[k] >= domain->boxhi[k]) avx[k] -= box[k];
          if (avx[k] < domain->boxlo[k]) avx[k] += box[k];

          avx[k] -= domain->boxlo[k];
          avx[k] /= box[k];
        }
        fprintf(pos, "\t%.8f \t%.8f \t%.8f \t%.8f", avq, avx[0], avx[1], avx[2]);
      }
      fprintf(pos, "\n");
    }
  }
  if (comm->me == 0 && !multipos) fprintf(pos, "#\n");
  memory->destroy(Nameall);
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::DeleteSpecies(int Nmole, int Nspec)
{
  int ndeletions;
  int headroom = -1;
  if (delete_Nsteps > 0) {
    if (delete_Tcount[delete_Nsteps - 1] == -1) return;
    ndeletions = delete_Tcount[0] - delete_Tcount[delete_Nsteps - 1];
    if (delete_Nlimit_varid > -1)
      delete_Nlimit = input->variable->compute_equal(delete_Nlimit_varid);
    headroom = MAX(0, delete_Nlimit - ndeletions);
    if (headroom == 0) return;
  }

  int i, j, m, n, itype, cid;
  int ndel, ndelone, count, count_tmp;
  int *Nameall;
  int *mask = atom->mask;
  double localmass, totalmass;
  std::string species_str;

  AtomVec *avec = atom->avec;

  mark = nullptr;
  memory->create(mark, nlocal, "reaxff/species:mark");
  for (i = 0; i < nlocal; i++) mark[i] = 0;

  Nameall = nullptr;
  memory->create(Nameall, ntypes, "reaxff/species:Nameall");

  int ndelcomm;
  if (masslimitflag)
    ndelcomm = Nspec;
  else
    ndelcomm = ndelspec;

  double *deletecount;
  memory->create(deletecount, ndelcomm, "reaxff/species:deletecount");
  for (i = 0; i < ndelcomm; i++) deletecount[i] = 0;

  int nmarklist;
  int *marklist;
  memory->create(marklist, nlocal, "reaxff/species:marklist");

  std::random_device rnd;
  std::minstd_rand park_rng(rnd());
  int *molrange;
  memory->create(molrange, Nmole, "reaxff/species:molrange");
  for (m = 0; m < Nmole; m++) molrange[m] = m + 1;
  if (delete_Nsteps > 0) {
    // shuffle index when using rate_limit, in case order is biased
    if (comm->me == 0) std::shuffle(&molrange[0], &molrange[Nmole], park_rng);
    MPI_Bcast(&molrange[0], Nmole, MPI_INT, 0, world);
  }

  int this_delete_Tcount = 0;
  for (int mm = 0; mm < Nmole; mm++) {
    if (this_delete_Tcount == headroom) break;
    m = molrange[mm];
    localmass = totalmass = count = nmarklist = 0;
    for (n = 0; n < ntypes; n++) Name[n] = 0;

    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      cid = nint(clusterID[i]);
      if (cid == m) {
        itype = atom->type[i] - 1;
        Name[itype]++;
        count++;
        marklist[nmarklist++] = i;
        localmass += atom->mass[atom->type[i]];
      }
    }

    MPI_Allreduce(&count, &count_tmp, 1, MPI_INT, MPI_SUM, world);
    count = count_tmp;

    MPI_Allreduce(Name, Nameall, ntypes, MPI_INT, MPI_SUM, world);
    for (n = 0; n < ntypes; n++) Name[n] = Nameall[n];

    MPI_Allreduce(&localmass, &totalmass, 1, MPI_DOUBLE, MPI_SUM, world);

    species_str = "";
    for (j = 0; j < ntypes; j++) {
      if (Name[j] != 0) {
        if (eletype)
          species_str += eletype[j];
        else
          species_str += ele[j];
        if (Name[j] != 1) species_str += fmt::format("{}", Name[j]);
      }
    }

    if (masslimitflag) {

      // find corresponding moltype

      if (totalmass > massmin && totalmass < massmax) {
        this_delete_Tcount++;
        for (j = 0; j < nmarklist; j++) {
          mark[marklist[j]] = 1;
          deletecount[Mol2Spec[m - 1]] += 1.0 / (double) count;
        }
      }
    } else {
      if (count > 0) {
        for (i = 0; i < ndelspec; i++) {
          if (del_species[i] == species_str) {
            this_delete_Tcount++;
            for (j = 0; j < nmarklist; j++) {
              mark[marklist[j]] = 1;
              deletecount[i] += 1.0 / (double) count;
            }
            break;
          }
        }
      }
    }
  }

  // delete atoms. loop in reverse order to avoid copying marked atoms

  ndel = ndelone = 0;
  for (i = atom->nlocal - 1; i >= 0; i--) {
    if (mark[i] == 1) {
      avec->copy(atom->nlocal - 1, i, 1);
      atom->nlocal--;
      ndelone++;
    }
  }

  MPI_Allreduce(&ndelone, &ndel, 1, MPI_INT, MPI_SUM, world);

  atom->natoms -= ndel;

  if (comm->me == 0)
    MPI_Reduce(MPI_IN_PLACE, deletecount, ndelcomm, MPI_DOUBLE, MPI_SUM, 0, world);
  else
    MPI_Reduce(deletecount, deletecount, ndelcomm, MPI_DOUBLE, MPI_SUM, 0, world);

  if (comm->me == 0) {
    if (masslimitflag) {
      int printflag = 0;
      for (int m = 0; m < Nspec; m++) {
        if (deletecount[m] > 0) {
          if (printflag == 0) {
            fmt::print(fdel, "Timestep {}", update->ntimestep);
            printflag = 1;
          }
          fprintf(fdel, " %g ", deletecount[m]);
          for (j = 0; j < ntypes; j++) {
            int itemp = MolName[ntypes * m + j];
            if (itemp != 0) {
              if (eletype)
                fprintf(fdel, "%s", eletype[j]);
              else
                fprintf(fdel, "%c", ele[j]);
              if (itemp != 1) fprintf(fdel, "%d", itemp);
            }
          }
        }
      }
      if (printflag) {
        fprintf(fdel, "\n");
        fflush(fdel);
      }
    } else {
      int writeflag = 0;
      for (i = 0; i < ndelspec; i++)
        if (deletecount[i]) writeflag = 1;

      if (writeflag) {
        fmt::print(fdel, "{}", update->ntimestep);
        for (i = 0; i < ndelspec; i++) { fprintf(fdel, "\t%g", deletecount[i]); }
        fprintf(fdel, "\n");
        fflush(fdel);
      }
    }
  }

  // push back delete_Tcount on every step
  if (delete_Nsteps > 0) {
    for (i = delete_Nsteps - 1; i > 0; i--) delete_Tcount[i] = delete_Tcount[i - 1];
    delete_Tcount[0] += this_delete_Tcount;
  }

  if (ndel && (atom->map_style != Atom::MAP_NONE)) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  next_reneighbor = update->ntimestep;

  memory->destroy(Nameall);
  memory->destroy(marklist);
  memory->destroy(mark);
  memory->destroy(deletecount);
  memory->destroy(molrange);
}

/* ---------------------------------------------------------------------- */

double FixReaxFFSpecies::compute_vector(int n)
{
  if (n == 0) return vector_nmole;
  if (n == 1) return vector_nspec;
  return 0.0;
}

/* ---------------------------------------------------------------------- */

int FixReaxFFSpecies::nint(const double &r)
{
  int i = 0;
  if (r > 0.0)
    i = static_cast<int>(r + 0.5);
  else if (r < 0.0)
    i = static_cast<int>(r - 0.5);
  return i;
}

/* ---------------------------------------------------------------------- */

int FixReaxFFSpecies::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m] = clusterID[j];
    buf[m + 1] = x0[j].x;
    buf[m + 2] = x0[j].y;
    buf[m + 3] = x0[j].z;
    m += 4;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFSpecies::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    clusterID[i] = buf[m];
    x0[i].x = buf[m + 1];
    x0[i].y = buf[m + 2];
    x0[i].z = buf[m + 3];
    m += 4;
  }
}

/* ---------------------------------------------------------------------- */

double FixReaxFFSpecies::memory_usage()
{
  double bytes;

  bytes = 4 * nmax * sizeof(double);    // clusterID + x0

  return bytes;
}

/* ---------------------------------------------------------------------- */
