// clang-format off
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
   Contributing author: Greg Wagner (SNL)
------------------------------------------------------------------------- */

#include "pair_meam.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "meam.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <algorithm>
#include <cstring>
#include <memory>

using namespace LAMMPS_NS;

#define MAXLINE 1024

static const int nkeywords = 22;
static const char *keywords[] = {
  "Ec","alpha","rho0","delta","lattce",
  "attrac","repuls","nn2","Cmin","Cmax","rc","delr",
  "augt1","gsmooth_factor","re","ialloy",
  "mixture_ref_t","erose_form","zbl",
  "emb_lin_neg","bkgd_dyn", "theta"};

/* ---------------------------------------------------------------------- */

PairMEAM::PairMEAM(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  allocated = 0;

  nlibelements = 0;

  meam_inst = new MEAM(memory);
  meam_inst->msmeamflag = msmeamflag = 0;
  myname = "meam";

  scale = nullptr;
}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMEAM::~PairMEAM()
{
  if (copymode) return;

  if (meam_inst)
    delete meam_inst;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void PairMEAM::compute(int eflag, int vflag)
{
  int i,ii,n,inum_half,errorflag;
  int *ilist_half,*numneigh_half,**firstneigh_half;
  int *numneigh_full,**firstneigh_full;
  ev_init(eflag,vflag);

  // neighbor list info

  inum_half = listhalf->inum;
  ilist_half = listhalf->ilist;
  numneigh_half = listhalf->numneigh;
  firstneigh_half = listhalf->firstneigh;
  numneigh_full = listfull->numneigh;
  firstneigh_full = listfull->firstneigh;

  // strip neighbor lists of any special bond flags before using with MEAM
  // necessary before doing neigh_f2c and neigh_c2f conversions each step

  if (neighbor->ago == 0) {
    neigh_strip(inum_half,ilist_half,numneigh_half,firstneigh_half);
    neigh_strip(inum_half,ilist_half,numneigh_full,firstneigh_full);
  }

  // check size of scrfcn based on half neighbor list

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  n = 0;
  for (ii = 0; ii < inum_half; ii++) n += numneigh_half[ilist_half[ii]];

  meam_inst->meam_dens_setup(atom->nmax, nall, n);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int ntype = atom->ntypes;

  // 3 stages of MEAM calculation
  // loop over my atoms followed by communication

  int offset = 0;
  errorflag = 0;
  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    meam_inst->meam_dens_init(i,ntype,type,map,x,
                    numneigh_half[i],firstneigh_half[i],
                    numneigh_full[i],firstneigh_full[i],
                    offset);
    offset += numneigh_half[i];
  }
  comm->reverse_comm(this);
  meam_inst->meam_dens_final(nlocal,eflag_either,eflag_global,eflag_atom,
                   &eng_vdwl,eatom,ntype,type,map,scale,errorflag);
  if (errorflag)
    error->one(FLERR,"MEAM library error {}",errorflag);

  comm->forward_comm(this);

  offset = 0;

  // vptr is first value in vatom if it will be used by meam_force()
  // else vatom may not exist, so pass dummy ptr

  double **vptr = nullptr;
  if (vflag_atom) vptr = vatom;
  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    meam_inst->meam_force(i,eflag_global,eflag_atom,vflag_global,
                          vflag_atom,&eng_vdwl,eatom,ntype,type,map,scale,x,
                          numneigh_half[i],firstneigh_half[i],
                          numneigh_full[i],firstneigh_full[i],
                          offset,f,vptr,virial);
    offset += numneigh_half[i];
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairMEAM::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(scale,n+1,n+1,"pair:scale");

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMEAM::settings(int narg, char ** /*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style {} command", myname);

  // set comm size needed by this Pair

  if (msmeamflag) {
    comm_forward = 38+23; // plus 23 for msmeam
    comm_reverse = 30+23; // plus 23 for msmeam
  } else {
    comm_forward = 38;
    comm_reverse = 30;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMEAM::coeff(int narg, char **arg)
{
  int m,n;

  if (!allocated) allocate();

  if (narg < 6) error->all(FLERR,"Incorrect args for pair style {} coefficients", myname);

  // check for presence of first meam file

  std::string lib_file = utils::get_potential_file_path(arg[2]);
  if (lib_file.empty())
    error->all(FLERR,"Cannot open MEAM library file {}",lib_file);

  // find meam parameter file in arguments:
  // first word that is a file or "NULL" after the MEAM library file
  // we need to extract at least one element, so start from index 4

  int paridx=-1;
  std::string par_file;
  for (int i = 4; i < narg; ++i) {
    if (strcmp(arg[i],"NULL") == 0) {
      par_file = "NULL";
      paridx = i;
      break;
    }
    par_file = utils::get_potential_file_path(arg[i]);
    if (!par_file.empty()) {
      paridx=i;
      break;
    }
  }
  if (paridx < 0) error->all(FLERR,"No MEAM parameter file in pair coefficients");
  if ((narg - paridx - 1) != atom->ntypes)
    error->all(FLERR,"Incorrect args for pair style {} coefficients", myname);

  // MEAM element names between 2 filenames
  // nlibelements = # of MEAM elements
  // elements = list of unique element names

  if (nlibelements) {
    libelements.clear();
    mass.clear();
  }

  nlibelements = paridx - 3;
  if (nlibelements < 1) error->all(FLERR,"Incorrect args for pair coefficients");
  if (nlibelements > maxelt)
    error->all(FLERR,"Too many elements extracted from MEAM library (current limit: {}). "
               "Increase 'maxelt' in meam.h and recompile.", maxelt);

  for (int i = 0; i < nlibelements; i++) {
    if (std::any_of(libelements.begin(), libelements.end(),
                    [&](const std::string &elem) { return elem == arg[i+3]; }))
      error->all(FLERR,"Must not extract the same element ({}) from MEAM library twice. ", arg[i+3]);

    libelements.emplace_back(arg[i+3]);
    mass.push_back(0.0);
  }

  // read MEAM library and parameter files
  // pass all parameters to MEAM package
  // tell MEAM package that setup is done

  read_files(lib_file,par_file);
  meam_inst->meam_setup_done(&cutmax);

  // read args that map atom types to MEAM elements
  // map[i] = which element the Ith atom type is, -1 if not mapped

  for (int i = 4 + nlibelements; i < narg; i++) {
    m = i - (4+nlibelements) + 1;
    int j;
    for (j = 0; j < nlibelements; j++)
      if (libelements[j] == arg[i]) break;
    if (j < nlibelements) map[m] = j;
    else if (strcmp(arg[i],"NULL") == 0) map[m] = -1;
    else error->all(FLERR,"Incorrect args for pair style {} coefficients", myname);
  }

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass for i,i in atom class

  int count = 0;
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR,i,mass[map[i]]);
        count++;
      }
      scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair style {} coefficients", myname);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMEAM::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style {} requires newton pair on", myname);

  // need a full and a half neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL)->set_id(1);
  neighbor->add_request(this)->set_id(2);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   half or full
------------------------------------------------------------------------- */

void PairMEAM::init_list(int id, NeighList *ptr)
{
  if (id == 1) listfull = ptr;
  else if (id == 2) listhalf = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMEAM::init_one(int i, int j)
{
  if (setflag[i][j] == 0) scale[i][j] = 1.0;
  scale[j][i] = scale[i][j];
  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairMEAM::read_files(const std::string &globalfile,
                           const std::string &userfile)
{
  read_global_meam_file(globalfile);
  read_user_meam_file(userfile);
}

/* ---------------------------------------------------------------------- */

void PairMEAM::read_global_meam_file(const std::string &globalfile)
{

  // allocate parameter arrays

  std::vector<lattice_t> lat(nlibelements);
  std::vector<int> ielement(nlibelements);
  std::vector<int> ibar(nlibelements);
  std::vector<double> z(nlibelements);
  std::vector<double> atwt(nlibelements);
  std::vector<double> alpha(nlibelements);
  std::vector<double> b0(nlibelements);
  std::vector<double> b1(nlibelements);
  std::vector<double> b2(nlibelements);
  std::vector<double> b3(nlibelements);
  std::vector<double> alat(nlibelements);
  std::vector<double> esub(nlibelements);
  std::vector<double> asub(nlibelements);
  std::vector<double> t0(nlibelements);
  std::vector<double> t1(nlibelements);
  std::vector<double> t2(nlibelements);
  std::vector<double> t3(nlibelements);
  std::vector<double> rozero(nlibelements);
  std::vector<bool> found(nlibelements, false);

  // allocate 6 extra arrays for msmeam

  std::vector<double> b1m(nlibelements);
  std::vector<double> b2m(nlibelements);
  std::vector<double> b3m(nlibelements);
  std::vector<double> t1m(nlibelements);
  std::vector<double> t2m(nlibelements);
  std::vector<double> t3m(nlibelements);

  // open global meamf file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, globalfile, "MEAM", " library");
    char * line;

    constexpr int params_per_line = 19;
    int nset = 0;

    while ((line = reader.next_line(params_per_line))) {
      try {
        ValueTokenizer values(line, "' \t\n\r\f");

        // read each set of params from global MEAM file
        // one set of params can span multiple lines
        // store params if element name is in element list
        // if element name appears multiple times, only store 1st entry
        std::string element = values.next_string();

        // skip if element name isn't in element list

        int index;
        for (index = 0; index < nlibelements; index++)
          if (libelements[index] == element) break;
        if (index == nlibelements) continue;

        // skip if element already appeared (technically error in library file, but always ignored)

        if (found[index]) continue;
        found[index] = true;

        // map lat string to an integer
        std::string lattice_type = values.next_string();

        if (!MEAM::str_to_lat(lattice_type, true, lat[index]))
          error->one(FLERR,"Unrecognized lattice type in MEAM library file: {}", lattice_type);

        // store parameters

        z[index] = values.next_double();
        ielement[index] = values.next_int();
        atwt[index] = values.next_double();
        alpha[index] = values.next_double();
        b0[index] = values.next_double();
        b1[index] = values.next_double();
        b2[index] = values.next_double();
        b3[index] = values.next_double();
        if (msmeamflag) {
          b1m[index] = values.next_double();
          b2m[index] = values.next_double();
          b3m[index] = values.next_double();
        }
        alat[index] = values.next_double();
        esub[index] = values.next_double();
        asub[index] = values.next_double();
        t0[index] = values.next_double();
        t1[index] = values.next_double();
        t2[index] = values.next_double();
        t3[index] = values.next_double();
        if (msmeamflag) {
          t1m[index] = values.next_double();
          t2m[index] = values.next_double();
          t3m[index] = values.next_double();
        }
        rozero[index] = values.next_double();
        ibar[index] = values.next_int();

        if (!isone(t0[index]))
          error->one(FLERR,"Unsupported parameter in MEAM library file: t0 != 1");

        // z given is ignored: if this is mismatched, we definitely won't do what the user said -> fatal error
        if (z[index] != MEAM::get_Zij(lat[index]))
          error->one(FLERR,"Mismatched parameter in MEAM library file: z != lat");

        nset++;
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }
    }

    // error if didn't find all elements in file

    if (nset != nlibelements) {
      std::string msg = "Did not find all elements in MEAM library file, missing:";
      for (int i = 0; i < nlibelements; i++)
        if (!found[i]) {
          msg += " ";
          msg += libelements[i];
        }
      error->one(FLERR,msg);
    }
  }

  // distribute complete parameter sets
  MPI_Bcast(lat.data(), nlibelements, MPI_INT, 0, world);
  MPI_Bcast(ielement.data(), nlibelements, MPI_INT, 0, world);
  MPI_Bcast(ibar.data(), nlibelements, MPI_INT, 0, world);
  MPI_Bcast(z.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(atwt.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(alpha.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(b0.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(b1.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(b2.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(b3.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(alat.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(esub.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(asub.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(t0.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(t1.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(t2.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(t3.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(rozero.data(), nlibelements, MPI_DOUBLE, 0, world);
  // distribute msmeam parameter sets
  MPI_Bcast(b1m.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(b2m.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(b3m.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(t1m.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(t2m.data(), nlibelements, MPI_DOUBLE, 0, world);
  MPI_Bcast(t3m.data(), nlibelements, MPI_DOUBLE, 0, world);

  // pass element parameters to MEAM package

  if (msmeamflag) {
    meam_inst->meam_setup_global(nlibelements, lat.data(), ielement.data(), atwt.data(),
                                alpha.data(), b0.data(), b1.data(), b2.data(), b3.data(),
                                alat.data(), esub.data(), asub.data(), t0.data(), t1.data(),
                                t2.data(), t3.data(), rozero.data(), ibar.data(), b1m.data(),
                                b2m.data(), b3m.data(), t1m.data(), t2m.data(), t3m.data());
  } else {
    meam_inst->meam_setup_global(nlibelements, lat.data(), ielement.data(), atwt.data(),
                                alpha.data(), b0.data(), b1.data(), b2.data(), b3.data(),
                                alat.data(), esub.data(), asub.data(), t0.data(), t1.data(),
                                t2.data(), t3.data(), rozero.data(), ibar.data(), nullptr,
                                nullptr, nullptr, nullptr, nullptr, nullptr);
  }

  // set element masses

  for (int i = 0; i < nlibelements; i++) mass[i] = atwt[i];
}

/* ---------------------------------------------------------------------- */

void PairMEAM::read_user_meam_file(const std::string &userfile)
{
  // done if user param file is "NULL"

  if (userfile == "NULL") return;

  // open user param file on proc 0

  std::shared_ptr<PotentialFileReader> reader;

  if (comm->me == 0) {
    reader = std::make_shared<PotentialFileReader>(lmp, userfile, "MEAM");
  }

  // read settings
  // pass them one at a time to MEAM package
  // match strings to list of corresponding ints
  char * line = nullptr;
  char buffer[MAXLINE];

  while (true) {
    int which;
    int nindex, index[3];
    double value;
    int nline;
    if (comm->me == 0) {
      line = reader->next_line();
      if (line == nullptr) {
        nline = -1;
      } else nline = strlen(line) + 1;
    } else {
      line = buffer;
    }

    MPI_Bcast(&nline,1,MPI_INT,0,world);
    if (nline<0) break;
    MPI_Bcast(line,nline,MPI_CHAR,0,world);

    ValueTokenizer values(line, "=(), '\t\n\r\f");
    int nparams = values.count();
    std::string keyword = values.next_string();

    for (which = 0; which < nkeywords; which++)
      if (keyword == keywords[which]) break;
    if (which == nkeywords)
      error->all(FLERR,"Keyword {} in MEAM parameter file not recognized", keyword);

    nindex = nparams - 2;
    for (int i = 0; i < nindex; i++) index[i] = values.next_int() - 1;

    // map lattce_meam value to an integer
    if (which == 4) {
      std::string lattice_type = values.next_string();
      lattice_t latt;
      if (!MEAM::str_to_lat(lattice_type, false, latt))
        error->all(FLERR, "Unrecognized lattice type in MEAM parameter file: {}", lattice_type);
      value = latt;
    }
    else value = values.next_double();

    // pass single setting to MEAM package

    int errorflag = 0;
    meam_inst->meam_setup_param(which,value,nindex,index,&errorflag);
    if (errorflag) {
      const char *descr[] = { "has an unknown error",
              "is out of range (please report a bug)",
              "expected more indices",
              "has out of range element index"};
      if ((errorflag < 0) || (errorflag > 3)) errorflag = 0;
      error->all(FLERR,"Error in MEAM parameter file: keyword {} {}", keyword, descr[errorflag]);
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairMEAM::pack_forward_comm(int n, int *list, double *buf,
                                int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = meam_inst->rho0[j];
    buf[m++] = meam_inst->rho1[j];
    buf[m++] = meam_inst->rho2[j];
    buf[m++] = meam_inst->rho3[j];
    buf[m++] = meam_inst->frhop[j];
    buf[m++] = meam_inst->gamma[j];
    buf[m++] = meam_inst->dgamma1[j];
    buf[m++] = meam_inst->dgamma2[j];
    buf[m++] = meam_inst->dgamma3[j];
    buf[m++] = meam_inst->arho2b[j];
    buf[m++] = meam_inst->arho1[j][0];
    buf[m++] = meam_inst->arho1[j][1];
    buf[m++] = meam_inst->arho1[j][2];
    buf[m++] = meam_inst->arho2[j][0];
    buf[m++] = meam_inst->arho2[j][1];
    buf[m++] = meam_inst->arho2[j][2];
    buf[m++] = meam_inst->arho2[j][3];
    buf[m++] = meam_inst->arho2[j][4];
    buf[m++] = meam_inst->arho2[j][5];
    for (k = 0; k < 10; k++) buf[m++] = meam_inst->arho3[j][k];
    buf[m++] = meam_inst->arho3b[j][0];
    buf[m++] = meam_inst->arho3b[j][1];
    buf[m++] = meam_inst->arho3b[j][2];
    buf[m++] = meam_inst->t_ave[j][0];
    buf[m++] = meam_inst->t_ave[j][1];
    buf[m++] = meam_inst->t_ave[j][2];
    buf[m++] = meam_inst->tsq_ave[j][0];
    buf[m++] = meam_inst->tsq_ave[j][1];
    buf[m++] = meam_inst->tsq_ave[j][2];
    if (msmeamflag) {
      buf[m++] = meam_inst->arho2mb[j];
      buf[m++] = meam_inst->arho1m[j][0];
      buf[m++] = meam_inst->arho1m[j][1];
      buf[m++] = meam_inst->arho1m[j][2];
      buf[m++] = meam_inst->arho2m[j][0];
      buf[m++] = meam_inst->arho2m[j][1];
      buf[m++] = meam_inst->arho2m[j][2];
      buf[m++] = meam_inst->arho2m[j][3];
      buf[m++] = meam_inst->arho2m[j][4];
      buf[m++] = meam_inst->arho2m[j][5];
      for (k = 0; k < 10; k++) buf[m++] = meam_inst->arho3m[j][k];
      buf[m++] = meam_inst->arho3mb[j][0];
      buf[m++] = meam_inst->arho3mb[j][1];
      buf[m++] = meam_inst->arho3mb[j][2];
    }

  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairMEAM::unpack_forward_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    meam_inst->rho0[i] = buf[m++];
    meam_inst->rho1[i] = buf[m++];
    meam_inst->rho2[i] = buf[m++];
    meam_inst->rho3[i] = buf[m++];
    meam_inst->frhop[i] = buf[m++];
    meam_inst->gamma[i] = buf[m++];
    meam_inst->dgamma1[i] = buf[m++];
    meam_inst->dgamma2[i] = buf[m++];
    meam_inst->dgamma3[i] = buf[m++];
    meam_inst->arho2b[i] = buf[m++];
    meam_inst->arho1[i][0] = buf[m++];
    meam_inst->arho1[i][1] = buf[m++];
    meam_inst->arho1[i][2] = buf[m++];
    meam_inst->arho2[i][0] = buf[m++];
    meam_inst->arho2[i][1] = buf[m++];
    meam_inst->arho2[i][2] = buf[m++];
    meam_inst->arho2[i][3] = buf[m++];
    meam_inst->arho2[i][4] = buf[m++];
    meam_inst->arho2[i][5] = buf[m++];
    for (k = 0; k < 10; k++) meam_inst->arho3[i][k] = buf[m++];
    meam_inst->arho3b[i][0] = buf[m++];
    meam_inst->arho3b[i][1] = buf[m++];
    meam_inst->arho3b[i][2] = buf[m++];
    meam_inst->t_ave[i][0] = buf[m++];
    meam_inst->t_ave[i][1] = buf[m++];
    meam_inst->t_ave[i][2] = buf[m++];
    meam_inst->tsq_ave[i][0] = buf[m++];
    meam_inst->tsq_ave[i][1] = buf[m++];
    meam_inst->tsq_ave[i][2] = buf[m++];
    if (msmeamflag) {
      meam_inst->arho2mb[i] = buf[m++];
      meam_inst->arho1m[i][0] = buf[m++];
      meam_inst->arho1m[i][1] = buf[m++];
      meam_inst->arho1m[i][2] = buf[m++];
      meam_inst->arho2m[i][0] = buf[m++];
      meam_inst->arho2m[i][1] = buf[m++];
      meam_inst->arho2m[i][2] = buf[m++];
      meam_inst->arho2m[i][3] = buf[m++];
      meam_inst->arho2m[i][4] = buf[m++];
      meam_inst->arho2m[i][5] = buf[m++];
      for (k = 0; k < 10; k++) meam_inst->arho3m[i][k] = buf[m++];
      meam_inst->arho3mb[i][0] = buf[m++];
      meam_inst->arho3mb[i][1] = buf[m++];
      meam_inst->arho3mb[i][2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairMEAM::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = meam_inst->rho0[i];
    buf[m++] = meam_inst->arho2b[i];
    buf[m++] = meam_inst->arho1[i][0];
    buf[m++] = meam_inst->arho1[i][1];
    buf[m++] = meam_inst->arho1[i][2];
    buf[m++] = meam_inst->arho2[i][0];
    buf[m++] = meam_inst->arho2[i][1];
    buf[m++] = meam_inst->arho2[i][2];
    buf[m++] = meam_inst->arho2[i][3];
    buf[m++] = meam_inst->arho2[i][4];
    buf[m++] = meam_inst->arho2[i][5];
    for (k = 0; k < 10; k++) buf[m++] = meam_inst->arho3[i][k];
    buf[m++] = meam_inst->arho3b[i][0];
    buf[m++] = meam_inst->arho3b[i][1];
    buf[m++] = meam_inst->arho3b[i][2];
    buf[m++] = meam_inst->t_ave[i][0];
    buf[m++] = meam_inst->t_ave[i][1];
    buf[m++] = meam_inst->t_ave[i][2];
    buf[m++] = meam_inst->tsq_ave[i][0];
    buf[m++] = meam_inst->tsq_ave[i][1];
    buf[m++] = meam_inst->tsq_ave[i][2];
    if (msmeamflag) {
      buf[m++] = meam_inst->arho2mb[i];
      buf[m++] = meam_inst->arho1m[i][0];
      buf[m++] = meam_inst->arho1m[i][1];
      buf[m++] = meam_inst->arho1m[i][2];
      buf[m++] = meam_inst->arho2m[i][0];
      buf[m++] = meam_inst->arho2m[i][1];
      buf[m++] = meam_inst->arho2m[i][2];
      buf[m++] = meam_inst->arho2m[i][3];
      buf[m++] = meam_inst->arho2m[i][4];
      buf[m++] = meam_inst->arho2m[i][5];
      for (k = 0; k < 10; k++) buf[m++] = meam_inst->arho3m[i][k];
      buf[m++] = meam_inst->arho3mb[i][0];
      buf[m++] = meam_inst->arho3mb[i][1];
      buf[m++] = meam_inst->arho3mb[i][2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairMEAM::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    meam_inst->rho0[j] += buf[m++];
    meam_inst->arho2b[j] += buf[m++];
    meam_inst->arho1[j][0] += buf[m++];
    meam_inst->arho1[j][1] += buf[m++];
    meam_inst->arho1[j][2] += buf[m++];
    meam_inst->arho2[j][0] += buf[m++];
    meam_inst->arho2[j][1] += buf[m++];
    meam_inst->arho2[j][2] += buf[m++];
    meam_inst->arho2[j][3] += buf[m++];
    meam_inst->arho2[j][4] += buf[m++];
    meam_inst->arho2[j][5] += buf[m++];
    for (k = 0; k < 10; k++) meam_inst->arho3[j][k] += buf[m++];
    meam_inst->arho3b[j][0] += buf[m++];
    meam_inst->arho3b[j][1] += buf[m++];
    meam_inst->arho3b[j][2] += buf[m++];
    meam_inst->t_ave[j][0] += buf[m++];
    meam_inst->t_ave[j][1] += buf[m++];
    meam_inst->t_ave[j][2] += buf[m++];
    meam_inst->tsq_ave[j][0] += buf[m++];
    meam_inst->tsq_ave[j][1] += buf[m++];
    meam_inst->tsq_ave[j][2] += buf[m++];
    if (msmeamflag) {
      meam_inst->arho2mb[j] += buf[m++];
      meam_inst->arho1m[j][0] += buf[m++];
      meam_inst->arho1m[j][1] += buf[m++];
      meam_inst->arho1m[j][2] += buf[m++];
      meam_inst->arho2m[j][0] += buf[m++];
      meam_inst->arho2m[j][1] += buf[m++];
      meam_inst->arho2m[j][2] += buf[m++];
      meam_inst->arho2m[j][3] += buf[m++];
      meam_inst->arho2m[j][4] += buf[m++];
      meam_inst->arho2m[j][5] += buf[m++];
      for (k = 0; k < 10; k++) meam_inst->arho3m[j][k] += buf[m++];
      meam_inst->arho3mb[j][0] += buf[m++];
      meam_inst->arho3mb[j][1] += buf[m++];
      meam_inst->arho3mb[j][2] += buf[m++];
    }
  }


}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairMEAM::memory_usage()
{
  double bytes = 11 * meam_inst->nmax * sizeof(double);
  bytes += (double)(3 + 6 + 10 + 3 + 3 + 3) * meam_inst->nmax * sizeof(double);
  bytes += (double)3 * meam_inst->maxneigh * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   strip special bond flags from neighbor list entries
   are not used with MEAM
   need to do here so Fortran lib doesn't see them
   done once per reneighbor so that neigh_f2c and neigh_c2f don't see them
------------------------------------------------------------------------- */

void PairMEAM::neigh_strip(int inum, int *ilist,
                           int *numneigh, int **firstneigh)
{
  int i,j,ii,jnum;
  int *jlist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (j = 0; j < jnum; j++) jlist[j] &= NEIGHMASK;
  }
}

/* ---------------------------------------------------------------------- */

void *PairMEAM::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return nullptr;
}
