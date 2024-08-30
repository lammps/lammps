// clang-format off
// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

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
   Contributing author:  Axel Kohlmeyer (Temple U)
   Currently maintained by:  Giacomo Fiorin (NIH)
------------------------------------------------------------------------- */

#include "fix_colvars.h"
#include "inthash.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "respa.h"
#include "universe.h"
#include "update.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_lammps.h"
#include "colvars_memstream.h"
#include "colvarscript.h"

/* struct for packed data communication of coordinates and forces. */
struct LAMMPS_NS::commdata {
  int tag, type;
  double x, y, z, m, q;
};

/***************************************************************/

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace IntHash_NS;

// initialize static class members
int FixColvars::instances = 0;

/***************************************************************
 create class and parse arguments in LAMMPS script. Syntax:

 fix ID group-ID colvars <config_file> [optional flags...]

 optional keyword value pairs:

  input   <input prefix>    (for restarting/continuing, defaults to
                             nullptr, but set to <output prefix> at end)
  output  <output prefix>   (defaults to 'out')
  seed    <integer>         (seed for RNG, defaults to '1966')
  tstat   <fix label>       (label of thermostatting fix)

 ***************************************************************/

FixColvars::FixColvars(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix colvars command: too few arguments");

  if (instances > 0)
    error->all(FLERR,"Only one colvars fix can be active at a time");
  ++instances;

  scalar_flag = 1;
  global_freq = 1;
  nevery = 1;
  extscalar = 1;
  restart_global = 1;
  energy_global_flag = 1;

  root2root = MPI_COMM_NULL;
  proxy = nullptr;

  if (strcmp(arg[3], "none") == 0) {
    conf_file = nullptr;
  } else {
    conf_file = utils::strdup(arg[3]);
  }

  rng_seed = 1966;
  unwrap_flag = 1;

  inp_name = nullptr;
  out_name = nullptr;
  tfix_name = nullptr;

  /* initialize various state variables. */
  energy = 0.0;
  nlevels_respa = 0;
  init_flag = 0;
  num_coords = 0;
  comm_buf = nullptr;
  taglist = nullptr;
  force_buf = nullptr;
  idmap = nullptr;

  script_args[0] = reinterpret_cast<unsigned char *>(utils::strdup("fix_modify"));

  parse_fix_arguments(narg, arg, true);

  if (!out_name) out_name = utils::strdup("out");

  if (comm->me == 0) {
#ifdef LAMMPS_BIGBIG
    utils::logmesg(lmp, "colvars: Warning: cannot handle atom ids > 2147483647\n");
#endif
    proxy = new colvarproxy_lammps(lmp);
    proxy->init();
    proxy->set_random_seed(rng_seed);
    proxy->set_target_temperature(t_target);
    if (conf_file) {
      proxy->add_config("configfile", conf_file);
    }
  }

  /* storage required to communicate a single coordinate or force. */
  size_one = sizeof(struct commdata);
}


int FixColvars::parse_fix_arguments(int narg, char **arg, bool fix_constructor)
{
  int const iarg_start = fix_constructor ? 4 : 0;
  int iarg = iarg_start;
  while (iarg < narg) {

    bool is_fix_keyword = false;

    if (0 == strcmp(arg[iarg], "input")) {
      inp_name = utils::strdup(arg[iarg+1]);
      // input prefix is set in FixColvars::setup()
      is_fix_keyword = true;
    } else if (0 == strcmp(arg[iarg], "output")) {
      out_name = utils::strdup(arg[iarg+1]);
      // output prefix is set in FixColvars::setup()
      is_fix_keyword = true;
    } else if (0 == strcmp(arg[iarg], "seed")) {
      rng_seed = utils::inumeric(FLERR, arg[iarg+1], false, lmp);
      is_fix_keyword = true;
    } else if (0 == strcmp(arg[iarg], "unwrap")) {
      unwrap_flag = utils::logical(FLERR, arg[iarg+1], false, lmp);
      is_fix_keyword = true;
    } else if (0 == strcmp(arg[iarg], "tstat")) {
      tfix_name = utils::strdup(arg[iarg+1]);
      if (comm->me == 0) set_thermostat_temperature();
      is_fix_keyword = true;
    }

    if (is_fix_keyword) {

      // Valid LAMMPS fix keyword: raise error if it has no argument
      if (iarg + 1 == narg) {
        if (fix_constructor) {
          error->all(FLERR, ("Missing argument to keyword \""+
                             std::string(arg[iarg]) +"\""));
        } else {
          // Error code consistent with Fix::modify_param()
          return 0;
        }
      }

    } else {

      if (fix_constructor) {
        error->all(FLERR, "Unrecognized fix colvars argument: please note that "
                   "Colvars script commands are not allowed until after the "
                   "fix is created");
      } else {
        if (iarg > iarg_start) {
          error->all(FLERR,
                     "Unrecognized fix colvars argument: please note that "
                     "you cannot combine fix colvars keywords and Colvars "
                     "script commands in the same line");
        } else {
          // Return negative error code to try the Colvars script commands
          return -1;
        }
      }
    }

    iarg += 2;
  }

  return iarg;
}


FixColvars::~FixColvars()
{
  delete[] conf_file;
  delete[] inp_name;
  delete[] out_name;
  delete[] tfix_name;
  delete[] script_args[0];

  memory->sfree(comm_buf);

  if (proxy) delete proxy;

  if (idmap) {
    inthash_destroy(idmap);
    delete idmap;
  }

  if (root2root != MPI_COMM_NULL)
    MPI_Comm_free(&root2root);

  --instances;
}

/* ---------------------------------------------------------------------- */

int FixColvars::setmask()
{
  int mask = 0;
  mask |= MIN_POST_FORCE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  mask |= POST_RUN;
  return mask;
}


void FixColvars::init()
{
  const auto me = comm->me;
  if (atom->tag_enable == 0)
    error->all(FLERR, "Cannot use fix colvars without atom IDs");

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix colvars requires an atom map, see atom_modify");

  if ((me == 0) && (update->whichflag == 2))
    error->warning(FLERR, "Using fix colvars with minimization");

  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if (init_flag) return;
  init_flag = 1;

  if (universe->nworlds > 1) {
    // create inter root communicator
    int color = 1;
    if (me == 0) {
      color = 0;
    }
    MPI_Comm_split(universe->uworld, color, universe->iworld, &root2root);
    if (me == 0) {
      proxy->set_replicas_communicator(root2root);
    }
  }
}


void FixColvars::set_thermostat_temperature()
{
  if (comm->me == 0) {
    if (tfix_name) {
      if (strcmp(tfix_name, "NULL") != 0) {
        Fix *tstat_fix = modify->get_fix_by_id(tfix_name);
        if (!tstat_fix) {
          error->one(FLERR, "Could not find thermostat fix ID {}", tfix_name);
        }
        int tmp = 0;
        double *tt = reinterpret_cast<double *>(tstat_fix->extract("t_target", tmp));
        if (tt) {
          t_target = *tt;
        } else {
          error->one(FLERR, "Fix ID {} is not a thermostat fix", tfix_name);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixColvars::init_taglist()
{
  int new_taglist_size = -1;
  const auto me = comm->me;

  if (me == 0) {

    // Number of atoms requested by Colvars
    num_coords = static_cast<int>(proxy->modify_atom_positions()->size());

    if (proxy->modified_atom_list()) {
      new_taglist_size = num_coords;
      proxy->reset_modified_atom_list();
    } else {
      new_taglist_size = -1;
    }
  }

  // Broadcast number of colvar atoms; negative means no updates
  MPI_Bcast(&new_taglist_size, 1, MPI_INT, 0, world);

  if (new_taglist_size < 0) {
    return;
  }

  num_coords = new_taglist_size;

  if (taglist) {
    memory->destroy(taglist);
    memory->destroy(force_buf);
  }
  memory->create(taglist, num_coords, "colvars:taglist");
  memory->create(force_buf, 3*num_coords, "colvars:force_buf");

  if (me == 0) {

    // Initialize and build hashtable on MPI rank 0

    std::vector<int> const &tl = *(proxy->get_atom_ids());

    if (idmap) {
      delete idmap;
      idmap = nullptr;
    }

    idmap = new inthash_t;
    inthash_init(idmap, num_coords);
    for (int i = 0; i < num_coords; ++i) {
      taglist[i] = tl[i];
      inthash_insert(idmap, tl[i], i);
    }
  }

  // Broadcast colvar atom ID list
  MPI_Bcast(taglist, num_coords, MPI_LMP_TAGINT, 0, world);
}


int FixColvars::modify_param(int narg, char **arg)
{
  if (narg >= 100) {
    error->one(FLERR, "Too many arguments for fix_modify command");
    return 2;
  }

  // Parse arguments to fix colvars
  int return_code = parse_fix_arguments(narg, arg, false);

  if (return_code >= 0) {
    // A fix colvars argument was detected, return directly
    return return_code;
  }

  // Any unknown arguments will go through the Colvars scripting interface
  if (comm->me == 0) {
    int error_code = COLVARSCRIPT_OK;
    colvarscript *script = proxy->script;
    script->set_cmdline_main_cmd("fix_modify " + std::string(id));

    for (int i = 0; i < narg; i++) {

      // Substitute LAMMPS variables

      char *new_arg = arg[i];
      int ncopy = strlen(new_arg) + 1;
      int nwork = ncopy;
      auto *copy = (char *) memory->smalloc(ncopy * sizeof(char), "fix/colvar:copy");
      auto *work = (char *) memory->smalloc(ncopy * sizeof(char), "fix/colvar:work");
      strncpy(copy, new_arg, ncopy);
      lmp->input->substitute(copy,work,ncopy,nwork,0);
      memory->sfree(work);
      new_arg = copy;

      script_args[i+1] = reinterpret_cast<unsigned char *>(new_arg);
    }

    // Run the command through Colvars
    error_code |= script->run(narg+1, script_args);

    std::string const result = proxy->get_error_msgs() + script->str_result();
    if (result.size()) utils::logmesg(lmp, result);

    // free allocated memory
    for (int i = 0; i < narg; i++) memory->sfree(script_args[i+1]);
    return (error_code == COLVARSCRIPT_OK) ? narg : 0;

  } else {

    // Return without error, don't block Fix::modify_params()
    return narg;
  }

  return 0;
}


void FixColvars::setup_io()
{
  if (comm->me == 0) {
    proxy->set_input_prefix(std::string(inp_name ? inp_name : ""));
    if (proxy->input_prefix().size() > 0) {
      proxy->log("Will read input state from file \""+
                 proxy->input_prefix()+".colvars.state\"");
    }

    proxy->set_output_prefix(std::string(out_name ? out_name : ""));

    // Try to extract a restart prefix from a potential restart command
    LAMMPS_NS::Output *outp = lmp->output;
    if ((outp->restart_every_single > 0) &&
        (outp->restart1 != nullptr)) {

      proxy->set_default_restart_frequency(outp->restart_every_single);
      proxy->set_restart_output_prefix(std::string(outp->restart1));

    } else if ((outp->restart_every_double > 0) &&
               (outp->restart2a != nullptr)) {

      proxy->set_default_restart_frequency(outp->restart_every_double);
      proxy->set_restart_output_prefix(std::string(outp->restart2a));
    }
  }
}


void FixColvars::setup(int vflag)
{
  const tagint * const tag  = atom->tag;
  const int * const type = atom->type;
  int i,nme,tmp,ndata;
  const auto nlocal = atom->nlocal;
  const auto me = comm->me;

  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    setup_io();
    proxy->parse_module_config();
  }

  init_taglist();

  // determine size of comm buffer
  nme=0;
  for (i=0; i < num_coords; ++i) {
    const tagint k = atom->map(taglist[i]);
    if ((k >= 0) && (k < nlocal))
      ++nme;
  }

  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  memory->create(comm_buf,nmax,"colvars:comm_buf");

  const double * const * const x = atom->x;
  const imageint * const image = atom->image;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;

  if (me == 0) {

    std::vector<int>     const &id = *(proxy->get_atom_ids());
    std::vector<int>           &tp = *(proxy->modify_atom_types());
    std::vector<cvm::atom_pos> &cd = *(proxy->modify_atom_positions());
    std::vector<cvm::rvector>  &of = *(proxy->modify_atom_total_forces());
    std::vector<cvm::real>     &m  = *(proxy->modify_atom_masses());
    std::vector<cvm::real>     &q  = *(proxy->modify_atom_charges());

    // store coordinate data in holding array, clear old forces


    for (i=0; i<num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal)) {

        of[i].x = of[i].y = of[i].z = 0.0;

        if (unwrap_flag) {
          const int ix = (image[k] & IMGMASK) - IMGMAX;
          const int iy = (image[k] >> IMGBITS & IMGMASK) - IMGMAX;
          const int iz = (image[k] >> IMG2BITS) - IMGMAX;

          cd[i].x = x[k][0] + ix * xprd + iy * xy + iz * xz;
          cd[i].y = x[k][1] + iy * yprd + iz * yz;
          cd[i].z = x[k][2] + iz * zprd;
        } else {
          cd[i].x = x[k][0];
          cd[i].y = x[k][1];
          cd[i].z = x[k][2];
        }
        if (atom->rmass_flag) {
          m[i] = atom->rmass[k];
        } else {
          m[i] = atom->mass[type[k]];
        }
        if (atom->q_flag) {
          q[i] = atom->q[k];
        }
      }
    }

    // loop over procs to receive and apply remote data

    for (i=1; i < comm->nprocs; ++i) {
      int maxbuf = nmax*size_one;
      MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (int k=0; k<ndata; ++k) {

        const int j = inthash_lookup(idmap, comm_buf[k].tag);

        if (j != HASH_FAIL) {

          tp[j] = comm_buf[k].type;

          cd[j].x = comm_buf[k].x;
          cd[j].y = comm_buf[k].y;
          cd[j].z = comm_buf[k].z;

          m[j] = comm_buf[k].m;
          q[j] = comm_buf[k].q;

          of[j].x = of[j].y = of[j].z = 0.0;
        }
      }
    }
  } else { // me != 0

    // copy coordinate data into communication buffer

    nme = 0;
    for (i=0; i<num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal)) {

        comm_buf[nme].tag  = tag[k];
        comm_buf[nme].type = type[k];

        if (unwrap_flag) {
          const int ix = (image[k] & IMGMASK) - IMGMAX;
          const int iy = (image[k] >> IMGBITS & IMGMASK) - IMGMAX;
          const int iz = (image[k] >> IMG2BITS) - IMGMAX;

          comm_buf[nme].x = x[k][0] + ix * xprd + iy * xy + iz * xz;
          comm_buf[nme].y = x[k][1] + iy * yprd + iz * yz;
          comm_buf[nme].z = x[k][2] + iz * zprd;
        } else {
          comm_buf[nme].x = x[k][0];
          comm_buf[nme].y = x[k][1];
          comm_buf[nme].z = x[k][2];
        }

        if (atom->rmass_flag) {
          comm_buf[nme].m = atom->rmass[k];
        } else {
          comm_buf[nme].m = atom->mass[type[k]];
        }

        if (atom->q_flag) {
          comm_buf[nme].q = atom->q[k];
        }

        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  // run pre-run setup in colvarproxy
  if (me == 0)
    proxy->setup();

  // initialize forces
  if (utils::strmatch(update->integrate_style,"^verlet") || (update->whichflag == 2))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */
/* Main colvars handler:
 * Send coodinates and add colvar forces to atoms. */
void FixColvars::post_force(int /*vflag*/)
{
  const auto me = comm->me;

  // some housekeeping: update status of the proxy as needed.
  if (me == 0) {
    if (proxy->want_exit())
      error->one(FLERR,"Run aborted on request from colvars module.\n");
  }

  const tagint * const tag = atom->tag;
  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const imageint * const image = atom->image;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;
  const int nlocal = atom->nlocal;

  /* check and potentially grow local communication buffers. */
  int i,nmax_new,nme=0;
  for (i=0; i < num_coords; ++i) {
    const tagint k = atom->map(taglist[i]);
    if ((k >= 0) && (k < nlocal))
      ++nme;
  }

  MPI_Allreduce(&nme,&nmax_new,1,MPI_INT,MPI_MAX,world);
  if (nmax_new > nmax) {
    nmax = nmax_new;
    memory->grow(comm_buf,nmax,"colvars:comm_buf");
  }

  MPI_Status status;
  MPI_Request request;
  int tmp, ndata;

  if (me == 0) {

    std::vector<cvm::atom_pos> &cd = *(proxy->modify_atom_positions());

    // store coordinate data

    for (i=0; i<num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal)) {

        if (unwrap_flag) {
          const int ix = (image[k] & IMGMASK) - IMGMAX;
          const int iy = (image[k] >> IMGBITS & IMGMASK) - IMGMAX;
          const int iz = (image[k] >> IMG2BITS) - IMGMAX;

          cd[i].x = x[k][0] + ix * xprd + iy * xy + iz * xz;
          cd[i].y = x[k][1] + iy * yprd + iz * yz;
          cd[i].z = x[k][2] + iz * zprd;
        } else {
          cd[i].x = x[k][0];
          cd[i].y = x[k][1];
          cd[i].z = x[k][2];
        }
      }
    }

    /* loop over procs to receive remote data */
    for (i=1; i < comm->nprocs; ++i) {
      int maxbuf = nmax*size_one;
      MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (int k=0; k<ndata; ++k) {
        const int j = inthash_lookup(idmap, comm_buf[k].tag);
        if (j != HASH_FAIL) {
          cd[j].x = comm_buf[k].x;
          cd[j].y = comm_buf[k].y;
          cd[j].z = comm_buf[k].z;
        }
      }
    }

  } else { // me != 0
    /* copy coordinate data into communication buffer */
    nme = 0;
    for (i=0; i<num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal)) {
        comm_buf[nme].tag = tag[k];

        if (unwrap_flag) {
          const int ix = (image[k] & IMGMASK) - IMGMAX;
          const int iy = (image[k] >> IMGBITS & IMGMASK) - IMGMAX;
          const int iz = (image[k] >> IMG2BITS) - IMGMAX;

          comm_buf[nme].x = x[k][0] + ix * xprd + iy * xy + iz * xz;
          comm_buf[nme].y = x[k][1] + iy * yprd + iz * yz;
          comm_buf[nme].z = x[k][2] + iz * zprd;
        } else {
          comm_buf[nme].x = x[k][0];
          comm_buf[nme].y = x[k][1];
          comm_buf[nme].z = x[k][2];
        }

        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  ////////////////////////////////////////////////////////////////////////
  // call our workhorse and retrieve additional information.
  if (me == 0) {
    energy = proxy->compute();
    store_forces = proxy->total_forces_enabled();
  }
  ////////////////////////////////////////////////////////////////////////

  // broadcast store_forces flag and energy data to all processors
  MPI_Bcast(&energy, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&store_forces, 1, MPI_INT, 0, world);

  // broadcast and apply biasing forces

  if (me == 0) {

    std::vector<cvm::rvector> &fo = *(proxy->modify_atom_applied_forces());

    double *fbuf = force_buf;
    for (int j=0; j < num_coords; ++j) {
      *fbuf++ = fo[j].x;
      *fbuf++ = fo[j].y;
      *fbuf++ = fo[j].z;
    }
  }
  MPI_Bcast(force_buf, 3*num_coords, MPI_DOUBLE, 0, world);

  for (int i=0; i < num_coords; ++i) {
    const tagint k = atom->map(taglist[i]);
    if ((k >= 0) && (k < nlocal)) {
      f[k][0] += force_buf[3*i+0];
      f[k][1] += force_buf[3*i+1];
      f[k][2] += force_buf[3*i+2];
    }
  }
}

/* ---------------------------------------------------------------------- */
void FixColvars::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */
void FixColvars::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  /* only process colvar forces on the outmost RESPA level. */
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */
void FixColvars::end_of_step()
{
  if (store_forces) {

    const tagint * const tag = atom->tag;
    double * const * const f = atom->f;
    const int nlocal = atom->nlocal;

    /* check and potentially grow local communication buffers. */
    int i,nmax_new,nme=0;
    for (i=0; i < num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal))
        ++nme;
    }

    MPI_Allreduce(&nme,&nmax_new,1,MPI_INT,MPI_MAX,world);
    if (nmax_new > nmax) {
      nmax = nmax_new;
      memory->grow(comm_buf,nmax,"colvars:comm_buf");
    }

    MPI_Status status;
    MPI_Request request;
    int tmp, ndata;

    if (comm->me == 0) {

      // store old force data
      std::vector<cvm::rvector> &of = *(proxy->modify_atom_total_forces());

      for (i=0; i<num_coords; ++i) {
        const tagint k = atom->map(taglist[i]);
        if ((k >= 0) && (k < nlocal)) {

          const int j = inthash_lookup(idmap, tag[k]);
          if (j != HASH_FAIL) {
            of[j].x = f[k][0];
            of[j].y = f[k][1];
            of[j].z = f[k][2];
          }
        }
      }

      /* loop over procs to receive remote data */
      for (i=1; i < comm->nprocs; ++i) {
        int maxbuf = nmax*size_one;
        MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_BYTE, &ndata);
        ndata /= size_one;

        for (int k=0; k<ndata; ++k) {
          const int j = inthash_lookup(idmap, comm_buf[k].tag);
          if (j != HASH_FAIL) {
            of[j].x = comm_buf[k].x;
            of[j].y = comm_buf[k].y;
            of[j].z = comm_buf[k].z;
          }
        }
      }

    } else { // me != 0
      /* copy total force data into communication buffer */
      nme = 0;
      for (i=0; i<num_coords; ++i) {
        const tagint k = atom->map(taglist[i]);
        if ((k >= 0) && (k < nlocal)) {
          comm_buf[nme].tag  = tag[k];
          comm_buf[nme].x    = f[k][0];
          comm_buf[nme].y    = f[k][1];
          comm_buf[nme].z    = f[k][2];
          ++nme;
        }
      }
      /* blocking receive to wait until it is our turn to send data. */
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixColvars::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    cvm::memory_stream ms;
    if (proxy->colvars->write_state(ms)) {
      int len_cv_state = ms.length();
      // Will write the buffer's length twice, so that the fix can read it later, too
      int len = len_cv_state + sizeof(int);
      fwrite(&len, sizeof(int), 1, fp);
      fwrite(&len, sizeof(int), 1, fp);
      fwrite(ms.output_buffer(), 1, len_cv_state, fp);
    } else {
      error->all(FLERR, "Failed to write Colvars state to binary file");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixColvars::restart(char *buf)
{
  if (comm->me == 0) {
    // Read the buffer's length, then load it into Colvars starting right past that location
    int length = *(reinterpret_cast<int *>(buf));
    unsigned char *colvars_state_buffer = reinterpret_cast<unsigned char *>(buf + sizeof(int));
    if (proxy->colvars->set_input_state_buffer(length, colvars_state_buffer) != COLVARS_OK) {
      error->all(FLERR, "Failed to set the Colvars input state from string buffer");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixColvars::post_run()
{
  if (comm->me == 0) {
    proxy->post_run();
    if (lmp->citeme) {
      lmp->citeme->add(proxy->colvars->feature_report(1));
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixColvars::compute_scalar()
{
  return energy;
}

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixColvars::memory_usage()
{
  double bytes = (double) (num_coords * (2*sizeof(int)+3*sizeof(double)));
  bytes += (double)(double) (nmax*size_one) + sizeof(this);
  return bytes;
}
