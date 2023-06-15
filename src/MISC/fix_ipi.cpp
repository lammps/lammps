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
   Contributing author: Michele Ceriotti (EPFL), Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_ipi.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "irregular.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/******************************************************************************************
 * A fix to interface LAMMPS with i-PI - A Python interface for path integral molecular dynamics
 * Michele Ceriotti, EPFL (2014)
 * Please cite:
 * Ceriotti, M., More, J., & Manolopoulos, D. E. (2014).
 * i-PI: A Python interface for ab initio path integral molecular dynamics simulations.
 * Computer Physics Communications, 185, 1019–1026. doi:10.1016/j.cpc.2013.10.027
 * And see [https://github.com/i-pi/i-pi] to download a version of i-PI
 ******************************************************************************************/

// socket interface
#ifndef _WIN32
#include <netdb.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/un.h>
#include <unistd.h>
#else
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <io.h>
#include <windows.h>
#endif

#define MSGLEN 12

/* Utility functions to simplify the interface with POSIX sockets */

static void open_socket(int &sockfd, int inet, int port, char *host, Error *error)
/* Opens a socket.

   Args:
   sockfd: The id of the socket that will be created.
   inet: An integer that determines whether the socket will be an inet or unix
   domain socket. Gives unix if 0, inet otherwise.
   port: The port number for the socket to be created. Low numbers are often
   reserved for important channels, so use of numbers of 4 or more digits is
   recommended.
   host: The name of the host server.
   error: pointer to a LAMMPS Error object
*/
{
  int ai_err;

#ifdef _WIN32
  error->one(FLERR, "i-PI socket implementation requires UNIX environment");
#else
  if (inet > 0) {    // creates an internet socket

    // fetches information on the host
    struct addrinfo hints, *res;

    memset(&hints, 0, sizeof(hints));
    hints.ai_socktype = SOCK_STREAM;
    hints.ai_family = AF_UNSPEC;
    hints.ai_flags = AI_PASSIVE;

    ai_err = getaddrinfo(host, std::to_string(port).c_str(), &hints, &res);
    if (ai_err != 0) error->one(FLERR, "Error fetching host data. Wrong host name?");

    // creates socket
    sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
    if (sockfd < 0) error->one(FLERR, "Error opening socket");

    // makes connection
    if (connect(sockfd, res->ai_addr, res->ai_addrlen) < 0)
      error->one(FLERR, "Error opening INET socket: wrong port or server unreachable");
    freeaddrinfo(res);

  } else {    // creates a unix socket
    struct sockaddr_un serv_addr;

    // fills up details of the socket address
    memset(&serv_addr, 0, sizeof(serv_addr));
    serv_addr.sun_family = AF_UNIX;
    strcpy(serv_addr.sun_path, "/tmp/ipi_");
    strcpy(serv_addr.sun_path + 9, host);

    // creates the socket
    sockfd = socket(AF_UNIX, SOCK_STREAM, 0);

    // connects
    if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0)
      error->one(FLERR,
                 "Error opening UNIX socket: server may not be running "
                 "or the path to the socket unavailable");
  }
#endif
}

static void writebuffer(int sockfd, const char *data, int len, Error *error)
/* Writes to a socket.

   Args:
   sockfd: The id of the socket that will be written to.
   data: The data to be written to the socket.
   len: The length of the data in bytes.
*/
{
  int n;

  n = write(sockfd, data, len);
  if (n < 0) error->one(FLERR, "Error writing to socket: broken connection");
}

static void readbuffer(int sockfd, char *data, int len, Error *error)
/* Reads from a socket.

   Args:
   sockfd: The id of the socket that will be read from.
   data: The storage array for data read from the socket.
   len: The length of the data in bytes.
*/
{
  int n, nr;

  n = nr = read(sockfd, data, len);

  while (nr > 0 && n < len) {
    nr = read(sockfd, &data[n], len - n);
    n += nr;
  }

  if (n == 0) error->one(FLERR, "Error reading from socket: broken connection");
}

/* ---------------------------------------------------------------------- */

FixIPI::FixIPI(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), irregular(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix ipi", error);

  if (atom->tag_enable == 0) error->all(FLERR, "Cannot use fix ipi without atom IDs");
  if (atom->tag_consecutive() == 0) error->all(FLERR, "Fix ipi requires consecutive atom IDs");
  if (strcmp(update->unit_style, "lj") == 0) error->all(FLERR, "Fix ipi does not support lj units");

  if ((strcmp(arg[1], "all") != 0) && (comm->me == 0))
    error->warning(FLERR, "Not using group 'all' with fix ipi can result in undefined behavior");

  host = strdup(arg[3]);
  port = utils::inumeric(FLERR, arg[4], false, lmp);

  master = (comm->me == 0) ? 1 : 0;
  inet = 1;
  reset_flag = 0;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "unix") == 0) {
      inet = 0;
      ++iarg;
    } else if (strcmp(arg[iarg], "reset") == 0) {
      reset_flag = 1;
      ++iarg;
    } else {
      error->all(FLERR, "Unknown fix ipi keyword: {}", arg[iarg]);
    }
  }

  // sanity check
  if (inet && ((port <= 1024) || (port > 65536)))
    error->all(FLERR, "Invalid port for fix ipi: {}", port);

  hasdata = bsize = 0;

  // creates a temperature compute for all atoms
  modify->add_compute("IPI_TEMP all temp");

  // creates a  pressure compute to extract the virial
  modify->add_compute("IPI_PRESS all pressure IPI_TEMP virial");

  // create instance of Irregular class
  irregular = new Irregular(lmp);

  // yet, we have not assigned a socket
  socketflag = 0;
}

/* ---------------------------------------------------------------------- */

FixIPI::~FixIPI()
{
  if (bsize) delete[] buffer;
  free(host);
  modify->delete_compute("IPI_TEMP");
  modify->delete_compute("IPI_PRESS");
  delete irregular;
}

/* ---------------------------------------------------------------------- */

int FixIPI::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIPI::init()
{
  //only opens socket on master process
  if (master) {
    if (!socketflag) open_socket(ipisock, inet, port, host, error);
  } else
    ipisock = 0;
  // TODO: should check for success in socket opening,
  // but the current open_socket routine dies brutally if unsuccessful
  // tell lammps we have assigned a socket
  socketflag = 1;

  // asks for evaluation of PE at first step
  modify->compute[modify->find_compute("thermo_pe")]->invoked_scalar = -1;
  modify->addstep_compute_all(update->ntimestep + 1);

  kspace_flag = (force->kspace) ? 1 : 0;

  // makes sure that neighbor lists are re-built at each step
  // (cannot make assumptions when cycling over beads!)
  neighbor->delay = 0;
  neighbor->every = 1;
}

// clang-format off
void FixIPI::initial_integrate(int /*vflag*/)
{
  /* This is called at the beginning of the integration loop,
   * and will be used to read positions from the socket. Then,
   * everything should be updated, since there is no guarantee
   * that successive snapshots will be close together (think
   * of parallel tempering for instance) */

  char header[MSGLEN+1];

  if (hasdata)
    error->all(FLERR, "i-PI got out of sync in initial_integrate and will die!");

  double cellh[9], cellih[9];
  int nat;
  if (master)   { // only read positions on master

    // wait until something happens
    while (true) {
      // while i-PI just asks for status, signal we are ready and wait
      readbuffer(ipisock, header, MSGLEN, error); header[MSGLEN]=0;

      if (strcmp(header,"STATUS      ") == 0 )
        writebuffer(ipisock,"READY       ",MSGLEN, error);
      else break;
    }

    if (strcmp(header,"EXIT        ") == 0 )
      error->one(FLERR, "Got EXIT message from i-PI. Now leaving!");

    // when i-PI signals it has positions to evaluate new forces,
    // read positions and cell data
    if (strcmp(header,"POSDATA     ") == 0 )  {
      readbuffer(ipisock, (char*) cellh, 9*8, error);
      readbuffer(ipisock, (char*) cellih, 9*8, error);
      readbuffer(ipisock, (char*) &nat, 4, error);

      // allocate buffer, but only do this once.
      if (bsize==0) {
        bsize=3*nat;
        buffer = new double[bsize];
      } else if (bsize != 3*nat)
        error->one(FLERR, "Number of atoms changed along the way.");

      // finally read position data into buffer
      readbuffer(ipisock, (char*) buffer, 8*bsize, error);
    } else
      error->one(FLERR, "Wrapper did not send positions, I will now die!");
  }

  // shares the atomic coordinates with everyone
  MPI_Bcast(&nat,1,MPI_INT,0,world);
  // must also allocate the buffer on the non-head nodes
  if (bsize==0) {
    bsize=3*nat;
    buffer = new double[bsize];
  }
  MPI_Bcast(cellh,9,MPI_DOUBLE,0,world);
  MPI_Bcast(cellih,9,MPI_DOUBLE,0,world);
  MPI_Bcast(buffer,bsize,MPI_DOUBLE,0,world);

  //updates atomic coordinates and cell based on the data received
  double *boxhi = domain->boxhi;
  double *boxlo = domain->boxlo;
  double posconv;
  posconv=0.52917721*force->angstrom;
  boxlo[0] = -0.5*cellh[0]*posconv;
  boxlo[1] = -0.5*cellh[4]*posconv;
  boxlo[2] = -0.5*cellh[8]*posconv;
  boxhi[0] = -boxlo[0];
  boxhi[1] = -boxlo[1];
  boxhi[2] = -boxlo[2];
  domain->xy = cellh[1]*posconv;
  domain->xz = cellh[2]*posconv;
  domain->yz = cellh[5]*posconv;

  // do error checks on simulation box and set small for triclinic boxes
  domain->set_initial_box();
  // reset global and local box using the new box dimensions
  domain->reset_box();
  // signal that the box has (or may have) changed
  domain->box_change = 1;

  // picks local atoms from the buffer
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0]=buffer[3*(atom->tag[i]-1)+0]*posconv;
      x[i][1]=buffer[3*(atom->tag[i]-1)+1]*posconv;
      x[i][2]=buffer[3*(atom->tag[i]-1)+2]*posconv;
    }
  }

  // ensure atoms are in current box & update box via shrink-wrap
  // has to be be done before invoking Irregular::migrate_atoms()
  //   since it requires atoms be inside simulation box

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // move atoms to new processors via irregular()
  // only needed if migrate_check() says an atom moves to far
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (irregular->migrate_check()) irregular->migrate_atoms();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if kspace solver is used
  if (reset_flag && kspace_flag) {
    // reset kspace, pair, angles, ... b/c simulation box might have changed.
    //   kspace->setup() is in some cases not enough since, e.g., g_ewald needs
    //   to be reestimated due to changes in box dimensions.
    force->init();
    // reset_grid() is necessary for pppm since init() is not calling
    //   setup() nor reset_grid() upon calling init().
    if (force->kspace->pppmflag) force->kspace->reset_grid();
    // other kspace styles might need too another setup()?
  } else if (!reset_flag && kspace_flag) {
    // original version
    force->kspace->setup();
  }

  // compute PE. makes sure that it will be evaluated at next step
  modify->compute[modify->find_compute("thermo_pe")]->invoked_scalar = -1;
  modify->addstep_compute_all(update->ntimestep+1);

  hasdata=1;
}

void FixIPI::final_integrate()
{
  /* This is called after forces and energy have been computed. Now we only need to
   * communicate them back to i-PI so that the integration can continue. */
  char header[MSGLEN+1];
  double vir[9], pot=0.0;
  double forceconv, potconv, posconv, pressconv, posconv3;
  char retstr[1024];

  // conversions from LAMMPS units to atomic units, which are used by i-PI
  potconv=3.1668152e-06/force->boltz;
  posconv=0.52917721*force->angstrom;
  posconv3=posconv*posconv*posconv;
  forceconv=potconv*posconv;
  pressconv=1/force->nktv2p*potconv*posconv3;

  // compute for potential energy
  pot=modify->compute[modify->find_compute("thermo_pe")]->compute_scalar();
  pot*=potconv;

  // probably useless check
  if (!hasdata)
    error->all(FLERR, "i-PI got out of sync in final_integrate and will die!");

  int nat=bsize/3;
  double **f= atom->f;
  auto lbuf = new double[bsize];

  // reassembles the force vector from the local arrays
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  for (int i = 0; i < bsize; ++i) lbuf[i]=0.0;
  for (int i = 0; i < nlocal; i++) {
    lbuf[3*(atom->tag[i]-1)+0]=f[i][0]*forceconv;
    lbuf[3*(atom->tag[i]-1)+1]=f[i][1]*forceconv;
    lbuf[3*(atom->tag[i]-1)+2]=f[i][2]*forceconv;
  }
  MPI_Allreduce(lbuf,buffer,bsize,MPI_DOUBLE,MPI_SUM,world);
  delete[] lbuf;

  for (int i = 0; i < 9; ++i) vir[i]=0.0;

  int press_id = modify->find_compute("IPI_PRESS");
  Compute* comp_p = modify->compute[press_id];
  comp_p->compute_vector();
  double myvol = domain->xprd*domain->yprd*domain->zprd/posconv3;

  vir[0] = comp_p->vector[0]*pressconv*myvol;
  vir[4] = comp_p->vector[1]*pressconv*myvol;
  vir[8] = comp_p->vector[2]*pressconv*myvol;
  vir[1] = comp_p->vector[3]*pressconv*myvol;
  vir[2] = comp_p->vector[4]*pressconv*myvol;
  vir[5] = comp_p->vector[5]*pressconv*myvol;
  retstr[0]=0;

  if (master) {
    while (true) {
      readbuffer(ipisock, header, MSGLEN, error); header[MSGLEN]=0;

      if (strcmp(header,"STATUS      ") == 0 )
        writebuffer(ipisock,"HAVEDATA    ",MSGLEN, error);
      else break;
    }

    if (strcmp(header,"EXIT        ") == 0 )
      error->one(FLERR, "Got EXIT message from i-PI. Now leaving!");

    if (strcmp(header,"GETFORCE    ") == 0 )  {
      writebuffer(ipisock,"FORCEREADY  ",MSGLEN, error);
      writebuffer(ipisock,(char*) &pot,8, error);
      writebuffer(ipisock,(char*) &nat,4, error);
      writebuffer(ipisock,(char*) buffer, bsize*8, error);
      writebuffer(ipisock,(char*) vir,9*8, error);
      nat=strlen(retstr);  writebuffer(ipisock,(char*) &nat,4, error);
      writebuffer(ipisock,(char*) retstr, nat, error);
    }
    else
      error->one(FLERR, "Wrapper did not ask for forces, I will now die!");
  }

  hasdata=0;
}


