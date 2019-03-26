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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)

   Please cite the related publication:
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_neb_spin.h"
#include "universe.h"
#include "update.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{SINGLE_PROC_DIRECT,SINGLE_PROC_MAP,MULTI_PROC};

#define BUFSIZE 8

/* ---------------------------------------------------------------------- */

FixNEB_spin::FixNEB_spin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_pe(NULL), pe(NULL), nlenall(NULL), xprev(NULL), 
  xnext(NULL), fnext(NULL), spprev(NULL), spnext(NULL), fmnext(NULL), springF(NULL), 
  tangent(NULL), xsend(NULL), xrecv(NULL), fsend(NULL), frecv(NULL), spsend(NULL), 
  sprecv(NULL), fmsend(NULL), fmrecv(NULL), tagsend(NULL), tagrecv(NULL), 
  xsendall(NULL), xrecvall(NULL), fsendall(NULL), frecvall(NULL), spsendall(NULL), 
  sprecvall(NULL), fmsendall(NULL), fmrecvall(NULL), tagsendall(NULL), tagrecvall(NULL), 
  counts(NULL), displacements(NULL)
{

  if (narg < 4) error->all(FLERR,"Illegal fix neb_spin command");

  kspring = force->numeric(FLERR,arg[3]);
  if (kspring <= 0.0) error->all(FLERR,"Illegal fix neb command");

  // optional params

  NEBLongRange = false; // see if needed (comb. with pppm/spin?)
  StandardNEB = true;   // only option for now
  PerpSpring = FreeEndIni = FreeEndFinal = false; 
  FreeEndFinalWithRespToEIni = FinalAndInterWithRespToEIni = false;
  kspringPerp = 0.0;
  kspringIni = 1.0;
  kspringFinal = 1.0;
  SpinLattice = false;	// no spin-lattice neb for now

  // no available fix neb/spin options for now
  
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"parallel") == 0) {
      error->all(FLERR,"Illegal fix neb command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"perp") == 0) {
      error->all(FLERR,"Illegal fix neb command");
      iarg += 2;
    } else if (strcmp (arg[iarg],"end") == 0) {
      iarg += 3;
    } else if (strcmp (arg[iarg],"lattice") == 0) {
      iarg += 2;
    } else error->all(FLERR,"Illegal fix neb command");
  }
  
  // nreplica = number of partitions
  // ireplica = which world I am in universe
  // nprocs_universe = # of procs in all replicase
  // procprev,procnext = root proc in adjacent replicas

  me = comm->me;
  nprocs = comm->nprocs;

  nprocs_universe = universe->nprocs;
  nreplica = universe->nworlds;
  ireplica = universe->iworld;

  if (ireplica > 0) procprev = universe->root_proc[ireplica-1];
  else procprev = -1;
  if (ireplica < nreplica-1) procnext = universe->root_proc[ireplica+1];
  else procnext = -1;

  uworld = universe->uworld;
  int *iroots = new int[nreplica];
  MPI_Group uworldgroup,rootgroup;
  if (NEBLongRange) {
    for (int i=0; i<nreplica; i++)
      iroots[i] = universe->root_proc[i];
    MPI_Comm_group(uworld, &uworldgroup);
    MPI_Group_incl(uworldgroup, nreplica, iroots, &rootgroup);
    MPI_Comm_create(uworld, rootgroup, &rootworld);
  }
  delete [] iroots;

  // create a new compute pe style
  // id = fix-ID + pe, compute group = all

  int n = strlen(id) + 4;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  char **newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  delete [] newarg;

  // initialize local storage

  maxlocal = -1;
  ntotal = -1;
}

/* ---------------------------------------------------------------------- */

FixNEB_spin::~FixNEB_spin()
{
  modify->delete_compute(id_pe);
  delete [] id_pe;

  // memory destroy of all spin and lattice arrays

  memory->destroy(xprev);
  memory->destroy(xnext);
  memory->destroy(tangent);
  memory->destroy(fnext);
  memory->destroy(spprev);
  memory->destroy(spnext);
  memory->destroy(fmnext);
  memory->destroy(springF);
  memory->destroy(xsend);
  memory->destroy(xrecv);
  memory->destroy(fsend);
  memory->destroy(frecv);
  memory->destroy(spsend);
  memory->destroy(sprecv);
  memory->destroy(fmsend);
  memory->destroy(fmrecv);
  memory->destroy(tagsend);
  memory->destroy(tagrecv);

  memory->destroy(xsendall);
  memory->destroy(xrecvall);
  memory->destroy(fsendall);
  memory->destroy(frecvall);
  memory->destroy(spsendall);
  memory->destroy(sprecvall);
  memory->destroy(fmsendall);
  memory->destroy(fmrecvall);
  memory->destroy(tagsendall);
  memory->destroy(tagrecvall);

  memory->destroy(counts);
  memory->destroy(displacements);

  if (NEBLongRange) {
    if (rootworld != MPI_COMM_NULL) MPI_Comm_free(&rootworld);
    memory->destroy(nlenall);
  }
}

/* ---------------------------------------------------------------------- */

int FixNEB_spin::setmask()
{
  int mask = 0;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNEB_spin::init()
{
  int icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"Potential energy ID for fix neb does not exist");
  pe = modify->compute[icompute];

  // turn off climbing mode, NEB command turns it on after init()

  rclimber = -1;

  // nebatoms = # of atoms in fix group = atoms with inter-replica forces

  bigint count = group->count(igroup);
  if (count > MAXSMALLINT) error->all(FLERR,"Too many active NEB atoms");
  nebatoms = count;

  // comm mode for inter-replica exchange of coords

  if (nreplica == nprocs_universe &&
      nebatoms == atom->natoms && atom->sortfreq == 0)
    cmode = SINGLE_PROC_DIRECT;
  else if (nreplica == nprocs_universe) cmode = SINGLE_PROC_MAP;
  else cmode = MULTI_PROC;

  // ntotal = total # of atoms in system, NEB atoms or not

  if (atom->natoms > MAXSMALLINT) error->all(FLERR,"Too many atoms for NEB");
  ntotal = atom->natoms;

  if (atom->nmax > maxlocal) reallocate();

  if (MULTI_PROC && counts == NULL) {
    memory->create(xsendall,ntotal,3,"neb:xsendall");
    memory->create(xrecvall,ntotal,3,"neb:xrecvall");
    memory->create(fsendall,ntotal,3,"neb:fsendall");
    memory->create(frecvall,ntotal,3,"neb:frecvall");
    memory->create(tagsendall,ntotal,"neb:tagsendall");
    memory->create(tagrecvall,ntotal,"neb:tagrecvall");
    memory->create(spsendall,ntotal,3,"neb:xsendall");
    memory->create(sprecvall,ntotal,3,"neb:xrecvall");
    memory->create(fmsendall,ntotal,3,"neb:fsendall");
    memory->create(fmrecvall,ntotal,3,"neb:frecvall");
    memory->create(counts,nprocs,"neb:counts");
    memory->create(displacements,nprocs,"neb:displacements");
  }
}

/* ---------------------------------------------------------------------- */

void FixNEB_spin::min_setup(int vflag)
{
  min_post_force(vflag);

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixNEB_spin::min_post_force(int /*vflag*/)
{
  double vprev,vnext;
  double delspxp,delspyp,delspzp;
  double delspxn,delspyn,delspzn;
  double templen;
  double vIni=0.0;
  double spi[3],spj[3];

  vprev = vnext = veng = pe->compute_scalar();

  if (ireplica < nreplica-1 && me == 0)
    MPI_Send(&veng,1,MPI_DOUBLE,procnext,0,uworld);
  if (ireplica > 0 && me == 0)
    MPI_Recv(&vprev,1,MPI_DOUBLE,procprev,0,uworld,MPI_STATUS_IGNORE);

  if (ireplica > 0 && me == 0)
    MPI_Send(&veng,1,MPI_DOUBLE,procprev,0,uworld);
  if (ireplica < nreplica-1 && me == 0)
    MPI_Recv(&vnext,1,MPI_DOUBLE,procnext,0,uworld,MPI_STATUS_IGNORE);

  if (cmode == MULTI_PROC) {
    MPI_Bcast(&vprev,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&vnext,1,MPI_DOUBLE,0,world);
  }

  if (FreeEndFinal && ireplica == nreplica-1 && (update->ntimestep == 0)) 
    error->all(FLERR,"NEB_spin Free End option not yet active");

  if (ireplica == 0) vIni=veng;

  if (FreeEndFinalWithRespToEIni) 
    error->all(FLERR,"NEB_spin Free End option not yet active");

  if (FreeEndIni && ireplica == 0 && (update->ntimestep == 0))
    error->all(FLERR,"NEB_spin Free End option not yet active");


  // communicate atoms to/from adjacent replicas to fill xprev,xnext

  inter_replica_comm();

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **sp = atom->sp;
  double dot = 0.0;
  double prefactor = 0.0;
  double **fm = atom->fm;

  //calculating separation between images

  plen = 0.0;
  nlen = 0.0;
  double tlen = 0.0;
  double gradnextlen = 0.0;
  double delndots, delpdots;

  dotgrad = gradlen = dotpath = dottangrad = 0.0;

  // computation of the tangent vector

  if (ireplica == nreplica-1) {		// final replica
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	
	// tangent vector
	
	delspxp = sp[i][0] - spprev[i][0];
	delspyp = sp[i][1] - spprev[i][1];
	delspzp = sp[i][2] - spprev[i][2];
        
	// project delp vector on tangent space
	
	delpdots = delspxp*sp[i][0]+delspyp*sp[i][1]+delspzp*sp[i][2];
	delspxp -= delpdots*sp[i][0];
	delspyp -= delpdots*sp[i][1];
	delspzp -= delpdots*sp[i][2];
        
	// calc. geodesic length
	
	spi[0] = sp[i][0];
	spi[1] = sp[i][1];
	spi[2] = sp[i][2];
	spj[0] = spprev[i][0];
	spj[1] = spprev[i][1];
	spj[2] = spprev[i][2];
	templen = geodesic_distance(spi,spj);
	plen += templen*templen;
	dottangrad += delspxp*fm[i][0]+ delspyp*fm[i][1]+delspzp*fm[i][2];
        gradlen += fm[i][0]*fm[i][0] + fm[i][1]*fm[i][1] + fm[i][2]*fm[i][2];
	
	// no free end option for now 
        
	if (FreeEndFinal||FreeEndFinalWithRespToEIni) 
	  error->all(FLERR,"Free End option not yet active");
      
      }
  } else if (ireplica == 0) {		// initial replica
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

	// tangent vector
	
	delspxn = spnext[i][0]- sp[i][0];
	delspyn = spnext[i][1]- sp[i][1];
	delspzn = spnext[i][2]- sp[i][2];

        // project deln vector on tangent space
	
	delndots = delspxn*sp[i][0]+delspyn*sp[i][1]+delspzn*sp[i][2];
	delspxn -= delndots*sp[i][0];
	delspyn -= delndots*sp[i][1];
	delspzn -= delndots*sp[i][2];
	
	// calc. geodesic length
	
	spi[0]=sp[i][0];
	spi[1]=sp[i][1];
	spi[2]=sp[i][2];
	spj[0]=spnext[i][0];
	spj[1]=spnext[i][1];
	spj[2]=spnext[i][2];
	templen = geodesic_distance(spi,spj);
	nlen += templen*templen;
	dottangrad += delspxn*fm[i][0] + delspyn*fm[i][1] + delspzn*fm[i][2];
        gradlen += fm[i][0]*fm[i][0] + fm[i][1]*fm[i][1] + fm[i][2]*fm[i][2];
	
	// no free end option for now 
        
        if (FreeEndIni) 
	  error->all(FLERR,"Free End option not yet active");

      }
  } else {			// intermediate replica

    double vmax = MAX(fabs(vnext-veng),fabs(vprev-veng));
    double vmin = MIN(fabs(vnext-veng),fabs(vprev-veng));

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        
	// calc. delp vector

	delspxp = sp[i][0] - spprev[i][0];
        delspyp = sp[i][1] - spprev[i][1];
        delspzp = sp[i][2] - spprev[i][2];
        
	// project delp vector on tangent space

	delndots = delspxp*sp[i][0]+delspyp*sp[i][1]+delspzp*sp[i][2];
	delspxp -= delpdots*sp[i][0];
	delspyp -= delpdots*sp[i][1];
	delspzp -= delpdots*sp[i][2];

	// calc. prev. geodesic length

	spi[0]=sp[i][0];
	spi[1]=sp[i][1];
	spi[2]=sp[i][2];
	spj[0]=spprev[i][0];
	spj[1]=spprev[i][1];
	spj[2]=spprev[i][2];
	templen = geodesic_distance(spi, spj);
	plen += templen*templen;

	// calc. deln vector

	delspxn = spnext[i][0] - sp[i][0];
        delspyn = spnext[i][1] - sp[i][1];
        delspzn = spnext[i][2] - sp[i][2];
        
	// project deln vector on tangent space

	delndots = delspxn*sp[i][0]+delspyn*sp[i][1]+delspzn*sp[i][2];
	delspxn -= delndots*sp[i][0];
	delspyn -= delndots*sp[i][1];
	delspzn -= delndots*sp[i][2];

	// evaluate best path tangent

        if (vnext > veng && veng > vprev) {
          tangent[i][0] = delspxn;
          tangent[i][1] = delspyn;
          tangent[i][2] = delspzn;
        } else if (vnext < veng && veng < vprev) {
          tangent[i][0] = delspxp;
          tangent[i][1] = delspyp;
          tangent[i][2] = delspzp;
        } else {
          if (vnext > vprev) {
            tangent[i][0] = vmax*delspxn + vmin*delspxp;
            tangent[i][1] = vmax*delspyn + vmin*delspyp;
            tangent[i][2] = vmax*delspzn + vmin*delspzp;
          } else if (vnext < vprev) {
            tangent[i][0] = vmin*delspxn + vmax*delspxp;
            tangent[i][1] = vmin*delspyn + vmax*delspyp;
            tangent[i][2] = vmin*delspzn + vmax*delspzp;
          } else { // vnext == vprev, e.g. for potentials that do not compute an energy
            tangent[i][0] = delspxn + delspxp;
            tangent[i][1] = delspyn + delspyp;
            tangent[i][2] = delspzn + delspzp;
          }
        }

	// calc. next geodesic length
	
	spi[0]=sp[i][0];
	spi[1]=sp[i][1];
	spi[2]=sp[i][2];
	spj[0]=spnext[i][0];
	spj[1]=spnext[i][1];
	spj[2]=spnext[i][2];
	templen = geodesic_distance(spi, spj);
	nlen += templen*templen;
    
        tlen += tangent[i][0]*tangent[i][0] + tangent[i][1]*tangent[i][1] + 
	  tangent[i][2]*tangent[i][2];
        gradlen += fm[i][0]*fm[i][0] + fm[i][1]*fm[i][1] + fm[i][2]*fm[i][2];
        dotpath += delspxp*delspxn + delspyp*delspyn + delspzp*delspzn;
        dottangrad += tangent[i][0]*fm[i][0] + tangent[i][1]*fm[i][1] + 
	  tangent[i][2]*fm[i][2];
        gradnextlen += fnext[i][0]*fnext[i][0] + fnext[i][1]*fnext[i][1] +
	  fnext[i][2]*fnext[i][2];
        dotgrad += fm[i][0]*fnext[i][0] + fm[i][1]*fnext[i][1] +
          fm[i][2]*fnext[i][2];

        // no Perpendicular nudging force option active yet
        
	if (kspringPerp != 0.0) 
	  error->all(FLERR,"NEB_spin Perpendicular nudging force not yet active");

      }
  }

  // MPI reduce if more than one proc for world

  double bufin[BUFSIZE], bufout[BUFSIZE];
  bufin[0] = nlen;
  bufin[1] = plen;
  bufin[2] = tlen;
  bufin[3] = gradlen;
  bufin[4] = gradnextlen;
  bufin[5] = dotpath;
  bufin[6] = dottangrad;
  bufin[7] = dotgrad;
  MPI_Allreduce(bufin,bufout,BUFSIZE,MPI_DOUBLE,MPI_SUM,world);
  nlen = sqrt(bufout[0]);
  plen = sqrt(bufout[1]);
  tlen = sqrt(bufout[2]);
  gradlen = sqrt(bufout[3]);
  gradnextlen = sqrt(bufout[4]);
  dotpath = bufout[5];
  dottangrad = bufout[6];
  dotgrad = bufout[7];

  // project tangent vector on tangent space and normalize it

  double buftan[3];
  double tandots;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      tandots = tangent[i][0]*sp[i][0]+tangent[i][1]*sp[i][1]+
	tangent[i][2]*sp[i][2];
      buftan[0] = tangent[i][0]-tandots*sp[i][0];
      buftan[1] = tangent[i][1]-tandots*sp[i][1];
      buftan[2] = tangent[i][2]-tandots*sp[i][2];
      tangent[i][0] = buftan[0];
      tangent[i][1] = buftan[1];
      tangent[i][2] = buftan[2];

      if (tlen > 0.0) {
	double tleninv = 1.0/tlen;
        tangent[i][0] *= tleninv;
        tangent[i][1] *= tleninv;
        tangent[i][2] *= tleninv;
      }
    }

  // first or last replica has no change to forces, just return

  if (ireplica > 0 && ireplica < nreplica-1)
    dottangrad = dottangrad/(tlen*gradlen);
  if (ireplica == 0)
    dottangrad = dottangrad/(nlen*gradlen);
  if (ireplica == nreplica-1)
    dottangrad = dottangrad/(plen*gradlen);
  if (ireplica < nreplica-1)
    dotgrad = dotgrad /(gradlen*gradnextlen);

  // no Free End options active yet
  
  if (FreeEndIni && ireplica == 0) 
    error->all(FLERR,"NEB_spin Free End option not yet active");
  if (FreeEndFinal && ireplica == nreplica -1) 
    error->all(FLERR,"NEB_spin Free End option not yet active");
  if (FreeEndFinalWithRespToEIni&&ireplica == nreplica -1) 
    error->all(FLERR,"NEB_spin Free End option not yet active");

  // no NEB_spin long range option 
  
  if (NEBLongRange) 
    error->all(FLERR,"NEB_spin long range option not yet active");

  // test output length
  
  //printf("testi irep / plen: %d %g \n",ireplica,nlen);

  // exit calc. if first or last replica (no gneb force)

  if (ireplica == 0 || ireplica == nreplica-1) return ;

  double AngularContr;
  dotpath = dotpath/(plen*nlen);
  AngularContr = 0.5 *(1+cos(MY_PI * dotpath));


  for (int i = 0; i < nlocal; i++) 
    if (mask[i] & groupbit) {
      dot += fm[i][0]*tangent[i][0] + fm[i][1]*tangent[i][1] +
        fm[i][2]*tangent[i][2];
    }

  // gather all dot for this replica
  
  double dotall;
  MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);
  dot=dotall;

  // for intermediate replica
  // calc. GNEB force prefactor

  if (ireplica == rclimber) prefactor = -2.0*dot;	// for climbing replica
  else {
    if (NEBLongRange) {	
      error->all(FLERR,"Long Range NEB_spin climber option not yet active");
    } else if (StandardNEB) {
      prefactor = -dot + kspring*(nlen-plen);
    }

    if (FinalAndInterWithRespToEIni && veng<vIni) {
      for (int i = 0; i < nlocal; i++)
       if (mask[i] & groupbit) {
          fm[i][0] = 0;
          fm[i][1] = 0;
          fm[i][2] = 0;
        }
      prefactor =  kspring*(nlen-plen);
      AngularContr=0;
    }
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fm[i][0] += prefactor*tangent[i][0];
      fm[i][1] += prefactor*tangent[i][1];
      fm[i][2] += prefactor*tangent[i][2];
    }
  
  // project GNEB force on tangent space

  double fdots;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      fdots = fm[i][0]*sp[i][0] + fm[i][1]*sp[i][1] +
	fm[i][2]*sp[i][2];
      fm[i][0] -= fdots*sp[i][0];
      fm[i][1] -= fdots*sp[i][1];
      fm[i][2] -= fdots*sp[i][2];
    }

}

/* ----------------------------------------------------------------------
   geodesic distance calculation (Vincenty's formula)
------------------------------------------------------------------------- */

double FixNEB_spin::geodesic_distance(double spi[3], double spj[3])
{
  double dist;
  double crossx,crossy,crossz;
  double dotx,doty,dotz;
  double normcross,dots;

  crossx = spi[1]*spj[2]-spi[2]*spj[1];
  crossy = spi[2]*spj[0]-spi[0]*spj[2];
  crossz = spi[0]*spj[1]-spi[1]*spj[0];
  normcross = sqrt(crossx*crossx + crossy*crossy + crossz*crossz);
  
  dotx = spi[0]*spj[0];
  doty = spi[1]*spj[1];
  dotz = spi[2]*spj[2];
  dots = dotx+doty+dotz;

  if (normcross == 0.0 && dots == 0.0) 
    error->all(FLERR,"Incorrect calc. of geodesic_distance in Fix NEB/spin");
  
    dist = atan2(normcross,dots);

  return dist;
}

/* ----------------------------------------------------------------------
   send/recv NEB atoms to/from adjacent replicas
   received atoms matching my local atoms are stored in xprev,xnext
   replicas 0 and N-1 send but do not receive any atoms
------------------------------------------------------------------------- */

void FixNEB_spin::inter_replica_comm()
{
  int i,m;
  MPI_Request request;
  MPI_Request requests[2];
  MPI_Status statuses[2];

  // reallocate memory if necessary

  if (atom->nmax > maxlocal) reallocate();

  double **x = atom->x;
  double **f = atom->f;
  double **sp = atom->sp;
  double **fm = atom->fm;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // -----------------------------------------------------
  // 3 cases: two for single proc per replica
  //          one for multiple procs per replica
  // -----------------------------------------------------


  // case 1 => to be done

  // single proc per replica
  // all atoms are NEB atoms and no atom sorting
  // direct comm of x -> xprev and x -> xnext

  if (cmode == SINGLE_PROC_DIRECT) {
    if (ireplica > 0)
      MPI_Irecv(xprev[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld,&request);
      MPI_Irecv(spprev[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld,&request);
    if (ireplica < nreplica-1)
      MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld);
      MPI_Send(sp[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld);
    if (ireplica > 0) MPI_Wait(&request,MPI_STATUS_IGNORE);
    if (ireplica < nreplica-1)
      MPI_Irecv(xnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
      MPI_Irecv(spnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
    if (ireplica > 0)
      MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(sp[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
    if (ireplica < nreplica-1) MPI_Wait(&request,MPI_STATUS_IGNORE);

    if (ireplica < nreplica-1)
      MPI_Irecv(fnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
      MPI_Irecv(fmnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
    if (ireplica > 0)
      MPI_Send(f[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(fm[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
    if (ireplica < nreplica-1) MPI_Wait(&request,MPI_STATUS_IGNORE);

    return;
  }

  // single proc per replica
  // but only some atoms are NEB atoms or atom sorting is enabled
  // send atom IDs and coords of only NEB atoms to prev/next proc
  // recv procs use atom->map() to match received coords to owned atoms

  if (cmode == SINGLE_PROC_MAP) {
    m = 0;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        tagsend[m] = tag[i];
        xsend[m][0] = x[i][0];
        xsend[m][1] = x[i][1];
        xsend[m][2] = x[i][2];
        fsend[m][0] = f[i][0];
        fsend[m][1] = f[i][1];
        fsend[m][2] = f[i][2];
        spsend[m][0] = sp[i][0];
        spsend[m][1] = sp[i][1];
        spsend[m][2] = sp[i][2];
        fmsend[m][0] = fm[i][0];
        fmsend[m][1] = fm[i][1];
        fmsend[m][2] = fm[i][2];
        m++;
      }

    if (ireplica > 0) {
      MPI_Irecv(xrecv[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld,&requests[0]);
      MPI_Irecv(sprecv[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld,&requests[0]);
      MPI_Irecv(tagrecv,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld,&requests[1]);
    }
    if (ireplica < nreplica-1) {
      MPI_Send(xsend[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld);
      MPI_Send(spsend[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld);
      MPI_Send(tagsend,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld);
    }

    if (ireplica > 0) {
      MPI_Waitall(2,requests,statuses);
      for (i = 0; i < nebatoms; i++) {
        m = atom->map(tagrecv[i]);
        xprev[m][0] = xrecv[i][0];
        xprev[m][1] = xrecv[i][1];
        xprev[m][2] = xrecv[i][2];
        spprev[m][0] = sprecv[i][0];
        spprev[m][1] = sprecv[i][1];
        spprev[m][2] = sprecv[i][2];
      }
    }
    if (ireplica < nreplica-1) {
      MPI_Irecv(xrecv[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
      MPI_Irecv(frecv[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
      MPI_Irecv(sprecv[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
      MPI_Irecv(fmrecv[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
      MPI_Irecv(tagrecv,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld,&requests[1]);
    }
    if (ireplica > 0) {
      MPI_Send(xsend[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(fsend[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(spsend[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(fmsend[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(tagsend,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld);
    }

    if (ireplica < nreplica-1) {
      MPI_Waitall(2,requests,statuses);
      for (i = 0; i < nebatoms; i++) {
        m = atom->map(tagrecv[i]);
        xnext[m][0] = xrecv[i][0];
        xnext[m][1] = xrecv[i][1];
        xnext[m][2] = xrecv[i][2];
        fnext[m][0] = frecv[i][0];
        fnext[m][1] = frecv[i][1];
        fnext[m][2] = frecv[i][2];
        spnext[m][0] = sprecv[i][0];
        spnext[m][1] = sprecv[i][1];
        spnext[m][2] = sprecv[i][2];
        fmnext[m][0] = fmrecv[i][0];
        fmnext[m][1] = fmrecv[i][1];
        fmnext[m][2] = fmrecv[i][2];
      }
    }

    return;
  }

  // multiple procs per replica
  // MPI_Gather all coords and atom IDs to root proc of each replica
  // send to root of adjacent replicas
  // bcast within each replica
  // each proc extracts info for atoms it owns via atom->map()

  m = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      tagsend[m] = tag[i];
      xsend[m][0] = x[i][0];
      xsend[m][1] = x[i][1];
      xsend[m][2] = x[i][2];
      fsend[m][0] = f[i][0];
      fsend[m][1] = f[i][1];
      fsend[m][2] = f[i][2];
      spsend[m][0] = sp[i][0];
      spsend[m][1] = sp[i][1];
      spsend[m][2] = sp[i][2];
      fmsend[m][0] = fm[i][0];
      fmsend[m][1] = fm[i][1];
      fmsend[m][2] = fm[i][2];
      m++;
    }

  MPI_Gather(&m,1,MPI_INT,counts,1,MPI_INT,0,world);
  displacements[0] = 0;
  for (i = 0; i < nprocs-1; i++)
    displacements[i+1] = displacements[i] + counts[i];
  MPI_Gatherv(tagsend,m,MPI_LMP_TAGINT,
              tagsendall,counts,displacements,MPI_LMP_TAGINT,0,world);
  for (i = 0; i < nprocs; i++) counts[i] *= 3;
  for (i = 0; i < nprocs-1; i++)
    displacements[i+1] = displacements[i] + counts[i];
  if (xsend) {
    MPI_Gatherv(xsend[0],3*m,MPI_DOUBLE,
                xsendall[0],counts,displacements,MPI_DOUBLE,0,world);
    MPI_Gatherv(fsend[0],3*m,MPI_DOUBLE,
                fsendall[0],counts,displacements,MPI_DOUBLE,0,world);
  } else {
    MPI_Gatherv(NULL,3*m,MPI_DOUBLE,
                xsendall[0],counts,displacements,MPI_DOUBLE,0,world);
    MPI_Gatherv(NULL,3*m,MPI_DOUBLE,
                fsendall[0],counts,displacements,MPI_DOUBLE,0,world);
  }
  if (spsend) {
    MPI_Gatherv(spsend[0],3*m,MPI_DOUBLE,
                spsendall[0],counts,displacements,MPI_DOUBLE,0,world);
    MPI_Gatherv(fmsend[0],3*m,MPI_DOUBLE,
                fmsendall[0],counts,displacements,MPI_DOUBLE,0,world);
  } else {
    MPI_Gatherv(NULL,3*m,MPI_DOUBLE,
                spsendall[0],counts,displacements,MPI_DOUBLE,0,world);
    MPI_Gatherv(NULL,3*m,MPI_DOUBLE,
                fmsendall[0],counts,displacements,MPI_DOUBLE,0,world);
  }

  if (ireplica > 0 && me == 0) {
    MPI_Irecv(xrecvall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld,&requests[0]);
    MPI_Irecv(sprecvall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld,&requests[0]);
    MPI_Irecv(tagrecvall,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld,
              &requests[1]);
  }
  if (ireplica < nreplica-1 && me == 0) {
    MPI_Send(xsendall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld);
    MPI_Send(spsendall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld);
    MPI_Send(tagsendall,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld);
  }

  if (ireplica > 0) {
    if (me == 0) MPI_Waitall(2,requests,statuses);

    MPI_Bcast(tagrecvall,nebatoms,MPI_INT,0,world);
    MPI_Bcast(xrecvall[0],3*nebatoms,MPI_DOUBLE,0,world);
    MPI_Bcast(sprecvall[0],3*nebatoms,MPI_DOUBLE,0,world);

    for (i = 0; i < nebatoms; i++) {
      m = atom->map(tagrecvall[i]);
      if (m < 0 || m >= nlocal) continue;
      xprev[m][0] = xrecvall[i][0];
      xprev[m][1] = xrecvall[i][1];
      xprev[m][2] = xrecvall[i][2];
      spprev[m][0] = sprecvall[i][0];
      spprev[m][1] = sprecvall[i][1];
      spprev[m][2] = sprecvall[i][2];
    }
  }

  if (ireplica < nreplica-1 && me == 0) {
    MPI_Irecv(xrecvall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
    MPI_Irecv(frecvall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
    MPI_Irecv(sprecvall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
    MPI_Irecv(sprecvall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
    MPI_Irecv(tagrecvall,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld,
              &requests[1]);
  }
  if (ireplica > 0 && me == 0) {
    MPI_Send(xsendall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
    MPI_Send(fsendall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
    MPI_Send(spsendall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
    MPI_Send(fmsendall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
    MPI_Send(tagsendall,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld);
  }

  if (ireplica < nreplica-1) {
    if (me == 0) MPI_Waitall(2,requests,statuses);

    MPI_Bcast(tagrecvall,nebatoms,MPI_INT,0,world);
    MPI_Bcast(xrecvall[0],3*nebatoms,MPI_DOUBLE,0,world);
    MPI_Bcast(frecvall[0],3*nebatoms,MPI_DOUBLE,0,world);
    MPI_Bcast(sprecvall[0],3*nebatoms,MPI_DOUBLE,0,world);
    MPI_Bcast(fmrecvall[0],3*nebatoms,MPI_DOUBLE,0,world);

    for (i = 0; i < nebatoms; i++) {
      m = atom->map(tagrecvall[i]);
      if (m < 0 || m >= nlocal) continue;
      xnext[m][0] = xrecvall[i][0];
      xnext[m][1] = xrecvall[i][1];
      xnext[m][2] = xrecvall[i][2];
      fnext[m][0] = frecvall[i][0];
      fnext[m][1] = frecvall[i][1];
      fnext[m][2] = frecvall[i][2];
      spnext[m][0] = sprecvall[i][0];
      spnext[m][1] = sprecvall[i][1];
      spnext[m][2] = sprecvall[i][2];
      fmnext[m][0] = fmrecvall[i][0];
      fmnext[m][1] = fmrecvall[i][1];
      fmnext[m][2] = fmrecvall[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate xprev,xnext,tangent arrays if necessary
   reallocate communication arrays if necessary
------------------------------------------------------------------------- */

void FixNEB_spin::reallocate()
{
  maxlocal = atom->nmax;

  memory->destroy(xprev);
  memory->destroy(xnext);
  memory->destroy(tangent);
  memory->destroy(fnext);
  memory->destroy(springF);
  memory->destroy(spprev);
  memory->destroy(spnext);
  memory->destroy(fmnext);

  memory->create(xprev,maxlocal,3,"neb:xprev");
  memory->create(xnext,maxlocal,3,"neb:xnext");
  memory->create(tangent,maxlocal,3,"neb:tangent");
  memory->create(fnext,maxlocal,3,"neb:fnext");
  memory->create(springF,maxlocal,3,"neb:springF");
  memory->create(spprev,maxlocal,3,"neb:xprev");
  memory->create(spnext,maxlocal,3,"neb:xnext");
  memory->create(fmnext,maxlocal,3,"neb:fnext");

  if (cmode != SINGLE_PROC_DIRECT) {
    memory->destroy(xsend);
    memory->destroy(fsend);
    memory->destroy(xrecv);
    memory->destroy(frecv);
    memory->destroy(spsend);
    memory->destroy(fmsend);
    memory->destroy(sprecv);
    memory->destroy(fmrecv);
    memory->destroy(tagsend);
    memory->destroy(tagrecv);
    memory->create(xsend,maxlocal,3,"neb:xsend");
    memory->create(fsend,maxlocal,3,"neb:fsend");
    memory->create(xrecv,maxlocal,3,"neb:xrecv");
    memory->create(frecv,maxlocal,3,"neb:frecv");
    memory->create(spsend,maxlocal,3,"neb:xsend");
    memory->create(fmsend,maxlocal,3,"neb:fsend");
    memory->create(sprecv,maxlocal,3,"neb:xrecv");
    memory->create(fmrecv,maxlocal,3,"neb:frecv");
    memory->create(tagsend,maxlocal,"neb:tagsend");
    memory->create(tagrecv,maxlocal,"neb:tagrecv");
  }

  if (NEBLongRange) {
    memory->destroy(nlenall);
    memory->create(nlenall,nreplica,"neb:nlenall");
  }
}
