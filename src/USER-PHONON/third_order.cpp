//
// Created by charlie sievers on 7/5/18.
//


#include <mpi.h>
#include <cstdlib>
#include "third_order.h"
#include "atom.h"
#include "complex"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "force.h"
#include "math_extra.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "update.h"
#include "neighbor.h"
#include "pair.h"
#include "timer.h"
#include "finish.h"


using namespace LAMMPS_NS;
enum{REGULAR,BALLISTICO};

/* ---------------------------------------------------------------------- */

ThirdOrder::ThirdOrder(LAMMPS *lmp) : Pointers(lmp), fp(NULL)
{
    external_force_clear = 1;
}

/* ---------------------------------------------------------------------- */

ThirdOrder::~ThirdOrder()
{
    if (fp) fclose(fp);
    fp = NULL;
}

/* ----------------------------------------------------------------------
   setup without output or one-time post-init setup
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void ThirdOrder::setup()
{
    // setup domain, communication and neighboring
    // acquire ghosts
    // build neighbor lists
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    domain->image_check();
    domain->box_too_small_check();
    neighbor->build(1);
    neighbor->ncalls = 0;
    neighbor->every = 3;                       // build every this many steps
    neighbor->delay = 1;
    neighbor->ago = 0;
    neighbor->ndanger = 0;

    // compute all forces
    force_clear();
    external_force_clear = 0;

    eflag=0;
    vflag=0;
    if (pair_compute_flag) force->pair->compute(eflag,vflag);
    else if (force->pair) force->pair->compute_dummy(eflag,vflag);

    if (atom->molecular) {
        if (force->bond) force->bond->compute(eflag,vflag);
        if (force->angle) force->angle->compute(eflag,vflag);
        if (force->dihedral) force->dihedral->compute(eflag,vflag);
        if (force->improper) force->improper->compute(eflag,vflag);
    }

    if (force->kspace) {
        force->kspace->setup();
        if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
        else force->kspace->compute_dummy(eflag,vflag);
    }

    if (force->newton) comm->reverse_comm();
}

/* ---------------------------------------------------------------------- */

void ThirdOrder::command(int narg, char **arg)
{
    MPI_Comm_rank(world,&me);

    if (domain->box_exist == 0)
        error->all(FLERR,"Dynamical_matrix command before simulation box is defined");
    if (narg < 2) error->all(FLERR,"Illegal dynamical_matrix command");

    lmp->init();

    // orthogonal vs triclinic simulation box

    triclinic = domain->triclinic;

    if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
    else pair_compute_flag = 0;
    if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
    else kspace_compute_flag = 0;

    // group and style

    igroup = group->find(arg[0]);
    if (igroup == -1) error->all(FLERR,"Could not find dynamical matrix group ID");
    groupbit = group->bitmask[igroup];
    update->setupflag = 1;

    int style = -1;
    if (strcmp(arg[1],"regular") == 0) style = REGULAR;
    else if (strcmp(arg[1],"ballistico") == 0) style = BALLISTICO;
    else error->all(FLERR,"Illegal Dynamical Matrix command");

    // set option defaults

    binaryflag = 0;
    scaleflag = 0;
    compressed = 0;
    file_flag = 0;
    file_opened = 0;
    conversion = 1;

    // read options from end of input line
    if (style == REGULAR) options(narg-3,&arg[3]);  //COME BACK
    else if (style == BALLISTICO) options(narg-3,&arg[3]); //COME BACK
    else if (comm->me == 0 && screen) fprintf(screen,"Illegal Dynamical Matrix command\n");
    del = force->numeric(FLERR, arg[2]);

    // move atoms by 3-vector or specified variable(s)

    if (style == REGULAR) {
        setup();
        calculateMatrix();
    }

    if (style == BALLISTICO) {
        setup();
        convert_units(update->unit_style);
        conversion = conv_energy/conv_distance/conv_distance;
        calculateMatrix();
    }

    Finish finish(lmp);
    finish.end(1);
}

/* ----------------------------------------------------------------------
   parse optional parameters
------------------------------------------------------------------------- */

void ThirdOrder::options(int narg, char **arg)
{
    if (narg < 0) error->all(FLERR,"Illegal dynamical_matrix command");
    int iarg = 0;
    const char *filename = "third_order.txt";
    std::stringstream fss;

    while (iarg < narg) {
        if (strcmp(arg[iarg],"file") == 0) {
            if (iarg+2 > narg) error->all(FLERR, "Illegal dynamical_matrix command");
            fss << arg[iarg + 1] << me;
            filename = fss.str().c_str();
            file_flag = 1;
            iarg += 2;
        }
        else if (strcmp(arg[iarg],"binary") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal dynamical_matrix command");
            if (strcmp(arg[iarg+1],"gzip") == 0) {
                compressed = 1;
            }
            else if (strcmp(arg[iarg+1],"yes") == 0) {
                binaryflag = 1;
            }
            iarg += 2;
        } else error->all(FLERR,"Illegal dynamical_matrix command");
    }
    if (file_flag == 1) {
        openfile(filename);
    }
}

/* ----------------------------------------------------------------------
   generic opening of a file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void ThirdOrder::openfile(const char* filename)
{
    // if file already opened, return
    if (file_opened) return;

    if (compressed) {
#ifdef LAMMPS_GZIP
        char gzip[128];
    sprintf(gzip,"gzip -6 > %s",filename);
#ifdef _WIN32
    fp = _popen(gzip,"wb");
#else
    fp = popen(gzip,"w");
#endif
#else
        error->one(FLERR,"Cannot open gzipped file");
#endif
    } else if (binaryflag) {
        fp = fopen(filename,"wb");
    } else {
        fp = fopen(filename,"w");
    }

    if (fp == NULL) error->one(FLERR,"Cannot open dump file");

    file_opened = 1;
}

/* ----------------------------------------------------------------------
   create dynamical matrix
------------------------------------------------------------------------- */

void ThirdOrder::calculateMatrix()
{
    int local_idx; // local index
    int local_jdx; // second local index
    int local_kdx; // third local index
    int natoms = atom->natoms;
    int *mask = atom->mask;
    double **f = atom->f;

    energy_force(0);

    if (comm->me == 0 && screen) fprintf(screen,"Calculating Anharmonic Dynamical Matrix...\n");

    for (int i=1; i<=natoms; i++){
        local_idx = atom->map(i);
        if (local_idx >= 0 && mask[local_idx] && groupbit){
            for (int alpha=0; alpha<3; alpha++){
                displace_atom(local_idx, alpha, 1);
                for (int j=1; j<=natoms; j++){
                    local_jdx = atom->map(j);
                    if (local_jdx >= 0&& mask[local_jdx] && groupbit){
                        for (int beta=0; beta<3; beta++){
                        }
                    }
                }
                displace_atom(local_idx,alpha,-2);
                for (int j=1; j<=natoms; j++){
                    local_jdx = atom->map(j);
                    if (local_jdx >= 0 && mask[local_jdx] && groupbit){
                        for (int beta=0; beta<3; beta++){

                        }
                    }
                }
                displace_atom(local_idx,alpha,1);
            }
        }
    }

    //for (int proc1=0; proc1 < comm->nprocs; proc1++) {
    //    plocal1 = atom->nlocal;  // 1 proc nlocal = 8
    //    MPI_Bcast(&plocal1, 1, MPI_INT, proc1, MPI_COMM_WORLD); // plocal1 = 8
    //    for (int i = 0; i < plocal1; i++) {
    //        if (me==proc1 & mask[i] & groupbit)
	//	group_flag_1 = 1;
    //        MPI_Bcast(&group_flag_1, 1, MPI_INT, proc1, MPI_COMM_WORLD);
    //        if (group_flag_1) {
    //            if (me == proc1) id1 = aid[i];
    //            MPI_Bcast(&id1, 1, MPI_INT, proc1, MPI_COMM_WORLD);
    //            for (int alpha = 0; alpha < 3; alpha++) {
    //                for (int proc2 = 0; proc2 < comm->nprocs; proc2++) {
    //                    plocal2 = atom->nlocal;
    //                    MPI_Bcast(&plocal2, 1, MPI_INT, proc2, MPI_COMM_WORLD);
    //                    for (int j = 0; j < plocal2; j++) {
    //                        if (me==proc2 & mask[j] & groupbit)
	//			group_flag_2 = 1;
    //                        MPI_Bcast(&group_flag_2, 1, MPI_INT, proc2, MPI_COMM_WORLD);
    //                        if (mask[j] & groupbit) {
    //                            if (me == proc2) id2 = aid[j];
    //                            MPI_Bcast(&id2, 1, MPI_INT, proc2, MPI_COMM_WORLD);
    //                            for (int beta = 0; beta < 3; beta++) {
//
    //                                if (me == proc1) x[i][alpha] += del;
//
    //                                if (me == proc2) x[j][beta] += del;
    //                                energy_force(0);
////
    //                                for (int gamma = 0; gamma < 3; gamma++) {
    //                                    for (int k = 0; k < nlocal; k++)
    //                                        if (mask[k] & groupbit) {
    //                                            first_derv[k*3+gamma] = f[k][gamma];
    //                                        }
    //                                }
//
    //                                if (me == proc2) x[j][beta] -= 2 * del;
    //                                energy_force(0);
////
    //                                for (int gamma = 0; gamma < 3; gamma++) {
    //                                    for (int k = 0; k < nlocal; k++)
    //                                        if (mask[k] & groupbit) {
    //                                            first_derv[k*3+gamma] -= f[k][gamma];
    //                                        }
    //                                }
//
    //                                if (me == proc2) x[j][beta] += 2 * del;
//
    //                                if (me == proc1) x[i][alpha] -= 2 * del;
//
    //                                energy_force(0);
////
    //                                for (int gamma = 0; gamma < 3; gamma++) {
    //                                    for (int k = 0; k < nlocal; k++)
    //                                        if (mask[k] & groupbit) {
    //                                            first_derv[k*3+gamma] -= f[k][gamma];
    //                                        }
    //                                }
////
    //                                if (me == proc2) x[j][beta] -= 2 * del;
    //                                energy_force(0);
////
    //                                for (int k = 0; k < nlocal; k++)
    //                                    if (mask[k] & groupbit) {
    //                                        for (int gamma = 0; gamma < 3; gamma++) {
    //                                            first_derv[k*3+gamma] += f[k][gamma];
    //                                            first_derv[k*3+gamma] /= -4*del*del;
    //                                        }
    //                                        double norm = pow(first_derv[k*3], 2)
    //                                                      + pow(first_derv[k*3+1], 2)
    //                                                      + pow(first_derv[k+3+2], 2);
    //                                        if (fp && norm > 1.0e-16)
    //                                            fprintf(fp,
    //                                                    "%d %d %d %d %d %7.8f %7.8f %7.8f\n",
    //                                                    id1, alpha + 1, id2, beta + 1, aid[k],
    //                                                    first_derv[k*3] * conversion,
    //                                                    first_derv[k*3+1] * conversion,
    //                                                    first_derv[k*3+2] * conversion);
    //                                    }
////
    //                                if (me == proc2) x[j][beta] += del;
//
    //                                if (me == proc1) x[i][alpha] += del;
    //                            }
    //                        }
    //                    }
    //                }
////
    //            }
    //        }
    //    }
    //}
//

    if (screen && me ==0 ) fprintf(screen,"Finished Calculating Third Order Tensor\n");

}

/* ----------------------------------------------------------------------
  Displace atoms
   ---------------------------------------------------------------------- */

void ThirdOrder::displace_atom(int local_idx, int direction, int magnitude)
{
    double **x = atom->x;
    int *sametag = atom->sametag;
    int j = local_idx;

    x[local_idx][direction] += del*magnitude;

    while (sametag[j] >= 0){
        j = sametag[j];
        x[j][direction] += del*magnitude;
    }
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms due to reneighboring
   return new energy, which should include nextra_global dof
   return negative gradient stored in atom->f
   return negative gradient for nextra_global dof in fextra
------------------------------------------------------------------------- */

void ThirdOrder::energy_force(int resetflag)
{
    // check for reneighboring
    // always communicate since atoms move
    int nflag = neighbor->decide();

    if (nflag == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
    } else {
        if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
            domain->reset_box();
            comm->setup();
            if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();

        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(Timer::COMM);

        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
    }

    force_clear();

    timer->stamp();

    if (pair_compute_flag) {
        force->pair->compute(eflag,vflag);
        timer->stamp(Timer::PAIR);
    }

    if (atom->molecular) {
        if (force->bond) force->bond->compute(eflag,vflag);
        if (force->angle) force->angle->compute(eflag,vflag);
        if (force->dihedral) force->dihedral->compute(eflag,vflag);
        if (force->improper) force->improper->compute(eflag,vflag);
        timer->stamp(Timer::BOND);
    }

    if (kspace_compute_flag) {
        force->kspace->compute(eflag,vflag);
        timer->stamp(Timer::KSPACE);
    }

    if (force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
    }
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void ThirdOrder::force_clear()
{
    if (external_force_clear) return;

    // clear global force array
    // if either newton flag is set, also include ghosts

    size_t nbytes = sizeof(double) * atom->nlocal;
    if (force->newton) nbytes += sizeof(double) * atom->nghost;

    if (nbytes) {
        memset(&atom->f[0][0],0,3*nbytes);
    }
}

/* ---------------------------------------------------------------------- */

void ThirdOrder::convert_units(const char *style)
{
    // physical constants from:
    // http://physics.nist.gov/cuu/Constants/Table/allascii.txt
    // using thermochemical calorie = 4.184 J

    if (strcmp(style,"lj") == 0) {
        error->all(FLERR,"Conversion Not Set");
        //conversion = 1; // lj -> 10 J/mol

    } else if (strcmp(style,"real") == 0) {
        conv_energy = 418.4; // kcal/mol -> 10 J/mol
        conv_mass = 1; // g/mol -> g/mol
        conv_distance = 1; // angstrom -> angstrom

    } else if (strcmp(style,"metal") == 0) {
        conv_energy = 9648.5; // eV -> 10 J/mol
        conv_mass = 1; // g/mol -> g/mol
        conv_distance = 1; // angstrom -> angstrom

    } else if (strcmp(style,"si") == 0) {
        if (comm->me) error->warning(FLERR,"Conversion Warning: Multiplication by Large Float");
        conv_energy = 6.022E22; // J -> 10 J/mol
        conv_mass = 6.022E26; // kg -> g/mol
        conv_distance = 1E-10; // meter -> angstrom

    } else if (strcmp(style,"cgs") == 0) {
        if (comm->me) error->warning(FLERR,"Conversion Warning: Multiplication by Large Float");
        conv_energy = 6.022E12; // Erg -> 10 J/mol
        conv_mass = 6.022E23; // g -> g/mol
        conv_distance = 1E-7; // centimeter -> angstrom

    } else if (strcmp(style,"electron") == 0) {
        conv_energy = 262550; // Hartree -> 10 J/mol
        conv_mass = 1; // amu -> g/mol
        conv_distance = 0.529177249; // bohr -> angstrom

    } else if (strcmp(style,"micro") == 0) {
        if (comm->me) error->warning(FLERR,"Conversion Warning: Untested Conversion");
        conv_energy = 6.022E10; // picogram-micrometer^2/microsecond^2 -> 10 J/mol
        conv_mass = 6.022E11; // pg -> g/mol
        conv_distance = 1E-4; // micrometer -> angstrom

    } else if (strcmp(style,"nano") == 0) {
        if (comm->me) error->warning(FLERR,"Conversion Warning: Untested Conversion");
        conv_energy = 6.022E4; // attogram-nanometer^2/nanosecond^2 -> 10 J/mol
        conv_mass = 6.022E5; // ag -> g/mol
        conv_distance = 0.1; // angstrom -> angstrom

    } else error->all(FLERR,"Units Type Conversion Not Found");

}
