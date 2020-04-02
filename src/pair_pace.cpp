//
// Created by lysogy36 on 27.02.20.
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_pace.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"


#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4

//added YL
int elements_num = 104;
char const *const elements[104] = {"X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
                                   "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
                                   "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
                                   "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
                                   "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                                   "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                                   "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                                   "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
};

int AtomicNumberByName(char *elname) {
    int i = 0;
    for (i = 1; i < elements_num; i++)
        if (strcmp(elname, elements[i]) == 0)
            return i;
    return -1;
}

/* ---------------------------------------------------------------------- */
PairPACE::PairPACE(LAMMPS *lmp) : Pair(lmp) {
    //single_enable = 0;
    restartinfo = 0;
    one_coeff = 1;
    manybody_flag = 1;

    nelements = 0;

    ace = NULL;
    potential_file_name = NULL;
    elements = NULL;
    map = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairPACE::~PairPACE() {
    if (copymode) return;

    if (elements)
        for (int i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;


    delete[] potential_file_name;

    delete basis_set;
    delete ace;

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(map);
    }
}

/* ---------------------------------------------------------------------- */

void PairPACE::compute(int eflag, int vflag) {
    //printf("--> PairPACE::compute\n");
    int i, j, ii, jj, inum, jnum, k;
    double delx, dely, delz, evdwl, rsq;
    double fij[3];
    int *ilist, *jlist, *numneigh, **firstneigh;
    evdwl = 0.0;

    ev_init(eflag, vflag);

    // downwards modified by YL

    double **x = atom->x;
    double **f = atom->f;
    tagint *tag = atom->tag;
    int *type = atom->type;
    // number of atoms in cell
    int nlocal = atom->nlocal;
    //printf("nlocal=%d\n",nlocal);
    int newton_pair = force->newton_pair;
    //printf("newton_pair=%d\n",newton_pair);

    // number of atoms including ghost atoms
    int nall = nlocal + atom->nghost;
    //printf("nall=%d\n",nall);

    // inum: length of the neighborlists list
    inum = list->inum;
    //printf("inum=%d\n",inum);
    // ilist: list of "i" atoms for which neighbor lists exist
    ilist = list->ilist;
    //numneigh: the length of each these neigbor list
    numneigh = list->numneigh;
    // the pointer to the list of neighbors of "i"
    firstneigh = list->firstneigh;

    if (inum != nlocal) {
        printf("inum: %d nlocal: %d are different.\n", inum, nlocal);
        exit(0);
    }


    // Aidan Thompson told RD (26 July 2019) that practically always holds:
    // inum = nlocal
    // i = ilist(ii) < inum
    // j = jlist(jj) < nall
    // neighborlist contains neighbor atoms plus skin atoms,
    //       skin atoms can be removed by setting skin to zero but here
    //       they are disregarded anyway

    // TODO: need to get elements properly from LAMMPS
    // bad: need proper string assignment
    // - make memory allocation more efficient and avoid regular malloc/free operations
    // - allocate, populate and distribute three vectors

    // loop all atoms and neighbors

    //from ACE.compute
    ace->e_atom = 0;

    //determine the maximum numer of neighbours
    int max_jnum = -1;
    int nei = 0;
    // printf("list->inum=%d\n",list->inum);
    for (ii = 0; ii < list->inum; ii++) {
        i = ilist[ii];
        //   printf("ilist[%d]=%d\n",ii,i);
        jnum = numneigh[i];
        // printf("numneigh[%d]=%d\n",i,jnum);
        nei = nei + jnum;
        //printf("firstneigh[%d]/jlist=\n",i);
        jlist = firstneigh[i];
        // loop all neighbors of atom i
//        for (jj = 0; jj < jnum; jj++) {
//                printf("%d ", jlist[jj]);
//        }
//        printf("\n");

        if (jnum > max_jnum)
            max_jnum = jnum;
    }

    //printf("Max num of neighbours = %d\n", max_jnum);

    ace->resize_neighbours_cache(max_jnum);
    //printf("Caches resized to %d\n", max_jnum);

    //loop over atoms
    for (ii = 0; ii < list->inum; ii++) {
        i = list->ilist[ii];
        k = map[type[i]];

        // this to be removed once fully integrated in LAMMPS
//        if (k != 0) {
//            printf("Error: Cannot handle multiple species.\n");
//            printf("Stopping.\n");
//            exit(0);
//        }

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];

        const int itype = type[i];
        const int ielem = map[itype];


        //	printf("Position of at.%d=(%f,%f,%f)\n",i,x[i][0],x[i][1],x[i][2]);
        jlist = firstneigh[i];
        jnum = numneigh[i];

        ace->compute_atom(i, x, type, jnum, jlist);

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
            delx = x[j][0] - xtmp;
            dely = x[j][1] - ytmp;
            delz = x[j][2] - ztmp;

            fij[0] = ace->neighbours_forces(jj, 0);
            fij[1] = ace->neighbours_forces(jj, 1);
            fij[2] = ace->neighbours_forces(jj, 2);

            //printf("f_ij(i=%d, j=%d) = (%f, %f, %f)\n",i,j,fij[0], fij[1], fij[2]);

            f[i][0] += fij[0];
            f[i][1] += fij[1];
            f[i][2] += fij[2];
            f[j][0] -= fij[0];
            f[j][1] -= fij[1];
            f[j][2] -= fij[2];

            // tally per-atom virial contribution
            if (vflag)
                ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0,
                             fij[0], fij[1], fij[2],
                             -delx, -dely, -delz);
        }

        // tally energy contribution

        if (eflag) {

            // evdwl = energy of atom I

            evdwl = ace->e_atom;

            ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
    }


    //printf("vflag_fdotr=%d\n", vflag_fdotr);
    if (vflag_fdotr) virial_fdotr_compute();


    // end modifications YL
    //-------------------------------------------------------
}

/* ---------------------------------------------------------------------- */

void PairPACE::allocate() {
    //   printf("--> PairPACE::allocate\n");
    allocated = 1;
    int n = atom->ntypes;

    memory->create(setflag, n + 1, n + 1, "pair:setflag");
    memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
    memory->create(map, n + 1, "pair:map");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPACE::settings(int narg, char **arg) {
    printf("--> PairPACE::settings\n");
    if (narg > 0)
        error->all(FLERR,
                   "Illegal pair_style command. Correct form:\npair_style pace");

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairPACE::coeff(int narg, char **arg) {
    //TODO: FOR DEBUG
    for (int i = 0; i < narg; i++) {
        printf("PairPACE::coeff: arg[%d] = %s\n", i, arg[i]);
    }

    if (narg < 4)
        error->all(FLERR,
                   "Incorrect args for pair coefficients. Correct form:\npair_coeff * * <potential.ace> elem1 elem2 ...");

    if (!allocated) allocate();

    //number of provided elements in pair_coeff line
    int ntypes_coeff = narg - 3;

    if (ntypes_coeff != atom->ntypes) {
        char error_message[1024];
        sprintf(error_message,
                "Incorrect args for pair coefficients. You provided %d elements in pair_coeff, but structure has %d atom types",
                ntypes_coeff, atom->ntypes);
        error->all(FLERR, error_message);
    }

    char *type1 = arg[0];
    char *type2 = arg[1];
    char *potential_file_name = arg[2];
    char **elemtypes = &arg[3];

    // insure I,J args are * *

    if (strcmp(type1, "*") != 0 || strcmp(type2, "*") != 0)
        error->all(FLERR, "Incorrect args for pair coefficients");

    printf("Create basis set \n");
    basis_set = new ACECTildeBasisSet();
    printf("Loading %s\n", potential_file_name);
    basis_set->load(potential_file_name);
    printf("Loaded \n");

    // read args that map atom types to pACE elements
    // map[i] = which element the Ith atom type is, -1 if not mapped
    // map[0] is not used

    for (int i = 0; i < atom->ntypes; i++) {
        char *elemname = elemtypes[i];
        if (AtomicNumberByName(elemname) == -1) {
            printf("String '%s' is not a valid element\n", elemname);
            error->all(FLERR, "Incorrect args for pair coefficients");
        }
        map[i] = i;
    }

    // clear setflag since coeff() called once with I,J = * *

    int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) {
            setflag[i][j] = 1;
        }
    }

    // set setflag i,j for type pairs where both are mapped to elements

    //TODO: wroaround why serflag doesn't work here
    int count = 1;
    for (int i = 1; i <= n; i++)
        for (int j = i; j <= n; j++)
            if (map[i] >= 0 && map[j] >= 0) {
                setflag[i][j] = 1;
                count++;
            }

    if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

    printf("Create ACE\n");
    ace = new ACECTildeEvaluator();
    printf("ACE.set basis\n");
    ace->set_basis(*basis_set);

    //TODO: adapt
//    if (ncoeff != snaptr->ncoeff) {
//        if (comm->me == 0)
//            printf("ncoeff = %d snancoeff = %d \n",ncoeff,snaptr->ncoeff);
//        error->all(FLERR,"Incorrect SNAP parameter file");
//    }

    // Calculate maximum cutoff for all elements

//    rcutmax = 0.0;
//    for (int ielem = 0; ielem < nelements; ielem++)
//        rcutmax = MAX(2.0*radelem[ielem]*rcutfac,rcutmax);


    //TODO: debug
    printf("PairPACE::coeff done\n");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPACE::init_style() {
    if (atom->tag_enable == 0)
        error->all(FLERR, "Pair style pACE requires atom IDs");
    if (force->newton_pair == 0)
        error->all(FLERR, "Pair style pACE requires newton pair on");

    // need a full neighbor list

    int irequest = neighbor->request(this, instance_me);
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPACE::init_one(int i, int j) {
    if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");
    printf("init_ine i=%d, j=%d: ", i, j);
    //TODO: adapt
    printf("%f\n", basis_set->radial_functions.cut(i - 1, j - 1));
    return basis_set->radial_functions.cut(i - 1, j - 1);
    //return basis_set->cutoffmax;
}

/* ---------------------------------------------------------------------- */

