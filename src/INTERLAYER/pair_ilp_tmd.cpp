/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Wengen Ouyang (Wuhan University)
   e-mail: w.g.ouyang at gmail dot com

   This is a full version of the potential described in
   [Ouyang et al., J. Chem. Theory Comput. 17, 7237 (2021).]
------------------------------------------------------------------------- */

#include "pair_ilp_tmd.h"

#include "atom.h"
#include "citeme.h"
#include "error.h"
#include "force.h"
#include "interlayer_taper.h"
#include "memory.h"
#include "my_page.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace InterLayer;

static const char cite_ilp_tmd[] =
    "ilp/tmd potential doi:10.1021/acs.jctc.1c00782\n"
    "@Article{Ouyang2021\n"
    "  author = {W. Ouyang and R. Sofer and X. Gao and J. Hermann and\n"
    "    A. Tkatchenko and L. Kronik and M. Urbakh and O. Hod},\n"
    "  title = {Anisotropic Interlayer Force Field for Transition\n"
    "    Metal Dichalcogenides: The Case of Molybdenum Disulfide},\n"
    "  journal = {J.~Chem.\\ Theory Comput.},\n"
    " volume   = 17,\n"
    " pages    = {7237--7245}\n"
    " year     = 2021,\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

PairILPTMD::PairILPTMD(LAMMPS *lmp) : PairILPGrapheneHBN(lmp)
{
  variant = ILP_TMD;
  single_enable = 0;

  // for TMD, each atom have six neighbors
  Nnei = 6;

  if (lmp->citeme) lmp->citeme->add(cite_ilp_tmd);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairILPTMD::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR, "Illegal pair_style command");
  if (!utils::strmatch(force->pair_style, "^hybrid/overlay"))
    error->all(FLERR, "Pair style ilp/tmd must be used as sub-style with hybrid/overlay");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (narg == 2) tap_flag = utils::numeric(FLERR, arg[1], false, lmp);
}

/* ----------------------------------------------------------------------
   Repulsive forces and energy
------------------------------------------------------------------------- */

void PairILPTMD::calc_FRep(int eflag, int /* vflag */)
{
  int i, j, ii, jj, inum, jnum, itype, jtype, k, kk;
  double prodnorm1, fkcx, fkcy, fkcz;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair, fpair1;
  double rsq, r, Rcut, rhosq1, exp0, exp1, Tap, dTap, Vilp;
  double frho1, Erep, fsum, rdsq1;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *ILP_neighs_i;

  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double delki[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  //calculate exp(-lambda*(r-z0))*[epsilon/2 + f(rho_ij)]
  // loop over neighbors of owned atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      // only include the interation between different layers
      if (rsq < cutsq[itype][jtype] && atom->molecule[i] != atom->molecule[j]) {

        int iparam_ij = elem2param[map[itype]][map[jtype]];
        Param &p = params[iparam_ij];

        r = sqrt(rsq);
        // turn on/off taper function
        if (tap_flag) {
          Rcut = sqrt(cutsq[itype][jtype]);
          Tap = calc_Tap(r, Rcut);
          dTap = calc_dTap(r, Rcut);
        } else {
          Tap = 1.0;
          dTap = 0.0;
        }

        // Calculate the transverse distance
        // note that rho_ij does not equal to rho_ji except when normals are all along z
        prodnorm1 = normal[i][0] * delx + normal[i][1] * dely + normal[i][2] * delz;
        rhosq1 = rsq - prodnorm1 * prodnorm1;    // rho_ij
        rdsq1 = rhosq1 * p.delta2inv;            // (rho_ij/delta)^2

        // store exponents
        exp0 = exp(-p.lambda * (r - p.z0));
        exp1 = exp(-rdsq1);

        frho1 = exp1 * p.C;
        Erep = 0.5 * p.epsilon + frho1;
        Vilp = exp0 * Erep;

        // derivatives
        fpair = p.lambda * exp0 / r * Erep;
        fpair1 = 2.0 * exp0 * frho1 * p.delta2inv;
        fsum = fpair + fpair1;
        // derivatives of the product of rij and ni, the result is a vector
        dprodnorm1[0] =
            dnormdri[i][0][0] * delx + dnormdri[i][1][0] * dely + dnormdri[i][2][0] * delz;
        dprodnorm1[1] =
            dnormdri[i][0][1] * delx + dnormdri[i][1][1] * dely + dnormdri[i][2][1] * delz;
        dprodnorm1[2] =
            dnormdri[i][0][2] * delx + dnormdri[i][1][2] * dely + dnormdri[i][2][2] * delz;
        fp1[0] = prodnorm1 * normal[i][0] * fpair1;
        fp1[1] = prodnorm1 * normal[i][1] * fpair1;
        fp1[2] = prodnorm1 * normal[i][2] * fpair1;
        fprod1[0] = prodnorm1 * dprodnorm1[0] * fpair1;
        fprod1[1] = prodnorm1 * dprodnorm1[1] * fpair1;
        fprod1[2] = prodnorm1 * dprodnorm1[2] * fpair1;

        fkcx = (delx * fsum - fp1[0]) * Tap - Vilp * dTap * delx / r;
        fkcy = (dely * fsum - fp1[1]) * Tap - Vilp * dTap * dely / r;
        fkcz = (delz * fsum - fp1[2]) * Tap - Vilp * dTap * delz / r;

        f[i][0] += fkcx - fprod1[0] * Tap;
        f[i][1] += fkcy - fprod1[1] * Tap;
        f[i][2] += fkcz - fprod1[2] * Tap;
        f[j][0] -= fkcx;
        f[j][1] -= fkcy;
        f[j][2] -= fkcz;

        // calculate the forces acted on the neighbors of atom i from atom j
        ILP_neighs_i = ILP_firstneigh[i];
        for (kk = 0; kk < ILP_numneigh[i]; kk++) {
          k = ILP_neighs_i[kk];
          if (k == i) continue;
          // derivatives of the product of rij and ni respect to rk, k=0,1,2, where atom k is the neighbors of atom i
          dprodnorm1[0] = dnormal[i][0][kk][0] * delx + dnormal[i][1][kk][0] * dely +
              dnormal[i][2][kk][0] * delz;
          dprodnorm1[1] = dnormal[i][0][kk][1] * delx + dnormal[i][1][kk][1] * dely +
              dnormal[i][2][kk][1] * delz;
          dprodnorm1[2] = dnormal[i][0][kk][2] * delx + dnormal[i][1][kk][2] * dely +
              dnormal[i][2][kk][2] * delz;
          fk[0] = (-prodnorm1 * dprodnorm1[0] * fpair1) * Tap;
          fk[1] = (-prodnorm1 * dprodnorm1[1] * fpair1) * Tap;
          fk[2] = (-prodnorm1 * dprodnorm1[2] * fpair1) * Tap;
          f[k][0] += fk[0];
          f[k][1] += fk[1];
          f[k][2] += fk[2];
          delki[0] = x[k][0] - x[i][0];
          delki[1] = x[k][1] - x[i][1];
          delki[2] = x[k][2] - x[i][2];
          if (evflag)
            ev_tally_xyz(k, i, nlocal, newton_pair, 0.0, 0.0, fk[0], fk[1], fk[2], delki[0],
                         delki[1], delki[2]);
        }

        if (eflag) pvector[1] += evdwl = Tap * Vilp;
        if (evflag)
          ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fkcx, fkcy, fkcz, delx, dely, delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create ILP neighbor list from main neighbor list to calculate normals
------------------------------------------------------------------------- */

void PairILPTMD::ILP_neigh()
{
  int i, j, l, ii, jj, ll, n, inum, jnum, itype, jtype, ltype, imol, jmol, count;
  double xtmp, ytmp, ztmp, delx, dely, delz, deljx, deljy, deljz, rsq, rsqlj;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *neighsort;
  int neighptr[10], check[10];

  double **x = atom->x;
  int *type = atom->type;

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(ILP_numneigh);
    memory->sfree(ILP_firstneigh);
    memory->create(ILP_numneigh, maxlocal, "ILPTMD:numneigh");
    ILP_firstneigh = (int **) memory->smalloc(maxlocal * sizeof(int *), "ILPTMD:firstneigh");
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // store all ILP neighs of owned and ghost atoms
  // scan full neighbor list of I

  ipage->reset();

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    //initialize varibles
    n = 0;
    neighsort = ipage->vget();
    for (ll = 0; ll < 10; ll++) {
      neighptr[ll] = -1;
      check[ll] = -1;
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    imol = atom->molecule[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      jmol = atom->molecule[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      // check if the atom i is a TMD atom, i.e., Mo/S/W/Se
      if (strcmp(elements[itype], "Mo") == 0 || strcmp(elements[itype], "W")  == 0 ||
          strcmp(elements[itype], "S")  == 0 || strcmp(elements[itype], "Se") == 0 ||
          strcmp(elements[itype], "Te") == 0) {
        if (rsq != 0 && rsq < cutILPsq[itype][jtype] && imol == jmol && type[i] == type[j]) {
          neighptr[n++] = j;
        }
      } else {    // atom i can be P, C, B, N or H.
        if (rsq != 0 && rsq < cutILPsq[itype][jtype] && imol == jmol) { neighptr[n++] = j; }
      }
    }    // loop over jj

    // if atom i is Mo/W/S/Se/Te, then sorting the orders of neighbors
    if (strcmp(elements[itype], "Mo") == 0 || strcmp(elements[itype], "W")  == 0 ||
        strcmp(elements[itype], "S")  == 0 || strcmp(elements[itype], "Se") == 0 ||
        strcmp(elements[itype], "Te") == 0) {
      // initialize neighsort
      for (ll = 0; ll < n; ll++) {
        neighsort[ll] = neighptr[ll];
        check[ll] = neighptr[ll];
      }

      // select the first neighbor of atomi
      if (n == Nnei) {
        neighsort[0] = neighptr[0];
        check[0] = -1;
      } else if (n < Nnei && n > 0) {
        for (jj = 0; jj < n; jj++) {    //identify the first neighbor
          j = neighptr[jj];
          jtype = map[type[j]];
          count = 0;
          for (ll = 0; ll < n; ll++) {
            l = neighptr[ll];
            ltype = map[type[l]];
            if (l == j) continue;
            deljx = x[l][0] - x[j][0];
            deljy = x[l][1] - x[j][1];
            deljz = x[l][2] - x[j][2];
            rsqlj = deljx * deljx + deljy * deljy + deljz * deljz;
            if (rsqlj != 0 && rsqlj < cutILPsq[ltype][jtype]) { count++; }
          }
          if (count == 1) {
            neighsort[0] = neighptr[jj];
            check[jj] = -1;
            break;
          }
        }    // end of idenfying the first neighbor
      } else if (n > Nnei) {
        error->one(FLERR,
                   "There are too many neighbors for TMD atoms, please check your configuration");
      }

      // sort the order of neighbors of atomi
      for (jj = 0; jj < n; jj++) {
        j = neighsort[jj];
        jtype = map[type[j]];
        ll = 0;
        while (ll < n) {
          l = neighptr[ll];
          if (check[ll] == -1) {
            ll++;
            continue;
          }
          ltype = map[type[l]];
          deljx = x[l][0] - x[j][0];
          deljy = x[l][1] - x[j][1];
          deljz = x[l][2] - x[j][2];
          rsqlj = deljx * deljx + deljy * deljy + deljz * deljz;
          if (rsqlj != 0 && rsqlj < cutILPsq[ltype][jtype]) {
            neighsort[jj + 1] = l;
            check[ll] = -1;
            break;
          }
          ll++;
        }
      }         // end of sorting the order of neighbors
    } else {    // for P/B/N/C/H atoms
      if (n > 3)
        error->one(
            FLERR,
            "There are too many neighbors for P/B/N/C/H atoms, please check your configuration");
      for (ll = 0; ll < n; ll++) { neighsort[ll] = neighptr[ll]; }
    }

    ILP_firstneigh[i] = neighsort;
    ILP_numneigh[i] = n;

    ipage->vgot(n);
    if (ipage->status()) error->one(FLERR, "Neighbor list overflow, boost neigh_modify one");
  }
}

/* ----------------------------------------------------------------------
   Calculate the normals for each atom
------------------------------------------------------------------------- */
void PairILPTMD::calc_normal()
{
  int i, j, ii, jj, inum, jnum;
  int cont, id, ip, m, k, itype;
  int *ilist, *jlist;
  int iH1,iH2,jH1,jH2;
  double Nave[3], dni[3], dpvdri[3][3];
  double nn, xtp, ytp, ztp, delx, dely, delz, nn2;

  double **x = atom->x;
  int *type = atom->type;
  tagint *tag = atom->tag;

  memory->destroy(dnn);
  memory->destroy(vect);
  memory->destroy(pvet);
  memory->destroy(dpvet1);
  memory->destroy(dpvet2);
  memory->destroy(dNave);
  memory->create(dnn, Nnei, 3, "ILPTMD:dnn");
  memory->create(vect, Nnei, 3, "ILPTMD:vect");
  memory->create(pvet, Nnei, 3, "ILPTMD:pvet");
  memory->create(dpvet1, Nnei, 3, 3, "ILPTMD:dpvet1");
  memory->create(dpvet2, Nnei, 3, 3, "ILPTMD:dpvet2");
  memory->create(dNave, 3, Nnei, 3, "ILPTMD:dNave");

  // grow normal array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(normal);
    memory->destroy(dnormal);
    memory->destroy(dnormdri);
    nmax = atom->nmax;
    memory->create(normal, nmax, 3, "ILPTMD:normal");
    memory->create(dnormdri, nmax, 3, 3, "ILPTMD:dnormdri");
    memory->create(dnormal, nmax, 3, Nnei, 3, "ILPTMD:dnormal");
  }

  inum = list->inum;
  ilist = list->ilist;
  //Calculate normals
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = map[type[i]];
    xtp = x[i][0];
    ytp = x[i][1];
    ztp = x[i][2];

    //   Initialize the arrays
    for (id = 0; id < 3; id++) {
      Nave[id] = 0.0;
      dni[id] = 0.0;
      normal[i][id] = 0.0;
      for (ip = 0; ip < 3; ip++) {
        dpvdri[ip][id] = 0.0;
        dnormdri[i][ip][id] = 0.0;
        for (m = 0; m < Nnei; m++) {
          dnn[m][id] = 0.0;
          vect[m][id] = 0.0;
          pvet[m][id] = 0.0;
          dpvet1[m][ip][id] = 0.0;
          dpvet2[m][ip][id] = 0.0;
          dNave[id][m][ip] = 0.0;
          dnormal[i][id][m][ip] = 0.0;
        }
      }
    }

    cont = 0;
    jlist = ILP_firstneigh[i];
    jnum = ILP_numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = x[j][0] - xtp;
      dely = x[j][1] - ytp;
      delz = x[j][2] - ztp;
      vect[cont][0] = delx;
      vect[cont][1] = dely;
      vect[cont][2] = delz;
      cont++;
    }

    //############################ For the dangling atoms ############################
    if (cont <= 1) {
      normal[i][0] = 0.0;
      normal[i][1] = 0.0;
      normal[i][2] = 1.0;
      for (id = 0; id < 3; id++) {
        for (ip = 0; ip < 3; ip++) {
          dnormdri[i][id][ip] = 0.0;
          for (m = 0; m < Nnei; m++) { dnormal[i][id][m][ip] = 0.0; }
        }
      }
      // for hydrogen in water molecule
      if (strcmp(elements[itype], "Hw") == 0) {
        if(cont == 0) {
          jH1 = atom->map(tag[i] - 1);
          jH2 = atom->map(tag[i] - 2);
          iH1 = map[type[jH1]];
          iH2 = map[type[jH2]];
          if (strcmp(elements[iH1], "Ow") == 0) {
            vect[0][0] = x[jH1][0] - xtp;
            vect[0][1] = x[jH1][1] - ytp;
            vect[0][2] = x[jH1][2] - ztp;
          } else if (strcmp(elements[iH2], "Ow") == 0) {
            vect[0][0] = x[jH2][0] - xtp;
            vect[0][1] = x[jH2][1] - ytp;
            vect[0][2] = x[jH2][2] - ztp;
          } else {
            error->one(FLERR, "The order of atoms in water molecule should be O H H !");
          }
        }
        Nave[0] = vect[0][0];
        Nave[1] = vect[0][1];
        Nave[2] = vect[0][2];
        // the magnitude of the normal vector
        nn2 = Nave[0] * Nave[0] + Nave[1] * Nave[1] + Nave[2] * Nave[2];
        nn = sqrt(nn2);
        if (nn == 0) error->one(FLERR, "The magnitude of the normal vector is zero");
        // the unit normal vector
        normal[i][0] = Nave[0] / nn;
        normal[i][1] = Nave[1] / nn;
        normal[i][2] = Nave[2] / nn;

        // Calculte dNave/dri, defined as dpvdri
        for (id = 0; id < 3; id++) {
          for (ip = 0; ip < 3; ip++) {
            if (ip == id) { dpvdri[id][ip] = -1.0;}
            else {dpvdri[id][ip] = 0.0;}
          }
        }

        // derivatives of nn, dnn:3x1 vector
        dni[0] = (Nave[0] * dpvdri[0][0] + Nave[1] * dpvdri[1][0] + Nave[2] * dpvdri[2][0]) / nn;
        dni[1] = (Nave[0] * dpvdri[0][1] + Nave[1] * dpvdri[1][1] + Nave[2] * dpvdri[2][1]) / nn;
        dni[2] = (Nave[0] * dpvdri[0][2] + Nave[1] * dpvdri[1][2] + Nave[2] * dpvdri[2][2]) / nn;
        // derivatives of unit vector ni respect to ri, the result is 3x3 matrix
        for (id = 0; id < 3; id++) {
          for (ip = 0; ip < 3; ip++) {
            dnormdri[i][id][ip] = dpvdri[id][ip] / nn - Nave[id] * dni[ip] / nn2;
            dnormal[i][id][0][ip] = -dnormdri[i][id][ip];
          }
        }
      }
    }
    //############################ For the edge atoms of TMD ################################
    else if (cont > 1 && cont < Nnei) {
      if (strcmp(elements[itype], "Mo") == 0 || strcmp(elements[itype], "W") == 0 ||
          strcmp(elements[itype], "S") == 0 || strcmp(elements[itype], "Se") == 0) {
        // derivatives of Ni[l] respect to the cont neighbors
        for (k = 0; k < cont - 1; k++) {
          for (ip = 0; ip < 3; ip++) {
            pvet[k][ip] = vect[k][modulo(ip + 1, 3)] * vect[k + 1][modulo(ip + 2, 3)] -
                vect[k][modulo(ip + 2, 3)] * vect[k + 1][modulo(ip + 1, 3)];
          }
          // dpvet1[k][l][ip]: the derivatve of the k (=0,...cont-1)th Nik respect to the ip component of atom l
          // derivatives respect to atom l
          // dNik,x/drl
          dpvet1[k][0][0] = 0.0;
          dpvet1[k][0][1] = vect[modulo(k + 1, Nnei)][2];
          dpvet1[k][0][2] = -vect[modulo(k + 1, Nnei)][1];
          // dNik,y/drl
          dpvet1[k][1][0] = -vect[modulo(k + 1, Nnei)][2];
          dpvet1[k][1][1] = 0.0;
          dpvet1[k][1][2] = vect[modulo(k + 1, Nnei)][0];
          // dNik,z/drl
          dpvet1[k][2][0] = vect[modulo(k + 1, Nnei)][1];
          dpvet1[k][2][1] = -vect[modulo(k + 1, Nnei)][0];
          dpvet1[k][2][2] = 0.0;

          // dpvet2[k][l][ip]: the derivatve of the k (=0,...cont-1)th Nik respect to the ip component of atom l+1
          // derivatives respect to atom l+1
          // dNik,x/drl+1
          dpvet2[k][0][0] = 0.0;
          dpvet2[k][0][1] = -vect[modulo(k, Nnei)][2];
          dpvet2[k][0][2] = vect[modulo(k, Nnei)][1];
          // dNik,y/drl+1
          dpvet2[k][1][0] = vect[modulo(k, Nnei)][2];
          dpvet2[k][1][1] = 0.0;
          dpvet2[k][1][2] = -vect[modulo(k, Nnei)][0];
          // dNik,z/drl+1
          dpvet2[k][2][0] = -vect[modulo(k, Nnei)][1];
          dpvet2[k][2][1] = vect[modulo(k, Nnei)][0];
          dpvet2[k][2][2] = 0.0;
        }

        // average the normal vectors by using the Nnei neighboring planes
        for (ip = 0; ip < 3; ip++) {
          Nave[ip] = 0.0;
          for (k = 0; k < cont - 1; k++) { Nave[ip] += pvet[k][ip]; }
          Nave[ip] /= (cont - 1);
        }
        // the magnitude of the normal vector
        nn2 = Nave[0] * Nave[0] + Nave[1] * Nave[1] + Nave[2] * Nave[2];
        nn = sqrt(nn2);
        if (nn == 0) error->one(FLERR, "The magnitude of the normal vector is zero");
        // the unit normal vector
        normal[i][0] = Nave[0] / nn;
        normal[i][1] = Nave[1] / nn;
        normal[i][2] = Nave[2] / nn;

        // derivatives of non-normalized normal vector, dNave:3xcontx3 array
        // dNave[id][m][ip]: the derivatve of the id component of Nave respect to the ip component of atom m
        for (id = 0; id < 3; id++) {
          for (ip = 0; ip < 3; ip++) {
            for (m = 0; m < cont; m++) {
              if (m == 0) {
                dNave[id][m][ip] = dpvet1[m][id][ip] / (cont - 1);
              } else if (m == cont - 1) {
                dNave[id][m][ip] = dpvet2[m - 1][id][ip] / (cont - 1);
              } else {    // sum of the derivatives of the mth and (m-1)th normal vector respect to the atom m
                dNave[id][m][ip] = (dpvet1[m][id][ip] + dpvet2[m - 1][id][ip]) / (cont - 1);
              }
            }
          }
        }
        // derivatives of nn, dnn:contx3 vector
        // dnn[m][id]: the derivative of nn respect to r[m][id], m=0,...Nnei-1; id=0,1,2
        // r[m][id]: the id's component of atom m
        for (m = 0; m < cont; m++) {
          for (id = 0; id < 3; id++) {
            dnn[m][id] = (Nave[0] * dNave[0][m][id] + Nave[1] * dNave[1][m][id] +
                          Nave[2] * dNave[2][m][id]) / nn;
          }
        }
        // dnormal[i][id][m][ip]: the derivative of normal[i][id] respect to r[m][ip], id,ip=0,1,2.
        // for atom m, which is a neighbor atom of atom i, m = 0,...,Nnei-1
        for (m = 0; m < cont; m++) {
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) {
              dnormal[i][id][m][ip] = dNave[id][m][ip] / nn - Nave[id] * dnn[m][ip] / nn2;
            }
          }
        }
        // Calculte dNave/dri, defined as dpvdri
        for (id = 0; id < 3; id++) {
          for (ip = 0; ip < 3; ip++) {
            dpvdri[id][ip] = 0.0;
            for (k = 0; k < cont; k++) { dpvdri[id][ip] -= dNave[id][k][ip]; }
          }
        }

        // derivatives of nn, dnn:3x1 vector
        dni[0] = (Nave[0] * dpvdri[0][0] + Nave[1] * dpvdri[1][0] + Nave[2] * dpvdri[2][0]) / nn;
        dni[1] = (Nave[0] * dpvdri[0][1] + Nave[1] * dpvdri[1][1] + Nave[2] * dpvdri[2][1]) / nn;
        dni[2] = (Nave[0] * dpvdri[0][2] + Nave[1] * dpvdri[1][2] + Nave[2] * dpvdri[2][2]) / nn;
        // derivatives of unit vector ni respect to ri, the result is 3x3 matrix
        for (id = 0; id < 3; id++) {
          for (ip = 0; ip < 3; ip++) {
            dnormdri[i][id][ip] = dpvdri[id][ip] / nn - Nave[id] * dni[ip] / nn2;
          }
        }
      }    // for TMD
      //############################ For Oxygen in the water molecule #######################
      else if (strcmp(elements[itype], "Ow") == 0) {
        if(cont == 0) {
          jH1 = atom->map(tag[i] + 1);
          jH2 = atom->map(tag[i] + 2);
          iH1 = map[type[jH1]];
          iH2 = map[type[jH2]];
          if (strcmp(elements[iH1], "Hw") == 0 && strcmp(elements[iH2], "Hw") == 0) {
            vect[0][0] = x[jH1][0] - xtp;
            vect[0][1] = x[jH1][1] - ytp;
            vect[0][2] = x[jH1][2] - ztp;

            vect[1][0] = x[jH2][0] - xtp;
            vect[1][1] = x[jH2][1] - ytp;
            vect[1][2] = x[jH2][2] - ztp;

            cont = 2;
          } else {
            error->one(FLERR, "The order of atoms in water molecule should be O H H !");
          }
        }
        if (cont == 2) {
          Nave[0] = (vect[0][0] + vect[1][0])/cont;
          Nave[1] = (vect[0][1] + vect[1][1])/cont;
          Nave[2] = (vect[0][2] + vect[1][2])/cont;
          // the magnitude of the normal vector
          nn2 = Nave[0] * Nave[0] + Nave[1] * Nave[1] + Nave[2] * Nave[2];
          nn = sqrt(nn2);
          if (nn == 0) error->one(FLERR, "The magnitude of the normal vector is zero");
          // the unit normal vector
          normal[i][0] = Nave[0] / nn;
          normal[i][1] = Nave[1] / nn;
          normal[i][2] = Nave[2] / nn;

          // derivatives of non-normalized normal vector, dNave:3xcontx3 array
          // dNave[id][m][ip]: the derivatve of the id component of Nave
          // respect to the ip component of atom m
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) {
              for (m = 0; m < cont; m++) {
                if (ip == id) { dNave[id][m][ip] = 0.5;}
                else {dNave[id][m][ip] = 0.0;}
              }
            }
          }
          // derivatives of nn, dnn:contx3 vector
          // dnn[m][id]: the derivative of nn respect to r[m][id], m=0,...Nnei-1; id=0,1,2
          // r[m][id]: the id's component of atom m
          for (m = 0; m < cont; m++) {
            for (id = 0; id < 3; id++) {
              dnn[m][id] = (Nave[0] * dNave[0][m][id] + Nave[1] * dNave[1][m][id] +
                            Nave[2] * dNave[2][m][id]) / nn;
            }
          }
          // dnormal[i][id][m][ip]: the derivative of normal[i][id] respect to r[m][ip], id,ip=0,1,2.
          // for atom m, which is a neighbor atom of atom i, m = 0,...,Nnei-1
          for (m = 0; m < cont; m++) {
            for (id = 0; id < 3; id++) {
              for (ip = 0; ip < 3; ip++) {
                dnormal[i][id][m][ip] = dNave[id][m][ip] / nn - Nave[id] * dnn[m][ip] / nn2;
              }
            }
          }
          // Calculte dNave/dri, defined as dpvdri
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) {
              dpvdri[id][ip] = 0.0;
              for (k = 0; k < cont; k++) { dpvdri[id][ip] -= dNave[id][k][ip]; }
            }
          }

          // derivatives of nn, dnn:3x1 vector
          dni[0] = (Nave[0] * dpvdri[0][0] + Nave[1] * dpvdri[1][0] + Nave[2] * dpvdri[2][0]) / nn;
          dni[1] = (Nave[0] * dpvdri[0][1] + Nave[1] * dpvdri[1][1] + Nave[2] * dpvdri[2][1]) / nn;
          dni[2] = (Nave[0] * dpvdri[0][2] + Nave[1] * dpvdri[1][2] + Nave[2] * dpvdri[2][2]) / nn;
          // derivatives of unit vector ni respect to ri, the result is 3x3 matrix
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) {
              dnormdri[i][id][ip] = dpvdri[id][ip] / nn - Nave[id] * dni[ip] / nn2;
            }
          }
        }
        else if (cont >= 3) {
          error->one(FLERR,
                     "There are too many neighbors for calculating normals of water molecules");
        }
      }
      //############################ For the edge & bulk atoms of GrhBN ################################
      else {
        if (cont == 2) {
          for (ip = 0; ip < 3; ip++) {
            pvet[0][ip] = vect[0][modulo(ip + 1, 3)] * vect[1][modulo(ip + 2, 3)] -
                          vect[0][modulo(ip + 2, 3)] * vect[1][modulo(ip + 1, 3)];
          }
          // dpvet1[k][l][ip]: the derivatve of the k (=0,...cont-1)th Nik respect to the ip component of atom l
          // derivatives respect to atom l
          // dNik,x/drl
          dpvet1[0][0][0] = 0.0;
          dpvet1[0][0][1] = vect[1][2];
          dpvet1[0][0][2] = -vect[1][1];
          // dNi0,y/drl
          dpvet1[0][1][0] = -vect[1][2];
          dpvet1[0][1][1] = 0.0;
          dpvet1[0][1][2] = vect[1][0];
          // dNi0,z/drl
          dpvet1[0][2][0] = vect[1][1];
          dpvet1[0][2][1] = -vect[1][0];
          dpvet1[0][2][2] = 0.0;

          // dpvet2[0][l][ip]: the derivatve of the 0 (=0,...cont-1)th Ni0 respect to the ip component of atom l+1
          // derivatives respect to atom l+1
          // dNi0,x/drl+1
          dpvet2[0][0][0] = 0.0;
          dpvet2[0][0][1] = -vect[0][2];
          dpvet2[0][0][2] = vect[0][1];
          // dNi0,y/drl+1
          dpvet2[0][1][0] = vect[0][2];
          dpvet2[0][1][1] = 0.0;
          dpvet2[0][1][2] = -vect[0][0];
          // dNi0,z/drl+1
          dpvet2[0][2][0] = -vect[0][1];
          dpvet2[0][2][1] = vect[0][0];
          dpvet2[0][2][2] = 0.0;

          for (ip = 0; ip < 3; ip++) { Nave[ip] += pvet[0][ip]; }
          // the magnitude of the normal vector
          nn2 = Nave[0] * Nave[0] + Nave[1] * Nave[1] + Nave[2] * Nave[2];
          nn = sqrt(nn2);
          if (nn == 0) error->one(FLERR, "The magnitude of the normal vector is zero");
          // the unit normal vector
          normal[i][0] = Nave[0] / nn;
          normal[i][1] = Nave[1] / nn;
          normal[i][2] = Nave[2] / nn;

          // derivatives of non-normalized normal vector, dNave:3xcontx3 array
          // dNave[id][m][ip]: the derivatve of the id component of Nave respect to the ip component of atom m
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) {
              for (m = 0; m < cont; m++) {
                if (m == 0) {
                  dNave[id][m][ip] = dpvet1[m][id][ip] / (cont - 1);
                } else if (m == cont - 1) {
                  dNave[id][m][ip] = dpvet2[m - 1][id][ip] / (cont - 1);
                } else {    // sum of the derivatives of the mth and (m-1)th normal vector respect to the atom m
                  dNave[id][m][ip] = (dpvet1[m][id][ip] + dpvet2[m - 1][id][ip]) / (cont - 1);
                }
              }
            }
          }
          // derivatives of nn, dnn:contx3 vector
          // dnn[m][id]: the derivative of nn respect to r[m][id], m=0,...Nnei-1; id=0,1,2
          // r[m][id]: the id's component of atom m
          for (m = 0; m < cont; m++) {
            for (id = 0; id < 3; id++) {
              dnn[m][id] = (Nave[0] * dNave[0][m][id] + Nave[1] * dNave[1][m][id] +
                            Nave[2] * dNave[2][m][id]) / nn;
            }
          }
          // dnormal[i][id][m][ip]: the derivative of normal[i][id] respect to r[m][ip], id,ip=0,1,2.
          // for atom m, which is a neighbor atom of atom i, m = 0,...,Nnei-1
          for (m = 0; m < cont; m++) {
            for (id = 0; id < 3; id++) {
              for (ip = 0; ip < 3; ip++) {
                dnormal[i][id][m][ip] = dNave[id][m][ip] / nn - Nave[id] * dnn[m][ip] / nn2;
              }
            }
          }
          // Calculte dNave/dri, defined as dpvdri
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) {
              dpvdri[id][ip] = 0.0;
              for (k = 0; k < cont; k++) { dpvdri[id][ip] -= dNave[id][k][ip]; }
            }
          }
          // derivatives of nn, dnn:3x1 vector
          dni[0] = (Nave[0] * dpvdri[0][0] + Nave[1] * dpvdri[1][0] + Nave[2] * dpvdri[2][0]) / nn;
          dni[1] = (Nave[0] * dpvdri[0][1] + Nave[1] * dpvdri[1][1] + Nave[2] * dpvdri[2][1]) / nn;
          dni[2] = (Nave[0] * dpvdri[0][2] + Nave[1] * dpvdri[1][2] + Nave[2] * dpvdri[2][2]) / nn;
          // derivatives of unit vector ni respect to ri, the result is 3x3 matrix
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) {
              dnormdri[i][id][ip] = dpvdri[id][ip] / nn - Nave[id] * dni[ip] / nn2;
            }
          }
        }    // end of cont == 2
        else if (cont == 3) {
          // derivatives of Ni[l] respect to the 3 neighbors
          for (k = 0; k < 3; k++) {
            for (ip = 0; ip < 3; ip++) {
              pvet[k][ip] = vect[modulo(k, 3)][modulo(ip + 1, 3)] *
                      vect[modulo(k + 1, 3)][modulo(ip + 2, 3)] -
                  vect[modulo(k, 3)][modulo(ip + 2, 3)] * vect[modulo(k + 1, 3)][modulo(ip + 1, 3)];
            }
            // dpvet1[k][l][ip]: the derivatve of the k (=0,...cont-1)th Nik respect to the ip component of atom l
            // derivatives respect to atom l
            // dNik,x/drl
            dpvet1[k][0][0] = 0.0;
            dpvet1[k][0][1] = vect[modulo(k + 1, 3)][2];
            dpvet1[k][0][2] = -vect[modulo(k + 1, 3)][1];
            // dNik,y/drl
            dpvet1[k][1][0] = -vect[modulo(k + 1, 3)][2];
            dpvet1[k][1][1] = 0.0;
            dpvet1[k][1][2] = vect[modulo(k + 1, 3)][0];
            // dNik,z/drl
            dpvet1[k][2][0] = vect[modulo(k + 1, 3)][1];
            dpvet1[k][2][1] = -vect[modulo(k + 1, 3)][0];
            dpvet1[k][2][2] = 0.0;

            // dpvet2[k][l][ip]: the derivatve of the k (=0,...cont-1)th Nik respect to the ip component of atom l+1
            // derivatives respect to atom l+1
            // dNik,x/drl+1
            dpvet2[k][0][0] = 0.0;
            dpvet2[k][0][1] = -vect[modulo(k, 3)][2];
            dpvet2[k][0][2] = vect[modulo(k, 3)][1];
            // dNik,y/drl+1
            dpvet2[k][1][0] = vect[modulo(k, 3)][2];
            dpvet2[k][1][1] = 0.0;
            dpvet2[k][1][2] = -vect[modulo(k, 3)][0];
            // dNik,z/drl+1
            dpvet2[k][2][0] = -vect[modulo(k, 3)][1];
            dpvet2[k][2][1] = vect[modulo(k, 3)][0];
            dpvet2[k][2][2] = 0.0;
          }

          // average the normal vectors by using the 3 neighboring planes
          for (ip = 0; ip < 3; ip++) {
            Nave[ip] = 0.0;
            for (k = 0; k < 3; k++) { Nave[ip] += pvet[k][ip]; }
            Nave[ip] /= 3;
          }
          // the magnitude of the normal vector
          nn2 = Nave[0] * Nave[0] + Nave[1] * Nave[1] + Nave[2] * Nave[2];
          nn = sqrt(nn2);
          if (nn == 0) error->one(FLERR, "The magnitude of the normal vector is zero");
          // the unit normal vector
          normal[i][0] = Nave[0] / nn;
          normal[i][1] = Nave[1] / nn;
          normal[i][2] = Nave[2] / nn;

          // for the central atoms, dnormdri is always zero
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) { dnormdri[i][id][ip] = 0.0; }
          }

          // derivatives of non-normalized normal vector, dNave:3x3x3 array
          // dNave[id][m][ip]: the derivatve of the id component of Nave respect to the ip component of atom m
          for (id = 0; id < 3; id++) {
            for (ip = 0; ip < 3; ip++) {
              for (
                  m = 0; m < 3;
                  m++) {    // sum of the derivatives of the mth and (m-1)th normal vector respect to the atom m
                dNave[id][m][ip] =
                    (dpvet1[modulo(m, 3)][id][ip] + dpvet2[modulo(m - 1, 3)][id][ip]) / 3;
              }
            }
          }
          // derivatives of nn, dnn:3x3 vector
          // dnn[m][id]: the derivative of nn respect to r[m][id], m=0,...3-1; id=0,1,2
          // r[m][id]: the id's component of atom m
          for (m = 0; m < 3; m++) {
            for (id = 0; id < 3; id++) {
              dnn[m][id] = (Nave[0] * dNave[0][m][id] + Nave[1] * dNave[1][m][id] +
                            Nave[2] * dNave[2][m][id]) /
                  nn;
            }
          }
          // dnormal[i][id][m][ip]: the derivative of normal[i][id] respect to r[m][ip], id,ip=0,1,2.
          // for atom m, which is a neighbor atom of atom i, m = 0,...,3-1
          for (m = 0; m < 3; m++) {
            for (id = 0; id < 3; id++) {
              for (ip = 0; ip < 3; ip++) {
                dnormal[i][id][m][ip] = dNave[id][m][ip] / nn - Nave[id] * dnn[m][ip] / nn2;
              }
            }
          }
        }    // end of cont == 3
        else
          error->one(FLERR,
                     "There are too many neighbors for calculating normals of B/N/C/H atoms");
      }    // for B/N/C/H
    }      // end of if(cont<Nnei)
           //############################ For the bulk atoms ################################
    else if (cont == Nnei) {
      // derivatives of Ni[l] respect to the Nnei neighbors
      for (k = 0; k < Nnei; k++) {
        for (ip = 0; ip < 3; ip++) {
          pvet[k][ip] = vect[modulo(k, Nnei)][modulo(ip + 1, 3)] *
                  vect[modulo(k + 1, Nnei)][modulo(ip + 2, 3)] -
              vect[modulo(k, Nnei)][modulo(ip + 2, 3)] *
                  vect[modulo(k + 1, Nnei)][modulo(ip + 1, 3)];
        }
        // dpvet1[k][l][ip]: the derivatve of the k (=0,...cont-1)th Nik respect to the ip component of atom l
        // derivatives respect to atom l
        // dNik,x/drl
        dpvet1[k][0][0] = 0.0;
        dpvet1[k][0][1] = vect[modulo(k + 1, Nnei)][2];
        dpvet1[k][0][2] = -vect[modulo(k + 1, Nnei)][1];
        // dNik,y/drl
        dpvet1[k][1][0] = -vect[modulo(k + 1, Nnei)][2];
        dpvet1[k][1][1] = 0.0;
        dpvet1[k][1][2] = vect[modulo(k + 1, Nnei)][0];
        // dNik,z/drl
        dpvet1[k][2][0] = vect[modulo(k + 1, Nnei)][1];
        dpvet1[k][2][1] = -vect[modulo(k + 1, Nnei)][0];
        dpvet1[k][2][2] = 0.0;

        // dpvet2[k][l][ip]: the derivatve of the k (=0,...cont-1)th Nik respect to the ip component of atom l+1
        // derivatives respect to atom l+1
        // dNik,x/drl+1
        dpvet2[k][0][0] = 0.0;
        dpvet2[k][0][1] = -vect[modulo(k, Nnei)][2];
        dpvet2[k][0][2] = vect[modulo(k, Nnei)][1];
        // dNik,y/drl+1
        dpvet2[k][1][0] = vect[modulo(k, Nnei)][2];
        dpvet2[k][1][1] = 0.0;
        dpvet2[k][1][2] = -vect[modulo(k, Nnei)][0];
        // dNik,z/drl+1
        dpvet2[k][2][0] = -vect[modulo(k, Nnei)][1];
        dpvet2[k][2][1] = vect[modulo(k, Nnei)][0];
        dpvet2[k][2][2] = 0.0;
      }

      // average the normal vectors by using the Nnei neighboring planes
      for (ip = 0; ip < 3; ip++) {
        Nave[ip] = 0.0;
        for (k = 0; k < Nnei; k++) { Nave[ip] += pvet[k][ip]; }
        Nave[ip] /= Nnei;
      }
      // the magnitude of the normal vector
      nn2 = Nave[0] * Nave[0] + Nave[1] * Nave[1] + Nave[2] * Nave[2];
      nn = sqrt(nn2);
      if (nn == 0.0) error->one(FLERR, "The magnitude of the normal vector is zero");
      // the unit normal vector
      normal[i][0] = Nave[0] / nn;
      normal[i][1] = Nave[1] / nn;
      normal[i][2] = Nave[2] / nn;

      // for the central atoms, dnormdri is always zero
      for (id = 0; id < 3; id++) {
        for (ip = 0; ip < 3; ip++) { dnormdri[i][id][ip] = 0.0; }
      }

      // derivatives of non-normalized normal vector, dNave:3xNneix3 array
      // dNave[id][m][ip]: the derivatve of the id component of Nave respect to the ip component of atom m
      for (id = 0; id < 3; id++) {
        for (ip = 0; ip < 3; ip++) {
          for (
              m = 0; m < Nnei;
              m++) {    // sum of the derivatives of the mth and (m-1)th normal vector respect to the atom m
            dNave[id][m][ip] =
                (dpvet1[modulo(m, Nnei)][id][ip] + dpvet2[modulo(m - 1, Nnei)][id][ip]) / Nnei;
          }
        }
      }
      // derivatives of nn, dnn:Nneix3 vector
      // dnn[m][id]: the derivative of nn respect to r[m][id], m=0,...Nnei-1; id=0,1,2
      // r[m][id]: the id's component of atom m
      for (m = 0; m < Nnei; m++) {
        for (id = 0; id < 3; id++) {
          dnn[m][id] =
              (Nave[0] * dNave[0][m][id] + Nave[1] * dNave[1][m][id] + Nave[2] * dNave[2][m][id]) /
              nn;
        }
      }
      // dnormal[i][id][m][ip]: the derivative of normal[i][id] respect to r[m][ip], id,ip=0,1,2.
      // for atom m, which is a neighbor atom of atom i, m = 0,...,Nnei-1
      for (m = 0; m < Nnei; m++) {
        for (id = 0; id < 3; id++) {
          for (ip = 0; ip < 3; ip++) {
            dnormal[i][id][m][ip] = dNave[id][m][ip] / nn - Nave[id] * dnn[m][ip] / nn2;
          }
        }
      }
    } else {
      error->one(FLERR, "There are too many neighbors for calculating normals of TMD atoms");
    }    // end of four cases of cont
  }      // end of i loop
}
