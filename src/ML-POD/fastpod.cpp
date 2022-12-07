/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ngoc Cuong Nguyen (MIT)
------------------------------------------------------------------------- */

// POD header file
#include "fastpod.h"

// LAMMPS header files

#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "tokenizer.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PI;
using MathSpecial::cube;
using MathSpecial::powint;

#define MAXLINE 1024

// constructor
FASTPOD::FASTPOD(LAMMPS *_lmp, const std::string &pod_file, const std::string &coeff_file) :
    Pointers(_lmp), elemindex(nullptr), Phi(nullptr), Lambda(nullptr), coeff(nullptr),
    newcoeff(nullptr), tmpmem(nullptr), tmpint(nullptr), pn3(nullptr), pq3(nullptr),
    pc3(nullptr), pq4(nullptr), pa4(nullptr), pb4(nullptr), pc4(nullptr), ind23(nullptr),
    ind32(nullptr), ind33(nullptr), ind34(nullptr), ind43(nullptr), ind44(nullptr)
{
  rin = 0.5;
  rcut = 5.0;
  nelements = 1;
  onebody = 1;
  besseldegree = 4;
  inversedegree = 8;
  nbesselpars = 3;
  ns = nbesselpars*besseldegree + inversedegree;
  Njmax = 100;
  nrbf2 = 6;
  nrbf3 = 5;
  nrbf4 = 4;
  nabf3 = 5;
  nabf4 = 4;
  nrbf23 = 0;
  nabf23 = 0;
  nrbf33 = 0;
  nabf33 = 0;
  nrbf34 = 0;
  nabf34 = 0;
  nabf43 = 0;
  nrbf44 = 0;
  nabf44 = 0;
  P3 = 4;
  P4 = 3;
  P23 = 0;
  P33 = 0;
  P34 = 0;
  P44 = 0;
  pdegree[0] = besseldegree;
  pdegree[1] = inversedegree;
  pbc[0] = 1;
  pbc[1] = 1;
  pbc[2] = 1;
  besselparams[0] = 1e-3;
  besselparams[1] = 2.0;
  besselparams[2] = 4.0;

  // read pod input file to podstruct
  read_pod_file(pod_file);

  // read pod coefficient file to podstruct
  if (coeff_file != "") {
    ncoeff = read_coeff_file(coeff_file);
    mknewcoeff();
  }
}

// destructor
FASTPOD::~FASTPOD()
{
  memory->destroy(elemindex);
  memory->destroy(Phi);
  memory->destroy(Lambda);
  memory->destroy(coeff);
  memory->destroy(newcoeff);
  memory->destroy(tmpmem);
  memory->destroy(tmpint);
  memory->destroy(pn3);
  memory->destroy(pq3);
  memory->destroy(pc3);
  memory->destroy(pa4);
  memory->destroy(pb4);
  memory->destroy(pc4);
  memory->destroy(pq4);
  memory->destroy(ind23);
  memory->destroy(ind32);
  memory->destroy(ind33);
  memory->destroy(ind34);
  memory->destroy(ind43);
  memory->destroy(ind44);
}

void FASTPOD::read_pod_file(std::string pod_file)
{
  std::string podfilename = pod_file;
  FILE *fppod;
  if (comm->me == 0) {

    fppod = utils::open_potential(podfilename,lmp,nullptr);
    if (fppod == nullptr)
      error->one(FLERR,"Cannot open POD coefficient file {}: ",
                                   podfilename, utils::getsyserror());
  }

  // loop through lines of POD file and parse keywords

  char line[MAXLINE],*ptr;
  int eof = 0;

  while (true) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fppod);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fppod);
      }
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    // words = ptrs to all words in line
    // strip single and double quotes from words

    std::vector<std::string> words;
    try {
      words = Tokenizer(utils::trim_comment(line),"\"' \t\n\r\f").as_vector();
    } catch (TokenizerException &) {
      // ignore
    }

    if (words.size() == 0) continue;

    auto keywd = words[0];

    if (keywd == "species") {
      nelements = words.size()-1;
      for (int ielem = 1; ielem <= nelements; ielem++) {
        species.push_back(words[ielem]);
      }
    }

    if (keywd == "pbc") {
      if (words.size() != 4)
        error->one(FLERR,"Improper POD file.", utils::getsyserror());
      pbc[0] = utils::inumeric(FLERR,words[1],false,lmp);
      pbc[1] = utils::inumeric(FLERR,words[2],false,lmp);
      pbc[2] = utils::inumeric(FLERR,words[3],false,lmp);
    }

    if ((keywd != "#") && (keywd != "species") && (keywd != "pbc")) {

      if (words.size() != 2)
        error->one(FLERR,"Improper POD file.", utils::getsyserror());

      if (keywd == "rin") rin = utils::numeric(FLERR,words[1],false,lmp);
      if (keywd == "rcut") rcut = utils::numeric(FLERR,words[1],false,lmp);
      if (keywd == "bessel_polynomial_degree")
        besseldegree = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "inverse_polynomial_degree")
        inversedegree = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "onebody") onebody = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "twobody_number_radial_basis_functions")
        nrbf2 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "threebody_number_radial_basis_functions")
        nrbf3 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "threebody_angular_degree")
        P3 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "fourbody_number_radial_basis_functions")
        nrbf4 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "fourbody_angular_degree")
        P4 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "fivebody_number_radial_basis_functions")
        nrbf33 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "fivebody_angular_degree")
        P33 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "sixbody_number_radial_basis_functions")
        nrbf34 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "sixbody_angular_degree")
        P34 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "sevenbody_number_radial_basis_functions")
        nrbf44 = utils::inumeric(FLERR,words[1],false,lmp);
      if (keywd == "sevenbody_angular_degree")
        P44 = utils::inumeric(FLERR,words[1],false,lmp);
    }
  }
  if (nrbf2 < nrbf3) error->all(FLERR,"number of three-body radial basis functions must be equal or less than number of two-body radial basis functions");
  if (nrbf3 < nrbf4) error->all(FLERR,"number of four-body radial basis functions must be equal or less than number of three-body radial basis functions");
  if (nrbf4 < nrbf33) error->all(FLERR,"number of five-body radial basis functions must be equal or less than number of four-body radial basis functions");
  if (nrbf33 < nrbf34) error->all(FLERR,"number of six-body radial basis functions must be equal or less than number of five-body radial basis functions");
  if (nrbf34 < nrbf44) error->all(FLERR,"number of seven-body radial basis functions must be equal or less than number of six-body radial basis functions");

  if (P3 < P4) error->all(FLERR,"four-body angular degree must be equal or less than three-body angular degree");
  if (P4 < P33) error->all(FLERR,"five-body angular degree must be equal or less than four-body angular degree");
  if (P33 < P34) error->all(FLERR,"six-body angular degree must be equal or less than five-body angular degree");
  if (P34 < P44) error->all(FLERR,"seven-body angular degree must be equal or less than six-body angular degree");

  if (P3 > 12) error->all(FLERR,"three-body angular degree must be equal or less than 12");
  if (P34 > 6) error->all(FLERR,"six-body angular degree must be equal or less than 6");
  if (P44 > 6) error->all(FLERR,"seven-body angular degree must be equal or less than 6");

  // four-body potential
  if ((nrbf4 > 0) && (nrbf33 == 0)) {
    if (P4 > 6) {
      nrbf23 = nrbf4;
      P23 = P4;
      nrbf4 = 0;
      P4 = 0;
    }
    else {
      nrbf23 = 0;
      P23 = 0;
    }
  }

  // five-body potential
  if ((nrbf33 > 0) && (nrbf34 == 0)) {
    nrbf23 = nrbf4;
    P23 = P4;
    nrbf4 = 0;
    P4 = 0;
  }

  // six-body potential or seven-body potential
  if (nrbf34 > 0) {
    nrbf23 = nrbf4;
    P23 = P4;
    nrbf4 = nrbf34;
    P4 = P34;
  }

  int Ne = nelements;

  memory->create(elemindex, Ne*Ne, "elemindex");
  int k = 0;
  for (int i1 = 0; i1<Ne; i1++)
    for (int i2 = i1; i2<Ne; i2++) {
      elemindex[i2 + Ne*i1] = k;
      elemindex[i1 + Ne*i2] = k;
      k += 1;
    }

  init2body();
  init3body(P3);
  init4body(P4);

  int nb[] = {1,     2,     4,     7,    11,    16,    23};
  nabf23 = P23+1;
  nabf33 = P33+1;
  nabf34 = P34+1;
  nabf43 = nb[P34];
  nabf44 = nb[P44];

  if (onebody==0)
    nd1 = 0;
  else
    nd1 = Ne;

  nl1 = nd1/Ne;
  nl2 = nrbf2*Ne;
  nl3 = nabf3*nrbf3*Ne*(Ne+1)/2;
  nl4 = nabf4*nrbf4*Ne*(Ne+1)*(Ne+2)/6;

  nd2 = nrbf2*Ne*(Ne+1)/2;
  nd3 = nabf3*nrbf3*Ne*Ne*(Ne+1)/2;
  nd4 = nabf4*nrbf4*Ne*Ne*(Ne+1)*(Ne+2)/6;

  n23 = nrbf23*Ne;
  n32 = nabf23*nrbf23*Ne*(Ne+1)/2;
  n33 = nabf33*nrbf33*Ne*(Ne+1)/2;
  n34 = nabf34*nrbf34*Ne*(Ne+1)/2;
  n43 = nabf43*nrbf34*Ne*(Ne+1)*(Ne+2)/6;
  n44 = nabf44*nrbf44*Ne*(Ne+1)*(Ne+2)/6;

  nl23 = n23*n32;
  nl34 = n34*n43;
  nl33 = n33*(n33+1)/2;
  nl44 = n44*(n44+1)/2;

  nd23 = nl23*Ne;
  nd33 = nl33*Ne;
  nd34 = nl34*Ne;
  nd44 = nl44*Ne;

  nl = nl1 + nl2 + nl3 + nl4 + nl23 + nl33 + nl34 + nl44;
  nd = nd1 + nd2 + nd3 + nd4 + nd23 + nd33 + nd34 + nd44;

  memory->create(ind23, n23, "ind23");
  memory->create(ind32, n32, "ind32");
  memory->create(ind33, n33, "ind33");
  memory->create(ind34, n34, "ind34");
  memory->create(ind43, n43, "ind43");
  memory->create(ind44, n44, "ind44");

  indexmap3(ind23, 1, nrbf23, Ne, 1, nrbf2);
  indexmap3(ind32, nabf23, nrbf23, Ne*(Ne+1)/2, nabf3, nrbf3);
  indexmap3(ind33, nabf33, nrbf33, Ne*(Ne+1)/2, nabf3, nrbf3);
  indexmap3(ind34, nabf34, nrbf34, Ne*(Ne+1)/2, nabf3, nrbf3);
  indexmap3(ind43, nabf43, nrbf34, Ne*(Ne+1)*(Ne+2)/6, nabf4, nrbf4);
  indexmap3(ind44, nabf44, nrbf44, Ne*(Ne+1)*(Ne+2)/6, nabf4, nrbf4);

  estimate_memory(Njmax);
  memory->create(tmpmem, ndblmem, "tmpmem");
  memory->create(tmpint, nintmem, "tmpint");

  if (comm->me == 0) {
    utils::logmesg(lmp, "**************** Begin of POD Potentials ****************\n");
    utils::logmesg(lmp, "species: ");
    for (int i=0; i<nelements; i++)
      utils::logmesg(lmp, "{} ", species[i]);
    utils::logmesg(lmp, "\n");
    utils::logmesg(lmp, "periodic boundary conditions: {} {} {}\n", pbc[0], pbc[1], pbc[2]);
    utils::logmesg(lmp, "inner cut-off radius: {}\n", rin);
    utils::logmesg(lmp, "outer cut-off radius: {}\n", rcut);
    utils::logmesg(lmp, "bessel polynomial degree: {}\n", besseldegree);
    utils::logmesg(lmp, "inverse polynomial degree: {}\n",inversedegree);
    utils::logmesg(lmp, "one-body potential: {}\n", onebody);
    utils::logmesg(lmp, "two-body radial basis functions: {}\n", nrbf2);
    utils::logmesg(lmp, "three-body radial basis functions: {}\n", nrbf3);
    utils::logmesg(lmp, "three-body angular degree: {}\n", P3);
    if (P23 < P4) {
      utils::logmesg(lmp, "four-body radial basis functions: {}\n", nrbf4);
      utils::logmesg(lmp, "four-body angular degree: {}\n", P4);
    }
    else {
      utils::logmesg(lmp, "four-body radial basis functions: {}\n", nrbf23);
      utils::logmesg(lmp, "four-body angular degree: {}\n", P23);
    }
    utils::logmesg(lmp, "five-body radial basis functions: {}\n", nrbf33);
    utils::logmesg(lmp, "five-body angular degree: {}\n", P33);
    utils::logmesg(lmp, "six-body radial basis functions: {}\n", nrbf34);
    utils::logmesg(lmp, "six-body angular degree: {}\n", P34);
    utils::logmesg(lmp, "seven-body radial basis functions: {}\n", nrbf44);
    utils::logmesg(lmp, "seven-body angular degree: {}\n", P44);
    utils::logmesg(lmp, "number of descriptors for one-body potential: {}\n", nd1);
    utils::logmesg(lmp, "number of descriptors for two-body potential: {}\n", nd2);
    utils::logmesg(lmp, "number of descriptors for three-body potential: {}\n", nd3);
    utils::logmesg(lmp, "number of descriptors for four-body potential: {}\n", nd4+nd23);
    utils::logmesg(lmp, "number of descriptors for five-body potential: {}\n", nd33);
    utils::logmesg(lmp, "number of descriptors for six-body potential: {}\n", nd34);
    utils::logmesg(lmp, "number of descriptors for seven-body potential: {}\n", nd44);
    utils::logmesg(lmp, "total number of descriptors for all potentials: {}\n", nd);
    utils::logmesg(lmp, "**************** End of POD Potentials ****************\n\n");
  }
}

int FASTPOD::read_coeff_file(std::string coeff_file)
{
  std::string coefffilename = coeff_file;
  FILE *fpcoeff;
  if (comm->me == 0) {

    fpcoeff = utils::open_potential(coefffilename,lmp,nullptr);
    if (fpcoeff == nullptr)
      error->one(FLERR,"Cannot open POD coefficient file {}: ",
                                   coefffilename, utils::getsyserror());
  }

  // check format for first line of file

  char line[MAXLINE],*ptr;
  int eof = 0;
  int nwords = 0;
  while (nwords == 0) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpcoeff);
      }
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    nwords = utils::count_words(utils::trim_comment(line));
  }

  if (nwords != 2)
    error->all(FLERR,"Incorrect format in POD coefficient file");

  // strip single and double quotes from words

  int ncoeffall;
  std::string tmp_str;
  try {
    ValueTokenizer words(utils::trim_comment(line),"\"' \t\n\r\f");
    tmp_str = words.next_string();
    ncoeffall = words.next_int();
  } catch (TokenizerException &e) {
    error->all(FLERR,"Incorrect format in POD coefficient file: {}", e.what());
  }

  // loop over single block of coefficients and insert values in coeff

  memory->create(coeff, ncoeffall, "pod:pod_coeff");

  for (int icoeff = 0; icoeff < ncoeffall; icoeff++) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fpcoeff);
      if (ptr == nullptr) {
        eof = 1;
        fclose(fpcoeff);
      }
    }

    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof)
      error->all(FLERR,"Incorrect format in POD coefficient file");
    MPI_Bcast(line,MAXLINE,MPI_CHAR,0,world);

    try {
      ValueTokenizer cff(utils::trim_comment(line));
      if (cff.count() != 1)
        error->all(FLERR,"Incorrect format in POD coefficient file");

      coeff[icoeff] = cff.next_double();
    } catch (TokenizerException &e) {
      error->all(FLERR,"Incorrect format in POD coefficient file: {}", e.what());
    }
  }
  if (comm->me == 0) {
    if (!eof) fclose(fpcoeff);
  }

  if (comm->me == 0) {
    utils::logmesg(lmp, "**************** Begin of POD Coefficients ****************\n");
    utils::logmesg(lmp, "total number of coefficients for POD potential: {}\n", ncoeffall);
    utils::logmesg(lmp, "**************** End of POD Potentials ****************\n\n");
  }

  return ncoeffall;
}

double FASTPOD::peratomenergyforce(double *fij, double *rij, double *temp,
        int *ti, int *tj, int Nj)
{
  double *coeff1 = &newcoeff[0];
  double *coeff2 = &newcoeff[nl1*nelements];
  double *coeff3 = &newcoeff[(nl1 + nl2)*nelements];
  double *coeff4 = &newcoeff[(nl1 + nl2 + nl3)*nelements];
  double *coeff23 = &newcoeff[(nl1 + nl2 + nl3 + nl4)*nelements];
  double *coeff33 = &newcoeff[(nl1 + nl2 + nl3 + nl4 + nl23)*nelements];
  double *coeff34 = &newcoeff[(nl1 + nl2 + nl3 + nl4 + nl23 + nl33)*nelements];
  double *coeff44 = &newcoeff[(nl1 + nl2 + nl3 + nl4 + nl23 + nl33 + nl34)*nelements];

  int t0 = ti[0]-1;
  int n1 = Nj*K3*nrbf3;
  int n2 = Nj*nrbf2;
  int n3 = Nj*ns;
  int n4 = Nj*K3;
  int n5 = K3*nrbf3*nelements;

  double *U = &temp[0]; // Nj*K3*nrbf3
  double *Ux = &temp[n1]; // Nj*K3*nrbf3
  double *Uy = &temp[2*n1]; // Nj*K3*nrbf3
  double *Uz = &temp[3*n1]; // Nj*K3*nrbf3
  double *sumU = &temp[4*n1]; // K3*nrbf3*nelements

  double *rbf = &temp[4*n1 + n5]; // Nj*nrbf2
  double *rbfx = &temp[4*n1 + n5 + n2]; // Nj*nrbf2
  double *rbfy = &temp[4*n1 + n5 + 2*n2]; // Nj*nrbf2
  double *rbfz = &temp[4*n1 + n5 + 3*n2]; // Nj*nrbf2

  double *rbft = &temp[4*n1 + n5 + 4*n2]; // Nj*ns
  double *rbfxt = &temp[4*n1 + n5 + 4*n2 + n3]; // Nj*ns
  double *rbfyt = &temp[4*n1 + n5 + 4*n2 + 2*n3]; // Nj*ns
  double *rbfzt = &temp[4*n1 + n5 + 4*n2 + 3*n3]; // Nj*ns

  // orthogonal radial basis functions
  radialbasis(rbft, rbfxt, rbfyt, rbfzt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nj);
  MatMul(rbf, rbft, Phi, Nj, ns, nrbf2);
  MatMul(rbfx, rbfxt, Phi, Nj, ns, nrbf2);
  MatMul(rbfy, rbfyt, Phi, Nj, ns, nrbf2);
  MatMul(rbfz, rbfzt, Phi, Nj, ns, nrbf2);

  for (int j=0; j<3*Nj; j++) fij[j] = 0.0;

  double e1=0, e2=0, e3=0, e4=0, e23=0, e33=0, e34=0, e44=0;

  e1 = coeff1[t0];
  e2 = tallytwobodylocalforce(fij, &coeff2[nl2*t0], rbf, rbfx, rbfy, rbfz, tj, nrbf2, Nj);

  if (nd3 > 0) {
    double *abf = &temp[4*n1 + n5 + 4*n2]; // Nj*K3
    double *abfx = &temp[4*n1 + n5 + 4*n2 + n4]; // Nj*K3
    double *abfy = &temp[4*n1 + n5 + 4*n2 + 2*n4]; // Nj*K3
    double *abfz = &temp[4*n1 + n5 + 4*n2 + 3*n4]; // Nj*K3
    double *tm = &temp[4*n1 + n5 + 4*n2 + 4*n4]; // 4*K3

    angularbasis(abf, abfx, abfy, abfz, rij, tm, pq3, Nj, K3);

    radialangularbasis(U, Ux, Uy, Uz, rbf, rbfx, rbfy, rbfz, abf, abfx, abfy, abfz, Nj, K3, nrbf3);

    sumradialangularfunctions(sumU, U, tj, Nj, K3, nrbf3, nelements);

    double *d2 =  &temp[4*n1 + n5 + 4*n2]; // nl2
    double *dd2 = &temp[4*n1 + n5 + 4*n2 + nl2]; // 3*Nj*nl2
    double *d3 =  &temp[4*n1 + n5 + 4*n2 + nl2 + 3*Nj*nl2]; // nl3
    double *dd3 = &temp[4*n1 + n5 + 4*n2 + nl2 + 3*Nj*nl2 + nl3]; // 3*Nj*nl3
    double *d4 =  &temp[4*n1 + n5 + 4*n2 + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3]; // nl4
    double *dd4 = &temp[4*n1 + n5 + 4*n2 + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3 + nl4]; // 3*Nj*nl4

    if (nd23>0) {
      twobodydescderiv(d2, dd2, rbf, rbfx, rbfy, rbfz, tj, Nj);
    }

    if ((nd23>0) || (nd33>0) || (nd34>0)) {
      threebodydesc(d3, sumU, Nj);
      threebodydescderiv(dd3, sumU, Ux, Uy, Uz, tj, Nj);
    }

    double *cU = &temp[4*n1 + n5 + 4*n2 + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3 + nl4 + 3*Nj*nl4];
    e3 = threebodycoeff(cU, &coeff3[nl3*t0], sumU, Nj);
    tallylocalforce(fij, cU, Ux, Uy, Uz, tj, Nj, K3, nrbf3, nelements);

    if (nd23>0) {
      double *d23 = &temp[0];
      fourbodydesc23(d23, d2, d3);
      e23 = dotproduct(&coeff23[nl23*t0], d23, nl23);
      fourbodyfij23(fij, temp, &coeff23[nl23*t0], d2, d3, dd2, dd3, 3*Nj);
    }

    if (nd33>0) {
      double *d33 = &temp[0];
      fivebodydesc33(d33, d3);
      e33 = dotproduct(&coeff33[nl33*t0], d33, nl33);
      fivebodyfij33(fij, temp, &coeff33[nl33*t0], d3, dd3, 3*Nj);
    }

    if (nd4 > 0) {
      if (K4 < K3) {
        for (int m=0; m<nrbf4; m++)
          for (int k=0; k<K4; k++)
            for (int i=0; i<nelements; i++)
              sumU[i + nelements*k + nelements*K4*m] = sumU[i + nelements*k + nelements*K3*m];

        for (int m=0; m<nrbf4; m++)
          for (int k=0; k<K4; k++)
            for (int i=0; i<Nj; i++) {
              int ii = i + Nj*k + Nj*K4*m;
              int jj = i + Nj*k + Nj*K3*m;
              Ux[ii] = Ux[jj];
              Uy[ii] = Uy[jj];
              Uz[ii] = Uz[jj];
            }
      }

      if ((nd34>0) || (nd44>0)) {
        fourbodydescderiv(d4, dd4, sumU, Ux, Uy, Uz, tj, Nj);
      }

      e4 = fourbodycoeff(cU, sumU, &coeff4[nl4*t0], Nj);
      tallylocalforce(fij, cU, Ux, Uy, Uz, tj, Nj, K4, nrbf4, nelements);

      if (nd34>0) {
        double *d34 = &temp[0];
        sixbodydesc34(d34, d3, d4);
        e34 = dotproduct(&coeff34[nl34*t0], d34, nl34);
        sixbodyfij34(fij, temp, &coeff34[nl34*t0], d3, d4, dd3, dd4, 3*Nj);
      }

      if (nd44>0) {
        double *d44 = &temp[0];
        sevenbodydesc44(d44, d4);
        e44 = dotproduct(&coeff44[nl44*t0], d44, nl44);
        sevenbodyfij44(fij, temp, &coeff44[nl44*t0], d4, dd4, 3*Nj);
      }
    }
  }

  return (e1+e2+e3+e4+e23+e33+e34+e44);
}

double FASTPOD::energyforce(double *force, double *x, int *atomtype, int *alist,
          int *jlist, int *pairnumsum, int natom)
{
  double etot = 0.0;
  for (int i=0; i<3*natom; i++) force[i] = 0.0;

  for (int i=0; i<natom; i++) {
    int Nj = pairnumsum[i+1] - pairnumsum[i]; // # neighbors around atom i

    // reallocate temporary memory
    if (Nj>Njmax) {
      Njmax = Nj;
      int nmem = estimate_memory(Njmax);
      memory->destroy(tmpmem);
      memory->destroy(tmpint);
      memory->create(tmpmem, nmem, "tmpmem");
      memory->create(tmpint, 4*Njmax, "tmpint");
    }

    double *rij = &tmpmem[0];    // 3*Nj
    double *fij = &tmpmem[3*Nj]; // 3*Nj
    int *ai = &tmpint[0];        // Nj
    int *aj = &tmpint[Nj];       // Nj
    int *ti = &tmpint[2*Nj];     // Nj
    int *tj = &tmpint[3*Nj];     // Nj

    myneighbors(rij, x, ai, aj, ti, tj, jlist, pairnumsum, atomtype, alist, i);

    etot += peratomenergyforce(fij, rij, &tmpmem[6*Nj], ti, tj, Nj);

    tallyforce(force, fij, ai, aj, Nj);
  }

  return etot;
}

void FASTPOD::descriptors(double *gd, double *gdd, double *x, int *atomtype, int *alist,
          int *jlist, int *pairnumsum, int natom)
{
  for (int i=0; i<nd; i++) gd[i] = 0.0;
  for (int i=0; i<3*natom*nd; i++) gdd[i] = 0.0;

  double *gd1 = &gd[0];
  double *gd2 = &gd[nd1];
  double *gd3 = &gd[nd1+nd2];
  double *gd4 = &gd[nd1+nd2+nd3];
  double *gd23 = &gd[nd1+nd2+nd3+nd4];
  double *gd33 = &gd[nd1+nd2+nd3+nd4+nd23];
  double *gd34 = &gd[nd1+nd2+nd3+nd4+nd23+nd33];
  double *gd44 = &gd[nd1+nd2+nd3+nd4+nd23+nd33+nd34];

  int N = 3*natom;
  double *gdd1 = &gdd[0];
  double *gdd2 = &gdd[N*(nd1)];
  double *gdd3 = &gdd[N*(nd1+nd2)];
  double *gdd4 = &gdd[N*(nd1+nd2+nd3)];
  double *gdd23 = &gdd[N*(nd1+nd2+nd3+nd4)];
  double *gdd33 = &gdd[N*(nd1+nd2+nd3+nd4+nd23)];
  double *gdd34 = &gdd[N*(nd1+nd2+nd3+nd4+nd23+nd33)];
  double *gdd44 = &gdd[N*(nd1+nd2+nd3+nd4+nd23+nd33+nd34)];

  for (int i=0; i<natom; i++) {
    int Nj = pairnumsum[i+1] - pairnumsum[i]; // # neighbors around atom i

    // reallocate temporary memory
    if (Nj>Njmax) {
      Njmax = Nj;
      int nmem = estimate_memory(Njmax);
      memory->destroy(tmpmem);
      memory->destroy(tmpint);
      memory->create(tmpmem, nmem, "tmpmem");
      memory->create(tmpint, 4*Njmax, "tmpint");
    }

    double *rij = &tmpmem[0]; // 3*Nj
    int *ai = &tmpint[0];     // Nj
    int *aj = &tmpint[Nj];   // Nj
    int *ti = &tmpint[2*Nj]; // Nj
    int *tj = &tmpint[3*Nj]; // Nj

    myneighbors(rij, x, ai, aj, ti, tj, jlist, pairnumsum, atomtype, alist, i);

    if (nd1>0) onebodydescriptors(gd1, gdd1, ti, natom, i);

    double *d2= &tmpmem[3*Nj];
    double *dd2= &tmpmem[3*Nj + nl2];
    double *tmp = &tmpmem[3*Nj + nl2 + 3*Nj*nl2];

    for (int j=0; j<nl2; j++) d2[j] = 0.0;
    for (int j=0; j<3*Nj*nl2; j++) dd2[j] = 0.0;

    if (nd2>0) twobodydescriptors(gd2, gdd2, d2, dd2, rij, tmp, ai, aj, ti, tj, Nj, natom);

    double *d3= &tmpmem[3*Nj + nl2 + 3*Nj*nl2];
    double *dd3= &tmpmem[3*Nj + nl2 + 3*Nj*nl2 + nl3];
    tmp = &tmpmem[3*Nj + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3];

    for (int j=0; j<nl3; j++) d3[j] = 0.0;
    for (int j=0; j<3*Nj*nl3; j++) dd3[j] = 0.0;

    if (nd3>0) threebodydescriptors(gd3, gdd3, d3, dd3, rij, tmp, ai, aj, ti, tj, Nj, natom);

    double *d4= &tmpmem[3*Nj + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3];
    double *dd4= &tmpmem[3*Nj + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3 + nl4];
    tmp = &tmpmem[3*Nj + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3 + nl4 + 3*Nj*nl4];

    for (int j=0; j<nl4; j++) d4[j] = 0.0;
    for (int j=0; j<3*Nj*nl4; j++) dd4[j] = 0.0;

    if (nd4>0) fourbodydescriptors(gd4, gdd4, d4, dd4, rij, tmp, ai, aj, ti, tj, Nj, natom);

    int nld = 3*Nj + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3 + nl4 + 3*Nj*nl4;

    double *d23 = &tmpmem[nld];
    double *dd23 = &tmpmem[nld + nl23];
    if (nd23>0) fourbodydescriptors23(gd23, gdd23, d23, dd23, d2, d3, dd2, dd3,
            ai, aj, ti, tj, Nj, natom);

    double *d33 = &tmpmem[nld];
    double *dd33 = &tmpmem[nld + nl33];
    if (nd33>0) fivebodydescriptors33(gd33, gdd33, d33, dd33, d3, dd3,
            ai, aj, ti, tj, Nj, natom);

    double *d34 = &tmpmem[nld];
    double *dd34 = &tmpmem[nld + nl34];
    if (nd34>0) sixbodydescriptors34(gd34, gdd34, d34, dd34, d3, d4, dd3, dd4,
            ai, aj, ti, tj, Nj, natom);

    double *d44 = &tmpmem[nld];
    double *dd44 = &tmpmem[nld + nl44];
    if (nd44>0) sevenbodydescriptors44(gd44, gdd44, d44, dd44, d4, dd4,
            ai, aj, ti, tj, Nj, natom);
  }
}

void FASTPOD::fourbodydesc23(double *d23, double *d2, double *d3)
{
  for (int j = 0; j<n32; j++)
    for (int i = 0; i<n23; i++)
      d23[i + n23*j] = d2[ind23[i]]*d3[ind32[j]];
}

void FASTPOD::fourbodydescderiv23(double* dd23, double *d2, double *d3, double *dd2, double *dd3, int N)
{
  for (int j = 0; j<n32; j++)
    for (int i = 0; i<n23; i++)
      for (int n=0; n<N; n++)
        dd23[n + N*i + N*n23*j] = d2[ind23[i]]*dd3[n + N*ind32[j]] + dd2[n + N*ind23[i]]*d3[ind32[j]];
}

void FASTPOD::fourbodyfij23(double *fij, double *cf, double *coeff23, double *d2, double *d3,
        double *dd2, double *dd3, int N)
{
  for (int j = 0; j<n32; j++) {
    cf[j] = 0.0;
    for (int i = 0; i<n23; i++)
      cf[j] += d2[ind23[i]]*coeff23[i + n23*j];
  }

  for (int j = 0; j<n32; j++)
    for (int n=0; n<N; n++)
      fij[n] += cf[j]*dd3[n + N*ind32[j]];

  for (int i = 0; i<n23; i++) {
    cf[i] = 0.0;
    for (int j = 0; j<n32; j++)
      cf[i] += d3[ind32[j]]*coeff23[i + n23*j];
  }

  for (int i = 0; i<n23; i++)
    for (int n=0; n<N; n++)
      fij[n] += cf[i]*dd2[n + N*ind23[i]];
}

void FASTPOD::fivebodydesc33(double *d33, double *d3)
{
  int k = 0;
  for (int j = 0; j<n33; j++)
    for (int i = j; i<n33; i++) {
      d33[k] = d3[ind33[i]]*d3[ind33[j]];
      k += 1;
    }
}

void FASTPOD::fivebodydescderiv33(double *dd33, double *d3, double *dd3, int N)
{
  int k = 0;
  for (int j = 0; j<n33; j++)
    for (int i = j; i<n33; i++) {
      for (int n=0; n<N; n++)
        dd33[n + N*k] = d3[ind33[i]]*dd3[n + N*ind33[j]] + dd3[n + N*ind33[i]]*d3[ind33[j]];
      k += 1;
    }
}

void FASTPOD::fivebodyfij33(double *fij, double *cf, double *coeff33,
        double *d3, double *dd3, int N)
{
  int k = 0;
  for (int j = 0; j<n33; j++) {
    cf[j] = 0.0;
    for (int i = j; i<n33; i++) {
      cf[j] += d3[ind33[i]]*coeff33[k];
      k += 1;
    }
  }

  for (int j = 0; j<n33; j++)
    for (int n=0; n<N; n++)
      fij[n] += cf[j]*dd3[n + N*ind33[j]];

  for (int i = 0; i<n33; i++)cf[i] = 0.0;

  k = 0;
  for (int j = 0; j<n33; j++)
    for (int i = j; i<n33; i++) {
      cf[i] += d3[ind33[j]]*coeff33[k];
      k += 1;
    }

  for (int i=0; i<n33; i++)
    for (int n=0; n<N; n++)
      fij[n] += cf[i]*dd3[n + N*ind33[i]];

}

void FASTPOD::sixbodydesc34(double *d34, double *d3, double *d4)
{
  for (int j = 0; j<n43; j++)
    for (int i = 0; i<n34; i++)
      d34[i + n34*j] = d3[ind34[i]]*d4[ind43[j]];
}

void FASTPOD::sixbodydescderiv34(double *dd34, double *d3, double *d4, double *dd3, double *dd4, int N)
{
  for (int j = 0; j<n43; j++)
    for (int i = 0; i<n34; i++)
      for (int n=0; n<N; n++)
        dd34[n + N*i + N*n34*j] = d3[ind34[i]]*dd4[n + N*ind43[j]] + dd3[n + N*ind34[i]]*d4[ind43[j]];
}

void FASTPOD::sixbodyfij34(double *fij, double *cf, double *coeff34, double *d3, double *d4,
        double *dd3, double *dd4, int N)
{
  for (int j = 0; j<n43; j++) {
    cf[j] = 0.0;
    for (int i = 0; i<n34; i++)
      cf[j] += d3[ind34[i]]*coeff34[i + n34*j];
  }

  for (int j = 0; j<n43; j++)
    for (int n=0; n<N; n++)
      fij[n] += cf[j]*dd4[n + N*ind43[j]];

  for (int i = 0; i<n34; i++) {
    cf[i] = 0.0;
    for (int j = 0; j<n43; j++)
      cf[i] += d4[ind43[j]]*coeff34[i + n34*j];
  }

  for (int i = 0; i<n34; i++)
    for (int n=0; n<N; n++)
      fij[n] += cf[i]*dd3[n + N*ind34[i]];
}

void FASTPOD::sevenbodydesc44(double* d44, double *d4)
{
  int k = 0;
  for (int j = 0; j<n44; j++)
    for (int i = j; i<n44; i++) {
      d44[k] = d4[ind44[i]]*d4[ind44[j]];
      k += 1;
    }
}

void FASTPOD::sevenbodydescderiv44(double* dd44, double *d4, double *dd4, int N)
{
  int k = 0;
  for (int j = 0; j<n44; j++)
    for (int i = j; i<n44; i++) {
      for (int n=0; n<N; n++)
        dd44[n + N*k] = d4[ind44[i]]*dd4[n + N*ind44[j]] + dd4[n + N*ind44[i]]*d4[ind44[j]];
      k += 1;
    }
}

void FASTPOD::sevenbodyfij44(double *fij, double *cf, double *coeff44,
        double *d4, double *dd4, int N)
{
  int k = 0;
  for (int j = 0; j<n44; j++) {
    cf[j] = 0.0;
    for (int i = j; i<n44; i++) {
      cf[j] += d4[ind44[i]]*coeff44[k];
      k += 1;
    }
  }

  for (int j = 0; j<n44; j++)
    for (int n=0; n<N; n++)
      fij[n] += cf[j]*dd4[n + N*ind44[j]];

  for (int i = 0; i<n44; i++)cf[i] = 0.0;

  k = 0;
  for (int j = 0; j<n44; j++)
    for (int i = j; i<n44; i++) {
      cf[i] += d4[ind44[j]]*coeff44[k];
      k += 1;
    }

  for (int i=0; i<n44; i++)
    for (int n=0; n<N; n++)
      fij[n] += cf[i]*dd4[n + N*ind44[i]];
}

void FASTPOD::fourbodydescriptors23(double *gd23, double *gdd23, double *d23, double *dd23,
        double* d2, double *d3, double* dd2, double *dd3, int *ai, int *aj, int *ti, int *tj,
        int Nj, int natom)
{
  fourbodydesc23(d23, d2, d3);
  fourbodydescderiv23(dd23, d2, d3, dd2, dd3, 3*Nj);
  tallyglobdesc(gd23, d23, nl23, ti[0]-1);
  tallyglobdescderiv(gdd23, dd23,  ai, aj, natom, Nj, nl23, ti[0]-1);
}

void FASTPOD::fivebodydescriptors33(double *gd33, double *gdd33, double *d33, double *dd33,
        double *d3, double *dd3, int *ai, int *aj, int *ti, int *tj, int Nj, int natom)
{
  fivebodydesc33(d33, d3);
  fivebodydescderiv33(dd33, d3, dd3,  3*Nj);
  tallyglobdesc(gd33, d33, nl33, ti[0]-1);
  tallyglobdescderiv(gdd33, dd33,  ai, aj, natom, Nj, nl33, ti[0]-1);
}

void FASTPOD::sixbodydescriptors34(double *gd34, double *gdd34, double *d34, double *dd34,
        double* d3, double *d4, double* dd3, double *dd4, int *ai, int *aj, int *ti, int *tj,
        int Nj, int natom)
{
  sixbodydesc34(d34, d3, d4);
  sixbodydescderiv34(dd34, d3, d4, dd3, dd4, 3*Nj);
  tallyglobdesc(gd34, d34, nl34, ti[0]-1);
  tallyglobdescderiv(gdd34, dd34,  ai, aj, natom, Nj, nl34, ti[0]-1);
}

void FASTPOD::sevenbodydescriptors44(double *gd44, double *gdd44, double *d44, double *dd44,
        double *d4, double *dd4, int *ai, int *aj, int *ti, int *tj, int Nj, int natom)
{
  sevenbodydesc44(d44, d4);
  sevenbodydescderiv44(dd44, d4, dd4,  3*Nj);
  tallyglobdesc(gd44, d44, nl44, ti[0]-1);
  tallyglobdescderiv(gdd44, dd44,  ai, aj, natom, Nj, nl44, ti[0]-1);
}

void FASTPOD::myneighbors(double *rij, double *x, int *ai, int *aj, int *ti, int *tj,
        int *jlist, int *pairnumsum, int *atomtype, int *alist, int i)
{
  int itype = atomtype[i];
  int start = pairnumsum[i];
  int m = pairnumsum[i+1] - start; // number of neighbors around i
  for (int l=0; l<m ; l++) {   // loop over each atom around atom i
    int j = jlist[l + start];  // atom j
    ai[l]        = i;
    aj[l]        = alist[j];
    ti[l]        = itype;
    tj[l]        = atomtype[alist[j]];
    rij[0 + 3*l]   = x[0 + 3*j] -  x[0 + 3*i];
    rij[1 + 3*l]   = x[1 + 3*j] -  x[1 + 3*i];
    rij[2 + 3*l]   = x[2 + 3*j] -  x[2 + 3*i];
  }
}

void FASTPOD::fourbodydescderiv(double *d4, double *dd4, double *sumU, double *Ux, double *Uy,
        double *Uz, int *atomtype, int N)
{
  int Me = nelements*(nelements+1)*(nelements+2)/6; //count4(nelements);
  for (int m=0; m<nabf4*nrbf4*Me; m++)
    d4[m] = 0.0;

  for (int m=0; m<3*N*nabf4*nrbf4*Me; m++)
    dd4[m] = 0.0;

  int Q = pa4[nabf4];
  for (int m=0; m<nrbf4; m++)
    for (int p=0; p<nabf4; p++) {
      int n1 = pa4[p];
      int n2 = pa4[p+1];
      int nn = n2 - n1;
      for (int q=0; q<nn; q++) {
        int c = pc4[n1+q];
        int j1 = pb4[n1+q];
        int j2 = pb4[n1+q + Q];
        int j3 = pb4[n1+q + 2*Q];
        int k = 0;
        for (int i1=0; i1<nelements; i1++) {
          double c1 = c*sumU[i1 + nelements*j1 + nelements*K4*m];
          for (int i2=i1; i2<nelements; i2++) {
            double c2 = c*sumU[i2 + nelements*j2 + nelements*K4*m];
            double t12 = c1*sumU[i2 + nelements*j2 + nelements*K4*m];
            for (int i3=i2; i3<nelements; i3++) {
              double c3 = sumU[i3 + nelements*j3 + nelements*K4*m];
              double t13 = c1*c3;
              double t23 = c2*c3;
              int kk = p + nabf4*m + nabf4*nrbf4*k;
              int ii = 3*N*(p + nabf4*m + nabf4*nrbf4*k);
              d4[kk] += t12*c3;
              for (int j=0; j<N; j++) {
                int tj = atomtype[j]-1;
                if (tj==i3) {
                  int jj = j + N*j3 + N*K4*m;
                  dd4[0 + 3*j + ii] += t12*Ux[jj];
                  dd4[1 + 3*j + ii] += t12*Uy[jj];
                  dd4[2 + 3*j + ii] += t12*Uz[jj];
                }
                if (tj==i2) {
                  int jj = j + N*j2 + N*K4*m;
                  dd4[0 + 3*j + ii] += t13*Ux[jj];
                  dd4[1 + 3*j + ii] += t13*Uy[jj];
                  dd4[2 + 3*j + ii] += t13*Uz[jj];
                }
                if (tj==i1) {
                  int jj = j + N*j1 + N*K4*m;
                  dd4[0 + 3*j + ii] += t23*Ux[jj];
                  dd4[1 + 3*j + ii] += t23*Uy[jj];
                  dd4[2 + 3*j + ii] += t23*Uz[jj];
                }
              }
              k += 1;
            }
          }
        }
      }
    }
}

void FASTPOD::fourbodydescriptors(double *gd4, double *gdd4, double *d4, double *dd4, double *rij,
        double *temp, int *ai, int *aj, int *ti, int *tj, int Nj, int natom)
{
  int n1 = Nj*K4*nrbf4;
  int Me = nelements*(nelements+1)*(nelements+2)/6; //count4(nelements);
  int ndesc = nrbf4*nabf4*Me;

  double *U = &temp[0]; // Nj*K4*nrbf4
  double *Ux = &temp[n1]; // Nj*K4*nrbf4
  double *Uy = &temp[2*n1]; // Nj*K4*nrbf4
  double *Uz = &temp[3*n1]; // Nj*K4*nrbf4
  double *sumU = &temp[4*n1]; // K4*nrbf4*nelements

  unifiedbasis(U, Ux, Uy, Uz, sumU, rij, Phi, besselparams, temp, rin, rcut, pdegree,
        tj, pq4, nbesselpars, nrbf4, K4, nelements, Nj);

  fourbodydescderiv(d4, dd4, sumU, Ux, Uy, Uz, tj, Nj);

  tallyglobdesc(gd4, d4, ndesc, ti[0]-1);

  tallyglobdescderiv(gdd4, dd4,  ai, aj, natom, Nj, ndesc, ti[0]-1);
}

void FASTPOD::fourbodydescriptors(double *d4, double *dd4, double *rij, double *temp, int *tj, int Nj)
{
  int n1 = Nj*K4*nrbf4;

  double *U = &temp[0]; // Nj*K4*nrbf4
  double *Ux = &temp[n1]; // Nj*K4*nrbf4
  double *Uy = &temp[2*n1]; // Nj*K4*nrbf4
  double *Uz = &temp[3*n1]; // Nj*K4*nrbf4
  double *sumU = &temp[4*n1]; // K4*nrbf4*nelements

  unifiedbasis(U, Ux, Uy, Uz, sumU, rij, Phi, besselparams, temp, rin, rcut, pdegree,
        tj, pq4, nbesselpars, nrbf4, K4, nelements, Nj);

  fourbodydescderiv(d4, dd4, sumU, Ux, Uy, Uz, tj, Nj);
}

double FASTPOD::fourbodycoeff(double *cU, double *sumU, double *coeff4, int N)
{
  for (int m=0; m<nelements*K4*nrbf4; m++)
    cU[m] = 0.0;

  int Q = pa4[nabf4];
  double e = 0.0;

  for (int m=0; m<nrbf4; m++)
    for (int p=0; p<nabf4; p++) {
      int n1 = pa4[p];
      int n2 = pa4[p+1];
      int nn = n2 - n1;
      for (int q=0; q<nn; q++) {
        int c = pc4[n1+q];
        int j1 = pb4[n1+q];
        int j2 = pb4[n1+q + Q];
        int j3 = pb4[n1+q + 2*Q];
        int k = 0;
        for (int i1=0; i1<nelements; i1++) {
          double c1 = c*sumU[i1 + nelements*j1 + nelements*K4*m];
          for (int i2=i1; i2<nelements; i2++) {
            double c0 = sumU[i2 + nelements*j2 + nelements*K4*m];
            double c2 = c*c0;
            double c4 = c1*c0;
            for (int i3=i2; i3<nelements; i3++) {
              double c5 = coeff4[p + nabf4*m + nabf4*nrbf4*k];
              double c6 = c5*sumU[i3 + nelements*j3 + nelements*K4*m];
              e += c4*c6;
              cU[i3 + nelements*j3 + nelements*K4*m] += c5*c4;
              cU[i2 + nelements*j2 + nelements*K4*m] += c6*c1;
              cU[i1 + nelements*j1 + nelements*K4*m] += c6*c2;
              k += 1;
            }
          }
        }
      }
    }

  return e;
}

void FASTPOD::threebodydesc(double *d3, double *sumU, int N)
{
  int Me = nelements*(nelements+1)/2;
  for (int m=0; m<nabf3*nrbf3*Me; m++)
    d3[m] = 0.0;

  for (int m=0; m<nrbf3; m++)
    for (int p=0; p<nabf3; p++) {
      int n1 = pn3[p];
      int n2 = pn3[p+1];
      int nn = n2 - n1;
      for (int q=0; q<nn; q++) {
        int k = 0;
        for (int i1=0; i1<nelements; i1++) {
          double t1 = pc3[n1+q]*sumU[i1 + nelements*(n1+q) + nelements*K3*m];
          for (int i2=i1; i2<nelements; i2++) {
            d3[p + nabf3*m + nabf3*nrbf3*k] += t1*sumU[i2 + nelements*(n1+q) + nelements*K3*m];
            k += 1;
          }
        }
      }
    }
}

void FASTPOD::threebodydescderiv(double *dd3, double *sumU, double *Ux, double *Uy, double *Uz,
        int *atomtype, int N)
{
  int Me = nelements*(nelements+1)/2;
  for (int m=0; m<3*N*nabf3*nrbf3*Me; m++)
    dd3[m] = 0.0;

  for (int m=0; m<nrbf3; m++)
    for (int p=0; p<nabf3; p++) {
      int n1 = pn3[p];
      int n2 = pn3[p+1];
      int nn = n2 - n1;
      for (int q=0; q<nn; q++) {
        for (int i1=0; i1<nelements; i1++) {
          double t1 = pc3[n1+q]*sumU[i1 + nelements*(n1+q) + nelements*K3*m];
          for (int j=0; j<N; j++) {
            int i2 = atomtype[j]-1;
            int k = elemindex[i2 + nelements*i1];
            double f = (i1==i2) ? 2.0*t1 : t1;
            int ii = 3*j + 3*N*(p + nabf3*m + nabf3*nrbf3*k);
            int jj = j + N*(n1+q) + N*K3*m;
            dd3[0 + ii] += f*Ux[jj];
            dd3[1 + ii] += f*Uy[jj];
            dd3[2 + ii] += f*Uz[jj];
          }
        }
      }
    }
}

void FASTPOD::threebodydescriptors(double *d3, double *dd3, double *rij, double *temp, int *tj, int Nj)
{
  int n1 = Nj*K3*nrbf3;

  double *U = &temp[0]; // Nj*K*nrbf
  double *Ux = &temp[n1]; // Nj*K*nrbf
  double *Uy = &temp[2*n1]; // Nj*K*nrbf
  double *Uz = &temp[3*n1]; // Nj*K*nrbf
  double *sumU = &temp[4*n1]; // K*nrbf*nelements

  unifiedbasis(U, Ux, Uy, Uz, sumU, rij, Phi, besselparams, temp, rin, rcut, pdegree,
        tj, pq3, nbesselpars, nrbf3, K3, nelements, Nj);

  threebodydesc(d3, sumU, Nj);

  threebodydescderiv(dd3, sumU, Ux, Uy, Uz, tj, Nj);
}

void FASTPOD::threebodydescriptors(double *gd3, double *gdd3, double *d3, double *dd3, double *rij,
        double *temp, int *ai, int *aj, int *ti, int *tj, int Nj, int natom)
{
  int ndesc = nrbf3*nabf3*nelements*(1+nelements)/2;
  int n1 = Nj*K3*nrbf3;
  //int n5 = K3*nrbf3*nelements;

  double *U = &temp[0]; // Nj*K3*nrbf3
  double *Ux = &temp[n1]; // Nj*K3*nrbf3
  double *Uy = &temp[2*n1]; // Nj*K3*nrbf3
  double *Uz = &temp[3*n1]; // Nj*K3*nrbf3
  double *sumU = &temp[4*n1]; // K3*nrbf3*nelements

  unifiedbasis(U, Ux, Uy, Uz, sumU, rij, Phi, besselparams, temp, rin, rcut, pdegree,
        tj, pq3, nbesselpars, nrbf3, K3, nelements, Nj);

  threebodydesc(d3, sumU, Nj);

  threebodydescderiv(dd3, sumU, Ux, Uy, Uz, tj, Nj);

  tallyglobdesc(gd3, d3, ndesc, ti[0]-1);

  tallyglobdescderiv(gdd3, dd3,  ai, aj, natom, Nj, ndesc, ti[0]-1);
}

double FASTPOD::threebodycoeff(double *cU, double *coeff3, double *sumU, int N)
{
  for (int m=0; m<nrbf3*K3*nelements; m++)
    cU[m] = 0.0;

  double e=0.0;
  for (int m=0; m<nrbf3; m++)
    for (int p=0; p<nabf3; p++) {
      int n1 = pn3[p];
      int n2 = pn3[p+1];
      int nn = n2 - n1;
      for (int q=0; q<nn; q++) {
        int k = 0;
        for (int i1=0; i1<nelements; i1++) {
          double t1 = pc3[n1+q]*sumU[i1 + nelements*(n1+q) + nelements*K3*m];
          for (int i2=i1; i2<nelements; i2++) {
            double c2 = sumU[i2 + nelements*(n1+q) + nelements*K3*m];
            double c3 = coeff3[p + nabf3*m + nabf3*nrbf3*k];
            double t2 = c3*t1;
            e += t2*c2;
            cU[i2 + nelements*(n1+q) + nelements*K3*m] += t2;
            cU[i1 + nelements*(n1+q) + nelements*K3*m] += pc3[n1+q]*c2*c3;
            k += 1;
          }
        }
      }
    }
  return e;
}

void FASTPOD::twobodydescderiv(double *d2, double *dd2, double *rbf, double *rbfx,
        double *rbfy, double *rbfz, int *tj, int N)
{
  for (int m=0; m<nl2; m++)
    d2[m] = 0.0;
  for (int m=0; m<3*N*nl2; m++)
    dd2[m] = 0.0;

  for (int m=0; m<nrbf2; m++) {
    for (int n=0; n<N; n++) {
      int i2 = n + N*m;
      int i1 = n + N*m + N*nrbf2*(tj[n]-1);
      d2[m + nrbf2*(tj[n]-1)] += rbf[i2];
      dd2[0 + 3*i1] += rbfx[i2];
      dd2[1 + 3*i1] += rbfy[i2];
      dd2[2 + 3*i1] += rbfz[i2];
    }
  }
}

void FASTPOD::twobodydescriptors(double *d2, double *dd2, double *rij, double *temp, int *tj, int Nj)
{
  orthogonalradialbasis(temp, rij, Phi, besselparams, rin, rcut-rin,
          pdegree[0], pdegree[1], nbesselpars, nrbf2, Nj);

  int n2 = Nj*nrbf2;
  twobodydescderiv(d2, dd2, temp, &temp[n2], &temp[2*n2], &temp[3*n2], tj, Nj);
}

void FASTPOD::twobodydescriptors(double *gd2, double *gdd2, double *d2, double *dd2, double *rij,
        double *temp, int *ai, int *aj, int *ti, int *tj, int Nj, int natom)
{
  orthogonalradialbasis(temp, rij, Phi, besselparams, rin, rcut-rin,
          pdegree[0], pdegree[1], nbesselpars, nrbf2, Nj);

  int n2 = Nj*nrbf2;
  twobodydescderiv(d2, dd2, temp, &temp[n2], &temp[2*n2], &temp[3*n2], tj, Nj);

  tallytwobodyglobdesc(gd2, d2, elemindex, nrbf2, nelements, ti[0]-1);

  tallytwobodyglobdescderiv(gdd2, dd2, ai, aj, ti, tj, elemindex, nrbf2, nelements, natom, Nj);
}


void FASTPOD::onebodydescriptors(double *gd1, double *gdd1, int *ti, int natom, int i)
{
  for (int m=0; m<nelements; m++) {
    double onehot = (ti[0] == (m+1)) ? 1.0 : 0.0;
    gd1[m] += onehot;
    gdd1[0 + 3*(i + natom*m)] = 0.0;
    gdd1[1 + 3*(i + natom*m)] = 0.0;
    gdd1[2 + 3*(i + natom*m)] = 0.0;
  }
}

double FASTPOD::onebodyenergy(double *coeff1, int *ti)
{
  double e1 = 0.0;
  for (int m=0; m<nelements; m++)
    e1 += (ti[0] == (m+1)) ? coeff1[m] : 0.0;

  return e1;
}

void FASTPOD::radialfunctions(double *rbf, double *rij, double *besselparams, double rin,
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N)
{
  for (int n=0; n<N; n++) {
    double rij1 = rij[0+3*n];
    double rij2 = rij[1+3*n];
    double rij3 = rij[2+3*n];

    double dij = pow(rij1*rij1 + rij2*rij2 + rij3*rij3, 0.5);
    double r = dij - rin;
    double y = r/rmax;
    double y2 = y*y;
    double y3 = 1.0 - y2*y;
    double y4 = y3*y3 + 1e-6;
    double y5 = pow(y4, 0.5);
    double y6 = exp(-1.0/y5);
    double fcut = y6/exp(-1.0);

    for (int j=0; j<nbesselpars; j++) {
      double x =  (1.0 - exp(-besselparams[j]*r/rmax))/(1.0-exp(-besselparams[j]));
      for (int i=0; i<besseldegree; i++)
        rbf[n + N*i + N*besseldegree*j] = ((sqrt(2.0/(rmax))/(i+1)))*fcut*sin((i+1)*M_PI*x)/r;
    }

    for (int i=0; i<inversedegree; i++) {
      int p = besseldegree*nbesselpars + i;
      double a = pow(dij, (double) (i+1.0));
      rbf[n + N*p] = fcut/a;
    }
  }
}

void FASTPOD::radialbasis(double *rbf, double *rbfx, double *rbfy, double *rbfz, double *rij, double *besselparams, double rin,
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N)
{
  for (int n=0; n<N; n++) {
    double rij1 = rij[0+3*n];
    double rij2 = rij[1+3*n];
    double rij3 = rij[2+3*n];

    double dij = pow(rij1*rij1 + rij2*rij2 + rij3*rij3, 0.5);
    double dr1 = rij1/dij;
    double dr2 = rij2/dij;
    double dr3 = rij3/dij;

    double r = dij - rin;
    double y = r/rmax;
    double y2 = y*y;
    double y3 = 1.0 - y2*y;
    double y4 = y3*y3 + 1e-6;
    double y5 = pow(y4, 0.5);
    double y6 = exp(-1.0/y5);
    double y7 = pow(y4, 1.5);
    double fcut = y6/exp(-1.0);
    double dfcut = ((3.0/(rmax*exp(-1.0)))*(y2)*y6*(y*y2 - 1.0))/y7;

    for (int j=0; j<nbesselpars; j++) {
      double alpha = besselparams[j];
      if (fabs(alpha) <= 1.0e-6) alpha = 1e-3;
      double x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));
      double dx = (alpha/rmax)*exp(-(alpha*r/rmax))/(1.0 - exp(-alpha));

      for (int i=0; i<besseldegree; i++) {
        double a = (i+1)*M_PI;
        double b = (sqrt(2.0/(rmax))/(i+1));
        int nij = n + N*i + N*besseldegree*j;
        rbf[nij] = b*fcut*sin(a*x)/r;
        double drbfdr = b*(dfcut*sin(a*x)/r - fcut*sin(a*x)/(r*r) + a*cos(a*x)*fcut*dx/r);
        rbfx[nij] = drbfdr*dr1;
        rbfy[nij] = drbfdr*dr2;
        rbfz[nij] = drbfdr*dr3;
      }
    }

    for (int i=0; i<inversedegree; i++) {
      int p = besseldegree*nbesselpars + i;
      int nij = n + N*p;
      double a = pow(dij, (double) (i+1.0));
      rbf[nij] = fcut/a;
      double drbfdr = dfcut/a - (i+1.0)*fcut/(a*dij);
      rbfx[nij] = drbfdr*dr1;
      rbfy[nij] = drbfdr*dr2;
      rbfz[nij] = drbfdr*dr3;
    }
  }
}

void FASTPOD::orthogonalradialbasis(double *orthorbf, double *rij, double *Phi, double *besselparams,
        double rin, double rmax, int besseldegree, int inversedegree, int nbesselpars, int nrbf2, int N)
{
  int ns = besseldegree*nbesselpars + inversedegree;
  int n2 = N*nrbf2;
  int n3 = N*ns;

  double *rbf = &orthorbf[0]; // Nj*nrbf2
  double *rbfx = &orthorbf[n2]; // Nj*nrbf2
  double *rbfy = &orthorbf[2*n2]; // Nj*nrbf2
  double *rbfz = &orthorbf[3*n2]; // Nj*nrbf2

  double *rbft = &orthorbf[4*n2]; // Nj*ns
  double *rbfxt = &orthorbf[4*n2 + n3]; // Nj*ns
  double *rbfyt = &orthorbf[4*n2 + 2*n3]; // Nj*ns
  double *rbfzt = &orthorbf[4*n2 + 3*n3]; // Nj*ns

  // orthogonal radial basis functions
  radialbasis(rbft, rbfxt, rbfyt, rbfzt, rij, besselparams, rin, rmax,
          besseldegree, inversedegree, nbesselpars, N);
  MatMul(rbf, rbft, Phi, N, ns, nrbf2);
  MatMul(rbfx, rbfxt, Phi, N, ns, nrbf2);
  MatMul(rbfy, rbfyt, Phi, N, ns, nrbf2);
  MatMul(rbfz, rbfzt, Phi, N, ns, nrbf2);
}

void FASTPOD::angularfunctions(double *abf, double *rij, double *tm, int *pq, int N, int K)
{
  tm[0] = 1.0;
  for (int j=0; j<N; j++) {
    double x = rij[0+3*j];
    double y = rij[1+3*j];
    double z = rij[2+3*j];

    double xx = x*x;
    double yy = y*y;
    double zz = z*z;

    double dij = sqrt(xx + yy + zz);
    double u = x/dij;
    double v = y/dij;
    double w = z/dij;

    abf[j] = tm[0];
    for (int n=1; n<K; n++) {
      int m = pq[n]-1;
      int d = pq[n + K];
      if (d==1)
        tm[n] = tm[m]*u;
      else if (d==2)
        tm[n] = tm[m]*v;
      else if (d==3)
        tm[n] = tm[m]*w;
      abf[j + N*n] = tm[n];
    }
  }
}

void FASTPOD::angularbasis(double *abf, double *abfx, double *abfy, double *abfz, double *rij, double *tm, int *pq, int N, int K)
{
  double *tmu = &tm[K];
  double *tmv = &tm[2*K];
  double *tmw = &tm[3*K];

  tm[0] = 1.0;
  tmu[0] = 0.0;
  tmv[0] = 0.0;
  tmw[0] = 0.0;
  for (int j=0; j<N; j++) {
    double x = rij[0+3*j];
    double y = rij[1+3*j];
    double z = rij[2+3*j];

    double xx = x*x;
    double yy = y*y;
    double zz = z*z;
    double xy = x*y;
    double xz = x*z;
    double yz = y*z;

    double dij = sqrt(xx + yy + zz);
    double u = x/dij;
    double v = y/dij;
    double w = z/dij;

    double dij3 = dij*dij*dij;
    double dudx = (yy+zz)/dij3;
    double dudy = -xy/dij3;
    double dudz = -xz/dij3;

    double dvdx = -xy/dij3;
    double dvdy = (xx+zz)/dij3;
    double dvdz = -yz/dij3;

    double dwdx = -xz/dij3;
    double dwdy = -yz/dij3;
    double dwdz = (xx+yy)/dij3;

    abf[j] = tm[0];
    abfx[j] = 0.0;
    abfy[j] = 0.0;
    abfz[j] = 0.0;
    for (int n=1; n<K; n++) {
      int m = pq[n]-1;
      int d = pq[n + K];
      if (d==1) {
        tm[n] = tm[m]*u;
        tmu[n] = tmu[m]*u + tm[m];
        tmv[n] = tmv[m]*u;
        tmw[n] = tmw[m]*u;
      }
      else if (d==2) {
        tm[n] = tm[m]*v;
        tmu[n] = tmu[m]*v;
        tmv[n] = tmv[m]*v + tm[m];
        tmw[n] = tmw[m]*v;
      }
      else if (d==3) {
        tm[n] = tm[m]*w;
        tmu[n] = tmu[m]*w;
        tmv[n] = tmv[m]*w;
        tmw[n] = tmw[m]*w + tm[m];
      }
      abf[j + N*n] = tm[n];
      abfx[j + N*n] = tmu[n]*dudx + tmv[n]*dvdx + tmw[n]*dwdx;
      abfy[j + N*n] = tmu[n]*dudy + tmv[n]*dvdy + tmw[n]*dwdy;
      abfz[j + N*n] = tmu[n]*dudz + tmv[n]*dvdz + tmw[n]*dwdz;
    }
  }
}

void FASTPOD::radialangularfunctions(double *U, double *rbf, double *abf, int N, int K, int M)
{
  for (int m=0; m<M; m++)
    for (int k=0; k<K; k++)
      for (int n=0; n<N; n++)
        U[n + N*k + N*K*m] = abf[n + N*k]*rbf[n + N*m];
}

void FASTPOD::radialangularbasis(double *U, double *Ux, double *Uy, double *Uz,
        double *rbf, double *rbfx, double *rbfy, double *rbfz, double *abf,
        double *abfx, double *abfy, double *abfz, int N, int K, int M)
{
  for (int m=0; m<M; m++)
    for (int k=0; k<K; k++)
      for (int n=0; n<N; n++) {
        int ia = n + N*k;
        int ib = n + N*m;
        int ii = ia + N*K*m;
        double c1 = rbf[ib];
        double c2 = abf[ia];
        U[ii] = c1*c2;
        Ux[ii] = abfx[ia]*c1 + c2*rbfx[ib];
        Uy[ii] = abfy[ia]*c1 + c2*rbfy[ib];
        Uz[ii] = abfz[ia]*c1 + c2*rbfz[ib];
      }
}

void FASTPOD::sumradialangularfunctions(double *sumU, double *U, int *atomtype, int N, int K, int M, int Ne)
{
  for (int m=0; m<Ne*K*M; m++)
    sumU[m] = 0.0;

  for (int m=0; m<M; m++)
    for (int k=0; k<K; k++)
      for (int n=0; n<N; n++)
        sumU[atomtype[n]-1 + Ne*k + Ne*K*m] += U[n + N*k + N*K*m];

}

void FASTPOD::unifiedbasis(double *U, double *Ux, double *Uy, double *Uz, double *sumU, double *rij,
        double *Phi, double *besselparams, double *tmpmem, double rin, double rcut, int *pdegree,
        int *tj, int *pq, int nbesselpars, int nrbf, int K, int nelements, int Nj)
{
  int ns = pdegree[0]*nbesselpars + pdegree[1];
  int n1 = Nj*K*nrbf;
  int n2 = Nj*nrbf;
  int n3 = Nj*ns;
  int n4 = Nj*K;

  double *rbf = &tmpmem[4*n1]; // Nj*nrbf
  double *rbfx = &tmpmem[4*n1 + n2]; // Nj*nrbf
  double *rbfy = &tmpmem[4*n1 + 2*n2]; // Nj*nrbf
  double *rbfz = &tmpmem[4*n1 + 3*n2]; // Nj*nrbf

  double *rbft = &tmpmem[4*n1 + 4*n2]; // Nj*ns
  double *rbfxt = &tmpmem[4*n1 + 4*n2 + n3]; // Nj*ns
  double *rbfyt = &tmpmem[4*n1 + 4*n2 + 2*n3]; // Nj*ns
  double *rbfzt = &tmpmem[4*n1 + 4*n2 + 3*n3]; // Nj*ns

  // orthogonal radial basis functions
  radialbasis(rbft, rbfxt, rbfyt, rbfzt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nj);
  MatMul(rbf, rbft, Phi, Nj, ns, nrbf);
  MatMul(rbfx, rbfxt, Phi, Nj, ns, nrbf);
  MatMul(rbfy, rbfyt, Phi, Nj, ns, nrbf);
  MatMul(rbfz, rbfzt, Phi, Nj, ns, nrbf);

  double *abf = &tmpmem[4*n1 + 4*n2]; // Nj*K
  double *abfx = &tmpmem[4*n1 + 4*n2 + n4]; // Nj*K
  double *abfy = &tmpmem[4*n1 + 4*n2 + 2*n4]; // Nj*K
  double *abfz = &tmpmem[4*n1 + 4*n2 + 3*n4]; // Nj*K
  double *tm = &tmpmem[4*n1 + 4*n2 + 4*n4]; // 4*K

  angularbasis(abf, abfx, abfy, abfz, rij, tm, pq, Nj, K);

  radialangularbasis(U, Ux, Uy, Uz, rbf, rbfx, rbfy, rbfz, abf, abfx, abfy, abfz, Nj, K, nrbf);

  sumradialangularfunctions(sumU, U, tj, Nj, K, nrbf, nelements);
}

void FASTPOD::tallytwobodyglobdesc(double *gd, double *d, int *elemindex, int nrbf, int nelements, int ti)
{
  for (int m = 0; m<nelements; m++) {
    int i1 = nrbf*elemindex[m + nelements*ti];
    int i2 = nrbf*m;
    for (int k = 0; k<nrbf; k++)
      gd[k + i1] += d[k + i2];
  }
}

void FASTPOD::tallytwobodyglobdescderiv(double *gdd, double *dd, int *ai, int *aj, int *ti, int *tj,
        int *elemindex, int nrbf, int nelements, int natom, int N)
{
    for (int k = 0; k<nrbf; k++)
      for (int n=0; n<N; n++) {
        int i1 = nrbf*elemindex[(tj[n]-1) + nelements*(ti[n]-1)];
        int i2 = nrbf*(tj[n]-1);
        int im =  3*ai[n] + 3*natom*(k+i1);
        int jm =  3*aj[n] + 3*natom*(k+i1);
        int nm = 3*n + 3*N*(k+i2);
        gdd[0 + im] += dd[0 + nm];
        gdd[1 + im] += dd[1 + nm];
        gdd[2 + im] += dd[2 + nm];
        gdd[0 + jm] -= dd[0 + nm];
        gdd[1 + jm] -= dd[1 + nm];
        gdd[2 + jm] -= dd[2 + nm];
      }
}

void FASTPOD::tallyglobdesc(double *gd, double *d, int ndesc, int ti)
{
  int k = ndesc*ti;
  for (int n = 0; n<ndesc; n++)
    gd[n + k] += d[n];
}

void FASTPOD::tallyglobdescderiv(double *gdd, double *dd,  int *ai, int *aj, int natom, int N, int ndesc, int ti)
{
  int k = 3*natom*ndesc*ti;
  for (int m=0; m<ndesc; m++)
    for (int n=0; n<N; n++)
    {
      int im =  3*ai[n] + 3*natom*m + k;
      int jm =  3*aj[n] + 3*natom*m + k;
      int nm = 3*n + 3*N*m;
      gdd[0 + im] += dd[0 + nm];
      gdd[1 + im] += dd[1 + nm];
      gdd[2 + im] += dd[2 + nm];
      gdd[0 + jm] -= dd[0 + nm];
      gdd[1 + jm] -= dd[1 + nm];
      gdd[2 + jm] -= dd[2 + nm];
    }
}

double FASTPOD::tallytwobodylocalforce(double *fij, double *coeff2,  double *rbf, double *rbfx,
        double *rbfy, double *rbfz, int *tj, int nbf, int N)
{
  double e = 0.0;
  for (int m=0; m<nbf; m++)
    for (int n=0; n<N; n++) {
      int nm = n + N*m;
      double c = coeff2[m + nbf*(tj[n]-1)];
      e += c*rbf[nm];
      fij[0 + 3*n] += c*rbfx[nm];
      fij[1 + 3*n] += c*rbfy[nm];
      fij[2 + 3*n] += c*rbfz[nm];
    }
  return e;
}

void FASTPOD::tallylocalforce(double *fij, double *cU, double *Ux, double *Uy, double *Uz,
        int *atomtype, int N, int K, int M, int Ne)
{
  for (int m=0; m<M; m++)
    for (int k=0; k<K; k++)
      for (int j=0; j<N; j++) {
        int i2 = atomtype[j]-1;
        int ii = 3*j;
        int jj = j + N*k + N*K*m;
        double c = cU[i2 + Ne*k + Ne*K*m];
        fij[0 + ii] += c*Ux[jj];
        fij[1 + ii] += c*Uy[jj];
        fij[2 + ii] += c*Uz[jj];
      }
}

void FASTPOD::tallyforce(double *force, double *fij,  int *ai, int *aj, int N)
{
  for (int n=0; n<N; n++) {
    int im =  3*ai[n];
    int jm =  3*aj[n];
    int nm = 3*n;
    force[0 + im] += fij[0 + nm];
    force[1 + im] += fij[1 + nm];
    force[2 + im] += fij[2 + nm];
    force[0 + jm] -= fij[0 + nm];
    force[1 + jm] -= fij[1 + nm];
    force[2 + jm] -= fij[2 + nm];
  }
}

void FASTPOD::tallyforce(double **force, double *fij,  int *ai, int *aj, int N)
{
  for (int n=0; n<N; n++) {
    int im =  ai[n];
    int jm =  aj[n];
    int nm = 3*n;
    force[im][0] += fij[0 + nm];
    force[im][1] += fij[1 + nm];
    force[im][2] += fij[2 + nm];
    force[jm][0] -= fij[0 + nm];
    force[jm][1] -= fij[1 + nm];
    force[jm][2] -= fij[2 + nm];
  }
}

void FASTPOD::twobodycoeff(double *newcoeff2, double *coeff2)
{
  for (int i=0; i<nelements; i++)
    for (int j=0; j<nelements; j++)
      for (int m=0; m<nrbf2; m++)
        newcoeff2[m + nrbf2*j + nrbf2*nelements*i] = coeff2[m + nrbf2*elemindex[i + nelements*j]];
}

void FASTPOD::mknewcoeff()
{
  memory->create(newcoeff, nl*nelements, "newcoeff");

  for (int n=0; n<nd1; n++)
    newcoeff[n] = coeff[n];

  if (nd2>0) twobodycoeff(&newcoeff[nd1], &coeff[nd1]);

  for (int n=0; n<(nd3+nd4+nd23+nd33+nd34+nd44); n++)
    newcoeff[(nl1 + nl2)*nelements + n] = coeff[nd1+nd2+n];
}

void FASTPOD::mknewcoeff(double *c)
{
  memory->create(newcoeff, nl*nelements, "newcoeff");

  for (int n=0; n<nd1; n++)
    newcoeff[n] = c[n];

  if (nd2>0) twobodycoeff(&newcoeff[nd1], &c[nd1]);

  for (int n=0; n<(nd3+nd4+nd23+nd33+nd34+nd44); n++)
    newcoeff[(nl1 + nl2)*nelements + n] = c[nd1+nd2+n];
}

void FASTPOD::snapshots(double *rbf, double *xij, int N)
{
  double rmax = rcut-rin;
  for (int n=0; n<N; n++) {
    double dij = xij[n];

    double r = dij - rin;
    double y = r/rmax;
    double y2 = y*y;
    double y3 = 1.0 - y2*y;
    double y4 = y3*y3 + 1e-6;
    double y5 = pow(y4, 0.5);
    double y6 = exp(-1.0/y5);
    double fcut = y6/exp(-1.0);

    for (int j=0; j<nbesselpars; j++) {
      double alpha = besselparams[j];
      if (fabs(alpha) <= 1.0e-6) alpha = 1e-3;
      double x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));

      for (int i=0; i<besseldegree; i++) {
        double a = (i+1)*M_PI;
        double b = (sqrt(2.0/(rmax))/(i+1));
        int nij = n + N*i + N*besseldegree*j;
        rbf[nij] = b*fcut*sin(a*x)/r;
      }
    }

    for (int i=0; i<inversedegree; i++) {
      int p = besseldegree*nbesselpars + i;
      int nij = n + N*p;
      double a = pow(dij, (double) (i+1.0));
      rbf[nij] = fcut/a;
    }
  }
}

void FASTPOD::eigenvaluedecomposition(double *Phi, double *Lambda, int N)
{
  int ns = besseldegree*nbesselpars + inversedegree;

  double *xij = (double *) malloc(N*sizeof(double));
  double *S = (double *) malloc(N*ns*sizeof(double));
  double *Q = (double *) malloc(N*ns*sizeof(double));
  double *A = (double *) malloc(ns*ns*sizeof(double));
  double *b = (double *) malloc(ns*sizeof(double));

  for (int i=0; i<N; i++)
    xij[i] = (rin+1e-6) + (rcut-rin-1e-6)*(i*1.0/(N-1));

  snapshots(S, xij, N);

  char chn = 'N';
  char cht = 'T';
  char chv = 'V';
  char chu = 'U';
  double alpha = 1.0, beta = 0.0;
  DGEMM(&cht, &chn, &ns, &ns, &N, &alpha, S, &N, S, &N, &beta, A, &ns);

  for (int i=0; i<ns*ns; i++)
    A[i] = A[i]*(1.0/N);

  int lwork = ns * ns;  // the length of the array work, lwork >= max(1,3*N-1)
  int info = 1;     // = 0:  successful exit
  double work[ns*ns];
  DSYEV(&chv, &chu, &ns, A, &ns, b, work, &lwork, &info);

  // order eigenvalues and eigenvectors from largest to smallest

  for (int j=0; j<ns; j++)
    for (int i=0; i<ns; i++)
      Phi[i + ns*(ns-j-1)] = A[i + ns*j];

  for (int i=0; i<ns; i++)
    Lambda[(ns-i-1)] = b[i];

  DGEMM(&chn, &chn, &N, &ns, &ns, &alpha, S, &N, Phi, &ns, &beta, Q, &N);
  for (int i=0; i<(N-1); i++)
    xij[i] = xij[i+1] - xij[i];
  double area;
  for (int m=0; m<ns; m++) {
    area = 0.0;
    for (int i=0; i<(N-1); i++)
      area += 0.5*xij[i]*(Q[i + N*m]*Q[i + N*m] + Q[i+1 + N*m]*Q[i+1 + N*m]);
    for (int i=0; i<ns; i++)
      Phi[i + ns*m] = Phi[i + ns*m]/sqrt(area);
  }

  free(xij); free(S); free(A); free(b); free(Q);
}

void FASTPOD::init2body()
{
  pdegree[0] = besseldegree;
  pdegree[1] = inversedegree;
  ns = nbesselpars*pdegree[0] + pdegree[1];

  memory->create(Phi, ns*ns, "Phi");
  memory->create(Lambda, ns, "Lambda");

  eigenvaluedecomposition(Phi, Lambda, 2000);
}

void FASTPOD::init3body(int Pa3)
{
  int npa[] = {0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455};

  nabf3 = Pa3+1;
  K3 = npa[nabf3];
  P3 = nabf3-1;

  memory->create(pn3, nabf3+1, "pn3");
  memory->create(pq3, K3*2, "pq3");
  memory->create(pc3, K3, "pc3");

  init3bodyarray(pn3, pq3, pc3, nabf3-1);
}

void FASTPOD::init4body(int Pa4)
{
  int npa[] = {0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455};
  int nb[] = {1,     2,     4,     7,    11,    16,    23};
  int ns[] = {0, 1, 4, 10, 19, 29, 47, 74, 89, 119, 155, 209, 230, 275, 335, 425, 533, 561, 624, 714, 849, 949, 1129, 1345};

  P4 = Pa4;
  K4 = npa[Pa4+1];

  int *pn4, *tm4;
  memory->create(pn4, Pa4+2, "pn4");
  memory->create(pq4, K4*2, "pq4");
  memory->create(tm4, K4, "tm4");

  init3bodyarray(pn4, pq4, tm4, Pa4);

  nabf4 = nb[Pa4];
  Q4 = ns[nabf4];
  memory->create(pa4, nabf4+1, "pa4");
  memory->create(pb4, Q4*3, "pb4");
  memory->create(pc4, Q4, "pc4");

  init4bodyarray(pa4, pb4, pc4, Pa4);

  memory->destroy(pn4);
  memory->destroy(tm4);
}

int FASTPOD::estimate_memory(int Nj)
{
  int Kmax = (K3 > K4) ? K3 : K4;
  int nrbf34 = (nrbf3 > nrbf4) ? nrbf3 : nrbf4;
  int nrbfmax = (nrbf2 > nrbf34) ? nrbf2 : nrbf34;
  int Knrbf34 = (K3*nrbf3 > K4*nrbf4) ? K3*nrbf3 : K4*nrbf4;

  int nld = (nl23 > nl33) ? nl23 : nl33;
  nld = (nld > nl34) ? nld : nl34;
  nld = (nld > nl44) ? nld : nl44;

  // rij, fij, and d2, dd2, d3, dd3, d4, dd4
  int nmax1 = 6*Nj + nl2 + 3*Nj*nl2 + nl3 + 3*Nj*nl3 + nl4 + 3*Nj*nl4 + nld + 3*Nj*nld;

  // U, Ux, Uy, Uz
  int nmax2 = 4*Nj*Knrbf34;

  // sumU and cU
  int nmax3 = 2*nelements*Knrbf34;

  // rbf, rbfx, rbfy, rbfz
  int nmax4 = 4*Nj*nrbfmax;

  // rbft, rbfxt, rbfyt, rbfzt
  int nmax5 = 4*Nj*ns;

  // abf, abfx, abfy, abfz
  int nmax6 = 4*(Nj+1)*Kmax;

  // U, Ux, Uy, Uz, sumU, cU, rbf, rbfx, rbfy, rbfz, abf, abfx, abfy, abfz
  int nmax7 = (nmax5 > nmax6) ? nmax5 : nmax6;
  int nmax8 = nmax2 + nmax3 + nmax4 + nmax7;

  // all double memory
  ndblmem = (nmax1 + nmax8);

  // ai, aj, ti, tj
  nintmem = 4*Nj;

  return ndblmem;
}

int FASTPOD::indexmap3(int *indx, int n1, int n2, int n3, int N1, int N2)
{
  int k = 0;
  for (int i3=0; i3<n3; i3++)
    for (int i2=0; i2<n2; i2++)
      for (int i1=0; i1<n1; i1++)
       {
        indx[k] = i1 + N1*i2 + N1*N2*i3;
        k += 1;
      }

  return k;
}

void FASTPOD::print_matrix(const char* desc, int m, int n, int* a, int lda )
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ )
    {
        for( j = 0; j < n; j++ ) printf( " %d", a[i+j*lda] );
        printf( "\n" );
    }
}

void FASTPOD::print_matrix(const char* desc, int m, int n, double* a, int lda )
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ )
    {
        for( j = 0; j < n; j++ ) printf( " %6.12f", a[i+j*lda] );
        printf( "\n" );
    }
}

void FASTPOD::scalarproduct(double *d, double c, int N)
{
  for (int n=0; n<N; n++)
    d[n] = d[n]*c;
}

double FASTPOD::dotproduct(double *c, double *d, int ndesc)
{
  double e = 0.0;
  for (int n = 0; n<ndesc; n++)
    e += d[n]*c[n];
  return e;
}

void FASTPOD::mvproduct(double *fij, double *c, double *dd, int N, int ndesc)
{
  for (int m=0; m<ndesc; m++)
    for (int n=0; n<N; n++)
      fij[n] += dd[n + N*m]*c[m];
}

void FASTPOD::MatMul(double *c, double *a, double *b, int r1, int c1, int c2)
{
  int i, j, k;

  for(j = 0; j < c2; j++)
    for(i = 0; i < r1; i++)
      c[i + r1*j] = 0.0;

  for(j = 0; j < c2; j++)
    for(i = 0; i < r1; i++)
      for(k = 0; k < c1; k++)
        c[i + r1*j] += a[i + r1*k] * b[k + c1*j];
}

void FASTPOD::init3bodyarray(int *np, int *pq, int *pc, int Pa)
{
  int npa[] = {0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455};

  int poly[455][6] =
  {
    {0, 0, 0, 0, 0, 1},
    {1, 0, 0, 1, 1, 1},
    {0, 1, 0, 1, 2, 1},
    {0, 0, 1, 1, 3, 1},
    {2, 0, 0, 2, 1, 1},
    {1, 1, 0, 2, 2, 2},
    {0, 2, 0, 3, 2, 1},
    {1, 0, 1, 2, 3, 2},
    {0, 1, 1, 3, 3, 2},
    {0, 0, 2, 4, 3, 1},
    {3, 0, 0, 5, 1, 1},
    {2, 1, 0, 5, 2, 3},
    {1, 2, 0, 6, 2, 3},
    {0, 3, 0, 7, 2, 1},
    {2, 0, 1, 5, 3, 3},
    {1, 1, 1, 6, 3, 6},
    {0, 2, 1, 7, 3, 3},
    {1, 0, 2, 8, 3, 3},
    {0, 1, 2, 9, 3, 3},
    {0, 0, 3, 10, 3, 1},
    {4, 0, 0, 11, 1, 1},
    {3, 1, 0, 11, 2, 4},
    {2, 2, 0, 12, 2, 6},
    {1, 3, 0, 13, 2, 4},
    {0, 4, 0, 14, 2, 1},
    {3, 0, 1, 11, 3, 4},
    {2, 1, 1, 12, 3, 12},
    {1, 2, 1, 13, 3, 12},
    {0, 3, 1, 14, 3, 4},
    {2, 0, 2, 15, 3, 6},
    {1, 1, 2, 16, 3, 12},
    {0, 2, 2, 17, 3, 6},
    {1, 0, 3, 18, 3, 4},
    {0, 1, 3, 19, 3, 4},
    {0, 0, 4, 20, 3, 1},
    {5, 0, 0, 21, 1, 1},
    {4, 1, 0, 21, 2, 5},
    {3, 2, 0, 22, 2, 10},
    {2, 3, 0, 23, 2, 10},
    {1, 4, 0, 24, 2, 5},
    {0, 5, 0, 25, 2, 1},
    {4, 0, 1, 21, 3, 5},
    {3, 1, 1, 22, 3, 20},
    {2, 2, 1, 23, 3, 30},
    {1, 3, 1, 24, 3, 20},
    {0, 4, 1, 25, 3, 5},
    {3, 0, 2, 26, 3, 10},
    {2, 1, 2, 27, 3, 30},
    {1, 2, 2, 28, 3, 30},
    {0, 3, 2, 29, 3, 10},
    {2, 0, 3, 30, 3, 10},
    {1, 1, 3, 31, 3, 20},
    {0, 2, 3, 32, 3, 10},
    {1, 0, 4, 33, 3, 5},
    {0, 1, 4, 34, 3, 5},
    {0, 0, 5, 35, 3, 1},
    {6, 0, 0, 36, 1, 1},
    {5, 1, 0, 36, 2, 6},
    {4, 2, 0, 37, 2, 15},
    {3, 3, 0, 38, 2, 20},
    {2, 4, 0, 39, 2, 15},
    {1, 5, 0, 40, 2, 6},
    {0, 6, 0, 41, 2, 1},
    {5, 0, 1, 36, 3, 6},
    {4, 1, 1, 37, 3, 30},
    {3, 2, 1, 38, 3, 60},
    {2, 3, 1, 39, 3, 60},
    {1, 4, 1, 40, 3, 30},
    {0, 5, 1, 41, 3, 6},
    {4, 0, 2, 42, 3, 15},
    {3, 1, 2, 43, 3, 60},
    {2, 2, 2, 44, 3, 90},
    {1, 3, 2, 45, 3, 60},
    {0, 4, 2, 46, 3, 15},
    {3, 0, 3, 47, 3, 20},
    {2, 1, 3, 48, 3, 60},
    {1, 2, 3, 49, 3, 60},
    {0, 3, 3, 50, 3, 20},
    {2, 0, 4, 51, 3, 15},
    {1, 1, 4, 52, 3, 30},
    {0, 2, 4, 53, 3, 15},
    {1, 0, 5, 54, 3, 6},
    {0, 1, 5, 55, 3, 6},
    {0, 0, 6, 56, 3, 1},
    {7, 0, 0, 57, 1, 1},
    {6, 1, 0, 57, 2, 7},
    {5, 2, 0, 58, 2, 21},
    {4, 3, 0, 59, 2, 35},
    {3, 4, 0, 60, 2, 35},
    {2, 5, 0, 61, 2, 21},
    {1, 6, 0, 62, 2, 7},
    {0, 7, 0, 63, 2, 1},
    {6, 0, 1, 57, 3, 7},
    {5, 1, 1, 58, 3, 42},
    {4, 2, 1, 59, 3, 105},
    {3, 3, 1, 60, 3, 140},
    {2, 4, 1, 61, 3, 105},
    {1, 5, 1, 62, 3, 42},
    {0, 6, 1, 63, 3, 7},
    {5, 0, 2, 64, 3, 21},
    {4, 1, 2, 65, 3, 105},
    {3, 2, 2, 66, 3, 210},
    {2, 3, 2, 67, 3, 210},
    {1, 4, 2, 68, 3, 105},
    {0, 5, 2, 69, 3, 21},
    {4, 0, 3, 70, 3, 35},
    {3, 1, 3, 71, 3, 140},
    {2, 2, 3, 72, 3, 210},
    {1, 3, 3, 73, 3, 140},
    {0, 4, 3, 74, 3, 35},
    {3, 0, 4, 75, 3, 35},
    {2, 1, 4, 76, 3, 105},
    {1, 2, 4, 77, 3, 105},
    {0, 3, 4, 78, 3, 35},
    {2, 0, 5, 79, 3, 21},
    {1, 1, 5, 80, 3, 42},
    {0, 2, 5, 81, 3, 21},
    {1, 0, 6, 82, 3, 7},
    {0, 1, 6, 83, 3, 7},
    {0, 0, 7, 84, 3, 1},
    {8, 0, 0, 85, 1, 1},
    {7, 1, 0, 85, 2, 8},
    {6, 2, 0, 86, 2, 28},
    {5, 3, 0, 87, 2, 56},
    {4, 4, 0, 88, 2, 70},
    {3, 5, 0, 89, 2, 56},
    {2, 6, 0, 90, 2, 28},
    {1, 7, 0, 91, 2, 8},
    {0, 8, 0, 92, 2, 1},
    {7, 0, 1, 85, 3, 8},
    {6, 1, 1, 86, 3, 56},
    {5, 2, 1, 87, 3, 168},
    {4, 3, 1, 88, 3, 280},
    {3, 4, 1, 89, 3, 280},
    {2, 5, 1, 90, 3, 168},
    {1, 6, 1, 91, 3, 56},
    {0, 7, 1, 92, 3, 8},
    {6, 0, 2, 93, 3, 28},
    {5, 1, 2, 94, 3, 168},
    {4, 2, 2, 95, 3, 420},
    {3, 3, 2, 96, 3, 560},
    {2, 4, 2, 97, 3, 420},
    {1, 5, 2, 98, 3, 168},
    {0, 6, 2, 99, 3, 28},
    {5, 0, 3, 100, 3, 56},
    {4, 1, 3, 101, 3, 280},
    {3, 2, 3, 102, 3, 560},
    {2, 3, 3, 103, 3, 560},
    {1, 4, 3, 104, 3, 280},
    {0, 5, 3, 105, 3, 56},
    {4, 0, 4, 106, 3, 70},
    {3, 1, 4, 107, 3, 280},
    {2, 2, 4, 108, 3, 420},
    {1, 3, 4, 109, 3, 280},
    {0, 4, 4, 110, 3, 70},
    {3, 0, 5, 111, 3, 56},
    {2, 1, 5, 112, 3, 168},
    {1, 2, 5, 113, 3, 168},
    {0, 3, 5, 114, 3, 56},
    {2, 0, 6, 115, 3, 28},
    {1, 1, 6, 116, 3, 56},
    {0, 2, 6, 117, 3, 28},
    {1, 0, 7, 118, 3, 8},
    {0, 1, 7, 119, 3, 8},
    {0, 0, 8, 120, 3, 1},
    {9, 0, 0, 121, 1, 1},
    {8, 1, 0, 121, 2, 9},
    {7, 2, 0, 122, 2, 36},
    {6, 3, 0, 123, 2, 84},
    {5, 4, 0, 124, 2, 126},
    {4, 5, 0, 125, 2, 126},
    {3, 6, 0, 126, 2, 84},
    {2, 7, 0, 127, 2, 36},
    {1, 8, 0, 128, 2, 9},
    {0, 9, 0, 129, 2, 1},
    {8, 0, 1, 121, 3, 9},
    {7, 1, 1, 122, 3, 72},
    {6, 2, 1, 123, 3, 252},
    {5, 3, 1, 124, 3, 504},
    {4, 4, 1, 125, 3, 630},
    {3, 5, 1, 126, 3, 504},
    {2, 6, 1, 127, 3, 252},
    {1, 7, 1, 128, 3, 72},
    {0, 8, 1, 129, 3, 9},
    {7, 0, 2, 130, 3, 36},
    {6, 1, 2, 131, 3, 252},
    {5, 2, 2, 132, 3, 756},
    {4, 3, 2, 133, 3, 1260},
    {3, 4, 2, 134, 3, 1260},
    {2, 5, 2, 135, 3, 756},
    {1, 6, 2, 136, 3, 252},
    {0, 7, 2, 137, 3, 36},
    {6, 0, 3, 138, 3, 84},
    {5, 1, 3, 139, 3, 504},
    {4, 2, 3, 140, 3, 1260},
    {3, 3, 3, 141, 3, 1680},
    {2, 4, 3, 142, 3, 1260},
    {1, 5, 3, 143, 3, 504},
    {0, 6, 3, 144, 3, 84},
    {5, 0, 4, 145, 3, 126},
    {4, 1, 4, 146, 3, 630},
    {3, 2, 4, 147, 3, 1260},
    {2, 3, 4, 148, 3, 1260},
    {1, 4, 4, 149, 3, 630},
    {0, 5, 4, 150, 3, 126},
    {4, 0, 5, 151, 3, 126},
    {3, 1, 5, 152, 3, 504},
    {2, 2, 5, 153, 3, 756},
    {1, 3, 5, 154, 3, 504},
    {0, 4, 5, 155, 3, 126},
    {3, 0, 6, 156, 3, 84},
    {2, 1, 6, 157, 3, 252},
    {1, 2, 6, 158, 3, 252},
    {0, 3, 6, 159, 3, 84},
    {2, 0, 7, 160, 3, 36},
    {1, 1, 7, 161, 3, 72},
    {0, 2, 7, 162, 3, 36},
    {1, 0, 8, 163, 3, 9},
    {0, 1, 8, 164, 3, 9},
    {0, 0, 9, 165, 3, 1},
    {10, 0, 0, 166, 1, 1},
    {9, 1, 0, 166, 2, 10},
    {8, 2, 0, 167, 2, 45},
    {7, 3, 0, 168, 2, 120},
    {6, 4, 0, 169, 2, 210},
    {5, 5, 0, 170, 2, 252},
    {4, 6, 0, 171, 2, 210},
    {3, 7, 0, 172, 2, 120},
    {2, 8, 0, 173, 2, 45},
    {1, 9, 0, 174, 2, 10},
    {0, 10, 0, 175, 2, 1},
    {9, 0, 1, 166, 3, 10},
    {8, 1, 1, 167, 3, 90},
    {7, 2, 1, 168, 3, 360},
    {6, 3, 1, 169, 3, 840},
    {5, 4, 1, 170, 3, 1260},
    {4, 5, 1, 171, 3, 1260},
    {3, 6, 1, 172, 3, 840},
    {2, 7, 1, 173, 3, 360},
    {1, 8, 1, 174, 3, 90},
    {0, 9, 1, 175, 3, 10},
    {8, 0, 2, 176, 3, 45},
    {7, 1, 2, 177, 3, 360},
    {6, 2, 2, 178, 3, 1260},
    {5, 3, 2, 179, 3, 2520},
    {4, 4, 2, 180, 3, 3150},
    {3, 5, 2, 181, 3, 2520},
    {2, 6, 2, 182, 3, 1260},
    {1, 7, 2, 183, 3, 360},
    {0, 8, 2, 184, 3, 45},
    {7, 0, 3, 185, 3, 120},
    {6, 1, 3, 186, 3, 840},
    {5, 2, 3, 187, 3, 2520},
    {4, 3, 3, 188, 3, 4200},
    {3, 4, 3, 189, 3, 4200},
    {2, 5, 3, 190, 3, 2520},
    {1, 6, 3, 191, 3, 840},
    {0, 7, 3, 192, 3, 120},
    {6, 0, 4, 193, 3, 210},
    {5, 1, 4, 194, 3, 1260},
    {4, 2, 4, 195, 3, 3150},
    {3, 3, 4, 196, 3, 4200},
    {2, 4, 4, 197, 3, 3150},
    {1, 5, 4, 198, 3, 1260},
    {0, 6, 4, 199, 3, 210},
    {5, 0, 5, 200, 3, 252},
    {4, 1, 5, 201, 3, 1260},
    {3, 2, 5, 202, 3, 2520},
    {2, 3, 5, 203, 3, 2520},
    {1, 4, 5, 204, 3, 1260},
    {0, 5, 5, 205, 3, 252},
    {4, 0, 6, 206, 3, 210},
    {3, 1, 6, 207, 3, 840},
    {2, 2, 6, 208, 3, 1260},
    {1, 3, 6, 209, 3, 840},
    {0, 4, 6, 210, 3, 210},
    {3, 0, 7, 211, 3, 120},
    {2, 1, 7, 212, 3, 360},
    {1, 2, 7, 213, 3, 360},
    {0, 3, 7, 214, 3, 120},
    {2, 0, 8, 215, 3, 45},
    {1, 1, 8, 216, 3, 90},
    {0, 2, 8, 217, 3, 45},
    {1, 0, 9, 218, 3, 10},
    {0, 1, 9, 219, 3, 10},
    {0, 0, 10, 220, 3, 1},
    {11, 0, 0, 221, 1, 1},
    {10, 1, 0, 221, 2, 11},
    {9, 2, 0, 222, 2, 55},
    {8, 3, 0, 223, 2, 165},
    {7, 4, 0, 224, 2, 330},
    {6, 5, 0, 225, 2, 462},
    {5, 6, 0, 226, 2, 462},
    {4, 7, 0, 227, 2, 330},
    {3, 8, 0, 228, 2, 165},
    {2, 9, 0, 229, 2, 55},
    {1, 10, 0, 230, 2, 11},
    {0, 11, 0, 231, 2, 1},
    {10, 0, 1, 221, 3, 11},
    {9, 1, 1, 222, 3, 110},
    {8, 2, 1, 223, 3, 495},
    {7, 3, 1, 224, 3, 1320},
    {6, 4, 1, 225, 3, 2310},
    {5, 5, 1, 226, 3, 2772},
    {4, 6, 1, 227, 3, 2310},
    {3, 7, 1, 228, 3, 1320},
    {2, 8, 1, 229, 3, 495},
    {1, 9, 1, 230, 3, 110},
    {0, 10, 1, 231, 3, 11},
    {9, 0, 2, 232, 3, 55},
    {8, 1, 2, 233, 3, 495},
    {7, 2, 2, 234, 3, 1980},
    {6, 3, 2, 235, 3, 4620},
    {5, 4, 2, 236, 3, 6930},
    {4, 5, 2, 237, 3, 6930},
    {3, 6, 2, 238, 3, 4620},
    {2, 7, 2, 239, 3, 1980},
    {1, 8, 2, 240, 3, 495},
    {0, 9, 2, 241, 3, 55},
    {8, 0, 3, 242, 3, 165},
    {7, 1, 3, 243, 3, 1320},
    {6, 2, 3, 244, 3, 4620},
    {5, 3, 3, 245, 3, 9240},
    {4, 4, 3, 246, 3, 11550},
    {3, 5, 3, 247, 3, 9240},
    {2, 6, 3, 248, 3, 4620},
    {1, 7, 3, 249, 3, 1320},
    {0, 8, 3, 250, 3, 165},
    {7, 0, 4, 251, 3, 330},
    {6, 1, 4, 252, 3, 2310},
    {5, 2, 4, 253, 3, 6930},
    {4, 3, 4, 254, 3, 11550},
    {3, 4, 4, 255, 3, 11550},
    {2, 5, 4, 256, 3, 6930},
    {1, 6, 4, 257, 3, 2310},
    {0, 7, 4, 258, 3, 330},
    {6, 0, 5, 259, 3, 462},
    {5, 1, 5, 260, 3, 2772},
    {4, 2, 5, 261, 3, 6930},
    {3, 3, 5, 262, 3, 9240},
    {2, 4, 5, 263, 3, 6930},
    {1, 5, 5, 264, 3, 2772},
    {0, 6, 5, 265, 3, 462},
    {5, 0, 6, 266, 3, 462},
    {4, 1, 6, 267, 3, 2310},
    {3, 2, 6, 268, 3, 4620},
    {2, 3, 6, 269, 3, 4620},
    {1, 4, 6, 270, 3, 2310},
    {0, 5, 6, 271, 3, 462},
    {4, 0, 7, 272, 3, 330},
    {3, 1, 7, 273, 3, 1320},
    {2, 2, 7, 274, 3, 1980},
    {1, 3, 7, 275, 3, 1320},
    {0, 4, 7, 276, 3, 330},
    {3, 0, 8, 277, 3, 165},
    {2, 1, 8, 278, 3, 495},
    {1, 2, 8, 279, 3, 495},
    {0, 3, 8, 280, 3, 165},
    {2, 0, 9, 281, 3, 55},
    {1, 1, 9, 282, 3, 110},
    {0, 2, 9, 283, 3, 55},
    {1, 0, 10, 284, 3, 11},
    {0, 1, 10, 285, 3, 11},
    {0, 0, 11, 286, 3, 1},
    {12, 0, 0, 287, 1, 1},
    {11, 1, 0, 287, 2, 12},
    {10, 2, 0, 288, 2, 66},
    {9, 3, 0, 289, 2, 220},
    {8, 4, 0, 290, 2, 495},
    {7, 5, 0, 291, 2, 792},
    {6, 6, 0, 292, 2, 924},
    {5, 7, 0, 293, 2, 792},
    {4, 8, 0, 294, 2, 495},
    {3, 9, 0, 295, 2, 220},
    {2, 10, 0, 296, 2, 66},
    {1, 11, 0, 297, 2, 12},
    {0, 12, 0, 298, 2, 1},
    {11, 0, 1, 287, 3, 12},
    {10, 1, 1, 288, 3, 132},
    {9, 2, 1, 289, 3, 660},
    {8, 3, 1, 290, 3, 1980},
    {7, 4, 1, 291, 3, 3960},
    {6, 5, 1, 292, 3, 5544},
    {5, 6, 1, 293, 3, 5544},
    {4, 7, 1, 294, 3, 3960},
    {3, 8, 1, 295, 3, 1980},
    {2, 9, 1, 296, 3, 660},
    {1, 10, 1, 297, 3, 132},
    {0, 11, 1, 298, 3, 12},
    {10, 0, 2, 299, 3, 66},
    {9, 1, 2, 300, 3, 660},
    {8, 2, 2, 301, 3, 2970},
    {7, 3, 2, 302, 3, 7920},
    {6, 4, 2, 303, 3, 13860},
    {5, 5, 2, 304, 3, 16632},
    {4, 6, 2, 305, 3, 13860},
    {3, 7, 2, 306, 3, 7920},
    {2, 8, 2, 307, 3, 2970},
    {1, 9, 2, 308, 3, 660},
    {0, 10, 2, 309, 3, 66},
    {9, 0, 3, 310, 3, 220},
    {8, 1, 3, 311, 3, 1980},
    {7, 2, 3, 312, 3, 7920},
    {6, 3, 3, 313, 3, 18480},
    {5, 4, 3, 314, 3, 27720},
    {4, 5, 3, 315, 3, 27720},
    {3, 6, 3, 316, 3, 18480},
    {2, 7, 3, 317, 3, 7920},
    {1, 8, 3, 318, 3, 1980},
    {0, 9, 3, 319, 3, 220},
    {8, 0, 4, 320, 3, 495},
    {7, 1, 4, 321, 3, 3960},
    {6, 2, 4, 322, 3, 13860},
    {5, 3, 4, 323, 3, 27720},
    {4, 4, 4, 324, 3, 34650},
    {3, 5, 4, 325, 3, 27720},
    {2, 6, 4, 326, 3, 13860},
    {1, 7, 4, 327, 3, 3960},
    {0, 8, 4, 328, 3, 495},
    {7, 0, 5, 329, 3, 792},
    {6, 1, 5, 330, 3, 5544},
    {5, 2, 5, 331, 3, 16632},
    {4, 3, 5, 332, 3, 27720},
    {3, 4, 5, 333, 3, 27720},
    {2, 5, 5, 334, 3, 16632},
    {1, 6, 5, 335, 3, 5544},
    {0, 7, 5, 336, 3, 792},
    {6, 0, 6, 337, 3, 924},
    {5, 1, 6, 338, 3, 5544},
    {4, 2, 6, 339, 3, 13860},
    {3, 3, 6, 340, 3, 18480},
    {2, 4, 6, 341, 3, 13860},
    {1, 5, 6, 342, 3, 5544},
    {0, 6, 6, 343, 3, 924},
    {5, 0, 7, 344, 3, 792},
    {4, 1, 7, 345, 3, 3960},
    {3, 2, 7, 346, 3, 7920},
    {2, 3, 7, 347, 3, 7920},
    {1, 4, 7, 348, 3, 3960},
    {0, 5, 7, 349, 3, 792},
    {4, 0, 8, 350, 3, 495},
    {3, 1, 8, 351, 3, 1980},
    {2, 2, 8, 352, 3, 2970},
    {1, 3, 8, 353, 3, 1980},
    {0, 4, 8, 354, 3, 495},
    {3, 0, 9, 355, 3, 220},
    {2, 1, 9, 356, 3, 660},
    {1, 2, 9, 357, 3, 660},
    {0, 3, 9, 358, 3, 220},
    {2, 0, 10, 359, 3, 66},
    {1, 1, 10, 360, 3, 132},
    {0, 2, 10, 361, 3, 66},
    {1, 0, 11, 362, 3, 12},
    {0, 1, 11, 363, 3, 12},
    {0, 0, 12, 364, 3, 1}
  };

  for (int i = 0; i<= Pa+1; i++)
    np[i] = npa[i];

  int nmax = np[Pa+1];

  for (int i=0; i<nmax; i++) {
    pq[i]        = poly[i][3];
    pq[i+nmax]   = poly[i][4];
    pc[i]        = poly[i][5];
  }
}

void FASTPOD::init4bodyarray(int *ns4, int *pb4, int *pc4, int Pa)
{
  int nb[] = {1,     2,     4,     7,    11,    16,    23};

  int ns[] =
  {
    0,
    1,
    4,
    10,
    19,
    29,
    47,
    74,
    89,
    119,
    155,
    209,
    230,
    275,
    335,
    425,
    533,
    561,
    624,
    714,
    849,
    949,
    1129,
    1345
  };

  int poly[1345][4] =
  {
    {1, 1, 1, 1},
    {2, 2, 1, 1},
    {3, 3, 1, 1},
    {4, 4, 1, 1},
    {5, 5, 1, 1},
    {6, 6, 1, 2},
    {8, 8, 1, 2},
    {7, 7, 1, 1},
    {9, 9, 1, 2},
    {10, 10, 1, 1},
    {5, 2, 2, 1},
    {6, 2, 3, 1},
    {6, 3, 2, 1},
    {8, 2, 4, 1},
    {8, 4, 2, 1},
    {7, 3, 3, 1},
    {9, 3, 4, 1},
    {9, 4, 3, 1},
    {10, 4, 4, 1},
    {11, 11, 1, 1},
    {12, 12, 1, 3},
    {15, 15, 1, 3},
    {13, 13, 1, 3},
    {16, 16, 1, 6},
    {18, 18, 1, 3},
    {14, 14, 1, 1},
    {17, 17, 1, 3},
    {19, 19, 1, 3},
    {20, 20, 1, 1},
    {11, 5, 2, 1},
    {12, 5, 3, 1},
    {12, 6, 2, 2},
    {15, 5, 4, 1},
    {15, 8, 2, 2},
    {13, 6, 3, 2},
    {13, 7, 2, 1},
    {16, 6, 4, 2},
    {16, 8, 3, 2},
    {16, 9, 2, 2},
    {18, 8, 4, 2},
    {18, 10, 2, 1},
    {14, 7, 3, 1},
    {17, 7, 4, 1},
    {17, 9, 3, 2},
    {19, 9, 4, 2},
    {19, 10, 3, 1},
    {20, 10, 4, 1},
    {5, 5, 5, 1},
    {5, 6, 6, 1},
    {5, 8, 8, 1},
    {6, 5, 6, 1},
    {6, 6, 5, 1},
    {6, 6, 7, 1},
    {6, 8, 9, 1},
    {6, 7, 6, 1},
    {6, 9, 8, 1},
    {8, 5, 8, 1},
    {8, 6, 9, 1},
    {8, 8, 5, 1},
    {8, 8, 10, 1},
    {8, 9, 6, 1},
    {8, 10, 8, 1},
    {7, 6, 6, 1},
    {7, 7, 7, 1},
    {7, 9, 9, 1},
    {9, 6, 8, 1},
    {9, 8, 6, 1},
    {9, 7, 9, 1},
    {9, 9, 7, 1},
    {9, 9, 10, 1},
    {9, 10, 9, 1},
    {10, 8, 8, 1},
    {10, 9, 9, 1},
    {10, 10, 10, 1},
    {21, 21, 1, 1},
    {22, 22, 1, 4},
    {26, 26, 1, 4},
    {23, 23, 1, 6},
    {27, 27, 1, 12},
    {30, 30, 1, 6},
    {24, 24, 1, 4},
    {28, 28, 1, 12},
    {31, 31, 1, 12},
    {33, 33, 1, 4},
    {25, 25, 1, 1},
    {29, 29, 1, 4},
    {32, 32, 1, 6},
    {34, 34, 1, 4},
    {35, 35, 1, 1},
    {21, 11, 2, 1},
    {22, 11, 3, 1},
    {22, 12, 2, 3},
    {26, 11, 4, 1},
    {26, 15, 2, 3},
    {23, 12, 3, 3},
    {23, 13, 2, 3},
    {27, 12, 4, 3},
    {27, 15, 3, 3},
    {27, 16, 2, 6},
    {30, 15, 4, 3},
    {30, 18, 2, 3},
    {24, 13, 3, 3},
    {24, 14, 2, 1},
    {28, 13, 4, 3},
    {28, 16, 3, 6},
    {28, 17, 2, 3},
    {31, 16, 4, 6},
    {31, 18, 3, 3},
    {31, 19, 2, 3},
    {33, 18, 4, 3},
    {33, 20, 2, 1},
    {25, 14, 3, 1},
    {29, 14, 4, 1},
    {29, 17, 3, 3},
    {32, 17, 4, 3},
    {32, 19, 3, 3},
    {34, 19, 4, 3},
    {34, 20, 3, 1},
    {35, 20, 4, 1},
    {21, 5, 5, 1},
    {22, 5, 6, 2},
    {22, 6, 5, 2},
    {26, 5, 8, 2},
    {26, 8, 5, 2},
    {23, 5, 7, 1},
    {23, 6, 6, 4},
    {23, 7, 5, 1},
    {27, 5, 9, 2},
    {27, 6, 8, 4},
    {27, 8, 6, 4},
    {27, 9, 5, 2},
    {30, 5, 10, 1},
    {30, 8, 8, 4},
    {30, 10, 5, 1},
    {24, 6, 7, 2},
    {24, 7, 6, 2},
    {28, 6, 9, 4},
    {28, 8, 7, 2},
    {28, 7, 8, 2},
    {28, 9, 6, 4},
    {31, 6, 10, 2},
    {31, 8, 9, 4},
    {31, 9, 8, 4},
    {31, 10, 6, 2},
    {33, 8, 10, 2},
    {33, 10, 8, 2},
    {25, 7, 7, 1},
    {29, 7, 9, 2},
    {29, 9, 7, 2},
    {32, 7, 10, 1},
    {32, 9, 9, 4},
    {32, 10, 7, 1},
    {34, 9, 10, 2},
    {34, 10, 9, 2},
    {35, 10, 10, 1},
    {11, 11, 5, 1},
    {11, 12, 6, 1},
    {11, 15, 8, 1},
    {12, 11, 6, 1},
    {12, 12, 5, 2},
    {12, 12, 7, 1},
    {12, 15, 9, 1},
    {12, 13, 6, 2},
    {12, 16, 8, 2},
    {15, 11, 8, 1},
    {15, 12, 9, 1},
    {15, 15, 5, 2},
    {15, 15, 10, 1},
    {15, 16, 6, 2},
    {15, 18, 8, 2},
    {13, 12, 6, 2},
    {13, 13, 5, 1},
    {13, 13, 7, 2},
    {13, 16, 9, 2},
    {13, 14, 6, 1},
    {13, 17, 8, 1},
    {16, 12, 8, 2},
    {16, 15, 6, 2},
    {16, 13, 9, 2},
    {16, 16, 5, 2},
    {16, 16, 7, 2},
    {16, 16, 10, 2},
    {16, 18, 9, 2},
    {16, 17, 6, 2},
    {16, 19, 8, 2},
    {18, 15, 8, 2},
    {18, 16, 9, 2},
    {18, 18, 5, 1},
    {18, 18, 10, 2},
    {18, 19, 6, 1},
    {18, 20, 8, 1},
    {14, 13, 6, 1},
    {14, 14, 7, 1},
    {14, 17, 9, 1},
    {17, 13, 8, 1},
    {17, 16, 6, 2},
    {17, 14, 9, 1},
    {17, 17, 7, 2},
    {17, 17, 10, 1},
    {17, 19, 9, 2},
    {19, 16, 8, 2},
    {19, 18, 6, 1},
    {19, 17, 9, 2},
    {19, 19, 7, 1},
    {19, 19, 10, 2},
    {19, 20, 9, 1},
    {20, 18, 8, 1},
    {20, 19, 9, 1},
    {20, 20, 10, 1},
    {36, 36, 1, 1},
    {37, 37, 1, 5},
    {42, 42, 1, 5},
    {38, 38, 1, 10},
    {43, 43, 1, 20},
    {47, 47, 1, 10},
    {39, 39, 1, 10},
    {44, 44, 1, 30},
    {48, 48, 1, 30},
    {51, 51, 1, 10},
    {40, 40, 1, 5},
    {45, 45, 1, 20},
    {49, 49, 1, 30},
    {52, 52, 1, 20},
    {54, 54, 1, 5},
    {41, 41, 1, 1},
    {46, 46, 1, 5},
    {50, 50, 1, 10},
    {53, 53, 1, 10},
    {55, 55, 1, 5},
    {56, 56, 1, 1},
    {36, 21, 2, 1},
    {37, 21, 3, 1},
    {37, 22, 2, 4},
    {42, 21, 4, 1},
    {42, 26, 2, 4},
    {38, 22, 3, 4},
    {38, 23, 2, 6},
    {43, 22, 4, 4},
    {43, 26, 3, 4},
    {43, 27, 2, 12},
    {47, 26, 4, 4},
    {47, 30, 2, 6},
    {39, 23, 3, 6},
    {39, 24, 2, 4},
    {44, 23, 4, 6},
    {44, 27, 3, 12},
    {44, 28, 2, 12},
    {48, 27, 4, 12},
    {48, 30, 3, 6},
    {48, 31, 2, 12},
    {51, 30, 4, 6},
    {51, 33, 2, 4},
    {40, 24, 3, 4},
    {40, 25, 2, 1},
    {45, 24, 4, 4},
    {45, 28, 3, 12},
    {45, 29, 2, 4},
    {49, 28, 4, 12},
    {49, 31, 3, 12},
    {49, 32, 2, 6},
    {52, 31, 4, 12},
    {52, 33, 3, 4},
    {52, 34, 2, 4},
    {54, 33, 4, 4},
    {54, 35, 2, 1},
    {41, 25, 3, 1},
    {46, 25, 4, 1},
    {46, 29, 3, 4},
    {50, 29, 4, 4},
    {50, 32, 3, 6},
    {53, 32, 4, 6},
    {53, 34, 3, 4},
    {55, 34, 4, 4},
    {55, 35, 3, 1},
    {56, 35, 4, 1},
    {36, 11, 5, 1},
    {37, 11, 6, 2},
    {37, 12, 5, 3},
    {42, 11, 8, 2},
    {42, 15, 5, 3},
    {38, 11, 7, 1},
    {38, 12, 6, 6},
    {38, 13, 5, 3},
    {43, 11, 9, 2},
    {43, 12, 8, 6},
    {43, 15, 6, 6},
    {43, 16, 5, 6},
    {47, 11, 10, 1},
    {47, 15, 8, 6},
    {47, 18, 5, 3},
    {39, 12, 7, 3},
    {39, 13, 6, 6},
    {39, 14, 5, 1},
    {44, 12, 9, 6},
    {44, 15, 7, 3},
    {44, 13, 8, 6},
    {44, 16, 6, 12},
    {44, 17, 5, 3},
    {48, 12, 10, 3},
    {48, 15, 9, 6},
    {48, 16, 8, 12},
    {48, 18, 6, 6},
    {48, 19, 5, 3},
    {51, 15, 10, 3},
    {51, 18, 8, 6},
    {51, 20, 5, 1},
    {40, 13, 7, 3},
    {40, 14, 6, 2},
    {45, 13, 9, 6},
    {45, 16, 7, 6},
    {45, 14, 8, 2},
    {45, 17, 6, 6},
    {49, 13, 10, 3},
    {49, 16, 9, 12},
    {49, 18, 7, 3},
    {49, 17, 8, 6},
    {49, 19, 6, 6},
    {52, 16, 10, 6},
    {52, 18, 9, 6},
    {52, 19, 8, 6},
    {52, 20, 6, 2},
    {54, 18, 10, 3},
    {54, 20, 8, 2},
    {41, 14, 7, 1},
    {46, 14, 9, 2},
    {46, 17, 7, 3},
    {50, 14, 10, 1},
    {50, 17, 9, 6},
    {50, 19, 7, 3},
    {53, 17, 10, 3},
    {53, 19, 9, 6},
    {53, 20, 7, 1},
    {55, 19, 10, 3},
    {55, 20, 9, 2},
    {56, 20, 10, 1},
    {21, 21, 5, 1},
    {21, 22, 6, 1},
    {21, 26, 8, 1},
    {22, 21, 6, 1},
    {22, 22, 5, 3},
    {22, 22, 7, 1},
    {22, 26, 9, 1},
    {22, 23, 6, 3},
    {22, 27, 8, 3},
    {26, 21, 8, 1},
    {26, 22, 9, 1},
    {26, 26, 5, 3},
    {26, 26, 10, 1},
    {26, 27, 6, 3},
    {26, 30, 8, 3},
    {23, 22, 6, 3},
    {23, 23, 5, 3},
    {23, 23, 7, 3},
    {23, 27, 9, 3},
    {23, 24, 6, 3},
    {23, 28, 8, 3},
    {27, 22, 8, 3},
    {27, 26, 6, 3},
    {27, 23, 9, 3},
    {27, 27, 5, 6},
    {27, 27, 7, 3},
    {27, 27, 10, 3},
    {27, 30, 9, 3},
    {27, 28, 6, 6},
    {27, 31, 8, 6},
    {30, 26, 8, 3},
    {30, 27, 9, 3},
    {30, 30, 5, 3},
    {30, 30, 10, 3},
    {30, 31, 6, 3},
    {30, 33, 8, 3},
    {24, 23, 6, 3},
    {24, 24, 5, 1},
    {24, 24, 7, 3},
    {24, 28, 9, 3},
    {24, 25, 6, 1},
    {24, 29, 8, 1},
    {28, 23, 8, 3},
    {28, 27, 6, 6},
    {28, 24, 9, 3},
    {28, 28, 5, 3},
    {28, 28, 7, 6},
    {28, 28, 10, 3},
    {28, 31, 9, 6},
    {28, 29, 6, 3},
    {28, 32, 8, 3},
    {31, 27, 8, 6},
    {31, 30, 6, 3},
    {31, 28, 9, 6},
    {31, 31, 5, 3},
    {31, 31, 7, 3},
    {31, 31, 10, 6},
    {31, 33, 9, 3},
    {31, 32, 6, 3},
    {31, 34, 8, 3},
    {33, 30, 8, 3},
    {33, 31, 9, 3},
    {33, 33, 5, 1},
    {33, 33, 10, 3},
    {33, 34, 6, 1},
    {33, 35, 8, 1},
    {25, 24, 6, 1},
    {25, 25, 7, 1},
    {25, 29, 9, 1},
    {29, 24, 8, 1},
    {29, 28, 6, 3},
    {29, 25, 9, 1},
    {29, 29, 7, 3},
    {29, 29, 10, 1},
    {29, 32, 9, 3},
    {32, 28, 8, 3},
    {32, 31, 6, 3},
    {32, 29, 9, 3},
    {32, 32, 7, 3},
    {32, 32, 10, 3},
    {32, 34, 9, 3},
    {34, 31, 8, 3},
    {34, 33, 6, 1},
    {34, 32, 9, 3},
    {34, 34, 7, 1},
    {34, 34, 10, 3},
    {34, 35, 9, 1},
    {35, 33, 8, 1},
    {35, 34, 9, 1},
    {35, 35, 10, 1},
    {21, 11, 11, 1},
    {21, 12, 12, 1},
    {21, 15, 15, 1},
    {22, 11, 12, 2},
    {22, 12, 11, 2},
    {22, 12, 13, 2},
    {22, 15, 16, 2},
    {22, 13, 12, 2},
    {22, 16, 15, 2},
    {26, 11, 15, 2},
    {26, 12, 16, 2},
    {26, 15, 11, 2},
    {26, 15, 18, 2},
    {26, 16, 12, 2},
    {26, 18, 15, 2},
    {23, 11, 13, 1},
    {23, 12, 12, 4},
    {23, 12, 14, 1},
    {23, 15, 17, 1},
    {23, 13, 11, 1},
    {23, 13, 13, 4},
    {23, 16, 16, 4},
    {23, 14, 12, 1},
    {23, 17, 15, 1},
    {27, 11, 16, 2},
    {27, 12, 15, 4},
    {27, 12, 17, 2},
    {27, 15, 12, 4},
    {27, 15, 19, 2},
    {27, 13, 16, 4},
    {27, 16, 11, 2},
    {27, 16, 13, 4},
    {27, 16, 18, 4},
    {27, 18, 16, 4},
    {27, 17, 12, 2},
    {27, 19, 15, 2},
    {30, 11, 18, 1},
    {30, 12, 19, 1},
    {30, 15, 15, 4},
    {30, 15, 20, 1},
    {30, 16, 16, 4},
    {30, 18, 11, 1},
    {30, 18, 18, 4},
    {30, 19, 12, 1},
    {30, 20, 15, 1},
    {24, 12, 13, 2},
    {24, 13, 12, 2},
    {24, 13, 14, 2},
    {24, 16, 17, 2},
    {24, 14, 13, 2},
    {24, 17, 16, 2},
    {28, 12, 16, 4},
    {28, 15, 13, 2},
    {28, 13, 15, 2},
    {28, 13, 17, 4},
    {28, 16, 12, 4},
    {28, 16, 14, 2},
    {28, 16, 19, 4},
    {28, 18, 17, 2},
    {28, 14, 16, 2},
    {28, 17, 13, 4},
    {28, 17, 18, 2},
    {28, 19, 16, 4},
    {31, 12, 18, 2},
    {31, 15, 16, 4},
    {31, 13, 19, 2},
    {31, 16, 15, 4},
    {31, 16, 17, 4},
    {31, 16, 20, 2},
    {31, 18, 12, 2},
    {31, 18, 19, 4},
    {31, 17, 16, 4},
    {31, 19, 13, 2},
    {31, 19, 18, 4},
    {31, 20, 16, 2},
    {33, 15, 18, 2},
    {33, 16, 19, 2},
    {33, 18, 15, 2},
    {33, 18, 20, 2},
    {33, 19, 16, 2},
    {33, 20, 18, 2},
    {25, 13, 13, 1},
    {25, 14, 14, 1},
    {25, 17, 17, 1},
    {29, 13, 16, 2},
    {29, 16, 13, 2},
    {29, 14, 17, 2},
    {29, 17, 14, 2},
    {29, 17, 19, 2},
    {29, 19, 17, 2},
    {32, 13, 18, 1},
    {32, 16, 16, 4},
    {32, 18, 13, 1},
    {32, 14, 19, 1},
    {32, 17, 17, 4},
    {32, 17, 20, 1},
    {32, 19, 14, 1},
    {32, 19, 19, 4},
    {32, 20, 17, 1},
    {34, 16, 18, 2},
    {34, 18, 16, 2},
    {34, 17, 19, 2},
    {34, 19, 17, 2},
    {34, 19, 20, 2},
    {34, 20, 19, 2},
    {35, 18, 18, 1},
    {35, 19, 19, 1},
    {35, 20, 20, 1},
    {57, 57, 1, 1},
    {58, 58, 1, 6},
    {64, 64, 1, 6},
    {59, 59, 1, 15},
    {65, 65, 1, 30},
    {70, 70, 1, 15},
    {60, 60, 1, 20},
    {66, 66, 1, 60},
    {71, 71, 1, 60},
    {75, 75, 1, 20},
    {61, 61, 1, 15},
    {67, 67, 1, 60},
    {72, 72, 1, 90},
    {76, 76, 1, 60},
    {79, 79, 1, 15},
    {62, 62, 1, 6},
    {68, 68, 1, 30},
    {73, 73, 1, 60},
    {77, 77, 1, 60},
    {80, 80, 1, 30},
    {82, 82, 1, 6},
    {63, 63, 1, 1},
    {69, 69, 1, 6},
    {74, 74, 1, 15},
    {78, 78, 1, 20},
    {81, 81, 1, 15},
    {83, 83, 1, 6},
    {84, 84, 1, 1},
    {57, 36, 2, 1},
    {58, 36, 3, 1},
    {58, 37, 2, 5},
    {64, 36, 4, 1},
    {64, 42, 2, 5},
    {59, 37, 3, 5},
    {59, 38, 2, 10},
    {65, 37, 4, 5},
    {65, 42, 3, 5},
    {65, 43, 2, 20},
    {70, 42, 4, 5},
    {70, 47, 2, 10},
    {60, 38, 3, 10},
    {60, 39, 2, 10},
    {66, 38, 4, 10},
    {66, 43, 3, 20},
    {66, 44, 2, 30},
    {71, 43, 4, 20},
    {71, 47, 3, 10},
    {71, 48, 2, 30},
    {75, 47, 4, 10},
    {75, 51, 2, 10},
    {61, 39, 3, 10},
    {61, 40, 2, 5},
    {67, 39, 4, 10},
    {67, 44, 3, 30},
    {67, 45, 2, 20},
    {72, 44, 4, 30},
    {72, 48, 3, 30},
    {72, 49, 2, 30},
    {76, 48, 4, 30},
    {76, 51, 3, 10},
    {76, 52, 2, 20},
    {79, 51, 4, 10},
    {79, 54, 2, 5},
    {62, 40, 3, 5},
    {62, 41, 2, 1},
    {68, 40, 4, 5},
    {68, 45, 3, 20},
    {68, 46, 2, 5},
    {73, 45, 4, 20},
    {73, 49, 3, 30},
    {73, 50, 2, 10},
    {77, 49, 4, 30},
    {77, 52, 3, 20},
    {77, 53, 2, 10},
    {80, 52, 4, 20},
    {80, 54, 3, 5},
    {80, 55, 2, 5},
    {82, 54, 4, 5},
    {82, 56, 2, 1},
    {63, 41, 3, 1},
    {69, 41, 4, 1},
    {69, 46, 3, 5},
    {74, 46, 4, 5},
    {74, 50, 3, 10},
    {78, 50, 4, 10},
    {78, 53, 3, 10},
    {81, 53, 4, 10},
    {81, 55, 3, 5},
    {83, 55, 4, 5},
    {83, 56, 3, 1},
    {84, 56, 4, 1},
    {57, 21, 5, 1},
    {58, 21, 6, 2},
    {58, 22, 5, 4},
    {64, 21, 8, 2},
    {64, 26, 5, 4},
    {59, 21, 7, 1},
    {59, 22, 6, 8},
    {59, 23, 5, 6},
    {65, 21, 9, 2},
    {65, 22, 8, 8},
    {65, 26, 6, 8},
    {65, 27, 5, 12},
    {70, 21, 10, 1},
    {70, 26, 8, 8},
    {70, 30, 5, 6},
    {60, 22, 7, 4},
    {60, 23, 6, 12},
    {60, 24, 5, 4},
    {66, 22, 9, 8},
    {66, 26, 7, 4},
    {66, 23, 8, 12},
    {66, 27, 6, 24},
    {66, 28, 5, 12},
    {71, 22, 10, 4},
    {71, 26, 9, 8},
    {71, 27, 8, 24},
    {71, 30, 6, 12},
    {71, 31, 5, 12},
    {75, 26, 10, 4},
    {75, 30, 8, 12},
    {75, 33, 5, 4},
    {61, 23, 7, 6},
    {61, 24, 6, 8},
    {61, 25, 5, 1},
    {67, 23, 9, 12},
    {67, 27, 7, 12},
    {67, 24, 8, 8},
    {67, 28, 6, 24},
    {67, 29, 5, 4},
    {72, 23, 10, 6},
    {72, 27, 9, 24},
    {72, 30, 7, 6},
    {72, 28, 8, 24},
    {72, 31, 6, 24},
    {72, 32, 5, 6},
    {76, 27, 10, 12},
    {76, 30, 9, 12},
    {76, 31, 8, 24},
    {76, 33, 6, 8},
    {76, 34, 5, 4},
    {79, 30, 10, 6},
    {79, 33, 8, 8},
    {79, 35, 5, 1},
    {62, 24, 7, 4},
    {62, 25, 6, 2},
    {68, 24, 9, 8},
    {68, 28, 7, 12},
    {68, 25, 8, 2},
    {68, 29, 6, 8},
    {73, 24, 10, 4},
    {73, 28, 9, 24},
    {73, 31, 7, 12},
    {73, 29, 8, 8},
    {73, 32, 6, 12},
    {77, 28, 10, 12},
    {77, 31, 9, 24},
    {77, 33, 7, 4},
    {77, 32, 8, 12},
    {77, 34, 6, 8},
    {80, 31, 10, 12},
    {80, 33, 9, 8},
    {80, 34, 8, 8},
    {80, 35, 6, 2},
    {82, 33, 10, 4},
    {82, 35, 8, 2},
    {63, 25, 7, 1},
    {69, 25, 9, 2},
    {69, 29, 7, 4},
    {74, 25, 10, 1},
    {74, 29, 9, 8},
    {74, 32, 7, 6},
    {78, 29, 10, 4},
    {78, 32, 9, 12},
    {78, 34, 7, 4},
    {81, 32, 10, 6},
    {81, 34, 9, 8},
    {81, 35, 7, 1},
    {83, 34, 10, 4},
    {83, 35, 9, 2},
    {84, 35, 10, 1},
    {36, 36, 5, 1},
    {36, 37, 6, 1},
    {36, 42, 8, 1},
    {37, 36, 6, 1},
    {37, 37, 5, 4},
    {37, 37, 7, 1},
    {37, 42, 9, 1},
    {37, 38, 6, 4},
    {37, 43, 8, 4},
    {42, 36, 8, 1},
    {42, 37, 9, 1},
    {42, 42, 5, 4},
    {42, 42, 10, 1},
    {42, 43, 6, 4},
    {42, 47, 8, 4},
    {38, 37, 6, 4},
    {38, 38, 5, 6},
    {38, 38, 7, 4},
    {38, 43, 9, 4},
    {38, 39, 6, 6},
    {38, 44, 8, 6},
    {43, 37, 8, 4},
    {43, 42, 6, 4},
    {43, 38, 9, 4},
    {43, 43, 5, 12},
    {43, 43, 7, 4},
    {43, 43, 10, 4},
    {43, 47, 9, 4},
    {43, 44, 6, 12},
    {43, 48, 8, 12},
    {47, 42, 8, 4},
    {47, 43, 9, 4},
    {47, 47, 5, 6},
    {47, 47, 10, 4},
    {47, 48, 6, 6},
    {47, 51, 8, 6},
    {39, 38, 6, 6},
    {39, 39, 5, 4},
    {39, 39, 7, 6},
    {39, 44, 9, 6},
    {39, 40, 6, 4},
    {39, 45, 8, 4},
    {44, 38, 8, 6},
    {44, 43, 6, 12},
    {44, 39, 9, 6},
    {44, 44, 5, 12},
    {44, 44, 7, 12},
    {44, 44, 10, 6},
    {44, 48, 9, 12},
    {44, 45, 6, 12},
    {44, 49, 8, 12},
    {48, 43, 8, 12},
    {48, 47, 6, 6},
    {48, 44, 9, 12},
    {48, 48, 5, 12},
    {48, 48, 7, 6},
    {48, 48, 10, 12},
    {48, 51, 9, 6},
    {48, 49, 6, 12},
    {48, 52, 8, 12},
    {51, 47, 8, 6},
    {51, 48, 9, 6},
    {51, 51, 5, 4},
    {51, 51, 10, 6},
    {51, 52, 6, 4},
    {51, 54, 8, 4},
    {40, 39, 6, 4},
    {40, 40, 5, 1},
    {40, 40, 7, 4},
    {40, 45, 9, 4},
    {40, 41, 6, 1},
    {40, 46, 8, 1},
    {45, 39, 8, 4},
    {45, 44, 6, 12},
    {45, 40, 9, 4},
    {45, 45, 5, 4},
    {45, 45, 7, 12},
    {45, 45, 10, 4},
    {45, 49, 9, 12},
    {45, 46, 6, 4},
    {45, 50, 8, 4},
    {49, 44, 8, 12},
    {49, 48, 6, 12},
    {49, 45, 9, 12},
    {49, 49, 5, 6},
    {49, 49, 7, 12},
    {49, 49, 10, 12},
    {49, 52, 9, 12},
    {49, 50, 6, 6},
    {49, 53, 8, 6},
    {52, 48, 8, 12},
    {52, 51, 6, 4},
    {52, 49, 9, 12},
    {52, 52, 5, 4},
    {52, 52, 7, 4},
    {52, 52, 10, 12},
    {52, 54, 9, 4},
    {52, 53, 6, 4},
    {52, 55, 8, 4},
    {54, 51, 8, 4},
    {54, 52, 9, 4},
    {54, 54, 5, 1},
    {54, 54, 10, 4},
    {54, 55, 6, 1},
    {54, 56, 8, 1},
    {41, 40, 6, 1},
    {41, 41, 7, 1},
    {41, 46, 9, 1},
    {46, 40, 8, 1},
    {46, 45, 6, 4},
    {46, 41, 9, 1},
    {46, 46, 7, 4},
    {46, 46, 10, 1},
    {46, 50, 9, 4},
    {50, 45, 8, 4},
    {50, 49, 6, 6},
    {50, 46, 9, 4},
    {50, 50, 7, 6},
    {50, 50, 10, 4},
    {50, 53, 9, 6},
    {53, 49, 8, 6},
    {53, 52, 6, 4},
    {53, 50, 9, 6},
    {53, 53, 7, 4},
    {53, 53, 10, 6},
    {53, 55, 9, 4},
    {55, 52, 8, 4},
    {55, 54, 6, 1},
    {55, 53, 9, 4},
    {55, 55, 7, 1},
    {55, 55, 10, 4},
    {55, 56, 9, 1},
    {56, 54, 8, 1},
    {56, 55, 9, 1},
    {56, 56, 10, 1},
    {57, 11, 11, 1},
    {58, 11, 12, 3},
    {58, 12, 11, 3},
    {64, 11, 15, 3},
    {64, 15, 11, 3},
    {59, 11, 13, 3},
    {59, 12, 12, 9},
    {59, 13, 11, 3},
    {65, 11, 16, 6},
    {65, 12, 15, 9},
    {65, 15, 12, 9},
    {65, 16, 11, 6},
    {70, 11, 18, 3},
    {70, 15, 15, 9},
    {70, 18, 11, 3},
    {60, 11, 14, 1},
    {60, 12, 13, 9},
    {60, 13, 12, 9},
    {60, 14, 11, 1},
    {66, 11, 17, 3},
    {66, 12, 16, 18},
    {66, 15, 13, 9},
    {66, 13, 15, 9},
    {66, 16, 12, 18},
    {66, 17, 11, 3},
    {71, 11, 19, 3},
    {71, 12, 18, 9},
    {71, 15, 16, 18},
    {71, 16, 15, 18},
    {71, 18, 12, 9},
    {71, 19, 11, 3},
    {75, 11, 20, 1},
    {75, 15, 18, 9},
    {75, 18, 15, 9},
    {75, 20, 11, 1},
    {61, 12, 14, 3},
    {61, 13, 13, 9},
    {61, 14, 12, 3},
    {67, 12, 17, 9},
    {67, 15, 14, 3},
    {67, 13, 16, 18},
    {67, 16, 13, 18},
    {67, 14, 15, 3},
    {67, 17, 12, 9},
    {72, 12, 19, 9},
    {72, 15, 17, 9},
    {72, 13, 18, 9},
    {72, 16, 16, 36},
    {72, 18, 13, 9},
    {72, 17, 15, 9},
    {72, 19, 12, 9},
    {76, 12, 20, 3},
    {76, 15, 19, 9},
    {76, 16, 18, 18},
    {76, 18, 16, 18},
    {76, 19, 15, 9},
    {76, 20, 12, 3},
    {79, 15, 20, 3},
    {79, 18, 18, 9},
    {79, 20, 15, 3},
    {62, 13, 14, 3},
    {62, 14, 13, 3},
    {68, 13, 17, 9},
    {68, 16, 14, 6},
    {68, 14, 16, 6},
    {68, 17, 13, 9},
    {73, 13, 19, 9},
    {73, 16, 17, 18},
    {73, 18, 14, 3},
    {73, 14, 18, 3},
    {73, 17, 16, 18},
    {73, 19, 13, 9},
    {77, 13, 20, 3},
    {77, 16, 19, 18},
    {77, 18, 17, 9},
    {77, 17, 18, 9},
    {77, 19, 16, 18},
    {77, 20, 13, 3},
    {80, 16, 20, 6},
    {80, 18, 19, 9},
    {80, 19, 18, 9},
    {80, 20, 16, 6},
    {82, 18, 20, 3},
    {82, 20, 18, 3},
    {63, 14, 14, 1},
    {69, 14, 17, 3},
    {69, 17, 14, 3},
    {74, 14, 19, 3},
    {74, 17, 17, 9},
    {74, 19, 14, 3},
    {78, 14, 20, 1},
    {78, 17, 19, 9},
    {78, 19, 17, 9},
    {78, 20, 14, 1},
    {81, 17, 20, 3},
    {81, 19, 19, 9},
    {81, 20, 17, 3},
    {83, 19, 20, 3},
    {83, 20, 19, 3},
    {84, 20, 20, 1},
    {36, 21, 11, 1},
    {36, 22, 12, 1},
    {36, 26, 15, 1},
    {37, 21, 12, 2},
    {37, 22, 11, 3},
    {37, 22, 13, 2},
    {37, 26, 16, 2},
    {37, 23, 12, 3},
    {37, 27, 15, 3},
    {42, 21, 15, 2},
    {42, 22, 16, 2},
    {42, 26, 11, 3},
    {42, 26, 18, 2},
    {42, 27, 12, 3},
    {42, 30, 15, 3},
    {38, 21, 13, 1},
    {38, 22, 12, 6},
    {38, 22, 14, 1},
    {38, 26, 17, 1},
    {38, 23, 11, 3},
    {38, 23, 13, 6},
    {38, 27, 16, 6},
    {38, 24, 12, 3},
    {38, 28, 15, 3},
    {43, 21, 16, 2},
    {43, 22, 15, 6},
    {43, 22, 17, 2},
    {43, 26, 12, 6},
    {43, 26, 19, 2},
    {43, 23, 16, 6},
    {43, 27, 11, 6},
    {43, 27, 13, 6},
    {43, 27, 18, 6},
    {43, 30, 16, 6},
    {43, 28, 12, 6},
    {43, 31, 15, 6},
    {47, 21, 18, 1},
    {47, 22, 19, 1},
    {47, 26, 15, 6},
    {47, 26, 20, 1},
    {47, 27, 16, 6},
    {47, 30, 11, 3},
    {47, 30, 18, 6},
    {47, 31, 12, 3},
    {47, 33, 15, 3},
    {39, 22, 13, 3},
    {39, 23, 12, 6},
    {39, 23, 14, 3},
    {39, 27, 17, 3},
    {39, 24, 11, 1},
    {39, 24, 13, 6},
    {39, 28, 16, 6},
    {39, 25, 12, 1},
    {39, 29, 15, 1},
    {44, 22, 16, 6},
    {44, 26, 13, 3},
    {44, 23, 15, 6},
    {44, 23, 17, 6},
    {44, 27, 12, 12},
    {44, 27, 14, 3},
    {44, 27, 19, 6},
    {44, 30, 17, 3},
    {44, 24, 16, 6},
    {44, 28, 11, 3},
    {44, 28, 13, 12},
    {44, 28, 18, 6},
    {44, 31, 16, 12},
    {44, 29, 12, 3},
    {44, 32, 15, 3},
    {48, 22, 18, 3},
    {48, 26, 16, 6},
    {48, 23, 19, 3},
    {48, 27, 15, 12},
    {48, 27, 17, 6},
    {48, 27, 20, 3},
    {48, 30, 12, 6},
    {48, 30, 19, 6},
    {48, 28, 16, 12},
    {48, 31, 11, 3},
    {48, 31, 13, 6},
    {48, 31, 18, 12},
    {48, 33, 16, 6},
    {48, 32, 12, 3},
    {48, 34, 15, 3},
    {51, 26, 18, 3},
    {51, 27, 19, 3},
    {51, 30, 15, 6},
    {51, 30, 20, 3},
    {51, 31, 16, 6},
    {51, 33, 11, 1},
    {51, 33, 18, 6},
    {51, 34, 12, 1},
    {51, 35, 15, 1},
    {40, 23, 13, 3},
    {40, 24, 12, 2},
    {40, 24, 14, 3},
    {40, 28, 17, 3},
    {40, 25, 13, 2},
    {40, 29, 16, 2},
    {45, 23, 16, 6},
    {45, 27, 13, 6},
    {45, 24, 15, 2},
    {45, 24, 17, 6},
    {45, 28, 12, 6},
    {45, 28, 14, 6},
    {45, 28, 19, 6},
    {45, 31, 17, 6},
    {45, 25, 16, 2},
    {45, 29, 13, 6},
    {45, 29, 18, 2},
    {45, 32, 16, 6},
    {49, 23, 18, 3},
    {49, 27, 16, 12},
    {49, 30, 13, 3},
    {49, 24, 19, 3},
    {49, 28, 15, 6},
    {49, 28, 17, 12},
    {49, 28, 20, 3},
    {49, 31, 12, 6},
    {49, 31, 14, 3},
    {49, 31, 19, 12},
    {49, 33, 17, 3},
    {49, 29, 16, 6},
    {49, 32, 13, 6},
    {49, 32, 18, 6},
    {49, 34, 16, 6},
    {52, 27, 18, 6},
    {52, 30, 16, 6},
    {52, 28, 19, 6},
    {52, 31, 15, 6},
    {52, 31, 17, 6},
    {52, 31, 20, 6},
    {52, 33, 12, 2},
    {52, 33, 19, 6},
    {52, 32, 16, 6},
    {52, 34, 13, 2},
    {52, 34, 18, 6},
    {52, 35, 16, 2},
    {54, 30, 18, 3},
    {54, 31, 19, 3},
    {54, 33, 15, 2},
    {54, 33, 20, 3},
    {54, 34, 16, 2},
    {54, 35, 18, 2},
    {41, 24, 13, 1},
    {41, 25, 14, 1},
    {41, 29, 17, 1},
    {46, 24, 16, 2},
    {46, 28, 13, 3},
    {46, 25, 17, 2},
    {46, 29, 14, 3},
    {46, 29, 19, 2},
    {46, 32, 17, 3},
    {50, 24, 18, 1},
    {50, 28, 16, 6},
    {50, 31, 13, 3},
    {50, 25, 19, 1},
    {50, 29, 17, 6},
    {50, 29, 20, 1},
    {50, 32, 14, 3},
    {50, 32, 19, 6},
    {50, 34, 17, 3},
    {53, 28, 18, 3},
    {53, 31, 16, 6},
    {53, 33, 13, 1},
    {53, 29, 19, 3},
    {53, 32, 17, 6},
    {53, 32, 20, 3},
    {53, 34, 14, 1},
    {53, 34, 19, 6},
    {53, 35, 17, 1},
    {55, 31, 18, 3},
    {55, 33, 16, 2},
    {55, 32, 19, 3},
    {55, 34, 17, 2},
    {55, 34, 20, 3},
    {55, 35, 19, 2},
    {56, 33, 18, 1},
    {56, 34, 19, 1},
    {56, 35, 20, 1},
    {21, 21, 21, 1},
    {21, 22, 22, 2},
    {21, 26, 26, 2},
    {21, 23, 23, 1},
    {21, 27, 27, 2},
    {21, 30, 30, 1},
    {22, 21, 22, 2},
    {22, 22, 21, 2},
    {22, 22, 23, 4},
    {22, 26, 27, 4},
    {22, 23, 22, 4},
    {22, 23, 24, 2},
    {22, 27, 26, 4},
    {22, 27, 28, 4},
    {22, 30, 31, 2},
    {22, 24, 23, 2},
    {22, 28, 27, 4},
    {22, 31, 30, 2},
    {26, 21, 26, 2},
    {26, 22, 27, 4},
    {26, 26, 21, 2},
    {26, 26, 30, 4},
    {26, 23, 28, 2},
    {26, 27, 22, 4},
    {26, 27, 31, 4},
    {26, 30, 26, 4},
    {26, 30, 33, 2},
    {26, 28, 23, 2},
    {26, 31, 27, 4},
    {26, 33, 30, 2},
    {23, 21, 23, 1},
    {23, 22, 22, 4},
    {23, 22, 24, 2},
    {23, 26, 28, 2},
    {23, 23, 21, 1},
    {23, 23, 23, 8},
    {23, 23, 25, 1},
    {23, 27, 27, 8},
    {23, 27, 29, 2},
    {23, 30, 32, 1},
    {23, 24, 22, 2},
    {23, 24, 24, 4},
    {23, 28, 26, 2},
    {23, 28, 28, 8},
    {23, 31, 31, 4},
    {23, 25, 23, 1},
    {23, 29, 27, 2},
    {23, 32, 30, 1},
    {27, 21, 27, 2},
    {27, 22, 26, 4},
    {27, 22, 28, 4},
    {27, 26, 22, 4},
    {27, 26, 31, 4},
    {27, 23, 27, 8},
    {27, 23, 29, 2},
    {27, 27, 21, 2},
    {27, 27, 23, 8},
    {27, 27, 30, 8},
    {27, 27, 32, 4},
    {27, 30, 27, 8},
    {27, 30, 34, 2},
    {27, 24, 28, 4},
    {27, 28, 22, 4},
    {27, 28, 24, 4},
    {27, 28, 31, 8},
    {27, 31, 26, 4},
    {27, 31, 28, 8},
    {27, 31, 33, 4},
    {27, 33, 31, 4},
    {27, 29, 23, 2},
    {27, 32, 27, 4},
    {27, 34, 30, 2},
    {30, 21, 30, 1},
    {30, 22, 31, 2},
    {30, 26, 26, 4},
    {30, 26, 33, 2},
    {30, 23, 32, 1},
    {30, 27, 27, 8},
    {30, 27, 34, 2},
    {30, 30, 21, 1},
    {30, 30, 30, 8},
    {30, 30, 35, 1},
    {30, 28, 28, 4},
    {30, 31, 22, 2},
    {30, 31, 31, 8},
    {30, 33, 26, 2},
    {30, 33, 33, 4},
    {30, 32, 23, 1},
    {30, 34, 27, 2},
    {30, 35, 30, 1},
    {24, 22, 23, 2},
    {24, 23, 22, 2},
    {24, 23, 24, 4},
    {24, 27, 28, 4},
    {24, 24, 23, 4},
    {24, 24, 25, 2},
    {24, 28, 27, 4},
    {24, 28, 29, 4},
    {24, 31, 32, 2},
    {24, 25, 24, 2},
    {24, 29, 28, 4},
    {24, 32, 31, 2},
    {28, 22, 27, 4},
    {28, 26, 23, 2},
    {28, 23, 26, 2},
    {28, 23, 28, 8},
    {28, 27, 22, 4},
    {28, 27, 24, 4},
    {28, 27, 31, 8},
    {28, 30, 28, 4},
    {28, 24, 27, 4},
    {28, 24, 29, 4},
    {28, 28, 23, 8},
    {28, 28, 30, 4},
    {28, 28, 25, 2},
    {28, 28, 32, 8},
    {28, 31, 27, 8},
    {28, 31, 29, 4},
    {28, 31, 34, 4},
    {28, 33, 32, 2},
    {28, 25, 28, 2},
    {28, 29, 24, 4},
    {28, 29, 31, 4},
    {28, 32, 28, 8},
    {28, 32, 33, 2},
    {28, 34, 31, 4},
    {31, 22, 30, 2},
    {31, 26, 27, 4},
    {31, 23, 31, 4},
    {31, 27, 26, 4},
    {31, 27, 28, 8},
    {31, 27, 33, 4},
    {31, 30, 22, 2},
    {31, 30, 31, 8},
    {31, 24, 32, 2},
    {31, 28, 27, 8},
    {31, 28, 29, 4},
    {31, 28, 34, 4},
    {31, 31, 23, 4},
    {31, 31, 30, 8},
    {31, 31, 32, 8},
    {31, 31, 35, 2},
    {31, 33, 27, 4},
    {31, 33, 34, 4},
    {31, 29, 28, 4},
    {31, 32, 24, 2},
    {31, 32, 31, 8},
    {31, 34, 28, 4},
    {31, 34, 33, 4},
    {31, 35, 31, 2},
    {33, 26, 30, 2},
    {33, 27, 31, 4},
    {33, 30, 26, 2},
    {33, 30, 33, 4},
    {33, 28, 32, 2},
    {33, 31, 27, 4},
    {33, 31, 34, 4},
    {33, 33, 30, 4},
    {33, 33, 35, 2},
    {33, 32, 28, 2},
    {33, 34, 31, 4},
    {33, 35, 33, 2},
    {25, 23, 23, 1},
    {25, 24, 24, 2},
    {25, 28, 28, 2},
    {25, 25, 25, 1},
    {25, 29, 29, 2},
    {25, 32, 32, 1},
    {29, 23, 27, 2},
    {29, 27, 23, 2},
    {29, 24, 28, 4},
    {29, 28, 24, 4},
    {29, 28, 31, 4},
    {29, 31, 28, 4},
    {29, 25, 29, 2},
    {29, 29, 25, 2},
    {29, 29, 32, 4},
    {29, 32, 29, 4},
    {29, 32, 34, 2},
    {29, 34, 32, 2},
    {32, 23, 30, 1},
    {32, 27, 27, 4},
    {32, 30, 23, 1},
    {32, 24, 31, 2},
    {32, 28, 28, 8},
    {32, 28, 33, 2},
    {32, 31, 24, 2},
    {32, 31, 31, 8},
    {32, 33, 28, 2},
    {32, 25, 32, 1},
    {32, 29, 29, 4},
    {32, 29, 34, 2},
    {32, 32, 25, 1},
    {32, 32, 32, 8},
    {32, 32, 35, 1},
    {32, 34, 29, 2},
    {32, 34, 34, 4},
    {32, 35, 32, 1},
    {34, 27, 30, 2},
    {34, 30, 27, 2},
    {34, 28, 31, 4},
    {34, 31, 28, 4},
    {34, 31, 33, 4},
    {34, 33, 31, 4},
    {34, 29, 32, 2},
    {34, 32, 29, 2},
    {34, 32, 34, 4},
    {34, 34, 32, 4},
    {34, 34, 35, 2},
    {34, 35, 34, 2},
    {35, 30, 30, 1},
    {35, 31, 31, 2},
    {35, 33, 33, 2},
    {35, 32, 32, 1},
    {35, 34, 34, 2},
    {35, 35, 35, 1}
  };

  int nb4 = nb[Pa];
  for (int i = 0; i<= nb4; i++)
    ns4[i] = ns[i];

  int nmax = ns4[nb4];
  for (int i=0; i<nmax; i++) {
    pb4[i]        = poly[i][0]-1;
    pb4[i+nmax]   = poly[i][1]-1;
    pb4[i+2*nmax] = poly[i][2]-1;
    pc4[i]        = poly[i][3];
  }
}
