// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/*  ----------------------------------------------------------------------
   Contributing authors: Christopher Barrett (MSU) barrett@me.msstate.edu
                              Doyl Dickel (MSU) doyl@me.msstate.edu
    ----------------------------------------------------------------------*/
/*
“The research described and the resulting data presented herein, unless
otherwise noted, was funded under PE 0602784A, Project T53 "Military
Engineering Applied Research", Task 002 under Contract No. W56HZV-17-C-0095,
managed by the U.S. Army Combat Capabilities Development Command (CCDC) and
the Engineer Research and Development Center (ERDC).  The work described in
this document was conducted at CAVS, MSU.  Permission was granted by ERDC
to publish this information. Any opinions, findings and conclusions or
recommendations expressed in this material are those of the author(s) and
do not necessarily reflect the views of the United States Army.​”

DISTRIBUTION A. Approved for public release; distribution unlimited. OPSEC#4918
 */

#include "pair_rann.h"

#include "atom.h"
#include "citeme.h"
#include "error.h"
#include "math_special.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "tokenizer.h"
#include "update.h"

#include <cmath>
#include <cstring>

#include "rann_activation_linear.h"
#include "rann_activation_sig_i.h"
#include "rann_fingerprint_bond.h"
#include "rann_fingerprint_bondscreened.h"
#include "rann_fingerprint_bondscreenedspin.h"
#include "rann_fingerprint_bondspin.h"
#include "rann_fingerprint_radial.h"
#include "rann_fingerprint_radialscreened.h"
#include "rann_fingerprint_radialscreenedspin.h"
#include "rann_fingerprint_radialspin.h"

#define MAXLINE 1024

using namespace LAMMPS_NS;

static const char cite_ml_rann_package[] =
  "ML-RANN package: doi:10.1016/j.commatsci.2020.110207\n\n"
  "@Article{Nitol2021,\n"
  " author = {Nitol, Mashroor S and Dickel, Doyl E and Barrett, Christopher D},\n"
  " title = {Artificial Neural Network Potential for Pure Zinc},\n"
  " journal = {Computational Materials Science},\n"
  " year =    2021,\n"
  " volume =  188,\n"
  " pages =   110207\n"
  "}\n\n";


PairRANN::PairRANN(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_ml_rann_package);

  //initialize ints and bools
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  allocated = 0;
  nelements = -1;
  nelementsp = -1;
  comm_forward = 0;
  comm_reverse = 0;
  res = 10000;
  cutmax = 0;
  dospin = false;
  memguess = 0;
  nmax1 = 0;
  nmax2 = 0;
  fmax = 0;
  fnmax = 0;
  //at least one of the following two flags will change during fingerprint definition:
  doscreen = false;
  allscreen = true;

  //null init for arrays with sizes not yet determined.
  elements = nullptr;
  mass = nullptr;
  elementsp = nullptr;
  map  = nullptr;
  fingerprintcount = nullptr;
  fingerprintlength = nullptr;
  fingerprintperelement = nullptr;
  screening_min = nullptr;
  screening_max = nullptr;
  weightdefined = nullptr;
  biasdefined = nullptr;
  xn = nullptr;
  yn = nullptr;
  zn = nullptr;
  Sik = nullptr;
  dSikx = nullptr;
  dSiky = nullptr;
  dSikz = nullptr;
  dSijkx = nullptr;
  dSijky = nullptr;
  dSijkz = nullptr;
  sx = nullptr;
  sy = nullptr;
  sz = nullptr;
  dSijkxc = nullptr;
  dSijkyc = nullptr;
  dSijkzc = nullptr;
  dfeaturesx = nullptr;
  dfeaturesy = nullptr;
  dfeaturesz = nullptr;
  features = nullptr;
  layer = nullptr;
  sum = nullptr;
  sum1 = nullptr;
  dlayerx = nullptr;
  dlayery = nullptr;
  dlayerz = nullptr;
  dlayersumx = nullptr;
  dlayersumy = nullptr;
  dlayersumz = nullptr;
  dsx = nullptr;
  dsy = nullptr;
  dsz = nullptr;
  dssumx = nullptr;
  dssumy = nullptr;
  dssumz = nullptr;
  tn = nullptr;
  jl = nullptr;
  Bij = nullptr;
  sims = nullptr;
  net = nullptr;
  activation = nullptr;
  fingerprints = nullptr;
}

PairRANN::~PairRANN()
{
  deallocate();
}

void PairRANN::deallocate()
{
  //clear memory
  delete[] mass;
  for (int i=0;i<nelements;i++) {delete [] elements[i];}
  delete[] elements;
  for (int i=0;i<nelementsp;i++) {delete [] elementsp[i];}
  delete[] elementsp;
  for (int i=0;i<=nelements;i++) {
    if (net[i].layers>0) {
      for (int j=0;j<net[i].layers-1;j++) {
        delete[] net[i].Weights[j];
        delete[] net[i].Biases[j];
        delete activation[i][j];
      }
      delete[] activation[i];
      delete[] net[i].Weights;
      delete[] net[i].Biases;
      delete[] weightdefined[i];
      delete[] biasdefined[i];
    }
    delete[] net[i].dimensions;
  }
  delete[] net;
  delete[] map;
  for (int i=0;i<nelementsp;i++) {
    if (fingerprintlength[i]>0) {
      for (int j=0;j<fingerprintperelement[i];j++) {
        delete fingerprints[i][j];
      }
      delete[] fingerprints[i];
    }
  }
  delete[] fingerprints;
  delete[] activation;
  delete[] weightdefined;
  delete[] biasdefined;
  delete[] fingerprintcount;
  delete[] fingerprintperelement;
  delete[] fingerprintlength;
  delete[] screening_min;
  delete[] screening_max;
  memory->destroy(xn);
  memory->destroy(yn);
  memory->destroy(zn);
  memory->destroy(tn);
  memory->destroy(jl);
  memory->destroy(features);
  memory->destroy(dfeaturesx);
  memory->destroy(dfeaturesy);
  memory->destroy(dfeaturesz);
  memory->destroy(layer);
  memory->destroy(sum);
  memory->destroy(sum1);
  memory->destroy(dlayerx);
  memory->destroy(dlayery);
  memory->destroy(dlayerz);
  memory->destroy(dlayersumx);
  memory->destroy(dlayersumy);
  memory->destroy(dlayersumz);
  memory->destroy(Sik);
  memory->destroy(Bij);
  memory->destroy(dSikx);
  memory->destroy(dSiky);
  memory->destroy(dSikz);
  memory->destroy(dSijkx);
  memory->destroy(dSijky);
  memory->destroy(dSijkz);
  memory->destroy(dSijkxc);
  memory->destroy(dSijkyc);
  memory->destroy(dSijkzc);
  memory->destroy(sx);
  memory->destroy(sy);
  memory->destroy(sz);
  memory->destroy(dsx);
  memory->destroy(dsy);
  memory->destroy(dsz);
  memory->destroy(dssumx);
  memory->destroy(dssumy);
  memory->destroy(dssumz);
  memory->destroy(setflag);
  memory->destroy(cutsq);
}

void PairRANN::allocate(const std::vector<std::string> &elementwords)
{
  int i,n;
  n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  cutmax = 0;
  nmax1 = 100;
  nmax2 = 20;
  fmax = 0;
  fnmax = 0;
  nelementsp=nelements+1;
  //initialize arrays
  elements = new char *[nelements];
  elementsp = new char *[nelementsp];//elements + 'all'
  mass = new double[nelements];
  net = new NNarchitecture[nelementsp];
  weightdefined = new bool*[nelementsp];
  biasdefined = new bool *[nelementsp];
  activation = new RANN::Activation**[nelementsp];
  fingerprints = new RANN::Fingerprint**[nelementsp];
  fingerprintlength = new int[nelementsp];
  fingerprintperelement = new int[nelementsp];
  fingerprintcount = new int[nelementsp];
  screening_min = new double[nelements*nelements*nelements];
  screening_max = new double[nelements*nelements*nelements];
  for (i=0;i<nelements;i++) {
    for (int j =0;j<nelements;j++) {
      for (int k=0;k<nelements;k++) {
        screening_min[i*nelements*nelements+j*nelements+k] = 0.8;//default values. Custom values may be read from potential file later.
        screening_max[i*nelements*nelements+j*nelements+k] = 2.8;//default values. Custom values may be read from potential file later.
      }
    }
  }
  for (i=0;i<=nelements;i++) {
    fingerprintlength[i]=0;
    fingerprintperelement[i] = -1;
    fingerprintcount[i] = 0;
    if (i<nelements) {
      mass[i]=-1.0;
      elements[i]= utils::strdup(elementwords[i]);
    }
    elementsp[i] = utils::strdup(elementwords[i]);
    net[i].layers = 0;
    net[i].dimensions = new int[1];
    net[i].dimensions[0]=0;
  }
}

void PairRANN::settings(int narg, char ** /*arg*/)
{
  //read pair_style command in input file
  if (narg > 0) error->one(FLERR,"Illegal pair_style command");
}

void PairRANN::coeff(int narg, char **arg)
{
  int i,j;
  deallocate();//clear allocation from any previous coeff
  map = new int[atom->ntypes+1];
  if (narg != 3 + atom->ntypes) error->one(FLERR,"Incorrect args for pair coefficients");
  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0) error->one(FLERR,"Incorrect args for pair coefficients");
  nelements = -1;
  read_file(arg[2]);
  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++) {
      if (strcmp(arg[i],elements[j]) == 0) break;
    }
    if (j < nelements) map[i-2] = j;
    else error->one(FLERR,"No matching element in NN potential file");
  }
  // clear setflag since coeff() called once with I,J = * *
  int n = atom->ntypes;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      setflag[i][j] = 0;
    }
  }
  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j
  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR,i,mass[map[i]]);
        count++;
      }
    }
  }
  if (count == 0) error->one(FLERR,"Incorrect args for pair coefficients");
  for (i=0;i<nelementsp;i++) {
    for (j=0;j<fingerprintperelement[i];j++) {
      fingerprints[i][j]->allocate();
    }
  }
  allocated=1;
}

void PairRANN::read_file(char *filename)
{
  FILE *fp;
  int eof = 0;
  std::string line,line1;
  const int longline = 4096;
  int linenum=0;
  char linetemp[longline];
  std::string strtemp;
  char *ptr;
  std::vector<std::string> linev,line1v;
  fp = utils::open_potential(filename,lmp,nullptr);
  if (fp == nullptr) {error->one(FLERR,"Cannot open RANN potential file");}
  ptr=fgets(linetemp,longline,fp);
  linenum++;
  strtemp=utils::trim_comment(linetemp);
  while (strtemp.empty()) {
          ptr=fgets(linetemp,longline,fp);
          strtemp=utils::trim_comment(linetemp);
          linenum++;
  }
  line=strtemp;
  while (eof == 0) {
    ptr=fgets(linetemp,longline,fp);
    linenum++;
    if (ptr == nullptr) {
      fclose(fp);
      if (check_potential()) {
        error->one(FLERR,"Invalid syntax in potential file, values are inconsistent or missing");
      }
      else {
        eof = 1;
        break;
      }
    }
    strtemp=utils::trim_comment(linetemp);
    while (strtemp.empty()) {
        ptr=fgets(linetemp,longline,fp);
        strtemp=utils::trim_comment(linetemp);
        linenum++;
    }
    line1=linetemp;
    Tokenizer values = Tokenizer(line,": ,\t_\n");
    Tokenizer values1 = Tokenizer(line1,": ,\t_\n");
    linev = values.as_vector();
    line1v = values1.as_vector();
    if (linev[0]=="atomtypes") read_atom_types(line1v,filename,linenum);
    else if (linev[0]=="mass") read_mass(linev,line1v,filename,linenum);
    else if (linev[0]=="fingerprintsperelement") read_fpe(linev,line1v,filename,linenum);
    else if (linev[0]=="fingerprints") read_fingerprints(linev,line1v,filename,linenum);
    else if (linev[0]=="fingerprintconstants") read_fingerprint_constants(linev,line1v,filename,linenum);
    else if (linev[0]=="networklayers") read_network_layers(linev,line1v,filename,linenum);
    else if (linev[0]=="layersize") read_layer_size(linev,line1v,filename,linenum);
    else if (linev[0]=="weight") read_weight(linev,line1v,fp,filename,&linenum);
    else if (linev[0]=="bias") read_bias(linev,line1v,fp,filename,&linenum);
    else if (linev[0]=="activationfunctions") read_activation_functions(linev,line1v,filename,linenum);
    else if (linev[0]=="screening") read_screening(linev,line1v,filename,linenum);
    else if (linev[0]=="calibrationparameters") continue;//information on how the network was trained
    else error->one(FLERR,"Could not understand file syntax: unknown keyword");
    ptr=fgets(linetemp,longline,fp);
    linenum++;
    strtemp=utils::trim_comment(linetemp);
    while (strtemp.empty()) {
        ptr=fgets(linetemp,longline,fp);
        strtemp=utils::trim_comment(linetemp);
        linenum++;
    }
    if (ptr == nullptr) {
      if (check_potential()) {
        error->one(FLERR,"Invalid syntax in potential file, values are inconsistent or missing");
      }
      else {
        eof = 1;
        break;
      }
    }
    line=linetemp;
  }
}

void PairRANN::read_atom_types(std::vector<std::string> line,char *filename,int linenum) {
  int nwords = line.size();
  if (nwords < 1) error->one(filename,linenum,"Incorrect syntax for atom types");
  nelements = nwords;
  line.emplace_back("all");
  allocate(line);
}

void PairRANN::read_mass(const std::vector<std::string> &line1, const std::vector<std::string> &line2, const char *filename,int linenum) {
  if (nelements == -1)error->one(filename,linenum-1,"atom types must be defined before mass in potential file.");
  for (int i=0;i<nelements;i++) {
    if (line1[1].compare(elements[i])==0) {
      mass[i]=utils::numeric(filename,linenum,line2[0],true,lmp);
      return;
    }
  }
  error->one(filename,linenum-1,"mass element not found in atom types.");
}

void PairRANN::read_fpe(std::vector<std::string> line,std::vector<std::string> line1,char *filename,int linenum) {
  int i;
  if (nelements == -1)error->one(filename,linenum-1,"atom types must be defined before fingerprints per element in potential file.");
  for (i=0;i<nelementsp;i++) {
    if (line[1].compare(elementsp[i])==0) {
      fingerprintperelement[i] = utils::inumeric(filename,linenum,line1[0],true,lmp);
      fingerprints[i] = new RANN::Fingerprint *[fingerprintperelement[i]];
      for (int j=0;j<fingerprintperelement[i];j++) {
        fingerprints[i][j]=new RANN::Fingerprint(this);
      }
      return;
    }
  }
  error->one(filename,linenum-1,"fingerprint-per-element element not found in atom types");
}

void PairRANN::read_fingerprints(std::vector<std::string> line,std::vector<std::string> line1,char *filename,int linenum) {
  int nwords1,nwords,i,j,k,i1,*atomtypes;
  bool found;
  nwords1 = line1.size();
  nwords = line.size();
  if (nelements == -1)error->one(filename,linenum-1,"atom types must be defined before fingerprints in potential file.");
  atomtypes = new int[nwords-1];
  for (i=1;i<nwords;i++) {
    found = false;
    for (j=0;j<nelementsp;j++) {
      if (line[i].compare(elementsp[j])==0) {
        atomtypes[i-1]=j;
        found = true;
        break;
      }
    }
    if (!found) {error->one(filename,linenum-1,"fingerprint element not found in atom types");}
  }
  i = atomtypes[0];
  k = 0;
  if (fingerprintperelement[i]==-1) {error->one(filename,linenum-1,"fingerprint per element must be defined before fingerprints");}
  while (k<nwords1) {
    i1 = fingerprintcount[i];
    if (i1>=fingerprintperelement[i]) {error->one(filename,linenum,"more fingerprints defined than fingerprint per element");}
    delete fingerprints[i][i1];
    fingerprints[i][i1] = create_fingerprint(line1[k].c_str());
    if (fingerprints[i][i1]->n_body_type!=nwords-1) {error->one(filename,linenum,"invalid fingerprint for element combination");}
    k++;
    fingerprints[i][i1]->init(atomtypes,utils::inumeric(filename,linenum,line1[k++],true,lmp));
    fingerprintcount[i]++;
  }
  delete[] atomtypes;
}

void PairRANN::read_fingerprint_constants(std::vector<std::string> line,std::vector<std::string> line1,char *filename,int linenum) {
  int i,j,k,i1,*atomtypes;
  bool found;
  int nwords = line.size();
  if (nelements == -1)error->one(filename,linenum-1,"atom types must be defined before fingerprints in potential file.");
  int n_body_type = nwords-4;
  atomtypes = new int[n_body_type];
  for (i=1;i<=n_body_type;i++) {
    found = false;
    for (j=0;j<nelementsp;j++) {
      if (line[i].compare(elementsp[j])==0) {
        atomtypes[i-1]=j;
        found = true;
        break;
      }
    }
    if (!found) {error->one(filename,linenum-1,"fingerprint element not found in atom types");}
  }
  i = atomtypes[0];
  found = false;
  for (k=0;k<fingerprintperelement[i];k++) {
    if (fingerprints[i][k]->empty) {continue;}
    if (n_body_type!=fingerprints[i][k]->n_body_type) {continue;}
    for (j=0;j<n_body_type;j++) {
      if (fingerprints[i][k]->atomtypes[j]!=atomtypes[j]) {break;}
      if (j==n_body_type-1) {
        if (line[nwords-3].compare(fingerprints[i][k]->style)==0 && utils::inumeric(filename,linenum,line[nwords-2],true,lmp)==fingerprints[i][k]->id) {
          found=true;
          i1 = k;
          break;
        }
      }
    }
    if (found) {break;}
  }
  if (!found) {error->one(filename,linenum-1,"cannot define constants for unknown fingerprint");}
  fingerprints[i][i1]->fullydefined=fingerprints[i][i1]->parse_values(line[nwords-1],line1);
  delete[] atomtypes;
}

void PairRANN::read_network_layers(std::vector<std::string> line,std::vector<std::string> line1,char *filename,int linenum) {
  int i,j;
  if (nelements == -1)error->one(filename,linenum-1,"atom types must be defined before network layers in potential file.");
  for (i=0;i<nelements;i++) {
    if (line[1].compare(elements[i])==0) {
      net[i].layers = utils::inumeric(filename,linenum,line1[0],true,lmp);
      if (net[i].layers < 1)error->one(filename,linenum,"invalid number of network layers");
      delete[] net[i].dimensions;
      weightdefined[i] = new bool [net[i].layers];
      biasdefined[i] = new bool [net[i].layers];
      net[i].dimensions = new int[net[i].layers];
      net[i].Weights = new double * [net[i].layers-1];
      net[i].Biases = new double * [net[i].layers-1];
      for (j=0;j<net[i].layers;j++) {
        net[i].dimensions[j]=0;
        weightdefined[i][j] = false;
        biasdefined[i][j] = false;
      }
      activation[i]=new RANN::Activation* [net[i].layers-1];
      for (int j=0;j<net[i].layers-1;j++) {
        activation[i][j]= new RANN::Activation(this);
      }
      return;
    }
  }
  error->one(filename,linenum-1,"network layers element not found in atom types");
}

void PairRANN::read_layer_size(std::vector<std::string> line,std::vector<std::string> line1,char *filename,int linenum) {
  int i;
  for (i=0;i<nelements;i++) {
    if (line[1].compare(elements[i])==0) {
      if (net[i].layers==0)error->one(filename,linenum-1,"networklayers for each atom type must be defined before the corresponding layer sizes.");
      int j = utils::inumeric(filename,linenum,line[2],true,lmp);
      if (j>=net[i].layers || j<0) {error->one(filename,linenum,"invalid layer in layer size definition");};
      net[i].dimensions[j]= utils::inumeric(filename,linenum,line1[0],true,lmp);
      return;
    }
  }
  error->one(filename,linenum-1,"layer size element not found in atom types");
}

void PairRANN::read_weight(std::vector<std::string> line,std::vector<std::string> line1,FILE* fp,char *filename,int *linenum) {
  int i,j,k,l,nwords;
  char *ptr;
  const int longline = 4096;
  char linetemp [longline];
  for (l=0;l<nelements;l++) {
    if (line[1].compare(elements[l])==0) {
      if (net[l].layers==0)error->one(filename,*linenum-1,"networklayers must be defined before weights.");
      i=utils::inumeric(filename,*linenum,line[2],true,lmp);
      if (i>=net[l].layers || i<0)error->one(filename,*linenum-1,"invalid weight layer");
      if (net[l].dimensions[i]==0 || net[l].dimensions[i+1]==0) error->one(filename,*linenum-1,"network layer sizes must be defined before corresponding weight");
      net[l].Weights[i] = new double[net[l].dimensions[i]*net[l].dimensions[i+1]];
      weightdefined[l][i] = true;
      nwords = line1.size();
      if (nwords != net[l].dimensions[i])error->one(filename,*linenum,"invalid weights per line");
      for (k=0;k<net[l].dimensions[i];k++) {
        net[l].Weights[i][k] = utils::numeric(filename,*linenum,line1[k],true,lmp);
      }
      for (j=1;j<net[l].dimensions[i+1];j++) {
        ptr = fgets(linetemp,longline,fp);
        (*linenum)++;
        Tokenizer values1 = Tokenizer(linetemp,": ,\t_\n");
        line1 = values1.as_vector();
        if (ptr==nullptr)error->one(filename,*linenum,"unexpected end of potential file!");
        nwords = line1.size();
        if (nwords != net[l].dimensions[i])error->one(filename,*linenum,"invalid weights per line");
        for (k=0;k<net[l].dimensions[i];k++) {
          net[l].Weights[i][j*net[l].dimensions[i]+k] = utils::numeric(filename,*linenum,line1[k],true,lmp);
        }
      }
      return;
    }
  }
  error->one(filename,*linenum-1,"weight element not found in atom types");
}

void PairRANN::read_bias(std::vector<std::string> line,std::vector<std::string> line1,FILE* fp,char *filename,int *linenum) {
  int i,j,l;
  char linetemp[MAXLINE],*ptr;
  for (l=0;l<nelements;l++) {
    if (line[1].compare(elements[l])==0) {
      if (net[l].layers==0)error->one(filename,*linenum-1,"networklayers must be defined before biases.");
      i=utils::inumeric(filename,*linenum,line[2],true,lmp);
      if (i>=net[l].layers || i<0)error->one(filename,*linenum-1,"invalid bias layer");
      if (net[l].dimensions[i]==0) error->one(filename,*linenum-1,"network layer sizes must be defined before corresponding bias");
      biasdefined[l][i] = true;
      net[l].Biases[i] = new double[net[l].dimensions[i+1]];
      net[l].Biases[i][0] = utils::numeric(filename,*linenum,line1[0],true,lmp);
      for (j=1;j<net[l].dimensions[i+1];j++) {
        ptr=fgets(linetemp,MAXLINE,fp);
        if (ptr==nullptr)error->one(filename,*linenum,"unexpected end of potential file!");
        (*linenum)++;
        Tokenizer values1 = Tokenizer(linetemp,": ,\t_\n");
        line1 = values1.as_vector();
        net[l].Biases[i][j] = utils::numeric(filename,*linenum,line1[0],true,lmp);
      }
      return;
    }
  }
  error->one(filename,*linenum-1,"bias element not found in atom types");
}

void PairRANN::read_activation_functions(std::vector<std::string> line,std::vector<std::string> line1,char *filename,int linenum) {
  int i,l;
  for (l=0;l<nelements;l++) {
    if (line[1].compare(elements[l])==0) {
      if (net[l].layers==0)error->one(filename,linenum-1,"networklayers must be defined before activation functions.");
      i = strtol(line[2].c_str(),nullptr,10);
      if (i>=net[l].layers || i<0)error->one(filename,linenum-1,"invalid activation layer");
      delete activation[l][i];
      activation[l][i]=create_activation(line1[0].c_str());
      return;
    }
  }
  error->one(filename,linenum-1,"activation function element not found in atom types");
}

void PairRANN::read_screening(std::vector<std::string> line,std::vector<std::string> line1,char *filename,int linenum) {
  int i,j,k,*atomtypes;
  bool found;
  int nwords = line.size();
  if (nelements == -1)error->one(filename,linenum-1,"atom types must be defined before fingerprints in potential file.");
  if (nwords!=5)error->one(filename,linenum-1,"invalid screening command");
  int n_body_type = 3;
  atomtypes = new int[n_body_type];
  for (i=1;i<=n_body_type;i++) {
    found = false;
    for (j=0;j<nelementsp;j++) {
      if (line[i].compare(elementsp[j])==0) {
        atomtypes[i-1]=j;
        found = true;
        break;
      }
    }
    if (!found) {error->one(filename,linenum-1,"fingerprint element not found in atom types");}
  }
  i = atomtypes[0];
  j = atomtypes[1];
  k = atomtypes[2];
  int index = i*nelements*nelements+j*nelements+k;
  if (line[4].compare("Cmin")==0)  {
    screening_min[index] = utils::numeric(filename,linenum,line1[0],true,lmp);
  }
  else if (line[4].compare("Cmax")==0) {
    screening_max[index] = utils::numeric(filename,linenum,line1[0],true,lmp);
  }
  else error->one(filename,linenum-1,"unrecognized screening keyword");
  delete[] atomtypes;
}

//Called after finishing reading the potential file to make sure it is complete. True is bad.
//also does the rest of the memory allocation.
bool PairRANN::check_potential() {
  int i,j,k,l;
  if (nelements==-1) {return true;}
  for (i=0;i<=nelements;i++) {
    if (i<nelements) {
      if (mass[i]<0)return true;//uninitialized mass
    }
    if (net[i].layers==0)break;//no definitions for this starting element, not considered an error.
    net[i].maxlayer=0;
    for (j=0;j<net[i].layers;j++) {
      if (net[i].dimensions[j]==0)return true;//incomplete network definition
      if (net[i].dimensions[j]>net[i].maxlayer)net[i].maxlayer = net[i].dimensions[j];
    }
    if (net[i].maxlayer>fnmax) {fnmax = net[i].maxlayer;}
    if (net[i].dimensions[net[i].layers-1]!=1)return true;//output layer must have single neuron (the energy)
    if (net[i].dimensions[0]>fmax)fmax=net[i].dimensions[0];
    for (j=0;j<net[i].layers-1;j++) {
      if (!weightdefined[i][j])return true;//undefined weights
      if (!biasdefined[i][j])return true;//undefined biases
      if (activation[i][j]->empty)return true;//undefined activations
      for (k=0;k<net[i].dimensions[j+1];k++) {
        for (l=0;l<net[i].dimensions[j];l++) {
          if (net[i].Weights[j][k*net[i].dimensions[j]+l]==0)return true;//undefined weights
        }
        if (net[i].Biases[j][k]==0)return true;//undefined biases
      }
    }
    for (j=0;j<fingerprintperelement[i];j++) {
      if (fingerprints[i][j]->fullydefined==false)return true;
      fingerprints[i][j]->startingneuron = fingerprintlength[i];
      fingerprintlength[i] +=fingerprints[i][j]->get_length();
      if (fingerprints[i][j]->rc>cutmax) {cutmax = fingerprints[i][j]->rc;}
    }
    if (net[i].dimensions[0]!=fingerprintlength[i])return true;
  }
  memory->create(xn,nmax1,"pair:xn");
  memory->create(yn,nmax1,"pair:yn");
  memory->create(zn,nmax1,"pair:zn");
  memory->create(tn,nmax1,"pair:tn");
  memory->create(jl,nmax1,"pair:jl");
  memory->create(features,fmax,"pair:features");
  memory->create(dfeaturesx,fmax*nmax2,"pair:dfeaturesx");
  memory->create(dfeaturesy,fmax*nmax2,"pair:dfeaturesy");
  memory->create(dfeaturesz,fmax*nmax2,"pair:dfeaturesz");
  memory->create(layer,fnmax,"pair:layer");
  memory->create(sum,fnmax,"pair:sum");
  memory->create(sum1,fnmax,"pair:sum1");
  memory->create(dlayerx,nmax2,fnmax,"pair:dlayerx");
  memory->create(dlayery,nmax2,fnmax,"pair:dlayery");
  memory->create(dlayerz,nmax2,fnmax,"pair:dlayerz");
  memory->create(dlayersumx,nmax2,fnmax,"pair:dlayersumx");
  memory->create(dlayersumy,nmax2,fnmax,"pair:dlayersumy");
  memory->create(dlayersumz,nmax2,fnmax,"pair:dlayersumz");
  if (doscreen) {
    memory->create(Sik,nmax2,"pair:Sik");
    memory->create(Bij,nmax2,"pair:Bij");
    memory->create(dSikx,nmax2,"pair:dSikx");
    memory->create(dSiky,nmax2,"pair:dSiky");
    memory->create(dSikz,nmax2,"pair:dSikz");
    memory->create(dSijkx,nmax2*nmax2,"pair:dSijkx");
    memory->create(dSijky,nmax2*nmax2,"pair:dSijky");
    memory->create(dSijkz,nmax2*nmax2,"pair:dSijkz");
    memory->create(dSijkxc,nmax2,nmax2,"pair:dSijkxc");
    memory->create(dSijkyc,nmax2,nmax2,"pair:dSijkyc");
    memory->create(dSijkzc,nmax2,nmax2,"pair:dSijkzc");
  }
  if (dospin) {
    memory->create(sx,fmax*nmax2,"pair:sx");
    memory->create(sy,fmax*nmax2,"pair:sy");
    memory->create(sz,fmax*nmax2,"pair:sz");
    memory->create(dsx,nmax2,fnmax,"pair:dsx");
    memory->create(dsy,nmax2,fnmax,"pair:dsy");
    memory->create(dsz,nmax2,fnmax,"pair:dsz");
    memory->create(dssumx,nmax2,fnmax,"pair:dssumx");
    memory->create(dssumy,nmax2,fnmax,"pair:dssumy");
    memory->create(dssumz,nmax2,fnmax,"pair:dssumz");
  }
  return false;//everything looks good
}

void PairRANN::compute(int eflag, int vflag)
{
  //perform force/energy computation_
  if (dospin) {
    if (strcmp(update->unit_style,"metal") != 0)
      error->one(FLERR,"Spin pair styles require metal units");
    if (!atom->sp_flag)
        error->one(FLERR,"Spin pair styles requires atom/spin style");
  }
  ev_init(eflag,vflag);

  // only global virial via fdotr is supported by this pair style

  if (vflag_atom)
    error->all(FLERR,"Pair style rann does not support computing per-atom stress");
  if (vflag && !vflag_fdotr)
    error->all(FLERR,"Pair style rann does not support 'pair_modify nofdotr'");

  int ii,i,j;
  int nn = 0;
  sims = new Simulation[1];
  sims->inum = listfull->inum;
  sims->ilist=listfull->ilist;
  sims->id = listfull->ilist;
  sims->type = atom->type;
  sims->x = atom->x;
  sims->numneigh = listfull->numneigh;
  sims->firstneigh = listfull->firstneigh;
  if (dospin) {
    sims->s = atom->sp;
  }
  int itype,f,jnum,len;
  if (eflag_global) {eng_vdwl=0;eng_coul=0;}
  double energy=0;
  double **force = atom->f;
  double **fm = atom->fm;

  //loop over atoms
  for (ii=0;ii<sims->inum;ii++) {
      i = sims->ilist[ii];
      itype = map[sims->type[i]];
      f = net[itype].dimensions[0];
      jnum = sims->numneigh[i];
      if (jnum>nmax1) {
        nmax1 = jnum;
        memory->grow(xn,nmax1,"pair:xn");
        memory->grow(yn,nmax1,"pair:yn");
        memory->grow(zn,nmax1,"pair:zn");
        memory->grow(tn,nmax1,"pair:tn");
        memory->grow(jl,nmax1,"pair:jl");
      }
      cull_neighbor_list(&jnum,i,0);
      if (jnum>nmax2) {
        nmax2=jnum;
        memory->grow(dfeaturesx,fmax*nmax2,"pair:dfeaturesx");
        memory->grow(dfeaturesy,fmax*nmax2,"pair:dfeaturesy");
        memory->grow(dfeaturesz,fmax*nmax2,"pair:dfeaturesz");
        memory->grow(layer,fnmax,"pair:layer");
        memory->grow(sum,fnmax,"pair:sum");
        memory->grow(sum1,fnmax,"pair:sum1");
        memory->grow(dlayerx,nmax2,fnmax,"pair:dlayerx");
        memory->grow(dlayery,nmax2,fnmax,"pair:dlayery");
        memory->grow(dlayerz,nmax2,fnmax,"pair:dlayerz");
        memory->grow(dlayersumx,nmax2,fnmax,"pair:dlayersumx");
        memory->grow(dlayersumy,nmax2,fnmax,"pair:dlayersumy");
        memory->grow(dlayersumz,nmax2,fnmax,"pair:dlayersumz");
        if (doscreen) {
          memory->grow(Sik,nmax2,"pair:Sik");
          memory->grow(Bij,nmax2,"pair:Bij");
          memory->grow(dSikx,nmax2,"pair:dSikx");
          memory->grow(dSiky,nmax2,"pair:dSiky");
          memory->grow(dSikz,nmax2,"pair:dSikz");
          memory->grow(dSijkx,nmax2*nmax2,"pair:dSijkx");
          memory->grow(dSijky,nmax2*nmax2,"pair:dSijky");
          memory->grow(dSijkz,nmax2*nmax2,"pair:dSijkz");
          memory->destroy(dSijkxc);
          memory->destroy(dSijkyc);
          memory->destroy(dSijkzc);
          memory->create(dSijkxc,nmax2,nmax2,"pair:dSijkxc");
          memory->create(dSijkyc,nmax2,nmax2,"pair:dSijkyc");
          memory->create(dSijkzc,nmax2,nmax2,"pair:dSijkzc");
        }
        if (dospin) {
          memory->grow(sx,fmax*nmax2,"pair:sx");
          memory->grow(sy,fmax*nmax2,"pair:sy");
          memory->grow(sz,fmax*nmax2,"pair:sz");
          memory->grow(dsx,nmax2,fnmax,"pair:dsx");
          memory->grow(dsy,nmax2,fnmax,"pair:dsy");
          memory->grow(dsz,nmax2,fnmax,"pair:dsz");
          memory->grow(dssumx,nmax2,fnmax,"pair:dssumx");
          memory->grow(dssumy,nmax2,fnmax,"pair:dssumy");
          memory->grow(dssumz,nmax2,fnmax,"pair:dssumz");
        }
      }
      for (j=0;j<f;j++) {
        features[j]=0;
      }
      for (j=0;j<f*jnum;j++) {
        dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
      }
      //screening is calculated once for all atoms if any fingerprint uses it.
      if (dospin) {
        for (j=0;j<f*jnum;j++) {
          sx[j]=sy[j]=sz[j]=0;
        }
      }
      if (doscreen) screening(ii,0,jnum-1);
      if (allscreen) screen_neighbor_list(&jnum);
      //do fingerprints for atom type
      len = fingerprintperelement[itype];
      for (j=0;j<len;j++) {
        if      (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
        else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
        else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
        else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
      }
      itype = nelements;
      //do fingerprints for type "all"
      len = fingerprintperelement[itype];
      for (j=0;j<len;j++) {
        if      (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
        else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
        else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
        else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
      }
      //run fingerprints through network
      if (dospin) {
        propagateforwardspin(energy,force,fm,ii,jnum);
      } else {
        propagateforward(energy,force,ii,jnum);
      }
  }
  if (vflag_fdotr) virial_fdotr_compute();
  delete[] sims;
}

void PairRANN::cull_neighbor_list(int* jnum,int i,int sn) {
  int *jlist,j,count,jj,*type,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double **x = sims[sn].x;
  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];
  type = sims[sn].type;
  jlist = sims[sn].firstneigh[i];
  count = 0;
  for (jj=0;jj<jnum[0];jj++) {
    j = jlist[jj];
    j &= NEIGHMASK;
    jtype = map[type[j]];
    delx = xtmp - x[j][0];
    dely = ytmp - x[j][1];
    delz = ztmp - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq>cutmax*cutmax) {
      continue;
    }
    xn[count]=delx;
    yn[count]=dely;
    zn[count]=delz;
    tn[count]=jtype;
    jl[count]=j;
    count++;
  }
  jnum[0]=count+1;
}

void PairRANN::screen_neighbor_list(int *jnum) {
  int jj,kk,count,count1;
  count = 0;
  for (jj=0;jj<jnum[0]-1;jj++) {
    if (Bij[jj]) {
      count1 = 0;
      if (jj!=count) {
      xn[count]=xn[jj];
      yn[count]=yn[jj];
      zn[count]=zn[jj];
      tn[count]=tn[jj];
      jl[count]=jl[jj];
      Sik[count]=Sik[jj];
      dSikx[count]=dSikx[jj];
      dSiky[count]=dSiky[jj];
      dSikz[count]=dSikz[jj];
      }
      for (kk=0;kk<jnum[0]-1;kk++) {
        if (Bij[kk]) {
          dSijkxc[count][count1] = dSijkx[jj*(jnum[0]-1)+kk];
          dSijkyc[count][count1] = dSijky[jj*(jnum[0]-1)+kk];
          dSijkzc[count][count1] = dSijkz[jj*(jnum[0]-1)+kk];
          count1++;
        }
      }
      count++;
    }
  }
  jnum[0]=count+1;
  for (jj=0;jj<count;jj++) {
    Bij[jj]=true;
    for (kk=0;kk<count;kk++) {
      dSijkx[jj*count+kk] = dSijkxc[jj][kk];
      dSijky[jj*count+kk] = dSijkyc[jj][kk];
      dSijkz[jj*count+kk] = dSijkzc[jj][kk];
    }
  }
}


void PairRANN::screening(int ii,int sid,int jnum)
{
  //see Baskes, Materials Chemistry and Physics 50 (1997) 152-1.58
  int i,jj,kk,itype,jtype,ktype;
  double Sijk,Cijk,Cn,Cd,Dij,Dik,Djk,C,dfc,dC;
  PairRANN::Simulation *sim = &sims[sid];
  double delx,dely,delz,rij,delx2,dely2,delz2,rik,delx3,dely3,delz3,rjk;
  i = sim->ilist[ii];
  itype = map[sim->type[i]];
  for (int jj=0;jj<jnum;jj++) {
    Sik[jj]=1;
    Bij[jj]=true;
    dSikx[jj]=0;
    dSiky[jj]=0;
    dSikz[jj]=0;
    for (kk=0;kk<jnum;kk++) {
      dSijkx[jj*jnum+kk]=0;
      dSijky[jj*jnum+kk]=0;
      dSijkz[jj*jnum+kk]=0;
    }
  }
  for (kk=0;kk<jnum;kk++) {//outer sum over k in accordance with source, some others reorder to outer sum over jj
    if (Bij[kk]==false) {continue;}
    ktype = tn[kk];
    delx2 = xn[kk];
    dely2 = yn[kk];
    delz2 = zn[kk];
    rik = delx2*delx2+dely2*dely2+delz2*delz2;
    if (rik>cutmax*cutmax) {
      Bij[kk]= false;
      continue;
    }
    for (jj=0;jj<jnum;jj++) {
      if (jj==kk) {continue;}
      if (Bij[jj]==false) {continue;}
      jtype = tn[jj];
      delx = xn[jj];
      dely = yn[jj];
      delz = zn[jj];
      rij = delx*delx+dely*dely+delz*delz;
      if (rij>cutmax*cutmax) {
        Bij[jj] = false;
        continue;
      }
      delx3 = delx2-delx;
      dely3 = dely2-dely;
      delz3 = delz2-delz;
      rjk = delx3*delx3+dely3*dely3+delz3*delz3;
      if (rik+rjk<=rij) {continue;}//bond angle > 90 degrees
      if (rik+rij<=rjk) {continue;}//bond angle > 90 degrees
      double Cmax = screening_max[itype*nelements*nelements+jtype*nelements+ktype];
      double Cmin = screening_min[itype*nelements*nelements+jtype*nelements+ktype];
      double temp1 = rij-rik+rjk;
      Cn = temp1*temp1-4*rij*rjk;
      //a few expanded computations provided for clarity:
      //Cn = (rij-rik+rjk)*(rij-rik+rjk)-4*rij*rjk;
      //Cd = (rij-rjk)*(rij-rjk)-rik*rik;
      //Cijk = 1+2*(rik*rij+rik*rjk-rik*rik)/(rik*rik-(rij-rjk)*(rij-rjk));
      temp1 = rij-rjk;
      Cd = temp1*temp1-rik*rik;
      Cijk = Cn/Cd;
      C = (Cijk-Cmin)/(Cmax-Cmin);
      if (C>=1) {continue;}
      else if (C<=0) {
        Bij[kk]=false;
        break;
      }
      dC = Cmax-Cmin;
      dC *= dC;
      dC *= dC;
      temp1 = 1-C;
      temp1 *= temp1;
      temp1 *= temp1;
      Sijk = 1-temp1;
      Sijk *= Sijk;
      Dij = 4*rik*(Cn+4*rjk*(rij+rik-rjk))/Cd/Cd;
      Dik = -4*(rij*Cn+rjk*Cn+8*rij*rik*rjk)/Cd/Cd;
      Djk = 4*rik*(Cn+4*rij*(rik-rij+rjk))/Cd/Cd;
      temp1 = Cijk-Cmax;
      double temp2 = temp1*temp1;
      dfc = 8*temp1*temp2/(temp2*temp2-dC);
      Sik[kk] *= Sijk;
      dSijkx[kk*jnum+jj] = dfc*(delx*Dij-delx3*Djk);
      dSikx[kk] += dfc*(delx2*Dik+delx3*Djk);
      dSijky[kk*jnum+jj] = dfc*(dely*Dij-dely3*Djk);
      dSiky[kk] += dfc*(dely2*Dik+dely3*Djk);
      dSijkz[kk*jnum+jj] = dfc*(delz*Dij-delz3*Djk);
      dSikz[kk] += dfc*(delz2*Dik+delz3*Djk);
    }
  }
}


//Called by getproperties. Propagate features and dfeatures through network. Updates force and energy
void PairRANN::propagateforward(double &energy,double **force,int ii,int jnum) {
  int i,j,k,jj,j1,itype,i1;
  int *ilist;
  ilist = listfull->ilist;
  int *type = atom->type;
  i1=ilist[ii];
  itype = map[type[i1]];
  NNarchitecture net1 = net[itype];
  int L = net1.layers-1;
  //energy output with forces from analytical derivatives
  double dsum1;
  int f = net1.dimensions[0];
  for (i=0;i<net1.layers-1;i++) {
    for (j=0;j<net1.dimensions[i+1];j++) {
      //energy forward propagation
      sum[j]=0;
      for (k=0;k<net1.dimensions[i];k++) {
        if (i==0&&j==0) {
          layer[k]=features[k];
        }
        sum[j] += net1.Weights[i][j*net1.dimensions[i]+k]*layer[k];
      }
      sum[j] += net1.Biases[i][j];
      dsum1 = activation[itype][i]->dactivation_function(sum[j]);
      sum[j] = activation[itype][i]->activation_function(sum[j]);
      if (i==L-1) {
        energy = sum[j];
        if (eflag_atom) eatom[i1]=sum[j];
        if (eflag_global) eng_vdwl +=sum[j];
      }
      //force propagation
      for (jj=0;jj<jnum;jj++) {
        dlayersumx[jj][j]=0;
        dlayersumy[jj][j]=0;
        dlayersumz[jj][j]=0;
        for (k=0;k<net1.dimensions[i];k++) {
          if (i==0&&j==0) {
            dlayerx[jj][k]=dfeaturesx[jj*f+k];
            dlayery[jj][k]=dfeaturesy[jj*f+k];
            dlayerz[jj][k]=dfeaturesz[jj*f+k];
          }
          double w1 = net1.Weights[i][j*net1.dimensions[i]+k];
          dlayersumx[jj][j] += w1*dlayerx[jj][k];
          dlayersumy[jj][j] += w1*dlayery[jj][k];
          dlayersumz[jj][j] += w1*dlayerz[jj][k];
        }
        dlayersumx[jj][j] = dsum1*dlayersumx[jj][j];
        dlayersumy[jj][j] = dsum1*dlayersumy[jj][j];
        dlayersumz[jj][j] = dsum1*dlayersumz[jj][j];
        if (i==L-1 && jj < (jnum-1)) {
          j1 = jl[jj];
          force[j1][0]+=dlayersumx[jj][j];
          force[j1][1]+=dlayersumy[jj][j];
          force[j1][2]+=dlayersumz[jj][j];
        }
      }
      if (i==L-1) {
        j1 = ilist[ii];
        jj = jnum-1;
        force[j1][0]+=dlayersumx[jj][j];
        force[j1][1]+=dlayersumy[jj][j];
        force[j1][2]+=dlayersumz[jj][j];
      }
    }
    //update values for next iteration
    for (j=0;j<net1.dimensions[i+1];j++) {
      layer[j]=sum[j];
      for (jj=0;jj<jnum;jj++) {
        dlayerx[jj][j] = dlayersumx[jj][j];
        dlayery[jj][j] = dlayersumy[jj][j];
        dlayerz[jj][j] = dlayersumz[jj][j];
      }
    }
  }
}

//Called by getproperties. Propagate features and dfeatures through network. Updates force and energy
void PairRANN::propagateforwardspin(double &energy,double **force,double **fm,int ii,int jnum) {
  int i,j,k,jj,j1,itype,i1;
  int *ilist;
  ilist = listfull->ilist;
  int *type = atom->type;
  i1=ilist[ii];
  itype = map[type[i1]];
  NNarchitecture net1 = net[itype];
  int L = net1.layers-1;
  //energy output with forces from analytical derivatives
  double dsum1;
  int f = net1.dimensions[0];
  for (i=0;i<net1.layers-1;i++) {
    for (j=0;j<net1.dimensions[i+1];j++) {
      //energy forward propagation
      sum[j]=0;
      for (k=0;k<net1.dimensions[i];k++) {
        if (i==0&&j==0) {
          layer[k]=features[k];
        }
        sum[j] += net1.Weights[i][j*net1.dimensions[i]+k]*layer[k];
      }
      sum[j] += net1.Biases[i][j];
      dsum1 = activation[itype][i]->dactivation_function(sum[j]);
      sum[j] = activation[itype][i]->activation_function(sum[j]);
      if (i==L-1) {
        energy = sum[j];
        if (eflag_atom) eatom[i1]=sum[j];
        if (eflag_global) eng_vdwl +=sum[j];
      }
      //force propagation
      for (jj=0;jj<jnum;jj++) {
        dlayersumx[jj][j]=0;
        dlayersumy[jj][j]=0;
        dlayersumz[jj][j]=0;
        dssumx[jj][j]=0;
        dssumy[jj][j]=0;
        dssumz[jj][j]=0;
        for (k=0;k<net1.dimensions[i];k++) {
          if (i==0&&j==0) {
            dlayerx[jj][k]=dfeaturesx[jj*f+k];
            dlayery[jj][k]=dfeaturesy[jj*f+k];
            dlayerz[jj][k]=dfeaturesz[jj*f+k];
            dsx[jj][k]=-sx[jj*f+k];
            dsy[jj][k]=-sy[jj*f+k];
            dsz[jj][k]=-sz[jj*f+k];
          }
          double w1 = net1.Weights[i][j*net1.dimensions[i]+k];
          dlayersumx[jj][j] += w1*dlayerx[jj][k];
          dlayersumy[jj][j] += w1*dlayery[jj][k];
          dlayersumz[jj][j] += w1*dlayerz[jj][k];
          dssumx[jj][j] += w1*dsx[jj][k];
          dssumy[jj][j] += w1*dsy[jj][k];
          dssumz[jj][j] += w1*dsz[jj][k];
        }
        dlayersumx[jj][j] = dsum1*dlayersumx[jj][j];
        dlayersumy[jj][j] = dsum1*dlayersumy[jj][j];
        dlayersumz[jj][j] = dsum1*dlayersumz[jj][j];
        dssumx[jj][j] *= dsum1;
        dssumy[jj][j] *= dsum1;
        dssumz[jj][j] *= dsum1;
        if (i==L-1 && jj < (jnum-1)) {
          j1 = jl[jj];
          force[j1][0]+=dlayersumx[jj][j];
          force[j1][1]+=dlayersumy[jj][j];
          force[j1][2]+=dlayersumz[jj][j];
          fm[j1][0]+=dssumx[jj][j];
          fm[j1][1]+=dssumy[jj][j];
          fm[j1][2]+=dssumz[jj][j];
        }
      }
      if (i==L-1) {
        j1 = ilist[ii];
        jj = jnum-1;
        force[j1][0]+=dlayersumx[jj][j];
        force[j1][1]+=dlayersumy[jj][j];
        force[j1][2]+=dlayersumz[jj][j];
        fm[j1][0]+=dssumx[jj][j];
        fm[j1][1]+=dssumy[jj][j];
        fm[j1][2]+=dssumz[jj][j];
      }
    }
    //update values for next iteration
    for (j=0;j<net1.dimensions[i+1];j++) {
      layer[j]=sum[j];
      for (jj=0;jj<jnum;jj++) {
        dlayerx[jj][j] = dlayersumx[jj][j];
        dlayery[jj][j] = dlayersumy[jj][j];
        dlayerz[jj][j] = dlayersumz[jj][j];
        dsx[jj][j] = dssumx[jj][j];
        dsy[jj][j] = dssumy[jj][j];
        dsz[jj][j] = dssumz[jj][j];
      }
    }
  }
}

void PairRANN::init_list(int /*which*/, NeighList *ptr)
{
  listfull = ptr;
}

void PairRANN::init_style()
{
  neighbor->add_request(this, NeighConst::REQ_FULL);
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairRANN::init_one(int /*i*/, int /*j*/)
{
  return cutmax;
}

void PairRANN::errorf(const char *file, int line, const char * message) {
  error->one(file,line,message);
}

int PairRANN::factorial(int n) {
  return round(MathSpecial::factorial(n));
}

RANN::Fingerprint *PairRANN::create_fingerprint(const char *style)
{
  if (strcmp(style,"radial")==0) {
    return new RANN::Fingerprint_radial(this);
  }
  else if (strcmp(style,"radialscreened")==0) {
    return new RANN::Fingerprint_radialscreened(this);
  }
  else if (strcmp(style,"radialscreenedspin")==0) {
    return new RANN::Fingerprint_radialscreenedspin(this);
  }
  else if (strcmp(style,"radialspin")==0) {
    return new RANN::Fingerprint_radialspin(this);
  }
  else if (strcmp(style,"bond")==0) {
    return new RANN::Fingerprint_bond(this);
  }
  else if (strcmp(style,"bondscreened")==0) {
    return new RANN::Fingerprint_bondscreened(this);
  }
  else if (strcmp(style,"bondscreenedspin")==0) {
    return new RANN::Fingerprint_bondscreenedspin(this);
  }
  else if (strcmp(style,"bondspin")==0) {
    return new RANN::Fingerprint_bondspin(this);
  }
  error->one(FLERR,"Unknown fingerprint style {}",style);
  return nullptr;
}


RANN::Activation *PairRANN::create_activation(const char *style)
{
  if (strcmp(style,"linear")==0) {
    return new RANN::Activation_linear(this);
  }
  else if (strcmp(style,"sigI")==0) {
    return new RANN::Activation_sigI(this);
  }
  error->one(FLERR,"Unknown activation style {}",style);
  return nullptr;
}

