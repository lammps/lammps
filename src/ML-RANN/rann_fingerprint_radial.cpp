// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#include "rann_fingerprint_radial.h"
#include "pair_rann.h"

#include <cmath>

using namespace LAMMPS_NS::RANN;

Fingerprint_radial::Fingerprint_radial(PairRANN *_pair) : Fingerprint(_pair)
{
  n_body_type = 2;
  dr = 0;
  re = 0;
  rc = 0;
  alpha = new double[1];
  alpha[0] = -1;
  nmax = 0;
  omin = 0;
  id = -1;
  style = "radial";
  atomtypes = new int[n_body_type];
  empty = true;
  fullydefined = false;
  _pair->allscreen = false;
}

Fingerprint_radial::~Fingerprint_radial()
{
  delete[] atomtypes;
  delete[] radialtable;
  delete[] alpha;
  delete[] dfctable;
  delete[] rinvsqrttable;
}

bool Fingerprint_radial::parse_values(std::string constant,std::vector<std::string> line1) {
  int l;
  int nwords=line1.size();
  if (constant.compare("re")==0) {
    re = strtod(line1[0].c_str(),nullptr);
  }
  else if (constant.compare("rc")==0) {
    rc = strtod(line1[0].c_str(),nullptr);
  }
  else if (constant.compare("alpha")==0) {
    delete[] alpha;
    alpha = new double[nwords];
    for (l=0;l<nwords;l++) {
      alpha[l]=strtod(line1[l].c_str(),nullptr);
    }
  }
  else if (constant.compare("dr")==0) {
    dr = strtod(line1[0].c_str(),nullptr);
  }
  else if (constant.compare("n")==0) {
    nmax = strtol(line1[0].c_str(),nullptr,10);
  }
  else if (constant.compare("o")==0) {
    omin = strtol(line1[0].c_str(),nullptr,10);
  }
  else pair->errorf(FLERR,"Undefined value for radial power");
  //code will run with default o=0 if o is never specified. All other values must be defined in potential file.
  if (re!=0 && rc!=0 && alpha!=nullptr && dr!=0 && nmax!=0)return true;
  return false;
}

void Fingerprint_radial::write_values(FILE *fid) {
  int i;
  fprintf(fid,"fingerprintconstants:");
  fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
  for (i=1;i<n_body_type;i++) {
    fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
  }
  fprintf(fid,":%s_%d:re:\n",style,id);
  fprintf(fid,"%f\n",re);
  fprintf(fid,"fingerprintconstants:");
  fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
  for (i=1;i<n_body_type;i++) {
    fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
  }
  fprintf(fid,":%s_%d:rc:\n",style,id);
  fprintf(fid,"%f\n",rc);
  fprintf(fid,"fingerprintconstants:");
  fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
  for (i=1;i<n_body_type;i++) {
    fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
  }
  fprintf(fid,":%s_%d:alpha:\n",style,id);
  for (i=0;i<(nmax-omin+1);i++) {
    fprintf(fid,"%f ",alpha[i]);
  }
  fprintf(fid,"\n");
  fprintf(fid,"fingerprintconstants:");
  fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
  for (i=1;i<n_body_type;i++) {
    fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
  }
  fprintf(fid,":%s_%d:dr:\n",style,id);
  fprintf(fid,"%f\n",dr);
  fprintf(fid,"fingerprintconstants:");
  fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
  for (i=1;i<n_body_type;i++) {
    fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
  }
  fprintf(fid,":%s_%d:o:\n",style,id);
  fprintf(fid,"%d\n",omin);
  fprintf(fid,"fingerprintconstants:");
  fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
  for (i=1;i<n_body_type;i++) {
    fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
  }
  fprintf(fid,":%s_%d:n:\n",style,id);
  fprintf(fid,"%d\n",nmax);
}

//called after fingerprint is fully defined and tables can be computed.
void Fingerprint_radial::allocate()
{
  int k,m;
  double r1;
  int buf = 5;
  int res = pair->res;
  double cutmax = pair->cutmax;
  radialtable = new double[(res+buf)*get_length()];
  dfctable = new double[res+buf];
  for (k=0;k<(res+buf);k++) {
    r1 = cutmax*cutmax*(double)(k)/(double)(res);
    for (m=0;m<=(nmax-omin);m++) {
      radialtable[k*(nmax-omin+1)+m]=pow(sqrt(r1)/re,m+omin)*exp(-alpha[m]*sqrt(r1)/re)*cutofffunction(sqrt(r1),rc,dr);
    }
    if (sqrt(r1)>=rc || sqrt(r1) <= (rc-dr)) {
        dfctable[k]=0;
    }
    else{
    dfctable[k]=-8*pow(1-(rc-sqrt(r1))/dr,3)/dr/(1-pow(1-(rc-sqrt(r1))/dr,4));
    }
  }
  generate_rinvssqrttable();
}

//called after fingerprint is declared for i-j type, but before its parameters are read.
void Fingerprint_radial::init(int *i,int _id)
{
  empty = false;
  for (int j=0;j<n_body_type;j++) {atomtypes[j] = i[j];}
  id = _id;
}

void Fingerprint_radial::compute_fingerprint(double * features,double * dfeaturesx,double *dfeaturesy,double *dfeaturesz,int ii,int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int * /*jl*/)
{
  int nelements = pair->nelements;
  int res = pair->res;
  int i,jj,itype,jtype,l;
  double delx,dely,delz,rsq;
  int *ilist;
  //
  PairRANN::Simulation *sim = &pair->sims[sid];
  int count=0;
//  double **x = sim->x;
  int *type = sim->type;
  ilist = sim->ilist;
  double cutmax = pair->cutmax;
  i = ilist[ii];
  itype = pair->map[type[i]];
  int f = pair->net[itype].dimensions[0];
  double cutinv2 = 1/cutmax/cutmax;
  //loop over neighbors
  for (jj = 0; jj < jnum; jj++) {
    jtype =tn[jj];
    if (atomtypes[1] != nelements && atomtypes[1] != jtype)continue;
    delx = xn[jj];
    dely = yn[jj];
    delz = zn[jj];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq > rc*rc)continue;
    count = startingneuron;
    double r1 = (rsq*((double)res)*cutinv2);
    int m1 = (int)r1;
    if (m1>res || m1<1) {pair->errorf(FLERR,"invalid neighbor radius!");}
    if (radialtable[m1]==0) {continue;}
    //cubic interpolation from tables
    double *p1 = &radialtable[m1*(nmax-omin+1)];
    double *p2 = &radialtable[(m1+1)*(nmax-omin+1)];
    double *p3 = &radialtable[(m1+2)*(nmax-omin+1)];
    double *p0 = &radialtable[(m1-1)*(nmax-omin+1)];
    double *q = &dfctable[m1-1];
    double *rinvs = &rinvsqrttable[m1-1];
    r1 = r1-trunc(r1);
    double dfc = q[1] + 0.5 * r1*(q[2] - q[0] + r1*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + r1*(3.0*(q[1] - q[2]) + q[3] - q[0])));
    double ri = rinvs[1] + 0.5 * r1*(rinvs[2] - rinvs[0] + r1*(2.0*rinvs[0] - 5.0*rinvs[1] + 4.0*rinvs[2] - rinvs[3] + r1*(3.0*(rinvs[1] - rinvs[2]) + rinvs[3] - rinvs[0])));
    for (l=0;l<=(nmax-omin);l++) {
      double rt = p1[l]+0.5*r1*(p2[l]-p0[l]+r1*(2.0*p0[l]-5.0*p1[l]+4.0*p2[l]-p3[l]+r1*(3.0*(p1[l]-p2[l])+p3[l]-p0[l])));
      features[count]+=rt;
      rt *= (l+omin)/rsq+(-alpha[l]/re+dfc)*ri;
      //update neighbor's features
      dfeaturesx[jj*f+count]+=rt*delx;
      dfeaturesy[jj*f+count]+=rt*dely;
      dfeaturesz[jj*f+count]+=rt*delz;
      //update atom's features
      dfeaturesx[jnum*f+count]-=rt*delx;
      dfeaturesy[jnum*f+count]-=rt*dely;
      dfeaturesz[jnum*f+count]-=rt*delz;
      count++;
    }
  }
}

int Fingerprint_radial::get_length()
{
  return nmax-omin+1;
}
