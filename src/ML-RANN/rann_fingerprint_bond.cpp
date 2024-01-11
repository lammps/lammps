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

#include "rann_fingerprint_bond.h"
#include "pair_rann.h"

#include <cmath>

using namespace LAMMPS_NS::RANN;

Fingerprint_bond::Fingerprint_bond(PairRANN *_pair) : Fingerprint(_pair)
{
  n_body_type = 3;
  dr = 0;
  re = 0;
  rc = 0;
  alpha_k = new double[1];
  alpha_k[0] = -1;
  kmax = 0;
  mlength = 0;
  id = -1;
  style = "bond";
  atomtypes = new int[n_body_type];
  empty = true;
  _pair->allscreen = false;
}

Fingerprint_bond::~Fingerprint_bond() {
  delete[] alpha_k;
  delete[] atomtypes;
  delete[] expcuttable;
  delete[] dfctable;
  for (int i=0;i<(mlength*(mlength+1))>>1;i++) {
    delete[] coeff[i];
    delete[] coeffx[i];
    delete[] coeffy[i];
    delete[] coeffz[i];
    delete[] Mf[i];
  }
  delete[] coeff;
  delete[] coeffx;
  delete[] coeffy;
  delete[] coeffz;
  delete[] Mf;
  delete[] rinvsqrttable;
}

bool Fingerprint_bond::parse_values(std::string constant,std::vector<std::string> line1) {
  int nwords,l;
  nwords=line1.size();
  if (constant.compare("re")==0) {
    re = strtod(line1[0].c_str(),nullptr);
  }
  else if (constant.compare("rc")==0) {
    rc = strtod(line1[0].c_str(),nullptr);
  }
  else if (constant.compare("alphak")==0) {
    delete[] alpha_k;
    alpha_k = new double[nwords];
    for (l=0;l<nwords;l++) {
      alpha_k[l]=strtod(line1[l].c_str(),nullptr);
    }
  }
  else if (constant.compare("dr")==0) {
    dr = strtod(line1[0].c_str(),nullptr);
  }
  else if (constant.compare("k")==0) {
    kmax = strtol(line1[0].c_str(),nullptr,10);
  }
  else if (constant.compare("m")==0) {
    mlength = strtol(line1[0].c_str(),nullptr,10);
  }
  else pair->errorf(FLERR,"Undefined value for bond power");
  if (re!=0.0 && rc!=0.0 && alpha_k[0]!=-1 && dr!=0.0 && mlength!=0 && kmax!=0)return true;
  return false;
}

void Fingerprint_bond::write_values(FILE *fid) {
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
  fprintf(fid,":%s_%d:alphak:\n",style,id);
  for (i=0;i<kmax;i++) {
    fprintf(fid,"%f ",alpha_k[i]);
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
  fprintf(fid,":%s_%d:k:\n",style,id);
  fprintf(fid,"%d\n",kmax);
  fprintf(fid,"fingerprintconstants:");
  fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
  for (i=1;i<n_body_type;i++) {
    fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
  }
  fprintf(fid,":%s_%d:m:\n",style,id);
  fprintf(fid,"%d\n",mlength);
}

void Fingerprint_bond::init(int *i,int _id) {
  for (int j=0;j<n_body_type;j++) {atomtypes[j] = i[j];}
  re = 0;
  rc = 0;
  mlength = 0;
  kmax = 0;
  delete[] alpha_k;
  alpha_k = new double[1];
  alpha_k[0]=-1;
  empty = false;
  id = _id;
}

//number of neurons defined by this fingerprint
int Fingerprint_bond::get_length() {
  return mlength*kmax;
}

void Fingerprint_bond::allocate() {
  generate_exp_cut_table();
  generate_coefficients();
  generate_rinvssqrttable();
}

//Generate table of complex functions for quick reference during compute. Used by do3bodyfeatureset_singleneighborloop and do3bodyfeatureset_doubleneighborloop.
void Fingerprint_bond::generate_exp_cut_table() {
  int m,n;
  double r1;
  int buf = 5;
  int res = pair->res;
  double cutmax = pair->cutmax;
  expcuttable = new double[(res+buf)*(kmax)];
  dfctable = new double[res+buf];
  for (m=0;m<(res+buf);m++) {
    r1 = cutmax*cutmax*(double)(m)/(double)(res);
    for (n=0;n<(kmax);n++) {
      expcuttable[n+m*(kmax)] = exp(-alpha_k[n]/re*sqrt(r1))*cutofffunction(sqrt(r1),rc,dr);
    }
    if (sqrt(r1)>=rc || sqrt(r1) <= (rc-dr)) {
      dfctable[m]=0;
    }
    else{
      dfctable[m]=-8*pow(1-(rc-sqrt(r1))/dr,3)/dr/(1-pow(1-(rc-sqrt(r1))/dr,4));
    }
  }
}

//Generate table of complex functions for quick reference during compute. Used by do3bodyfeatureset_singleneighborloop.
void Fingerprint_bond::generate_coefficients() {      //calculates multinomial coefficient for each term
  int p,mb,mc;
  int n,p1,i1;
  mb = mlength;
  mc=(mb*(mb+1))>>1;
  coeff = new int *[mc];
  coeffx = new int *[mc];
  coeffy = new int *[mc];
  coeffz = new int *[mc];
  for (p=0;p<mc;p++) {
    coeff[p]=new int [mb];
    coeffx[p]=new int [mb];
    coeffy[p]=new int [mb];
    coeffz[p]=new int [mb];
  }
  Mf = new int*[mc];
  int *M = new int[mlength+1];
  for (p=0;p<mlength+1;p++) {
    M[p]=0;
  }
  for (p1=0;p1<mc;p1++) {
    Mf[p1] = new int[mlength+1];
    for (p=0;p<mlength+1;p++) {
      Mf[p1][p]=0;
    }
  }
  M[0] = 2;
  Mf[0][0] = 2;
  n = 1;
  int m1 = 1;
  bool go = true;
  bool broke = false;
  while (go) {
    broke = false;
    for (i1=0;i1<mlength-1;i1++) {
      if (M[i1+1] == 0) {
        M[i1+1]=M[i1+1]+1;
        for (p1=0;p1<mlength+1;p1++) {
          Mf[n][p1] = M[p1];
        }
        n = n+1;
        broke = true;
        break;
      }
    }
    if (m1<mlength && !broke) {
      M[m1]=M[m1]+1;
      for (p1=m1+1;p1<mlength+1;p1++) {
        M[p1]=0;
      }
      for (p1=0;p1<mlength+1;p1++) {
        Mf[n][p1]=M[p1];
      }
      n=n+1;
      broke = true;
      m1 = m1+1;
    }
    if (!broke) {
      go = false;
    }
  }
  for (p=0;p<mb;p++) {
    for (p1=0;p1<mc;p1++) {
      if (p==0) {
        coeffx[p1][p]=0;
        coeffy[p1][p]=0;
        coeffz[p1][p]=0;
      }
      else{
        coeffx[p1][p]=coeffx[p1][p-1];
        if (Mf[p1][p]==0) {
          coeffx[p1][p]++;
        }
        coeffy[p1][p]=coeffy[p1][p-1];
        if (Mf[p1][p]==1) {
          coeffy[p1][p]++;
        }
        coeffz[p1][p]=coeffz[p1][p-1];
        if (Mf[p1][p]==2) {
          coeffz[p1][p]++;
        }
      }
      coeff[p1][p] = pair->factorial(p)/pair->factorial(coeffx[p1][p])/pair->factorial(coeffy[p1][p])/pair->factorial(coeffz[p1][p]);
    }
  }
  delete[] M;
}


//Called by getproperties. Gets 3-body features and dfeatures
void Fingerprint_bond::compute_fingerprint(double * features,double * dfeaturesx,double *dfeaturesy,double *dfeaturesz, int ii,int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int *jl) {
  //select the more efficient algorithm for this particular potential and environment.
  if (jnum*2>(mlength+1)*mlength*20) {
    do3bodyfeatureset_singleneighborloop(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,sid,xn,yn,zn,tn,jnum,jl);
  }
  else{
    do3bodyfeatureset_doubleneighborloop(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,sid,xn,yn,zn,tn,jnum,jl);

  }
}

//Called by do3bodyfeatureset. Algorithm for high neighbor numbers and small series of bond angle powers
void Fingerprint_bond::do3bodyfeatureset_singleneighborloop(double * features,double * dfeaturesx,double *dfeaturesy,double *dfeaturesz,int ii,int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int * /*jl*/) {
  int i,jj,itype,jtype,kk,m,n,mcount,a,a1,a2,ai;
  double delx,dely,delz,rsq;
  int *ilist;
  int count=0;
  PairRANN::Simulation *sim = &pair->sims[sid];
  int *type = sim->type;
  double cutmax = pair->cutmax;
  int res = pair->res;
  double cutinv2 = 1/cutmax/cutmax;
  ilist = sim->ilist;
  int nelements=pair->nelements;
  i = ilist[ii];
  itype = pair->map[type[i]];
  int f = pair->net[itype].dimensions[0];
  double expr[jnum][kmax+12];
  int p = kmax;
  int countmb=((mlength)*(mlength+1))>>1;
  // calculate interpolation expr, rinvs and dfc, for each neighbor
  for (jj = 0; jj < jnum; jj++) {
    jtype = tn[jj];
    if (atomtypes[1] != nelements && atomtypes[1] != jtype && atomtypes[2] != nelements && atomtypes[2] != jtype) {
      expr[jj][0]=0;
      continue;
    }
    delx = xn[jj];
    dely = yn[jj];
    delz = zn[jj];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq>rc*rc) {
    expr[jj][0]=0;
    continue;
    }
    double r1 = (rsq*((double)res)*cutinv2);
    int m1 = (int)r1;
    r1 = r1-trunc(r1);
    double *p0 = &expcuttable[(m1-1)*kmax];
    double *p1 = &expcuttable[m1*kmax];
    double *p2 = &expcuttable[(m1+1)*kmax];
    double *p3 = &expcuttable[(m1+2)*kmax];
    for (kk=0;kk<kmax;kk++) {
      expr[jj][kk] = p1[kk]+0.5*r1*(p2[kk]-p0[kk]+r1*(2.0*p0[kk]-5.0*p1[kk]+4.0*p2[kk]-p3[kk]+r1*(3.0*(p1[kk]-p2[kk])+p3[kk]-p0[kk])));
    }
    double* q = &dfctable[m1-1];
    double* ri = &rinvsqrttable[m1-1];
    double dfc = q[1] + 0.5 * r1*(q[2] - q[0] + r1*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + r1*(3.0*(q[1] - q[2]) + q[3] - q[0])));
    double rinvs = ri[1] + 0.5 * r1*(ri[2] - ri[0] + r1*(2.0*ri[0] - 5.0*ri[1] + 4.0*ri[2] - ri[3] + r1*(3.0*(ri[1] - ri[2]) + ri[3] - ri[0])));

    expr[jj][p]=delx*rinvs;
    expr[jj][p+1]=dely*rinvs;
    expr[jj][p+2]=delz*rinvs;
    //Hack to avoid nan when x y or z component of radial vector is exactly 0. Shouldn't affect accuracy.
    if (expr[jj][p]*expr[jj][p]<0.000000000001) {
      expr[jj][p] = 0.000001;
    }
    if (expr[jj][p+1]*expr[jj][p+1]<0.000000000001) {
      expr[jj][p+1] = 0.000001;
    }
    if (expr[jj][p+2]*expr[jj][p+2]<0.000000000001) {
      expr[jj][p+2] = 0.000001;
    }
    expr[jj][p+3] = -dfc*expr[jj][p];
    expr[jj][p+4] = rinvs/expr[jj][p];
    expr[jj][p+5] = rinvs*expr[jj][p];
    expr[jj][p+6] = -dfc*expr[jj][p+1];
    expr[jj][p+7] = rinvs/expr[jj][p+1];
    expr[jj][p+8] = rinvs*expr[jj][p+1];
    expr[jj][p+9] = -dfc*expr[jj][p+2];
    expr[jj][p+10] = rinvs/expr[jj][p+2];
    expr[jj][p+11] = rinvs*expr[jj][p+2];
  }

  int kb = kmax;
  int mb = mlength;
  count = startingneuron;
  double Bb[mb];
  double dBbx;
  double dBby;
  double dBbz;
  double yprod;
  for (mcount=0;mcount<countmb;mcount++) {
    int *M = Mf[mcount];
    int *_coeffx = coeffx[mcount];
    int *_coeffy = coeffy[mcount];
    int *_coeffz = coeffz[mcount];
    int *_coeff = coeff[mcount];
    a = mb+1;
    for (a1=0;a1<mb;a1++) {
      if (Mf[mcount][a1+1]==0) {
        a = a1;
        break;
      }
    }
    for (n=0;n<kb;n++) {
      for (a1=0;a1<mb;a1++) {
        Bb[a1]=0;
      }
      ai = n;
      double y1 = alpha_k[ai]/re;
      //loop over jtype to get Bb
      for (jj=0;jj<jnum;jj++) {
        if (expr[jj][0]==0) {continue;}
        jtype = tn[jj];
        if (atomtypes[1] != nelements && atomtypes[1] != jtype) {
          continue;
        }
        double yprod = expr[jj][ai];
        double *y4 = &expr[jj][p];
        for (a2=0;a2<a;a2++) {
          yprod *= y4[M[a2+1]];
        }
        for (a2=a;a2<mb;a2++) {
          Bb[a2]=Bb[a2]+yprod;
          yprod *= y4[M[a2+1]];
        }
      }
      if (atomtypes[1]!=atomtypes[2]) {//Bb!=Bg
        double Bg[mb];
        for (a1=0;a1<mb;a1++) {
          Bg[a1]=0;
        }
        ai = n;
        double y1 = alpha_k[ai]/re;
        //loop over ktype to get Bg
        for (jj=0;jj<jnum;jj++) {
          if (expr[jj][0]==0) {continue;}
          jtype = tn[jj];
          if (atomtypes[2] != nelements && atomtypes[2] != jtype) {
            continue;
          }
          double yprod = expr[jj][ai];
          double *y4 = &expr[jj][p];
          for (a2=0;a2<a;a2++) {
            yprod *= y4[M[a2+1]];
          }
          for (a2=a;a2<mb;a2++) {
            Bg[a2]=Bg[a2]+yprod;
            yprod *= y4[M[a2+1]];
          }
        }
        double B1;
        //loop over ktype to get dBg*Bb
        for (jj=0;jj<jnum;jj++) {
          if (expr[jj][0]==0) {continue;}
          jtype = tn[jj];
          if (atomtypes[2] != nelements && atomtypes[2] != jtype) {
            continue;
          }
          double *y3 = &expr[jj][p+3];
          double *y4 = &expr[jj][p];
          ai = n;
          yprod = expr[jj][ai];
          for (a2=0;a2<a;a2++) {
            yprod *= y4[M[a2+1]];
          }
          ai = n*(mb)+a+count+jj*f;
          for (a2=a;a2<mb;a2++) {
            B1 = Bb[a2]*_coeff[a2]*yprod;
            dBbx = -B1*(y1*y4[0]+y3[0]-_coeffx[a2]*y3[1]+a2*y3[2]);
            dBby = -B1*(y1*y4[1]+y3[3]-_coeffy[a2]*y3[4]+a2*y3[5]);
            dBbz = -B1*(y1*y4[2]+y3[6]-_coeffz[a2]*y3[7]+a2*y3[8]);
            dfeaturesx[ai] += dBbx;
            dfeaturesy[ai] += dBby;
            dfeaturesz[ai] += dBbz;
            yprod *= y4[M[a2+1]];
            ai++;
          }
        }
        //loop over jtype to get dBb*Bg
        for (jj=0;jj<jnum;jj++) {
          if (expr[jj][0]==0) {continue;}
          jtype = tn[jj];
          if (atomtypes[1] != nelements && atomtypes[1] != jtype) {
            continue;
          }
          double *y3 = &expr[jj][p+3];
          double *y4 = &expr[jj][p];
          ai = n;
          yprod = expr[jj][ai];
          for (a2=0;a2<a;a2++) {
            yprod *= y4[M[a2+1]];
          }
          ai = n*(mb)+a+count+jj*f;
          for (a2=a;a2<mb;a2++) {
            B1 = Bg[a2]*_coeff[a2]*yprod;
            dBbx = -B1*(y1*y4[0]+y3[0]-_coeffx[a2]*y3[1]+a2*y3[2]);
            dBby = -B1*(y1*y4[1]+y3[3]-_coeffy[a2]*y3[4]+a2*y3[5]);
            dBbz = -B1*(y1*y4[2]+y3[6]-_coeffz[a2]*y3[7]+a2*y3[8]);
            dfeaturesx[ai] += dBbx;
            dfeaturesy[ai] += dBby;
            dfeaturesz[ai] += dBbz;
            yprod *= y4[M[a2+1]];
            ai++;
          }
        }
        //central atom derivative
        for (a2=a;a2<mb;a2++) {
          ai = n*(mb)+a2+count+jnum*f;
          features[ai-jnum*f] += Bb[a2]*Bg[a2]*_coeff[a2];
        }
      }
      else{//Bb=Bg
        double B1;
        //loop over jtype to get 2*Bb*dBb
        for (jj=0;jj<jnum;jj++) {
          if (expr[jj][0]==0) {continue;}
          jtype = tn[jj];
          if (atomtypes[1] != nelements && atomtypes[1] != jtype) {
            continue;
          }
          double *y3 = &expr[jj][p+3];
          double *y4 = &expr[jj][p];
          ai = n;
          yprod = expr[jj][ai];
          for (a2=0;a2<a;a2++) {
            yprod *= y4[M[a2+1]];
          }
          ai = n*(mb)+a+count+jj*f;
          for (a2=a;a2<mb;a2++) {
            B1 = 2*Bb[a2]*_coeff[a2]*yprod;
            dBbx = -B1*(y1*y4[0]+y3[0]-_coeffx[a2]*y3[1]+a2*y3[2]);
            dBby = -B1*(y1*y4[1]+y3[3]-_coeffy[a2]*y3[4]+a2*y3[5]);
            dBbz = -B1*(y1*y4[2]+y3[6]-_coeffz[a2]*y3[7]+a2*y3[8]);
            dfeaturesx[ai] += dBbx;
            dfeaturesy[ai] += dBby;
            dfeaturesz[ai] += dBbz;
            yprod *= y4[M[a2+1]];
            ai++;
          }
        }
        //central atom derivative
        for (a2=a;a2<mb;a2++) {
          ai = n*(mb)+a2+count+jnum*f;
          features[ai-jnum*f] += Bb[a2]*Bb[a2]*_coeff[a2];
        }
      }
    }
  }
  for (jj=0;jj<jnum;jj++) {
    if (expr[jj][0]==0) {continue;}
    count = startingneuron;
    for (n=0;n<kb;n++) {
      for (m=0;m<mb;m++) {
        dfeaturesx[jnum*f+count]-=dfeaturesx[jj*f+count];
        dfeaturesy[jnum*f+count]-=dfeaturesy[jj*f+count];
        dfeaturesz[jnum*f+count]-=dfeaturesz[jj*f+count];
        count++;
      }
    }
  }
}

//Called by do3bodyfeatureset. Algorithm for low neighbor numbers and large series of bond angle powers
  int *ilist,*jlist,*numneigh,**firstneigh;
void Fingerprint_bond::do3bodyfeatureset_doubleneighborloop(double * features,double * dfeaturesx,double *dfeaturesy,double *dfeaturesz,int ii, int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int * /*jl*/) {
  int i,jj,itype,jtype,ktype,kk,m,n;
  double delx,dely,delz,rsq;
  int jtypes = atomtypes[1];
  int ktypes = atomtypes[2];
  int count=0;
  PairRANN::Simulation *sim = &pair->sims[sid];
  int *type = sim->type;
  int nelements = pair->nelements;
  int res = pair->res;
  double cutmax = pair->cutmax;
  double cutinv2 = 1/cutmax/cutmax;
  ilist = sim->ilist;
  i = ilist[ii];
  itype = pair->map[type[i]];
  int f = pair->net[itype].dimensions[0];
  double expr[jnum][kmax];
  double y[jnum][3];
  double ri[jnum];
  double dfc[jnum];
  int kb = kmax;
  int mb = mlength;
  double c41[kmax];
  double c51[kmax];
  double c61[kmax];
  double ct[kmax];
  for (jj = 0; jj < jnum; jj++) {
    jtype = tn[jj];
    if (jtypes != nelements && jtypes != jtype && ktypes != nelements && ktypes != jtype) {
      expr[jj][0]=0;
      continue;
    }
    delx = xn[jj];
    dely = yn[jj];
    delz = zn[jj];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq>rc*rc) {
        expr[jj][0]=0;
        continue;
    }
    double r1 = (rsq*((double)res)*cutinv2);
    int m1 = (int)r1;
    if (!(m1>=1 && m1 <= res))pair->errorf(FLERR,"Neighbor list is invalid.");//usually results from nan somewhere.
    r1 = r1-trunc(r1);
    double *p0 = &expcuttable[(m1-1)*kmax];
    double *p1 = &expcuttable[m1*kmax];
    double *p2 = &expcuttable[(m1+1)*kmax];
    double *p3 = &expcuttable[(m1+2)*kmax];
    for (kk=0;kk<kmax;kk++) {
        expr[jj][kk] = p1[kk]+0.5*r1*(p2[kk]-p0[kk]+r1*(2.0*p0[kk]-5.0*p1[kk]+4.0*p2[kk]-p3[kk]+r1*(3.0*(p1[kk]-p2[kk])+p3[kk]-p0[kk])));
    }
    double* q = &dfctable[m1-1];
    double* r2 = &rinvsqrttable[m1-1];
    dfc[jj] = q[1] + 0.5 * r1*(q[2] - q[0] + r1*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + r1*(3.0*(q[1] - q[2]) + q[3] - q[0])));
    ri[jj] = r2[1] + 0.5 * r1*(r2[2] - r2[0] + r1*(2.0*r2[0] - 5.0*r2[1] + 4.0*r2[2] - r2[3] + r1*(3.0*(r2[1] - r2[2]) + r2[3] - r2[0])));
    y[jj][0]=delx*ri[jj];
    y[jj][1]=dely*ri[jj];
    y[jj][2]=delz*ri[jj];
  }
  for (jj = 0; jj < jnum; jj++) {
    if (expr[jj][0]==0)continue;
    jtype = tn[jj];
    if (jtypes != nelements && jtypes != jtype) {
      continue;
    }
    for (n = 0;n<kmax;n++) {
      ct[n] = 2*alpha_k[n]/re;
      c41[n]=(-ct[n]+2*dfc[jj])*y[jj][0];
      c51[n]=(-ct[n]+2*dfc[jj])*y[jj][1];
      c61[n]= (-ct[n]+2*dfc[jj])*y[jj][2];
    }
    if (jtypes==ktypes) {
      for (kk=jj+1;kk< jnum; kk++) {
        if (expr[kk][0]==0)continue;
        ktype = tn[kk];
        if (ktypes != nelements && ktypes != ktype) {
          continue;
        }
        count = startingneuron;
        double dot = (y[jj][0]*y[kk][0]+y[jj][1]*y[kk][1]+y[jj][2]*y[kk][2]);
        double c1  = 2*ri[jj]*(y[kk][0]-dot*y[jj][0]);
        double c2  = 2*ri[jj]*(y[kk][1]-dot*y[jj][1]);
        double c3  = 2*ri[jj]*(y[kk][2]-dot*y[jj][2]);
        double c10 = 2*ri[kk]*(y[jj][0]-dot*y[kk][0]);
        double c11 = 2*ri[kk]*(y[jj][1]-dot*y[kk][1]);
        double c12 = 2*ri[kk]*(y[jj][2]-dot*y[kk][2]);
        //alternate formulation:
//        double c1 = 2*ri[jj]*y[kk][0]*(1-y[jj][0]*y[jj][0]);
//        double c2 = 2*ri[jj]*y[kk][1]*(1-y[jj][1]*y[jj][1]);
//        double c3 = 2*ri[jj]*y[kk][2]*(1-y[jj][2]*y[jj][2]);
//        double c10 = 2*ri[kk]*y[jj][0]*(1-y[kk][0]*y[kk][0]);
//        double c11 = 2*ri[kk]*y[jj][1]*(1-y[kk][1]*y[kk][1]);
//        double c12 = 2*ri[kk]*y[jj][2]*(1-y[kk][2]*y[kk][2]);
        for (n=0;n<kb;n++) {
          double dot1=expr[jj][n]*expr[kk][n];
          double c4 = c41[n];
          double c5 = c51[n];
          double c6 = c61[n];
          double ct2 = -ct[n]+2*dfc[kk];
          double c42 = ct2*y[kk][0];
          double c52 = ct2*y[kk][1];
          double c62 = ct2*y[kk][2];
          //m=0
          features[count]+=2*dot1;
          dfeaturesx[jj*f+count]+=dot1*c4;
          dfeaturesy[jj*f+count]+=dot1*c5;
          dfeaturesz[jj*f+count]+=dot1*c6;
          dfeaturesx[kk*f+count]+=dot1*c42;
          dfeaturesy[kk*f+count]+=dot1*c52;
          dfeaturesz[kk*f+count]+=dot1*c62;
          c4*=dot;
          c5*=dot;
          c6*=dot;
          c42*=dot;
          c52*=dot;
          c62*=dot;
          count++;
          for (m=1;m<mb;m++) {
            double c7 = dot1*(m*c1+c4);
            double c8 = dot1*(m*c2+c5);
            double c9 = dot1*(m*c3+c6);
            dfeaturesx[jj*f+count]+=c7;
            dfeaturesy[jj*f+count]+=c8;
            dfeaturesz[jj*f+count]+=c9;
            dfeaturesx[kk*f+count]+=dot1*(m*c10+c42);
            dfeaturesy[kk*f+count]+=dot1*(m*c11+c52);
            dfeaturesz[kk*f+count]+=dot1*(m*c12+c62);
            dot1*=dot;
            features[count++]+=2*dot1;
          }
        }
      }
      kk=jj;
      if (ktypes == nelements || ktypes == jtype) {
        count = startingneuron;
        double dot = (y[jj][0]*y[kk][0]+y[jj][1]*y[kk][1]+y[jj][2]*y[kk][2]);
        double c1 = 2*ri[jj]*(y[kk][0]-dot*y[jj][0]);
        double c2 = 2*ri[jj]*(y[kk][1]-dot*y[jj][1]);
        double c3 = 2*ri[jj]*(y[kk][2]-dot*y[jj][2]);
        //alternate formulation:
//        double c1 = 2*ri[jj]*y[kk][0]*(1-y[jj][0]*y[jj][0]);
//        double c2 = 2*ri[jj]*y[kk][1]*(1-y[jj][1]*y[jj][1]);
//        double c3 = 2*ri[jj]*y[kk][2]*(1-y[jj][2]*y[jj][2]);
        for (n=0;n<kb;n++) {
          double dot1=expr[jj][n]*expr[kk][n];
          double c4 = c41[n];
          double c5 = c51[n];
          double c6 = c61[n];
          //m=0
          features[count]+=dot1;
          dfeaturesx[jj*f+count]+=dot1*c4;
          dfeaturesy[jj*f+count]+=dot1*c5;
          dfeaturesz[jj*f+count]+=dot1*c6;
          c4*=dot;
          c5*=dot;
          c6*=dot;
          count++;
          for (m=1;m<mb;m++) {
            double c7 = dot1*(m*c1+c4);
            double c8 = dot1*(m*c2+c5);
            double c9 = dot1*(m*c3+c6);
            dfeaturesx[jj*f+count]+=c7;
            dfeaturesy[jj*f+count]+=c8;
            dfeaturesz[jj*f+count]+=c9;
            dot1*=dot;
            features[count++]+=dot1;
          }
        }
      }
    }
    else {
      for (kk=0;kk<jnum; kk++) {
        if (expr[kk][0]==0)continue;
        ktype = tn[kk];
        if (ktypes != nelements && ktypes != ktype) {
          continue;
        }
        count = startingneuron;
        double dot = (y[jj][0]*y[kk][0]+y[jj][1]*y[kk][1]+y[jj][2]*y[kk][2]);
        double c1  = ri[jj]*(y[kk][0]-dot*y[jj][0]);
        double c2  = ri[jj]*(y[kk][1]-dot*y[jj][1]);
        double c3  = ri[jj]*(y[kk][2]-dot*y[jj][2]);
        double c10 = ri[kk]*(y[jj][0]-dot*y[kk][0]);
        double c11 = ri[kk]*(y[jj][1]-dot*y[kk][1]);
        double c12 = ri[kk]*(y[jj][2]-dot*y[kk][2]);
        //alternate formulation:
//        double c1 = 2*ri[jj]*y[kk][0]*(1-y[jj][0]*y[jj][0]);
//        double c2 = 2*ri[jj]*y[kk][1]*(1-y[jj][1]*y[jj][1]);
//        double c3 = 2*ri[jj]*y[kk][2]*(1-y[jj][2]*y[jj][2]);
//        double c10 = 2*ri[kk]*y[jj][0]*(1-y[kk][0]*y[kk][0]);
//        double c11 = 2*ri[kk]*y[jj][1]*(1-y[kk][1]*y[kk][1]);
//        double c12 = 2*ri[kk]*y[jj][2]*(1-y[kk][2]*y[kk][2]);
        for (n=0;n<kb;n++) {
          double dot1=expr[jj][n]*expr[kk][n];
          double c4 = c41[n]/2;
          double c5 = c51[n]/2;
          double c6 = c61[n]/2;
          double ct2 = -ct[n]/2+dfc[kk];
          double c42 = ct2*y[kk][0];
          double c52 = ct2*y[kk][1];
          double c62 = ct2*y[kk][2];
          //m=0
          features[count]+=dot1;
          dfeaturesx[jj*f+count]+=dot1*c4;
          dfeaturesy[jj*f+count]+=dot1*c5;
          dfeaturesz[jj*f+count]+=dot1*c6;
          dfeaturesx[kk*f+count]+=dot1*c42;
          dfeaturesy[kk*f+count]+=dot1*c52;
          dfeaturesz[kk*f+count]+=dot1*c62;
          c4*=dot;
          c5*=dot;
          c6*=dot;
          c42*=dot;
          c52*=dot;
          c62*=dot;
          count++;
          for (m=1;m<mb;m++) {
            double c7 = dot1*(m*c1+c4);
            double c8 = dot1*(m*c2+c5);
            double c9 = dot1*(m*c3+c6);
            dfeaturesx[jj*f+count]+=c7;
            dfeaturesy[jj*f+count]+=c8;
            dfeaturesz[jj*f+count]+=c9;
            dfeaturesx[kk*f+count]+=dot1*(m*c10+c42);
            dfeaturesy[kk*f+count]+=dot1*(m*c11+c52);
            dfeaturesz[kk*f+count]+=dot1*(m*c12+c62);
            dot1*=dot;
            features[count++]+=dot1;
          }
        }
      }
    }
  }
  for (jj=0;jj<jnum;jj++) {
    if (expr[jj][0]==0) {continue;}
    count = startingneuron;
    for (n=0;n<kb;n++) {
      for (m=0;m<mb;m++) {
        dfeaturesx[jnum*f+count]-=dfeaturesx[jj*f+count];
        dfeaturesy[jnum*f+count]-=dfeaturesy[jj*f+count];
        dfeaturesz[jnum*f+count]-=dfeaturesz[jj*f+count];
        count++;
      }
    }
  }
}


