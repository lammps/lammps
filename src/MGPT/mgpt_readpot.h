// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This file is part of the MGPT implementation. See further comments
   in pair_mgpt.cpp and pair_mgpt.h.
------------------------------------------------------------------------- */

#ifndef READPOT__
#define READPOT__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "mgpt_splinetab.h"

struct potdata {
  double va,vb,vc,vd,ve,p1,al,rp,r00,pn;
  double dva,dvb,dvc,dvd,dve,dp1,dal,drp,dr00;
  double evol0,devol0;
  double (*vpair_spline)[4],(*dvpair_spline)[4];
  double r0,r1;
  int nr;

  double mass,rcrit,rmax;

  int lang,lmax;
  double anorm3,anorm4;
  double ddl[5];

  int ipot,mode;
  char metal[80];

  double input_vol;

  void readpot(const char *parmin_file,const char *potin_file,double vol);

  void eval_pot(double r,double *e_p,double *f_p) {
    double d2y;
    evalspline(nr-1,r0,r1,vpair_spline,r,e_p,f_p,&d2y);
  }

  void eval_vir(double r,double *v_p) {
    double dy,d2y;
    evalspline(nr-1,r0,r1,dvpair_spline,r,v_p,&dy,&d2y);
  }
};


struct potdata2 {
  typedef double (*spline)[4];
  spline va,vb,vc,vd,ve,p1,al,rp,r00;
  spline dva,dvb,dvc,dvd,dve,dp1,dal,drp,dr00;
  spline evol0,devol0;
  double (*vpair)[4][4],(*dvpair)[4][4];
  double r0,r1,T0,T1;
  int nr,nt;

  potdata *potlist;

  double mass,rcrit,rmax;

  int lang,lmax;
  spline ddl[5];

  int ipot,mode;
  char metal[80];

  double input_vol;

  /* Functions to retrieve temperature dependent parameters */

  double eval_tdep(spline y,double T) {
    double f,df,d2f;

    if(0) if(T != 3000.0)
      printf("%s:%d: Error, T = %.3f\n",__FILE__,__LINE__,T);
    evalspline(nt-1,T0,T1,y,T,&f,&df,&d2f);
    return f;
  }
  double eval_tdepderiv(spline y,double T) {
    double f,df,d2f;

    if(0) if(T != 3000.0)
      printf("%s:%d: Error, T = %.3f\n",__FILE__,__LINE__,T);
    evalspline(nt-1,T0,T1,y,T,&f,&df,&d2f);
    return df;
  }


#define make_get(param) \
  double get_##param(double T) { return eval_tdep(param,T); } \
  double get_d##param(double T) { return eval_tdep(d##param,T); } \
  double get_##param##_Tderiv(double T) { return eval_tdepderiv(param,T); }

  /*
#define make_get(param) \
  double get_##param(double T) { if(T != 3000.0) printf("%s:%d: Error, T = %.3f\n",__FILE__,__LINE__,T); return potlist[3].param; } \
  double get_d##param(double T) {  if(T != 3000.0) printf("%s:%d: Error, T = %.3f\n",__FILE__,__LINE__,T); return potlist[3].d##param; } \
  double get_##param##_Tderiv(double T) {  if(T != 3000.0) printf("%s:%d: Error, T = %.3f\n",__FILE__,__LINE__,T); return 0.0; }
  */

  make_get(va)  make_get(vb)  make_get(vc)  make_get(vd)  make_get(ve)
  make_get(p1)  make_get(al)  make_get(rp)  make_get(r00) make_get(evol0)
#undef make_get

  void get_anorm34(double T,double anorm3_p[1],double anorm4_p[1]) {
    double
      s = eval_tdep(ddl[1],T),
      p = eval_tdep(ddl[2],T),
      d = eval_tdep(ddl[3],T),
      f = eval_tdep(ddl[4],T);
    double ss = s*s, pp = p*p, dd = d*d, ff = f*f;
    anorm3_p[0] =  s*ss + 2.0*( p*pp +  d*dd +  f*ff);
    anorm4_p[0] = ss*ss + 2.0*(pp*pp + dd*dd + ff*ff);
  }
  double get_anorm3(double T) {
    double a3,a4;
    get_anorm34(T,&a3,&a4);
    return a3;
  }
  double get_anorm4(double T) {
    double a3,a4;
    get_anorm34(T,&a3,&a4);
    return a4;
  }
  void get_anorm34_Tderiv(double T,double danorm3_p[1],double danorm4_p[1]) {
    double
      s = eval_tdep(ddl[1],T),      ds = eval_tdepderiv(ddl[1],T),
      p = eval_tdep(ddl[2],T),      dp = eval_tdepderiv(ddl[2],T),
      d = eval_tdep(ddl[3],T),      d_d = eval_tdepderiv(ddl[3],T),
      f = eval_tdep(ddl[4],T),      df = eval_tdepderiv(ddl[4],T);
    double ss = s*s, pp = p*p, dd = d*d, ff = f*f;
    danorm3_p[0] =  3.0*ds*ss + 6.0*( dp*pp +  d_d*dd +  df*ff);
    danorm4_p[0] = 4.0*ds*s*ss + 8.0*(dp*p*pp + d_d*d*dd + df*f*ff);
  }
  double get_anorm3_Tderiv(double T) {
    double da3,da4;
    get_anorm34_Tderiv(T,&da3,&da4);
    return da3;
  }
  double get_anorm4_Tderiv(double T) {
    double da3,da4;
    get_anorm34_Tderiv(T,&da3,&da4);
    return da4;
  }

  /* ... */



  char * parsefname(const char *nametemplate,int *i0,int *i1,int *stride) {
    char *s,*p;

    if(0) {
      s = new char[strlen(nametemplate)+1];
    } else {
      int len = 0;
      while(nametemplate[len] != '\0') len = len + 1;
      s = new char[len+1];
    }
    strcpy(s,nametemplate);

    p = strchr(s,'{');
    if(p != nullptr) {
      if(sscanf(p+1,"%d:%d:%d",i0,stride,i1) != 3) {
        fprintf(stderr,"Error in template (\'%s\'), can not parse range.\n",nametemplate);
        exit(1);
      }
      *p = '\0';
    } else {
      *i0 = -1;
      *i1 = -1;
      *stride = 1;
    }
    return s;
  }

  spline maketempspline(int n,potdata data[],double *ptr) {
    int stride = &(data[1].va) - &(data[0].va);
    spline s = new double[n-1][4];

    makespline(n,stride,ptr,s);
    return s;
  }

  void readpot2(const char *parmin_template,const char *potin_template,double vol) {
    int i0,i1,stride,i0x,i1x,stridex;
    char *parmin_file = parsefname(parmin_template,&i0 ,&i1 ,&stride );
    char *potin_file  = parsefname( potin_template,&i0x,&i1x,&stridex);
    int ntemp;

    potdata2 &tdeppot = *this;

    if(i0x != i0 || i1x != i1 || stridex != stride) {
      fprintf(stderr,"Inconsistent templates. parmin_template=\'%s\', potin_template=\'%s\'\n",
              parmin_template,potin_template);
      exit(1);
    }
    if(i0 < 0 || i1 < i0 || stride <= 0 || (i1-i0)/stride+1 < 4) {
      fprintf(stderr,"Improper temperature range. Need at least 4 temperature samples. "
              "i0=%d,i1=%d,stride=%d,basename=\'%s\'\n",
              i0,i1,stride,parmin_file);
      exit(1);
    }

    const char *parmin_suffix = strchr(parmin_template,'}')+1;
    const char * potin_suffix = strchr( potin_template,'}')+1;

    if(parmin_suffix-1 == nullptr) {
      fprintf(stderr,"No closing }. parmin_template=\'%s\'\n",
              parmin_template);
      exit(1);
    }
    if(potin_suffix-1 == nullptr) {
      fprintf(stderr,"No closing }. potin_template=\'%s\'\n",
              potin_template);
      exit(1);
    }

    printf("parmin_template = %s\n"
           "parmin_file = %s\n"
           "parmin_suffix = %s\n"
           "T0=%d , T1=%d , stride=%d\n",
           parmin_template,parmin_file,parmin_suffix,i0,i1,stride);

    ntemp = (i1-i0)/stride + 1;
    /*potdata **/potlist = new potdata[ntemp];
    char *parend = parmin_file+strlen(parmin_file);
    char *potend = potin_file +strlen( potin_file);
    for(int k=0; k<ntemp; k++) {
      sprintf(parend,"%d%s",i0+k*stride,parmin_suffix);
      sprintf(potend,"%d%s",i0+k*stride,potin_suffix);


      printf("Calling readpot(%s,%s,%.3f)\n",
             parmin_file,potin_file,vol);
      potlist[k].readpot(parmin_file,potin_file,vol);

      if(k > 0) {
        if(potlist[k].nr != potlist[k-1].nr) {
          fprintf(stderr,"nr differs between file %d and %d. Exiting.\n",
                  k,k-1);
          exit(1);
        }

        if(potlist[k].r0 != potlist[k-1].r0) {
          fprintf(stderr,"r0 differs between file %d and %d. Exiting.\n",
                  k,k-1);
          exit(1);
        }

        if(potlist[k].r1 != potlist[k-1].r1) {
          fprintf(stderr,"r1 differs between file %d and %d. Exiting.\n",
                  k,k-1);
          exit(1);
        }
      }
    }
    tdeppot.r0 = potlist[0].r0;
    tdeppot.r1 = potlist[0].r1;
    tdeppot.nr = potlist[0].nr;
    tdeppot.T0 = i0;
    tdeppot.T1 = i1;
    tdeppot.nt = ntemp;

    tdeppot.mass = potlist[0].mass;
    tdeppot.rcrit = potlist[0].rcrit;
    tdeppot.rmax = potlist[0].rmax;

    tdeppot.lang = potlist[0].lang;
    tdeppot.lmax = potlist[0].lmax;
    tdeppot.ipot = potlist[0].ipot;
    tdeppot.mode = potlist[0].mode;
    tdeppot.input_vol = potlist[0].input_vol;

    strncpy(tdeppot.metal,potlist[0].metal,sizeof(tdeppot.metal)/sizeof(char));
    tdeppot.metal[sizeof(tdeppot.metal)/sizeof(char) - 1] = '\0';

    delete[] parmin_file;
    delete[] potin_file;

    // Base parameters
    tdeppot.va = maketempspline(ntemp,potlist,&(potlist[0].va));
    tdeppot.vb = maketempspline(ntemp,potlist,&(potlist[0].vb));
    tdeppot.vc = maketempspline(ntemp,potlist,&(potlist[0].vc));
    tdeppot.vd = maketempspline(ntemp,potlist,&(potlist[0].vd));
    tdeppot.ve = maketempspline(ntemp,potlist,&(potlist[0].ve));

    tdeppot.p1 = maketempspline(ntemp,potlist,&(potlist[0].p1));
    tdeppot.al = maketempspline(ntemp,potlist,&(potlist[0].al));
    tdeppot.rp = maketempspline(ntemp,potlist,&(potlist[0].rp));
    tdeppot.r00 = maketempspline(ntemp,potlist,&(potlist[0].r00));

    tdeppot.evol0 = maketempspline(ntemp,potlist,&(potlist[0].evol0));

    // Volume derivatives of base parameters
    tdeppot.dva = maketempspline(ntemp,potlist,&(potlist[0].dva));
    tdeppot.dvb = maketempspline(ntemp,potlist,&(potlist[0].dvb));
    tdeppot.dvc = maketempspline(ntemp,potlist,&(potlist[0].dvc));
    tdeppot.dvd = maketempspline(ntemp,potlist,&(potlist[0].dvd));
    tdeppot.dve = maketempspline(ntemp,potlist,&(potlist[0].dve));

    tdeppot.dp1 = maketempspline(ntemp,potlist,&(potlist[0].dp1));
    tdeppot.dal = maketempspline(ntemp,potlist,&(potlist[0].dal));
    tdeppot.drp = maketempspline(ntemp,potlist,&(potlist[0].drp));
    tdeppot.dr00 = maketempspline(ntemp,potlist,&(potlist[0].dr00));

    tdeppot.devol0 = maketempspline(ntemp,potlist,&(potlist[0].devol0));


    tdeppot.ddl[0] = 0;
    for(int k = 1; k<=4; k++)
      tdeppot.ddl[k] = maketempspline(ntemp,potlist,&(potlist[0].ddl[k]));

    {
      double *v = new double[ntemp];
      double (*C)[4] = new double[ntemp-1][4];

      int sz = (nr-1)*(ntemp-1);
      //printf("Allocation:: nr=%d  ntemp=%d  size=%d\n",nr,ntemp,sz);
      tdeppot.vpair = new double[sz][4][4];
      tdeppot.dvpair = new double[sz][4][4];
      /*
      printf("vpair = %llx , dvpair = %llx",
             (unsigned long long int) tdeppot.vpair,
             (unsigned long long int) tdeppot.dvpair);
      printf("   @@@@@@@@@@@@@@ nr = %d\n",nr);
      */
      for(int i = 0; i<nr-1; i++)
        for(int j = 0; j<4; j++) {
          /*
          if(j == 5)
            printf("    ############### i=%d\n",i);
          */

          /* Make pair interaction interpolation functions */
          for(int k = 0; k<ntemp; k++) {
            if(0) if(i >= potlist[k].nr-1)
              printf("Index error, local_nr=%d, k=%d, i=%d, nr=%d\n",nr,k,i,potlist[k].nr);
            v[k] = potlist[k].vpair_spline[i][j];
          }
          makespline(ntemp,1,v,C);

          for(int k = 0; k<ntemp-1; k++)
            for(int m = 0; m<4; m++)
              tdeppot.vpair[k*(nr-1) + i][j][m] = C[k][m];

          /* Make pair virial interpolation functions */
          for(int k = 0; k<ntemp; k++)
            v[k] = potlist[k].dvpair_spline[i][j];
          makespline(ntemp,1,v,C);

          for(int k = 0; k<ntemp-1; k++)
            for(int m = 0; m<4; m++)
              tdeppot.dvpair[k*(nr-1) + i][j][m] = C[k][m];

        }

      delete[] C;
      delete[] v;
    }

  }

  void eval2Dspline(double fun[][4][4],double r,double T,double *e_p,double *f_p,double *dedT_p) {
    double C[4],dC[4],dd;
    double Tfrac = (T-T0)/(T1-T0) * (nt-1);
    double rfrac = (r-r0)/(r1-r0) * (nr-1);
    int k = (int) Tfrac;
    int i = (int) rfrac;
    int j;

    if(k < 0) k = 0;
    else if(k > nt-2) k = nt-2;

    if(i < 0) i = 0;
    else if(i > nr-2) i = nr-2;

    /*
    printf("eval_pot nr=%d  nt=%d\n",nr,nt);
    printf("eval_pot r=%.3f  T=%.3f  k=%d  i=%d\n",
           r,T,k,i);
    printf("Tfrac=%.3f  Tfrac-k=%.3f\n",Tfrac,Tfrac-k);
    printf("rfrac=%.3f  rfrac-i=%.3f\n",rfrac,rfrac-i);
    */
    for(j = 0; j<4; j++) {
      evalcubic(fun[k*(nr-1) + i][j],Tfrac-k,&C[j],&dC[j],&dd);
      dC[j] = dC[j] * ((nt-1) / (T1-T0));
    }
    /*
    printf("C coeff: %.3e  %.3e  %.3e  %.3e\n",
           C[0],C[1],C[2],C[3]);
    */
    evalcubic(C,rfrac-i,e_p,f_p,&dd);
    evalcubic(dC,rfrac-i,dedT_p,&dd,&dd);
    *f_p *= (nr-1) / (r1-r0);
  }
  void eval_pot(double r,double T,double *e_p,double *f_p,double *dedT_p) {
    eval2Dspline(vpair,r,T,e_p,f_p,dedT_p);
  }
  void eval_vir(double r,double T,double *v_p) {
    double vf,dvdT;
    eval2Dspline(dvpair,r,T,v_p,&vf,&dvdT);
  }
};


#endif
