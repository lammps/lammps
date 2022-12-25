// clang-format off
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
   This file is part of the MGPT implementation. See further comments
   in pair_mgpt.cpp and pair_mgpt.h.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "mgpt_splinetab.h"

#include "mgpt_readpot.h"

static double fgauss(double x,double al) {
  return exp(-al * pow(x-1.0, 2));
}
static double hgauss(double x,double al) {
  return (1.0 + al * pow(x-1.0, 2)) * exp(-al * pow(x-1.0, 2));
}
static double fl(double r,int mode,double rp,double p1,double al,double r0,double pn)
{
  double term;
  //double pn=1.0;
  if (mode <= 4)
    term = pow(rp/r, p1);
  else
    term = exp(-p1*(pow(r/rp, pn) - 1.0)/pn);

  if (r <= r0) return term;
  double quan = al*(r/r0 - 1.0)*(r/r0 - 1.0);
  if (mode <= 2)
    return term*exp(-quan);
  else
    return term*(1.0 + quan)*exp(-quan);
}

static int cmp_double(const void *ap,const void *bp) {
  double a = *((const double *) ap);
  double b = *((const double *) bp);
  if (a < b)
    return -1;
  else if (a > b)
    return 1;
  else
    return 0;
}
static void getparmindata(const char *potin_file,int nvol[1],double vol0[1],double x0[1],double x1[1]) {
  int n,vsize;
  double *volarr;
  char metal[80],metalx[80];
  int ipot,ipotx,mode,modex;
  FILE *in = fopen(potin_file,"r");
  char line[1024];

  if (in == nullptr) {
    fprintf(stderr,"@%s:%d: Error reading potin file. Can not open file \'%s\'.\n",
            __FILE__,__LINE__,potin_file);
    exit(1);
  }

  vsize = 10;
  volarr = (double *) malloc(sizeof(double) * vsize);
  n = 0;
  while (fgets(line,sizeof(line),in) != nullptr) {
    double zval,ivol,rws,mass;
    double r0x,r1x,drx;
    int nrx,i;

    if (line[strspn(line," \t")] == '#') continue;

    if (n == 0) {
      metal[0] = 0;
      if (sscanf(line,"%s %d %d",metal,&ipot,&mode) != 3) {
        fprintf(stderr,"@%s:%d: Error on potin file. line = %s\n",
                __FILE__,__LINE__,line);
        exit(1);
      }
    } else {
      metalx[0] = 0;
      if (sscanf(line,"%s %d %d",metalx,&ipotx,&modex) != 3) {
        fprintf(stderr,"@%s:%d: Error on potin file. line = %s\n",
                __FILE__,__LINE__,line);
        exit(1);
      } else if (strcmp(metal,metalx) != 0 || ipot != ipotx || mode != modex) {
        fprintf(stderr,"@%s:%d: Error on potin file, parameter mismatch:\n"
                "    metal  = \'%s\'    ipot  = %d    mode  = %d\n"
                "    metalx = \'%s\'    ipotx = %d    modex = %d\n",
                __FILE__,__LINE__,
                metal,ipot,mode,
                metalx,ipotx,modex);
        exit(1);
      }
    }

    fgets(line,sizeof(line),in);
    sscanf(line,"%lf %lf %lf %lf",&zval,&ivol,&rws,&mass);
    if (n >= vsize) {
      vsize = 2*vsize;
      volarr = (double *) realloc(volarr,sizeof(double) * vsize);
    }
    volarr[n] = ivol;
    n = n + 1;

    for (i = 0; i<5; i++)
      fgets(line,sizeof(line),in);
    sscanf(line,"%lf %lf %lf",&r0x,&r1x,&drx);
    nrx = (int) ((r1x-r0x)/drx + 1.1); /* Really: 1+round((r1-r0)/dr) */
    for (i = 0; i<nrx; i++)
      fgets(line,sizeof(line),in);
  }
  fclose(in);

  if (n == 0) {
    fprintf(stderr,"@%s:%d: Invalid potin file \'%s\', no volume records.\n",
            __FILE__,__LINE__,potin_file);
    exit(1);
  }

  if (false) {
    printf("Before sort:\n");
    for (int i = 0; i<n; i++)
      printf("%3d :: %.3f%s",i,volarr[i], (i%5==4) ? "\n" : "  ");
    printf("\n\n");
  }
  qsort(volarr,n,sizeof(double),cmp_double);

  if (false) {
    printf("After sort:\n");
    for (int i = 0; i<n; i++)
      printf("%3d :: %.3f%s",i,volarr[i], (i%5==4) ? "\n" : "  ");
    printf("\n\n");
  }

  nvol[0] = n;
  vol0[0] = volarr[n/2];
  x0[0] = pow(volarr[0]/vol0[0],1.0/3.0);
  x1[0] = pow(volarr[n-1]/vol0[0],1.0/3.0);

  free(volarr);
}

void potdata::readpot(const char *parmin_file,const char *potin_file,const double vol) {
  FILE *in;
  double x0,x1,dx,dr;
  int nx;

  double r0x,r1x,drx;
  int nrx;

  char metalx[80];
  int ipotx,modex; double pnx;
  double vol0;

  double *vatab,*vbtab,*vctab,*vdtab,*vetab,*p1tab,*altab,*vpairtab = nullptr;
  double *r0rwstab,*evol0tab;
  double (*C)[4];
  double *y,*dy;
  double x,dxdv;
  double unused;

  double zval,rws,ivol,r0rws,rcrws,rmrws;//,mass

  int i,j;
  int L;

  char line[1024];

  input_vol = vol;

  /* Read potential data */
  in = fopen(parmin_file,"r");
  do {
    fgets(line,sizeof(line),in);
  } while (line[strspn(line," \t")] == '#');

  /* Test to see whether this is a one-line or two-line version of parmin */
  if (sscanf(line,"%lf %lf %lf %lf %d",&ddl[1],&ddl[2],&ddl[3],&ddl[4],&L) == 5) {
    /* One-line version, call getparmindata to figure out volume table spacing. */
    int nvol;
    getparmindata(potin_file,&nvol,&vol0,&x0,&x1);
    dx = (x1-x0)/(nvol-1);
    if (false) {
      printf("getparmindata() ==> nvol = %d, vol0 = %.6f, x0= %.6f, x1 = %.6f, dx = %.6f\n",
             nvol,vol0,x0,x1,dx);
    }
  } else {
    /* Two-line version, reparse this line, and read second line */
    sscanf(line,"%lf %lf %lf %lf",&x0,&x1,&dx,&vol0);
    fgets(line,sizeof(line),in);
    sscanf(line,"%lf %lf %lf %lf %d",&ddl[1],&ddl[2],&ddl[3],&ddl[4],&L);

  }
  double rws_scale = pow(3.0*vol0/(16.0*atan(1.0)),1.0/3.0);
  fclose(in);

  lang = L+1;
  lmax = 2*L+1;
  double s = ddl[1],p = ddl[2],d = ddl[3],f = ddl[4];
  double ss = s*s, pp = p*p, dd = d*d, ff = f*f;
  anorm3 =  s*ss + 2.0*( p*pp +  d*dd +  f*ff);
  anorm4 = ss*ss + 2.0*(pp*pp + dd*dd + ff*ff);
  /*
    for (i = 1; i<=lmax; i++) {
      for (j = 1; j<=lmax; j++)
        del0.m[i][j] = 0.0;
      del0[i][i] = 1.0;
    }
    Matrix::sz = lmax;
  */
  nx = (int) ((x1-x0)/dx + 1.1); /* Really: 1+round((x1-x0)/dx) */
  vatab = new double[nx];
  vbtab = new double[nx];
  vctab = new double[nx];
  vdtab = new double[nx];
  vetab = new double[nx];

  p1tab = new double[nx];
  altab = new double[nx];

  r0rwstab = new double[nx];
  evol0tab = new double[nx];

  in = fopen(potin_file,"r");

  int *tag = new int[nx];
  for (i = 0; i<nx; i++) tag[i] = 0;

  int ii;
  for (ii = 0; ii<nx; ii++) {

    do {
      fgets(line,sizeof(line),in);
    } while (line[strspn(line," \t")] == '#');

    metalx[0] = 0;

    /* Read element type, mode, and pn parameter */ {
      int nf = sscanf(line,"%s %d %d %lf",metalx,&ipotx,&modex,&pnx);
      if (nf < 3) {
        printf("Error in %s() @ %s:%d: Inconsistency in potential input file (%s) "
               "at record %d:\n"
               "  Expected at least three fields. Number of fields = %d\n",
               __func__,__FILE__,__LINE__,potin_file,ii,
               nf);
        exit(1);
      }
      if (modex <= 4) {
        pnx = 1.0;
      } else if (modex <= 6) {
        if (nf != 4) {
          printf("Error in %s() @ %s:%d: Inconsistency in potential input file (%s) "
                 "at record %d:\n"
                 "  mode = %d, number of fields = %d\n",
                 __func__,__FILE__,__LINE__,potin_file,ii,
                 modex,nf);
          exit(1);
        }
      } else {
          printf("Error in %s() @ %s:%d: Inconsistency in potential input file (%s): "
                 "at record %d\n"
                 "  Invalid mode. mode = %d\n",
                 __func__,__FILE__,__LINE__,potin_file,ii,
                 modex);
      }
    }

    if (ii == 0) {
      sscanf(line,"%s %d %d %lf",metal,&ipot,&mode,&pn);
      if (modex <= 4) pn = pnx;
    } else {
      /* Check that {metal,ipot,mode}x == {metal,ipot,mode} */
      if (strcmp(metal,metalx) != 0 ||
         ipotx != ipot ||
         modex != mode ||
         pnx != pn) {
        printf("Error in %s() @ %s:%d: Inconsistency in potential input file (%s) "
               "at record %d:\n"
               "metalx != metal (%s != %s) or\n"
               "ipotx  != ipot  (%d != %d) or\n"
               "modex  != mode  (%d != %d) or\n"
               "pnx    != pn    (%.3f != %.3f).\n",
               __func__,__FILE__,__LINE__,potin_file,ii,
               metalx,metal,
               ipotx,ipot,
               modex,mode,
               pnx,pn);
        exit(1);
      }
    }
    //printf("LINE: %s\n",line);
    //printf("metal = \'%s\'  ipot = %d  mode = %d\n",metalx,ipotx,modex);
    fgets(line,sizeof(line),in);
    sscanf(line,"%lf %lf %lf %lf",&zval,&ivol,&rws,&mass);
    /*{
      double xi = x0 + i/((double) (nx-1)) * (x1-x0);
      double volguess = vol0 * xi*xi*xi;
      if (fabs(volguess/ivol - 1.0) > 1e-3)
        printf("Wrong volume guess,  i=%d  volgues=%15.5e  ivol=%15.5e\n",
               i,volguess,ivol);
     }*/

    double ifrac = (pow(ivol/vol0,1.0/3.0) - x0)/((x1-x0)/(nx-1));
    i = (int) (ifrac + 0.1);
    if (fabs(i - ifrac) > 0.01) {
      printf("Volume point not in table... ii=%d i=%d ifrac=%15.5e  vol=%15.5e\n",
             ii,i,ifrac,ivol);
      printf("vol0 = %15.5e  zval = %15.5e  mass = %15.5e\n",vol0,zval,mass);
      exit(1);
    } else if (tag[i] == 1) {
      printf("Duplicate volume point in table.... ii=%d i=%d ifrac=%15.5e  vol=%15.5e\n",
             ii,i,ifrac,ivol);
      exit(1);
    } else tag[i] = 1;

    fgets(line,sizeof(line),in);
    sscanf(line,"%lf %lf %lf %lf",&r0rwstab[i],&altab[i],&rcrws,&rmrws);
    fgets(line,sizeof(line),in);
    sscanf(line,"%lf %lf %lf",&p1tab[i],&unused,&evol0tab[i]);

    fgets(line,sizeof(line),in);
    sscanf(line,"%lf %lf %lf %lf %lf",
           &vatab[i],&vbtab[i],&vctab[i],&vdtab[i],&vetab[i]);
    if (ipot == 1) {
      vatab[i] *= vdtab[i];
      vctab[i] *= vctab[i];
      vetab[i] *= vetab[i];
    }

    fgets(line,sizeof(line),in);

    fgets(line,sizeof(line),in);
    sscanf(line,"%lf %lf %lf",&r0x,&r1x,&drx);
    nrx = (int) ((r1x-r0x)/drx + 1.1); /* Really: 1+round((r1-r0)/dr) */

    if (ii == 0) {
      r0 = r0x; r1 = r1x; dr = drx; nr = nrx;
      vpairtab = new double[nx*nr];
    } else {
      /* Check that {r0,r1,dr,nr}x == {r0,r1,dr,nr} */
    }

    for (j = 0; j<nr; j++) {
      double rj,ktan,dvdvol;
      fgets(line,sizeof(line),in);
      sscanf(line,"%lf %lf %lf %lf",
             &rj,&vpairtab[i*nr+j],&ktan,&dvdvol);

      { /* Add screening and fl() part to pair energy table */

        double al = altab[i];
        double p1 = p1tab[i];

        int bscreen = (al > 0.0);

        double xi = x0 + i/((double) (nx-1)) * (x1-x0);
        double rws = rws_scale * xi;

        double r0rws = r0rwstab[i];
        double r00 = r0rws*rws,rp = 1.8*rws;
        if (bscreen == 0) r0rws = 10.0;
        double alp = al,alm = al;
        if (mode == 2 || mode == 4 || mode == 6) alm = 125.0;
        al = alp;

        double r = r0 + j*(r1-r0)/(nr-1);

        double rrws = r/rws;
        //double rsqr = r*r;
        // double fl(double r,int mode,double rp,double p1,double al,double r0)
        double flr = fl(r,mode,rp,p1,al,r00,pn);
        double fl2 = flr*flr;
        double v2a = vatab[i]*fl2*fl2;
        double v2b = vbtab[i]*fl2;
        double fscr = 1.0;

        if (bscreen == 1 && rrws >= r0rws) {
          double arg = rrws/r0rwstab[i];
          double arg1 = arg - 1.0;
          double arg12 = arg1*arg1;
          double f,dp;
          if (mode <= 2) {
            f = fgauss(arg,al);
            dp=2.*al*arg*arg1;
          }
          else {
            f = hgauss(arg,al);
            double arg13 = arg1*arg12;
            dp=2.0*al*al*arg*arg13/(1.+al*arg12);
          }
          fscr = f*f;
        }

        double vpair_tmp = vpairtab[i*nr+j];
        vpairtab[i*nr+j] = vpairtab[i*nr+j]*fscr + v2a - v2b;

        if (false) if (fabs(vol-ivol) < 0.01) {
          static FILE *xfile = nullptr;
          if (j == 0) {
            xfile = fopen("mgpt5-pot.dat","w");
            fprintf(xfile,"%%%%  vol = %15.5e  ivol = %15.5e i = %d  ii = %d\n",
                    vol,ivol,i,ii);
          }
          fprintf(xfile,"%15.5e %15.5e %15.5e %15.5e %15.5e %20.10e\n",
                  r,vpair_tmp,fscr,v2a,v2b,flr);
          if (j == nr-1) fclose(xfile);
        }


      }

    }
  }
  fclose(in);

  for (i = 0; i<nx; i++)
    if (tag[i] == 0) {
      printf("Volume point missing in table. i = %d\n",i);
      exit(1);
    }

  /* Make table */
  x = pow(vol/vol0,1.0/3.0);
  dxdv = 1.0/(3.0*vol0*x*x);

  C = new double[(nr > nx) ? nr : nx][4];
  makespline(nx,1,vatab,C);
  evalspline(nx-1,x0,x1,C,x,&va,&dva,&unused);
  dva *= dxdv;

  makespline(nx,1,vbtab,C);
  evalspline(nx-1,x0,x1,C,x,&vb,&dvb,&unused);
  dvb *= dxdv;

  makespline(nx,1,vctab,C);
  evalspline(nx-1,x0,x1,C,x,&vc,&dvc,&unused);
  dvc *= dxdv;

  makespline(nx,1,vdtab,C);
  evalspline(nx-1,x0,x1,C,x,&vd,&dvd,&unused);
  dvd *= dxdv;

  makespline(nx,1,vetab,C);
  evalspline(nx-1,x0,x1,C,x,&ve,&dve,&unused);
  dve *= dxdv;

  makespline(nx,1,p1tab,C);
  evalspline(nx-1,x0,x1,C,x,&p1,&dp1,&unused);
  dp1 *= dxdv;

  makespline(nx,1,altab,C);
  evalspline(nx-1,x0,x1,C,x,&al,&dal,&unused);
  dal *= dxdv;
  if (mode == 2 || mode == 4 || mode == 6) {
    al = 125.0;
    dal = 0.0;
  }


  {
    double dr0rws;
    makespline(nx,1,r0rwstab,C);
    evalspline(nx-1,x0,x1,C,x,&r0rws,&dr0rws,&unused);
    dr0rws *= dxdv;
    rws = rws_scale*x;
    r00 = r0rws * rws;
    dr00 = dr0rws*rws + r0rws*rws_scale*dxdv;
    rp = 1.8 * rws;
    drp = 1.8 * rws_scale*dxdv;
  }

  makespline(nx,1,evol0tab,C);
  evalspline(nx-1,x0,x1,C,x,&evol0,&devol0,&unused);
  devol0 *= dxdv;

  if (true) {
    printf("%% READPOT PARAMETERS:\n");

    printf("%% ddl = %15.5e  %15.5e  %15.5e  %15.5e\n",ddl[1],ddl[2],ddl[3],ddl[4]);
    printf("%% anorm3 = %15.5e  anorm4 = %15.5e\n",anorm3,anorm4);

    printf("%% x  = %15.5e    pn  = %15.5e\n",x,pn);
    printf("%% va = %15.5e    dva = %15.5e\n",va,dva);
    printf("%% vb = %15.5e    dvb = %15.5e\n",vb,dvb);
    printf("%% vc = %15.5e    dvc = %15.5e\n",vc,dvc);
    printf("%% vd = %15.5e    dvd = %15.5e\n",vd,dvd);
    printf("%% ve = %15.5e    dve = %15.5e\n",ve,dve);
    printf("%% p1 = %15.5e    dp1 = %15.5e\n",p1,dp1);
    printf("%% al = %15.5e    dal = %15.5e\n",al,dal);
    printf("%% rp = %15.5e    drp = %15.5e\n",rp,drp);
    printf("%% r00= %15.5e    dr00= %15.5e\n",r00,dr00);
    printf("\n");
  }

  y = new double[nr];
  dy = new double[nr];

  for (j = 0; j<nr; j++) {
    double d2y;
    makespline(nx,nr,&vpairtab[j],C);
    evalspline(nx-1,x0,x1,C,x,&y[j],&dy[j],&d2y);
    dy[j] *= dxdv;
  }
  vpair_spline = new double[nr-1][4];
  dvpair_spline = new double[nr-1][4];
  makespline(nr,1,y,vpair_spline);
  makespline(nr,1,dy,dvpair_spline);


  rcrit = rcrws * rws;
  rmax = rmrws * rws;

  delete[] dy;
  delete[] y;
  delete[] C;
  delete[] evol0tab;
  delete[] r0rwstab;
  delete[] altab;
  delete[] p1tab;
  delete[] vetab;
  delete[] vdtab;
  delete[] vctab;
  delete[] vbtab;
  delete[] vatab;
}

/*
int main(int argc,char *argv[]) {
  double vol = atof(argv[3]);
  int n = 25,i;

  printf("%% parmin = %s\n%% potin = %s\n%% vol = %15.5e\n",
         argv[1],argv[2],vol);

  readpot(argv[1],argv[2],vol);

  for (i = 0; i<n; i++) {
    double x,u,f,vir,dy,d2y;

    x = r0 + i*(r1-r0)/(n-1);
    evalspline(nr-1,r0,r1,vpair_spline,x,&u,&f,&d2y);
    evalspline(nr-1,r0,r1,dvpair_spline,x,&vir,&dy,&d2y);
    printf("  %15.5e  %15.5e  %15.5e  %15.5e\n",
           x,u,f,vir);
  }


  delete[] dvpair_spline;
  delete[] vpair_spline;

  return 0;
}
*/
