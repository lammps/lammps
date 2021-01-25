/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*  ----------------------------------------------------------------------
   Contributing author: Christopher Barrett (MSU) barrett@me.msstate.edu
    ----------------------------------------------------------------------*/


#include "fingerprint_radialscreened.h"
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>


using namespace LAMMPS_NS;

Fingerprint_radialscreened::Fingerprint_radialscreened(PairRANN *pair) : Fingerprint(pair)
{
	n_body_type = 2;
	dr = 0;
	re = 0;
	rc = 0;
	alpha = new double [1];
	alpha[0] = -1;
	n = 0;
	o = 0;
	id = -1;
	style = "radialscreened";
	atomtypes = new int [n_body_type];
	empty = true;
	fullydefined = false;
	pair->doscreen = true;
	screen = true;
}

Fingerprint_radialscreened::~Fingerprint_radialscreened()
{
	delete [] atomtypes;
	delete [] radialtable;
	delete [] alpha;
	delete [] dfctable;
	delete [] rinvsqrttable;
}

bool Fingerprint_radialscreened::parse_values(char* constant,char * line1){
	int l;
	char **words=new char *[MAXLINE];
	int nwords;
	nwords=0;
	words[nwords++] = strtok(line1,": ,\t\n");
	while ((words[nwords++] = strtok(NULL,": ,\t\n"))) continue;
	nwords -= 1;
	if (strcmp(constant,"re")==0){
		re = strtod(line1,NULL);
	}
	else if (strcmp(constant,"rc")==0){
		rc = strtod(line1,NULL);
	}
	else if (strcmp(constant,"alpha")==0){
		delete [] alpha;
		alpha = new double [nwords];
		for (l=0;l<nwords;l++){
			alpha[l]=strtod(words[l],NULL);
		}
	}
	else if (strcmp(constant,"dr")==0){
		dr = strtod(line1,NULL);
	}
	else if (strcmp(constant,"n")==0){
		n = strtol(line1,NULL,10);
	}
	else if (strcmp(constant,"o")==0){
		o = strtol(line1,NULL,10);
	}
	else pair->errorf("Undefined value for radial power");
	//code will run with default o=0 if o is never specified. All other values must be defined in potential file.
	delete [] words;
	if (re!=0 && rc!=0 && alpha!=0 && dr!=0 && n!=0)return true;
	return false;
}

void Fingerprint_radialscreened::write_values(FILE *fid){
	int i;
	fprintf(fid,"fingerprintconstants:");
	fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
	for (i=1;i<n_body_type;i++){
		fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
	}
	fprintf(fid,":%s_%d:re:\n",style,id);
	fprintf(fid,"%f\n",re);
	fprintf(fid,"fingerprintconstants:");
	fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
	for (i=1;i<n_body_type;i++){
		fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
	}
	fprintf(fid,":%s_%d:rc:\n",style,id);
	fprintf(fid,"%f\n",rc);
	fprintf(fid,"fingerprintconstants:");
	fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
	for (i=1;i<n_body_type;i++){
		fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
	}
	fprintf(fid,":%s_%d:alpha:\n",style,id);
	for (i=0;i<(n-o+1);i++){
		fprintf(fid,"%f ",alpha[i]);
	}
	fprintf(fid,"\n");
	fprintf(fid,"fingerprintconstants:");
	fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
	for (i=1;i<n_body_type;i++){
		fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
	}
	fprintf(fid,":%s_%d:dr:\n",style,id);
	fprintf(fid,"%f\n",dr);
	fprintf(fid,"fingerprintconstants:");
	fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
	for (i=1;i<n_body_type;i++){
		fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
	}
	fprintf(fid,":%s_%d:o:\n",style,id);
	fprintf(fid,"%d\n",o);
	fprintf(fid,"fingerprintconstants:");
	fprintf(fid,"%s",pair->elementsp[atomtypes[0]]);
	for (i=1;i<n_body_type;i++){
		fprintf(fid,"_%s",pair->elementsp[atomtypes[i]]);
	}
	fprintf(fid,":%s_%d:n:\n",style,id);
	fprintf(fid,"%d\n",n);
}

//called after fingerprint is fully defined and tables can be computed.
void Fingerprint_radialscreened::allocate()
{
	int k,m;
	double r1;
	int buf = 5;
	int res = pair->res;
	double cutmax = pair->cutmax;
	radialtable = new double [(res+buf)*get_length()];
	dfctable = new double [res+buf];
	for (k=0;k<(res+buf);k++){
		r1 = cutmax*cutmax*(double)(k)/(double)(res);
		for (m=0;m<=(n-o);m++){
			radialtable[k*(n-o+1)+m]=pow(sqrt(r1)/re,m+o)*exp(-alpha[m]*sqrt(r1)/re)*cutofffunction(sqrt(r1),rc,dr);
		}
		if (sqrt(r1)>=rc || sqrt(r1) <= (rc-dr)){
				dfctable[k]=0;
		}
		else{
		dfctable[k]=-8*pow(1-(rc-sqrt(r1))/dr,3)/dr/(1-pow(1-(rc-sqrt(r1))/dr,4));
		}
	}
	generate_rinvssqrttable();
}

//called after fingerprint is declared for i-j type, but before its parameters are read.
void Fingerprint_radialscreened::init(int *i,int id)
{
	empty = false;
	for (int j=0;j<n_body_type;j++){atomtypes[j] = i[j];}
	this->id = id;
}

void Fingerprint_radialscreened::compute_fingerprint(double * features,double * dfeaturesx,double *dfeaturesy,double *dfeaturesz,double *Sik, double *dSikx, double*dSiky, double *dSikz, double *dSijkx, double *dSijky, double *dSijkz, bool *Bij,int ii,int sid,double *xn,double *yn,double*zn,int *tn,int jnum,int *jl)
{
		int nelements = pair->nelements;
		int res = pair->res;
		int i,j,jj,itype,jtype,l,kk;
		double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
		int *ilist,*jlist,*numneigh,**firstneigh;
		//
		PairRANN::Simulation *sim = &pair->sims[sid];
		int count=0;
		double **x = sim->x;
		int *type = sim->type;
		ilist = sim->ilist;
		double cutmax = pair->cutmax;
		i = ilist[ii];
		itype = pair->map[type[i]];
		int f = pair->net[itype].dimensions[0];
		double cutinv2 = 1/cutmax/cutmax;
//		numneigh = sim->numneigh;
//		firstneigh = sim->firstneigh;
//		xtmp = x[i][0];
//		ytmp = x[i][1];
//		ztmp = x[i][2];
//		jlist = firstneigh[i];
//		jnum = numneigh[i];
		//loop over neighbors
		for (jj = 0; jj < jnum; jj++) {
			if (Bij[jj]==false){continue;}
//			j = jlist[jj];
//			j &= NEIGHMASK;
//			jtype = pair->map[type[j]];
			jtype = tn[jj];
			if (this->atomtypes[1] != nelements && this->atomtypes[1] != jtype)continue;
//			delx = xtmp - x[j][0];
//			dely = ytmp - x[j][1];
//			delz = ztmp - x[j][2];
			delx = xn[jj];
			dely = yn[jj];
			delz = zn[jj];
			rsq = delx*delx + dely*dely + delz*delz;
			if (rsq > rc*rc)continue;
			count = startingneuron;
			double r1 = (rsq*((double)res)*cutinv2);
			int m1 = (int)r1;
			if (m1>res || m1<1){pair->errorf("invalid neighbor radius!");}
			if (radialtable[m1]==0){continue;}
			//cubic interpolation from tables
			double *p1 = &radialtable[m1*(n-o+1)];
			double *p2 = &radialtable[(m1+1)*(n-o+1)];
			double *p3 = &radialtable[(m1+2)*(n-o+1)];
			double *p0 = &radialtable[(m1-1)*(n-o+1)];
			double *q = &dfctable[m1-1];
			double *rinvs = &rinvsqrttable[m1-1];
			r1 = r1-trunc(r1);
			double dfc = q[1] + 0.5 * r1*(q[2] - q[0] + r1*(2.0*q[0] - 5.0*q[1] + 4.0*q[2] - q[3] + r1*(3.0*(q[1] - q[2]) + q[3] - q[0])));
			double ri = rinvs[1] + 0.5 * r1*(rinvs[2] - rinvs[0] + r1*(2.0*rinvs[0] - 5.0*rinvs[1] + 4.0*rinvs[2] - rinvs[3] + r1*(3.0*(rinvs[1] - rinvs[2]) + rinvs[3] - rinvs[0])));
			for (l=0;l<=(n-o);l++){
				double rt = Sik[jj]*(p1[l]+0.5*r1*(p2[l]-p0[l]+r1*(2.0*p0[l]-5.0*p1[l]+4.0*p2[l]-p3[l]+r1*(3.0*(p1[l]-p2[l])+p3[l]-p0[l]))));
				features[count]+=rt;
				double rt1 = rt*((l+o)/rsq+(-alpha[l]/re+dfc)*ri);
				//update neighbor's features
				dfeaturesx[jj*f+count]+=rt1*delx+rt*dSikx[jj];
				dfeaturesy[jj*f+count]+=rt1*dely+rt*dSiky[jj];
				dfeaturesz[jj*f+count]+=rt1*delz+rt*dSikz[jj];
				for (kk=0;kk<jnum;kk++){
					if (Bij[kk]==false){continue;}
					dfeaturesx[kk*f+count]+=rt*dSijkx[jj*jnum+kk];
					dfeaturesy[kk*f+count]+=rt*dSijky[jj*jnum+kk];
					dfeaturesz[kk*f+count]+=rt*dSijkz[jj*jnum+kk];
				}
				count++;
			}
		}
		for (jj=0;jj<jnum;jj++){
			if (Bij[jj]==false){continue;}
			count = startingneuron;
			for (l=0;l<=(n-o);l++){
				dfeaturesx[jnum*f+count]-=dfeaturesx[jj*f+count];
				dfeaturesy[jnum*f+count]-=dfeaturesy[jj*f+count];
				dfeaturesz[jnum*f+count]-=dfeaturesz[jj*f+count];
				count++;
			}
		}

	}

int Fingerprint_radialscreened::get_length()
{
	return n-o+1;
}
