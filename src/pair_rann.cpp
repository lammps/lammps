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
   Contributing authors: Christopher Barrett (MSU) barrett@me.msstate.edu
   	   	   	   	   	     Doyl Dickel (MSU) doyl@cavs.msstate.edu
    ----------------------------------------------------------------------*/

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "atom.h"
#include "style_fingerprint.h"
#include "style_activation.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "pair_rann.h"



using namespace LAMMPS_NS;

PairRANN::PairRANN(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  allocated = 0;

  nelements = -1;
  elements = NULL;
  mass = NULL;

  // set comm size needed by this Pair
  // comm unused for now.

  comm_forward = 38;
  comm_reverse = 30;
  res = 10000;
  cutmax = 0;
  //at least one of the following will change during fingerprint definition:
  doscreen = false;
  allscreen = true;
  dospin = false;

	fingerprint_map = new FingerprintCreatorMap();

	#define FINGERPRINT_CLASS
	#define FingerprintStyle(key,Class) \
	  (*fingerprint_map)[#key] = &fingerprint_creator<Class>;
	#include "style_fingerprint.h"
	#undef FingerprintStyle
	#undef FINGERPRINT_CLASS

	activation_map = new ActivationCreatorMap();

	#define ACTIVATION_CLASS
	#define ActivationStyle(key,Class) \
	  (*activation_map)[#key] = &activation_creator<Class>;
	#include "style_activation.h"
	#undef ActivationStyle
	#undef ACTIVATION_CLASS
}

PairRANN::~PairRANN()
{
	//clear memory
	delete [] mass;
	for (int i=0;i<nelements;i++){delete [] elements[i];}
	delete [] elements;
	for (int i=0;i<nelementsp;i++){delete [] elementsp[i];}
	delete [] elementsp;
	for (int i=0;i<=nelements;i++){
		if (net[i].layers>0){
			for (int j=0;j<net[i].layers-1;j++){
				delete [] net[i].Weights[j];
				delete [] net[i].Biases[j];
			}
			delete [] net[i].dimensions;
			delete [] net[i].Weights;
			delete [] net[i].Biases;
		}
	}
	delete [] net;
	delete [] map;
	for (int i=0;i<nelementsp;i++){
		if (fingerprintlength[i]>0){
			delete [] fingerprints[i];
			delete [] activation[i];
		}
	}
	delete [] fingerprints;
	delete [] activation;
	delete [] fingerprintcount;
	delete [] fingerprintperelement;
	delete [] fingerprintlength;
	delete [] screening_min;
	delete [] screening_max;
}



void PairRANN::allocate(char **elementword)
{
	int i,j,k,l,n;
	n = atom->ntypes;
	memory->create(setflag,n+1,n+1,"pair:setflag");
	memory->create(cutsq,n+1,n+1,"pair:cutsq");
	cutmax = 0;
	nelementsp=nelements+1;
	//initialize arrays
	elements = new char *[nelements];
	elementsp = new char *[nelementsp];//elements + 'all'
	mass = new double[nelements];
	net = new NNarchitecture[nelementsp];
	weightdefined = new bool*[nelementsp];
	biasdefined = new bool *[nelementsp];
	activation = new Activation**[nelementsp];
	fingerprints = new Fingerprint**[nelementsp];
	fingerprintlength = new int[nelementsp];
	fingerprintperelement = new int [nelementsp];
	fingerprintcount = new int[nelementsp];
	screening_min = new double [nelements*nelements*nelements];
	screening_max = new double [nelements*nelements*nelements];
	for (i=0;i<nelements;i++){
		for (int j =0;j<nelements;j++){
			for (int k=0;k<nelements;k++){
				screening_min[i*nelements*nelements+j*nelements+k] = 0.8;//default values. Custom values may be read from potential file later.
				screening_max[i*nelements*nelements+j*nelements+k] = 2.8;//default values. Custom values may be read from potential file later.
			}
		}
	}
	for (i=0;i<=nelements;i++){
		n = strlen(elementword[i])+1;
		fingerprintlength[i]=0;
		fingerprintperelement[i] = -1;
		fingerprintcount[i] = 0;
		if (i<nelements){
			mass[i]=-1.0;
			elements[i]= new char[n];
			strcpy(elements[i],elementword[i]);
		}
		elementsp[i] = new char[n];
		strcpy(elementsp[i],elementword[i]);

		net[i].layers = 0;
		net[i].dimensions = new int[1];
		net[i].dimensions[0]=0;
	}
}

void PairRANN::settings(int narg, char **arg)
{
	//read pair_style command in input file
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

void PairRANN::coeff(int narg, char **arg)
{
	int i,j;
	map = new int [atom->ntypes+1];


	  if (narg != 3 + atom->ntypes)
	    error->all(FLERR,"Incorrect args for pair coefficients");


	  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
	    error->all(FLERR,"Incorrect args for pair coefficients");
	  nelements = -1;
	  read_file(arg[2]);

	  // read args that map atom types to elements in potential file
	  // map[i] = which element the Ith atom type is, -1 if NULL

	  for (i = 3; i < narg; i++) {
	    if (strcmp(arg[i],"NULL") == 0) {
	      map[i-2] = -1;
	      continue;
	    }
	    for (j = 0; j < nelements; j++)
	    	{
	    	if (strcmp(arg[i],elements[j]) == 0) break;
	    	}
	    if (j < nelements) map[i-2] = j;
	    else error->all(FLERR,"No matching element in NN potential file");
	  }
	  // clear setflag since coeff() called once with I,J = * *

	  int n = atom->ntypes;
	  for (i = 1; i <= n; i++)
	    for (j = i; j <= n; j++)
	      setflag[i][j] = 0;

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
	  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
	  for (i=0;i<nelementsp;i++){
		  for (j=0;j<fingerprintperelement[i];j++){
			  fingerprints[i][j]->allocate();
		  }
	  }
	  allocated=1;
}

void PairRANN::read_file(char *filename)
{
	FILE *fp;
	int eof = 0,i,j,k,l;
	int n,nwords;
	char line [MAXLINE],line1[MAXLINE];

	char *ptr;
	bool comment;
	char str[128];
	fp = force->open_potential(filename);
	if (fp == NULL) {
	  sprintf(str,"Cannot open NN potential file %s",filename);
	  error->one(FLERR,str);
	}
	while (eof == 0){
		  ptr = fgets(line,MAXLINE,fp);
		  if (ptr == NULL) {
			  if (check_potential()) {//looks to see if everything needed appears to be defined
				  error->one(FLERR,"Invalid syntax in potential file, values are inconsistent or missing");
			  }
			  else{
				    update_stack_size();
					fclose(fp);
					eof = 1;
					break;
			  }
		  }
		  else n = strlen(line) + 1;
		// strip comment, skip line if blank
		if ((ptr = strchr(line,'#'))) *ptr = '\0';
		nwords = count_words(line);
		char **words = new char* [strlen(line)];
		if (nwords == 0) continue;
		comment = true;
		while (comment==true){
			ptr = fgets(line1,MAXLINE,fp);
			if (ptr==NULL)error->one(FLERR,"Unexpected end of file");
			if ((ptr = strchr(line1,'#'))) *ptr = '\0';
			nwords = count_words(line1);
			if (nwords == 0) continue;
			comment = false;
		}
		nwords = 0;
		words[nwords++] = strtok(line,": ,\t_\n");
		while ((words[nwords++] = strtok(NULL,": ,\t_\n"))) continue;
		if (strcmp(words[0],"atomtypes")==0)read_atom_types(words,line1);
		else if (strcmp(words[0],"mass")==0)read_mass(words,line1);
		else if (strcmp(words[0],"fingerprintsperelement")==0)read_fpe(words,line1);
		else if (strcmp(words[0],"fingerprints")==0)read_fingerprints(words,nwords-1,line1);
		else if (strcmp(words[0],"fingerprintconstants")==0)read_fingerprint_constants(words,nwords-1,line1);
		else if (strcmp(words[0],"networklayers")==0)read_network_layers(words,line1);
		else if (strcmp(words[0],"layersize")==0)read_layer_size(words,line1);
		else if (strcmp(words[0],"weight")==0)read_weight(words,line1,fp);
		else if (strcmp(words[0],"bias")==0)read_bias(words,line1,fp);
		else if (strcmp(words[0],"activationfunctions")==0)read_activation_functions(words,line1);
		else if (strcmp(words[0],"calibrationparameters")==0)continue;//information on how the network was trained
		else if (strcmp(words[0],"screening")==0)read_screening(words,nwords-1,line1);
		else error->one(FLERR,"Could not understand file syntax: unknown keyword");
		delete [] words;
	}
}

void PairRANN::read_atom_types(char **words,char * line1){
	int nwords = 0;
	int t = count_words(line1)+1;
	char **elementword = new char *[t];
	elementword[nwords++] = strtok(line1," ,\t:_\n");
	while ((elementword[nwords++] = strtok(NULL," ,\t:_\n"))) continue;
	if (nwords < 1) errorf("Incorrect syntax for atom types");
	elementword[nwords-1] = new char [strlen("all")+1];
	char elt [] = "all";
	strcpy(elementword[nwords-1],elt);
	nelements = nwords-1;
	allocate(elementword);
}

void PairRANN::read_mass(char **words,char * line1){
	if (nelements == -1)error->one(FLERR,"atom types must be defined before mass in potential file.");
	int nwords = 0,i;
	for (i=0;i<nelements;i++){
		if (strcmp(words[1],elements[i])==0){
			mass[i]=strtod(line1,NULL);
			return;
		}
	}
	error->one(FLERR,"mass element not found in atom types.");
}

void PairRANN::read_fpe(char **words,char * line1){
	int i,j;
	if (nelements == -1)error->one(FLERR,"atom types must be defined before fingerprints per element in potential file.");
	for (i=0;i<nelementsp;i++){
		if (strcmp(words[1],elementsp[i])==0){
			fingerprintperelement[i] = strtol(line1,NULL,10);
			fingerprints[i] = new Fingerprint *[fingerprintperelement[i]];
			for (int j=0;j<fingerprintperelement[i];j++){
				fingerprints[i][j]=new Fingerprint(this);
			}
			return;
		}
	}
	error->one(FLERR,"fingerprint-per-element element not found in atom types");
}

void PairRANN::read_fingerprints(char **words,int nwords,char * line1){
	int nwords1=0,i,j,k,l,m,i1;
	bool found;
	char str[MAXLINE];
	char **words1 = new char * [count_words(line1)+1];
	words1[nwords1++] = strtok(line1," ,\t:_\n");
	while ((words1[nwords1++] = strtok(NULL," ,\t:_\n"))) continue;
	nwords1 -= 1;
	if (nelements == -1)error->one(FLERR,"atom types must be defined before fingerprints in potential file.");
	int atomtypes[nwords-1];
	for (i=1;i<nwords;i++){
		found = false;
		for (j=0;j<nelementsp;j++){
			if (strcmp(words[i],elementsp[j])==0){
				atomtypes[i-1]=j;
				found = true;
				break;
			}
		}
		if (!found){error->one(FLERR,"fingerprint element not found in atom types");}
	}
	i = atomtypes[0];
	k = 0;
	if (fingerprintperelement[i]==-1){error->one(FLERR,"fingerprint per element must be defined before fingerprints");}
	while (k<nwords1){
		i1 = fingerprintcount[i];
		delete fingerprints[i][i1];
		fingerprints[i][i1] = create_fingerprint(words1[k]);
		sprintf(str,"%d %d\n",nwords-1,fingerprints[i][i1]->n_body_type);
		std::cout<<str;
		if (fingerprints[i][i1]->n_body_type!=nwords-1){error->one(FLERR,"invalid fingerprint for element combination");}
		k++;
		fingerprints[i][i1]->init(atomtypes,strtol(words1[k++],NULL,10));
		fingerprintcount[i]++;
	}
	delete [] words1;
}


void PairRANN::read_fingerprint_constants(char **words,int nwords,char * line1){
	int i,j,k,l,m,i1;
	bool found;
	char str [128];
	if (nelements == -1)error->one(FLERR,"atom types must be defined before fingerprints in potential file.");
	int n_body_type = nwords-4;
	int atomtypes[n_body_type];
	for (i=1;i<=n_body_type;i++){
		found = false;
		for (j=0;j<nelementsp;j++){
			if (strcmp(words[i],elementsp[j])==0){
				atomtypes[i-1]=j;
				found = true;
				break;
			}
		}
		if (!found){error->one(FLERR,"fingerprint element not found in atom types");}
	}
	i = atomtypes[0];
	found = false;
	for (k=0;k<fingerprintperelement[i];k++){
		if (fingerprints[i][k]->empty){continue;}
		if (n_body_type!=fingerprints[i][k]->n_body_type){continue;}
		for (j=0;j<n_body_type;j++){
			if (fingerprints[i][k]->atomtypes[j]!=atomtypes[j]){break;}
			if (j==n_body_type-1){
				if (strcmp(words[nwords-3],fingerprints[i][k]->style)==0 && strtol(words[nwords-2],NULL,10)==fingerprints[i][k]->id){
					found=true;
					i1 = k;
					break;
				}
			}
		}
		if (found){break;}
	}
	if (!found){error->one(FLERR,"cannot define constants for unknown fingerprint");}
	fingerprints[i][i1]->fullydefined=fingerprints[i][i1]->parse_values(words[nwords-1],line1);
}

void PairRANN::read_network_layers(char **words,char *line1){
	int i,j;
	if (nelements == -1)error->one(FLERR,"atom types must be defined before network layers in potential file.");
	for (i=0;i<nelements;i++){
		if (strcmp(words[1],elements[i])==0){
			net[i].layers = strtol(line1,NULL,10);
			if (net[i].layers < 1)error->one(FLERR,"invalid number of network layers");
			delete [] net[i].dimensions;
			net[i].dimensions = new int [net[i].layers];
			weightdefined[i] = new bool [net[i].layers];
			biasdefined[i] = new bool [net[i].layers];
			net[i].Weights = new double * [net[i].layers-1];
			net[i].Biases = new double * [net[i].layers-1];
			for (j=0;j<net[i].layers;j++){
				net[i].dimensions[j]=0;
				weightdefined[i][j] = false;
				biasdefined[i][j] = false;
			}
			activation[i]=new Activation* [net[i].layers-1];
			for (int j=0;j<net[i].layers-1;j++){
				activation[i][j]= new Activation(this);
			}
			return;
		}
	}
	error->one(FLERR,"network layers element not found in atom types");
}

void PairRANN::read_layer_size(char **words,char* line1){
	int i;
	for (i=0;i<nelements;i++){
		if (strcmp(words[1],elements[i])==0){
			if (net[i].layers==0)error->one(FLERR,"networklayers for each atom type must be defined before the corresponding layer sizes.");
			int j = strtol(words[2],NULL,10);
			if (j>=net[i].layers || j<0){errorf("invalid layer in layer size definition");};
			net[i].dimensions[j]= strtol(line1,NULL,10);
//			net[i].Weights[j] = new double [1];
//			net[i].Weights[j][0]=0;
//			net[i].Biases[j] = new double [1];
//			net[i].Biases[j][0] = 0;
			return;
		}
	}
	error->one(FLERR,"layer size element not found in atom types");
}

void PairRANN::read_weight(char **words,char* line1,FILE* fp){
	int i,j,k,l,nwords;
	char *ptr;
	char **words1;
	for (l=0;l<nelements;l++){
		if (strcmp(words[1],elements[l])==0){
			if (net[l].layers==0)error->one(FLERR,"networklayers must be defined before weights.");
			i=strtol(words[2],NULL,10);
			if (i>=net[l].layers || i<0)error->one(FLERR,"invalid weight layer");
			if (net[l].dimensions[i]==0 || net[l].dimensions[i+1]==0) error->one(FLERR,"network layer sizes must be defined before corresponding weight");
//			delete [] net[l].Weights[i];
			net[l].Weights[i] = new double [net[l].dimensions[i]*net[l].dimensions[i+1]];
			weightdefined[l][i] = true;
			int n = count_words(line1)+1;
			words1 = new char* [n];
			nwords=0;
			words1[nwords++] = strtok(line1," ,\t:_\n");
			while ((words1[nwords++] = strtok(NULL," ,\t:_\n"))) continue;
			nwords -= 1;
			if (nwords != net[l].dimensions[i])error->one(FLERR,"invalid weights per line");
			for (k=0;k<net[l].dimensions[i];k++){
				net[l].Weights[i][k] = strtod(words1[k],NULL);
			}
			for (j=1;j<net[l].dimensions[i+1];j++){
				ptr = fgets(line1,MAXLINE,fp);
				nwords=0;
				words1[nwords++] = strtok(line1," ,\t:_\n");
				while ((words1[nwords++] = strtok(NULL," ,\t:_\n"))) continue;
				nwords -= 1;
				if (nwords != net[l].dimensions[i])error->one(FLERR,"invalid weights per line");
				for (k=0;k<net[l].dimensions[i];k++){
					net[l].Weights[i][j*net[l].dimensions[i]+k] = strtod(words1[k],NULL);
				}
			}
			delete [] words1;
			return;
		}
	}
	error->one(FLERR,"weight element not found in atom types");
}

void PairRANN::read_bias(char **words,char* line1,FILE* fp){
	int i,j,l,nwords;
	char *ptr;
	for (l=0;l<nelements;l++){
		if (strcmp(words[1],elements[l])==0){
			if (net[l].layers==0)error->one(FLERR,"networklayers must be defined before biases.");
			i=strtol(words[2],NULL,10);
			if (i>=net[l].layers || i<0)error->one(FLERR,"invalid bias layer");
			if (net[l].dimensions[i]==0) error->one(FLERR,"network layer sizes must be defined before corresponding bias");
//			delete [] net[l].Biases[i];
			biasdefined[l][i] = true;
			net[l].Biases[i] = new double [net[l].dimensions[i+1]];
			words[0] = strtok(line1," ,\t:_\n");
			net[l].Biases[i][0] = strtod(words[0],NULL);
			for (j=1;j<net[l].dimensions[i+1];j++){
				ptr = fgets(line1,MAXLINE,fp);
				words[0] = strtok(line1," ,\t:_\n");
				net[l].Biases[i][j] = strtod(words[0],NULL);
			}
			return;
		}
	}
	error->one(FLERR,"bias element not found in atom types");
}

void PairRANN::read_activation_functions(char** words,char * line1){
	int i,j,l,nwords;
	int *ptr;
	for (l=0;l<nelements;l++){
		if (strcmp(words[1],elements[l])==0){
			if (net[l].layers==0)error->one(FLERR,"networklayers must be defined before activation functions.");
			i = strtol(words[2],NULL,10);
			if (i>=net[l].layers || i<0)error->one(FLERR,"invalid activation layer");
			nwords=0;
			words[nwords++] = strtok(line1," ,\t:_\n");
			delete activation[l][i];
			activation[l][i]=create_activation(line1);
			return;
		}
	}
	error->one(FLERR,"activation function element not found in atom types");
}

void PairRANN::read_screening(char** words,int nwords,char *line1){
	int i,j,k;
	bool found;
//	char str[MAXLINE];
//	sprintf(str,"%d\n",nwords);
//	for (i=0;i<nwords;i++){
//		std::cout<<words[i];
//		std::cout<<"\n";
//	}
//	std::cout<<str;
	if (nelements == -1)errorf("atom types must be defined before fingerprints in potential file.");
	if (nwords!=5)errorf("invalid screening command");
	int n_body_type = 3;
	int atomtypes[n_body_type];
	for (i=1;i<=n_body_type;i++){
		found = false;
		for (j=0;j<nelementsp;j++){
			if (strcmp(words[i],elementsp[j])==0){
				atomtypes[i-1]=j;
				found = true;
				break;
			}
		}
		if (!found){errorf("fingerprint element not found in atom types");}
	}
	i = atomtypes[0];
	j = atomtypes[1];
	k = atomtypes[2];
	int index = i*nelements*nelements+j*nelements+k;
	int index1 = i*nelements*nelements+k*nelements+j;
	if (strcmp(words[4],"Cmin")==0)	{
		screening_min[index] = strtod(line1,NULL);
		screening_min[index1] = screening_min[index];
	}
	else if (strcmp(words[4],"Cmax")==0) {
		screening_max[index] = strtod(line1,NULL);
		screening_max[index1] = screening_max[index];
	}
	else errorf("unrecognized screening keyword");
}

//Called after finishing reading the potential file to make sure it is complete. True is bad.
//also allocates maxlayer and fingerprintlength.
bool PairRANN::check_potential(){
  int i,j,k,l;
  if (nelements==-1){return true;}
  for (i=0;i<=nelements;i++){
	  if (i<nelements){
		  if (mass[i]<0)return true;//uninitialized mass
	  }
	  if (net[i].layers==0)break;//no definitions for this starting element, not considered an error.
	  net[i].maxlayer=0;
	  for (j=0;j<net[i].layers;j++){
		  if (net[i].dimensions[j]==0)return true;//incomplete network definition
		  if (net[i].dimensions[j]>net[i].maxlayer)net[i].maxlayer = net[i].dimensions[j];
	  }
	  if (net[i].dimensions[net[i].layers-1]!=1)return true;//output layer must have single neuron (the energy)
	  for (j=0;j<net[i].layers-1;j++){
		  if (!weightdefined[i][j])return true;//undefined weights
		  if (!biasdefined[i][j])return true;//undefined biases
		  if (activation[i][j]->empty)return true;//undefined activations
		  for (k=0;k<net[i].dimensions[j+1];k++){
			  for (l=0;l<net[i].dimensions[j];l++){
				  if (net[i].Weights[j][k*net[i].dimensions[j]+l]==0)return true;//undefined weights
			  }
			  if (net[i].Biases[j][k]==0)return true;//undefined biases
		  }
	  }
	  for (j=0;j<fingerprintperelement[i];j++){
		  if (fingerprints[i][j]->fullydefined==false)return true;
		  fingerprints[i][j]->startingneuron = fingerprintlength[i];
		  fingerprintlength[i] +=fingerprints[i][j]->get_length();
		  if (fingerprints[i][j]->rc>cutmax){cutmax = fingerprints[i][j]->rc;}
	  }
	  if (net[i].dimensions[0]!=fingerprintlength[i])return true;
  }
  return false;//everything looks good
}

void PairRANN::compute(int eflag, int vflag)
{
		//perform force/energy computation_
	if (dospin){
	  if (strcmp(update->unit_style,"metal") != 0)
	    error->all(FLERR,"Spin pair styles require metal units");
	  if (!atom->sp_flag)
	      error->all(FLERR,"Spin pair styles requires atom/spin style");
	}
	if (eflag || vflag) ev_setup(eflag,vflag);
	else evflag = vflag_fdotr = vflag_atom = 0;
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
	if (dospin){
	  sims->s = atom->sp;
	}
	int itype,f,jnum,len;
	if (eflag || vflag) ev_setup(eflag,vflag);
	else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;
	if (eflag_global){eng_vdwl=0;eng_coul=0;}
	double energy=0;
	double **force = atom->f;
	double **fm = atom->fm;
	double **virial = vatom;
	char str[MAXLINE];
	//loop over atoms
	for (ii=0;ii<sims->inum;ii++){
		  i = sims->ilist[ii];
		  itype = map[sims->type[i]];
		  f = net[itype].dimensions[0];
		  jnum = sims->numneigh[i];
		  double xn[jnum];
		  double yn[jnum];
		  double zn[jnum];
		  int tn[jnum];
		  int jl[jnum];
		  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
		  double features [f];
		  double dfeaturesx[f*jnum];
		  double dfeaturesy[f*jnum];
		  double dfeaturesz[f*jnum];
		  for (j=0;j<f;j++){
			  features[j]=0;
		  }
		  for (j=0;j<f*jnum;j++){
			  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
		  }
		  //screening is calculated once for all atoms if any fingerprint uses it.
		  double Sik[jnum];
		  double dSikx[jnum];
		  double dSiky[jnum];
		  double dSikz[jnum];
		  double dSijkx[jnum*jnum];
		  double dSijky[jnum*jnum];
		  double dSijkz[jnum*jnum];
		  bool Bij[jnum];
		  double sx[f*jnum];
		  double sy[f*jnum];
		  double sz[f*jnum];
		  if (dospin){
			  for (j=0;j<f*jnum;j++){
				  sx[j]=sy[j]=sz[j]=0;
			  }
		  }
	      clock_t t1 = clock();
		  if (doscreen){
				screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
		  }
		  if (allscreen){
			  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
		  }
		  clock_t t2 = clock();
		  //do fingerprints for atom type
		  len = fingerprintperelement[itype];
		  for (j=0;j<len;j++){
					   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		  }
		  itype = nelements;
		  //do fingerprints for type "all"
		  len = fingerprintperelement[itype];
		  for (j=0;j<len;j++){
			  	  	   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  	  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  	  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  	  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		  }
		  clock_t t3 = clock();
		  //run fingerprints through network
		  if (dospin){
			  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&energy,force,fm,virial,ii,jnum,jl);
		  }
		  else {
			  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&energy,force,virial,ii,jnum,jl);
		  }
		  //testdfeatures(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii);
		  clock_t t4 = clock();
		  double ts = (double) (t2-t1) / CLOCKS_PER_SEC * 1000.0;
		  double tf = (double) (t3-t2) / CLOCKS_PER_SEC * 1000.0;
		  double tp = (double) (t4-t3) / CLOCKS_PER_SEC * 1000.0;
		  sprintf(str,"screen time: %f, fingerprint time: %f, propagation time: %f\n",ts,tf,tp);
//		  std::cout<<str;
	}
//	testdenergy();
	if (vflag_fdotr) virial_fdotr_compute();
}

void PairRANN::cull_neighbor_list(double *xn,double *yn, double *zn,int *tn, int* jnum,int *jl,int i,int sn){
	int *jlist,j,count,jj,*type,jtype;
	double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
	double **x = sims[sn].x;
	xtmp = x[i][0];
	ytmp = x[i][1];
	ztmp = x[i][2];
	type = sims[sn].type;
	jlist = sims[sn].firstneigh[i];
	count = 0;
	for (jj=0;jj<jnum[0];jj++){
		j = jlist[jj];
		j &= NEIGHMASK;
		jtype = map[type[j]];
		delx = xtmp - x[j][0];
		dely = ytmp - x[j][1];
		delz = ztmp - x[j][2];
		rsq = delx*delx + dely*dely + delz*delz;
		if (rsq>cutmax*cutmax){
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

void PairRANN::screen_neighbor_list(double *xn,double *yn, double *zn,int *tn, int* jnum,int *jl,int i,int sn,bool *Bij,double *Sik, double *dSikx, double*dSiky, double *dSikz, double *dSijkx, double *dSijky, double *dSijkz){
	double xnc[jnum[0]],ync[jnum[0]],znc[jnum[0]];
	double Sikc[jnum[0]];
	double dSikxc[jnum[0]];
	double dSikyc[jnum[0]];
	double dSikzc[jnum[0]];
	double dSijkxc[jnum[0]][jnum[0]];
	double dSijkyc[jnum[0]][jnum[0]];
	double dSijkzc[jnum[0]][jnum[0]];
	int jj,kk,count,count1,tnc[jnum[0]],jlc[jnum[0]];
	count = 0;
	for (jj=0;jj<jnum[0]-1;jj++){
		if (Bij[jj]){
			count1 = 0;
			xnc[count]=xn[jj];
			ync[count]=yn[jj];
			znc[count]=zn[jj];
			tnc[count]=tn[jj];
			jlc[count]=jl[jj];
			Sikc[count]=Sik[jj];
			dSikxc[count]=dSikx[jj];
			dSikyc[count]=dSiky[jj];
			dSikzc[count]=dSikz[jj];
			for (kk=0;kk<jnum[0]-1;kk++){
				if (Bij[kk]){
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
	for (jj=0;jj<count;jj++){
		xn[jj]=xnc[jj];
		yn[jj]=ync[jj];
		zn[jj]=znc[jj];
		tn[jj]=tnc[jj];
		jl[jj]=jlc[jj];
		Bij[jj] = true;
		Sik[jj]=Sikc[jj];
		dSikx[jj]=dSikxc[jj];
		dSiky[jj]=dSikyc[jj];
		dSikz[jj]=dSikzc[jj];
		for (kk=0;kk<count;kk++){
			dSijkx[jj*count+kk] = dSijkxc[jj][kk];
			dSijky[jj*count+kk] = dSijkyc[jj][kk];
			dSijkz[jj*count+kk] = dSijkzc[jj][kk];
		}
	}
}

void PairRANN::testdfeatures(double *features,double *dfeaturesx,double *dfeaturesy,double *dfeaturesz, double *dspinx,double *dspiny, double *dspinz,int ii){
	char str[MAXLINE];
	char str1[128];
	double del = 0.00005;
	int i,j,k,jj,inum,jnum,itype,jtype,ktype,l,kk,m,p,jj1;
	double xtmp,ytmp,ztmp,delx,dely,delz,rsq,rp,rbp,delx1,dely1,delz1,rsq1;
	int *ilist,*jlist,*numneigh,**firstneigh,len;
	int maxindex = nelementsp*nelementsp-1;
	//

	int nn = 0;
	int count=0;
	int *type = atom->type;
	inum = sims->inum;
	ilist = sims->ilist;
	int f = net[itype].dimensions[0];
	double cutinv2 = 1/cutmax/cutmax;
	i = sims->ilist[ii];
	numneigh = sims->numneigh;
	firstneigh = sims->firstneigh;
	jlist = firstneigh[i];
	itype = map[sims->type[i]];
	f = net[itype].dimensions[0];
	jnum = sims->numneigh[i];
	double xn[jnum];
	double yn[jnum];
	double zn[jnum];
	int tn[jnum];
	int jl[jnum];
	double features1x [f];
	double features2x [f];
	double features1y [f];
	double features2y [f];
	double features1z [f];
	double features2z [f];
	double features1sx [f];
	double features2sx [f];
	double features1sy [f];
	double features2sy [f];
	double features1sz [f];
	double features2sz [f];
	double dfeaturesx1[f*jnum];
	double dfeaturesy1[f*jnum];
	double dfeaturesz1[f*jnum];
	double dspinx1[f*jnum];
	double dspiny1[f*jnum];
	double dspinz1[f*jnum];
	double dfxtest[f*jnum];
	double dfytest[f*jnum];
	double dfztest[f*jnum];
	double dfstestx[f*jnum];
	double dfstesty[f*jnum];
	double dfstestz[f*jnum];
	for (j=0;j<f*jnum;j++){
	  dfeaturesx1[j]=dfeaturesy1[j]=dfeaturesz1[j]=0;
	}
	//screening is calculated once for all atoms if any fingerprint uses it.
	double Sik[jnum];
	double dSikx[jnum];
	double dSiky[jnum];
	double dSikz[jnum];
	double dSijkx[jnum*jnum];
	double dSijky[jnum*jnum];
	double dSijkz[jnum*jnum];
	bool Bij[jnum];
	double sx[jnum*f];
	double sy[jnum*f];
	double sz[jnum*f];
	itype = map[type[i]];
	cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
	//do fingerprints for atom type
	len = fingerprintperelement[itype];
	if (doscreen){
		screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
	}
	if (allscreen){
	  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
	}
	for (jj=0;jj<jnum-1;jj++){
		for (l=0;l<f;l++){
		  features1x[l]=0;
		  features2x[l]=0;
		  features1y[l]=0;
		  features2y[l]=0;
		  features1z[l]=0;
		  features2z[l]=0;
		  features1sx[l]=0;
		  features2sx[l]=0;
		  features1sy[l]=0;
		  features2sy[l]=0;
		  features1sz[l]=0;
		  features2sz[l]=0;
		}
		for (l=0;l<f*jnum;l++){
		  dfeaturesx1[l]=dfeaturesy1[l]=dfeaturesz1[l]=sx[l]=sy[l]=sz[l]=0;
		}
		xn[jj] = xn[jj]+del;
		if (doscreen){
			screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
		}
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1x,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1x,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1x,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1x,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		xn[jj] = xn[jj]-2*del;
		if (doscreen){
			screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
		}
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2x,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2x,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2x,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2x,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		for (k=0;k<f;k++){
			dfxtest[jj*f+k] = (features1x[k]-features2x[k])/2/del;
		}
		xn[jj] = xn[jj]+del;
		yn[jj] = yn[jj]+del;
		if (doscreen){
			screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
		}
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1y,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1y,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1y,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1y,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		yn[jj] = yn[jj]-2*del;
		if (doscreen){
			screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
		}
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2y,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2y,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2y,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2y,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		for (k=0;k<f;k++){
			dfytest[jj*f+k] = (features1y[k]-features2y[k])/2/del;
		}
		yn[jj] = yn[jj]+del;
		zn[jj] = zn[jj]+del;
		if (doscreen){
			screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
		}
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1z,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1z,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1z,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1z,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		zn[jj] = zn[jj]-2*del;
		if (doscreen){
			screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
		}
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2z,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2z,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2z,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2z,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		for (k=0;k<f;k++){
			dfztest[jj*f+k] = (features1z[k]-features2z[k])/2/del;
		}
		zn[jj] = zn[jj]+del;
		jj1 = jl[jj];
		sims->s[jj1][0]+=del;
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1sx,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1sx,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1sx,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1sx,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		sims->s[jj1][0]-=2*del;
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2sx,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2sx,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2sx,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2sx,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		for (k=0;k<f;k++){
			dfstestx[jj*f+k] = (features1sx[k]-features2sx[k])/2/del;
		}
		sims->s[jj1][0]+=del;
		sims->s[jj1][1]+=del;
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1sy,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1sy,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1sy,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1sy,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		sims->s[jj1][1]-=2*del;
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2sy,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2sy,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2sy,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2sy,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		for (k=0;k<f;k++){
			dfstesty[jj*f+k] = (features1sy[k]-features2sy[k])/2/del;
		}
		sims->s[jj1][1]+=del;
		sims->s[jj1][2]+=del;
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1sz,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1sz,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features1sz,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features1sz,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		sims->s[jj1][2]-=2*del;
		for (j=0;j<len;j++){
				   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2sz,dfeaturesx1,dfeaturesy1,dfeaturesz1,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2sz,dfeaturesx1,dfeaturesy1,dfeaturesz1,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features2sz,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features2sz,dfeaturesx1,dfeaturesy1,dfeaturesz1,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
		}
		for (k=0;k<f;k++){
			dfstestz[jj*f+k] = (features1sz[k]-features2sz[k])/2/del;
		}
		sims->s[jj1][2]+=del;
		for (k=0;k<f;k++){
			sprintf(str,"ii: %d, jj: %d, k: %d, f:%.10f, fx1: %.10f, fx2: %.10f, dfxa: %.10f, dfxn: %.10f, dfya: %.10f, dfyn: %.10f, dfza: %.10f, dfzn: %.10f,dfsxa: %.10f, dfsxn: %.10f, dfsya: %.10f, dfsyn: %.10f, dfsza: %.10f, dfszn: %.10f\n",ii,jj,k,features[k],features1x[k],features2x[k],dfeaturesx[jj*f+k],dfxtest[jj*f+k],dfeaturesy[jj*f+k],dfytest[jj*f+k],dfeaturesz[jj*f+k],dfztest[jj*f+k],dspinx[jj*f+k],dfstestx[jj*f+k],dspiny[jj*f+k],dfstesty[jj*f+k],dspinz[jj*f+k],dfstestz[jj*f+k]);
			std::cout<<str;
		}
		sprintf(str,"got here %d %d %d %d\n",ii,f,jnum,jj);
		std::cout<<str;
	}
	if (ii==inum-1){
		error->one(FLERR,"terminate");
	}
}

void PairRANN::testdenergy(){
	double **force = atom->f;
	double **fm = atom->fm;
	double del = 0.00002;
	int ii,i,j;
	int nn = 0;
	double ex1,ex2,ey1,ey2,ez1,ez2,esx1,esx2,esy1,esy2,esz1,esz2;
	ex1=ex2=ey1=ey2=ez1=ez2=esx1=esx2=esy1=esy2=esz1=esz2=0.0;
	double **force1 = new double *[listfull->maxatom];
	double **fm1 = new double *[listfull->maxatom];
	double **force1n = new double *[listfull->maxatom];
	double **fm1n = new double *[listfull->maxatom];
	for (i=0;i<listfull->maxatom;i++){
		force1[i]=new double [3];
		fm1[i] = new double [3];
		force1n[i]=new double [3];
		fm1n[i] = new double [3];
	}
	for (int n=0;n<sims->inum;n++){
		int itype,f,jnum,len;
		double **virial = vatom;
		char str[MAXLINE];
		eng_vdwl=0;
		//loop over atoms
		sims->x[n][0]+=del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
				  	  	   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  	  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  	  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  	  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&ex1,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&ex1,force1,virial,ii,jnum,jl);
			  }
			  ex1=eng_vdwl;
		}
		eng_vdwl=0;
		//loop over atoms
		sims->x[n][0]-=2*del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&ex2,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&ex2,force1,virial,ii,jnum,jl);
			  }
			  ex2=eng_vdwl;
		}
		force1n[n][0]=(ex1-ex2)/2/del;
		eng_vdwl=0;
		//loop over atoms
		sims->x[n][0]+=del;
		sims->x[n][1]+=del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&ey1,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&ey1,force1,virial,ii,jnum,jl);
			  }
			  ey1=eng_vdwl;
		}
		eng_vdwl=0;
		//loop over atoms
		sims->x[n][1]-=2*del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&ey2,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&ey2,force1,virial,ii,jnum,jl);
			  }
			  ey2=eng_vdwl;
		}
		eng_vdwl=0;
		force1n[n][1]=(ey1-ey2)/2/del;
		//loop over atoms
		sims->x[n][1]+=del;
		sims->x[n][2]+=del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&ez1,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&ez1,force1,virial,ii,jnum,jl);
			  }
			  ez1=eng_vdwl;
		}
		eng_vdwl=0;
		//loop over atoms
		sims->x[n][2]-=2*del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&ez2,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&ez2,force1,virial,ii,jnum,jl);
			  }
			  ez2=eng_vdwl;
		}
		eng_vdwl=0;
		sims->x[n][2]+=del;
		force1n[n][2]=(ez1-ez2)/2/del;

		//loop over atoms
		sims->s[n][0]+=del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
				  	  	   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  	  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  	  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
				  	  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&esx1,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&esx1,force1,virial,ii,jnum,jl);
			  }
			  esx1=eng_vdwl;
		}
		eng_vdwl=0;
		//loop over atoms
		sims->s[n][0]-=2*del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&esx2,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&esx2,force1,virial,ii,jnum,jl);
			  }
			  esx2=eng_vdwl;
		}
		eng_vdwl=0;
		fm1n[n][0]=(esx1-esx2)/2/del;
		//loop over atoms
		sims->s[n][0]+=del;
		sims->s[n][1]+=del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&esy1,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&esy1,force1,virial,ii,jnum,jl);
			  }
			  esy1=eng_vdwl;
		}
		eng_vdwl=0;
		//loop over atoms
		sims->s[n][1]-=2*del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&esy2,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&esy2,force1,virial,ii,jnum,jl);
			  }
			  esy2=eng_vdwl;
		}
		eng_vdwl=0;
		fm1n[n][1]=(esy1-esy2)/2/del;
		//loop over atoms
		sims->x[n][1]+=del;
		sims->x[n][2]+=del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&esz1,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&esz1,force1,virial,ii,jnum,jl);
			  }
			  esz1=eng_vdwl;
		}
		eng_vdwl=0;
		//loop over atoms
		sims->s[n][2]-=2*del;
		for (ii=0;ii<sims->inum;ii++){
			  i = sims->ilist[ii];
			  itype = map[sims->type[i]];
			  f = net[itype].dimensions[0];
			  jnum = sims->numneigh[i];
			  double xn[jnum];
			  double yn[jnum];
			  double zn[jnum];
			  int tn[jnum];
			  int jl[jnum];
			  cull_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0);
			  double features [f];
			  double dfeaturesx[f*jnum];
			  double dfeaturesy[f*jnum];
			  double dfeaturesz[f*jnum];
			  for (j=0;j<f;j++){
				  features[j]=0;
			  }
			  for (j=0;j<f*jnum;j++){
				  dfeaturesx[j]=dfeaturesy[j]=dfeaturesz[j]=0;
			  }
			  //screening is calculated once for all atoms if any fingerprint uses it.
			  double Sik[jnum];
			  double dSikx[jnum];
			  double dSiky[jnum];
			  double dSikz[jnum];
			  double dSijkx[jnum*jnum];
			  double dSijky[jnum*jnum];
			  double dSijkz[jnum*jnum];
			  bool Bij[jnum];
			  double sx[f*jnum];
			  double sy[f*jnum];
			  double sz[f*jnum];
			  if (dospin){
				  for (j=0;j<f*jnum;j++){
					  sx[j]=sy[j]=sz[j]=0;
				  }
			  }
			  if (doscreen){
					screen(Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,0,xn,yn,zn,tn,jnum-1);
			  }
			  if (allscreen){
				  screen_neighbor_list(xn,yn,zn,tn,&jnum,jl,i,0,Bij,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz);
			  }
			  //do fingerprints for atom type
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  itype = nelements;
			  //do fingerprints for type "all"
			  len = fingerprintperelement[itype];
			  for (j=0;j<len;j++){
						   if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==false && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==false)fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,ii,nn,xn,yn,zn,tn,jnum-1,jl);
					  else if (fingerprints[itype][j]->spin==true  && fingerprints[itype][j]->screen==true) fingerprints[itype][j]->compute_fingerprint(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,Sik,dSikx,dSiky,dSikz,dSijkx,dSijky,dSijkz,Bij,ii,nn,xn,yn,zn,tn,jnum-1,jl);
			  }
			  //run fingerprints through network
			  if (dospin){
				  propagateforwardspin(features,dfeaturesx,dfeaturesy,dfeaturesz,sx,sy,sz,&esz2,force1,fm1,virial,ii,jnum,jl);
			  }
			  else {
				  propagateforward(features,dfeaturesx,dfeaturesy,dfeaturesz,&esz2,force1,virial,ii,jnum,jl);
			  }
			  esz2=eng_vdwl;
		}
		eng_vdwl=0;
		sims->s[n][2]+=del;
		fm1n[n][2]=(esz1-esz2)/2/del;
		sprintf(str,"atom: %d fx: %f fxn: %f fy: %f fyn: %f fz: %f fzn: %f fmx: %f fmxn: %f fmy: %f fmyn: %f fmz: %f fmzn: %f\n",n,force[n][0],force1n[n][0],force[n][1],force1n[n][1],force[n][2],force1n[n][2],fm[n][0],fm1n[n][0],fm[n][1],fm1n[n][1],fm[n][2],fm1n[n][2]);
		std::cout<<str;
	}
}

void PairRANN::screen(double *Sik, double *dSikx, double*dSiky, double *dSikz, double *dSijkx, double *dSijky, double *dSijkz, bool *Bij, int ii,int sid,double *xn,double *yn,double *zn,int *tn,int jnum)
{
	//see Baskes, Materials Chemistry and Physics 50 (1997) 152-1.58
	int i,*jlist,jj,j,kk,k,itype,jtype,ktype;
	double Sijk,Cijk,Cn,Cd,Dij,Dik,Djk,C,dfc,dC,**x;
	PairRANN::Simulation *sim = &sims[sid];
//	x = sim->x;
	double xtmp,ytmp,ztmp,delx,dely,delz,rij,delx2,dely2,delz2,rik,delx3,dely3,delz3,rjk;
	i = sim->ilist[ii];
	itype = map[sim->type[i]];
//	jnum = sim->numneigh[i];
//	jlist = sim->firstneigh[i];
//	xtmp = x[i][0];
//	ytmp = x[i][1];
//	ztmp = x[i][2];
	for (int jj=0;jj<jnum;jj++){
		Sik[jj]=1;
		Bij[jj]=true;
		dSikx[jj]=0;
		dSiky[jj]=0;
		dSikz[jj]=0;
		for (kk=0;kk<jnum;kk++){
			dSijkx[jj*jnum+kk]=0;
			dSijky[jj*jnum+kk]=0;
			dSijkz[jj*jnum+kk]=0;
		}
	}
	for (kk=0;kk<jnum;kk++){//outer sum over k in accordance with source, some others reorder to outer sum over jj
		if (Bij[kk]==false){continue;}
//		k = jlist[kk];
//		k &= NEIGHMASK;
//		ktype = map[sim->type[k]];
		ktype = tn[kk];
//		delx2 = xtmp - x[k][0];
//		dely2 = ytmp - x[k][1];
//		delz2 = ztmp - x[k][2];
		delx2 = xn[kk];
		dely2 = yn[kk];
		delz2 = zn[kk];
		rik = delx2*delx2+dely2*dely2+delz2*delz2;
		if (rik>cutmax*cutmax){
			Bij[kk]= false;
			continue;
		}
		for (jj=0;jj<jnum;jj++){
			if (jj==kk){continue;}
			if (Bij[jj]==false){continue;}
//			j = jlist[jj];
//			j &= NEIGHMASK;
//			jtype = map[sim->type[j]];
//			delx = xtmp - x[j][0];
//			dely = ytmp - x[j][1];
//			delz = ztmp - x[j][2];
			jtype = tn[jj];
			delx = xn[jj];
			dely = yn[jj];
			delz = zn[jj];
			rij = delx*delx+dely*dely+delz*delz;
			if (rij>cutmax*cutmax){
				Bij[jj] = false;
				continue;
			}
//			delx3 = x[j][0]-x[k][0];
//			dely3 = x[j][1]-x[k][1];
//			delz3 = x[j][2]-x[k][2];
			delx3 = delx2-delx;
			dely3 = dely2-dely;
			delz3 = delz2-delz;
			rjk = delx3*delx3+dely3*dely3+delz3*delz3;
			if (rik+rjk<=rij){continue;}//bond angle > 90 degrees
			if (rik+rij<=rjk){continue;}//bond angle > 90 degrees
			double Cmax = screening_max[itype*nelements*nelements+jtype*nelements+ktype];
			double Cmin = screening_min[itype*nelements*nelements+jtype*nelements+ktype];
			double temp1 = rij-rik+rjk;
			Cn = temp1*temp1-4*rij*rjk;
			//Cn = (rij-rik+rjk)*(rij-rik+rjk)-4*rij*rjk;
			temp1 = rij-rjk;
			Cd = temp1*temp1-rik*rik;
			//Cd = (rij-rjk)*(rij-rjk)-rik*rik;
			Cijk = Cn/Cd;
			//Cijk = 1+2*(rik*rij+rik*rjk-rik*rik)/(rik*rik-(rij-rjk)*(rij-rjk));
			C = (Cijk-Cmin)/(Cmax-Cmin);
			if (C>=1){continue;}
			else if (C<=0){
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
void PairRANN::propagateforward(double *features,double *dfeaturesx,double *dfeaturesy,double *dfeaturesz, double * energy,double **force,double **virial, int ii,int jnum,int *jl){
	int i,j,k,jj,j1,itype,i1;
	int *ilist,*numneigh;
	ilist = listfull->ilist;
	int inum = listfull->inum;
	int *type = atom->type;
	i1=ilist[ii];
	itype = map[type[i1]];
	NNarchitecture net1 = net[itype];
//	numneigh = listfull->numneigh;
//	jnum = numneigh[ilist[ii]]+1;//extra value on the end of the array is the self term.
//	firstneigh = listfull->firstneigh;
//	jlist = firstneigh[i1];
	int L = net1.layers-1;
	double layer[net1.maxlayer];
	double sum[net1.maxlayer];
	double sum1[net1.maxlayer];
	double dlayerx[jnum][net1.maxlayer];
	double dlayersumx[jnum][net1.maxlayer];
	double dlayery[jnum][net1.maxlayer];
	double dlayersumy[jnum][net1.maxlayer];
	double dlayerz[jnum][net1.maxlayer];
	double dlayersumz[jnum][net1.maxlayer];
	//energy output with forces from analytical derivatives
	double dsum1;
	int f = net1.dimensions[0];
	for (i=0;i<net1.layers-1;i++){
		for (j=0;j<net1.dimensions[i+1];j++){
			//energy forward propagation
			sum[j]=0;
			for (k=0;k<net1.dimensions[i];k++){
				if (i==0&&j==0){
					layer[k]=features[k];
				}
				sum[j] += net1.Weights[i][j*net1.dimensions[i]+k]*layer[k];
			}
			sum[j] += net1.Biases[i][j];
			dsum1 = activation[itype][i]->dactivation_function(sum[j]);
			sum[j] = activation[itype][i]->activation_function(sum[j]);
			if (i==L-1){
				energy[j] = sum[j];
				if (eflag_atom)eatom[i1]=sum[j];
				if (eflag_global){eng_vdwl +=sum[j];}
			}
			//force propagation
			for (jj=0;jj<jnum;jj++){
//				if (Bij[jj]==false)continue;
				dlayersumx[jj][j]=0;
				dlayersumy[jj][j]=0;
				dlayersumz[jj][j]=0;
				for (k=0;k<net1.dimensions[i];k++){
					if (i==0&&j==0){
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
				if (i==L-1 && jj < (jnum-1)){
//					j1 = jlist[jj];
//					j1 &= NEIGHMASK;
					j1 = jl[jj];
					force[j1][0]+=dlayersumx[jj][j];
					force[j1][1]+=dlayersumy[jj][j];
					force[j1][2]+=dlayersumz[jj][j];
				}
			}
			if (i==L-1){
				j1 = ilist[ii];
				jj = jnum-1;
				force[j1][0]+=dlayersumx[jj][j];
				force[j1][1]+=dlayersumy[jj][j];
				force[j1][2]+=dlayersumz[jj][j];
			}
		}
		//update values for next iteration
		for (j=0;j<net1.dimensions[i+1];j++){
			layer[j]=sum[j];
			for (jj=0;jj<jnum;jj++){
//				if (Bij[jj]==false)continue;
				dlayerx[jj][j] = dlayersumx[jj][j];
				dlayery[jj][j] = dlayersumy[jj][j];
				dlayerz[jj][j] = dlayersumz[jj][j];
			}
		}
	}
}

//Called by getproperties. Propagate features and dfeatures through network. Updates force and energy
void PairRANN::propagateforwardspin(double *features,double *dfeaturesx,double *dfeaturesy,double *dfeaturesz, double *sx, double *sy, double *sz, double * energy,double **force,double **fm,double **virial, int ii,int jnum,int *jl){
	int i,j,k,jj,j1,itype,i1;
	int *ilist,*numneigh;
	ilist = listfull->ilist;
	int inum = listfull->inum;
	int *type = atom->type;
	i1=ilist[ii];
	itype = map[type[i1]];
	NNarchitecture net1 = net[itype];
//	numneigh = listfull->numneigh;
//	jnum = numneigh[ilist[ii]]+1;//extra value on the end of the array is the self term.
//	firstneigh = listfull->firstneigh;
//	jlist = firstneigh[i1];
	int L = net1.layers-1;
	double layer[net1.maxlayer];
	double sum[net1.maxlayer];
	double sum1[net1.maxlayer];
	double dlayerx[jnum][net1.maxlayer];
	double dlayersumx[jnum][net1.maxlayer];
	double dlayery[jnum][net1.maxlayer];
	double dlayersumy[jnum][net1.maxlayer];
	double dlayerz[jnum][net1.maxlayer];
	double dlayersumz[jnum][net1.maxlayer];
	double dsx[jnum][net1.maxlayer];
	double dssumx[jnum][net1.maxlayer];
	double dsy[jnum][net1.maxlayer];
	double dssumy[jnum][net1.maxlayer];
	double dsz[jnum][net1.maxlayer];
	double dssumz[jnum][net1.maxlayer];
	//energy output with forces from analytical derivatives
	double dsum1;
	int f = net1.dimensions[0];
	for (i=0;i<net1.layers-1;i++){
		for (j=0;j<net1.dimensions[i+1];j++){
			//energy forward propagation
			sum[j]=0;
			for (k=0;k<net1.dimensions[i];k++){
				if (i==0&&j==0){
					layer[k]=features[k];
				}
				sum[j] += net1.Weights[i][j*net1.dimensions[i]+k]*layer[k];
			}
			sum[j] += net1.Biases[i][j];
			dsum1 = activation[itype][i]->dactivation_function(sum[j]);
			sum[j] = activation[itype][i]->activation_function(sum[j]);
			if (i==L-1){
				energy[j] = sum[j];
				if (eflag_atom)eatom[i1]=sum[j];
				if (eflag_global){eng_vdwl +=sum[j];}
			}
			//force propagation
			for (jj=0;jj<jnum;jj++){
//				if (Bij[jj]==false)continue;
				dlayersumx[jj][j]=0;
				dlayersumy[jj][j]=0;
				dlayersumz[jj][j]=0;
				dssumx[jj][j]=0;
				dssumy[jj][j]=0;
				dssumz[jj][j]=0;
				for (k=0;k<net1.dimensions[i];k++){
					if (i==0&&j==0){
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
				if (i==L-1 && jj < (jnum-1)){
//					j1 = jlist[jj];
//					j1 &= NEIGHMASK;
					j1 = jl[jj];
					force[j1][0]+=dlayersumx[jj][j];
					force[j1][1]+=dlayersumy[jj][j];
					force[j1][2]+=dlayersumz[jj][j];
					fm[j1][0]+=dssumx[jj][j];
					fm[j1][1]+=dssumy[jj][j];
					fm[j1][2]+=dssumz[jj][j];
				}
			}
			if (i==L-1){
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
		for (j=0;j<net1.dimensions[i+1];j++){
			layer[j]=sum[j];
			for (jj=0;jj<jnum;jj++){
//				if (Bij[jj]==false)continue;
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

void PairRANN::update_stack_size(){
	//get very rough guess of memory usage

	//neighborlist memory use:
	memguess = 0;
	for (int i=0;i<nelementsp;i++){
		memguess+=8*net[i].dimensions[0]*20*3;
	}
	memguess+=8*20*12;
	memguess+=8*20*20*3;
	//generous buffer:
	memguess *= 16;
	const rlim_t kStackSize = memguess;
	struct rlimit rl;
	int result;
	result = getrlimit(RLIMIT_STACK, &rl);
	if (result == 0)
	{
		if (rl.rlim_cur < kStackSize)
		{
			rl.rlim_cur += kStackSize;
			result = setrlimit(RLIMIT_STACK, &rl);
			if (result != 0)
			{
				fprintf(stderr, "setrlimit returned result = %d\n", result);
			}
		}
	}
}

//treats # as starting a comment to be ignored.
int PairRANN::count_words(char *line){
	int n = strlen(line) + 1;
	char copy[n];
	strncpy(copy,line,n);
	char *ptr;
	if ((ptr = strchr(copy,'#'))) *ptr = '\0';
	if (strtok(copy," ,\t:_\n") == NULL) {
		return 0;
	}
	n=1;
	while ((strtok(NULL," ,\t:_\n"))) n++;
	return n;
}

void PairRANN::init_list(int which, NeighList *ptr)
{
  listfull = ptr;
}


void PairRANN::init_style()
{
	  int irequest_full = neighbor->request(this,instance_me);
	  neighbor->requests[irequest_full]->id = 1;
	  neighbor->requests[irequest_full]->half = 0;
	  neighbor->requests[irequest_full]->full = 1;
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairRANN::init_one(int i, int j)
{
  return cutmax;
}

void PairRANN::errorf(char * message){
	this->error->all(FLERR,message);
}

template <typename T>
Fingerprint *PairRANN::fingerprint_creator(PairRANN* pair)
{
  return new T(pair);
}

Fingerprint *PairRANN::create_fingerprint(const char *style)
{
	if (fingerprint_map->find(style) != fingerprint_map->end()) {
		FingerprintCreator fingerprint_creator = (*fingerprint_map)[style];
		return fingerprint_creator(this);
	}
	char str[128];
	sprintf(str,"Unknown fingerprint style %s",style);
	error->all(FLERR,str);
	return NULL;
}

template <typename T>
Activation *PairRANN::activation_creator(PairRANN* pair)
{
  return new T(pair);
}

Activation *PairRANN::create_activation(const char *style)
{
	if (activation_map->find(style) != activation_map->end()) {
		ActivationCreator activation_creator = (*activation_map)[style];
		return activation_creator(this);
	}
	char str[128];
	sprintf(str,"Unknown activation style %s",style);
	error->all(FLERR,str);
	return NULL;
}

