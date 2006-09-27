/* replicate tool

   read LAMMPS data file as input
   create new LAMMPS data file as output

   new file can be translated, dilated, mapped to periodic box
     and/or replicated to create a larger system
   create new molecule tags as needed

   Syntax:  replicate [options] < infile > outfile

     Options:

     -style atom_style

          use the LAMMPS atom style that corresponds to this file
	    e.g. atomic or bond or angle or full or granular
          will be used for reading and writing files
	  this option must be specified
	  for atom_style = hybrid, append subsequent styles
	    e.g. -style hybrid bond charge

     -repeat nx ny nz

          use the input data as a unit cell and replicate the unit
	    cell by (nx,ny,nz) in the +x, +y, +z directions
	  this creates new atoms, bonds, angles, etc
	  if not specified, nx = ny = nz = 1

     -inbox xlo xhi ylo yhi zlo zhi

          use these values for the bounding box of the input system
	    and for the unit cell used for replicating purposes
	  if not specified use the box bounds read in from infile

     -outbox xlo xhi ylo yhi zlo zhi

          use these values for the bounding box of the output system
	  if not specified the output bounding box is computed from the input
	    bounding box by holding the "lo" values constant and creating 
	    new "hi" values by multiplying by the replicating factors

     -unmap

          unmap all input atom positions using input image flags
            and inbox in a periodic sense
	  this enables molecular bonds to be replicated

     -map

          map all final atom positions into the output bounding box
	    in a periodic sense
	  works no matter how far away from the box the atoms are
	  this operation is independent of unmap option

     -imageout

          image flags will be written to outfile, else they won't be
	  if unmapping/mapping is done, image flags are modified accordingly

   Notes:
      
      works with "Velocities" entry by giving same velocity to
        all replicated images of an atom, which is probably not ideal

      can perform system dilation/contraction by using "-repeat 1 1 1"
        and specifying a slightly larger/smaller output box than input
        box.

      can add/delete image flags by using -imageout
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct box {
  double xlo,xhi,xsize;
  double ylo,yhi,ysize;
  double zlo,zhi,zsize;
};

#define MAXLINE 1000

main(int argc, char *argv[])

{
  char line[MAXLINE];        /* strings for reading/parsing input file */
  char linetoken[MAXLINE];
  char *token,*ptr;             
                             /* numbers of various quantities as read in */

  int natoms,nbonds,nangles,ndihedrals,nimpropers;
  int ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes;

  int nmolecules;            /* total # of mols specified by mol tags */

                             /* vectors of various quantities as read in */

  int *tag,*molecule,*type;
  double *q,*radius,*density,*x,*y,*z,*mux,*muy,*muz;
  int *imagex,*imagey,*imagez;
  int *vtag;
  double *vx,*vy,*vz,*vphix,*vphiy,*vphiz;
  int *btype,*atype,*dtype,*itype;
  int *bond1,*bond2,*angle1,*angle2,*angle3;
  int *dihed1,*dihed2,*dihed3,*dihed4,*impro1,*impro2,*impro3,*impro4;

  int atomstyle;             /* atom styles and settings */
  int style_angle,style_atomic,style_bond,style_charge,style_dipole;
  int style_dpd,style_full,style_granular,style_molecular;
  int charge,molecular,dipole;

  int nx,ny,nz;              /* replication factors */
  int inbox,outbox;          /* flags for whether in/out boxes were explicitly
				specified */
  int unmap;                 /* 0 = no unmapping, 1 = unmap from input box */
  int map;                   /* 0 = no mapping, 1 = map to output box */
  int imageout;              /* 0/1 = no/yes output of image flags */
  int imageflag;             /* 0/1 = no/yes if image flags exist in input */

  struct box in,out;         /* input and output bounding boxes */
  int ntotal;                /* total replication factor = nx*ny*nz */
  double newx,newy,newz;     /* replicated atom position */
  int size_atom_valid;       /* values that should be on atom line */
  int size_atom_actual;      /* values that are on atom line */

  int i,j,k,m,n,tmp,del,tag_offset,mol_offset;   /* temp variables */
  char **values;

  char *gettoken(char *, int);    /* function defs */
  int count_words(char *); 
  void copyline(int); 
  void scale(struct box, int, int, int, 
	     struct box, double *, double *, double *);
  void pbc(struct box, double *, double *, double *, int *, int *, int *);

  /* default input values */

  atomstyle = 0;
  style_angle = style_atomic = style_bond = style_charge = style_dipole = 
    style_dpd = style_full = style_granular = style_molecular = 0;
  nx = ny = nz = 1;
  inbox = outbox = 0;
  unmap = 0;
  map = 0;
  imageout = 0;
  
  /* read input options from command line
     should do more error checking for missing args */

  i = 1;
  while (i < argc) {

    if (!strcmp(argv[i],"-style")) {
      atomstyle = 1;
      if (strcmp(argv[i+1],"angle") == 0) style_angle = 1;
      else if (strcmp(argv[i+1],"atomic") == 0) style_atomic = 1;
      else if (strcmp(argv[i+1],"bond") == 0) style_bond = 1;
      else if (strcmp(argv[i+1],"charge") == 0) style_charge = 1;
      else if (strcmp(argv[i+1],"dipole") == 0) style_dipole = 1;
      else if (strcmp(argv[i+1],"dpd") == 0) style_dpd = 1;
      else if (strcmp(argv[i+1],"full") == 0) style_full = 1;
      else if (strcmp(argv[i+1],"granular") == 0) style_granular = 1;
      else if (strcmp(argv[i+1],"molecular") == 0) style_molecular = 1;
      else if (strcmp(argv[i+1],"hybrid") == 0) {
	n = i + 2;
	while (n < argc && argv[n][0] != '-') {
	  if (strcmp(argv[n],"angle") == 0) style_angle = 1;
	  else if (strcmp(argv[n],"atomic") == 0) style_atomic = 1;
	  else if (strcmp(argv[n],"bond") == 0) style_bond = 1;
	  else if (strcmp(argv[n],"charge") == 0) style_charge = 1;
	  else if (strcmp(argv[n],"dipole") == 0) style_dipole = 1;
	  else if (strcmp(argv[n],"dpd") == 0) style_dpd = 1;
	  else if (strcmp(argv[n],"full") == 0) style_full = 1;
	  else if (strcmp(argv[n],"granular") == 0) style_granular = 1;
	  else if (strcmp(argv[n],"molecular") == 0) style_molecular = 1;
	  else {
	    fprintf(stderr,"Error with command-line arg style\n");
	    exit(1);
	  }
	  n++;
	}
	if (i == n - 2) {
	  fprintf(stderr,"Error with command-line arg style\n");
	  exit(1);
	}
	i = n - 2;
      } else {
	fprintf(stderr,"Error with command-line arg style\n");
	exit(1);
      }
      i += 2;
    } else if (!strcmp(argv[i],"-repeat")) {
      sscanf(argv[i+1],"%d",&nx);
      sscanf(argv[i+2],"%d",&ny);
      sscanf(argv[i+3],"%d",&nz);
      i += 4;
    } else if (!strcmp(argv[i],"-inbox")) {
      inbox = 1;
      sscanf(argv[i+1],"%lg",&in.xlo);
      sscanf(argv[i+2],"%lg",&in.xhi);
      sscanf(argv[i+3],"%lg",&in.ylo);
      sscanf(argv[i+4],"%lg",&in.yhi);
      sscanf(argv[i+5],"%lg",&in.zlo);
      sscanf(argv[i+6],"%lg",&in.zhi);
      i += 7;
    } else if (!strcmp(argv[i],"-outbox")) {
      outbox = 1;
      sscanf(argv[i+1],"%lg",&out.xlo);
      sscanf(argv[i+2],"%lg",&out.xhi);
      sscanf(argv[i+3],"%lg",&out.ylo);
      sscanf(argv[i+4],"%lg",&out.yhi);
      sscanf(argv[i+5],"%lg",&out.zlo);
      sscanf(argv[i+6],"%lg",&out.zhi);
      i += 7;
    } else if (!strcmp(argv[i],"-unmap")) {
      unmap = 1;
      i += 1;
    } else if (!strcmp(argv[i],"-map")) {
      map = 1;
      i += 1;
    } else if (!strcmp(argv[i],"-imageout")) {
      imageout = 1;
      i += 1;
    } else {
      fprintf(stderr,"Syntax error: replicate [options] < infile > outfile\n");
      exit(1);
    }
  }

  /* set low-level info from style flags */

  charge = 0;
  if (style_charge || style_full || style_dipole) charge = 1;

  molecular = 0;
  if (style_angle || style_bond || style_full || style_molecular)
    molecular = 1;

  dipole = 0;
  if (style_dipole) dipole = 1;

  size_atom_valid = 5;
  if (charge) size_atom_valid += 1;
  if (style_dipole) size_atom_valid += 3;
  if (style_granular) size_atom_valid += 2;
  if (molecular) size_atom_valid += 1;

  /* error checks */

  if (atomstyle == 0) {
    fprintf(stderr,"ERROR: must use -style to set atom style\n");
    exit(1);
  }

  /* ntotal = total replication factor */

  ntotal = nx*ny*nz;

  /* defaults */

  natoms = nbonds = nangles = ndihedrals = nimpropers = 0;
  ntypes = nbondtypes = nangletypes = ndihedtypes = nimprotypes = 0;

  /* read/write header */

  copyline(1);   // echo 1st line

  while (1) {

    fgets(line,MAXLINE,stdin);

    // trim anything from '#' onward
    // if line is blank, continue

    if (ptr = strchr(line,'#')) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) {
      printf("\n");
      continue;
    }

    if (strstr(line,"atoms")) {
      sscanf(line,"%d",&natoms);
      printf("%d atoms\n",natoms*ntotal);
    } else if (strstr(line,"bonds")) {
      sscanf(line,"%d",&nbonds);
      printf("%d bonds\n",nbonds*ntotal);
    } else if (strstr(line,"angles")) {
      sscanf(line,"%d",&nangles);
      printf("%d angles\n",nangles*ntotal);
    } else if (strstr(line,"dihedrals")) {
      sscanf(line,"%d",&ndihedrals);
      printf("%d dihedrals\n",ndihedrals*ntotal);
    } else if (strstr(line,"impropers")) {
      sscanf(line,"%d",&nimpropers);
      printf("%d impropers\n",nimpropers*ntotal);
    } else if (strstr(line,"atom types")) {
      sscanf(line,"%d",&ntypes);
      printf("%d atom types\n",ntypes);
    } else if (strstr(line,"bond types")) {
      sscanf(line,"%d",&nbondtypes);
      printf("%d bond types\n",nbondtypes);
    } else if (strstr(line,"angle types")) {
      sscanf(line,"%d",&nangletypes);
      printf("%d angle types\n",nangletypes);
    } else if (strstr(line,"dihedral types")) {
      sscanf(line,"%d",&ndihedtypes);
      printf("%d dihedral types\n",ndihedtypes);
    } else if (strstr(line,"improper types")) {
      sscanf(line,"%d",&nimprotypes);
      printf("%d improper types\n",nimprotypes);
    } else if (strstr(line,"xlo xhi")) {
      if (!inbox) sscanf(line,"%lg %lg",&in.xlo,&in.xhi);
      in.xsize = in.xhi - in.xlo;
      if (!outbox) {
	out.xlo = in.xlo;
	if (nx > 1) out.xhi = out.xlo + (in.xhi-in.xlo)*nx;
	else out.xhi = in.xhi;
      }
      out.xsize = out.xhi - out.xlo;
      printf("%g %g xlo xhi\n",out.xlo,out.xhi);
    } else if (strstr(line,"ylo yhi")) {
      if (!inbox) sscanf(line,"%lg %lg",&in.ylo,&in.yhi);
      in.ysize = in.yhi - in.ylo;
      if (!outbox) {
	out.ylo = in.ylo;
	if (ny > 1) out.yhi = out.ylo + (in.yhi-in.ylo)*ny;
	else out.yhi = in.yhi;
      }
      out.ysize = out.yhi - out.ylo;
      printf("%g %g ylo yhi\n",out.ylo,out.yhi);
    } else if (strstr(line,"zlo zhi")) {
      if (!inbox) sscanf(line,"%lg %lg",&in.zlo,&in.zhi);
      in.zsize = in.zhi - in.zlo;
      if (!outbox) {
	out.zlo = in.zlo;
	if (nz > 1) out.zhi = out.zlo + (in.zhi-in.zlo)*nz;
	else out.zhi = in.zhi;
      }
      out.zsize = out.zhi - out.zlo;
      printf("%g %g zlo zhi\n",out.zlo,out.zhi);
    } else break;
  }

  /* malloc space for input atoms and topology */

  tag = (int *) malloc(natoms*sizeof(int));
  type = (int *) malloc(natoms*sizeof(int));
  x = (double *) malloc(natoms*sizeof(double));
  y = (double *) malloc(natoms*sizeof(double));
  z = (double *) malloc(natoms*sizeof(double));
  imagex = (int *) malloc(natoms*sizeof(int));
  imagey = (int *) malloc(natoms*sizeof(int));
  imagez = (int *) malloc(natoms*sizeof(int));

  vtag = (int *) malloc(natoms*sizeof(int));
  vx = (double *) malloc(natoms*sizeof(double));
  vy = (double *) malloc(natoms*sizeof(double));
  vz = (double *) malloc(natoms*sizeof(double));

  if (tag == NULL || type == NULL ||
      x == NULL || y == NULL || z == NULL ||
      imagex == NULL || imagey == NULL || imagez == NULL ||
      vtag == NULL || vx == NULL || vy == NULL || vz == NULL) {
    fprintf(stderr,"Error in atom malloc - no space\n");
    exit(1);
  }

  if (charge) {
    q = (double *) malloc(natoms*sizeof(double));
    if (q == NULL) {
      fprintf(stderr,"Error in atom malloc - no space\n");
      exit(1);
    }
  }

  if (molecular) {
    molecule = (int *) malloc(natoms*sizeof(int));
    if (molecule == NULL) {
      fprintf(stderr,"Error in atom malloc - no space\n");
      exit(1);
    }
  }

  if (dipole) {
    mux = (double *) malloc(natoms*sizeof(double));
    muy = (double *) malloc(natoms*sizeof(double));
    muz = (double *) malloc(natoms*sizeof(double));
    if (mux == NULL || muy == NULL || muz == NULL) {
      fprintf(stderr,"Error in atom malloc - no space\n");
      exit(1);
    }
  }

  if (style_granular) {
    radius = (double *) malloc(natoms*sizeof(double));
    density = (double *) malloc(natoms*sizeof(double));
    vphix = (double *) malloc(natoms*sizeof(double));
    vphiy = (double *) malloc(natoms*sizeof(double));
    vphiz = (double *) malloc(natoms*sizeof(double));
    if (radius == NULL || density == NULL || 
	vphix == NULL || vphiy == NULL || vphiz == NULL) {
      fprintf(stderr,"Error in atom malloc - no space\n");
      exit(1);
    }
  }

  if (nbonds) {
    btype = (int *) malloc(nbonds*sizeof(int));
    bond1 = (int *) malloc(nbonds*sizeof(int));
    bond2 = (int *) malloc(nbonds*sizeof(int));
    
    if (btype == NULL || bond1 == NULL || bond2 == NULL) {
      fprintf(stderr,"Error in bond malloc - no space\n");
      exit(1);
    }
  }

  if (nangles) {
    atype = (int *) malloc(nangles*sizeof(int));
    angle1 = (int *) malloc(nangles*sizeof(int));
    angle2 = (int *) malloc(nangles*sizeof(int));
    angle3 = (int *) malloc(nangles*sizeof(int));
    
    if (atype == NULL || angle1 == NULL || angle2 == NULL || angle3 == NULL) {
      fprintf(stderr,"Error in angle malloc - no space\n");
      exit(1);
    }
  }

  if (ndihedrals) {
    dtype = (int *) malloc(ndihedrals*sizeof(int));
    dihed1 = (int *) malloc(ndihedrals*sizeof(int));
    dihed2 = (int *) malloc(ndihedrals*sizeof(int));
    dihed3 = (int *) malloc(ndihedrals*sizeof(int));
    dihed4 = (int *) malloc(ndihedrals*sizeof(int));
    
    if (dtype == NULL || dihed1 == NULL || dihed2 == NULL || 
	dihed3 == NULL || dihed4 == NULL) {
      fprintf(stderr,"Error in dihedral malloc - no space\n");
      exit(1);
    }
  }

  if (nimpropers) {
    itype = (int *) malloc(nimpropers*sizeof(int));
    impro1 = (int *) malloc(nimpropers*sizeof(int));
    impro2 = (int *) malloc(nimpropers*sizeof(int));
    impro3 = (int *) malloc(nimpropers*sizeof(int));
    impro4 = (int *) malloc(nimpropers*sizeof(int));
    
    if (itype == NULL || impro1 == NULL || impro2 == NULL || 
	impro3 == NULL || impro4 == NULL) {
      fprintf(stderr,"Error in improper malloc - no space\n");
      exit(1);
    }
  }

  /* read identifier strings one by one in free-form part of data file */

  token = gettoken(line,1);

  while (token) {

    /* read atoms */

    if (!strcmp(token,"Atoms")) {

      for (i = 0; i < natoms; i++) {
	fgets(line,MAXLINE,stdin);

	if (i == 0) {
	  size_atom_actual = count_words(line);
	  if (size_atom_actual == size_atom_valid) imageflag = 0;
	  else if (size_atom_actual == size_atom_valid + 3) imageflag = 1;
	  else {
	    fprintf(stderr,"Error in atom lines\n");
	    exit(1);
	  }
	  if (unmap == 1 && imageflag == 0) {
	    fprintf(stderr,
		    "If unmap is requested, file must have image flags\n");
	    exit(1);
	  }
	  values = (char **) malloc(size_atom_actual*sizeof(char *));
	}

	values[0] = strtok(line," \t\n");
	for (m = 1; m < size_atom_actual; m++)
	  values[m] = strtok(NULL," \t\n");

	m = 0;
	tag[i] = atoi(values[m++]);
	if (molecular) molecule[i] = atoi(values[m++]);
	type[i] = atoi(values[m++]);
	if (charge) q[i] = atof(values[m++]);
	if (style_granular) {
	  radius[i] = 0.5 * atof(values[m++]);
	  density[i] = atof(values[m++]);
	}
	x[i] = atof(values[m++]);
	y[i] = atof(values[m++]);
	z[i] = atof(values[m++]);
	if (style_dipole) {
	  mux[i] = atof(values[m++]);
	  muy[i] = atof(values[m++]);
	  muz[i] = atof(values[m++]);
	}
	if (imageflag) {
	  imagex[i] = atoi(values[m++]);
	  imagey[i] = atoi(values[m++]);
	  imagez[i] = atoi(values[m++]);
	} else imagex[i] = imagey[i] = imagez[i] = 0;
      }
      
      /* unmap atoms if requested */

      if (unmap) {
	for (i = 0; i < natoms; i++) {
	  x[i] = x[i] + imagex[i]*in.xsize;
	  y[i] = y[i] + imagey[i]*in.ysize;
	  z[i] = z[i] + imagez[i]*in.zsize;
	  imagex[i] = imagey[i] = imagez[i] = 0;
	}
      }

      /* find number of molecules */

      nmolecules = 0;
      if (molecular)
	for (i = 0; i < natoms; i++)
	  if (molecule[i] > nmolecules) nmolecules = molecule[i];

      /* replicate set of N atoms as many times as requested
	 generate new tags, mol number, coords, image flags as needed */

      for (k = 0; k < nz; k++) {
	for (j = 0; j < ny; j++) {
	  for (i = 0; i < nx; i++) {
	    for (m = 0; m < natoms; m++) {
	      newx = x[m] + i*in.xsize;
	      newy = y[m] + j*in.ysize;
	      newz = z[m] + k*in.zsize;
	      if (outbox) scale(in,nx,ny,nz,out,&newx,&newy,&newz);
	      if (map) pbc(out,&newx,&newy,&newz,
			   &imagex[m],&imagey[m],&imagez[m]);
              tag_offset = k*ny*nx*natoms + j*nx*natoms + i*natoms;
              mol_offset = k*ny*nx*nmolecules + j*nx*nmolecules + i*nmolecules;

	      line[0] = '\0';
	      sprintf(&line[strlen(line)],"%d",tag[m]+tag_offset);
	      if (molecular) sprintf(&line[strlen(line)]," %d",
				     molecule[m]+mol_offset);
	      sprintf(&line[strlen(line)]," %d",type[m]);
	      if (charge) sprintf(&line[strlen(line)]," %g",q[m]);
	      if (style_granular) sprintf(&line[strlen(line)]," %g %g",
					  radius[m],density[m]);
	      sprintf(&line[strlen(line)]," %g %g %g",newx,newy,newz);
	      if (style_dipole) sprintf(&line[strlen(line)]," %g %g %g",
					mux[m],muy[m],muz[m]);
	      if (imageout) sprintf(&line[strlen(line)]," %d %d %d",
					imagex[m],imagey[m],imagez[m]);
	      printf("%s\n",line);
	    }
	  }
	}
      }

    }

    /* read velocities and replicate */
    
    else if (!strcmp(token,"Velocities")) {

      for (i = 0; i < natoms; i++) {
	fgets(line,MAXLINE,stdin);
	if (!style_granular) 
	  sscanf(line,"%d %lg %lg %lg",&vtag[i],&vx[i],&vy[i],&vz[i]);
	else 
	  sscanf(line,"%d %lg %lg %lg %lg %lg %lg",&vtag[i],
		 &vx[i],&vy[i],&vz[i],&vphix[i],&vphiy[i],&vphiz[i]);
      }
      
      for (k = 0; k < nz; k++) {
	for (j = 0; j < ny; j++) {
	  for (i = 0; i < nx; i++) {
	    for (m = 0; m < natoms; m++) {
              tag_offset = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	      if (!style_granular) 
		printf("%d %g %g %g\n",vtag[m]+tag_offset,vx[m],vy[m],vz[m]);
	      else 
		printf("%d %g %g %g %g %g %g\n",vtag[m]+tag_offset,
		       vx[m],vy[m],vz[m],vphix[m],vphiy[m],vphiz[m]);
	    }
	  }
	}
      }

    }

    /* read bonds and replicate */

    else if (!strcmp(token,"Bonds")) {

      for (i = 0; i < nbonds; i++) {
	fgets(line,MAXLINE,stdin);
	sscanf(line,"%d %d %d %d",&n,&btype[i],&bond1[i],&bond2[i]);
      }
    
      n = 0;
      for (k = 0; k < nz; k++)
	for (j = 0; j < ny; j++)
	  for (i = 0; i < nx; i++)
	    for (m = 0; m < nbonds; m++) {
	      n++;
	      del = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	      printf("%d %d %d %d\n",n,btype[m],bond1[m]+del,bond2[m]+del);
	    }
    }

    /* read angles and replicate */

    else if (!strcmp(token,"Angles")) {

      for (i = 0; i < nangles; i++) {
	fgets(line,MAXLINE,stdin);
	sscanf(line,"%d %d %d %d %d",&n,&atype[i],
	       &angle1[i],&angle2[i],&angle3[i]);
      }
    
      n = 0;
      for (k = 0; k < nz; k++)
	for (j = 0; j < ny; j++)
	  for (i = 0; i < nx; i++)
	    for (m = 0; m < nangles; m++) {
	      n++;
	      del = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	      printf("%d %d %d %d %d\n",n,atype[m],
		     angle1[m]+del,angle2[m]+del,angle3[m]+del);
	    }
    }
    
    /* read dihedrals and replicate */
    
    else if (!strcmp(token,"Dihedrals")) {

      for (i = 0; i < ndihedrals; i++) {
	fgets(line,MAXLINE,stdin);
	sscanf(line,"%d %d %d %d %d %d",&n,&dtype[i],&dihed1[i],&dihed2[i],
	       &dihed3[i],&dihed4[i]);
      }

      n = 0;
      for (k = 0; k < nz; k++)
	for (j = 0; j < ny; j++)
	  for (i = 0; i < nx; i++)
	    for (m = 0; m < ndihedrals; m++) {
	      n++;
	      del = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	      printf("%d %d %d %d %d %d\n",n,dtype[m],
		     dihed1[m]+del,dihed2[m]+del,dihed3[m]+del,dihed4[m]+del);
	    }
    }

    /* read impropers and replicate */

    else if (!strcmp(token,"Impropers")) {

      for (i = 0; i < nimpropers; i++) {
	fgets(line,MAXLINE,stdin);
	sscanf(line,"%d %d %d %d %d %d",&n,&itype[i],&impro1[i],&impro2[i],
	     &impro3[i],&impro4[i]);
      }

      n = 0;
      for (k = 0; k < nz; k++)
	for (j = 0; j < ny; j++)
	  for (i = 0; i < nx; i++)
	    for (m = 0; m < nimpropers; m++) {
	      n++;
	      del = k*ny*nx*natoms + j*nx*natoms + i*natoms;
	      printf("%d %d %d %d %d %d\n",n,itype[m],
		     impro1[m]+del,impro2[m]+del,impro3[m]+del,impro4[m]+del);
	    }
    }
    
    /* non-replicated sections
       just copy proper number of lines from in to out */

    else if (!strcmp(token,"Masses"))
      copyline(ntypes);
    else if (!strcmp(token,"Dipoles"))
      copyline(ntypes);
    else if (!strcmp(token,"Pair Coeffs"))
      copyline(ntypes);
    else if (!strcmp(token,"Bond Coeffs"))
      copyline(nbondtypes);
    else if (!strcmp(token,"Angle Coeffs"))
      copyline(nangletypes);
    else if (!strcmp(token,"Dihedral Coeffs"))
      copyline(ndihedtypes);
    else if (!strcmp(token,"Improper Coeffs"))
      copyline(nimprotypes);
    else if (!strcmp(token,"BondBond Coeffs"))
      copyline(nangletypes);
    else if (!strcmp(token,"BondAngle Coeffs"))
      copyline(nangletypes);
    else if (!strcmp(token,"MiddleBondTorsion Coeffs"))
      copyline(ndihedtypes);
    else if (!strcmp(token,"EndBondTorsion Coeffs"))
      copyline(ndihedtypes);
    else if (!strcmp(token,"AngleTorsion Coeffs"))
      copyline(ndihedtypes);
    else if (!strcmp(token,"AngleAngleTorsion Coeffs"))
      copyline(ndihedtypes);
    else if (!strcmp(token,"BondBond13 Coeffs"))
      copyline(ndihedtypes);
    else if (!strcmp(token,"AngleAngle Coeffs"))
      copyline(nimprotypes);
    else {
      fprintf(stderr,
	      "Error in input data file - unknown identifier %s\n",token);
      exit(1);
    }

    token = gettoken(line,0);
  }
  
}

/* ------------------------------------------------------------------- */

/* return a LAMMPS keyword from data file
   if first is 1, non-blank line with token is in line
   else read until find a non-blank line
   keyword is all text on line w/out leading & trailing white space
   read one additional line after non-blank line (assumed blank)
   return ptr to keyword
   return NULL if end of file */

char *gettoken(char *line, int first)
{
  char buffer[MAXLINE];

  // read upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  int eof = 0;

  if (!first) if (fgets(line,MAXLINE,stdin) == NULL) eof = 1;
  while (eof == 0 && strspn(line," \t\n\r") == strlen(line))
    if (fgets(line,MAXLINE,stdin) == NULL) eof = 1;
  if (fgets(buffer,MAXLINE,stdin) == NULL) eof = 1;

  // if eof, return NULL

  if (eof) return NULL;

  // bracket non-whitespace portion of line

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t' 
	 || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';

  // print keyword into output file

  if (first) printf("%s\n\n",&line[start]);
  else printf("\n%s\n\n",&line[start]);

  // return ptr to keyword

  return &line[start];
}

/* count and return words in a single line
   make copy of line before using strtok so as not to change line
   trim anything from '#' onward */

int count_words(char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) malloc(n*sizeof(char));
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n") == NULL) return 0;
  n = 1;
  while (strtok(NULL," \t\n")) n++;

  free(copy);
  return n;
}

/* copy n lines from stdin to stdout */

void copyline(int n)
{
  char line[1000];

  while (n) {
    fgets(line,MAXLINE,stdin);
    printf("%s",line);
    n--;
  }
}

/* point (x,y,z) is somewhere in (or near) the replicated box "in"
   rescale (x,y,z) so it is in the same relative location in box "out"
   do this by 1st converting to 0-1 coord system in "in" space */
   
void scale(struct box in, int nx, int ny, int nz,
	   struct box out, double *x, double *y, double *z)
{
  double rel;

  rel = (*x-in.xlo) / (nx*in.xsize);
  *x = out.xlo + rel*out.xsize;
  rel = (*y-in.ylo) / (ny*in.ysize);
  *y = out.ylo + rel*out.ysize;
  rel = (*z-in.zlo) / (nz*in.zsize);
  *z = out.zlo + rel*out.zsize;
}

/* remap the point (x,y,z) so it is guaranteed to be inside the box 
   using standard periodic imaging, update image flags */

void pbc(struct box box, double *x, double *y, double *z, 
	 int *imagex, int *imagey, int *imagez)
{
  while (*x < box.xlo) {
    *x += box.xsize;
    *imagex -= 1;
  }
  while (*x >= box.xhi) {
    *x -= box.xsize;
    *imagex += 1;
  }
  while (*y < box.ylo) {
    *y += box.ysize;
    *imagey -= 1;
  }
  while (*y >= box.yhi) {
    *y -= box.ysize;
    *imagey += 1;
  }
  while (*z < box.zlo) {
    *z += box.zsize;
    *imagez -= 1;
  }
  while (*z >= box.zhi) {
    *z -= box.zsize;
    *imagez += 1;
  }
}

