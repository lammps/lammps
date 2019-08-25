
#include "msi2lmp.h"

#include <stdlib.h>
#include <string.h>

static int count_bonds();
static int count_angles();
static int count_dihedrals();
static int count_oops();
static int count_angle_angles();

static void build_bonds_list();
static void build_angles_list();
static void build_dihedrals_list();
static void build_oops_list();
static void build_angleangles_list();

static void build_atomtypes_list();
static void build_bondtypes_list();
static void build_angletypes_list();
static void build_dihedraltypes_list();
static void build_ooptypes_list();
static void build_angleangletypes_list();

static void swap_ints(int *,int *);
static void bubble_sort(int, int *, int *);

void MakeLists()
{

  total_no_bonds = count_bonds();
  total_no_angles = count_angles();
  total_no_dihedrals = count_dihedrals();
  total_no_oops = count_oops();
  total_no_angle_angles = count_angle_angles();


  atomtypes = (struct AtomTypeList *)calloc(MAX_ATOM_TYPES,
                                            sizeof(struct AtomTypeList));
  if (atomtypes == NULL) {
    fprintf(stderr,"Trouble allocating memory for atomtypes list - Exiting\n");
    exit(1);
  }

  build_atomtypes_list();


  if (total_no_bonds > 0) {
    bonds = (struct BondList *)calloc(total_no_bonds,sizeof(struct BondList));
    if (bonds == NULL) {
      fprintf(stderr,"Trouble allocating memory for bonds list - Exiting\n");
      exit(1);
    }

    build_bonds_list();

    bondtypes = (struct BondTypeList *)calloc(MAX_BOND_TYPES,
                                              sizeof(struct BondTypeList));
    if (bondtypes == NULL) {
      fprintf(stderr,"Trouble allocating memory for bondtypes list - Exiting\n");
      exit(1);
    }

    build_bondtypes_list();
  }

  if (total_no_angles > 0) {
    angles = (struct AngleList *)calloc(total_no_angles,
                                        sizeof(struct AngleList));
    if (angles == NULL) {
      fprintf(stderr,"Trouble allocating memory for angles list - Exiting\n");
      exit(1);
    }

    build_angles_list();

    angletypes = (struct AngleTypeList *)calloc(MAX_ANGLE_TYPES,
                                                sizeof(struct AngleTypeList));
    if (angletypes == NULL) {
      fprintf(stderr,"Trouble allocating memory for angletypes list - Exiting\n");
      exit(1);
    }

    build_angletypes_list();
  }

  if (total_no_dihedrals > 0) {

    dihedrals = (struct DihedralList *)calloc(total_no_dihedrals,
                                              sizeof(struct DihedralList));
    if (dihedrals == NULL) {
      fprintf(stderr,"Trouble allocating memory for dihedrals list - Exiting\n");
      exit(1);
    }

    build_dihedrals_list();

    dihedraltypes = (struct DihedralTypeList *)calloc(MAX_DIHEDRAL_TYPES,
                                                      sizeof(struct DihedralTypeList));
    if (dihedraltypes == NULL) {
      fprintf(stderr,"Trouble allocating memory for dihedraltypes list - Exiting\n");
      exit(1);
    }

    build_dihedraltypes_list();
  }

  if (total_no_oops > 0) {
    oops = (struct OOPList *)calloc(total_no_oops,sizeof(struct OOPList));
    if (oops == NULL) {
      fprintf(stderr,"Trouble allocating memory for oops list - Exiting\n");
      exit(1);
    }
    build_oops_list();

    ooptypes = (struct OOPTypeList *)calloc(MAX_OOP_TYPES,
                                            sizeof(struct OOPTypeList));
    if (ooptypes == NULL) {
      fprintf(stderr,"Trouble allocating memory for ooptypes list - Exiting\n");
      exit(1);
    }

    build_ooptypes_list();
  }

  if ((forcefield & FF_TYPE_CLASS2) && (total_no_angle_angles > 0)) {

    angleangles = (struct AngleAngleList *)calloc(total_no_angle_angles,
                                                  sizeof(struct AngleAngleList));

    if (angleangles == NULL) {
      fprintf(stderr,"Trouble allocating memory for angleangles list - Exiting\n");
      exit(1);
    }
    build_angleangles_list();

    angleangletypes = (struct AngleAngleTypeList *)calloc(MAX_ANGLEANGLE_TYPES,
                                                          sizeof(struct AngleAngleTypeList));
    if (angleangletypes == NULL) {
      fprintf(stderr,"Trouble allocating memory for angleangletypes list - Exiting\n");
      exit(1);
    }
    build_angleangletypes_list();
  }


  if (pflag > 2) {
    int i;
    fprintf(stderr,"Atom Types\n N Potential\n");
    for (i=0; i < no_atom_types; i++) {
      fprintf(stderr," %d %s\n",i,atomtypes[i].potential);
    }

    fprintf(stderr,"Atoms\n");
    for (i=0; i < total_no_atoms; i++) {
      fprintf(stderr,"Atom %3d %2d %-5s  %7.4f  %9.6f %9.6f %9.6f\n",
              i,atoms[i].type,atoms[i].potential,atoms[i].q,
              atoms[i].x[0],atoms[i].x[1],atoms[i].x[2]);
    }

    if (total_no_bonds > 0) {
      fprintf(stderr,"Bond Types\n");
      for (i=0; i < no_bond_types; i++) {
        fprintf(stderr," %d  %d %d  %-s %-s\n",i,bondtypes[i].types[0],
                bondtypes[i].types[1],
                atomtypes[bondtypes[i].types[0]].potential,
                atomtypes[bondtypes[i].types[1]].potential);
      }

      fprintf(stderr,"Bonds\n N  Type  I J\n");
      for (i=0; i < total_no_bonds; i++) {
        fprintf(stderr," %d %d %d %d\n",i,bonds[i].type,
                bonds[i].members[0],
                bonds[i].members[1]);
      }
    }

    if (total_no_angles > 0) {
      fprintf(stderr,"Angle Types\n");
      for (i=0; i < no_angle_types; i++) {
        fprintf(stderr," %d  %d %d %d  %-s %-s %-s\n",i,angletypes[i].types[0],
                angletypes[i].types[1],angletypes[i].types[2],
                atomtypes[angletypes[i].types[0]].potential,
                atomtypes[angletypes[i].types[1]].potential,
                atomtypes[angletypes[i].types[2]].potential);
      }

      fprintf(stderr,"Angles\n N  Type  I  J  K\n");
      for (i=0; i < total_no_angles; i++) {
        fprintf(stderr," %d %d %d %d %d\n",i,angles[i].type,
                angles[i].members[0],
                angles[i].members[1],
                angles[i].members[2]);
      }
    }

    if (total_no_dihedrals > 0) {
      fprintf(stderr,"Dihedral Types\n");
      for (i=0; i < no_dihedral_types; i++) {
        fprintf(stderr," %d  %d %d %d %d %-s %-s %-s %-s\n",i,
                dihedraltypes[i].types[0],
                dihedraltypes[i].types[1],
                dihedraltypes[i].types[2],
                dihedraltypes[i].types[3],
                atomtypes[dihedraltypes[i].types[0]].potential,
                atomtypes[dihedraltypes[i].types[1]].potential,
                atomtypes[dihedraltypes[i].types[2]].potential,
                atomtypes[dihedraltypes[i].types[3]].potential);
      }

      fprintf(stderr,"Dihedrals\n N  Type  I  J  K  L\n");
      for (i=0; i < total_no_dihedrals; i++) {
        fprintf(stderr," %d %d %d %d %d %d\n",i,dihedrals[i].type,
                dihedrals[i].members[0],
                dihedrals[i].members[1],
                dihedrals[i].members[2],
                dihedrals[i].members[3]);
      }
    }

    if (total_no_oops > 0) {
      fprintf(stderr,"Oop Types\n");
      for (i=0; i < no_oop_types; i++) {
        fprintf(stderr," %d  %d %d %d %d %-s %-s %-s %-s\n",i,
                ooptypes[i].types[0],
                ooptypes[i].types[1],
                ooptypes[i].types[2],
                ooptypes[i].types[3],
                atomtypes[ooptypes[i].types[0]].potential,
                atomtypes[ooptypes[i].types[1]].potential,
                atomtypes[ooptypes[i].types[2]].potential,
                atomtypes[ooptypes[i].types[3]].potential);
      }

      fprintf(stderr,"Oops\n N  Type  I  J  K  L\n");
      for (i=0; i < total_no_oops; i++) {
        fprintf(stderr," %d %d %d %d %d %d\n",i,oops[i].type,
                oops[i].members[0],
                oops[i].members[1],
                oops[i].members[2],
                oops[i].members[3]);
      }
    }

    if ((forcefield & FF_TYPE_CLASS2) & (total_no_angle_angles > 0)) {

      fprintf(stderr,"Angleangle Types\n");
      for (i=0; i < no_angleangle_types; i++) {
        fprintf(stderr," %d  %d %d %d %d %-s %-s %-s %-s\n",i,
                angleangletypes[i].types[0],
                angleangletypes[i].types[1],
                angleangletypes[i].types[2],
                angleangletypes[i].types[3],
                atomtypes[angleangletypes[i].types[0]].potential,
                atomtypes[angleangletypes[i].types[1]].potential,
                atomtypes[angleangletypes[i].types[2]].potential,
                atomtypes[angleangletypes[i].types[3]].potential);
      }
      fprintf(stderr,"AngleAngles\n N  Type  I  J  K  L\n");
      for (i=0; i < total_no_angle_angles; i++) {
        fprintf(stderr," %d %d %d %d %d %d\n",i,angleangles[i].type,
                angleangles[i].members[0],
                angleangles[i].members[1],
                angleangles[i].members[2],
                angleangles[i].members[3]);
      }
    }
  }

  if (pflag > 1) {
    fprintf(stderr,"\n");
    fprintf(stderr," Number of bonds, types = %7d %3d\n",
            total_no_bonds,no_bond_types);
    fprintf(stderr," Number of angles, types = %7d %3d\n",
            total_no_angles, no_angle_types);
    fprintf(stderr," Number of dihedrals, types = %7d %3d\n",
            total_no_dihedrals, no_dihedral_types);
    fprintf(stderr," Number of out-of-planes, types = %7d %3d\n",
            total_no_oops, no_oop_types);
    if (forcefield & FF_TYPE_CLASS2)
      fprintf(stderr," Number of Angle Angle Terms, types =   %7d %3d\n",
              total_no_angle_angles, no_angleangle_types);
  }
}

int count_bonds()
{
  int i,j,n;

  for (n=0,i=0; i < total_no_atoms; i++) {
    for (j=0; j < atoms[i].no_connect; j++) {
      if (i < atoms[i].conn_no[j]) n++;
    }
  }
  return n;
}

void build_bonds_list()
{
  int i,j,n;

  for (n=0,i=0; i < total_no_atoms; i++) {
    for (j=0; j < atoms[i].no_connect; j++) {
      if (i < atoms[i].conn_no[j]) {
        bonds[n  ].type = 0;
        bonds[n  ].members[0] = i;
        bonds[n++].members[1] = atoms[i].conn_no[j];
      }
    }
  }
  return;
}

int count_angles()
{
  int i,j,k,n;

  for (n=0,j=0; j < total_no_atoms; j++) {
    if (atoms[j].no_connect > 1) {
      for (i=0; i < atoms[j].no_connect-1; i++) {
        for (k=i+1; k < atoms[j].no_connect; k++) {
          n++;
        }
      }
    }
  }
  return n;
}

void build_angles_list()
{
  int i,j,k,n;

  for (n=0,j=0; j < total_no_atoms; j++) {
    if (atoms[j].no_connect > 1) {
      for (i=0; i < atoms[j].no_connect-1; i++) {
        for (k=i+1; k < atoms[j].no_connect; k++) {
          angles[n  ].type = 0;
          angles[n  ].members[0] = atoms[j].conn_no[i];
          angles[n  ].members[1] = j;
          angles[n++].members[2] = atoms[j].conn_no[k];
        }
      }
    }
  }
  return;
}

int count_dihedrals()
{
  int i,j,k,l,n;
  int ii,kk,ll;

  for (n=0,j=0; j < total_no_atoms; j++) {
    if (atoms[j].no_connect > 1) {
      for (kk=0; kk < atoms[j].no_connect; kk++) {
        k = atoms[j].conn_no[kk];
        if (atoms[k].no_connect > 1) {
          if (j < k) {
            for (ii=0; ii < atoms[j].no_connect; ii++) {
              i = atoms[j].conn_no[ii];
              if (i != k) {
                for (ll=0; ll < atoms[k].no_connect; ll++) {
                  l = atoms[k].conn_no[ll];
                  if ((l != j) && (i != l)) n++;
                }
              }
            }
          }
        }
      }
    }
  }
  return n;
}

void build_dihedrals_list()
{
  int i,j,k,l,n;
  int ii,kk,ll;

  for (n=0,j=0; j < total_no_atoms; j++) {
    if (atoms[j].no_connect > 1) {
      for (kk=0; kk < atoms[j].no_connect; kk++) {
        k = atoms[j].conn_no[kk];
        if (atoms[k].no_connect > 1) {
          if (j < k) {
            for (ii=0; ii < atoms[j].no_connect; ii++) {
              i = atoms[j].conn_no[ii];
              if (i != k) {
                for (ll=0; ll < atoms[k].no_connect; ll++) {
                  l = atoms[k].conn_no[ll];
                  if ((l != j) && (i != l)) {
                    dihedrals[n  ].type = 0;
                    dihedrals[n  ].members[0] = i;
                    dihedrals[n  ].members[1] = j;
                    dihedrals[n  ].members[2] = k;
                    dihedrals[n++].members[3] = l;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}

int
count_oops()
{
  int j,n;

  for (n=0,j=0; j < total_no_atoms; j++) {
    if (atoms[j].no_connect == 3) n++;
  }
  return n;
}

void build_oops_list()
{
  int j,n;

  for (n=0,j=0; j < total_no_atoms; j++) {

    if (atoms[j].no_connect == 3) {
      oops[n  ].type = 0;
      oops[n  ].members[0] = atoms[j].conn_no[0];
      oops[n  ].members[1] = j;
      oops[n  ].members[2] = atoms[j].conn_no[1];
      oops[n++].members[3] = atoms[j].conn_no[2];
    }
  }
  return;
}

int count_angle_angles()
{
  int num_triples[11] = {0,0,0,0,4,10,20,35,56,84,120};
  int j,n;

  for (n=0,j=0; j < total_no_atoms; j++) {
    n += num_triples[atoms[j].no_connect];
  }
  return n;
}

void build_angleangles_list()
{
  int i,j,k,l,nc,n;

  for (n=0,j=0; j < total_no_atoms; j++) {
    nc = atoms[j].no_connect;
    if (nc > 3) {
      for (i=0; i < nc-2; i++) {
        for (k=i+1; k < nc-1; k++) {
          for (l=k+1; l < nc; l++) {
            angleangles[n].type = 0;
            angleangles[n  ].members[0] = atoms[j].conn_no[i];
            angleangles[n  ].members[1] = j;
            angleangles[n  ].members[2] = atoms[j].conn_no[k];
            angleangles[n++].members[3] = atoms[j].conn_no[l];
          }
        }
      }
    }
  }
  return;
}


void build_atomtypes_list()
{
  int j,k,n,match,atom_type=0;

  strncpy(atomtypes[0].potential,atoms[0].potential,5);
  atoms[0].type = 0;

  atomtypes[0].no_connect = atoms[0].no_connect;

  for (n=1,j=1; j < total_no_atoms; j++) {
    match = 0;
    k = 0;
    while (!match && (k < n)) {
      if (strncmp(atomtypes[k].potential,atoms[j].potential,5) == 0) {
        match = 1;
        atom_type = k;
        if (atomtypes[k].no_connect != atoms[j].no_connect) {
          if (pflag > 0) fprintf(stderr," WARNING inconsistent # of connects on atom %d type %s\n",j,
                                 atomtypes[k].potential);
        }
      } else k++;
    }
    if (match == 0) {
      atom_type = n;
      atomtypes[n].no_connect = atoms[j].no_connect;
      strncpy(atomtypes[n++].potential,atoms[j].potential,5);
    }
    if (n >= MAX_ATOM_TYPES) {
      fprintf(stderr,"Too many atom types (> 100) - error\n");
      exit(1);
    }
    atoms[j].type = atom_type;
  }
  no_atom_types = n;
  return;
}

void build_bondtypes_list() {
  int j,k,n,match,bond_type=0;
  int typei,typej;

  for (n=0,j=0; j < total_no_bonds; j++) {
    typei = atoms[bonds[j].members[0]].type;
    typej = atoms[bonds[j].members[1]].type;
    if (typej < typei) {
      swap_ints(&typei,&typej);
      swap_ints(&bonds[j].members[0],&bonds[j].members[1]);
    }

    match = 0;
    k = 0;
    while (!match && (k < n)) {
      if ((typei == bondtypes[k].types[0]) &&
          (typej == bondtypes[k].types[1])) {
        match = 1;
        bond_type = k;
      } else k++;
    }
    if (match == 0) {
      bond_type = n;
      bondtypes[n  ].types[0] = typei;
      bondtypes[n++].types[1] = typej;
    }
    if (n >= MAX_BOND_TYPES) {
      fprintf(stderr,"Too many bond types (> 200) - error\n");
      exit(1);
    }

    bonds[j].type = bond_type;
  }
  no_bond_types = n;
  return;
}

void build_angletypes_list()
{
  int j,k,n,match,angle_type=0;
  int typei,typej,typek;

  for (n=0,j=0; j < total_no_angles; j++) {
    typei = atoms[angles[j].members[0]].type;
    typej = atoms[angles[j].members[1]].type;
    typek = atoms[angles[j].members[2]].type;
    if (typek < typei) {
      swap_ints(&typei,&typek);
      swap_ints(&angles[j].members[0],&angles[j].members[2]);
    }

    match = 0;
    k = 0;
    while (!match && (k < n)) {
      if ((typei == angletypes[k].types[0]) &&
          (typej == angletypes[k].types[1]) &&
          (typek == angletypes[k].types[2])) {
        match = 1;
        angle_type = k;
      } else k++;
    }
    if (match == 0) {
      angle_type = n;
      angletypes[n  ].types[0] = typei;
      angletypes[n  ].types[1] = typej;
      angletypes[n++].types[2] = typek;
    }
    if (n >= MAX_ANGLE_TYPES) {
      fprintf(stderr,"Too many angle types (> 300) - error\n");
      exit(1);
    }
    angles[j].type = angle_type;
  }
  no_angle_types = n;
  return;
}

void build_dihedraltypes_list()
{
  int j,k,n,match,dihedral_type=0;
  int typei,typej,typek,typel;

  for (n=0,j=0; j < total_no_dihedrals; j++) {
    typei = atoms[dihedrals[j].members[0]].type;
    typej = atoms[dihedrals[j].members[1]].type;
    typek = atoms[dihedrals[j].members[2]].type;
    typel = atoms[dihedrals[j].members[3]].type;
    if ((typek < typej) || ((typej == typek) && (typel < typei))) {
      swap_ints(&typej,&typek);
      swap_ints(&dihedrals[j].members[1],&dihedrals[j].members[2]);
      swap_ints(&typei,&typel);
      swap_ints(&dihedrals[j].members[0],&dihedrals[j].members[3]);
    }

    match = 0;
    k = 0;
    while (!match && (k < n)) {
      if ((typei == dihedraltypes[k].types[0]) &&
          (typej == dihedraltypes[k].types[1]) &&
          (typek == dihedraltypes[k].types[2]) &&
          (typel == dihedraltypes[k].types[3])) {
        match = 1;
        dihedral_type = k;
      } else k++;
    }
    if (match == 0) {
      dihedral_type = n;
      dihedraltypes[n  ].types[0] = typei;
      dihedraltypes[n  ].types[1] = typej;
      dihedraltypes[n  ].types[2] = typek;
      dihedraltypes[n++].types[3] = typel;
    }
    if (n >= MAX_DIHEDRAL_TYPES) {
      fprintf(stderr,"Too many dihedral types (> 400) - error\n");
      exit(1);
    }
    dihedrals[j].type = dihedral_type;
  }
  no_dihedral_types = n;
  return;
}

void build_ooptypes_list()
{
  int j,k,n,match,oop_type=0;
  int temp_types[3],temp_pos[3];
  int typei,typej,typek,typel;

  for (n=0,j=0; j < total_no_oops; j++) {
    typei = atoms[oops[j].members[0]].type;
    typej = atoms[oops[j].members[1]].type;
    typek = atoms[oops[j].members[2]].type;
    typel = atoms[oops[j].members[3]].type;
    temp_types[0] = typei;
    temp_types[1] = typek;
    temp_types[2] = typel;

    bubble_sort(3,temp_types,temp_pos);

    typei = temp_types[0];
    typek = temp_types[1];
    typel = temp_types[2];
    temp_types[0] = oops[j].members[0];
    temp_types[1] = oops[j].members[2];
    temp_types[2] = oops[j].members[3];
    oops[j].members[0] = temp_types[temp_pos[0]];
    oops[j].members[2] = temp_types[temp_pos[1]];
    oops[j].members[3] = temp_types[temp_pos[2]];

    match = 0;
    k = 0;
    while (!match && (k < n)) {
      if ((typei == ooptypes[k].types[0]) &&
          (typej == ooptypes[k].types[1]) &&
          (typek == ooptypes[k].types[2]) &&
          (typel == ooptypes[k].types[3])) {
        match = 1;
        oop_type = k;
      } else k++;
    }
    if (match == 0) {

      oop_type = n;
      ooptypes[n  ].types[0] = typei;
      ooptypes[n  ].types[1] = typej;
      ooptypes[n  ].types[2] = typek;
      ooptypes[n++].types[3] = typel;
    }
    if (n >= MAX_OOP_TYPES) {
      fprintf(stderr,"Too many oop types (> 400) - error\n");
      exit(1);
    }
    oops[j].type = oop_type;
  }
  no_oop_types = n;
  return;
}

void build_angleangletypes_list()
{
  int j,k,n,match,angleangle_type=0;
  int temp_types[3],temp_pos[3];
  int typei,typej,typek,typel;

  for (n=0,j=0; j < total_no_angle_angles; j++) {

    typei = atoms[angleangles[j].members[0]].type;
    typej = atoms[angleangles[j].members[1]].type;
    typek = atoms[angleangles[j].members[2]].type;
    typel = atoms[angleangles[j].members[3]].type;

    temp_types[0] = typei;
    temp_types[1] = typek;
    temp_types[2] = typel;

    bubble_sort(3,temp_types,temp_pos);

    typei = temp_types[0];
    typek = temp_types[1];
    typel = temp_types[2];

    temp_types[0] = angleangles[j].members[0];
    temp_types[1] = angleangles[j].members[2];
    temp_types[2] = angleangles[j].members[3];

    angleangles[j].members[0] = temp_types[temp_pos[0]];
    angleangles[j].members[2] = temp_types[temp_pos[1]];
    angleangles[j].members[3] = temp_types[temp_pos[2]];

    match = 0;
    k = 0;
    while (!match && (k < n)) {
      if ((typei == angleangletypes[k].types[0]) &&
          (typej == angleangletypes[k].types[1]) &&
          (typek == angleangletypes[k].types[2]) &&
          (typel == angleangletypes[k].types[3])) {
        match = 1;
        angleangle_type = k;
      } else k++;
    }
    if (match == 0) {
      angleangle_type = n;
      angleangletypes[n  ].types[0] = typei;
      angleangletypes[n  ].types[1] = typej;
      angleangletypes[n  ].types[2] = typek;
      angleangletypes[n++].types[3] = typel;
    }

    if (n >= MAX_ANGLEANGLE_TYPES) {
      fprintf(stderr,"Too many angleangle types (> 400) - error\n");
      exit(1);
    }
    angleangles[j].type = angleangle_type;
  }
  no_angleangle_types = n;
  return;
}

void swap_ints(int *i, int *j)
{
  int temp;

  temp = *i;
  *i = *j;
  *j = temp;

  return;
}

void bubble_sort(int n, int *val, int *pos)
{
  int i,j;

  for (i=0; i < n; i++) pos[i] = i;
  for (i=0; i < n-1; i++) {
    for (j=1; j < n; j++) {
      if (val[j] < val[i]) {
        swap_ints(&val[i],&val[j]);
        swap_ints(&pos[i],&pos[j]);
      }
    }
  }
}
