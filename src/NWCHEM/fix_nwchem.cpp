    qm2lmp_force = HARTREE_TO_KCAL_MOLE * ANGSTROM_TO_BOHR;
  }

  // initializations

  nqm = 0;
  qmIDs = nullptr;
  xqm = fqm = nullptr;
  qpotential = qqm = nullptr;
  qm2lmp = nullptr;

  qmenergy = 0.0;
}

/* ---------------------------------------------------------------------- */

FixNWChem::~FixNWChem()
{
  delete [] id_pe;
  memory->destroy(qmIDs);
  memory->destroy(xqm);
  memory->destroy(qpotential);
  memory->destroy(fqm);
  memory->destroy(qqm);
  memory->destroy(qm2lmp);
}

/* ---------------------------------------------------------------------- */

int FixNWChem::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNWChem::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix nwchem requires 3d problem");

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix nwchem requires fully periodic or "
                  "fully non-periodic system");
  
  if (atom->q_flag == 0 || force->pair == nullptr || force->kspace == nullptr)
    error->all(FLERR,"Fix nwchem cannot compute Coulomb potential");

  // c_pe = instance of compute pe/atom

  int ipe = modify->find_compute(id_pe);
  if (ipe < 0) error->all(FLERR,"Could not find fix nwchem compute ID");
  c_pe = modify->compute[ipe];

  // pair_coul = hybrid pair style that computes short-range Coulombics
  // NOTE: could be another coul like coul/msm ?
  
  if (!force->pair) error->all(FLERR,"Fix nwchem requires a pair style");

  pair_coul = force->pair_match("hybrid/overlay",1,0);
  if (!pair_coul) 
    error->all(FLERR,"Fix nwchem requires pair style hybrid/overlay");

  pair_coul = force->pair_match("coul/long",1,0);
  if (!pair_coul) 
    error->all(FLERR,"Fix nwchem requires pair sub-style coul/long");

  // fix group = QM atoms
  // one-time initialization of qmIDs
  // NOTE: need to sort in ascending order to match NWChem ?
  // NOTE: make nqm an int, check that group count is not 0 or > MAXBIGINT
  
  if (nqm == 0) {
    nqm = group->count(igroup);
    memory->create(qmIDs,nqm,"nwchem:qmIDs");
    memory->create(xqm,nqm,3,"nwchem:xqm");
    memory->create(qpotential,nqm,"nwchem:qpotential");
    memory->create(fqm,nqm,3,"nwchem:fqm");
    memory->create(qqm,nqm,"nwchem:qqm");
    memory->create(qm2lmp,nqm,"nwchem:qm2lmp");

    tagint *tag = atom->tag;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    nqm = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) qmIDs[nqm++] = tag[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixNWChem::setup(int vflag)
{
  //newsystem = 1;
  post_force(vflag);
  //newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixNWChem::min_setup(int vflag)
{
  //newsystem = 1;
  post_force(vflag);
  //newsystem = 0;
}

/* ---------------------------------------------------------------------- */

void FixNWChem::pre_force(int vflag)
{
  int ilocal,jlocal;
  double delx,dely,delz,rsq;

  // invoke pair hybrid sub-style pair coul/long and Kspace directly
  // set eflag = 2 so they calculate per-atom energy
  // NOTE: need to comm ghost per-atom energy (serial or parallel)

  pair_coul->compute(2,0);
  double *eatom_pair = pair_coul->eatom;

  double *eatom_kspace = nullptr;
  if (force->kspace) {
    kspace->compute(2,0);
    eatom_kspace = force->kspace->eatom;
  }

  // on reneigh step, reset qm2lmp vector = indices of QM atoms

  if (neighbor->ago == 0)
    for (int i = 0; i < nqm; i++)
      qm2lmp[i] = atom->map(qmIDs[i]);

  // create 2 NWChem inputs: xqm and qpotential (Coulomb potential)
  // qpotential[i] = (eatom[i] from pair_coul + kspace) / Qi
  // subtract off Qj/Rij for QM I interacting with all other QM J atoms
 
  double **x = atom->x;
  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];

    // NOTE: what if LAMMPS atom moves via PBC to other side of LAMMPS box ?

    xqm[i][0] = x[ilocal][0];
    xqm[i][1] = x[ilocal][1];
    xqm[i][2] = x[ilocal][2];

    if (q[i] == 0.0) qpotential[i] = 0.0;
    else if (!force->kspace)
      qpotential[i] = eatom_pair[ilocal] / q[ilocal];
    else
      qpotential[i] = (eatom_pair[ilocal] + eatom_kspace[ilocal]) / q[ilocal];

    for (int j = 0; j < nqm; j++) {
      if (j == i) continue;
      jlocal = qm2lmp[j];
      // NOTE: apply PBC
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      qpotential[i] -= qqrd2e * q[j] / sqrt(rsq);
    }
  }

  // unit conversion from LAMMPS to NWChem

  for (int i = 0; i < nqm; i++) {
    xqm[i][0] *= lmp2qm_distance;
    xqm[i][1] *= lmp2qm_distance;
    xqm[i][2] *= lmp2qm_distance;
    qpotential[i] *= lmp2qm_energy;
  }

  // call to NWChem
  // input/outputs are only for QM atoms
  // inputs = xqm and qpotential
  // outputs = fqm and qqm

  int nwerr = pwdft::pspw_minimizer(world,&xqm[0][0],qpotential,&fqm[0][0],qqm);
  if (nwerr) error->all(FLERR,"Internal NWChem error");

  //latte(flags,&natoms,coords,type,&ntypes,mass,boxlo,boxhi,&domain->xy,
  //      &domain->xz,&domain->yz,forces,&maxiter,&latte_energy,
  //      &atom->v[0][0],&update->dt,virial,&newsystem,&latteerror);

  // unit conversion from NWChem to LAMMPS

  for (int i = 0; i < nqm; i++) {
    fqm[i][0] *= qm2lmp_force;
    fqm[i][1] *= qm2lmp_force;
    fqm[i][2] *= qm2lmp_force;
  }

  // reset Q of QM atoms

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];
    q[ilocal] = qqm[i];
  }

  // reset LAMMPS forces to zero
  // NOTE: do this via Verlet or rRESPA force_clear() ?
  // NOTE: should rRESPA be not allowed for this fix ?
  // NOTE: what about min
  
  // update->integrate->force_clear();  // what is external_force_clear = OPENMP ?
}

/* ---------------------------------------------------------------------- */

void FixNWChem::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::post_force(int vflag)
{
  // add NWChem QM forces to QM atoms

  double **f = atom->f;

  for (int i = 0; i < nqm; i++) {
    ilocal = qm2lmp[i];
    f[ilocal][0] += fqm[i][0];
    f[ilocal][1] += fqm[i][1];
    f[ilocal][2] += fqm[i][2];
  }

  // NOTE: what is qmenergy for contrib of this fix to system eng

  // trigger per-atom energy computation on next step by pair/kspace

  c_pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixNWChem::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   QM energy from NWChem
------------------------------------------------------------------------- */

double FixNWChem::compute_scalar()
{
  return qmenergy;
}
