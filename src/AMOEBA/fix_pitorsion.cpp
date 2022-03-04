// clang-format off
/* ----------------------------------------------------------------------
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
   Contributing authors:
   Xiaohu Hu, CMB/ORNL (hux2@ornl.gov)
   David Hyde-Volpe, Tigran Abramyan, and Robert A. Latour (Clemson University)
   Chris Lorenz (Kings College-London)

   Implementation of the CHARMM CMAP; adds an extra energy term for the
   peptide backbone dihedrals.  The tools/ch2lmp/charmm2lammps.pl
   conversion script, which generates an extra section in the LAMMPS data
   file, is needed in order to generate the info used by this fix style.

   References:
   - MacKerell et al., J. Am. Chem. Soc. 126(2004):698-699.
   - MacKerell et al., J. Comput. Chem. 25(2004):1400-1415.
------------------------------------------------------------------------- */

// read_data data.ubiquitin fix pitorsion pitorsions PiTorsions


#include "fix_pitorsion.h"

#include <cmath>

#include <cstring>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "force.h"
#include "comm.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define MAXLINE 256
#define LISTDELTA 10000
#define LB_FACTOR 1.5

/* ---------------------------------------------------------------------- */

FixPiTorsion::FixPiTorsion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  crosstermlist(nullptr), num_crossterm(nullptr), crossterm_type(nullptr),
  crossterm_atom1(nullptr), crossterm_atom2(nullptr), crossterm_atom3(nullptr),
  crossterm_atom4(nullptr), crossterm_atom5(nullptr),
  g_axis(nullptr), cmapgrid(nullptr), d1cmapgrid(nullptr), d2cmapgrid(nullptr),
  d12cmapgrid(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal fix cmap command");

  restart_global = 1;
  restart_peratom = 1;
  energy_global_flag = energy_peratom_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_energy = thermo_virial = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  peratom_freq = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  wd_header = 1;
  wd_section = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // allocate memory for CMAP data

  memory->create(g_axis,CMAPDIM,"cmap:g_axis");
  memory->create(cmapgrid,6,CMAPDIM,CMAPDIM,"cmap:grid");
  memory->create(d1cmapgrid,6,CMAPDIM,CMAPDIM,"cmap:d1grid");
  memory->create(d2cmapgrid,6,CMAPDIM,CMAPDIM,"cmap:d2grid");
  memory->create(d12cmapgrid,6,CMAPDIM,CMAPDIM,"cmap:d12grid");

  // read and setup CMAP data

  read_grid_map(arg[3]);

  // perform initial allocation of atom-based arrays
  // register with Atom class

  num_crossterm = nullptr;
  crossterm_type = nullptr;
  crossterm_atom1 = nullptr;
  crossterm_atom2 = nullptr;
  crossterm_atom3 = nullptr;
  crossterm_atom4 = nullptr;
  crossterm_atom5 = nullptr;

  nmax_previous = 0;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // local list of crossterms

  ncmap = 0;
  maxcrossterm = 0;
  crosstermlist = nullptr;
}

/* --------------------------------------------------------------------- */

FixCMAP::~FixCMAP()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);

  memory->destroy(g_axis);
  memory->destroy(cmapgrid);
  memory->destroy(d1cmapgrid);
  memory->destroy(d2cmapgrid);
  memory->destroy(d12cmapgrid);

  memory->destroy(crosstermlist);

  memory->destroy(num_crossterm);
  memory->destroy(crossterm_type);
  memory->destroy(crossterm_atom1);
  memory->destroy(crossterm_atom2);
  memory->destroy(crossterm_atom3);
  memory->destroy(crossterm_atom4);
  memory->destroy(crossterm_atom5);
}

/* ---------------------------------------------------------------------- */

int FixCMAP::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCMAP::init()
{
  int i;
  double angle;

  i = 0;
  angle = -180.0;
  while (angle < 180.0) {
    g_axis[i] = angle;
    angle += CMAPDX;
    i++;
  }

  // pre-compute the derivatives of the maps

  for (i = 0; i < 6; i++)
    set_map_derivatives(cmapgrid[i],d1cmapgrid[i],d2cmapgrid[i],d12cmapgrid[i]);

  // define newton_bond here in case restart file was read (not data file)

  newton_bond = force->newton_bond;

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* --------------------------------------------------------------------- */

void FixCMAP::setup(int vflag)
{
  pre_neighbor();

  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* --------------------------------------------------------------------- */

void FixCMAP::setup_pre_neighbor()
{
  pre_neighbor();
}

/* --------------------------------------------------------------------- */

void FixCMAP::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* --------------------------------------------------------------------- */

void FixCMAP::min_setup(int vflag)
{
  pre_neighbor();
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   store local neighbor list as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

void FixCMAP::pre_neighbor()
{
  int i,m,atom1,atom2,atom3,atom4,atom5;

  // guesstimate initial length of local crossterm list
  // if ncmap was not set (due to read_restart, no read_data),
  //   then list will grow by LISTDELTA chunks

  if (maxcrossterm == 0) {
    if (nprocs == 1) maxcrossterm = ncmap;
    else maxcrossterm = static_cast<int> (LB_FACTOR*ncmap/nprocs);
    memory->create(crosstermlist,maxcrossterm,6,"cmap:crosstermlist");
  }

  int nlocal = atom->nlocal;

  ncrosstermlist = 0;

  for (i = 0; i < nlocal; i++) {
    for (m = 0; m < num_crossterm[i]; m++) {
      atom1 = atom->map(crossterm_atom1[i][m]);
      atom2 = atom->map(crossterm_atom2[i][m]);
      atom3 = atom->map(crossterm_atom3[i][m]);
      atom4 = atom->map(crossterm_atom4[i][m]);
      atom5 = atom->map(crossterm_atom5[i][m]);

      if (atom1 == -1 || atom2 == -1 || atom3 == -1 ||
          atom4 == -1 || atom5 == -1)
        error->one(FLERR,"CMAP atoms {} {} {} {} {} missing on "
                                     "proc {} at step {}",
                                     crossterm_atom1[i][m],crossterm_atom2[i][m],
                                     crossterm_atom3[i][m],crossterm_atom4[i][m],
                                     crossterm_atom5[i][m],me,update->ntimestep);
      atom1 = domain->closest_image(i,atom1);
      atom2 = domain->closest_image(i,atom2);
      atom3 = domain->closest_image(i,atom3);
      atom4 = domain->closest_image(i,atom4);
      atom5 = domain->closest_image(i,atom5);

      if (i <= atom1 && i <= atom2 && i <= atom3 &&
          i <= atom4 && i <= atom5) {
        if (ncrosstermlist == maxcrossterm) {
          maxcrossterm += LISTDELTA;
          memory->grow(crosstermlist,maxcrossterm,6,"cmap:crosstermlist");
        }
        crosstermlist[ncrosstermlist][0] = atom1;
        crosstermlist[ncrosstermlist][1] = atom2;
        crosstermlist[ncrosstermlist][2] = atom3;
        crosstermlist[ncrosstermlist][3] = atom4;
        crosstermlist[ncrosstermlist][4] = atom5;
        crosstermlist[ncrosstermlist][5] = crossterm_type[i][m];
        ncrosstermlist++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixCMAP::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
}

/* ----------------------------------------------------------------------
   compute PiTorsion terms as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

void FixCMAP::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;

  for (int n = 0; n < npitors; n++) {
    ia = ipit[n][0];
    ib = ipit[n][1];
    ic = ipit[n][2];
    id = ipit[n][3];
    ie = ipit[n][4];
    ig = ipit[n][5];

    // compute the value of the pi-system torsion angle

    xia = x[ia][0];
    yia = x[ia][1];
    zia = x[ia][2];
    xib = x[ib][0];
    yib = x[ib][1];
    zib = x[ib][2];
    xic = x[ic][0];
    yic = x[ic][1];
    zic = x[ic][2];
    xid = x[id][0];
    yid = x[id][1];
    zid = x[id][2];
    xie = x[ie][0];
    yie = x[ie][1];
    zie = x[ie][2];
    xig = x[ig][0];
    yig = x[ig][1];
    zig = x[ig][2];

    xad = xia - xid;
    yad = yia - yid;
    zad = zia - zid;
    xbd = xib - xid;
    ybd = yib - yid;
    zbd = zib - zid;
    xec = xie - xic;
    yec = yie - yic;
    zec = zie - zic;
    xgc = xig - xic;
    ygc = yig - yic;
    zgc = zig - zic;

    xip = yad*zbd - ybd*zad + xic;
    yip = zad*xbd - zbd*xad + yic;
    zip = xad*ybd - xbd*yad + zic;
    xiq = yec*zgc - ygc*zec + xid;
    yiq = zec*xgc - zgc*xec + yid;
    ziq = xec*ygc - xgc*yec + zid;
    xcp = xic - xip;
    ycp = yic - yip;
    zcp = zic - zip;
    xdc = xid - xic;
    ydc = yid - yic;
    zdc = zid - zic;
    xqd = xiq - xid;
    yqd = yiq - yid;
    zqd = ziq - zid;

    xt = ycp*zdc - ydc*zcp;
    yt = zcp*xdc - zdc*xcp;
    zt = xcp*ydc - xdc*ycp;
    xu = ydc*zqd - yqd*zdc;
    yu = zdc*xqd - zqd*xdc;
    zu = xdc*yqd - xqd*ydc;
    xtu = yt*zu - yu*zt;
    ytu = zt*xu - zu*xt;
    ztu = xt*yu - xu*yt;
    rt2 = xt*xt + yt*yt + zt*zt;
    ru2 = xu*xu + yu*yu + zu*zu;
    rtru = sqrt(rt2*ru2);

    if (rtru <= 0.0) continue;

    rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc);
    cosine = (xt*xu + yt*yu + zt*zu) / rtru;
    sine = (xdc*xtu + ydc*ytu + zdc*ztu) / (rdc*rtru);

    // set the pi-system torsion parameters for this angle
        
    v2 = kpit(i);
    c2 = -1.0;
    s2 = 0.0;

    // compute the multiple angle trigonometry and the phase terms

    cosine2 = cosine*cosine - sine*sine;
    sine2 = 2.0 * cosine * sine;
    phi2 = 1.0 + (cosine2*c2 + sine2*s2);
    dphi2 = 2.0 * (cosine2*s2 - sine2*c2);

    // calculate pi-system torsion energy and master chain rule term

    e = ptorunit * v2 * phi2;
    dedphi = ptorunit * v2 * dphi2;

    // chain rule terms for first derivative components

    xdp = xid - xip;
    ydp = yid - yip;
    zdp = zid - zip;
    xqc = xiq - xic;
    yqc = yiq - yic;
    zqc = ziq - zic;
    dedxt = dedphi * (yt*zdc - ydc*zt) / (rt2*rdc);
    dedyt = dedphi * (zt*xdc - zdc*xt) / (rt2*rdc);
    dedzt = dedphi * (xt*ydc - xdc*yt) / (rt2*rdc);
    dedxu = -dedphi * (yu*zdc - ydc*zu) / (ru2*rdc);
    dedyu = -dedphi * (zu*xdc - zdc*xu) / (ru2*rdc);
    dedzu = -dedphi * (xu*ydc - xdc*yu) / (ru2*rdc);

    // compute first derivative components for pi-system angle

    dedxip = zdc*dedyt - ydc*dedzt;
    dedyip = xdc*dedzt - zdc*dedxt;
    dedzip = ydc*dedxt - xdc*dedyt;
    dedxic = ydp*dedzt - zdp*dedyt + zqd*dedyu - yqd*dedzu;
    dedyic = zdp*dedxt - xdp*dedzt + xqd*dedzu - zqd*dedxu;
    dedzic = xdp*dedyt - ydp*dedxt + yqd*dedxu - xqd*dedyu;
    dedxid = zcp*dedyt - ycp*dedzt + yqc*dedzu - zqc*dedyu;
    dedyid = xcp*dedzt - zcp*dedxt + zqc*dedxu - xqc*dedzu;
    dedzid = ycp*dedxt - xcp*dedyt + xqc*dedyu - yqc*dedxu;
    dedxiq = zdc*dedyu - ydc*dedzu;
    dedyiq = xdc*dedzu - zdc*dedxu;
    dedziq = ydc*dedxu - xdc*dedyu;

    // compute first derivative components for individual atoms

    dedxia = ybd*dedzip - zbd*dedyip;
    dedyia = zbd*dedxip - xbd*dedzip;
    dedzia = xbd*dedyip - ybd*dedxip;
    dedxib = zad*dedyip - yad*dedzip;
    dedyib = xad*dedzip - zad*dedxip;
    dedzib = yad*dedxip - xad*dedyip;
    dedxie = ygc*dedziq - zgc*dedyiq;
    dedyie = zgc*dedxiq - xgc*dedziq;
    dedzie = xgc*dedyiq - ygc*dedxiq;
    dedxig = zec*dedyiq - yec*dedziq;
    dedyig = xec*dedziq - zec*dedxiq;
    dedzig = yec*dedxiq - xec*dedyiq;
    dedxic = dedxic + dedxip - dedxie - dedxig;
    dedyic = dedyic + dedyip - dedyie - dedyig;
    dedzic = dedzic + dedzip - dedzie - dedzig;
    dedxid = dedxid + dedxiq - dedxia - dedxib;
    dedyid = dedyid + dedyiq - dedyia - dedyib;
    dedzid = dedzid + dedziq - dedzia - dedzib;

    // increment the total pi-system torsion energy and gradient

    ept += e;

    f[ia][0] += dedxia;
    f[ia][1] += dedyia;
    f[ia][2] += dedzia;
    f[ib][0] += dedxib;
    f[ib][1] += dedyib;
    f[ib][2] += dedzib;
    f[ic][0] += dedxic;
    f[ic][1] += dedyic;
    f[ic][2] += dedzic;
    f[id][0] += dedxid;
    f[id][1] += dedyid;
    f[id][2] += dedzid;
    f[ie][0] += dedxie;
    f[ie][1] += dedyie;
    f[ie][2] += dedzie;
    f[ig][0] += dedxig;
    f[ig][1] += dedyig;
    f[ig][2] += dedzig;

    // increment the internal virial tensor components

    vxterm = dedxid + dedxia + dedxib;
    vyterm = dedyid + dedyia + dedyib;
    vzterm = dedzid + dedzia + dedzib;
    vxx = xdc*vxterm + xcp*dedxip - xqd*dedxiq;
    vyx = ydc*vxterm + ycp*dedxip - yqd*dedxiq;
    vzx = zdc*vxterm + zcp*dedxip - zqd*dedxiq;
    vyy = ydc*vyterm + ycp*dedyip - yqd*dedyiq;
    vzy = zdc*vyterm + zcp*dedyip - zqd*dedyiq;
    vzz = zdc*vzterm + zcp*dedzip - zqd*dedziq;

    virial[0] += vxx;
    virial[1] += vyy;
    virial[2] += vzz;
    virial[3] += vyx;
    virial[4] += vzx;
    virial[5] += vzy;
  }
}

/* ---------------------------------------------------------------------- */

void FixCMAP::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCMAP::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of CMAP term
------------------------------------------------------------------------- */

double FixCMAP::compute_scalar()
{
  double all;
  MPI_Allreduce(&ecmap,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods to read and write data file
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void FixCMAP::read_data_header(char *line)
{
  if (strstr(line,"pitorsions")) {
    sscanf(line,BIGINT_FORMAT,&ncmap);
  } else error->all(FLERR,"Invalid read data header line for fix pitorsion");

  // didn't set in constructor because this fix could be defined
  // before newton command

  newton_bond = force->newton_bond;
}

/* ----------------------------------------------------------------------
   unpack N lines in buf from section of data file labeled by keyword
   id_offset is applied to atomID fields if multiple data files are read
   store CMAP interactions as if newton_bond = OFF, even if actually ON
------------------------------------------------------------------------- */

void FixCMAP::read_data_section(char *keyword, int n, char *buf,
                                 tagint id_offset)
{
  int m,tmp,itype;
  tagint atom1,atom2,atom3,atom4,atom5;
  char *next;

  next = strchr(buf,'\n');
  *next = '\0';
  int nwords = utils::count_words(utils::trim_comment(buf));
  *next = '\n';

  if (nwords != 7)
    error->all(FLERR,"Incorrect {} format in data file",keyword);

  // loop over lines of PiTorsions
  // tokenize the line into values
  // add crossterm to one of my atoms, depending on newton_bond

  for (int i = 0; i < n; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    sscanf(buf,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT
           " " TAGINT_FORMAT " " TAGINT_FORMAT,
           &tmp,&itype,&atom1,&atom2,&atom3,&atom4,&atom5);

    atom1 += id_offset;
    atom2 += id_offset;
    atom3 += id_offset;
    atom4 += id_offset;
    atom5 += id_offset;

    if ((m = atom->map(atom1)) >= 0) {
      if (num_crossterm[m] == CMAPMAX)
        error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    if ((m = atom->map(atom2)) >= 0) {
      if (num_crossterm[m] == CMAPMAX)
        error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    if ((m = atom->map(atom3)) >= 0) {
      if (num_crossterm[m] == CMAPMAX)
        error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    if ((m = atom->map(atom4)) >= 0) {
      if (num_crossterm[m] == CMAPMAX)
        error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    if ((m = atom->map(atom5)) >= 0) {
      if (num_crossterm[m] == CMAPMAX)
        error->one(FLERR,"Too many CMAP crossterms for one atom");
      crossterm_type[m][num_crossterm[m]] = itype;
      crossterm_atom1[m][num_crossterm[m]] = atom1;
      crossterm_atom2[m][num_crossterm[m]] = atom2;
      crossterm_atom3[m][num_crossterm[m]] = atom3;
      crossterm_atom4[m][num_crossterm[m]] = atom4;
      crossterm_atom5[m][num_crossterm[m]] = atom5;
      num_crossterm[m]++;
    }

    buf = next + 1;
  }
}

/* ---------------------------------------------------------------------- */

bigint FixCMAP::read_data_skip_lines(char * /*keyword*/)
{
  return ncmap;
}

/* ----------------------------------------------------------------------
   write Mth header line to file
   only called by proc 0
------------------------------------------------------------------------- */

void FixCMAP::write_data_header(FILE *fp, int /*mth*/)
{
  fprintf(fp,BIGINT_FORMAT " cmap crossterms\n",ncmap);
}

/* ----------------------------------------------------------------------
   return size I own for Mth data section
   # of data sections = 1 for this fix
   nx = # of crossterms owned by my local atoms
     if newton_bond off, atom only owns crossterm if it is atom3
   ny = columns = type + 5 atom IDs
------------------------------------------------------------------------- */

void FixCMAP::write_data_section_size(int /*mth*/, int &nx, int &ny)
{
  int i,m;

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  nx = 0;
  for (i = 0; i < nlocal; i++)
    for (m = 0; m < num_crossterm[i]; m++)
      if (crossterm_atom3[i][m] == tag[i]) nx++;

  ny = 6;
}

/* ----------------------------------------------------------------------
   pack values for Mth data section into 2d buf
   buf allocated by caller as owned crossterms by 6
------------------------------------------------------------------------- */

void FixCMAP::write_data_section_pack(int /*mth*/, double **buf)
{
  int i,m;

  // 1st column = CMAP type
  // 2nd-6th columns = 5 atom IDs

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int n = 0;
  for (i = 0; i < nlocal; i++) {
    for (m = 0; m < num_crossterm[i]; m++) {
      if (crossterm_atom3[i][m] != tag[i]) continue;
      buf[n][0] = ubuf(crossterm_type[i][m]).d;
      buf[n][1] = ubuf(crossterm_atom1[i][m]).d;
      buf[n][2] = ubuf(crossterm_atom2[i][m]).d;
      buf[n][3] = ubuf(crossterm_atom3[i][m]).d;
      buf[n][4] = ubuf(crossterm_atom4[i][m]).d;
      buf[n][5] = ubuf(crossterm_atom5[i][m]).d;
      n++;
    }
  }
}

/* ----------------------------------------------------------------------
   write section keyword for Mth data section to file
   use Molecules or Charges if that is only field, else use fix ID
   only called by proc 0
------------------------------------------------------------------------- */

void FixCMAP::write_data_section_keyword(int /*mth*/, FILE *fp)
{
  fprintf(fp,"\nCMAP\n\n");
}

/* ----------------------------------------------------------------------
   write N lines from buf to file
   convert buf fields to int or double depending on styles
   index can be used to prepend global numbering
   only called by proc 0
------------------------------------------------------------------------- */

void FixCMAP::write_data_section(int /*mth*/, FILE *fp,
                                  int n, double **buf, int index)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT
            " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT "\n",
            index+i,(int) ubuf(buf[i][0]).i,(tagint) ubuf(buf[i][1]).i,
            (tagint) ubuf(buf[i][2]).i,(tagint) ubuf(buf[i][3]).i,
            (tagint) ubuf(buf[i][4]).i,(tagint) ubuf(buf[i][5]).i);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// methods for restart and communication
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixCMAP::write_restart(FILE *fp)
{
  if (comm->me == 0) {
    int size = sizeof(bigint);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&ncmap,sizeof(bigint),1,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixCMAP::restart(char *buf)
{
  ncmap = *((bigint *) buf);
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixCMAP::pack_restart(int i, double *buf)
{
  int n = 1;
  for (int m = 0; m < num_crossterm[i]; m++) {
    buf[n++] = ubuf(MAX(crossterm_type[i][m],-crossterm_type[i][m])).d;
    buf[n++] = ubuf(crossterm_atom1[i][m]).d;
    buf[n++] = ubuf(crossterm_atom2[i][m]).d;
    buf[n++] = ubuf(crossterm_atom3[i][m]).d;
    buf[n++] = ubuf(crossterm_atom4[i][m]).d;
    buf[n++] = ubuf(crossterm_atom5[i][m]).d;
  }
  // pack buf[0] this way because other fixes unpack it
  buf[0] = n;

  return n;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixCMAP::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

   int n = 0;
   for (int i = 0; i < nth; i++) n += static_cast<int> (extra[nlocal][n]);

   int count = static_cast<int> (extra[nlocal][n++]);
   num_crossterm[nlocal] = (count-1)/6;

   for (int m = 0; m < num_crossterm[nlocal]; m++) {
     crossterm_type[nlocal][m] = (int) ubuf(extra[nlocal][n++]).i;
     crossterm_atom1[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     crossterm_atom2[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     crossterm_atom3[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     crossterm_atom4[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
     crossterm_atom5[nlocal][m] = (tagint) ubuf(extra[nlocal][n++]).i;
   }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixCMAP::maxsize_restart()
{
  return 1 + CMAPMAX*6;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixCMAP::size_restart(int nlocal)
{
  return 1 + num_crossterm[nlocal]*6;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixCMAP::grow_arrays(int nmax)
{
  num_crossterm = memory->grow(num_crossterm,nmax,"cmap:num_crossterm");
  crossterm_type = memory->grow(crossterm_type,nmax,CMAPMAX,
                                "cmap:crossterm_type");
  crossterm_atom1 = memory->grow(crossterm_atom1,nmax,CMAPMAX,
                                 "cmap:crossterm_atom1");
  crossterm_atom2 = memory->grow(crossterm_atom2,nmax,CMAPMAX,
                                 "cmap:crossterm_atom2");
  crossterm_atom3 = memory->grow(crossterm_atom3,nmax,CMAPMAX,
                                 "cmap:crossterm_atom3");
  crossterm_atom4 = memory->grow(crossterm_atom4,nmax,CMAPMAX,
                                 "cmap:crossterm_atom4");
  crossterm_atom5 = memory->grow(crossterm_atom5,nmax,CMAPMAX,
                                 "cmap:crossterm_atom5");

  // must initialize num_crossterm to 0 for added atoms
  // may never be set for some atoms when data file is read

  for (int i = nmax_previous; i < nmax; i++) num_crossterm[i] = 0;
  nmax_previous = nmax;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixCMAP::copy_arrays(int i, int j, int /*delflag*/)
{
  num_crossterm[j] = num_crossterm[i];

  for (int k = 0; k < num_crossterm[j]; k++) {
    crossterm_type[j][k] = crossterm_type[i][k];
    crossterm_atom1[j][k] = crossterm_atom1[i][k];
    crossterm_atom2[j][k] = crossterm_atom2[i][k];
    crossterm_atom3[j][k] = crossterm_atom3[i][k];
    crossterm_atom4[j][k] = crossterm_atom4[i][k];
    crossterm_atom5[j][k] = crossterm_atom5[i][k];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixCMAP::set_arrays(int i)
{
  num_crossterm[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixCMAP::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = ubuf(num_crossterm[i]).d;
  for (int m = 0; m < num_crossterm[i]; m++) {
    buf[n++] = ubuf(crossterm_type[i][m]).d;
    buf[n++] = ubuf(crossterm_atom1[i][m]).d;
    buf[n++] = ubuf(crossterm_atom2[i][m]).d;
    buf[n++] = ubuf(crossterm_atom3[i][m]).d;
    buf[n++] = ubuf(crossterm_atom4[i][m]).d;
    buf[n++] = ubuf(crossterm_atom5[i][m]).d;
  }
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixCMAP::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  num_crossterm[nlocal] = (int) ubuf(buf[n++]).i;
  for (int m = 0; m < num_crossterm[nlocal]; m++) {
    crossterm_type[nlocal][m] = (int) ubuf(buf[n++]).i;
    crossterm_atom1[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    crossterm_atom2[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    crossterm_atom3[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    crossterm_atom4[nlocal][m] = (tagint) ubuf(buf[n++]).i;
    crossterm_atom5[nlocal][m] = (tagint) ubuf(buf[n++]).i;
  }
  return n;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixCMAP::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = (double)nmax * sizeof(int);        // num_crossterm
  bytes += (double)nmax*CMAPMAX * sizeof(int);      // crossterm_type
  bytes += (double)5*nmax*CMAPMAX * sizeof(int);    // crossterm_atom 12345
  bytes += (double)maxcrossterm*6 * sizeof(int);    // crosstermlist
  return bytes;
}
