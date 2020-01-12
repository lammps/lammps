#include "manifold_thylakoid.h"
#include <cmath>
#include "manifold_thylakoid_shared.h"
#include "comm.h"
#include "domain.h" // For some checks regarding the simulation box.
#include "error.h"


#define MANIFOLD_THYLAKOID_DEBUG

using namespace LAMMPS_NS;
using namespace user_manifold;


manifold_thylakoid::manifold_thylakoid( LAMMPS *lmp, int /*narg*/, char ** /*arg*/)
  : manifold(lmp)
{
  // You can NOT depend on proper construction of the domains in
  // the constructor, because the params are set by the function that
  // calls the factory-construction method. Instead, the constructing
  // fix should call post_param_init();
}



manifold_thylakoid::~manifold_thylakoid()
{
  for( std::size_t i = 0; i < parts.size(); ++i ){
    delete parts[i];
  }
}



void manifold_thylakoid::post_param_init()
{
  // Set coefficients:
  lT  = 3;  // Radius of cylinder edges of stacks
  LT  = 15; // Size of faces of cylinder stacks
  pad = 3;  // Padding to prevent improper look-ups
  wB  = 3.0;
  LB = 10.0;
  lB = 3.0;

  wB = params[0];
  LB = params[1];
  lB = params[2];

  if (comm->me == 0) {
    fprintf(screen,"My params are now: lT = %f, LT = %f, pad = %f, "
            "wB = %f, LB = %f, lB = %f\n", lT, LT, pad, wB, LB, lB );
    fprintf(screen,"Calling init_domains() from post_param_init().\n");
  }
  init_domains();
  checkup();
}

void manifold_thylakoid::checkup()
{
  if (comm->me == 0 ) {
    fprintf(screen,"This is checkup of thylakoid %p\n", this);
    fprintf(screen,"I have %ld parts. They are:\n", parts.size());
    for( int i = 0; i < parts.size(); ++i ){
      fprintf(screen, "[%f, %f] x [%f, %f] x [%f, %f]\n",
              parts[i]->xlo, parts[i]->xhi,
              parts[i]->ylo, parts[i]->yhi,
              parts[i]->zlo, parts[i]->zhi );
    }
    fprintf(screen,"My params are:\n");
    for( int i = 0; i < NPARAMS; ++i ){
      fprintf(screen,"%f\n", params[i]);
    }
  }
}


double manifold_thylakoid::g( const double *x )
{
  int err = 0;
  std::size_t idx;
  thyla_part *p = get_thyla_part(x,&err,&idx);
  if(err){
    char msg[2048];
    sprintf(msg,"Error getting thyla_part for x = (%f, %f, %f)",x[0],x[1],x[2]);
    error->one(FLERR,msg);
  }
  double con_val = p->g(x);
  if (std::isfinite(con_val)) {
    return con_val;
  } else {
    char msg[2048];
    sprintf(msg,"Error, thyla_part of type %d returned %f as constraint val!",
            p->type, con_val);
    error->one(FLERR,msg);
    return 0;
  }
}

void   manifold_thylakoid::n( const double *x, double *n )
{
  int err = 0;
  std::size_t idx;
  thyla_part *p = get_thyla_part(x,&err,&idx);
  if(err){
    char msg[2048];
    sprintf(msg,"Error getting thyla_part for x = (%f, %f, %f)",x[0],x[1],x[2]);
    error->one(FLERR,msg);
  }
  p->n(x,n);
  if (std::isfinite(n[0]) && std::isfinite(n[1]) && std::isfinite(n[2])) {
    return;
  } else {
    char msg[2048];
    sprintf(msg,"Error, thyla_part of type %d returned (%f,%f,%f) as gradient!",
            p->type, n[0], n[1], n[2]);
    error->one(FLERR,msg);
  }
}

thyla_part *manifold_thylakoid::get_thyla_part( const double *x, int * /*err_flag*/, std::size_t *idx )
{

  for( std::size_t i = 0; i < parts.size(); ++i ){
    thyla_part *p = parts[i];
    if (is_in_domain(p,x)) {
      if (idx != NULL) *idx = i;
      return p;
    }
  }
  char msg[2048];
  sprintf(msg,"Could not find thyla_part for x = (%f,%f,%f)", x[0],x[1],x[2]);
  error->one(FLERR,msg);
  return NULL;
}






void manifold_thylakoid::init_domains()
{
  if (wB + 2*lB > LT) {
    char msg[2048];
    sprintf(msg,"LT = %f not large enough to accomodate bridge with "
            "wB = %f and lB = %f! %f > %f\n", LT, wB, lB, wB + 2*lB, LT);
    error->one(FLERR,msg);
  }

  // Determine some constant coordinates:
  x0 = -( 0.5*LB + lB + lT + LT + lT + pad);
  y0 = -( 0.5*LT + lT + pad );
  z0 = -15;
#ifdef  MANIFOLD_THYLAKOID_DEBUG
  if (comm->me == 0) {
    fprintf(screen,"x0, y0, z0 = %f, %f, %f\n",x0,y0,z0);
  }
#endif // MANIFOLD_THYLAKOID_DEBUG

#ifndef USE_PHONY_LAMMPS
  if (x0 < domain->boxlo[0]) {
    char msg[2048];
    sprintf(msg,"Thylakoid expects xlo of at most %f, but found %f",
            x0, domain->boxlo[0]);
    error->one(FLERR,msg);
  }
  if (y0 < domain->boxlo[1]) {
    char msg[2048];
    sprintf(msg,"Thylakoid expects ylo of at most %f, but found %f",
            y0, domain->boxlo[1]);
    error->one(FLERR,msg);
  }
#endif

  // Add some padding to prevent improper lookups.
  z0 -= pad;

  x1 = -x0;
  y1 = -y0;
  z1 = -z0;

  Lx = x1 - x0;
  Ly = y1 - y0;
  Lz = z1 - z0;

#ifndef USE_PHONY_LAMMPS
  char msg[2048];
  if(x1 > domain->boxhi[0]){
    sprintf(msg,"Expected xhi larger than current box has: %f > %f",
            x1, domain->boxhi[0]);
    error->one(FLERR,msg);
  }
  if(y1 > domain->boxhi[1]){
    sprintf(msg,"Expected yhi larger than current box has: %f > %f",
            y1, domain->boxhi[1]);
    error->one(FLERR,msg);
  }
  // if(z1 > domain->boxhi[2]){
  //   sprintf(msg,"Expected zhi larger than current box has: %f > %f",
  //           z1, domain->boxhi[2]);
  //   error->one(FLERR,msg);
  // }
#endif

  // Create and add the manifold parts to the array.
  thyla_part *p;


  // Determine coordinates of domain boundaries and centres of "mass":
  thyla_part_geom cllb, cllt, clrb, clrt;      // Left thylakoid cylinder parts
  thyla_part_geom pll, plb, plt, plr;          // Left thylakoid plane parts

  thyla_part_geom crlb, crlt, crrb, crrt;      // Right thylakoid cylinder parts
  thyla_part_geom prl, prb, prt, prr;          // Right thylakoid plane parts

  thyla_part_geom bl, br, bc; // Bridge left, right connectors and cylinder.





  // The bridge is three parts.
  // 1. A connector between the right face of the left  grana and a cylinder
  // 2. A connector between the left  face of the right grana and a cylinder
  // 3. The aforementioned cylinder.

  // 1.:
  // Args: pt, X0, R0, R, x0, y0, z0, sign
  double rB = 0.5*wB;
  double Reff = rB + lB;
  double safety_fac = 1.0;
  bl.pt[0] = -0.5*LB;
  bl.pt[1] = 0;
  bl.pt[2] = 0;

  bl.lo[0] = bl.pt[0] - safety_fac*lB;
  bl.lo[1] = -(1.0 + safety_fac) * Reff;
  bl.lo[2] = -(1.0 + safety_fac) * Reff;

  bl.hi[0] = bl.pt[0];
  bl.hi[1] = (1.0 + safety_fac) * Reff;
  bl.hi[2] = (1.0 + safety_fac) * Reff;

  // double X0, double R0, double R, double s,
#ifdef MANIFOLD_THYLAKOID_DEBUG
  if (comm->me == 0) {
    fprintf(screen,"x0, r0, R = %f, %f, %f\n", bl.pt[0], rB, lB);
  }
#endif // MANIFOLD_THYLAKOID_DEBUG
  p = make_cyl_to_plane_part(bl.pt[0], rB, lB, -1, bl.pt);
  set_domain(p, bl.lo, bl.hi);
  parts.push_back(p);

  // 2.:
  br.pt[0] = 0.5*LB;
  br.pt[1] = 0;
  br.pt[2] = 0;

  br.lo[0] = br.pt[0];
  br.lo[1] = -(1.0 + safety_fac) * Reff;
  br.lo[2] = -(1.0 + safety_fac) * Reff;

  br.hi[0] = br.pt[0] + safety_fac*lB;
  br.hi[1] = (1.0 + safety_fac) * Reff;
  br.hi[2] = (1.0 + safety_fac) * Reff;

  // double X0, double R0, double R, double s,
#ifdef MANIFOLD_THYLAKOID_DEBUG
  if (comm->me == 0) {
    fprintf(screen,"x0, r0, R = %f, %f, %f\n", br.pt[0], rB, lB);
  }
#endif // MANIFOLD_THYLAKOID_DEBUG
  p = make_cyl_to_plane_part(br.pt[0], rB, lB,  1, br.pt);
  set_domain(p, br.lo, br.hi);
  parts.push_back(p);



  // 3.:
  // Cylinder in between:
  bc.pt[0] = 0;
  bc.pt[1] = 0;
  bc.pt[2] = 0;

  bc.lo[0] = bl.pt[0];
  bc.lo[1] = -Reff;
  bc.lo[2] = -Reff;

  bc.hi[0] = br.pt[0];
  bc.hi[1] = Reff;
  bc.hi[2] = Reff;

  p = make_cyl_part( 0, 1, 1, bc.pt, rB );
  set_domain( p, bc.lo, bc.hi );
#ifdef MANIFOLD_THYLAKOID_DEBUG
  if (comm->me == 0) {
    fprintf(screen,"Cylinder lives on [ %f x %f ] x [ %f x %f ] x [ %f x %f]\n",
            bc.lo[0], bc.hi[0], bc.lo[1], bc.hi[1], bc.lo[2], bc.hi[2]);
  }
#endif // MANIFOLD_THYLAKOID_DEBUG

  parts.push_back(p);


  // The stack on the left:
  cllb.lo[0] = x0;
  cllb.lo[1] = y0;
  cllb.lo[2] = z0;
  cllb.pt[0] = x0 + pad + lT;
  cllb.pt[1] = y0 + pad + lT;
  cllb.pt[2] = 0;
  cllb.hi[0] = cllb.pt[0];
  cllb.hi[1] = cllb.pt[1];
  cllb.hi[2] = z1;

  p = make_cyl_part(1,1,0,cllb.pt,lT);
  set_domain(p, cllb.lo, cllb.hi);
  parts.push_back(p);

  // left left top cylinder
  cllt = cllb;
  cllt.lo[1] = y1 - pad - lT;
  cllt.hi[1] = y1;
  cllt.pt[1] = cllb.pt[1] + LT;

  p = make_cyl_part(1,1,0,cllt.pt,lT);
  set_domain(p, cllt.lo, cllt.hi);
  parts.push_back(p);

  // left right bottom cylinder
  clrb = cllb;
  clrb.pt[0] += LT;
  clrb.lo[0] = clrb.pt[0];
  clrb.hi[0] = clrb.lo[0] + lT + lB;

  p = make_cyl_part(1,1,0,clrb.pt,lT);
  set_domain(p, clrb.lo, clrb.hi);
  parts.push_back(p);

  // left right top cylinder
  clrt = clrb;
  clrt.pt[1] += LT;
  clrt.lo[1] = y1 - pad - lT;
  clrt.hi[1] = y1;

  p = make_cyl_part(1,1,0,clrt.pt,lT);
  set_domain(p, clrt.lo, clrt.hi);
  parts.push_back(p);


  // left left plane
  pll.pt[0] = x0 + pad;
  pll.pt[1] = 0;
  pll.pt[2] = 0;
  pll.lo[0] = x0;
  pll.lo[1] = cllb.pt[1];
  pll.lo[2] = z0;
  pll.hi[0] = pll.lo[0] + pad + lT;
  pll.hi[1] = pll.lo[1] + LT;
  pll.hi[2] = z1;

  p = make_plane_part(1,0,0,pll.pt);
  set_domain(p, pll.lo, pll.hi);
  parts.push_back(p);

  // left bottom plane
  plb.pt[0] = x0 + pad + lT + 0.5*LT;
  plb.pt[1] = y0 + pad;
  plb.pt[2] = 0;
  plb.lo[0] = x0 + pad + lT;
  plb.lo[1] = y0;
  plb.lo[2] = z0;
  plb.hi[0] = plb.lo[0] + LT;
  plb.hi[1] = plb.lo[1] + pad + lT;
  plb.hi[2] = z1;

  p = make_plane_part(0,1,0,plb.pt);
  set_domain(p, plb.lo, plb.hi);
  parts.push_back(p);

  // left top plane
  plt = plb;
  plt.lo[1] = cllb.pt[1] + LT;
  plt.hi[1] = y1;
  plt.pt[1] = y1 - pad;

  p = make_plane_part(0,1,0,plt.pt);
  set_domain(p, plt.lo, plt.hi);
  parts.push_back(p);

  // left right plane
  plr = pll;
  plr.lo[0] = bl.lo[0] - 0.5;
  plr.lo[1] = y0 - pad;
  plr.hi[0] = bl.lo[0] + 0.5;
  plr.hi[1] = y1 + pad;
  plr.pt[0] = bl.pt[0] - lB;
  plr.pt[1] = 0.0;
  plr.pt[2] = 0.0;
  plr.hi[2] = z1 + pad;
  plr.lo[2] = z0 - pad;

  p = make_plane_part(1,0,0,plr.pt);
  set_domain(p, plr.lo, plr.hi);
  parts.push_back(p);

  // Check if this plane lines up with bl:
  if (fabs(plr.pt[0] - bl.pt[0] + lB) > 1e-8) {
    char msg[2048];
    sprintf(msg,"Origins of plane left right and bridge left misaligned! %f != %f!\n",
            plr.pt[0], bl.pt[0] - lB );
    error->one(FLERR,msg);
  }

  // Now, for the right stack, you can mirror the other...
  // To mirror them you need to invert lo[0] and hi[0] and flip their sign.

  thyla_part_geom::mirror(thyla_part_geom::DIR_X, &crlb, &cllb);
  p = make_cyl_part(1,1,0,crlb.pt,lT);
  set_domain(p, crlb.lo, crlb.hi);
  parts.push_back(p);

  thyla_part_geom::mirror(thyla_part_geom::DIR_X, &crlt, &cllt);
  p = make_cyl_part(1,1,0,crlt.pt,lT);
  set_domain(p, crlt.lo, crlt.hi);
  parts.push_back(p);

  thyla_part_geom::mirror(thyla_part_geom::DIR_X, &crrb, &clrb);
  p = make_cyl_part(1,1,0,crrb.pt,lT);
  set_domain(p, crrb.lo, crrb.hi);
  parts.push_back(p);

  thyla_part_geom::mirror(thyla_part_geom::DIR_X, &crrt, &clrt);
  p = make_cyl_part(1,1,0,crrt.pt,lT);
  set_domain(p, crrt.lo, crrt.hi);
  parts.push_back(p);



  thyla_part_geom::mirror(thyla_part_geom::DIR_X, &prl , &pll );
  p = make_plane_part(1,0,0,prl.pt);
  set_domain(p, prl.lo, prl.hi);
  parts.push_back(p);

  thyla_part_geom::mirror(thyla_part_geom::DIR_X, &prb , &plb );
  p = make_plane_part(0,1,0,prb.pt);
  set_domain(p, prb.lo, prb.hi);
  parts.push_back(p);

  thyla_part_geom::mirror(thyla_part_geom::DIR_X, &prt , &plt );
  p = make_plane_part(0,1,0,prt.pt);
  set_domain(p, prt.lo, prt.hi);
  parts.push_back(p);

  // Careful, this one is wrongly named.
  thyla_part_geom::mirror(thyla_part_geom::DIR_X, &prr, &plr);
  p = make_plane_part(1,0,0,prr.pt);
  set_domain(p, prr.lo, prr.hi);
  parts.push_back(p);

  if (fabs(prr.pt[0] - br.pt[0] - lB) > 1e-8) {
    char msg[2048];
    sprintf(msg,"Origins of plane left right and bridge left misaligned! %f != %f!\n",
            prr.pt[0], br.pt[0] + lB);
    error->one(FLERR,msg);
  }

  // For debugging, print the domains and coms:
#ifdef MANIFOLD_THYLAKOID_DEBUG
  if (comm->me == 0) {
    FILE *fp_doms = fopen("test_doms.dat","w");
    FILE *fp_coms = fopen("test_coms.dat","w");
    print_part_data(fp_doms, fp_coms);
    fclose(fp_doms);
    fclose(fp_coms);
  }
#endif // MANIFOLD_THYLAKOID_DEBUG
}


void manifold_thylakoid::set_domain( thyla_part *p, const std::vector<double> &lo,
                                     const std::vector<double> &hi )
{
#ifdef MANIFOLD_THYLAKOID_DEBUG
  if (comm->me == 0) {
    fprintf(screen,"Adding part with domain [%f, %f] x [%f, %f] x [%f, %f]\n",
            lo[0],hi[0],lo[1],hi[1],lo[2],hi[2] );
  }
#endif // MANIFOLD_THYLAKOID_DEBUG
  if (lo[0] >= hi[0]) {
    char msg[2048];
    sprintf(msg,"xlo >= xhi (%f >= %f)",lo[0],hi[0]);
    error->one(FLERR,msg);
  } else if (lo[1] >= hi[1]) {
    char msg[2048];
    sprintf(msg,"ylo >= yhi (%f >= %f)",lo[1],hi[1]);
    error->one(FLERR,msg);
  } else if (lo[2] >= hi[2]) {
    char msg[2048];
    sprintf(msg,"zlo >= zhi (%f >= %f)",lo[2],hi[2]);
    error->one(FLERR,msg);
  }
  p->xlo = lo[0];
  p->ylo = lo[1];
  p->zlo = lo[2];

  p->xhi = hi[0];
  p->yhi = hi[1];
  p->zhi = hi[2];
}

int manifold_thylakoid::is_in_domain( thyla_part *part, const double *x )
{
  bool domain_ok = (x[0] >= part->xlo) && (x[0] <= part->xhi) &&
          (x[1] >= part->ylo) && (x[1] <= part->yhi) &&
          (x[2] >= part->zlo) && (x[2] <= part->zhi);

  if (!domain_ok) return false;

  // From here on out, domain is ok.

  if (part->type == thyla_part::THYLA_TYPE_CYL_TO_PLANE) {

    double R0 = part->params[1];
    double R  = part->params[2];
    double y = x[1];
    double z = x[2];
    double dist2 = y*y + z*z;
    double RR = R0+R;
    double RR2 = RR*RR;


    if (dist2 < RR2) {
      return true;
    } else {
      // Domain was ok, but radius not.
      return false;
    }
  } else {
    return true;
  }
}


thyla_part *manifold_thylakoid::make_plane_part (double a, double b, double c,
                                                 const std::vector<double> &pt )
{
  double args[7];
  args[0] = a;
  args[1] = b;
  args[2] = c;
  args[3] = pt[0];
  args[4] = pt[1];
  args[5] = pt[2];
  thyla_part *p = new thyla_part(thyla_part::THYLA_TYPE_PLANE,args,0,0,0,0,0,0);
  return p;
}

thyla_part *manifold_thylakoid::make_cyl_part   (double a, double b, double c,
                                                 const std::vector<double> &pt, double R)
{
  double args[7];
  args[0] = a;
  args[1] = b;
  args[2] = c;
  args[3] = pt[0];
  args[4] = pt[1];
  args[5] = pt[2];
  args[6] = R;
  thyla_part *p = new thyla_part(thyla_part::THYLA_TYPE_CYL,args,0,0,0,0,0,0);
  return p;
}


thyla_part *manifold_thylakoid::make_sphere_part(const std::vector<double> &pt, double R)
{
  double args[7];
  args[0] = R;
  args[1] = pt[0];
  args[2] = pt[1];
  args[3] = pt[2];
  thyla_part *p = new thyla_part(thyla_part::THYLA_TYPE_SPHERE,args,0,0,0,0,0,0);
  return p;
}


thyla_part *manifold_thylakoid::make_cyl_to_plane_part(double X0, double R0, double R,
                                                       double s, const std::vector<double> &pt )
{
  double args[7];
  args[0] = X0;
  args[1] = R0;
  args[2] = R;
  args[3] = pt[0];
  args[4] = pt[1];
  args[5] = pt[2];
  args[6] = s;
  thyla_part *p = new thyla_part(thyla_part::THYLA_TYPE_CYL_TO_PLANE,args,0,0,0,0,0,0);
  return p;
}




void manifold_thylakoid::print_part_data( FILE *fp_doms, FILE *fp_coms )
{
  for( std::size_t i = 0; i < parts.size(); ++i ){
    thyla_part *p = parts[i];
    fprintf(fp_doms, "%f   %f\n",  p->xlo, p->ylo);
    fprintf(fp_doms, "%f   %f\n",  p->xlo, p->yhi);
    fprintf(fp_doms, "%f   %f\n",  p->xhi, p->yhi);
    fprintf(fp_doms, "%f   %f\n",  p->xhi, p->ylo);
    fprintf(fp_doms, "%f   %f\n\n",p->xlo, p->ylo);
    fprintf(fp_coms, "%f   %f\n",  p->x0, p->y0 );
  }
}


