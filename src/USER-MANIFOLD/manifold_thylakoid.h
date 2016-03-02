#ifndef LMP_MANIFOLD_THYLAKOID_H
#define LMP_MANIFOLD_THYLAKOID_H

#include "manifold.h"

#include <cassert>
#include <vector>

#include "manifold_thylakoid_shared.h"

namespace LAMMPS_NS {

  using LAMMPS_NS::thyla_part;
  using LAMMPS_NS::thyla_part_geom;



  class manifold_thylakoid : public manifold, protected Pointers {
   public:
    enum { NPARAMS = 3 };
    manifold_thylakoid( LAMMPS *lmp, int, char ** );
    virtual ~manifold_thylakoid();
    
    virtual double g( const double *x );
    virtual void   n( const double *x, double *n );

    static const char* type(){ return "thylakoid"; }
    virtual const char *id(){ return type(); }
    static int expected_argc(){ return NPARAMS; }	
    virtual int nparams(){ return NPARAMS; }

    
    virtual void post_param_init();
    virtual void checkup(); // Some diagnostics...
   private:
    void init_domains();

    thyla_part *get_thyla_part( const double *x, int *err_flag, std::size_t *idx = NULL );
    int is_in_domain( thyla_part *p, const double *x );
    void check_overlap();
    std::vector<thyla_part*> parts;

    thyla_part *make_plane_part (double a, double b, double c,
                                 const std::vector<double> &pt);
    thyla_part *make_cyl_part   (double a, double b, double c,
                                 const std::vector<double> &pt, double R);
    thyla_part *make_sphere_part(const std::vector<double> &pt, double R);
    thyla_part *make_cyl_to_plane_part(double X0, double R0, double R, double s,
                                       const std::vector<double> &pt );


    void set_domain( thyla_part *p, const std::vector<double> &lo,
                     const std::vector<double> &hi );

    void add_part( thyla_part *part, double x0, double y0,
                   double z0, double x1, double y1, double z1 )
    {
      int me = -1;
      MPI_Comm_rank(world, &me);

      std::vector<double> lo(3), hi(3);
      lo[0] = x0; lo[1] = y0; lo[2] = z0;
      hi[0] = x1; hi[1] = y1; hi[2] = z1;
      
      set_domain(part, lo, hi);
      parts.push_back( part );
      
      fprintf(my_screen,"This is proc %d: Added part %p of type %d, "
              "size is now %lu.\n", me, part, part->type, parts.size());
      fflush(my_screen);
    }
	  
    
    void print_part_data( FILE *fp_doms, FILE *fp_coms );

    void assert_part_ok( const thyla_part *part ){
      if( !part ) error->one(FLERR,"Error! Part was NULL!");
    }

	  

    // Coefficients for the thylakoid model. At the moment it is just
    // a cylinder, we slowly expand it.
    double pad;  // Padding to make sure periodic images are mapped back properly.
    FILE *my_screen;


    double LB, lT, lB, wB, LT;

    // Domain size:
    double x0, y0, z0;
    double x1, y1, z1;
    double Lx, Ly, Lz;
  };



} // namespace LAMMPS_NS


#endif // LMP_MANIFOLD_THYLAKOID_H
