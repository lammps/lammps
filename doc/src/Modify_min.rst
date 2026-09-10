Minimization styles
===================

Classes that perform energy minimization derived from the Min class.
New styles can be created to add new minimization algorithms to
LAMMPS.

``src/min_cg.cpp`` is an example of conjugate gradient minimization.

Here is a brief description of methods you define in your new derived
class.  See ``src/min.h`` for details.

+-------------------+----------------------------------------------------------+
| Required          | "pure" methods that *must* be overridden                 |
+===================+==========================================================+
| setup_style       | initialize minimizer-specific data structures            |
+-------------------+----------------------------------------------------------+
| reset_vectors     | connect internal xvec/fvec pointers to atom coordinates  |
+-------------------+----------------------------------------------------------+
| iterate           | perform the minimization iterations; returns stop-code   |
+-------------------+----------------------------------------------------------+

+-------------------+----------------------------------------------------------+
| Optional          | methods that have a default or base-class implementation |
+===================+==========================================================+
| init              | initialize the minimization before a run                 |
+-------------------+----------------------------------------------------------+
| setup             | called before minimization starts                        |
+-------------------+----------------------------------------------------------+
| setup_minimal     | minimal setup variant used for restarts                  |
+-------------------+----------------------------------------------------------+
| init_style        | check style-specific conditions at init time             |
+-------------------+----------------------------------------------------------+
| run               | called by LAMMPS to drive the full minimization          |
+-------------------+----------------------------------------------------------+
| force_clear       | clear forces; override for GPU/Kokkos variants           |
+-------------------+----------------------------------------------------------+
| energy_force      | evaluate energy and forces; override for Kokkos          |
+-------------------+----------------------------------------------------------+
| fnorm_sqr         | return the squared force norm; override for spin styles  |
+-------------------+----------------------------------------------------------+
| fnorm_inf         | return the inf-norm of forces; override for spin styles  |
+-------------------+----------------------------------------------------------+
| fnorm_max         | return the max force norm; override for spin styles      |
+-------------------+----------------------------------------------------------+
| modify_param      | handle arguments from the min_modify command             |
+-------------------+----------------------------------------------------------+
| memory_usage      | tally of memory usage                                    |
+-------------------+----------------------------------------------------------+

The ``iterate()`` method must return one of the following integer
stop-condition codes defined as an enum in ``src/min.h``:
``MAXITER``, ``MAXEVAL``, ``ETOL``, ``FTOL``, ``DOWNHILL``,
``ZEROALPHA``, ``ZEROFORCE``, ``ZEROQUAD``, ``TRSMALL``,
``INTERROR``, ``TIMEOUT``, or ``MAXVDOTF``.

The line-search style used by the minimizer is controlled by the
``linestyle`` member, which takes values from the enum
``{ BACKTRACK, QUADRATIC, FORCEZERO, SPIN_CUBIC, SPIN_NONE }``.
