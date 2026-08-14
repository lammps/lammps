Kspace styles
=============

Classes that compute long-range Coulombic interactions via K-space
representations (Ewald, PPPM) are derived from the KSpace class.  New
styles can be created to add new K-space options to LAMMPS.

``src/ewald.cpp`` is an example of computing K-space interactions.

Here is a brief description of methods you define in your new derived
class.  See ``src/kspace.h`` for details.

+------------------------------+------------------------------------------------------------------+
| Required                     | "pure" methods that *must* be overridden in a derived class      |
+==============================+==================================================================+
| init                         | initialize the calculation before a run                          |
+------------------------------+------------------------------------------------------------------+
| setup                        | computation before the first timestep of a run                   |
+------------------------------+------------------------------------------------------------------+
| compute                      | every-timestep computation                                       |
+------------------------------+------------------------------------------------------------------+

+------------------------------+------------------------------------------------------------------+
| Optional                     | methods that have a default or empty implementation              |
+==============================+==================================================================+
| settings                     | process the arguments to the kspace_style command                |
+------------------------------+------------------------------------------------------------------+
| compute_group_group          | compute energy/force between two groups of atoms                 |
+------------------------------+------------------------------------------------------------------+
| qsum_qsq                     | compute total charge and charge-squared sum                      |
+------------------------------+------------------------------------------------------------------+
| reset_grid                   | called when the FFT grid size changes                            |
+------------------------------+------------------------------------------------------------------+
| pack_forward_grid            | pack FFT grid data into a buffer for forward communication       |
+------------------------------+------------------------------------------------------------------+
| unpack_forward_grid          | unpack FFT grid data from a forward communication buffer         |
+------------------------------+------------------------------------------------------------------+
| pack_reverse_grid            | pack FFT grid data into a buffer for reverse communication       |
+------------------------------+------------------------------------------------------------------+
| unpack_reverse_grid          | unpack FFT grid data from a reverse communication buffer         |
+------------------------------+------------------------------------------------------------------+
| timing                       | report timing information (used by kspace_modify)                |
+------------------------------+------------------------------------------------------------------+
| timing_1d                    | report 1d FFT timing                                             |
+------------------------------+------------------------------------------------------------------+
| timing_3d                    | report 3d FFT timing                                             |
+------------------------------+------------------------------------------------------------------+
| modify_param                 | handle arguments from the kspace_modify command                  |
+------------------------------+------------------------------------------------------------------+
| memory_usage                 | tally of memory usage                                            |
+------------------------------+------------------------------------------------------------------+

Here is a list of flags or settings that should be set in the
constructor of the derived class when they differ from the default
setting.

+---------------------------------+-------------------------------------------------------------+---------+
| Name of flag                    | Description                                                 | default |
+=================================+=============================================================+=========+
| ewaldflag                       | 1 if this is an Ewald solver                                | 0       |
+---------------------------------+-------------------------------------------------------------+---------+
| pppmflag                        | 1 if this is a PPPM solver                                  | 0       |
+---------------------------------+-------------------------------------------------------------+---------+
| msmflag                         | 1 if this is an MSM solver                                  | 0       |
+---------------------------------+-------------------------------------------------------------+---------+
| dispersionflag                  | 1 if this is a LJ/dispersion solver                         | 0       |
+---------------------------------+-------------------------------------------------------------+---------+
| tip4pflag                       | 1 if this is a TIP4P solver                                 | 0       |
+---------------------------------+-------------------------------------------------------------+---------+
| dipoleflag                      | 1 if this is a dipole solver                                | 0       |
+---------------------------------+-------------------------------------------------------------+---------+
| spinflag                        | 1 if this is a spin solver                                  | 0       |
+---------------------------------+-------------------------------------------------------------+---------+
| triclinic_support               | 1 if the style supports triclinic boxes                     | 1       |
+---------------------------------+-------------------------------------------------------------+---------+
| group_group_enable              | 1 if compute_group_group() is implemented                   | 0       |
+---------------------------------+-------------------------------------------------------------+---------+
| centroidstressflag              | CENTROID_SAME/CENTROID_AVAIL/CENTROID_NOTAVAIL              | NOTAV.  |
+---------------------------------+-------------------------------------------------------------+---------+

The pair-style flags ``ewaldflag``, ``pppmflag``, ``msmflag``,
``dispersionflag``, ``tip4pflag``, ``dipoleflag``, and ``spinflag``
must match the corresponding flags set in the pair style that is used
together with the kspace style.  The default value of
``centroidstressflag`` is ``CENTROID_NOTAVAIL`` for kspace styles, as
the centroid stress for long-range interactions is generally not equal
to the two-body virial.
