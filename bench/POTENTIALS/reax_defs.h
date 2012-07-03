#define PORTABLECOMMENTFLAG
#ifndef PORTABLECOMMENTFLAG
// This is just a way to have portable comments 
// for both C++ and FORTRAN preprocessing.
 /* ///:EOH~                                                                 */
 /*                                                                          */
 /* This file contains array dimension parameters for all the main           */
 /* ReaxFF data structures, some of which need to be directly accessed       */
 /* by Grasp C++ functions. If they are set too small, the calculation       */
 /* will run out of allocated memory. If they are set too big, the machine   */
 /* will not be able to allocate enough memory.                              */
 /*                                                                          */
 
 /*     NNEIGHMAXDEF =   Max number of neighbors / NATDEF                    */
 /*     NATDEF =         Max number of atoms                                 */
 /*     NATTOTDEF =      Max number of global atoms                          */
 /*     NSORTDEF =       Max number of atom types                            */
 /*     MBONDDEF =       Max number of bonds connected to one atom           */
 /*     NAVIBDEF =       for 2nd derivatives                                 */
 /*     NBOTYMDEF =      Max number of bond types                            */
 /*     NVATYMDEF =      Max number of valency angle types                   */
 /*     NTOTYMDEF =      Max number of torsion angle types                   */
 /*     NHBTYMDEF =      Max number of hydrogen bond types                   */
 /*     NODMTYMDEF =     Max number of off-diagonal Morse types              */
 /*     NBOALLMAXDEF =   Max number of all bonds                             */
 /*     NBOMAXDEF =      Max number of bonds                                 */
 /*     NHBMAXDEF =      Max number of  hydrogen bonds                       */
 /*     NVAMAXDEF =      Max number of valency angles                        */
 /*     NOPMAXDEF =      Max number of out of plane angles                   */
 /*     NTOMAXDEF =      Max number of torsion angles                        */
 /*     NPAMAXDEF =      Max number of general parameters in force field     */
 /*     NMOLMAXDEF =     Max number of molecules in system                   */
 /*     NMOLSETDEF =     Max number of molecules in training set             */
 /*     MRESTRADEF =     Max number of restraints                            */
 /*     MTREGDEF =       Max number of temperature regimes                   */
 /*     MTZONEDEF =      Max number of temperature zones                     */
 /*     MVREGDEF =       Max number of volume regimes                        */
 /*     MVZONEDEF =      Max number of volume zones                          */
 /*     MEREGDEF =       Max number of electric field regimes                */
 /*     MEZONEDEF =      Max number of electric field zones                  */
#endif

#define NNEIGHMAXDEF 200
#define NATDEF 50000
#define NATTOTDEF 1
#define NSORTDEF 20
#define MBONDDEF 20
#define NAVIBDEF 50
#define NBOTYMDEF 200
#define NVATYMDEF 200
#define NTOTYMDEF 200
#define NHBTYMDEF 200
#define NODMTYMDEF 20
#define NBOALLMAXDEF 180000
#define NBOMAXDEF 90000
#define NHBMAXDEF 400000
#define NVAMAXDEF 300000
#define NOPMAXDEF 00010
#define NTOMAXDEF 200000
#define NPAMAXDEF 50
#define NMOLMAXDEF 2000
#define NMOLSETDEF 1500
#define MRESTRADEF 100
#define MTREGDEF 100
#define MTZONEDEF 5
#define MVREGDEF 100
#define MVZONEDEF 6
#define MEREGDEF 100
#define MEZONEDEF 3
