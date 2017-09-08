This directory contains other versions of the same molecule
(with the same force-field), implemented in different ways.

charmm/1beadProtSci2010.lt   <-- This applies multiple "charmm" dihedral angle
                                 forces to the same quartet of atoms to create
                                 a Fourier series. (No packages needed.)

                                 NOTE: You must run moltemplate this way:

                                 moltemplate.sh -overlay-dihdedrals system.lt

class2/1beadProtSci2010.lt   <-- This uses the "class2" dihedral angles forces
                                 (You must build LAMMPS with the CLASS2 package)
