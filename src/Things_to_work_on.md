# TILD Package to LAMMPS by Christian Tabedzki and Zachariah Vicars

## Things to be implemented

* Insuring the fftw and PPPM scheme have the same points so we can do a
  convolution
  * Figuring out how to properly break up the grid so that the same amount of
    points are calculated for both.
* Determine how grid point (0 to N-1) corresponds to the xyz coordinates to
  properly map the points for FFTW
* Figure out kspace_settings to figure out how to do the values

## Implementation Questions

* How does one go about setting up multiple field interactions?
  * Should there be two different types of parameter lines?
* What function oconverts between the values that the current processor stores
  from 0 to N-1 and the grid values that correspond to it?
