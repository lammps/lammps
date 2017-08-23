A Menger cube is a fractal composed of subunits that resemble
a Rubik's-cubes with a hollow interior:
      _________
     /        /|
     ######### |
     # ## ## # |
     ######### |
     ###   ### |
     # #   # # |
     ###   ### |
     ######### |
     # ## ## #/
     #########

There are several ways to build them in moltemplate:
1) You can define each cube as a 3x3x3 array of smaller cubes, and then delete
   the 7 interior cubes using the "delete" command.  (Each smaller cube is a
   similar structure containing an array of 3x3x3 even smaller cubes...)
2) You can define each cube as a list of 20 smaller cubes corresponding to the
   cubes that would have remained after deleting the 7 interior cubes.

Method 1 is a little bit simpler, but method 2 is much more efficient because
it never has to create sub-cubes which will be deleted later.

This example uses method 1.

If you are running out of memory, or if moltemplate is taking too long use
method2. It is located in the "memory_efficient_but_ugly_version/" subdirectory.
