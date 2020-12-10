The input script in2d.selfpropel illustrates how to add an active force
of magnitude fp=400.0 to each particle in a 2D geometry.
There is no time-stepping involved here, so one needs to add their
own integrating fix. The contribution to the pressure is
calculated as fp*<mu dot r>/(2*A), and will return the ideal
active swim pressure in a system of no interactions (with
e.g. overdamped brownian dynamics).


The input script in3d_no_active_pressure.selfpropel
illustrates how to add an active force
of magnitude fp=10.0 to each particle in a 3D geometry.
There is no time-stepping involved here, so one needs to add their
own integrating fix. The pressure contribution from the
active force is turned off using fix_modify.

To confirm that these fixes are working correctly, just add
a dump computation in the input scripts and directly calculate
the forces that should be applied (fx = mux*fp, fy = muy*fp,
fz = muz*fp), and compare them to the forces calculated via
the simulation.
