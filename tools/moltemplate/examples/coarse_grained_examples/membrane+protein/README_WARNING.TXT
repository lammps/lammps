WARNING:
   This is not a realistic simulation of a protein in a lipid membrane.

   When I released this example in 2013, the simulations ran without crashing and produced no obvious visual flaws.  However recently (2015) I've noticed that there are some serious problems:  Currently, the protein no-longer remains perpendicular to the membrane during the simulations.  (In some cases, it also unfolds.)  The original version of this example did not behave this way.  I'm not sure what changed since then.

   This example was created to demonstrate how to combine totally different kinds of coarse-grained molecules (with different kinds of force-fields) together in the same simulation in LAMMPS.  Tuning the force-field parameters to get realistic results was not the goal.

   One day I might tinker with the shape of the protein and the force-field parameters to try and avoid this unsightly outcome.  (Probably shortening the length of the helices would help.)  Either way, I doubt it would make the simulations any more accurate.

   I apologize that the simulations no longer produce the same pretty images shown on the web page.

Andrew 2015-8-02
