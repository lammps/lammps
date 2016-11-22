 WARNING:  A reader has emailed me to point out:

    THIS IS NOT A REALISTIC MODEL OF A GRAPHENE-NANOTUBE JUNCTION.
    A real junction would likely be curved near the boundary, 
    not a 90 degree junction.  (Although both graphene and nanotubes
    consist of hexagons of carbon atoms, you would need 6 heptagons
    near the junction between the nanotube and the graphene 
    to account for the negative Gaussian curvature there).

 To solve this problem:
    Moltemplate allows you to move, add, customize or delete individual
    atoms near the boundary.  You can move atoms by overwriting their 
    coordinates using additional write("Data Atoms") statements (after
    the walls and tube are created).  You can also adjust their partial charge.

 Alternately, you could start with the structure provided here, add or delete
    atoms if necessary, and relax/minimize the coordinates of the carbon 
    atoms using LAMMPS.  You could also run a high temperature annealing
    simulation to relax their positions.  If it helps, the AIREBO 
    force-field has used in LAMMPS to simulate carbon nanotubes breaking:
     http://scitation.aip.org/content/aip/journal/jcp/134/20/10.1063/1.3594197
     http://lammps.sandia.gov/pictures.html#cnt
