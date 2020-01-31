.. index:: fix edpd/source

fix edpd/source command
=======================

fix tdpd/source command
=======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID edpd/source keyword values ...
   fix ID group-ID tdpd/source cc_index keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* edpd/source or tdpd/source = style name of this fix command
* index (only specified for tdpd/source) = index of chemical species (1 to Nspecies)
* keyword = *sphere* or *cuboid*
  
  .. parsed-literal::
  
       *sphere* values = cx,cy,cz,radius,source
         cx,cy,cz = x,y,z center of spherical domain (distance units)
         radius = radius of a spherical domain (distance units)
         source = heat source or concentration source (flux units, see below)
       *cuboid* values = cx,cy,cz,dLx,dLy,dLz,source
         cx,cy,cz = x,y,z lower left corner of a cuboid domain (distance units)
         dLx,dLy,dLz = x,y,z side length of a cuboid domain (distance units)
         source = heat source or concentration source (flux units, see below)



Examples
""""""""


.. parsed-literal::

   fix 1 all edpd/source sphere 0.0 0.0 0.0 5.0 0.01
   fix 1 all edpd/source cuboid 0.0 0.0 0.0 20.0 10.0 10.0 -0.01
   fix 1 all tdpd/source 1 sphere 5.0 0.0 0.0 5.0 0.01
   fix 1 all tdpd/source 2 cuboid 0.0 0.0 0.0 20.0 10.0 10.0 0.01

Description
"""""""""""

Fix *edpd/source* adds a heat source as an external heat flux to each
atom in a spherical or cuboid domain, where the *source* is in units
of energy/time.  Fix *tdpd/source* adds an external concentration
source of the chemical species specified by *index* as an external
concentration flux for each atom in a spherical or cuboid domain,
where the *source* is in units of mole/volume/time.

This command can be used to give an additional heat/concentration
source term to atoms in a simulation, such as for a simulation of a
heat conduction with a source term (see Fig.12 in :ref:`(Li2014) <Li2014b>`)
or diffusion with a source term (see Fig.1 in :ref:`(Li2015) <Li2015b>`), as
an analog of a periodic Poiseuille flow problem.

If the *sphere* keyword is used, the *cx,cy,cz,radius* defines a
spherical domain to apply the source flux to.

If the *cuboid* keyword is used, the *cx,cy,cz,dLx,dLy,dLz* defines a
cuboid domain to apply the source flux to.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the USER-MESO package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

Fix *edpd/source* must be used with the :doc:`pair_style edpd <pair_meso>` command.  Fix *tdpd/source* must be used with the
:doc:`pair_style tdpd <pair_meso>` command.

Related commands
""""""""""""""""

:doc:`pair_style edpd <pair_meso>`, :doc:`pair_style tdpd <pair_meso>`,
:doc:`compute edpd/temp/atom <compute_edpd_temp_atom>`, :doc:`compute tdpd/cc/atom <compute_tdpd_cc_atom>`

**Default:** none


----------


.. _Li2014b:



**(Li2014)** Z. Li, Y.-H. Tang, H. Lei, B. Caswell and G.E. Karniadakis,
"Energy-conserving dissipative particle dynamics with
temperature-dependent properties", J. Comput. Phys., 265: 113-127
(2014). DOI: 10.1016/j.jcp.2014.02.003

.. _Li2015b:



**(Li2015)** Z. Li, A. Yazdani, A. Tartakovsky and G.E. Karniadakis,
"Transport dissipative particle dynamics model for mesoscopic
advection-diffusion-reaction problems", J. Chem. Phys., 143: 014101
(2015).  DOI: 10.1063/1.4923254


