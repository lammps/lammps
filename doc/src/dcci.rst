.. index:: dcci

dcci command
===================

Syntax
""""""

.. parsed-literal::

   dcci Tcoex Pcoex lambda fix-ID barostat-ID t_sc keyword values

* Tcoex = intial temperature of coexistence for the two phases 
* Pcoex = intial pressure of coexistence for the two phases 
* lambda = initial scaling variable of this command
* fix-ID = ID of *fix adapt/dcci*
* barostat-ID = ID of the barostat that controls the pressure
* t_sc = number of step for the dynamical Clausius-Clapeyron integration
* keyword = *temp* or *press*
  
  .. parsed-literal::

       *temp* values = Tstart Tstop =  desired temperature of coexistence at start/end of run (temperature units)
       *press* values = Pstart Pstop =  desired pressure of coexistence at start/end of run (pressure units) 

Examples
""""""""

.. code-block:: LAMMPS

   dcci 0.7 1.0 ${lambda} fxADAPT-DCCI fxnpt 50000 press 1.0 174.0
   dcci 0.7 1.0 ${lambda} fxADAPT-DCCI fxnpt 50000 temp 0.7 7.8
   dcci ${Tcoex} ${Pcoex} ${lambda} fxADAPTDCCI fxnph ${t_sc} press ${Pi} ${Pf}
   dcci ${Tcoex} ${Pcoex} ${lambda} fxADAPTDCCI fxnph ${t_sc} temp ${Ti} ${Tf}   

Description
"""""""""""

This command allows you to compute the phase boundaries, knowing only a single
coexistence point between phases, by performing a nonequilibrium dynamical
Clausius-Clapeyron integration of both phases 
simultaneously :ref:`(deKoning) <deKoning>`. A detailed explanation of how to use this
command and choose its parameters for optimal performance and accuracy is
given in the paper by :ref:`(Cajahuaringa) <Cajahuaringa>`. This is done
by run only two systems in parallel (one for each phase). See the :doc:`Howto replica <Howto_replica>` doc for details of how to run replicas on multiple partitions of one or more processors.

Instead of variables mention above, this command need scaling the forces in
both systems in parallel by the :doc:`fix adapt/dcci <fix_adapt_dcci>`. The scaling parameter *lambda* is controlled by the dcci command, which is typically a variable previously set in the input script, so that each system is assigned the same scaling parameter and scaling pressure, see the paper :ref:`(Cajahuaringa) <Cajahuaringa>` for more details. For example:

.. code-block:: LAMMPS

   variable phase world solid liquid
   ...
   read_data ${phase}/phase.lammps
   ...
   variable lambda equal 1.0
   fix fxnpt all npt temp ${Tcoex} ${Tcoex} ${Tdamp} iso ${Pcoex} ${Pcoex} ${Pdamp}
   fix fxADAPT-DCCI all adapt/dcci ${lambda} pair lj/cut fscale * *
   dcci ${Tcoex} ${Pcoex} ${lambda} fxADAPT-DCCI fxnpt ${t_sc} press ${Pi} ${Pf}

in the similar way of :doc:`temper <temper>` command is necessary especify the path of each phase, after that it is define a initial *lambda* value, who is pass to the :doc:`fix adapt/dcci <fix_adapt_dcci>` for scaling the forces in each phase and the dcci command calculates the new value of the *lambda* variable and scaling pressure to perform the simulation of both phase in the next coexistence condition, this process is repeated *t_sc* steps until achevies final temperature (pressure).

The last argument *temp* or *press* keyword allows use two different control to computed the phase bounadaries. If is select the *temp* keyword is calculated the coxistence pressure as function of the temperature from *Tstart* to *Tstop*, otherwise it is select the *press* keyword is calculated the coxistence temperature as function of the pressure from *Pstart* to *Pstop*.

As a dcci run proceeds, two log files and screen output files are created, one per phase. By default these files are named log.lammps.M and screen.M where *M* is the phase number from 0 to 1. See the :doc:`-log and -screen command-line swiches <Run_options>` for info on how to change these names.

The main screen and log file (log.lammps) will list information about dcci dynamics assigned to each phase at each coxistence condition output timestep. E.g. for a simulation of two phases:

.. parsed-literal::

   Running on 2 partitions of processors
   Step Tcoex  Pcoex  lambda Pcoex_rs pe1  pe2  vol1  vol2
   50000    0.770000    1.000000    1.00000000    1.00000000    -3657.45830    -3049.38969    508.9074    576.0835
   50001    0.770265    1.003460    0.99965590    1.00311471    -3658.25936    -3048.49558    508.9249    576.0910
   50002    0.770529    1.006920    0.99931295    1.00622820    -3658.60164    -3048.24218    508.9424    576.0989
   50003    0.770793    1.010380    0.99897058    1.00933989    -3658.79853    -3048.69379    508.9599    576.1074
   50004    0.771058    1.013840    0.99862835    1.01244937    -3658.54627    -3049.93109    508.9774    576.1165
   50005    0.771322    1.017300    0.99828565    1.01555599    -3658.06732    -3051.29644    508.9948    576.1260
   50006    0.771588    1.020760    0.99794228    1.01865956    -3657.88431    -3052.18231    509.0121    576.1360
   50007    0.771854    1.024220    0.99759863    1.02176047    -3657.73010    -3053.17874    509.0294    576.1465
   ...

The second and third columns rotuled as Tcoex and Pcoex, are the temperature and the pressure coexistences conditions calculated by the dcci method, the fourth and fifth columns rotuled as lambda and Pcoex_rs, are the lambda scaling parameter and the scaling pressure, the pe1 and pe2 are the instantaneus potential energy of each phase and the vol1 and vol2 are the instantaneus volume each phase.

If you need calculated the phase boundaries between isotropic and anisotropic phases, for example: the phase boundaries between cubic diamond and beta-tin phases of silicon, this can be done using the follow protocol:

.. code-block:: LAMMPS

   variable phase world Si-I Si-II
   variable ensemble world iso aniso
   ...
   read_data ${phase}/phase.lammps
   ...
   variable lambda equal 1.0
   fix f1 all npt temp ${Tcoex} ${Tcoex} ${Tdamp} ${ensemble} ${Pcoex} ${Pcoex} ${Pdamp}
   fix f2 all adapt/dcci ${lambda} pair sw fscale * *
   dcci ${Tcoex} ${Pcoex} ${lambda} f1 f2 ${t_sc} temp ${Pi} ${Pf}
   
   
.. note::

   As described in :ref:`(Cajahuaringa) <Cajahuaringa>`, for the calculation of the phase  
   boundaries with small slope (in absolute value), which is common in solid-solid coexistence 
   curves, it is more convenient control of the coexistence temperature and if, on other hand, the 
   phase boundaries is expected to have a large slope (in absolute value),which is common in 
   solid-fluid  and fluid-fluid coexistence curves, for these case is more convenient control the  
   coexistence pressure.

----------

Restrictions
""""""""""""

This command must be used with :doc:`fix adapt/dcci <fix_adapt_dcci>`.

This command only works with atomistic systems.

This command should be used with a fix that maintains the isothermal-isobaric (NPT) ensemble.

This command can only be used if LAMMPS was built with the REPLICA
package.  See the :doc:`Build package <Build_package>` doc
page for more info.


Related commands
""""""""""""""""

:doc:`fix adapt/dcci <fix_adapt_dcci>`, :doc:`temper <temper>`, :doc:`prd <prd>`

Default
"""""""

none

----------

.. _Cajahuaringa:

**(Cajahuaringa)** Cajahuaringa and Antonelli, Comput. Mater. Sci., 111275, 207 (2022).

.. _deKoning:

**(deKoning)** de Koning, Antonelli and Sidney, J Chem Phys 115, 11025 (2001).
