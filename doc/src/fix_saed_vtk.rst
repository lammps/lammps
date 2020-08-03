.. index:: fix saed/vtk

fix saed/vtk command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID saed/vtk Nevery Nrepeat Nfreak c_ID attribute args ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* saed/vtk = style name of this fix command
* Nevery = use input values every this many timesteps
* Nrepeat = # of times to use input values for calculating averages
* Nfreq = calculate averages every this many timesteps
* c_ID = saed compute ID

  .. parsed-literal::

     keyword = *file* or *ave* or *start* or *file* or *overwrite*\ :l
       *ave* args = *one* or *running* or *window M*
         one = output a new average value every Nfreq steps
         running = output cumulative average of all previous Nfreq steps
         window M = output average of M most recent Nfreq steps
       *start* args = Nstart
         Nstart = start averaging on this timestep
       *file* arg = filename
         filename = name of file to output time averages to
       *overwrite* arg = none = overwrite output file with only latest output

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all saed 0.0251 Al O Kmax 1.70 Zone 0 0 1 dR_Ewald 0.01 c 0.5 0.5 0.5
   compute 2 all saed 0.0251 Ni Kmax 1.70 Zone 0 0 0 c 0.05 0.05 0.05 manual echo

   fix 1 all saed/vtk 1 1 1 c_1 file Al2O3_001.saed
   fix 2 all saed/vtk 1 1 1 c_2 file Ni_000.saed

Description
"""""""""""

Time average computed intensities from :doc:`compute saed <compute_saed>` and
write output to a file in the third generation vtk image data format for
visualization directly in parallelized visualization software packages
like ParaView and VisIt. Note that if no time averaging is done, this
command can be used as a convenient way to simply output diffraction
intensities at a single snapshot.

To produce output in the image data vtk format ghost data is added
outside the *Kmax* range assigned in the compute saed. The ghost data is
assigned a value of -1 and can be removed setting a minimum isovolume
of 0 within the visualization software. SAED images can be created by
visualizing a spherical slice of the data that is centered at
R_Ewald\*[h k l]/norm([h k l]), where R_Ewald=1/lambda.

The group specified within this command is ignored. However, note that
specified values may represent calculations performed by saed computes
which store their own "group" definitions.

Fix saed/vtk is designed to work only with :doc:`compute saed <compute_saed>`
values, e.g.

.. code-block:: LAMMPS

   compute 3 top saed 0.0251 Al O
   fix saed/vtk 1 1 1 c_3 file Al2O3_001.saed

----------

The *Nevery*\ , *Nrepeat*\ , and *Nfreq* arguments specify on what
timesteps the input values will be used in order to contribute to the
average.  The final averaged quantities are generated on timesteps
that are a multiple of *Nfreq*\ .  The average is over *Nrepeat*
quantities, computed in the preceding portion of the simulation every
*Nevery* timesteps.  *Nfreq* must be a multiple of *Nevery* and
*Nevery* must be non-zero even if *Nrepeat* is 1.
Also, the timesteps
contributing to the average value cannot overlap,
i.e. Nrepeat\*Nevery can not exceed Nfreq.

For example, if Nevery=2, Nrepeat=6, and Nfreq=100, then values on
timesteps 90,92,94,96,98,100 will be used to compute the final average
on timestep 100.  Similarly for timesteps 190,192,194,196,198,200 on
timestep 200, etc.  If Nrepeat=1 and Nfreq = 100, then no time
averaging is done; values are simply generated on timesteps
100,200,etc.

----------

The output for fix ave/time/saed is a file written with the third generation
vtk image data formatting.  The filename assigned by the *file* keyword is
appended with _N.vtk where N is an index (0,1,2...) to account for multiple
diffraction intensity outputs.

By default the header contains the following information (with example data):

.. parsed-literal::

   # vtk DataFile Version 3.0 c_SAED
   Image data set
   ASCII
   DATASET STRUCTURED_POINTS
   DIMENSIONS 337 219 209
   ASPECT_RATIO 0.00507953 0.00785161 0.00821458
   ORIGIN -0.853361 -0.855826 -0.854316
   POINT_DATA 15424827
   SCALARS intensity float
   LOOKUP_TABLE default
   ...data

In this example, kspace is sampled across a 337 x 219 x 209 point mesh
where the mesh spacing is approximately 0.005, 0.007, and 0.008
inv(length) units in the k1, k2, and k3 directions, respectively.
The data is shifted by -0.85, -0.85, -0.85 inv(length) units so that
the origin will lie at 0, 0, 0.   Here, 15,424,827 kspace points are
sampled in total.

----------

Additional optional keywords also affect the operation of this fix.

The *ave* keyword determines how the values produced every *Nfreq*
steps are averaged with values produced on previous steps that were
multiples of *Nfreq*\ , before they are accessed by another output
command or written to a file.

If the *ave* setting is *one*\ , then the values produced on timesteps
that are multiples of *Nfreq* are independent of each other; they are
output as-is without further averaging.

If the *ave* setting is *running*\ , then the values produced on
timesteps that are multiples of *Nfreq* are summed and averaged in a
cumulative sense before being output.  Each output value is thus the
average of the value produced on that timestep with all preceding
values.  This running average begins when the fix is defined; it can
only be restarted by deleting the fix via the :doc:`unfix <unfix>`
command, or by re-defining the fix by re-specifying it.

If the *ave* setting is *window*\ , then the values produced on
timesteps that are multiples of *Nfreq* are summed and averaged within
a moving "window" of time, so that the last M values are used to
produce the output.  E.g. if M = 3 and Nfreq = 1000, then the output
on step 10000 will be the average of the individual values on steps
8000,9000,10000.  Outputs on early steps will average over less than M
values if they are not available.

The *start* keyword specifies what timestep averaging will begin on.
The default is step 0.  Often input values can be 0.0 at time 0, so
setting *start* to a larger value can avoid including a 0.0 in a
running or windowed average.

The *file* keyword allows a filename to be specified.  Every *Nfreq*
steps, the vector of saed intensity data is written to a new file using
the third generation vtk format.  The base of each file is assigned by
the *file* keyword and this string is appended with _N.vtk where N is
an index (0,1,2...) to account for situations with multiple diffraction
intensity outputs.

The *overwrite* keyword will continuously overwrite the output file
with the latest output, so that it only contains one timestep worth of
output.  This option can only be used with the *ave running* setting.

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The attributes for fix_saed_vtk must match the values assigned in the
associated :doc:`compute_saed <compute_saed>` command.

Related commands
""""""""""""""""

:doc:`compute_saed <compute_saed>`

Default
"""""""

The option defaults are ave = one, start = 0, no file output.

----------

.. _Coleman:

**(Coleman)** Coleman, Spearot, Capolungo, MSMSE, 21, 055020
(2013).
