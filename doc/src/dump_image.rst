.. index:: dump image

dump image command
==================

dump movie command
==================

Syntax
""""""


.. parsed-literal::

   dump ID group-ID style N file color diameter keyword value ...

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be imaged
* style = *image* or *movie* = style of dump command (other styles *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write image to
* color = atom attribute that determines color of each atom
* diameter = atom attribute that determines size of each atom
* zero or more keyword/value pairs may be appended
* keyword = *atom* or *adiam* or *bond* or *line* or *tri* or *body* or *fix* or *size* or *view* or *center* or *up* or *zoom* or *persp* or *box* or *axes* or *subbox* or *shiny* or *ssao*
  
  .. parsed-literal::
  
       *atom* = yes/no = do or do not draw atoms
       *adiam* size = numeric value for atom diameter (distance units)
       *bond* values = color width = color and width of bonds
         color = *atom* or *type* or *none*
         width = number or *atom* or *type* or *none*
           number = numeric value for bond width (distance units)
       *line* = color width
         color = *type*
         width = numeric value for line width (distance units)
       *tri* = color tflag width
         color = *type*
         tflag = 1 for just triangle, 2 for just tri edges, 3 for both
         width = numeric value for tringle edge width (distance units)
       *body* = color bflag1 bflag2
         color = *type*
         bflag1,bflag2 = 2 numeric flags to affect how bodies are drawn
       *fix* = fixID color fflag1 fflag2
         fixID = ID of fix that generates objects to dray
         color = *type*
         fflag1,fflag2 = 2 numeric flags to affect how fix objects are drawn
       *size* values = width height = size of images
         width = width of image in # of pixels
         height = height of image in # of pixels
       *view* values = theta phi = view of simulation box
         theta = view angle from +z axis (degrees)
         phi = azimuthal view angle (degrees)
         theta or phi can be a variable (see below)
       *center* values = flag Cx Cy Cz = center point of image
         flag = "s" for static, "d" for dynamic
         Cx,Cy,Cz = center point of image as fraction of box dimension (0.5 = center of box)
         Cx,Cy,Cz can be variables (see below)
       *up* values = Ux Uy Uz = direction that is "up" in image
         Ux,Uy,Uz = components of up vector
         Ux,Uy,Uz can be variables (see below)
       *zoom* value = zfactor = size that simulation box appears in image
         zfactor = scale image size by factor > 1 to enlarge, factor < 1 to shrink
         zfactor can be a variable (see below)
       *persp* value = pfactor = amount of "perspective" in image
         pfactor = amount of perspective (0 = none, < 1 = some, > 1 = highly skewed)
         pfactor can be a variable (see below)
       *box* values = yes/no diam = draw outline of simulation box
         yes/no = do or do not draw simulation box lines
         diam = diameter of box lines as fraction of shortest box length
       *axes* values = yes/no length diam = draw xyz axes
         yes/no = do or do not draw xyz axes lines next to simulation box
         length = length of axes lines as fraction of respective box lengths
         diam = diameter of axes lines as fraction of shortest box length
       *subbox* values = yes/no diam = draw outline of processor sub-domains
         yes/no = do or do not draw sub-domain lines
         diam = diameter of sub-domain lines as fraction of shortest box length
       *shiny* value = sfactor = shinyness of spheres and cylinders
         sfactor = shinyness of spheres and cylinders from 0.0 to 1.0
       *ssao* value = yes/no seed dfactor = SSAO depth shading
         yes/no = turn depth shading on/off
         seed = random # seed (positive integer)
         dfactor = strength of shading from 0.0 to 1.0



Examples
""""""""


.. code-block:: LAMMPS

   dump d0 all image 100 dump.\*.jpg type type
   dump d1 mobile image 500 snap.\*.png element element ssao yes 4539 0.6
   dump d2 all image 200 img-\*.ppm type type zoom 2.5 adiam 1.5 size 1280 720
   dump m0 all movie 1000 movie.mpg type type size 640 480
   dump m1 all movie 1000 movie.avi type type size 640 480
   dump m2 all movie 100 movie.m4v type type zoom 1.8 adiam v_value size 1280 720

Description
"""""""""""

Dump a high-quality rendered image of the atom configuration every N
timesteps and save the images either as a sequence of JPEG or PNG or
PPM files, or as a single movie file.  The options for this command as
well as the :doc:`dump_modify <dump_modify>` command control what is
included in the image or movie and how it appears.  A series of such
images can easily be manually converted into an animated movie of your
simulation or the process can be automated without writing the
intermediate files using the dump movie style; see further details
below.  Other dump styles store snapshots of numerical data associated
with atoms in various formats, as discussed on the :doc:`dump <dump>`
doc page.

Note that a set of images or a movie can be made after a simulation
has been run, using the :doc:`rerun <rerun>` command to read snapshots
from an existing dump file, and using these dump commands in the rerun
script to generate the images/movie.

Here are two sample images, rendered as 1024x1024 JPEG files.  Click
to see the full-size images:

.. image:: JPG/dump1_small.jpg
   :target: JPG/dump1.jpg
   :width: 48%
.. image:: JPG/dump2_small.jpg
   :target: JPG/dump2.jpg
   :width: 48%

Only atoms in the specified group are rendered in the image.  The
:doc:`dump_modify region and thresh <dump_modify>` commands can also
alter what atoms are included in the image.
The filename suffix determines whether a JPEG, PNG, or PPM file is
created with the *image* dump style.  If the suffix is ".jpg" or
".jpeg", then a JPEG format file is created, if the suffix is ".png",
then a PNG format is created, else a PPM (aka NETPBM) format file is
created.  The JPEG and PNG files are binary; PPM has a text mode
header followed by binary data. JPEG images have lossy compression;
PNG has lossless compression; and PPM files are uncompressed but can
be compressed with gzip, if LAMMPS has been compiled with
-DLAMMPS\_GZIP and a ".gz" suffix is used.

Similarly, the format of the resulting movie is chosen with the
*movie* dump style. This is handled by the underlying FFmpeg converter
and thus details have to be looked up in the FFmpeg documentation.
Typical examples are: .avi, .mpg, .m4v, .mp4, .mkv, .flv, .mov, .gif
Additional settings of the movie compression like bitrate and
framerate can be set using the :doc:`dump_modify <dump_modify>` command.

To write out JPEG and PNG format files, you must build LAMMPS with
support for the corresponding JPEG or PNG library. To convert images
into movies, LAMMPS has to be compiled with the -DLAMMPS\_FFMPEG
flag. See the :doc:`Build settings <Build_settings>` doc page for
details.

.. note::

   Because periodic boundary conditions are enforced only on
   timesteps when neighbor lists are rebuilt, the coordinates of an atom
   in the image may be slightly outside the simulation box.


----------


Dumps are performed on timesteps that are a multiple of N (including
timestep 0) and on the last timestep of a minimization if the
minimization converges.  Note that this means a dump will not be
performed on the initial timestep after the dump command is invoked,
if the current timestep is not a multiple of N.  This behavior can be
changed via the :doc:`dump_modify first <dump_modify>` command, which
can be useful if the dump command is invoked after a minimization
ended on an arbitrary timestep.  N can be changed between runs by
using the :doc:`dump_modify every <dump_modify>` command.

Dump *image* filenames must contain a wildcard character "\*", so that
one image file per snapshot is written.  The "\*" character is replaced
with the timestep value.  For example, tmp.dump.\*.jpg becomes
tmp.dump.0.jpg, tmp.dump.10000.jpg, tmp.dump.20000.jpg, etc.  Note
that the :doc:`dump_modify pad <dump_modify>` command can be used to
insure all timestep numbers are the same length (e.g. 00010), which
can make it easier to convert a series of images into a movie in the
correct ordering.

Dump *movie* filenames on the other hand, must not have any wildcard
character since only one file combining all images into a single
movie will be written by the movie encoder.


----------


The *color* and *diameter* settings determine the color and size of
atoms rendered in the image.  They can be any atom attribute defined
for the :doc:`dump custom <dump>` command, including *type* and
*element*\ .  This includes per-atom quantities calculated by a
:doc:`compute <compute>`, :doc:`fix <fix>`, or :doc:`variable <variable>`,
which are prefixed by "c\_", "f\_", or "v\_" respectively.  Note that the
*diameter* setting can be overridden with a numeric value applied to
all atoms by the optional *adiam* keyword.

If *type* is specified for the *color* setting, then the color of each
atom is determined by its atom type.  By default the mapping of types
to colors is as follows:

* type 1 = red
* type 2 = green
* type 3 = blue
* type 4 = yellow
* type 5 = aqua
* type 6 = cyan

and repeats itself for types > 6.  This mapping can be changed by the
:doc:`dump_modify acolor <dump_modify>` command.

If *type* is specified for the *diameter* setting then the diameter of
each atom is determined by its atom type.  By default all types have
diameter 1.0.  This mapping can be changed by the :doc:`dump_modify adiam <dump_modify>` command.

If *element* is specified for the *color* and/or *diameter* setting,
then the color and/or diameter of each atom is determined by which
element it is, which in turn is specified by the element-to-type
mapping specified by the "dump\_modify element" command.  By default
every atom type is C (carbon).  Every element has a color and diameter
associated with it, which is the same as the colors and sizes used by
the `AtomEye <atomeye_>`_ visualization package.

.. _atomeye: http://mt.seas.upenn.edu/Archive/Graphics/A

If other atom attributes are used for the *color* or *diameter*
settings, they are interpreted in the following way.

If "vx", for example, is used as the *color* setting, then the color
of the atom will depend on the x-component of its velocity.  The
association of a per-atom value with a specific color is determined by
a "color map", which can be specified via the
:doc:`dump_modify <dump_modify>` command.  The basic idea is that the
atom-attribute will be within a range of values, and every value
within the range is mapped to a specific color.  Depending on how the
color map is defined, that mapping can take place via interpolation so
that a value of -3.2 is halfway between "red" and "blue", or
discretely so that the value of -3.2 is "orange".

If "vx", for example, is used as the *diameter* setting, then the atom
will be rendered using the x-component of its velocity as the
diameter.  If the per-atom value <= 0.0, them the atom will not be
drawn.  Note that finite-size spherical particles, as defined by
:doc:`atom_style sphere <atom_style>` define a per-particle radius or
diameter, which can be used as the *diameter* setting.

----------

The various keywords listed above control how the image is rendered.
As listed below, all of the keywords have defaults, most of which you
will likely not need to change.  The :doc:`dump modify <dump_modify>`
also has options specific to the dump image style, particularly for
assigning colors to atoms, bonds, and other image features.

----------

The *atom* keyword allow you to turn off the drawing of all atoms, if
the specified value is *no*\ .  Note that this will not turn off the
drawing of particles that are represented as lines, triangles, or
bodies, as discussed below.  These particles can be drawn separately
if the *line*\ , *tri*\ , or *body* keywords are used.

The *adiam* keyword allows you to override the *diameter* setting to
set a single numeric *size*\ .  All atoms will be drawn with that
diameter, e.g. 1.5, which is in whatever distance :doc:`units <units>`
the input script defines, e.g. Angstroms.

----------

The *bond* keyword allows to you to alter how bonds are drawn.  A bond
is only drawn if both atoms in the bond are being drawn due to being
in the specified group and due to other selection criteria
(e.g. region, threshold settings of the
:doc:`dump_modify <dump_modify>` command).  By default, bonds are drawn
if they are defined in the input data file as read by the
:doc:`read_data <read_data>` command.  Using *none* for both the bond
*color* and *width* value will turn off the drawing of all bonds.

If *atom* is specified for the bond *color* value, then each bond is
drawn in 2 halves, with the color of each half being the color of the
atom at that end of the bond.

If *type* is specified for the *color* value, then the color of each
bond is determined by its bond type.  By default the mapping of bond
types to colors is as follows:

* type 1 = red
* type 2 = green
* type 3 = blue
* type 4 = yellow
* type 5 = aqua
* type 6 = cyan

and repeats itself for bond types > 6.  This mapping can be changed by
the :doc:`dump_modify bcolor <dump_modify>` command.

The bond *width* value can be a numeric value or *atom* or *type* (or
*none* as indicated above).

If a numeric value is specified, then all bonds will be drawn as
cylinders with that diameter, e.g. 1.0, which is in whatever distance
:doc:`units <units>` the input script defines, e.g. Angstroms.

If *atom* is specified for the *width* value, then each bond
will be drawn with a width corresponding to the minimum diameter
of the 2 atoms in the bond.

If *type* is specified for the *width* value then the diameter of each
bond is determined by its bond type.  By default all types have
diameter 0.5.  This mapping can be changed by the :doc:`dump_modify bdiam <dump_modify>` command.

----------

The *line* keyword can be used when :doc:`atom_style line <atom_style>`
is used to define particles as line segments, and will draw them as
lines.  If this keyword is not used, such particles will be drawn as
spheres, the same as if they were regular atoms.  The only setting
currently allowed for the *color* value is *type*\ , which will color
the lines according to the atom type of the particle.  By default the
mapping of types to colors is as follows:

* type 1 = red
* type 2 = green
* type 3 = blue
* type 4 = yellow
* type 5 = aqua
* type 6 = cyan

and repeats itself for types > 6.  There is not yet an option to
change this via the :doc:`dump_modify <dump_modify>` command.

The line *width* can only be a numeric value, which specifies that all
lines will be drawn as cylinders with that diameter, e.g. 1.0, which
is in whatever distance :doc:`units <units>` the input script defines,
e.g. Angstroms.


----------


The *tri* keyword can be used when :doc:`atom_style tri <atom_style>` is
used to define particles as triangles, and will draw them as triangles
or edges (3 lines) or both, depending on the setting for *tflag*\ .  If
edges are drawn, the *width* setting determines the diameters of the
line segments.  If this keyword is not used, triangle particles will
be drawn as spheres, the same as if they were regular atoms.  The only
setting currently allowed for the *color* value is *type*\ , which will
color the triangles according to the atom type of the particle.  By
default the mapping of types to colors is as follows:

* type 1 = red
* type 2 = green
* type 3 = blue
* type 4 = yellow
* type 5 = aqua
* type 6 = cyan

and repeats itself for types > 6.  There is not yet an option to
change this via the :doc:`dump_modify <dump_modify>` command.


----------


The *body* keyword can be used when :doc:`atom_style body <atom_style>`
is used to define body particles with internal state
(e.g. sub-particles), and will drawn them in a manner specific to the
body style.  If this keyword is not used, such particles will be drawn
as spheres, the same as if they were regular atoms.

The :doc:`Howto body <Howto_body>` doc page describes the body styles
LAMMPS currently supports, and provides more details as to the kind of
body particles they represent and how they are drawn by this dump
image command.  For all the body styles, individual atoms can be
either a body particle or a usual point (non-body) particle.  Non-body
particles will be drawn the same way they would be as a regular atom.
The *bflag1* and *bflag2* settings are numerical values which are
passed to the body style to affect how the drawing of a body particle
is done.  See the :doc:`Howto body <Howto_body>` doc page for a
description of what these parameters mean for each body style.

The only setting currently allowed for the *color* value is *type*\ ,
which will color the body particles according to the atom type of the
particle.  By default the mapping of types to colors is as follows:

* type 1 = red
* type 2 = green
* type 3 = blue
* type 4 = yellow
* type 5 = aqua
* type 6 = cyan

and repeats itself for types > 6.  There is not yet an option to
change this via the :doc:`dump_modify <dump_modify>` command.


----------


The *fix* keyword can be used with a :doc:`fix <fix>` that produces
objects to be drawn.

The *fflag1* and *fflag2* settings are numerical values which are
passed to the fix to affect how the drawing of its objects is done.
See the individual fix doc page for a description of what these
parameters mean for a particular fix.

The only setting currently allowed for the *color* value is *type*\ ,
which will color the fix objects according to their type.  By default
the mapping of types to colors is as follows:

* type 1 = red
* type 2 = green
* type 3 = blue
* type 4 = yellow
* type 5 = aqua
* type 6 = cyan

and repeats itself for types > 6.  There is not yet an option to
change this via the :doc:`dump_modify <dump_modify>` command.


----------


The *size* keyword sets the width and height of the created images,
i.e. the number of pixels in each direction.


----------


The *view*\ , *center*\ , *up*\ , *zoom*\ , and *persp* values determine how
3d simulation space is mapped to the 2d plane of the image.  Basically
they control how the simulation box appears in the image.

All of the *view*\ , *center*\ , *up*\ , *zoom*\ , and *persp* values can be
specified as numeric quantities, whose meaning is explained below.
Any of them can also be specified as an :doc:`equal-style variable <variable>`, by using v\_name as the value, where "name" is
the variable name.  In this case the variable will be evaluated on the
timestep each image is created to create a new value.  If the
equal-style variable is time-dependent, this is a means of changing
the way the simulation box appears from image to image, effectively
doing a pan or fly-by view of your simulation.

The *view* keyword determines the viewpoint from which the simulation
box is viewed, looking towards the *center* point.  The *theta* value
is the vertical angle from the +z axis, and must be an angle from 0 to
180 degrees.  The *phi* value is an azimuthal angle around the z axis
and can be positive or negative.  A value of 0.0 is a view along the
+x axis, towards the *center* point.  If *theta* or *phi* are
specified via variables, then the variable values should be in
degrees.

The *center* keyword determines the point in simulation space that
will be at the center of the image.  *Cx*\ , *Cy*\ , and *Cz* are
specified as fractions of the box dimensions, so that (0.5,0.5,0.5) is
the center of the simulation box.  These values do not have to be
between 0.0 and 1.0, if you want the simulation box to be offset from
the center of the image.  Note, however, that if you choose strange
values for *Cx*\ , *Cy*\ , or *Cz* you may get a blank image.  Internally,
*Cx*\ , *Cy*\ , and *Cz* are converted into a point in simulation space.
If *flag* is set to "s" for static, then this conversion is done once,
at the time the dump command is issued.  If *flag* is set to "d" for
dynamic then the conversion is performed every time a new image is
created.  If the box size or shape is changing, this will adjust the
center point in simulation space.

The *up* keyword determines what direction in simulation space will be
"up" in the image.  Internally it is stored as a vector that is in the
plane perpendicular to the view vector implied by the *theta* and
*pni* values, and which is also in the plane defined by the view
vector and user-specified up vector.  Thus this internal vector is
computed from the user-specified *up* vector as


.. parsed-literal::

   up_internal = view cross (up cross view)

This means the only restriction on the specified *up* vector is that
it cannot be parallel to the *view* vector, implied by the *theta* and
*phi* values.

The *zoom* keyword scales the size of the simulation box as it appears
in the image.  The default *zfactor* value of 1 should display an
image mostly filled by the atoms in the simulation box.  A *zfactor* >
1 will make the simulation box larger; a *zfactor* < 1 will make it
smaller.  *Zfactor* must be a value > 0.0.

The *persp* keyword determines how much depth perspective is present
in the image.  Depth perspective makes lines that are parallel in
simulation space appear non-parallel in the image.  A *pfactor* value
of 0.0 means that parallel lines will meet at infinity (1.0/pfactor),
which is an orthographic rendering with no perspective.  A *pfactor*
value between 0.0 and 1.0 will introduce more perspective.  A *pfactor*
value > 1 will create a highly skewed image with a large amount of
perspective.

.. note::

   The *persp* keyword is not yet supported as an option.


----------


The *box* keyword determines if and how the simulation box boundaries
are rendered as thin cylinders in the image.  If *no* is set, then the
box boundaries are not drawn and the *diam* setting is ignored.  If
*yes* is set, the 12 edges of the box are drawn, with a diameter that
is a fraction of the shortest box length in x,y,z (for 3d) or x,y (for
2d).  The color of the box boundaries can be set with the :doc:`dump_modify boxcolor <dump_modify>` command.

The *axes* keyword determines if and how the coordinate axes are
rendered as thin cylinders in the image.  If *no* is set, then the
axes are not drawn and the *length* and *diam* settings are ignored.
If *yes* is set, 3 thin cylinders are drawn to represent the x,y,z
axes in colors red,green,blue.  The origin of these cylinders will be
offset from the lower left corner of the box by 10%.  The *length*
setting determines how long the cylinders will be as a fraction of the
respective box lengths.  The *diam* setting determines their thickness
as a fraction of the shortest box length in x,y,z (for 3d) or x,y (for
2d).

The *subbox* keyword determines if and how processor sub-domain
boundaries are rendered as thin cylinders in the image.  If *no* is
set (default), then the sub-domain boundaries are not drawn and the
*diam* setting is ignored.  If *yes* is set, the 12 edges of each
processor sub-domain are drawn, with a diameter that is a fraction of
the shortest box length in x,y,z (for 3d) or x,y (for 2d).  The color
of the sub-domain boundaries can be set with the :doc:`dump_modify boxcolor <dump_modify>` command.


----------


The *shiny* keyword determines how shiny the objects rendered in the
image will appear.  The *sfactor* value must be a value 0.0 <=
*sfactor* <= 1.0, where *sfactor* = 1 is a highly reflective surface
and *sfactor* = 0 is a rough non-shiny surface.

The *ssao* keyword turns on/off a screen space ambient occlusion
(SSAO) model for depth shading.  If *yes* is set, then atoms further
away from the viewer are darkened via a randomized process, which is
perceived as depth.  The calculation of this effect can increase the
cost of computing the image by roughly 2x.  The strength of the effect
can be scaled by the *dfactor* parameter.  If *no* is set, no depth
shading is performed.


----------


A series of JPEG, PNG, or PPM images can be converted into a movie
file and then played as a movie using commonly available tools. Using
dump style *movie* automates this step and avoids the intermediate
step of writing (many) image snapshot file. But LAMMPS has to be
compiled with -DLAMMPS\_FFMPEG and an FFmpeg executable have to be
installed.

To manually convert JPEG, PNG or PPM files into an animated GIF or
MPEG or other movie file you can use:

* a) Use the ImageMagick convert program.
  
  .. code-block:: bash
  
     % convert *.jpg foo.gif
     % convert -loop 1 *.ppm foo.mpg


  Animated GIF files from ImageMagick are not optimized. You can use
  a program like gifsicle to optimize and thus massively shrink them.
  MPEG files created by ImageMagick are in MPEG-1 format with a rather
  inefficient compression and low quality compared to more modern
  compression styles like MPEG-4, H.264, VP8, VP9, H.265 and so on.

* b) Use QuickTime.

  Select "Open Image Sequence" under the File menu Load the images into
  QuickTime to animate them Select "Export" under the File menu Save the
  movie as a QuickTime movie (\*.mov) or in another format.  QuickTime
  can generate very high quality and efficiently compressed movie
  files. Some of the supported formats require to buy a license and some
  are not readable on all platforms until specific runtime libraries are
  installed.

* c) Use FFmpeg

  FFmpeg is a command line tool that is available on many platforms and
  allows extremely flexible encoding and decoding of movies.

  
  .. code-block:: bash
  
     cat snap.*.jpg | ffmpeg -y -f image2pipe -c:v mjpeg -i - -b:v 2000k movie.m4v
     cat snap.*.ppm | ffmpeg -y -f image2pipe -c:v ppm -i - -b:v 2400k movie.avi


  Front ends for FFmpeg exist for multiple platforms. For more
  information see the `FFmpeg homepage <http://www.ffmpeg.org/>`_

----------


Play the movie:

* a) Use your browser to view an animated GIF movie.

  Select "Open File" under the File menu
  Load the animated GIF file

* b) Use the freely available mplayer or ffplay tool to view a
  movie. Both are available for multiple OSes and support a large
  variety of file formats and decoders.
  
  .. code-block:: bash
  
     % mplayer foo.mpg
     % ffplay bar.avi

* c) Use the `Pizza.py <http://www.sandia.gov/~sjplimp/pizza.html>`_
  `animate tool <http://www.sandia.gov/~sjplimp/pizza/doc/animate.html>`_,
  which works directly on a series of image files.
  
  .. code-block:: python
  
     a = animate("foo*.jpg")

* d) QuickTime and other Windows- or MacOS-based media players can
  obviously play movie files directly. Similarly for corresponding tools
  bundled with Linux desktop environments.  However, due to licensing
  issues with some file formats, the formats may require installing
  additional libraries, purchasing a license, or may not be
  supported.

----------

See the :doc:`Modify <Modify>` doc page for information on how to add
new compute and fix styles to LAMMPS to calculate per-atom quantities
which could then be output into dump files.


----------


Restrictions
""""""""""""


To write JPEG images, you must use the -DLAMMPS\_JPEG switch when
building LAMMPS and link with a JPEG library. To write PNG images, you
must use the -DLAMMPS\_PNG switch when building LAMMPS and link with a
PNG library.

To write *movie* dumps, you must use the -DLAMMPS\_FFMPEG switch when
building LAMMPS and have the FFmpeg executable available on the
machine where LAMMPS is being run.  Typically it's name is lowercase,
i.e. ffmpeg.

See the :doc:`Build settings <Build_settings>` doc page for details.

Note that since FFmpeg is run as an external program via a pipe,
LAMMPS has limited control over its execution and no knowledge about
errors and warnings printed by it. Those warnings and error messages
will be printed to the screen only. Due to the way image data is
communicated to FFmpeg, it will often print the message


.. parsed-literal::

   pipe:: Input/output error

which can be safely ignored. Other warnings
and errors have to be addressed according to the FFmpeg documentation.
One known issue is that certain movie file formats (e.g. MPEG level 1
and 2 format streams) have video bandwidth limits that can be crossed
when rendering too large of image sizes. Typical warnings look like
this:


.. parsed-literal::

   [mpeg @ 0x98b5e0] packet too large, ignoring buffer limits to mux it
   [mpeg @ 0x98b5e0] buffer underflow st=0 bufi=281407 size=285018
   [mpeg @ 0x98b5e0] buffer underflow st=0 bufi=283448 size=285018

In this case it is recommended to either reduce the size of the image
or encode in a different format that is also supported by your copy of
FFmpeg, and which does not have this limitation (e.g. .avi, .mkv,
mp4).

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`

Default
"""""""

The defaults for the keywords are as follows:

* adiam = not specified (use diameter setting)
* atom = yes
* bond = none none (if no bonds in system)
* bond = atom 0.5 (if bonds in system)
* size = 512 512
* view = 60 30 (for 3d)
* view = 0 0 (for 2d)
* center = s 0.5 0.5 0.5
* up = 0 0 1 (for 3d)
* up = 0 1 0 (for 2d)
* zoom = 1.0
* persp = 0.0
* box = yes 0.02
* axes = no 0.0 0.0
* subbox no 0.0
* shiny = 1.0
* ssao = no
