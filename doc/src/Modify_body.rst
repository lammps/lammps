Body styles
===========

Classes that define body particles are derived from the Body class.
Body particles can represent complex entities, such as surface meshes
of discrete points, collections of sub-particles, deformable objects,
etc.

See the :doc:`Howto body <Howto_body>` page for an overview of using
body particles and the various body styles LAMMPS supports.  New
styles can be created to add new kinds of body particles to LAMMPS.

Body_nparticle.cpp is an example of a body particle that is treated as
a rigid body containing N sub-particles.

Here is a brief description of methods you define in your new derived
class.  See ``src/body.h`` for details.

+----------------------+-----------------------------------------------------------+
| Required             | "pure" methods that *must* be overridden                  |
+======================+===========================================================+
| data_body            | process a line from the Bodies section of a data file     |
+----------------------+-----------------------------------------------------------+
| pack_data_body       | pack body data into a buffer for writing to a data file   |
+----------------------+-----------------------------------------------------------+
| write_data_body      | write body data to a data file                            |
+----------------------+-----------------------------------------------------------+
| noutrow              | number of sub-particles output is generated for           |
+----------------------+-----------------------------------------------------------+
| noutcol              | number of values per-sub-particle output is generated for |
+----------------------+-----------------------------------------------------------+
| output               | output values for the Mth sub-particle                    |
+----------------------+-----------------------------------------------------------+
| image                | provide geometry for dump image rendering                 |
+----------------------+-----------------------------------------------------------+

+----------------------+-----------------------------------------------------------+
| Optional             | methods that have a default or empty implementation       |
+======================+===========================================================+
| pack_comm_body       | body attributes to communicate every timestep             |
+----------------------+-----------------------------------------------------------+
| unpack_comm_body     | unpacking of those attributes                             |
+----------------------+-----------------------------------------------------------+
| pack_border_body     | body attributes to communicate when reneighboring is done |
+----------------------+-----------------------------------------------------------+
| unpack_border_body   | unpacking of those attributes                             |
+----------------------+-----------------------------------------------------------+
| radius_body          | return effective radius for use by granular pair styles   |
+----------------------+-----------------------------------------------------------+
