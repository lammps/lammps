 ------- A note on building the graphene sheet in VMD: ------

Probably you can ignore these instructions.
These instructions are not necessary for this example to run.

This example contains several pyramid shaped objects resting on a surface
made of graphene.  The instructions in this file explain how to build the
graphene (representing the "ground") using VMD instead of with moltemplate.
  Why do this?
VMD can create graphene sheets with bonds connecting neighboring carbon atoms,
(which looks more pretty).  However, as of 2013-4-29, moltemplate currently
can not generate these bonds.  It does not matter physically in this case,
because the graphene sheet used here does not move.  It is only used as
scenery, to graphically represent the ground surface.

Select "Extensions"->"Modeling"->"Carbon Nanotube Builder"
     Build a graphene sheet of size 39.8 x 39.8 (units: nm)
     400.3358398  399.876008
     (try to use a size compatible with the periodic boundaries)
Select "Extensions"->"Tk Console", and type
     display backgroundgradient on

Note: If you want to do this, before you run moltemplate, you may want to delete
      the sections of the "system.lt" file (located in "moltemplate_files")
      which define the graphene wall.  Instead create the graphene data file
      in VMD.  You will have to manually merge the data file for graphene
      with the data file for the pyramids created by moltemplate,
      (taking care to avoid overlapping atom-id numbers).
