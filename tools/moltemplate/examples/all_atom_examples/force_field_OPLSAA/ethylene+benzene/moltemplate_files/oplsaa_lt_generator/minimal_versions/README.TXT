MOST USERS SHOULD IGNORE THIS DIRECTORY.

This directory contains versions of the oplsaa_subset.prm file
which nearly all of the OPLSAA force-field information removed.
However for the "ethylene+benzene" example, all of the essential
parameters are contained in these files.  You can use oplsaa_moltemplate.py
with either of these files and the physics should be the same.

However there is no reason to do this.
When you download the "oplsaa.prm" file from:
http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm
(also http://dasher.wustl.edu/ffe/distribution/params/oplsaa.prm)
...just remove the lines beginning with "atom" for atoms you don't need.
You don't have to delete all the other irrelevant interactions.
(In fact, it is hard to do that without making a mistake.
 I recommend that you leave the rest of the oplsaa.prm file alone.)
