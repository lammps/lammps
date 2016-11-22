# This example shows an alternative way to setup the
# aluminum crystal loading simulation described here:
# http://icme.hpc.msstate.edu/mediawiki/index.php/Uniaxial_Compression
# by Mark Tschopp and Nathan R. Rhodes
# For additional backgroumd information, please consult that web page.
#
# In this example, I use moltemplate to build a "DATA" file for this system.
# (I can't think of a compelling reason to do this for simple simulations like
# this. But this approach might be useful if you want to artificially create 
# unusual structures out of aluminum crystals, or mix them with other molecules.
# I created this example in response to a user request.)
#
# Use these commands to generate the LAMMPS input script and data file:

moltemplate.sh system.lt

# This will generate system.data, system.in.init, system.in.settings.
# In addition to will need to download "Al99.eam.alloy" file.
# (It was not included in this directory because if its large size.)
# As of 2012-11, I was able to obtain it here:
# http://www.ctcms.nist.gov/~cbecker/Download/Al-YM/Al99.eam.alloy

