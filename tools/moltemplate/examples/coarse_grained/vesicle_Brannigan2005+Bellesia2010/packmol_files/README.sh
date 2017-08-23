# Here we generate the starting coordinates of the simulation
# using PACKMOL.



# You must run each packmol commend one after the other.
# NOTE: If PACKMOL gets stuck in an endless loop, then edit the corresponding
# "inp" file.  This should not happen.  You can also usually interrupt
# packmol after 30 minutes, and the solution at that point (an .xyz file)
# should be good enough for use.
packmol < step1_proteins.inp   # This step determines the protein's location
                               # It takes ~20 minutes (on an intel i7)
packmol < step2_innerlayer.inp # this step builds the inner monolayer
                               # It takes ~90 minutes
packmol < step3_outerlayer.inp # this step builds the outer monolayer
                               # It takes ~4 hours


# NOTE: PLEASE USE "packmol", NOT "ppackmol".  ("ppackmol" is the
#       parallel-version of packmol using OpemMP.  This example has NOT been
#       tested with "ppackmol".  Our impression was that the "ppackmol"
#       version is more likely to get stuck in an infinite loop. -Andrew 2015-8)


# Step3 creates a file named "step3_outerlayer.xyz" containing the coordinates
# in all the atoms of all the molecules.  Later we will run moltemplate.sh
# using the "-xyz ../system.xyz" command line argument.  That will instruct
# moltemplate to look for a file named "system.xyz" in the parent directory.
# So I rename the "step3_outerlayer.xyz" file to "system.xyz", and move it
# to this directory so that later moltemplate.sh can find it.

mv -f step3_outerlayer.xyz ../system.xyz

