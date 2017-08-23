  # First, rescale and interpolate the positions
  # where the monomers will be located.  (This step is not needed
  # if the coords_orig.raw file already has correct coordinates.)
  # The first argument, 32768 is the number of atoms in the desired file.
  # The second argument, 1.6198059 = 4.25^(1/3), tells interpolate_coords.py
  # to multiply all the coordinates (scale them up) by 1.6198059.

  ./interpolate_coords.py 32768 1.6198059 < coords_orig.raw > coords.raw

  # Then, build the "system.lt" file

  ./generate_system_lt.py 32768 51 < coords.raw > system.lt

  # 32768 is the number of monomers in the polymer
  # (which may be different from the number of coordinates
  # in the "coords_orig.raw" file)  This number will vary
  # depending on how long you want the polymer to be.
  # The second argument "51" is the average interval between
  # condensin anchors (IE the "loop size" in monomers.)
