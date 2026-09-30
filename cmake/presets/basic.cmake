# preset that turns on just a few, frequently used packages
# this will be compiled quickly and handle a lot of common inputs.

set(ALL_PACKAGES GRAPHICS KSPACE MANYBODY MOLECULE RIGID)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
