# Install/unInstall package files in LAMMPS
# edit Makefile.package to include/exclude CUDA library
# do not copy child files if parent does not exist

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
      sed -i -e '/include ..\/..\/lib\/cuda\/Makefile.common/d' ../Makefile.package
      sed -i -e 's/-llammpscuda -lcuda -lcudart -lrt //' ../Makefile.package
      sed -i -e 's/-I..\/..\/lib\/cuda -I$(CUDA_INSTALL_PATH)\/include //' ../Makefile.package
      sed -i -e 's/-L..\/..\/lib\/cuda -L$(CUDA_INSTALL_PATH)\/lib64 -L$(CUDA_INSTALL_PATH)\/lib $(USRLIB_CONDITIONAL) -DLMP_USER_CUDA //' ../Makefile.package
      sed -i '1 i include ..\/..\/lib\/cuda\/Makefile.common' ../Makefile.package
      sed -i -e 's|^PKG_INC =[ \t]*|&-I..\/..\/lib\/cuda -I$(CUDA_INSTALL_PATH)\/include |' ../Makefile.package
      sed -i -e 's|^PKG_PATH =[ \t]*|&-L..\/..\/lib\/cuda -L$(CUDA_INSTALL_PATH)\/lib64 -L$(CUDA_INSTALL_PATH)\/lib $(USRLIB_CONDITIONAL) |' ../Makefile.package
      sed -i -e 's|^PKG_LIB =[ \t]*|&-llammpscuda -lcuda -lcudart -lrt |' ../Makefile.package
  fi

  cp comm_cuda.cpp ..
  cp domain_cuda.cpp ..
  cp modify_cuda.cpp ..
  cp neighbor_cuda.cpp ..
  cp neigh_full_cuda.cpp ..
  cp verlet_cuda.cpp ..

  cp cuda.cpp ..
  cp cuda_neigh_list.cpp ..

  cp comm_cuda.h ..
  cp domain_cuda.h ..
  cp modify_cuda.h ..
  cp neighbor_cuda.h ..
  cp verlet_cuda.h ..

  cp cuda.h ..
  cp cuda_common.h ..
  cp cuda_data.h ..
  cp cuda_modify_flags.h ..
  cp cuda_neigh_list.h ..
  cp cuda_precision.h ..
  cp cuda_shared.h ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e '/include ..\/..\/lib\/cuda\/Makefile.common/d' ../Makefile.package
    sed -i -e 's/-llammpscuda -lcuda -lcudart -lrt //' ../Makefile.package
    sed -i -e 's/-I..\/..\/lib\/cuda -I$(CUDA_INSTALL_PATH)\/include //' ../Makefile.package
    sed -i -e 's/-L..\/..\/lib\/cuda -L$(CUDA_INSTALL_PATH)\/lib64 -L$(CUDA_INSTALL_PATH)\/lib $(USRLIB_CONDITIONAL) -DLMP_USER_CUDA //' ../Makefile.package
  fi

  rm ../comm_cuda.cpp
  rm ../domain_cuda.cpp
  rm ../modify_cuda.cpp
  rm ../neighbor_cuda.cpp
  rm ../neigh_full_cuda.cpp
  rm ../verlet_cuda.cpp

  rm ../cuda.cpp
  rm ../cuda_neigh_list.cpp

  rm ../comm_cuda.h
  rm ../domain_cuda.h
  rm ../modify_cuda.h
  rm ../neighbor_cuda.h
  rm ../verlet_cuda.h

  rm ../cuda.h
  rm ../cuda_common.h
  rm ../cuda_data.h
  rm ../cuda_modify_flags.h
  rm ../cuda_neigh_list.h
  rm ../cuda_precision.h
  rm ../cuda_shared.h
fi
