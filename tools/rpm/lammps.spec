# verified on Fedora 19     / x86_64 / 2016-05-11
# verified on Fedora 19     / i386   / 2016-05-11
# verified on Fedora 20     / x86_64 / 2016-05-11
# verified on Fedora 20     / i386   / 2016-05-11
# verified on Fedora 21     / x86_64 / 2016-05-11
# verified on Fedora 22     / x86_64 / 2016-05-11
# verified on Fedora 23     / x86_64 / 2016-05-11
# verified on Fedora 24     / x86_64 / 2016-05-11
# verified on CentOS 6.6    / x86_64 / 2016-05-11
# verified on CentOS 6.6    / i386   / 2016-05-11
# verified on CentOS 7.0    / x86_64 / 2016-05-11
# verified on OpenSuSE 12.3 / x86_64 / 2016-05-11
# verified on OpenSuSE 12.3 / i586   / 2016-05-11
# verified on OpenSuSE 13.1 / x86_64 / 2016-05-11
# verified on OpenSuSE 13.1 / i586   / 2016-05-11
# verified on OpenSuSE 13.2 / x86_64 / 2016-05-11
# verified on OpenSuSE 42.1 / x86_64 / 2016-05-11

%ifnarch s390 s390x
%global with_openmpi 1
%else
%global with_openmpi 0
%endif

%ifarch x86_64
%global bigintsize -DLAMMPS_SMALLBIG
%else
%global bigintsize -DLAMMPS_SMALLSMALL
%endif

%if %{defined suse_version}
%global with_suse 1
%global with_mpich2 0
%global _openmpi_load :
%global _openmpi_unload :
%else
%global with_suse 0
%if %{defined _mpich2_load}
%global with_mpich2 1
%else
%global with_mpich2 0
%endif
%endif

# to find the proper location for installing the lammps.py* files
%{!?python_sitearch: %global python_sitearch %(%{__python} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1))")}

Name:           lammps
Version:        20160511
Release:        14%{?dist}
Summary:        LAMMPS Molecular Dynamics Simulator
Group:          Applications/Engineering

License:        GPLv2
URL:            http://lammps.sandia.gov
Source0:        lammps-current.tar.gz
Source1:        lammps.sh
Source2:        lammps.csh

BuildRequires:  gcc-c++
BuildRequires:  fftw-devel
BuildRequires:  libpng-devel
BuildRequires:  python-devel
%if %{with_suse}
BuildRequires:  gcc-fortran
BuildRequires:  libjpeg8-devel
%else
BuildRequires:  gcc-gfortran
BuildRequires:  libjpeg-devel
%endif
Requires:       lammps-common = %{version}-%{release}

%global lammps_desc \
LAMMPS is an acronym for Large-scale Atomic/Molecular Massively Parallel\
Simulator. LAMMPS has potentials for soft materials (biomolecules, polymers)\
and solid-state materials (metals, semiconductors) and coarse-grained or\
mesoscopic systems. It can be used to model atoms or, more generically,\
as a parallel particle simulator at the atomic, meso, or continuum scale.\
LAMMPS runs on single processors or in parallel using message-passing\
techniques and a spatial-decomposition of the simulation domain, as well\
as using multi-threading via OpenMP inside each simulation domain.

%description
%{lammps_desc}

This package contains a LAMMPS executable compiled without MPI support.

%package common
Summary:        LAMMPS utilities and potentials
Group:          Applications/Engineering
Requires:       python

%description common
%{lammps_desc}

This package contains common utilities and potential files

%package doc
Summary:        LAMMPS documentation, example inputs, and benchmark tests
Group:          Applications/Engineering

%description doc
%{lammps_desc}

This package contains the LAMMPS manual and addtional documentation as
well as example and benchmark test inputs.

%if %{with_openmpi}
%package openmpi
Summary:        LAMMPS OpenMPI executable
Group:          Applications/Engineering
Requires:       lammps-common
Requires:       openmpi
BuildRequires:  openmpi-devel

%description openmpi
%{lammps_desc}

This package contains a parallel LAMMPS executable for OpenMPI.
%endif

%if %{with_suse}
# no out-of-the-box support for MPICH2
%else
# older Fedora and RHEL/CenOS use MPICH2
%if %{with_mpich2}
%package mpich2
Summary:        LAMMPS MPICH2 executable
Group:          Applications/Engineering
Requires:       lammps-common
Requires:       mpich2
BuildRequires:  mpich2-devel

%description mpich2
%{lammps_desc}

This package contains a parallel LAMMPS executable for MPICH2.
%else
# package names and files for MPICH3 went back to plain MPICH
%package mpich
Summary:        LAMMPS MPICH executable
Group:          Applications/Engineering
Requires:       lammps-common
Requires:       mpich
BuildRequires:  mpich-devel >= 3.0.2

%description mpich
%{lammps_desc}

This package contains a parallel LAMMPS executable for MPICH.
%endif
%endif

%package python
Summary:        LAMMPS Python module
Group:          Applications/Engineering
Requires:       lammps-common
Requires:       python

%description python
%{lammps_desc}

This package contains the LAMMPS Python module

%prep
%setup -q -n lammps-current

%build
# build supporting libraries for MPI stubs
cd lib/atc
make -f Makefile.g++ CC=g++ CCFLAGS="-fPIC -I../../src -I../../src/STUBS  ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.linalg
cd ../awpmd
make -f Makefile.mpicc CC=g++ CCFLAGS="-fPIC -Isystems/interact/TCP/ -Isystems/interact -Iivutils/include ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.linalg
cd ../colvars
make -f Makefile.g++ CXX=g++ CXXFLAGS="-fPIC ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.empty
cd ../linalg
make -f Makefile.gfortran FC=gfortran FFLAGS="-fPIC ${RPM_OPT_FLAGS}" FFLAGS0="${RPM_OPT_FLAGS} -O0 -fPIC"
cd ../meam
make -f Makefile.gfortran F90=gfortran F90LAGS="-fPIC ${RPM_OPT_FLAGS}" EXTRAMAKE=Makefile.lammps.gfortran
cd ../poems
make -f Makefile.g++ CC=g++ CCFLAGS="-fPIC ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.empty
cd ../voronoi
make -f Makefile.g++ CXX=g++ CXXFLAGS="-fPIC ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.empty

# now build in main source directory
cd ../../src

# install packages
# fortran reax is obsolete, no GPU support, QM/MM requires a Q-E library, USER-INTEL requires Intel Compiler, USER-LB and MPIIO require MPI-IO, user-vtk requires vtk.
make yes-all no-kokkos no-kim no-gpu no-user-cuda no-user-qmmm no-user-intel no-user-quip no-user-lb no-mpiio no-user-h5md no-user-vtk

make -C STUBS

make g++ CC=g++ CCFLAGS="${RPM_OPT_FLAGS} -fopenmp -fPIC" LINK=g++ LINKFLAGS="${RPM_LD_FLAGS} -fopenmp" LMP_INC="-DLAMMPS_GZIP -DLAMMPS_FFMPEG -DLAMMPS_JPEG -DLAMMPS_PNG -DLAMMPS_MEMALIGN=64 %{bigintsize} " MPI_INC="-I../STUBS" MPI_PATH="-L../STUBS" MPI_LIB="-lmpi_stubs -lrt" FFT_INC=-DFFT_FFTW3 FFT_LIB=-lfftw3 JPG_LIB="-ljpeg -lpng" mode=shexe OMP=yes LIBRT=yes

# build shared library for python bindings
make g++ CC=g++ CCFLAGS="${RPM_OPT_FLAGS} -fopenmp -fPIC" LINK=g++ LINKFLAGS="${RPM_LD_FLAGS} -fopenmp" LMP_INC="-DLAMMPS_GZIP -DLAMMPS_FFMPEG -DLAMMPS_JPEG -DLAMMPS_PNG -DLAMMPS_MEMALIGN=64 %{bigintsize} " MPI_INC="-I../STUBS" MPI_PATH="-L../STUBS" MPI_LIB="-lmpi_stubs -lrt" FFT_INC=-DFFT_FFTW3 FFT_LIB=-lfftw3 JPG_LIB="-ljpeg -lpng" mode=shlib OMP=yes LIBRT=yes

# stash executable and shared lib away
cd ../
mkdir serial
mv src/lmp_g++ serial/
mv src/liblammps_g++.so serial/liblammps.so

# byte compile python script wrapper
cd python
ln -s ../serial/liblammps.so
cat > dummy.py <<EOF
from lammps import lammps
lmp = lammps()
EOF
LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH} %{__python} -O dummy.py
cd ../

# now build some tools
g++ -o serial/restart2data ${RPM_OPT_FLAGS} %{bigintsize} tools/restart2data.cpp
g++ -o serial/binary2txt ${RPM_OPT_FLAGS} %{bigintsize} tools/binary2txt.cpp
gfortran -o serial/chain.x ${RPM_OPT_FLAGS} tools/chain.f
gfortran -o serial/createatoms ${RPM_OPT_FLAGS} tools/createatoms/createAtoms.f
g++ -o serial/lammpsplot ${RPM_OPT_FLAGS} tools/xmgrace/lammpsplot.cpp
make -C tools/msi2lmp/src CC=gcc CCFLAGS="${RPM_OPT_FLAGS} %{bigintsize}"
mv tools/msi2lmp/src/msi2lmp.exe serial/msi2lmp
make -C tools/colvars CXX=g++ CXXFLAGS="${RPM_OPT_FLAGS}" EXT=
mv tools/colvars/abf_integrate serial/abf_integrate
cp tools/msi2lmp/README tools/msi2lmp/README.msi2lmp
cp tools/createatoms/Manual.pdf tools/createatoms/createatoms.pdf

# build OpenMPI parallel version, if supported
%if %{with_openmpi}
%{_openmpi_load}
# need to rebuild lib/atc
cd lib/atc
make -f Makefile.g++ clean
make -f Makefile.g++ CC=mpicxx CCFLAGS="-fPIC -I../../src -DMPICH_IGNORE_CXX_SEEK -DOMPI_SKIP_MPICXX=1 ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.linalg

# back to main source dir
cd ../../src
make clean-g++

# enable USER-LB and MPIIO since we have a full MPI library now
make yes-user-lb yes-mpiio
make g++ CC=mpicxx CCFLAGS="${RPM_OPT_FLAGS} -fopenmp" LINK=mpicxx LINKFLAGS="${RPM_LD_FLAGS} -fopenmp" LMP_INC="-DLAMMPS_GZIP -DLAMMPS_FFMPEG -DLAMMPS_JPEG -DLAMMPS_PNG -DLAMMPS_MEMALIGN=64 %{bigintsize} " MPI_INC="" MPI_PATH="" MPI_LIB=-lrt FFT_INC=-DFFT_FFTW3 FFT_LIB=-lfftw3 JPG_LIB="-ljpeg -lpng" OMP=yes LIBRT=yes

# and save the executable
cd ../
mkdir openmpi
mv src/lmp_g++ openmpi/
%{_openmpi_unload}
%endif

%if %{with_suse}
# no MPICH2 build with SuSE
%else
# build MPICH parallel version
%if %{with_mpich2}
%{_mpich2_load}
%else
%{_mpich_load}
%endif
# need to rebuild lib/atc
cd lib/atc
make -f Makefile.g++ clean
make -f Makefile.g++ CC=mpicxx CCFLAGS="-fPIC -I../../src -DMPICH_IGNORE_CXX_SEEK -DOMPI_SKIP_MPICXX=1 ${RPM_OPT_FLAGS} %{bigintsize}" EXTRAMAKE=Makefile.lammps.linalg

# back to main source dir
cd ../../src
make clean-g++
# enable USER-LB since we have a full MPI library now
make yes-user-lb yes-mpiio

make g++ CC=mpicxx CCFLAGS="${RPM_OPT_FLAGS} -fopenmp" LINK=mpicxx LINKFLAGS="${RPM_LD_FLAGS} -fopenmp" LMP_INC="-DLAMMPS_GZIP -DLAMMPS_FFMPEG -DLAMMPS_JPEG -DLAMMPS_PNG -DLAMMPS_MEMALIGN=64 %{bigintsize} " MPI_INC="" MPI_PATH="" MPI_LIB="-lrt" FFT_INC=-DFFT_FFTW3 FFT_LIB=-lfftw3 JPG_LIB="-ljpeg -lpng" OMP=yes LIBRT=yes

# and save the executable
cd ../
mkdir mpich
mv src/lmp_g++ mpich/
%if %{with_mpich2}
%{_mpich2_unload}
%else
%{_mpich_unload}
%endif
%endif

# build done (so far)

%install
rm -rf %{buildroot}
mkdir -p %{buildroot}%{_bindir}
install -p -m 755 serial/lmp_g++ %{buildroot}%{_bindir}
install -p -m 755 serial/restart2data %{buildroot}%{_bindir}
install -p -m 755 serial/binary2txt %{buildroot}%{_bindir}
install -p -m 755 serial/chain.x %{buildroot}%{_bindir}
install -p -m 755 serial/createatoms %{buildroot}%{_bindir}
install -p -m 755 serial/lammpsplot %{buildroot}%{_bindir}
install -p -m 755 serial/msi2lmp %{buildroot}%{_bindir}
install -p -m 755 serial/abf_integrate %{buildroot}%{_bindir}

%if %{with_openmpi}
%if %{with_suse}
mkdir -p $RPM_BUILD_ROOT/%{_libdir}/mpi/gcc/openmpi/bin
install -p -m 755 openmpi/lmp_g++ $RPM_BUILD_ROOT/%{_libdir}/mpi/gcc/openmpi/bin/
%else
%{_openmpi_load}
mkdir -p $RPM_BUILD_ROOT/%{_libdir}/openmpi/bin
install -p -m 755 openmpi/lmp_g++ $RPM_BUILD_ROOT/%{_libdir}/openmpi/bin/
%{_openmpi_unload}
%endif
%endif

%if %{with_suse}
# no MPICH2 support in SuSE
%else
%if %{with_mpich2}
%{_mpich2_load}
mkdir -p $RPM_BUILD_ROOT/%{_libdir}/mpich2/bin
install -p -m 755 mpich/lmp_g++ $RPM_BUILD_ROOT/%{_libdir}/mpich2/bin/
%{_mpich2_unload}
%else
%{_mpich_load}
mkdir -p $RPM_BUILD_ROOT/%{_libdir}/mpich/bin
install -p -m 755 mpich/lmp_g++ $RPM_BUILD_ROOT/%{_libdir}/mpich/bin/
%{_mpich_unload}
%endif
%endif

mkdir -p $RPM_BUILD_ROOT/%{python_sitearch}
mkdir -p $RPM_BUILD_ROOT/%{_libdir}
mkdir -p $RPM_BUILD_ROOT/%{_datadir}/lammps
mkdir -p $RPM_BUILD_ROOT/%{_sysconfdir}/profile.d/
cp python/lammps.py* $RPM_BUILD_ROOT/%{python_sitearch}
cp serial/liblammps.so $RPM_BUILD_ROOT/%{_libdir}
sed -e s,@DATADIR@,%{_datadir}/lammps, < %{SOURCE1} > $RPM_BUILD_ROOT/%{_sysconfdir}/profile.d/lammps.sh
sed -e s,@DATADIR@,%{_datadir}/lammps, < %{SOURCE2} > $RPM_BUILD_ROOT/%{_sysconfdir}/profile.d/lammps.csh
cp -arp potentials $RPM_BUILD_ROOT/%{_datadir}/lammps/
cp -arp tools/msi2lmp/frc_files $RPM_BUILD_ROOT/%{_datadir}/lammps/
cp -arp bench $RPM_BUILD_ROOT/%{_datadir}/lammps/
cp -arp examples $RPM_BUILD_ROOT/%{_datadir}/lammps/
mkdir -p $RPM_BUILD_ROOT/%{_datadir}/lammps/examples/msi2lmp
cp -ap tools/msi2lmp/test/*.{mdf,car} $RPM_BUILD_ROOT/%{_datadir}/lammps/examples/msi2lmp/
cp -ap tools/msi2lmp/test/reference/*.data $RPM_BUILD_ROOT/%{_datadir}/lammps/examples/msi2lmp/
cp -ap tools/msi2lmp/test/in.* $RPM_BUILD_ROOT/%{_datadir}/lammps/examples/msi2lmp
mkdir -p $RPM_BUILD_ROOT/%{_datadir}/lammps/examples/createatoms
cp -ap tools/createatoms/create.input $RPM_BUILD_ROOT/%{_datadir}/lammps/examples/createatoms

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%{_bindir}/lmp_g++

%files common
%defattr(-,root,root,-)
%{_bindir}/restart2data
%{_bindir}/binary2txt
%{_bindir}/chain.x
%{_bindir}/createatoms
%{_bindir}/lammpsplot
%{_bindir}/msi2lmp
%{_bindir}/abf_integrate
%{_sysconfdir}/profile.d/lammps.*sh
%{_datadir}/lammps/potentials
%{_datadir}/lammps/frc_files
%doc README LICENSE
%doc tools/createatoms/createatoms.pdf
%doc tools/msi2lmp/README.msi2lmp
%doc tools/xmgrace/lammpsplot.pdf

%files doc
%defattr(-,root,root,-)
%{_datadir}/lammps/bench
%{_datadir}/lammps/examples
%doc doc/Manual.pdf
%doc doc/src/PDF/*.pdf

%if %{with_openmpi}
%files openmpi
%defattr(-,root,root,-)
%if %{with_suse}
%{_libdir}/mpi/gcc/openmpi/bin/lmp_g++
%else
%{_libdir}/openmpi/bin/lmp_g++
%endif
%endif

%if %{with_suse}
# no MPICH2 package for suse
%else
%if %{with_mpich2}
%files mpich2
%defattr(-,root,root,-)
%{_libdir}/mpich2/bin/lmp_g++
%else
%files mpich
%defattr(-,root,root,-)
%{_libdir}/mpich/bin/lmp_g++
%endif
%endif

%files python
%defattr(-,root,root,-)
%doc python/README python/examples
%{python_sitearch}/*
%{_libdir}/liblammps.so


%changelog
* Wed May 11 2016 Axel Kohlmeyer <akohlmey@gmail.com> - 20160511-14
- Update for new documentation tree. addition of Fedora 24

* Fri Oct 16 2015 Axel Kohlmeyer <akohlmey@gmail.com> - 20151016-12
- Update for removal of Fedora 17/18 and addition of Fedora 23

* Thu Jan 23 2014 Axel Kohlmeyer <akohlmey@gmail.com> - 20140123-7
- Handle dependencies on MPI-IO for USER-LB and MPIIO. 

* Mon Dec 16 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20131120-6
- Include abf_integrate as tool for colvars

* Wed Nov 20 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20131120-5
- Update for Fedora 20 and OpenSuSE 13.1

* Fri Aug  2 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130802-4
- Rename $BIOSYM_LIBRARY to $MSI2LMP_LIBRARY and biosym_frc_files to frc_files

* Sun Jul  7 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130707-3
- included support for compiling and distributing msi2lmp and included support for setting LAMMPS_POTENTIALS and BIOSYM_LIBRARY environment variables

* Thu Jul  4 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130615-2
- Added flags to compile with high resolution timers

* Sat Jun 15 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130615-1
- added proper installation directory for SuSE OpenMPI version

* Fri Jun 14 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130614-1
- Added bundling (most) example inputs with the doc subpackage

* Fri Jun 14 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130609-1
- Included more documentation pdfs and move potentials and benchmarks to _datadir/lammps

* Thu Jun 13 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130609-1
- Split off manual and benchmarks from common into doc package

* Sun Jun  9 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130609-1
- Added subpackage for python wrapper

* Sun Jun  9 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130609-1
- Initial Fedora/RedHat style SPEC file

