Summary: LAMMPS Molecular Dynamics Simulator
Name: lammps
# XXX: set this to the major release
Version: 15_Jan_2010
# XXX: set this to the patch level
Release: 20022010
License: GPL
Group: Applications/Education
Source: %{name}.tar.gz
URL: http://lammps.sandia.gov
Prefix: /usr
Buildroot: %{_tmppath}/%{name}-root
Packager: Axel Kohlmeyer <akohlmey@gmail.com>
BuildRequires: gcc-gfortran, fftw2-devel, lapack-devel, blas-devel

# description take from the website (see URL)
%description
LAMMPS is an acronym for Large-scale Atomic/Molecular Massively Parallel 
Simulator. LAMMPS has potentials for soft materials (biomolecules, polymers)
and solid-state materials (metals, semiconductors) and coarse-grained or 
mesoscopic systems. It can be used to model atoms or, more generically, 
as a parallel particle simulator at the atomic, meso, or continuum scale. 
LAMMPS runs on single processors or in parallel using message-passing 
techniques and a spatial-decomposition of the simulation domain. The code 
is designed to be easy to modify or extend with new functionality. 

%clean
case "$RPM_BUILD_ROOT" in *-root) rm -rf $RPM_BUILD_ROOT ;; esac

# XXX: adjust the -n argument to the directory name in the tar file
%prep
%setup -n lammps-20Feb10
# XXX: install all packages except GPU (for now, since there is no
#      fedora/rpmfusion sanctioned CUDA rpm package yet)
cd src ; make yes-all no-gpu

# serial build
%build
# first build support libs
make -C lib/atc -f Makefile.serial CCFLAGS="${RPM_OPT_FLAGS} -fpermissive -I../../src -I../../src/STUBS"
make -C lib/meam -f Makefile.gfortran F90FLAGS="${RPM_OPT_FLAGS} -fno-second-underscore"
make -C lib/poems -f Makefile.g++ CCFLAGS="${RPM_OPT_FLAGS}"
make -C lib/reax -f Makefile.gfortran F90FLAGS="${RPM_OPT_FLAGS} -fno-second-underscore"
make -C src/STUBS CCFLAGS="${RPM_OPT_FLAGS}"
make -C src serial  CCFLAGS="${RPM_OPT_FLAGS}" \
        FFT_INC=-DFFT_FFTW FFT_LIB=-lfftw \
        meam_SYSLIB=-lgfortran reax_SYSLIB= \
        user-atc_SYSLIB="-lblas -llapack"


%install
case "$RPM_BUILD_ROOT" in *-root) rm -rf $RPM_BUILD_ROOT ;; esac
mkdir -p $RPM_BUILD_ROOT%{prefix}/bin
cp src/lmp_serial  $RPM_BUILD_ROOT%{prefix}/bin/lmp_serial

%files
%defattr(-,root,root,-)
%doc README LICENSE
/usr/bin/lmp_serial

%changelog
* Thu Feb 25 2010    Axel Kohlmeyer <akohlemy@gmail.com> 
- Initial RPM specfile 
