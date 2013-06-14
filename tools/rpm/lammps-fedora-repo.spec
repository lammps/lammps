Name:           lammps-fedora-repo
Version:        1
Release:        1
Summary:        LAMMPS-ICMS Snapshot Repository Configuration

Group:          System Environment/Base
License:        BSD
URL:            https://sites.google.com/site/akohlmey/software/lammps-icms
Source0:        lammps-fedora.repo

BuildArch:      noarch
ExclusiveArch:  noarch

%description
This package contains the Yum package manager configuration files for the
LAMMPS-ICMS snapshot repository of precompiled LAMMPS binaries for Fedora
Linux installations

%prep
# no setup

%build
# no build

%install
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/etc/yum.repos.d/
install -m 0644 %{SOURCE0} $RPM_BUILD_ROOT/%{_sysconfdir}/yum.repos.d/

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,-)
%{_sysconfdir}/yum.repos.d/*.repo

%changelog

* Wed Jun  12 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 20130612-1
- Initial Fedora style SPEC file


