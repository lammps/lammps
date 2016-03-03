Name:           lammps-centos-rhel-repo
Version:        1
Release:        3
Summary:        LAMMPS-ICMS Repository Configuration

Group:          System Environment/Base
License:        BSD
URL:            https://sites.google.com/site/akohlmey/software/lammps-icms
Source0:        lammps-centos-rhel.repo

BuildArch:      noarch
ExclusiveArch:  noarch

%description
This package contains the Yum package manager configuration files for the
LAMMPS-ICMS repository of precompiled LAMMPS binaries for CentOS
and RedHat Enterprise Linux installations.

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
* Sat Jun 15 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 1-2
- Updated repo file to not include source packages. No need.


* Wed Jun  12 2013 Axel Kohlmeyer <akohlmey@gmail.com> - 1-1
- Initial RHEL style SPEC file


