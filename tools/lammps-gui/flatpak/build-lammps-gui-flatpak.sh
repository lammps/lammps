#!/bin/sh
flatpak-builder --force-clean --user --install-deps-from=flathub --repo=repo --install builddir org.lammps.lammps-gui.yml
flatpak build-bundle repo lammps-gui.flatpak org.lammps.lammps-gui --runtime-repo=https://flathub.org/repo/flathub.flatpakrepo
