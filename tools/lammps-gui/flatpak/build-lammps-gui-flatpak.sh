#!/bin/sh
flatpak-builder --gpg-sign=EEA103764C6C633EDC8AC428D9B44E93BF0C375A --force-clean --user --install-deps-from=flathub --repo=repo --install builddir org.lammps.lammps-gui.yml
flatpak build-bundle repo lammps-gui.flatpak org.lammps.lammps-gui --runtime-repo=https://flathub.org/repo/flathub.flatpakrepo
