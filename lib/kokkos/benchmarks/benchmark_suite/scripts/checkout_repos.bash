#!/bin/bash

# Kokkos
if [ ! -d "kokkos" ]; then
  git clone https://github.com/kokkos/kokkos
fi
cd kokkos
git checkout develop
git pull
cd ..

# KokkosKernels
if [ ! -d "kokkos-kernels" ]; then
git clone https://github.com/kokkos/kokkos-kernels
fi
cd kokkos-kernels
git pull
cd ..

# MiniMD
if [ ! -d "miniMD" ]; then
  git clone https://github.com/mantevo/miniMD
fi
cd miniMD
git pull
cd ..

# MiniFE
if [ ! -d "miniFE" ]; then
  git clone https://github.com/mantevo/miniFE
fi
cd miniFE
git pull
cd ..



