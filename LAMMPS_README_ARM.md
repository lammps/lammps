# LAMMPS 5 Jun 2019 Porting Guide (CentOS 7.6)
## Introduction
Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) is a molecular dynamics program from Sandia National Laboratories. LAMMPS makes use of Message Passing Interface (MPI) for parallel communication and is free and open-source software, distributed under the terms of the GNU General Public License. The latest LAMMPS version supports CUDA and OpenCL-based GPU computing. LAMMPS was originally developed under a Cooperative Research and Development Agreement (CRADA) between the United States Department of Energy and three laboratories from private sector firms.The Weather Research and Forecasting (WRF) Model can be used for fine-scale weather simulation and forecasting, which is one of the important application scenarios of high-performance computing (HPC).

For more information about LAMMPS, visit the official LAMMPS website(https://lammps.sandia.gov/download.html).

Programming language: C++

Brief description: open-source molecular dynamics program

Open-source license: GPL 2.0

##### Recommended Version
The recommended version is LAMMPS June 5, 2019
## Environment Requirements
### Hardware Requirements
| Item  | Description |
| ------| ----------- |
| CPU   | kunpeng 920 |
### Software Requirements
| Item  | Version  |  Download Address |
| ----- | ---------|  ---------------- |
| LAMMPS | June 5, 2019 | https://lammps.sandia.gov/tars/ |
| FFTW | 3.3.8 | http://www.fftw.org/fftw-3.3.8.tar.gz |
| Test computing instance | in.lj | lammps-stable_5Jun2019/bench |
### OS Requirements
| Item  | Version  | How to Obtain  |
| ------------ | ------------ | ------------ |
|  CentOS | 7.6  |  https://www.centos.org/download/ |
| Kernel  | 4.14.0-115  |  Included in the OS image. |
## Configuring the Compilation Environment
### Installing dependencies


    yum install time -y
    yum install curl* -y
    yum install csh -y
    yum install zlib* -y
### Installing GNU 9.3


    yum install -y centos-release-scl
    yum install -y devtoolset-9-gcc
    yum install -y devtoolset-9-gcc-c++
    yum install -y devtoolset-9-binutils
    scl enable devtoolset-9 bash
    echo "souce /opt/rh/devtoolset-9/enable" >> /etc/profile
### Installing Open MPI
1. Run the following command to install the system dependency package:


    yum install libxml2* systemd-devel.aarch64 numa* -y
2. Run the following commands to install Open MPI:


    wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.3.tar.gz
    
    tar -zxvf openmpi-4.0.3.tar.gz
    
    cd openmpi-4.0.3
    
    ./configure --prefix=/path/to/OPENMPI --enable-pretty-print-stacktrace --enable-orterun-prefix-by-default  --with-cma --enable-mpi1-compatibility
    
    make -j 16
    
    make install
3. Configure environment variables:

    export PATH=/path/to/OPENMPI/bin:$PATH
    
    export LD_LIBRARY_PATH=/path/to/OPENMPI/lib:$LD_LIBRARY_PATH
### Installing FFTW
1. Run the following command to decompress the FFTW installation package:
   
   tar -xvf fftw-3.3.8.tar.gz
2. Run the following command to switch to the directory generated after decompression:

   cd fftw-3.3.8

3. Run the following command to perform configuration:
   
   ./configure --prefix=/path/to/FFTW --enable-shared --enable-static --enable-fma --enable-neon

4. Run the following commands to compile and install FFTW.

   make -j 96

   make install

5. Run the following commands to load environment variables:
   
   export PATH=/path/to/FFTW/bin:$PATH

   export LD_LIBRARY_PATH=/path/to/FFTW/lib:$LD_LIBRARY_PATH

## Obtaining the Source Code
1. Download the LAMMPS installation package lammps-5Jun19.tar.gz from:https://lammps.sandia.gov/tars/

2. Use the SFTP tool to upload the LAMMPS installation package to the /path/to/LAMMPS directory on the server.

## Compiling and Installing LAMMPS
1. Run the following command to switch to the directory, in which the LAMMPS installation package is stored:
   
   cd /path/to/LAMMPS

2. Run the following command to decompress the LAMMPS installation package:
   
   tar -xvf lammps-5Jun19.tar.gz

3. Run the following command to switch to the directory generated after decompression:
   
   cd lammps-5Jun2019

4. Run the following command to switch to the src directory:
   
   cd src

5. Perform the following operations to modify the MAKE/OPTIONS/Makefile.g++_openmpi file:
   
   a. Run the vi MAKE/OPTIONS/Makefile.g++_openmpi command.

   b. Press I to enter the insert mode and modify lines 54 to 56 in the MAKE/OPTIONS/Makefile.g++_openmpi file. Pay attention to the information in bold.

   ```shell
      FFT_INC =    -DFFT_FFTW -I/path/to/FFTW/include
      FFT_PATH =   -L/path/to/FFTW/lib
      FFT_LIB =    -lfftw3
   ```

   c. Press Esc, type :wq!, and press Enter to save the file and exit.

6. Run the following commands to perform compilation:
   
   make yes-std
   
   make no-lib

   make -j 96 g++_openmpi

## Running and Verifying LAMMPS
1. Run the following command to create the working directory:

   mkdir -p path/to/CASE

2. Perform the following commands to switch to the working directory and copy the computing instance and binary files to the working directory.

   cd /path/to/CASE

   cp /path/to/LAMMPS/lammps-stable_5Jun2019/bench/in.lj ./

   cp /path/to/LAMMPS/lammps-stable_5Jun2019/src/lmp_g++_openmpi ./

3. Run the following command:
   
   mpirun --allow-run-as-root -np 96 --mca btl ^openib ./lmp_g++_openmpi -in in.lj

   Check the value (in the unit of timesteps/s) of Performance in the log.lammps log file. A higher value indicates better performance.

   The following is an example of the test result.

   ```shell
      Performance: 593019.799 tau/day, 1372.731 timesteps/s
      98.1% CPU use with 96 MPI tasks x no OpenMP threads
      MPI task timing breakdown:
      Section |  min time  |  avg time  |  max time  |%varavg| %total
      ---------------------------------------------------------------
      Pair    | 0.028215   | 0.040367   | 0.051175   |   3.4 | 55.41
      Neigh   | 0.0036562  | 0.0049894  | 0.0061901  |   1.2 |  6.85
      Comm    | 0.013968   | 0.026072   | 0.039602   |   4.8 | 35.79
      Output  | 7.905e-05  | 0.000128   | 0.00023195  |   0.0 |  0.18
      Modify  | 0.00056847 | 0.00080538  | 0.00099514 |  0.0 |  1.11
      Other   |            | 0.0004857  |          |      |  0.67
   ```






    






