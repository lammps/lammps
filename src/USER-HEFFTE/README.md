# Prototype kernels for using heFFTe in LAMMPS

[HeFFTe](https://bitbucket.org/icl/heffte/) provides a C++11 implementation of 2D and 3D FFT computation on hybrid architectures. For experimentation with LAMMPS, we have created a package `USER-HEFFTE`, which contains modifications of `PPPM`routines within KSPACE package available by default in LAMMPS.

#### Requirements

A one dimensional FFT library, heFFTe currently supports:

- CPU backends: [FFTW](http://www.fftw.org/), [MKL](https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/fourier-transform-functions.html) 

- GPU backends: [CUFFT](https://developer.nvidia.com/cufft), [RocFFT](https://github.com/ROCmSoftwarePlatform/rocFFT)


#### Enabling heFFTe package

We have added a file called `Makefile.heffte` for compilation. This follows the same standards as other packages when using the installation with GNU Make according to LAMMPS manual, copy the Makefile into the source folder:

```
    cd lammps/src/
    cp MAKE/OPTIONS/Makefile.heffte MAKE/.
```

Install and compile accordingly to LAMMPS procedure for GNU Make:
```
    make yes-user-heffte
    make -j heffte
```

This will produce the executable `lmp_heffte`.

Another approach would be using CMake, for which we have added heFFTe options to file `KSPACE.cmake`, located in the `cmake/Modules/Packages` folder. This also allows an easier way to benchmark on CPUs (with FFTW by default), or GPUs (with CUFFT by default), via enabling or disabling the flag `Heffte_ENABLE_GPU`.

#### Basic Usage

We have target FFT calls from KSPACE package, in particular the routine `poisson_ik`. We have modified the following files:

-  pppm.cpp
-  pppm_dipole.cpp
-  pppm_dipole_spin.cpp
-  pppm_disp.cpp        
-  pppm.h
-  pppm_disp.h

Where, for the case of `pppm.cpp` we have showed how to enable GPU computations. This can be considered as a prototype, and further modification must and will be performed according to LAMMPS requirements.

Note that heFFTe also provides more and different flags and tuning parameters as those available in FFTMPI; for example, the use of slab decompositions or handling strided data. This can be set manually in the code, or read at runtime from the input file, where the latter approach would require to further modify LAMMPS kernels for input file lecture.

#### Testing

For testing purposes, we have created an input file, `in.heffte` , which is based on the *Rhodopsin* protein benchmark) and is available within the `bench` directory. It provides a 3D FFT benchmark at output. Also, a better analysis of the impact and speedup achieved for FFTs, can be achieved by changing the geometry (`data.rhodo`).

Running examples (Using Summit supercomputer):

- On CPUs:
    ```
    jsrun -n4 -g1 -c1 -a1 ./lmp_heffte -in in.heffte
    ```

- On GPUs:
    ```
    jsrun --smpiargs="-gpu" -n2 -c1 -g1 -a1 ./lmp_heffte -in in.heffte
    ```