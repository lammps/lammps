# Harmonic-surface-mapping-algorithm-in-LAMMPS
[HSMA3D](https://aip.scitation.org/doi/10.1063/1.5044438) and [HSMA2D](https://aip.scitation.org/doi/10.1063/5.0003293) (with planar dielectric interfaces) have been implemented into LAMMPS as a k-space module. This module is written via C++ and is paralleled via MPI + OpenMP. We recommend user install 'user-omp' package in Lammps. Fewer MPIs and more OpenMPs are the most efficient choice. We suggest not to use pure MPI. With optimal choice of parameters, the speed of this package is comparable to PPPM (with Intel optimization) in LAMMPS, or even faster than it. 

## Installation
For employing HSMA3D after download this full package, the first thing is to include the USER-HSMA package and other appropriate packages in your LAMMPS (cd ./src catalogue):
```
make yes-user-hsma yes-molecule yes-manybody yes-kspace
```

We recommend to install 'USER-OMP' package for the better performance of the evaluation of the LJ potential:
```
make yes-user-omp
```

Finally, compile the LAMMPS under ./src catalogue using
```
make intel_cpu_intelmpi
```

Note that if you want to use HSMA with your own LAMMPS, please just copy /src/USER-HSMA to your lammps/src catalogue, and then refer to the same installation procedure listed in this part.

## User guide
For employing HSMA3D, the only thing needed is to change the k-space solver in your Lammps in-file, just as 
```
kspace_style HSMA3D 1.0e-3 1.3 8 128 55.0 89.0 0 0
```

The first two parameters don't need to modify. Just keep them as `"kspace_style HSMA".`

1.0e-3 : the tolerance of FMM or directly compute.  

1.3 : If the dimensions of simulation box are Lx=Ly=Lz=100, then this parameter means that one find all the images within a sphere with radius 1.3 * sqrt(LxLx+LyLy+LzLz)
.   

8 : means that the number of basis is 8*8=64.  

128 : the number of monitoring points.  

55.0 and 89.0 : the parameters of Fibonacci quadrature. One can set these two numbers as two adjacent Fibonacci numbers, like "89.0 144.0" or "144.0 233.0".

1 : Indicate if one employs FMM (O(N)) to evaluate the potential of near-field. "0" indicates directly computing which has O(N^2) complexity.

0 : Indicate if one employs FMM (O(N)) to evaluate the potential of far-field. "0" indicates directly computing which has O(N^2) complexity.

For employing HSMA2D, the only thing needed is to change the k-space solver in your Lammps in-file, just as 
```
kspace_style HSMA2D 1.0e-3 1.5 0.0 6 40 16 55.0 89.0 1 0
```  
The first two parameters don't need to modify. Just keep them as `"kspace_style HSMA2D".`  

1.0e-3 : the tolerance of FMM or directly compute.  

1.5 : If the dimensions of simulation box are Lx=Ly=Lz=100, then this parameter means that one find all the images within a sphere with radius 1.3 * sqrt(Lx*Lx+Ly*Ly+Lz*Lz)
.  

0.00 : the dielectric mismatch. The range of mismatch is [-1,1]
  
6 : means that the number of basis is 6*6=36. 

40 : the number of Gaussian quadratures.  

16 : the parameter "w" in our paper (relative to the 2D dilation quadrature).

55.0 and 89.0 : the parameters of Fibonacci quadrature. One can set these two numbers as two adjacent Fibonacci numbers, like "89.0 144.0" or "144.0 233.0".

1 : Indicate if one employs FMM (O(N)) to evaluate the potential of near-field. "0" indicates directly computing which has O(N^2) complexity.

0 : Indicate if one employs FMM (O(N)) to evaluate the potential of far-field. "0" indicates directly computing which has O(N^2) complexity.

For more details of parameter setting, please refer to our SCI papers which contain the set of parameter within given accuracy. 

Note that the 'pair style' should be set as 'lj/cut' (or lj/cut/omp, we recommend using user-omp package in Lammps) if you want to evaluate LJ potential. Please do not use pair styles which are combined with the near part of a Coulomb solver, such as'lj/cut/coul/long', etc. 

## Examples
Some examples of in-file are provided in the folder `HSMA_Example`.
```
salt_1-1.in : 1:1 electrolyte solution

salt_2-1.in : 2:1 electrolyte solution

salt_3-1.in : 3:1 electrolyte solution

SPCE_Water.spce-bulk-nvt : the SPCE water system (17496 atoms)
```

To set the number of OpenMP threads per MPI, please type
```
export OMP_NUM_THREADS=10
```
in the command line or dynamiclly set in the code (not recommand).

To run the input script, an example is cd ./HSMA_Example and then type
```
srun --mpi=pmi2 ../src/lmp_intel_cpu_intelmpi -i salt_3-1.in
```
in the command line and hit Enter.

## Citation
In this package, we utilize an efficient implementation of FMM ([FMM3D](https://github.com/flatironinstitute/FMM3D)) which is developed by Greengard's group. 

If you use this package in your work and feel that this package is helpful to you, please cite one (or more) of the following papers in your citation：

1. [Q. Zhao, J. Liang, and Z. Xu, J. Chem. Phys. 149, 084111 (2018).](https://aip.scitation.org/doi/10.1063/1.5044438)
2. [J. Liang, J. Yuan, E. Luijten, and Z. Xu, J. Chem. Phys. 152, 134109 (2020).](https://aip.scitation.org/doi/10.1063/5.0003293)
3. [J. Liang, J. Yuan, and Z. Xu, preprint (2021).](https://arxiv.org/abs/2104.05260)

This version still need optimization. If you have any questions and suggestions, please send a email to liangjiuyang@sjtu.edu.cn (both Chinese and English are OK).

Good luck to all of you!

## `Introduction of HSMA3D and HSMA2D (with planar dielectric interfaces)`  
Harmonic Surface Mapping Algorithm for 3D periodic systems, first described in paper [Harmonic surface mapping algorithm for fast electrostatic sums](https://aip.scitation.org/doi/10.1063/1.5044438) published by The Journal of Chemical Physics, is an efficient implementation for electrostatic pairwise sums of an infinite number of images accelerated by Fast Multiple method(FMM) and graphics processing units(GPU) (cuda codes are available in another repository). Numerical calculations of the Madelung constant, electrostatic energy of ions in a metallic cavity, and the time performance for large-scale systems show that the HSMA is accurate and fast, and thus is attractive for many applications.

Our recent work, [Harmonic Surface Mapping Algorithm for 2D periodic systems with planar dielectric interfaces](https://aip.scitation.org/doi/10.1063/5.0003293), is published in The Journal of Chemical Physics. We have developed an accurate and efficient method for molecular dynamics simulations of charged particles confined by planar dielectric interfaces. The algorithm combines the image-charge method for near field with the harmonic surface mapping, which converts the contribution of infinite far-field charges into a finite number of charges on an auxiliary spherical surface. We approximate the electrostatic potential of far-field charges via spherical harmonic expansion and determine the coefficients by fitting the Dirichlet-to-Neumann boundary condition, which only requires the potential within the simulation cell. Instead of performing the direct evaluation of spherical harmonic series expansion, we use Green’s second identity to transform the series expansion into a spherical integral, which can be accurately represented by discrete charges on the sphere. Therefore, the fast multipole method can be readily employed to sum over all charges within and on the sphere, achieving truly linear O(N) complexity. Our algorithm can be applied to a broad range of charged complex fluids under dielectric confinement.

Now, we implement two recently developed fast Coulomb solvers, HSMA3D [J. Chem. Phys. 149 (8) (2018) 084111] and HSMA2D [J. Chem. Phys. 152 (13) (2020) 134109], into a new user package [HSMA](https://arxiv.org/abs/2104.05260) for molecular dynamics simulation engine LAMMPS. The HSMA package is designed for efficient and accurate modeling of electrostatic interactions in 3D and 2D periodic systems with dielectric effects at the O(N) cost. The implementation is hybrid MPI and OpenMP parallelized and compatible with existing LAMMPS functionalities. The vectorization technique following AVX512 instructions is adopted for acceleration. To establish the validity of our implementation, we have presented extensive comparisons to the widely used particle-particle particle-mesh (PPPM) algorithm in LAMMPS and other dielectric solvers. With the proper choice of algorithm parameters and parallelization setup, the package enables calculations of electrostatic interactions that outperform the standard PPPM in speed for a wide range of particle numbers.

```
                           Jiuyang Liang
                           Ph.D candidate
                           School of Mathematical Science and Institute of Natural Science
                           Shang Hai Jiao Tong University
                           [Homepage in Github](https://github.com/LiangJiuyang/)
                           [Homepage in Researchgate](https://www.researchgate.net/profile/Liang-Jiuyang)
```
