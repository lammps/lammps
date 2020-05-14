**Summary**

New version of the code that compiles and runs with mixtures of gaussian and erfc particles.
Hangs and segmentation faults are fixed.
Interaction potentials and energies agree with previous TILD code.


**Author(s)**

Andrew Santos, Sandia National Labs, asanto@sandia.gov

**Implementation Notes**

The code has a more pair-style like format, but is still very flexible.  Most parts of the 

As it is now it is still all in one TILD file. In the future interaction-specific routines could/should be in child classes.  That way there would be little work to add new potentials. As of now `initialize_potential` is the main interaction-specific routine.

I still have questions about how the TILD energy contribution, pressure tensor and the forces have and should be calculated, especially concerning the Helfand compressibility.

Should we be able to reproduce pppm of point charges with Coulombic interactions?

**Tasks**

As of now I see future work, prioritized in decreasing order, including:
0. Fix force calculation, may be an issue with `grad_potent_hat`
1. Nail down how the Helfand  `- rho0` contribution is implemented to the energy, pressure tensor and forces
2. Compare speed one single to multiple cores with previous TILD implementations, and try to improve speed/performance.
3. Timing tests on mulitple nodes
4. If we can, reproduce pppm of point charges with Coulombic interactions
5. Add triclinic box support
6. Add documentation

