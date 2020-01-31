Long-range dispersion settings
==============================

The PPPM method computes interactions by splitting the pair potential
into two parts, one of which is computed in a normal pairwise fashion,
the so-called real-space part, and one of which is computed using the
Fourier transform, the so called reciprocal-space or kspace part.  For
both parts, the potential is not computed exactly but is approximated.
Thus, there is an error in both parts of the computation, the
real-space and the kspace error. The just mentioned facts are true
both for the PPPM for Coulomb as well as dispersion interactions. The
deciding difference - and also the reason why the parameters for
pppm/disp have to be selected with more care - is the impact of the
errors on the results: The kspace error of the PPPM for Coulomb and
dispersion interaction and the real-space error of the PPPM for
Coulomb interaction have the character of noise. In contrast, the
real-space error of the PPPM for dispersion has a clear physical
interpretation: the underprediction of cohesion. As a consequence, the
real-space error has a much stronger effect than the kspace error on
simulation results for pppm/disp.  Parameters must thus be chosen in a
way that this error is much smaller than the kspace error.

When using pppm/disp and not making any specifications on the PPPM
parameters via the kspace modify command, parameters will be tuned
such that the real-space error and the kspace error are equal.  This
will result in simulations that are either inaccurate or slow, both of
which is not desirable. For selecting parameters for the pppm/disp
that provide fast and accurate simulations, there are two approaches,
which both have their up- and downsides.

The first approach is to set desired real-space an kspace accuracies
via the *kspace\_modify force/disp/real* and *kspace\_modify
force/disp/kspace* commands. Note that the accuracies have to be
specified in force units and are thus dependent on the chosen unit
settings. For real units, 0.0001 and 0.002 seem to provide reasonable
accurate and efficient computations for the real-space and kspace
accuracies.  0.002 and 0.05 work well for most systems using lj
units. PPPM parameters will be generated based on the desired
accuracies. The upside of this approach is that it usually provides a
good set of parameters and will work for both the *kspace\_modify diff
ad* and *kspace\_modify diff ik* options.  The downside of the method
is that setting the PPPM parameters will take some time during the
initialization of the simulation.

The second approach is to set the parameters for the pppm/disp
explicitly using the *kspace\_modify mesh/disp*, *kspace\_modify
order/disp*, and *kspace\_modify gewald/disp* commands. This approach
requires a more experienced user who understands well the impact of
the choice of parameters on the simulation accuracy and
performance. This approach provides a fast initialization of the
simulation. However, it is sensitive to errors: A combination of
parameters that will perform well for one system might result in
far-from-optimal conditions for other simulations. For example,
parameters that provide accurate and fast computations for
all-atomistic force fields can provide insufficient accuracy or
united-atomistic force fields (which is related to that the latter
typically have larger dispersion coefficients).

To avoid inaccurate or inefficient simulations, the pppm/disp stops
simulations with an error message if no action is taken to control the
PPPM parameters. If the automatic parameter generation is desired and
real-space and kspace accuracies are desired to be equal, this error
message can be suppressed using the *kspace\_modify disp/auto yes*
command.

A reasonable approach that combines the upsides of both methods is to
make the first run using the *kspace\_modify force/disp/real* and
*kspace\_modify force/disp/kspace* commands, write down the PPPM
parameters from the output, and specify these parameters using the
second approach in subsequent runs (which have the same composition,
force field, and approximately the same volume).

Concerning the performance of the pppm/disp there are two more things
to consider. The first is that when using the pppm/disp, the cutoff
parameter does no longer affect the accuracy of the simulation
(subject to that gewald/disp is adjusted when changing the cutoff).
The performance can thus be increased by examining different values
for the cutoff parameter. A lower bound for the cutoff is only set by
the truncation error of the repulsive term of pair potentials.

The second is that the mixing rule of the pair style has an impact on
the computation time when using the pppm/disp. Fastest computations
are achieved when using the geometric mixing rule. Using the
arithmetic mixing rule substantially increases the computational cost.
The computational overhead can be reduced using the *kspace\_modify
mix/disp geom* and *kspace\_modify splittol* commands. The first
command simply enforces geometric mixing of the dispersion
coefficients in kspace computations.  This introduces some error in
the computations but will also significantly speed-up the
simulations. The second keyword sets the accuracy with which the
dispersion coefficients are approximated using a matrix factorization
approach.  This may result in better accuracy then using the first
command, but will usually also not provide an equally good increase of
efficiency.

Finally, pppm/disp can also be used when no mixing rules apply.
This can be achieved using the *kspace\_modify mix/disp none* command.
Note that the code does not check automatically whether any mixing
rule is fulfilled. If mixing rules do not apply, the user will have
to specify this command explicitly.


