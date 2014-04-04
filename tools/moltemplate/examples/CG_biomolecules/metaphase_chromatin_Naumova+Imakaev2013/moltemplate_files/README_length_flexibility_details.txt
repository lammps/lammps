---- Andrew's comments ----

The two-stage model at the end of Naumova et al Science 2013 uses the "30nm-fiber" model, whose details are (somewhat vaguely) described in the supplemental materials section.

For the 10nm model, 
    n=128000, 
    L=200, 
    U(alpha)=5*(1 - cos(alpha))
    bond_length=1.0     (=10nm)
    sigma=1.0 (particle radius = 10nm)

30nm-fiber model details:
"The 30nm-like fiber was modeled by increasing the volume of each monomer and the amount of DNA represented by each monomer by a factor of 4.25, while keeping other parameters the same at the monomer level."

I interpret this to mean that, for the 30nm model, 
    n=128000/4.25~=30117
    L=200/4.25~=47
    U(alpha)=1.17647*(1 - cos(alpha))    (5/4.25=1.17647)
To increase the volume by a factor o 4.25, originally I thought I should increase the "sigma" parameter from 1.0 to 4.25^(1/3)~=1.6198.  But I suspect that the bond-lengths between monomers should be fixed at 1.0.  If that is the case, then, perhaps I should increase "sigma" from 1.0 to 4.25^(1/2)~=2.061552, and keep the bond-length fixed at 1.0 (which in the units used by thsi paper, corresponds to 10nm).  (That would increase the volume of a cylinder of radius "sigma" and length="bond-length" by a factor of 4.25)
    bond_length=1.0 (10nm again)
    sigma=2.061552  (Yes, this is less than 3.0<-->30nm. See below.)




--- Excerpts from the Supplemental section of Naumova et al Science 2013 ---

From p. 18 of the supplemental materials section of Naumova et al Science 2013.

    (This section was probably written by Maxim Imakaev.)

    In vivo, the structure of the chromatin fiber can be complicated and many details remain unknown, particularly in metaphase. Given this uncertainty, we simulated chromatin as a homogeneous “beads-on-a-string” polymer fiber. We consider a 10nm fiber, as the pervasiveness of the 30nm fiber in vivo has become increasingly contested. In our simulations, 77Mb is represented by a densely-packed 10nm fiber of 128,000 monomers. Each monomer represents a 10nm-sized DNA-histone complex containing 3 nucleosomes (around 600bp). The fiber has a persistence length of 4 monomers (~2.4Kb), which is based on earlier estimates of 5-10 nucleosomes for interphase (14). Those estimates arise from the assumption that 5-10 linker DNA fragments, each of 20-40bp, can collectively provide flexibility equal to that of the 150bp persistence length of DNA. Binding of proteins to the linker DNA (e.g. histone H1) and interactions between neighboring nucleosomes can further constrain dynamics, requiring more linkers to provide the persistence length. Due to the tight packing of nucleosomes in metaphase, we use the upper limit of this range, i.e. 12 nucleosomes.
    For the consecutive loops on a scaffold model (the final folded state model with the best agreement with Hi-C data), we also performed simulations with a more flexible 10nm fiber, or with a 30nm fiber, and found similar results. The more flexible 10nm fiber was modeled by decreasing persistence length to 1.8 monomers. The 30nm-like fiber was modeled by increasing the volume of each monomer and the amount of DNA represented by each monomer by a factor of 4.25, while keeping other parameters the same at the monomer level. We note that a classic model of a 30nm fiber is much less dense than a compact metaphase chromosome. A textbook model of a 30nm fiber assumes packing of about 6 nucleosomes per 10nm of fiber length. This model predicts that only 28% of the volume of the fiber (a 30nm-diameter cylinder) is occupied by nucleosomes, assuming a nucleosome shell volume of 328 nm^3. This is much less than the estimated 30-50% density of nucleosomes in a metaphase chromosome, assuming a diameter of 600nm, a packing density of 50-70 Mb/um, and the same nucleosome volume. (See also (15), which gives an estimate of 0.14-0.18 pg/μm for DNA only, and would give about twice the density if DNA is counted with nucleosomes). As follows, these fibers would have to interdigitate, and fill in gaps within each other. We account for this overlap by assuming the effective diameter of the fiber to be less than 30nm. The effective diameter was chosen to make the volume of the fiber equal the total volume of all the nucleosomes.

   We accounted for topoisomerase II activity by allowing chromatin fibers to pass through each other, while still having excluded volume interactions. This was achieved by using a soft-core Lennard-Jones potential with 1kT energy cost for monomer overlap (see below). This allows for changes in the topological state of a chromosome that are known to occur during compaction in vivo.

Our simulations of a two-step folding process show that Hi-C data for mitotic chromosomes is consistent with a linearly compressed array of consecutive chromatin loops. Whereas mechanisms for formation of consecutive chromatin loops have been proposed, the process of axial compression is less understood. Chromatid compression cannot be accomplished by increased chromatin-chromatin affinity alone, as this would lead to condensation into a globular geometry (14, 16, 17). However, mechanisms which locally compress the fiber of loop bases naturally allow for anisotropic compression into a shorter and thicker fiber, with the same width regardless of chromosome length (18). Differences in the duration or efficiency of the first and second stages of chromosomal condensation provide a natural mechanism for condensation-related proteins to separately affect mitotic chromosome length and width (19). We also note that the axis of loop-bases in our two-stage model does not necessarily form a continuous and rigid scaffold (Figure S26). As follows, we remain agnostic about the molecular details of the chromosomal scaffold, which might for example be formed by a network consisting of protein-protein and/or protein-DNA interactions (20).

   1. Polymer simulations

    To perform Langevin dynamics polymer simulations we used OpenMM, a high-performance GPU-assisted molecular dynamics API (21, 22). To represent chromatin fibers as polymers, we used a sequence of spherical monomers of 1 unit of length in diameter. Here and below all distances are measured in monomer sizes, set to be 10nm unless specified otherwise.  Neighboring monomers are connected by harmonic bonds, with a potential U = 100*(r - 1)^2 (here and below in units of kT). Polymer stiffness is modeled with a three point interaction term, with the potential U = 5*(1 - cos(alpha)), where alpha is the angle between neighboring bonds.  All monomers interact via either a shifted Lennard-Jones (LJ) repulsive potential, or an attractive Lennard-Jones potential. At high densities in a confined volume, the details of the inter-monomer interactions become negligible due to screening (23), and we therefore used the computationally efficient shifted LJ potential. The shifted LJ potential allows for a short-range repulsion by truncating the LJ potential at its minimum and shifting the minimum to zero: U = 4 * (1/r^12 - 1/r^6) + 1, for r<2^(1/6); U=0 for r > 2^(1/6). The shifted LJ potential is one of the most computationally efficient repulsive potentials due to a very short cutoff radius.

    To allow chain passing, which represents activity of topoisomerase II, we softened the shifted LJ potential by truncating the interaction energy at Ecutoff = 1 kT. At energies more than 0.5 Ecutoff, the LJ potential was softened via: Usoftened = 0.5 * Ecutoff * (1 + tanh(2*U/Ecutoff - 1)). To avoid numerical 19instabilities in the calculation of U at r ~ 0, the interaction radius r was truncated at r=0.3 via: rtruncated = (r^10 + (0.3)^10)^0.1, which introduced negligible shift in a final softened potential. For an attractive LJ potential, we used: U = 4 * e * (1/r^12 - 1/r^6), with e = 0.46 kT, slightly below the theta-temperature.  The attractive potential was similarly softened at 2 kT and cut off at r=2.5. Unless noted, we used a softened shifted LJ repulsive potential.

    Polymer models were visualized using Pymol and Rasmol. For images with loop bases highlighted, a base of each loop and 3 monomers surrounding it in each direction were labeled in red.

   SECTIONS 2-5 SKIPPED

6. Two-stage process: linear compaction - axial compression

To simulate the two-stage process of metaphase chromosome folding, we used the 30nm fiber representation described above for its computational efficiency. Simulations were initialized from 30000 monomer fractal globule conformations; fractal globule is a model for interphase chromatin organization. First, random consecutive loops with L=100 monomers (see above) were introduced, and anchors of neighboring loops were brought together using harmonic springs with a potential U = k * (r – r0)2; r0=0.5. To avoid abrupt motion of the loop anchors, the force was gradually turned on over the first
300000 timesteps, with k linearly increasing in time from 0 to 10 kT. We used softened shifted repulsive LJ potential for inter-monomer interaction.

Upon completion of linear compaction, axial compression was initiated. This involves following changes: the repulsive LJ force is replaced with an attractive LJ force for all monomers, and the chromosomal core of loop anchors is homogeneously compressed. To achieve the latter, all anchor pairs separated by less than 30 anchors were attracted via a potential U = step(d-3) * abs(d-3) * 10 kT, which implements a constant attractive force between two anchors if they are separated by a distance larger than 3. The interactions between neighboring loop anchors were kept throughout this process.

To obtain the contact map from this simulation, 50 independent runs of 1.5e7 timesteps were performed, and 250 conformations were collected from the second half of each run. The contact map was calculated from all conformations of all runs at a 30-monomer resolution, and was further averaged over three 10000-monomer blocks along the diagonal of the heatmap. The latter was done to show contact map at a relevant length scale (0 to 25 Mb), and to achieve a better averaging of the contact map.

