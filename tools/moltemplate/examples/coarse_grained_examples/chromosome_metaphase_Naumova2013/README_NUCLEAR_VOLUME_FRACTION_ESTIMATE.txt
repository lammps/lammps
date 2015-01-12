---------------
The average diameter of a mammalian cell nucleus is is 6 micrometers (μm),
which occupies about 10% of the total cell volume.

(See "Molecular Biology of the Cell, Chapter 4, pages 191–234 (4th ed.)", by
 Bruce Alberts, Alexander Johnson, Julian Lewis, Martin Raff, Keith Roberts, Peter Walter, 2002)

... of that, 25% of it is occupied by the nucleolus
http://en.wikipedia.org/wiki/Nucleolus
("citation needed")
---------------

From the supplemental material for the original HiC paper
(Lieberman-Aiden et al., Science 2009)

Appendix 1.
   Estimate of the Volume Fraction of Chromatin in Human Cells
In the simulations we sought to obtain an ensemble of structures that, in their statistical properties, resemble some of the features of chromatin arrangement in the cell.  Below we demonstrate that chromatin occupies a significant fraction of cell volume, a property that we reproduced in simulations. Taking the nuclear diameter of a tissue culture cell to be 5-10um, and assuming close to a spherical shape we obtain the volume in the range 50-500 um^3, with a (geometric) mean of ~160 um^3. If we assume that the chromatin is built of DNA wrapped around nucleosomes, then we have 6x10^9bp/200bp=3x10^7 nucleosomes. Each may be approximated as a cylinder ~10nm in diameter and ~5nm in height, suggesting a volume of about 500nm3 each. The linker DNA after each nucleosome is about 50bps long, suggesting a volume of about 50*0.34nm*3.14*1nm^2=50nm^3. Thus the total volume of chromatin = 550x3x10^7 =16 um^3, or ~10% (3-23%) of the nuclear volume. This strikingly large volume fraction is itself a significant underestimate, since we ignored, among other things, all other DNA-bound proteins. Note that any further packing or localization of chromatin inside the nucleus will increase local density.
---- This next section mostly only justifies why they   ----
---- they did not stop the simulation when the globules ----
---- were fully crumpled (ie with uniform density)      ----
   In our simulations, the radius of the final crumpled globule was R≈12.5 and the volume V≈8000 cubic units. The total volume of the 4000 monomers, 1 unit in diameters each, is V≈2000. This implies a volume fraction of about 25%, which is consistent with the volume fraction estimated above.
 ----  ----

Appendix 2.
   Monomer length in base pairs
Each monomer of the chain corresponds to a fragment of chromatin that equals the Kuhn length of the chromatin fiber, i.e. approximately twice the persistence length of the fiber. Although the persistence length of the chromatin fiber is unknown it can be estimated using the following arguments. DNA is packed into nucleosomes, where 150 bps are wrapped around the histone core and do not contribute to flexibility of the fiber. The linker DNA of about 50 bps that connects consecutive nucleosomes is bendable, and is the source of flexibility in the fiber. Since the persistence length of double-stranded DNA is 150 bps, an equally flexible region of the nucleosomal DNA should contain 3 linkers, i.e. 3 consecutive nucleosomes packing about 600 bps of DNA. The excluded volume of the nucleosomes, nucleosome interactions, and other DNA-bound proteins can make the fiber less flexible or prohibit certain conformation and may tend to increase the persistence length of the fiber. Using this estimated lower bound estimate for the persistence length, we obtain the Kuhn length of the equivalent freely-jointed chain to be 6 nucleosomes, or ~ 1200bp. A simulated chain of 4000 monomers corresponds to 4.8Mb of packed DNA. The size of each monomer was chosen such that its volume is equal to (or slightly above) that of 6 nucleosomes (V=6 x 600 nm^3); thus the radius of the spherical monomer is R=10nm. The diameter of each globule shown in Figure 4 is about 200 nm.


