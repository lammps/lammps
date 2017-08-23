README!

v0.2

This is only part of the SDK force field and is to be used for lipids only. Only parameters from:

Shinoda et al. J. Phys. Chem. B, Vol. 114, No. 20, 2010
http://dx.doi.org/10.1021/jp9107206

are used.

NOTE: We extracted the parameters from that publication from the files 
distributed with the "EMC" tool.  If you use these files, please also cite:
P. J. in ‘t Veld and G. C. Rutledge, Macromolecules 2003, 36, 7358.


This works for any topology built using the following types:

Name    Structure         Charge
NC      -CH2CH2-N-(CH3)3  +1
NH      -CH2CH2-NH3       +1
PH      -PO4-             -1
PHE     -PO4- (PE lipid)  -1
GL      -CH2CH-CH2-
EST1    -CH2CO2-
EST2    -H2CO2-
CMD2    -HC=CH- (cis)
CM      -CH2CH2CH2-
CT      CH3CH2CH2-
CT2     CH3CH2-
W       (H2O)3

This coarse-grainng allows for design of a wide variety of lipids. 

NEW! in v0.2:
  SDK Cholesterol model has been added! Using parameters from:

MacDermaid et al. J. Chem. Phys, 143(24), 243144, 2015
http://dx.doi.org/10.1063/1.4937153

are used. The following types are defined specifically for cholesterol:

Name    Structure       Location
C2T     -CH-(CH3)2      Tail
CM2     -CH2-CH2-       Tail
CM2R    -CH2-CH2-       Ring A
CMDB    -CH2-C=CH-      Ring A/B
CMB     -CH2-CH-CH-     Ring B/C
CMR     -CH-CH2-CH2-    Ring B/C
CMR5    -CH2-CH2-CH-    Ring D
CTB     -CH2-CH-CH3-    Tail
CTBA    -C-CH3          Ring A/B
CTBB    -C-CH3          Ring C/D
OAB     -CH-OH          Ring A

See the provided reference for details on the CG cholesterol topology. 
A 5-10 timestep is used when using cholesterol.

Several limiations, due to missing parameters:
-use of cholesterol with type "NH" is not possible.
-use of cholesterol with type "PHE" is not possible.

---
