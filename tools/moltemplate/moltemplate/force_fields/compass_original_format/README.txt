This directoroy contains the COMPASS force field parameters in the original
MSI file format ("compass_published.frc" which is distributed with "msi2lmp").
It can be converted into moltemplate format using the following command:

msifrc2lt.py -name COMPASS < compass_published.frc > compass_published.lt

--- Credits: ----
   This is an incomplete version of the COMPASS force field based on available
public sources.  Parameters for some common chemical groups are missing
(for example, the NH2 amine group).  The commercial version of COMPASS is
much larger and presumably includes more up to date parameters and parameters
for a wider range of atom types and molecule types.  (However files
containing those force field parameters are not publicly available.)

   This file has been graciously made available by Materials Design Inc.

   Here is a comment from "compass_published.frc":

   "This file created by Materials Design, Inc. (www.materialsdesign.com)
Please realize that we neither support this version, nor make any warranty
as to the correctness of the parameters.  We have checked the numbers against
the literature, but of course there may still be errors, including errors of
interpretation. Also, the current version of COMPASS may well be different 
that that originally published.
    If you have comments or suggestions, feel free to email Paul Saxe
at psaxe (at) materialsdesign.com"

(Note: This file predates moltemplate and was intended for use with other
       software. Paul Saxe cannot be expected to answer questions related to
       moltemplate.)
