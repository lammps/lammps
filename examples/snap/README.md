See `compute_snap_dgrad.py` for a test that compares the dBi/dRj from compute snap (`dbirjflag=1`) to the sum of dBi/dRj from usual compute snap (`dbirjflag=0`).

The format of the global array from `dbirjflag=1` is as follows.

The first N rows belong to bispectrum components for each atom, if we use `bikflag=1`. The first `K` columns correspond to each bispectrum coefficient. The final 3 columns contain reference force components for each atom.

The rows after the first N rows contain dBi/dRj values for all pairs. These values are arranged in row chunks for each atom `j`, where all the rows in a chunk are associated with the neighbors `i` of `j`, as well as the self-terms where `i=j`. So for atom `j`, the number of rows is equal to the number of atoms within the SNAP cutoff, plus 1 for the `i=j` terms, times 3 for each Cartesian component. The total number of dBi/dRj rows is therefore equal to `N*(Nneigh+1)*3`, and `Nneigh` may be different for each atom. To facilitate with determining which row belong to which atom pair `ij`, the last 3 columns contain indices; the 3rd to last column contains global indices of atoms `i` (the neighbors), the 2nd to last column contains global indices of atoms `j`, and the last column contains an index 0,1,2 for the Cartesian component. Like the `bik` rows, the first `K` columns correspond to each bispectrum coefficient.

Finally, the first column of the last row contains the reference energy.
