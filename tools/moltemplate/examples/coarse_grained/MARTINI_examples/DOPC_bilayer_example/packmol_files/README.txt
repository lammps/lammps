You can use packmol to create a file containing the atomic coordinates
for a system of coarse-grained lipids mixed with water using this command:

If it takes too long for packmol to run, try lowering the tolerance.
(tolerance 2.0 should work)

packmol < mix_lipids+water.inp

