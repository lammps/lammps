# create set of randomly oriented tri boxes 

p = patch(0.05)

p.seed = 54321
p.dim = 3
p.extratype = 1
p.lattice = [5,5,5]

p.build(125,"tribox",0.5,2,0.5,2,0.5,2,1)

p.write("data.tri.srd")

print "all done ... type CTRL-D to exit Pizza.py"
