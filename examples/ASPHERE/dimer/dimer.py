# create set of randomly oriented dimers

p = patch(0.3)
p.dim = 2
p.extratype = 1
p.style = "sphere"
p.extra = "Molecules"
p.seed = 54321
p.build(100,"dimer",0.8,1)
p.write("data.dimer")
