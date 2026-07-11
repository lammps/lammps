# create set of randomly oriented line boxes and triangles

p = patch(0.05)

p.seed = 54321
p.dim = 2
p.lattice = [10,10,1]

p.build(50,"linebox",0.5,2,0.5,2,1)
p.build(50,"linetri",0.5,1.5,0.5,1.5,1)

p.write("data.line")

print "all done ... type CTRL-D to exit Pizza.py"
