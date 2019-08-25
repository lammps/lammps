# create set of randomly oriented line boxes 

p = patch(0.05)

p.seed = 7878382
p.dim = 2
p.extratype = 1
p.lattice = [10,10,1]

p.build(100,"linebox",0.5,2,0.5,2,1)

p.write("data.line.srd")

print "all done ... type CTRL-D to exit Pizza.py"
