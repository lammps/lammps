# create set of randomly oriented boxes

p = patch(0.15)
p.dim = 2
p.extratype = 1
p.style = "sphere"
p.extra = "Molecules"
p.seed = 596982
p.build(30,"box2d",6,3,0.8,1)
p.write("data.box")
