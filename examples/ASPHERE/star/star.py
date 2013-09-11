# create set of randomly oriented stars

p = patch(0.2)
p.dim = 2
p.extratype = 1
p.style = "sphere"
p.extra = "Molecules"
p.seed = 5932894
#p.seed = 54321
p.build(30,"star2d",5,0.8,1)
p.write("data.star")
