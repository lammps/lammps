# Pizza.py viz of line output

d = dump("dump1.atom dump2.atom")
ld = ldump("dump1.line dump2.line")
ld.map(1,"id",2,"type",3,"end1x",4,"end1y",5,"end2x",6,"end2y")
d.extra(ld)

g = gl(d)
g.arad(1,0.2)

g.lrad(1,5)
g.lcol(1,"yellow")
v = vcr(g)
