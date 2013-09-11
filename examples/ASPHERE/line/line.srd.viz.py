# Pizza.py viz of line + SRD output

d = dump("dump1.atom.srd dump2.atom.srd")
ld = ldump("dump1.line.srd dump2.line.srd")
ld.map(1,"id",2,"type",3,"end1x",4,"end1y",5,"end2x",6,"end2y")
d.extra(ld)

g = gl(d)
g.arad(1,0.2)

g.arad(2,0.05)
g.acol(2,"green")

g.lrad(1,5)
g.lcol(1,"yellow")
v = vcr(g)
