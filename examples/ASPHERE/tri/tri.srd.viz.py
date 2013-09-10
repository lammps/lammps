# Pizza.py viz of triangle + SRD output

d = dump("dump1.atom.srd dump2.atom.srd")
t = tdump("dump1.tri.srd dump2.tri.srd")
t.map(1,"id",2,"type",
      3,"corner1x",4,"corner1y",5,"corner1z",
      6,"corner2x",7,"corner2y",8,"corner2z",
      9,"corner3x",10,"corner3y",11,"corner3z")
d.extra(t)

g = gl(d)
g.arad(1,0.02)
g.acol(1,"green")

g.arad(2,0.05)
g.acol(2,"green")

v = vcr(g)
