
from math import sqrt

lines = """
1 1 1 -0.8476   0.13513   -0.05627    0.09445
2 1 2  0.4238   0.79058    0.37197   -0.48187
3 1 2  0.4238  -0.33549   -0.65246   -0.50563
4 2 1 -0.8476   0.02105    0.48666    2.81228
5 2 2  0.4238  -0.73085    1.07786    2.68447
6 2 2  0.4238   0.21349    0.19973    1.90016
"""

qqrd2e = 332.06371

q = []
x = []

lines = lines.split('\n')
for line in lines[1:7]:
  print line
  id,imol,itype,qone,xone,yone,zone = line.split()
  q.append(float(qone))
  x.append((float(xone),float(yone),float(zone)))

qp4 = 0.0
qp5 = 0.0
qp6 = 0.0

for i in range(6):
  for j in range(6):
    if i == j: continue
    dx = x[i][0] - x[j][0]
    dy = x[i][1] - x[j][1]
    dz = x[i][2] - x[j][2]
    r = sqrt(dx*dx + dy*dy + dz*dz)
    eng = qqrd2e * q[i]*q[j] / r
    print "Eng of atoms %d-%d = %g",i+1,j+1,eng
    if i+1 == 4: qp4 += eng/q[i]
    if i+1 == 5: qp5 += eng/q[i]
    if i+1 == 6: qp6 += eng/q[i]

print "QPOT full: 4 %g 5 %g 6 %g" % (qp4,qp5,qp6)
