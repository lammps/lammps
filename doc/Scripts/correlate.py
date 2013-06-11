#!/usr/bin/env python
"""
  function: 
    parse the block of thermo data in a lammps logfile and perform auto- and 
    cross correlation of the specified column data.  The total sum of the 
    correlation is also computed which can be converted to an integral by 
    multiplying by the timestep. 
  output:
    standard output contains column data for the auto- & cross correlations 
    plus the total sum of each. Note, only the upper triangle of the 
    correlation matrix is computed.
  usage: 
    correlate.py [-c col] <-c col2> <-s max_correlation_time>  [logfile] 
"""
import sys
import re
import array

# parse command line

maxCorrelationTime = 0
cols = array.array("I")
nCols = 0
args = sys.argv[1:]
index = 0
while index < len(args):
  arg = args[index]
  index += 1
  if   (arg == "-c"):
    cols.append(int(args[index])-1)
    nCols += 1
    index += 1
  elif (arg == "-s"):
    maxCorrelationTime = int(args[index])
    index += 1
  else :
    filename = arg
if (nCols < 1): raise RuntimeError, 'no data columns requested'
data = [array.array("d")]
for s in range(1,nCols)  : data.append( array.array("d") )

# read data block from log file

start = False
input = open(filename)
nSamples = 0
pattern = re.compile('\d')
line = input.readline()
while line :
  columns = line.split()
  if (columns and pattern.match(columns[0])) :
    for i in range(nCols): 
      data[i].append( float(columns[cols[i]]) )
    nSamples += 1
    start = True
  else :
     if (start) : break
  line = input.readline()
print "# read :",nSamples," samples of ", nCols," data"
if( maxCorrelationTime < 1): maxCorrelationTime = int(nSamples/2);
 
# correlate and integrate

correlationPairs = []
for i in range(0,nCols):
  for j in range(i,nCols): # note only upper triangle of the correlation matrix
    correlationPairs.append([i,j])
header = "# "
for k in range(len(correlationPairs)):
  i = str(correlationPairs[k][0]+1)
  j = str(correlationPairs[k][1]+1)
  header += " C"+i+j+" sum_C"+i+j
print header
nCorrelationPairs = len(correlationPairs)
sum = [0.0] * nCorrelationPairs
for s in range(maxCorrelationTime)  :
  correlation = [0.0] * nCorrelationPairs
  nt = nSamples-s
  for t in range(0,nt)  :
    for p in range(nCorrelationPairs):
      i = correlationPairs[p][0]
      j = correlationPairs[p][1]
      correlation[p] += data[i][t]*data[j][s+t]
  output = ""
  for p in range(0,nCorrelationPairs):
    correlation[p] /= nt
    sum[p] += correlation[p]
    output += str(correlation[p]) + " " + str(sum[p]) + " "
  print output
