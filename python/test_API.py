#test pyAVL wrapper
import numpy as np

import pyAVL


case = pyAVL.avlAnalysis(geo_file='aircraft.txt', mass_file ='aircraft.mass')

# case.addConstraint('alpha', 3)
case.addConstraint('elevator', 0.00)
case.addConstraint('rudder', 0.00)

# case.addTrimCondition('CL', 1.0)
# case.addConstraint('alpha', 0.00)

case.executeRun()

print '----------------- Neutral Point ----------------'
case.calcNP()
print case.NP
case.clearVals()


case.alphaSweep(0, 10)

print '----------------- alpha sweep ----------------'
print 'Angle      Cl         Cd         Cm'
for i in xrange(len(case.alpha)):
    print '%8f   %8f   %8f   %8f   '%(case.alpha[i]*(180/np.pi),case.CL[i],case.CD[i],case.CM[i])


case.clearVals()

case.CLSweep(0.6, 1.6)

print '----------------- CL sweep ----------------'
print 'Angle      Cl         Cd         Cm'
for i in xrange(len(case.alpha)):
    print '%8f   %8f   %8f   %8f   '%(case.alpha[i]*(180/np.pi),case.CL[i],case.CD[i],case.CM[i])


   
 # correct output

#     ----------------- alpha sweep ----------------
# Angle      Cl         Cd         Cm
# 0.000000   0.938616   0.063291   0.000000   
# 1.000000   1.011160   0.071276   0.000000   
# 2.000000   1.083098   0.079784   -0.000000   
# 3.000000   1.154380   0.088796   0.000000   
# 4.000000   1.224957   0.098291   0.000000   
# 5.000000   1.294782   0.108248   -0.000000   
# 6.000000   1.363809   0.118642   -0.000000   
# 7.000000   1.431993   0.129449   -0.000000   
# 8.000000   1.499289   0.140643   0.000000   
# 9.000000   1.565657   0.152196   0.000000   
# 10.000000   1.631056   0.164079   -0.000000   
# ----------------- CL sweep ----------------
# Angle      Cl         Cd         Cm
# -4.580877   0.600000   0.033839   -0.000000   
# -3.240347   0.700000   0.041206   0.000000   
# -1.890446   0.800000   0.049685   -0.000000   
# -0.529164   0.900000   0.059283   -0.000000   
# 0.845640   1.000000   0.070009   -0.000000   
# 2.236265   1.100000   0.081868   -0.000000   
# 3.645202   1.200000   0.094868   -0.000000   
# 5.075178   1.300000   0.109014   -0.000000   
# 6.529203   1.400000   0.124311   -0.000000   
# 8.010633   1.500000   0.140764   -0.000000   
# 9.523246   1.600000   0.158375   -0.000000   