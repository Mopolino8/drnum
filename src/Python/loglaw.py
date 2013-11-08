from math import *

y_plus = 1
y_plus_min = 20
y_wall = 0.00002
u = 100
u_tau = 0
for i in range(10):
  y_plus = max(y_plus_min, u_tau*y_wall/1.5e-5)
  u_tau = u/(log(y_plus)/0.41 + 5)
  print "y+ = " + str(y_plus) + "   u_tau = " + str(u_tau)


