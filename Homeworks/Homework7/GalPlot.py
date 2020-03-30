# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# read in data from files

MW = np.genfromtxt('Orbit_MW.txt')
M31 = np.genfromtxt('Orbit_M31.txt')
M33 = np.genfromtxt('Orbit_M33.txt')

M33a = np.genfromtxt('M33_Analytical.txt')

# cols: t x y z vx vy vz

def VecDiff(x1, y1, z1, x2, y2, z2):
    # function to find the resultant vector of two subtracted vectors and returns the magnitude
    # inputs:
    # (x1,y1,z1) are the components of the first vector
    # (x2,y2,z2) are the components of the second vector
    # returns magnitude of resultant vector

    return np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

MW_M31_Pdiff = VecDiff(MW[:,1],MW[:,2],MW[:,3],M31[:,1],M31[:,2],M31[:,3])
M31_M33_Pdiff = VecDiff(M31[:,1],M31[:,2],M31[:,3],M33[:,1],M33[:,2],M33[:,3])
MW_M31_Vdiff = VecDiff(MW[:,4],MW[:,5],MW[:,6],M31[:,4],M31[:,5],M31[:,6])
M31_M33_Vdiff = VecDiff(M31[:,4],M31[:,5],M31[:,6],M33[:,4],M33[:,5],M33[:,6])

M33_M31_analytic = np.sqrt(M33a[:,1]**2 + M33a[:,2]**2 + M33a[:,3]**2)

M33_M31_analytic_v = np.sqrt(M33a[:,4]**2 + M33a[:,5]**2 + M33a[:,6]**2)



"""
plt.plot(MW[:,0], M31_M33_Pdiff, label = "M31-M33", color ='b')
plt.plot(M33a[:,0],M33_M31_analytic, label = "M31-M33 Analytical", color ='r')
plt.title("M31-M33 Separation as a Function of Time")
"""

plt.plot(MW[:,0], M31_M33_Vdiff, label = "M31-M33", color ='b')
plt.plot(M33a[:,0],M33_M31_analytic_v, label = "M31-M33 Analytical", color ='r')
plt.title("M31-M33 Velocities as a Function of Time")

plt.legend()
plt.xlabel("Time in Myr")
#plt.ylabel("Seperation in  kpc")
plt.ylabel("Relative Velocity in kpc/gyr")

#check for the merger
#plt.semilogy()

plt.savefig("M31_M33_Velocity_Comparison.png")
plt.show()
