import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy.units as u
from MassCalc import MassProfile

# declare class, radii, scale length
MProf = MassProfile("MW",0)
r = np.arange(0.25,30.5,1.5)
a = 0.5

# get mass profiles
MPHalo = MProf.MassEnclosed(r,1)
MPDisk = MProf.MassEnclosed(r,2)
MPBulge = MProf.MassEnclosed(r,3) #comment out for M33
MPTotal = MProf.MassEnclosedTotal(r)
MPHern = MProf.HernquistProfile(r,a,MPHalo)

# get velocity profiles
VHalo = MProf.CircularVelocity(r,1)
VDisk = MProf.CircularVelocity(r,2)
VBulge = MProf.CircularVelocity(r,3) #comment out for M33
VTotal = MProf.CircularVelocityTotal(r)
VHern = MProf.HernquistVCirc(r,a,MPHalo)

fig,ax = plt.subplots(figsize=(10,8))

label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

#section for mass profile
"""
plt.plot(r, MPHalo,linewidth = 3, color= 'k', label='Halo Mass')
plt.plot(r, MPDisk,linewidth = 3, color ='b', label='Disk Mass')
#plt.plot(r, MPBulge,linewidth = 3, color ='r', label='Bulge Mass') #comment out for M33
plt.plot(r, MPTotal,linewidth = 5, color ='g', label='Total Mass')
plt.plot(r, MPHern,linewidth = 5, color ='orange', label='Hernquist Mass, a=0.3')
"""
#section for velocity profile
plt.plot(r, VHalo,linewidth = 3, color= 'k', label='Halo Velocity')
plt.plot(r, VDisk,linewidth = 3, color ='b', label='Disk Velocity')
plt.plot(r, VBulge,linewidth = 3, color ='r', label='Bulge Velocity') #comment out for M33
plt.plot(r, VTotal,linewidth = 5, color ='g', label='Total Mass Velocity Curve')
plt.plot(r, VHern,linewidth = 5, color ='orange', label='Hernquist Profile Velocity, a=0.5')

# Axes labels 
plt.xlabel('Radius kpc',fontsize=22) 
#plt.ylabel('Enclosed Mass Msun', fontsize=22) # for mass
plt.ylabel('Circular Velocity km/s', fontsize=22) # for velocity
#plt.semilogy # comment out for velocity

# Legend and title
plt.legend(loc='upper left',fontsize='large')
plt.title("MW Velocity Profile")


plt.savefig("MW_VelocityProfile.png")
plt.show()
