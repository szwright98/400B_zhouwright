import numpy as np
import astropy.units as u
from CenterOfMassMod import CenterOfMass
from ReadFile import Read
from MassCalc import MassProfile


"""
I want to map the tidal features of M31 and MW as they merge, with the ultimte goal of figuring out what happens to the tidal features after the merger is complete. Per Prof Besla's suggestion in my proposal, I will be finding the jacobi radius at each particle's position and picking out those outside the jacobi radius, and those that hav a high enough escape velocity to escape A) just their host galaxy and B) the combined gravitational attraction of both galaxies. The second part is to figure out how much mass the galaxy loses as a whole, to figure out whetherthe tidal f eatures are lost or not.  
"""

# We want to iterate over the entire time, once and only once
# need mass profile object for each step
# want to store at each step: particles that are outside the jacobi radius of their galaxy, particles that are escaping their own galaxy, particles tha are escaping the merged galaxies
# have to calculate jacobi radius for each particle, in addition to the escape velocity of each particle

for i in range(10):
    print i
    # declare all the mass profile objects
    MP = MassProfile('MW','M31',i)
    #MP = MassProfile('M31',i)
    #MP = MassProfile('MW','M31',i,True)
    #MP = MassProfile('M31',i,True,'MW')

    a = 0.5 # pretty much the same for both

    # get particles outside the jacobi radius
    
    # get enclosed mass at each particles radius, then get Jacobi radius at each particle
    Menc = MP.MassEnclosedTotal(MP.R1)
    JRadii = MP.JacobiRadius(MP.R1,MP.m1,Menc)
    
    # make index that only pulls out the outliers

    out_ind = np.where(MP.R1 > JRadii)



    # get particles that are escaping

    bulge = MP.MassEnclosed(MP.R1,1)
    disk = MP.MassEnclosed(MP.R1,2)
    halo = MP.MassEnclosed(MP.R1,3)

    EVel = MP.EVelocity(MP.R1,a,MP.z1,bulge,disk,halo)

    # make index that only pulls out escaping objects


    esc_ind = np.where(MP.V1 > EVel)

    # make all the arrays or whatever

    outliers  = np.zeros((1,8))
    time_outliers = np.full((len(MP.m1[out_ind]),1),i)
    outliers = np.vstack((outliers,[time_outliers,MP.m1[out_ind],MP.x1[out_ind],MP.y1[out_ind],\
                                    MP.z1[out_ind],MP.vx1[out_ind],MP.vy1[out_ind],\
                                    MP.vz1[out_ind]]))
    escapers  = np.zeros((1,8))
    time_escapers = np.full((len(MP.m1[esc_ind]),1),i)
    escapers = np.vstack((escapers,[time_escapers,MP.m1[esc_ind],MP.x1[esc_ind],MP.y1[esc_ind],\
                                    MP.z1[esc_ind],MP.vx1[esc_ind],MP.vy1[esc_ind],\
                                    MP.vz1[esc_ind]]))
    
    
outliers = np.delete(outliers,0,0)
escapers = np.delete(escapers,0,0)


np.savetxt("MW_Outliers.txt", outliers, fmt = "%11.3f"*8, comments='#',
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

np.savetxt("MW_Escapers.txt", escapers, fmt = "%11.3f"*8, comments='#',
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

