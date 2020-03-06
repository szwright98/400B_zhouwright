

# Homework 6 Template
# G. Besla & R. Li

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMassMod import CenterOfMass




def OrbitCom(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
    galaxy, string containing the galaxy name
    start, the number of the first snapshot to be read
    end, the number of the last snapshot to be read
    n, an integer indicating the intervals which will return the COM
    returns: time, position, and velocity of the COM
    """
    
    # compose the filename for output
    fileout = "Orbit_" + galaxy + ".txt"
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    d = 0.1
    VD = 4
    # for M33 that is stripped more, use different values for VolDec (4)
    print galaxy
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
   
    snap_ids = np.arange(start,end,n)
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros((len(snap_ids),7))
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000'+ str(snap_id)
        ilbl = ilbl[-3:]
        filename = "./VLowRes/%s_"%(galaxy)+ilbl+".txt"
        #print filename #troubleshooting
        #initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename,2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COMP = COM.COM_P(d,VD)
        COMV = COM.COM_V(COMP[0],COMP[1],COMP[2])
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value)
        #print COMV[0] #troubleshooting
        orbit[i][0] = COM.time.value/1000
        orbit[i][1] = COMP[0].value
        orbit[i][2] = COMP[1].value
        orbit[i][3] = COMP[2].value
        orbit[i][4] = COMV[0].value
        orbit[i][5] = COMV[1].value
        orbit[i][6] = COMV[2].value
        
        # note that you can store 
        # a[i] = var1, *tuple(array1)

        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 

OrbitCom("M33", 0, 800,5) 

#OrbitCom("M31", 0, 800,5) 

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt




# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  




# Determine the magnitude of the relative position and velocities 

# of MW and M31

# of M33 and M31




# Plot the Orbit of the galaxies 
#################################




# Plot the orbital velocities of the galaxies 
#################################

