import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename,particle_type,particle_num):
    #read in the text file 
    time, part_total, data =  Read(filename)
    #make an index to only select the particle type we want
    index = np.where(data['type'] == particle_type)
    #make new arrays only with the data from the desired particle type
    xnew = data['x'][index]
    ynew = data['y'][index]
    znew = data['z'][index]

    vxnew = data['vx'][index]
    vynew = data['vy'][index]
    vznew = data['vz'][index]

    mnew = data['m'][index]
    #pull the desired particle's information from the reduced arrays 
    xpos = xnew[particle_num]
    ypos = ynew[particle_num]
    zpos = znew[particle_num]
    xvel = vxnew[particle_num]
    yvel = vynew[particle_num]
    zvel = vznew[particle_num]
    mass = mnew[particle_num]
    #acquire the final3D position from the 3 listed positions
    pos = np.round((np.sqrt(xpos**2+ypos**2+zpos**2)),3)
    pos = pos * u.kpc 
    #acquire the final 3D velocity from the 3 listed components
    vel = np.round((np.sqrt(xvel**2+yvel**2+zvel**2)),3)
    vel = vel * u.km / u.s 
    #return the information as a tuple
    return pos, vel, mass
#Call function, and store output as tuple 100th particle is at index 99 
result = ParticleInfo("MW_000.txt",2,99)
#unpack tuple into independent variables
pos = result[0]
vel = result[1]
mass = result[2]
#print out 
print("3D Distance is: {}".format(pos))
print("3D Velocity is: {}".format(vel))
print("Mass is: {}".format(mass))
#convert and print
print("3D Distance in ly is: {}".format(pos.to(u.lyr)))
