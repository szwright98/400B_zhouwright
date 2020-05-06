import numpy as np
import astropy.units as u
import CenterOfMassMod
from ReadFile import Read
from MassCalc import MassProfile


"""
I want to map the tidal features of M31 and MW as they merge, with the ultimte goal of figuring out what happens to the tidal features after the merger is complete. Per Prof Besla's suggestion in my proposal, I will be finding the jacobi radius at each particle's position and picking out those outside the jacobi radius, and those that hav a high enough escape velocity to escape A) just their host galaxy and B) the combined gravitational attraction of both galaxies. The second part is to figure out how much mass the galaxy loses as a whole, to figure out whetherthe tidal f eatures are lost or not.  
"""

# We want to iterate over the entire time, once and only once
# need mass profile object for each step
# want to store at each step: particles that are outside the jacobi radius of their galaxy, particles that are escaping their own galaxy, particles tha are escaping the merged galaxies
# have to calculate jacobi radius for each particle, in addition to the escape velocity of each particle


# make all the arrays

   
outliers  = np.zeros((1,8))
escapers  = np.zeros((1,8))
    
i = 250
while i < 451:
    print i
    # declare all the mass profile objects
    MP = MassProfile('MW','M31',i)
    # need this oe to get the mass enclosed for jacobi radius (since it depends on the other galaxy's mass)
    MP2 = MassProfile('M31','MW',i)

    
    # need COM ojects to calculate Jacobi Radius
    ME_COM1 = CenterOfMassMod.CenterOfMass(MP.filename1,2) # COM found ONLY with one galaxy
    ME_COMP1 = ME_COM1.COM_P(0.1,2)
    ME_COMV1 = ME_COM1.COM_V(ME_COMP1[0].value,ME_COMP1[1].value,ME_COMP1[2].value)
    ME_COM2 = CenterOfMassMod.CenterOfMass(MP.filename2,2) # COM found ONLY with one galaxy
    ME_COMP2 = ME_COM2.COM_P(0.1,2) - ME_COMP1 # shift center of mass location of satellite's radius (since we only want the magnitude, it doesn't matter which frame we shift too)
    #radius between the galaxies at this snapshot
    GalR = np.sqrt(ME_COMP2[0]**2 + ME_COMP2[1]**2 + ME_COMP2[2]**2).value

    a = 0.5 # pretty much the same for both


    
    # get particles outside the jacobi radius
    # get enclosed mass at each particles radius, then get Jacobi radius at each particle
    Menc = MP2.MassEnclosedTotal(GalR) # mass of the host glaxy enclosed at the satellite's radius
    Msat = MP.MassEnclosedTotal(GalR) #enclosed mass of the satellite within GalR
    JRadii = MP.JacobiRadius(GalR,Msat,Menc)
    #print("JRadii Calculated")
  
    # convert outliers location to the COM frame

    x1New = MP.x1[MP.diskindex] - ME_COMP1[0].value
    y1New = MP.y1[MP.diskindex] - ME_COMP1[1].value
    z1New = MP.z1[MP.diskindex] - ME_COMP1[2].value
    vx1New = MP.vx1[MP.diskindex] - ME_COMV1[0].value
    vy1New = MP.vy1[MP.diskindex] - ME_COMV1[1].value
    vz1New = MP.vz1[MP.diskindex] - ME_COMV1[2].value
    #print zNew[0] #troubleshoot
    R1New = np.sqrt(x1New**2 + y1New**2 + z1New**2)
    V1New = np.sqrt(vx1New**2 + vy1New**2 + vz1New**2)

    # make index that only pulls out the outliers
    out_ind = np.where(R1New > JRadii)
    #print(R1New) #troubleshooting
    #print(JRadii) #troubleshooting
    #print ("Index Created")


    # get particles that are escaping
    #EVel = MP.EVelocity(R1New,a,z1New)
    #print(V1New)
    #print(EVel) #troubleshooting
    # make index that only pulls out escaping objects
    #esc_ind = np.where(V1New > EVel)



    

    # only do the above if that snapshot actually has outliers
    if len(MP.m1[out_ind]) > 0:
        # first organize the data of all the outlying particles into columns of mass, x, y, z, vx, vz, vy
        time_outliers = np.full(len(MP.m1[out_ind]),i)
        # make a temporary array storing all the relevant values for this snap
        outliers_temp = np.array([time_outliers, MP.m1[out_ind],x1New[out_ind],y1New[out_ind],\
                                        z1New[out_ind],vx1New[out_ind],vy1New[out_ind],\
                                        vz1New[out_ind]])
       
        
        outliers = np.vstack((outliers,outliers_temp.T))
    i += 10

    """ 
    # exact same process as above 
    if len(MP.m1[esc_ind]) > 0:
       
        time_escapers = np.full(len(MP.m1[esc_ind]),i)
    
        escape_temp = np.array([time_escapers,MP.m1[esc_ind],x1New[esc_ind],y1New[esc_ind],\
                                        z1New[esc_ind],vx1New[esc_ind],vy1New[esc_ind],\
                                        vz1New[esc_ind]])
        #print(np.shape(escape_temp)) #troubleshooting
        #print(len(MP.x1))
        
        escapers = np.vstack((escapers,escape_temp.T)) # need transverse of the array  
        #print(escapers) # troubleshooting
	""" 
    
   
	
outliers = np.delete(outliers,0,0)
#escapers = np.delete(escapers,0,0)


np.savetxt("MW_Outliers250to450.txt", outliers, fmt = "%11.3f"*8, comments='#',
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
"""
np.savetxt("MW_Escapers500.txt", escapers, fmt = "%11.3f"*8, comments='#',
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
"""
