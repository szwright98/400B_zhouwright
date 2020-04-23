import numpy as np
from astropy import units as u
from ReadFile import Read
import CenterOfMassMod
from astropy.constants import G

"""
Modified the MassProfile class to have functions capable of computing:
Jacobi Radius 
Escape Velocity (and the associated potentials)
Also the code now supports an optional (boolean) feature that allows for the input of two galaxies, which then allows the mass of the other galaxy's particles to be included in the calculation.  

"""

G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
G = G.value

# class to describe the masses and positions of masses in a galaxy
class MassProfile:
	def __init__(self,galaxy1,galaxy2,snap,IncludeCompanion = False):
		# initialize with:
		# galaxy1, string that contains the name of the galaxy USED TO FIND THE CENTER OF MASS i.e. "MW"
                #galaxy2, string that contains name of other galaxy
		# snap, integer that gives the snapshot number, corresponding to time
	        # IncludeCompanion, to develop a mass profile that includes the mass of both galaxies

                self.IncludeCompanion = IncludeCompanion
		# use galaxy and snap to reconstruct string name
		ilbl1 = '000' +str(snap)
		# remove all but three last digits
		ilbl1 = ilbl1[-3:]
		self.filename1 = "./VLowRes/%s_"%(galaxy1) + ilbl1 + '.txt'
	        #print (self.filename) #troubleshooting
		# usual file reading
		self.time1, self.total1, self.data1 = Read(self.filename1)
		self.gname1 = galaxy1 
		# pull data values. no need to distinguish by particle type (yet)
		self.m1 = self.data1['m']
		self.x1 = self.data1['x']
		self.y1 = self.data1['y']
		self.z1 = self.data1['z']
                self.R1 = np.sqrt(self.x1**2 + self.y1**2 +self.z1**2)
                self.vx1 = self.data1['vx']
		self.vy1 = self.data1['vy']
		self.vz1 = self.data1['vz']
                self.V1 = np.sqrt(self.vx1**2 + self.vy1**2 +self.vz1**2)
	        #print self.z #troubleshoot

                self.diskindex = np.where(self.data1['type'] == 2)
                
                # repeat above with second galaxy
		# use galaxy and snap to reconstruct string name
	        ilbl2 = '000' +str(snap)
		# remove all but three last digits
	        ilbl2 = ilbl2[-3:]
		self.filename2 = "./VLowRes/%s_"%(galaxy2) + ilbl2 + '.txt'
	        #print (self.filename) #troubleshooting
		# usual file reading
		self.time2, self.total2, self.data2 = Read(self.filename2)
		self.gname2 = galaxy2 


                if(self.IncludeCompanion is True):
                        # only do  this part if you actually need it 
		        # pull data values. no need to distinguish by particle type (yet)
		        self.m2 = self.data2['m']
		        self.x2 = self.data2['x']
		        self.y2 = self.data2['y']
		        self.z2 = self.data2['z']
                        self.R2 = np.sqrt(self.x2**2 + self.y2**2 +self.z2**2)
	                #print self.z #troubleshoot

                
	def MassEnclosed(self,radii,ptype):
		# function to calculate the mass enclosed out to varying radii for a certain particle type
		# inputs:
		# radii, array of radius magnitudes in kpc
		# ptype, particle type, integer
		# outputs: array of total masses within the radii specified from the center of mass, in Msun

		# need a COM object, to find COM position
		ME_COM = CenterOfMassMod.CenterOfMass(self.filename1, ptype) # COM found ONLY with one galaxy
		ME_COMP = ME_COM.COM_P(0.1,2)
	        #print ME_COMP #troubleshoot

                # initialize return array
                # length is an issue, since radii can be a single value
                # convert the single value to an array
                #print(radii) # troubleshooting
                try:
                        numberofradii = len(radii)
                except TypeError:
                        radii = np.array([radii])
                        numberofradii = len(radii)
		masses_enclosed = np.zeros(numberofradii)

                # give radii array proper units
                #radii *= u.kpc 

                # index to pick out particletype, find the radii of the particles in COM frame
		self.index = np.where(self.data1['type']==ptype)
		m1New = self.m1[self.index] 
		x1New = self.x1[self.index] - ME_COMP[0].value
		y1New = self.y1[self.index] - ME_COMP[1].value
		z1New = self.z1[self.index] - ME_COMP[2].value
                #print zNew[0] #troubleshoot
		R1New = np.sqrt(x1New**2 + y1New**2 + z1New**2)
                if (self.IncludeCompanion is True):
                        m2New = self.m2[self.index] 
		        x2New = self.x2[self.index] - ME_COMP[0].value
		        y2New = self.y2[self.index] - ME_COMP[1].value
		        z2New = self.z2[self.index] - ME_COMP[2].value
                        #print zNew[0] #troubleshoot
		        R2New = np.sqrt(x2New**2 + y2New**2 + z2New**2)

                # loop over the different radii 
		for i in range(numberofradii): 
		        MEnc = np.sum(m1New[np.where(R1New <= radii[i])])
                        if (self.IncludeCompanion is True):
                                MEnc += np.sum(m2New[np.where(R2New <= radii[i])])

			masses_enclosed[i] = MEnc

                return masses_enclosed

        def MassEnclosedTotal(self,radii):
                # function to calculate the total mass of all components out to varying radii
                # input is an array of radii in kpc
                # output is array of total ecnlosed masses at each radii in Msun

                # prepare for galaxies w/out a bulge (i.e. M33)
                if self.gname1 != "M33":
                        Mtot = self.MassEnclosed(radii,1)+ self.MassEnclosed(radii,2) + self.MassEnclosed(radii,3)
                else:
                        Mtot = self.MassEnclosed(radii,1)+ self.MassEnclosed(radii,2)
                return Mtot
        
        def HernquistProfile(self,r,a,Mhalo):
                # function that computes the mass within a radius using a theoretical profile
                # M = (Mhalo * r**2)/(a+r)**2
                # inputs:
                # radius in kpc
                # scale factor in kpc 
                # total halo mass Mhalo in Msun
                # returns enclosed mass predicted by hernquist profile
                return (Mhalo * r**2)/(a+r)**2
        
        def CircularVelocity(self,radii,ptype):
                # function to calculate circular speed of orbiters at a certain radius, assuming spherical symmetry
                # inputs:
                # particle type as integer
                # array of radii in kpc
                # returns array of velocities at each radius
                # Vcirc = sqrt(G*MEnc/r)
                Vcircs = np.sqrt(G*self.MassEnclosed(radii,ptype)/radii)
                return Vcircs
        
        def CircularVelocityTotal(self,radii):
                # function to find the circular speed of orbiters at a certain radius, assuming spherical symmetry
                # array of radii in kpc
                # returns array of velocities at each radius
                # Vcirc = sqrt(G*MEnc/r)
                VcircTotal = np.sqrt(G*self.MassEnclosedTotal(radii)/radii)
                return VcircTotal
        
        def HernquistVCirc(self,radii,a,Mhalo):
                # function to find the circular speed of orbiters at a certain radius, assuming spherical symmetry
                # inputs:
                # radii, array of radii in kpc
                # a scale radius in kpc
                # MHalo dark matter halo mass in Msun
                # returns array of velocities at each radius
                # Vcirc = sqrt(G*MEnc/r)
                return np.sqrt(G*self.HernquistProfile(radii,a,self.MassEnclosed(radii,1))/radii)
        
        def HernquistPotential(self,radii,a,M):
                # find the gravitational potential of the bulge/halo
                # inputs:
                # radii, an array of radii or a a single radius
                # a, scale length in kpc
                # M, mass of the component (bulge or halo)
                # phi = -GM/(r+a)
                return -G*M/(radii+a)
        
        def MiyamotoNagaiPotential(self,radii,a,z,Mdisk):
                # find gravitational potential of the disk
                # inputs:
                # radii, an array of radii or a a single radius
                # z, particle distance away from
                # a, scale length in kpc
                # M, mass of the component (bulge or halo)
                # phi = -GM/sqrt(R^2+(a+sqrt(b^2+z^2))^2)
                b = a/5.0
                B = a+np.sqrt((b**2)+(z**2))
                Denom = np.sqrt((radii**2) + (B**2))
                return -G*Mdisk/Denom
        
        def EVelocity(self,radii,a,z):
                # find the total escape velocity for all galactic components
                # inputs:
                # raddi an array of radii or a single radius
                # z, a particle's height
                # a, scale length
                # Mbulge,Mdisk,Mhalo, enclosed masses of the components of the galaxy
                DiskPhi = self.MiyamotoNagaiPotential(radii,a,z,self.MassEnclosed(radii,2))
                BulgePhi = self.HernquistPotential(radii,a,self.MassEnclosed(radii,3))
                HaloPhi = self.HernquistPotential(radii,a,self.MassEnclosed(radii,1))
                #print(HaloPhi) # troublshooting
                #print(self.MassEnclosed(radii,2)) # troublshooting


                # gravitational potentials obey the law of superposition
                return np.sqrt(2*np.abs((DiskPhi + BulgePhi + HaloPhi)))
        
        def JacobiRadius(self, R, Msat, Mhost):
                # calculate the jacobi radius, which we will say is the outer limit of the galaxy.
                # ideally tidal features will fall out of this radius
                # inputs:
                # R, distancce between the galaies
                # Msat, mass of orbiting galaxy (alternates in each case) 
                # MEnclosed, enclosed mass of the host galaxy (at radius of the satellite)


                return R*(Msat/(2*Mhost))**(1/3)

        """
        def JacobiOutliers(self,radii,m, MEnclosed):
                # function to find particles outside the jacobi radius annd store their kinematic values
                # inputs: 
                # takes an array of radii,
                # takes the mass of the particle
                # takes the mass enclosed at the particle's radius
                # returns mass, x, y, z
                ME = MassEnclosedTotal(self.R1)
                JRadii = JacobiRadius(self.R1,self.m1,ME) 
                self.outliers = np.where(self.R1 > JRadii)

                return [self.m1[self.outliers],self.x1[self.outliers],self.y1[self.outliers],\
                        self.z1[self.outliers]]
        def VelocityOutliers(self,radii,a,z,Mbulge, Mdisk, Mhalo):
                
                #function to get the 
        """

        
"""
#test prints
MW = MassProfile("MW",0)
r = np.arange(0.25,30.5,1.5)
ME = MW.MassEnclosed(r,1) 

print ME
"""
