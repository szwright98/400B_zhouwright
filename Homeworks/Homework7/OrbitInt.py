import numpy as np

G = 4.498768e-6 #in units of kpc, Msun, Gyr

class M33AnalyticOrbit:

    # class that will contain a series of functions to determine acceleration of M33, and predict is position and velocity

    def __init__(self,fileout):
        # initialize class, with the name of the output file.
        self.fileout = fileout
        # define starting positions, velocities in kpc and km/s with respect to M31 center of mass (HW4)

        self.x = 98.56
        self.y = 119.99
        self.z = 127.76

        self.vx = 28.43
        self.vy = -173.92
        self.vz = -93.23
        

        # get scale lengths and masses of M31 components in kpc and Msun

        self.rdisk = 5. # from HW7 pdf
        self.Mdisk = 1.2e11 # from HW3 solutions
        self.rbulge = 1. # from HW7 pdf
        self.Mbulge = 1.9e10 # from HW3 solutions
        self.rhalo = 0.5 # from HW5 
        self.Mhalo = 1.921e12 # from HW3 solutions

    # function to determine accelaration as predicted by the Hernquist profile (for approximately spherical distributions)
    # a = -(GM/(r(r_a+r)^2)) * r vector 
    def HernquistAccel(self, M, r_a, x, y, z):
        # inputs:
        # M, mass of the halo or bulge in Msun
        # r_a, scale length of component in kpc
        # x,y,z, current position of the satellite wrt center of mass, in kpc
        # returns acceleration vector

        # create tuple to store the position as a vector
        rvec = [x,y,z]
        r = np.sqrt(x**2+y**2+z**2)
        #print r # troubleshooting
        # calculate coeffecient out front
        C = -(G*M)/(r*(r_a+r)**2)
        #print np.multiply(C,rvec) # troubleshooting
        return np.multiply(C,rvec)

    # function to calculate the acceleration due to the gravity of the disk using the Miyamoto Nagai profile
    # potential = -GM/sqrt(R^2 +B^2)
    # R = sqrt(x^2 + y^2)
    # B = r_d + sqrt(z^2+z_d^2)
    # r_d is scale length of the disk, z_d is r_d/5.0 according to pdf
    # actual acceleration vector has slightly different vectors for each component, see HW7
    def MiyamotoNagaiAccel(self,M,rd,x,y,z):
        # inputs:
        # M, mass of the disk in Msun
        # rd, scale length of disk in kpc
        # x,y,z, current position of the satellite wrt center of mass, in kpc
        # returns acceleration vector
        R = np.sqrt(x**2 + y**2)
        #print R #troubleshooting
        zd = rd/5.0
        B = rd + np.sqrt(z**2 + zd**2)
        #print B #troubleshooting
        #print z #troubleshooting
        rvec = [x,y,z]
        avec = [0,0,0]
        #print G * M #troubleshooting
        top = -(G*M)
        D = (R**2+B**2)**1.5
        #print D # troubleshooting

        E = top/D
        #pick out z component, as it has a different form
        for i in range(len(rvec)):
            if i == 2:
                E *= B / np.sqrt(z**2+zd**2)
                #print E #troubleshooting
                avec[i] = E*rvec[i]
            else:
                #print E #troubleshooting
                avec[i] = E*rvec[i]
        #print avec # troubleshooting
        return avec

    #function that returns the summed acceleration vector
    def M31Accel(self,x,y,z):
        # x,y,z, current position of the satellite wrt center of mass, in kpc
        # returns summed acceleration vector
        DiskAcc = self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,x,y,z)
        HaloAcc = self.HernquistAccel(self.Mhalo, self.rhalo,x,y,z)
        BulgeAcc = self.HernquistAccel(self.Mbulge,self.rbulge,x,y,z)

        TotalAcc = np.add(DiskAcc,HaloAcc)
        TotalAcc = np.add(TotalAcc,BulgeAcc)
        #print TotalAcc #troubleshooting
        return TotalAcc


    # a function to take a step in an integration algorithm to calculate the orbits
    #uses leapfrof integration scheme
    def LeapFrog(self,dt,x,y,z,vx,vy,vz):
        # inputs:
        # dt, time integration step length
        # x,y,z, current position of the satellite wrt center of mass, in kpc
        # vx,vy,vz, current velocities of satellite
        # returns new positions and velocities

        rinit = [x,y,z]
        vinit = [vx,vy,vz]
        halfstep = dt/2.
        #predict 3d position vector at middle of the timestep
       
        vhalf = np.multiply(vinit,halfstep)
        rhalf = np.add(rinit,vhalf)
        #advance position and velocity by full timestep
        vstep = np.multiply(dt,self.M31Accel(rhalf[0],rhalf[1],rhalf[2]))
        vnext = np.add(vinit,vstep)
        rstep = np.multiply(vnext,halfstep)
        rnext = np.add(rhalf,rstep)

        #print vnext #troubleshooting
        return rnext,vnext

    # a function to actually run the integration
    def OrbitIntegrator(self,t_o,dt,tmax):
        # inputs:
        # starting time, t_o
        # time step dt
        # maximum time tmax
        # returns an array containing the values of velocity and position for M33

        #initialize array, store initial values
        orbit = np.zeros((1,7))
        orbit[0][0] = t_o
        orbit[0][1] = self.x
        orbit[0][2] = self.y
        orbit[0][3] = self.z
        orbit[0][4] = self.vx
        orbit[0][5] = self.vy
        orbit[0][6] = self.vz
        
        #initialize counter variables
        t = t_o
        i = 0
        
        while t < tmax:
            print(t)
            pos, vel = self.LeapFrog(dt,orbit[i][1],orbit[i][2],orbit[i][3],orbit[i][4],orbit[i][5],orbit[i][6])
            t += dt
            orbit = np.vstack((orbit, [t,pos[0],pos[1],pos[2],vel[0],vel[1],vel[2]]))
            i+=1
        print("Ding!")
        np.savetxt(self.fileout, orbit, fmt = "%11.3f"*7, comments='#',
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        return None

# using the function 

M33 = M33AnalyticOrbit('M33_Analytical.txt')

M33.OrbitIntegrator(0,.01,12)
