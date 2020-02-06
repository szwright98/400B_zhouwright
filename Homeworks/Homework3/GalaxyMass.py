import numpy as np
import astropy.units as u
from astropy.table import Table
from ReadFile import Read
from astropy.io import ascii
#Function to return the total mass of a galaxy component
def ComponentMass(filename,particletype):
	#inputs:
		#particletype is the location or type or particle in the galaxy(Halo, disk, bulge)
		#filename is the name of the file that contains the data
	#returns the summed mass of all the particles of a given type in 1e12 Msun (converted from 1e10 Msun)
	#read in file 
	time, total, data = Read(filename)
	#store indeices of a certain particle type
	index= np.where(data['type'] == particletype)
	#create a new array with the masses of each particle of the given type
	#mass is stored in units of 1e10 Msun
	mdata = data['m'][index]*1e10*u.Msun
	#get total mass
	mtotal = np.round(np.sum(mdata))
	#perform conversion
	mtotal.to(1e12*u.Msun)
	return mtotal
#call function to get values
MW_halo = ComponentMass("MW_000.txt",1)
MW_disk = ComponentMass("MW_000.txt",2)
MW_bulge = ComponentMass("MW_000.txt",3)
MW_total = MW_disk+MW_halo+MW_bulge
MW_fbar = (MW_disk+MW_bulge)/MW_total

M33_halo = ComponentMass("M33_000.txt",1)
M33_disk = ComponentMass("M33_000.txt",2)
M33_bulge = ComponentMass("M33_000.txt",3)
M33_total = M33_disk+M33_halo+M33_bulge
M33_fbar = (M33_disk+M33_bulge)/M33_total

M31_halo = ComponentMass("M31_000.txt",1)
M31_disk = ComponentMass("M31_000.txt",2)
M31_bulge = ComponentMass("M31_000.txt",3)
M31_total = M31_disk+M31_halo+M31_bulge
M31_fbar = (M31_disk+M31_bulge)/M31_total
#declare Qtable
t = Table(names = ('Galaxy Name', 'Halo Mass(Msun)','Disk Mass(Msun)','Bulge Mass(Msun)','Total(Msun)','fbar'), dtype=('str','f','f','f','f','f'))
#fill Qtable with function calls
t.add_row(('Milky Way',MW_halo,MW_disk,MW_bulge,MW_total,MW_fbar))
t.add_row(('M31',M31_halo,M31_disk,M31_bulge,M31_total,M31_fbar))
t.add_row(('M33',M33_halo,M33_disk,M33_bulge,M33_total,M33_fbar))
#export table in a format that is readable for LaTeX editors
ascii.write(table=t, output='fin_table', format='latex')

