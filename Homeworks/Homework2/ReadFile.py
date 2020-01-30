import numpy as np
import astropy.units as u


def Read(filename):
    #open file in read mode
    file = open(filename, 'r')
    #pull time from first line
    line1 = file.readline()
    label,value =  line1.split()
    time = float(value)*u.Myr
    #pull particle total from second line
    line2 = file.readline()
    label, value = line2.split()
    part_total = float(value)
    #close file when done
    file.close()
    #make arrays from data, ignore first 3 lines
    data = np.genfromtxt(filename, dtype=None,names=True, skip_header=3)
    #print(data['vx'][4]) #testing line, ensuredata is stored correctly
    #return information by tuple
    return time, part_total, data

#Read("MW_000.txt") #testing line, ensure function calls correctly 



