#!/usr/bin/python
from pylab import *
from os import popen

def loadData(inputFile):
    data = []
    for line in inputFile:
        if line.startswith("#"):
            continue
        data.append([float(v) for v in line.strip().split()])
    return data

def transpose(data):
        return [[data[j][i] for j in range(len(data))] for i in range(len(data[0]))]

def getSphData():
    data = transpose(loadData(open('./datums.txt')))
    mass_star = 1.1	#need to input these
    mass_disk = 0.3	#need to input these
    J_star = 0.
    J_disk = 0.
    J_total = 0.
    for i in range(len(data[1])):
	J_total = J_total + data[2][i]
	if data[1][i] < mass_star:
	   J_star = J_star + data[2][i]
	else:
	   J_disk = J_disk + data[2][i]
    #print J_star, J_disk, 'total: ', J_star + J_disk, J_total
    #plot(data[1],data[2],'ro')
    #show()
    return J_star, J_disk, mass_star, mass_disk

if __name__ == '__main__':
   J_star, J_disk, M_star, M_disk = getSphData()
   print '''
		*** Reading the out of the SPH data file *** 

	 '''
   Blist = [2e6,2e7,1e8,2e8,6e8]
   print '''
		*** writing the angular momentum data into data.in 
		*** and
		*** running spin_WD for B values of
	 ''', Blist, ' G'
   for B in Blist:
	spinfile = open('data.in','w')
   	line = '%1.3e %3.2e %3.4e %3.2e %3.4e'%(B,M_star,J_star,M_disk,J_disk)
   	spinfile.write(line)
   	spinfile.close()
   	popen('./spin_WD')
   	popen('mv star.out star_%1.1e.out'%B)
   	popen('mv disk.out disk_%1.1e.out'%B)
   print ''' finished '''
