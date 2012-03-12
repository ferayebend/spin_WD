#!/usr/bin/python
from pylab import *
from os import popen

def loadData(inputFile):
    data = []
    for line in inputFile:
        if '#' in line:
            continue
        data.append([float(v) for v in line.strip().split()])
    return data

def transpose(data):
        return [[data[j][i] for j in range(len(data))] for i in range(len(data[0]))]

def getSphData():
    data = transpose(loadData(open('./data/datums.txt')))
    mass_star = 1.1	#need to input these
    mass_disk = 0.28	#need to input these
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

def runSpin(Blist):
   J_star, J_disk, M_star, M_disk = getSphData()
   print '''
		*** Reading the out of the SPH data file *** 

	 '''

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

def arrayplot(data,style):
   if data:
   	x = transpose(data)[0]
   	y = transpose(data)[1]
   	loglog(x,y,style,linewidth=1.7)

def totalplot(Blist,T_final):
   #Blist = [2e6,2e7,1e8,2e8,6e8]
   #T_final = 30000.
   finaldata = []
   final = open('finaldata.txt','w')
   for B in Blist:
	data = loadData(open('star_%1.1e.out'%B))
   	spin = []
   	spindown = []
	for d in data:
	    if d[9] > T_final and d[0] > 1.:	
	       if d[8] == 1:
		  spin.append([d[0],d[4]])
  	       elif d[8] == 0:
		  spindown.append([d[0],d[4]])
	    elif d[9] <= T_final and d[0] > 1:
		mass = d[1]
		B = d[10]
		spin_final = d[4]
		t_final = d[0]
		finaldata.append([mass, B, spin_final])
		final.write('%3.3f %3.2f %3.3f \n'%(mass, log10(B*1e6), log10(spin_final)))
		break
	arrayplot(spin,'r')
	arrayplot(spindown,'g')
	xlabel('time / year',size=14)
	ylabel('period / seconds',size=14)
	text(t_final, spin_final,'%4.1f MG'%B)
   savefig('period_evo_Bs.png')
   show()


if __name__ == '__main__':
   Blist = [3.5e5,6e5,1.2e6,1e7,2e7,8e7]
   #files = []
   runSpin(Blist)
   totalplot(Blist,30000)
