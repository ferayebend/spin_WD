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

def getI(mass):
	I_star = 3.21601*1e50*(mass)**0.34158*(1.e0-(mass/(1.437))**1.2499)**1.43773 
	return I_star

def getSphData(index):
    files = ['L_0305.dat','L_0408.dat',	'L_0606.dat', 'L_0608.dat']
    M_center = [0.5,0.8,0.6,0.8]
    M_stars = [0.62,0.92,1.10,1.08]
    data = transpose(loadData(open('/home/dagon/kulebi/repositories/applications/spin_WD/data/Baybars/%s'%files[index])))
    mass_star = M_stars[index]
    mass_center = M_center[index]
    mass_disk = data[1][-1]-mass_star#-0.02 #chandrasekhar kutlesine inmemesi icin uc kagit
    J_star = 0.
    J_disk = 0.
    J_total = 0.
    for i in range(len(data[1])):
	J_total = J_total + data[2][i]
	if data[1][i] < mass_center:
	   '''omega ile hesap'''
	   J_star = J_star + data[2][i]
	elif data[1][i] > mass_star:
	   J_disk = J_disk + data[2][i]
    #print J_star, J_disk, 'total: ', J_star + J_disk, J_total
    #plot(data[1],data[2],'ro')
    #show()
    J_star = getI(mass_star)/getI(mass_center)*J_star
    return J_star, J_disk, mass_star, mass_disk

def getSphInput(j):
    mass_star = [1.198,1.298,0.798]
    mass_disk = [2e-3,2e-3,2e-3]
    #J_star = 2.2471e49 #calculate from critical rotation
    omega = [1.5,2.0,0.4]
    J_star = getI(mass_star[j])*omega[j]
    J_disk = [2.6475e50,2.5e50,2.64e50]
    return J_star, J_disk[j], mass_star[j], mass_disk[j]

def runSpin(Blist,index):
   if index > 3:
   	J_star, J_disk, M_star, M_disk = getSphInput(1)#getSphData(index)
   else: 
   	J_star, J_disk, M_star, M_disk = getSphData(index)
   print '''
		*** Reading the out of the SPH data file *** 

	 '''

   print '''
		*** writing the angular momentum data into data.in 
		*** and
		*** running spin_WD for B values of
	 ''', Blist, ' G'
   popen('rm MRI_turnoff_vs_B_%i.dat'%(index))
   for B in Blist:
	spinfile = open('data.in','w')
   	line = '%1.3e %4.3e %3.4e %3.2e %3.4e'%(B,M_star,J_star,M_disk,J_disk)
   	spinfile.write(line)
   	spinfile.close()
   	#popen('./spin_WD_break')
	popen('./spin_WD_Ts02')
   	popen('mv star.out ./output/starTs_%i_%1.1e.out'%(index,B))
   	popen('mv disk.out ./output/diskTs_%i_%1.1e.out'%(index,B))
	popen('mv spectrum.out ./output/spectrumTs_%i_%1.1e.out'%(index,B))
	popen('more MRI_turnoff.out >> MRI_turnoff_vs_B_%i.dat'%(index))
   print ''' finished '''

def arrayplot(data,style):
   if data:
   	x = transpose(data)[0]
   	y = transpose(data)[1]
   	loglog(x,y,style,linewidth=1.7)

def convertMdot(Mdotgs):
   Msun = 1.989e33
   day = 8.64e4
   year = 365.25*day
   MdotMsunyr = Mdotgs*year/Msun
   return MdotMsunyr
   
def totalplot(Blist,T_final,index):
   #Blist = [2e6,2e7,1e8,2e8,6e8]
   #T_final = 30000.
   finaldata = []
   final = open('finaldataTs%i.txt'%index,'w')
   for B in Blist:
	data = loadData(open('./output/starTs_%i_%1.1e.out'%(index,B)))
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
	text(t_final, spin_final,'%4.2f MG'%B, size=10)
   savefig('period_evo_Bs_40kK_Ts%i.eps'%index)
   show()


if __name__ == '__main__':
   #Blist = [2e5,6e5,1.5e6,3e6]#,7e6,1e7,2e7,5e7,1.8e8,3e8,5e8,6e8] #0305 icin
   Blist = [1e5,6e5,2e6,5e6,1e7,5e7,9e7,2e8,5e8] #normalde kullanilan
   #Blist = [5e8]
   #Blist = [1e5,2e5,4e5,6e5,1e6,5e6,1e7,3e7,9e7] #0608 icin
   #files = []
   #files = ['L_0305.dat','L_0408.dat',	'L_0606.dat', 'L_0608.dat']
   runSpin(Blist,123)
   totalplot(Blist,40000,123)
