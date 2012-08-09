#!/usr/bin/python
from pylab import *
from os import popen

def loadData(inputFile):
    data = []
    for line in inputFile:
        if "#" in line:
            continue
        data.append([float(v) for v in line.strip().split()])
    return data

def transpose(data):
        return [[data[j][i] for j in range(len(data))] for i in range(len(data[0]))]

def wcutoff():
   files = popen('ls ./output/diskTs_123_*').read().split() 
   print files

   for f in files:
	starfile = f.replace('disk','star')
	b = loadData(open(starfile))[-1][-1]
	print f
	switch = 0
       	data = loadData(open(f))
	for d in data:
	   if d[6] == 0:
		if switch == 0:
 		   text(d[0],d[1],'%s MG'%b, size = 10)
		   switch = 1
		d[1] = 1e-16
        dtrs = transpose(data)
       	loglog(dtrs[0],dtrs[1])
   xlim(1,1e8)
   xlabel('time / year',size=14)
   ylabel(r'$\dot{M}$ / $M_{\odot}\rm{year}^{-1}$',size=14)
   savefig('./paper/Mdot_evo_123.eps')
   show()

def wocutoff():
   files = popen('ls ./output/diskTs_133_*').read().split() 
   print files

   for f in [files[-1]]:
       	data = loadData(open(f))
        dtrs = transpose(data)
       	loglog(dtrs[0],dtrs[1],'r-')
   xlim(1,1e8)
   xlabel('time / year',size=14)
   ylabel(r'$\dot{M}$ / $M_{\odot}\rm{year}^{-1}$',size=14)
   savefig('./paper/Mdot_evo_130.eps')
   show()

if __name__ == '__main__':
   wcutoff()
   wocutoff()
