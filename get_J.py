#!/usr/bin/python
from pylab import *
from scipy import interpolate

def loadData(inputFile):
    data = []
    for line in inputFile:
        if line.startswith("#"):
            continue
        data.append([float(v) for v in line.strip().split()])
    return data

def transpose(data):
        return [[data[j][i] for j in range(len(data))] for i in range(len(data[0]))]

if __name__ == '__main__':
    data = transpose(loadData(open('datums.txt')))
    mass_star = 1.1
    J_star = 0.
    J_disk = 0.
    J_total = 0.
    for i in range(len(data[1])):
	J_total = J_total + data[2][i]
	if data[1][i] < mass_star:
	   J_star = J_star + data[2][i]
	else:
	   J_disk = J_disk + data[2][i]
    print J_star, J_disk, 'total: ', J_star + J_disk, J_total
    plot(data[1],data[2],'ro')
    show()
