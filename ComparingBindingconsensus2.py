#!/usr/bin/env python

import itertools

## This code counts for each peak, how many other peaks in the network set overlap at least 50% 
fname = raw_input('Enter network file: ')
#file1 = open('MTBchipnetwork031814.txt','r')
file1 = open('%s' %(fname),'r')

bindWhole = file1.read()
file1.close()

bindHalf = bindWhole.split('\n')
bindList = []
for item in bindHalf:
	x = item.split('\t')
	bindList.append(x)
	
bindList = bindList[:-1]

fout = open('ChIPpeakOverlapCounts031814.txt','w')
peakOverlapCount = [0] * len(bindList)

#ran1 = range(0,len(bindList))

for pair in itertools.combinations(range(len(bindList)),r=2):
	s1 = int(bindList[pair[0]][5])
	e1 = int(bindList[pair[0]][6])
	r1 = range(s1,e1+1)

	s2 = int(bindList[pair[1]][5])
	e2 = int(bindList[pair[1]][6])
	r2 = range(s2,e2+1)

	overlap = len(set(r1).intersection(set(r2)))
	if overlap > 0.5*len(r1):
		peakOverlapCount[pair[0]] = peakOverlapCount[pair[0]] + 1
		peakOverlapCount[pair[1]] = peakOverlapCount[pair[1]] + 1

for ix1 in range(len(bindList)):
	fout.write('%s\t%s\t%d\t%d\t%d\n' %(bindList[ix1][0],bindList[ix1][1],bindList[ix1][5],bindList[ix1][6],peakOverlapCount[ix1]))