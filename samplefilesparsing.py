#!/usr/bin/env python

fname = raw_input('Enter file with sample names: ')
#f1 = open('files.txt','r')
f1 = open('%s' %(fname),'r')
tb = f1.read()
tb = tb.split('\n')
tb = tb[:-1]
f1.close()


tb2 = []
for item in tb:
	tb2.append(item.split('_'))

TFdict = {}
for count in range(0,len(tb)):
	try:
		TFdict[tb2[count][0]].append(tb[count])
	except:
		TFdict[tb2[count][0]] = [tb[count]]

lrange = []
for key in TFdict:
	lrange.append(len(TFdict[key]))

lrange = list(set(lrange))

for item in lrange:
	fout = open('samplefiles%d.txt' %(item),'w')
	for key in TFdict:
		if len(TFdict[key]) == item:
			tmpline = '\t'.join(TFdict[key])
			fout.write('%s\n' %(tmpline))
	fout.close()

