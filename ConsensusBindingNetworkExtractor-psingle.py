#!/usr/bin/env python
import csv
import math
import re

#fname = raw_input('Enter file with sample names: ')
#file1 = open('%s' %(fname),'r')

file1 = open('samplefiles1.txt','r') # NonControlSamples.txt contains the names of control experiments
sbatches = file1.read()
sbatches = sbatches.split('\n')
print len(sbatches)
file1.close()

AllPeakFreq = [0]*4411532

fout = open('MTB_ChIPbindingnetworkinfo1.txt','w')
fgout = open('MTBChIPpeakoverlay1.txt','w')

pv = raw_input('Enter p-value significance threshold:')

for dir in sbatches:
	tmp = []
	print dir
	ft = open('./%s/%s.MTb.ChIPpeaks.Excel.csv' %(dir,dir),'r')
	reader = csv.reader(ft,delimiter=',',quotechar='"')
	
	x = dir.split('_')
	sTF = x[0]
	
	for row in reader:
		tmp.append(row) # tmp is a list containing all of the peak information for each experiment
	# tmp2 = tmp[1:]	# want to get rid of the first row containing row headers
	ft.close()
	
	for peak in tmp[1:]:
		Pval = float(peak[3])
		
		if Pval < pv:
			Target = peak[1]
			VPM = float(peak[2])
			Score = float(peak[4])
			Ccenter = float(peak[6])
			Fcenter = float(peak[14])
			Rcenter = float(peak[15])

			DNAseq = peak[50].split('\n')[-1]

			Fstart = peak[32]
			#if re.search('NA',Fstart):
			#	continue
			#else:
			Fstart = int(math.floor(float(Fstart)))
			#Fstart = int(math.floor(float(peak([32])))
			
			Rstop = peak[39]
			## Convert Rstop to an integer if it's not "NA"
			#if re.search('NA',Rstop):
			#	continue
			#else:
			Rstop = int(math.ceil(float(Rstop)))
			#Rstop = int(math.ceil(float(peak[39])))
			
			# for each peak in each experiment, find the span of the peak (from Fstart to Rstop)
			# Convert Fstart to an integer if it's not "NA"
			#Find all of the bases spanned by the peak
			prange = range(Fstart,Rstop+1) 
			
			#Tally this peak base range in the variable AllPeakFreq
			for base in prange:
				AllPeakFreq[base-1] = AllPeakFreq[base-1] + 1
			
			fout.write('%s\t%s\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\n' %(sTF,Target,Pval,Score,VPM,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq,1))
			
			#Clear the variables used
			del Fstart
			del Rstop
			del Target
			del VPM
			del Score
			del prange
			del Ccenter
			del Fcenter
			del Rcenter
			del DNAseq
			
		del Pval
	
	del reader
	del x
	del tmp
	del sTF
	

for basefreq in AllPeakFreq:
	fgout.write('%d\n' %(basefreq))
fgout.close()