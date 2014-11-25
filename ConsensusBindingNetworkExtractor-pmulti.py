#!/usr/bin/env python
import csv
import math
import re
import itertools

fname1 = raw_input('Enter Number of Replicates: ')
file3 = open('samplefiles%s.txt' %(fname1),'r') # NonControlSamples.txt contains the names of control experiments
apn = raw_input('Enter name of previous MTBChIPpeakoverlay file:')
#apf = open('MTBChIPpeakoverlay6.txt','r')
apf = open('%s' %(apn),'r')
fout = open('MTB_ChIPbindingnetworkinfo%s.txt' %(fname1),'w')
fgout = open('MTBChIPpeakoverlay%s.txt' %(fname1),'w')

#pv = 0.05
pv = raw_input('Enter p-value significance threshold:')

AllPeakFreq = apf.read()
AllPeakFreq = AllPeakFreq.split('\n')
AllPeakFreq = map(int,AllPeakFreq[:-1])

tb = file3.read()
tb = tb.split('\n')
tbatches = []
for item in tb:
	x = item.split('\t')
	tbatches.append(x)
file3.close()

#### Functions used for the main part of the code #####

# Reads the information from each .csv file
def readExcelChIP(samplename):
	ft1 = open('./%s/%s.MTb.ChIPpeaks.Excel.csv' %(samplename,samplename),'r')
	reader1 = csv.reader(ft1,delimiter=',',quotechar='"')
	chiparray = []
	for row in reader1:
		chiparray.append(row) # chiparray is a list containing all of the peak information for each experiment
	
	expname = dir[0].split('_')
	dTF = expname[0]
	ft1.close()
	return dTF, chiparray

# Extracts the relevant peak properties from each ChIP experiment	
def extractPeakParts(chiparray):
	Pval = []
	Target = []
	VPM = []
	Score = []
	Fstart = []
	Rstop = []
	Fcenter = []
	Ccenter = []
	Rcenter = []
	DNAseq = []

	for peak in chiparray[1:]:
		try:
			VPM.append(float(peak[2]))
		except:
			break
		Pval.append(float(peak[3]))
		Target.append(peak[1])
		Score.append(float(peak[4]))
		Ccenter.append(float(peak[6]))
		Fcenter.append(float(peak[14]))
		Rcenter.append(float(peak[15]))
		DNAseq.append(peak[50].split('\n')[-1])

			
		Fst = peak[32] # the start of the forward strand peak (start of the footprint)
		Fst = int(math.floor(float(Fst)))
		Fstart.append(Fst)	
		del Fst
			
		Rop = peak[39] # the end of the reverse strand peak (end of the footprint)
		Rop = int(math.ceil(float(Rop)))
		Rstop.append(Rop)
		del Rop
	return Pval,Target,VPM,Score,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq

# Outputs the peak properties to text file
def outputGoodPeak(TF,Target,Pval,Score,VPM,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq,AllPeakFreq,count):
	fout.write('%s\t%s\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%s\t%d\n' %(TF,Target,Pval,Score,VPM,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq,count))	
	prange = range(Fstart,Rstop+1) 
	#Tally this peak base range in the variable AllPeakFreq
	for base in prange:
		AllPeakFreq[base-1] = AllPeakFreq[base-1] + 1
	return AllPeakFreq

# For peaks that show up in only one experiment, outputs the peak if p-value < 0.05
def SinglePeakProcessing(Tonly,TF,Target,Pval,Score,VPM,Fstart,Rstop,Fcenter,Ccenter,Rcenter,DNAseq,AllPeakFreq):
	for tar in Tonly:
		ix = Target.index(tar)
		if Pval[ix] < pv:
			AllPeakFreq = outputGoodPeak(TF,Target[ix],Pval[ix],Score[ix],VPM[ix],Fstart[ix],Rstop[ix],Fcenter[ix],Ccenter[ix],Rcenter[ix],DNAseq[ix],AllPeakFreq,1)
		
		del ix
	return AllPeakFreq

# For peaks that show up in multiple experiments, outputs the peak with the lowest p-value if P < 0.05
def DuplicatePeakProcessing(dir,dTF1,Target1,Pval1,Score1,VPM1,Fstart1,Rstop1,Fcenter1,Ccenter1,Rcenter1,DNAseq1,AllPeakFreq):
	# commonTargets contain target genes found in at least 2 of the experiments
	combos = list(itertools.combinations(range(len(dir)),r=2))
	commonTargets = [[] for item in range(len(combos))]
	#print dTF1
	TFdict = {}
	
	for item in range(len(combos)):
		combo0 = combos[item][0]
		combo1 = combos[item][1]
		commonTargets[item] = list(set(Target1[combo0]) & set(Target1[combo1]))

		noncombos = list(set(range(len(dir))) - set(combos[item]))

		for ctar in commonTargets[item]:
			ix = [0]*len(dir)
			ix[combo0] = Target1[combo0].index(ctar)
			ix[combo1] = Target1[combo1].index(ctar)
			
			pvlist = [100]*len(dir)
			pvlist[combo0]=Pval1[combo0][ix[combo0]]
			pvlist[combo1]=Pval1[combo1][ix[combo1]]

			for noncount in noncombos:
				if set([ctar]) & set(Target1[noncount]):
					ix[noncount] = Target1[noncount].index(ctar)
					pvlist[noncount]=Pval1[noncount][ix[noncount]]

			pmin = min(pvlist)
			
			if pmin < pv:
				pix = pvlist.index(pmin)
				#print pix
				#print len(dTF1[pix])
				kstring = '_'.join([dTF1[pix],Target1[pix][ix[pix]],str(Fstart1[pix][ix[pix]]),str(Rstop1[pix][ix[pix]])])
				astring = '\t'.join([dTF1[pix],Target1[pix][ix[pix]],str(Pval1[pix][ix[pix]]),str(Score1[pix][ix[pix]]),str(VPM1[pix][ix[pix]]),str(Fstart1[pix][ix[pix]]),str(Rstop1[pix][ix[pix]]),str(Fcenter1[pix][ix[pix]]),str(Ccenter1[pix][ix[pix]]),str(Rcenter1[pix][ix[pix]]),DNAseq1[pix][ix[pix]]])
				TFdict[kstring] = astring
	
	for keys in TFdict:
		peakproperties = TFdict[keys].split('\t')
		AllPeakFreq = outputGoodPeak(peakproperties[0],peakproperties[1],float(peakproperties[2]),float(peakproperties[3]),float(peakproperties[4]),int(peakproperties[5]),int(peakproperties[6]),float(peakproperties[7]),float(peakproperties[8]),float(peakproperties[9]),peakproperties[10],AllPeakFreq,2)
	return AllPeakFreq

######################################

for dir in tbatches: #dir is separate for each experiment
	tmp1 = []
	dTF1 = []
	Pval1 = []
	Target1 = []
	VPM1 = []
	Score1 = []
	Fstart1 = []
	Rstop1 = []
	Fcenter1 = []
	Ccenter1 = []
	Rcenter1 = []
	DNAseq1 = []

	for fnum in range(len(dir)):
		##print(dir[fnum])
		tmp1.append([])
		dTF1.append([])
		Pval1.append([])
		Target1.append([])
		VPM1.append([])
		Score1.append([])
		Fstart1.append([])
		Rstop1.append([])
		Fcenter1.append([])
		Ccenter1.append([])
		Rcenter1.append([])
		DNAseq1.append([])

		dTF1[fnum],tmp1[fnum] = readExcelChIP(dir[fnum])

		Pval1[fnum],Target1[fnum],VPM1[fnum],Score1[fnum],Fstart1[fnum],Rstop1[fnum],Fcenter1[fnum],Ccenter1[fnum],Rcenter1[fnum],DNAseq1[fnum]=extractPeakParts(tmp1[fnum])
	
	intersects = list(itertools.combinations(range(len(dir)),r=len(dir)-1))
	##print(len(intersects))
	# bigIntersect contains the union of target genes in all but one of the experiments
	bigIntersect = [[] for item in range(len(dir))]
	for item in range(len(dir)):
		for element in intersects[item]:
			bigIntersect[item] = list(set(bigIntersect[item]) | set(Target1[element]))
	
	# uniqueTargets contains target genes found in only one of the experiments
	uniqueTargets = [[] for item in range(len(dir))]
	for item in range(len(dir)):		
		uniqueTargets[item] = list(set(Target1[item]).difference(bigIntersect[item-1]))
		AllPeakFreq = SinglePeakProcessing(uniqueTargets[item],dTF1[item],Target1[item],Pval1[item],Score1[item],VPM1[item],Fstart1[item],Rstop1[item],Fcenter1[item],Ccenter1[item],Rcenter1[item],DNAseq1[item],AllPeakFreq)

	AllPeakFreq = DuplicatePeakProcessing(dir,dTF1,Target1,Pval1,Score1,VPM1,Fstart1,Rstop1,Fcenter1,Ccenter1,Rcenter1,DNAseq1,AllPeakFreq)
			
	del Fstart1
	del Rstop1
	del Target1
	del VPM1
	del Score1
	del Pval1
	del Fcenter1
	del Ccenter1
	del Rcenter1
	del DNAseq1
	del uniqueTargets	

for basefreq in AllPeakFreq:
	fgout.write('%d\n' %(basefreq))
fgout.close()