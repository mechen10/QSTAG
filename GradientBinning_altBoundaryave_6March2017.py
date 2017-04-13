# BIOL548CD- Melissa Chen
# Final project

# In this script, I take the taxasummaries (output from QIIME) and metadata with salinity (manually created) and make two graphs:
	# 1. 'Distribution of Boundaries Across Salinity'
	# 2. 'Composition of Specialists across Salinity'
# Additionally, there are four text outputs:
	# 1. boundaries.txt (a list of boundaries)
	# 2. modelBoundaries_type.txt (list of taxa,type, and boundaires)
	# 3. taxa_abundances_across_salinity.txt (table with salinity as headers and abundances of each taxa)
	# 4. types_across_salinity.txt (table with salinity as headers and abundances of each 'type' at each salinity bin)
	
# What this script does is:
	# - For each OTU in the taxasummaries, decides whether it is a Low,Inter, or High specialist (or none of them)
		# It does this by comparing the mean relative abundance of three groups (groupA, groupB, groupC),
		# The three groups are separated by boundaries A and B.
		# The script iterates through all combinations of A and B and finds the combination that best 'fits' the data
		# Best 'fit' is determined by calculated squared differences between means and actual values
		# Then, it uses Welch's T-test to see if there are significant differences in the means between groupA, groupB, and groupC
		# eg. if groupA > groupC (significant), then it is classified as Low.
		# It uses A and B to calculate the places were each OTU seems to be 'turning over' (ie, changing in abundance significantly) 
	# - Then, it plots the distribution of ALL 'boundaries' (A and B) across salinity for all taxa
	# - Makes file outputs 
	# - Finally, it takes the classifications of each OTU and creates a stacked bar graph that shows the relative abundance of each 'type'.
	
# REQUIREMENTS:
	# You must have:
		# graphics
		# math
		# scipy
		# argparse
	# Also:
		# taxasummaries.txt file somewhere
		# metadata.txt file somewhere
		

#==========================================

import math 
import numpy # For standard deviation
	# Required for BINNING and BINNING SALINITY sections
from scipy import stats # need for welch's t-test
	# Required for Welch's t-test
# from graphics import *
import argparse
import os
import subprocess

#==========================================

# BEGIN

#==========================================
# FUNCTION TO LOAD FILES AND GET CORRECT FORMAT

def makeTaxaSummaries(taxaTablePWD,metadataPWD):
	global metadata_name
	print "Making and loading taxa table and metadata table..."
	os.system('biom convert -i ' + taxaTablePWD + ' --to-tsv --header-key taxonomy --table-type="OTU table" -o ./OTUTableText.txt')
# 	subprocess.call(['biom', 'convert', str('-i'+taxaTablePWD),'--to-tsv', '--header-key','taxonomy', '-o./OTUTableText.txt', '--table-type="OTU table"'])
	# Open up and read all lines in
	taxaOpen = open('./OTUTableText.txt', 'r')
	taxaOpenTemp = []
	for i in taxaOpen:	
		taxaOpenTemp += [i]
	taxaOpen.close()
	# Split each line into lists by stripping ends and splitting by 'tab'
	tempList =[]
	for j in taxaOpenTemp:
		tempLine = j.strip()
		tempList += [tempLine.split('\t')]
	del tempList[0]
	# Sort information in relavent piles
	taxaIDs = {}
	taxaCountsTemp = {}
	first = True
	for y in tempList:
		if y[0] == '#OTU ID':
			sites = y[1:len(y)-1]
		else:
			# Make dictionary of taxaIDs for later
			taxaIDs[y[0]] = y[len(y)-1]
			# Make file of 'total' counts for each site
			for x in range(len(y)):
				if (x != 0) and (x != (len(y)-1)):
					if first:
						taxaCountsTemp[str(x)] = float(y[x])
					else:
						taxaCountsTemp[str(x)] += float(y[x])
			first = False
	taxaTable = {}
	for i in range(len(tempList)):
		taxaTable[tempList[i][0]] = [[]]
		if tempList[i][0] == '#OTU ID':
			pass
		else:
			for j in range(1,len(tempList[i])-1):
				# Sum of all abundances to make relative abundances
				sumAbund = int(taxaCountsTemp[str(j)])
				value = tempList[i][j]
				taxaTable[tempList[i][0]][0].append(float(value)/float(sumAbund))
				# Save values as relative abundances instead of absolute ones
# 				taxaTable[tempList[i][0]] = [[float(x)/sumAbund for x in values]]
	metadataOpen = open(metadataPWD, 'r')
	metadataOpenTemp = []
	for i in metadataOpen:
		metadataOpenTemp += [i]
	metadataOpen.close()
	tempMeta =[]
	for j in metadataOpenTemp:
		tempLine = j.strip()
		tempMeta += [tempLine.split('\t')]
	positionSal = tempMeta[0].index(metadata_name)
	metadata = []
	for line in tempMeta:
		metadata.append([line[0],line[positionSal]])
	# Now, change key names so they're not longer sites; they're salinities
	for site in metadata:
		sites = [site[1] if x == site[0] else x for x in sites]
	# Make proper format, but with site names instead of numbers
	for taxa in taxaTable:
		taxaTable[taxa].append(sites)
	# Make abundance values integer as well
	for x in taxaTable:
		taxaTable[x] = [[float(strings) for strings in keys] for keys in taxaTable[x]]
	return taxaTable,taxaIDs

#==========================================
# Delete certain taxa based on total abundance; must be sufficiently abundant in order to work

def deleteLowAbund(taxasummary): # Change to absolute threshold, or get rid of entirely. Get rid of show ups in < 3 samples
	newTaxaSummary = {}
	for taxa in taxasummary:
		nonzeroCount = 0
		for i in taxasummary[taxa][0]:
			if i > 0.0:
				nonzeroCount += 1
		if nonzeroCount >= 3:
			newTaxaSummary[taxa] = taxasummary[taxa]
	return newTaxaSummary

#==========================================
# Msc funcitons

def average(listValues): # shortcut for 'average'
	# finds average for list of numbers
	if len(listValues) == 0:
		return None
	else:
		average = float(sum(listValues))/len(listValues)
		return average
		
def is_numeric(s):
	try:
		float(s)
		return True
	except ValueError:
		return False
		
#==========================================
# MAKING BINS FOR EACH SET OF BOUNDARIES-- USED BELOW

def sortBins(A, B, listAbund): # Sorts list of abundances into bins according to A and B and their associated salinities
	lista = []
	listb = []
	listc = []
	for i in range(len(listAbund[1])):
		if listAbund[1][i] < A:  # if A = B, value will be sorted into A
			lista.append(listAbund[0][i])
		if listAbund[1][i] >= A and listAbund[1][i] <= B:
			listb.append(listAbund[0][i])
		if listAbund[1][i] > B:
			listc.append(listAbund[0][i])
	binValues = {'lista':lista, 'listb':listb, 'listc':listc}
	return binValues # output is a dictionary with abundance observations sorted into lists according to what salinity they were found in

# def sortBinsFM(A,B,listAbund): # Same sorting abundances, except with one boundary only
# 	lista = []
# 	listb = []
# 	listc = []
# 	bValueGrad = []
# 	for i in range(len(listAbund[1])):
# 		if listAbund[1][i] < A:  # if A = B, value will be sorted into A
# 			lista.append(listAbund[0][i])
# 		if listAbund[1][i] >= A and listAbund[1][i] <= B:
# 			listb.append(listAbund[0][i])
# 			bValueGrad.append(listAbund[1][i])
# 		if listAbund[1][i] > B:
# 			listc.append(listAbund[0][i])
# 	binValues = {'lista':lista, 'listb':listb, 'listc':listc}
# 	return binValues,bValueGrad
#==========================================

# MAKE AND TEST A MODEL USING ERROR SQUARED


def makeAModel(A,B, binValues):
	aMean = average(binValues['lista'])
	bMean = average(binValues['listb'])
	cMean = average(binValues['listc'])
	# Make a combined 'piecewise' list that has each value and its corresponding mean
	abundance = binValues['lista'] + binValues['listb'] + binValues['listc']
	means = [aMean]*len(binValues['lista']) + [bMean]*len(binValues['listb']) + [cMean]*len(binValues['listc'])
	combined = {'abundance': abundance, 'means': means}
	return combined
	
# def bLineEquation(A,B,aMean,cMean,bValueGrad):
# 	m =(float(aMean) - cMean)/(float(A)-B) # Calculate slope
# 	b = aMean-m*A # Calculate y intercept
# 	bLine = []
# 	for i in bValueGrad:
# 		bLine.append(m*i+b) # Calculate y for equation
# 	return bLine,m,b

# def makeAModelFM(A,B,binValues,bValueGrad):
# 	aMean = average(binValues['lista'])
# 	cMean = average(binValues['listc'])
# 	bLine,m,b = bLineEquation(A,B,aMean,cMean,bValueGrad)
# 	# Make a combined 'piecewise' list that has each value and its corresponding mean
# 	abundance = binValues['lista'] + binValues['listb'] + binValues['listc']
# 	means = [aMean]*len(binValues['lista']) + bLine + [cMean]*len(binValues['listc'])
# 	combined = {'abundance': abundance, 'means': means}
# 	return combined
	
def meanErrorSquared(combined):
	# find the sum of squared difference for each value from the mean in that bin
	errorSquared = 0
	n = range(len(combined['abundance']))
	for i in n:
		error = combined['means'][i] - combined['abundance'][i]
		sqerror = error**2
		errorSquared += sqerror
	MES = float(errorSquared)/len(combined['abundance'])
	return MES

# def scaleBySig(A,B,pAB,pBC): # TBD: Should fix this
	# This function scales the location of the boundary if we have 2 boundaries that need to be combined into one. 
	# The way I chose to do this is by square rooting the inverse of the significance, such that a greater significance will yield a higher number
	# Then, I find the relative 'significance' of these two numbers. This relative significance is used to scale the final boundary position.
# 	adjpAB = (1/pAB)**0.5
# 	adjpBC = (1/pBC)**0.5
# 	totalSig = adjpAB + adjpBC
# 	aScale = adjpAB/totalSig
# 	bScale = adjpBC/totalSig
# 	scaledBoundary = round((A*aScale + B*bScale),1)
# 	scaledBoundary = average([A,B])
# 	return scaledBoundary # Output is the scaled Boundary based on significance
	
def scaleByDiff(A,B,meanA,meanB,meanC):
	if meanA > meanC:
		diffA = meanA-meanB
		diffB = meanB-meanC
		diffAB = B-A
		if diffA <= 0:
			finalBoundary = B
		elif diffB <= 0:
			finalBoundary = A
		else:
			scaleFactor = diffA/(diffA+diffB)
			diffscaled = diffAB*scaleFactor
			finalBoundary = B-diffscaled
	if meanC > meanA:
		diffA = meanB-meanA
		diffB = meanC-meanB
		diffAB = B-A
		if diffB <= 0:
			finalBoundary = A
		elif diffA <= 0:
			finalBoundary = B
		else:
			scaleFactor = diffB/(diffA+diffB)
			diffscaled = diffAB*scaleFactor
			finalBoundary = A+diffscaled
	if meanA == meanC:
		print "ERROR"
	return finalBoundary
# def typeTaxaForM(A,B,listAbund):
# 	global threshold
# 	global Low,Inter,High
# 	binValues,bValueGrad = sortBinsFM(A,B,listAbund)
# 	groupA = binValues['lista']
# 	groupB = binValues['listb']
# 	groupC = binValues['listc']
# 	meanA = average(groupA)
# 	meanC = average(groupC)
# 	lineB,m,b = bLineEquation(A,B,meanA,meanC,bValueGrad)
# 	# Set threshold
# 	if threshold[0] == True:
# 		maxAbund = max(groupA + groupB + groupC)
# 		thresh = maxAbund*threshold[1]
# 	else:
# 		thresh = threshold[1]
# 	pAC = stats.ttest_ind(groupA,groupC, equal_var = False)[1]
# 	pvalue = 0.05
# 	sigAC = pAC < pvalue
# 	bloom = None
# 	noclassTF = False
# 	if meanA > meanC and sigAC:
# 		if meanC < thresh:
# 			typeFM = Low
# 			boundaries = [average([A,B])]
# 		else:
# 			typeFM = 'half'+Low
# 			boundaries = [average([A,B])]
# 	if meanA < meanC and sigAC:
# 		if meanA < thresh:
# 			typeFM = High
# 			boundaries = [average([A,B])]
# 		else:
# 			typeFM = 'half'+High
# 			boundaries = [average([A,B])]
# 	if meanA == meanC or not sigAC:
# 		typeFM = 'noclass'
# 		Avalue = None
# 		boundaries = []
# 		bloom = 'TBA'
# 		noclassTF = True
# 	return typeFM,boundaries,meanA,meanC,[m,b],pAC,bloom, noclassTF

	# 
# def notBrackishIteration(listAbund):
# 	global minA,maxB,ABdiff
# 	global threshold
# 	bestfitAB = []
# 	modelDiffList = []
# 	currentbest = 0
# 	for i in range(minA,(maxB-ABdiff)):
# 		for j in range((i+ABdiff),maxB):
# 			A = i
# 			B = j
# 			binValues,bValueGrad = sortBinsFM(A,B,listAbund)
# 			if len(binValues['lista']) <= 1 or len(binValues['listb']) <= 1 or len(binValues['listc']) <= 1: # if the bins have nothing in them, then don't use that bin combination
# 				pass
# 			else:
# 				combined = makeAModelFM(A,B,binValues,bValueGrad)
# 				MES = meanErrorSquared(combined)
# 				if MES == 0:
# 					pass
# 				else:
# 					invMES = 1/MES
# # 					binValues,bValueGrad = sortBinsFM(A,B,listAbund)
# # 					combined = makeAModelFM(A,B,binValues,bValueGrad)
# 					modelDiff = max(combined['means'])-min(combined['means'])
# 					if invMES > currentbest:
# 						bestfitAB = [[A,B]]
# 						modelDiffList = [modelDiff]
# 						currentbest = invMES
# 					elif invMES == currentbest:
# 						bestfitAB.append([A,B]) # Assumes they are due to some lack of data, not that they are exactly the same
# 						modelDiffList.append([modelDiff])
# # 	if len(bestfitAB) > 1: # Is there more than 1 boundary?
# # 		differencesA = numpy.diff(bestfitAB) # Check if it is consecutive or not
# # 	else:
# # 		differencesA = [False]
# # 	if (True in [True for i in differencesA if i>1]): # if consecutive numbers, then just average them. if not, then print warning but average them anyway.
# # 		print "WARNING: TWO OR MORE IDENTICAL MODELS-ForM"
# # # 		print bestfitA
# 	maxDifferenceFound = [modelDiffList.index(i) for i in modelDiffList if i == max(modelDiffList)]
# 	if len(maxDifferenceFound) > 1:
# 		print "WARNING: SAME MODELDIFF FOR IDENTICAL MODELS"
# 	inputA = float(bestfitAB[maxDifferenceFound[0]][0])# input into typeTaxaForM needs to be float or int, not list
# 	inputB = float(bestfitAB[maxDifferenceFound[0]][1])
# # 	inputA = float(Avalue[0]) # input into typeTaxaForM needs to be float or int, not list
# 	typeFM,boundaries,meanA,meanC,bEquation,pAC,bloom,noclassTF = typeTaxaForM(inputA,inputB,listAbund)
# 	return typeFM,boundaries,pAC,meanA,meanC,bEquation,inputA,inputB,bloom,noclassTF
# 	

def typeTaxa(A, B, listAbund): # Uses Welch's t-test and bins above to classify taxa as marine, fresh, or brackish (or other)
	global threshold
	global Low,Inter,High
	typeOutput = {'boundaries':[], 'type': '', 'typeSimple':'', 'meanA': [], 'meanB': [], 'meanC': [], 'A': [], 'B': [], 'sigAB': [], 'sigBC': [], 'sigAC': [], 'bloom': 'No'} # Empty dictionary
	binValues = sortBins(A,B,listAbund) # Use function above to create dictionary with abundance and salinity information
	# Find out whether someting is marine-,fresh-, or brackish- specific
	groupA = binValues['lista']
	groupB = binValues['listb'] # might be 1 in length
	groupC = binValues['listc']
	meanA = average(groupA)
	meanB = average(groupB) # might be 'None'
	meanC = average(groupC)
	# Find threshold by using proportion of max, if necessariy
	if threshold[0] == True:
		maxAbund = max(groupA + groupB + groupC)
		thresh = maxAbund*threshold[1]
	else:
		thresh = threshold[1]
# 	typeOutput['meanA'] = meanA
# 	typeOutput['meanB'] = meanB
# 	typeOutput['meanC'] = meanC
	bonef = 0.05/2 # boneferri correction; for each taxa, we are repeatedly comparing means (comparing twice; AB-BC OR AB-AC OR AC-BC)
	# Calculate significance
# 	if [True for x in [groupA,groupB,groupC] if numpy.var(x) == 0]:
# 		print "NO VARIANCE IN ONE OF GROUPS"
	# Calculate variance, but first test if each combination of groups has 0 variance.
	# If variance is 0, set p** to 1, which is maximum
	if average(groupA) == average(groupC) and numpy.var(groupA) == 0 and numpy.var(groupC) == 0:
		pAC = 1
	else:
		pAC = stats.ttest_ind(groupA,groupC, equal_var = False)[1] # p-value of A vs C
	if average(groupA) == average(groupB) and numpy.var(groupA) == 0 and numpy.var(groupB) == 0:
		pAB = 1
	else:
		pAB = stats.ttest_ind(groupA,groupB, equal_var = False)[1] # p-value of A vs B
	if average(groupB) == average(groupC) and numpy.var(groupB) == 0 and numpy.var(groupC) == 0:
		pBC = 1
	else:
		pBC = stats.ttest_ind(groupB,groupC, equal_var = False)[1] # p-value of B vs C
# 	print pAC,pAB,pBC
	sigAB = pAB < bonef # True if significant
	sigBC = pBC < bonef # True if significant
	sigAC = pAC < bonef # True if significant
# 	typeOutput['sigAB'] = pAB # Save in output file
# 	typeOutput['sigBC'] = pBC # Save in output file
# 	typeOutput['sigAC'] = pAC # Save in output file
	Stda = numpy.std(groupA) # For 'bloom' test
	Stdb = numpy.std(groupB)
	Stdc = numpy.std(groupC)
	typeOutput['sigAB'] = pAB
	typeOutput['sigBC'] = pBC
	typeOutput['sigAC'] = pAC
	typeOutput['meanA'] = meanA
	typeOutput['meanB'] = meanB
	typeOutput['meanC'] = meanC
	typeOutput['A'] = A
	typeOutput['B'] = B
	isBrack = False	# See if there is a brackish community or not; if there isn't, then I compare just A and C.
	if sigAB and sigBC: # When the middle group is REAL (and not just 1 number), and it is significantly different than both flanking groups
		if meanB > meanA and meanB > meanC: # brackish water species
			isBrack = True
			typeOutput['boundaries'] = [A,B]
			typeOutput['sigAB'] = pAB
			typeOutput['sigBC'] = pBC
			typeOutput['sigAC'] = pBC 
 			typeOutput['meanA'] = meanA
 			typeOutput['meanB'] = meanB
 			typeOutput['meanC'] = meanC
 			typeOutput['A'] = A
			typeOutput['B'] = B
			if meanA < threshold and meanC < thresh: # "very" brackish; that is, the abundances are basically 0 on either side # TBD: Values are arbitrary
				typeOutput['type'] = Inter
				typeOutput['typeSimple'] = Inter
			elif meanA < thresh and meanC > thresh: # leaning towards being marine
				typeOutput['type'] = 'hi'+Inter
				typeOutput['typeSimple'] = Inter
			elif meanA > thresh and meanC < thresh: # leaning towards being fresh
				typeOutput['type'] = 'lo'+Inter
				typeOutput['typeSimple'] = Inter
			else: # both meanA and meanC are larger than 0
				typeOutput['type'] = 'half'+Inter
				typeOutput['typeSimple'] = Inter
		elif meanB < meanA and meanB < meanC: # inv-brackish water-- shouldn't exist according to hypothesis, but I put it in as a fail-safe
			isBrack = True
			typeOutput['type'] = 'inv'+Inter
			typeOutput['typeSimple'] = 'noclass'
			typeOutput['boundaries'] = [A,B]
			typeOutput['sigAB'] = pAB
			typeOutput['sigBC'] = pBC
			typeOutput['sigAC'] = pBC 
 			typeOutput['meanA'] = meanA
 			typeOutput['meanB'] = meanB
 			typeOutput['meanC'] = meanC
 			typeOutput['A'] = A
			typeOutput['B'] = B
		elif Stdb >= 2*Stda and Stdb >= 2*Stdc and meanC < thresh and meanA < thresh and meanB > meanC and meanB > meanA:
	 		isBrack = True
	 		typeOutput['type'] = 'noclass'
	 		typeOutput['bloom'] = Inter
	 		typeOutput['typeSimple'] = Inter
	 		typeOutput['boundaries'] = [A,B]
		else:
 			pass	# All other combos mean the brackish do NOT exist, so we can just compare A and C
	if isBrack: # If the group is already classified as either brackish or inv-brack, then we skip the next loop. If it is not classified, we continue.
		pass
# 	elif not sigAB and not sigBC and not sigAC:
# 		typeOutput['type'] = 'noclass' # A 'catch-all' for things that have NO significant differences between any of the groups
# 		Cv = numpy.std(listAbund[0])/average(listAbund[0])
# 		typeOutput['bloom'] = Cv
# 		typeOutput['boundaries'] = [A,B]
# 		typeOutput['sigAB'] = pAB
# 		typeOutput['sigBC'] = pBC
# 		typeOutput['sigAC'] = pBC 
#  		typeOutput['meanA'] = meanA
#  		typeOutput['meanB'] = meanB
#  		typeOutput['meanC'] = meanC
# 		typeOutput['A'] = A
# 		typeOutput['B'] = B


# 		print numpy.std(listAbund[0])
# 		print Cv
# 		maxMean = max([meanA,meanB,meanC])
# 		pos = [i for i, j in enumerate([meanA,meanB,meanC]) if j == maxMean]
# 		if len(pos) != 1:
# 			pass
# 		else:
# 			print pos
# 			print [groupA,groupB,groupC][int(pos)]
# 			StD = numpy.std([groupA,groupB,groupC][pos[0]]) # DOES NOT WORK WELL.
# 			uBar = [meanA,meanB,meanC][pos[0]]
# 			CvA = StD/uBar # StD/mean
# 			if CvA >= 3.0:
# 				typeOutput['bloom'] = 'yes'
# 	else:
# 		typeFM,boundaries,pAC,meanA,meanC,bEquation,inputA,inputB,bloom,noclassTF = notBrackishIteration(listAbund)
# 		if not noclassTF:
# 			typeOutput['type'] = typeFM
# 			typeOutput['boundaries'] = boundaries
# 			typeOutput['sigAC'] = pAC
# 			typeOutput['sigAB'] = None
# 			typeOutput['sigBC'] = None
# 			typeOutput['meanA'] = meanA
# 			typeOutput['meanC'] = meanC
# 			typeOutput['meanB'] = bEquation
# 			typeOutput['A'] = inputA
# 			typeOutput['B'] = inputB
# 			typeOutput['bloom'] = 'NA'
# 		else:
# 			typeOutput['type'] = 'noclass' # A 'catch-all' for things that have NO significant differences between any of the groups
# #  			Cv = numpy.std(listAbund[0])/average(listAbund[0])
# #  			typeOutput['bloom'] = Cv
#  			typeOutput['boundaries'] = [A,B]
#  			typeOutput['sigAB'] = pAB
# 	 		typeOutput['sigBC'] = pBC
# 	 		typeOutput['sigAC'] = pBC 
# 	  		typeOutput['meanA'] = meanA
# 	  		typeOutput['meanB'] = meanB
# 	  		typeOutput['meanC'] = meanC
# 	 		typeOutput['A'] = A
# 	 		typeOutput['B'] = B
# 	 		if Stdc >= 2*Stda and meanC > meanA and meanC > meanB and meanA < thresh:
# 	 			typeOutput['bloom'] = High
# 	 		elif Stda >= 2*Stdc and meanA > meanC and meanA > meanB and meanC < thresh:
# 	 			typeOutput['bloom'] =  Low
# 	 		else:
# 	 			typeOutput['bloom'] = 'No'
	 		
	 		
	elif (meanA > meanC and sigAC): #or (meanB > meanC and sigBC and (meanB-meanC) > (meanB-meanA)): # More in fresh water and it's significant
		# Above, you can have EITHER A>C or B>C but they must be significant, and if it's B, the distance between B and C must be greater than the distance between A and B (to prevent brackish-looking ones)
		# Note that meanB-meanC should ALWAYS be greater than meanB-meanA because if meanB-meanA is negative, it means it's truly fresh!
		if (sigAB and sigBC) or (not sigAB and not sigBC): # If there are actually significant difference between all three groups, then...
			typeOutput['boundaries'] = [scaleByDiff(A,B,meanA,meanB,meanC)] # We scale by how significant each difference is. (eg. If a-b is very significant but b-c is not very significant, then the 'true' boundary is approximated to be closer to a-b than to b-c. 
			# See function above for details
		elif sigAB and not sigBC: # If there are ONLY significance differences between a-b, then we choose this as the boundary
			typeOutput['boundaries'] = [scaleByDiff(A,B,meanA,meanB,meanC)]
		elif sigBC and not sigAC: # If there are ONLY significance difference between b-c, then we choose this as the boundary.
			typeOutput['boundaries'] = [scaleByDiff(A,B,meanA,meanB,meanC)]
		if meanC < thresh:
			typeOutput['type'] = Low
			typeOutput['typeSimple'] = Low
		else:
			typeOutput['type'] = 'half'+Low
			typeOutput['typeSimple'] = Low
	elif (meanC > meanA and sigAC): #or (meanB > meanA and sigAB and (meanB-meanA) > (meanB-meanC)): # This is same process as above, except for marine samples
		if (sigAB and sigBC) or (not sigAB and not sigBC):
			typeOutput['boundaries'] = [scaleByDiff(A,B,meanA,meanB,meanC)] # See above
		elif sigAB and not sigBC:
			typeOutput['boundaries'] = [scaleByDiff(A,B,meanA,meanB,meanC)]
		elif sigBC and not sigAB:
			typeOutput['boundaries'] = [scaleByDiff(A,B,meanA,meanB,meanC)]
		if meanA < thresh:
			typeOutput['type'] = High
			typeOutput['typeSimple'] = High
		else:
			typeOutput['type'] = 'half'+High 
			typeOutput['typeSimple'] = High
	else:
		typeOutput['type'] = 'noclass' # A 'catch-all' for things that have NO significant differences between any of the groups
#  			Cv = numpy.std(listAbund[0])/average(listAbund[0])
#  			typeOutput['bloom'] = Cv
		typeOutput['typeSimple'] = 'noclass'
		typeOutput['boundaries'] = ['',''] 
		typeOutput['sigAB'] = pAB
		typeOutput['sigBC'] = pBC
		typeOutput['sigAC'] = pBC 
		typeOutput['meanA'] = meanA
		typeOutput['meanB'] = meanB
		typeOutput['meanC'] = meanC
		typeOutput['A'] = A
		typeOutput['B'] = B
		if Stdc >= 2*Stda and meanC > meanA and meanC > meanB and meanA < thresh:
			typeOutput['bloom'] = High
			typeOutput['boundaries'] = [scaleByDiff(A,B,meanA,meanB,meanC)]
			typeOutput['typeSimple'] = High
		elif Stda >= 2*Stdc and meanA > meanC and meanA > meanB and meanC < thresh:
			typeOutput['bloom'] =  Low
			typeOutput['boundaries'] = [scaleByDiff(A,B,meanA,meanB,meanC)]
			typeOutput['typeSimple'] = Low
		else:
			typeOutput['bloom'] = 'No'

	return typeOutput # Output is a dictionary, in each taxa name there is 'type' and 'boundaries'; also, meanA, meanB, meanC, A, B-- this is for downstream stuff

#==========================================
# Take 'boundaries' that were significantly different and plot them to see how boundaries are distributed across salinities

def summaryBoundaryTypes(taxaInfo):
	transitionsList = []
	transitions = []
	for taxa in taxaInfo:
		transitionsList = taxaInfo[taxa]['boundaries']
		if isinstance(transitionsList, list):
			for i in transitionsList:
				if i == '':
					pass
				else:
					transitions.append(i)
		else:
			if transitionsList == '':
				pass
			else:
				transitions.append(transitionsList)	
# 	transitionsFiltered = [x for x in transitions if x is not None]
	return transitions # Output is summary of all boundaries in all taxa, excluding 'Nones'

#==========================================
# BINNING

def sequenceGenerator(min,max,binSize): # Make sequence from max value, min value, and binsize. List includes min and max bin values. Integers.
	if max < min:
		print "Error, max is less than min"
		return
	breadth = max-min
	nbins = math.ceil(float(breadth)/float(binSize))
	nbins = int(nbins) # Make into integer from float
	seq = []
	for i in range(nbins+1):
		current = float(min) + float(binSize)*i
		seq.append(current)
	return seq
# 	return nbins,seq

# def binIt(upper,lower,list,binSize): # List is a list of all values to be sorted into bins
# 	numBins,binRanges = sequenceGenerator(lower,upper,binSize)
# 	for i in range(len(binRanges)-1):
# 		binCount=[0]*numBins
# 		for x in list:
# 			if x==upper:		# treats last bin different
# 				binCount[numBins-1]+=1 # Put in last bin
# 			else:
# 				location=int((x-lower)*numBins/(upper-lower))		#puts each value into a bin
# 				binCount[location]+=1
# 	return binCount,binRanges

#==========================================
# Determine composition (fresh, marine, brackish etc) within each salinity bin

class compositionAtSalinity():
	def __init__(self):
		self.lowGrad = 0
		self.halfLowGrad = 0
		self.interGrad = 0
		self.hiInterGrad = 0
		self.loInterGrad = 0
		self.halfInterGrad = 0
		self.invInterGrad = 0
		self.highGrad = 0
		self.halfHighGrad = 0
		self.noclass = 0
		self.bloomLowGrad = 0
		self.bloomHighGrad = 0
		self.bloomInterGrad = 0
		
	def addValue(self,type,value,bloom): # Add 'raw' relative abundances from original taxasummaries
		if type == Low:
			self.lowGrad += value
		if type == 'half'+Low:
			self.halfLowGrad += value
		if type == Inter:
			self.interGrad += value
		if type == 'hi'+Inter:
			self.hiInterGrad += value
		if type == 'lo'+Inter:
			self.loInterGrad += value
		if type == 'half'+Inter:
			self.halfInterGrad += value
		if type == High:
			self.highGrad += value
		if type == 'half'+High:
			self.halfHighGrad += value
		if type == 'noclass':
			if bloom == 'No':
				self.noclass += value
			if bloom == Low:
				self.bloomLowGrad += value
			if bloom == High:
				self.bloomHighGrad += value
			if bloom == Inter:
				self.bloomInterGrad += value
		if type == 'inv-brackish':
			self.halfInterGrad += value

	def reportComp(self):
		total = self.lowGrad + self.halfLowGrad + self.interGrad + self.hiInterGrad + self.loInterGrad + self.halfInterGrad + self.highGrad + self.halfHighGrad + self.noclass + self.invInterGrad + self.bloomLowGrad + self.bloomHighGrad + self.bloomInterGrad
		if total == 0: # Aka, there are no samples that exist in the data
			summary = 'Pass' 
		else: # By standardizing with "total", it should be a new 'average' relative abundance for each salinity bin.
			lowGrad = float(self.lowGrad)/total
			halfLowGrad = float(self.halfLowGrad)/total
			interGrad = float(self.interGrad)/total
			hiInterGrad = float(self.hiInterGrad)/total
			loInterGrad = float(self.loInterGrad)/total
			halfInterGrad = float(self.halfInterGrad)/total
			highGrad = float(self.highGrad)/total
			halfHighGrad = float(self.halfHighGrad)/total
			noclass = float(self.noclass)/total
			invInterGrad = float(self.invInterGrad)/total
			bloomLowGrad = float(self.bloomLowGrad)/total
			bloomHighGrad = float(self.bloomHighGrad)/total
			bloomInterGrad = float(self.bloomInterGrad)/total
			summary = [lowGrad,bloomLowGrad,halfLowGrad,loInterGrad,interGrad,bloomInterGrad,halfInterGrad,hiInterGrad,halfHighGrad,highGrad,bloomHighGrad,noclass,invInterGrad]
		return summary # Output is list in order listed above.


def countListEvenLevels(listToCount,AllOptions):
	returnList = {}
	for i in AllOptions:
		returnList[i] = 0
		for j in listToCount:
			if j == i:
				returnList[i] += 1
			else:
				pass
	return returnList

#=========================================================


#=========================================================

parser = argparse.ArgumentParser(
	description="Bins and classifies OTUs according to salinity specialization")
parser.add_argument(
	'-t',
	'--taxaTable',
	help = "Taxa table file from QIIME-- biom format",
	required = True,)
parser.add_argument(
	'-m',
	'--metadata',
	help = 'File where first column are site name and second column is gradient values. No header required',
	required = True)
parser.add_argument(
	'-M',
	'--metadata_name',
	help = 'Header name of column in metadata file that specifies the variable you want',
	required = True)
parser.add_argument(
	'-o',
	'--output_dir',
	help = 'Output directory [default: GradientProcessing]',
	required = False,
	default = 'GradientProcessing')
parser.add_argument(
	'-g',
	'--gradient',
	help = 'Names for three points of gradient-- must have low, intermediate, and high name (Eg. Fresh, Brackish, Marine). Comma separated. [default: Low, Inter, High]',
	required = False,
	type = str,
	default = 'Low,Inter,High')
parser.add_argument(
	'--min_threshold_proportion',
	help = 'Proportional threshold for OTU abundances to be considered "Major" [default 0.10]',
	required = False,
	type = float,
	default = 0.10)
parser.add_argument(
	'--min_threshold_constant',
	help = 'Constant threshold for OTU abundances to be considered "Major" [default None]',
	required = False,
	default = None)
parser.add_argument(
	'--minA',
	help = "Minimum boundary A can be [default is minVal]",
	required = False,
	default = 'Check')
parser.add_argument(
	'--maxB',
	help = 'Maximum boundary B can be [default is maxVal]',
	required = False,
	default = 'Check')
parser.add_argument(
	'--ABdiff',
	help = 'Difference between A ad B [default 2 units]',
	required = False,
	type = float,
	default = 2)
parser.add_argument(
	'--division_Size',
	help = 'Division size for bins for graphing bar graph [default 2]',
	required = False,
	default = 2)
parser.add_argument(
	'--unit_Size',
	help = 'Unit size to iterate through [default 1]',
	required = False,
	default = 1)
parser.add_argument(
	'-T',
	'--tree',
	help = 'Pathway to tree, required',
	required = True)
parser.add_argument(
	'-a',
	'--allBins',
	action =  'store_false',
	help = 'Include this flag if you DO NOT want to make generate all bins. [Default: True]',
	required = False,
	default = True)
parser.add_argument(
	'-c',
	'--condensedBins',
	action =  'store_false',
	help = 'Include this flag if you DO NOT want to make condensed bins. [Default: True]',
	required = False,
	default = True)
parser.add_argument(
	'-R',
	'--R_script',
	help = 'If given, will graph outputs through R script. Need pathway to Rscript file. If not provided, will not graph.',
	required = False,
	default = False)
	
args = parser.parse_args()

taxaTable = args.taxaTable
metadata = args.metadata
metadata_name = args.metadata_name
gradient = args.gradient
minA = args.minA
maxB = args.maxB
ABdiff = args.ABdiff
min_threshold_proportion = args.min_threshold_proportion
min_threshold_constant = args.min_threshold_constant
divisionSize = args.division_Size
unitSize = args.unit_Size
output_dir = args.output_dir
allBins = args.allBins
condensedBins = args.condensedBins
tree = args.tree
R_script = args.R_script


#==========================================
# Prestep: Check to see if macqiime is loaded

# How to do this??

# Step one: load files and set variables

gradient = gradient.split(',')
Low = gradient[0]
Inter = gradient[1]
High = gradient[2]

ABdiffunit = ABdiff*unitSize

# Make taxa summary using function
taxasummariesRaw,taxaIDs = makeTaxaSummaries(taxaTable, metadata)
taxasummaries = deleteLowAbund(taxasummariesRaw)
# if type(taxasummaries.get('Unassigned;Other;Other;Other;Other;Other')) == list:
# 	del taxasummaries['Unassigned;Other;Other;Other;Other;Other']
# if type(taxasummaries.get('Unassigned')) == list:
# 	del taxasummaries['Unassigned'] # FIX THIS
	

maxVal = max(taxasummaries[taxasummaries.keys()[1]][1])
minVal = min(taxasummaries[taxasummaries.keys()[1]][1])

if str(minA) == 'Check':
	minA = int(minVal)
else:
	minA = int(minA)
if str(maxB) == 'Check':
	maxB = int(maxVal)
else:
	maxB = int(maxB)

# If 'notes' flag is true, then you will be given a chance to make an input

LOG = open("LOG.txt", 'wr')
LOG.write(output_dir + '\n')
# if notes == True:
# 	LOG.write(raw_input("Type short description or comment, then press enter:\n") + '\n')
LOG.write("Settings: \n"+ "minA =" + str(minA) + "\nmaxB =" + str(maxB) + "\nABdiff =" + str(ABdiff) + "\nminVal =" + str(minVal) + "\nmaxVal =" + str(maxVal) + "\nDivSize =" + str(divisionSize) + "\nunitSize = " + str(unitSize))
LOG.write("\nthreshold =" + str(min_threshold_proportion))
LOG.write("\nFull PWD: \ntaxasummaries = " + taxaTable + "\nmetadata = " + metadata + "\ntree = " + tree)
LOG.write("\nLow: "+ Low)
LOG.write("\nInter: "+ Inter)
LOG.write("\nHigh: "+ High)
LOG.write("\nGradient Header: "+ metadata_name)
LOG.write("\nallBins,condensedBins: "+ str(allBins) + ","+ str(condensedBins))
LOG.close()


# Check if they want an even threshold or a proportional threshold
threshold = [True,min_threshold_proportion]
if min_threshold_constant == None:
	threshold = [True,min_threshold_proportion]
else:
	threshold = [False,min_threshold_constant]

#==========================================
# Step two: Make data


# Dictionary that will have taxa as keys and a short list of paired [A,B]
bestfitAB = {}
taxaInfo = {}
# modelDiffList = {}

print "Iterating through all combinations..."
# Go through each taxa and iterate through all combinations of boundaries A and B
for taxa in taxasummaries:
	currentbest = 0 # going to compare 1/n of error squared because I know it can't be lower than 0. Conversely, not sure what maximum of n will be.
	modelDiffList = []
	listAbund = taxasummaries[taxa]
	for i in sequenceGenerator(minA,(maxB-ABdiffunit),unitSize):
		A = i
		for j in sequenceGenerator((A+ABdiffunit),maxB,unitSize):
			B = j
			binValues = sortBins(A,B,listAbund)
			if len(binValues['lista']) <= 3 or len(binValues['listb']) <= 3 or len(binValues['listc']) <= 3: # if the bins have nothing in them, then don't use that bin combination
				pass
			else:
				combined = makeAModel(A,B,binValues)
				MES = meanErrorSquared(combined)
				if MES == 0: # If the meanErrorSquared is EXACTLY 0, there is probably something wrong, OR it is in such low abundance it's not worth lookinga t
					print "MES IS ZERO"
					pass
				else:
					invMES = 1/MES
					modelDiffAsq = (average(binValues['lista']) - average(binValues['listb']))**2
					modelDiffBsq = (average(binValues['listb']) - average(binValues['lista']))**2
					modelDiff = (modelDiffAsq + modelDiffBsq)
					if invMES > currentbest:
						bestfitAB[taxa] = [[A,B]]
						currentbest = invMES
						modelDiffList = [modelDiff]
					elif invMES == currentbest:
						bestfitAB[taxa].append([A,B])
						modelDiffList.append([modelDiff])
	firstBoundary = []
	secondBoundary = []
	for i in bestfitAB[taxa]:
		firstBoundary.append(i[0])
		secondBoundary.append(i[1])
# 	overallDiff = []
# 	for i in modelDiffList:
# 		overallDiff.append(i)
	maxDifferenceFoundPosition = [modelDiffList.index(i) for i in modelDiffList if i == max(modelDiffList)]
	differencesA = numpy.diff(firstBoundary)
	differencesB = numpy.diff(secondBoundary)
	differencesAB = list(differencesA) + list(differencesB)
	if len(maxDifferenceFoundPosition) > 1 and (True in [True for i in differencesAB if i > 1]): # If consecutive numbers AND same diff
		print "WARNING: SAME MODELDIFF FOR IDENTICAL MODELS-Br"
		print taxa
		print firstBoundary
		print secondBoundary
	firstBoundaries = []
	secondBoundaries = []
	for i in maxDifferenceFoundPosition:
		firstBoundaries.append(firstBoundary[i])
		secondBoundaries.append(secondBoundary[i])
	aveFirst = average(firstBoundaries)
	aveSecond = average(secondBoundaries)
	# Here, I make the assumption that all 'equal' bestfit boundaries are essentially equivalent. Aka, they are actually salinities that are lacking data and thus will repeat multiple times.
	# There is a possibility that the error is EXACTLY equal with two totally different boundaries, but I'm assuming they aren't.
# 	if len(firstBoundary) > 1 or len(secondBoundary) > 1:
# 		differencesA = numpy.diff(firstBoundary)
# 		differencesB = numpy.diff(secondBoundary)
# 		differencesAB = list(differencesA) + list(differencesB)
# 	else:
# 		differencesAB = [False]
# 	if (True in [True for i in differencesAB if i > 1]): # if consecutive numbers, then just average them. if not, then print warning but average them anyway.
# 		print "WARNING: TWO OR MORE IDENTICAL MODELS-Br"
# 	aveFirst = average(firstBoundary)
# 	aveSecond = average(secondBoundary)
	typeInfo = typeTaxa(aveFirst,aveSecond,listAbund) # type is going to be list of type and boundaries
	taxaInfo[taxa] = typeInfo


transitions = summaryBoundaryTypes(taxaInfo) # list of all boundaries
binRanges = sequenceGenerator(minVal,maxVal,divisionSize)
# binCount, binRanges = binIt(maxVal,minVal,transitions,divisionSize) # input for histogram; don't actually need binCount but need binRanges

#==========================================
# Step three: make list of percent composition of salinities

# Make a file out of the taxasummaries file where you sum up the relative abundances at each salinity for each organism
# maxValinity = max(binRanges)

if allBins == True:
	compositionAtGradAll = [0]*(len(binRanges)-1) # Empty matrix
	for i in range(len(binRanges)-1): # For each bin
		compositionAtGradAll[i] = compositionAtSalinity()
		for x in taxasummaries: # For each taxa
			type = taxaInfo[x]['type']
			bloom = taxaInfo[x]['bloom']
			for y in range(len(taxasummaries[x][1])): # Add the relative abundance to the class, "compositionAtSalinity"
				if taxasummaries[x][1][y] >= binRanges[i] and taxasummaries[x][1][y] < binRanges[i+1]:
					compositionAtGradAll[i].addValue(type,taxasummaries[x][0][y],bloom)

	composition = [0]*len(compositionAtGradAll)
	for i in range(len(compositionAtGradAll)):
		composition[i] = compositionAtGradAll[i].reportComp() # Report the raw numbers for composition at each salinity
		# Order is: fresh,halffresh,frebrackish,brackish,halfbrackish,marbrackish,halfmarine,marine,noclass,invbrackish
if condensedBins == True:
	compositionAtGradCondensed = [0]*(len(binRanges)-1) # Empty matrix
	for i in range(len(binRanges)-1): # For each bin
		compositionAtGradCondensed[i] = compositionAtSalinity()
		for x in taxasummaries: # For each taxa
			type = taxaInfo[x]['typeSimple']
			bloom = taxaInfo[x]['bloom']
			for y in range(len(taxasummaries[x][1])): # Add the relative abundance to the class, "compositionAtSalinity"
				if taxasummaries[x][1][y] >= binRanges[i] and taxasummaries[x][1][y] < binRanges[i+1]:
					compositionAtGradCondensed[i].addValue(type,taxasummaries[x][0][y],bloom)

	compositionCondensed = [0]*len(compositionAtGradCondensed)
	for i in range(len(compositionAtGradCondensed)):
		compositionCondensed[i] = compositionAtGradCondensed[i].reportComp() # Report the raw numbers for composition at each salinity
		# Order is: fresh,halffresh,frebrackish,brackish,halfbrackish,marbrackish,halfmarine,marine,noclass,invbrackish

			
#==========================================
# Making OTU presence/absence table

allTypes = []
for i in taxaInfo:
	if taxaInfo[i]['type'] == 'noclass':
		if taxaInfo[i]['bloom'] == 'No':
			allTypes.append('noclass')
		else:
			allTypes.append(taxaInfo[i]['bloom'])
	else:
		allTypes.append(taxaInfo[i]['type'])
allTypesUnique = list(set(allTypes))

OTUtableByType = {}
uniqueOTUList = []
for type in allTypesUnique: 
	OTUtableByType[type] = []
	for OTU in taxaInfo.keys():
		uniqueOTUList.append(OTU)
		if taxaInfo[OTU]['type'] == type:
				OTUtableByType[type].append(OTU)
		else: 
			pass
for type in allTypesUnique:
	tempdict = countListEvenLevels(OTUtableByType[type], uniqueOTUList)
	OTUtableByType[type] = tempdict
	
#==========================================




#==========================================
# Step four: Output the numerical data so we can look at it manually.

# combined = result from makeAModel
# transitions = list of all boundaries; can do histogram with output
# bestfitAB = dictionary, where indices are 'taxa' and then it lists their boundaries
# listAbund = [abundance,salinity]
print "Printing and saving..."

os.system('mkdir ' + output_dir)
os.chdir(output_dir)

if allBins == True:
	typeAbundance = open('types_across_gradient_all.txt', 'wr')
	toWrite = 'Gradient' + '\t'
	for i in range(len(binRanges)-1):
		toWrite += str(binRanges[i+1]) + '\t'
	toWrite = toWrite.strip()
	toWrite += '\r'
	typeAbundance.write(toWrite)
	categories = [Low,'bloom'+Low,'half'+Low,'lo'+Inter,Inter,'bloom'+Inter,'half'+Inter,'hi'+Inter,'half'+High,High,'bloom'+High,'noclass','inv'+Inter]
	for i in range(len(categories)):
		toWrite = categories[i] + '\t'
		for j in composition:
			if j == 'Pass':
				toWrite += '\t'
			else:
				toWrite += str(j[i]) + '\t'
		toWrite = toWrite.strip()
		toWrite += '\r'
		typeAbundance.write(toWrite)
	typeAbundance.close()

if condensedBins == True: # TESTING
	typeAbundance = open('types_across_gradient_condensed.txt', 'wr')
	toWrite = 'Gradient' + '\t'
	for i in range(len(binRanges)-1):
		toWrite += str(binRanges[i+1]) + '\t'
	toWrite = toWrite.strip()
	toWrite += '\r'
	typeAbundance.write(toWrite)
	categories = [Low,'bloom'+Low,'half'+Low,'lo'+Inter,Inter,'bloom'+Inter,'half'+Inter,'hi'+Inter,'half'+High,High,'bloom'+High,'noclass','inv'+Inter]
	positionCategories = [0,4,9,11]
	for i in positionCategories:
		toWrite = categories[i] + '\t'
		for j in compositionCondensed:
			if j == 'Pass':
				toWrite += '\t'
			else:
				toWrite += str(j[i]) + '\t'
		toWrite = toWrite.strip()
		toWrite += '\r'
		typeAbundance.write(toWrite)
	typeAbundance.close()
	

taxaFile = open('taxa_abundances_across_gradient.txt', 'wr')
firstLine = 'Gradient\t'
for x in taxasummaries.values()[0][1]:
	firstLine += str(x) + '\t'
firstLine = firstLine.strip()
firstLine += '\r'
taxaFile.write(firstLine)
for taxa in taxasummaries:
	toWrite = taxa + '\t'
	for i in taxasummaries[taxa][0]:
		toWrite += str(i) + '\t'
	toWrite = toWrite.strip()
	toWrite += '\r'
	taxaFile.write(toWrite)
taxaFile.close()

listBoundaries = open('boundaries.txt','wr')
for x in transitions:
	listBoundaries.write(str(x) + '\r')
listBoundaries.close()

bestfitBoundaries = open('modelBoundaries_type.txt', 'wr')
header = 'taxa\ttype\ttypeSimple\tmeanA\tmeanB\tmeanC\tA\tB\tsigAB\tsigBC\tsigAC\tbloom\tboundaries\tboundariestwo\r'
bestfitBoundaries.write(header)
for taxa in taxaInfo:
	toWrite = taxa + '\t'
	toWrite += taxaInfo[taxa]['type'] + '\t'
	toWrite += taxaInfo[taxa]['typeSimple'] + '\t'
	toWrite += str(taxaInfo[taxa]['meanA']) + '\t'
	toWrite += str(taxaInfo[taxa]['meanB']) + '\t'
	toWrite += str(taxaInfo[taxa]['meanC']) + '\t'
	toWrite += str(taxaInfo[taxa]['A']) + '\t'
	toWrite += str(taxaInfo[taxa]['B']) + '\t'
	toWrite += str(taxaInfo[taxa]['sigAB']) + '\t'
	toWrite += str(taxaInfo[taxa]['sigBC']) + '\t'
	toWrite += str(taxaInfo[taxa]['sigAC']) + '\t'
	toWrite += str(taxaInfo[taxa]['bloom']) + '\t'
	for i in taxaInfo[taxa]['boundaries']:
		toWrite += str(i) + '\t'
	toWrite = toWrite.strip()
	toWrite += '\r'
	bestfitBoundaries.write(toWrite)
bestfitBoundaries.close()

taxaIDsprint = open('taxaIDLegend.txt','wr')
toWrite = []
for taxa in taxaIDs:
	toWrite = str(taxa) + '\t' + str(taxaIDs[taxa]) + '\r'
	taxaIDsprint.write(toWrite)
taxaIDsprint.close()

gradientTemp = open("gradient.txt", 'w')
gradientTemp.write(Low+','+Inter+','+High+','+metadata_name+'\r')
gradientTemp.close()

# Print OTU presence/absence table 
OTUtabletoPrint = open('OTUTablebyType.txt','w')
firstDone = False
for type in OTUtableByType:
	lineToPrint = ''
	if not firstDone:
		lineToPrint += "#OTUTable \t"
		for i in OTUtableByType[type].keys():
			lineToPrint += str(i) + "\t"
		lineToPrint = lineToPrint.strip()
		lineToPrint += "\r"
		firstDone = True
	lineToPrint += type +"\t"
	for i in OTUtableByType[type]:
		lineToPrint += str(OTUtableByType[type][i]) + "\t"
	lineToPrint = lineToPrint.strip()
	lineToPrint += "\r"
	OTUtabletoPrint.write(lineToPrint)
OTUtabletoPrint.close()

uniqueListOfOTUs = open('OTUs.txt','w')
toPrint = ''
OTUs = set(uniqueOTUList)
for taxa in taxaInfo:
	toPrint += str(taxa) + '\r'
toPrint = toPrint.strip()
uniqueListOfOTUs.write(toPrint)
uniqueListOfOTUs.close()

os.system('mv ../LOG.txt ./')
os.system('mv ../OTUTableText.txt ./')

#==========================================
# Step five: Draw the histogram and the other graphs
# Using R
print "Drawing graphs..."

if isinstance(R_script,str):
	subprocess.call(['Rscript',str(R_script)])

print 'DONE'

