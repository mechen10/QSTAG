#!/bin/bash

# Script is for figuring out whether my binning script's classifications have taxonomic patterns
#============================================

import argparse
import os
import numpy

#============================================

parser = argparse.ArgumentParser(
	description = "This program determines which taxonomies are binned as what type")
parser.add_argument(
	'-m',
	'--modelboundaries',
	help = "Model Boundaries output from salinity binning script",
	required = True)
parser.add_argument(
	'-t',
	'--taxaID',
	help = 'Taxa ID output from salinity binning script',
	required = True)
parser.add_argument(
	'-o',
	'--output',
	help = 'Output folder',
	required = False,
	default = 'ClassificationTaxonomy')
parser.add_argument(
	'-d',
	'--delimiter',
	help = 'The delimiter for taxonomy legen. Default: [" __; "]',
	required = False,
	default = ' __; ')
args = parser.parse_args()

modelboundaries = args.modelboundaries
taxaID = args.taxaID
output = args.output
delimiter = args.delimiter

#============================================

os.system("mkdir " + output)
	
# Load file
filemb = open(modelboundaries,'r')
for i in filemb:
	newFile = i.strip()
	newFile = newFile.split("\r")
filemb.close()

# Make file
modelBound = []
first = True
for line in newFile:
	newLine = line
	sepLine = newLine.split("\t")
	if first:
		headerNames = sepLine
		first = False
	else:
		modelBound.append(sepLine)

# Load file	
mbDict = {}
for r in modelBound:
	mbDict[str(r[0])] = {}
	for n in range(len(r)):
		if n == 0:
			continue
		else:
			mbDict[str(r[0])][headerNames[n]] = r[n]
			
filet = open(taxaID,'r')
for i in filet:
	newFile = i.strip()
	newFile = newFile.split("\r")
filet.close()

# Make taxaID file
taxaIDs = {}
for line in newFile:
	newLine = line
	sepLine = newLine.split("\t")
	taxaIDs[sepLine[0]] = sepLine[1]
	
for taxa in taxaIDs:
	newTaxa = taxaIDs[taxa].split(delimiter)
	if newTaxa == 'Unassigned':
		taxaIDs[taxa] = ['Unassigned']
	else:
		taxaIDs[taxa] = newTaxa

#============================================
# FUNCTIONS

def countList(listToCount): # Given a list, it will count how many of each element there are
	unique = list(set(listToCount))
	returnList = {}
	for i in unique:
		count = 0
		for j in listToCount:
			if j == i:
				count += 1
			else:
				pass
		returnList[i] = count
	return returnList

def countListEvenLevels(listToCount,AllOptions):
	returnList = {}
	# For sanity check, sum of this should equal listToCount
	totalCount = 0
	for i in AllOptions:
		returnList[i] = listToCount.count(i)
		totalCount += listToCount.count(i)
		# if i in listToCount:
# 			returnList[i] += 1
# 			totalCount += 1
	if totalCount == len(listToCount):
		return returnList
	else:
		print "ERROR-- check countListEvenLevels function. Not all items counted."
			
		# for j in listToCount:
# 			if j == i:
# 				returnList[i] += 1
# 			else:
# 				pass
# 	return returnList
			
def countTaxa(listToCount): # tally's up the fresh, brack, etc within each taxa. Input is a 2-element list of [taxa,type]
	 returnDict = {}
	 allTaxa = []
	 allTypes = []
	 for i in listToCount:
	 	allTaxa.append(i[0])
	 	allTypes.append(i[1])
	 allTaxaSetList = list(set(allTaxa))
	 allTypesSetList = list(set(allTypes))
	 for i in allTaxaSetList:
	 	tempTypeList = []
	 	for j in listToCount:
	 		if i == j[0]:
	 			tempTypeList.append(j[1])
	 		else: 
	 			pass
	 	tempTypeListCounted = countListEvenLevels(tempTypeList, allTypesSetList)
	 	returnDict[i] = tempTypeListCounted
	 return returnDict
	
	 
def percentageCalc(threeLayerDictionary): # Three layers dictionary: [Level][taxa][type] OR [Type][Level][taxa]
	percentV = threeLayerDictionary
	for i in percentV:
		for j in percentV[i]:
			totalCount = 0
			for k in percentV[i][j]:
				totalCount += percentV[i][j][k]
			if totalCount == 0:
				pass
			else:
				for k in percentV[i][j]:
					percentV[i][j][k] = float(percentV[i][j][k])/float(totalCount)
	return percentV
	
def printDictionary(threeLayerDictionary,foldername):
	os.system('mkdir ' + output + '/' + foldername)
	for i in threeLayerDictionary: # Each Level = separate file
		fileopen = open(str(output+'/' + foldername + '/'+ i + '.txt'), 'w')
		toPrint = ''
		headerDone = False
		for j in threeLayerDictionary[i]: # Each taxa or type
			newLine = ''
			if not headerDone:
				for headerkeys in threeLayerDictionary[i][j].keys(): # Take first one and make headers
					toPrint += headerkeys + '\t'
				toPrint = toPrint.strip()
				toPrint = '\t'+toPrint+'\r'
				headerDone = True
			newLine = j + '\t'
			for k in threeLayerDictionary[i][j]:
				newLine += str(threeLayerDictionary[i][j][k]) + '\t'
			newLine = newLine.strip()
			newLine += '\r'
			toPrint += newLine
		fileopen.write(toPrint)
		fileopen.close()

def weightedAverage(twoLayerList):
	weightedList = []
	for i in twoLayerList:
		value = i[0]
		weight = int(i[1]**(1/2))
		for i in range(weight):
			weightedList.append(value)
	weightedAve = sum(weightedList)/len(weightedList)
	return weightedAve
	
def dictToTable(dictOne,dictTwo):
	OTUID = []
	firstCol = []
	secondCol = []
	for i in dictOne:
		if i == 'unassigned':
			pass
		else:
			OTUID.append(i)
			firstCol.append(dictOne[i])
	for i in dictTwo:
		if i == 'unassigned':
			pass
		else:
			secondCol.append(dictTwo[i])
	table = numpy.array([firstCol,secondCol])
	return table	
	
def bray_curtis_distance(table, sample1_id, sample2_id):
    numerator = 0
    denominator = 0
    sample1_counts = table[sample1_id]
    sample2_counts = table[sample2_id]
    for sample1_count, sample2_count in zip(sample1_counts, sample2_counts):
        numerator += abs(sample1_count - sample2_count)
        denominator += sample1_count + sample2_count
    return float(numerator) / denominator
			
#============================================
# STATS

# # Change types of noclass to 'bloom'
# for taxa in mbDict:
# 	if mbDict[taxa]['type'] == 'noclass':
# 		if mbDict[taxa]['bloom'] == 'No':
# 			pass
# 		else:
# 			mbDict[taxa]['type'] = 'bloom'+mbDict[taxa]['bloom']
# 	else:
# 		pass
# 		


# Get all types
allTypes = []
for taxa in mbDict:
	allTypes.append(mbDict[taxa]['type'])
allTypesUnique = list(set(allTypes))

# Make list of all types
byType = {}

for i in allTypesUnique:
	byType[i] = []
for taxa in mbDict:
	tempType = mbDict[taxa]['type']
	byType[tempType].append(taxa)

# ============================================
# Get 3-layer dictionaries

taxonomyList = ['Domain','Phylum','Class','Order','Family','Genus','Species']
byTypeList = {}
for level in range(len(taxonomyList)):
	byTypeList[taxonomyList[level]] = {}
	typesListTemp = {}
	allTaxaList = []
	for type in allTypesUnique:
		typesListTemp[type] = []
# 		byTypeList[taxonomyList[level]][type] = []
		for OTU in mbDict.keys():
			if mbDict[OTU]['type'] == type:
				if len(taxaIDs[OTU]) > level:
					if level == 6:
						typesListTemp[type].append(str(taxaIDs[OTU][5]+'_'+taxaIDs[OTU][6]))
						allTaxaList.append(str(taxaIDs[OTU][5]+'_'+taxaIDs[OTU][6]))
					else:
						typesListTemp[type].append(taxaIDs[OTU][level])
						allTaxaList.append(taxaIDs[OTU][level])
				else:
					name = taxaIDs[OTU][len(taxaIDs[OTU])-1]
					nUnass = level - (len(taxaIDs[OTU])-1)
					for nU in range(0,nUnass):
						name += '_unassigned'
					typesListTemp[type].append(name)
					allTaxaList.append(name)
					# name = ''
# 					for j in range(len(taxaIDs[OTU])):
# 						name += taxaIDs[OTU][j] + '_'
# 					name += 'unassigned'
# 					typesListTemp[type].append(name)
# 					allTaxaList.append(name)
			else: 
				pass
	uniqueTaxaList = list(set(allTaxaList))
	for type in allTypesUnique:
		tempdict = countListEvenLevels(typesListTemp[type], uniqueTaxaList)
		byTypeList[taxonomyList[level]][type] = tempdict

######------------------------------------------
# Version where OTU IDs are retained; to be used downstream for R Unifrac	
uniqueOTUList = []
for OTU in mbDict.keys():
	uniqueOTUList.append(OTU)

OTUtableByType = {}
for type in allTypesUnique:
	OTUtableByType[type] = []
	for OTU in mbDict.keys():
# 		uniqueOTUList.append(OTU)
		if mbDict[OTU]['type'] == type:
				OTUtableByType[type].append(OTU)
		else: 
			pass
for type in allTypesUnique:
	tempdict = countListEvenLevels(OTUtableByType[type], uniqueOTUList)
	OTUtableByType[type] = tempdict
# Print OTU table
OTUtabletoPrint = open(output + '/OTUTablebyType.txt','w')
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
		
		
######------------------------------------------
# Make dictionary: first set of keys are 'Levels' [Domain, Phylum, etc]. Second layer is taxa (bacteria, archaea, etc). Third is composition of types
# Make list of taxonomies
byTaxaList = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
byTaxaListDict = {}

for i in range(len(byTaxaList)):
	taxaListTemp = []
	for taxa in mbDict:
	# Determine the type of each taxa
		currentType = mbDict[taxa]['type']
# 		if mbDict[taxa]['type'] == 'noclass':
# 			currentType = mbDict[taxa]['bloom']
# 			if currentType == 'No':
# 				currentType = 'noclass'
# 		else: 
# 			currentType = mbDict[taxa]['type']
	# Then take this and input it into list
		if len(taxaIDs[taxa]) > i: # Get all taxonomy IDs and make set list
			taxaListTemp.append([taxaIDs[taxa][i],currentType])
		else:
			name = taxaIDs[taxa][len(taxaIDs[taxa])-1]
			nUnass = i - (len(taxaIDs[taxa])-1)
			for nU in range(0,nUnass):
				name += '_unassigned'
			taxaListTemp.append([name, currentType])
	taxaListTemp = countTaxa(taxaListTemp)
	byTaxaListDict[byTaxaList[i]] = taxaListTemp

# ====================================
# Calculate 'uneveness' between two groups. 

# For each 'taxa', you want to give points for either evenness or unevennes
# BUT you also want to weight evenly between 'taxa'
# Each taxa unevenness has a scale of 0-1; 1 means totally absent in one
# 0 means totally even in abundance
# then average all 'uneveness' so that the number of OTUs doesn't matter
# Should weight by count?

# Ignore both 0 ones
# Take ratio. multiply by total raw observations--> this weights by overall 'sample size'.
# Then, average this across all taxa
	

unevenness = {}
brayCurtis = {}
for level in byTypeList:
	unevenness[level] = {}
	brayCurtis[level] = {}
	for firstGroup in range(len(byTypeList[level].keys())):
		firstGroupName = byTypeList[level].keys()[firstGroup]
		BCdictOne = byTypeList[level][firstGroupName]
		for secondGroup in range(firstGroup+1,len(byTypeList[level].keys())):
			secondGroupName = byTypeList[level].keys()[secondGroup]
			colName = firstGroupName + 'vs' + secondGroupName
			unevenness[level][colName] = []
			brayCurtis[level][colName] = []
			allDiff = []
			# BRAY CURTIS
			BCdictTwo = byTypeList[level][secondGroupName]
			tempTable = dictToTable(BCdictOne,BCdictTwo)
			BCdistance = bray_curtis_distance(tempTable,0,1)
			# UNEVENNESS
			for i in byTypeList[level][firstGroupName].keys():
				if i == 'unassigned':
					pass
				else:
					firstVal = byTypeList[level][firstGroupName][i]
					secondVal = byTypeList[level][secondGroupName][i]
					if (firstVal == 0 and secondVal == 0):
						pass
					else:
						allDiff.append([float(abs(firstVal-secondVal))/(firstVal+secondVal),(firstVal + secondVal)])
			unevenness[level][colName] = weightedAverage(allDiff)
			brayCurtis[level][colName] = BCdistance
			
# Print the unevenness data and brayCurtis data
unevennessToPrint = open(str(output+'/UnevennessData.txt'), 'w')
header = False
for i in taxonomyList:
	lineToPrint = i + '\t'
	if not header:
		headerToPrint = '#TYPE \t'
		for comparison in unevenness[i].keys():
			headerToPrint += comparison + '\t'
		headerToPrint = headerToPrint.strip()
		headerToPrint += '\r'
		unevennessToPrint.write(headerToPrint)
		header = True
	for j in unevenness[i]:
		lineToPrint += str(unevenness[i][j]) + '\t'
	lineToPrint = lineToPrint.strip()
	lineToPrint += '\r'
	unevennessToPrint.write(lineToPrint)
unevennessToPrint.close()

brayCurtisToPrint = open(str(output+'/brayCurtisData.txt'), 'w')
header = False
for i in taxonomyList:
	lineToPrint = i + '\t'
	if not header:
		headerToPrint = '#TYPE \t'
		for comparison in brayCurtis[i].keys():
			headerToPrint += comparison + '\t'
		headerToPrint = headerToPrint.strip()
		headerToPrint += '\r'
		brayCurtisToPrint.write(headerToPrint)
		header = True
	for j in brayCurtis[i]:
		lineToPrint += str(brayCurtis[i][j]) + '\t'
	lineToPrint = lineToPrint.strip()
	lineToPrint += '\r'
	brayCurtisToPrint.write(lineToPrint)
brayCurtisToPrint.close()

# ====================================

# Calculate percentage of each taxa that is a certain 'type', then print
printDictionary(byTypeList, 'TypeCounts')
printDictionary(byTaxaListDict, 'TaxaCounts')

byTypePercentage = percentageCalc(byTypeList)
byTaxaPercentage = percentageCalc(byTaxaListDict)

printDictionary(byTypePercentage,'TypePercentage')
printDictionary(byTaxaPercentage,'TaxaPercentage')

	
	

	
	
	
	