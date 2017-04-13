#!/bin/bash/python

# This script takes the modelBoundaries from two datasets, the taxamap, and names to make a comparison between shared OTUs
# the outputs are:
# metadataShareOTUs.txt
# SharedOTUs.txt
# SharedOTUsCond.txt
# SharedOTUsNames.txt



# Do optparse here
import argparse
import numpy as np


parser = argparse.ArgumentParser(
	description="Making simple table summary")
parser.add_argument(
	'-n',
	'--namesData',
	help = "Names of data set; comma separated",
	type = str,
	required = True,)
parser.add_argument(
	'-m',
	'--modelBoundariesFP',
	help = 'File Paths to modelBoundary files; must be in same order as namesData',
	type = str,
	required = True)
parser.add_argument(
	'-N',
	'--gradientNames',
	help = 'Gradient file; must be identical between two datasets',
	required = True)
parser.add_argument(
	'-t',
	'--taxamap',
	help = 'Taxa map for IDs',
	required = True)
	
		
args = parser.parse_args()

namesDatatmp = args.namesData
namesData = namesDatatmp.split(",")

modelBoundariestmp = args.modelBoundariesFP
PWDData = modelBoundariestmp.split(',')

gradientFP = args.gradientNames
gradienttemp = open(gradientFP, 'r')
tempfile = ''
for i in gradienttemp:
	tempfile += i
gradient = tempfile.split(',')

taxamapFP = args.taxamap


# Things I have to load
# namesData = ['Baltic','Fraser']
# PWDData = ['/Users/melissachen/Documents/BETATEST_13march2017/16sBaltic/modelBoundaries_type.txt','/Users/melissachen/Documents/BETATEST_13march2017/16sFraser/modelBoundaries_type.txt']
# gradient = ['fresh','brackish','marine','SalinityEnviron']
#  

#########################################################



def loadMBIntoDictionary(PWD):
	tempFile = open(PWD, 'r')
	tempContent = ''
	for line in tempFile:
		tempContent += line
	tempContent = tempContent.strip()
	tempContent = tempContent.split('\r')
	tempLine = []
	outputDictionary = {}
	first = True
	for line in tempContent:
		tempLine = line.split('\t')
		if first:
			headers = tempLine[1:]
			first = False
		else:
			outputDictionary[tempLine[0]] = {}
			otherLines = tempLine[1:]
			for n in range(len(headers)):
				if n > (len(otherLines)-1):
					outputDictionary[tempLine[0]][str(headers[n])] = []
				else:
					outputDictionary[tempLine[0]][str(headers[n])] = otherLines[n]
	return outputDictionary

def loadTaxaMap(PWD):
	tempFile = open(PWD, 'r')
	tempContent = ''
	for line in tempFile:
		tempContent += line
	tempContent = tempContent.strip()
	tempContent = tempContent.split('\n')
	tempLine = []
	outputTaxaMap = {}
	for line in tempContent:
		tempLine = line.split('\t')
		outputTaxaMap[tempLine[0]] = tempLine[1]
	return outputTaxaMap

def changeBloom(dictionary):
	for taxa in dictionary:
		if dictionary[taxa]['bloom'] != 'No':
			dictionary[taxa]['type'] = 'bloom' + dictionary[taxa]['bloom']
	return dictionary

def countType(dictionary, dataset, type):
	totalCount = 0
	for taxa in dictionary[dataset]:
		if dictionary[dataset][taxa]['type'] == type:
			totalCount += 1
	return totalCount

def countObs(dictionary, dataset):
	relAbundData = dictionary[dataset].copy()
	totalCount = 0
	for type in dictionary[dataset]:
		totalCount += dictionary[dataset][type]
	return totalCount
	
def makeRelAbund(dictionary, dataset):
	newDictionary = dictionary[dataset].copy()
	totalCount = dictionary[dataset]['TOTAL']
	for type in dictionary[dataset]:
		if type == 'TOTAL':
			pass
		else:
			newDictionary[type] = dictionary[dataset][type]/float(totalCount)
	return newDictionary

def printTable(dictionary,name):
	global listTypes
	first = True
	toWrite = ''
	printingTable = open(name,'w')
	for i in dictionary:
		if first == True:
			headers = dictionary[i].keys()
			toWrite += "TABLE\t"
			for j in listTypes:
				if headers.count(j) == 0:
					toWrite += "\t"
				else:
					toWrite += j + " (%)\t"
			toWrite += "TOTAL # OTUs"
			toWrite = toWrite.strip()
			toWrite += "\r"
			first = False
		toWrite += i + "\t" 
		for k in listTypes:
			if k == "SPACE":
				toWrite += "\t"
			elif listTypes.count(k) == 0:
				pass
			else:
				tempNumber = round(dictionary[i][k]*100,1)
				toWrite += str(tempNumber) + "\t"
		toWrite += str(dictionary[i]['TOTAL'])
		toWrite = toWrite.strip()
		toWrite += "\r"
	printingTable.write(toWrite)
	printingTable.close()
	return "TABLE PRINTED"
	
def printMatrix(x,y,matrix,headers,totalcounts,name):
	xTotalCounts = totalcounts[x]
	yTotalCounts = totalcounts[y]
	xTotalSum = 0
	yTotalSum = 0
	for a in xTotalCounts.keys():
		xTotalSum += xTotalCounts[a]
		yTotalSum += yTotalCounts[a]
	toWrite = ''
	toWrite += x+"vs"+y+"\t"
	# First row is going to be Fraser (y)
	for i in headers:
		n = yTotalCounts[i]
		toWrite += i + " ("+str(n)+")"+ "\t"
	toWrite += "TOTAL (" + str(yTotalSum) + ")"+ "\r"
	for j in range(len(headers)):
		headerTemp = headers[j]
		toWrite += headerTemp + " (" + str(xTotalCounts[headerTemp]) + ")" + "\t"
		for k in range(len(headers)+1):
			toWrite += str(matrix[j,k]) + "\t"
		toWrite = toWrite.strip()
		toWrite += "\r"
	toWrite += "TOTAL ("+ str(xTotalSum) + ")" + "\t"
	for l in range(len(headers)+1):
		toWrite += str(matrix[len(headers),l]) + "\t"
	toWrite = toWrite.strip()
	matrixToPrint = open(name, 'w')
	matrixToPrint.write(toWrite)
	matrixToPrint.close()
	return "TABLE PRINTED"

def printMatrixNames(x,y,matrix,headers,name):
	toWrite = ''
	toWrite = str(x) + "vs" + str(y)
	for i in headers:
		toWrite += "\t" + i
	toWrite += "\r"
	for i in range(len(matrix)):
		toWrite += headers[i] 
		for j in range(len(matrix[i])):
			toWrite += "\t" + matrix[i][j] 
		toWrite += "\r"
	matrixNamesToPrint = open(name, "w")
	matrixNamesToPrint.write(toWrite)
	matrixNamesToPrint.close()
	return "TABLE PRINTED"

def printMetadata(dict, name):
	toWrite = ''
	first = True
	for i in dict.keys():
		if first:
			toWrite += "TABLE"
			for j in dict[i].keys():
				toWrite += "\t" + j
			toWrite += "\r"
			first = False
		toWrite += i
		for k in dict[i].keys():
			toWrite += "\t" + dict[i][k]
		toWrite += "\r"
	metadataToPrint = open(name, "w")
	metadataToPrint.write(toWrite)
	metadataToPrint.close()
	return "TABLE PRINTED"

	
	
#########################################################
import sys


#########################################################
# BEGIN

if len(namesData) != len(PWDData):
	sys.exit("Different number of names of datasets and datasets themselves")
	
comboDictionary = {}
for i in range(len(namesData)):
	comboDictionary[namesData[i]] = changeBloom(loadMBIntoDictionary(PWDData[i]))
	
# Get list of types to make sure it's in the listTypes later

listTypesCheck = []
for dataset in comboDictionary:
	for taxa in comboDictionary[dataset]:
		listTypesCheck.append(comboDictionary[dataset][taxa]['type'])
listTypesCheckUnique = list(set(listTypesCheck))

listTypes = [gradient[0],gradient[1],gradient[2],'SPACE','half'+gradient[0],'lo'+gradient[1],'half'+gradient[1],'hi'+gradient[1],'half'+gradient[2],'SPACE','bloom'+gradient[0],'bloom'+gradient[1],'bloom'+gradient[2],'SPACE','noclass','inv'+gradient[1]]

# Check to make sure there aren't any weird names that the program doesn't recognize
if False in [True for x in listTypesCheckUnique if x in listTypes]:
	sys.exit("There are type names that are not recognized; check that the gradient file is correct")


taxamap = loadTaxaMap(taxamapFP)

print "DONE LOADING"

#########################################################

# Make table that shows counts/percentages of each site and type

# Make count table
binningTableCounts = {}
for dataset in comboDictionary:
	binningTableCounts[dataset]= {}
	for type in listTypes:
		if type == 'SPACE':
			pass
		else:
			binningTableCounts[dataset][type] = countType(comboDictionary,dataset,type)
			
# Count total for each and add as column
for dataset in binningTableCounts:
	binningTableCounts[dataset]['TOTAL'] = countObs(binningTableCounts, dataset)

# Make percentage version

binningTableRelAbund = binningTableCounts.copy()
for dataset in binningTableCounts:
	binningTableRelAbund[dataset] = makeRelAbund(binningTableCounts, dataset)
	
printTable(binningTableRelAbund, 'BinningSummary.txt')

print "DONE MAKING BINNING SUMMARY"

#########################################################

# Now, find out how many OTUs they have in common. ONLY WORKS WITH 2 SETS AT A TIME
# Make it do it twice; once with true OTUs and once with OTU names
# 
# 
# # Find all OTUs
# allCommonOTUs = list(set(comboDictionary[comboDictionary.keys()[0]]) & set(comboDictionary[comboDictionary.keys()[1]]))
# 
# # Find all common names
# namesAll = {}
# for dataset in namesData:
# 	namesAll[dataset] = []
# 	for i in comboDictionary[dataset].keys():
# 		keyPos = taxamap.keys().index(i)
# 		fullName = taxamap[taxamap.keys()[keyPos]]
# 		namesAll[dataset].append((fullName))
# 
# print namesAll


######################################################### OTU version

allCommonOTUs = list(set(comboDictionary[comboDictionary.keys()[0]]) & set(comboDictionary[comboDictionary.keys()[1]]))
# Find out total list of types, without spaces
listTypesNospace = [x for x in listTypes if x != "SPACE"]

# make internal dictionary for total counts
tempIntDict = {}
for i in listTypesNospace:
	tempIntDict[i]= 0

# Count total counts of each
dictTotalTypeCounts = {}
for i in comboDictionary.keys():
	dictTotalTypeCounts[i] = tempIntDict.copy()
	for j in comboDictionary[i].keys():
		typeTemp = listTypesNospace[listTypesNospace.index(comboDictionary[i][j]['type'])]
		dictTotalTypeCounts[i][typeTemp] += 1


# Make an empty table to fill with all overlapping OTUs
commonMat = np.zeros((len(listTypesNospace)+1, len(listTypesNospace) +1))


for i in allCommonOTUs:
	firstDictTypeTemp = comboDictionary[comboDictionary.keys()[0]][i]['type']
	secondDictTypeTemp = comboDictionary[comboDictionary.keys()[1]][i]['type']
	xPosTemp = listTypesNospace.index(firstDictTypeTemp)
	yPosTemp = listTypesNospace.index(secondDictTypeTemp)
	commonMat[xPosTemp,yPosTemp] += 1

for i in range(len(listTypesNospace)+1):
	commonMat[i,len(listTypesNospace)] = sum(commonMat[i,0:len(listTypesNospace)])
	commonMat[len(listTypesNospace),i] = sum(commonMat[0:len(listTypesNospace),i])

print "Printing SharedOTUs.txt"
printMatrix(comboDictionary.keys()[0],comboDictionary.keys()[1],commonMat, listTypesNospace, dictTotalTypeCounts, "SharedOTUs.txt")


#### DO WITH CONDENSED

listTypesCond = [gradient[0],gradient[1],gradient[2],"noclass"]

# Make temp dict for counts
tempIntDictCond = {}
for i in listTypesCond:
	tempIntDictCond[i] = 0

# Get number in each dataset of total types, to compared with number shared
dictTotalTypeCountsCond = {}
for i in comboDictionary.keys():
	dictTotalTypeCountsCond[i] = tempIntDictCond.copy()
	for j in comboDictionary[i].keys():
		typeTemp = listTypesCond[listTypesCond.index(comboDictionary[i][j]['typeSimple'])]
		dictTotalTypeCountsCond[i][typeTemp] += 1

commonMatCond = np.zeros((len(listTypesCond)+1, len(listTypesCond) +1))

for i in allCommonOTUs:
	firstDictTypeTemp = comboDictionary[comboDictionary.keys()[0]][i]['typeSimple']
	secondDictTypeTemp = comboDictionary[comboDictionary.keys()[1]][i]['typeSimple']
	xPosTemp = listTypesCond.index(firstDictTypeTemp)
	yPosTemp = listTypesCond.index(secondDictTypeTemp)
	commonMatCond[xPosTemp,yPosTemp] += 1

for i in range(len(listTypesCond)+1):
	commonMatCond[i,len(listTypesCond)] = sum(commonMatCond[i,0:len(listTypesCond)])
	commonMatCond[len(listTypesCond),i] = sum(commonMatCond[0:len(listTypesCond),i])

print "Printing SharedOTUsCond.txt"

printMatrix(comboDictionary.keys()[0],comboDictionary.keys()[1],commonMatCond, listTypesCond, dictTotalTypeCountsCond, "SharedOTUsCond.txt")


#### MAKE MATRIX WITH OTU NAMES- ALL

# Make an empty table to fill with OTU names
# commonMatNames = np.empty([len(listTypesNospace)+1, len(listTypesNospace) +1], dtype=object)
# for i in range(len(commonMatNames)):
# 	commonMatNames[i] = []

commonMatNames = []
for i in range(len(listTypesNospace)):
	commonMatNames.append([])
	for j in range(len(listTypesNospace)):
		commonMatNames[i].append('')
	
# Find all OTUs
# already called allCommonOTUs

for i in allCommonOTUs:
	firstDictTypeTemp = comboDictionary[comboDictionary.keys()[0]][i]['type']
	secondDictTypeTemp = comboDictionary[comboDictionary.keys()[1]][i]['type']
	xPosTemp = listTypesNospace.index(firstDictTypeTemp)
	yPosTemp = listTypesNospace.index(secondDictTypeTemp)
	if commonMatNames[xPosTemp][yPosTemp] == '':
		commonMatNames[xPosTemp][yPosTemp] = str(i)
	else:
		commonMatNames[xPosTemp][yPosTemp] += "," + str(i)

print "Printing SharedOTUsNames.txt"

printMatrixNames(comboDictionary.keys()[0],comboDictionary.keys()[1],commonMatNames, listTypesNospace, "SharedOTUsNames.txt")


#### Make list of shared OTUs, but with percent representation as well

metadataCommonOTUs = {}

for i in allCommonOTUs:
	metadataCommonOTUs[i] = {}
	firstDataSet = comboDictionary.keys()[0]
	secondDataSet = comboDictionary.keys()[1]
	metadataCommonOTUs[i][firstDataSet] = comboDictionary[firstDataSet][i]['type']
	metadataCommonOTUs[i][secondDataSet] = comboDictionary[firstDataSet][i]['type']
	taxamap[i]
	metadataCommonOTUs[i]['taxaID'] = taxamap[i]

print "Printing metadataShareOTUs.txt"

printMetadata(metadataCommonOTUs, 'metadataShareOTUs.txt')
	


######################################################### Names version

# Find all common names
# namesAll = {}
# for dataset in namesData:
# 	namesAll[dataset] = []
# 	for i in comboDictionary[dataset].keys():
# 		keyPos = taxamap.keys().index(i)
# 		fullName = taxamap[taxamap.keys()[keyPos]]
# 		namesAll[dataset].append(set(fullName))
# 		
# allCommonNames = list(namesAll[namesData[0]] & namesAll[namesData[1]])
# print(len(allCommonNames))
# 
# 		
# commonMatNames = np.zeros((len(listTypesNospace)+1, len(listTypesNospace) +1))
# 
# for i in allCommonNames:
# 	ti
# 	firstDictTypeTemp = comboDictionary[comboDictionary.keys()[0]][i]['type']
# 	secondDictTypeTemp = comboDictionary[comboDictionary.keys()[1]][i]['type']
# 	xPosTemp = listTypesNospace.index(firstDictTypeTemp)
# 	yPosTemp = listTypesNospace.index(secondDictTypeTemp)
# 	commonMat[xPosTemp,yPosTemp] += 1
# 
# for i in range(len(listTypesNospace)+1):
# 	commonMat[i,len(listTypesNospace)] = sum(commonMat[i,0:len(listTypesNospace)])
# 	commonMat[len(listTypesNospace),i] = sum(commonMat[0:len(listTypesNospace),i])

