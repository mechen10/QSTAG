#!bin/bash/python

# Temp script to extract OTU ids with other info

modelBoundariesFP = "modelBoundaries_type.txt"
taxaIDLegendFP = "taxaIDLegend.txt"

################################## Functions

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
	tempContent = tempContent.split('\r')
	tempLine = []
	outputTaxaMap = {}
	for line in tempContent:
		tempLine = line.split('\t')
		outputTaxaMap[tempLine[0]] = tempLine[1]
	return outputTaxaMap
	
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
	
def changeBloom(dictionary):
	for taxa in dictionary:
		if dictionary[taxa]['bloom'] != 'No':
			dictionary[taxa]['type'] = 'bloom' + dictionary[taxa]['bloom']
	return dictionary
	
################################## Functions

	
modelBoundaries = loadMBIntoDictionary(modelBoundariesFP)
modelBoundaries = changeBloom(modelBoundaries)
taxaIDLegend = loadTaxaMap(taxaIDLegendFP)

OTUs = []
metadata = {}
for i in modelBoundaries.keys():
	OTUs.append(i)
	metadata[i] = {}
	metadata[i]['type'] = modelBoundaries[i]['type']
	metadata[i]['taxonomy'] = taxaIDLegend[i]
	
	
	

# tempFile = ''
# for line in modelBoundaries:
# 	tempFile += line
# 	
# tempFile = tempFile.strip()
# tempFile = tempFile.split('\r')
# 
# taxaIDtemp = ''
# for line in taxaIDLegend:
# 	taxaIDtemp += lin
# 
# taxaIDtemp = taxaIDtemp.strip()
# taxaIDtemp = taxaIDtemp.split('\r')
# 
# taxaIDdict = {}
# for line in taxaIDtemp:
# 	tempLine = line.split("\t")
# 	taxaIDdict[tempLine[0]] = tempLine[1]
# 
# OTUs = []
# for line in tempFile[1:]:
# 	splitLine = line.split('\t')
# 	OTUs.append(splitLine[0], splitLine)

toPrint = open('OTUsTaxaType.txt', 'w')
toWrite = ''
for i in OTUs:
	toWrite += i + '\r'
toWrite = toWrite.strip()
toPrint.write(toWrite)
toPrint.close()

printMetadata(metadata, 'metadataIDTypeTaxonomy.txt')