import statistics
import regex
import sys
import os
import re
from Bio import SeqIO
from Bio import Entrez
from Bio.Align.Applications import ClustalOmegaCommandline
import numpy as np
import matplotlib.pyplot as plt


# Opens FASTA file and stores them in headerAndSequence, key = header value = sequence
# edit: Great spot for sequence length edit
def fastaFormatFileInput(multiFASTA):
    currentHeader = ''

    # Opens FASTA file and stores them in headerAndSequence
    with open(multiFASTA) as a_file:

        for line in a_file:

            strippedLine = line.strip()

            if len(strippedLine) == 0:  # skips empty lines

                continue

            # Uses sequence header as key
            elif strippedLine.startswith(tuple(headerIdentifiers)):

                headerAndSequence[strippedLine] = ''

                currentHeader = strippedLine

            # if strippedLine is not empty or header it is sequence and will be appended
            else:

                headerAndSequence[currentHeader] += strippedLine

    for key, value in headerAndSequence.items():
        # print(len(value))
        if endPositionofSurvey > len(str(value)):
            # print("End position is greater than sequence length. Please choose a different value")
            break

        elif startPositionofSurvey > 0 and endPositionofSurvey == 0:
            headerAndSequence[key] = [value[startPositionofSurvey:]]


        elif startPositionofSurvey and endPositionofSurvey > 0:

            headerAndSequence[key] = [value[startPositionofSurvey:endPositionofSurvey]]

        else:
            headerAndSequence[key] = [value]

        # used to ensure .gb file is read correctly
        headerAndSequenceFullLength[key] = [value]

    # print(headerAndSequence)
    # outputComparison['ENTRY COUNT'] = len(headerAndSequence.keys())
    a_file.close()


# Links .gb files with sequences in FASTA using sequence as identifier
def genBankFileLink(dict, path, emptyDict):
    fastaDict = dict
    keyAndAccessionDict = {}

    # first parameter is only used in myco analysis, operational without as long as .gb files are provided
    # with open(direc) as b_file:
    #
    #     for line in b_file:
    #         lineAsList = (line.strip()).split("\t")
    #
    #         genBankFileName = lineAsList[0].upper()
    #
    #         accessionNumber = lineAsList[1]
    #
    #         keyAndAccessionDict[genBankFileName] = accessionNumber
    #
    # b_file.close()

    for file in os.listdir(path):

        with open(os.path.join(path, file)) as gbFile:
            fullGBText = SeqIO.read(gbFile, "genbank")
            # print(type(fullGBText.seq))
            # print(type(fastaSeq))
            # print(len(fastaSeq))
            # print(len(fullGBText.seq))
            # print(str(fullGBText.seq))
            # print(str(fastaSeq))
            # print(str(fastaSeq))

        for k, v in fastaDict.items():
            # print(v)
            fastaName = str(k)

            if fastaName in headerAndSequenceFullLength.keys():

                fastaSeq = headerAndSequenceFullLength[fastaName]

                fastaSeq = str(fastaSeq)[2:-2]

                # print(fastaSeq)
                # print(fullGBText.seq)

            if str(fastaSeq) == str(fullGBText.seq):
                # print("yes")
                genBankDict[fastaName] = [fullGBText.description, fullGBText]


# usefull to change set() to list[] while keeping variable order
def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]


# variables required for function, DO NOT EDIT.
overLapDictionary = {}
overLapDictionaryJustNames = {}
overLapDictionaryNumbers = {}
motifOccurancesDict = {}
histogramMotifPos = {}
headerAndSequence = {}
headerAndSequenceAllEntries = {}
headerIdentifiers = ['>', '.', '+']
sequencesWithoutMotif = set()
genBankDict = {}
totalSequenceCount = 0
countOfSequencesMatchingMotif = 0
headerAndSequenceFullLength = {}


#used to debug
currentCluster = 'A12'


######OPTIONAL USER INPUT VARIABLES###################
#############################################

#If a motif is abundant and users do not want it inclued they can place in in this list and seqSearch will skip.
#This is helpful for repeats "ie AAAA", but by default no motifs are ommited
#Additionally, if users want to explore variable motifs excluding a few specific instances they can be placed here
motifsToIgnore = []

#If required, users can choose to only explore indicated regions of the inputer sequences.
#If nothing is entered the full sequence of each entry will be surveyed
#Default value of both is 0

startPositionofSurvey = 0
endPositionofSurvey = 0

if endPositionofSurvey > startPositionofSurvey:
    print("Start position must be GREATER than end position.")

#Do the entered motifs form distinct groups?
#By default, this is selected as 'N', but this parameter reports at the top of the summary file
#Which sequences can be placed in multiple motif groups

motifsFormDistinctGroups = 'Y'

######REQUIRED USER INPUT VARIABLES###################
#############################################

# enter all motifs of interest here. This is a required function, unlike in motifSurvey where it is optional.
regexList = ['TTACGA[AG]TCA', 'GTGCGATGTCAA', 'GTCCGTTGTCA',
             'ACGTAGTCAA', 'GGGAT[AT]GTCAA', '[ATCG][ATCG]TG[ATCG][ATCG]TGTCAA', 'TTC[ATCG][ATCG]TGTCAA',
             'GGGG[GT]ATGTCAA', 'TACGGTGTCCAA', 'GTAC[AG]GGGTCAA']

# regexList = ['TGA[ATCG]TCGTAA','TGATTCGTAA']

# regexList = ['TTGACA[ATCG][ATCG]CA[ATCG][ATCG]', 'AACTGT[ATCG][ATCG]GT[ATCG][ATCG]']

regexList = ['TTGACATCGCAC', 'AACTGTAGCGTG']

# #MPME #1 & #2 test
# regexList = ['TTATC[AT]GGGGT', 'ACCTCAGATAA']

#used to store current file name for comparison
nameOfPath = ''

# enter a name for summary file
fileName = nameOfPath + '_output.txt'

# enter a prefex for additional outputs
filePrefix = nameOfPath + 'FS'

# directory to output file
outputDirectory = 'C:/Users/Jason Thompson/Desktop/seqSearch/finalOutputs/scrapCheck/debugMotif/'

# path to FASTA file(s), function will detect whether it is a file or directory.
# If this is a directory, please ensure only files pertaining to analysis are present
# This will be changed in future seqSearch updates.
pathToFasta = 'C:/Users/Jason Thompson/Desktop/seqSearch/debugFiles/A4OutputFile.txt'

# path to GenBank file directory
# Please ensure only files pertaining to analysis are present
# This will be changed in future seqSearch updates.
pathToGenBankDirectory = 'C:/Users/Jason Thompson/Desktop/seqSearch/finalOutputs/scrapCheck/allGBFiles'

# minimum Motifs to be inclued in a motif group
minimumMotifs = 4

# input max sequence length of dataset here
maxSeqLength = 60000

# input desired bin size here
binSize = 1000

# this is required for alignments. Please provide the path to the clustal omega .exe file here.
pathToClustalOmega = 'C:/Users/Jason Thompson/Desktop/seqSearch/clustal-omega-1.2.2-win64/clustalo.exe'


# Variables For Bar Graph

xAxisLabel = 'Sequence position (bp)'

yAxisLabel = 'Number Of Motifs Found At Sequence Position On X-Axis'

#############################################
#############################################


# All motifs in regexList

nameOfPath = os.path.basename(pathToFasta)
print(nameOfPath)
for entry in regexList:
    overLapDictionary[entry] = []
    overLapDictionaryJustNames[entry] = []
    overLapDictionaryNumbers[entry] = []
    # Makes regex expression into a re.Pattern object, used later to search using re.finditer
    r = regex.compile('%s' % entry)
    motifName = entry
    print(entry)

    if os.path.isfile(pathToFasta):

        numberOfMatches = 0
        fastaFormatFileInput(pathToFasta)
        genBankFileLink(headerAndSequence, pathToGenBankDirectory, genBankDict)
        nameOfPath = os.path.basename(pathToFasta)
        ifMatchMadeEqualsOne = 0

        # Searches every sequence for sites that match motif
        for h, v in headerAndSequence.items():

            listOfPositions = []
            sequencesInDict = len(headerAndSequence.keys())
            sequence = str(v)

            # test section, used for phylo analysis not part of motifCheck functionality
            #######################

            # finds position of all matches in sequence
            for match in regex.finditer(r, sequence):

                if match in motifsToIgnore:
                    continue
                else:
                    s = match.start()

                    e = match.end()

                    subSeqKey = v[s:e]

                    matchAsRange = "(" + str(s) + "," + str(e) + ")"

                    listOfPositions.append(matchAsRange)

            nameForInput = h + "\t" + str(len(listOfPositions)) + "\t" + str(listOfPositions)

            # If value minimumMotifs or more sites are detected, sequences is included as a match
            if len(listOfPositions) > minimumMotifs:
                numberOfMatches = numberOfMatches + 1
                overLapDictionary[entry].append(str(nameForInput))
                overLapDictionaryJustNames[entry].append(str(h))
                sequencesWithoutMotif.add(str(h))

        inputForSplitFunction = nameOfPath + ":" + str((numberOfMatches)) + ":" + str((sequencesInDict))
        countOfSequencesMatchingMotif = countOfSequencesMatchingMotif + numberOfMatches
        # Stores (num. Sequences with Matches)/(total sequences) for each cluster
        overLapDictionaryNumbers[entry].append(inputForSplitFunction)
        # A dictionary containing all sequences in all clusters
        headerAndSequenceAllEntries.update(headerAndSequence)
        # resets dictionary for next group
        headerAndSequence = {}



    elif os.path.isdir(pathToFasta):
        for file in os.listdir(pathToFasta):
            print(file)
            pathToFile = pathToFasta + file
            nameOfPath = str(file)
            numberOfMatches = 0
            fastaFormatFileInput(pathToFile)
            genBankFileLink(headerAndSequence, pathToGenBankDirectory, genBankDict)
            ifMatchMadeEqualsOne = 0

            # Searches every sequence for sites that match motif
            for h, v in headerAndSequence.items():

                listOfPositions = []
                sequencesInDict = len(headerAndSequence.keys())
                sequence = str(v)

                # finds position of all matches in sequence
                for match in regex.finditer(r, sequence):
                    s = match.start()

                    e = match.end()

                    subSeqKey = v[s:e]

                    matchAsRange = "(" + str(s) + "," + str(e) + ")"

                    listOfPositions.append(matchAsRange)

                nameForInput = h + "\t" + str(len(listOfPositions)) + "\t" + str(listOfPositions)

                # If value minimumMotifs or more sites are detected, sequences is included as a match
                if len(listOfPositions) > minimumMotifs:
                    numberOfMatches = numberOfMatches + 1
                    overLapDictionary[entry].append(str(nameForInput))
                    overLapDictionaryJustNames[entry].append(str(h))
                    sequencesWithoutMotif.add(str(h))

            inputForSplitFunction = nameOfPath + ":" + str((numberOfMatches)) + ":" + str((sequencesInDict))
            countOfSequencesMatchingMotif = countOfSequencesMatchingMotif + numberOfMatches
            # Stores (num. Sequences with Matches)/(total sequences) for each cluster
            overLapDictionaryNumbers[entry].append(inputForSplitFunction)
            # A dictionary containing all sequences in all clusters
            headerAndSequenceAllEntries.update(headerAndSequence)
            # resets dictionary for next group
            headerAndSequence = {}




original_stdout = sys.stdout


#This section handles all file output, file names bellow can be changed as users see fit.

#Displays the sequences and positions each motif was found
with open(outputDirectory + "motifPositionsInSequence.txt",'w') as q:

    sys.stdout = q

    for k,v in sorted(overLapDictionary.items()):

        motifStatsCount = []
        q = set(v)
        value = sorted(list(q))
        print("Motif: " + k + ":\t" + "SEQUENCE WITH MOTIF: " + str(len(value)))
        print("____________________________________")
        print()
        for item in value:
            itemAsList = item.split("\t")
            #print(itemAsList[1])
            motifStatsCount.append(itemAsList[1])
            # print(len(itemAsList))
            # testList = itemAsList[2]
            # testList = testList.lstrip("['")
            # testList = testList.rstrip("']")
            # testList = testList.replace("'","")
            # finalList = [testList]
            # print(type(finalList))
            # print(finalList)
            print(itemAsList[0])
            print("TOTAL OCCURRENCE IN SEQUENCE: " + str(itemAsList[1]))
            # print(itemAsList[2])
            # print(type(itemAsList[2]))
            positionForOutput = itemAsList[2].split("',")
            # print(type(positionForOutput))
            # print(positionForOutput)

            for numValue in positionForOutput:
                testName = numValue
                #print(testName)
                justNum = re.findall(r'\d+', testName)
                justNum = [int(x)+startPositionofSurvey for x in justNum]
                print(justNum)


            print()
            motifOccurancesDict[k] = motifStatsCount
        print()
        print("###END OF MOTIF###")
        print()


#Provides a summary of all motifs. How many sequences contained motifs, which files contained them,
# average found per sequence etc.
with open(outputDirectory + "summaryOfMotifsFound.txt",'w') as f:

    sys.stdout = f
    finalOverlapDict = {}

    #makes sure fraction isn't > 100%
    if countOfSequencesMatchingMotif < 0:

        matchVsTotal = str(countOfSequencesMatchingMotif) + "/" + str(len(headerAndSequenceAllEntries.keys())) + \
                   " SEQUENCES WERE FOUND TO HAVE INPUT MOTIFS"

    elif countOfSequencesMatchingMotif >= len(headerAndSequenceAllEntries.keys()):
        matchVsTotal = str(len(headerAndSequenceAllEntries.keys())) + "/" + \
                       str(len(headerAndSequenceAllEntries.keys())) + \
                       " SEQUENCES WERE FOUND TO HAVE INPUT MOTIFS"

        print(matchVsTotal)
        print()

    #Compares dictionary for sequences that are found in >1 of the input motifs
    for y, v in overLapDictionaryJustNames.items():

        for j, t in overLapDictionaryJustNames.items():

            overlappingSequences = set(v).intersection(set(t))
            #print(len(overlappingSequences))

            duplicateCheck = j + ":" + y

            #skips same key
            if y == j:

                continue

            #as dict. comapres to self, >1 is needed to avoid its own inclusion
            elif len(overlappingSequences) > 1 and duplicateCheck not in finalOverlapDict.keys():

                finalOverlapDict[y + ":" + j] = overlappingSequences

    if motifsFormDistinctGroups == 'Y':
        print("BELOW ARE ALL SEQUENCES THAT CAN BE PLACED IN TWO OR MORE GROUPS")
        for k, v in finalOverlapDict.items():
            print(k + "\t" + str(len(v)))
            for item in v:
                print(item)

    print()
    print("#####################")
    print()
    print(str())
    entryInFile = 0
    # print (num. Sequences with Matches)/(total sequences) for each cluster
    for k, v in overLapDictionaryNumbers.items():
        entryInFile += 1
        print(str(entryInFile) + ")" + "\t" + k)
        # print(v)
        countOfSequences = 0
        for entry in v:
            # print(entry)
            valueAsList = entry.split(":")
            # print(valueAsList)
            if int(valueAsList[1]) > 0:
                print(valueAsList[0] + ":\t" + valueAsList[1] + "/" + valueAsList[2])
                # Total matches
                countOfSequences = countOfSequences + int(valueAsList[1])
        print("____________________________")
        print("SEQUENCE COUNT:" + str(countOfSequences))
        print("____________________________")

        print()
        print("#####################")
        print()


    for k,v in motifOccurancesDict.items():
        print(k)
        print("SEQUENCES WITH MOTIF:\t" + str(len(v)))
        valueAsInt = []
        for item in v:
            valueAsInt.append(int(item))
        #print(len(valueAsInt))
        print("____________")
        print("MEAN: " + str(statistics.fmean(valueAsInt)))
        print("MEDIAN: " + str(statistics.median(valueAsInt)))
        print("STDEV: " + str(statistics.pstdev(valueAsInt)))
        print("MIN: " + str(min(valueAsInt)))
        print("MAX: " + str(max(valueAsInt)))
        print("\n")


    print("SEQUENCES NOT MATCHING ANY MOTIF")

    sequencesNotIncluded = set(sequencesWithoutMotif).symmetric_difference(set(headerAndSequenceAllEntries))
    for i in sequencesNotIncluded:
        print(i)


#Uses user input max sequence length and bin size to view motif matches in each bin
#Each bin has its own attributes such as '% in NC regions', 'Found in x/Total Sequences' etc.
#Using the same bins as ‘motifPositionsInBins’, for each bin all genes’ motifs were found in are output in file
# *_genesWithMotif.
with open(outputDirectory + "motifPositionsInBins.txt",'w') as x:

    sys.stdout = x

    #edit: move to function input, keep first value and move 60k and 1k.
    bins = list(range(0,maxSeqLength,binSize))
    #print(len(bins))




    for k,v in sorted(overLapDictionary.items()):

        # print(k)
        # print(v)

        motifStatsCount = []
        binDict = {}
        binDictSequences = {}
        binDictPos = {}
        motifsWithGeneOutput = {}
        motifsWithGeneOutputCount = {}
        motifName = str(k)
        sequencesWithGB = set()
        tv = set(v)

        #creates bins necessary for comparison, each sequences positions are compared to input bin values
        for value in range(len(bins) - 1):
            binRangeAsKey = str(int(bins[value]) + 1) + "," + str(bins[value + 1])
            binDict[binRangeAsKey] = 0
            binDictSequences[binRangeAsKey] = set()
            binDictPos[binRangeAsKey] = 0
            motifsWithGeneOutput[binRangeAsKey] = []

        valueEdit = sorted(list(tv))
        print("Motif: " + k + ":\t" + "SEQUENCES WITH MOTIF: " + str(len(valueEdit)))
        print("____________________________________")
        print("BIN\tTot. Occur\tSeq w/ Motif \t% Motif in Non-Coding")
        print()
        countOfNonCodingInverse = 0
        genBankNoEntry = 0
        sequencesWithoutGB = set()

        #the sequences each motif is found in is used here. Each item is a sequence and its list of attributes
        for item in valueEdit:
            itemAsList = item.split("\t")
            # print(itemAsList[1])
            motifStatsCount.append(itemAsList[1])
            testList = itemAsList[2]
            sequenceName = itemAsList[0]
            testList = testList.lstrip("['")
            testList = testList.rstrip("']")
            testList = testList.replace("'", "")
            testList = testList.replace("(", "")
            testList = testList.replace(")", "")
            finalListRaw = testList.split(",")
            # Reads list in from file
            # print(type(finalList))
            # print(finalList)
            # print(type(finalList))
            # print(item + "\n")
            motifOccurancesDict[k] = motifStatsCount
            # print(sequenceName)

            # print(sequenceName)
            # print(finalListRaw)
            finalList = []
            #only uses the first value os each motif range
            for item in range(len(finalListRaw)):
                if item % 2 != 0:
                    continue
                else:
                    finalList.append(finalListRaw[item])
            # print(finalList)
            for pos in finalList:
                # print(str(pos))
                #if users enter own start and end points of sequence to survey, the if/else here checks that
                if startPositionofSurvey == 0 and endPositionofSurvey == 0:
                    for k,v in binDict.items():
                        keyAsInt = k.split(',')
                        rangeValOne = int(keyAsInt[0])
                        rangeValTwo = int(keyAsInt[1])
                        if int(pos) in range(rangeValOne, rangeValTwo):
                            # print('1')
                            if sequenceName in genBankDict:
                                # print('1')
                                binDict[k] += 1
                                binDictSequences[k].add(sequenceName)
                                sequenceRecord = genBankDict[sequenceName][1]
                                matchesIngenBank = 0
                                seqRecordPosAsList = []
                                inverseSeqRecordPos = []
                                sequencesWithGB.add(sequenceName)
                                #searches .gb record for positions
                                for record in range(len(sequenceRecord.features)):

                                    locationRangeString = str(sequenceRecord.features[record].location)[1:-4]

                                    rangeAsList = locationRangeString.split(":")

                                    if re.search('^(oin)', str(rangeAsList[0])):

                                        continue

                                    elif int(rangeAsList[0]) == 0:

                                        continue

                                    if len(sequenceRecord.features[record].qualifiers) > 2 and int(rangeAsList[0]) != 0:
                                        seqRecordPosAsList.append(str(rangeAsList[0] + "," + rangeAsList[1]))

                                #used to make note of regions not included in .gb records (non-coding/non-transcribable)
                                for item in range(len(seqRecordPosAsList) - 1):
                                    indexPos = seqRecordPosAsList[item].split(',')
                                    indexPosNextIter = seqRecordPosAsList[item + 1].split(',')

                                    firstPosOfGap = int(indexPos[1]) + 1
                                    secondPosOfGap = int(indexPosNextIter[0]) - 1

                                    inverseSeqRecordPos.append(str(firstPosOfGap) + ',' + str(secondPosOfGap))



                                for record in range(len(sequenceRecord.features)):

                                    locationRangeString = str(sequenceRecord.features[record].location)[1:-4]

                                    rangeAsList = locationRangeString.split(":")

                                    currentProduct = str(sequenceRecord.features[record].qualifiers.get('product', ''))

                                    currentProduct = currentProduct[2:-2].strip()

                                    if re.search('^(oin)', str(rangeAsList[0])):

                                        continue

                                    elif int(rangeAsList[0]) == 0:

                                        continue

                                    #Checks each .gb record to see if motif position occurs in range of record
                                    elif int(pos) in range(int(rangeAsList[0]), int(rangeAsList[1])) and len(
                                            sequenceRecord.features[record].qualifiers) > 2:
                                        binDictPos[k] += 1
                                        recordTest = sequenceRecord.features[record]
                                        motifsWithGeneOutput[k].append(recordTest)


                                            #Keeps track of which positions occur in gene
                                        if currentProduct not in motifsWithGeneOutputCount.keys():

                                            motifsWithGeneOutputCount[currentProduct] = 1

                                        elif currentProduct in motifsWithGeneOutputCount.keys():

                                            motifsWithGeneOutputCount[currentProduct] += 1






                                for r in inverseSeqRecordPos:
                                    rangeSplit = r.split(',')
                                    if int(pos) in range(int(rangeSplit[0]), int(rangeSplit[1])):
                                        countOfNonCodingInverse += 1

                            elif sequenceName not in genBankDict:
                                    sequencesWithoutGB.add(sequenceName)

                elif startPositionofSurvey > 0:
                    posEdit = int(pos) + startPositionofSurvey
                    for k,v in binDict.items():
                        keyAsInt = k.split(',')
                        rangeValOne = int(keyAsInt[0])
                        rangeValTwo = int(keyAsInt[1])

                        if int(posEdit) in range(rangeValOne, rangeValTwo):
                            if itemAsList[0] in genBankDict:
                                binDict[k] += 1
                                binDictSequences[k].add(itemAsList[0])
                                sequenceRecord = genBankDict[sequenceName][1]
                                matchesIngenBank = 0
                                seqRecordPosAsList = []
                                inverseSeqRecordPos = []
                                sequencesWithGB.add(sequenceName)
                                #searches .gb record for positions
                                for record in range(len(sequenceRecord.features)):

                                    locationRangeString = str(sequenceRecord.features[record].location)[1:-4]

                                    rangeAsList = locationRangeString.split(":")

                                    if re.search('^(oin)', str(rangeAsList[0])):

                                        continue

                                    elif int(rangeAsList[0]) == 0:

                                        continue

                                    if len(sequenceRecord.features[record].qualifiers) > 2 and int(rangeAsList[0]) != 0:
                                        seqRecordPosAsList.append(str(rangeAsList[0] + "," + rangeAsList[1]))

                                #used to make note of regions not included in .gb records (non-coding/non-transcribable)
                                for item in range(len(seqRecordPosAsList) - 1):
                                    indexPos = seqRecordPosAsList[item].split(',')
                                    indexPosNextIter = seqRecordPosAsList[item + 1].split(',')

                                    firstPosOfGap = int(indexPos[1]) + 1
                                    secondPosOfGap = int(indexPosNextIter[0]) - 1

                                    inverseSeqRecordPos.append(str(firstPosOfGap) + ',' + str(secondPosOfGap))



                                for record in range(len(sequenceRecord.features)):

                                    locationRangeString = str(sequenceRecord.features[record].location)[1:-4]

                                    rangeAsList = locationRangeString.split(":")

                                    currentProduct = str(sequenceRecord.features[record].qualifiers.get('product', ''))

                                    currentProduct = currentProduct[2:-2].strip()

                                    if re.search('^(oin)', str(rangeAsList[0])):

                                        continue

                                    elif int(rangeAsList[0]) == 0:

                                        continue

                                    #Checks each .gb record to see if motif position occurs in range of record
                                    elif int(posEdit) in range(int(rangeAsList[0]), int(rangeAsList[1])) and len(
                                            sequenceRecord.features[record].qualifiers) > 2:
                                        binDictPos[k] += 1
                                        recordTest = sequenceRecord.features[record]
                                        motifsWithGeneOutput[k].append(recordTest)


                                            #Keeps track of which positions occur in gene
                                        if currentProduct not in motifsWithGeneOutputCount.keys():

                                            motifsWithGeneOutputCount[currentProduct] = 1

                                        elif currentProduct in motifsWithGeneOutputCount.keys():

                                            motifsWithGeneOutputCount[currentProduct] += 1






                                for r in inverseSeqRecordPos:
                                    rangeSplit = r.split(',')
                                    if int(posEdit) in range(int(rangeSplit[0]), int(rangeSplit[1])):
                                        countOfNonCodingInverse += 1

                            elif itemAsList[0] not in genBankDict:
                                    sequencesWithoutGB.add(sequenceName)

        totalMotifOccur = 0
        totalBins = []
        # print(str(binDict))
        #Only print bins with values to keep output concise
        for k,v in binDict.items():

            # print(str(k))
            # print(str(v))

            if int(v) > 0:
                # print(v)
                NCbinValue = round((1-int(binDictPos[k])/int(v))*100, 2)
                print(str(k) + "\t" + str(v) + "\t" + str(len(binDictSequences[k])) + "/" + str(len(valueEdit) - len(sequencesWithoutGB))+ "\t" +
                      str(NCbinValue))
                totalMotifOccur += v
                totalBins.append(str(k))

        print()
        print("TOTAL MOTIFS:" + str(totalMotifOccur))
        print("TOTAL BINS:" + str(len(totalBins)))

        if totalMotifOccur > 0:
            percentMotifsNC = countOfNonCodingInverse/totalMotifOccur
            # print(countOfNonCodingInverse)
              # print(totalMotifOccur)
            print("PERCENT MOTIFS FOUND IN NC: "+ str(percentMotifsNC*100))

        if genBankNoEntry > 0:
            print("MOTIFS FOUND OUT OF RANGE: " + str(round(int(genBankNoEntry)/int(totalMotifOccur), 2)))

        if len(sequencesWithoutGB) > 0:
            print("SEQUENCES WITHOUT .gb FILE: " + str(len(sequencesWithoutGB)))


        for k,v in binDict.items():
            binDict[k] = 0

        for k,v in binDictSequences.items():
            binDictSequences[k] = set()

        for k,v in binDictPos.items():
            binDictPos[k] = 0

        outFile = open(outputDirectory + str(motifName) +
                      "_genesWithMotifs.txt", 'w')
        #do count of genes here

        for k,v in motifsWithGeneOutput.items():
            if len(v) > 0:
                outFile.write("BIN:" + k + "\n\n")
                for item in v:
                    outFile.write(str(item) + "\n")

        outFile.write("\n\n" + "GENES ENTRIES WERE FOUND IN: \n\n")

        totalGeneProducts = 0
        totalDNABinding = 0
        totalHypoProteins = 0

        # print(str(motifsWithGeneOutputCount))

        for k, v in motifsWithGeneOutputCount.items():

            outFile.write(str(k) + ": " + str(v) + "\n")

            totalGeneProducts += v
            if k == 'DNA binding protein':
                totalDNABinding = v
            if k == 'hypothetical protein':
                totalHypoProteins = v

        totalWithoutHypo = totalGeneProducts - totalHypoProteins

        # outFile.write("RAW VALUES" + "\n")
        # outFile.write(str(totalWithoutHypo) + "\n")
        # outFile.write(str(totalHypoProteins) + "\n")
        # outFile.write(str(totalDNABinding) + "\n")
        outFile.write("\n")

        outFile.write("TOTAL GENE PRODUCTS: " + str(totalGeneProducts) + "\n")
        outFile.write("SEQUENCES WITH .gb FILE: " + str(len(sequencesWithGB)) + "\n")
        outFile.write("\n")

        if totalGeneProducts > 0:
            # outFile.write("WHEN ALL PROTEINS ARE CONSIDERED, 'DNA binding protein' is: " + str(
            #     round(totalDNABinding / totalGeneProducts, 3) * 100) + "% OF TOTAL GENE ENTRIES")
            # outFile.write("\n")
            if totalWithoutHypo > 0:
                outFile.write("WHEN ONLY NON-HYPOTHETICAL PROTEINS ARE CONSIDERED, 'DNA binding protein' is: " + str(
                round(totalDNABinding / totalWithoutHypo, 3) * 100) + "% OF TOTAL GENE ENTRIES")
                outFile.write("\n")

        outFile.close()
        print()
        print("###END OF MOTIF###")
        print()
        # p.close()




    x.close()


#For each input motif, the average proportion of each sequence that is non-transcribable is output
with open(outputDirectory + "nonCodingofGB.txt",'w') as z:

    sys.stdout = z
    #key is motif, value is a list of all sequences and positions found in
    for k,v in sorted(overLapDictionary.items()):

        motifStatsCount = []
        # nonCodingWithMotifs = []
        # nonCodingWithoutMotifs = []
        allNCDict = {}
        currentMotif = k
        tv = set(v)
        value = sorted(list(tv))
        #Creates a section for each motif
        print("Motif: " + k + ":\t" + "SEQUENCE WITH MOTIF: " + str(len(value)))
        print("____________________________________")
        print()

        #Iterates through each sequence motif was found in, using position of occurence to search .gb file
        #Also uses positional information to find commonalities between all positions with motif

        for item in value:

            itemAsList = item.split("\t")
            sequenceName = itemAsList[0]
            finalList = []
            gBDuplicate = 0
            inverseSeqRecordPos = set()

            for item in range(len(finalListRaw)):
                if item % 2 != 0:
                    continue
                else:
                    finalList.append(finalListRaw[item])

            if sequenceName in genBankDict:

                sequenceRecord = genBankDict[sequenceName][1]
                gBDuplicate = sequenceRecord
                seqRecordPosAsList = []
                matchesIngenBank = 0
                #If .gb is found, searches .gb file for area where position occurs
                for record in range(len(sequenceRecord.features)):

                    locationRangeString = str(sequenceRecord.features[record].location)[1:-4]
                    rangeAsList = locationRangeString.split(":")

                    if re.search('^(oin)', str(rangeAsList[0])):

                        continue

                    elif int(rangeAsList[0]) == 0:

                        continue
                    #Makes list of all start/end positions in .gb file
                    if len(sequenceRecord.features[record].qualifiers) > 2 and int(rangeAsList[0]) != 0:
                        seqRecordPosAsList.append(str(rangeAsList[0] + "," + rangeAsList[1]))

                #Uses list of .gb gene start/end positions and uses it to look at the inverse (areas without genes)
                for item in range(len(seqRecordPosAsList)-1):

                    indexPos = seqRecordPosAsList[item].split(',')
                    indexPosNextIter = seqRecordPosAsList[item + 1].split(',')

                    firstPosOfGap = int(indexPos[1]) + 1
                    secondPosOfGap = int(indexPosNextIter[0]) - 1

                    inverseSeqRecordPos.add(str(firstPosOfGap) + ',' + str(secondPosOfGap))


                #test, ONLY needed for phylo analysis.

                #testSection this section is only for phylo analysis:


                #####

                # for r in inverseSeqRecordPos:
                #     rangeSplit = r.split(',')
                # print(sequenceName)
                listInOrder = sorted(inverseSeqRecordPos, key=lambda i : int(i.split(',')[0]))
                allNCDict[sequenceName] = listInOrder

        # for k,v in allNCDict.items():
        #     print("SEQUENCE NAME: " + str(k))
        #     print()
        #     print(v[0:10])
        #     print()

        allNCDictOneSubstring = {}
        allNCDictWithoutNCOccurences = {}

        for k,v in allNCDict.items():
            sequenceJustNC = ''
            currentSequence = str(headerAndSequenceAllEntries[k])
            # print(len(currentSequence))
            for i in v:
                rangeSplit = i.split(',')
                posOne = int(rangeSplit[0])
                posTwo = int(rangeSplit[1])
                if len(range(posOne,posTwo)) > 0:
                    sequenceJustNC += currentSequence[posOne: posTwo]

            allNCDictOneSubstring[k] = sequenceJustNC


        for k,v in allNCDict.items():
            sequenceJustNC = ''
            currentSequence = str(headerAndSequenceAllEntries[k])
            # print(len(currentSequence))
            for i in v:
                rangeSplit = i.split(',')
                posOne = int(rangeSplit[0])
                posTwo = int(rangeSplit[1])
                if len(range(posOne,posTwo)) > 0:
                    currentNCRegion = currentSequence[posOne: posTwo]
                    ifMotifFound = 0

                    for motif in regexList:
                        if re.search(motif, currentNCRegion):
                            ifMotifFound = 1

                    if ifMotifFound == 0:
                        sequenceJustNC += currentNCRegion

            allNCDictWithoutNCOccurences[k] = sequenceJustNC

        # for k,v in allNCDict.items():
        #     # print("Full")
        #     # print(len(allNCDictOneSubstring[k]))
        #     # print("W/O")
        #     # print(len(allNCDictWithoutNCOccurences[k]))




        listForNCCalc = []

        for k,v in allNCDictOneSubstring.items():
            # print(k + ":" + str(len(v)))
            currentSequence = str(headerAndSequenceAllEntries[k])
            percNCInGenome = int(len(v))/int(len(currentSequence))
            listForNCCalc.append(percNCInGenome)

        if len(listForNCCalc) > 0:
            print("AVERAGE PORTION OF EACH SEQUENCE IN NON-CODING REGION")
            print()
            print("MEAN: " + str(round(statistics.fmean(listForNCCalc),3)*100) + "%")
            print("MEDIAN: " + str(round(statistics.median(listForNCCalc),3)*100))
            print("STDEV: " + str(round(statistics.pstdev(listForNCCalc),2)*100))
            print("MIN: " + str(round(min(listForNCCalc),3)*100))
            print("MAX: " + str(round(max(listForNCCalc),3)*100))
            print()


sys.stdout = original_stdout

print(outputDirectory)

#This section handles bar graph generation

for k, v in sorted(overLapDictionary.items()):

    histogramMotifPos[k] = []
    q = set(v)
    value = sorted(list(q))
    for item in value:
        itemAsList = item.split("\t")
        #print(itemAsList[0])
        positionInBrackets = itemAsList[2].split(',')
        #print(positionInBrackets)
        positionOfString = 0
        for num in positionInBrackets:
            valueOnlyDigits =  re.sub("[^0-9]", "", num)

            if positionOfString%2 == 0:
                #print(valueOnlyDigits)
                histogramMotifPos[k].append(valueOnlyDigits)


for k,v in histogramMotifPos.items():

    x = []

    for value in v:
        x.append(int(value))


    weights = np.ones(len(x)) / len(x)
    plt.hist(x, 10,
             histtype='bar',
             color='green',  lw=3, stacked=True)  # , weights=weights)
    # plt.hist(x)

    plt.xlabel(xAxisLabel)
    plt.ylabel(yAxisLabel)
    plt.title(str(k) + '\n', fontsize =12)
    plt.savefig(outputDirectory + "/" + str(k) + '_histogram.jpg')
    plt.clf()