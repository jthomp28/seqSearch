import os
import statistics
from Bio import SeqIO
from Bio import Entrez
from Bio.Align.Applications import ClustalOmegaCommandline
import re
import regex
import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt



#Used to create object to store information about each sequence
class inputSeqCapsule:

    #Used on object initiation, creates multiple variables from list input)
    def __init__(self, list):

        self.sequence = list[0][0]

        # print(self.sequence)

        self.GCcontent = list[0][1]

        # print(self.GCcontent)

        self.gapSegements = list[0][2]

        # print(self.gapSegements)

        self.sequenceLength = list[0][3]

        # print(self.sequenceLength)

        self.countA = list[0][4]

        # print(self.countA)

        self.countT = list[0][5]

        # print(self.countT)

        self.countC = list[0][6]

        # print(self.countT)

        self.countG = list[0][7]

        # print(self.countG)

        self.name = list[0][8]

        self.subsequenceDict = {}

        self.genBankFile = ''

    #Prints output to screen
    def objectInformation(self, l="N", ):

        a = "Sequence Name:\t" + self.name
        b = "Sequence Length:\t" + str(self.sequenceLength)
        c = "GC%:\t" + str(self.GCcontent)
        d = "Gaps:\t" + str(self.gapSegements)
        e = "Adenosine:\t" + str(self.countA)
        f = "Thymine:\t" + str(self.countT)
        g = "Cytosine:\t" + str(self.countC)
        h = "Guanine:\t" + str(self.countG)

        if l == "N":

            return a + "\n" + b + "\n" + c + "\n"

        elif l == "Y":

            return a + "\n" + b + "\n" + c + "\n" + d + "\n" + e + "\n" + f + "\n" + g + "\n" + h + "\n"

    #Main function called by each object, creates a dictionary of substrings that occur more than once in each sequence
    #Records positions of matches using (x,y) notation, stored as set() in each dictionary value
    def seqSearch(self, subsequence="", window=0):


        #If a nucletotide is repeated more than 10 (default value) times, it is stored here and output to ""
        #Planned update, currently seqSearch handles motifs >3 nt in length
        repeatNtCheck = []

        #Checks if there is input, if "" then using sliding window to match substring
        if subsequence == "":

            windowStart = 0
            windowStop = windowStart + window

            #If sequence length < window returns sequence
            if (len(self.sequence) <= windowStop):

                return (self.sequence)



            else:

                #Uses sliding window, subtracts end one position to account starting position
                for i in range(len(self.sequence) - (windowStop + 1)):

                    s1 = self.sequence

                    #Sets search pattern
                    pattern = (self.sequence[i:i + windowStop])

                    if pattern in motifsToIgnore:

                        continue

                    else:

                        #Index is used here to denote the exact substring position
                        subStringStart = ((self.sequence).index(self.sequence[i:i + windowStop]))

                        subStringEnd = ((self.sequence).index(self.sequence[i:i + windowStop]) + window)

                        #%s is string formatting, this allows use of 'r' in other regex function
                        r = regex.compile('(%s)' % pattern, re.IGNORECASE)


                        #Finds match and then
                        for match in regex.finditer(r, s1):

                            s = match.start()

                            e = match.end()

                            subSeqKey = s1[s:e]

                            matchAsRange = "(" + str(s) + "," + str(e) + ")"

                            originalStringAsRange = "(" + str(subStringStart) + "," + str(subStringEnd) + ")"

                            #Does not include itself as a match
                            if s == subStringStart and e == subStringEnd:

                                next

                            #If subSeqKey is not in dictionary, it creates an entry
                            elif subSeqKey not in self.subsequenceDict.keys():

                                self.subsequenceDict[subSeqKey] = set()

                                self.subsequenceDict[subSeqKey].add(originalStringAsRange)

                            #If subSeqKey is in dictionary it prints output to screen of match
                            else:

                                print('Motif match "%s" at %d:%d' % (s1[s:e], s, e))
                                print('Original motif is at ' + str(subStringStart) + ": " + str(subStringEnd))
                                self.subsequenceDict[subSeqKey].add(matchAsRange)

#Opens FASTA file and stores them in headerAndSequence
#edit: great spot for seq. length options
def fastaFormatFileInput(multiFASTA):

    currentHeader = ''

    #Opens FASTA file and stores them in headerAndSequence
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
    #edit: could alter what region of sequence to search here

    for key, value in headerAndSequence.items():

        if endPositionofSurvey > len(str(value)):
            print("End position is greater than sequence length. Please choose a different value")
            break

        elif startPositionofSurvey > 0 and endPositionofSurvey == 0:
            headerAndSequence[key] = [value[startPositionofSurvey:]]


        elif startPositionofSurvey and endPositionofSurvey > 0:

            headerAndSequence[key] = [value[startPositionofSurvey:endPositionofSurvey]]

        else:
            headerAndSequence[key] = [value]

        #used to ensure .gb file is read correctly
        headerAndSequenceFullLength[key] = [value]

    # print(headerAndSequence)
    outputComparison['ENTRY COUNT'] = len(headerAndSequence.keys())
    a_file.close()

#Appends information about genome to previous entry in headerAndSequence using updatedValue
#updatedValue = [value + sequenceInformation + [(key)]]
#headerAndSequence[key] = updatedValue
#updatedValue is also appended to inputSeqMasterList
#inputSeqMasterList is used to initialize
def sequenceInformation(dictObj):
    # Displays preliminary sequence information

    seqInfoDict = dictObj
    entryNumber = 1

    for key, value in seqInfoDict.items():

        GCcount = 0
        gapSegments = 0
        countA, countT, countC, countG = (0,) * 4
        #edit: use count function here instead of list
        for i in str(value):

            if i == 'A':
                countA += 1
            elif i == 'T':
                countT += 1
            elif i == 'C':
                countC += 1
            elif i == 'G':
                countG += 1

            elif i == '-':

                gapSegments += 1

            if i == 'G' or i == 'C':
                GCcount += 1

        GCpercentage = GCcount / len(str(value)[2:-2])

        sequenceLength = len(str(value)[2:-2])



        sequenceInformation = [GCpercentage, gapSegments, sequenceLength, countA, countT, countC, countG]

        fullLengthOfSequence = str(headerAndSequenceFullLength[key])



        updatedValue = [value + sequenceInformation + [(key)]]



        headerAndSequence[key] = updatedValue


        # print(str(entryNumber) + ')')
        # print(key)
        # print('GC Percetage:\t' + str(GCpercentage))
        # print('Gaps:\t' + str(gapSegments))
        # print('Sequence Length:\t' + str(sequenceLength))
        # print('Adenosine:\t' + str(countA))
        # print('Thymine:\t' + str(countT))
        # print('Cytosine:\t' + str(countC))
        # print('Guanine:\t' + str(countG))
        # print('\n')

        entryNumber += 1
        inputSeqMasterList.append(updatedValue)
        # print(inputSeqMasterList)

    # print(inputSeqMasterList)
    # print(inputSeqMasterList)
    # print(headerAndSequence)

#Iterates through genBank file previously associated with seqObj, attempts to find position in range of entries
def genBankLocationObj(seqObj, arrayOne, windowSize, notFoundList, substring, pos=0):

    matchesIngenBank = 0

    sequenceRecord = seqObj

    #Iterates through all GenBank features
    for record in range(len(sequenceRecord.features)):

        locationRangeString = str(sequenceRecord.features[record].location)[1:-4]

        rangeAsList = locationRangeString.split(":")

        currentProduct = (sequenceRecord.features[record].qualifiers.get('product', ''))

        positionForLabel = " (" + str(pos) + "-" + str(pos + windowSize) + ") "

        labelForNoEntry = str(sequenceRecord.description) + positionForLabel + substring


        # print(rangeAsList[0])

        if re.search('^(oin)', str(rangeAsList[0])):

            continue

        elif int(rangeAsList[0]) == 0:

            continue


        elif (pos) in range(int(rangeAsList[0]), int(rangeAsList[1])) and len(
                sequenceRecord.features[record].qualifiers) > 2:

            print("START POSITION OF MOTIF:" + str(pos) + "\n")
            print(sequenceRecord.features[record])

            arrayOne.append(currentProduct)

            matchesIngenBank = 1


    #If not matches are found, this is the output
    if matchesIngenBank == 0:

        print(str(pos) + " WAS NOT FOUND IN THE RANGE OF ANY RECORDED FEATURES \n\n")
        #Takes note of positions not found in recorded features
        notFoundList.append(labelForNoEntry)

#this function links .gb files with input FASTA sequences using the sequence itself as identifier
def genBankFileLink(dict, path):

    fastaDict = dict
    keyAndAccessionDict = {}

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
        #print(file)
        with open(os.path.join(path, file)) as gbFile:
            fullGBText = SeqIO.read(gbFile, "genbank")

        for k, v in fastaDict.items():

            seqObjLabel = str(k)

            gbSeq = ''.join(str(fullGBText.seq))

            fastaName = fastaDict[seqObjLabel].name

            if fastaName in headerAndSequenceFullLength.keys():

                fastaSeq = headerAndSequenceFullLength[fastaName]

                fastaSeq = ''.join(fastaSeq)

            # print(str(len(fastaSeq)))
            # print(str(len(gbSeq)))

            if str(fastaSeq) == str(gbSeq):
                fastaDict[seqObjLabel].genBankFile = fullGBText

#Calculates the % of the genome not contained in gb Records
#Gives perspective if a substring occurs in an area more often than expected by chance
def genomePercentNotCoveredInGB(seqObj):

    overLapAdjustment = 0
    valueForOverlapAdjustment = 0

    if type(seqObj) == str:

        return

    if seqObj is None:

        return

    else:
        #seqObj is gbFile saved to object
        sequenceRecord = seqObj
        dictionaryCheck = {}
        lengthOfSequenceNoSpace = 0
        itemCheck = set()
        itemtest = []

        #Iterates through each record in gbFile
        for record in range(len(sequenceRecord.features)):

            locationRangeString = str(sequenceRecord.features[record].location)[1:-4]
            #print(locationRangeString)
            rangeAsList = locationRangeString.split(":")
            #print(rangeAsList)


            #print((rangeAsList[1]))
            #print(len(sequenceRecord.seq))

            if re.search('^(oin)', str(rangeAsList[0])):

                continue

            #Adds the length of the record to lengthOfSequenceNoSpace
            else:

                    rangeInputLength = len(range(int(rangeAsList[0]), int(rangeAsList[1])))
                    #print(rangeInputLength)
                    rangeInput = range(int(rangeAsList[0]), int(rangeAsList[1]))
                    #print(rangeInput)




                    if int(rangeAsList[0]) == 0:

                        continue

                    elif rangeInput in itemCheck:

                        continue

                    else:

                        if valueForOverlapAdjustment > int(rangeAsList[0]):
                            adjustment = valueForOverlapAdjustment - int(rangeAsList[0])
                            overLapAdjustment = overLapAdjustment + adjustment


                        #print(rangeAsList)
                        valueForOverlapAdjustment = int(rangeAsList[1])
                        #print(valueForOverlapAdjustment)
                        itemCheck.add(rangeInput)
                        #print(rangeInputLength)
                        #print(rangeInput)
                        lengthOfSequenceNoSpace = lengthOfSequenceNoSpace + rangeInputLength

        # print('TESTING')
        # print(lengthOfSequenceNoSpace)
        # print(overLapAdjustment)
        # print(len(sequenceRecord.seq))
        # print(lengthOfSequenceNoSpace/len(sequenceRecord.seq))
        # print((lengthOfSequenceNoSpace - overLapAdjustment) / len(sequenceRecord.seq))

        #print("\n")

        #divides lengthOfSequenceNoSpace by full sequence length to determine difference in length
        #this sequence length is taken directly from the gb file
        return (lengthOfSequenceNoSpace - overLapAdjustment) / len(sequenceRecord.seq)

#Used to create masterSequenceDictionary, this dictionary assigns each header and sequence its own object
#SO,S1,S2,S3..... is the key
#The value will be the class object, this is neccessary for substring comparisons between the sequence objects
#This object contains sequence information, a seqSearch function and variables needed for downstream analysis
#This dictionary is one of the values reset after each run (A1, A2, A3 etc.)
def classCreation(list):
    objectNameDict = {}

    allSequenceInformation = list

    # print((len(allSequenceInformation)))
    # print((len(headerAndSequence.values())))

    for x in range(len(headerAndSequence.values())):
        objectName = "S" + str(x)

        objectNameDict[objectName] = ""

        # print(objectNameDict)

    numberLabel = 0

    #inputSeqCapsule is class creation function. Takes list of all sequences and gives each its own object. numberLabel
    #is used to ensure all sequences get an object of their own
    for k, v in objectNameDict.items():
        objectNameDict[k] = inputSeqCapsule(allSequenceInformation[numberLabel])

        numberLabel += 1

        # print(objectNameDict[k].name)

        # print(objectNameDict[k].subsequenceDict)

        print(objectNameDict[k].objectInformation("Y"))

    return (objectNameDict)

#Using masterSequenceDictionary as, the seqSearch function within each object is called
#Heavily commented
#Each seqSearch input is saved in comparisonDict, which is used for file output
#Substring = key, position of occurrence in that ONE sequence = value (subseqDictionary)
#Substring = key, position and specific sequence = value (comparisonDict)
#edit: optional parameter for sequence length
def seqSearch(dInput, searchWindow, summaryFileName, maxSequenceLength=0, filePrefix='', numberOfAlignments=0,
              outputDirectory=''):
    #Contains masterSequenceDictionary, S#: classObject containing information about 1 specific sequence
    #Each sequence from FASTA input file gets its own object
    totalSeqDict = dInput

    substringPositionsForRange = []

    productAsList = []

    translationAsDict = {}

    #used to set ranges for later alignments, observers alignments in a 10000 bp range
    binSize = 0

    originalBinSize = 0

    seqWithoutGenBank = []

    #stores paths for alignments using ClustalO
    fastaPathsForAlignment = []

    substringAlignmentList = []

    listOfNoGBEntry = []

    gbGenevSpace = []

    #edit: Used to determine string length, better method to be added
    inputMessage = 0

    masterTextList = []

    averageForSegmentOccurances = []

    #Used to set the range for alignments to be made. ie/ 60000 bp with 6 alignments is 6 10000 bp bins
    if maxSequenceLength > 0 and numberOfAlignments > 0:
        binSize = int(maxSequenceLength / numberOfAlignments)
        originalBinSize = binSize

        # print(binSize)

    #Calls seqSearch in each object of masterSequenceDictionary
    for k, v in totalSeqDict.items():
        print(totalSeqDict[k].name + "\n")

        totalSeqDict[k].seqSearch(window=searchWindow)

        #subSequenceDict is the dictionary within each object for seqSearch output. Substring = key, position and sequence = value
        print(totalSeqDict[k].subsequenceDict.items())
        print(len(totalSeqDict[k].subsequenceDict.keys()))

        print("END OF ENTRY")
        print("_________________\n\n")

    #Dictionary that compares sequences (object by object). Makes note
    comparisonDict = {}
    uniqueSubSeqs = set()

    #Handles case where only one sequence is in FASTA
    #Compares masterSequenceDictionary to itself, ignoring exact comparisons. Matching substrings are passed to comparisonDict
    #comparisonDict[sharedKey] = []
    #The list value stores the sequence and positions in the sequence where the match occured
    #Thus, a dictionary with substrings as the key will be created. The values will be all the sequences and respective positions
    #where matches to that substring occured
    if len(totalSeqDict.items()) == 1:

        # print("LENGTH == 1")

        for k, v in totalSeqDict.items():

            # print(k)
            # print(v)

            for j, w in totalSeqDict.items():
                # point of intersection between sequeunces, if there is a match that substrng is used as a key and the two sequences
                # as well as where the match occured will be stored in the list contained in the value

                subSeqIntersection = list(
                    totalSeqDict[k].subsequenceDict.keys() & totalSeqDict[j].subsequenceDict.keys())
                # print(subSeqIntersection)
                # print(totalSeqDict[k].name)
                # print(totalSeqDict[j].name)


                if len(subSeqIntersection) > 0:

                    for x in range(len(subSeqIntersection)):

                        sharedKey = subSeqIntersection[x]
                        # print(sharedKey)
                        # x will be key in both dictionaries, can use to access them

                        outerLoopMatch = [totalSeqDict[k].name, totalSeqDict[k].subsequenceDict[sharedKey], k]

                        innerLoopMatch = [totalSeqDict[j].name, totalSeqDict[j].subsequenceDict[sharedKey], j]

                        if subSeqIntersection[x] not in comparisonDict.keys():

                            comparisonDict[sharedKey] = []

                            comparisonDict[sharedKey].append(outerLoopMatch)

                            comparisonDict[sharedKey].append(innerLoopMatch)

                        elif subSeqIntersection[x] in comparisonDict.keys():

                            comparisonDict[sharedKey].append(outerLoopMatch)

                            comparisonDict[sharedKey].append(innerLoopMatch)


        # print(comparisonDict)
        comparisonOfAllInput.append(comparisonDict)

        original_stdout = sys.stdout

        with open(outputDirectory + "/" + summaryFileName, 'w') as f:

            sys.stdout = f

            entryNumber = 1

            comparisonDictOutput = {}

            #removes duplicate entries from dictionary
            for key, value in comparisonDict.items():
                NoDupOutput = []

                for item in comparisonDict[key]:
                    if item not in NoDupOutput:
                        NoDupOutput.append(item)
                comparisonDictOutput[key] = NoDupOutput

            comparisonDictMotifEdit = {}

            #edit: Add filter variable to start
            #Uses regexList to filter only potential motifs, without this all entries are kept
            for x in regexList:

                for k,v in comparisonDict.items():

                    if re.search(x,k):

                        comparisonDictMotifEdit[k] = v

                    else:

                        continue

            #print(comparisonDictOutput)
            # print(len(comparisonDictOutput))

            #edit: motif Edit
            #Saves to original variable name
            comparisonDictOutput = comparisonDictMotifEdit

            #print(comparisonDictOutput)
            # print(len(comparisonDictOutput))

            #edit: Change this to reflect all motif occurence at top of page
            print("TOTAL SEQUENCES IN FILE: " + str(len(totalSeqDict)))
            print("MOTIF AND NUMBER OF SEQUENCES LOCATED IN\n")

            #prints output to file using comparisonDict and string manipulations.
            #Lists all motifs shared

            totalMotifs = 0
            for y in sorted(comparisonDictOutput, key=lambda y: len(list(comparisonDictOutput[y])), reverse=True):

                print(str(y) + ":" + str(len(list(comparisonDictOutput[y]))))
                substringAlignmentList.append(str(y))
                totalMotifs += len(list(comparisonDictOutput[y]))

            print("\n")
            print('TOTAL MOTIFS: ' + str(totalMotifs))
            print("\n")
            print("#" * 40)
            print("#" * 40)

            totalSubstringMatches = 0

            #section dedicated to each motif
            for k in sorted(comparisonDictOutput, key=lambda k: len(comparisonDictOutput[k]), reverse=True):

                similarKeys = comparisonDict[k]
                # print(len(similarKeys))

                numberOfMatches = 0

                # print("LOOK:"+ str(similarKeys))

                similiarKeysNoDupicates = []

                for item in similarKeys:
                    if item not in similiarKeysNoDupicates:
                        similiarKeysNoDupicates.append(item)

                similarKeys = similiarKeysNoDupicates

                # print("LOOK:" + str(similarKeys) + "\n")

                print("_" * 40)
                print("_" * 40)
                print("_" * 40)

                print("ENTRY " + "(" + str(entryNumber) + ")")

                print("MOTIF: " + str(k))

                print("SEQUENCES WITH MOTIF: " + str(len(similarKeys)))

                print("_" * 40)

                #This is the updated input, similar keys are the entries for the current motif iteration
                for x in sorted(similarKeys, key=len, reverse=True):

                    name = x[0]

                    #uses a function to sort for simplicity
                    positionInSequence = sorted(list(x[1]), key = motifPositionSort)

                    # print(positionInSequence)

                    #All individuals matches are calculated addition this each iteration
                    numberOfSubstringMatches = len(positionInSequence)

                    averageForSegmentOccurances.append(numberOfSubstringMatches)

                    # print(str(positionInSequence))


                    pattern = re.compile(r"[0-9]+")

                    #adds only numbers to list
                    posAsList = pattern.findall(str(positionInSequence))

                    seqObjLabel = x[2]

                    gbFile = totalSeqDict[seqObjLabel].genBankFile

                    sequenceLength = len(totalSeqDict[seqObjLabel].sequence)

                    # print(str(gbFile))

                    currentName = totalSeqDict[seqObjLabel].name

                    print("_" * 40)

                    print("SEQUENCE NAME: " + str(name))

                    print("POSITION(S) OF OCCURRENCE IN ORIGINAL SEQUENCE: " + str((positionInSequence)))

                    print("TOTAL OCCURRENCES IN SEQUENCE: " + str(numberOfSubstringMatches))

                    totalSubstringMatches = totalSubstringMatches + numberOfSubstringMatches

                    print("_" * 40)
                    print("_" * 40)
                    # print(str(posAsList) + "\n")

                    for value in range(len(posAsList)):
                        if type(gbFile) == '':

                            # print(currentName + " NO GENBANK FILE PROVIDED")
                            seqWithoutGenBank.append(currentName)

                            continue

                        elif type(gbFile) == str:

                            continue

                        #List input is [0,1,2,3,4], this ensures only the start position is used to find sequence
                        elif value % 2 == 0:

                            currentPos = posAsList[value]
                            totalLength = len(gbFile.seq)
                            adjustedPosition = startPositionofSurvey + int(currentPos)
                            print(adjustedPosition)

                            if startPositionofSurvey == 0 and endPositionofSurvey == 0:

                                print("POSITION OF MOTIF: " + currentPos + "\n")
                                genBankLocationObj(gbFile, productAsList, searchWindow,
                                                   listOfNoGBEntry, str(k),
                                                   int(currentPos))

                                substringPositionsForRange.append(currentPos)

                            elif startPositionofSurvey > 0:

                                print("POSITION OF MOTIF: " + str(adjustedPosition) + "\n")
                                # print(sequenceLength)
                                print(totalLength)
                                print(startPositionofSurvey)
                                print(adjustedPosition)

                                genBankLocationObj(gbFile, productAsList, searchWindow,
                                                   listOfNoGBEntry, str(k),
                                                   int(adjustedPosition))

                                substringPositionsForRange.append(adjustedPosition)

                            # print(returnTest)

                entryNumber += 1

                print("END OF ENTRY")
                print("#" * 40)
                print("#" * 40)
                print("\n")

            print("TOTAL OCCURRENCES OF ALL MOTIFS: " + str(totalSubstringMatches))
            masterTextList.append("TOTAL OCCURRENCES OF ALL MOTIFS: " + str(totalSubstringMatches))


            #if entry not found in range of gb it is appended to listOfNoGBEntry
            print(str(len(listOfNoGBEntry)) + " POSITIONS WERE NOT FOUND IN RANGE OF GB RECORDS")
            #print(comparisonDictOutput)
            print("\n")

            if totalSubstringMatches > 0:
                print(str(100*round(len(listOfNoGBEntry) / totalSubstringMatches,
                                3)) + "% OF DETECTED MOTIFS WERE FOUND OUT OF RANGE OF GB RECORDS\n")
                masterTextList.append(str(100*round(len(listOfNoGBEntry) / totalSubstringMatches,
                                3)) + "% OF DETECTED MOTIFS WERE FOUND OUT OF RANGE OF GB RECORDS\n")


            #edit: Put this higher in notes, change function if needed
            substringAlignmentListEdit = set(substringAlignmentList)
            # print(substringAlignmentListEdit)
            print("TOTAL MOTIF COUNT: " + str(len(substringAlignmentListEdit)) + "\n")
            masterTextList.append("TOTAL MOTIF COUNT: " + str(len(substringAlignmentListEdit)) + "\n")


            for k, v in totalSeqDict.items():

                # print(k)
                # print(v)

                currentGenBank = totalSeqDict[k].genBankFile

                valueForGBCalc = genomePercentNotCoveredInGB(currentGenBank)

                if valueForGBCalc is None or valueForGBCalc < 0:

                    continue

                else:
                    gbGenevSpace.append(valueForGBCalc)

            print(gbGenevSpace)
            if len(gbGenevSpace) > 0:
                print("MEAN: " + str(statistics.fmean(gbGenevSpace)))
                print("MEDIAN: " + str(statistics.median(gbGenevSpace)))
                print("STDEV: " + str(statistics.pstdev(gbGenevSpace)))
                print("MIN: " + str(min(gbGenevSpace)))
                print("MAX: " + str(max(gbGenevSpace)))
                noRecord = 1 - statistics.fmean(gbGenevSpace)

                print("ON AVERAGE " + str(100*round(1 - statistics.fmean(gbGenevSpace),
                                                3)) + "% OF EACH SEQUENCE WAS FOUND OUT OF RANGE OF LISTED GB RECORDS")
                masterTextList.append("ON AVERAGE " + str(100*round(1 - statistics.fmean(gbGenevSpace),
                                                3)) + "% OF EACH SEQUENCE WAS FOUND OUT OF RANGE OF LISTED GB RECORDS")

                print()

            fileOutput = open(outputDirectory + "/" + str(fileName.split('.')[0]) + "_SubstringAlignmentFASTA.txt", 'w')

            fileInputAsString = outputDirectory + "/" + str(fileName.split('.')[0]) + "_SubstringAlignmentFASTA.txt"

            fileOutputAsString = outputDirectory + "/" + str(fileName.split('.')[0]) + "_SubstringAlignmentALN.txt"

            #Aligns all substring using ClustalO for similarities between them

            sequencesInFile = 0

            if len(substringAlignmentList) > 1:
                for item in range(len(substringAlignmentList)):

                    sequencesInFile += 1

                    substringEntry = str(substringAlignmentList[item])

                    fastaLabel = "MOTIF: " + str(substringAlignmentList[item])

                    fileOutput.write(">" + fastaLabel + "\n" + substringEntry + "\n")




                fileOutput.close()


            clustalo_exe = pathToClustalEXE
            # print(fileInputAsString)
            if len(substringAlignmentList) > 1:
                clustalomega_cline = ClustalOmegaCommandline(clustalo_exe, infile=fileInputAsString,
                                                         outfile=fileOutputAsString, force=True)

                clustalomega_cline()

            if len(seqWithoutGenBank) > 0:
                print("THE FOLLOWING ENTRIES HAVE NO GENBANK FILE PROVIDED")
                for entry in seqWithoutGenBank:
                    print(entry)

            print("\n")

            if len(substringPositionsForRange) > 0:
                rangeLadder(substringPositionsForRange, 60000, 6)
                matplotHistogramList[fileName] = substringPositionsForRange

            print("\n")

            productOutputDict = {}

            for entry in productAsList:

                # print(str(entry))

                if str(entry) in productOutputDict:

                    productOutputDict[str(entry)] += 1

                else:

                    productOutputDict[str(entry)] = 1

            # print(productOutputDict)

            for k, v in sorted(productOutputDict.items(), key=lambda item: item[1], reverse=True):
                print(str(k) + ": " + str(v))


            print('\n')
            if len(averageForSegmentOccurances) > 0:
                masterTextList.append(["STATS OF OCCURANCE GIVEN EACH SITE", "MEAN: " + str(statistics.fmean(averageForSegmentOccurances)),
                                       "MEDIAN: " + str(statistics.median(averageForSegmentOccurances)), "STDEV: " + str(statistics.pstdev(averageForSegmentOccurances)),
                                       "MIN: " + str(min(averageForSegmentOccurances)),"MAX: " + str(max(averageForSegmentOccurances))])
                outputComparison[fileName] = masterTextList
                outputComparison[fileName].append(outputComparison['ENTRY COUNT'])
                del outputComparison['ENTRY COUNT']


            # for item in masterTextList:
            #
            #     print(str(item) + "\n")


            # print("STATS OF OCCURANCE GIVEN EACH SITE")
            # print("MEAN: " + str(statistics.fmean(averageForSegmentOccurances)))
            # print("MEDIAN: " + str(statistics.median(averageForSegmentOccurances)))
            # print("STDEV: " + str(statistics.pstdev(averageForSegmentOccurances)))
            # print("MIN: " + str(min(averageForSegmentOccurances)))
            # print("MAX: " + str(max(averageForSegmentOccurances)))

            sys.stdout = original_stdout


        f.close()


    else:

        for k, v in totalSeqDict.items():

            # print(k)
            # print(v)

            for j, w in totalSeqDict.items():
                #point of intersection between sequeunces, if there is a match that substrng is used as a key and the two sequences
                #as well as where the match occured will be stored in the list contained in the value
                subSeqIntersection = list(totalSeqDict[k].subsequenceDict.keys() & totalSeqDict[j].subsequenceDict.keys())
                # print(subSeqIntersection)
                # print(totalSeqDict[k].name)
                # print(totalSeqDict[j].name)

                if totalSeqDict[k].subsequenceDict.keys() == totalSeqDict[j].subsequenceDict.keys() and totalSeqDict[k].name == totalSeqDict[j].name:

                    next


                elif len(subSeqIntersection) > 0:

                    for x in range(len(subSeqIntersection)):

                        sharedKey = subSeqIntersection[x]
                        # print(sharedKey)
                        # x will be key in both dictionaries, can use to access them

                        #subsequenceDict match is the position of substring matches
                        outerLoopMatch = [totalSeqDict[k].name, totalSeqDict[k].subsequenceDict[sharedKey], k]

                        innerLoopMatch = [totalSeqDict[j].name, totalSeqDict[j].subsequenceDict[sharedKey], j]

                        if subSeqIntersection[x] not in comparisonDict.keys():

                            comparisonDict[sharedKey] = []

                            comparisonDict[sharedKey].append(outerLoopMatch)

                            comparisonDict[sharedKey].append(innerLoopMatch)

                        elif subSeqIntersection[x] in comparisonDict.keys():

                            comparisonDict[sharedKey].append(outerLoopMatch)


                            comparisonDict[sharedKey].append(innerLoopMatch)

        #matches gb files with respective sequence object sharing a name
        # customGenBankAssociation('C:/Users/Jason Thompson/Desktop/seqSearch/PhagesDB_accession_numbers.txt', 2)

        # print(comparisonDict)
        #Used to compare all entries as a whole (A1, A2, A3) not individual sequences
        comparisonOfAllInput.append(comparisonDict)

        original_stdout = sys.stdout


        with open(outputDirectory + "/" + str(summaryFileName), 'w') as f:

            sys.stdout = f

            entryNumber = 1

            comparisonDictOutput = {}

            for key, value in comparisonDict.items():
                NoDupOutput = []

                for item in comparisonDict[key]:
                    if item not in NoDupOutput:
                        NoDupOutput.append(item)
                comparisonDictOutput[key] = NoDupOutput

            comparisonDictMotifEdit = {}

            for x in regexList:

                for k, v in comparisonDictOutput.items():

                    if re.search(x, k):

                        comparisonDictMotifEdit[k] = v

                    else:

                        continue


            # print(len(comparisonDictOutput))
            comparisonDictOutput = comparisonDictMotifEdit
            #print(comparisonDictOutput)
            # print(len(comparisonDictOutput))


            print("MOTIF AND NUMBER OF SEQUENCES LOCATED IN\n")
            print("TOTAL SEQUENCES IN FILE: " + str(len(totalSeqDict)))

            totalMotifs = 0
            for y in sorted(comparisonDictOutput, key=lambda y: len((comparisonDictOutput[y])), reverse=True):

                print(str(y) + ":" + str(len((comparisonDictOutput[y]))))
                totalMotifs += len(list(comparisonDictOutput[y]))

                substringAlignmentList.append(str(y))

            print("\n")
            print('TOTAL MOTIFS: ' + str(totalMotifs))
            print("\n")
            print("#" * 40)
            print("#" * 40)

            totalSubstringMatches = 0

            for k in sorted(comparisonDictOutput, key=lambda k: len(comparisonDictOutput[k]), reverse=True):

                similarKeys = comparisonDict[k]
                #print(len(similarKeys))

                numberOfMatches = 0

                # print("LOOK:"+ str(similarKeys))

                similiarKeysNoDupicates = []

                for item in similarKeys:
                    if item not in similiarKeysNoDupicates:
                        similiarKeysNoDupicates.append(item)

                similarKeys = similiarKeysNoDupicates

                # print("LOOK:" + str(similarKeys) + "\n")

                print("_" * 40)
                print("_" * 40)
                print("_" * 40)

                print("ENTRY " + "(" + str(entryNumber) + ")")

                print("MOTIF: " + str(k))

                print("SEQUENCES WITH MOTIF: " + str(len(similarKeys)))

                print("_" * 40)

                for x in sorted(similarKeys, key=len, reverse=True):

                    name = x[0]

                    #uses a function to sort for simplicity
                    positionInSequence = sorted(list(x[1]), key = motifPositionSort)

                    numberOfSubstringMatches = len(positionInSequence)
                    averageForSegmentOccurances.append(numberOfSubstringMatches)

                    # print(str(positionInSequence))

                    pattern = re.compile(r"[0-9]+")

                    posAsList = pattern.findall(str(positionInSequence))


                    seqObjLabel = x[2]

                    gbFile = totalSeqDict[seqObjLabel].genBankFile

                    sequenceLength = len(totalSeqDict[seqObjLabel].sequence)

                    # print(str(gbFile))

                    currentName = totalSeqDict[seqObjLabel].name

                    print("_" * 40)

                    print("SEQUENCE NAME: " + str(name))

                    print("POSITION(S) OF OCCURRENCE IN SEQUENCE: " + str(positionInSequence))

                    print("TOTAL OCCURRNCES IN SEQUENCE: " + str(numberOfSubstringMatches))

                    totalSubstringMatches = totalSubstringMatches + numberOfSubstringMatches

                    print("_" * 40)
                    print("_" * 40)
                    # print(str(posAsList) + "\n")

                    for value in range(len(posAsList)):


                        if type(gbFile) == '':

                            # print(currentName + " NO GENBANK FILE PROVIDED")
                            seqWithoutGenBank.append(currentName)

                            continue

                        elif type(gbFile) == str:

                            continue

                        elif value % 2 == 0:


                            currentPos = posAsList[value]
                            differenceInLength = len(gbFile.seq) - sequenceLength
                            adjustedPosition = startPositionofSurvey + int(currentPos)

                            if startPositionofSurvey == 0 and endPositionofSurvey == 0:

                                print("POSITION OF SUBSTRING: " + currentPos + "\n")
                                genBankLocationObj(gbFile, productAsList, searchWindow, listOfNoGBEntry, str(k),
                                                   int(currentPos))

                                substringPositionsForRange.append(currentPos)

                            elif startPositionofSurvey > 0:

                                print("POSITION OF SUBSTRING: " + str(adjustedPosition) + "\n")

                                # print(sequenceLength)
                                # print(differenceInLength)
                                # print(adjustedPosition)

                                genBankLocationObj(gbFile, productAsList, searchWindow, listOfNoGBEntry, str(k),
                                                   int(adjustedPosition))

                                substringPositionsForRange.append(adjustedPosition)

                            # print(returnTest)

                entryNumber += 1

                print("END OF ENTRY")
                print("#" * 40)
                print("#" * 40)
                print("\n")

            # print(len(substringPositionsForRange))
            # print(substringPositionsForRange)
            # print(productAsList)
            # print(translationAsDict)

            print("TOTAL OCCURRENCES OF ALL MOTIFS: " + str(totalSubstringMatches))
            masterTextList.append("TOTAL OCCURENCES OF ALL MOTIFS: " + str(totalSubstringMatches))




            print(str(len(listOfNoGBEntry)) + " POSITIONS WERE NOT FOUND IN RANGE OF GB RECORDS")
            # print(comparisonDictOutput)
            if totalSubstringMatches > 0:
                print(str(100*round(len(listOfNoGBEntry)/totalSubstringMatches, 3)) + "% OF DETECTED MOTIFS WERE FOUND OUT OF RANGE OF GB RECORDS\n")
                masterTextList.append(str(round(len(listOfNoGBEntry) / totalSubstringMatches,
                                                3)) + "% OF DETECTED MOTIFS WERE FOUND OUT OF RANGE OF GB RECORDS\n")

            listForSubstringPosComp = []


            print()

            substringAlignmentListEdit = set(substringAlignmentList)

            print("TOTAL MOTIFS COUNT: " + str(len(substringAlignmentListEdit)) + "\n")
            masterTextList.append("TOTAL MOTIFS COUNT: " + str(len(substringAlignmentListEdit)) + "\n")


            print()
            print()

            for k,v in totalSeqDict.items():

                #print(k)
                #print(v)

                currentGenBank = totalSeqDict[k].genBankFile

                valueForGBCalc = genomePercentNotCoveredInGB(currentGenBank)

                if valueForGBCalc is None:

                    continue

                else:
                    gbGenevSpace.append(valueForGBCalc)

            if len(gbGenevSpace) > 0:
                print("THE FOLLOWING VALUES DESCRIBE WHAT PERCENTAGE OF THE INPUT SEQUENCE(S) CONSIST OF GENE CONTENT")
                print("MEAN: " + str(statistics.fmean(gbGenevSpace)))
                print("MEDIAN: " + str(statistics.median(gbGenevSpace)))
                print("STDEV: " + str(statistics.pstdev(gbGenevSpace)))
                print("MIN: " + str(min(gbGenevSpace)))
                print("MAX: " + str(max(gbGenevSpace)))
                print()
                print("ON AVERAGE " + str(100*round(1 - statistics.fmean(gbGenevSpace), 3)) + "% OF EACH SEQUENCE WAS FOUND OUT OF RANGE OF LISTED GB RECORDS")
                masterTextList.append("ON AVERAGE " + str(round(1 - statistics.fmean(gbGenevSpace), 3)) + "% OF EACH SEQUENCE WAS FOUND OUT OF RANGE OF LISTED GB RECORDS")

                print()


            #this section handles substring alignment
            fileOutput = open(outputDirectory + "/" + str(fileName.split('.')[0]) + "_SubstringAlignmentFASTA.txt", 'w')

            fileInputAsString = (outputDirectory + "/" + str(fileName.split('.')[0]) + "_SubstringAlignmentFASTA.txt")

            fileOutputAsString = (outputDirectory + "/" + str(fileName.split('.')[0]) + "_SubstringAlignmentALN.txt")


            sequencesInFile = 0
            if len(substringAlignmentList) > 1:
                for item in range(len(substringAlignmentList)):

                    sequencesInFile += 1
                    substringEntry = str(substringAlignmentList[item])

                    fastaLabel = "MOTIF: " + str(substringAlignmentList[item])

                    fileOutput.write(">" + fastaLabel + "\n" + substringEntry + "\n")

                fileOutput.close()

            clustalo_exe = pathToClustalEXE
            # print(os.path.getsize(fileInputAsString))
            if len(substringAlignmentList) > 1:
                clustalomega_cline = ClustalOmegaCommandline(clustalo_exe, infile=fileInputAsString,
                                                         outfile=fileOutputAsString, force=True)
                clustalomega_cline()


            if len(seqWithoutGenBank) > 0:
                print("THE FOLLOWING ENTRIES HAVE NO GENBANK FILE PROVIDED")
                for entry in seqWithoutGenBank:
                    print(entry)




            print("\n")
            #edit: No parameters in function
            print("MOTIFS WERE FOUND IN THE FOLLOWING BIN RANGES")
            rangeLadder(substringPositionsForRange, 60000, 6)
            matplotHistogramList[fileName] = substringPositionsForRange

            print("\n")



            productOutputDict = {}

            for entry in productAsList:

                # print(str(entry))

                if str(entry) in productOutputDict:

                    productOutputDict[str(entry)] += 1

                else:

                    productOutputDict[str(entry)] = 1

            # print(productOutputDict)

            #last entry of output file
            for k, v in sorted(productOutputDict.items(), key = lambda item:item[1], reverse = True):
                print(str(k) + ": " + str(v))

            if len(averageForSegmentOccurances) > 0:
                masterTextList.append(
                    ["STATS OF OCCURANCE GIVEN EACH SITE", "MEAN: " + str(statistics.fmean(averageForSegmentOccurances)),
                     "MEDIAN: " + str(statistics.median(averageForSegmentOccurances)),
                     "STDEV: " + str(statistics.pstdev(averageForSegmentOccurances)),
                     "MIN: " + str(min(averageForSegmentOccurances)), "MAX: " + str(max(averageForSegmentOccurances))])
                outputComparison[fileName] = masterTextList
                outputComparison[fileName].append(outputComparison['ENTRY COUNT'])
                del outputComparison['ENTRY COUNT']

                # for item in masterTextList:
                #     print(str(item))



                # print("MEAN: " + str(statistics.fmean(averageForSegmentOccurances)))
                # print("MEDIAN: " + str(statistics.median(averageForSegmentOccurances)))
                # print("STDEV: " + str(statistics.pstdev(averageForSegmentOccurances)))
                # print("MIN: " + str(min(averageForSegmentOccurances)))
                # print("MAX: " + str(max(averageForSegmentOccurances)))

            sys.stdout = original_stdout




        f.close()

#Used to make note of areas where substrings appear, bins are decided by user.
#ie/ 60000 bp/ 6 bins = 10000 bp a bin or occurrences per 10000bp
#edit: potential revision for cleaner function
def rangeLadder(numList, length, bins):
    valuesToCheck = numList
    sequenceLength = length
    numberOfBins = bins
    rangeOutputDictionary = {}

    z = 0

    for entry in range(numberOfBins):
        keyName = str(z)

        rangeOutputDictionary[keyName] = []

        z = z + 1

    # print(rangeOutputDictionary)

    for value in valuesToCheck:

        # print(value)

        sizeOfBin = int(sequenceLength / numberOfBins)
        valueToAdd = int(sequenceLength / numberOfBins)
        startString = 0

        # print(str(sizeOfBin), str(valueToAdd))

        for x in range(numberOfBins):

            # print(str(x))
            # print(str(startString) + ":" + str(sizeOfBin))

            if int(value) in range(startString, sizeOfBin):
                # print("YES IT IS")

                rangeOutputDictionary[str(x)].append(value)

            startString = startString + valueToAdd

            sizeOfBin = sizeOfBin + valueToAdd

    sizeOfBin = int(sequenceLength / numberOfBins)
    valueToAdd = int(sequenceLength / numberOfBins)
    startString = 0

    for j, k in rangeOutputDictionary.items():
        print(str(startString) + " - " + str(sizeOfBin) + ":\t" + str(len(rangeOutputDictionary[str(int(j))])))

        startString = startString + valueToAdd

        sizeOfBin = sizeOfBin + valueToAdd

    # print(rangeOutputDictionary)

#Extra function
def variableWipe(listOne, listTwo, listThree, DictOne, DictTwo):
    genBankFileNames = []

    inputSeqMasterList = []

    headerAndSequence = {}

    headerIdentifiers = ['>', '.', '+']

    masterSequenceDictionary = {}

    print(genBankFileNames, inputSeqMasterList, headerAndSequence, headerIdentifiers, masterSequenceDictionary)

#compares masterSequenceDictionary with regexList to determine how motifs are shared between groups
#masterSeuqneceDict key = S# value = classObject
#Only uses keys of each masterSequenceDict
#Those with matching motif to any in regexList will be reported
def substringCompareCoverage(list, dictInput):

    listInput = list
    sequenceDict = dictInput

    for x in listInput:

        dictForOutput = {}

        for k, v in sequenceDict.items():

            name = sequenceDict[k].name

            subseqDict = sequenceDict[k].subsequenceDict.keys()



            for item in subseqDict:

                if re.search(x, item):

                    if name in dictForOutput.keys():

                        dictForOutput[name].append(item)

                    else:

                        dictForOutput[name] = [item]

        print("SUBSTRING: " + x)
        print("LENGTH OF ALL ENTRIES: " + str(len(sequenceDict)))
        print("LENGTH OF ENTRIES CONTAINING POT. MOTIF: " + str(len(dictForOutput)))
        print("_______________________________________")

        for k, v in dictForOutput.items():
            print(str(k) + "\t" + str(v) + "\t" + str(len(v)))

        print("END OF ENTRY")
        print()

#This function is used to sort entries of the output file in the 3' to 5' direction (0 - Sequence End)
#Passed directly to sorted() in seqSearch function to ensure positions are in ascending order
def motifPositionSort(ele):
    positionInput = ele
    positionInput = positionInput.split(',')
    startingPoint = positionInput[0]
    return int(startingPoint[1:])

#Variables required for functions, do not edit.

headerAndSequence = {}

headerIdentifiers = ['>', '.', '+']

genBankFileNames = []

inputSeqMasterList = []

inputSeqMasterListFull = []

comparisonOfAllInput = []

outputComparison = {}

matplotHistogramList = {}

headerAndSequenceFullLength = {}

#used for debugging
currentCluster = 'A12'

###### OPTIONAL USER INPUT VARIABLES###################
#######################################################

#This is an optional variable, this allows users to filter sequences using any motifs of interest
#Without this input, all motifs matching window size will be captured
regexList = ['TTACGA[AG]TCA', 'GGTG[ATCG][ATCG]TGTCAA', 'TTCTCTGTCAA', 'GTGCGATGTCAA', 'GTCCGTTGTCA', 'TTCGGTGTCAA',
             'ACGTAGTCAA', 'GATGAGTGTCAA', 'AGTGGGTGTCAA', 'TTGGGTGTCAA', 'GGGAT[AT]GTCAA', 'GTTGGGTGTCAA',
             '[ATCG][ATCG]TG[ATCG][ATCG]TGTCAA', 'TTC[ATCG][ATCG]TGTCAA', 'TTGTGTGTCTA',
             'GGGG[GT]ATGTCAA', 'GTTCGATGTCAA', 'TACGGTGTCCAA', 'GTAC[AG]GGGTCAA']


#If a motif is abundant and users do not want it inclued they can place in in this list and seqSearch will skip.
#This is helpful for repeats "ie AAAA", but by default no motifs are ommited
#Motifs here should match  chosen window size

motifsToIgnore = []

#If required, users can choose to only explore indicated regions of the inputer sequences.
#If nothing is entered the full sequence of each entry will be surveyed

startPositionofSurvey = 0
endPositionofSurvey = 0

###### REQUIRED USER INPUT VARIABLES###################
#######################################################


# #enter a name for summary file
# fileName = currentCluster + '_output.txt'

#enter a prefix for accessory file output
filePrefix = currentCluster

#enter a directory to output file
outputDirectory = 'C:/Users/Jason Thompson/Desktop/seqSearch/finalOutputs/autoOutput/filter/FullSequence/debugOutput'

# path to FASTA file(s), function will detect whether it is a file or directory.
#If this is a directory, please ensure only files pertaining to analysis are present
#This will be changed in future seqSearch updates.
pathToFasta = 'C:/Users/Jason Thompson/Desktop/seqSearch/debugFiles/A1OutputFile.txt'

#path to GenBank file directory
# Please ensure only files pertaining to analysis are present
# This will be changed in future seqSearch updates.
pathToGenBankDirectory = 'C:/Users/Jason Thompson/Desktop/seqSearch/finalOutputs/scrapCheck/allGBFiles'

#input the window size (k) here
windowSizeInput = 16

#input max sequence length of dataset here
maxSeqLength = 60000

#input desired bin size here
binSize = 1000

#path to clustalOmega.exe
pathToClustalEXE = 'C:/Users/Jason Thompson/Desktop/seqSearch/clustal-omega-1.2.2-win64/clustalo.exe'

#Variables For Bar Graph

xAxisLabel = 'Sequence position (bp)'

yAxisLabel = 'Number Of Motifs Found At Sequence Position On X-Axis'

#############################################
#############################################

if os.path.exists(pathToFasta):
    if os.path.isfile(pathToFasta):

        fileName = os.path.basename(pathToFasta)

        #Goes through each file, putting header/dict into headerAndSequence
        fastaFormatFileInput(pathToFasta)

        #adds GC% and other sequence info to headerAndSequence, edits dictionary not a new instance
        #also appends sequence to inputSeqMasterList for following function
        sequenceInformation(headerAndSequence)

        #using classCreation function, inputSeqMasterList containing all sequence information is used
        #to give each sequence its own object with its information stored in local variable
        masterSequenceDictionary = classCreation(inputSeqMasterList)


        # Adds .gb file to each object if it is the input directory (parameter 3)
        #if statement handles scenario where only a portion of each sequence is to be surveyed

        genBankFileLink(masterSequenceDictionary,
                        pathToGenBankDirectory)

        #the seqSearch function in each object is called, searching sequences and storing motif in a local dict.
        #variable stored in each object
        seqSearch(masterSequenceDictionary, windowSizeInput, fileName, maxSeqLength, filePrefix, binSize, outputDirectory)

        # resets variables
        genBankFileNames = []

        inputSeqMasterList = []

        headerAndSequence = {}

        headerIdentifiers = ['>', '.', '+']

        masterSequenceDictionary = {}

        outputComparison = {}


    elif os.path.isdir(pathToFasta):
        for file in os.listdir(pathToFasta):
            fileName = file
            filePrefix = 'one'
            print(file)
            pathToFile = pathToFasta + file
            titleOfGraph = file

            # Goes through each file, putting header/dict into headerAndSequence
            fastaFormatFileInput(pathToFile)

            # adds GC% and other sequence info to headerAndSequence, edits dictionary not a new instance
            # also appends sequence to inputSeqMasterList for following function
            sequenceInformation(headerAndSequence)

            # using classCreation function, inputSeqMasterList containing all sequence information is used
            # to give each sequence its own object with its information stored in local variable
            masterSequenceDictionary = classCreation(inputSeqMasterList)

            # Adds .gb file to each object if it is the input directory (parameter 3)
            # if statement handles scenario where only a portion of each sequence is to be surveyed

            genBankFileLink(masterSequenceDictionary,
                            pathToGenBankDirectory)

            # the seqSearch function in each object is called, searching sequences and storing motif in a local dict.
            # variable stored in each object
            seqSearch(masterSequenceDictionary, windowSizeInput, fileName, maxSeqLength, filePrefix, binSize,
                      outputDirectory)

            # resets variables
            genBankFileNames = []

            inputSeqMasterList = []

            headerAndSequence = {}

            headerIdentifiers = ['>', '.', '+']

            masterSequenceDictionary = {}

            outputComparison = {}


for k,v in matplotHistogramList.items():

    x = []

    for value in v:
        x.append(int(value))

    weights = np.ones(len(x)) / len(x)
    plt.hist(x, 10,
             histtype='bar',
             color='green', lw=3, stacked=True)
    # plt.hist(x)

    titleNoSuffix = str(k.split('.')[0])
    plt.xlabel(xAxisLabel)
    plt.ylabel(yAxisLabel)
    plt.title(str(k) + '\n', fontsize =10)
    plt.savefig(outputDirectory + "/" + titleNoSuffix + '.jpg')
    plt.clf()


print(outputDirectory)
print('\n')

