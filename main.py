from Bio import SeqIO
import os
import typing
import glob
import math

class VSM:

    def getFasta(self):
        sequence = []
        listOfFastaFiles = []
        os.chdir("..\\FASTA")
        currentDirectory = os.getcwd()
        # Description: currentDirectory = F:\University\University_Term8\Course BioInformation\AssignMent\AssignMent2\FASTA
        # with glob you can specify a type of file to list like this
        # Method 1
        # listOfFastaFiles = [f for f in glob.glob("*.fasta")]
        # Method 2
        for file in glob.glob("*.fasta"):
            listOfFastaFiles.append(file)

        for eachFastaFile in listOfFastaFiles:

            gene = currentDirectory + "\\" + eachFastaFile
            records = SeqIO.parse(gene, 'fasta')
            for record in records:
                sequence.append(record.seq)

        return sequence
        # print(listOfFastaFiles)
        # gene = currentDirectory + "\\PPARD.fasta"
        # records = SeqIO.parse(gene, 'fasta')

        for record in records:
            sequence.append(record.seq)

        return sequence

    def calFrequency(self, sequenceList):

        listOfFrequency = []

        A = {
            "f1": 0,
            "f2": 0,
            "f3": 0,
            "p1": 0.0,
            "p2": 0.0,
            "p3": 0.0
        }
        C = {
            "f1": 0,
            "f2": 0,
            "f3": 0,
            "p1": 0.0,
            "p2": 0.0,
            "p3": 0.0
        }
        G = {
            "f1": 0,
            "f2": 0,
            "f3": 0,
            "p1": 0.0,
            "p2": 0.0,
            "p3": 0.0
        }
        T = {
            "f1": 0,
            "f2": 0,
            "f3": 0,
            "p1": 0.0,
            "p2": 0.0,
            "p3": 0.0
        }


        for index, eachSequence in enumerate(sequenceList):
            for eachCharacter in eachSequence:
                if eachCharacter == 'A':
                    if index == 0:
                        A['f1'] += 1
                    elif index == 1:
                        A['f2'] += 1
                    else:
                        A['f3'] += 1
                elif eachCharacter == 'C':
                    if index == 0:
                        C['f1'] += 1
                    elif index == 1:
                        C['f2'] += 1
                    else:
                        C['f3'] += 1
                elif eachCharacter == "G":
                    if index == 0:
                        G['f1'] += 1
                    elif index == 1:
                        G['f2'] += 1
                    else:
                        G['f3'] += 1
                elif eachCharacter == "T":
                    if index == 0:
                        T['f1'] += 1
                    elif index == 1:
                        T['f2'] += 1
                    else:
                        T['f3'] += 1
                else:
                    print("[-] Error: \"Not a Valid Character\"")
                    break

        listOfFrequency.append(A)
        listOfFrequency.append(C)
        listOfFrequency.append(G)
        listOfFrequency.append(T)
        return listOfFrequency

    def calTotalFrequency(self, listOfFrequency):
        A = {}
        C = {}
        G = {}
        T = {}
        total = {
            "f1":0,
            "f2":0,
            "f3":0,
            "p1":0.0,
            "p2":0.0,
            "p3":0.0,
        }


        for index, value in enumerate(listOfFrequency):
            if index == 0:
                A = value
            elif index == 1:
                C = value
            elif index == 2:
                G = value
            else:
                T = value

        total["f1"] = A["f1"] + C["f1"] + G["f1"] + T["f1"]
        total["f2"] = A["f2"] + C["f2"] + G["f2"] + T["f2"]
        total["f3"] = A["f3"] + C["f3"] + G["f3"] + T["f3"]
        return total

    def calProbability(self, listOfFrequency, totalFrequency):

        A = {}
        C = {}
        G = {}
        T = {}

        for index, value in enumerate(listOfFrequency):
            if index == 0:
                A = value
            elif index == 1:
                C = value
            elif index == 2:
                G = value
            else:
                T = value

        A["p1"] = (A["f1"] + 1)/ (totalFrequency["f1"] + 4)
        C["p1"] = (C["f1"] + 1)/ (totalFrequency["f1"] + 4)
        G["p1"] = (G["f1"] + 1)/ (totalFrequency["f1"] + 4)
        T["p1"] = (T["f1"] + 1)/ (totalFrequency["f1"] + 4)

        A["p2"] = (A["f2"] + 1)/ (totalFrequency["f2"] + 4)
        C["p2"] = (C["f2"] + 1)/ (totalFrequency["f2"] + 4)
        G["p2"] = (G["f2"] + 1)/ (totalFrequency["f2"] + 4)
        T["p2"] = (T["f2"] + 1)/ (totalFrequency["f2"] + 4)

        A["p3"] = (A["f3"] + 1)/ (totalFrequency["f3"] + 4)
        C["p3"] = (C["f3"] + 1)/ (totalFrequency["f3"] + 4)
        G["p3"] = (G["f3"] + 1)/ (totalFrequency["f3"] + 4)
        T["p3"] = (T["f3"] + 1)/ (totalFrequency["f3"] + 4)

        listOfFrequency = []
        listOfFrequency.append(A)
        listOfFrequency.append(C)
        listOfFrequency.append(T)
        listOfFrequency.append(G)
        return listOfFrequency

    def calVectorSize(self, listOfFrequencyWithProbability):
        A = {}
        C = {}
        G = {}
        T = {}

        vectorSize = []

        for index, value in enumerate(listOfFrequencyWithProbability):
            if index == 0:
                A = value
            elif index == 1:
                C = value
            elif index == 2:
                G = value
            else:
                T = value

        p1 = math.sqrt(math.pow(A["p1"], 2) + math.pow(C["p1"], 2) + math.pow(G["p1"], 2) + math.pow(T["p1"], 2))
        p2 = math.sqrt(math.pow(A["p2"], 2) + math.pow(C["p2"], 2) + math.pow(G["p2"], 2) + math.pow(T["p2"], 2))
        p3 = math.sqrt(math.pow(A["p3"], 2) + math.pow(C["p3"], 2) + math.pow(G["p3"], 2) + math.pow(T["p3"], 2))

        print("|P1|=" + str(p1))
        print("|P2|=" + str(p2))
        print("|P3|=" + str(p3))

        vectorSize.append(p1)
        vectorSize.append(p2)
        vectorSize.append(p3)


        return vectorSize

    def cosineSimilarity(self, listOfFrequencyWithProbability, vectorSize):
        A = {}
        C = {}
        G = {}
        T = {}

        for index, value in enumerate(listOfFrequencyWithProbability):

            similarityP1P2 = 0.0
            similarityP1P3 = 0.0
            similarityP2P3 = 0.0

            similarityList = []


            if index == 0:
                A = value
            elif index == 1:
                C = value
            elif index == 2:
                G = value
            else:
                T = value


        similarityP1P2 = (A["p1"] * A["p2"] + C["p1"] * C["p2"] + G["p1"]*G["p2"] + T["p1"]*T["p2"]) / (vectorSize[0] * vectorSize[1])
        similarityList.append( "cos(p1, p2)=" + str(similarityP1P2))
        similarityP1P3 = (A["p1"] * A["p3"] + C["p1"] * C["p3"] + G["p1"]*G["p3"] + T["p1"]*T["p3"]) / (vectorSize[0] * vectorSize[2])
        similarityList.append("cos(p1, p3)=" + str(similarityP1P3))
        similarityP2P3 = (A["p2"] * A["p3"] + C["p2"] * C["p3"] + G["p2"]*G["p3"] + T["p2"]*T["p3"]) / (vectorSize[1] * vectorSize[2])
        similarityList.append("cos(p2, p3)=" + str(similarityP2P3))


        return similarityList

if __name__ == '__main__':
    # object instantiation
    vsm = VSM()
    # get sequences from .FASTA files
    sequence = vsm.getFasta()
    # list of frequencies for each Nucleotide sequence
    listOfFrequncy = vsm.calFrequency(sequence)
    # list of total frequency for each Nucleotide sequence
    totalFreqeuncy = vsm.calTotalFrequency(listOfFrequncy)
    # list of total frequencies with probability for each Nucleotide sequence
    listOfFrequencyWithProbability = vsm.calProbability(listOfFrequncy, totalFreqeuncy)
    # vectorSize example (|P1|=0.5008889011201644 # |P2|=0.5062315613026167 # |P3|=0.5025505211612074)
    vectorSize = vsm.calVectorSize(listOfFrequencyWithProbability)
    # calculate Cosine Similarity
    result = vsm.cosineSimilarity(listOfFrequencyWithProbability, vectorSize)
    [print(val) for val in result]