# Oragnism object
# It allocates the full data structure
import numpy as np
import random
class OrganismObject():

    # Organism constructor
    def __init__(self, ID, conf, rootNode = None):
        self.ID = ID
        self.rootNode = rootNode
        self.bestFitnessScore = None
        self.CUMULATIVE_FIT_METHOD = conf["CUMULATIVE_FIT_METHOD"]
    
    # Setters an getters
    def setRootNode(self, rootNode):
        self.rootNode = rootNode
    
    def getID(self):
        return self.ID
    
    def setID(self, iD):
        self.ID = iD

    # Mutates a part of the organism
    def mutate(self, orgFactory):
        #print("Mutating organism {}".format(self.ID))
        nNodes = self.countNodes()
        randomNode = random.randint(0, nNodes - 1)
        mutatedNode = self.getNode(randomNode)
        mutatedNode.mutate(orgFactory)
    
    # Return the fitness of the organism for a given DNA sequence
    def getBestAllFitness(self, sDNA):
        
        sequenceLength = len(sDNA)

        # get all PSSM recognizers
        aPssmObjects = self.rootNode.getAllPssm([])

        #check position where it fits (position)
        pssmPositionScoreTable = []

        # Go through all PSSM objects to check where it fits
        for pssm in aPssmObjects:
            pssmLength = pssm.length
            maxScore = float("-inf")    
            position = 0

            # Check every position possible on the sequence
            for pos in range(sequenceLength - pssmLength):
                score = pssm.getScore(sDNA[pos:pos+pssmLength])

                # Update max score if the actual score is acctually better
                # Also check that the position is not overlapping other pssm object
                if score > maxScore and not self.isOverlapping(pos, pssmLength,
                        pssmPositionScoreTable):
                    
                    maxScore = score
                    position = pos

            # Add to a table the ID, maxScore and position of a pssm object
            pssmPositionScoreTable.append((pssm.ID, maxScore, position, pssmLength))

        # This is to see how its going... not definitive
        ''' 
        verbose = random.random()<0.001
        if verbose:
            print("\n{}".format(sDNA))
            mapPositions = " "*len(sDNA)
            
            for ids, sc, pos, length in pssmPositionScoreTable:
                strId = str(ids)
                while len(strId) < length:
                    strId += "*"
                mapPositions = mapPositions[0:pos] + strId + mapPositions[pos+length:]
                
            print(mapPositions + "\n")
        '''
        # call recursively to get the total fitness of the organism
        finalScore, distance = self.rootNode.getBestAll(pssmPositionScoreTable)           
        
        # return score in that dataset
        return finalScore
    
    # Given a table of positioned pssms, check that new pssm is not overlapping
    def isOverlapping(self, positionN, lengthN, table):
        
        for idP, scoreP, positionP, lengthP in table:

            # On the left side
            if positionN < positionP:
                if positionN + lengthN > positionP:
                    return True

            # On the right side
            elif positionN > positionP:
                if positionP + lengthP > positionN:
                    return True

            # Same position
            else:
                return True




        return False

    # Return the total Fitness for an array of DNA sequences and the fitness method (for now only one) 
    def getScore(self, aDNA):
        
        score = 0
        # sum method returns the sum of all fitness to DNA sequences
        if self.CUMULATIVE_FIT_METHOD == "sum":

            for sDNA in aDNA:
                score+=self.getBestAllFitness(sDNA.lower())

        # mean method returns the mean of all fitness to SNA sequences
        elif self.CUMULATIVE_FIT_METHOD == "mean":

            scores = []
            for sDNA in aDNA:
                scores.append(self.getBestAllFitness(sDNA.lower()))
            score = np.mean(scores)

        return score

    # Returns a node Number N based on in-order search. [0-N)
    def getNode(self, objective):
        nodeCount = 0
        return self.rootNode.getNode(objective, nodeCount)

    # Set positions a given a node and ID where is has to be
    def setNode(self, node, ID):
        
        # CANNOT BE BY ID BC SECOND IS OVERWRITED

        print("node.ID = {} ID to change {}".format(node.ID, ID))
        if self.rootNode.ID == ID:
            self.rootNode = node
        else:
            self.rootNode.setNode(node, ID)

    # Get the parent node of a given ID and if it is the left child
    def getParent(self, ID):
        
        if self.rootNode.ID == ID:
            return [True]

        return self.rootNode.getParent(ID)

        '''
        if node1.ID == ID:
            # Return itself and specify child is on left side
            return [self, True]
        elif node2.ID == ID:
            # Return itself and specify child is on right side
            return [self, False]
        checkNode = None
        checkNode = node1.getParent(ID)
        if checkNode == None:
            checkNode = node2.getParent(ID)
        return checkNode
        '''

    # Returns the number of nodes of the organism
    def countNodes(self):
        return self.rootNode.countNodes()


    # Reset IDs of the full organism
    def resetIDs(self):
        firstID = 0
        self.rootNode.resetID(firstID)

    # Prints the whole tree data structure
    def print(self):
        print()
        print("***** Organism {} *****".format(self.ID))
        self.rootNode.print(1)

    # Exports the whole tree data structure
    def export(self, filename):
        organismFile = open(filename, "w+")
        organismFile.write("***** Organism {} *****".format(self.ID))

        self.rootNode.export(organismFile, 0)
        organismFile.write("\n")
        organismFile.close()

    # Exports all DNA sequences organism binding to a file
    def exportResults(self, aDNA, filename):

        #Sort the array, so its always shown in the same order
        aDNA.sort()

        resultsFile = open(filename, "w+")

        for sDNA in aDNA:
            
            sDNA = sDNA.lower()
            sequenceLength = len(sDNA)

            # get all PSSM recognizers
            aPssmObjects = self.rootNode.getAllPssm([])

            #check position where it fits (position)
            pssmPositionScoreTable = []

            # Go through all PSSM objects to check where it fits
            for pssm in aPssmObjects:
                pssmLength = pssm.length
                maxScore = float("-inf")    
                position = 0

                # Check every position possible on the sequence
                for pos in range(sequenceLength - pssmLength):
                    score = pssm.getScore(sDNA[pos:pos+pssmLength])

                    # Update max score if the actual score is acctually better
                    # Also check that the position is not overlapping other pssm object
                    if score > maxScore and not self.isOverlapping(pos, pssmLength,
                        pssmPositionScoreTable):
                    
                        maxScore = score
                        position = pos

                # Add to a table the ID, maxScore and position of a pssm object
                pssmPositionScoreTable.append((pssm.ID, maxScore, position, pssmLength))
            
            resultsFile.write("\n{}\n".format(sDNA))
            mapPositions = " "*len(sDNA)
            
            for ids, sc, pos, length in pssmPositionScoreTable:
                strId = str(ids)
                while len(strId) < length:
                    strId += "*"
                mapPositions = mapPositions[0:pos] + strId + mapPositions[pos+length:]
                
            resultsFile.write(mapPositions + "\n")

        resultsFile.close()
                


