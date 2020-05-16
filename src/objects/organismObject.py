# Oragnism object
# It allocates the full data structure
import numpy as np
import random


class OrganismObject:

    # Organism constructor
    def __init__(self, ID, conf, rootNode=None):
        self.ID = ID
        self.rootNode = rootNode
        self.meanFitnessScore = 0  # Mean fitness score on all datasets computed
        self.numSequencesScored = 0  # Number of datasets computed
        self.numNodes = 0  # Number of nodes of the oranism
        self.CUMULATIVE_FIT_METHOD = conf["CUMULATIVE_FIT_METHOD"]
        self.MUTATE_PROBABILITY_SUBSTITUTE_PSSM = conf[
            "MUTATE_PROBABILITY_SUBSTITUTE_PSSM"
        ]
        self.MUTATE_PROBABILITY_RISE_CHILD = conf["MUTATE_PROBABILITY_RISE_CHILD"]
        self.MUTATE_PROBABILITY_SUNK_CHILD = conf["MUTATE_PROBABILITY_SUNK_CHILD"]
        self.MUTATE_PROBABILITY_NODE_MUTATION = conf["MUTATE_PROBABILITY_NODE_MUTATION"]
        self.isTracked = False

    # Setters an getters
    def setRootNode(self, rootNode):
        self.rootNode = rootNode
        # Compute nodes on the organism
        self.numNodes = 0

    def getID(self):
        return self.ID

    def setID(self, iD):

        self.ID = iD

    def setIsTracked(self, newTracked):
        self.isTracked = newTracked

    # Mutates a part of the organism
    def mutate(self, orgFactory):

        # Substitute a random node by a random PSSM
        if random.random() < self.MUTATE_PROBABILITY_SUBSTITUTE_PSSM:

            nNodes = self.countNodes()
            randomNode = random.randint(0, nNodes - 1)
            substitutedNode = self.getNode(randomNode)
            parentNode = self.getParent(substitutedNode.ID)
            # It can be improved so you can not pass the length by parameter jsjs
            newNode = orgFactory.createPSSM(orgFactory.PWM_LENGTH)

            if parentNode["isRootNode"]:

                self.rootNode = newNode
            else:

                if parentNode["isLeftSide"]:

                    parentNode["self"].setNode1(newNode)
                else:

                    parentNode["self"].setNode2(newNode)

            self.resetIDs()

        # Rise the level of a child
        if random.random() < self.MUTATE_PROBABILITY_RISE_CHILD:

            nNodes = self.countNodes()
            randomNode = random.randint(0, nNodes - 1)
            risedNode = self.getNode(randomNode)
            parentNode = self.getParent(risedNode.ID)

            #
            if not risedNode.ID == self.rootNode.ID:

                parent1 = self.getParent(risedNode.ID)
                parent2 = self.getParent(parent1["self"].ID)

                if parent2["isRootNode"]:
                    self.rootNode = risedNode
                else:
                    if parent2["isLeftSide"]:
                        parent2["self"].setNode1(risedNode)
                    else:
                        parent2["self"].setNode2(risedNode)
                self.resetIDs()

        # Add complexity to the organism
        if random.random() < self.MUTATE_PROBABILITY_SUNK_CHILD:
            nNodes = self.countNodes()
            randomNode = random.randint(0, nNodes - 1)
            sukenNode = self.getNode(randomNode)
            parentNode = self.getParent(sukenNode.ID)

            newNode = orgFactory.createConnection(0)

            if parentNode["isRootNode"]:
                self.rootNode = newNode
            else:
                if parentNode["isLeftSide"]:

                    parentNode["self"].setNode1(newNode)
                else:
                    parentNode["self"].setNode1(newNode)

            if random.random() < 0.5:
                newNode.setNode1(sukenNode)
            else:
                newNode.setNode2(sukenNode)
            self.resetIDs()

        # Mutate a random node
        if random.random() < self.MUTATE_PROBABILITY_NODE_MUTATION:

            nNodes = self.countNodes()
            randomNode = random.randint(0, nNodes - 1)
            mutatedNode = self.getNode(randomNode)
            mutatedNode.mutate(orgFactory)

    # Returns the complexity of the organism
    def getComplexity(self, meanNodes, meanFitness):
        # Complexity is calculed as:
        # meanFitnessScore * # nodes / meanNodes

        return meanFitness * self.numNodes / meanNodes

    # Return the fitness of the organism for a given DNA sequence
    def getBestAllFitness(self, sDNA):

        sequenceLength = len(sDNA)

        # get all PSSM recognizers
        aPssmObjects = self.rootNode.getAllPssm([])

        # check position where it fits (position)
        pssmPositionScoreTable = []

        # Go through all PSSM objects to check where it fits
        for pssm in aPssmObjects:
            pssmLength = pssm.length
            maxScore = float("-inf")
            position = 0

            # Possible max scores of the secuences are in aPositions
            aPositions = []
            # Check every position possible on the sequence
            for pos in range(sequenceLength - pssmLength):
                score = pssm.getScore(sDNA[pos : pos + pssmLength])

                # Update max score if the actual score is acctually better
                # Also check that the position is not overlapping other pssm object
                if not self.isOverlapping(pos, pssmLength, pssmPositionScoreTable):

                    if score > maxScore:  # First ocurrence of the max score
                        maxScore = score
                        aPositions = [pos]

                    elif (
                        score == maxScore
                    ):  # Another ocurrence of the max score and possible binding site
                        aPositions.append(pos)

            # Here we checked the whole sequence, but now we have to select the
            # best pssm based on pssm memory
            position = min(aPositions, key=lambda x: abs(x - pssm.positionMemory))

            # And tell th pssm to remember it for the future
            pssm.rememberPosition(pos)
            # Add to a table the ID, maxScore and position of a pssm object
            pssmPositionScoreTable.append((pssm.ID, maxScore, position, pssmLength))

        # call recursively to get the total fitness of the organism
        finalScore, distance = self.rootNode.getBestAll(pssmPositionScoreTable)

        # return score in that dataset
        return finalScore

    # Given a table of positioned pssms, check that new pssm is not overlapping
    def isOverlapping(self, positionN, lengthN, table):
        # Compare the new position with already positioned pssms
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

    # Return the total Fitness for an array of DNA sequences and the fitness method
    # It also updates mean sequences score
    def getScore(self, aDNA):

        score = 0
        aDNALength = len(aDNA)

        # Set the pssm memory to the middle of the secuence
        self.rootNode.setPositionMemory(len(aDNA[0]) / 2.0)

        # sum method returns the sum of all fitness to DNA sequences
        if self.CUMULATIVE_FIT_METHOD == "sum":

            for sDNA in aDNA:
                score += self.getBestAllFitness(sDNA.lower())

            # Assign values if its the first time
            if self.numSequencesScored == 0:

                self.numSequencesScored = aDNALength
                self.meanFitnessScore = score / aDNALength
            # Calculate new values else
            else:
                self.meanFitnessScore = (
                    self.meanFitnessScore * self.numSequencesScored + score
                ) / (aDNALength + self.numSequencesScored)
                self.numSequencesScored += aDNALength

        # mean method returns the mean of all fitness to SNA sequences
        elif self.CUMULATIVE_FIT_METHOD == "mean":

            scores = []
            for sDNA in aDNA:
                scores.append(self.getBestAllFitness(sDNA.lower()))
            score = np.mean(scores)

            # Assign values if its the first time
            if self.numSequencesScored == 0:

                self.numSequencesScored = aDNALength
                self.meanFitnessScore = score
            # Calculate new values else
            else:
                self.meanFitnessScore = (
                    self.meanFitnessScore * self.numSequencesScored + np.sum(scores)
                ) / (aDNALength + self.numSequencesScored)
                self.numSequencesScored += aDNALength

        return score

    # Returns a node Number N based on in-order search. [0-N)
    def getNode(self, objective):
        nodeCount = 0
        return self.rootNode.getNode(objective, nodeCount)

    # Set positions a given a node and ID where is has to be
    def setNode(self, node, ID):

        print("node.ID = {} ID to change {}".format(node.ID, ID))
        if self.rootNode.ID == ID:
            self.rootNode = node
        else:
            self.rootNode.setNode(node, ID)

    # Get the parent node of a given ID and if it is the left child
    def getParent(self, ID):

        if self.rootNode.ID == ID:
            return {"isRootNode": True}

        return self.rootNode.getParent(ID)

        """
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
        """

    # Returns the number of nodes of the organism
    def countNodes(self):
        self.numNodes = self.rootNode.countNodes()

        return self.numNodes

    # Reset IDs of the full organism
    def resetIDs(self):
        firstID = 0
        self.rootNode.resetID(firstID)

    # Reset mean sequences scores and mean sequences scored
    def resetScores(self):
        self.meanFitnessScore = 0
        self.numSequencesScored = 0

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

        # Sort the array, so its always shown in the same order
        aDNA.sort()

        resultsFile = open(filename, "w+")

        for sDNA in aDNA:

            sDNA = sDNA.lower()
            sequenceLength = len(sDNA)

            # get all PSSM recognizers
            aPssmObjects = self.rootNode.getAllPssm([])

            # check position where it fits (position)
            pssmPositionScoreTable = []

            # Go through all PSSM objects to check where it fits
            for pssm in aPssmObjects:
                pssmLength = pssm.length
                maxScore = float("-inf")
                position = 0

                aPositions = []

                # Check every position possible on the sequence
                for pos in range(sequenceLength - pssmLength):
                    score = pssm.getScore(sDNA[pos : pos + pssmLength])

                    # Update max score if the actual score is acctually better
                    # Also check that the position is not overlapping other pssm object
                    if not self.isOverlapping(pos, pssmLength, pssmPositionScoreTable):
                        if score > maxScore:
                            maxScore = score
                            aPositions = [pos]

                        elif score == maxScore:
                            aPositions.append(pos)

                position = min(aPositions, key=lambda x: abs(x - pssm.positionMemory))

                pssm.rememberPosition(position)

                # Add to a table the ID, maxScore and position of a pssm object
                pssmPositionScoreTable.append((pssm.ID, maxScore, position, pssmLength))

            resultsFile.write("\n{}\n".format(sDNA))
            mapPositions = " " * len(sDNA)

            for ids, sc, pos, length in pssmPositionScoreTable:
                strId = str(ids)
                while len(strId) < length:
                    strId += "*"
                mapPositions = (
                    mapPositions[0:pos] + strId + mapPositions[pos + length :]
                )

            resultsFile.write(mapPositions + "\n")

        resultsFile.close()

    def printResult(self, sDNA):

        sDNA = sDNA.lower()
        sequenceLength = len(sDNA)

        # get all PSSM recognizers
        aPssmObjects = self.rootNode.getAllPssm([])

        # check position where it fits (position)
        pssmPositionScoreTable = []

        # Go through all PSSM objects to check where it fits
        for pssm in aPssmObjects:
            pssmLength = pssm.length
            maxScore = float("-inf")
            position = 0

            aPositions = []
            # Update max score if the actual score is acctually better
            # Also check that the position is not overlapping other pssm object

            # Check every position possible on the sequence
            for pos in range(sequenceLength - pssmLength):
                score = pssm.getScore(sDNA[pos : pos + pssmLength])

                if not self.isOverlapping(pos, pssmLength, pssmPositionScoreTable):
                    if score > maxScore:
                        maxScore = score
                        aPositions = [pos]

                    elif score == maxScore:
                        aPositions.append(pos)

            position = min(aPositions, key=lambda x: abs(x - pssm.positionMemory))

            pssm.rememberPosition(position)

            # Add to a table the ID, maxScore and position of a pssm object
            pssmPositionScoreTable.append((pssm.ID, maxScore, position, pssmLength))

        mapPositions = " " * len(sDNA)

        for ids, sc, pos, length in pssmPositionScoreTable:
            strId = str(ids)
            while len(strId) < length:
                strId += "*"
            mapPositions = mapPositions[0:pos] + strId + mapPositions[pos + length :]

        return "{}\n{}".format(sDNA, mapPositions)
