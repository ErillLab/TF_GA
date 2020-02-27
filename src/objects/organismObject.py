# Oragnism object
# It allocates the full data structure

class OrganismObject():

    # Organism constructor
    def __init__(self, ID = None, rootNode = None):
        self.ID = ID
        self.rootNode = rootNode
        self.bestFitnessScore = None

    def setRootNode(self, rootNode):
        self.rootNode = rootNode
    
    def getID(self):
        return self.ID


    # Mutates a part of the organism
    def mutate(self):
        print("Mutating organism {}".format(self.ID))
    
    # Return the fitness of the organism for a given DNA sequence
    def getBestAllFitness(self, sDNA):

        return 0

    # Return the total Fitness for an array of DNA sequences and the fitness method 
    def getScore(self, aDNA):
        
        return 0

    def print(self):
        print("Organism " + str(self.ID))
        print()
        self.rootNode.print(1)
