# Oragnism object
# It allocates the full data structure

class OrganismObject():

    # Organism constructor
    def __init__(self, ID = None, rootNode = None):
        self.ID = ID
        self.rootNode = rootNode
        self.bestFitnessScore = None
    
    # Setters an getters
    def setRootNode(self, rootNode):
        self.rootNode = rootNode
    
    def getID(self):
        return self.ID
    
    def setID(self, iD):
        self.ID = iD

    # TODO: Mutates a part of the organism
    def mutate(self):
        print("Mutating organism {}".format(self.ID))
    
    # TODO: Return the fitness of the organism for a given DNA sequence
    def getBestAllFitness(self, sDNA):

        return 0

    # TODO: Return the total Fitness for an array of DNA sequences and the fitness method 
    def getScore(self, aDNA):
        
        return 0

    # Returns a node Number N based on in-order search. [0-N)
    def getNode(self, objective):
        nodeCount = 0
        return self.rootNode.getNode(objective, nodeCount)

    # Set positions a given a node and ID where is has to be
    def setNode(self, node, ID):
        # CANNOT BE ID BC SECOND IS OVERWRITED

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

    # Prints the whole tree data structure
    def print(self):
        print()
        print("Organism " + str(self.ID))
        self.rootNode.print(1)
