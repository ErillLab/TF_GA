# C object
# Connects two nodes at an specific distance

from objects.nodeObject import Node
import numpy as np
import random 
class ConnectorObject(Node):

    # Connector constructor gets mu, sigma and can get the two initial nodes.
    def __init__(self, mu, sigma, config, node1 = None, node2 = None):
        Node.__init__(self)
        self.mu = mu            # Mean discance
        self.sigma = sigma      # Variance between elements
        self.MUTATE_PROBABILITY_SIGMA = config["MUTATE_PROBABILITY_SIGMA"]
        self.MUTATE_PROBABILITY_MU = config["MUTATE_PROBABILITY_MU"]
        self.MUTATE_PROBABILITY_SWAP = config["MUTATE_PROBABILITY_SWAP"]
        self.TAU = config["TAU"]
        self.MUTATE_VARIANCE_SIGMA = config["MUTATE_VARIANCE_SIGMA"]
        self.MUTATE_VARIANCE_MU = config["MUTATE_VARIANCE_MU"]

        self.node1 = node1
        self.node2 = node2

    # Setters
    def setMu(self, mu):
        self.mu = mu

    def setSigma(self, sigma):
        self.sigma = sigma

    def setNode1(self, node1):
        self.node1 = node1

    def setNode2(self, node2):
        self.node2 = node2

    # Counts the number of nodes below the node (including itself)
    def countNodes(self):
        return self.node1.countNodes() + self.node2.countNodes() + 1
    
    # Get a specific node based on a count and the objective node
    def getNode(self, objective, nodeCount):
        leftNodes = self.node1.countNodes()
        nodeNumber = leftNodes + nodeCount

        if(nodeNumber == objective):
            # We are on the node we are searching
            return self
        elif(nodeNumber < objective):
            #The node is on the right side
            return self.node2.getNode(objective, nodeNumber + 1)
        elif(nodeNumber > objective):
            # The node is on the left side
            return self.node1.getNode(objective, nodeCount)

    # Get the parent node of a given ID and if it is the left child
    def getParent(self, ID):

        if self.node1.ID == ID:
            # Return itself and specify child is on left side
            return { "isRootNode": False, "self": self, "isLeftSide": True }

        elif self.node2.ID == ID:
            # Return itself and specify child is on right side
            return { "isRootNode": False, "self": self, "isLeftSide": False }

        checkNode = None
        checkNode = self.node1.getParent(ID)
        if checkNode == None:
            checkNode = self.node2.getParent(ID)
        return checkNode

    # returns an array of all pssm objects of the organism 
    def getAllPssm(self, aPssms):
        aPssms = self.node1.getAllPssm(aPssms)
        aPssms = self.node2.getAllPssm(aPssms)
        return aPssms

    # calculate best all strategy for pssm position based on a table of pssm scores.
    def getBestAll(self, table):
        # This tau shows how much value we give to the connector fit 
        tau = self.TAU
        
        eNode1, distance1 = self.node1.getBestAll(table)
        eNode2, distance2 = self.node2.getBestAll(table)
        
        numerator = (self.mu - abs(distance1-distance2))**2
        exponent = - numerator / (1 + 2 * (self.sigma ** 2))
        eConnector = (tau/np.log10(10 + self.sigma ** 2)) * np.exp(exponent)
        
        #energy = (eNode1 + eNode2) + eConnector
        #print("N1:{} N2:{} C:{} SIGMA:{} MU {}\n".format(eNode1, eNode2, eConnector, self.sigma, self.mu))
        #energy = (eNode1 + eNode2) * eConnector
        energy = max(eNode1*eConnector +eNode2, eNode2*eConnector +eNode1)
        distance = (distance1 + distance2) / 2

        return energy, distance 

    # Sets the node on a given ID
    def setNode(self, node, ID):
        
        if self.node1.ID == ID:
            self.node1 = node
        elif self.node2 == ID:
            self.node2 = node
        else:
            self.node1.setNode(node, ID)    
            self.node2.setNode(node, ID)    

    def resetID(self, newID):
        newID = self.node1.resetID(newID)
        self.ID = newID
        newID = self.node2.resetID(newID + 1)
        return newID

    # mutation for a connector
    def mutate(self, orgFactory):
        #print("Mutating Connector..." + str(self.ID))
        if random.random() < self.MUTATE_PROBABILITY_SIGMA:
            # Alter sigma
            self.sigma += random.randint(-self.MUTATE_VARIANCE_SIGMA, self.MUTATE_VARIANCE_SIGMA)

        if random.random() < self.MUTATE_PROBABILITY_MU:
            # Alter mu
            self.mu += random.randint(-self.MUTATE_VARIANCE_MU, self.MUTATE_VARIANCE_MU)

        if random.random() < self.MUTATE_PROBABILITY_SWAP:
            # Swap connectors
            tmpNode = self.node1
            self.node1 = self.node2
            self.node2 = tmpNode

    # It prints the connector mu, sigma values and its children values in tree structure
    def print(self, distance):
        self.node1.print(distance + 1)
        print("----    "*distance +"Connection: "+str(self.ID)+" mu: {} sigma: {}".format(self.mu, self.sigma))
        self.node2.print(distance + 1)


    # Exports Connector data to the given file
    def export(self, exportFile, level):
        self.node1.export(exportFile, level + 1)
        exportFile.write("\n"+"----    "*level +"Connection "+str(self.ID)+": mu: {} sigma: {}".format(self.mu, self.sigma))
        self.node2.export(exportFile, level + 1)
    
    def isConnector(self): 
        return True

    def isPSSM(self):
        return False
