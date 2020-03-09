# C object
# Connects two nodes at an specific distance

from objects.nodeObject import Node

class ConnectorObject(Node):

    # Connector constructor gets mu, sigma and can get the two initial nodes.
    def __init__(self, mu, sigma, node1 = None, node2 = None):
        Node.__init__(self)
        self.mu = mu            # Mean discance
        self.sigma = sigma      # Variance between elements
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
            return [False, self, True]
        elif self.node2.ID == ID:
            # Return itself and specify child is on right side
            return [False, self, False]
        checkNode = None
        checkNode = self.node1.getParent(ID)
        if checkNode == None:
            checkNode = self.node2.getParent(ID)
        return checkNode




    # Sets the node on a given ID
    def setNode(self, node, ID):
        
        if self.node1.ID == ID:
            self.node1 = node
        elif self.node2 == ID:
            self.node2 = node
        else:
            self.node1.setNode(node, ID)    
            self.node2.setNode(node, ID)    


    # TODO: implement
    def mutate(self):
        print("Mutating Connector..." + str(self.ID))
        self.node1.mutate()
        self.node2.mutate()

    # It prints the connector mu, sigma values and its children values in tree structure
    def print(self, distance):
        self.node1.print(distance + 1)
        print("----"*distance +"Connection: "+str(self.ID)) #+str(self.mu) + " " + str(self.sigma))
        self.node2.print(distance + 1)
