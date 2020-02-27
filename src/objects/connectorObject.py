# C object
# Connects two nodes at an specific distance

from objects.nodeObject import Node

class ConnectorObject(Node):

    # Connector constructor
    def __init__(self, mu, sigma, node1 = None, node2 = None):
        Node.__init__(self)
        self.mu = mu # Mean discance
        self.sigma = sigma # Variance between elements
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
    
    #Operational functions
    def doNothing(self):
        print("Doing nothing... from P Object.")

    def mutate(self):
        print("Mutating Connector..." + str(self.ID))
        self.node1.mutate()
        self.node2.mutate()
    def print(self, distance):
        print("--"*distance +"Cnnection: "+str(self.mu) + " " + str(self.sigma))
        self.node1.print(distance + 1)
        self.node2.print(distance+ 1)
