# P object
# Saves al specific PSSM data structure

from objects.nodeObject import Node

class PssmObject(Node):

        
    # PSSM Constructor
    def __init__(self, pwm):
        Node.__init__(self)
        self.length = len(pwm) #length of the numpy array
        self.pwm = pwm # numpy array of dictionaries
        self.pssm = None

        # It first calculates PSSM Matrix based on  pwm
        self.recalculatePSSM()

    # returns itself as a node
    def countNodes(self):
        return 1

    # return itself if its the node, otherwise return None
    def getNode(self, objective, nodeCount):
        return self if objective == nodeCount else None

    #  pssm objects cannot be parents
    def getParent(self, ID):
        return None

    # TODO: Mutate PSSM object
    def mutate(self):
        print("Mutating..." + str(self.ID))
        return 0

    # TODO: Calculate self.pssm based on self.pwm
    def recalculatePSSM(self):
        self.pssm = self.pssm

    # Nodes cannot be setted from recognizer objects
    def setNode(self, node, ID):
        return 0

    # TODO: Print PSSM object (Matrix)
    def print(self, distance):
        print("----"*distance + "Node "+str(self.ID))#+ str(self.pwm))

