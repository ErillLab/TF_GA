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
        self.recalculatePSSM()

    def mutate(self):
        print("Mutating..." + str(self.ID))
        return 0

    def recalculatePSSM(self):
        self.pssm = self.pssm

    def doNothing(self):
        print("Doing nothing... from P Object.")

    def print(self, distance):
        print("--"*distance + "Node "+ str(self.pwm))

