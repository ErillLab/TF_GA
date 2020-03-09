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
    # TODO: Mutate PSSM object
    def mutate(self):
        print("Mutating..." + str(self.ID))
        return 0

    # TODO: Calculate self.pssm based on self.pwm
    def recalculatePSSM(self):
        self.pssm = self.pssm

    # TODO: Print PSSM object (Matrix)
    def print(self, distance):
        print("--"*distance + "Node "+ str(self.pwm))

