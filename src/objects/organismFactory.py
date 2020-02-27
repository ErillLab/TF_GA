from objects.organismObject import OrganismObject
from objects.connectorObject import ConnectorObject
from objects.pssmObject import PssmObject

import random
import numpy

INITIAL_PROBABILITY = 0.9
PROBABILITY_REDUCED_FACTOR = 0.6

MIN_MU = 0
MAX_MU = 50

MIN_SIGMA = 0
MAX_SIGMA = 25

PWM_LENGTH = 4

STEP = 5.0 # It should be a 100 divisor Ex: 1, 2, 4, 5, 10, 25...
BASE_PROBABILITY = 100.0


class OrganismFactory:

    def __init__(self):
        self.ID = 0
        
    def getID(self):
        self.ID +=1
        return self.ID
    
    def getOrganism(self):
        
        newOrganism = OrganismObject(self.getID())
        
        rootNode = None
        if(random.random() < INITIAL_PROBABILITY):
            rootNode = self.createConnection(INITIAL_PROBABILITY*PROBABILITY_REDUCED_FACTOR)
        else:
            rootNode = self.createPSSM(PWM_LENGTH)

        newOrganism.setRootNode(rootNode)

        return newOrganism

    def createConnection(self, connectionProbability):
        # Create the new connection
        mu = random.randint(MIN_MU,MAX_MU)
        sigma = random.randint(MIN_SIGMA, MAX_SIGMA)

        newConnection = ConnectorObject(mu, sigma)

        # Set the connection nodes
        node1 = None
        if (random.random() < connectionProbability):
            node1 = self.createConnection(connectionProbability*PROBABILITY_REDUCED_FACTOR)
        else:
            node1 = self.createPSSM(PWM_LENGTH)

        newConnection.setNode1(node1)

        node2 = None
        if (random.random() < connectionProbability):
            node2 = self.createConnection(connectionProbability*PROBABILITY_REDUCED_FACTOR)
        else:
            node2 = self.createPSSM(PWM_LENGTH)

        newConnection.setNode2(node2)

        return newConnection

    def createPSSM(self, length):

        pwm = []
        # Generate as much 
        for i in range(length):
            pwm.append(self.getPwmColumn())
            
        return PssmObject(numpy.array(pwm))
    
    def getPwmColumn(self):
        
        initialProbability = BASE_PROBABILITY / STEP
        probabilities = []
        leftProbability = initialProbability
        minProbability = 0
        maxProbability = 4
        decimals = 2
        # Generate 4 random probabilities out of 100
        
        while(leftProbability > minProbability and len(probabilities) < maxProbability -1):
            newProbability = random.randint(0, leftProbability)
            probabilities.append(float(newProbability))
            leftProbability -= newProbability

        if(leftProbability > 0):
            probabilities.append(initialProbability - sum(probabilities))
        else:
            while(len(probabilities)< maxProbability):
                probabilities.append(0.0)

        random.shuffle(probabilities)
        numpyProbabilities = numpy.array(probabilities)*STEP*(1/BASE_PROBABILITY)
        probabilities = numpyProbabilities.tolist()
        return {
                "a":round(probabilities[0],decimals),
                "g":round(probabilities[1],decimals),
                "c":round(probabilities[2],decimals),
                "t":round(probabilities[3],decimals)
                } 

