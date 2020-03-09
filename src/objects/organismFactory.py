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
# The number of positions, PSSM object can recognize
PWM_LENGTH = 4


STEP = 5.0 # It should be a BASE_PROBABILITY divisor Ex: 1, 2, 4, 5, 10, 25...
BASE_PROBABILITY = 100.0


class OrganismFactory:

    def __init__(self):
        self.ID = 0
        
    # This should be a function so all the count of IDs, including assigned
    # outside the class, keep consistency between all organisms
    def getID(self):
        self.ID +=1
        return self.ID
    
    
    # It returns a full organism datastructure
    def getOrganism(self):
        
        newOrganism = OrganismObject(self.getID())
        rootNode = None

        # Based on a random probability, we assign a connector or a PSSM object
        # to the root Node
        if(random.random() < INITIAL_PROBABILITY):
            rootNode = self.createConnection(INITIAL_PROBABILITY*PROBABILITY_REDUCED_FACTOR)
        else:
            rootNode = self.createPSSM(PWM_LENGTH)

        newOrganism.setRootNode(rootNode)

        return newOrganism
    
    
    # It returns a connector object with its nodes assigned depending on connectionProbability
    def createConnection(self, connectionProbability):

        # Assign a random value to mu and sigma
        mu = random.randint(MIN_MU,MAX_MU)
        sigma = random.randint(MIN_SIGMA, MAX_SIGMA)

        # Create the new connection
        newConnection = ConnectorObject(mu, sigma)

        # Set the connection node to a connector or PSSM object, based on a 
        # random probability. If connector object is selected, probability of 
        # getting another connector is reduced by PROBABILITY_REDUCED_FACTOR
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
    # It return a PSSM object with a specific length
    def createPSSM(self, length):

        pwm = []
        # Generate as much as needed
        for i in range(length):
            pwm.append(self.getPwmColumn())
            
        return PssmObject(numpy.array(pwm))
    
    def getPwmColumn(self):
        
        initialProbability = BASE_PROBABILITY / STEP
        probabilities = []
        # Left probability is
        leftProbability = initialProbability
        # Minimum and maximum number of probabilities to be generated 
        minProbability = 0
        maxProbability = 4
        # Number of decimals on the probability
        decimals = 2
        # Generate 4 random probabilities out of initialProbability, one for each base
        
        # Add a probability while we have less than 3 and and total probability is not 1
        while(leftProbability > minProbability and len(probabilities) < maxProbability -1):
            newProbability = random.randint(0, leftProbability)
            probabilities.append(float(newProbability))
            leftProbability -= newProbability
        # Add the last probability or fill with 0 probability
        if(leftProbability > 0):
            probabilities.append(initialProbability - sum(probabilities))
        else:
            while(len(probabilities)< maxProbability):
                probabilities.append(0.0)


        # Shuffle the array is needed so high probability is not always on first positions
        random.shuffle(probabilities)

        # Transform probabilities array from integer [0-(BASE_PROBABILITY / STEP)] to complementary float probabilities [0.0-1.0]
        numpyProbabilities = numpy.array(probabilities)*STEP*(1/BASE_PROBABILITY)
        probabilities = numpyProbabilities.tolist()
        
        # Return object with "decimals" decimals probability to each base 
        return {
                "a":round(probabilities[0],decimals),
                "g":round(probabilities[1],decimals),
                "c":round(probabilities[2],decimals),
                "t":round(probabilities[3],decimals)
                } 

