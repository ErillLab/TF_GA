from objects.organismObject import OrganismObject
from objects.connectorObject import ConnectorObject
from objects.pssmObject import PssmObject

import random
import numpy

class OrganismFactory:

    def __init__(self, confOrg, confOrgFac, confCon, confPssm):

        self.ID = 0
        self.INITIAL_CONNECTOR_PROBABILITY = confOrgFac["INITIAL_CONNECTOR_PROBABILITY"]
        self.REDUCER_PROBABILITY_FACTOR = confOrgFac["REDUCER_PROBABILITY_FACTOR"]

        self.MIN_MU = confOrgFac["MIN_MU"]
        self.MAX_MU = confOrgFac["MAX_MU"]

        self.MIN_SIGMA = confOrgFac["MIN_SIGMA"]
        self.MAX_SIGMA = confOrgFac["MAX_SIGMA"]
        # The number of positions, PSSM object can recognize
        
        self.PWM_LENGTH = confOrgFac["PWM_LENGTH"]
        
        self.PWM_PROBABILITY_STEP = confOrgFac["PWM_PROBABILITY_STEP"] # It should be a BASE_PROBABILITY divisor Ex: 1, 2, 4, 5, 10, 25...
        self.PWM_PROBABILITY_BASE = confOrgFac["PWM_PROBABILITY_BASE"]
        self.PWM_PROBABILITY_DECIMALS = confOrgFac["PWM_PROBABILITY_DECIMALS"]

        self.confOrg = confOrg
        self.confCon = confCon
        self.confPssm = confPssm

    # This should be a function so all the count of IDs, including assigned
    # outside the class, keep consistency between all organisms
    def getID(self):
        self.ID +=1
        return self.ID
    
    
    # It returns a full organism datastructure
    def getOrganism(self):
        
        newOrganism = OrganismObject(self.getID(), self.confOrg)
        rootNode = None

        # Based on a random probability, we assign a connector or a PSSM object
        # to the root Node
        if(random.random() < self.INITIAL_CONNECTOR_PROBABILITY):
            rootNode = self.createConnection(self.INITIAL_CONNECTOR_PROBABILITY*self.REDUCER_PROBABILITY_FACTOR)
        else:
            rootNode = self.createPSSM(self.PWM_LENGTH)

        newOrganism.setRootNode(rootNode)

        return newOrganism
    
    
    # It returns a connector object with its nodes assigned depending on connectionProbability
    def createConnection(self, connectionProbability):

        # Assign a random value to mu and sigma
        mu = random.randint(self.MIN_MU,self.MAX_MU)
        sigma = random.randint(self.MIN_SIGMA, self.MAX_SIGMA)

        # Create the new connection
        newConnection = ConnectorObject(mu, sigma, self.confCon)

        # Set the connection node to a connector or PSSM object, based on a 
        # random probability. If connector object is selected, probability of 
        # getting another connector is reduced by PROBABILITY_REDUCED_FACTOR
        node1 = None
        if (random.random() < connectionProbability):
            node1 = self.createConnection(connectionProbability*self.REDUCER_PROBABILITY_FACTOR)
        else:
            node1 = self.createPSSM(self.PWM_LENGTH)

        newConnection.setNode1(node1)

        node2 = None
        if (random.random() < connectionProbability):
            node2 = self.createConnection(connectionProbability*self.REDUCER_PROBABILITY_FACTOR)
        else:
            node2 = self.createPSSM(self.PWM_LENGTH)

        newConnection.setNode2(node2)

        return newConnection
    # It return a PSSM object with a specific length
    def createPSSM(self, length):

        pwm = []
        # Generate as much as needed
        for i in range(length):
            pwm.append(self.getPwmColumn())
            
        return PssmObject(numpy.array(pwm), self.confPssm)
    
    def getPwmColumn(self):
        
        initialProbability = self.PWM_PROBABILITY_BASE / self.PWM_PROBABILITY_STEP
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
        numpyProbabilities = numpy.array(probabilities)*self.PWM_PROBABILITY_STEP*(1/self.PWM_PROBABILITY_BASE)
        probabilities = numpyProbabilities.tolist()
        
        # Return object with "decimals" decimals probability to each base 
        return {
                "a":round(probabilities[0],decimals),
                "g":round(probabilities[1],decimals),
                "c":round(probabilities[2],decimals),
                "t":round(probabilities[3],decimals)
                } 

