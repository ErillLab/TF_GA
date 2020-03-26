from objects.organismObject import OrganismObject
from objects.connectorObject import ConnectorObject
from objects.pssmObject import PssmObject

import random
import numpy
import json
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
    #Import Organism from file
    def importOrganisms(self, fileName):
        organismList = []
        organismJSON = {}
        with open(fileName) as json_file:
            organismJSON = json.load(json_file)
        
        for organism in organismJSON:



            newOrganism = OrganismObject(self.getID(), self.confOrg)
            rootNode = None

            if organism["rootNode"]["objectType"] == "pssm":
                rootNode = self.importPSSM(organism["rootNode"])
            else:
                rootNode = self.importConnector(organism["rootNode"])
            
            newOrganism.setRootNode(rootNode)

            organismList.append(newOrganism)

        return organismList

    # Import Connector from JSON object
    def importConnector(self, connector):
        newConnector = ConnectorObject(connector["mu"], connector["sigma"], self.confCon)
        
        node1 = None
        if connector["node1"]["objectType"] == "pssm":
            node1 = self.importPSSM(connector["node1"])
        else:
            node1 = self.importConnector(connector["node1"])

        newConnector.setNode1(node1)

        node2 = None
        if connector["node2"]["objectType"] == "pssm":
            node2 = self.importPSSM(connector["node2"])
        else:
            node2 = self.importConnector(connector["node2"])

        newConnector.setNode2(node2)

        return newConnector


    # Import PSSM from JSON object
    def importPSSM(self, pssm):
        return PssmObject(numpy.array(pssm["pwm"]), self.confPssm)


    # Export a list of organisms to JSON format
    def exportOrganisms(self, aOrganisms, filename):
        listJsonOrganisms = []
        for oOrganism in aOrganisms:
            organism = {}
            if oOrganism.rootNode.isConnector():
                organism["rootNode"] = self.exportConnector(oOrganism.rootNode)
            else:
                organism["rootNode"] = self.exportPSSM(oOrganism.rootNode)
            listJsonOrganisms.append(organism)


        with open(filename, "w+") as json_file:
            json.dump(listJsonOrganisms, json_file, indent=2)
    

    #Export connector object
    def exportConnector(self, oConnector):
        connector = {}
        connector["objectType"] = "connector"
        connector["mu"] = oConnector.mu
        connector["sigma"] = oConnector.sigma
        # Check if its pssm
        if oConnector.node1.isConnector():
            connector["node1"] = self.exportConnector(oConnector.node1)
        else:
            connector["node1"] = self.exportPSSM(oConnector.node1)
        

        if oConnector.node2.isConnector():
            connector["node2"] = self.exportConnector(oConnector.node2)
        else:
            connector["node2"] = self.exportPSSM(oConnector.node2)

        return connector


    # Export PSSM object
    def exportPSSM(self, oPssm):
        pssm = {}
        pssm["objectType"] = "pssm"
        pssm["pwm"] = oPssm.pwm.tolist()
        return pssm
