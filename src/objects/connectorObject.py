# C object
# Connects two nodes at an specific distance

from objects.nodeObject import Node
import numpy as np
import random


class ConnectorObject(Node):

    # Connector constructor gets mu, sigma and can get the two initial nodes.
    def __init__(self, mu, sigma, config, node1=None, node2=None):
        Node.__init__(self)
        self.mu = mu  # Mean discance
        self.sigma = sigma  # Variance between elements
        self.MUTATE_PROBABILITY_SIGMA = config["MUTATE_PROBABILITY_SIGMA"]
        self.MUTATE_PROBABILITY_MU = config["MUTATE_PROBABILITY_MU"]
        self.MUTATE_PROBABILITY_SWAP = config["MUTATE_PROBABILITY_SWAP"]
        self.TAU = config["TAU"]
        self.MUTATE_VARIANCE_SIGMA = config["MUTATE_VARIANCE_SIGMA"]
        self.MUTATE_VARIANCE_MU = config["MUTATE_VARIANCE_MU"]
        self.PLACEMENT_OPTIONS = config["PLACEMENT_OPTIONS"]

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

    def setPositionMemory(self, pos):
        self.node1.setPositionMemory(pos)
        self.node2.setPositionMemory(pos)

    # Counts the number of nodes below the node (including itself)
    def countNodes(self):
        return self.node1.countNodes() + self.node2.countNodes() + 1

    # Get a specific node based on a count and the objective node
    def getNode(self, objective, nodeCount):
        leftNodes = self.node1.countNodes()
        nodeNumber = leftNodes + nodeCount

        if nodeNumber == objective:
            # We are on the node we are searching
            return self
        elif nodeNumber < objective:
            # The node is on the right side
            return self.node2.getNode(objective, nodeNumber + 1)
        elif nodeNumber > objective:
            # The node is on the left side
            return self.node1.getNode(objective, nodeCount)

    # Get the parent node of a given ID and if it is the left child
    def getParent(self, ID):

        if self.node1.ID == ID:
            # Return itself and specify child is on left side
            return {"isRootNode": False, "self": self, "isLeftSide": True}

        elif self.node2.ID == ID:
            # Return itself and specify child is on right side
            return {"isRootNode": False, "self": self, "isLeftSide": False}

        checkNode = None
        checkNode = self.node1.getParent(ID)
        if checkNode is None:
            checkNode = self.node2.getParent(ID)
        return checkNode

    # returns an array of all pssm objects of the organism
    def getAllPssm(self):
        return self.node1.getAllPssm() + self.node2.getAllPssm()

    # calculate organism pacement on a sequence 
    """This is a position/score propagation method, defined for connector 
       objects.
       It is invoked by the placement method in the organism, for the root
       connector object, and calls itself recursively.
       The function calls itself until reaching terminal connectors, which call
       onto PSSM objects. 
       At that point, the call to the getPlacement function in PSSM nodes leads
       to the evaluation of the PSSM node across all the sequence, and it
       returns the score/position pairs, sorted by descending score.
       
       The connector function then propagates this up, taking the middle 
       position between both PSSMs and adding the connector energy to the 
       energies provided by the PSSMs.
       
       The connector determines (i.e. freezes) the PSSM locations, adding them
       to the block list that is passed as a parameter.
       
       Further connector objects proceed in the same manner, computing middle
       distance and adding their energy contribution, until the energy reaches
       the root node, and is returned as the fitness for the organism.

       The energy contribution of each connector is: EN1 + EN2 + Tau * EC, 
       where EN1 is the energy of its daugher element 1, and EN2 that of 
       daughter element 2. The EC connector energy component is an exponential 
       function controlled by the difference in the observed distance of the 
       elements  of the connector with respect to an ideal mean distance (mu), 
       and modulated by a dispersion parameter (sigma). Tau controls the 
       "weight" of the connector contribution to energy.
    """
    def getPlacement(self, sDNA, sDNAlen, blocks, blockers):
        # This tau shows how much value we give to the connector fit
        tau = self.TAU
        
        #ask daughter nodes what their placement is
        #placement call will return a vector of positions (for PSSMs) and their
        #scores (sorted descending by score) , as well as an udpated 
        #block/blocker vector
        node1 = self.node1.getPlacement(sDNA, sDNAlen, blocks, blockers)
        node2 = self.node2.getPlacement(sDNA, sDNAlen, blocks, blockers)

        #precompute connector energy term (not dependent on PSSM placement)
        logterm = np.log10(10 + self.sigma ** 2)

        maxenergy=-np.inf
        maxposition=0
        max1=0
        max2=0
        placecnt = 0
        placeopt = self.PLACEMENT_OPTIONS
        confset=False
        #iterate over all possible configurations of sub-node placements
        #and determine the optimal one. this goes on for a _minimum_ number of 
        #iterations (user-defined), but continues if a satisfactory 
        #configuration has not been found (for instance, because
        #all tried configurations generated conflicts with blocked positions)
        #the two for loops iterate over possible subset of configurations
        #if a valid solution is not found, they iterate over next subset
        #this will take place when at least one of the nodes is a PSSM and has
        #returned more than one possible position (although in theory
        #connectors could also return more than one position)
        while not confset:
            #set range according to availability of positions
            if len(node1['pspairs']) > placeopt:
                range1st=placecnt
                range1nd=placecnt+placeopt
            else:
                range1st=0
                range1nd=1
            if len(node2['pspairs']) > placeopt:
                range2st=placecnt
                range2nd=placecnt+placeopt
            else:
                range2st=0
                range2nd=1
                
            for n1count in range(range1st,range1nd):
                for n2count in range(range2st,range2nd):
                    #compute connector energy terms that depend on PSSM placement
                    numerator = (self.mu - (node2['pspairs'][n2count]['pos'] \
                                 - node1['pspairs'][n1count]['pos'])) ** 2
                    exponent = -1.0 * numerator / (1 + 2 * (self.sigma ** 2))
                    expterm = np.exp(exponent)
                    #compute additive connector energy term
                    eConnector = (tau / logterm) * expterm
                    #submodel energy: additive
                    energy = (node1['pspairs'][n1count]['energy'] \
                              + node2['pspairs'][n2count]['energy']) \
                              + eConnector
                    #submodel position: average of daughter positions
                    position = (node1['pspairs'][n1count]['pos'] \
                                + node2['pspairs'][n2count]['pos']) / 2
                        
                    
                    #if BOTH nodes are pssms, avoid pssm overlapping 
                    if self.node1.isPSSM() and self.node2.isPSSM():
                        #determine largest PSSM
                        p_len = self.node1.getLength() \
                            if self.node1.getLength() > self.node2.getLength() \
                                else self.node2.getLength() 
                        #determine if overlap exists, skip combo if there is
                        if abs(round(node2['pspairs'][n2count]['pos']) \
                               - round(node1['pspairs'][n1count]['pos'])) < p_len:
                            continue
                    
                    #if ONE of the nodes is a PSSM, make sure its position does
                    #not overlap with blocked positions
                    if self.node1.isPSSM():
                        p_len = self.node1.getLength()
                        startpos1=round(node1['pspairs'][n1count]['pos'])\
                                  - round(p_len/2)
                        blocked=False
                        for jpos in range(startpos1,startpos1+p_len):
                            if jpos in blocks:
                                blocked=True
                                break
                        #if position has been blocked, skip this combo
                        if blocked:
                            continue
                    if self.node2.isPSSM():
                        p_len = self.node2.getLength()
                        startpos2=round(node2['pspairs'][n2count]['pos'])\
                                  - round(p_len/2)
                        blocked=False
                        for jpos in range(startpos2,startpos2+p_len):
                            if jpos in blocks:
                                blocked=True
                                break
                        #if position has been blocked, skip this combo
                        if blocked:
                            continue                
                    #if this energy is the best so far, annotate it
                    #set confset to true, since we have obtained a valid
                    #configuration
                    if energy > maxenergy:
                        maxenergy = energy
                        maxposition = position
                        max1 = n1count
                        max2 = n2count
                        confset=True
            
                #increase base configuration counter
                placecnt=placecnt+1        
            
            # print('Max1: ',node1['pspairs'][max1]['pos'])        
            # print('Max2: ',node2['pspairs'][max2]['pos'])        
        
        #block the positions of this connector's PSSMs
        if self.node1.isPSSM():
            n1length=self.node1.getLength()
            blockstartpos1=round(node1['pspairs'][max1]['pos'])\
                           - round(n1length/2)
            for blockade in range(n1length):
                blocks.append(blockstartpos1+blockade)
                blockers.append(self.node1.ID)        
        if self.node2.isPSSM():
            n2length=self.node2.getLength()
            blockstartpos2=round(node2['pspairs'][max2]['pos'])\
                           - round(n2length/2)
            for blockade in range(n2length):
                blocks.append(blockstartpos2+blockade)
                blockers.append(self.node2.ID)        
        
        pair={'pos' : maxposition, 'energy' : maxenergy}

        return {'pspairs': [pair], 'blocked' : blocks, 'blocker' : blockers}

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
        # print("Mutating Connector..." + str(self.ID))
        if random.random() < self.MUTATE_PROBABILITY_SIGMA:
            # Alter sigma
            self.sigma += random.randint(
                -self.MUTATE_VARIANCE_SIGMA, self.MUTATE_VARIANCE_SIGMA
            )

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
        print(
            "   |" * distance
            + " - C"
            + str(self.ID)
            + " m: {} s: {}".format(self.mu, self.sigma)
        )
        self.node1.print(distance + 1)
        self.node2.print(distance + 1)

    # Exports Connector data to the given file
    def export(self, exportFile, level):
        exportFile.write(
            "\n"
            + "   |" * level
            + " - C"
            + str(self.ID)
            + " m: {} s: {}".format(self.mu, self.sigma)
        )
        self.node1.export(exportFile, level + 1)
        self.node2.export(exportFile, level + 1)

    def isConnector(self):
        return True

    def isPSSM(self):
        return False
