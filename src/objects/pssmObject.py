# P object
# Saves al specific PSSM data structure

from objects.nodeObject import Node
import numpy as np
import random

class PssmObject(Node):

        
    # PSSM Constructor
    def __init__(self, pwm, config):
        Node.__init__(self)
        self.length = len(pwm) #length of the numpy array
        self.pwm = pwm # numpy array of dictionaries
        self.pssm = None
        self.MUTATE_PROBABILITY_RANDOM_COL = config["MUTATE_PROBABILITY_RANDOM_COL"]
        self.MUTATE_PROBABILITY_FLIP_COLS = config["MUTATE_PROBABILITY_FLIP_COL"]
        self.MUTATE_PROBABILITY_FLIP_ROWS = config["MUTATE_PROBABILITY_FLIP_ROW"]
        self.MUTATE_PROBABILITY_SHIFT_LEFT = config["MUTATE_PROBABILITY_SHIFT_LEFT"]
        self.MUTATE_PROBABILITY_SHIFT_RIGHT = config["MUTATE_PROBABILITY_SHIFT_RIGHT"]
        self.PSEUDO_COUNT = config["PSEUDO_COUNT"]
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

    # Mutate PSSM object
    def mutate(self, orgFactory):

        #print("Mutating PSSM..." + str(self.ID))
        randNum = random.random()
        #print("Start: randum = {}\n".format(randNum)+str(self.pwm))

        if random.random() < self.MUTATE_PROBABILITY_RANDOM_COL:
            
            #Randomize PSSM column
            newCol = orgFactory.getPwmColumn()
            # Select a random col in self.pwm
            columnToUpdate = random.randint(0, self.length - 1)
            # Insert it in that position
            self.pwm[columnToUpdate] = newCol

        if random.random() < self.MUTATE_PROBABILITY_FLIP_COLS:
            # Flip two columns

            col1, col2 = random.sample(range(self.length),2)
            # Select two random columns and swap it
            tmpCol = self.pwm[col1]
            self.pwm[col1] = self.pwm[col2]
            self.pwm[col2] = tmpCol
            

        if random.random() < self.MUTATE_PROBABILITY_FLIP_ROWS:
            # Flip two rows
            # Save values of two rows and swap it (its gonna be one by one)
            bases = ['a','c','g','t']
            random.shuffle(bases)
            base1, base2 = bases[:2]
        
            # Swap rows
            for i in range(self.length):
                tmpBase = self.pwm[i][base1]
                self.pwm[i][base1] = self.pwm[i][base2]
                self.pwm[i][base2] = tmpBase

        if random.random() < self.MUTATE_PROBABILITY_SHIFT_LEFT:
            #Shift to right/left
            self.pwm = np.roll(self.pwm, -1)

        if random.random() < self.MUTATE_PROBABILITY_SHIFT_RIGHT:
            #Shift to right/left
            self.pwm = np.roll(self.pwm, 1)

        self.recalculatePSSM()

    # Calculate self.pssm based on self.pwm
    def recalculatePSSM(self):
        tmpPSSM = []
        for column in self.pwm:
            # From pwm to pssm
            # log2(base/0.25) = log2(4.0*base)
            decimals = 2
            tmpBases = []
            # cast to float so round function does not become crazy
            tmpBases.append(float(np.log2(4.0 * column["c"] + self.PSEUDO_COUNT)))
            tmpBases.append(float(np.log2(4.0 * column["t"] + self.PSEUDO_COUNT)))
            tmpBases.append(float(np.log2(4.0 * column["g"] + self.PSEUDO_COUNT)))
            tmpBases.append(float(np.log2(4.0 * column["a"] + self.PSEUDO_COUNT)))

            tmpPSSM.append({
                "c":round(tmpBases[0], decimals), 
                "t":round(tmpBases[1], decimals), 
                "g":round(tmpBases[2], decimals), 
                "a":round(tmpBases[3], decimals)
                })
        self.pssm = np.array(tmpPSSM)

    # Searchs himself on the table and returns position and score
    def getBestAll(self, table):
        for ID, score, position, length in table:
            if self.ID == ID:
                #Maybe add half the length so the position is centered 
                return score, float(position) 
        return 0, 0

    # Adds himself as a pssm recognizer    
    def getAllPssm(self, aPssms):
        aPssms.append(self)
        return aPssms

    # returns a score to that DNA secuence
    def getScore(self, sDNA):
        # gets a score from pssm
        score = 0
        if len(sDNA) != len(self.pssm):
            print("Not a valid length!")
            return -3000

        for i in range(len(sDNA)):
            score += self.pssm[i][sDNA[i]]
        return score

    # Nodes cannot be setted from recognizer objects
    def setNode(self, node, ID):
        return 0

    def resetID(self, newID):
        self.ID = newID
        return newID + 1

    # Print PSSM object (Matrix)
    def print(self, distance):
        print("----    "*distance + "Node "+str(self.ID))
        print(str(self.pssm))
    
    # Exports pssm to a file
    def export(self, exportFile, level):
        exportFile.write("\n"+"----    "*level + "Node "+str(self.ID))
        exportFile.write("\n"+str(self.pssm))
