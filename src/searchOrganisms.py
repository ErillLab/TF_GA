# Main execution

from Bio import SeqIO
import time
from objects.organismFactory import OrganismFactory
import random
import copy
import json


POPULATION_LENGTH = 0
DATASET_BASE_PATH_DIR = '' 
RESULT_BASE_PATH_DIR = '' 
POSITIVE_FILENAME  = ''
NEGATIVE_FILENAME  = ''
JSON_CONFIG_FILENAME = "config.json"
configOrganism = {}
configOrganismFactory = {}
configConnector = {}
configPssm = {}

MAX_SEQUENCES_TO_FIT = 0
MIN_ITERATIONS = 0
MIN_SCORE = 0

organismPopulation = []


positiveDataset = []
negativeDataset = []

THRESHOLD = 0.0
MUTATION_PROBABILITY = 0.0

def main():

    
    positiveDataset = readFastaFile(DATASET_BASE_PATH_DIR + POSITIVE_FILENAME)
    negativeDataset = readFastaFile(DATASET_BASE_PATH_DIR + NEGATIVE_FILENAME)


    # Generate initial population
    organismFactory = OrganismFactory(configOrganism, configOrganismFactory, configConnector, configPssm)

    for i in range(POPULATION_LENGTH):
        newOrganism = organismFactory.getOrganism()
        organismPopulation.append(newOrganism)

    iterations = 0

    maxScore = float("-inf")
    lastMaxScore = 0.0 
    bestOrganism = (None, 0.0)

    timeformat = "%Y-%m-%d--%H-%M-%S"
    # Main loop, it iterates until organisms do not get a significant change 
    # or MIN_ITERATIONS is reached.
    # it can be done untill you get a score over a value
    try:
        while not isFinished(END_WHILE_METHOD, iterations, maxScore, lastMaxScore):
        
            # Shuffle population & datasets
            random.shuffle(organismPopulation)
            random.shuffle(negativeDataset)
            #random.shuffle(positiveDataset) (HAS no effect)

            # Reset maxScore
            lastMaxScore = maxScore    
            maxScore = float("-inf")
            changedBestScore = False
            # Iterate over pairs of organisms
            for i in range(0, len(organismPopulation) - 1, 2):
                org1 = organismPopulation[i]
                org2 = organismPopulation[i+1]
            
                # Cross parents to get childs
                children = combineOrganisms(org1, org2, organismFactory)
            
                child1 = children["child1"]["child"]
                child2 = children["child2"]["child"]

                # Match with its closest child
                # There are 2 possible combinations p1-c1, p2-c2 & p1-c2, p2-c1
                # We select a combination based on a sum of similarities in combinations
                combination1 = children["child1"]["simOrg1"] + children["child2"]["simOrg2"] # Match the first parent to first child and second parent to second child 
                combination2 = children["child1"]["simOrg2"] + children["child2"]["simOrg1"] # Match the first parent to second child and second parent to first child

                pairChildren = []
            
                # Mutate children and parents with a probability pm

                if random.random() < MUTATION_PROBABILITY:
                    child1.mutate(organismFactory)
                if random.random() < MUTATION_PROBABILITY:
                    child2.mutate(organismFactory)


                if combination1 > combination2:
                    pairChildren.append((org1, child1))
                    pairChildren.append((org2, child2))
                else:
                    pairChildren.append((org1, child2))
                    pairChildren.append((org2, child1))

                # "Fight" two organisms (DO NOT take two organisms directly 
                # because it uses j index to reinsert the winner organisms 
                # into the global population)
                for j in range(len(pairChildren)):

                    firstOrganism = pairChildren[j][0]
                    secondOrganism = pairChildren[j][1]
                    p1 = firstOrganism.getScore(positiveDataset[:MAX_SEQUENCES_TO_FIT]) 
                    n1 = firstOrganism.getScore(negativeDataset[:MAX_SEQUENCES_TO_FIT])
                    #c1 = firstOrganism.getComplexity()

                    p2 = secondOrganism.getScore(positiveDataset[:MAX_SEQUENCES_TO_FIT])
                    n2 = secondOrganism.getScore(negativeDataset[:MAX_SEQUENCES_TO_FIT])
                    #c2 = secondOrganism.getComplexity()
                    
                    # Check negative scores (horrible matching)
                    #if(p1 < 0 or n1 < 0):
                    #    print("Negative values on {}!!Pos: {} Neg:{}".format(firstOrganism.ID, p1, n1))

                    #if(p2 < 0 or n2 < 0):
                    #    print("Negative values on {}!!POS: {} Neg:{}".format(secondOrganism.ID, p2, n2))


                    score1 = (p1 - n1)

                    score2 = (p2 - n2)

                
                    # print("ID1: {} Score1:{}/{} =  {} ID2: {} Score2: {}/{} = {}".format(firstOrganism.ID, p1, n1, score1, secondOrganism.ID, p2, n2, score2))
                
                    if(score1 > score2): # The first organism wins
                        # Set it back to the population
                        organismPopulation[i+j] = firstOrganism

                        # Check if its the max score in that iteration
                        if score1 > maxScore:
                            maxScore = score1

                        #Check if its the max score in the program
                        if score1 > bestOrganism[1]:
                            bestOrganism = (firstOrganism, score1)
                            changedBestScore = True

                    else: # The second organism wins
                        # Set it back to the population
                        organismPopulation[i+j] = secondOrganism

                        # Check if its the max score in that iteration
                        if score2 > maxScore:
                            maxScore = score2

                        #Check if its the max score in the program
                        if score2 > bestOrganism[1]:
                            bestOrganism = (secondOrganism, score2)
                            changedBestScore = True

                     
                    #END FOR j
            
                # END FOR i

            # Show IDs of final array
            #print("-"*10)
            print("Iter: {} Max Score: {} -  BestOrg: {} Score: {}".format(iterations, maxScore, bestOrganism[0].ID, bestOrganism[1]))
            if(changedBestScore):
                filename = "{}_{}".format(time.strftime(timeformat), bestOrganism[0].ID)
                exportOrganism(bestOrganism[0], positiveDataset, filename)
            #print("-"*10)
            iterations += 1
            # END WHILE

    except Exception as e:
        print("Exited program: \n{}\n".format(e))

    finally:
        print()
        print("-"*10)
        print("Best Organism: {}".format(bestOrganism[1]))
        bestOrganism[0].print()
        print("-"*10)

# Checks if main while loop is finished
# methods: 'Iterations', 'minScore', 'Threshold'
def isFinished(method, iterations, maxScore, lastMaxScore):

    if method.lower() == 'iterations':
        return iterations >= MIN_ITERATIONS

    elif method.lower() == 'minscore':
        return maxScore >= MIN_SCORE

    elif method.lower() == 'threshold':
        return abs(lastMaxScore - maxScore) <= THRESHOLD

    return True

def exportOrganism(organism, dataset, filename):
    
    organismFile = "{}{}_organism.txt".format(RESULT_BASE_PATH_DIR, filename)
    resultsFile = "{}{}_results.txt".format(RESULT_BASE_PATH_DIR, filename)

    organism.export(organismFile)
    organism.exportResults(dataset, resultsFile)



# Gets 2 organisms, and returns 2 children with format (child, similarity to parent 1, similarity to parent 2)
def combineOrganisms(organism1, organism2, organismFactory):
    # Save the number of nodes from the parents
    nNodesOrg1 = organism1.countNodes()
    nNodesOrg2 = organism2.countNodes()
    
    # Create the 2 childs and assign new IDs
    child1 = copy.deepcopy(organism1)
    child2 = copy.deepcopy(organism2)

    child1.setID(organismFactory.getID())
    child2.setID(organismFactory.getID())

    # Select random nodes from each child
    randomNode1 = random.randint(0,nNodesOrg1 - 1)
    randomNode2 = random.randint(0,nNodesOrg2 - 1)
    node1 = child1.getNode(randomNode1)
    node2 = child2.getNode(randomNode2)

    # Save the number of nodes taken from each  child
    nNodesFromOrg1 = node1.countNodes()
    nNodesFromOrg2 = node2.countNodes()

    # Search is done by ID, so this should be checked bc 
    # you could find duplicated IDs!!

    parentNode1 = child1.getParent(node1.ID)
    parentNode2 = child2.getParent(node2.ID)
    
    # Swap nodes
    # Set nodes in oposite children 
    if parentNode1["isRootNode"]:
        # Its the root node of child 1
        child1.setRootNode(node2)
    else:
        if parentNode1["isLeftSide"]:
            # Child on left side
            parentNode1["self"].setNode1(node2)
        else:
            # Child on right side
            parentNode1["self"].setNode2(node2)

    if parentNode2["isRootNode"]:
        # Its the root node of child 2
        child2.setRootNode(node1)
    else: 
        if parentNode2["isLeftSide"]:
            # Child on left side
            parentNode2["self"].setNode1(node1)
        else:
            # Child on right side
            parentNode2["self"].setNode2(node1)

    nNodesChild1 = child1.countNodes()
    nNodesChild2 = child2.countNodes()

    # Reset childs IDs
    child1.resetIDs()
    child2.resetIDs()

    
    # dictionary with an organism and similarities to each parent
    child1Similarities = {
            "simOrg1": (nNodesChild1 - nNodesFromOrg2) / nNodesChild1,
            "simOrg2": nNodesFromOrg2 / nNodesChild1,
            "child": child1
            }
    

    child2Similarities = {
            "simOrg1": nNodesFromOrg1 / nNodesChild2,
            "simOrg2": (nNodesChild2 - nNodesFromOrg1) / nNodesChild2,
            "child": child2
            }
    return {"child1": child1Similarities, "child2": child2Similarities}


# Reads configuration file and sets up all program variables
def setUp():

    # specify as global variable so it can be accesed in local contexts outside setUp

    global END_WHILE_METHOD
    global POPULATION_LENGTH
    global DATASET_BASE_PATH_DIR 
    global RESULT_BASE_PATH_DIR
    global POSITIVE_FILENAME 
    global NEGATIVE_FILENAME 
    global RESULT_PATH_PATH_DIR
    global MAX_SEQUENCES_TO_FIT
    global MIN_ITERATIONS
    global MIN_SCORE
    global THRESHOLD
    global MUTATION_PROBABILITY
    
    # Config data
    global configOrganism
    global configOrganismFactory
    global configConnector
    global configPssm

    config = readJsonFile(JSON_CONFIG_FILENAME)
    POPULATION_LENGTH = config["main"]["POPULATION_LENGTH"]
    DATASET_BASE_PATH_DIR = config["main"]["DATASET_BASE_PATH_DIR"]
    RESULT_BASE_PATH_DIR = config["main"]["RESULT_BASE_PATH_DIR"]
    POSITIVE_FILENAME = config["main"]["POSITIVE_FILENAME"]
    NEGATIVE_FILENAME = config["main"]["NEGATIVE_FILENAME"]
    MAX_SEQUENCES_TO_FIT = config["main"]["MAX_SEQUENCES_TO_FIT"]
    MIN_ITERATIONS = config["main"]["MIN_ITERATIONS"]
    MIN_SCORE = config["main"]["MIN_SCORE"]
    THRESHOLD = config["main"]["THRESHOLD"]
    MUTATION_PROBABILITY = config["main"]["MUTATION_PROBABILITY"]
    END_WHILE_METHOD = config["main"]["END_WHILE_METHOD"]


    configOrganism = config["organism"]
    configOrganismFactory = config["organismFactory"]
    configConnector = config["connector"]
    configPssm = config["pssm"]


# Reads a fasta file and returns an array of DNA sequences (strings)
def readFastaFile(filename):
    dataset = []

    fasta_sequences = SeqIO.parse(open(filename), 'fasta')    

    for fasta in fasta_sequences:
        dataset.append(str(fasta.seq))

    return dataset


# Reads a JSON file and returns a dictionary
def readJsonFile(filename):
    
    with open(filename) as json_content:

        return json.load(json_content)


# Entry point to app execution
# It calculates the time, but could include other app stats
if __name__ == '__main__':
    initial = time.time()
    setUp()
    main()


    print("\n")
    print("-"*50)
    print(" "*20+"STATS")
    print("-"*50)
    print("Time: {}".format(time.time()-initial))
    print("-"*50)



