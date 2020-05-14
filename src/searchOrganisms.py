# Main execution

from Bio import SeqIO
import time
from objects.organismFactory import OrganismFactory
import random
import copy
import json
import numpy as np
import os


POPULATION_LENGTH = 0
DATASET_BASE_PATH_DIR = ""
RESULT_BASE_PATH_DIR = ""
POSITIVE_FILENAME = ""
NEGATIVE_FILENAME = ""
POPULATION_ORIGIN = ""
POPULATION_FILL_TYPE = ""
INPUT_FILENAME = ""
OUTPUT_FILENAME = ""

JSON_CONFIG_FILENAME = "config.json"
configOrganism = {}
configOrganismFactory = {}
configConnector = {}
configPssm = {}

MAX_SEQUENCES_TO_FIT_POS = 0
MAX_SEQUENCES_TO_FIT_NEG = 0
MIN_ITERATIONS = 0
MIN_FITNESS = 0

RECOMBINATION_PROBABILITY = 0.0

organismPopulation = []
# meanNodes are the mean nodes of the population organisms lately
# used to calculate organism complexity
meanNodes = 0
meanFitness = 0


positiveDataset = []
negativeDataset = []

THRESHOLD = 0.0


def main():

    print("Loading parameters...")
    positiveDataset = readFastaFile(DATASET_BASE_PATH_DIR + POSITIVE_FILENAME)
    negativeDataset = readFastaFile(DATASET_BASE_PATH_DIR + NEGATIVE_FILENAME)

    meanNodes = 0
    meanFitness = 0
    # Generate initial population
    organismFactory = OrganismFactory(
        configOrganism, configOrganismFactory, configConnector, configPssm
    )
    organismPopulation = []

    if POPULATION_ORIGIN.lower() == "random":
        # For a random origin, we can generate #POPULATION_LENGTH organisms.

        for i in range(POPULATION_LENGTH):
            newOrganism = organismFactory.getOrganism()
            organismPopulation.append(newOrganism)
            meanNodes += newOrganism.countNodes()

    elif POPULATION_ORIGIN.lower() == "file":
        # Set the file organisms and fill with random/same organisms
        # POPULATION_LENGTH must be >= len(fileOrganisms)
        fileOrganisms = organismFactory.importOrganisms(INPUT_FILENAME)
        remainingOrganisms = POPULATION_LENGTH - len(fileOrganisms)
        fillOrganismPopulation = []

        if POPULATION_FILL_TYPE.lower() == "random":
            # FILL WITH RANDOM

            for i in range(remainingOrganisms):
                newOrganism = organismFactory.getOrganism()
                fillOrganismPopulation.append(newOrganism)

        elif POPULATION_FILL_TYPE.lower() == "uniform":
            # FILL WITH SAME

            for i in range(remainingOrganisms):
                newOrganism = copy.deepcopy(fileOrganisms[i % len(fileOrganisms)])
                fillOrganismPopulation.append(newOrganism)
                newOrganism.setID(organismFactory.getID())

        # join & calculate mean nodes
        organismPopulation = fileOrganisms + fillOrganismPopulation

        for org in organismPopulation:
            meanNodes += org.countNodes()

    else:
        print("Not a valid population origin, check the configuration file.")

        return -1

    meanNodes /= POPULATION_LENGTH
    # print("mean: {}".format(meanNodes))
    iterations = 0

    maxScore = float("-inf")
    lastMaxScore = 0.0
    bestOrganism = (None, 0.0, 0, 0.0)

    timeformat = "%Y-%m-%d--%H-%M-%S"
    print("Starting execution...")
    # Main loop, it iterates until organisms do not get a significant change
    # or MIN_ITERATIONS is reached.
    # it can be done untill you get a score over a value
    try:
        while not isFinished(END_WHILE_METHOD, iterations, maxScore, lastMaxScore):

            # Shuffle population & datasets
            random.shuffle(organismPopulation)
            random.shuffle(negativeDataset)
            random.shuffle(positiveDataset)

            # Reset maxScore
            lastMaxScore = maxScore
            maxScore = float("-inf")
            changedBestScore = False
            initial = time.time()

            aFitness = []
            aNodes = []

            # Iterate over pairs of organisms

            for i in range(0, len(organismPopulation) - 1, 2):
                org1 = organismPopulation[i]
                org2 = organismPopulation[i + 1]

                # Cross parents to get childs
                children = combineOrganisms(org1, org2, organismFactory)

                child1 = children["child1"]["child"]
                child2 = children["child2"]["child"]

                # Match with its closest child
                # There are 2 possible combinations p1-c1, p2-c2 & p1-c2, p2-c1
                # We select a combination based on a sum of similarities in combinations
                combination1 = (
                    children["child1"]["simOrg1"] + children["child2"]["simOrg2"]
                )  # Match the first parent to first child and second parent to second child
                combination2 = (
                    children["child1"]["simOrg2"] + children["child2"]["simOrg1"]
                )  # Match the first parent to second child and second parent to first child

                pairChildren = []

                # Mutate children

                child1.mutate(organismFactory)
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

                    firstOrganism = pairChildren[j][0]  # Parent Organism
                    secondOrganism = pairChildren[j][1]  # Chid Organism
                    p1 = firstOrganism.getScore(
                        positiveDataset[:MAX_SEQUENCES_TO_FIT_POS]
                    )
                    n1 = firstOrganism.getScore(
                        negativeDataset[:MAX_SEQUENCES_TO_FIT_NEG]
                    )
                    # Compute complexity after gettig the score
                    c1 = firstOrganism.getComplexity(meanNodes, meanFitness)

                    p2 = secondOrganism.getScore(
                        positiveDataset[:MAX_SEQUENCES_TO_FIT_POS]
                    )
                    n2 = secondOrganism.getScore(
                        negativeDataset[:MAX_SEQUENCES_TO_FIT_NEG]
                    )
                    # Compute complexity after gettig the score
                    c2 = secondOrganism.getComplexity(meanNodes, meanFitness)

                    # Check negative scores (horrible matching)
                    # if(p1 < 0 or n1 < 0):
                    #    print("Negative values on {}!!Pos: {} Neg:{}".format(firstOrganism.ID, p1, n1))

                    # if(p2 < 0 or n2 < 0):
                    #    print("Negative values on {}!!POS: {} Neg:{}".format(secondOrganism.ID, p2, n2))

                    fitness1 = p1 - n1
                    effectiveFitness1 = fitness1 - COMPLEXITY_FACTOR * c1

                    fitness2 = p2 - n2
                    effectiveFitness2 = fitness2 - COMPLEXITY_FACTOR * c2

                    # print("ID1: {} EFitness1:{:.2f}-{:.2f}-{:.2f} =  {:.2f} \nID2: {} EFitness2: {:.2f}-{:.2f}-{:.2f} = {:.2f}".format(firstOrganism.ID, p1, n1, c1, effectiveFitness1, secondOrganism.ID, p2, n2, c2, effectiveFitness2))

                    if effectiveFitness1 > effectiveFitness2:  # The first organism wins
                        # Set it back to the population and save fitness for next iteration
                        organismPopulation[i + j] = firstOrganism
                        aFitness.append(fitness1)
                        # If the parent wins, meanNodes don't change
                        aNodes.append(firstOrganism.countNodes())

                        # Check if its the max score in that iteration

                        if effectiveFitness1 > maxScore:
                            maxScore = effectiveFitness1
                            maxScoreP = p1
                            maxPenalty = c1
                            maxNodes = firstOrganism.countNodes()

                        # Check if its the max score in the program

                        if effectiveFitness1 > bestOrganism[1]:
                            # ID, EF, Nodes, Penalty applied
                            bestOrganism = (
                                firstOrganism,
                                effectiveFitness1,
                                firstOrganism.countNodes(),
                                c1,
                            )
                            changedBestScore = True

                    else:  # The second organism wins
                        # Set it back to the population and save fitness for next iteration
                        organismPopulation[i + j] = secondOrganism
                        aFitness.append(fitness2)
                        # If the child wins, update meanNodes
                        # meanNodes = ((meanNodes * POPULATION_LENGTH) + secondOrganism.countNodes() - firstOrganism.countNodes()) / POPULATION_LENGTH
                        aNodes.append(secondOrganism.countNodes())

                        # Set the tracking to the new child...
                        secondOrganism.setIsTracked(firstOrganism.isTracked)

                        if secondOrganism.isTracked:
                            # Export it If its being tracked
                            println(
                                "Evolution {}->{}".format(
                                    firstOrganism.ID, secondOrganism.ID
                                ),
                                RESULT_BASE_PATH_DIR + "evolution.txt",
                            )
                            filename = "tr{}_{}".format(
                                time.strftime(timeformat), secondOrganism.ID
                            )
                            exportOrganism(
                                secondOrganism,
                                positiveDataset,
                                filename,
                                organismFactory,
                            )

                        # Check if its the max score in that iteration

                        if effectiveFitness2 > maxScore:
                            maxScore = effectiveFitness2
                            maxScoreP = p2
                            maxPenalty = c2
                            maxNodes = secondOrganism.countNodes()

                        # Check if its the max score in the program

                        if effectiveFitness2 > bestOrganism[1]:
                            # ID, EF, Nodes, Penalty applied
                            bestOrganism = (
                                secondOrganism,
                                effectiveFitness2,
                                secondOrganism.countNodes(),
                                c2,
                            )
                            changedBestScore = True

                    # END FOR j

                # END FOR i

            # Compute mean fitness of the organisms
            meanFitness = np.mean(aFitness)
            meanNodes = np.mean(aNodes)

            # Show IDs of final array
            # print("-"*10)
            m, s = divmod((time.time() - initial), 60)
            h, m = divmod(m, 60)
            sTime = "{}h:{}m:{:.2f}s".format(int(h), int(m), s)
            println(
                "Iter: {} AN:{:.2f} AF:{:.2f} - MF: {:.2f} MSP: {:.2f} MP: {:.2f} MN: {} -  BO: {} BF: {:.2f} BN: {} BP: {:.2f} Time: {}".format(
                    iterations,
                    meanNodes,
                    meanFitness,
                    maxScore,
                    maxScoreP,
                    maxPenalty,
                    maxNodes,
                    bestOrganism[0].ID,
                    bestOrganism[1],
                    bestOrganism[2],
                    bestOrganism[3],
                    sTime,
                ),
                RESULT_BASE_PATH_DIR + OUTPUT_FILENAME,
            )

            # Print agains a random positive secuence
            random.shuffle(positiveDataset)
            print(bestOrganism[0].printResult(positiveDataset[0]))

            if changedBestScore or iterations % PERIODIC_EXPORT == 0:
                filename = "{}_{}".format(time.strftime(timeformat), bestOrganism[0].ID)
                exportOrganism(
                    bestOrganism[0], positiveDataset, filename, organismFactory
                )

            # print("-"*10)
            iterations += 1
            # END WHILE

    except Exception as e:
        print("Exited program: \n{}\n".format(e))

    finally:
        print()
        print("-" * 10)
        print("Best Organism {}: {}".format(bestOrganism[0].ID, bestOrganism[1]))


# Checks if main while loop is finished
# methods: 'Iterations', 'minScore', 'Threshold'
def isFinished(method, iterations, maxScore, lastMaxScore):

    if method.lower() == "iterations":
        return iterations >= MIN_ITERATIONS

    elif method.lower() == "fitness":
        return maxScore >= MIN_FITNESS

    elif method.lower() == "threshold":
        return abs(lastMaxScore - maxScore) <= THRESHOLD

    return True


def exportOrganism(organism, dataset, filename, factory):

    organismFile = "{}{}_organism.txt".format(RESULT_BASE_PATH_DIR, filename)
    organismFileJSON = "{}{}_organism.json".format(RESULT_BASE_PATH_DIR, filename)
    resultsFile = "{}{}_results.txt".format(RESULT_BASE_PATH_DIR, filename)

    organism.export(organismFile)
    organism.exportResults(dataset, resultsFile)
    factory.exportOrganisms([organism], organismFileJSON)


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

    # Combine parents with a probability p

    if random.random() < RECOMBINATION_PROBABILITY:

        # Select random nodes from each child
        randomNode1 = random.randint(0, nNodesOrg1 - 1)
        randomNode2 = random.randint(0, nNodesOrg2 - 1)
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

        # Reset fitness Scores
        child1.resetScores()
        child2.resetScores()

        # dictionary with an organism and similarities to each parent
        child1Similarities = {
            "simOrg1": (nNodesChild1 - nNodesFromOrg2) / nNodesChild1,
            "simOrg2": nNodesFromOrg2 / nNodesChild1,
            "child": child1,
        }

        child2Similarities = {
            "simOrg1": nNodesFromOrg1 / nNodesChild2,
            "simOrg2": (nNodesChild2 - nNodesFromOrg1) / nNodesChild2,
            "child": child2,
        }
    else:

        # Reset fitness Scores
        child1.resetScores()
        child2.resetScores()

        # If childs are not recombined, return the same organisms and their similarities
        child1Similarities = {
            "simOrg1": 1,  # Equal to organism 1
            "simOrg2": 0,
            "child": child1,
        }

        child2Similarities = {
            "simOrg1": 0,
            "simOrg2": 1,  # Equal to organism2
            "child": child2,
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
    global MAX_SEQUENCES_TO_FIT_POS
    global MAX_SEQUENCES_TO_FIT_NEG
    global MIN_ITERATIONS
    global MIN_FITNESS
    global THRESHOLD
    global COMPLEXITY_FACTOR
    global POPULATION_ORIGIN
    global POPULATION_FILL_TYPE
    global INPUT_FILENAME
    global OUTPUT_FILENAME
    global RECOMBINATION_PROBABILITY
    global PERIODIC_EXPORT

    # Config data
    global configOrganism
    global configOrganismFactory
    global configConnector
    global configPssm

    config = readJsonFile(JSON_CONFIG_FILENAME)
    POPULATION_LENGTH = config["main"]["POPULATION_LENGTH"]
    DATASET_BASE_PATH_DIR = config["main"]["DATASET_BASE_PATH_DIR"]
    RESULT_BASE_PATH_DIR = (
        config["main"]["RESULT_BASE_PATH_DIR"] + time.strftime("%Y%m%d%H%M%S") + "/"
    )
    POSITIVE_FILENAME = config["main"]["POSITIVE_FILENAME"]
    NEGATIVE_FILENAME = config["main"]["NEGATIVE_FILENAME"]
    MAX_SEQUENCES_TO_FIT_POS = config["main"]["MAX_SEQUENCES_TO_FIT_POS"]
    MAX_SEQUENCES_TO_FIT_NEG = config["main"]["MAX_SEQUENCES_TO_FIT_NEG"]
    MIN_ITERATIONS = config["main"]["MIN_ITERATIONS"]
    MIN_FITNESS = config["main"]["MIN_FITNESS"]
    THRESHOLD = config["main"]["THRESHOLD"]
    END_WHILE_METHOD = config["main"]["END_WHILE_METHOD"]
    COMPLEXITY_FACTOR = config["main"]["COMPLEXITY_FACTOR"]
    POPULATION_ORIGIN = config["main"]["POPULATION_ORIGIN"]
    POPULATION_FILL_TYPE = config["main"]["POPULATION_FILL_TYPE"]
    INPUT_FILENAME = config["main"]["INPUT_FILENAME"]
    OUTPUT_FILENAME = config["main"]["OUTPUT_FILENAME"]
    RECOMBINATION_PROBABILITY = config["main"]["RECOMBINATION_PROBABILITY"]
    PERIODIC_EXPORT = config["main"]["PERIODIC_EXPORT"]

    # Create directory where the output and results will be stored
    os.mkdir(RESULT_BASE_PATH_DIR)

    # Store Config into variables to use later
    configOrganism = config["organism"]
    configOrganismFactory = config["organismFactory"]
    configConnector = config["connector"]
    configPssm = config["pssm"]

    # Throw config on a file
    parametersPath = RESULT_BASE_PATH_DIR + "parameters.txt"
    println("-" * 50, parametersPath)
    println(" " * 20 + "PARAMETERS", parametersPath)
    println("-" * 50, parametersPath)
    printConfigJSON(config["main"], "Main Config", parametersPath)
    printConfigJSON(configOrganism, "Organism Config", parametersPath)
    printConfigJSON(configOrganismFactory, "Organism Factory Config", parametersPath)
    printConfigJSON(configConnector, "Connector Config", parametersPath)
    printConfigJSON(configPssm, "PSSM Config", parametersPath)
    println("-" * 50, parametersPath)


# Reads a fasta file and returns an array of DNA sequences (strings)
def readFastaFile(filename):
    dataset = []

    fasta_sequences = SeqIO.parse(open(filename), "fasta")

    for fasta in fasta_sequences:
        dataset.append(str(fasta.seq))

    return dataset


# Reads a JSON file and returns a dictionary
def readJsonFile(filename):

    with open(filename) as json_content:

        return json.load(json_content)


def printConfigJSON(config, name, path):
    println("{}:".format(name), path)

    for key in config.keys():
        println("{}: {}".format(key, config[key]), path)
    println("\n", path)


# Shows the string on stdout and write it to a file
def println(string, nameFile):

    print(string)

    # Here we are sure file exists
    f = open(nameFile, "a+")
    f.write(string + "\n")
    f.close()


# Entry point to app execution
# It calculates the time, but could include other app stats

if __name__ == "__main__":
    initial = time.time()
    setUp()
    main()

    m, s = divmod((time.time() - initial), 60)
    h, m = divmod(m, 60)
    println(
        "Time: {}h:{}m:{:.2f}s".format(int(h), int(m), s),
        RESULT_BASE_PATH_DIR + "parameters.txt",
    )
