# Main execution

from Bio import SeqIO
import time
from objects.organismFactory import OrganismFactory
import random
import copy

POPULATION_LENGTH = 10
DATASETS_BASE_DIR = 'datasets/'
POSITIVE_FILENAME = 'myFakeDataset1.fa'
NEGATIVE_FILENAME = 'myFakeDataset2.fa'

MAX_FIT_SEQUENCES = 50

organismPopulation = []

positiveDataset = []
negativeDataset = []

threshold = 5
pm = 0.2

def main():
    
    # TEST READ  DATASET
    # mySeq = readFastaFile(POSITIVE_FILENAME)
    # print(len(mySeq))
    positiveDataset = readFastaFile(POSITIVE_FILENAME)
    negativeDataset = readFastaFile(NEGATIVE_FILENAME)

    #print(positiveDataset)
    #print(negativeDataset)
    # Generate initial population
    organismFactory = OrganismFactory()
    for i in range(POPULATION_LENGTH):
        organismPopulation.append(organismFactory.getOrganism())
    print("Total population: "+str(len(organismPopulation)))
    maxScore = float("-inf")
    lastMaxScore = 0.0 
    # Main loop, it iterates untill organisms do not get a significant change.
    while abs(lastMaxScore - maxScore) > threshold:
        
        # Shuffle population
        random.shuffle(organismPopulation)
        random.shuffle(negativeDataset)
        random.shuffle(positiveDataset)
        lastMaxScore = maxScore    
        maxScore = float("-inf")
        # Iterate over pairs of organisms
        for i in range(0, len(organismPopulation) - 1, 2):
            org1 = organismPopulation[i]
            org2 = organismPopulation[i+1]
            print("Parents: {} {}".format(org1.ID, org2.ID))
            # Cross parents to get childs
            children = combineOrganisms(org1, org2, organismFactory)
            
            # Match with its closest child
            # There are 2 possible combinations p1-c1, p2-c2 & p1-c2, p2-c1
            # We select a combination based on a sum of similarities in combinations
            combination1 = children[0]["simOrg1"] + children[1]["simOrg2"] # Match the first parent to first child and second parent to second child 
            combination2 = children[0]["simOrg2"] + children[1]["simOrg1"] # Match the first parent to second child and second parent to first child
            pairChildren = []
            
            #TODO: Mutate children and parents with a probability pm

            if combination1 > combination2:
                pairChildren.append((org1, children[0]["child"]))
                pairChildren.append((org2, children[1]["child"]))
            else:
                pairChildren.append((org1, children[1]["child"]))
                pairChildren.append((org2, children[0]["child"]))
            for p, c in pairChildren:
                print("Fighting: {} {}".format(p.ID, c.ID))
            # "Fight" two organisms
            for j in range(len(pairChildren)):

                firstOrganism = pairChildren[j][0]
                secondOrganism = pairChildren[j][1]

                score1 = firstOrganism.getScore(positiveDataset[:MAX_FIT_SEQUENCES]) / firstOrganism.getScore(negativeDataset[:MAX_FIT_SEQUENCES])
                score2 = secondOrganism.getScore(positiveDataset[:MAX_FIT_SEQUENCES]) / secondOrganism.getScore(negativeDataset[:MAX_FIT_SEQUENCES])
                
                if(score1 > score2):
                    organismPopulation[i+j] = firstOrganism
                    print("Wins {} - {}".format(firstOrganism.ID, score1))
                    if score1 > maxScore:
                        maxScore = score1
                else:
                    organismPopulation[i+j] = secondOrganism
                    print("Wins {} - {}".format(secondOrganism.ID, score2))
                    if score2 > maxScore:
                        maxScore = score2
                #END FOR j
            
            # END FOR i

        # Show IDs of final array
        print("-"*10)
        for i in range(len(organismPopulation)):
            print(str(organismPopulation[i].ID))
        print("-"*10)
        # END WHILE

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

    # parentNodeX has an specific scructure:
    # parentNodeX[0] = isRootNode
    # parentNodeX[1] = parentNode
    # parentNodeX[2] = isLeftSideNode
    # 
    # Search is done by ID, so this should be checked bc 
    # you could find duplicated IDs!!

    parentNode1 = child1.getParent(node1.ID)
    parentNode2 = child2.getParent(node2.ID)
    
    # Swap nodes
    # Set nodes in oposite children 
    if parentNode1[0]:
        # Its the root node of child 1
        child1.setRootNode(node2)
    else:
        if parentNode1[2]:
            # Child on left side
            parentNode1[1].setNode1(node2)
        else:
            # Child on left side
            parentNode1[1].setNode2(node2)

    if parentNode2[0]:
        # Its the root node of child 2
        child2.setRootNode(node1)
    else: 
        if parentNode2[2]:
            # Child on left side
            parentNode2[1].setNode1(node1)
        else:
            # Child on left side
            parentNode2[1].setNode2(node1)

    nNodesChild1 = child1.countNodes()
    nNodesChild2 = child2.countNodes()
    
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
    return [child1Similarities, child2Similarities]


    
# Reads a fasta files and return an array of DNA sequences (strings)
def readFastaFile(filename):
    dataset = []

    fasta_sequences = SeqIO.parse(open(DATASETS_BASE_DIR + filename), 'fasta')    

    for fasta in fasta_sequences:
        dataset.append(str(fasta.seq))

    return dataset

# Entry point to app execution
# It calculates the time, but could include other app stats
if __name__ == '__main__':
    initial = time.time()
    main()
    print("\n")
    print("-"*50)
    print(" "*20+"STATS")
    print("-"*50)
    print("Time: {}".format(time.time()-initial))
    print("-"*50)



